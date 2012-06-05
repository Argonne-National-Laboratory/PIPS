/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_cgmip.c
 * @brief  Chvatal-Gomory cuts computed via a sub-MIP
 * @author Marc Pfetsch
 *
 * Separate Chv&aacute;tal-Gomory cuts using a sub-MIP. The approach is based on the following papers.
 *
 * M. Fischetti and A. Lodi@n
 * Optimizing over the first Chv&aacute;tal closure,@n
 * in: M. J&uuml;nger and V. Kaibel (eds.) Integer Programming and Combinatorial Optimization IPCO 2005,@n
 * LNCS 3509, pp. 12-22. Springer, Berlin Heidelberg New York (2005)
 *
 * P. Bonami, G. Cornu&eacute;jols, S. Dash, M. Fischetti, and A. Lodi@n
 * Projected Chv&aacute;tal-Gomory cuts for mixed integer linear programs,@n
 * Mathematical Programming 113, No. 2 (2008)
 *
 *
 * There are two versions to generate the final cut:
 *
 * - The CMIR-routines of SCIP can be used (if usecmir is true). One can determine which bound is
 *   used in the rounding operation (if cmirownbounds is true) or let SCIP choose the best. This
 *   version is generally numerically more stable than the second one.
 * - One can directly generate the CG-cut as computed (usecmir = false). The cut is not take from
 *   the solution of the MIP, but is recomputed, and some care (but not as much as in the first
 *   version) has been taken to create a valid cut.
 *
 * The computation time of the separation MIP is limited as follows:
 * - There is a node limit (parameter @a nodelimit).
 * - If paramter @a earlyterm is true, the separation is run until the first cut that is violated is
 *   found. (Note that these cuts are not necessarily added to the LP, because here also the norm of
 *   the cuts are taken into account - which cannot easily be included into the separation subscip.) 
 *   Then the solution is continued for a certain number of nodes.
 *
 * @todo Check whether one can weaken the conditions on the continuous variables.
 *
 * @warning This plugin is not yet fully tested.
 * @warning This separator should be used carefully - it may require a long separation time.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_cgmip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/pub_misc.h"


#define SEPA_NAME              "cgmip"
#define SEPA_DESC              "Chvatal-Gomory cuts via MIPs separator"
#define SEPA_PRIORITY             -1000
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP           TRUE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_MAXROUNDS             5 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        50 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXDEPTH             -1 /**< maximal depth at which the separator is applied */
#define DEFAULT_DYNAMICCUTS        TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_TIMELIMIT          1e20 /**< time limit for sub-MIP */
#define DEFAULT_MEMORYLIMIT        1e20 /**< memory limit for sub-MIP */
#define DEFAULT_NODELIMIT       10000LL /**< node limit for sub-MIP */
#define DEFAULT_OBJWEIGHT         1e-03 /**< objective weight for artificial variables */
#define DEFAULT_USECMIR            TRUE /**< use CMIR-generator (otherwise add cut directly) */
#define DEFAULT_CMIROWNBOUNDS     FALSE /**< tell CMIR-generator which bounds to used in rounding? */
#define DEFAULT_ALLOWLOCAL        FALSE /**< allow to generate local cuts */
#define DEFAULT_ONLYINTVARS       FALSE /**< generate cuts for problems with only integer variables? */
#define DEFAULT_ONLYACTIVEROWS    FALSE /**< use only active rows to generate cuts? */
#define DEFAULT_MAXROWAGE            -1 /**< maximal age of rows to consider if onlyactiverows is false */
#define DEFAULT_USECUTPOOL         TRUE /**< use cutpool to store CG-cuts even if the are not efficient? */
#define DEFAULT_PRIMALSEPARATION   TRUE /**< only separate cuts that are tight for the best feasible solution? */
#define DEFAULT_ONLYRANKONE       FALSE /**< whether only rank 1 inequalities should be separated */
#define DEFAULT_EARLYTERM          TRUE /**< terminate separation if a violated (but possibly sub-optimal) cut has been found? */
#define DEFAULT_ADDVIOLATIONCONS   TRUE /**< add constraint to subscip that only allows violated cuts? */
#define DEFAULT_ADDVIOLCONSHDLR   FALSE /**< add constraint handler to filter out violated cuts? */
#define DEFAULT_CONSHDLRUSENORM    TRUE /**< should the violation constraint handler use the norm of a cut to check for feasibility? */
#define DEFAULT_OBJLONE           FALSE /**< should the objective of the sub-MIP minimize the l1-norm of the multipliers? */
#define DEFAULT_CONTCONVERT        TRUE /**< convert some integral variables to be continuous to reduce the size of the sub-MIP */
#define DEFAULT_CONTCONVFRAC        0.5 /**< fraction of integral variables converted to be continuous (if contconvert) */
#define DEFAULT_CONTCONVMIN         100 /**< minimum number of integral variables before some are converted to be continuous */

#define NROWSTOOSMALL                 5 /**< only separate if the number of rows is larger than this number */
#define NCOLSTOOSMALL                 5 /**< only separate if the number of columns is larger than this number */

#define EPSILONVALUE              1e-03 /**< epsilon value needed to model strict-inequalities */
#define BETAEPSILONVALUE          1e-02 /**< epsilon value for fracbeta - is larger than EPSILONVALUE for numerical stability */
#define CUTCOEFBND               1000.0 /**< bounds on the values of the coefficients in the CG-cut */
#define STALLNODELIMIT           1000LL /**< number of stalling nodes if earlyterm is true */
#define MINEFFICACY                0.05 /**< minimum efficacy of a cut - compare set.c */

/* parameters used for CMIR-generation (taken from sepa_gomory) */
#define BOUNDSWITCH              0.9999
#define USEVBDS                    TRUE
#define MINFRAC                  0.0009 /* to allow a deviation of the same size as EPSILONVALUE */
#define MAXFRAC                  0.9991 /* to allow a deviation of the same size as EPSILONVALUE */
#define FIXINTEGRALRHS            FALSE
#define MAKECONTINTEGRAL          FALSE
#define MAXWEIGHTRANGE            1e+05 /**< maximal valid range max(|weights|)/min(|weights|) of row weights */

#if 0
#define MAXAGGRLEN(nvars)          (0.1*(nvars)+1000) /**< maximal length of base inequality */
#endif
#define MAXAGGRLEN(nvars)          nvars              /**< currently very large to allow any generation */


/** separator data */
struct SCIP_SepaData
{
   int                   maxrounds;          /**< maximal number of separation rounds per node (-1: unlimited) */
   int                   maxroundsroot;      /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxdepth;           /**< maximal depth at which the separator is applied */
   SCIP_Real             objweight;          /*<< objective weight for artificial variables */
   SCIP_Bool             dynamiccuts;        /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Real             timelimit;          /**< time limit for subscip */
   SCIP_Real             memorylimit;        /**< memory limit for subscip */
   SCIP_Longint          nodelimit;          /**< node limit for subscip */
   SCIP_Bool             usecmir;            /**< use CMIR-generator (otherwise add cut directly) */
   SCIP_Bool             cmirownbounds;      /**< tell CMIR-generator which bounds to used in rounding? */
   SCIP_Bool             allowlocal;         /**< allow local cuts */
   SCIP_Bool             onlyintvars;        /**< generate cuts for problems with only integer variables? */
   SCIP_Bool             onlyactiverows;     /**< use only active rows to generate cuts? */
   int                   maxrowage;          /**< maximal age of rows to consider if onlyactiverows is false */
   SCIP_Bool             usecutpool;         /**< use cutpool to store CG-cuts even if the are not efficient? */
   SCIP_Bool             primalseparation;   /**< only separate cuts that are tight for the best feasible solution? */
   SCIP_Bool             onlyrankone;        /**< whether only rank 1 inequalities should be separated */
   SCIP_Bool             earlyterm;          /**< terminate separation if a violated (but possibly sub-optimal) cut has been found? */
   SCIP_Bool             addviolationcons;   /**< add constraint to subscip that only allows violated cuts? */
   SCIP_Bool             addviolconshdlr;    /**< add constraint handler to filter out violated cuts? */
   SCIP_Bool             conshdlrusenorm;    /**< should the violation constraint handler use the cut-norm to check for feasibility? */
   SCIP_Bool             objlone;            /**< should the objective of the sub-MIP minimize the l1-norm of the multipliers? */
   SCIP_Bool             contconvert;        /**< convert some integral variables to be continuous to reduce the size of the sub-MIP */
   SCIP_Real             contconvfrac;       /**< fraction of integral variables converted to be continuous (if contconvert) */
   int                   contconvmin;        /**< minimum number of integral variables before some are converted to be continuous */
};


/** what happens for columns in the LP */
enum CGMIP_ColType
{
   colPresent    = 0,    /**< column is present in the separating MIP */
   colContinuous = 1,    /**< column corresponds to a continuous variable */
   colConverted  = 2,    /**< column is converted to be continuous */
   colAtUb       = 3,    /**< variable corresponding to column was at it's upper bound and was complemented */
   colAtLb       = 4     /**< variable corresponding to column was at it's lower bound (possibly complemented) */
};
typedef enum CGMIP_ColType CGMIP_COLTYPE;


/** data for the sub-MIP */
struct CGMIP_MIPData
{
   SCIP*                 subscip;            /**< pointer to (sub)SCIP data structure containing the auxiliary IP */
   unsigned int          m;                  /**< number of constraints of subscip */
   unsigned int          n;                  /**< number of variables of subscip */
   unsigned int          nrows;              /**< number of rows of original LP */
   unsigned int          ncols;              /**< number of columns of original LP */

   SCIP_VAR**            alpha;              /**< cut coefficient variable (NULL if not in separating MIP) */
   SCIP_VAR*             beta;               /**< rhs of cut */
   SCIP_VAR**            fracalpha;          /**< fractional part of lhs of cut (NULL if not present) */
   SCIP_VAR*             fracbeta;           /**< fractional part of rhs of cut */
   CGMIP_COLTYPE*        coltype;            /**< type for the columns */
   SCIP_Bool*            iscomplemented;     /**< whether the variable was complemented */
   SCIP_Bool*            isshifted;          /**< whether the variable was shifted to have 0 lower bound */

   SCIP_VAR**            ylhs;               /**< auxiliary row variables for lhs (NULL if not present) */
   SCIP_VAR**            yrhs;               /**< auxiliary row variables for rhs (NULL if not present) */

   SCIP_VAR**            z;                  /**< auxiliary variables for upper bounds (NULL if not present) */

   SCIP_Bool             conshdlrusenorm;    /**< copy from sepadata */
};
typedef struct CGMIP_MIPData CGMIP_MIPDATA;



/*
 * constraint handler to filter out violated cuts
 */

/* constraint handler properties */
#define CONSHDLR_NAME          "violatedCuts"
#define CONSHDLR_DESC          "only allow solutions corresponding to violated cuts"

/** constraint handler data */
struct SCIP_ConshdlrData
{
   CGMIP_MIPDATA*        mipdata;            /**< data of separating sub-MIP */
};


/** check whether cut corresponding to solution is violated */
static
SCIP_Bool solCutIsViolated(
   CGMIP_MIPDATA*        mipdata,            /**< data of separating sub-MIP */
   SCIP_SOL*             sol                 /**< solution to be checked */
   )
{
   SCIP* subscip;
   SCIP_Real act;
   SCIP_Real norm;
   SCIP_Real val;
   SCIP_VAR* var;
   SCIP_Real rhs;
   unsigned int j;

   assert( mipdata != NULL );
   subscip = mipdata->subscip;
   assert( subscip != NULL );

   /* initialize activity and norm */
   act = 0.0;
   norm = 1.0;

   /* compute activity and norm  */
   if ( mipdata->conshdlrusenorm )
   {
      SCIP_Real cutsqrnorm = 0.0;
      for (j = 0; j < mipdata->ncols; ++j)
      {
         var = mipdata->alpha[j];
         if ( var == NULL )
            continue;

         val = SCIPgetSolVal(subscip, sol, var);
         if ( !SCIPisZero(subscip, val) )
         {
            act += SCIPvarGetObj(var) * val;
            cutsqrnorm += SQR(val);
         }
      }
      norm = SQRT(cutsqrnorm);

      /* if norm is 0, the cut is trivial */
      if ( SCIPisZero(subscip, norm) )
         return FALSE;
   }
   else
   {
      for (j = 0; j < mipdata->ncols; ++j)
      {
         var = mipdata->alpha[j];
         if ( var == NULL )
            continue;

         val = SCIPgetSolVal(subscip, sol, var);
         if ( !SCIPisZero(subscip, val) )
            act += SCIPvarGetObj(var) * val;
      }
   }

   /* get rhs */
   rhs = SCIPgetSolVal(subscip, sol, mipdata->beta);

#ifdef SCIP_DEBUG
   if ( SCIPisEfficacious(subscip, (act - rhs)/norm) )
   {
      SCIPdebugMessage("Violated cut from solution - act: %f, rhs: %f, norm: %f, eff.: %f\n", act, rhs, norm, (act-rhs)/norm);
   }
#endif

   return SCIPisEfficacious(subscip, (act - rhs)/norm);   
}


/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeViolatedCuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpViolatedCuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( result != NULL );

   assert( SCIPgetNLPBranchCands(scip) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   violated = solCutIsViolated(conshdlrdata->mipdata, NULL);

   if ( violated )
      *result = SCIP_FEASIBLE;
   else
      *result = SCIP_CUTOFF;  /* cutoff, since all integer variables are integer, but the solution is not feasible */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsViolatedCuts)
{  /*lint --e{715}*/
   assert( result != NULL );

   /* this function should better not be called, since we need an LP solution for the sub-MIP to
      make sense, because of the multiplier variables. We therefore return SCIP_FEASIBLE. */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckViolatedCuts)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool violated;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( sol != NULL );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   violated = solCutIsViolated(conshdlrdata->mipdata, sol);

   if ( violated )
      *result = SCIP_FEASIBLE;
   else
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockViolatedCuts)
{  /*lint --e{715}*/
   /* do not lock variables */
   return SCIP_OKAY;
}


/** creates the violated CG-cut constraint handler and includes it in SCIP */
static
SCIP_RETCODE SCIPincludeConshdlrViolatedCut(
   SCIP*                 scip,               /**< SCIP data structure */
   CGMIP_MIPDATA*        mipdata             /**< data of separating sub-MIP */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->mipdata = mipdata;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         -1000000, -1000000, -1000000, -1, -1, 100, 0, FALSE, FALSE, FALSE, FALSE,
         SCIP_PROPTIMING_BEFORELP,
         NULL, consFreeViolatedCuts, NULL, NULL,
         NULL, NULL, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL,
         consEnfolpViolatedCuts, consEnfopsViolatedCuts, consCheckViolatedCuts, 
         NULL, NULL, NULL, consLockViolatedCuts,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
         conshdlrdata) );

   return SCIP_OKAY;
}





/*
 * local methods
 */


/** stores nonzero elements of dense coefficient vector as sparse vector, and calculates activity and norm
 *
 *  copied from sepa_gomory.c
 */
static
SCIP_RETCODE storeCutInArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars,               /**< problem variables */
   SCIP_Real*            cutcoefs,           /**< dense coefficient vector */
   SCIP_Real*            varsolvals,         /**< dense variable LP solution vector */
   char                  normtype,           /**< type of norm to use for efficacy norm calculation */
   SCIP_VAR**            cutvars,            /**< array to store variables of sparse cut vector */
   SCIP_Real*            cutvals,            /**< array to store coefficients of sparse cut vector */
   int*                  cutlen,             /**< pointer to store number of nonzero entries in cut */
   SCIP_Real*            cutact,             /**< pointer to store activity of cut */
   SCIP_Real*            cutnorm             /**< pointer to store norm of cut vector */
   )
{
   SCIP_Real val;
   SCIP_Real cutsqrnorm;
   SCIP_Real act;
   SCIP_Real norm;
   int len;
   int v;

   assert( nvars == 0 || cutcoefs != NULL );
   assert( nvars == 0 || varsolvals != NULL );
   assert( cutvars != NULL );
   assert( cutvals != NULL );
   assert( cutlen != NULL );
   assert( cutact != NULL );
   assert( cutnorm != NULL );

   len = 0;
   act = 0.0;
   norm = 0.0;
   switch ( normtype )
   {
   case 'e':
      cutsqrnorm = 0.0;
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutsqrnorm += SQR(val);
            cutvars[len] = vars[v];
            cutvals[len++] = val;
         }
      }
      norm = SQRT(cutsqrnorm);
      break;
   case 'm':
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            if ( REALABS(val) > norm )
               norm = REALABS(val);
            cutvars[len] = vars[v];
            cutvals[len++] = val;
         }
      }
      break;
   case 's':
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            norm += REALABS(val);
            cutvars[len] = vars[v];
            cutvals[len++] = val;
         }
      }
      break;
   case 'd':
      for (v = 0; v < nvars; ++v)
      {
         val = cutcoefs[v];
         if ( !SCIPisZero(scip, val) )
         {
            act += val * varsolvals[v];
            cutvars[len] = vars[v];
            cutvals[len++] = val;
         }
      }
      if ( len > 0 )
         norm = 1.0;
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", normtype);
      return SCIP_INVALIDDATA;
   }

   *cutlen = len;
   *cutact = act;
   *cutnorm = norm;

   return SCIP_OKAY;
}


/** Compute lhs/rhs for transformed column
 *
 *  Consider a variable \f$x_j\f$ and some row of the original system:
 *  \f[
 *       \gamma \leq a^T x \leq \delta, \quad \ell_j \leq x_j \leq u_j.
 *  \f]
 *  We perform the transformation
 *  \f[
 *       x_i' = \left\{
 *       \begin{array}{ll}
 *         s + \frac{1}{\sigma}\, x_j & \mbox{if }i = j\\
 *         x_i              & \mbox{otherwise},
 *       \end{array}
 *       \right.
 *  \f]
 *  where \f$s\f$ is the offset value and \f$\sigma\f$ is a scaling factor. The new system is
 *  \f[
 *     \gamma + \sigma\, a_j\,s \leq \sum_{i \neq j} a_i\, x_i' + \sigma a_j\, x_j' \leq \delta + \sigma\, a_j\, s
 *  \f]
 *  with bounds
 *  \f[
 *     \frac{1}{\sigma} \ell_j + s \leq x_j' \leq \frac{1}{\sigma} u_j + s, \qquad \mbox{ if }\sigma > 0
 *  \f]
 *  and
 *  \f[
 *     \frac{1}{\sigma} u_j + s \leq x_j' \leq \frac{1}{\sigma} \ell_j + s, \qquad \mbox{ if }\sigma < 0.
 *  \f]
 *
 *  This can be used as follows:
 *
 *  - If \f$x_j \geq \ell_j\f$ has a (nonzero) lower bound, one can use \f$s = -\ell_j\f$, \f$\sigma = 1\f$,
 *    and obtain \f$\gamma - a_j\,\ell_j \leq a^T x' \leq \delta - a_j\,\ell_j\f$, \f$0 \leq x_j' \leq u_j - \ell_j\f$.
 *
 *  - If \f$x_j \leq u_j\f$ has a (nonzero) upper bound, one can use \f$s = u_j\f$, \f$\sigma = -1\f$,
 *    and obtain \f$\gamma - a_j\,u_j \leq \sum_{i \neq j} a_i\, x_i' - a_j\, x_j' \leq \delta - a_j\, u_j\f$,
 *    \f$0 \leq x_j' \leq u_j - \ell_j\f$.
 */
static
SCIP_RETCODE transformColumn(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_COL*             col,                /**< column that should be complemented */
   SCIP_Real             offset,             /**< offset by which column should be shifted */
   SCIP_Real             sigma,              /**< scaling factor */
   SCIP_Real*            lhs,                /**< array of lhs of rows */
   SCIP_Real*            rhs,                /**< array rhs of rows */
   SCIP_Real*            lb,                 /**< pointer to lb of column */
   SCIP_Real*            ub,                 /**< pointer to ub of column */
   SCIP_Real*            primsol             /**< pointer to solution value */
   )
{
   SCIP_ROW** colrows;
   SCIP_Real* colvals;
   int pos, i;

   assert( scip != NULL );
   assert( lhs != NULL );
   assert( rhs != NULL );
   assert( col != NULL );

   colrows = SCIPcolGetRows(col);
   colvals = SCIPcolGetVals(col);
   assert( SCIPcolGetNLPNonz(col) == 0 || colrows != NULL );
   assert( SCIPcolGetNLPNonz(col) == 0 || colvals != NULL );
   assert( ! SCIPisZero(scip, sigma) );

   /* loop through rows that contain column */
   for (i = 0; i < SCIPcolGetNLPNonz(col); ++i)
   {
      SCIP_ROW* row;

      row = colrows[i];
      assert( row != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
         continue;

      pos = SCIProwGetLPPos(row);
      assert( 0 <= pos && pos < (int) mipdata->nrows );

      assert( ! SCIPisInfinity(scip, lhs[pos]) );
      if ( ! SCIPisInfinity(scip, -lhs[pos]) )
         lhs[pos] += sigma * colvals[i] * offset;

      assert( ! SCIPisInfinity(scip, -rhs[pos]) );
      if ( ! SCIPisInfinity(scip, rhs[pos]) )
         rhs[pos] += sigma * colvals[i] * offset;
   }

   /* correct lower and upper bounds and solution */
   if ( SCIPisNegative(scip, sigma) )
   {
      SCIP_Real l;

      assert( ! SCIPisInfinity(scip, -*ub) );
      if ( ! SCIPisInfinity(scip, *ub) )
         l = *ub/sigma + offset;
      else
         l = -SCIPinfinity(scip);

      assert( ! SCIPisInfinity(scip, *lb) );
      if ( ! SCIPisInfinity(scip, -*lb) )
         *ub = *lb/sigma + offset;
      else
         *ub = SCIPinfinity(scip);
      *lb = l;
   }
   else
   {
      assert( ! SCIPisInfinity(scip, *lb) );
      if ( ! SCIPisInfinity(scip, -*lb) )
         *lb = *lb/sigma + offset;
      assert( ! SCIPisInfinity(scip, -*ub) );
      if ( ! SCIPisInfinity(scip, *ub) )
         *ub = *ub/sigma + offset;
   }
   *primsol = *primsol/sigma + offset;

   return SCIP_OKAY;
}


/** Creates a subscip representing the separating MIP.
 *
 *  Let the constraints of the original MIP be of the following form:
 *  \f[
 *    \begin{array}{l@{\;}ll}
 *      a \leq A x + & C r & \leq b\\
 *      \ell \leq x & & \leq u\\
 *      c \leq & r & \leq d\\
 *      x \in Z^n.
 *    \end{array}
 *  \f]
 *  Here, some of the bounds may have value \f$\infty\f$ or \f$-\infty\f$.  Written in
 *  \f$\leq\f$-form this becomes:
 *  \f[
 *    \begin{array}{r@{\;}l}
 *      \tilde{A} x + \tilde{C} r & \leq \tilde{b}\\
 *      -x & \leq -\ell\\
 *      x & \leq u\\
 *      -r & \leq -c\\
 *      r & \leq d\\
 *      x \in Z^n,
 *    \end{array}
 *  \f]
 *  where we use
 *  \f[
 *    \tilde{A} =
 *    \left[
 *    \begin{array}{r}
 *      -A \\
 *      A 
 *    \end{array}
 *    \right],
 *    \quad
 *    \tilde{C} =
 *    \left[
 *    \begin{array}{r}
 *      - C\\
 *      C
 *    \end{array}
 *    \right]
 *    \qquad\mbox{ and }\qquad
 *    \tilde{b} =
 *    \left[
 *    \begin{array}{r}
 *      -a\\
 *      b
 *    \end{array}
 *    \right].
 *  \f]
 *  For the moment we assume that \f$c = 0\f$, i.e., the lower bounds on the continuous variables
 *  are 0.  To obtain a Chv&aacute;tal-Gomory cut we have to find nonnegative multipliers \f$y\f$,
 *  \f$\underline{z}\f$, and \f$\overline{z}\f$ such that
 *  \f[
 *      y^T \tilde{A} - \underline{z}^T + \overline{z}^T  \in Z \qquad\mbox{ and }\qquad
 *      y^T \tilde{C} \geq 0.
 *  \f]
 *  Note that we use zero multipliers for the bounds on the continuous variables \f$r\f$. Moreover,
 *  if some bounds are infinity, the corresponding multipliers are assumed to be 0. From these
 *  conditions, we obtain
 *  \f[
 *      (y^T \tilde{A} - \underline{z}^T + \overline{z}^T)\, x +
 *      y^T \tilde{C} \, r \leq
 *      y^T \tilde{b} - \underline{z}^T \ell + \overline{z}^T u.
 *  \f]
 *  Because \f$r \geq 0\f$, we can ignore the term \f$y^T \tilde{C} \, r \geq 0\f$ and obtain the
 *  following cut:
 *  \f[
 *      (y^T \tilde{A} - \underline{z}^T + \overline{z}^T )\, x \leq
 *      \lfloor y^T \tilde{b} - \underline{z}^T \ell + \overline{z}^T u \rfloor.
 *  \f]
 *  Assume that \f$\ell = 0\f$ for the meantime. Then the cut can be written as:
 *  \f[
 *      \lfloor y^T \tilde{A} + \overline{z}^T \rfloor \, x \leq
 *      \lfloor y^T \tilde{b} + \overline{z}^T u \rfloor.
 *  \f]
 *
 *  Following Fischetti and Lodi [2005], let \f$(x^*,r^*)\f$ be a fractional solution of the above
 *  original system.  The separating MIP created below is
 *  \f[
 *    \begin{array}{rlr@{\;}l}
 *       \max & \multicolumn{2}{@{}l}{(x^*)^T \alpha - \beta - w^T y} &\\
 *            & f = & \tilde{A}^T y + \overline{z} - \alpha & \\
 *            & \tilde{f} = & \tilde{b}^T y + u^T \overline{z} - \beta &\\
 *            & & \tilde{C}^T y & \geq 0\\
 *            & & 0 \leq f & \leq 1 - \epsilon \\
 *            & & 0 \leq \tilde{f} & \leq 1 - \epsilon\\
 *            & & 0 \leq y, \overline{z} & \leq 1 - \epsilon.\\
 *            & & \alpha \in Z^m, \beta & \in Z.
 *    \end{array}
 *  \f]
 *  Here, \f$w\f$ is a weight vector; it's idea is to make the sum over all components of \f$y\f$ as
 *  small as possible, in order to generate sparse cuts.
 *
 *  We perform the following additional computations:
 *
 *  - If the lower bounds on \f$x_i\f$ or \f$r_j\f$ are finite, we shift the variable to have a zero
 *    lower bound, i.e., we replace it by \f$x_i - \ell_i\f$ (or \f$r_j - u_j\f$). This is helpful in
 *    several ways: As seen above, the resulting inequalities/formulations simplify. Moreover, it
 *    allows to drop a variable if \f$x^*_i = 0\f$, see the next comment. If the lower bounds are not
 *    finite, but the upper bounds are finite, we can complement the variable. If the variables are
 *    free, the above formulation changes as follows: For free continuous variables, we require
 *    \f$\tilde{C}^T y = 0\f$. For a free integer variable \f$x_j\f$ (which rarely occurs in
 *    practice), we require \f$f_j = 0\f$, i.e., we force that \f$(\tilde{A}^T y + \overline{z})_j =
 *    \alpha_j\f$.
 *
 *  - If \f$x^*_j = 0 = \ell_j\f$ (after the above preprocessing), we drop variable \f$\alpha_j\f$
 *    from the formulation. Let \f$(\alpha^*, \beta^*, y^*, \overline{z}^*)\f$ be an
 *    optimal solution to the separating MIP. Then we can compute \f$\alpha_j =
 *    \lfloor(\tilde{A}_j^T y^* + \overline{z}^*)\rfloor\f$.
 *
 *  - If \f$x^*_i = u_i\f$, we complement the variable and drop it from the formulation, since the
 *    lower bound is 0 afterwards.
 *
 *  - If a variable has been shifted or complemented, we have to recompute \f$\beta\f$ with the
 *    original lhs/rhs.
 *
 *  - If a continuous variable \f$r_j\f$ is free, we have to force equality for the corresponding components in
 *    \f$y^T \tilde{C} \, r \geq 0\f$.
 *
 *  - If an integer variable \f$x_i\f$ is free, we are not allowed to round the cut down. In this
 *    case, the combintation of rows and bounds has to be integral. We force this by requiring that
 *    \f$f_i = 0\f$.
 *
 *  - If @p contconvert is true some integral variables are randomly treated as if they were
 *    continuous. This has the effect that in the resulting cut the corresponding coefficient has
 *    value 0. This makes the cuts more sparse. Moreover, the separation problems should become
 *    easier.
 *
 *  - If required, i.e., parameter @p primalseparation is true, we force a primal separation step. For
 *    this we require that the cut is tight at the currently best solution. To get reliable solutions
 *    we relax equality by EPSILONVALUE.
 */
static
SCIP_RETCODE createSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata             /**< data for sub-MIP */
   )
{
   SCIP* subscip;
   SCIP_COL** cols;
   SCIP_ROW** rows;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* primsol;

   int ncols;
   int nrows;
   int i, j;
   unsigned int cnt, ucnt;
   unsigned int nshifted;
   unsigned int ncomplemented;
   unsigned int nconverted;
   unsigned int nlbounds;
   unsigned int nubounds;

   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   SCIP_CONS* cons;
   int nconsvars;
   char name[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( sepadata != NULL );

   assert( mipdata->subscip == NULL );

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert( ncols > 0 && nrows > 0 );

   mipdata->m = 0;
   mipdata->n = 0;
   mipdata->nrows = (unsigned int) nrows;
   mipdata->ncols = (unsigned int) ncols;

   /* copy value */
   mipdata->conshdlrusenorm = sepadata->conshdlrusenorm;

   /* create subscip */
   SCIP_CALL( SCIPcreate( &(mipdata->subscip) ) );
   subscip = mipdata->subscip;
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* add violation constraint handler if requested */
   if ( sepadata->addviolconshdlr )
   {
      SCIP_CALL( SCIPincludeConshdlrViolatedCut(subscip, mipdata) );
   }

   SCIP_CALL( SCIPcreateProb(subscip, "sepa_cgmip separating MIP", NULL, NULL , NULL , NULL , NULL , NULL , NULL) );
   SCIP_CALL( SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE) );

   /* alloc memory for subscipdata elements */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->alpha), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->fracalpha), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->coltype), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->iscomplemented), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->isshifted), ncols) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->ylhs), nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->yrhs), nrows) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(mipdata->z), 2*ncols) );

   /* get temporary storage */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &primsol, ncols) );

   /* store lhs/rhs for complementing (see below) */
   for (i = 0; i < nrows; ++i)
   {
      SCIP_Real val;
      assert( rows[i] != NULL );

      val = SCIProwGetLhs(rows[i]) - SCIProwGetConstant(rows[i]);
      if ( SCIProwIsIntegral(rows[i]) )
         val = SCIPfeasCeil(scip, val); /* row is integral: round left hand side up */
      lhs[i] = val;

      val = SCIProwGetRhs(rows[i]) - SCIProwGetConstant(rows[i]);
      if ( SCIProwIsIntegral(rows[i]) )
         val = SCIPfeasFloor(scip, val); /* row is integral: round right hand side down */
      rhs[i] = val;
   }

   /* store lb/ub for complementing and perform preprocessing */
   nshifted = 0;
   ncomplemented = 0;
   nconverted = 0;
   nlbounds = 0;
   nubounds = 0;
   for (j = 0; j < ncols; ++j)
   {
      SCIP_COL* col;
      SCIP_VAR* var;

      col = cols[j];
      assert( col != NULL );
      var = SCIPcolGetVar(col);
      assert( var != NULL );

      primsol[j] = SCIPcolGetPrimsol(col);
      assert(SCIPisEQ(scip, SCIPgetVarSol(scip, var), primsol[j]) );

      lb[j] = SCIPvarGetLbGlobal(var);
      assert( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPcolGetLb(col)) );

      /* if allowed, try to use stronger local bound */
      if ( sepadata->allowlocal && SCIPisGT(scip, SCIPvarGetLbLocal(var), lb[j]) )
         lb[j] = SCIPvarGetLbLocal(var);

      ub[j] = SCIPvarGetUbGlobal(var);
      assert( SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPcolGetUb(col)) );

      /* if allowed, try to use stronger local bound */
      if ( sepadata->allowlocal && SCIPisLT(scip, SCIPvarGetUbLocal(var), ub[j]) )
         ub[j] = SCIPvarGetUbLocal(var);

      mipdata->coltype[j] = colPresent;
      mipdata->iscomplemented[j] = FALSE;
      mipdata->isshifted[j] = FALSE;

      /* check status of column/variable */
      if ( SCIPcolIsIntegral(col) )
      {
         /* possibly convert integral variables to be continuous */
         if ( sepadata->contconvert && ncols >= sepadata->contconvmin )
         {
            /* randomly convert variables */
            if ( ((SCIP_Real) rand())/((SCIP_Real) RAND_MAX) >= sepadata->contconvfrac )
            {
               /* preprocessing is also performed for converted columns */
               mipdata->coltype[j] = colConverted;
               ++nconverted;
            }
         }
      }
      else
      {
         /* detect continuous variables, but perform preprocessing for them */
         mipdata->coltype[j] = colContinuous;
      }

      /* if integer variable is at its upper bound -> complementing (this also generates a 0 lower bound) */
      if ( mipdata->coltype[j] == colPresent && SCIPisFeasEQ(scip, primsol[j], ub[j]) )
      {
         assert( ! SCIPisInfinity(scip, ub[j]) );
         SCIP_CALL( transformColumn(scip, sepadata, mipdata, col, ub[j], -1.0, lhs, rhs, &(lb[j]), &(ub[j]), &(primsol[j])) );
         mipdata->iscomplemented[j] = TRUE;
         mipdata->coltype[j] = colAtUb;
         ++nubounds;
      }
      else
      {
         /* if a variable has a finite nonzero lower bound -> shift */
         if ( ! SCIPisInfinity(scip, -lb[j]) )
         {
            if ( ! SCIPisZero(scip, lb[j]) )
            {
               SCIP_CALL( transformColumn(scip, sepadata, mipdata, col, -lb[j], 1.0, lhs, rhs, &(lb[j]), &(ub[j]), &(primsol[j])) );
               assert( SCIPisZero(scip, lb[j]) );
               mipdata->isshifted[j] = TRUE;
               ++nshifted;
            }

            /* if integer variable is at its lower bound */
            if ( mipdata->coltype[j] == colPresent && SCIPisZero(scip, primsol[j]) )
            {
               mipdata->coltype[j] = colAtLb;
               ++nlbounds;
            }
         }
         else
         {
            /* lower bound is minus-infinity -> check whether upper bound is finite */
            if ( ! SCIPisInfinity(scip, ub[j]) )
            {
               /* complement variable */
               SCIP_CALL( transformColumn(scip, sepadata, mipdata, col, ub[j], -1.0, lhs, rhs, &(lb[j]), &(ub[j]), &(primsol[j])) );
               assert( SCIPisZero(scip, lb[j]) );
               mipdata->iscomplemented[j] = TRUE;
               ++ncomplemented;

               /* if integer variable is at its lower bound */
               if ( mipdata->coltype[j] == colPresent && SCIPisZero(scip, primsol[j]) )
               {
                  mipdata->coltype[j] = colAtLb;
                  ++nlbounds;
               }
            }
         }
      }

      assert( SCIPisFeasLE(scip, lb[j], primsol[j]) );
      assert( SCIPisFeasLE(scip, primsol[j], ub[j]) );
   }

#ifndef NDEBUG
   if ( sepadata->contconvert && ncols >= sepadata->contconvmin )
   {
      SCIPdebugMessage("Converted %u integral variables to be continuous.\n", nconverted);
   }
#endif

   /* create artificial variables for row combinations (y-variables) */
   cnt = 0;
   for (i = 0; i < nrows; ++i)
   {
      SCIP_ROW* row;

      row = rows[i];
      assert( row != NULL );

      mipdata->ylhs[i] = NULL;
      mipdata->yrhs[i] = NULL;

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
         continue;

      /* skip rows that not have been active for a longer time */
      if ( ! sepadata->onlyactiverows && sepadata->maxrowage > 0 && SCIProwGetAge(row) > sepadata->maxrowage )
         continue;

      /* check whether we want to skip cut produced by the CGMIP separator */
      if ( sepadata->onlyrankone )
      {
         const char* rowname;
         rowname = SCIProwGetName(row);
         if ( strlen(rowname) > 5 )
         {
            /* check for name "cgcut..." */
            if ( rowname[0] == 'c' && rowname[1] == 'g' && rowname[2] == 'c' && rowname[3] == 'u' && rowname[4] == 't' )
               continue;
         }
      }

      /* if we have an equation */
      if ( SCIPisEQ(scip, lhs[i], rhs[i]) )
      {
         assert( ! SCIPisInfinity(scip, rhs[i]) );

         if ( ! sepadata->onlyactiverows || SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row)) )
         {
            /* create two variables for each equation */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yeq1_%d", i);
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->ylhs[i]), name, 0.0, 1.0-EPSILONVALUE,
                  -sepadata->objweight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, mipdata->ylhs[i]) );
            ++cnt;

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yeq2_%d", i);
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->yrhs[i]), name, 0.0, 1.0-EPSILONVALUE,
                  -sepadata->objweight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, mipdata->yrhs[i]) );
            ++cnt;
         }
      }
      else
      {
         /* create variable for lhs of row if necessary */
         if ( ! SCIPisInfinity(scip, -lhs[i]) && 
            ( ! sepadata->onlyactiverows || SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetLhs(row))) )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ylhs_%d", i);
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->ylhs[i]), name, 0.0, 1.0-EPSILONVALUE,
                  -sepadata->objweight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, mipdata->ylhs[i]) );
            ++cnt;
         }

         /* create variable for rhs of row if necessary */
         if ( ! SCIPisInfinity(scip, rhs[i]) && 
            ( ! sepadata->onlyactiverows || SCIPisFeasEQ(scip, SCIPgetRowLPActivity(scip, row), SCIProwGetRhs(row))) )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "yrhs_%d", i);
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->yrhs[i]), name, 0.0, 1.0-EPSILONVALUE,
                  -sepadata->objweight, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, mipdata->yrhs[i]) );
            ++cnt;
         }
      }
   }
   assert( (int) cnt <= 2 * nrows );
   mipdata->n += cnt;

   /* create alpha, bound, and fractional variables */
   cnt = 0;
   ucnt = 0;
   for (j = 0; j < ncols; ++j)
   {
      mipdata->z[j] = NULL;
      mipdata->alpha[j] = NULL;
      mipdata->fracalpha[j] = NULL;

      if ( mipdata->coltype[j] == colPresent )
      {
         SCIP_Real obj;

         if ( sepadata->objlone )
            obj = 0.0;
         else
            obj = primsol[j];

         /* create alpha variables */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "alpha_%d", j);
         SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->alpha[j]), name, -CUTCOEFBND, CUTCOEFBND, obj,
               SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(subscip, mipdata->alpha[j]) );
         ++cnt;

         /* create fractional variables */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "f_%d", j);
         if ( SCIPisInfinity(scip, -lb[j]) && SCIPisInfinity(scip, ub[j]) )
         {
            /* fix fractional value to be zero for free original variables */
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->fracalpha[j]), name, 0.0, 0.0, 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         }
         else
         {
            /* fractional value in [0, 1) for variables with finite bounds */
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->fracalpha[j]), name, 0.0, 1.0-EPSILONVALUE, 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         }
         SCIP_CALL( SCIPaddVar(subscip, mipdata->fracalpha[j]) );
         ++cnt;

         /* create variables for upper bounds */
         if ( ! SCIPisInfinity(scip, ub[j]) )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "zub_%d", j);
            SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->z[j]), name, 0.0, 1.0-EPSILONVALUE,
                  0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(subscip, mipdata->z[j]) );
            ++ucnt;
         }
      }
   }
   assert( (int) cnt <= 2 * ncols );
   assert( (int) ucnt <= ncols );

   /* create variable for the rhs of the cut */
   if ( sepadata->objlone )
   {
      SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->beta), "beta", -CUTCOEFBND, CUTCOEFBND, 0.0,
            SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   }
   else
   {
      SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->beta), "beta", -CUTCOEFBND, CUTCOEFBND, -1.0,
            SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   }
   SCIP_CALL( SCIPaddVar(subscip, mipdata->beta) );

   /* create fractional variable for the rhs */
   SCIP_CALL( SCIPcreateVar(subscip, &(mipdata->fracbeta), "fracbeta", 0.0, 1.0-BETAEPSILONVALUE, 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(subscip, mipdata->fracbeta) );
   mipdata->n += cnt + ucnt + 2;

   /* get temporary storage */
   SCIP_CALL( SCIPallocBufferArray(scip, &consvals, (int) mipdata->n) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, (int) mipdata->n) );

   /* create constraints for alpha variables of CG-cut */
   cnt = 0;
   for (j = 0; j < ncols; ++j)
   {
      SCIP_ROW** colrows;
      SCIP_Real* colvals;

      /* create ordinary part for all selected variables */
      if ( mipdata->coltype[j] == colPresent )
      {
         SCIP_Real sigma;

         assert( cols[j] != NULL );
         colrows = SCIPcolGetRows(cols[j]);
         colvals = SCIPcolGetVals(cols[j]);
         nconsvars = 0;

         if ( mipdata->iscomplemented[j] )
            sigma = -1.0;
         else
            sigma = 1.0;

         /* add part for columns */
         for (i = 0; i < SCIPcolGetNLPNonz(cols[j]); ++i)
         {
            SCIP_ROW* row;
            int pos;

            row = colrows[i];
            assert( row != NULL );

            /* skip modifiable rows and local rows, unless allowed */
            if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
               continue;

            pos = SCIProwGetLPPos(row);
            assert( 0 <= pos && pos < nrows );

            if ( mipdata->ylhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->ylhs[pos];
               consvals[nconsvars] = -sigma * colvals[i];
               ++nconsvars;
            }
            if ( mipdata->yrhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->yrhs[pos];
               consvals[nconsvars] = sigma * colvals[i];
               ++nconsvars;
            }
            assert( nconsvars <= (int) mipdata->n );
         }
         /* add part for upper bounds */
         if ( mipdata->z[j] != NULL )
         {
            assert( ! SCIPisInfinity(scip, ub[j]) );
            consvars[nconsvars] = mipdata->z[j];
            consvals[nconsvars] = 1.0;
            ++nconsvars;
         }
         assert( nconsvars <= (int) mipdata->n );

         /* add alpha variable */
         consvars[nconsvars] = mipdata->alpha[j];
         consvals[nconsvars] = -1.0;
         ++nconsvars;
         assert( nconsvars <= (int) mipdata->n );

         /* add fractional-alpha variable */
         consvars[nconsvars] = mipdata->fracalpha[j];
         consvals[nconsvars] = -1.0;
         ++nconsvars;
         assert( nconsvars <= (int) mipdata->n );

         /* add linear constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "alpha_%d", j);
         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, nconsvars, consvars, consvals, 0.0, 0.0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         ++cnt;
      }
      /* generate part that makes sure that cut is valid for continuous variables */
      else if ( mipdata->coltype[j] == colContinuous || mipdata->coltype[j] == colConverted )
      {
         SCIP_Real sigma;
         SCIP_Real r;

         assert( cols[j] != NULL );
         colrows = SCIPcolGetRows(cols[j]);
         colvals = SCIPcolGetVals(cols[j]);
         nconsvars = 0;

         if ( mipdata->iscomplemented[j] )
            sigma = -1.0;
         else
            sigma = 1.0;

         /* add part for columns */
         for (i = 0; i < SCIPcolGetNLPNonz(cols[j]); ++i)
         {
            SCIP_ROW* row;
            int pos;

            row = colrows[i];
            assert( row != NULL );

            /* skip modifiable rows and local rows, unless allowed */
            if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
               continue;

            pos = SCIProwGetLPPos(row);
            assert( 0 <= pos && pos < nrows );

            if ( mipdata->ylhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->ylhs[pos];
               consvals[nconsvars] = -sigma * colvals[i];
               ++nconsvars;
            }
            if ( mipdata->yrhs[pos] != NULL )
            {
               consvars[nconsvars] = mipdata->yrhs[pos];
               consvals[nconsvars] = sigma * colvals[i];
               ++nconsvars;
            }
            assert( nconsvars <= (int) mipdata->n );
         }

         /* add linear constraint */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cont_%d", j);

         /* for free continuous variables require equality */
         r = SCIPinfinity(subscip);
         if ( SCIPisInfinity(scip, -lb[j]) && SCIPisInfinity(scip, ub[j]) )
            r = 0.0;
         else
            assert( SCIPisZero(scip, lb[j]) );

         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, nconsvars, consvars, consvals, 0.0, r,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         ++cnt;
      }
   }
   assert( (int) cnt <= ncols );
   mipdata->m += cnt;

   /* create constraints for rhs of cut */
   nconsvars = 0;

   /* first for the rows */
   for (i = 0; i < nrows; ++i)
   {
      assert( rows[i] != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(rows[i]) || (SCIProwIsLocal(rows[i]) && !sepadata->allowlocal) )
         continue;

      /* if lhs is there */
      if ( mipdata->ylhs[i] != NULL && ! SCIPisZero(scip, lhs[i]) )
      {
         assert( ! SCIPisInfinity(scip, -lhs[i]) );
         consvars[nconsvars] = mipdata->ylhs[i];
         consvals[nconsvars] = -lhs[i];
         ++nconsvars;
      }
      /* if rhs is there */
      if ( mipdata->yrhs[i] != NULL && ! SCIPisZero(scip, rhs[i]) )
      {
         assert( ! SCIPisInfinity(scip, rhs[i]) );
         consvars[nconsvars] = mipdata->yrhs[i];
         consvals[nconsvars] = rhs[i];
         ++nconsvars;
      }
      assert( nconsvars <= (int) mipdata->n );
   }
   /* next for the columns */
   for (j = 0; j < ncols; ++j)
   {
      /* if ub is there */
      if ( mipdata->z[j] != NULL && ! SCIPisZero(scip, ub[j]) )
      {
         assert( mipdata->coltype[j] == colPresent );
         assert( ! SCIPisInfinity(scip, ub[j]) );
         consvars[nconsvars] = mipdata->z[j];
         consvals[nconsvars] = ub[j];
         ++nconsvars;
         assert( nconsvars <= (int) mipdata->n );
      }
   }
   /* add beta variable */
   consvars[nconsvars] = mipdata->beta;
   consvals[nconsvars] = -1.0;
   ++nconsvars;

   /* add fractional-beta variable */
   consvars[nconsvars] = mipdata->fracbeta;
   consvals[nconsvars] = -1.0;
   ++nconsvars;
   assert( nconsvars <= (int) mipdata->n );

   /* add linear constraint */
   SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "beta", nconsvars, consvars, consvals, 0.0, 0.0,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(subscip, cons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   ++mipdata->m;

   /* add primal separation constraint if required */
   if ( sepadata->primalseparation )
   {
      SCIP_SOL* bestsol;
      bestsol = SCIPgetBestSol(scip);
      if ( bestsol != NULL )
      {
         nconsvars = 0;
         for (j = 0; j < ncols; ++j)
         {
            if ( mipdata->alpha[j] != NULL )
            {
               SCIP_Real val;
               assert( mipdata->coltype[j] == colPresent );
               
               val = SCIPgetSolVal(scip, bestsol, SCIPcolGetVar(cols[j]));
               consvars[nconsvars] = mipdata->alpha[j];
               consvals[nconsvars] = val;
               ++nconsvars;
               assert( nconsvars <= (int) mipdata->n );
            }
         }
         consvars[nconsvars] = mipdata->beta;
         consvals[nconsvars] = -1.0;
         ++nconsvars;

         /* add linear constraint - allow slight deviation from equality */
         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "primalseparation", nconsvars, consvars, consvals, -EPSILONVALUE, EPSILONVALUE,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         ++mipdata->m;
      }
   }

   /* add constraint to force violated cuts if required */
   if ( sepadata->addviolationcons )
   {
      nconsvars = 0;
      for (j = 0; j < ncols; ++j)
      {
         if ( mipdata->alpha[j] != NULL )
         {
            consvars[nconsvars] = mipdata->alpha[j];
            consvals[nconsvars] = primsol[j];
            ++nconsvars;
            assert( nconsvars <= (int) mipdata->n );
         }
      }
      consvars[nconsvars] = mipdata->beta;
      consvals[nconsvars] = -1.0;
      ++nconsvars;

      /* add linear constraint - allow slight deviation from equality */
      SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, "violationConstraint", nconsvars, consvars, consvals, MINEFFICACY, SCIPinfinity(subscip),
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
      ++mipdata->m;
   }

   SCIPdebugMessage("subscip has %u variables and %u constraints (%u shifted, %u complemented, %u at lb, %u at ub).\n",
      mipdata->n, mipdata->m, nshifted, ncomplemented, nlbounds, nubounds);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &consvals);

   SCIPfreeBufferArray(scip, &primsol);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);

#ifdef SCIP_OUTPUT
   /* SCIPdebug( SCIP_CALL( SCIPprintOrigProblem(subscip, NULL, NULL, FALSE) ) ); */
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "debug.lp", "lp", FALSE) );
#endif

   return SCIP_OKAY;
}


/** sets parameters for subscip */
static
SCIP_RETCODE subscipSetParams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_Bool*            success             /**< if setting was successful -> stop solution otherwise */
   )
{
   SCIP* subscip;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( success != NULL );

   *success = TRUE;

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   /* set other limits of subscip */
   /* SCIP_CALL( SCIPsetObjlimit(subscip, mipdata->objectivelimit) ); */
   /* SCIP_CALL( SCIPsetIntParam(subscip, "limits/solutions", sepadata->sollimit) ); */

   /* do not abort subscip on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* determine output to console */
#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 4) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1000) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/nsols/active", 2) );
#else
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 1000) );
#endif

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/crossover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/mutation/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/dins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/undercover/freq", -1) );

   /* set other heuristics */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/shifting/freq", 3) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/simplerounding/freq", 1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rounding/freq", 1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/oneopt/freq", 1) );

   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/pscostdiving/freq", 1) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/feaspump/freq", 3) ); */

   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/coefdiving/freq", -1) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/fracdiving/freq", -1) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/guideddiving/freq", -1) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/linesearchdiving/freq", -1) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/objpscostdiving/freq", -1) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/rootsoldiving/freq", -1) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/veclendiving/freq", -1) ); */

   /* disable cut separation in subscip */
   SCIP_CALL( SCIPsetIntParam(subscip, "separating/cgmip/freq", -1) );

   /* disable expensive presolving */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "presolving/probing/maxrounds", 0) ); */
   /*     SCIP_CALL( SCIPsetIntParam(subscip, "constraints/linear/maxpresolpairrounds", 0) ); */
   /*     SCIP_CALL( SCIPsetRealParam(subscip, "constraints/linear/maxaggrnormscale", 0.0) ); */

   /* disable conflict analysis */
   /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) ); */
   /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) ); */
   /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) ); */
   /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) ); */
   /*     SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) ); */

   return SCIP_OKAY;
}


/** solve subscip */
static
SCIP_RETCODE solveSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_Bool*            success             /**< if setting was successful -> stop solution otherwise */
   )
{
   SCIP* subscip;
   SCIP_STATUS status;
   SCIP_Real timelimit;
   SCIP_Real memorylimit;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( success != NULL );

   *success = TRUE;

   subscip = mipdata->subscip;

   /* determine timelimit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( ! SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   if ( sepadata->timelimit < timelimit )
      timelimit = sepadata->timelimit;
   if ( timelimit > 0.0 )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   }
   else
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* determine memorylimit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if ( sepadata->memorylimit < memorylimit )
      memorylimit = sepadata->memorylimit;
   if ( ! SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if ( memorylimit > 0.0 )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   }
   else
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   /* set nodelimit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", sepadata->nodelimit) );

   /* check whether we want a complete solve */
   if ( ! sepadata->earlyterm )
   {
      SCIP_CALL( SCIPsolve(subscip) );
      status = SCIPgetStatus(subscip);

#ifdef SCIP_OUTPUT
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      /* if the solution process was terminated or the problem is infeasible (can happen because of violation constraint) */
      if ( status == SCIP_STATUS_TIMELIMIT || status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_INFORUNBD )
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      /* all other statuses except optimal are invalid */
      if ( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_NODELIMIT )
      {
         SCIPerrorMessage("Solution of subscip for CG-separation returned with invalid status %d.\n", status);
         return SCIP_ERROR;
      }
   }
   else
   {
      /* otherwise we want a heuristic solve */

      /* -> solve until first solution is found */
      SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", 1) );
      SCIP_CALL( SCIPsolve(subscip) );
      SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", -1) );

      status = SCIPgetStatus(subscip);

      /* if the solution process was terminated or the problem is infeasible (can happen because of violation constraint) */
      if ( status == SCIP_STATUS_TIMELIMIT || status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_NODELIMIT || 
         status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_INFORUNBD)
      {
         *success = FALSE;
         return SCIP_OKAY;
      }

      /* all other statuses except optimal or bestsollimit are invalid - (problem cannot be infeasible) */
      if ( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_BESTSOLLIMIT )
      {
         SCIPerrorMessage("Solution of subscip for CG-separation returned with invalid status %d.\n", status);
         return SCIP_ERROR;
      }

      /* solve some more, if a feasible solution was found */
      if ( status == SCIP_STATUS_BESTSOLLIMIT )
      {
         SCIPdebugMessage("Continue solving separation problem ...\n");

         SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", STALLNODELIMIT) );
         SCIP_CALL( SCIPsolve(subscip) );
         SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", -1LL) );

         status = SCIPgetStatus(subscip);
         assert( status != SCIP_STATUS_BESTSOLLIMIT );

#ifdef SCIP_OUTPUT
         SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

         /* if the solution process was terminated */
         if ( status == SCIP_STATUS_TIMELIMIT || status == SCIP_STATUS_USERINTERRUPT )
         {
            *success = FALSE;
            return SCIP_OKAY;
         }

         /* all other statuses except optimal or bestsollimit are invalid */
         if ( status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_STALLNODELIMIT && status != SCIP_STATUS_NODELIMIT )
         {
            SCIPerrorMessage("Solution of subscip for CG-separation returned with invalid status %d.\n", status);
            return SCIP_ERROR;
         }
      }
   }

   return SCIP_OKAY;
}


/** Computes cut from the given multipliers 
 *
 *  Note that the cut computed here in general will not be the same as the one computed with the
 *  sub-MIP, because of numerical differences. Here, we only combine rows whose corresponding
 *  multiplier is positive w.r.t. the feasibility tolerance. In the sub-MIP, however, the rows are
 *  combined in any case. This makes a difference, if the coefficients in the matrix are large and
 *  hence yield a value that is larger than the tolerance.
 *
 *  Because of the transformations we have the following:
 * 
 *  If variable \f$x_j\f$ was complemented, we have \f$x'_j = u_j - x_j\f$. If in the transformed
 *  system the lower bound is used, its corresponding multiplier is \f$y^T A'_j - \lfloor y^T A'_j
 *  \rfloor\f$, which corresponds to 
 *  \f[
 *      y^T A'_j - \lfloor y^T A'_j \rfloor = - y^T A_j - \lfloor - y^T A_j \rfloor = - y^T A_j + \lceil y^T A_j \rceil
 *  \f]
 *  in the original system.
 *
 *  If such a variable was at its upper bound before the transformation, it is at its lower bound
 *  afterwards. Hence, its contribution to the cut is 0.
 *
 *  Note that if the original LP-solution does not satisfy some of the rows with equality the
 *  violation of the cut might be smaller than what is computed with the reduced sub-MIP.
 *
 *  Furthermore, note that if continuous variables have been shifted, the computed violated may be
 *  different as well, because the necessary changes in the lhs/rhs are not used here anymore.
 *
 *  @todo check if cut is correct if continuous variables have been shifted.
 */
static
SCIP_RETCODE computeCut(
   SCIP*                 scip,               /**< original scip */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   SCIP_SOL*             sol,                /**< current solution for sub-MIP */
   SCIP_Real*            cutcoefs,           /**< coefficients of the cut */
   SCIP_Real*            cutrhs,             /**< rhs of the cut */
   SCIP_Bool*            localrowsused,      /**< pointer to store whether local rows were used in summation */
   SCIP_Bool*            localboundsused,    /**< pointer to store whether local bounds were used in summation */
   SCIP_Bool*            success             /**< whether we produced a valid cut */
   )
{
   SCIP* subscip;
   int i, j;
   int nvars;
   int ncols;
   SCIP_VAR** vars;
   SCIP_ROW** rows;
   SCIP_COL** cols;
   SCIP_Real val;
   SCIP_Real maxabsweight;
   int nrows;

   assert( scip != NULL );
   assert( mipdata != NULL );
   assert( sepadata != NULL );
   assert( sol != NULL );
   assert( cutcoefs != NULL );
   assert( cutrhs != NULL );
   assert( localrowsused != NULL );
   assert( success != NULL );

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   /* get data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   assert( nrows == (int) mipdata->nrows );
   assert( ncols == (int) mipdata->ncols );

   /* initialize */
   *success = TRUE;
   *localrowsused = FALSE;
   *localboundsused = FALSE;
   BMSclearMemoryArray(cutcoefs, nvars);
   *cutrhs = 0.0;

   /* find maximal absolute weight */
   maxabsweight = 0.0;
   for (i = 0; i < nrows; ++i)
   {
#ifndef NDEBUG
      const char* rowname;
#endif
      SCIP_ROW* row;
      SCIP_Real weight;

      row = rows[i];
      assert( row != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
      {
         assert( mipdata->ylhs[i] == NULL && mipdata->yrhs[i] == NULL );
         continue;
      }

#ifndef NDEBUG
      rowname = SCIProwGetName(row);
#endif

      /* get weight from solution */
      weight = 0.0;
      if ( mipdata->ylhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[i]);
         if ( SCIPisFeasPositive(scip, val) )
            weight = -val;

         assert( ! sepadata->onlyrankone || strlen(rowname) <= 5 || 
            rowname[0] != 'c' || rowname[1] != 'g' || rowname[2] != 'c' || rowname[3] != 'u' || rowname[4] != 't' );
      }
      if ( mipdata->yrhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[i]);

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, ABS(weight)) )
            weight = val;

         assert( ! sepadata->onlyrankone || strlen(rowname) <= 5 || 
            rowname[0] != 'c' || rowname[1] != 'g' || rowname[2] != 'c' || rowname[3] != 'u' || rowname[4] != 't' );
      }

      weight = REALABS(weight);
      if ( weight > maxabsweight )
         maxabsweight = weight;
   }

   /* calculate the row summation */
   for (i = 0; i < nrows; ++i)
   {
      SCIP_ROW* row;
      SCIP_Real weight;
      SCIP_Real absweight;
      SCIP_Bool uselhs;

      row = rows[i];
      assert( row != NULL );

      /* skip modifiable rows and local rows, unless allowed */
      if ( SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !sepadata->allowlocal) )
      {
         assert( mipdata->ylhs[i] == NULL && mipdata->yrhs[i] == NULL );
         continue;
      }

      /* get weight from solution */
      weight = 0.0;
      uselhs = FALSE;
      if ( mipdata->ylhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[i]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         if ( SCIPisFeasPositive(scip, val) )
         {
            uselhs = TRUE;
            weight = -val;
         }
      }
      if ( mipdata->yrhs[i] != NULL )
      {
         val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[i]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
         if ( SCIPisFeasGT(scip, val, ABS(weight)) )
            weight = val;
      }

      /* add row if weight is nonzero and lies within range */
      absweight = REALABS(weight);
      if ( ! SCIPisSumZero(scip, weight) && absweight * MAXWEIGHTRANGE >= maxabsweight )
      {
         SCIP_COL** rowcols;
         SCIP_Real* rowvals;

         rowcols = SCIProwGetCols(row);
         rowvals = SCIProwGetVals(row);

         /* add the row coefficients to the sum */
         for (j = 0; j < SCIProwGetNLPNonz(row); ++j)
         {
            int idx;
            SCIP_VAR* var;

            assert( rowcols[j] != NULL );
            var = SCIPcolGetVar(rowcols[j]);

            assert( var != NULL );
            assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );
            assert( SCIPvarGetCol(var) == rowcols[j] );

            idx = SCIPvarGetProbindex(var);
            assert( 0 <= idx && idx < nvars );

            cutcoefs[idx] += weight * rowvals[j];
         }

         /* compute rhs */
         if ( uselhs )
         {
            assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );
            val = SCIProwGetLhs(row) - SCIProwGetConstant(row);
            if ( SCIProwIsIntegral(row) )
               val = SCIPfeasCeil(scip, val); /* row is integral: round left hand side up */
         }
         else
         {
            assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );
            val = SCIProwGetRhs(row) - SCIProwGetConstant(row);
            if ( SCIProwIsIntegral(row) )
               val = SCIPfeasFloor(scip, val); /* row is integral: round right hand side down */
         }
         (*cutrhs) += weight * val;

         *localrowsused = *localrowsused || SCIProwIsLocal(row);
      }
   }

   /* add upper bounds */
   for (j = 0; j < ncols; ++j)
   {
      assert( cols[j] != NULL );
      if ( mipdata->z[j] != NULL )
      {
         assert( mipdata->coltype[j] == colPresent );

         val = SCIPgetSolVal(subscip, sol, mipdata->z[j]);
         assert( ! SCIPisFeasNegative(subscip, val) );

         /* if a bound has been used */
         if ( SCIPisSumPositive(subscip, val) )
         {
            SCIP_VAR* var;
            int idx;

            var = SCIPcolGetVar(cols[j]);

            assert( var != NULL );
            assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );
            assert( SCIPvarIsIntegral(var) );
            assert( SCIPvarGetCol(var) == cols[j] );

            idx = SCIPvarGetProbindex(var);
            assert( 0 <= idx && idx < nvars );

            /* check whether variable is complemented */
            if ( mipdata->iscomplemented[j] )
            {
               SCIP_Real lbnd;
               lbnd = SCIPvarGetLbGlobal(var);
               assert( ! SCIPisInfinity(scip, -lbnd) );
               assert( SCIPisIntegral(scip, lbnd) );
               assert( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPcolGetLb(cols[j])) );

               /* variable should not be free */
               assert( ! SCIPisInfinity(scip, -lbnd) || ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetLbLocal(var) - 0.5 > lbnd )
               {
                  lbnd = SCIPvarGetLbLocal(var);
                  assert( SCIPisIntegral(scip, lbnd) );
                  *localboundsused = TRUE;
               }

               cutcoefs[idx] -= val;
               *cutrhs -= lbnd * val;
            }
            else
            {
               SCIP_Real ubnd;
               ubnd = SCIPvarGetUbGlobal(var);
               assert( ! SCIPisInfinity(scip, ubnd) );
               assert( SCIPisIntegral(scip, ubnd) );
               assert( SCIPisEQ(scip, SCIPvarGetUbLocal(var), SCIPcolGetUb(cols[j])) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetUbLocal(var) + 0.5 < ubnd )
               {
                  ubnd = SCIPvarGetUbLocal(var);
                  assert( SCIPisIntegral(scip, ubnd) );
                  *localboundsused = TRUE;
               }

               cutcoefs[idx] += val;
               *cutrhs += ubnd * val;
            }
         }
      }
   }

   /* check lower bounds for integral variables */
   for (j = 0; j < nvars; ++j)
   {
      SCIP_VAR* var;
      int pos;

      var = vars[j];
      assert( var != NULL );
      pos = SCIPcolGetLPPos(SCIPvarGetCol(var));

      /* a variable may have status COLUMN, but the corresponding column may not (yet) be in the LP */
      if ( pos >= 0 && mipdata->coltype[pos] != colContinuous && mipdata->coltype[pos] != colConverted )
      {
         assert( pos < ncols );
         assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );
         assert( SCIPvarIsIntegral(var) );

         /* check whether variable is complemented */
         if ( mipdata->iscomplemented[pos] )
         {
            assert( ! mipdata->isshifted[pos] );
            /* if the variable is complemented, the multiplier for the upper bound arises from the
               lower bound multiplier for the transformed problem - because of the minus-sign in the
               transformation this yields a round-up operation. */
            val = SCIPfeasCeil(scip, cutcoefs[j]) - cutcoefs[j];
            assert( ! SCIPisFeasNegative(scip, val) );

            /* only if variable needs to be rounded */
            if ( SCIPisSumPositive(scip, val) )
            {
               SCIP_Real ubnd;
               ubnd = SCIPvarGetUbGlobal(var);
               assert( ! SCIPisInfinity(scip, ubnd) );
               assert( SCIPisIntegral(scip, ubnd) );

               /* variable should not be free */
               assert( ! SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) || ! SCIPisInfinity(scip, ubnd) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetUbLocal(var) + 0.5 < ubnd )
               {
                  ubnd = SCIPvarGetUbLocal(var);
                  assert( SCIPisIntegral(scip, ubnd) );
                  *localboundsused = TRUE;
               }

               /* round cut coefficients, i.e., add val to cutcoefs[j] */
               cutcoefs[j] = SCIPfeasCeil(scip, cutcoefs[j]);

               /* correct rhs */
               if ( ! SCIPisSumZero(scip, ubnd) )
                  *cutrhs += ubnd * val;
            }
         }
         else
         {
            /* compute multiplier for lower bound: */
            val = cutcoefs[j] - SCIPfeasFloor(scip, cutcoefs[j]);
            assert( ! SCIPisFeasNegative(scip, val) );

            /* only if variable needs to be rounded */
            if ( SCIPisSumPositive(scip, val) )
            {
               SCIP_Real lbnd;
               lbnd = SCIPvarGetLbGlobal(var);
               assert( ! SCIPisInfinity(scip, -lbnd) );
               assert( SCIPisIntegral(scip, lbnd) );

               /* variable should not be free */
               assert( ! SCIPisInfinity(scip, -lbnd) || ! SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) );

               /* if allowed, try to use stronger local bound */
               if ( sepadata->allowlocal && SCIPvarGetLbLocal(var) - 0.5 > lbnd )
               {
                  lbnd = SCIPvarGetLbLocal(var);
                  assert( SCIPisIntegral(scip, lbnd) );
                  *localboundsused = TRUE;
               }

               /* round cut coefficients, i.e., subtract val from cutcoefs[j] */
               cutcoefs[j] = SCIPfeasFloor(scip, cutcoefs[j]);

               /* correct rhs */
               if ( ! SCIPisSumZero(scip, lbnd) )
                  *cutrhs -= lbnd * val;
            }
         }
      }
      else
      {
         /* force coefficients of all continuous variables or of variables not in the lp to zero */
         assert( pos == -1 || mipdata->coltype[pos] == colContinuous || mipdata->coltype[pos] == colConverted );

         /* check whether all coefficients for continuous or converted variables are nonnegative */
         if ( pos >= 0 )
         {
            if ( SCIPisNegative(scip, cutcoefs[j]) )
            {
               *success = FALSE;
               break;
            }
         }

         cutcoefs[j] = 0.0;
      }
   }

   /* round rhs */
   *cutrhs = SCIPfeasFloor(scip, *cutrhs);

   return SCIP_OKAY;
}


/** Create CG-cuts directly from solution of sub-MIP */
static
SCIP_RETCODE createCGCutsDirect(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   unsigned int*         ngen                /**< number of generated cuts */
   )
{
   SCIP* subscip;
   SCIP_STAGE stage;
   SCIP_SOL** sols;
   SCIP_ROW** prevrows;
   int nprevrows;
   int nsols;
   int k, s;
   char name[SCIP_MAXSTRLEN];

   SCIP_VAR** vars;
   SCIP_Real* cutcoefs;
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_Real* varsolvals;
   SCIP_Bool cutislocal;
   SCIP_Bool success;
   char normtype;
   int nvars;

   SCIP_Real cutrhs;
   SCIP_Real cutact;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( ngen != NULL );

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   SCIPdebugMessage("Trying to generate cuts directly ...\n");
   *ngen = 0;

   /* check if solving was successful and get solutions */
   stage = SCIPgetStage(subscip);
   if ( stage == SCIP_STAGE_SOLVING || stage == SCIP_STAGE_SOLVED )
      nsols = SCIPgetNSols(subscip);
   else
      nsols = 0;

   /* only if solutions have been found */
   if ( nsols == 0 )
      return SCIP_OKAY;

   sols = SCIPgetSols(subscip);

   /* get the type of norm to use for efficacy calculations */
   SCIP_CALL( SCIPgetCharParam(scip, "separating/efficacynorm", &normtype) );

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prevrows, nsols) );

   /* get solution values */
   for (k = 0; k < nvars; ++k)
   {
      if ( SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_COLUMN )
         varsolvals[k] = SCIPvarGetLPSol(vars[k]);
      else
         varsolvals[k] = 0.0;
   }

   /* loop through solutions found */
   nprevrows = 0;
   for (s = 0; s < nsols; ++s)
   {
      SCIP_SOL* sol;
      SCIP_Bool localrowsused;
      SCIP_Bool localboundsused;
      sol = sols[s];

      /* compute coefficients */
      SCIP_CALL( computeCut(scip, mipdata, sepadata, sol, cutcoefs, &cutrhs, &localrowsused, &localboundsused, &success) );
      cutislocal = localrowsused || localboundsused;

      /* take next solution if cut was not valid */
      if ( ! success )
      {
         SCIPdebugMessage("cut not valid - skipping ...\n");
         continue;
      }

      /* compute activity */
      cutact = 0.0;
      for (k = 0; k < nvars; ++k)
         cutact += cutcoefs[k] * varsolvals[k];

      /* the following test should be treated with care because of numerical differences - see computeCut() */ 
#if 0
      {
         /* check for correctness of computed values */
         SCIP_Real obj = 0.0;
         SCIP_Real val;
         SCIP_Bool contVarShifted = FALSE;
         unsigned int j;
         SCIP_COL** cols;
         int ncols;

         SCIP_CALL( SCIPprintSol(subscip, sol, NULL, FALSE) );

         SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
         for (j = 0; j < mipdata->ncols; ++j)
         {
            if ( mipdata->coltype[j] == colPresent )
            {
               int idx;
               assert( mipdata->alpha[j] != NULL );
               val = SCIPgetSolVal(subscip, sol, mipdata->alpha[j]);
               assert( SCIPisFeasIntegral(subscip, val) );
               idx = SCIPvarGetProbindex(SCIPcolGetVar(cols[j]));
               assert( SCIPisFeasEQ(scip, val, cutcoefs[idx]) );
               obj += val * SCIPvarGetObj(mipdata->alpha[j]);
            }
            else
            {
               if ( (mipdata->coltype[j] == colContinuous || mipdata->coltype[j] == colConverted) && mipdata->isshifted[j] )
                  contVarShifted = TRUE;
            }
         }
         assert( mipdata->beta != NULL );
         val = SCIPgetSolVal(subscip, sol, mipdata->beta);
         assert( SCIPisFeasIntegral(subscip, val) );
         obj += val * SCIPvarGetObj(mipdata->beta);
         assert( contVarShifted || SCIPisFeasEQ(scip, obj, cutact - cutrhs) );
      }
#endif

      /* if successful, convert dense cut into sparse row, and add the row as a cut */
      if ( SCIPisFeasGT(scip, cutact, cutrhs) )
      {
         SCIP_Real cutnorm;
         int cutlen;

         /* store the cut as sparse row, calculate activity and norm of cut */
         SCIP_CALL( storeCutInArrays(scip, nvars, vars, cutcoefs, varsolvals, normtype,
               cutvars, cutvals, &cutlen, &cutact, &cutnorm) );

         SCIPdebugMessage("act=%f, rhs=%f, norm=%f, eff=%f\n", cutact, cutrhs, cutnorm, (cutact - cutrhs)/cutnorm);

         /* if norm is 0, the cut is trivial */
         if ( SCIPisPositive(scip, cutnorm) )
         {
            SCIP_Bool violated = SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm);

            if ( violated || (sepadata->usecutpool && ! cutislocal ) )
            {
               SCIP_ROW* cut;

               /* create the cut */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cgcut%d_%u", SCIPgetNLPs(scip), *ngen);
               SCIP_CALL( SCIPcreateEmptyRow(scip, &cut, name, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );
               SCIP_CALL( SCIPaddVarsToRow(scip, cut, cutlen, cutvars, cutvals) );
               /*SCIPdebug(SCIPprintRow(scip, cut, NULL));*/

               /* add cut to pool */
               if ( ! cutislocal )
               {
                  assert( violated || sepadata->usecutpool );
                  SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               }

               /* add cut if it is violated */
               if ( violated )
               {
                  /* check whether cut has been found before - may happend due to projection */
                  for (k = 0; k < nprevrows; ++k)
                  {
                     SCIP_Real parval;

                     assert( prevrows[k] != NULL );
                     parval = SCIProwGetParallelism(cut, prevrows[k], 'e');
                     /* exit if row is parallel to existing cut and rhs is not better */
                     if ( SCIPisEQ(scip, parval, 1.0) && SCIPisGE(scip, cutrhs, SCIProwGetRhs(prevrows[k])) )
                        break;
                  }

                  /* if cut is new */
                  if ( k >= nprevrows )
                  {
                     prevrows[nprevrows++] = cut;
                     assert( nprevrows <= nsols );

                     SCIPdebugMessage(" -> CG-cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                        name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                        SCIPgetCutEfficacy(scip, NULL, cut),
                        SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                        SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
                     SCIPdebug(SCIPprintRow(scip, cut, NULL));
                     SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE) );
                     ++(*ngen);
                  }
                  else
                  {
                     SCIPdebugMessage("Cut already exists.\n");
                     /* release the row */
                     SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  }
               }
               else
               {
                  /* release the row */
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }
            }
         }
      }
   }
   assert( nprevrows <= nsols );

   /* release rows */
   for (k = 0; k < nprevrows; ++k)
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(prevrows[k])) );
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &prevrows);
   SCIPfreeBufferArray(scip, &cutvals);
   SCIPfreeBufferArray(scip, &cutvars);
   SCIPfreeBufferArrayNull(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &cutcoefs);

   return SCIP_OKAY;
}


/** create CG-cuts via CMIR-function */
static
SCIP_RETCODE createCGCutsCMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data */
   CGMIP_MIPDATA*        mipdata,            /**< data for sub-MIP */
   unsigned int*         ngen                /**< number of generated cuts */
   )
{
   SCIP* subscip;
   SCIP_STAGE stage;
   SCIP_Real* weights;
   SCIP_SOL** sols;
   SCIP_ROW** prevrows;
   int nprevrows;
   int nsols;
   int k, s;
   char name[SCIP_MAXSTRLEN];

   SCIP_VAR** vars;
   SCIP_ROW** rows;
   SCIP_Real* cutcoefs;
   SCIP_VAR** cutvars;
   SCIP_Real* cutvals;
   SCIP_Real* varsolvals;
   SCIP_Bool success;
   SCIP_Bool cutislocal;
   char normtype;
   int nvars;
   int nrows;

   SCIP_Real cutrhs;
   SCIP_Real cutact;
   SCIP_Real maxscale;
   SCIP_Longint maxdnom;

   int* boundsfortrans;
   SCIP_BOUNDTYPE* boundtypesfortrans;

   assert( scip != NULL );
   assert( sepadata != NULL );
   assert( mipdata != NULL );
   assert( ngen != NULL );

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   SCIPdebugMessage("Trying to generate cuts via CMIR routines (use own bounds: %u) ...\n", sepadata->cmirownbounds);
   *ngen = 0;

   /* check if solving was successful and get solutions */
   stage = SCIPgetStage(subscip);
   if ( stage == SCIP_STAGE_SOLVING || stage == SCIP_STAGE_SOLVED )
      nsols = SCIPgetNSols(subscip);
   else
      nsols = 0;

   /* only if solutions have been found */
   if ( nsols == 0 )
      return SCIP_OKAY;

   sols = SCIPgetSols(subscip);

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   assert( nrows > 0 );
   assert( (int) mipdata->nrows == nrows );

   /* @todo more advanced settings - compare sepa_gomory.c */
   maxdnom = (SCIP_Longint) CUTCOEFBND+1;
   maxscale = 10000.0;

   /* get the type of norm to use for efficacy calculations */
   SCIP_CALL( SCIPgetCharParam(scip, "separating/efficacynorm", &normtype) );

   /* get variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &prevrows, nsols) );

   /* get solution values */
   for (k = 0; k < nvars; ++k)
   {
      /* do not have solution values for variables that are not in the LP (varstatus != column) */
      if ( SCIPvarGetStatus(vars[k]) == SCIP_VARSTATUS_COLUMN )
         varsolvals[k] = SCIPvarGetLPSol(vars[k]);
      else
         varsolvals[k] = 0.0;
   }

   /* prepare arrays for bound information, if requested */
   if ( sepadata->cmirownbounds )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &boundsfortrans, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundtypesfortrans, nvars) );
   }
   else
   {
      boundsfortrans = NULL;
      boundtypesfortrans = NULL;
   }

   /* loop through solutions found */
   nprevrows = 0;
   for (s = 0; s < nsols; ++s)
   {
      SCIP_SOL* sol;
      sol = sols[s];

      /* generate weights */
      for (k = 0; k < nrows; ++k)
      {
         SCIP_Real val;

         weights[k] = 0;
         if ( mipdata->ylhs[k] != NULL )
         {
            assert( !SCIProwIsModifiable(rows[k]) && (!SCIProwIsLocal(rows[k]) || sepadata->allowlocal) );

            val = SCIPgetSolVal(subscip, sol, mipdata->ylhs[k]);
            assert( ! SCIPisFeasNegative(subscip, val) );

            if ( SCIPisFeasPositive(subscip, val) )
               weights[k] = -val;
         }
         if ( mipdata->yrhs[k] != NULL )
         {
            assert( !SCIProwIsModifiable(rows[k]) && (!SCIProwIsLocal(rows[k]) || sepadata->allowlocal) );

            val = SCIPgetSolVal(subscip, sol, mipdata->yrhs[k]);
            assert( ! SCIPisFeasNegative(subscip, val) );

            /* in a suboptimal solution both values may be positive - take the one with larger absolute value */
            if ( SCIPisFeasGT(scip, val, ABS(weights[k])) )
               weights[k] = val;
         }
      }

      /* set up data for bounds to use */
      if ( sepadata->cmirownbounds )
      {
         int typefortrans;

         assert( boundsfortrans != NULL );
         assert( boundtypesfortrans != NULL );

         if ( sepadata->allowlocal )
            typefortrans = -2;
         else
            typefortrans = -1;

         /* check all variables */
         for (k = 0; k < nvars; ++k)
         {
            int pos;
            SCIP_VAR* var;

            var = vars[k];
            assert( var != NULL );
            pos = SCIPcolGetLPPos(SCIPvarGetCol(var));

            if ( pos < 0 )
               continue;

            assert( pos < (int) mipdata->ncols );
            assert( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN );

            boundsfortrans[k] = typefortrans;
            boundtypesfortrans[k] = SCIP_BOUNDTYPE_LOWER;

            if ( mipdata->coltype[pos] == colContinuous || mipdata->coltype[pos] == colConverted )
            {
               assert( SCIPvarIsIntegral(var) || mipdata->coltype[pos] != colContinuous );
               continue;
            }

            /* check upper bound */
            if ( mipdata->z[pos] != NULL && SCIPisSumPositive(subscip, SCIPgetSolVal(subscip, sol, mipdata->z[pos])) )
            {
               /* check whether variable is complemented */
               if ( ! mipdata->iscomplemented[pos] )
                  boundtypesfortrans[k] = SCIP_BOUNDTYPE_UPPER;
               /* otherwise use lower bound */
            }
            else
            {
               /* check whether variable is complemented */
               if ( mipdata->iscomplemented[pos] )
                  boundtypesfortrans[k] = SCIP_BOUNDTYPE_UPPER;
               /* otherwise use lower bound */
            }
         }
      }

      /* create a MIR cut using the above calculated weights */
      cutact = -1.0;
      cutrhs = -1.0;
      SCIP_CALL( SCIPcalcMIR(scip, NULL, BOUNDSWITCH, USEVBDS, sepadata->allowlocal, FIXINTEGRALRHS, boundsfortrans, boundtypesfortrans,
            (int) MAXAGGRLEN(nvars), MAXWEIGHTRANGE, MINFRAC, MAXFRAC,
            weights, 1.0, NULL, NULL, cutcoefs, &cutrhs, &cutact, &success, &cutislocal) );
      assert( sepadata->allowlocal || !cutislocal );
      SCIPdebugMessage(" -> success=%u: %g <= %g\n", success, cutact, cutrhs);

      /* if successful, convert dense cut into sparse row, and add the row as a cut */
      if ( success && SCIPisFeasGT(scip, cutact, cutrhs) )
      {
         SCIP_Real cutnorm;
         int cutlen;

         /* store the cut as sparse row, calculate activity and norm of cut */
         SCIP_CALL( storeCutInArrays(scip, nvars, vars, cutcoefs, varsolvals, normtype,
               cutvars, cutvals, &cutlen, &cutact, &cutnorm) );

         /* only proceed if norm is positive - otherwise the cut is trivial */
         if ( SCIPisPositive(scip, cutnorm) )
         {
            SCIP_Bool violated;

            violated = SCIPisEfficacious(scip, (cutact - cutrhs)/cutnorm);

            /* only if the cut if violated - if it is not violated we might store non-local cuts in the pool */
            if ( violated || ( sepadata->usecutpool && ! cutislocal) )
            {
               SCIP_ROW* cut;

               /* create the cut */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cgcut%d_%u", SCIPgetNLPs(scip), *ngen);
               SCIP_CALL( SCIPcreateEmptyRow(scip, &cut, name, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, sepadata->dynamiccuts) );
               SCIP_CALL( SCIPaddVarsToRow(scip, cut, cutlen, cutvars, cutvals) );
               assert( success );

#ifdef SCIP_DEBUG
               SCIPdebug(SCIPprintRow(scip, cut, NULL));
#endif

               /* try to scale the cut to integral values */
               SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                     maxdnom, maxscale, MAKECONTINTEGRAL, &success) );

               /* if the cut could be made integral */
               if ( success )
               {
                  /* add cut to pool */
                  if ( !cutislocal )
                  {
                     assert( violated || sepadata->usecutpool );
                     SCIP_CALL( SCIPaddPoolCut(scip, cut) );
                  }

                  if ( !SCIPisCutEfficacious(scip, NULL, cut) )
                  {
                     SCIPdebugMessage(" -> CG-cut <%s> no longer efficacious: act=%f, rhs=%f, norm=%f, eff=%f\n",
                        name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                        SCIPgetCutEfficacy(scip, NULL, cut));
                     success = FALSE;

                     /* release the row */
                     SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                  }
                  else
                  {
                     /* check whether cut has been found before - may happend due to projection */
                     for (k = 0; k < nprevrows; ++k)
                     {
                        SCIP_Real parval;

                        assert( prevrows[k] != NULL );
                        parval = SCIProwGetParallelism(cut, prevrows[k], 'e');
                        /* exit if row is parallel to existing cut and rhs is not better */
                        if ( SCIPisEQ(scip, parval, 1.0) && SCIPisGE(scip, cutrhs, SCIProwGetRhs(prevrows[k])) )
                           break;
                     }

                     /* if cut is new */
                     if ( k >= nprevrows )
                     {
                        prevrows[nprevrows++] = cut;
                        assert( nprevrows <= nsols );

                        SCIPdebugMessage(" -> CG-cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                           name, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                           SCIPgetCutEfficacy(scip, NULL, cut),
                           SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                           SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));
#ifdef SCIP_OUTPUT
                        SCIPdebug(SCIPprintRow(scip, cut, NULL));
#endif
                        SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE) );
                        ++(*ngen);
                     }
                     else
                     {
                        SCIPdebugMessage("Cut already exists.\n");
                        /* release the row */
                        SCIP_CALL( SCIPreleaseRow(scip, &cut) );
                     }
                  }
               }
               else
               {
                  SCIPdebugMessage(" -> CG-cut <%s> could not be scaled to integral coefficients: act=%f, rhs=%f, norm=%f, eff=%f\n",
                     name, cutact, cutrhs, cutnorm, SCIPgetCutEfficacy(scip, NULL, cut));

                  /* release the row */
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }
            }
         }
      }
   }
   assert( nprevrows <= nsols );

   /* release rows */
   for (k = 0; k < nprevrows; ++k)
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(prevrows[k])) );
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &prevrows);
   SCIPfreeBufferArrayNull(scip, &boundsfortrans);
   SCIPfreeBufferArrayNull(scip, &boundtypesfortrans);

   SCIPfreeBufferArray(scip, &cutvals);
   SCIPfreeBufferArray(scip, &cutvars);
   SCIPfreeBufferArrayNull(scip, &varsolvals);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &weights);

   return SCIP_OKAY;
}


/** frees "subscip" data */
static
SCIP_RETCODE freeSubscip(
   SCIP*                 scip,               /**< SCIP data structure */
   CGMIP_MIPDATA*        mipdata             /**< data for sub-MIP */
   )
{
   unsigned int i, j;
   SCIP* subscip;

   assert( scip != NULL );
   assert( mipdata != NULL );

   SCIPdebugMessage("Freeing subscip ...\n");

   subscip = mipdata->subscip;
   assert( subscip != NULL );

   for (j = 0; j < mipdata->ncols; ++j)
   {
      if ( mipdata->coltype[j] == colPresent )
      {
         assert( mipdata->alpha[j] != NULL );
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->alpha[j])) );
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->fracalpha[j])) );
      }
   }
   SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->beta)) );
   SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->fracbeta)) );

   for (i = 0; i < mipdata->nrows; ++i)
   {
      if ( mipdata->ylhs[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->ylhs[i])) );
      }
      if ( mipdata->yrhs[i] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->yrhs[i])) );
      }
   }

   for (j = 0; j < mipdata->ncols; ++j)
   {
      if ( mipdata->z[j] != NULL )
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &(mipdata->z[j])) );
      }
   }

   if ( mipdata->subscip != NULL )
   {
      SCIP_CALL( SCIPfree(&(mipdata->subscip)) );
   }

   SCIPfreeBlockMemoryArray(scip, &(mipdata->z), 2*mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->yrhs), mipdata->nrows);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->ylhs), mipdata->nrows);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->isshifted), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->iscomplemented), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->coltype), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->fracalpha), mipdata->ncols);
   SCIPfreeBlockMemoryArray(scip, &(mipdata->alpha), mipdata->ncols);

   return SCIP_OKAY;
}



/*
 * Callback methods
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyCGMIP)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaCGMIP(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeCGMIP)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert( sepadata != NULL );

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called when problem solving starts) */
#define sepaInitCGMIP NULL


/** deinitialization method of separator (called when problem solving exits) */
#define sepaExitCGMIP NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolCGMIP NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolCGMIP NULL


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpCGMIP)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   CGMIP_MIPDATA* mipdata;

   int depth;
   int ncalls;
   int ncols;
   int nrows;
   unsigned int ngen;

   SCIP_Bool success;

   assert( scip != NULL );
   assert( sepa != NULL );
   assert( strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;
   ngen = 0;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   depth = SCIPgetDepth(scip);

   /* only call separator, if we are not close to terminating */
   if ( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call separator up to a maximum depth */
   if ( sepadata->maxdepth >= 0 && depth > sepadata->maxdepth )
      return SCIP_OKAY;

   /* only call separator a given number of times at each node */
   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   if ( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* skip separation if there are continuous variables, but only integers required */
   if ( SCIPgetNContVars(scip) > 0 && sepadata->onlyintvars )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if ( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* get LP data */
   ncols = SCIPgetNLPCols(scip);
   nrows = SCIPgetNLPRows(scip);
   if ( ncols <= NCOLSTOOSMALL || nrows <= NROWSTOOSMALL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("separating CG-cuts via sub-MIPs: %d cols, %d rows\n", ncols, nrows);

   /* prepare data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &mipdata) );
   mipdata->subscip = NULL;
   mipdata->alpha = NULL;
   mipdata->fracalpha = NULL;
   mipdata->beta = NULL;
   mipdata->fracbeta = NULL;
   mipdata->coltype = NULL;
   mipdata->iscomplemented = NULL;
   mipdata->isshifted = NULL;
   mipdata->ylhs = NULL;
   mipdata->yrhs = NULL;
   mipdata->z = NULL;

   /* create subscip */
   SCIP_CALL( createSubscip(scip, sepadata, mipdata) );

   /* set parameters */
   SCIP_CALL( subscipSetParams(scip, sepadata, mipdata, &success) );

   if ( success && !SCIPisStopped(scip) )
   {      
      /* solve subscip */
      SCIP_CALL( solveSubscip(scip, sepadata, mipdata, &success) );

#ifdef SCIP_OUPUT
      /* print statistics */
      SCIP_CALL( SCIPprintStatistics(mipdata->subscip, NULL) );
#endif

      /* preceed if solution was successful */
      if ( success && !SCIPisStopped(scip) )
      {
         if ( sepadata->usecmir )
         {
            SCIP_CALL( createCGCutsCMIR(scip, sepadata, mipdata, &ngen) );
         }
         else
         {
            SCIP_CALL( createCGCutsDirect(scip, sepadata, mipdata, &ngen) );
         }
      }
   }

   SCIP_CALL( freeSubscip(scip, mipdata) );
   SCIPfreeBlockMemory(scip, &mipdata);

   SCIPdebugMessage("Found %u CG-cuts.\n", ngen);

   if ( ngen > 0 )
      *result = SCIP_SEPARATED;

#ifdef SCIP_OUTPUT
   /* SCIP_CALL( SCIPwriteLP(scip, "cuts.lp") ); */
   /* SCIP_CALL( SCIPwriteMIP(scip, "cuts.lp", FALSE, TRUE) ); */
#endif

   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
#define sepaExecsolCGMIP NULL




/*
 * separator specific interface methods
 */

/** creates the CGMIP MIR cut separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaCGMIP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST, 
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyCGMIP, sepaFreeCGMIP, sepaInitCGMIP, sepaExitCGMIP,
         sepaInitsolCGMIP, sepaExitsolCGMIP, sepaExeclpCGMIP, sepaExecsolCGMIP,
         sepadata) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cgmip/maxrounds",
         "maximal number of cgmip separation rounds per node (-1: unlimited)",
         &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cgmip/maxroundsroot",
         "maximal number of cgmip separation rounds in the root node (-1: unlimited)",
         &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cgmip/maxdepth",
         "maximal depth at which the separator is applied (-1: unlimited)",
         &sepadata->maxdepth, FALSE, DEFAULT_MAXDEPTH, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/dynamiccuts",
         "should generated cuts be removed from the LP if they are no longer tight?",
         &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cgmip/timelimit",
         "time limit for sub-MIP",
         &sepadata->timelimit, TRUE, DEFAULT_TIMELIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cgmip/memorylimit",
         "memory limit for sub-MIP",
         &sepadata->memorylimit, TRUE, DEFAULT_MEMORYLIMIT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "separating/cgmip/nodelimit",
         "node limit for sub-MIP (-1: unlimited)",
         &sepadata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cgmip/objweight",
         "weight used for the row combination coefficient in the sub-MIP objective",
         &sepadata->objweight, TRUE, DEFAULT_OBJWEIGHT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/usecmir",
         "use CMIR-generator (otherwise add cut directly)?",
         &sepadata->usecmir, FALSE, DEFAULT_USECMIR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/cmirownbounds",
         "tell CMIR-generator which bounds to used in rounding?",
         &sepadata->cmirownbounds, FALSE, DEFAULT_CMIROWNBOUNDS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/allowlocal",
         "allow to generate local cuts?",
         &sepadata->allowlocal, FALSE, DEFAULT_ALLOWLOCAL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/onlyintvars",
         "generate cuts for problems with only integer variables?",
         &sepadata->onlyintvars, FALSE, DEFAULT_ONLYINTVARS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/onlyactiverows",
         "use only active rows to generate cuts?",
         &sepadata->onlyactiverows, FALSE, DEFAULT_ONLYACTIVEROWS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cgmip/maxrowage",
         "maximal age of rows to consider if onlyactiverows is false",
         &sepadata->maxrowage, FALSE, DEFAULT_MAXROWAGE, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/usecutpool",
         "use cutpool to store CG-cuts even if the are not efficient?",
         &sepadata->usecutpool, FALSE, DEFAULT_USECUTPOOL, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/primalseparation",
         "only separate cuts that are tight for the best feasible solution?",
         &sepadata->primalseparation, FALSE, DEFAULT_PRIMALSEPARATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/onlyrankone",
         "whether only rank 1 inequalities should be separated",
         &sepadata->onlyrankone, FALSE, DEFAULT_ONLYRANKONE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/earlyterm",
         "terminate separation if a violated (but possibly sub-optimal) cut has been found?",
         &sepadata->earlyterm, FALSE, DEFAULT_EARLYTERM, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/addviolationcons",
         "add constraint to subscip that only allows violated cuts?",
         &sepadata->addviolationcons, FALSE, DEFAULT_ADDVIOLATIONCONS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/addviolconshdlr",
         "add constraint handler to filter out violated cuts?",
         &sepadata->addviolconshdlr, FALSE, DEFAULT_ADDVIOLCONSHDLR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/conshdlrusenorm",
         "should the violation constraint handler use the norm of a cut to check for feasibility?",
         &sepadata->conshdlrusenorm, FALSE, DEFAULT_CONSHDLRUSENORM, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/objlone",
         "should the objective of the sub-MIP minimize the l1-norm of the multipliers?",
         &sepadata->objlone, FALSE, DEFAULT_OBJLONE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/cgmip/contconvert",
         "convert some integral variables to be continuous to reduce the size of the sub-MIP?",
         &sepadata->contconvert, FALSE, DEFAULT_CONTCONVERT, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "separating/cgmip/contconvfrac",
         "fraction of integral variables converted to be continuous (if contconvert)",
         &sepadata->contconvfrac, FALSE, DEFAULT_CONTCONVFRAC, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "separating/cgmip/contconvmin",
         "minimum number of integral variables before some are converted to be continuous",
         &sepadata->contconvmin, FALSE, DEFAULT_CONTCONVMIN, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
