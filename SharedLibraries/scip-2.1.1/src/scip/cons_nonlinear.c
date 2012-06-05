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

/**@file   cons_nonlinear.c
 * @brief  constraint handler for nonlinear constraints \f$\textrm{lhs} \leq \sum_{i=1}^n a_ix_i + \sum_{j=1}^m c_jf_j(x) \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>

/* for the MS compiler, the function finite(a) is named _finite(a) */
#ifdef _MSC_VER
#define finite(a) _finite(a)
#endif

/* on SunOS, the function finite(a) is declared in ieeefp.h
 * but this header does not exist on every system, so include only if __sun is defined
 */
#ifdef __sun
#include <ieeefp.h>
#endif

#include "scip/cons_nonlinear.h"
#include "scip/cons_linear.h"
#include "scip/heur_trysol.h"
#include "scip/heur_subnlp.h"
#include "nlpi/exprinterpret.h"
#include "scip/debug.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "nonlinear"
#define CONSHDLR_DESC          "constraint handler for nonlinear constraints"
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -60 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000010 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_PROP_TIMING SCIP_PROPTIMING_BEFORELP

#define INTERVALINFTY             1E+43 /**< value for infinity in interval operations */
#define BOUNDTIGHTENING_MINSTRENGTH 0.05/**< minimal required bound tightening strength in expression graph domain tightening for propagating bound change */
#define INITLPMAXVARVAL          1000.0 /**< maximal absolute value of variable for still generating a linearization cut at that point in initlp */

/*
 * Data structures
 */

/** event data for linear variable bound change events */
struct SCIP_EventData
{
   SCIP_CONSHDLRDATA*    conshdlrdata;       /**< the constraint handler data */
   SCIP_CONSDATA*        consdata;           /**< the constraint data */
   int                   varidx;             /**< the index of the linear variable which bound change is catched */
   int                   filterpos;          /**< position of eventdata in SCIP's event filter */
};

/** constraint data for nonlinear constraints */
struct SCIP_ConsData
{
   SCIP_Real             lhs;                /**< left hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */

   int                   nlinvars;           /**< number of linear variables */
   int                   linvarssize;        /**< length of linear variable arrays */
   SCIP_VAR**            linvars;            /**< linear variables */
   SCIP_Real*            lincoefs;           /**< coefficients of linear variables */
   SCIP_EVENTDATA**      lineventdata;       /**< eventdata for bound change of linear variable */

   int                   nexprtrees;         /**< number of expression trees */
   SCIP_Real*            nonlincoefs;        /**< coefficients of expression trees */
   SCIP_EXPRTREE**       exprtrees;          /**< nonlinear part of constraint */
   SCIP_EXPRCURV*        curvatures;         /**< curvature of each expression tree (taking nonlincoefs into account) */
   SCIP_EXPRGRAPHNODE*   exprgraphnode;      /**< node in expression graph corresponding to expression tree of this constraint */
   SCIP_EXPRCURV         curvature;          /**< curvature of complete nonlinear part, if checked */

   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */

   unsigned int          linvarssorted:1;    /**< are the linear variables already sorted? */
   unsigned int          linvarsmerged:1;    /**< are equal linear variables already merged? */

   unsigned int          iscurvchecked:1;    /**< is nonlinear function checked on convexity or concavity ? */
   unsigned int          isremovedfixingslin:1; /**< did we removed fixed/aggr/multiaggr variables in linear part? */
   unsigned int          ispropagated:1;     /**< did we propagate the current bounds of linear variables in this constraint? */
   unsigned int          ispresolved:1;      /**< did we checked for possibilities of upgrading or implicit integer variables ? */

   SCIP_Real             minlinactivity;     /**< sum of minimal activities of all linear terms with finite minimal activity */
   SCIP_Real             maxlinactivity;     /**< sum of maximal activities of all linear terms with finite maximal activity */
   int                   minlinactivityinf;  /**< number of linear terms with infinite minimal activity */
   int                   maxlinactivityinf;  /**< number of linear terms with infinity maximal activity */
   SCIP_Real             activity;           /**< activity of constraint function w.r.t. current solution */
   SCIP_Real             lhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of lower bound by current solution (used temporarily inside constraint handler) */

   int                   linvar_maydecrease; /**< index of a variable in linvars that may be decreased without making any other constraint infeasible, or -1 if none */
   int                   linvar_mayincrease; /**< index of a variable in linvars that may be increased without making any other constraint infeasible, or -1 if none */

   SCIP_Real             lincoefsmin;        /**< maximal absolute value of coefficients in linear part, only available in solving stage */
   SCIP_Real             lincoefsmax;        /**< minimal absolute value of coefficients in linear part, only available in solving stage */
};

/** nonlinear constraint update method */
struct SCIP_NlConsUpgrade
{
   SCIP_DECL_NONLINCONSUPGD((*nlconsupgd));    /**< method to call for upgrading nonlinear constraint */
   SCIP_DECL_EXPRGRAPHNODEREFORM((*nodereform));/**< method to call for reformulating an expression graph node */
   int                     priority;           /**< priority of upgrading method */
   SCIP_Bool               active;             /**< is upgrading enabled */
};
typedef struct SCIP_NlConsUpgrade SCIP_NLCONSUPGRADE;

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_EXPRINT*         exprinterpreter;    /**< expression interpreter to compute gradients */

   SCIP_Real             mincutefficacysepa; /**< minimal efficacy of a cut in order to add it to relaxation during separation */
   SCIP_Real             mincutefficacyenfofac;/**< minimal target efficacy of a cut in order to add it to relaxation during enforcement as factor of feasibility tolerance (may be ignored) */
   SCIP_Real             cutmaxrange;        /**< maximal range (maximal coef / minimal coef) of a cut in order to be added to LP */
   SCIP_Bool             linfeasshift;       /**< whether to make solutions in check feasible if possible */
   SCIP_Bool             checkconvexexpensive;/**< whether to apply expensive curvature checking methods */
   SCIP_Bool             assumeconvex;       /**< whether functions in inequalities should be assumed to be convex */
   int                   maxproprounds;      /**< limit on number of propagation rounds for a single constraint within one round of SCIP propagation */
   SCIP_Bool             reformulate;        /**< whether to reformulate expression graph */
   int                   maxexpansionexponent;/**< maximal exponent where still expanding non-monomial polynomials in expression simplification */
   SCIP_Real             sepanlpmincont;     /**< minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation */

   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subNLP heuristic, if available */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the TRYSOL heuristic, if available */
   SCIP_EVENTHDLR*       linvareventhdlr;    /**< our handler for linear variable bound change events */
   SCIP_EVENTHDLR*       nonlinvareventhdlr; /**< our handler for nonlinear variable bound change events */
   int                   newsoleventfilterpos;/**< filter position of new solution event handler, if catched */

   SCIP_NLCONSUPGRADE**  nlconsupgrades;     /**< nonlinear constraint upgrade methods for specializing nonlinear constraints */
   int                   nlconsupgradessize; /**< size of nlconsupgrade array */
   int                   nnlconsupgrades;    /**< number of nonlinear constraint upgrade methods */

   SCIP_EXPRGRAPH*       exprgraph;          /**< expression graph */
   SCIP*                 scip;               /**< SCIP pointer for use in expression graph callbacks */
   unsigned int          isremovedfixings:1; /**< have fixed variables been removed in the expression graph? */
   unsigned int          ispropagated:1;     /**< have current bounds of linear variables in constraints and variables in expression graph been propagated? */
   unsigned int          isreformulated:1;   /**< has expression graph been reformulated? */
   int                   naddedreformconss;  /**< number of constraints added via reformulation */
   SCIP_Bool             sepanlp;            /**< where linearization of the NLP relaxation solution added? */
   SCIP_NODE*            lastenfolpnode;     /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenfolprounds;      /**< counter on number of enforcement rounds for the current node */
};

/*
 * Local methods
 */

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/* catches variable bound change events on a linear variable in a nonlinear constraint */
static
SCIP_RETCODE catchLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EVENTDATA* eventdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsEnabled(cons));
   assert(SCIPconsIsTransformed(cons));

   assert(SCIPconsGetHdlr(cons) != NULL);
   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->linvareventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);

   SCIP_CALL( SCIPallocBlockMemory(scip, &eventdata) );

   eventtype = SCIP_EVENTTYPE_VARFIXED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   eventdata->conshdlrdata = conshdlrdata;
   eventdata->consdata = consdata;
   eventdata->varidx = linvarpos;
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->linvars[linvarpos], eventtype, conshdlrdata->linvareventhdlr, eventdata, &eventdata->filterpos) );

   /* ensure lineventdata array is existing */
   if( consdata->lineventdata == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize) );
   }

   consdata->lineventdata[linvarpos] = eventdata;

   /* since bound changes were not catched before, a possibly stored linear activity may have become outdated, so set to invalid */
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->ispropagated   = FALSE;

   return SCIP_OKAY;
}

/* drops variable bound change events on a linear variable in a nonlinear constraint */
static
SCIP_RETCODE dropLinearVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint for which to catch bound change events */
   int                   linvarpos           /**< position of variable in linear variables array */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   assert(SCIPconsGetHdlr(cons) != NULL);
   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->linvareventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(linvarpos >= 0);
   assert(linvarpos < consdata->nlinvars);
   assert(consdata->lineventdata != NULL);
   assert(consdata->lineventdata[linvarpos] != NULL);
   assert(consdata->lineventdata[linvarpos]->consdata == consdata);
   assert(consdata->lineventdata[linvarpos]->varidx == linvarpos);
   assert(consdata->lineventdata[linvarpos]->filterpos >= 0);

   eventtype = SCIP_EVENTTYPE_VARFIXED;
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* if right hand side is finite, then a tightening in the lower bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
   }
   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* if left hand side is finite, then a tightening in the upper bound of coef*linvar is of interest */
      if( consdata->lincoefs[linvarpos] > 0.0 )
         eventtype |= SCIP_EVENTTYPE_UBCHANGED;
      else
         eventtype |= SCIP_EVENTTYPE_LBCHANGED;
   }

   SCIP_CALL( SCIPdropVarEvent(scip, consdata->linvars[linvarpos], eventtype, conshdlrdata->linvareventhdlr, consdata->lineventdata[linvarpos], consdata->lineventdata[linvarpos]->filterpos) );

   SCIPfreeBlockMemory(scip, &consdata->lineventdata[linvarpos]);  /*lint !e866*/

   return SCIP_OKAY;
}

/** locks a linear variable in a constraint */
static
SCIP_RETCODE lockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to lock a variable */
   SCIP_VAR*             var,                /**< variable to lock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** unlocks a linear variable in a constraint */
static
SCIP_RETCODE unlockLinearVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where to unlock a variable */
   SCIP_VAR*             var,                /**< variable to unlock */
   SCIP_Real             coef                /**< coefficient of variable in constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var != NULL);
   assert(coef != 0.0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( coef > 0.0 )
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
   }
   else
   {
      SCIP_CALL( SCIPunlockVarCons(scip, var, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
   }

   return SCIP_OKAY;
}

/** computes the minimal and maximal activity for the linear part in a constraint data
 *  only sums up terms that contribute finite values
 *  gives the number of terms that contribute infinite values
 *  only computes those activities where the corresponding side of the constraint is finite
 */
static
void consdataUpdateLinearActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{  /*lint --e{666}*/
   SCIP_ROUNDMODE prevroundmode;
   int       i;
   SCIP_Real bnd;

   assert(scip != NULL);
   assert(consdata != NULL);

   if( consdata->minlinactivity != SCIP_INVALID && consdata->maxlinactivity != SCIP_INVALID )  /*lint !e777*/
   {
      /* activities should be uptodate */
      assert(consdata->minlinactivityinf >= 0);
      assert(consdata->maxlinactivityinf >= 0);
      return;
   }

   consdata->minlinactivityinf = 0;
   consdata->maxlinactivityinf = 0;

   /* if lhs is -infinite, then we do not compute a maximal activity, so we set it to  infinity
    * if rhs is  infinite, then we do not compute a minimal activity, so we set it to -infinity
    */
   consdata->minlinactivity = SCIPisInfinity(scip,  consdata->rhs) ? -INTERVALINFTY : 0.0;
   consdata->maxlinactivity = SCIPisInfinity(scip, -consdata->lhs) ?  INTERVALINFTY : 0.0;

   if( consdata->nlinvars == 0 )
      return;

   /* if the activities computed here should be still uptodate after boundchanges,
    * variable events need to be catched */
   assert(consdata->lineventdata != NULL);

   prevroundmode = SCIPintervalGetRoundingMode();

   if( !SCIPisInfinity(scip,  consdata->rhs) )
   {
      /* compute minimal activity only if there is a finite right hand side */
      SCIPintervalSetRoundingModeDownwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(consdata->lineventdata[i] != NULL);
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = SCIPvarGetLbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         else
         {
            bnd = SCIPvarGetUbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->minlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         consdata->minlinactivity += consdata->lincoefs[i] * bnd;
      }
   }

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* compute maximal activity only if there is a finite left hand side */
      SCIPintervalSetRoundingModeUpwards();

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         assert(consdata->lineventdata[i] != NULL);
         if( consdata->lincoefs[i] >= 0.0 )
         {
            bnd = SCIPvarGetUbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip,  bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, -bnd)); /* do not like variables that are fixed at -infinity */
         }
         else
         {
            bnd = SCIPvarGetLbLocal(consdata->linvars[i]);
            if( SCIPisInfinity(scip, -bnd) )
            {
               ++consdata->maxlinactivityinf;
               continue;
            }
            assert(!SCIPisInfinity(scip, bnd)); /* do not like variables that are fixed at +infinity */
         }
         consdata->maxlinactivity += consdata->lincoefs[i] * bnd;
      }
   }

   SCIPintervalSetRoundingMode(prevroundmode);
}

/** update the linear activities after a change in the lower bound of a variable */
static
void consdataUpdateLinearActivityLbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with lower bounds at infinity */
   assert(!SCIPisInfinity(scip, oldbnd));
   assert(!SCIPisInfinity(scip, newbnd));

   /* assume lhs <= a*x + y <= rhs, then the following boundchanges can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */

   if( coef > 0.0 )
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->minlinactivity > -INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->minlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->maxlinactivity < INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, -oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->maxlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, -newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** update the linear activities after a change in the upper bound of a variable */
static
void consdataUpdateLinearActivityUbChange(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             coef,               /**< coefficient of variable in constraint */
   SCIP_Real             oldbnd,             /**< previous lower bound of variable */
   SCIP_Real             newbnd              /**< new lower bound of variable */
   )
{
   SCIP_ROUNDMODE prevroundmode;

   assert(scip != NULL);
   assert(consdata != NULL);
   /* we can't deal with upper bounds at -infinity */
   assert(!SCIPisInfinity(scip, -oldbnd));
   assert(!SCIPisInfinity(scip, -newbnd));

   /* assume lhs <= a*x + y <= rhs, then the following boundchanges can be deduced:
    * a > 0:  y <= rhs - a*lb(x),  y >= lhs - a*ub(x)
    * a < 0:  y <= rhs - a*ub(x),  y >= lhs - a*lb(x)
    */

   if( coef > 0.0 )
   {
      /* we should only be called if lhs is finite */
      assert(!SCIPisInfinity(scip, -consdata->lhs));

      /* we have no max activities computed so far, so cannot update */
      if( consdata->maxlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->maxlinactivity < INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeUpwards();

      /* update max activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->maxlinactivityinf;
         assert(consdata->maxlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->maxlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->maxlinactivityinf;
      }
      else
      {
         consdata->maxlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
   else
   {
      /* we should only be called if rhs is finite */
      assert(!SCIPisInfinity(scip, consdata->rhs));

      /* we have no min activities computed so far, so cannot update */
      if( consdata->minlinactivity == SCIP_INVALID )  /*lint !e777*/
         return;

      assert(consdata->minlinactivity > -INTERVALINFTY);

      prevroundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      /* update min activity */
      if( SCIPisInfinity(scip, oldbnd) )
      {
         --consdata->minlinactivityinf;
         assert(consdata->minlinactivityinf >= 0);
      }
      else
      {
         SCIP_Real minuscoef;
         minuscoef = -coef;
         consdata->minlinactivity += minuscoef * oldbnd;
      }

      if( SCIPisInfinity(scip, newbnd) )
      {
         ++consdata->minlinactivityinf;
      }
      else
      {
         consdata->minlinactivity += coef * newbnd;
      }

      SCIPintervalSetRoundingMode(prevroundmode);
   }
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processLinearVarEvent)
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   consdata = eventdata->consdata;
   assert(consdata != NULL);
   assert(eventdata->varidx >= 0);
   assert(eventdata->varidx < consdata->nlinvars);

   eventtype = SCIPeventGetType(event);

   if( eventtype & SCIP_EVENTTYPE_VARFIXED )
   {
      consdata->isremovedfixingslin = FALSE;
   }

   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      /* update activity bounds for linear terms */
      if( eventtype & SCIP_EVENTTYPE_LBCHANGED )
         consdataUpdateLinearActivityLbChange(scip, consdata, consdata->lincoefs[eventdata->varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));
      else
         consdataUpdateLinearActivityUbChange(scip, consdata, consdata->lincoefs[eventdata->varidx], SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

      if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
      {
         assert(eventdata->conshdlrdata != NULL);
         eventdata->conshdlrdata->ispropagated = FALSE;
         consdata->ispropagated = FALSE;
      }
   }

   return SCIP_OKAY;
}

/** processes bound change events for variables in expression graph */
static
SCIP_DECL_EVENTEXEC(processNonlinearVarEvent)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   conshdlrdata = (SCIP_CONSHDLRDATA*)SCIPeventhdlrGetData(eventhdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   eventtype = SCIPeventGetType(event);
   assert( eventtype & (SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED) );

   if( eventtype & SCIP_EVENTTYPE_BOUNDCHANGED )
   {
      SCIP_Real newbd;

      SCIPdebugMessage("changed %s bound on expression graph variable <%s> from %g to %g\n",
         eventtype & SCIP_EVENTTYPE_LBCHANGED ? "lower" : "upper",
         SCIPvarGetName(SCIPeventGetVar(event)), SCIPeventGetOldbound(event), SCIPeventGetNewbound(event));

      if( eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED )
         conshdlrdata->ispropagated = FALSE;
      /* @todo a global bound tightening may yield in convex/concave curvatures, maybe the iscurvcheck flag should be reset? */

      /* update variable bound in expression graph, with epsilon added */
      if( eventtype & SCIP_EVENTTYPE_LBCHANGED )
      {
         newbd = -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -SCIPeventGetNewbound(event));  /*lint !e666*/
         /* if newbd in [0,eps], then relax to 0.0, otherwise relax by -epsilon
          * this is to ensure that an original positive variable does not get negative by this, which may have an adverse effect on convexity recoginition, for example */
         if( newbd >= 0.0 && newbd <= SCIPepsilon(scip) )
            newbd = 0.0;
         else
            newbd -= SCIPepsilon(scip);
         SCIPexprgraphSetVarNodeLb(conshdlrdata->exprgraph, (SCIP_EXPRGRAPHNODE*)eventdata, newbd);
      }
      else
      {
         newbd = +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  SCIPeventGetNewbound(event));  /*lint !e666*/
         /* if newbd in [-eps,0], then relax to 0.0, otherwise relax by +epsilon */
         if( newbd >= -SCIPepsilon(scip) && newbd <= 0.0 )
            newbd = 0.0;
         else
            newbd += SCIPepsilon(scip);
         SCIPexprgraphSetVarNodeUb(conshdlrdata->exprgraph, (SCIP_EXPRGRAPHNODE*)eventdata, newbd);
      }
   }
   else
   {
      assert(eventtype & SCIP_EVENTTYPE_VARFIXED);
      conshdlrdata->isremovedfixings = FALSE;
   }

   return SCIP_OKAY;
}

/** callback method for variable addition in expression graph */
static
SCIP_DECL_EXPRGRAPHVARADDED( exprgraphVarAdded )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_INTERVAL varbounds;
   SCIP_VAR* var_;

   assert(exprgraph != NULL);
   assert(var != NULL);
   assert(varnode != NULL);

   var_ = (SCIP_VAR*)var;

   conshdlrdata = (SCIP_CONSHDLRDATA*)userdata;
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph == exprgraph);

   /* catch variable bound change events */
   SCIP_CALL( SCIPcatchVarEvent(conshdlrdata->scip, (SCIP_VAR*)var, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, conshdlrdata->nonlinvareventhdlr, (SCIP_EVENTDATA*)varnode, NULL) );
   SCIPdebugMessage("catch boundchange events on new expression graph variable <%s>\n", SCIPvarGetName(var_));

   /* set current bounds in expression graph */
   SCIPintervalSetBounds(&varbounds,
      -infty2infty(SCIPinfinity(conshdlrdata->scip), INTERVALINFTY, -MIN(SCIPvarGetLbLocal(var_), SCIPvarGetUbLocal(var_))),  /*lint !e666*/
      +infty2infty(SCIPinfinity(conshdlrdata->scip), INTERVALINFTY,  MAX(SCIPvarGetLbLocal(var_), SCIPvarGetUbLocal(var_)))   /*lint !e666*/
      );
   SCIPexprgraphSetVarNodeBounds(exprgraph, varnode, varbounds);

   SCIP_CALL( SCIPaddVarLocks(conshdlrdata->scip, var_, 1, 1) );
   SCIPdebugMessage("increased up- and downlocks of variable <%s>\n", SCIPvarGetName(var_));

   SCIP_CALL( SCIPcaptureVar(conshdlrdata->scip, var_) );
   SCIPdebugMessage("capture variable <%s>\n", SCIPvarGetName(var_));

   conshdlrdata->isremovedfixings &= SCIPvarIsActive(var_);
   conshdlrdata->ispropagated = FALSE;

   return SCIP_OKAY;
}

/** callback method for variable removal in expression graph */
static
SCIP_DECL_EXPRGRAPHVARREMOVE( exprgraphVarRemove )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR* var_;

   assert(exprgraph != NULL);
   assert(var != NULL);
   assert(varnode != NULL);

   var_ = (SCIP_VAR*)var;

   conshdlrdata = (SCIP_CONSHDLRDATA*)userdata;
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph == exprgraph);

   SCIP_CALL( SCIPdropVarEvent(conshdlrdata->scip, var_, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, conshdlrdata->nonlinvareventhdlr, (SCIP_EVENTDATA*)varnode, -1) );
   SCIPdebugMessage("drop boundchange events on expression graph variable <%s>\n", SCIPvarGetName(var_));

   SCIP_CALL( SCIPaddVarLocks(conshdlrdata->scip, var_, -1, -1) );
   SCIPdebugMessage("decreased up- and downlocks of variable <%s>\n", SCIPvarGetName(var_));

   SCIPdebugMessage("release variable <%s>\n", SCIPvarGetName(var_));
   SCIP_CALL( SCIPreleaseVar(conshdlrdata->scip, &var_) );

   return SCIP_OKAY;
}

/* sets expression trees of constraints, clears previously ones */
static
SCIP_RETCODE consdataSetExprtrees(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< nonlinear constraint data */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            coefs,              /**< coefficients of expression trees, or NULL if all 1.0 */
   SCIP_Bool             copytrees           /**< whether trees should be copied or ownership should be assumed */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->exprtrees != NULL || consdata->nexprtrees == 0);

   /* clear previous expression trees */
   for( i = 0; i < consdata->nexprtrees; ++i )
   {
      assert(consdata->exprtrees[i] != NULL);
      SCIP_CALL( SCIPexprtreeFree(&consdata->exprtrees[i]) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispresolved  = FALSE;

   if( nexprtrees == 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->exprtrees,   consdata->nexprtrees);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->nonlincoefs, consdata->nexprtrees);
      SCIPfreeBlockMemoryArrayNull(scip, &consdata->curvatures,  consdata->nexprtrees);
      consdata->nexprtrees = 0;

      consdata->curvature = SCIP_EXPRCURV_LINEAR;
      consdata->iscurvchecked = TRUE;

      return SCIP_OKAY;
   }

   if( consdata->nexprtrees == 0 )
   {
      assert(consdata->exprtrees   == NULL);
      assert(consdata->nonlincoefs == NULL);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->exprtrees,   nexprtrees) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nonlincoefs, nexprtrees) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->curvatures,  nexprtrees) );
   }
   else
   {
      assert(consdata->exprtrees   != NULL);
      assert(consdata->nonlincoefs != NULL);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->exprtrees,   consdata->nexprtrees, nexprtrees) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nonlincoefs, consdata->nexprtrees, nexprtrees) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->curvatures,  consdata->nexprtrees, nexprtrees) );
   }
   consdata->nexprtrees = nexprtrees;

   for( i = 0; i < nexprtrees; ++i )
   {
      assert(exprtrees[i] != NULL);
      /* the expression tree need to have SCIP_VAR*'s stored */
      assert(SCIPexprtreeGetNVars(exprtrees[i]) == 0 || SCIPexprtreeGetVars(exprtrees[i]) != NULL);

      if( copytrees )
      {
         SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &consdata->exprtrees[i], exprtrees[i]) );
      }
      else
      {
         consdata->exprtrees[i] = exprtrees[i];
      }

      consdata->nonlincoefs[i] = (coefs != NULL ? coefs[i] : 1.0);
      consdata->curvatures[i] = SCIP_EXPRCURV_UNKNOWN;
   }

   consdata->curvature = SCIP_EXPRCURV_UNKNOWN;
   consdata->iscurvchecked = FALSE;

   return SCIP_OKAY;
}

/** ensures, that linear vars and coefs arrays can store at least num entries */
static
SCIP_RETCODE consdataEnsureLinearVarsSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< nonlinear constraint data */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);
   assert(consdata->nlinvars <= consdata->linvarssize);

   if( num > consdata->linvarssize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, num);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->linvars,  consdata->linvarssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lincoefs, consdata->linvarssize, newsize) );
      if( consdata->lineventdata != NULL )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->lineventdata, consdata->linvarssize, newsize) );
      }
      consdata->linvarssize = newsize;
   }
   assert(num <= consdata->linvarssize);

   return SCIP_OKAY;
}

/** creates constraint data structure */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< a buffer to store pointer to new constraint data */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< array of linear variables */
   SCIP_Real*            lincoefs,           /**< array of coefficients of linear variables */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees */
   SCIP_Real*            nonlincoefs,        /**< coefficients of expression trees */
   SCIP_Bool             capturevars         /**< whether we should capture variables */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);

   assert(nlinvars == 0 || linvars  != NULL);
   assert(nlinvars == 0 || lincoefs != NULL);
   assert(nexprtrees == 0 || exprtrees != NULL);
   assert(nexprtrees == 0 || nonlincoefs != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;

   if( nlinvars > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->linvars,  linvars,  nlinvars) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*consdata)->lincoefs, lincoefs, nlinvars) );
      (*consdata)->nlinvars = nlinvars;
      (*consdata)->linvarssize = nlinvars;

      if( capturevars )
         for( i = 0; i < nlinvars; ++i )
         {
            SCIP_CALL( SCIPcaptureVar(scip, linvars[i]) );
         }
   }
   else
   {
      (*consdata)->linvarssorted = TRUE;
      (*consdata)->linvarsmerged = TRUE;
   }

   SCIP_CALL( consdataSetExprtrees(scip, *consdata, nexprtrees, exprtrees, nonlincoefs, TRUE) );

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->activity = SCIP_INVALID;
   (*consdata)->lhsviol  = SCIPisInfinity(scip, -lhs) ? 0.0 : SCIP_INVALID;
   (*consdata)->rhsviol  = SCIPisInfinity(scip,  rhs) ? 0.0 : SCIP_INVALID;

   return SCIP_OKAY;
}

/** creates empty constraint data structure */
static
SCIP_RETCODE consdataCreateEmpty(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< a buffer to store pointer to new constraint data */
   )
{
   assert(scip != NULL);
   assert(consdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );
   BMSclearMemory(*consdata);

   (*consdata)->lhs = -SCIPinfinity(scip);
   (*consdata)->rhs =  SCIPinfinity(scip);

   (*consdata)->linvarssorted    = TRUE;
   (*consdata)->linvarsmerged    = TRUE;

   (*consdata)->isremovedfixingslin = TRUE;
   (*consdata)->ispropagated     = TRUE;

   (*consdata)->linvar_maydecrease = -1;
   (*consdata)->linvar_mayincrease = -1;

   (*consdata)->minlinactivity = SCIP_INVALID;
   (*consdata)->maxlinactivity = SCIP_INVALID;
   (*consdata)->minlinactivityinf = -1;
   (*consdata)->maxlinactivityinf = -1;

   (*consdata)->curvature = SCIP_EXPRCURV_LINEAR;
   (*consdata)->iscurvchecked = TRUE;

   return SCIP_OKAY;
}

/** frees constraint data structure */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to constraint data to free */
   )
{
   int i;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release linear variables and free linear part */
   if( (*consdata)->linvarssize > 0 )
   {
      for( i = 0; i < (*consdata)->nlinvars; ++i )
      {
         assert((*consdata)->lineventdata == NULL || (*consdata)->lineventdata[i] == NULL); /* variable events should have been dropped earlier */
         SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->linvars[i]) );
      }
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->linvars,  (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->lincoefs, (*consdata)->linvarssize);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->lineventdata, (*consdata)->linvarssize);
   }
   assert((*consdata)->linvars == NULL);
   assert((*consdata)->lincoefs == NULL);
   assert((*consdata)->lineventdata == NULL);

   if( (*consdata)->nexprtrees > 0 )
   {
      assert((*consdata)->exprtrees   != NULL);
      assert((*consdata)->nonlincoefs != NULL);
      assert((*consdata)->curvatures  != NULL);
      for( i = 0; i < (*consdata)->nexprtrees; ++i )
      {
         assert((*consdata)->exprtrees[i] != NULL);
         SCIP_CALL( SCIPexprtreeFree(&(*consdata)->exprtrees[i]) );
         assert((*consdata)->exprtrees[i] == NULL);
      }
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->exprtrees,   (*consdata)->nexprtrees);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->nonlincoefs, (*consdata)->nexprtrees);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->curvatures,  (*consdata)->nexprtrees);
   }
   assert((*consdata)->exprtrees   == NULL);
   assert((*consdata)->nonlincoefs == NULL);
   assert((*consdata)->curvatures  == NULL);

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);
   *consdata = NULL;

   return SCIP_OKAY;
}

/** sorts linear part of constraint data */
static
void consdataSortLinearVars(
   SCIP_CONSDATA*        consdata            /**< nonlinear constraint data */
   )
{
   assert(consdata != NULL);

   if( consdata->linvarssorted )
      return;

   if( consdata->nlinvars <= 1 )
   {
      consdata->linvarssorted = TRUE;
      return;
   }

   if( consdata->lineventdata == NULL )
   {
      SCIPsortPtrReal((void**)consdata->linvars, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);
   }
   else
   {
      int i;

      SCIPsortPtrPtrReal((void**)consdata->linvars, (void**)consdata->lineventdata, consdata->lincoefs, SCIPvarComp, consdata->nlinvars);

      /* update variable indices in event data */
      for( i = 0; i < consdata->nlinvars; ++i )
         if( consdata->lineventdata[i] != NULL )
            consdata->lineventdata[i]->varidx = i;
   }

   consdata->linvarssorted = TRUE;
}

/* this function is currently not needed, but also to nice to be deleted, so it is only deactivated */
#if 0
/** returns the position of variable in the linear coefficients array of a constraint, or -1 if not found */
static
int consdataFindLinearVar(
   SCIP_CONSDATA*        consdata,           /**< nonlinear constraint data */
   SCIP_VAR*             var                 /**< variable to search for */
   )
{
   int pos;

   assert(consdata != NULL);
   assert(var != NULL);

   if( consdata->nlinvars == 0 )
      return -1;

   consdataSortLinearVars(consdata);

   if( !SCIPsortedvecFindPtr((void**)consdata->linvars, SCIPvarComp, (void*)var, consdata->nlinvars, &pos) )
      pos = -1;

   return pos;
}
#endif

/** moves a linear variable from one position to another */
static
void consdataMoveLinearVar(
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   int                   oldpos,             /**< position of variable that shall be moved */
   int                   newpos              /**< new position of variable */
   )
{
   assert(consdata != NULL);
   assert(oldpos >= 0);
   assert(oldpos < consdata->nlinvars);
   assert(newpos >= 0);
   assert(newpos < consdata->linvarssize);

   if( newpos == oldpos )
      return;

   consdata->linvars [newpos] = consdata->linvars [oldpos];
   consdata->lincoefs[newpos] = consdata->lincoefs[oldpos];

   if( consdata->lineventdata != NULL )
   {
      assert(newpos >= consdata->nlinvars || consdata->lineventdata[newpos] == NULL);

      consdata->lineventdata[newpos] = consdata->lineventdata[oldpos];
      consdata->lineventdata[newpos]->varidx = newpos;

      consdata->lineventdata[oldpos] = NULL;
   }

   consdata->linvarssorted = FALSE;
}

/** adds linear coefficient in nonlinear constraint */
static
SCIP_RETCODE addLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_VAR*             var,                /**< variable of constraint entry */
   SCIP_Real             coef                /**< coefficient of constraint entry */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool transformed;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);

   /* ignore coefficient if it is nearly zero */
   if( SCIPisZero(scip, coef) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* are we in the transformed problem? */
   transformed = SCIPconsIsTransformed(cons);

   /* always use transformed variables in transformed constraints */
   if( transformed )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, var, &var) );
   }
   assert(var != NULL);
   assert(transformed == SCIPvarIsTransformed(var));

   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars+1) );
   consdata->linvars [consdata->nlinvars] = var;
   consdata->lincoefs[consdata->nlinvars] = coef;

   ++consdata->nlinvars;

   /* catch variable events */
   if( SCIPconsIsEnabled(cons) )
   {
      /* catch bound change events of variable */
      SCIP_CALL( catchLinearVarEvents(scip, cons, consdata->nlinvars-1) );
   }

   /* invalidate activity information */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* install rounding locks for new variable */
   SCIP_CALL( lockLinearVariable(scip, cons, var, coef) );

   /* capture new variable */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;
   consdata->isremovedfixingslin = consdata->isremovedfixingslin && SCIPvarIsActive(var);
   if( consdata->nlinvars == 1 )
      consdata->linvarssorted = TRUE;
   else
      consdata->linvarssorted = consdata->linvarssorted &&
         (SCIPvarCompare(consdata->linvars[consdata->nlinvars-2], consdata->linvars[consdata->nlinvars-1]) == -1);
   /* always set too FALSE since the new linear variable should be checked if already existing as quad var term */
   consdata->linvarsmerged = FALSE;

   return SCIP_OKAY;
}

/** deletes linear coefficient at given position from nonlinear constraint data */
static
SCIP_RETCODE delLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   int                   pos                 /**< position of coefficient to delete */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos && pos < consdata->nlinvars);

   var  = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);

   /* remove rounding locks for deleted variable */
   SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );

   /* if constraint is enabled, drop the events on the variable */
   if( SCIPconsIsEnabled(cons) )
   {
      /* drop bound change events of variable */
      SCIP_CALL( dropLinearVarEvents(scip, cons, pos) );
   }

   /* release variable */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->linvars[pos]) );

   /* move the last variable to the free slot */
   consdataMoveLinearVar(consdata, consdata->nlinvars-1, pos);

   --consdata->nlinvars;

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved  = FALSE;

   return SCIP_OKAY;
}

/** changes linear coefficient value at given position of nonlinear constraint */
static
SCIP_RETCODE chgLinearCoefPos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   int                   pos,                /**< position of linear coefficient to change */
   SCIP_Real             newcoef             /**< new value of linear coefficient */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisZero(scip, newcoef));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(0 <= pos);
   assert(pos < consdata->nlinvars);
   assert(!SCIPisZero(scip, newcoef));

   var = consdata->linvars[pos];
   coef = consdata->lincoefs[pos];
   assert(var != NULL);
   assert(SCIPconsIsTransformed(cons) == SCIPvarIsTransformed(var));

   /* invalidate activity */
   consdata->activity = SCIP_INVALID;
   consdata->minlinactivity = SCIP_INVALID;
   consdata->maxlinactivity = SCIP_INVALID;
   consdata->minlinactivityinf = -1;
   consdata->maxlinactivityinf = -1;

   /* invalidate nonlinear row */
   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   /* if necessary, remove the rounding locks and event catching of the variable */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         assert(SCIPconsIsTransformed(cons));

         /* remove rounding locks for variable with old coefficient */
         SCIP_CALL( unlockLinearVariable(scip, cons, var, coef) );
      }

      if( SCIPconsIsEnabled(cons) )
      {
         /* drop bound change events of variable */
         SCIP_CALL( dropLinearVarEvents(scip, cons, pos) );
      }
   }

   /* change the coefficient */
   consdata->lincoefs[pos] = newcoef;

   /* if necessary, install the rounding locks and event catching of the variable again */
   if( newcoef * coef < 0.0 )
   {
      if( SCIPconsIsLocked(cons) )
      {
         /* install rounding locks for variable with new coefficient */
         SCIP_CALL( lockLinearVariable(scip, cons, var, newcoef) );
      }

      if( SCIPconsIsEnabled(cons) )
      {
         /* catch bound change events of variable */
         SCIP_CALL( catchLinearVarEvents(scip, cons, pos) );
      }
   }

   consdata->ispropagated = FALSE;
   consdata->ispresolved = FALSE;

   return SCIP_OKAY;
}


/* merges entries with same linear variable into one entry and cleans up entries with coefficient 0.0 */
static
SCIP_RETCODE mergeAndCleanLinearVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real newcoef;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( consdata->linvarsmerged )
      return SCIP_OKAY;

   if( consdata->nlinvars == 0 )
   {
      consdata->linvarsmerged = TRUE;
      return SCIP_OKAY;
   }

   i = 0;
   while( i < consdata->nlinvars )
   {
      /* make sure linear variables are sorted (do this in every round, since we may move variables around) */
      consdataSortLinearVars(consdata);

      /* sum up coefficients that correspond to variable i */
      newcoef = consdata->lincoefs[i];
      for( j = i+1; j < consdata->nlinvars && consdata->linvars[i] == consdata->linvars[j]; ++j )
         newcoef += consdata->lincoefs[j];
      /* delete the additional variables in backward order */
      for( j = j-1; j > i; --j )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, j) );
      }

      /* delete also entry at position i, if it became zero (or was zero before) */
      if( SCIPisZero(scip, newcoef) )
      {
         SCIP_CALL( delLinearCoefPos(scip, cons, i) );
      }
      else
      {
         SCIP_CALL( chgLinearCoefPos(scip, cons, i, newcoef) );
         ++i;
      }
   }

   consdata->linvarsmerged = TRUE;

   return SCIP_OKAY;
}

/** removes fixes (or aggregated) linear variables from a nonlinear constraint */
static
SCIP_RETCODE removeFixedLinearVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinearconstraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real coef;
   SCIP_Real offset;
   SCIP_VAR* var;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);

   if( !consdata->isremovedfixingslin )
   {
      i = 0;
      while( i < consdata->nlinvars )
      {
         var = consdata->linvars[i];

         if( SCIPvarIsActive(var) )
         {
            ++i;
            continue;
         }

         coef = consdata->lincoefs[i];
         offset = 0.0;

         SCIP_CALL( SCIPvarGetProbvarSum(&var, &coef, &offset) );

         SCIPdebugMessage("  linear term %g*<%s> is replaced by %g * <%s> + %g\n", consdata->lincoefs[i], SCIPvarGetName(consdata->linvars[i]), coef, SCIPvarGetName(var), offset);

         /* delete previous variable (this will move another variable to position i) */
         SCIP_CALL( delLinearCoefPos(scip, cons, i) );

         /* put constant part into bounds */
         if( offset != 0.0 )
         {
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs -= offset;
            if( !SCIPisInfinity(scip,  consdata->rhs) )
               consdata->rhs -= offset;
         }

         /* nothing left to do if variable had been fixed */
         if( coef == 0.0 )
            continue;

         /* if GetProbvar gave a linear variable, just add it
          * if it's a multilinear variable, add it's disaggregated variables */
         if( SCIPvarIsActive(var) )
         {
            SCIP_CALL( addLinearCoef(scip, cons, var, coef) );
         }
         else
         {
            int        naggrs;
            SCIP_VAR** aggrvars;
            SCIP_Real* aggrscalars;
            SCIP_Real  aggrconstant;

            assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR);

            naggrs = SCIPvarGetMultaggrNVars(var);
            aggrvars = SCIPvarGetMultaggrVars(var);
            aggrscalars = SCIPvarGetMultaggrScalars(var);
            aggrconstant = SCIPvarGetMultaggrConstant(var);

            SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, consdata->nlinvars + naggrs) );

            for( j = 0; j < naggrs; ++j )
            {
               SCIP_CALL( addLinearCoef(scip, cons, aggrvars[j], coef * aggrscalars[j]) );
            }

            if( aggrconstant != 0.0 )
            {
               if( !SCIPisInfinity(scip, -consdata->lhs) )
                  consdata->lhs -= coef * aggrconstant;
               if( !SCIPisInfinity(scip,  consdata->rhs) )
                  consdata->rhs -= coef * aggrconstant;
            }
         }
      }

      SCIP_CALL( mergeAndCleanLinearVars(scip, cons) );

      consdata->isremovedfixingslin = TRUE;
   }

   SCIPdebugMessage("removed fixations of linear variables from <%s>\n  -> ", SCIPconsGetName(cons));
   SCIPdebug( SCIPprintCons(scip, cons, NULL) );

#ifndef NDEBUG
   for( i = 0; i < consdata->nlinvars; ++i )
      assert(SCIPvarIsActive(consdata->linvars[i]));
#endif

   return SCIP_OKAY;
}

/** removes fixed variables from expression graph */
static
SCIP_RETCODE removeFixedNonlinearVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   int nvars;
   int varssize;
   SCIP_Real constant;
   int i;
   int requsize;
   SCIPdebug( int j );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   if( conshdlrdata->isremovedfixings )
      return SCIP_OKAY;

   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, varssize) );

   i = 0;
   while( i < SCIPexprgraphGetNVars(conshdlrdata->exprgraph) )
   {
      var = (SCIP_VAR*)SCIPexprgraphGetVars(conshdlrdata->exprgraph)[i];
      if( SCIPvarIsActive(var) )
      {
         ++i;
         continue;
      }

      do
      {
         vars[0]  = var;
         coefs[0] = 1.0;
         constant = 0.0;
         nvars = 1;
         SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, coefs, &nvars, varssize, &constant, &requsize, TRUE) );

         if( requsize > varssize )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  requsize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, requsize) );
            varssize = requsize;
            continue;
         }

      } while( FALSE );

#ifdef SCIP_DEBUG
      SCIPdebugMessage("replace fixed variable <%s> by %g", SCIPvarGetName(var), constant);
      for( j = 0; j < nvars; ++j )
      {
         SCIPdebugPrintf(" %+g <%s>", coefs[j], SCIPvarGetName(vars[j]));
      }
      SCIPdebugPrintf("\n");
#endif

      SCIP_CALL( SCIPexprgraphReplaceVarByLinearSum(conshdlrdata->exprgraph, var, nvars, coefs, (void**)vars, constant) );

      i = 0;
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &coefs);

   conshdlrdata->isremovedfixings = TRUE;

   return SCIP_OKAY;
}

/** moves constant and linear part from expression graph node into constraint sides and linear part
 * frees expression graph node if expression is constant or linear */
static
SCIP_RETCODE splitOffLinearPart(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR** linvars;
   SCIP_Real* lincoefs;
   SCIP_Real constant;
   int linvarssize;
   int nlinvars;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->exprgraphnode == NULL )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   /* number of children of expression graph node is a good upper estimate on number of linear variables */
   linvarssize = SCIPexprgraphGetNodeNChildren(consdata->exprgraphnode);
   SCIP_CALL( SCIPallocBufferArray(scip, &linvars,  linvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, linvarssize) );

   /* get linear and constant part from expression graph node
    * releases expression graph node if not uses otherwise */
   SCIP_CALL( SCIPexprgraphNodeSplitOffLinear(conshdlrdata->exprgraph, &consdata->exprgraphnode, linvarssize, &nlinvars, (void**)linvars, lincoefs, &constant) );

   if( constant != 0.0 )
   {
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         consdata->lhs -= constant;
      if( !SCIPisInfinity(scip,  consdata->rhs) )
         consdata->rhs -= constant;
   }

   for( i = 0; i < nlinvars; ++i )
   {
      SCIP_CALL( addLinearCoef(scip, cons, linvars[i], lincoefs[i]) );
   }

   SCIPfreeBufferArray(scip, &linvars);
   SCIPfreeBufferArray(scip, &lincoefs);

   /* @todo linear variables that are also children of exprgraphnode could be moved into the expression graph for certain nonlinear operators (quadratic, polynomial), since that may allow better bound tightening */

   return SCIP_OKAY;
}

/** create a nonlinear row representation of the constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< nonlinear constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   if( consdata->nexprtrees == 0 )
   {
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            consdata->nlinvars, consdata->linvars, consdata->lincoefs,
            0, NULL, 0, NULL,
            NULL, consdata->lhs, consdata->rhs) );
   }
   else if( consdata->nexprtrees == 1 && consdata->nonlincoefs[0] == 1.0 )
   {
      assert(consdata->exprtrees[0] != NULL);
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            consdata->nlinvars, consdata->linvars, consdata->lincoefs,
            0, NULL, 0, NULL,
            consdata->exprtrees[0], consdata->lhs, consdata->rhs) );
   }
   else
   {
      /* since expression trees may share variable, we cannot easily sum them up,
       * but we can request a single expression tree from the expression graph
       */
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_EXPRTREE* exprtree;

      assert(consdata->exprgraphnode != NULL); /* since nexprtrees > 0 */
      conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
      assert(conshdlrdata != NULL);

      SCIP_CALL( SCIPexprgraphGetTree(conshdlrdata->exprgraph, consdata->exprgraphnode, &exprtree) );
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            consdata->nlinvars, consdata->linvars, consdata->lincoefs,
            0, NULL, 0, NULL,
            exprtree, consdata->lhs, consdata->rhs) );
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** tries to automatically convert a nonlinear constraint (or a part of it) into a more specific and more specialized constraint */
static
SCIP_RETCODE presolveUpgrade(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data structure */
   SCIP_CONS*            cons,               /**< source constraint to try to convert */
   SCIP_Bool*            upgraded,           /**< buffer to store whether constraint was upgraded */
   int*                  nupgdconss,         /**< buffer to increase if constraint was upgraded */
   int*                  naddconss           /**< buffer to increase with number of additional constraints created during upgrade */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS** upgdconss;
   int upgdconsssize;
   int nupgdconss_;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsModifiable(cons));
   assert(upgraded   != NULL);
   assert(nupgdconss != NULL);
   assert(naddconss  != NULL);

   *upgraded = FALSE;

   nupgdconss_ = 0;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* if there are no upgrade methods, we can stop */
   if( conshdlrdata->nnlconsupgrades == 0 )
      return SCIP_OKAY;

   /* set debug solution in expression graph and evaluate nodes, so reformulation methods can compute debug solution values for new auxiliary variables */
#ifdef SCIP_DEBUG_SOLUTION
   {
      SCIP_Real* varvals;

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(conshdlrdata->exprgraph)) );

      for( i = 0; i < SCIPexprgraphGetNVars(conshdlrdata->exprgraph); ++i )
         SCIP_CALL( SCIPdebugGetSolVal(scip, (SCIP_VAR*)SCIPexprgraphGetVars(conshdlrdata->exprgraph)[i], &varvals[i]) );

      SCIP_CALL( SCIPexprgraphEval(conshdlrdata->exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }
#endif

   upgdconsssize = 2;
   SCIP_CALL( SCIPallocBufferArray(scip, &upgdconss, upgdconsssize) );

   /* call the upgrading methods */

   SCIPdebugMessage("upgrading nonlinear constraint <%s> (up to %d upgrade methods):\n",
      SCIPconsGetName(cons), conshdlrdata->nnlconsupgrades);
   SCIPdebug( SCIPprintCons(scip, cons, NULL) );

   /* try all upgrading methods in priority order in case the upgrading step is enable  */
   for( i = 0; i < conshdlrdata->nnlconsupgrades; ++i )
   {
      if( !conshdlrdata->nlconsupgrades[i]->active )
         continue;
      if( conshdlrdata->nlconsupgrades[i]->nlconsupgd == NULL )
         continue;

      SCIP_CALL( conshdlrdata->nlconsupgrades[i]->nlconsupgd(scip, cons, &nupgdconss_, upgdconss, upgdconsssize) );

      while( nupgdconss_ < 0 )
      {
         /* upgrade function requires more memory: resize upgdconss and call again */
         assert(-nupgdconss_ > upgdconsssize);
         upgdconsssize = -nupgdconss_;
         SCIP_CALL( SCIPreallocBufferArray(scip, &upgdconss, -nupgdconss_) );

         SCIP_CALL( conshdlrdata->nlconsupgrades[i]->nlconsupgd(scip, cons, &nupgdconss_, upgdconss, upgdconsssize) );

         assert(nupgdconss_ != 0);
      }

      if( nupgdconss_ > 0 )
      {
         /* got upgrade */
         int j;

         SCIPdebugMessage(" -> upgraded to %d constraints:\n", nupgdconss_);

         /* add the upgraded constraints to the problem and forget them */
         for( j = 0; j < nupgdconss_; ++j )
         {
            SCIPdebugPrintf("\t");
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, upgdconss[j], NULL) ) );

            SCIP_CALL( SCIPaddCons(scip, upgdconss[j]) );      /*lint !e613*/
            SCIP_CALL( SCIPreleaseCons(scip, &upgdconss[j]) ); /*lint !e613*/
         }

         /* count the first upgrade constraint as constraint upgrade and the remaining ones as added constraints */
         *nupgdconss += 1;
         *naddconss += nupgdconss_ - 1;
         *upgraded = TRUE;

         /* delete upgraded constraint */
         SCIPdebugMessage("delete constraint <%s> after upgrade\n", SCIPconsGetName(cons));
         SCIP_CALL( SCIPdelCons(scip, cons) );

         break;
      }
   }

   SCIPfreeBufferArray(scip, &upgdconss);

   return SCIP_OKAY;
}

/** checks a nonlinear constraint for convexity and/or concavity */
static
SCIP_RETCODE checkCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool             expensivechecks,    /**< whether also expensive checks should be executed */
   SCIP_Bool             assumeconvex        /**< whether to assume convexity in inequalities */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_INTERVAL* varbounds;
   int varboundssize;
   SCIP_VAR* var;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->iscurvchecked )
      return SCIP_OKAY;

   SCIPdebugMessage("Checking curvature of constraint <%s>\n", SCIPconsGetName(cons));

   consdata->curvature = SCIP_EXPRCURV_LINEAR;
   consdata->iscurvchecked = TRUE;

   varbounds = NULL;
   varboundssize = 0;

   for( i = 0; i < consdata->nexprtrees; ++i )
   {
      assert(consdata->exprtrees[i] != NULL);
      assert(SCIPexprtreeGetNVars(consdata->exprtrees[i]) > 0 );

      if( assumeconvex )
      {
         /* for constraints a*f(x) <= rhs, we assume that it is convex */
         if( SCIPisInfinity(scip, -consdata->lhs) )
            consdata->curvatures[i] = SCIP_EXPRCURV_CONVEX;

         /* for constraints lhs <= a*f(x), we assume that it is concave */
         if( SCIPisInfinity(scip,  consdata->rhs) )
            consdata->curvatures[i] = SCIP_EXPRCURV_CONCAVE;
      }
      else
      {
         if( varboundssize == 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &varbounds, SCIPexprtreeGetNVars(consdata->exprtrees[i])) );
            varboundssize = SCIPexprtreeGetNVars(consdata->exprtrees[i]);
         }
         else if( varboundssize < SCIPexprtreeGetNVars(consdata->exprtrees[i]) )
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &varbounds, SCIPexprtreeGetNVars(consdata->exprtrees[i])) );
            varboundssize = SCIPexprtreeGetNVars(consdata->exprtrees[i]);
         }
         assert(varbounds != NULL);

         for( j = 0; j < SCIPexprtreeGetNVars(consdata->exprtrees[i]); ++j )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[i])[j];
            SCIPintervalSetBounds(&varbounds[j],
               -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))),    /*lint !e666*/
               +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var))) );  /*lint !e666*/
         }

         SCIP_CALL( SCIPexprtreeCheckCurvature(consdata->exprtrees[i], INTERVALINFTY, varbounds, &consdata->curvatures[i], NULL) );
         consdata->curvatures[i] = SCIPexprcurvMultiply(consdata->nonlincoefs[i], consdata->curvatures[i]);

         if( consdata->curvatures[i] == SCIP_EXPRCURV_UNKNOWN && SCIPconshdlrGetData(SCIPconsGetHdlr(cons))->isreformulated )
         {
            SCIPwarningMessage("indefinite expression tree in constraint <%s>\n", SCIPconsGetName(cons));
            SCIPdebug( SCIPexprtreePrintWithNames(consdata->exprtrees[i], NULL) );
            SCIPdebugPrintf("\n");
         }
      }

      /* @todo implement some more expensive checks */

      consdata->curvature = SCIPexprcurvAdd(consdata->curvature, consdata->curvatures[i]);

      SCIPdebugMessage("-> tree %d with coef %g is %s -> nonlinear part is %s\n", i, consdata->nonlincoefs[i], SCIPexprcurvGetName(consdata->curvatures[i]), SCIPexprcurvGetName(consdata->curvature));
   }

   SCIPfreeBufferArrayNull(scip, &varbounds);

   return SCIP_OKAY;
}  /*lint !e715*/

/* replaces a node by another node in expression graph
 * moves all parents of node to replacement
 * replaces all exprgraphnode's in constraints that are node by replacement
 * node may be freed, if captured only by given constraints
 */
static
SCIP_RETCODE reformReplaceNode(
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE**  node,               /**< pointer to node to be replaced in expression graph */
   SCIP_EXPRGRAPHNODE*   replacement,        /**< node which takes node's place */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int c;

   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(*node != NULL);
   assert(replacement != NULL);
   assert(conss != NULL || nconss == 0);

   SCIP_CALL( SCIPexprgraphMoveNodeParents(exprgraph, node, replacement) );

   /* if node still exists, then because it is captured by some constraint (it should not have parents anymore)
    * thus, look into the given constraints and replace their exprgraphnode by replacement
    * @todo may be expensive when this is done more often, esp. when we know that node will not be freed due to an added auxiliary constraint
    */
   assert(*node == NULL || SCIPexprgraphGetNodeNParents(*node) == 0);
   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( consdata->exprgraphnode == *node )
      {
         SCIP_CALL( SCIPexprgraphReleaseNode(exprgraph, &consdata->exprgraphnode) );
         consdata->exprgraphnode = replacement;
         SCIPexprgraphCaptureNode(replacement);
      }
   }
   *node = NULL;

   return SCIP_OKAY;
}

/** creates a new auxiliary variable and a new auxiliary nonlinear constraint connecting the var and a given node
 * node is replaced by new new auxiliary variables node in all parents of node in expression graph and in all constraints that use node
 */
static
SCIP_RETCODE reformNode2Var(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_CONS**           conss,              /**< constraints where to update exprgraphnode */
   int                   nconss,             /**< number of constraints */
   int*                  naddcons,           /**< counter to increase when constraint is added */
   SCIP_Bool             donotmultaggr       /**< whether to mark auxiliary variable as not to multiaggregate */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* auxvar;
   SCIP_CONS* auxcons;
   SCIP_EXPRGRAPHNODE* auxvarnode;
   SCIP_INTERVAL bounds;
   SCIP_Real minusone;
   SCIP_Bool cutoff;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(SCIPexprgraphGetNodeDepth(node) >= 1); /* do not turn vars or consts into new vars */

   bounds = SCIPexprgraphGetNodeBounds(node);
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);

   SCIPdebugMessage("add auxiliary variable and constraint %s for node %p(%d,%d)\n", name, (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));

   SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, SCIPintervalGetInf(bounds), SCIPintervalGetSup(bounds), 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, auxvar) );
   SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, &auxvarnode) );
#ifdef SCIP_DEBUG_SOLUTION
   /* store debug sol value of node as value for auxvar in debug solution and as value for auxvarnode */
   SCIPexprgraphSetVarNodeValue(auxvarnode, SCIPexprgraphGetNodeVal(node));
   SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(node)) );
#endif

   if( donotmultaggr )
   {
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, auxvar) );
   }

   /* set also bounds of auxvarnode to bounds, so it is available for new parent nodes (currently node->parents)
    * when updating their curvature information; avoid having to run domain propagation through exprgraph
    */
   SCIPexprgraphTightenNodeBounds(exprgraph, auxvarnode, bounds, BOUNDTIGHTENING_MINSTRENGTH, &cutoff);
   assert(!cutoff); /* we tightenend bounds from [-inf,+inf] to bounds, this should not be infeasible */

   /* add new constraint auxvar == node */
   minusone = -1.0;
   SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 1, &auxvar, &minusone, node, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
         FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddCons(scip, auxcons) );

   /* move parents of node in expression graph to auxvarnode
    * replace node by auxvarnode in constraints that use node */
   SCIP_CALL( reformReplaceNode(exprgraph, &node, auxvarnode, conss, nconss) );

   SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
   SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

   ++*naddcons;

   return SCIP_OKAY;
}

/** ensures that all children of a node have at least a given curvature by adding auxiliary variables */
static
SCIP_RETCODE reformEnsureChildrenMinCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   SCIP_EXPRGRAPHNODE*   node,               /**< expression graph node */
   SCIP_EXPRCURV         mincurv,            /**< minimal desired curvature */
   SCIP_CONS**           conss,              /**< constraints to check whether they use one of the replaced nodes */
   int                   nconss,             /**< number of constraints to check */
   int*                  naddcons            /**< counter to increase when constraint is added */
   )
{
   SCIP_EXPRGRAPHNODE* child;
   SCIP_Bool needupdate;

   int i;
   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(SCIPexprgraphGetNodeDepth(node) >= 1); /* do not turn vars or consts into new vars */
   assert(mincurv != SCIP_EXPRCURV_UNKNOWN); /* this is trivial and makes no sense */

   needupdate = FALSE; /* whether we need to update curvature of node */

   for( i = 0; i < SCIPexprgraphGetNodeNChildren(node); ++i )
   {
      child = SCIPexprgraphGetNodeChildren(node)[i];
      assert(child != NULL);

      if( (SCIPexprgraphGetNodeCurvature(child) & mincurv) != mincurv )
      {
         SCIPdebugMessage("add auxiliary variable for child %p(%d,%d) with curvature %s\n",
            (void*)child, SCIPexprgraphGetNodeDepth(child), SCIPexprgraphGetNodePosition(child), SCIPexprcurvGetName(SCIPexprgraphGetNodeCurvature(child)) );

         SCIP_CALL( reformNode2Var(scip, exprgraph, child, conss, nconss, naddcons, FALSE) );
         needupdate = TRUE;

         /* i'th child of node should now be a variable */
         assert(SCIPexprgraphGetNodeChildren(node)[i] != child);
         assert(SCIPexprgraphGetNodeOperator(SCIPexprgraphGetNodeChildren(node)[i]) == SCIP_EXPR_VARIDX);
      }

      assert(SCIPexprgraphGetNodeCurvature(SCIPexprgraphGetNodeChildren(node)[i]) & mincurv);
   }

   if( needupdate )
   {
      SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
      assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(node)));
   }

   return SCIP_OKAY;
}

/** reformulates a monomial by adding auxiliary variables and constraints for bilinear terms */
static
SCIP_RETCODE reformMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph */
   int                   nfactors,           /**< number of factors */
   SCIP_EXPRGRAPHNODE**  factors,            /**< factors */
   SCIP_Real*            exponents,          /**< exponents, or NULL if all 1.0 */
   SCIP_EXPRGRAPHNODE**  resultnode,         /**< buffer to store node which represents the reformulated monomial */
   SCIP_Bool             createauxcons,      /**< whether to create auxiliary var/cons */
   int*                  naddcons            /**< buffer to increase by number of added cons */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR* auxvar;
   SCIP_CONS* auxcons;
   SCIP_Real minusone;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(nfactors > 0);
   assert(factors != NULL);
   assert(resultnode != NULL);
   assert(naddcons != NULL);

   /* factors are just one node */
   if( nfactors == 1 && (exponents == NULL || exponents[0] == 1.0) )
   {
      *resultnode = factors[0];
      return SCIP_OKAY;
   }

   /* only one factor, but with exponent < 0.0 and factor has mixed sign, e.g., x^(-3)
    * reformulate as auxvar * factor^(-exponent) = 1 and return the node for auxvar in resultnode
    */
   if( nfactors == 1 && exponents[0] < 0.0 && SCIPexprgraphGetNodeBounds(factors[0]).inf < 0.0 && SCIPexprgraphGetNodeBounds(factors[0]).sup > 0.0 )  /*lint !e613*/
   {
      SCIP_EXPRGRAPHNODE* auxnode;
      SCIP_EXPRGRAPHNODE* reformfactors[2];
      SCIP_Real reformexp[2];

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);
      SCIPdebugMessage("add auxiliary variable and constraint %s\n", name);

      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );
      SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, resultnode) );

#ifdef SCIP_DEBUG_SOLUTION
      /* store debug sol value of node as value for auxvar in debug solution and as value for resultnode */
      {
         SCIP_Real debugval;
         debugval = pow(SCIPexprgraphGetNodeVal(factors[0]), exponents[0]);
         SCIPexprgraphSetVarNodeValue(*resultnode, debugval);
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
      }
#endif

      /* increase naddcons before next call to reformMonomial, to avoid duplicate variable names, which is bad for debugging */
      ++*naddcons;

      /* add reformulation for resultnode(=auxvar) * factor^(-exponent) = 1.0
       * if exponent != -1.0, then factor^(-exponent) should be moved into extra variable
       * finally one should get an EXPR_MUL node */
      reformfactors[0] = *resultnode;
      reformfactors[1] = factors[0];
      reformexp[0] = 1.0;
      reformexp[1] = -exponents[0];  /*lint !e613*/
      SCIP_CALL( reformMonomial(scip, exprgraph, 2, reformfactors, reformexp, &auxnode, FALSE, naddcons) );

      SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 0, NULL, NULL, auxnode, 1.0, 1.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, auxcons) );

      SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

      return SCIP_OKAY;
   }

   /* only one factor, but with exponent != 1.0 */
   if( nfactors == 1 )
   {
      /* create some power expression node, if not existing already */
      SCIP_EXPRGRAPHNODE* expnode;
      SCIP_EXPRGRAPHNODE* parent;
      int p;

      assert(exponents != NULL);

      /* check if there is already a node for factors[0]^exponents[0] */
      expnode = NULL;
      for( p = 0; p < SCIPexprgraphGetNodeNParents(factors[0]); ++p)
      {
         parent = SCIPexprgraphGetNodeParents(factors[0])[p];
         if( SCIPisIntegral(scip, exponents[0]) &&
            SCIPexprgraphGetNodeOperator(parent) == SCIP_EXPR_INTPOWER &&
            SCIPexprgraphGetNodeIntPowerExponent(parent) == (int)SCIPround(scip, exponents[0]) )
         {
            expnode = parent;
            break;
         }
         if( SCIPexprgraphGetNodeOperator(parent) == SCIP_EXPR_REALPOWER &&
            SCIPisEQ(scip, SCIPexprgraphGetNodeRealPowerExponent(parent), exponents[0]) )
         {
            expnode = parent;
         }
      }
      if( expnode == NULL )
      {
         if( SCIPisIntegral(scip, exponents[0]) )
            SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &expnode, SCIP_EXPR_INTPOWER, (int)SCIPround(scip, exponents[0])) );
         else
            SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &expnode, SCIP_EXPR_REALPOWER, exponents[0]) );

         SCIP_CALL( SCIPexprgraphAddNode(exprgraph, expnode, -1, 1, &factors[0]) );
         SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(expnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
         assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(expnode)));
      }

      if( createauxcons )
      {
         /* @todo if there was already a node for factors[0]^exponents[0], then there may have been also been already an auxiliary variable and constraint (-> ex7_3_4) */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);
         SCIPdebugMessage("add auxiliary variable and constraint %s\n", name);

         SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, auxvar) );
         SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, resultnode) );

#ifdef SCIP_DEBUG_SOLUTION
         SCIPexprgraphSetVarNodeValue(*resultnode, SCIPexprgraphGetNodeVal(expnode));
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(expnode)) );
#endif

         /* add new constraint resultnode(=auxvar) = expnode */
         minusone = -1.0;
         SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 1, &auxvar, &minusone, expnode, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, auxcons) );

         SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

         ++*naddcons;
      }
      else
      {
         *resultnode = expnode;
      }

      return SCIP_OKAY;
   }

   if( nfactors == 2 && exponents != NULL && exponents[0] != 1.0 && exponents[0] == exponents[1] )  /*lint !e777*/
   {
      /* factor0^exponent * factor1^exponent with exponent != 1.0, reform as (factor0*factor1)^exponent */
      SCIP_EXPRGRAPHNODE* productnode;

      /* create node for factor0*factor1 */
      SCIP_CALL( reformMonomial(scip, exprgraph, 2, factors, NULL, &productnode, TRUE, naddcons) );

      /* create node for productnode^exponents[0] by just calling this method again */
      SCIP_CALL( reformMonomial(scip, exprgraph, 1, &productnode, &exponents[0], resultnode, createauxcons, naddcons) );

      return SCIP_OKAY;
   }
   /* @todo if nfactors > 2, assemble groups of factors with same exponent and replace these by a single variable first */

   {
      /* at least two factors */
      /* create auxvar for left half (recursively) and auxvar for right half (recursively) and maybe new auxvar for product */
      /* @todo it may be enough to replace single factors in a monomial to get it convex or concave, see Westerlund et.al. */

      SCIP_EXPRGRAPHNODE* productnode;
      SCIP_EXPRGRAPHNODE* leftright[2]; /* {left, right} */
      SCIP_EXPRGRAPHNODE* parent;
      int half;
      int p;

      half = nfactors / 2;
      assert(half > 0);
      assert(half < nfactors);

      SCIP_CALL( reformMonomial(scip, exprgraph, half, factors, exponents, &leftright[0], TRUE, naddcons) );
      SCIP_CALL( reformMonomial(scip, exprgraph, nfactors-half, &factors[half], exponents != NULL ? &exponents[half] : NULL, &leftright[1], TRUE, naddcons) );  /*lint !e826*/

      /* check if there is already a node for left * right */
      productnode = NULL;
      for( p = 0; p < SCIPexprgraphGetNodeNParents(leftright[0]); ++p)
      {
         parent = SCIPexprgraphGetNodeParents(factors[0])[p];
         if( SCIPexprgraphGetNodeOperator(parent) != SCIP_EXPR_MUL )
            continue;

         assert(SCIPexprgraphGetNodeNChildren(parent) == 2);
         if( (SCIPexprgraphGetNodeChildren(parent)[0] == leftright[0] && SCIPexprgraphGetNodeChildren(parent)[1] == leftright[1]) ||
            ( SCIPexprgraphGetNodeChildren(parent)[0] == leftright[1] && SCIPexprgraphGetNodeChildren(parent)[1] == leftright[0]) )
         {
            productnode = parent;
            break;
         }
      }
      if( productnode == NULL )
      {
         /* create node for left * right */
         SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &productnode, SCIP_EXPR_MUL, NULL) );
         SCIP_CALL( SCIPexprgraphAddNode(exprgraph, productnode, -1, 2, leftright) );
         SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(productnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
         assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(productnode)));
      }

      if( createauxcons )
      {
         /* @todo if there was already a node for factors[0]^exponents[0], then there may have been also been already an auxiliary variable and constraint (-> ex7_3_4) */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);
         SCIPdebugMessage("add auxiliary variable and constraint %s\n", name);

         SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, auxvar) );
         SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, resultnode) );

#ifdef SCIP_DEBUG_SOLUTION
         SCIPexprgraphSetVarNodeValue(*resultnode, SCIPexprgraphGetNodeVal(productnode));
         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(productnode)) );
#endif

         /* add new constraint resultnode(= auxvar) == left * right */
         minusone = -1.0;
         SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 1, &auxvar, &minusone, productnode, 0.0, 0.0, TRUE, TRUE, TRUE, TRUE, TRUE,
               FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, auxcons) );

         SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
         SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

         ++*naddcons;
      }
      else
      {
         *resultnode = productnode;
      }
   }

   return SCIP_OKAY;
}

/** reformulates expression graph into a form so that for each node under- and overestimators could be computed
 * similar to factorable reformulation in other global solvers, but sometimes does not split up complex operands (like quadratic)
 */
static
SCIP_RETCODE reformulate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  naddcons            /**< buffer to increase by the number of added constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   SCIP_EXPRGRAPHNODE** children;
   SCIP_EXPRGRAPHNODE* reformnode;
   SCIP_Bool havenonlinparent;
   SCIP_Bool domainerror;
   int nchildren;
   int c;
   int d;
   int i;
   int u;
#ifndef NDEBUG
   int j;
#endif

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(naddcons != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
   assert(!SCIPinProbing(scip));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->isreformulated )
   {
      SCIPdebugMessage("skip reformulation, already done\n");
      return SCIP_OKAY;
   }

   exprgraph = conshdlrdata->exprgraph;

   /* make sure current variable bounds are variable nodes of exprgraph */
   SCIP_CALL( SCIPexprgraphPropagateVarBounds(exprgraph, INTERVALINFTY, FALSE, &domainerror) );
   assert(!domainerror); /* should have been found by domain propagation */

   /* set debug solution in expression graph and evaluate nodes, so we can compute debug solution values for auxiliary variables */
#ifdef SCIP_DEBUG_SOLUTION
   {
      SCIP_Real* varvals;

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(exprgraph)) );

      for( i = 0; i < SCIPexprgraphGetNVars(exprgraph); ++i )
         SCIP_CALL( SCIPdebugGetSolVal(scip, (SCIP_VAR*)SCIPexprgraphGetVars(exprgraph)[i], &varvals[i]) );

      SCIP_CALL( SCIPexprgraphEval(exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }
#endif

   for( d = 1; d < SCIPexprgraphGetDepth(exprgraph); ++d )
   {
      i = 0;
      while( i < SCIPexprgraphGetNNodes(exprgraph)[d] )
      {
         node = SCIPexprgraphGetNodes(exprgraph)[d][i];
         assert(node != NULL);

         /* skip disabled nodes, they should be removed soon */
         if( !SCIPexprgraphIsNodeEnabled(node) )
         {
            ++i;
            continue;
         }

         /* make sure bounds and curvature of node are uptodate */
         SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
         assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(node)));

         /* try external reformulation methods */
         for( u = 0; u < conshdlrdata->nnlconsupgrades; ++u )
         {
            if( conshdlrdata->nlconsupgrades[u]->nodereform != NULL && conshdlrdata->nlconsupgrades[u]->active )
            {
               SCIP_CALL( conshdlrdata->nlconsupgrades[u]->nodereform(scip, exprgraph, node, naddcons, &reformnode) );
               if( reformnode == NULL )
                  continue;

               SCIPdebugMessage("external nodereform reformulated node %p(%d,%d), replace by %p\n",
                  (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node), (void*)reformnode);

               SCIP_CALL( reformReplaceNode(exprgraph, &node, reformnode, conss, nconss) );
               SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(reformnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
               assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(reformnode)));

               break;
            }
         }
         /* if node has been reformulated, continue with next node without increasing i */
         if( u < conshdlrdata->nnlconsupgrades )
            continue;

         /* leave nodes that are known to be convex/concave/linear as they are */
         if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
         {
            SCIPdebugMessage("skip reformulating node %p(%d,%d) = ", (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));
            SCIPdebug( SCIPexprgraphPrintNode(node, NULL) );
            SCIPdebugPrintf(", curv = %s\n", SCIPexprcurvGetName(SCIPexprgraphGetNodeCurvature(node)));
            ++i;
            continue;
         }

         /* get flag whether node has a nonlinear parent
          * we want to know whether the current node will be at the top of the tree after the next simplification run
          * due to the tricky reformulation of polynomials below, this may not be the case yet
          */
         havenonlinparent = SCIPexprgraphHasNodeNonlinearAncestor(node);

         /* take action */
         assert(SCIPexprgraphGetNodeCurvature(node) == SCIP_EXPRCURV_UNKNOWN);
         SCIPdebugMessage("think about reformulating %s node %p(%d,%d) = ", SCIPexpropGetName(SCIPexprgraphGetNodeOperator(node)), (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));
         SCIPdebug( SCIPexprgraphPrintNode(node, NULL) );
         SCIPdebugPrintf("\n");

         children  = SCIPexprgraphGetNodeChildren(node);
         nchildren = SCIPexprgraphGetNodeNChildren(node);
         assert(children != NULL || nchildren == 0);

#ifndef NDEBUG
         /* at this place, all children nodes should have a known curvature, except if they only appear only linearly in constraints
          * the latter for cases where constraints with unknown curvature are upgraded to other constraint handler that can handle these (quadratic, signpower,...)
          */
         for( j = 0; j < nchildren; ++j )
         {
            assert(children[j] != NULL);  /*lint !e613*/
            if( havenonlinparent ||
               (  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_PLUS  &&
                  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_MINUS &&
                  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_SUM   &&
                  SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_LINEAR) )
               assert(SCIPexprgraphGetNodeCurvature(children[j]) != SCIP_EXPRCURV_UNKNOWN);  /*lint !e613*/
         }
#endif

         switch( SCIPexprgraphGetNodeOperator(node) )
         {
         case SCIP_EXPR_VARIDX:
         case SCIP_EXPR_PARAM:
         case SCIP_EXPR_CONST:
            SCIPerrorMessage("node with operator %d cannot have unknown curvature\n", SCIPexprgraphGetNodeOperator(node));
            SCIPABORT();
            break;

            /* linear operands */
         case SCIP_EXPR_PLUS:
         case SCIP_EXPR_MINUS:
         case SCIP_EXPR_SUM:
         case SCIP_EXPR_LINEAR:
         {
            /* children have conflicting curvature, we can handle such sums in cons_nonlinear
             * thus, turn node into variable, if it has nonlinear parents */
            if( havenonlinparent )
            {
               SCIP_CALL( reformNode2Var(scip, exprgraph, node, conss, nconss, naddcons, FALSE) );
               assert(node != NULL);
               assert(SCIPexprgraphGetNodeNParents(node) == 0); /* node should now be at top of graph */
            }
            ++i;
            break;
         }

         /* quadratic operands */
         case SCIP_EXPR_MUL:
         case SCIP_EXPR_QUADRATIC:
         {
            /* ensure all children are linear, so next simplifier run makes sure all children will be variables */
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_LINEAR, conss, nconss, naddcons) );
            if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
            {
               /* if curvature is now known then we are done */
               ++i;
               break;
            }

            /* if we have nonlinear parents or a sibling, then add add auxiliary variable for this node, so an upgrade to cons_quadratic should take place
             * we assume that siblings are non-linear and non-quadratic, which should be the case if simplifier was run, and also if this node was created during reformulating a polynomial
             * @todo we could also add auxvars for the sibling nodes, e.g., if there is only one
             * @todo if sibling nodes are quadratic (or even linear) due to reformulation, then we do not need to reform here... (-> nvs16)
             *       maybe this step should not be done here at all if havenonlinparent is FALSE? e.g., move into upgrade from quadratic?
             */
            if( havenonlinparent || SCIPexprgraphHasNodeSibling(node) )
            {
               SCIP_CALL( reformNode2Var(scip, exprgraph, node, conss, nconss, naddcons, FALSE) );
               assert(node != NULL);
               assert(SCIPexprgraphGetNodeNParents(node) == 0); /* node should now be at top of graph, so it can be upgraded by cons_quadratic */
               break;
            }

            ++i;
            break;
         }

         case SCIP_EXPR_DIV:
         {
            /* reformulate as bilinear term
             * note that in the reformulation, a zero in the denominator forces the nominator to be zero too, but the auxiliary can be arbitrary
             */
            SCIP_EXPRGRAPHNODE* auxvarnode;
            SCIP_EXPRGRAPHNODE* auxnode;
            SCIP_EXPRGRAPHNODE* auxchildren[3];
            SCIP_Real           lincoefs[3];
            SCIP_QUADELEM       quadelem;
            SCIP_VAR*           auxvar;
            SCIP_CONS*          auxcons;
            char                name[SCIP_MAXSTRLEN];
            SCIP_INTERVAL       bounds;

            bounds = SCIPexprgraphGetNodeBounds(node);
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%d", *naddcons);

            SCIPdebugMessage("add auxiliary variable %s for division in node %p(%d,%d)\n", name, (void*)node, SCIPexprgraphGetNodeDepth(node), SCIPexprgraphGetNodePosition(node));

            SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, SCIPintervalGetInf(bounds), SCIPintervalGetSup(bounds), 0.0,
                  SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
            SCIP_CALL( SCIPaddVar(scip, auxvar) );
            SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&auxvar, &auxvarnode) );

#ifdef SCIP_DEBUG_SOLUTION
            {
               SCIP_Real debugval;
               debugval = SCIPexprgraphGetNodeVal(children[0]) / SCIPexprgraphGetNodeVal(children[1]);
               SCIPexprgraphSetVarNodeValue(auxvarnode, debugval);
               SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
            }
#endif

            /* add new constraint auxvar * child[1] - child[0] == 0 */
            auxchildren[0] = children[0];  /*lint !e613*/
            auxchildren[1] = children[1];  /*lint !e613*/
            auxchildren[2] = auxvarnode;

            lincoefs[0] = -1.0;
            lincoefs[1] =  0.0;
            lincoefs[2] =  0.0;

            quadelem.idx1 = 1;
            quadelem.idx2 = 2;
            quadelem.coef = 1.0;

            SCIP_CALL( SCIPexprgraphCreateNodeQuadratic(SCIPblkmem(scip), &auxnode, 3, lincoefs, 1, &quadelem, 0.0) );
            SCIP_CALL( SCIPexprgraphAddNode(exprgraph, auxnode, -1, 3, auxchildren) );

            SCIP_CALL( SCIPcreateConsNonlinear2(scip, &auxcons, name, 0, NULL, NULL, auxnode, 0.0, 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, auxcons) );

            /* replace node by auxvarnode in graph and constraints that use it */
            SCIP_CALL( reformReplaceNode(exprgraph, &node, auxvarnode, conss, nconss) );

            SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
            SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

            ++*naddcons;

            /* do not increase i, since node was removed and not necessarily replaced here */
            break;
         }

         case SCIP_EXPR_MIN:
         {
            /* make sure that both children are concave, because min of concave functions is concave */
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_CONCAVE, conss, nconss, naddcons) );
            assert(SCIPexprgraphGetNodeCurvature(node) & SCIP_EXPRCURV_CONCAVE);
            ++i;
            break;
         }

         case SCIP_EXPR_MAX:
         {
            /* make sure that both children are convex, because max of convex functions is convex */
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_CONVEX, conss, nconss, naddcons) );
            assert(SCIPexprgraphGetNodeCurvature(node) & SCIP_EXPRCURV_CONVEX);
            ++i;
            break;

         }

         case SCIP_EXPR_INTPOWER:
         {
            assert(nchildren == 1);

            /* for an intpower with mixed sign in the base and negative exponent, we reformulate similar as for EXPR_DIV */
            if( SCIPexprgraphGetNodeIntPowerExponent(node) < 0 && SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(children[0])) < 0.0 && SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(children[0])) > 0.0 )  /*lint !e613*/
            {
               SCIP_EXPRGRAPHNODE* auxvarnode;
               SCIP_Real exponent;

               /* if we have something like x^(-3) with mixed sign for x, then add auxvar and reform as auxvar*x^3 = 1 via reformMonomial */
               exponent = (SCIP_Real)SCIPexprgraphGetNodeIntPowerExponent(node);
               SCIP_CALL( reformMonomial(scip, exprgraph, 1, children, &exponent, &auxvarnode, TRUE, naddcons) );
               /* replace node by auxvarnode */
               SCIP_CALL( reformReplaceNode(exprgraph, &node, auxvarnode, conss, nconss) );
               break;
            }

            /* otherwise, we continue as for other univariate operands */
         }   /*lint -fallthrough*/

         /* univariate operands where the child does not have bounds and curvature from which we can deduce curvature of this node,
          * but where we can do more if the child is linear
          * thus, turn child into auxiliary variable
          */
         case SCIP_EXPR_SQUARE:
         case SCIP_EXPR_SQRT:
         case SCIP_EXPR_EXP:
         case SCIP_EXPR_LOG:
         case SCIP_EXPR_ABS:
         case SCIP_EXPR_REALPOWER:
         case SCIP_EXPR_SIGNPOWER:
         {
            assert(nchildren == 1);

            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_LINEAR, conss, nconss, naddcons) );

            if( SCIPexprgraphGetNodeCurvature(node) == SCIP_EXPRCURV_UNKNOWN )
            {
               /* the only case where f(x) for a linear term x is indefinite here is if f is intpower or signpower and x has mixed sign */
               assert(SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_INTPOWER || SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_SIGNPOWER);
               assert(SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(children[0])) < 0.0);  /*lint !e613*/
               assert(SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(children[0])) > 0.0);  /*lint !e613*/
            }

            /* update curvature of node */
            SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
            assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(node)));

            if( SCIPexprgraphGetNodeCurvature(node) == SCIP_EXPRCURV_UNKNOWN )
            {
               /* if intpower and signpower with positive exponent and a mixed sign in the child bounds still does not give a curvature,
                * we can do more if we make this node the root of a nonlinear constraints expression node, so it can be upgraded by cons_signpower
                * of course, this is only required if the node is still intermediate
                *
                * an intpower with negative exponent should have been handled above
                * for signpower, we assume the exponent is > 1.0
                */
               assert(SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_INTPOWER || SCIPexprgraphGetNodeOperator(node) == SCIP_EXPR_SIGNPOWER);
               assert(SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_INTPOWER  || SCIPexprgraphGetNodeIntPowerExponent(node) > 1);
               assert(SCIPexprgraphGetNodeOperator(node) != SCIP_EXPR_SIGNPOWER || SCIPexprgraphGetNodeSignPowerExponent(node) > 1.0);
               if( havenonlinparent )
               {
                  SCIP_CALL( reformNode2Var(scip, exprgraph, node, conss, nconss, naddcons, FALSE) );
                  assert(node != NULL); /* it should be used by some auxiliary constraint now */
                  assert(SCIPexprgraphGetNodeNParents(node) == 0); /* node should now be at top of graph (and used by new auxiliary constraint) */
               }
            }
            ++i;

            break;
         }

         case SCIP_EXPR_SIN:
         case SCIP_EXPR_COS:
         case SCIP_EXPR_TAN:
         case SCIP_EXPR_SIGN:
            /* case SCIP_EXPR_ERF   : */
            /* case SCIP_EXPR_ERFI  : */
         {
            SCIPerrorMessage("no support for trigonometric or sign operands yet\n");
            return SCIP_ERROR;
         }

         case SCIP_EXPR_PRODUCT:
         {
            /* ensure all children are linear */
            SCIP_CALL( reformEnsureChildrenMinCurvature(scip, exprgraph, node, SCIP_EXPRCURV_LINEAR, conss, nconss, naddcons) );
            if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
            {
               ++i;
               break;
            }

            /* if curvature is still unknown (quite likely), then turn into a cascade of bilinear terms
             * if node has parents, then ensure that it has a known curvature, otherwise we are also fine with a node that is a product of two (aux)variables */
            SCIP_CALL( reformMonomial(scip, exprgraph, nchildren, children, NULL, &reformnode, havenonlinparent, naddcons) );

            /* replace node by reformnode in graph and in all constraints that use it */
            SCIP_CALL( reformReplaceNode(exprgraph, &node, reformnode, conss, nconss) );

            /* do not increase i, since node was removed and not necessarily replaced here */
            break;
         }
         case SCIP_EXPR_POLYNOMIAL:
         {
            /* if polynomial has several monomials, replace by a sum of nodes each having a single monomial and one that has all linear and quadratic monomials
             * if polynomial has only a single monomial, then reformulate that one
             */
            SCIP_EXPRDATA_MONOMIAL** monomials;
            SCIP_EXPRDATA_MONOMIAL* monomial;
            int nmonomials;
            SCIP_Real* exponents;
            SCIP_Real coef;
            int* childidxs;
            int nfactors;
            int f;
            SCIP_INTERVAL childbounds;
            SCIP_EXPRCURV childcurv;
            SCIP_Bool modified;

            monomials  = SCIPexprgraphGetNodePolynomialMonomials(node);
            nmonomials = SCIPexprgraphGetNodePolynomialNMonomials(node);
            assert(nmonomials >= 1); /* constant polynomials should have been simplified away */

            if( nmonomials > 1 )
            {
               SCIP_EXPRGRAPHNODE* sumnode;
               SCIP_Real constant;
               int nquadelems;
               SCIP_QUADELEM* quadelems;
               SCIP_Real* lincoefs;
               int nmonomialnodes;
               SCIP_EXPRGRAPHNODE** childrennew;
               SCIP_EXPRGRAPHNODE** monomialnodes;
               int m;

               /* @todo if a monomial is a factor of another monomial, then we could (and should?) replace it there by the node we create for it here -> ex7_2_1
                * @todo factorizing the polynomial could be beneficial
                */

               /* constant part of polynomials, to add to first monomialnode, if any, or quadratic or linear part */
               constant = SCIPexprgraphGetNodePolynomialConstant(node);

               /* coefficients from linear monomials */
               lincoefs = NULL;

               /* quadratic elements */
               nquadelems = 0;
               quadelems = NULL;

               /* expression graph nodes representing single higher-degree monomials, and single node with linear and/or quadratic monomials */
               nmonomialnodes = 0;
               SCIP_CALL( SCIPallocBufferArray(scip, &monomialnodes, nmonomials) );

               /* children of new monomial nodes that are setup */
               childrennew = NULL;

               for( m = 0; m < nmonomials; ++m )
               {
                  monomial = monomials[m];
                  assert(monomial != NULL);

                  coef = SCIPexprGetMonomialCoef(monomial);
                  exponents = SCIPexprGetMonomialExponents(monomial);
                  childidxs = SCIPexprGetMonomialChildIndices(monomial);
                  nfactors = SCIPexprGetMonomialNFactors(monomial);
                  assert(nfactors >= 1); /* constant monomials should have been simplified away */
                  assert(coef != 0.0);  /* zero-monomials should have been simplified away */

                  if( nfactors == 1 && exponents[0] == 1.0 )
                  {
                     /* linear monomial */
                     if( lincoefs == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, nchildren) );
                        BMSclearMemoryArray(lincoefs, nchildren);
                     }
                     assert(0 <= childidxs[0] && childidxs[0] < nchildren);
                     assert(lincoefs[childidxs[0]] == 0.0); /* monomials should have been merged */
                     lincoefs[childidxs[0]] = coef;
                  }
                  else if( nfactors == 1 && exponents[0] == 2.0 )
                  {
                     /* square monomial */
                     if( quadelems == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nmonomials) );
                     }
                     quadelems[nquadelems].idx1 = childidxs[0];
                     quadelems[nquadelems].idx2 = childidxs[0];
                     quadelems[nquadelems].coef = coef;
                     ++nquadelems;
                  }
                  else if( nfactors == 2 && exponents[0] == 1.0 && exponents[1] == 1.0 )
                  {
                     /* bilinear monomial */
                     if( quadelems == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &quadelems, nmonomials) );
                     }
                     assert(childidxs[0] < childidxs[1]);
                     quadelems[nquadelems].idx1 = childidxs[0];
                     quadelems[nquadelems].idx2 = childidxs[1];
                     quadelems[nquadelems].coef = coef;
                     ++nquadelems;
                  }
                  else
                  {
                     /* general monomial -> pass into separate expression graph node */
                     SCIP_EXPRDATA_MONOMIAL* monomialnew;

                     /* create new node for this monomial, children will be those associated with factors */
                     SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomialnew, coef, nfactors, NULL, exponents) );
                     SCIP_CALL( SCIPexprgraphCreateNodePolynomial(SCIPblkmem(scip), &monomialnodes[nmonomialnodes], 1, &monomialnew, constant, FALSE) );
                     constant = 0.0;

                     if( childrennew == NULL )
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &childrennew, nchildren) );
                     }
                     assert(nfactors <= nchildren);
                     for( f = 0; f < nfactors; ++f )
                        childrennew[f] = children[childidxs[f]];  /*lint !e613*/

                     /* add new node to same depth as this node, so we will reformulate it during this run
                      * no need to refresh bounds/curvature here, since that will be done when we reach this node next */
                     SCIP_CALL( SCIPexprgraphAddNode(exprgraph, monomialnodes[nmonomialnodes], SCIPexprgraphGetNodeDepth(node), nfactors, childrennew) );

                     ++nmonomialnodes;
                  }
               }
               /* should have had at least one linear, quadratic, or general monomial */
               assert(lincoefs != NULL || nquadelems > 0 || nmonomialnodes > 0);

               if( nquadelems > 0 )
               {
                  /* create and add additional node for quadratic and linear part, simplifier should take care of removing unused children later */
                  SCIP_CALL( SCIPexprgraphCreateNodeQuadratic(SCIPblkmem(scip), &monomialnodes[nmonomialnodes], nchildren, lincoefs, nquadelems, quadelems, constant) );
                  constant = 0.0;
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, monomialnodes[nmonomialnodes], SCIPexprgraphGetNodeDepth(node), nchildren, children) );
                  ++nmonomialnodes;
               }
               else if( lincoefs != NULL )
               {
                  /* create additional node for linear part, simplifier should take care of removing unused children later */
                  SCIP_CALL( SCIPexprgraphCreateNodeLinear(SCIPblkmem(scip), &monomialnodes[nmonomialnodes], nchildren, lincoefs, constant) );
                  constant = 0.0;
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, monomialnodes[nmonomialnodes], SCIPexprgraphGetNodeDepth(node), nchildren, children) );
                  ++nmonomialnodes;
               }
               assert(constant == 0.0); /* the constant should have been used somewhere */

               SCIPfreeBufferArrayNull(scip, &lincoefs);
               SCIPfreeBufferArrayNull(scip, &quadelems);
               SCIPfreeBufferArrayNull(scip, &childrennew);

               assert(nmonomialnodes > 0);
               if( nmonomialnodes > 1 )
               {
                  /* add node for sum of monomials to expression graph */
                  SCIP_CALL( SCIPexprgraphCreateNode(SCIPblkmem(scip), &sumnode, nmonomialnodes == 2 ? SCIP_EXPR_PLUS : SCIP_EXPR_SUM) );
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, sumnode, -1, nmonomialnodes, monomialnodes) );
               }
               else
               {
                  /* if only one monomial, then because polynomial was linear or quadratic... */
                  assert(SCIPexprgraphGetNodeOperator(monomialnodes[0]) == SCIP_EXPR_LINEAR || SCIPexprgraphGetNodeOperator(monomialnodes[0]) == SCIP_EXPR_QUADRATIC);
                  sumnode = monomialnodes[0];
               }
               SCIPfreeBufferArray(scip, &monomialnodes);

               /* replace node by sumnode, and we are done */
               SCIP_CALL( reformReplaceNode(exprgraph, &node, sumnode, conss, nconss) );

               SCIPdebugMessage("splitup polynomial into sum of %d nodes\n", nmonomialnodes);

               break;
            }

            /* reformulate a monomial such that it becomes convex or concave, if necessary */

            monomial = monomials[0];
            assert(monomial != NULL);

            coef = SCIPexprGetMonomialCoef(monomial);
            exponents = SCIPexprGetMonomialExponents(monomial);
            childidxs = SCIPexprGetMonomialChildIndices(monomial);
            nfactors = SCIPexprGetMonomialNFactors(monomial);
            assert(nfactors >= 1); /* constant monomials should have been simplified away */
            assert(coef != 0.0);  /* zero-monomials should have been simplified away */

            /* check if we make monomial convex or concave by making a child linear */
            modified = FALSE;
            if( nfactors == 1 )
            {
               /* ensure that the child of an univariate monomial is linear if its current (bounds,curvature) yields an unknown curvature for the monomial
                * and with linear child it had a known curvature (rules out x^a, a negative, x not linear) */
               childcurv = SCIPexprgraphGetNodeCurvature(children[childidxs[0]]);  /*lint !e613*/
               childbounds = SCIPexprgraphGetNodeBounds(children[childidxs[0]]);  /*lint !e613*/
               assert(SCIPexprcurvPower(childbounds, childcurv, exponents[0]) == SCIP_EXPRCURV_UNKNOWN); /* this is exactly the curvature of the node, which is unknown */

               /* if monomial were convex or concave if child were linear, then make child linear */
               if( SCIPexprcurvPower(childbounds, SCIP_EXPRCURV_LINEAR, exponents[0]) != SCIP_EXPRCURV_UNKNOWN )
               {
                  assert(childcurv != SCIP_EXPRCURV_LINEAR);
                  SCIPdebugMessage("reform child %d (univar. monomial) with curv %s into var\n", childidxs[0], SCIPexprcurvGetName(childcurv));
                  SCIP_CALL( reformNode2Var(scip, exprgraph, children[childidxs[0]], conss, nconss, naddcons, FALSE) );  /*lint !e613*/
                  modified = TRUE;
               }
            }
            else
            {
               /* check if the conditions on the exponents allow for a convex or concave monomial, assuming that the children are linear
                * if one of these conditions is fulfilled but a children curvature does not fit, then make these children linear
                */
               int nnegative;
               int npositive;
               SCIP_Real sum;
               SCIP_Bool expcurvpos;
               SCIP_Bool expcurvneg;
               SCIP_EXPRCURV desiredcurv;

               nnegative = 0; /* number of negative exponents */
               npositive = 0; /* number of positive exponents */
               sum = 0.0;     /* sum of exponents */
               expcurvpos = TRUE; /* whether exp_j * f_j''(x) >= 0 for all factors (assuming f_j >= 0) */
               expcurvneg = TRUE; /* whether exp_j * f_j''(x) <= 0 for all factors (assuming f_j >= 0) */

               for( f = 0; f < nfactors; ++f )
               {
                  childcurv = SCIPexprgraphGetNodeCurvature(children[childidxs[f]]);  /*lint !e613*/
                  assert(childcurv != SCIP_EXPRCURV_UNKNOWN);
                  childbounds = SCIPexprgraphGetNodeBounds(children[childidxs[f]]);  /*lint !e613*/
                  if( childbounds.inf < 0.0 && childbounds.sup > 0.0 )
                     break;

                  if( exponents[f] < 0.0 )
                     ++nnegative;
                  else
                     ++npositive;
                  sum += exponents[f];

                  /* negate curvature if factor is negative */
                  if( childbounds.inf < 0.0 )
                     childcurv = SCIPexprcurvNegate(childcurv);

                  /* check if exp_j * checkcurv is convex (>= 0) and/or concave */
                  childcurv = SCIPexprcurvMultiply(exponents[f], childcurv);
                  if( !(childcurv & SCIP_EXPRCURV_CONVEX) )
                     expcurvpos = FALSE;
                  if( !(childcurv & SCIP_EXPRCURV_CONCAVE) )
                     expcurvneg = FALSE;
               }

               /* if some child can be both positive and negative, then nothing we can do here to get the monomial convex or concave
                * otherwise (i.e., f == nfactors), look further */
               desiredcurv = SCIP_EXPRCURV_UNKNOWN;
               if( f == nfactors )
               {
                  /* if all factors are linear, then a product f_j^exp_j with f_j >= 0 is convex if
                   * - all exponents are negative, or
                   * - all except one exponent j* are negative and exp_j* >= 1 - sum_{j!=j*}exp_j, but the latter is equivalent to sum_j exp_j >= 1
                   * further, the product is concave if
                   * - all exponents are positive and the sum of exponents is <= 1.0
                   *
                   * if factors are nonlinear, then we require additionally, that for convexity
                   * - each factor is convex if exp_j >= 0, or concave if exp_j <= 0, i.e., exp_j*f_j'' >= 0
                   * and for concavity, we require that
                   * - all factors are concave, i.e., exp_j*f_j'' <= 0
                   */

                  if( nnegative == nfactors || (nnegative == nfactors-1 && SCIPisGE(scip, sum, 1.0)) )
                  {
                     /* if exponents are such that we can be convex, but children curvature does not fit, make some children linear */
                     SCIPdebugMessage("%d-variate monomial is convex (modulo sign), child curv fits = %u\n", nfactors, expcurvpos);
                     /* since current node curvature is set to unknown, there must be such a child, since otherwise the node curvature had to be convex */
                     assert(!expcurvpos);
                     desiredcurv = SCIP_EXPRCURV_CONVEX;
                  }
                  else if( npositive == nfactors && SCIPisLE(scip, sum, 1.0) )
                  {
                     /* if exponents are such that we can be concave, but children curvature does not fit, make some children linear */
                     SCIPdebugMessage("%d-variate monomial is concave (modulo sign), child curv fits = %u\n", nfactors, expcurvneg);
                     /* since current node curvature is set to unknown, there must be such a child, since otherwise the node curvature had to be concave */
                     assert(!expcurvneg);
                     desiredcurv = SCIP_EXPRCURV_CONCAVE;
                  }
                  else
                  {
                     /* exponents are such that monomial is neither convex nor concave even if children were linear
                      * thus, reformulate monomial below
                      */
                  }
               }

               if( desiredcurv != SCIP_EXPRCURV_UNKNOWN )
               {
                  for( f = 0; f < nfactors; ++f )
                  {
                     childcurv = SCIPexprgraphGetNodeCurvature(children[childidxs[f]]);  /*lint !e613*/
                     assert(childcurv != SCIP_EXPRCURV_UNKNOWN);
                     childbounds = SCIPexprgraphGetNodeBounds(children[childidxs[f]]);  /*lint !e613*/
                     assert(childbounds.inf >= 0.0 || childbounds.sup <= 0.0);

                     /* negate curvature if factor is negative */
                     if( childbounds.inf < 0.0 )
                        childcurv = SCIPexprcurvNegate(childcurv);

                     /* check if exp_j * checkcurv is convex (>= 0) and/or concave */
                     childcurv = SCIPexprcurvMultiply(SCIPexprGetMonomialExponents(monomial)[f], childcurv);
                     if( (desiredcurv == SCIP_EXPRCURV_CONVEX  && !(childcurv & SCIP_EXPRCURV_CONVEX )) ||
                        (desiredcurv == SCIP_EXPRCURV_CONCAVE && !(childcurv & SCIP_EXPRCURV_CONCAVE)) )
                     {
                        SCIPdebugMessage("reform child %d (factor %d) with curv %s into var\n",
                           childidxs[f], f, SCIPexprcurvGetName(SCIPexprgraphGetNodeCurvature(children[childidxs[f]])));  /*lint !e613*/
                        SCIP_CALL( reformNode2Var(scip, exprgraph, children[childidxs[f]], conss, nconss, naddcons, FALSE) );  /*lint !e613*/
                        modified = TRUE;
                     }
                  }
               }
            }

            if( modified )
            {
               /* refresh curvature information in node, since we changed children, it should be convex or concave now */
               SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(node, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
               assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(node)));
               assert(SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN);

               /* we are done and can proceed with the next node */
               ++i;
               break;
            }

            /* monomial can only have unknown curvature here, if it has several factors
             * or is of form x^a with x both negative and positive and a an odd or negative integer (-> INTPOWER expression)
             */
            assert(nfactors > 1 ||
               (SCIPexprgraphGetNodeBounds(children[childidxs[0]]).inf < 0.0 && SCIPexprgraphGetNodeBounds(children[childidxs[0]]).sup > 0.0 &&
                  SCIPisIntegral(scip, exponents[0]) && (exponents[0] < 0.0 || ((int)SCIPround(scip, exponents[0]) % 2 != 0)))
               );  /*lint !e613*/

            /* bilinear monomials should not come up here, since simplifier should have turned them into quadratic expression nodes */
            assert(!(nfactors == 2 && exponents[0] == 1.0 && exponents[1] == 1.0));

            /* reform monomial if it is a product, or we need it to be on the top of the graph, or if it of the form x^a with a < 0.0 (and thus x having mixed sign, see assert above)
             * thus, in the case x^a with a an odd positive integer we assume that cons_signpower will do something */
            if( nfactors > 1 || havenonlinparent || exponents[0] < 0.0 )
            {
               SCIP_EXPRGRAPHNODE* auxnode;
               SCIP_EXPRGRAPHNODE** factors;

               if( nfactors > 1 )
               {
                  SCIP_CALL( SCIPallocBufferArray(scip, &factors, nfactors) );
                  for( f = 0; f < nfactors; ++f )
                     factors[f] = children[childidxs[f]];  /*lint !e613*/
               }
               else
                  factors = &children[childidxs[0]];  /*lint !e613*/

               SCIPdebugMessage("reform monomial node, create auxvar = %u\n", havenonlinparent);
               /* get new auxnode for monomial
                * if node has parents and monomial is of indefinite form x^a, then also create auxvar for it, since otherwise we create a auxnode with unknown curvature
                * note, that the case x^a with positive and odd a will still give an indefinite node (without parents), where we assume that signpower will pick it up at some point
                */
               SCIP_CALL( reformMonomial(scip, exprgraph, nfactors, factors, exponents, &auxnode, havenonlinparent, naddcons) );

               if( nfactors > 1 )
               {
                  SCIPfreeBufferArray(scip, &factors);
               }

               /* create node for monomialcoef * auxnode + monomialconstant, if not identical to auxnode */
               if( SCIPexprgraphGetNodePolynomialConstant(node) != 0.0 || coef != 1.0 )
               {
                  SCIP_EXPRGRAPHNODE* replnode;

                  SCIP_CALL( SCIPexprgraphCreateNodeLinear(SCIPblkmem(scip), &replnode, 1, &coef, SCIPexprgraphGetNodePolynomialConstant(node)) );
                  SCIP_CALL( SCIPexprgraphAddNode(exprgraph, replnode, -1, 1, &auxnode) );
                  auxnode = replnode;
               }

               /* replace node by auxnode and refresh its curvature */
               SCIP_CALL( reformReplaceNode(exprgraph, &node, auxnode, conss, nconss) );
               SCIP_CALL( SCIPexprgraphUpdateNodeBoundsCurvature(auxnode, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, TRUE) );
               assert(!SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(auxnode)));

               break;
            }
            else
            {
               SCIPdebugMessage("no reformulation of monomial node, assume signpower will take care of it\n");
            }

            ++i;
            break;
         }

         case SCIP_EXPR_LAST:
         default:
            SCIPerrorMessage("got expression with invalid operand\n");
         }
      }
   }

   /* for constraints with concave f(g(x)) with linear g:R^n -> R, n>1, reformulate to get a univariate concave function, since this is easier to underestimate
    * @todo this does not work yet for sums of functions, e.g., polynomials with more than one monomial
    */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_EXPRGRAPHNODE* multivarnode;
      SCIP_EXPRCURV curv;

      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( consdata->exprgraphnode == NULL )
         continue;

      curv = SCIPexprgraphGetNodeCurvature(consdata->exprgraphnode);

      /* if nothing concave, then continue */
      if( (SCIPisInfinity(scip,  consdata->rhs) || curv != SCIP_EXPRCURV_CONCAVE) &&
         ( SCIPisInfinity(scip, -consdata->lhs) || curv != SCIP_EXPRCURV_CONVEX) )
         continue;

      /* search for a descendant of node that has > 1 children
       * after simiplifier run, there should be no constant expressions left */
      multivarnode = consdata->exprgraphnode;
      while( SCIPexprgraphGetNodeNChildren(multivarnode) == 1 )
         multivarnode = SCIPexprgraphGetNodeChildren(multivarnode)[0];

      /* if node expression is obviously univariate, then continue */
      if( SCIPexprgraphGetNodeNChildren(multivarnode) == 0 )
      {
         assert(SCIPexprgraphGetNodeOperator(multivarnode) == SCIP_EXPR_CONST || SCIPexprgraphGetNodeOperator(multivarnode) == SCIP_EXPR_VARIDX);
         continue;
      }

      /* if node itself is multivariate, then continue */
      if( multivarnode == consdata->exprgraphnode )
         continue;

      /* if multivarnode is a linear expression, then replace this by an auxiliary variable/node
       * mark auxiliary variable as not to multiaggregate, so SCIP cannot undo what we just did */
      if( SCIPexprgraphGetNodeCurvature(multivarnode) == SCIP_EXPRCURV_LINEAR )
      {
         SCIPdebugMessage("replace linear multivariate node %p(%d,%d) in expression of cons <%s> by auxvar\n",
            (void*)multivarnode, SCIPexprgraphGetNodeDepth(multivarnode), SCIPexprgraphGetNodePosition(multivarnode), SCIPconsGetName(conss[c]));  /*lint !e613*/
         SCIPdebug( SCIPprintCons(scip, conss[c], NULL) );  /*lint !e613*/
         SCIP_CALL( reformNode2Var(scip, exprgraph, multivarnode, conss, nconss, naddcons, TRUE) );
      }
   }

   conshdlrdata->isreformulated = TRUE;

   return SCIP_OKAY;
}

/** gets maximal absolute element of gradient of nonlinear function */
static
SCIP_RETCODE getGradientMaxElement(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expressions interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Bool             newsol,             /**< have the expression trees been evaluated at sol before? */
   SCIP_Real*            maxelem             /**< buffer to store maximal element */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(maxelem != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprint != NULL);
   assert(consdata->nexprtrees != 0 || consdata->exprgraphnode == NULL);

   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
   {
      *maxelem = 0.0;
      for( i = 0; i < consdata->nlinvars; ++i )
         if( REALABS(consdata->lincoefs[i]) > *maxelem )
            *maxelem = REALABS(consdata->lincoefs[i]);
   }
   else
   {
      *maxelem = consdata->lincoefsmax;
   }

   for( j = 0; j < consdata->nexprtrees; ++j )
   {
      int nvars;
      SCIP_Real val;

      assert(consdata->exprtrees[j] != NULL);

      nvars = SCIPexprtreeGetNVars(consdata->exprtrees[j]);

      if( newsol )
      {
         /* compile expression tree, if not done before (can happen, if called from proposeFeasibleSolution) */
         if( SCIPexprtreeGetInterpreterData(consdata->exprtrees[j]) == NULL )
         {
            SCIP_CALL( SCIPexprintCompile(exprint, consdata->exprtrees[j]) );
         }

         if( nvars == 1 )
         {
            /* in the not so unusual case that an expression has only one variable, we do not extra alloc memory */
            double varval;
            SCIP_Real grad;

            varval = SCIPgetSolVal(scip, sol, SCIPexprtreeGetVars(consdata->exprtrees[j])[0]);
            SCIP_CALL( SCIPexprintGrad(exprint, consdata->exprtrees[j], &varval, TRUE, &val, &grad) );
            if( REALABS(grad) > *maxelem )
               *maxelem = REALABS(grad);
         }
         else
         {
            SCIP_Real* x;
            SCIP_Real* grad;

            SCIP_CALL( SCIPallocBufferArray(scip, &x, nvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &grad, nvars) );

            SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, SCIPexprtreeGetVars(consdata->exprtrees[j]), x) );
            SCIP_CALL( SCIPexprintGrad(exprint, consdata->exprtrees[j], x, TRUE, &val, grad) );

            for( i = 0; i < nvars; ++i )
            {
               grad[i] *= consdata->nonlincoefs[j];
               if( REALABS(grad[i]) > *maxelem )
                  *maxelem = REALABS(grad[i]);
            }

            SCIPfreeBufferArray(scip, &x);
            SCIPfreeBufferArray(scip, &grad);
         }
      }
      else
      {
         assert(SCIPexprtreeGetInterpreterData(consdata->exprtrees[j]) != NULL);

         if( nvars == 1 )
         {
            SCIP_Real grad;

            SCIP_CALL( SCIPexprintGrad(exprint, consdata->exprtrees[j], NULL, FALSE, &val, &grad) );
            if( REALABS(grad) > *maxelem )
               *maxelem = REALABS(grad);
         }
         else
         {
            SCIP_Real* grad;

            SCIP_CALL( SCIPallocBufferArray(scip, &grad, nvars) );

            SCIP_CALL( SCIPexprintGrad(exprint, consdata->exprtrees[j], NULL, FALSE, &val, grad) );

            for( i = 0; i < nvars; ++i )
            {
               grad[i] *= consdata->nonlincoefs[j];
               if( REALABS(grad[i]) > *maxelem )
                  *maxelem = REALABS(grad[i]);
            }

            SCIPfreeBufferArray(scip, &grad);
         }
      }
   }

   return SCIP_OKAY;
}

/** computes activity and violation of a constraint
 * during presolving and if the constraint is active, it is assumes that SCIPexprgraphEval has been called for sol before
 */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_SOL*             sol                 /**< solution or NULL if LP solution should be used */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real varval;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   consdata->activity = 0.0;
   varval = 0.0;

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      varval = SCIPgetSolVal(scip, sol, consdata->linvars[i]);
      if( SCIPisInfinity(scip, REALABS(varval)) )
      {
         consdata->activity = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      consdata->activity += consdata->lincoefs[i] * SCIPgetSolVal(scip, sol, consdata->linvars[i]);
   }

   for( i = 0; i < consdata->nexprtrees; ++i )
   {
      SCIP_Real val;
      int nvars;

      /* compile expression tree, if not done before */
      if( SCIPexprtreeGetInterpreterData(consdata->exprtrees[i]) == NULL )
      {
         SCIP_CALL( SCIPexprintCompile(exprint, consdata->exprtrees[i]) );
      }

      nvars = SCIPexprtreeGetNVars(consdata->exprtrees[i]);

      if( nvars == 1 )
      {
         /* in the not so unusual case that an expression has only one variable, we do not extra alloc memory */
         varval = SCIPgetSolVal(scip, sol, SCIPexprtreeGetVars(consdata->exprtrees[i])[0]);
         SCIP_CALL( SCIPexprintEval(exprint, consdata->exprtrees[i], &varval, &val) );
      }
      else
      {
         SCIP_Real* x;

         SCIP_CALL( SCIPallocBufferArray(scip, &x, nvars) );
         SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, SCIPexprtreeGetVars(consdata->exprtrees[i]), x) );
         SCIP_CALL( SCIPexprintEval(exprint, consdata->exprtrees[i], x, &val) );
         SCIPfreeBufferArray(scip, &x);
      }

      if( SCIPisInfinity(scip, REALABS(val)) || !finite(val) )
      {
         consdata->activity = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      consdata->activity += consdata->nonlincoefs[i] * val;
   }

   if( consdata->nexprtrees == 0 && consdata->exprgraphnode != NULL )
   {
      SCIP_Real val;

      assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);

      val = SCIPexprgraphGetNodeVal(consdata->exprgraphnode);
      assert(val != SCIP_INVALID);  /*lint !e777*/

      if( !finite(val) || SCIPisInfinity(scip, REALABS(val)) )
      {
         consdata->activity = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhsviol = SCIPinfinity(scip);
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhsviol = SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      consdata->activity += val;
   }

   if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisGT(scip, consdata->lhs - consdata->activity, SCIPfeastol(scip)) )
      consdata->lhsviol = consdata->lhs - consdata->activity;
   else
      consdata->lhsviol = 0.0;

   if( !SCIPisInfinity(scip,  consdata->rhs) && SCIPisGT(scip, consdata->activity - consdata->rhs, SCIPfeastol(scip)) )
      consdata->rhsviol = consdata->activity - consdata->rhs;
   else
      consdata->rhsviol = 0.0;

   /* scale violation by infinity norm of gradient, if violated
    * do only if we are linear or have expression trees, thus, not during presolve
    */
   if( (consdata->lhsviol > 0.0 || consdata->rhsviol > 0.0) && (consdata->exprgraphnode == NULL || consdata->nexprtrees > 0) )
   {
      SCIP_Real norm;

      SCIP_CALL( getGradientMaxElement(scip, exprint, cons, sol, FALSE, &norm) );
      if( norm > 1.0 && !SCIPisInfinity(scip, norm) )
      {
         consdata->lhsviol /= norm;
         consdata->rhsviol /= norm;
      }
   }

   return SCIP_OKAY;
}

/** computes violation of a set of constraints */
static
SCIP_RETCODE computeViolations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRGRAPH*       exprgraph,          /**< expression graph of constraint handler */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_CONS**           maxviolcon          /**< buffer to store constraint with largest violation, or NULL if solution is feasible */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      viol;
   SCIP_Real      maxviol;
   int            c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(maxviolcon != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_Real* varvals;

      assert(exprgraph != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(exprgraph)) );
      SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPexprgraphGetNVars(exprgraph), (SCIP_VAR**)SCIPexprgraphGetVars(exprgraph), varvals) );

      SCIP_CALL( SCIPexprgraphEval(exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }

   *maxviolcon = NULL;

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      SCIP_CALL( computeViolation(scip, exprint, conss[c], sol) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if( viol > maxviol && SCIPisGT(scip, viol, SCIPfeastol(scip)) )
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }

      /* SCIPdebugMessage("constraint <%s> violated by (%g, %g), activity = %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->activity); */
   }

   return SCIP_OKAY;
}

/** adds linearization of a constraints expression tree in reference point to a row */
static
SCIP_RETCODE addLinearization(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a linearization should be added */
   SCIP_Real*            x,                  /**< value of expression tree variables where to generate cut */
   SCIP_Bool             newx,               /**< whether the last evaluation of the expression with the expression interpreter was not at x */
   SCIP_ROW*             row,                /**< row where to add linearization */
   SCIP_Bool*            success             /**< buffer to store whether a linearization was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real treecoef;
   SCIP_Real val;
   SCIP_Real* grad;
   SCIP_Real constant;
   SCIP_Bool perturbedx;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x    != NULL);
   assert(row  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(newx || SCIPexprtreeGetInterpreterData(exprtree) != NULL);

   treecoef = consdata->nonlincoefs[exprtreeidx];

   *success = FALSE;

   /* compile expression if evaluated the first time; can only happen if newx is FALSE */
   if( newx && SCIPexprtreeGetInterpreterData(exprtree) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(exprint, exprtree) );
   }

   nvars = SCIPexprtreeGetNVars(exprtree);
   SCIP_CALL( SCIPallocBufferArray(scip, &grad, nvars) );

   perturbedx = FALSE;
   do
   {
      /* get value and gradient */
      SCIP_CALL( SCIPexprintGrad(exprint, exprtree, x, newx, &val, grad) );
      if( finite(val) && !SCIPisInfinity(scip, REALABS(val)) )
      {
         val *= treecoef;
         /* check gradient entries and compute constant f(refx) - grad * refx */
         constant = val;
         for( i = 0; i < nvars; ++i )
         {
            if( !finite(grad[i]) || SCIPisInfinity(scip, grad[i]) || SCIPisInfinity(scip, -grad[i]) )
               break;

            grad[i] *= treecoef;
            constant -= grad[i] * x[i];

            /* coefficients smaller than epsilon are rounded to 0.0 when added to row, this can be wrong if variable value is very large (bad numerics)
             * in this case, set gradient to 0.0 here, but modify constant so that cut is still valid (if possible)
             * i.e., estimate grad[i]*x >= grad[i] * bound(x) or grad[i]*x <= grad[i] * bound(x), depending on whether we compute an underestimator (convex) or an overestimator (concave)
             * if required bound of x is not finite, then give up
             */
            if( grad[i] != 0.0 && SCIPisZero(scip, grad[i]) )
            {
               SCIP_VAR* var;
               SCIP_Real xbnd;

               var = SCIPexprtreeGetVars(exprtree)[i];
               if( consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONVEX )
               {
                  xbnd = grad[i] > 0.0 ? SCIPvarGetLbGlobal(var) : SCIPvarGetUbGlobal(var);
               }
               else
               {
                  assert(consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONCAVE);
                  xbnd = grad[i] > 0.0 ? SCIPvarGetUbGlobal(var) : SCIPvarGetLbGlobal(var);
               }
               if( !SCIPisInfinity(scip, REALABS(xbnd)) )
               {
                  SCIPdebugMessage("var <%s> [%g,%g] has tiny gradient %g, replace coefficient by constant %g\n",
                     SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), grad[i], grad[i] * xbnd);
                  constant += grad[i] * xbnd;
                  grad[i] = 0.0;
               }
               else
               {
                  *success = FALSE;
                  SCIPdebugMessage("skipping linearization, var <%s> [%g,%g] has tiny gradient %g but no finite bound in this direction\n",
                     SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var), grad[i]);
                  SCIPfreeBufferArray(scip, &grad);
                  return SCIP_OKAY;
               }
            }
         }

         if( i == nvars )
            break;
      }

      SCIPdebugMessage("got nonfinite value in evaluation or gradient of <%s>: ", SCIPconsGetName(cons));
      if( !perturbedx )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         SCIPdebugPrintf("perturbing reference point and trying again\n");
         for( i = 0; i < nvars; ++i )
         {
            lb = SCIPvarGetLbGlobal(SCIPexprtreeGetVars(exprtree)[i]);
            ub = SCIPvarGetUbGlobal(SCIPexprtreeGetVars(exprtree)[i]);
            if( SCIPisEQ(scip, x[i], lb) )
               x[i] += MIN(0.9*(ub-lb), i*SCIPfeastol(scip));  /*lint !e666*/
            else if( SCIPisEQ(scip, x[i], ub) )
               x[i] -= MIN(0.9*(ub-lb), i*SCIPfeastol(scip));  /*lint !e666*/
            else
               x[i] += MIN3(0.9*(ub-x[i]), 0.9*(x[i]-lb), i*SCIPfeastol(scip)) * (i%2 != 0 ? -1.0 : 1.0);  /*lint !e666*/
         }
         newx = TRUE;
         perturbedx = TRUE;
      }
      else
      {
         SCIPdebugPrintf("skipping linearization\n");
         SCIPfreeBufferArray(scip, &grad);
         return SCIP_OKAY;
      }
   } while( TRUE );  /*lint !e506*/

   /* add linearization to SCIP row */
   if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, row, SCIProwGetLhs(row) - constant) );  /*lint !e644*/
   }
   if( !SCIPisInfinity(scip,  SCIProwGetRhs(row)) )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, row, SCIProwGetRhs(row) - constant) );
   }
   SCIP_CALL( SCIPaddVarsToRow(scip, row, nvars, SCIPexprtreeGetVars(exprtree), grad) );

   *success = TRUE;

   SCIPfreeBufferArray(scip, &grad);

   SCIPdebugMessage("added linearization for tree %d of constraint <%s>\n", exprtreeidx, SCIPconsGetName(cons));
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

   return SCIP_OKAY;
}

/** adds secant of a constraints univariate expression tree in reference point to a row */
static
SCIP_RETCODE addConcaveEstimatorUnivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_ROW*             row,                /**< row where to add secant */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real      treecoef;
   SCIP_VAR*      var;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      vallb;
   SCIP_Real      valub;
   SCIP_Real      slope;
   SCIP_Real      constant;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(SCIPexprtreeGetNVars(exprtree) == 1);

   treecoef = consdata->nonlincoefs[exprtreeidx];

   *success = FALSE;

   var = SCIPexprtreeGetVars(exprtree)[0];
   xlb = SCIPvarGetLbLocal(var);
   xub = SCIPvarGetUbLocal(var);

   /* if variable is fixed, then cannot really compute secant */
   if( SCIPisEQ(scip, xlb, xub) )
      return SCIP_OKAY;
   /* if variable is unbounded, then cannot really compute secant */
   if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
      return SCIP_OKAY;
   assert(SCIPisLT(scip, xlb, xub));

   SCIP_CALL( SCIPexprtreeEval(exprtree, &xlb, &vallb) );
   if( !finite(vallb) || SCIPisInfinity(scip, REALABS(vallb)) )
   {
      SCIPdebugMessage("skip secant for tree %d of constraint <%s> since function cannot be evaluated in lower bound\n", exprtreeidx, SCIPconsGetName(cons));
      return SCIP_OKAY;
   }
   vallb *= treecoef;

   SCIP_CALL( SCIPexprtreeEval(exprtree, &xub, &valub) );
   if( !finite(valub) || SCIPisInfinity(scip, REALABS(valub)) )
   {
      SCIPdebugMessage("skip secant for tree %d of constraint <%s> since function cannot be evaluated in lower bound\n", exprtreeidx, SCIPconsGetName(cons));
      return SCIP_OKAY;
   }
   valub *= treecoef;

   slope = (valub - vallb) / (xub - xlb);
   constant = vallb - slope * xlb;

   /* add secant to SCIP row */
   if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, row, SCIProwGetLhs(row) - constant) );
   }
   if( !SCIPisInfinity(scip,  SCIProwGetRhs(row)) )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, row, SCIProwGetRhs(row) - constant) );
   }
   SCIP_CALL( SCIPaddVarsToRow(scip, row, 1, &var, &slope) );

   *success = TRUE;

   SCIPdebugMessage("added secant for tree %d of constraint <%s>, slope = %g\n", exprtreeidx, SCIPconsGetName(cons), slope);
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

   return SCIP_OKAY;
}

/** given three points, constructs coefficient of equation for hyperplane generated by these three points
 * Three points a, b, and c are given.
 * Computes coefficients alpha, beta, gamma, and delta, such that a, b, and c, satisfy
 * alpha * x1 + beta * x2 + gamma * x3 = delta and gamma >= 0.0.
 */
static
void getAlphaBetaGammaDelta(
   SCIP_Real             a1,                 /* first coordinate of a */
   SCIP_Real             a2,                 /* second coordinate of a */
   SCIP_Real             a3,                 /* third coordinate of a */
   SCIP_Real             b1,                 /* first coordinate of b */
   SCIP_Real             b2,                 /* second coordinate of b */
   SCIP_Real             b3,                 /* third coordinate of b */
   SCIP_Real             c1,                 /* first coordinate of c */
   SCIP_Real             c2,                 /* second coordinate of c */
   SCIP_Real             c3,                 /* third coordinate of c */
   SCIP_Real*            alpha,              /* coefficient of first coordinate */
   SCIP_Real*            beta,               /* coefficient of second coordinate */
   SCIP_Real*            gamma_,             /* coefficient of third coordinate */
   SCIP_Real*            delta               /* constant right-hand side */
   )
{
   assert(alpha != NULL);
   assert(beta  != NULL);
   assert(gamma_ != NULL);
   assert(delta != NULL);

   *alpha  = -b3*c2 + a3*(-b2+c2) + a2*(b3-c3) + b2*c3;
   *beta   = -(-b3*c1 + a3*(-b1+c1) + a1*(b3-c3) + b1*c3);
   *gamma_ = -a2*b1 + a1*b2 + a2*c1 - b2*c1 - a1*c2 + b1*c2;
   *delta  = -a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3;

   if( *gamma_ < 0.0 )
   {
      *alpha  = -*alpha;
      *beta   = -*beta;
      *gamma_ = -*gamma_;
      *delta  = -*delta;
   }
}

/** adds estimator of a constraints bivariate expression tree to a row
 * a reference point is given to decide which hyperplane to choose
 */
static
SCIP_RETCODE addConcaveEstimatorBivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_Real*            ref,                /**< reference values of expression tree variables where to generate cut */
   SCIP_ROW*             row,                /**< row where to add secant */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real      treecoef;
   SCIP_VAR*      x;
   SCIP_VAR*      y;
   SCIP_Real      xlb;
   SCIP_Real      xub;
   SCIP_Real      ylb;
   SCIP_Real      yub;

   SCIP_Real      coefx;
   SCIP_Real      coefy;
   SCIP_Real      constant;

   SCIP_Real      p1[2];
   SCIP_Real      p2[2];
   SCIP_Real      p3[2];
   SCIP_Real      p4[2];
   SCIP_Real      p1val, p2val, p3val, p4val;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref  != NULL);
   assert(row  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(SCIPexprtreeGetNVars(exprtree) == 2);

   treecoef = consdata->nonlincoefs[exprtreeidx];

   *success = FALSE;

   x = SCIPexprtreeGetVars(exprtree)[0];
   y = SCIPexprtreeGetVars(exprtree)[1];
   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);
   ylb = SCIPvarGetLbLocal(y);
   yub = SCIPvarGetUbLocal(y);

   if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) || SCIPisInfinity(scip, -ylb) || SCIPisInfinity(scip, yub) )
   {
      SCIPdebugMessage("skip bivariate secant since <%s> or <%s> is unbounded\n", SCIPvarGetName(x), SCIPvarGetName(y));
      return SCIP_OKAY;
   }

   if( SCIPisEQ(scip, xlb, xub) && SCIPisEQ(scip, ylb, yub) )
   {
      SCIPdebugMessage("skip bivariate secant since both <%s> and <%s> are fixed\n", SCIPvarGetName(x), SCIPvarGetName(y));
      return SCIP_OKAY;
   }

   /* reference point should not be outside of bounds */
   assert(SCIPisFeasLE(scip, xlb, ref[0]));
   assert(SCIPisFeasGE(scip, xub, ref[0]));
   ref[0] = MIN(xub, MAX(xlb, ref[0]));
   assert(SCIPisFeasLE(scip, ylb, ref[1]));
   assert(SCIPisFeasGE(scip, yub, ref[1]));
   ref[1] = MIN(yub, MAX(ylb, ref[1]));

   /* lower left */
   p1[0] = xlb;
   p1[1] = ylb;

   /* lower right */
   p2[0] = xub;
   p2[1] = ylb;

   /* upper right */
   p3[0] = xub;
   p3[1] = yub;

   /* upper left */
   p4[0] = xlb;
   p4[1] = yub;

   if( SCIPisEQ(scip, xlb, xub) )
   {
      /* secant between p1 and p4: p1val + [(p4val - p1val) / (yub - ylb)] * (y - ylb) */
      assert(!SCIPisEQ(scip, ylb, yub));

      SCIP_CALL( SCIPexprtreeEval(exprtree, p1, &p1val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p4, &p4val) );
      if( !finite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !finite(p4val) || SCIPisInfinity(scip, REALABS(p4val)) )
      {
         SCIPdebugMessage("skip secant for tree %d of constraint <%s> since function cannot be evaluated\n", exprtreeidx, SCIPconsGetName(cons));
         return SCIP_OKAY;
      }
      p1val *= treecoef;
      p4val *= treecoef;

      coefx = 0.0;
      coefy = (p4val - p1val) / (yub - ylb);
      constant = p1val - coefy * ylb;
   }
   else if( SCIPisEQ(scip, ylb, yub) )
   {
      /* secant between p1 and p2: p1val + [(p2val - p1val) / (xub - xlb)] * (x - xlb) */
      assert(!SCIPisEQ(scip, xlb, xub));

      SCIP_CALL( SCIPexprtreeEval(exprtree, p1, &p1val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p2, &p2val) );
      if( !finite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !finite(p2val) || SCIPisInfinity(scip, REALABS(p2val)) )
      {
         SCIPdebugMessage("skip secant for tree %d of constraint <%s> since function cannot be evaluated\n", exprtreeidx, SCIPconsGetName(cons));
         return SCIP_OKAY;
      }

      p1val *= treecoef;
      p2val *= treecoef;

      coefx = (p2val - p1val) / (xub - xlb);
      coefy = 0.0;
      constant = p1val - coefx * xlb;
   }
   else
   {
      SCIP_Real alpha, beta, gamma_, delta;
      SCIP_Bool tryother;
      SCIP_Bool doover;

      /* if function is convex, then we want an overestimator, otherwise we want an underestimator */
      assert(consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONVEX || consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONCAVE);
      doover = (consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONVEX);  /*lint !e641*/

      SCIP_CALL( SCIPexprtreeEval(exprtree, p1, &p1val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p2, &p2val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p3, &p3val) );
      SCIP_CALL( SCIPexprtreeEval(exprtree, p4, &p4val) );
      if( !finite(p1val) || SCIPisInfinity(scip, REALABS(p1val)) || !finite(p2val) || SCIPisInfinity(scip, REALABS(p2val)) ||
         ! finite(p3val) || SCIPisInfinity(scip, REALABS(p3val)) || !finite(p4val) || SCIPisInfinity(scip, REALABS(p4val)) )
      {
         SCIPdebugMessage("skip secant for tree %d of constraint <%s> since function cannot be evaluated\n", exprtreeidx, SCIPconsGetName(cons));
         return SCIP_OKAY;
      }
      p1val *= treecoef;
      p2val *= treecoef;
      p3val *= treecoef;
      p4val *= treecoef;

      /* if we want an underestimator, flip f(x,y), i.e., do as if we compute an overestimator for -f(x,y) */
      if( !doover )
      {
         p1val = -p1val;
         p2val = -p2val;
         p3val = -p3val;
         p4val = -p4val;
      }

      SCIPdebugMessage("p1 = (%g, %g), f(p1) = %g\n", p1[0], p1[1], p1val);
      SCIPdebugMessage("p2 = (%g, %g), f(p2) = %g\n", p2[0], p2[1], p2val);
      SCIPdebugMessage("p3 = (%g, %g), f(p3) = %g\n", p3[0], p3[1], p3val);
      SCIPdebugMessage("p4 = (%g, %g), f(p4) = %g\n", p4[0], p4[1], p4val);

      /* Compute coefficients alpha, beta, gamma (>0), delta such that
       *   alpha*x + beta*y + gamma*z = delta
       * is satisfied by at least three of the corner points (p1,f(p1)), ..., (p4,f(p4)) and
       * the fourth corner point lies below this hyperplane.
       * Since we assume that f is convex, we then know that all points (x,y,f(x,y)) are below this hyperplane, i.e.,
       *    alpha*x + beta*y - delta <= -gamma * f(x,y),
       * or, equivalently,
       *   -alpha/gamma*x - beta/gamma*y + delta/gamma >= f(x,y).
       */

      tryother = FALSE;
      if( ref[1] <= ylb + (yub - ylb)/(xub - xlb) * (ref[0] - xlb) )
      {
         getAlphaBetaGammaDelta(p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val, &alpha, &beta, &gamma_, &delta);
         assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
         assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
         assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));

         /* if hyperplane through p1,p2,p3 does not overestimate f(p4), then it must be the other variant */
         if( alpha * p4[0] + beta * p4[1] + gamma_ * p4val > delta )
            tryother = TRUE;
         else if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
            (      !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
         {
            /* if numerically bad, take alternative hyperplane */
            getAlphaBetaGammaDelta(p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val, &alpha, &beta, &gamma_, &delta);
            assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
            assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

            /* if hyperplane through p1,p3,p4 does not overestimate f(p2), then it must be the other variant */
            if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val > delta )
               tryother = TRUE;
         }
      }
      else
      {
         getAlphaBetaGammaDelta(p1[0], p1[1], p1val, p3[0], p3[1], p3val, p4[0], p4[1], p4val, &alpha, &beta, &gamma_, &delta);
         assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
         assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
         assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

         /* if hyperplane through p1,p3,p4 does not overestimate f(p2), then it must be the other variant */
         if( alpha * p2[0] + beta * p2[1] + gamma_ * p2val > delta )
            tryother = TRUE;
         else if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
            (      !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
         {
            /* if numerically bad, take alternative */
            getAlphaBetaGammaDelta(p1[0], p1[1], p1val, p2[0], p2[1], p2val, p3[0], p3[1], p3val, &alpha, &beta, &gamma_, &delta);
            assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
            assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));

            /* if hyperplane through p1,p2,p3 does not overestimate f(p4), then it must be the other variant */
            if( alpha * p4[0] + beta * p4[1] + gamma_ * p4val > delta )
               tryother = TRUE;
         }
      }

      if( tryother )
      {
         if( ref[1] <= yub + (ylb - yub)/(xub - xlb) * (ref[0] - xlb) )
         {
            getAlphaBetaGammaDelta(p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1], p4val, &alpha, &beta, &gamma_, &delta);

            /* hyperplane should be above (p3,f(p3)) and other points should lie on hyperplane */
            assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
            assert(SCIPisRelLE(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
            assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

            if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
               ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
            {
               /* if numerically bad, take alternative */
               getAlphaBetaGammaDelta(p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1], p4val, &alpha, &beta, &gamma_, &delta);

               /* hyperplane should be above (p1,f(p1)) and other points should lie on hyperplane */
               assert(SCIPisRelLE(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
               assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
               assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
               assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));
            }
         }
         else
         {
            getAlphaBetaGammaDelta(p2[0], p2[1], p2val, p3[0], p3[1], p3val, p4[0], p4[1], p4val, &alpha, &beta, &gamma_, &delta);

            /* hyperplane should be above (p1,f(p1)) and other points should lie on hyperplane */
            assert(SCIPisRelLE(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
            assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
            assert(SCIPisRelEQ(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
            assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));

            if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, alpha/gamma_)) ||
               ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, beta /gamma_)) )
            {
               /* if numerically bad, take alternative */
               getAlphaBetaGammaDelta(p1[0], p1[1], p1val, p2[0], p2[1], p2val, p4[0], p4[1], p4val, &alpha, &beta, &gamma_, &delta);

               /* hyperplane should be above (p3,f(p3)) and other points should lie on hyperplane */
               assert(SCIPisRelEQ(scip, alpha * p1[0] + beta * p1[1] + gamma_ * p1val, delta));
               assert(SCIPisRelEQ(scip, alpha * p2[0] + beta * p2[1] + gamma_ * p2val, delta));
               assert(SCIPisRelLE(scip, alpha * p3[0] + beta * p3[1] + gamma_ * p3val, delta));
               assert(SCIPisRelEQ(scip, alpha * p4[0] + beta * p4[1] + gamma_ * p4val, delta));
            }
         }
      }

      SCIPdebugMessage("alpha = %g, beta = %g, gamma = %g, delta = %g\n", alpha, beta, gamma_, delta);

      /* check if bad luck: should not happen if xlb != xub and ylb != yub and numerics are fine */
      if( SCIPisZero(scip, gamma_) )
         return SCIP_OKAY;
      assert(!SCIPisNegative(scip, gamma_));

      /* flip hyperplane */
      if( !doover )
         gamma_ = -gamma_;

      coefx    = -alpha / gamma_;
      coefy    = -beta  / gamma_;
      constant =  delta / gamma_;

      /* if we loose coefficients because division by gamma makes them < SCIPepsilon(scip), then better not generate a cut here */
      if( (!SCIPisZero(scip, alpha) && SCIPisZero(scip, coefx)) ||
         ( !SCIPisZero(scip, beta)  && SCIPisZero(scip, coefy)) )
      {
         SCIPdebugMessage("skip bivar secant for <%s> tree %d due to bad numerics\n", SCIPconsGetName(cons), exprtreeidx);
         return SCIP_OKAY;
      }
   }

   /* add hyperplane coefs to SCIP row */
   if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, row, SCIProwGetLhs(row) - constant) );
   }
   if( !SCIPisInfinity(scip,  SCIProwGetRhs(row)) )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, row, SCIProwGetRhs(row) - constant) );
   }
   SCIP_CALL( SCIPaddVarsToRow(scip, row, 1, &x, &coefx) );
   SCIP_CALL( SCIPaddVarsToRow(scip, row, 1, &y, &coefy) );

   *success = TRUE;

   SCIPdebugMessage("added bivariate secant for tree %d of constraint <%s>\n", exprtreeidx, SCIPconsGetName(cons));
   SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );

   return SCIP_OKAY;
}

/** adds estimator of a constraints multivariate expression tree to a row
 * Given concave function f(x) and reference point ref.
 * Let (v_i: i=1,...,n) be corner points of current domain of x.
 * Find (coef,constant) such that <coef,v_i> + constant <= f(v_i) (cut validity) and
 * such that <coef, ref> + constant is maximized (cut efficacy).
 * Then <coef, x> + constant <= f(x) for all x in current domain.
 *
 * Similar to compute an overestimator for a convex function f(x).
 * Find (coef,constant) such that <coef,v_i> + constant >= f(v_i) and
 * such that <coef, ref> + constant is minimized.
 * Then <coef, x> + constant >= f(x) for all x in current domain.
 */
static
SCIP_RETCODE addConcaveEstimatorMultivariate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_Real*            ref,                /**< reference values of expression tree variables where to generate cut */
   SCIP_ROW*             row,                /**< row where to add secant */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real treecoef;
   SCIP_LPI* lpi;
   SCIP_Bool doupper;
   SCIP_Real funcval;
   SCIP_Real lpobj;
   SCIP_RETCODE lpret;

   SCIP_VAR** vars;
   int nvars;

   int ncols;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int nrows;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   int nnonz;
   int* beg;
   int* ind;
   SCIP_Real* val;

   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ref != NULL);
   assert(row != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);

   nvars = SCIPexprtreeGetNVars(exprtree);
   assert(nvars >= 2);

   *success = FALSE;

   /* size of LP is exponential in number of variables of tree, so do only for small trees */
   if( nvars > 10 )
   {
      SCIPwarningMessage("concave function in constraint <%s> too high-dimensional to compute underestimator\n", SCIPconsGetName(cons));
      return SCIP_OKAY;
   }

   treecoef = consdata->nonlincoefs[exprtreeidx];
   vars = SCIPexprtreeGetVars(exprtree);

   /* check whether bounds are finite
    * make sure reference point is strictly within bounds
    * otherwise we can easily get an unbounded LP below, e.g., with instances like ex6_2_* from GlobalLib
    */
   for( j = 0; j < nvars; ++j )
   {
      if( SCIPisInfinity(scip, -SCIPvarGetLbLocal(vars[j])) || SCIPisInfinity(scip, SCIPvarGetUbLocal(vars[j])) )
      {
         SCIPdebugMessage("cannot compute underestimator for concave because variable <%s> is unbounded\n", SCIPvarGetName(vars[j]));
         return SCIP_OKAY;
      }
      assert(SCIPisFeasLE(scip, SCIPvarGetLbLocal(vars[j]), ref[j]));
      assert(SCIPisFeasGE(scip, SCIPvarGetUbLocal(vars[j]), ref[j]));
      ref[j] = MIN(SCIPvarGetUbLocal(vars[j]), MAX(SCIPvarGetLbLocal(vars[j]), ref[j]));  /*lint !e666*/
   }

   assert(consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONVEX || consdata->curvatures[exprtreeidx] == SCIP_EXPRCURV_CONCAVE);
   doupper = (consdata->curvatures[exprtreeidx] & SCIP_EXPRCURV_CONVEX);  /*lint !e641*/

   lpi = NULL;

   /* columns are cut coefficients plus constant */
   ncols = nvars + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );

   /* one row for each corner of domain, i.e., 2^nvars many */
   nrows = (int)(1u << nvars);
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nrows) );

   /* dense coefficients matrix, i.e., ncols * nrows many potential nonzeros */
   nnonz = nrows * ncols;
   SCIP_CALL( SCIPallocBufferArray(scip, &beg, nrows+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ind, nnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nnonz) );

   /* setup LP data */
   for( i = 0; i < nrows; ++i )
   {
      beg[i] = i * ncols;
      /* assemble corner point */
      SCIPdebugMessage("f(");
      for( j = 0; j < nvars; ++j )
      {
         /* if j'th bit of row index i is set, then take upper bound on var j, otherwise lower bound var j
          * we check this by shifting i for j positions to the right and checking whether the j'th bit is set */
         if( ((unsigned int)i >> j) & 0x1 )
            val[i * ncols + j] = SCIPvarGetUbLocal(vars[j]);
         else
            val[i * ncols + j] = SCIPvarGetLbLocal(vars[j]);
         SCIPdebugPrintf("%g, ", val[i*ncols+j]);
         assert(!SCIPisInfinity(scip, REALABS(val[i*ncols+j])));

         ind[i * ncols + j] = j;
      }

      /* evaluate function in current corner */
      SCIP_CALL( SCIPexprtreeEval(exprtree, &val[i*ncols], &funcval) );
      SCIPdebugPrintf(") = %g\n", funcval);

      if( !finite(funcval) || SCIPisInfinity(scip, REALABS(funcval)) )
      {
         SCIPdebugMessage("cannot compute underestimator for concave because constaint <%s> cannot be evaluated\n", SCIPconsGetName(cons));
         goto TERMINATE;
      }

      funcval *= treecoef;

      if( !doupper )
      {
         lhs[i] = -SCIPlpiInfinity(lpi);
         rhs[i] = funcval;
      }
      else
      {
         lhs[i] = funcval;
         rhs[i] = SCIPlpiInfinity(lpi);
      }

      /* coefficient for constant is 1.0 */
      val[i * ncols + nvars] = 1.0;
      ind[i * ncols + nvars] = nvars;
   }
   beg[nrows] = nnonz;

   for( j = 0; j < ncols; ++j )
   {
      lb[j] = -SCIPlpiInfinity(lpi);
      ub[j] =  SCIPlpiInfinity(lpi);
   }

   /* objective coefficients are reference points, and an additional 1.0 for the constant */
   BMScopyMemoryArray(obj, ref, nvars);
   obj[nvars] = 1.0;

   /* get function value in reference point, so we can use this as a cutoff */
   SCIP_CALL( SCIPexprtreeEval(exprtree, ref, &funcval) );
   funcval *= treecoef;

   SCIP_CALL( SCIPlpiCreate(&lpi, "concaveunderest", doupper ? SCIP_OBJSEN_MINIMIZE : SCIP_OBJSEN_MAXIMIZE) );
   SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );
   SCIP_CALL( SCIPlpiAddRows(lpi, nrows, lhs, rhs, NULL, nnonz, beg, ind, val) );

   /* make use of this convenient features, since for us nrows >> ncols */
   /*SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_ROWREPSWITCH, 5.0) ); */
   /* get accurate coefficients */
   SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_FEASTOL, SCIPfeastol(scip)/100.0) );
   SCIP_CALL( SCIPlpiSetRealpar(lpi, doupper ? SCIP_LPPAR_LOBJLIM : SCIP_LPPAR_UOBJLIM, funcval) );
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, 10 * nvars) );
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_SCALING, 1) );
   SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FROMSCRATCH, 1) );

   /* SCIPdebug( SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPINFO, 1) ) ); */

   lpret = SCIPlpiSolveDual(lpi);
   if( lpret != SCIP_OKAY )
   {
      SCIPwarningMessage("solving auxiliary LP for underestimator of concave function returned %d\n", lpret);
      goto TERMINATE;
   }

   if( !SCIPlpiIsPrimalFeasible(lpi) )
   {
      SCIPdebugMessage("failed to find feasible solution for auxiliary LP for underestimator of concave function, iterlimexc = %u, cutoff = %u, unbounded = %u\n", SCIPlpiIsIterlimExc(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsPrimalUnbounded(lpi));
      goto TERMINATE;
   }
   /* should be either solved to optimality, or the objective or iteration limit be hit */
   assert(SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsIterlimExc(lpi));

   /* setup row coefficient, reuse obj array to store LP sol values */
   SCIP_CALL( SCIPlpiGetSol(lpi, &lpobj, obj, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVarsToRow(scip, row, nvars, vars, obj) );

   /* check that computed hyperplane is on right side of function in refpoint
    * if numerics is very bad (e.g., st_e32), then even this can happen */
   if( (!doupper && SCIPisFeasGT(scip, lpobj, funcval)) || (doupper && SCIPisFeasGT(scip, funcval, lpobj)) )
   {
      SCIPwarningMessage("computed cut does not underestimate concave function in refpoint\n");
      goto TERMINATE;
   }
   assert( doupper || SCIPisFeasLE(scip, lpobj, funcval) );
   assert(!doupper || SCIPisFeasLE(scip, funcval, lpobj) );

   /* substract constant from lhs or rhs */
   if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, row, SCIProwGetLhs(row) - obj[nvars]) );
   }
   if( !SCIPisInfinity(scip,  SCIProwGetRhs(row)) )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, row, SCIProwGetRhs(row) - obj[nvars]) );
   }

   *success = TRUE;

 TERMINATE:
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lhs);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &beg);
   SCIPfreeBufferArray(scip, &ind);
   SCIPfreeBufferArray(scip, &val);

   if( lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&lpi) );
   }

   return SCIP_OKAY;
}

/** adds estimator from interval gradient of a constraints univariate expression tree to a row
 * a reference point is used to decide in which corner to generate the cut
 */
static
SCIP_RETCODE addIntervalGradientEstimator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   int                   exprtreeidx,        /**< for which tree a secant should be added */
   SCIP_Real*            x,                  /**< value of expression tree variables where to generate cut */
   SCIP_Bool             newx,               /**< whether the last evaluation of the expression with the expression interpreter was not at x */
   SCIP_Bool             overestimate,       /**< whether to compute an overestimator instead of an underestimator */
   SCIP_ROW*             row,                /**< row where to add secant */
   SCIP_Bool*            success             /**< buffer to store whether a secant was succefully added to the row */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_Real treecoef;
   SCIP_Real* coefs;
   SCIP_Real constant;
   SCIP_Real val;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_INTERVAL* box;
   SCIP_INTERVAL* intgrad;
   SCIP_INTERVAL intval;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(x    != NULL);
   assert(row  != NULL);
   assert(success != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(exprtreeidx >= 0);
   assert(exprtreeidx < consdata->nexprtrees);
   assert(consdata->exprtrees != NULL);

   exprtree = consdata->exprtrees[exprtreeidx];
   assert(exprtree != NULL);
   assert(newx || SCIPexprtreeGetInterpreterData(exprtree) != NULL);

   *success = FALSE;

   /* skip interval gradient if expression interpreter cannot compute interval gradients */
   if( !(SCIPexprintGetCapability() & SCIP_EXPRINTCAPABILITY_INTGRADIENT) )
      return SCIP_OKAY;

   nvars = SCIPexprtreeGetNVars(exprtree);
   vars = SCIPexprtreeGetVars(exprtree);

   box = NULL;
   intgrad = NULL;
   coefs = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &box, nvars) );

   /* move reference point to bounds, setup box */
   for( i = 0; i < nvars; ++i )
   {
      lb  = SCIPvarGetLbLocal(vars[i]);
      ub  = SCIPvarGetUbLocal(vars[i]);
      if( SCIPisInfinity(scip, -lb) )
      {
         if( SCIPisInfinity(scip, ub) )
         {
            SCIPdebugMessage("skip interval gradient estimator for constraint <%s> because variable <%s> is still unbounded.\n", SCIPconsGetName(cons), SCIPvarGetName(vars[i]));
            goto INTGRADESTIMATOR_CLEANUP;
         }
         x[i] = ub;
      }
      else
      {
         if( SCIPisInfinity(scip, ub) )
            x[i] = lb;
         else
            x[i] = (2.0*x[i] < lb+ub) ? lb : ub;
      }
      SCIPintervalSetBounds(&box[i],
         -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -MIN(lb, ub)),
         +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  MAX(lb, ub)));
   }

   /* compile expression if evaluated the first time; can only happen if newx is FALSE */
   if( newx && SCIPexprtreeGetInterpreterData(exprtree) == NULL )
   {
      SCIP_CALL( SCIPexprintCompile(exprint, exprtree) );
   }

   /* evaluate in reference point */
   SCIP_CALL( SCIPexprintEval(exprint, exprtree, x, &val) );
   if( !finite(val) )
   {
      SCIPdebugMessage("Got nonfinite function value from evaluation of constraint %s tree %d. skipping interval gradient estimator.\n", SCIPconsGetName(cons), exprtreeidx);
      goto INTGRADESTIMATOR_CLEANUP;
   }

   treecoef = consdata->nonlincoefs[exprtreeidx];
   val *= treecoef;
   constant = val;

   /* compute interval gradient */
   SCIP_CALL( SCIPallocBufferArray(scip, &intgrad, nvars) );
   SCIP_CALL( SCIPexprintGradInt(exprint, exprtree, INTERVALINFTY, box, TRUE, &intval, intgrad) );
   SCIPintervalMulScalar(INTERVALINFTY, &intval, intval, treecoef);

   /* printf("nvars %d side %d xref = %g x = [%g,%g] intval = [%g,%g] intgrad = [%g,%g]\n", nvars, side, x[0],
      box[0].inf, box[0].sup, intval.inf, intval.sup, intgrad[0].inf, intgrad[0].sup); */

   /* compute coefficients and constant */
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      val = x[i];
      lb  = SCIPintervalGetInf(box[i]);
      ub  = SCIPintervalGetSup(box[i]);

      SCIPintervalMulScalar(INTERVALINFTY, &intgrad[i], intgrad[i], treecoef);

      if( (overestimate && val == ub) ||  /*lint !e777*/
         (!overestimate && val == lb) )   /*lint !e777*/
         coefs[i] = SCIPintervalGetInf(intgrad[i]);
      else
         coefs[i] = SCIPintervalGetSup(intgrad[i]);

      if( SCIPisZero(scip, coefs[i]) )
         continue;

      if( SCIPisInfinity(scip, -coefs[i]) || SCIPisInfinity(scip, coefs[i]) )
      {
         SCIPdebugMessage("skip intgrad estimator because of infinite interval bound\n");
         goto INTGRADESTIMATOR_CLEANUP;
      }

      constant -= coefs[i] * val;
   }

   /* add interval gradient estimator to row */
   if( !SCIPisInfinity(scip, -SCIProwGetLhs(row)) )
   {
      SCIP_CALL( SCIPchgRowLhs(scip, row, SCIProwGetLhs(row) - constant) );
   }
   if( !SCIPisInfinity(scip,  SCIProwGetRhs(row)) )
   {
      SCIP_CALL( SCIPchgRowRhs(scip, row, SCIProwGetRhs(row) - constant) );
   }
   SCIP_CALL( SCIPaddVarsToRow(scip, row, nvars, vars, coefs) );

 INTGRADESTIMATOR_CLEANUP:
   SCIPfreeBufferArrayNull(scip, &box);
   SCIPfreeBufferArrayNull(scip, &intgrad);
   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** generates a cut based on linearization (if convex), secant (if concave), or intervalgradient (if indefinite)
 */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRINT*         exprint,            /**< expression interpreter */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real**           ref,                /**< reference point for each exprtree, or NULL if sol should be used */
   SCIP_SOL*             sol,                /**< reference solution where cut should be generated, or NULL if LP solution should be used */
   SCIP_Bool             newsol,             /**< whether the last evaluation of the expression with the expression interpreter was not at sol */
   SCIP_SIDETYPE         side,               /**< for which side a cut should be generated */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Real             maxrange,           /**< maximal range allowed */
   SCIP_Bool             expensivecurvchecks,/**< whether also expensive checks should be executed */
   SCIP_Bool             assumeconvex        /**< whether to assume convexity in inequalities */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;
   SCIP_Real* x;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   SCIPdebugMessage("constructing cut for %s hand side of constraint <%s>\n", side == SCIP_SIDETYPE_LEFT ? "left" : "right", SCIPconsGetName(cons));

   SCIP_CALL( checkCurvature(scip, cons, expensivecurvchecks, assumeconvex) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nexprtrees == 0 )
   {
      /* if we are actually linear, add the constraint as row to the LP */
      SCIP_CALL( SCIPcreateEmptyRow(scip, row, SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPconsIsLocal(cons), FALSE , TRUE) );
      SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, SCIPconsGetName(cons),
         side == SCIP_SIDETYPE_LEFT  ? consdata->lhs : -SCIPinfinity(scip),
         side == SCIP_SIDETYPE_RIGHT ? consdata->rhs :  SCIPinfinity(scip),
         !(side == SCIP_SIDETYPE_LEFT  && (consdata->curvature & SCIP_EXPRCURV_CONCAVE)) &&
         !(side == SCIP_SIDETYPE_RIGHT && (consdata->curvature & SCIP_EXPRCURV_CONVEX )),
         FALSE, TRUE) );

   if( ref == NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &x, SCIPexprtreeGetNVars(consdata->exprtrees[0])) );
   }

   success = TRUE;
   for( i = 0; i < consdata->nexprtrees; ++i )
   {
      if( ref == NULL )
      {
         SCIP_CALL( SCIPreallocBufferArray(scip, &x, SCIPexprtreeGetNVars(consdata->exprtrees[i])) );
         SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPexprtreeGetNVars(consdata->exprtrees[i]), SCIPexprtreeGetVars(consdata->exprtrees[i]), x) );
      }
      else
      {
         x = ref[i];
      }

      if( (side == SCIP_SIDETYPE_LEFT && (consdata->curvatures[i] & SCIP_EXPRCURV_CONCAVE)) ||
         (side == SCIP_SIDETYPE_RIGHT && (consdata->curvatures[i] & SCIP_EXPRCURV_CONVEX )) )
      {
         SCIP_CALL( addLinearization(scip, exprint, cons, i, x, newsol, *row, &success) );
      }
      else if( (side == SCIP_SIDETYPE_LEFT  && (consdata->curvatures[i] & SCIP_EXPRCURV_CONVEX)) ||
         (      side == SCIP_SIDETYPE_RIGHT && (consdata->curvatures[i] & SCIP_EXPRCURV_CONCAVE)) )
      {
         switch( SCIPexprtreeGetNVars(consdata->exprtrees[i]) )
         {
         case 1:
         {
            SCIP_CALL( addConcaveEstimatorUnivariate(scip, cons, i, *row, &success) );
            break;
         }
         case 2:
         {
            SCIP_CALL( addConcaveEstimatorBivariate(scip, cons, i, x, *row, &success) );
            break;
         }
         default:
         {
            SCIP_CALL( addConcaveEstimatorMultivariate(scip, cons, i, x, *row, &success) );
            break;
         }
         }
         if( !success )
         {
            SCIPdebugMessage("failed to generate polyhedral estimator for concave function, fall back to intervalgradient cut\n");
            SCIP_CALL( addIntervalGradientEstimator(scip, exprint, cons, i, x, newsol, side == SCIP_SIDETYPE_LEFT, *row, &success) );
         }
      }
      else
      {
         SCIP_CALL( addIntervalGradientEstimator(scip, exprint, cons, i, x, newsol, side == SCIP_SIDETYPE_LEFT, *row, &success) );
      }

      if( !success )
         break;
   }

   if( ref == NULL )
   {
      SCIPfreeBufferArray(scip, &x);
   }

   /* check numerics */
   if( success )
   {
      SCIP_Real mincoef;
      SCIP_Real maxcoef;

      mincoef = SCIPgetRowMinCoef(scip, *row);
      maxcoef = SCIPgetRowMaxCoef(scip, *row);

      assert(SCIPgetStage(scip) == SCIP_STAGE_SOLVING);
      mincoef = MIN(mincoef, consdata->lincoefsmin);
      maxcoef = MAX(maxcoef, consdata->lincoefsmax);

      while( maxcoef / mincoef > maxrange )
      {
         SCIP_VAR* var;
         SCIP_Real coef;
         SCIP_Real constant;
         int j;

         /* if range of coefficients is bad, find very small coefficients (from nonlinear vars) and make them zero */
         SCIPdebugMessage("cut coefficients for constraint <%s> have very large range: mincoef = %g maxcoef = %g\n", SCIPconsGetName(cons), mincoef, maxcoef);

         /* if minimal coefficient is given by linear var, then give up (probably the maximal coefficient is the problem) */
         if( mincoef == consdata->lincoefsmin )  /*lint !e777*/
         {
            SCIPdebugMessage("could not eliminate small coefficient, since it comes from linear part\n");
            break;
         }

         constant = 0.0;
         for( j = 0; j < SCIProwGetNNonz(*row); ++j )
         {
            coef = SCIProwGetVals(*row)[j];
            if( !SCIPisEQ(scip, REALABS(coef), mincoef) )
               continue;

            var = SCIPcolGetVar(SCIProwGetCols(*row)[j]);
            assert(var != NULL);

            /* try to eliminate coefficient with minimal absolute value by weakening cut and try again */
            if( ((coef > 0.0 && side == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && side == SCIP_SIDETYPE_LEFT)) &&
               !SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
            {
               SCIPdebugMessage("eliminate coefficient %g for <%s> = %g [%g, %g]\n", coef, SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

               constant += coef * (SCIProwIsLocal(*row) ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var));
               SCIP_CALL( SCIPaddVarToRow(scip, *row, var, -coef) );
               continue;
            }

            if( ((coef < 0.0 && side == SCIP_SIDETYPE_RIGHT) || (coef > 0.0 && side == SCIP_SIDETYPE_LEFT)) &&
               !SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
            {
               SCIPdebugMessage("eliminate coefficient %g for <%s> = %g [%g, %g]\n", coef, SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));

               constant += coef * (SCIProwIsLocal(*row) ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var));
               SCIP_CALL( SCIPaddVarToRow(scip, *row, var, -coef) );
               continue;
            }

            break;
         }

         if( j < SCIProwGetNNonz(*row) )
         {
            SCIPdebugMessage("could not eliminate small coefficient\n");
            success = FALSE;
            break;
         }

         if( side == SCIP_SIDETYPE_LEFT )
         {
            SCIP_CALL( SCIPchgRowLhs(scip, *row, SCIProwGetLhs(*row) - constant) );
         }
         else
         {
            SCIP_CALL( SCIPchgRowRhs(scip, *row, SCIProwGetRhs(*row) - constant) );
         }

         /* update min/max coefficient */
         mincoef = SCIPgetRowMinCoef(scip, *row);
         maxcoef = SCIPgetRowMaxCoef(scip, *row);

         mincoef = MIN(mincoef, consdata->lincoefsmin);
         maxcoef = MAX(maxcoef, consdata->lincoefsmax);
      };

      /* avoid numerically very bad cuts */
      if( maxcoef / mincoef > maxrange )
      {
         SCIPdebugMessage("drop row for constraint <%s> because range of coefficients is too large: mincoef = %g, maxcoef = %g -> range = %g\n",
            SCIPconsGetName(cons), mincoef, maxcoef, maxcoef / mincoef);
         success = FALSE;
      }
   }

   if( success &&
      ((  side == SCIP_SIDETYPE_LEFT  && SCIPisInfinity(scip, -SCIProwGetLhs(*row))) ||
         (side == SCIP_SIDETYPE_RIGHT && SCIPisInfinity(scip,  SCIProwGetRhs(*row)))) )
   {
      SCIPdebugMessage("drop row for constraint <%s> because of very large side: %g\n", SCIPconsGetName(cons), side == SCIP_SIDETYPE_LEFT ? -SCIProwGetLhs(*row) : SCIProwGetRhs(*row));
      success = FALSE;
   }

   if( !success )
   {
      SCIP_CALL( SCIPreleaseRow(scip, row) );
      return SCIP_OKAY;
   }

   /* add coefficients for linear variables */
   SCIP_CALL( SCIPaddVarsToRow(scip, *row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );

   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 *
 *  assumes that constraint violations have been computed
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Bool             newsol,             /**< have the constraints just been evaluated at this point? */
   SCIP_Real             minefficacy,        /**< minimal efficacy of a cut if it should be added to the LP */
   SCIP_Bool             convexalways,       /**< whether to ignore minefficacy criteria for a convex constraint (and use feastol instead) */
   SCIP_RESULT*          result,             /**< result of separation */
   SCIP_Real*            bestefficacy        /**< buffer to store best efficacy of a cut that was added to the LP, if found; or NULL if not of interest */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          efficacy;
   SCIP_Real          feasibility;
   SCIP_Real          norm;
   SCIP_SIDETYPE      violside;
   int                c;
   SCIP_ROW*          row;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nusefulconss <= nconss);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( bestefficacy != NULL )
      *bestefficacy = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         /* we are not feasible anymore */
         if( *result == SCIP_FEASIBLE )
            *result = SCIP_DIDNOTFIND;

         violside = SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT;

         /* generate cut
          * if function is defined at sol (activity<infinity) and constraint is violated, then expression interpreter should have evaluated at sol to get gradient before
          */
         SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], NULL, sol, newsol || SCIPisInfinity(scip, consdata->activity), violside, &row, conshdlrdata->cutmaxrange, conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );

         if( row == NULL ) /* failed to generate cut */
            continue;

         if( sol == NULL )
            feasibility = SCIPgetRowLPFeasibility(scip, row);
         else
            feasibility = SCIPgetRowSolFeasibility(scip, row, sol);
         norm = SCIProwGetNorm(row);

         /* in difference to SCIPgetCutEfficacy, we scale by norm only if the norm is > 1.0
          * this avoid finding cuts efficiant which are only very slightly violated
          * CPLEX does not seem to scale row coefficients up too
          */
         if( norm > 1.0 )
            efficacy = -feasibility / norm;
         else
            efficacy = -feasibility;

         if( SCIPisGT(scip, efficacy, minefficacy) ||
            (convexalways &&
               ( ( violside == SCIP_SIDETYPE_RIGHT && (consdata->curvature & SCIP_EXPRCURV_CONVEX )) ||
                  (violside == SCIP_SIDETYPE_LEFT  && (consdata->curvature & SCIP_EXPRCURV_CONCAVE)) ) &&
               SCIPisGT(scip, efficacy, SCIPfeastol(scip))
               )
            )
         {
            /* cut cuts off solution */
            SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE /* forcecut */) );
            *result = SCIP_SEPARATED;
            SCIP_CALL( SCIPresetConsAge(scip, conss[c]) );
            SCIPdebugMessage("add cut with efficacy %g for constraint <%s> violated by %g\n", efficacy, SCIPconsGetName(conss[c]), MAX(consdata->lhsviol, consdata->rhsviol));
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            if( bestefficacy != NULL && efficacy > *bestefficacy )
               *bestefficacy = efficacy;
         }
         else
         {
            SCIPdebugMessage("drop cut since efficacy %g is too small (< %g)\n", efficacy, minefficacy);
         }

         SCIP_CALL( SCIPreleaseRow (scip, &row) );
      }

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */
      if( c >= nusefulconss && *result == SCIP_SEPARATED )
         break;
   }

   return SCIP_OKAY;
}

/** adds linearizations cuts for convex constraints w.r.t. a given reference point to cutpool and sepastore
 * if separatedlpsol is not NULL, then cuts that separate the LP solution are added to the sepastore too
 */
static
SCIP_RETCODE addLinearizationCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             ref,                /**< reference point where to linearize, or NULL for LP solution */
   SCIP_Bool*            separatedlpsol,     /**< buffer to store whether a cut that separates the current LP solution was found, or NULL if not of interest */
   SCIP_Real             minefficacy         /**< minimal efficacy of a cut when checking for separation of LP solution */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( separatedlpsol != NULL )
      *separatedlpsol = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      if( SCIPconsIsLocal(conss[c]) )  /*lint !e613*/
         continue;

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* if we cannot linearize, then skip constraint */
      if( (!(consdata->curvature & SCIP_EXPRCURV_CONVEX)  || SCIPisInfinity(scip,  consdata->rhs)) &&
         ( !(consdata->curvature & SCIP_EXPRCURV_CONCAVE) || SCIPisInfinity(scip, -consdata->lhs)) )
         continue;

      SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], NULL, ref, TRUE,
            consdata->curvature & SCIP_EXPRCURV_CONVEX ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT,
            &row, conshdlrdata->cutmaxrange, FALSE, FALSE) );  /*lint !e613*/

      if( row == NULL )
         continue;

      /* if caller wants, then check if cut separates LP solution and add to sepastore if so */
      if( separatedlpsol != NULL )
      {
         SCIP_Real feasibility;
         SCIP_Real norm;

         feasibility = SCIPgetRowLPFeasibility(scip, row);
         norm = SCIPgetRowMaxCoef(scip, row);

         if( -feasibility / MAX(1.0, norm) >= minefficacy )
         {
            *separatedlpsol = TRUE;
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );
         }
      }

      if( !SCIProwIsLocal(row) )
      {
         SCIP_CALL( SCIPaddPoolCut(scip, row) );
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** processes the event that a new primal solution has been found */
static
SCIP_DECL_EVENTEXEC(processNewSolutionEvent)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS**    conss;
   int            nconss;
   SCIP_SOL*      sol;

   assert(scip != NULL);
   assert(event != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   assert((SCIPeventGetType(event) & SCIP_EVENTTYPE_SOLFOUND) != 0);

   conshdlr = (SCIP_CONSHDLR*)eventdata;

   nconss = SCIPconshdlrGetNConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* we are only interested in solution coming from some heuristic other than trysol, but not from the tree
    * the reason for ignoring trysol solutions is that they may come from an NLP solve in sepalp, where we already added linearizations,
    * or are from the tree, but postprocessed via proposeFeasibleSolution
    */
   if( SCIPsolGetHeur(sol) == NULL || SCIPsolGetHeur(sol) == conshdlrdata->trysolheur )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   SCIPdebugMessage("catched new sol event %x from heur <%s>; have %d conss\n", SCIPeventGetType(event), SCIPheurGetName(SCIPsolGetHeur(sol)), nconss);

   SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, sol, NULL, 0.0) );

   return SCIP_OKAY;
}

/** registers unfixed variables in nonlinear terms of violated constraints as external branching candidates */
static
SCIP_RETCODE registerBranchingVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int c;
   int i;
   int j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->nexprtrees == 0 )
         continue;

      /* do not branch on violation of convex constraint */
      if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || (consdata->curvature & SCIP_EXPRCURV_CONCAVE)) &&
         ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || (consdata->curvature & SCIP_EXPRCURV_CONVEX )) )
         continue;
      SCIPdebugMessage("cons <%s> violation: %g %g  curvature: %s\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, SCIPexprcurvGetName(consdata->curvature));

      for( i = 0; i < consdata->nexprtrees; ++i )
      {
         /* skip convex summands */
         if( (!SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || (consdata->curvatures[i] & SCIP_EXPRCURV_CONCAVE)) &&
            ( !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) || (consdata->curvatures[i] & SCIP_EXPRCURV_CONVEX )) )
            continue;

         for( j = 0; j < SCIPexprtreeGetNVars(consdata->exprtrees[i]); ++j )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[i])[j];
            assert(var != NULL);

            if( SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            {
               SCIPdebugMessage("ignore fixed variable <%s>[%g, %g], width %g\n", SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetUbLocal(var) - SCIPvarGetLbLocal(var));
               continue;
            }

            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++*nnotify;
         }
      }
   }

   SCIPdebugMessage("registered %d branching candidates\n", *nnotify);

   return SCIP_OKAY;
}

/** registers a nonlinear variable from a violated constraint as branching candidate that has a large absolute value in the LP relaxation */
static
SCIP_RETCODE registerLargeLPValueVariableForBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_VAR**            brvar               /**< buffer to store branching variable */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_VAR*           var;
   SCIP_Real           val;
   SCIP_Real           brvarval;
   int i;
   int j;
   int c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   *brvar = NULL;
   brvarval = -1.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]); ++i )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[j])[i];
            /* do not propose fixed variables */
            if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
               continue;
            val = SCIPgetSolVal(scip, NULL, var);
            if( REALABS(val) > brvarval )
            {
               brvarval = ABS(val);
               *brvar = var;
            }
         }
      }
   }

   if( *brvar != NULL )
   {
      SCIP_CALL( SCIPaddExternBranchCand(scip, *brvar, brvarval, SCIP_INVALID) );
   }

   return SCIP_OKAY;
}

/** replaces violated nonlinear constraints where all nonlinear variables are fixed by linear constraints
 * only adds constraint if it is violated in current solution
 */
static
SCIP_RETCODE replaceViolatedByLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            addedcons           /**< buffer to store whether a linear constraint was added */
   )
{
   SCIP_CONS*          cons;
   SCIP_CONSDATA*      consdata;
   SCIP_Real           lhs;
   SCIP_Real           rhs;
   SCIP_RESULT         checkresult;
   int c;
   int i;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);
   assert(addedcons != NULL);

   *addedcons = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      lhs = consdata->lhs;
      rhs = consdata->rhs;

      for( i = 0; i < consdata->nexprtrees; ++i )
      {
         SCIP_INTERVAL nonlinactivity;

         SCIP_CALL( SCIPevalExprtreeLocalBounds(scip, consdata->exprtrees[i], INTERVALINFTY, &nonlinactivity) );
         assert(SCIPintervalGetInf(nonlinactivity) > -INTERVALINFTY);
         assert(SCIPintervalGetSup(nonlinactivity) <  INTERVALINFTY);
         SCIPintervalMulScalar(INTERVALINFTY, &nonlinactivity, nonlinactivity, consdata->nonlincoefs[i]);

         if( !SCIPisInfinity(scip, -lhs) )
            lhs -= SCIPintervalGetSup(nonlinactivity);

         if( !SCIPisInfinity(scip,  rhs) )
            rhs -= SCIPintervalGetInf(nonlinactivity);
      }

      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, SCIPconsGetName(conss[c]),
            consdata->nlinvars, consdata->linvars, consdata->lincoefs,
            lhs, rhs,
            SCIPconsIsInitial(conss[c]), SCIPconsIsSeparated(conss[c]), SCIPconsIsEnforced(conss[c]),
            SCIPconsIsChecked(conss[c]), SCIPconsIsPropagated(conss[c]),  TRUE,
            SCIPconsIsModifiable(conss[c]), SCIPconsIsDynamic(conss[c]), SCIPconsIsRemovable(conss[c]),
            SCIPconsIsStickingAtNode(conss[c])) );

      SCIPdebugMessage("replace violated nonlinear constraint <%s> by linear constraint after all nonlinear vars have been fixed\n", SCIPconsGetName(conss[c]) );
      SCIPdebug( SCIPprintCons(scip, conss[c], NULL) );
      SCIPdebug( SCIPprintCons(scip, cons, NULL) );

      SCIP_CALL( SCIPcheckCons(scip, cons, NULL, FALSE, FALSE, FALSE, &checkresult) );

      if( checkresult != SCIP_INFEASIBLE )
      {
         SCIPdebugMessage("linear constraint is feasible, thus do not add\n");
      }
      else
      {
         SCIP_CALL( SCIPaddConsLocal(scip, cons, NULL) );
         *addedcons = TRUE;
      }
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
   }

   return SCIP_OKAY;
}

/* tightens a lower bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where we currently propagate, or NULL if tightening is from expression graph */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new lower bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds            /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(bnd > -INTERVALINFTY);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   if( SCIPisInfinity(scip, bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }

   /* new lower bound is very low (between -INTERVALINFTY and -SCIPinfinity()) */
   if( SCIPisInfinity(scip, -bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarLb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarLb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMessage("%sfound constraint <%s> infeasible due to tightened lower bound %g for variable <%s>\n", SCIPinProbing(scip) ? "in probing " : "", cons != NULL ? SCIPconsGetName(cons) : "??", bnd, SCIPvarGetName(var));  /*lint !e585*/
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMessage("%stightened lower bound of variable <%s> in constraint <%s> to %.20g\n", SCIPinProbing(scip) ? "in probing " : "", SCIPvarGetName(var), cons != NULL ? SCIPconsGetName(cons) : "??", bnd);  /*lint !e585*/
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/* tightens an upper bound on a variable and checks the result */
static
SCIP_RETCODE propagateBoundsTightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint where we currently propagate, or NULL if tightening is from expression graph */
   SCIP_VAR*             var,                /**< variable which domain we might reduce */
   SCIP_Real             bnd,                /**< new upper bound for variable */
   SCIP_RESULT*          result,             /**< result to update if there was a tightening or cutoff */
   int*                  nchgbds             /**< counter to increase if a bound was tightened */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip != NULL);
   assert(bnd < INTERVALINFTY);
   assert(var != NULL);
   assert(result != NULL);
   assert(*result == SCIP_DIDNOTFIND || *result == SCIP_REDUCEDDOM);
   assert(nchgbds != NULL);

   if( SCIPisInfinity(scip, -bnd) )
   { /* domain will be outside [-infty, +infty] -> declare node infeasible */
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }

   /* new upper bound is very high (between SCIPinfinity() and INTERVALINFTY) */
   if( SCIPisInfinity(scip, bnd) )
      return SCIP_OKAY;

   bnd = SCIPadjustedVarUb(scip, var, bnd);
   SCIP_CALL( SCIPtightenVarUb(scip, var, bnd, FALSE, &infeas, &tightened) );
   if( infeas )
   {
      SCIPdebugMessage("%sfound constraint <%s> infeasible due to tightened upper bound %g for variable <%s>\n", SCIPinProbing(scip) ? "in probing " : "", cons != NULL ? SCIPconsGetName(cons) : "??", bnd, SCIPvarGetName(var));  /*lint !e585*/
      *result = SCIP_CUTOFF;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
      return SCIP_OKAY;
   }
   if( tightened )
   {
      SCIPdebugMessage("%stightened upper bound of variable <%s> in constraint <%s> to %g\n", SCIPinProbing(scip) ? "in probing " : "", SCIPvarGetName(var), cons != NULL ? SCIPconsGetName(cons) : "??", bnd);  /*lint !e585*/
      ++*nchgbds;
      *result = SCIP_REDUCEDDOM;
      if( cons != NULL )
      {
         SCIP_CALL( SCIPresetConsAge(scip, cons) );
      }
   }

   return SCIP_OKAY;
}

/** tightens bounds of linear variables in a single nonlinear constraint */
static
SCIP_RETCODE propagateBoundsCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation call */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   SCIP_Bool*            redundant           /**< buffer where to store whether constraint has been found to be redundant */
   )
{  /*lint --e{666}*/
   SCIP_CONSDATA*     consdata;
   SCIP_INTERVAL      consbounds;    /* lower and upper bounds of constraint */
   SCIP_INTERVAL      consactivity;  /* activity of linear plus nonlinear part */
   SCIP_VAR*          var;
   SCIP_INTERVAL      rhs;           /* right hand side of nonlinear equation */
   SCIP_ROUNDMODE     roundmode;
   SCIP_Real          bnd;
   int                i;
   SCIP_INTERVAL      nonlinactivity;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTRUN;
   *redundant = FALSE;

   if( consdata->ispropagated )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMessage("start linear vars domain propagation for constraint <%s>\n", SCIPconsGetName(cons));

   consdata->ispropagated = TRUE;

   /* make sure we have activity of linear term */
   consdataUpdateLinearActivity(scip, consdata);
   assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777*/
   assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777*/
   assert(consdata->minlinactivityinf >= 0);
   assert(consdata->maxlinactivityinf >= 0);
   assert(consdata->exprgraphnode != NULL || consdata->nexprtrees == 0);

   /* get activity of nonlinear part, should have been updated in propagateBounds */
   if( consdata->exprgraphnode != NULL )
   {
      nonlinactivity = SCIPexprgraphGetNodeBounds(consdata->exprgraphnode);
   }
   else
   {
      SCIPintervalSet(&nonlinactivity, 0.0);
   }
   assert(!SCIPintervalIsEmpty(nonlinactivity) );

   /* @todo adding SCIPepsilon may be sufficient? */
   SCIPintervalSetBounds(&consbounds,
      -infty2infty(SCIPinfinity(scip), INTERVALINFTY, -consdata->lhs + SCIPfeastol(scip)),
      +infty2infty(SCIPinfinity(scip), INTERVALINFTY,  consdata->rhs + SCIPfeastol(scip)));

   /* check redundancy and infeasibility */
   SCIPintervalSetBounds(&consactivity, consdata->minlinactivityinf > 0 ? -INTERVALINFTY : consdata->minlinactivity, consdata->maxlinactivityinf > 0 ? INTERVALINFTY : consdata->maxlinactivity);
   SCIPintervalAdd(INTERVALINFTY, &consactivity, consactivity, nonlinactivity);
   if( SCIPintervalIsSubsetEQ(INTERVALINFTY, consactivity, consbounds) )
   {
      SCIPdebugMessage("found constraint <%s> to be redundant: sides: [%g, %g], activity: [%g, %g]\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity));
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   if( SCIPintervalAreDisjoint(consbounds, consactivity) )
   {
      SCIPdebugMessage("found constraint <%s> to be infeasible; sides: [%g, %g], activity: [%g, %g], infeas: %.20g\n",
         SCIPconsGetName(cons), consdata->lhs, consdata->rhs, SCIPintervalGetInf(consactivity), SCIPintervalGetSup(consactivity),
         MAX(consdata->lhs - SCIPintervalGetSup(consactivity), SCIPintervalGetInf(consactivity) - consdata->rhs));
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* propagate linear part in rhs = consbounds - nonlinactivity (use the one from consdata, since that includes infinities) */
   SCIPintervalSub(INTERVALINFTY, &rhs, consbounds, nonlinactivity);
   if( !SCIPintervalIsEntire(INTERVALINFTY, rhs) )
   {
      SCIP_Real coef;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         coef = consdata->lincoefs[i];
         var  = consdata->linvars[i];

         /* skip fixed variables
          * @todo is that a good or a bad idea?
          *   we can't expect much more tightening, but may detect infeasiblity, but shouldn't the check on the constraints activity detect that?
          */
         if( SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            continue;

         if( coef > 0.0 )
         {
            if( SCIPintervalGetSup(rhs) < INTERVALINFTY )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the upper bound on var x */
               if( consdata->minlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
                  /* tighten upper bound on x to (rhs.sup - (minlinactivity - coef * xlb)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd += coef * SCIPvarGetLbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the minimal linear activity equal -infinity, so
                   * we tighten upper bound on x to just (rhs.sup - minlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = SCIPintervalGetSup(rhs);
                  bnd -= consdata->minlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetInf(rhs) > -INTERVALINFTY )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the lower bound on var x */
               if( consdata->maxlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* tighten lower bound on x to (rhs.inf - (maxlinactivity - coef * xub)) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd += coef * SCIPvarGetUbLocal(var);
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (rhs.inf - maxlinactivity) / coef */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = SCIPintervalGetInf(rhs);
                  bnd -= consdata->maxlinactivity;
                  bnd /= coef;
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is +infinity and x is not solely responsible for this */
            }
         }
         else
         {
            assert(coef < 0.0 );
            if( SCIPintervalGetInf(rhs) > -INTERVALINFTY )
            {
               assert(consdata->maxlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the upper bound on var x */
               if( consdata->maxlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetLbLocal(var)));
                  /* compute upper bound on x to (maxlinactivity - coef * xlb) - rhs.inf / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd += (-coef) * SCIPvarGetLbLocal(var);
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->maxlinactivityinf == 1 && SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal infinity, so
                   * we tighten upper bound on x to just (maxlinactivity - rhs.inf) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeUpwards();
                  bnd  = consdata->maxlinactivity;
                  bnd -= SCIPintervalGetInf(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarUb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the maximal activity is infinity and x is not solely responsible for this */
            }

            if( SCIPintervalGetSup(rhs) < INTERVALINFTY )
            {
               assert(consdata->minlinactivity != SCIP_INVALID);  /*lint !e777*/
               /* try to tighten the lower bound on var x */
               if( consdata->minlinactivityinf == 0 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  /* compute lower bound on x to (minlinactivity - coef * xub) - rhs.sup / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd += (-coef) * SCIPvarGetUbLocal(var);
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               else if( consdata->minlinactivityinf == 1 && SCIPisInfinity(scip, SCIPvarGetUbLocal(var)) )
               {
                  /* x was the variable that made the maximal linear activity equal -infinity, so
                   * we tighten lower bound on x to just (minlinactivity - rhs.sup) / (-coef) */
                  roundmode = SCIPintervalGetRoundingMode();
                  SCIPintervalSetRoundingModeDownwards();
                  bnd  = consdata->minlinactivity;
                  bnd -= SCIPintervalGetSup(rhs);
                  bnd /= (-coef);
                  SCIPintervalSetRoundingMode(roundmode);
                  SCIP_CALL( propagateBoundsTightenVarLb(scip, cons, var, bnd, result, nchgbds) );
                  if( *result == SCIP_CUTOFF )
                     break;
               }
               /* otherwise the minimal activity is -infinity and x is not solely responsible for this */
            }
         }
      }
      if( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** propagate constraints sides minus linear activity into nonlinear variables */
static
SCIP_RETCODE propagateConstraintSides(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation calls */
   int*                  nchgbds             /**< buffer where to add the number of changed bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int         nvars;
   SCIP_VAR**  vars;
   SCIP_EXPRGRAPHNODE** varnodes;
   SCIP_INTERVAL bounds;
   SCIP_Bool   cutoff;
   SCIP_ROUNDMODE roundmode;
   int         c;
   int         i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(result != NULL);
   assert(nchgbds != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   SCIPdebugMessage("start backward propagation in expression graph\n");

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_propconss1.dot", "w");
      SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, file, NULL) );
      fclose(file);
   }
#endif

   /* put constraint sides less linear activity into expression graph nodes
    * also add a [-feastol,feastol] range around constraint sides to cope with numerics */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->exprgraphnode == NULL )
         continue;

      /* skip (just) deleted or disabled constraints */
      if( SCIPconsIsDeleted(conss[c]) || !SCIPconsIsEnabled(conss[c]) )
         continue;

      roundmode = SCIPintervalGetRoundingMode();
      SCIPintervalSetRoundingModeDownwards();

      if( !SCIPisInfinity(scip, -consdata->lhs) && consdata->maxlinactivityinf == 0 )
         bounds.inf = consdata->lhs - consdata->maxlinactivity - SCIPfeastol(scip);
      else
         bounds.inf = -INTERVALINFTY;

      if( !SCIPisInfinity(scip,  consdata->rhs) && consdata->minlinactivityinf == 0 )
         bounds.sup = SCIPintervalNegateReal(consdata->minlinactivity - consdata->rhs - SCIPfeastol(scip));
      else
         bounds.sup =  INTERVALINFTY;

      SCIPintervalSetRoundingMode(roundmode);

      SCIPexprgraphTightenNodeBounds(conshdlrdata->exprgraph, consdata->exprgraphnode, bounds, BOUNDTIGHTENING_MINSTRENGTH, &cutoff);

      if( cutoff )
      {
         SCIPdebugMessage("found constraint <%s> infeasible%s\n", SCIPconsGetName(conss[c]), SCIPinProbing(scip) ? " in probing" : "");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   /* compute bound tightenings for nonlinear variables */
   SCIPexprgraphPropagateNodeBounds(conshdlrdata->exprgraph, INTERVALINFTY, BOUNDTIGHTENING_MINSTRENGTH, &cutoff);

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_propconss2.dot", "w");
      SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, file, NULL) );
      fclose(file);
   }
#endif

   if( cutoff )
   {
      SCIPdebugMessage("backward propagation found problem infeasible%s\n", SCIPinProbing(scip) ? " in probing" : "");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* put tighter bounds into variables */
   nvars = SCIPexprgraphGetNVars(conshdlrdata->exprgraph);
   vars  = (SCIP_VAR**)SCIPexprgraphGetVars(conshdlrdata->exprgraph);
   varnodes = SCIPexprgraphGetVarNodes(conshdlrdata->exprgraph);

   /* put back new bounds into SCIP variables */
   for( i = 0; i < nvars && *result != SCIP_CUTOFF; ++i )
   {
      if( !SCIPisInfinity(scip, -SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(varnodes[i]))) )
      {
         SCIP_CALL( propagateBoundsTightenVarLb(scip, NULL, vars[i], SCIPintervalGetInf(SCIPexprgraphGetNodeBounds(varnodes[i])), result, nchgbds) );
      }
      if( *result != SCIP_CUTOFF && !SCIPisInfinity(scip,  SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(varnodes[i]))) )
      {
         SCIP_CALL( propagateBoundsTightenVarUb(scip, NULL, vars[i], SCIPintervalGetSup(SCIPexprgraphGetNodeBounds(varnodes[i])), result, nchgbds) );
      }
   }

   return SCIP_OKAY;
}

/** calls domain propagation for a set of constraints */
static
SCIP_RETCODE propagateBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result of the propagation calls */
   int*                  nchgbds,            /**< buffer where to add the the number of changed bounds */
   int*                  ndelconss           /**< buffer where to increase if a constraint was deleted (locally) due to redundancy */
   )
{
#ifndef NDEBUG
   SCIP_CONSDATA* consdata;
#endif
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RESULT propresult;
   SCIP_Bool   domainerror;
   SCIP_Bool   redundant;
   int         roundnr;
   SCIP_Bool   success;
   int         c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   if( nconss == 0 || conshdlrdata->ispropagated )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   roundnr = 0;
   do
   {
      success = FALSE;

      SCIPdebugMessage("starting domain propagation round %d for %d constraints\n", roundnr, nconss);

      conshdlrdata->ispropagated = TRUE;

      /* propagate variable bounds through expression graph
       * roundnr == 0 clears remainings from a previous backward propagation
       * @todo could give FALSE if no linear variable in the constraints had been relaxed since last time
       */
      SCIP_CALL( SCIPexprgraphPropagateVarBounds(conshdlrdata->exprgraph, INTERVALINFTY, roundnr == 0, &domainerror) );

#ifdef SCIP_OUTPUT
      {
         FILE* file;
         file = fopen("exprgraph_propvars.dot", "w");
         SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, file, NULL) );
         fclose(file);
      }
#endif

      if( domainerror )
      {
         SCIPdebugMessage("current bounds out of domain for some expression, do cutoff\n");
         *result = SCIP_CUTOFF;
         break;
      }

      /* check for redundancy and infeasibility of constraints, tighten bounds on linear variables */
      for( c = 0; c < nconss && *result != SCIP_CUTOFF; ++c )
      {
         assert(conss != NULL);
         if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c]) )
            continue;

#ifndef NDEBUG
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);
         assert(consdata->exprgraphnode == NULL || !SCIPintervalIsEmpty(SCIPexprgraphGetNodeBounds(consdata->exprgraphnode)));
#endif

         SCIP_CALL( propagateBoundsCons(scip, conshdlr, conss[c], &propresult, nchgbds, &redundant) );
         if( propresult != SCIP_DIDNOTFIND && propresult != SCIP_DIDNOTRUN )
         {
            *result = propresult;
            success = TRUE;
         }
         if( redundant )
         {
            SCIPdebugMessage("delete redundant constraint <%s> locally\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
            ++*ndelconss;
         }
      }

      /* propagate backward through expression graph */
      if( *result != SCIP_CUTOFF )
      {
         propresult = SCIP_DIDNOTFIND;
         SCIP_CALL( propagateConstraintSides(scip, conshdlr, conss, nconss, &propresult, nchgbds) );

         if( propresult != SCIP_DIDNOTFIND )
         {
            *result = propresult;
            success = TRUE;
         }
      }

   } while( success && *result != SCIP_CUTOFF && ++roundnr < conshdlrdata->maxproprounds );

   return SCIP_OKAY;
}

/* checks for a linear variable that can be increase or decreased without harming feasibility */
static
void consdataFindUnlockedLinearVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   int i;
   int poslock;
   int neglock;

   consdata->linvar_maydecrease = -1;
   consdata->linvar_mayincrease = -1;

   /* check for a linear variable that can be increase or decreased without harming feasibility
    * setup lincoefsmin, lincoefsmax */
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      /* compute locks of i'th linear variable */
      assert(consdata->lincoefs[i] != 0.0);
      if( consdata->lincoefs[i] > 0.0 )
      {
         poslock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
      }
      else
      {
         poslock = !SCIPisInfinity(scip,  consdata->rhs) ? 1 : 0;
         neglock = !SCIPisInfinity(scip, -consdata->lhs) ? 1 : 0;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - neglock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can decrease x without harming other constraints */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_maydecrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_maydecrease]) / consdata->lincoefs[consdata->linvar_maydecrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_maydecrease = i;
      }

      if( SCIPvarGetNLocksDown(consdata->linvars[i]) - poslock == 0 )
      {
         /* for a*x + q(y) \in [lhs, rhs], we can increase x without harm */
         /* if we have already one candidate, then take the one where the loss in the objective function is less */
         if( (consdata->linvar_mayincrease < 0) ||
            (SCIPvarGetObj(consdata->linvars[consdata->linvar_mayincrease]) / consdata->lincoefs[consdata->linvar_mayincrease] > SCIPvarGetObj(consdata->linvars[i]) / consdata->lincoefs[i]) )
            consdata->linvar_mayincrease = i;
      }
   }

#ifdef SCIP_DEBUG
   if( consdata->linvar_mayincrease >= 0 )
   {
      SCIPdebugMessage("may increase <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_mayincrease]));
   }
   if( consdata->linvar_maydecrease >= 0 )
   {
      SCIPdebugMessage("may decrease <%s> to become feasible\n", SCIPvarGetName(consdata->linvars[consdata->linvar_maydecrease]));
   }
#endif
}

/** Given a solution where every nonlinear constraint is either feasible or can be made feasible by
 * moving a linear variable, construct the corresponding feasible solution and pass it to the trysol heuristic.
 * The method assumes that this is always possible and that not all constraints are feasible already.
 */
static
SCIP_RETCODE proposeFeasibleSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to process */
   SCIP_Bool*            success             /**< buffer to store whether we succeeded to construct a solution that satisfies all provided constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_SOL* newsol;
   SCIP_VAR* var;
   int c;
   SCIP_Real viol;
   SCIP_Real norm;
   SCIP_Real delta;
   SCIP_Real gap;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(success != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->trysolheur != NULL);

   *success = FALSE;

   if( sol != NULL )
   {
      SCIP_CALL( SCIPcreateSolCopy(scip, &newsol, sol) );
   }
   else
   {
      SCIP_CALL( SCIPcreateLPSol(scip, &newsol, NULL) );
   }
   SCIP_CALL( SCIPunlinkSol(scip, newsol) );

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* recompute violation of solution in case solution has changed
       * get absolution violation and sign */
      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conshdlrdata->exprinterpreter, conss[c], newsol) );  /*lint !e613*/
         viol = consdata->lhs - consdata->activity;
      }
      else if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         SCIP_CALL( computeViolation(scip, conshdlrdata->exprinterpreter, conss[c], newsol) );  /*lint !e613*/
         viol = consdata->rhs - consdata->activity;
      }
      else
         continue; /* constraint is satisfied */

      assert(viol != 0.0);
      if( consdata->linvar_mayincrease >= 0 &&
         ((  viol > 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) ||
            (viol < 0.0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0)) )
      {
         /* have variable where increasing makes the constraint less violated */
         var = consdata->linvars[consdata->linvar_mayincrease];
         /* compute how much we would like to increase var */
         delta = viol / consdata->lincoefs[consdata->linvar_mayincrease];
         assert(delta > 0.0);
         /* if var has an upper bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
         {
            gap = SCIPvarGetUbGlobal(var) - SCIPgetSolVal(scip, newsol, var);
            delta = MIN(MAX(0.0, gap), delta);
         }
         if( SCIPisPositive(scip, delta) )
         {
            /* if variable is integral, round delta up so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPceil(scip, delta);

            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMessage("increase <%s> by %g to %g\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->lincoefs[consdata->linvar_mayincrease] * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      assert(viol != 0.0);
      if( consdata->linvar_maydecrease >= 0 &&
         ((  viol > 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) ||
            (viol < 0.0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0)) )
      {
         /* have variable where decreasing makes constraint less violated */
         var = consdata->linvars[consdata->linvar_maydecrease];
         /* compute how much we would like to decrease var */
         delta = viol / consdata->lincoefs[consdata->linvar_maydecrease];
         assert(delta < 0.0);
         /* if var has a lower bound, may need to reduce delta */
         if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
         {
            gap = SCIPgetSolVal(scip, newsol, var) - SCIPvarGetLbGlobal(var);
            delta = MAX(MIN(0.0, gap), delta);
         }
         if( SCIPisNegative(scip, delta) )
         {
            /* if variable is integral, round delta down so that it will still have an integer value */
            if( SCIPvarIsIntegral(var) )
               delta = SCIPfloor(scip, delta);
            SCIP_CALL( SCIPincSolVal(scip, newsol, var, delta) );
            SCIPdebugMessage("increase <%s> by %g to %g\n", SCIPvarGetName(var), delta, SCIPgetSolVal(scip, newsol, var));

            /* adjust constraint violation, if satisfied go on to next constraint */
            viol -= consdata->lincoefs[consdata->linvar_maydecrease] * delta;
            if( SCIPisZero(scip, viol) )
               continue;
         }
      }

      /* still here... so maybe we could not make constraint feasible due to variable bounds
       * check if we are feasible w.r.t. (relative) feasibility tolerance */
      SCIP_CALL( getGradientMaxElement(scip, conshdlrdata->exprinterpreter, conss[c], newsol, TRUE, &norm) );  /*lint !e613*/
      if( norm > 1.0 )
         viol /= norm;
      /* if still violated, we give up */
      if( SCIPisGT(scip, REALABS(viol), SCIPfeastol(scip)) )
         break;

      /* if objective value is not better than current upper bound, we give up */
      if( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) && !SCIPisSumLT(scip, SCIPgetSolTransObj(scip, newsol), SCIPgetUpperbound(scip)) )
         break;
   }

   /* if we have a solution that should satisfy all nonlinear constraints and has a better objective than the current upper bound,
    * then pass it to the trysol heuristic */
   if( c == nconss )
   {
      SCIPdebugMessage("pass solution with objective value %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );
      *success = TRUE;
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyNonlinear)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   /* assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0); */

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprinterpreter != NULL);
   assert(conshdlrdata->exprgraph != NULL);
   assert(SCIPexprgraphGetNVars(conshdlrdata->exprgraph) == 0);

   /* free expression graph */
   SCIP_CALL( SCIPexprgraphFree(&conshdlrdata->exprgraph) );

   /* free upgrade functions */
   for( i = 0; i < conshdlrdata->nnlconsupgrades; ++i )
   {
      assert(conshdlrdata->nlconsupgrades[i] != NULL);
      SCIPfreeMemory(scip, &conshdlrdata->nlconsupgrades[i]);
   }
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->nlconsupgrades);

   /* free expressions interpreter */
   SCIP_CALL( SCIPexprintFree(&conshdlrdata->exprinterpreter) );

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");

   /* reset counter, since we have a new problem */
   conshdlrdata->naddedreformconss = 0;

#ifdef SCIP_OUTPUT
   {
      FILE* file;
      file = fopen("exprgraph_init.dot", "w");
      SCIP_CALL( SCIPexprgraphPrintDot(conshdlrdata->exprgraph, file, NULL) );
      fclose(file);
   }
#endif

   return SCIP_OKAY;
}  /*lint !e715*/

/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;

   return SCIP_OKAY;
}  /*lint !e715*/


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
static
SCIP_DECL_CONSINITPRE(consInitpreNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   *result = SCIP_FEASIBLE;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* forget expression trees */
      assert(consdata->nexprtrees == 0 || consdata->exprgraphnode != NULL);
      SCIP_CALL( consdataSetExprtrees(scip, consdata, 0, NULL, NULL, FALSE) );
   }

   return SCIP_OKAY;
}


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          havegraphchange;
   SCIP_Bool          havechange;
   SCIP_Bool          domainerror;
   int i;
   int j;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;
   havegraphchange = FALSE;

   if( !conshdlrdata->isremovedfixings )
   {
      SCIP_CALL( removeFixedNonlinearVariables(scip, conshdlr) );
      assert(conshdlrdata->isremovedfixings);

      havegraphchange = TRUE;
   }

   /* if undefined expressions in exprgraph, then declare problem as infeasible */
   SCIP_CALL( SCIPexprgraphSimplify(conshdlrdata->exprgraph, SCIPepsilon(scip), conshdlrdata->maxexpansionexponent, &havechange, &domainerror) );
   SCIPdebugMessage("expression graph simplifier found %schange, domain error = %u\n", havechange ? "" : "no ", domainerror);
   havegraphchange |= havechange;

   if( domainerror )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !consdata->isremovedfixingslin )
      {
         SCIP_CALL( removeFixedLinearVariables(scip, conss[c]) );
      }

      if( !consdata->ispresolved || havegraphchange )
      {
         SCIP_CALL( splitOffLinearPart(scip, conshdlr, conss[c]) );
      }

      SCIP_CALL( mergeAndCleanLinearVars(scip, conss[c]) );

      assert(consdata->isremovedfixingslin);
      assert(consdata->linvarsmerged);
#ifndef NDEBUG
      for( i = 0; i < consdata->nlinvars; ++i )
         assert(SCIPvarIsActive(consdata->linvars[i]));
#endif

      if( consdata->exprgraphnode != NULL )
      {
         /* get expression trees from expression graph */
         SCIP_EXPRTREE** exprtrees;
         SCIP_Real* coefs;
         int nexprtrees;
         int exprtreessize;

         exprtreessize = SCIPexprgraphGetSumTreesNSummands(consdata->exprgraphnode);

         SCIP_CALL( SCIPallocBufferArray(scip, &exprtrees, exprtreessize) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coefs,     exprtreessize) );

         SCIP_CALL( SCIPexprgraphGetSumTrees(conshdlrdata->exprgraph, consdata->exprgraphnode,
               exprtreessize, &nexprtrees, exprtrees, coefs) );
         assert(nexprtrees > 0);

         SCIP_CALL( consdataSetExprtrees(scip, consdata, nexprtrees, exprtrees, coefs, FALSE) );

         SCIPfreeBufferArray(scip, &exprtrees);
         SCIPfreeBufferArray(scip, &coefs);

         assert(consdata->nexprtrees > 0 );
#ifndef NDEBUG
         for( j = 0; j < consdata->nexprtrees; ++j )
            for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]); ++i )
               assert(SCIPvarIsActive(SCIPexprtreeGetVars(consdata->exprtrees[j])[i]));
#endif

         /* tell SCIP that we have something nonlinear */
         SCIPmarkNonlinearitiesPresent(scip);
         for( j = 0; !SCIPhasContinuousNonlinearitiesPresent(scip) && j < consdata->nexprtrees; ++j )
            for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]); ++i )
               if( SCIPvarGetType(SCIPexprtreeGetVars(consdata->exprtrees[j])[i]) >= SCIP_VARTYPE_CONTINUOUS )
               {
                  SCIPmarkContinuousNonlinearitiesPresent(scip);
                  break;
               }
      }
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int                c;
   int                i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check for a linear variable that can be increase or decreased without harming feasibility */
      consdataFindUnlockedLinearVar(scip, consdata);

      /* setup lincoefsmin, lincoefsmax */
      consdata->lincoefsmin = SCIPinfinity(scip);
      consdata->lincoefsmax = 0.0;
      for( i = 0; i < consdata->nlinvars; ++i )
      {
         consdata->lincoefsmin = MIN(consdata->lincoefsmin, REALABS(consdata->lincoefs[i]));  /*lint !e666*/
         consdata->lincoefsmax = MAX(consdata->lincoefsmax, REALABS(consdata->lincoefs[i]));  /*lint !e666*/
      }

      /* add nlrow respresentation to NLP, if NLP had been constructed */
      if( SCIPisNLPConstructed(scip) && SCIPconsIsChecked(conss[c]) )
      {
         if( consdata->nlrow == NULL )
         {
            SCIP_CALL( createNlRow(scip, conss[c]) );
            assert(consdata->nlrow != NULL);
         }
         SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
      }
   }

   conshdlrdata->newsoleventfilterpos = -1;
   if( nconss != 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, &conshdlrdata->newsoleventfilterpos) );
   }

   /* reset flags and counters */
   conshdlrdata->sepanlp = FALSE;
   conshdlrdata->lastenfolpnode = NULL;
   conshdlrdata->nenfolprounds = 0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->newsoleventfilterpos >= 0 )
   {
      SCIP_EVENTHDLR* eventhdlr;

      eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_newsolution");
      assert(eventhdlr != NULL);

      SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_SOLFOUND, eventhdlr, (SCIP_EVENTDATA*)conshdlr, conshdlrdata->newsoleventfilterpos) );
      conshdlrdata->newsoleventfilterpos = -1;
   }

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteNonlinear)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsActive(cons));
   assert(consdata != NULL);
   assert(SCIPconsGetData(cons) == *consdata);

   SCIPdebugMessage("consDelete for cons <%s>\n", SCIPconsGetName(cons));

   /* expression should have been removed from expression graph when constraint was deactivated */
   assert((*consdata)->exprgraphnode == NULL);

   SCIP_CALL( consdataFree(scip, consdata) );

   assert(*consdata == NULL);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransNonlinear)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int            i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->lhs, sourcedata->rhs,
         sourcedata->nlinvars, sourcedata->linvars, sourcedata->lincoefs,
         sourcedata->nexprtrees, sourcedata->exprtrees, sourcedata->nonlincoefs,
         FALSE) );

   /* copy information on curvature, if known in original constraint */
   if( sourcedata->iscurvchecked && sourcedata->nexprtrees > 0 )
   {
      BMScopyMemoryArray(targetdata->curvatures, sourcedata->curvatures, sourcedata->nexprtrees);
      targetdata->curvature = sourcedata->curvature;
      targetdata->iscurvchecked = TRUE;
   }

   for( i = 0; i < targetdata->nlinvars; ++i )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, targetdata->linvars[i], &targetdata->linvars[i]) );
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->linvars[i]) );
   }

   for( i = 0; i < targetdata->nexprtrees; ++i )
   {
      SCIP_CALL( SCIPgetExprtreeTransformedVars(scip, targetdata->exprtrees[i]) );
   }

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   SCIPdebugMessage("created transformed nonlinear constraint ");
   SCIPdebug( SCIPprintCons(scip, *targetcons, NULL) );

   return SCIP_OKAY;
}

/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_ROW*          row;
   int                c;
   SCIP_Real**        x;
   int                nvars;
   int                i;
   int                j;
   SCIP_VAR*          var;
   SCIP_Bool          haveunboundedvar;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      row = NULL;

      if( consdata->nexprtrees == 0 )
      {
         assert(consdata->exprgraphnode == NULL);
         /* if we are actually linear, add the constraint as row to the LP */
         SCIP_CALL( SCIPcreateEmptyRow(scip, &row, SCIPconsGetName(conss[c]), consdata->lhs, consdata->rhs, SCIPconsIsLocal(conss[c]), FALSE , TRUE) );  /*lint !e613*/
         SCIP_CALL( SCIPaddVarsToRow(scip, row, consdata->nlinvars, consdata->linvars, consdata->lincoefs) );
         SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE) );
         SCIP_CALL( SCIPreleaseRow (scip, &row) );
         continue;
      }

      /* setup reference points for each exprtree */
      SCIP_CALL( SCIPallocBufferArray(scip, &x, consdata->nexprtrees) );
      haveunboundedvar = FALSE;
      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         nvars = SCIPexprtreeGetNVars(consdata->exprtrees[j]);

         SCIP_CALL( SCIPallocBufferArray(scip, &x[j], nvars) );  /*lint !e866*/
         for( i = 0; i < nvars; ++i )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[j])[i];
            assert(var != NULL);
            /* use midpoint as reference value, if both bounds are finite
             * otherwise use 0.0, projected on bounds
             */
            if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
            {
               if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
               {
                  x[j][i] = 0.0;
                  haveunboundedvar = TRUE;
               }
               else
                  x[j][i] = MIN(0.0, SCIPvarGetUbGlobal(var));  /*lint !e666*/
            }
            else
            {
               if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
                  x[j][i] = MAX(0.0, SCIPvarGetLbGlobal(var));  /*lint !e666*/
               else
               {
                  x[j][i] = (SCIPvarGetLbGlobal(var) + SCIPvarGetUbGlobal(var)) / 2.0;
                  /* shift refpoint into [-INITLPMAXVARVAL, INITLPMAXVARVAL], if bounds allow */
                  if( x[j][i] < -INITLPMAXVARVAL && SCIPvarGetUbGlobal(var) >= -INITLPMAXVARVAL )
                     x[j][i] = -INITLPMAXVARVAL;
                  else if( x[j][i] > INITLPMAXVARVAL && SCIPvarGetLbGlobal(var) <= INITLPMAXVARVAL )
                     x[j][i] =  INITLPMAXVARVAL;
               }
            }
         }
      }

      /* for inequalities that are convex or that have bounded variables, try to generate a cut */
      if( !SCIPisInfinity(scip,  consdata->rhs) && ((consdata->curvature & SCIP_EXPRCURV_CONVEX)  || !haveunboundedvar) )
      {
         SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], x, NULL, TRUE, SCIP_SIDETYPE_RIGHT, &row, conshdlrdata->cutmaxrange, FALSE, FALSE) );  /*lint !e613*/

         if( row != NULL )
         {
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }
      }

      if( !SCIPisInfinity(scip, -consdata->lhs) && ((consdata->curvature & SCIP_EXPRCURV_CONCAVE) || !haveunboundedvar) )
      {
         SCIP_CALL( generateCut(scip, conshdlrdata->exprinterpreter, conss[c], x, NULL, TRUE, SCIP_SIDETYPE_LEFT, &row, conshdlrdata->cutmaxrange, FALSE, FALSE) );  /*lint !e613*/

         if( row != NULL )
         {
            SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );
            SCIPdebug( SCIP_CALL( SCIPprintRow(scip, row, NULL) ) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }
      }

      /* @todo could add more linearizations for convex or multivariate concave inequ. */

      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         SCIPfreeBufferArray(scip, &x[j]);
      }
      SCIPfreeBufferArray(scip, &x);
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          newsol;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlrdata->exprgraph, conshdlrdata->exprinterpreter, conss, nconss, NULL, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   newsol = FALSE;

   /* at root, check if we want to solve the NLP relaxation and use its solutions as reference point
    * if there is something convex, then linearizing in the solution of the NLP relaxation can be very useful
    */
   if( SCIPgetDepth(scip) == 0 && !conshdlrdata->sepanlp &&
      (SCIPgetNContVars(scip) >= conshdlrdata->sepanlpmincont * SCIPgetNVars(scip) || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY) &&
      SCIPisNLPConstructed(scip) && SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_CONSDATA* consdata;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Bool solvednlp;   /* whether we invoked an NLP solve here */
      int c;

      solstat = SCIPgetNLPSolstat(scip);
      solvednlp = FALSE;
      if( solstat == SCIP_NLPSOLSTAT_UNKNOWN )
      {
         /* NLP is not solved yet, so we might want to do this
          * but first check whether there is a violated constraint side which corresponds to a convex function
          * @todo put this check into initsol and update via consenable/consdisable
          */
         for( c = 0; c < nconss; ++c )
         {
            assert(conss[c] != NULL);  /*lint !e613*/

            consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
            assert(consdata != NULL);

            /* skip feasible constraints */
            if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
               continue;

            /* make sure curvature has been checked */
            SCIP_CALL( checkCurvature(scip, conss[c], conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );  /*lint !e613*/

            if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && (consdata->curvature & SCIP_EXPRCURV_CONVEX )) ||
               ( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && (consdata->curvature & SCIP_EXPRCURV_CONCAVE)) )
               break;
         }

         if( c < nconss )
         {
            /* try to solve NLP and update solstat */

            /* ensure linear conss are in NLP */
            if( conshdlrdata->subnlpheur != NULL )
            {
               SCIP_CALL( SCIPaddLinearConsToNlpHeurSubNlp(scip, conshdlrdata->subnlpheur, TRUE, TRUE) );
            }

            /* set LP solution as starting values, if available */
            if( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );
            }

            /* SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, 1) ); */
            SCIP_CALL( SCIPsolveNLP(scip) );

            solstat = SCIPgetNLPSolstat(scip);
            SCIPdebugMessage("solved NLP relax, solution status: %d\n", solstat);

            solvednlp = TRUE;
         }
      }

      conshdlrdata->sepanlp = TRUE;

      if( solstat == SCIP_NLPSOLSTAT_GLOBINFEASIBLE )
      {
         SCIPdebugMessage("NLP relaxation is globally infeasible, thus can cutoff node\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( solstat <= SCIP_NLPSOLSTAT_FEASIBLE )
      {
         /* if we have feasible NLP solution, generate linearization cuts there */
         SCIP_Bool lpsolseparated;
         SCIP_SOL* nlpsol;

         SCIP_CALL( SCIPcreateNLPSol(scip, &nlpsol, NULL) );
         assert(nlpsol != NULL);

         /* if we solved the NLP and solution is integral, then pass it to trysol heuristic */
         if( solvednlp && conshdlrdata->trysolheur != NULL )
         {
            int nfracvars;

            nfracvars = 0;
            if( SCIPgetNBinVars(scip) > 0 || SCIPgetNIntVars(scip) > 0 )
            {
               SCIP_CALL( SCIPgetNLPFracVars(scip, NULL, NULL, NULL, &nfracvars, NULL) );
            }

            if( nfracvars == 0 )
            {
               SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, nlpsol) );
            }
         }

         SCIP_CALL( addLinearizationCuts(scip, conshdlr, conss, nconss, nlpsol, &lpsolseparated, conshdlrdata->mincutefficacysepa) );
         newsol = TRUE;

         SCIP_CALL( SCIPfreeSol(scip, &nlpsol) );

         /* if a cut that separated the LP solution was added, then return, otherwise continue with usual separation in LP solution */
         if( lpsolseparated )
         {
            SCIPdebugMessage("linearization cuts separate LP solution\n");

            *result = SCIP_SEPARATED;

            return SCIP_OKAY;
         }
      }
   }
   /* if we do not want to try solving the NLP, or have no NLP, or have no NLP solver, or solving the NLP failed,
    * or separating with NLP solution as reference point failed, then try (again) with LP solution as reference point
    */

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, newsol, conshdlrdata->mincutefficacysepa, FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(sol != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conshdlrdata->exprgraph, conshdlrdata->exprinterpreter, conss, nconss, sol, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, FALSE, conshdlrdata->mincutefficacysepa, FALSE, result, NULL) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_CONS*         maxviolcons;
   SCIP_Real          maxviol;
   SCIP_RESULT        propresult;
   SCIP_RESULT        separateresult;
   int                dummy;
   int                nnotify;
   SCIP_Real          sepaefficacy;
   SCIP_Real          minefficacy;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlrdata->exprgraph, conshdlrdata->exprinterpreter, conss, nconss, NULL, &maxviolcons) );
   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   consdata = SCIPconsGetData(maxviolcons);
   assert(consdata != NULL);
   maxviol = consdata->lhsviol + consdata->rhsviol;
   assert(SCIPisGT(scip, maxviol, SCIPfeastol(scip)));

   SCIPdebugMessage("enfolp with max violation %g in cons <%s>\n", maxviol, SCIPconsGetName(maxviolcons));

   /* if we are above the 100'th enforcement round for this node, something is strange
    * (maybe the LP does not think that the cuts we add are violated, or we do ECP on a high-dimensional convex function)
    * in this case, check if some limit is hit or SCIP should stop for some other reason and terminate enforcement by creating a dummy node
    * (in optimized more, returning SCIP_INFEASIBLE in *result would be sufficient, but in debug mode this would give an assert in scip.c)
    * the reason to wait for 100 rounds is to avoid calls to SCIPisStopped in normal runs, which may be expensive
    * we only increment nenfolprounds until 101 to avoid an overflow
    */
   if( conshdlrdata->lastenfolpnode == SCIPgetCurrentNode(scip) )
   {
      if( conshdlrdata->nenfolprounds > 100 )
      {
         if( SCIPisStopped(scip) )
         {
            SCIP_NODE* child;

            SCIP_CALL( SCIPcreateChild(scip, &child, 1.0, SCIPnodeGetEstimate(SCIPgetCurrentNode(scip))) );
            *result = SCIP_BRANCHED;

            return SCIP_OKAY;
         }
      }
      else
         ++conshdlrdata->nenfolprounds;
   }
   else
   {
      conshdlrdata->lastenfolpnode = SCIPgetCurrentNode(scip);
      conshdlrdata->nenfolprounds = 0;
   }

   /* run domain propagation */
   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &dummy, &dummy) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we would like a cut that is efficient enough that it is not redundant in the LP (>feastol)
    * however, if the maximal violation is very small, also the best cut efficacy cannot be large
    * thus, in the latter case, we are also happy if the efficacy is at least, say, 75% of the maximal violation
    * but in any case we need an efficacy that is at least feastol
    */
   minefficacy = MIN(0.75*maxviol, conshdlrdata->mincutefficacyenfofac * SCIPfeastol(scip));  /*lint !e666*/
   minefficacy = MAX(minefficacy, SCIPfeastol(scip));  /*lint !e666*/
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, FALSE, minefficacy, TRUE, &separateresult, &sepaefficacy) );
   if( separateresult == SCIP_SEPARATED )
   {
      SCIPdebugMessage("separation succeeded (bestefficacy = %g, minefficacy = %g)\n", sepaefficacy, minefficacy);
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }

   /* we are not feasible, the whole node is not infeasible, and we cannot find a good cut
    * -> collect variables for branching
    */

   SCIPdebugMessage("separation failed (bestefficacy = %g < %g = minefficacy ); max viol: %g\n", sepaefficacy, minefficacy, maxviol);

   /* find branching candidates */
   SCIP_CALL( registerBranchingVariables(scip, conshdlr, conss, nconss, &nnotify) );

   if( nnotify == 0 && !solinfeasible && minefficacy > SCIPfeastol(scip) )
   {
      /* fallback 1: we also have no branching candidates, so try to find a weak cut */
      SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, FALSE, SCIPfeastol(scip), TRUE, &separateresult, &sepaefficacy) );
      if( separateresult == SCIP_SEPARATED )
      {
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   }

   if( nnotify == 0 && !solinfeasible )
   {
      /* fallback 2: separation probably failed because of numerical difficulties with a convex constraint;
       * if noone declared solution infeasible yet and we had not even found a weak cut, try to resolve by branching
       */
      SCIP_VAR* brvar = NULL;
      SCIP_CALL( registerLargeLPValueVariableForBranching(scip, conss, nconss, &brvar) );
      if( brvar == NULL )
      {
         /* fallback 3: all nonlinear variables in all violated constraints seem to be fixed -> replace by linear constraints */
         SCIP_Bool addedcons;

         SCIPdebugMessage("All nonlinear variables seem to be fixed. Replace remaining violated nonlinear constraints by linear constraints.\n");
         SCIP_CALL( replaceViolatedByLinearConstraints(scip, conss, nconss, &addedcons) );
         /* if the linear constraints are actually feasible, then adding them and returning SCIP_CONSADDED confuses SCIP when it enforces the new constraints again and nothing resolves the infeasiblity that we declare here
          * thus, we only add them if considered violated, and otherwise claim the solution is feasible (but print a warning) */
         if( addedcons )
         {
            *result = SCIP_CONSADDED;
         }
         else
         {
            *result = SCIP_FEASIBLE;
            SCIPwarningMessage("could not enforce feasibility by separating or branching; declaring solution with viol %g as feasible\n", maxviol);
         }
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage("Could not find any usual branching variable candidate. Proposed variable <%s> with LP value %g for branching.\n", SCIPvarGetName(brvar), SCIPgetSolVal(scip, NULL, brvar));
         nnotify = 1;
      }
   }

   assert(*result == SCIP_INFEASIBLE && (solinfeasible || nnotify > 0));
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcons;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        propresult;
   SCIP_VAR*          var;
   int                dummy;
   int                nnotify;
   int                c;
   int                i;
   int                j;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conshdlrdata->exprgraph, conshdlrdata->exprinterpreter, conss, nconss, NULL, &maxviolcons) );
   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   SCIPdebugMessage("enfops with max violation in cons <%s>\n", SCIPconsGetName(maxviolcons));

   /* run domain propagation */
   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, &dummy, &dummy) );
   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* we are not feasible and we cannot proof that the whole node is infeasible
    * -> collect all variables in violated constraints for branching
    */

   nnotify = 0;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      for( i = 0; i < consdata->nlinvars; ++i )
      {
         var = consdata->linvars[i];
         if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
            ++nnotify;
         }
      }

      for( j = 0; j < consdata->nexprtrees; ++j )
         for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]); ++i )
         {
            var = SCIPexprtreeGetVars(consdata->exprtrees[j])[i];
            if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
            {
               SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
               ++nnotify;
            }
         }
   }

   if( nnotify == 0 )
   {
      SCIPdebugMessage("All variables in violated constraints fixed (up to epsilon). Cannot find branching candidate. Forcing solution of LP.\n");
      *result = SCIP_SOLVELP;
   }

   assert(*result == SCIP_SOLVELP || (*result == SCIP_INFEASIBLE && nnotify > 0));
   return SCIP_OKAY;
}  /*lint !e715*/


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          maxviol;
   int                c;
   SCIP_Bool          maypropfeasible; /* whether we may be able to propose a feasible solution */

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(sol != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   /* during presolve, we do not have exprtrees in the constraints, but we can get values from the expression graph, if we have evaluated it */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
   {
      SCIP_Real* varvals;

      assert(conshdlrdata->exprgraph != NULL);

      SCIP_CALL( SCIPallocBufferArray(scip, &varvals, SCIPexprgraphGetNVars(conshdlrdata->exprgraph)) );
      SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPexprgraphGetNVars(conshdlrdata->exprgraph), (SCIP_VAR**)SCIPexprgraphGetVars(conshdlrdata->exprgraph), varvals) );

      SCIP_CALL( SCIPexprgraphEval(conshdlrdata->exprgraph, varvals) );

      SCIPfreeBufferArray(scip, &varvals);
   }

   /* @todo adapt proposeFeasibleSolution to function also during presolving */
   maxviol = 0.0;
   maypropfeasible = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL) &&
      SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED && SCIPgetStage(scip) <= SCIP_STAGE_SOLVING && SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( computeViolation(scip, conshdlrdata->exprinterpreter, conss[c], sol) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g (scaled: %.15g)\n", consdata->lhs - consdata->activity, consdata->lhsviol);
            }
            if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g (scaled: %.15g)\n", consdata->activity - consdata->rhs, consdata->rhsviol);
            }
         }
         if( (conshdlrdata->subnlpheur == NULL || sol == NULL) && !maypropfeasible )
            return SCIP_OKAY;
         if( consdata->lhsviol > maxviol || consdata->rhsviol > maxviol )
            maxviol = MAX(consdata->lhsviol, consdata->rhsviol);
         if( maypropfeasible )
         {
            /* update information on linear variables that may be in- or decreased */
            if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
               consdataFindUnlockedLinearVar(scip, consdata);

            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               /* check if there is a variable which may help to get the left hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] > 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] < 0.0) )
                  maypropfeasible = FALSE;
            }
            else
            {
               assert(SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)));
               /* check if there is a variable which may help to get the right hand side satisfied
                * if there is no such var, then we cannot get feasible */
               if( !(consdata->linvar_mayincrease >= 0 && consdata->lincoefs[consdata->linvar_mayincrease] < 0.0) &&
                  ! (consdata->linvar_maydecrease >= 0 && consdata->lincoefs[consdata->linvar_maydecrease] > 0.0) )
                  maypropfeasible = FALSE;
            }
         }
      }
      else
      {
         /* SCIPdebugMessage("constraint <%s> is feasible (%g, %g) in check, activity = %g, sides = [%g, %g]\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol, consdata->activity, consdata->lhs, consdata->rhs); */
      }
   }

   if( *result == SCIP_INFEASIBLE && maypropfeasible )
   {
      SCIP_Bool success;

      SCIP_CALL( proposeFeasibleSolution(scip, conshdlr, conss, nconss, sol, &success) );

      /* do not pass solution to NLP heuristic if we made it feasible this way */
      if( success )
         return SCIP_OKAY;
   }

   if( *result == SCIP_INFEASIBLE && conshdlrdata->subnlpheur != NULL && sol != NULL )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}  /*lint !e715*/


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropNonlinear)
{
   int dummy;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   dummy = 0;
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, result, &dummy, &dummy) );

   return SCIP_OKAY;
}  /*lint !e715*/

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolNonlinear)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_RESULT        propresult;
   SCIP_Bool          havechange;
   SCIP_Bool          domainerror;
   SCIP_Bool          havegraphchange;
   SCIP_Bool          doreformulations;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   /* if other presolvers did not find any changes (except for deleted conss) since last call,
    * then try the reformulations (replacing products with binaries, disaggregation, setting default variable bounds)
    * otherwise, we wait with these
    */
   doreformulations = (nrounds > 0 || SCIPconshdlrWasPresolvingDelayed(conshdlr)) &&
      nnewfixedvars == 0 && nnewaggrvars == 0 && nnewchgvartypes == 0 && nnewchgbds == 0 &&
      nnewholes == 0 && /* nnewdelconss == 0 && */ nnewaddconss == 0 && nnewupgdconss == 0 && nnewchgcoefs == 0 && nnewchgsides == 0;
   SCIPdebugMessage("presolving will %swait with reformulation\n", doreformulations ? "not " : "");

   havegraphchange = FALSE;

   if( !conshdlrdata->isremovedfixings )
   {
      SCIP_CALL( removeFixedNonlinearVariables(scip, conshdlr) );
      assert(conshdlrdata->isremovedfixings);

      havegraphchange = TRUE;
   }

   SCIP_CALL( SCIPexprgraphSimplify(conshdlrdata->exprgraph, SCIPepsilon(scip), conshdlrdata->maxexpansionexponent, &havechange, &domainerror) );
   SCIPdebugMessage("expression graph simplifier found %schange, domain error = %u\n", havechange ? "" : "no ", domainerror);
   havegraphchange |= havechange;

   /* if simplifier found some undefined expression, then declare problem as infeasible
    * usually, this should be discovered during domain propagation already, but since that is using interval arithmetics,
    *   it may overestimate in a way that actually undefined expressions still get a value assigned (e.g., 0^(-1) = [-inf,inf]) */
   if( domainerror )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   if( nrounds == 0 )
   {
      /* upgrade methods may look at expression graph bounds, which are not present in the first presolving round yet */
      SCIP_CALL( SCIPexprgraphPropagateVarBounds(conshdlrdata->exprgraph, INTERVALINFTY, TRUE, &domainerror) );

      if( domainerror )
      {
         SCIPdebugMessage("propagating variable bounds through expression graph found that some expressions cannot be evaluated w.r.t. current bounds, thus cutoff\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      SCIPdebugMessage("process constraint <%s>\n", SCIPconsGetName(conss[c]));
      SCIPdebug( SCIPprintCons(scip, conss[c], NULL) );

      havechange = FALSE;

      if( !consdata->isremovedfixingslin )
      {
         SCIP_CALL( removeFixedLinearVariables(scip, conss[c]) );
         assert(consdata->isremovedfixingslin);
         havechange = TRUE;
      }

      if( !consdata->ispresolved || havegraphchange )
      {
         SCIP_CALL( splitOffLinearPart(scip, conshdlr, conss[c]) );
      }

      if( consdata->nlinvars == 0 && consdata->exprgraphnode == NULL )
      {
         /* all variables fixed or removed, constraint function is 0.0 now */
         if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisFeasPositive(scip, consdata->lhs)) ||
            ( !SCIPisInfinity(scip,  consdata->rhs) && SCIPisFeasNegative(scip, consdata->rhs)) )
         {
            /* left hand side positive or right hand side negative */
            SCIPdebugMessage("constraint <%s> is constant and infeasible\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else
         {
            /* left and right hand side are consistent */
            SCIPdebugMessage("constraint <%s> is constant and feasible, deleting\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelCons(scip, conss[c]) );
            ++*ndelconss;
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      /* call upgrade methods if constraint was not presolved, has been changed, or the expression graph has changed */
      if( !consdata->ispresolved || havechange || havegraphchange )
      {
         SCIP_Bool upgraded;

         SCIP_CALL( presolveUpgrade(scip, conshdlr, conss[c], &upgraded, nupgdconss, naddconss) );
         if( upgraded )
         {
            *result = SCIP_SUCCESS;
            continue;
         }
      }

      consdata->ispresolved = TRUE;
   }

   /* run domain propagation */
   SCIP_CALL( propagateBounds(scip, conshdlr, conss, nconss, &propresult, nchgbds, ndelconss) );
   switch( propresult )
   {
   case SCIP_REDUCEDDOM:
      *result = SCIP_SUCCESS;
      break;
   case SCIP_CUTOFF:
      SCIPdebugMessage("propagation says problem is infeasible in presolve\n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   default:
      assert(propresult == SCIP_DIDNOTFIND || propresult == SCIP_DIDNOTRUN);
   }  /*lint !e788*/

   if( doreformulations && conshdlrdata->reformulate && !conshdlrdata->assumeconvex )
   {
      int naddconssbefore;

      naddconssbefore = conshdlrdata->naddedreformconss;
      SCIP_CALL( reformulate(scip, conshdlr, conss, nconss, &conshdlrdata->naddedreformconss) );

      if( conshdlrdata->naddedreformconss > naddconssbefore )
      {
         *result = SCIP_SUCCESS;
         *naddconss += conshdlrdata->naddedreformconss - naddconssbefore;

         /* if expression graph changed, ensure that we apply all presolving techniques (esp. upgrades) in next round again */
         for( c = 0; c < nconss; ++c )
         {
            assert(conss[c] != NULL);  /*lint !e794*/

            consdata = SCIPconsGetData(conss[c]);  /*lint !e794*/
            assert(consdata != NULL);

            consdata->ispresolved = FALSE;
         }
      }
   }

   /* if we did not try reformulations, ensure that presolving is called again even if there were only a few changes (< abortfac) */
   if( !doreformulations )
      *result = SCIP_DELAYED;

   return SCIP_OKAY;
}  /*lint !e715*/

/** propagation conflict resolving method of constraint handler */
#define consRespropNonlinear NULL

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockNonlinear)
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool      havelhs;
   SCIP_Bool      haverhs;
   int            i;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   havelhs = !SCIPisInfinity(scip, -consdata->lhs);
   haverhs = !SCIPisInfinity(scip,  consdata->rhs);

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      if( consdata->lincoefs[i] > 0 )
      {
         if( havelhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
         if( haverhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( havelhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlocksneg, nlockspos) );
         }
         if( haverhs )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->linvars[i], nlockspos, nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}  /*lint !e715*/

/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("activate cons <%s>\n", SCIPconsGetName(cons));

   if( consdata->nexprtrees > 0 )
   {
      SCIP_Bool exprtreeisnew;

      assert(consdata->exprgraphnode == NULL);

      /* add exprtrees to expression graph */
      SCIP_CALL( SCIPexprgraphAddExprtreeSum(conshdlrdata->exprgraph, consdata->nexprtrees, consdata->exprtrees, consdata->nonlincoefs, &consdata->exprgraphnode, &exprtreeisnew) );
      assert(consdata->exprgraphnode != NULL);
      /* @todo do something with exprtreeisnew? */

      /* if during presolving, then forget expression trees */
      if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      {
         SCIP_CALL( consdataSetExprtrees(scip, consdata, 0, NULL, NULL, FALSE) );
      }

      /* remember that we should run reformulation again */
      conshdlrdata->isreformulated = FALSE;
   }
   else if( consdata->exprgraphnode != NULL )
   {
      /* if constraint already comes with node in expression graph, then also remember that we should run reformulation again */
      conshdlrdata->isreformulated = FALSE;
   }

   return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->exprgraphnode != NULL || consdata->nexprtrees == 0);

   SCIPdebugMessage("deactivate cons <%s>\n", SCIPconsGetName(cons));

   if( consdata->exprgraphnode != NULL )
   {
      if( consdata->nexprtrees == 0 )
      {
         /* during presolving, the exprtrees in the constraint are removed, so put them back before releasing the exprgraphnode */
         SCIP_EXPRTREE* exprtree;

         /* if only presolve is run and problem is found infeasible there, then constraints may not be deactivated there, but in a later call to freeTransform */
         assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING || SCIPgetStage(scip) == SCIP_STAGE_FREETRANS);

         SCIP_CALL( SCIPexprgraphGetTree(conshdlrdata->exprgraph, consdata->exprgraphnode, &exprtree) );
         SCIP_CALL( consdataSetExprtrees(scip, consdata, 1, &exprtree, NULL, FALSE) );
      }

      SCIP_CALL( SCIPexprgraphReleaseNode(conshdlrdata->exprgraph, &consdata->exprgraphnode) );
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));
   assert(SCIPconsIsActive(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("enable cons <%s>\n", SCIPconsGetName(cons));

   if( consdata->exprgraphnode != NULL )
   {
      /* enable node of expression in expression graph */
      SCIPexprgraphEnableNode(conshdlrdata->exprgraph, consdata->exprgraphnode);
   }

   /* enable event catching for linear variables */
   consdata->isremovedfixingslin = TRUE;
   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( catchLinearVarEvents(scip, cons, i) );

      consdata->isremovedfixingslin = consdata->isremovedfixingslin && SCIPvarIsActive(consdata->linvars[i]);
   }

   return SCIP_OKAY;
}

/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableNonlinear)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(SCIPconsIsTransformed(cons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lineventdata != NULL || consdata->nlinvars == 0);

   SCIPdebugMessage("disable cons <%s>\n", SCIPconsGetName(cons));

   /* disable node of expression in expression graph */
   if( consdata->exprgraphnode != NULL )
   {
      SCIPexprgraphDisableNode(conshdlrdata->exprgraph, consdata->exprgraphnode);
   }

   for( i = 0; i < consdata->nlinvars; ++i )
   {
      SCIP_CALL( dropLinearVarEvents(scip, cons, i) );
   }

   return SCIP_OKAY;
}

/** variable deletion method of constraint handler */
#define consDelvarsNonlinear NULL

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintNonlinear)
{
   SCIP_CONSDATA* consdata;
   int            j;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   if( consdata->nlinvars == 0 && consdata->nexprtrees == 0 && consdata->exprgraphnode == 0 )
   {
      SCIPinfoMessage(scip, file, "0 ");
   }
   else
   {
      if( consdata->nexprtrees > 0 )
      {
         for( j = 0; j < consdata->nexprtrees; ++j )
         {
            if( j > 0 || consdata->nonlincoefs[j] != 1.0 )
               SCIPinfoMessage(scip, file, " %+.20g ", consdata->nonlincoefs[j]);
            SCIP_CALL( SCIPexprtreePrintWithNames(consdata->exprtrees[j], file) );
         }
      }
      else if( consdata->exprgraphnode != NULL )
      {
         SCIP_CONSHDLRDATA* conshdlrdata;
         SCIP_EXPRTREE* tree;

         conshdlrdata = SCIPconshdlrGetData(conshdlr);
         assert(conshdlrdata != NULL);
         SCIP_CALL( SCIPexprgraphGetTree(conshdlrdata->exprgraph, consdata->exprgraphnode, &tree) );

         SCIP_CALL( SCIPexprtreePrintWithNames(tree, file) );

         SCIP_CALL( SCIPexprtreeFree(&tree) );
      }

      for( j = 0; j < consdata->nlinvars; ++j )
      {
         SCIPinfoMessage(scip, file, "%+.15g<%s>[%c] ", consdata->lincoefs[j], SCIPvarGetName(consdata->linvars[j]),
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_BINARY ? 'B' :
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_INTEGER ? 'I' :
            SCIPvarGetType(consdata->linvars[j]) == SCIP_VARTYPE_IMPLINT ? 'I' : 'C');
      }
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   }
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   }
   else
   {
      SCIPinfoMessage(scip, file, " [free]");
   }

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyNonlinear)
{
   SCIP_CONSDATA*    consdata;
   SCIP_CONSDATA*    targetconsdata;
   SCIP_VAR**        linvars;
   SCIP_EXPRTREE**   exprtrees;
   int i;
   int j;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourceconshdlr != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);

   linvars = NULL;
   exprtrees = NULL;

   *valid = TRUE;

   if( consdata->nlinvars != 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &linvars, consdata->nlinvars) );
      for( i = 0; i < consdata->nlinvars && *valid; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->linvars[i], &linvars[i], varmap, consmap, global, valid) );
         assert(!*valid || linvars[i] != NULL);
      }
   }

   if( *valid && consdata->nexprtrees > 0 )
   {
      SCIP_VAR** nonlinvars;

      SCIP_CALL( SCIPallocBufferArray(sourcescip, &exprtrees, consdata->nexprtrees) );
      BMSclearMemoryArray(exprtrees, consdata->nexprtrees);
      SCIP_CALL( SCIPallocBufferArray(sourcescip, &nonlinvars, SCIPexprtreeGetNVars(consdata->exprtrees[0])) );

      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         SCIP_CALL( SCIPreallocBufferArray(sourcescip, &nonlinvars, SCIPexprtreeGetNVars(consdata->exprtrees[j])) );
         for( i = 0; i < SCIPexprtreeGetNVars(consdata->exprtrees[j]) && *valid; ++i )
         {
            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, SCIPexprtreeGetVars(consdata->exprtrees[j])[i], &nonlinvars[i], varmap, consmap, global, valid) );
            assert(!*valid || nonlinvars[i] != NULL);
         }

         if( *valid )
         {
            SCIP_CALL( SCIPexprtreeCopy(SCIPblkmem(scip), &exprtrees[j], consdata->exprtrees[j]) );
            SCIP_CALL( SCIPexprtreeSetVars(exprtrees[j], SCIPexprtreeGetNVars(consdata->exprtrees[j]), nonlinvars) );
         }
         else
            break;
      }

      SCIPfreeBufferArray(sourcescip, &nonlinvars);
   }

   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsNonlinear(scip, cons, name ? name : SCIPconsGetName(sourcecons),
            consdata->nlinvars, linvars, consdata->lincoefs,
            consdata->nexprtrees, exprtrees, consdata->nonlincoefs,
            consdata->lhs, consdata->rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

      /* copy information on curvature */
      targetconsdata = SCIPconsGetData(*cons);
      targetconsdata->curvature     = consdata->curvature;
      targetconsdata->iscurvchecked = consdata->iscurvchecked && global; /* if the copy is local, then curvature may change (get stronger) */
   }

   SCIPfreeBufferArrayNull(sourcescip, &linvars);
   if( exprtrees != NULL )
   {
      for( j = 0; j < consdata->nexprtrees; ++j )
      {
         if( exprtrees[j] != NULL )
         {
            SCIP_CALL( SCIPexprtreeFree(&exprtrees[j]) );
         }
      }
      SCIPfreeBufferArray(sourcescip, &exprtrees);
   }

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseNonlinear)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of nonlinear constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/** creates the handler for nonlinear constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrNonlinear(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create nonlinear constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyNonlinear,
         consFreeNonlinear, consInitNonlinear, consExitNonlinear,
         consInitpreNonlinear, consExitpreNonlinear, consInitsolNonlinear, consExitsolNonlinear,
         consDeleteNonlinear, consTransNonlinear, consInitlpNonlinear,
         consSepalpNonlinear, consSepasolNonlinear, consEnfolpNonlinear, consEnfopsNonlinear, consCheckNonlinear,
         consPropNonlinear, consPresolNonlinear, consRespropNonlinear, consLockNonlinear,
         consActiveNonlinear, consDeactiveNonlinear,
         consEnableNonlinear, consDisableNonlinear, consDelvarsNonlinear,
         consPrintNonlinear, consCopyNonlinear, consParseNonlinear,
         conshdlrdata) );

   /* add nonlinear constraint handler parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacysepa",
         "minimal efficacy for a cut to be added to the LP during separation; overwrites separating/efficacy",
         &conshdlrdata->mincutefficacysepa, TRUE, 0.0001, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacyenfofac",
         "minimal target efficacy of a cut in order to add it to relaxation during enforcement as a factor of the feasibility tolerance (may be ignored)",
         &conshdlrdata->mincutefficacyenfofac, TRUE, 2.0, 1.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/cutmaxrange",
         "maximal coef range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",
         &conshdlrdata->cutmaxrange, FALSE, 1e+7, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/linfeasshift",
         "whether to try to make solutions in check function feasible by shifting a linear variable (esp. useful if constraint was actually objective function)",
         &conshdlrdata->linfeasshift, FALSE, TRUE, NULL, NULL) );

#if 0 /* don't have any expensive checks yet, so we disable this parameter for now */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/checkconvexexpensive",
         "whether to apply expensive curvature checking methods",
         &conshdlrdata->checkconvexexpensive, FALSE, TRUE, NULL, NULL) );
#endif

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/assumeconvex",
         "whether to assume that nonlinear functions in inequalities (<=) are convex (disables reformulation)",
         &conshdlrdata->assumeconvex, TRUE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxproprounds",
         "limit on number of propagation rounds for a single constraint within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 1, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/reformulate",
         "whether to reformulate expression graph",
         &conshdlrdata->reformulate, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/maxexpansionexponent",
         "maximal exponent where still expanding non-monomial polynomials in expression simplification",
         &conshdlrdata->maxexpansionexponent, TRUE, 2, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/sepanlpmincont",
         "minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation",
         &conshdlrdata->sepanlpmincont, FALSE, 1.0, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_boundchange", "signals a bound change to a nonlinear constraint",
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, processLinearVarEvent, NULL) );
   conshdlrdata->linvareventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_boundchange");

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_boundchange2", "signals a bound change to a nonlinear constraint handler",
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, processNonlinearVarEvent, (SCIP_EVENTHDLRDATA*)conshdlrdata) );
   conshdlrdata->nonlinvareventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME"_boundchange2");

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_newsolution", "handles the event that a new primal solution has been found",
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, processNewSolutionEvent, NULL) );

   /* create expression interpreter */
   SCIP_CALL( SCIPexprintCreate(SCIPblkmem(scip), &conshdlrdata->exprinterpreter) );

   /* create expression graph */
   SCIP_CALL( SCIPexprgraphCreate(SCIPblkmem(scip), &conshdlrdata->exprgraph, -1, -1,
         exprgraphVarAdded, exprgraphVarRemove, NULL, (void*)conshdlrdata) );
   conshdlrdata->isremovedfixings = TRUE;
   conshdlrdata->ispropagated = TRUE;

   conshdlrdata->scip = scip;

   return SCIP_OKAY;
}

/** includes a nonlinear constraint upgrade method into the nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNonlinconsUpgrade(
   SCIP*                   scip,               /**< SCIP data structure */
   SCIP_DECL_NONLINCONSUPGD((*nonlinconsupgd)),/**< method to call for upgrading nonlinear constraint, or NULL */
   SCIP_DECL_EXPRGRAPHNODEREFORM((*nodereform)),/**< method to call for reformulating expression graph node, or NULL */
   int                     priority,           /**< priority of upgrading method */
   SCIP_Bool               active,             /**< should the upgrading method by active by default? */
   const char*             conshdlrname        /**< name of the constraint handler */
   )
{
   SCIP_CONSHDLR*        conshdlr;
   SCIP_CONSHDLRDATA*    conshdlrdata;
   SCIP_NLCONSUPGRADE*   nlconsupgrade;
   char                  paramname[SCIP_MAXSTRLEN];
   char                  paramdesc[SCIP_MAXSTRLEN];
   int                   i;

   assert(conshdlrname != NULL );

   /* ignore empty upgrade functions */
   if( nonlinconsupgd == NULL && nodereform == NULL )
      return SCIP_OKAY;

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether upgrade method exists already */
   for( i = conshdlrdata->nnlconsupgrades - 1; i >= 0; --i )
   {
      if( conshdlrdata->nlconsupgrades[i]->nlconsupgd == nonlinconsupgd && conshdlrdata->nlconsupgrades[i]->nodereform == nodereform)
      {
#ifdef SCIP_DEBUG
         SCIPwarningMessage("Try to add already known upgrade method pair (%p,%p) for constraint handler <%s>.\n", (void*)nonlinconsupgd, (void*)nodereform, conshdlrname);
#endif
         return SCIP_OKAY;
      }
   }

   /* create a nonlinear constraint upgrade data object */
   SCIP_CALL( SCIPallocMemory(scip, &nlconsupgrade) );
   nlconsupgrade->nlconsupgd = nonlinconsupgd;
   nlconsupgrade->nodereform = nodereform;
   nlconsupgrade->priority   = priority;
   nlconsupgrade->active     = active;

   /* insert nonlinear constraint upgrade method into constraint handler data */
   assert(conshdlrdata->nnlconsupgrades <= conshdlrdata->nlconsupgradessize);
   if( conshdlrdata->nnlconsupgrades+1 > conshdlrdata->nlconsupgradessize )
   {
      int newsize;

      newsize = SCIPcalcMemGrowSize(scip, conshdlrdata->nnlconsupgrades+1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->nlconsupgrades, newsize) );
      conshdlrdata->nlconsupgradessize = newsize;
   }
   assert(conshdlrdata->nnlconsupgrades+1 <= conshdlrdata->nlconsupgradessize);

   for( i = conshdlrdata->nnlconsupgrades; i > 0 && conshdlrdata->nlconsupgrades[i-1]->priority < nlconsupgrade->priority; --i )
      conshdlrdata->nlconsupgrades[i] = conshdlrdata->nlconsupgrades[i-1];
   assert(0 <= i && i <= conshdlrdata->nnlconsupgrades);
   conshdlrdata->nlconsupgrades[i] = nlconsupgrade;
   conshdlrdata->nnlconsupgrades++;

   /* adds parameter to turn on and off the upgrading step */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "constraints/"CONSHDLR_NAME"/upgrade/%s", conshdlrname);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "enable nonlinear upgrading for constraint handler <%s>", conshdlrname);
   SCIP_CALL( SCIPaddBoolParam(scip,
         paramname, paramdesc,
         &nlconsupgrade->active, FALSE, active, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 *  this variant takes expression trees as input
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   int                   nexprtrees,         /**< number of expression trees for nonlinear part of constraint */
   SCIP_EXPRTREE**       exprtrees,          /**< expression trees for nonlinear part of constraint */
   SCIP_Real*            nonlincoefs,        /**< coefficients for expression trees for nonlinear part, or NULL if all 1.0 */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;

   assert(linvars  != NULL || nlinvars == 0);
   assert(lincoefs != NULL || nlinvars == 0);
   assert(exprtrees   != NULL || nexprtrees == 0);
   assert(modifiable == FALSE); /* we do not support column generation */

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreateEmpty(scip, &consdata) );

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   /* add linear variables */
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, nlinvars) );
   for( i = 0; i < nlinvars; ++i )
   {
      if( SCIPisZero(scip, lincoefs[i]) )  /*lint !e613*/
         continue;

      SCIP_CALL( addLinearCoef(scip, *cons, linvars[i], lincoefs[i]) );  /*lint !e613*/
   }

   /* set expression trees */
   SCIP_CALL( consdataSetExprtrees(scip, consdata, nexprtrees, exprtrees, nonlincoefs, TRUE) );

   SCIPdebugMessage("created nonlinear constraint ");
   SCIPdebug( SCIPprintCons(scip, *cons, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a nonlinear constraint
 * this variant takes a node of the expression graph as input and can only be used during presolving
 * it is assumed that the nonlinear constraint will be added to the transformed problem short after creation
 * the given exprgraphnode is captured in this method
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsNonlinear2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nlinvars,           /**< number of linear variables in the constraint */
   SCIP_VAR**            linvars,            /**< array with linear variables of constraint entries */
   SCIP_Real*            lincoefs,           /**< array with coefficients of constraint linear entries */
   SCIP_EXPRGRAPHNODE*   exprgraphnode,      /**< expression graph node associated to nonlinear expression */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;

   assert(modifiable == FALSE); /* we do not support column generation */
   assert(SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);

   /* find the nonlinear constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("nonlinear constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreateEmpty(scip, &consdata) );

   consdata->lhs = lhs;
   consdata->rhs = rhs;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   /* add linear variables */
   SCIP_CALL( consdataEnsureLinearVarsSize(scip, consdata, nlinvars) );
   for( i = 0; i < nlinvars; ++i )
   {
      if( SCIPisZero(scip, lincoefs[i]) )
         continue;

      SCIP_CALL( addLinearCoef(scip, *cons, linvars[i], lincoefs[i]) );
   }

   /* set expression graph node */
   if( exprgraphnode != NULL )
   {
      consdata->exprgraphnode = exprgraphnode;
      consdata->curvature = SCIP_EXPRCURV_UNKNOWN;
      consdata->iscurvchecked = FALSE;
      consdata->activity = SCIP_INVALID;
      SCIPexprgraphCaptureNode(exprgraphnode);
   }

   SCIPdebugMessage("created nonlinear constraint ");
   SCIPdebug( SCIPprintCons(scip, *cons, NULL) );

   return SCIP_OKAY;
}

/** adds a linear variable with coefficient to a nonlinear constraint */
SCIP_RETCODE SCIPaddLinearVarNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< coefficient of variable */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(!SCIPisInfinity(scip, REALABS(coef)));

   SCIP_CALL( addLinearCoef(scip, cons, var, coef) );

   return SCIP_OKAY;
}

/** sets the expression trees in a nonlinear constraint
 * constraint must not be active yet
 */
SCIP_RETCODE SCIPsetExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nexprtrees,         /**< number of expression trees */
   SCIP_EXPRTREE**       exprtrees,          /**< new expression trees, or NULL if nexprtrees is 0 */
   SCIP_Real*            coefs               /**< coefficients of expression trees, or NULL if all 1.0 */
   )
{
   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPconsIsActive(cons));
   assert(SCIPconsGetData(cons) != NULL);
   assert(exprtrees != NULL || nexprtrees == 0);

   SCIP_CALL( consdataSetExprtrees(scip, SCIPconsGetData(cons), nexprtrees, exprtrees, coefs, TRUE) );

   return SCIP_OKAY;
}

/** gets the nonlinear constraint as a nonlinear row representation */
SCIP_RETCODE SCIPgetNlRowNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< a buffer where to store pointer to nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(nlrow != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow == NULL )
   {
      SCIP_CALL( createNlRow(scip, cons) );
   }
   assert(consdata->nlrow != NULL);
   *nlrow = consdata->nlrow;

   return SCIP_OKAY;
}

/** gets the number of variables in the linear term of a nonlinear constraint */
int SCIPgetNLinearVarsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->nlinvars;
}

/** gets the variables in the linear part of a nonlinear constraint */
SCIP_VAR** SCIPgetLinearVarsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->linvars;
}

/** gets the coefficients in the linear part of a nonlinear constraint */
SCIP_Real* SCIPgetLinearCoefsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lincoefs;
}

/** gets the number of expression trees of a nonlinear constraint */
int SCIPgetNExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING);

   return SCIPconsGetData(cons)->nexprtrees;
}

/** gets the expression trees of a nonlinear constraint */
SCIP_EXPRTREE** SCIPgetExprtreesNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING);

   return SCIPconsGetData(cons)->exprtrees;
}

/** gets the coefficients of the expression trees of a nonlinear constraint */
SCIP_Real* SCIPgetExprtreeCoefsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);
   assert(SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING);

   return SCIPconsGetData(cons)->nonlincoefs;
}

/** gets the expression graph node of a nonlinear constraint */
SCIP_EXPRGRAPHNODE* SCIPgetExprgraphNodeNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->exprgraphnode;
}

/** gets the left hand side of a nonlinear constraint */
SCIP_Real SCIPgetLhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->lhs;
}

/** gets the right hand side of a nonlinear constraint */
SCIP_Real SCIPgetRhsNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   assert(cons != NULL);
   assert(SCIPconsGetData(cons) != NULL);

   return SCIPconsGetData(cons)->rhs;
}

/** check the function of a nonlinear constraint for convexity/concavity, if not done yet */
SCIP_RETCODE SCIPcheckCurvatureNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(cons != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( checkCurvature(scip, cons, conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );

   return SCIP_OKAY;
}

/** gets the curvature of the nonlinear function of a nonlinear constraint */
SCIP_RETCODE SCIPgetCurvatureNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             checkcurv,          /**< whether to check constraint curvature, if not checked before */
   SCIP_EXPRCURV*        curvature           /**< buffer to store curvature of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(curvature != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( checkcurv && !consdata->iscurvchecked )
   {
      SCIP_CALL( checkCurvature(scip, cons, conshdlrdata->checkconvexexpensive, conshdlrdata->assumeconvex) );
   }

   *curvature = consdata->curvature;

   return SCIP_OKAY;
}

/** computes the violation of a nonlinear constraint by a solution */
SCIP_RETCODE SCIPgetViolationNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution which violation to calculate, or NULL for LP solution */
   SCIP_Real*            violation           /**< buffer to store violation of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(violation != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && SCIPconsIsActive(cons) )
   {
      /* @todo make available */
      SCIPwarningMessage("SCIPgetViolationNonlinear is not available for active constraints during presolve.\n");
      *violation = SCIP_INVALID;
      return SCIP_OKAY;
   }

   conshdlr = SCIPconsGetHdlr(cons);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolation(scip, conshdlrdata->exprinterpreter, cons, sol) );

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *violation = MAX(consdata->lhsviol, consdata->rhsviol);

   return SCIP_OKAY;
}

/** gets expression graph of nonlinear constraint handler */
SCIP_EXPRGRAPH* SCIPgetExprgraphNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< nonlinear constraint handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->exprgraph != NULL);

   return conshdlrdata->exprgraph;
}
