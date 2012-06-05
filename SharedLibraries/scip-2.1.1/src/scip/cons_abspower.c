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

/**@file   cons_abspower.c
 * @brief  Constraint handler for absolute power constraints \f$\textrm{lhs} \leq \textrm{sign}(x+a) |x+a|^n + c z \leq \textrm{rhs}\f$
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_abspower.h"
#include "scip/cons_nonlinear.h"
#include "scip/cons_indicator.h"
#include "scip/cons_quadratic.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "scip/intervalarith.h"
#include "scip/heur_subnlp.h"
#include "scip/heur_trysol.h"
#include "scip/debug.h"

/* constraint handler properties */
#define CONSHDLR_NAME          "abspower"
#define CONSHDLR_DESC          "constraint handler for absolute power constraints lhs <= sign(x+offset)abs(x+offset)^n + c*z <= rhs"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY       -30 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -3500000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_PROP_TIMING   SCIP_PROPTIMING_ALWAYS /**< when should the constraint handlers propagation routines be called? */

#define QUADCONSUPGD_PRIORITY     50000 /**< priority of the constraint handler for upgrading of quadratic constraints */
#define NONLINCONSUPGD_PRIORITY   50000 /**< priority of the constraint handler for upgrading of nonlinear constraints and reformulating expression graph nodes */

/*
 * Local defines
 */

#define PROPVARTOL    SCIPepsilon(scip) /**< tolerance to add to variable bounds in domain propagation */
#define PROPSIDETOL   SCIPepsilon(scip) /**< tolerance to add to constraint sides in domain propagation */
#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */
#define INITLPMAXVARVAL          1000.0 /**< maximal absolute value of variable for still generating a linearization cut at that point in initlp */

/**< power function type to be used by a constraint instead of the general pow */
#define DECL_MYPOW(x) SCIP_Real x (SCIP_Real base, SCIP_Real exponent)

/** sign of a value (-1 or +1)
 *
 * 0.0 has sign +1
 */
#define SIGN(x) ((x) >= 0.0 ? 1.0 : -1.0)

#define ROOTS_KNOWN 10                  /**< up to which (integer) exponents precomputed roots have been stored */

/*
 * Data structures
 */

/** The positive root of the polynom (n-1) y^n + n y^(n-1) - 1 is needed in separation.
 * Here we store these roots for small integer values of n.
 */
static
SCIP_Real roots[ROOTS_KNOWN+1] = {
   -1.0,                     /**< no root for n=0 */
   -1.0,                     /**< no root for n=1 */
   0.41421356237309504880,   /**< root for n=2 (-1+sqrt(2)) */
   0.5,                      /**< root for n=3 */
   0.56042566045031785945,   /**< root for n=4 */
   0.60582958618826802099,   /**< root for n=5 */
   0.64146546982884663257,   /**< root for n=6 */
   0.67033204760309682774,   /**< root for n=7 */
   0.69428385661425826738,   /**< root for n=8 */
   0.71453772716733489700,   /**< root for n=9 */
   0.73192937842370733350    /**< root for n=10 */
};

/** constraint data for absolute power constraints */
struct SCIP_ConsData
{
   SCIP_VAR*             x;                  /**< variable x in sign(x+offset)|x+offset|^n term */
   SCIP_VAR*             z;                  /**< linear variable in constraint */
   SCIP_Real             exponent;           /**< exponent n of |x+offset| */
   SCIP_Real             xoffset;            /**< offset in x+offset */
   SCIP_Real             zcoef;              /**< coefficient of linear variable z */
   SCIP_Real             lhs;                /**< left  hand side of constraint */
   SCIP_Real             rhs;                /**< right hand side of constraint */

   SCIP_Real             root;               /**< root of polynom */
   DECL_MYPOW            ((*pow));           /**< function for computing power*/

   SCIP_Real             lhsviol;            /**< current (scaled) violation of left  hand side */
   SCIP_Real             rhsviol;            /**< current (scaled) violation of right hand side */

   SCIP_Bool             isxpropagated;      /**< have all bound tightenings on x been propagated? */
   SCIP_Bool             iszpropagated;      /**< have all bound tightenings on z been propagated? */
   int                   xeventfilterpos;    /**< position of x var event in SCIP event filter */
   int                   zeventfilterpos;    /**< position of z var event in SCIP event filter */
   unsigned int          propvarbounds:1;    /**< have variable bounds been propagated? */

   SCIP_NLROW*           nlrow;              /**< nonlinear row representation of constraint */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Real             mincutefficacysepa; /**< minimal efficacy of a cut in order to add it to relaxation during separation */
   SCIP_Real             mincutefficacyenfofac;/**< minimal target efficacy of a cut in order to add it to relaxation during enforcement as factor of feasibility tolerance (may be ignored) */
   SCIP_Real             cutmaxrange;        /**< maximal coef range (maximal abs coef / minimal abs coef) of a cut in order to be added to LP */
   SCIP_Bool             projectrefpoint;    /**< whether to project the reference point when linearizing a absolute power constraint in a convex region */
   int                   preferzerobranch;   /**< how much we prefer to branch on 0.0 first */
   SCIP_Bool             branchminconverror; /**< whether to compute branching point such that the convexification error is minimized after branching on 0.0 */
   SCIP_Bool             addvarbounds;       /**< will variable bounds be added to the cutpool? */
   SCIP_Bool             linfeasshift;       /**< try linear feasibility shift heuristic in CONSCHECK */
   SCIP_Bool             sepainboundsonly;   /**< should tangents only be generated in variable bounds during separation? */
   SCIP_Real             sepanlpmincont;     /**< minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation */

   SCIP_HEUR*            subnlpheur;         /**< a pointer to the subnlp heuristic */
   SCIP_HEUR*            trysolheur;         /**< a pointer to the trysol heuristic */
   SCIP_EVENTHDLR*       eventhdlr;          /**< our handler for bound change events on variable x */
   SCIP_CONSHDLR*        conshdlrindicator;  /**< a pointer to the indicator constraint handler */
   int                   newsoleventfilterpos;/**< filter position of new solution event handler, if catched */
   SCIP_Bool             comparedpairwise;   /**< did we compare absolute power constraints pairwise in this run? */
   SCIP_Bool             sepanlp;            /**< where linearization of the NLP relaxation solution added? */
   SCIP_NODE*            lastenfolpnode;     /**< the node for which enforcement was called the last time (and some constraint was violated) */
   int                   nenfolprounds;      /**< counter on number of enforcement rounds for the current node */
};

/*
 * Propagation rules
 */

enum Proprule
   {
      PROPRULE_1,                          /**< left hand side and bounds on z -> lower bound on x */
      PROPRULE_2,                          /**< left hand side and upper bound on x -> bound on z */
      PROPRULE_3,                          /**< right hand side and bounds on z -> upper bound on x */
      PROPRULE_4,                          /**< right hand side and lower bound on x -> bound on z */
      PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
   };
typedef enum Proprule PROPRULE;

/*
 * Local methods
 */

/** power function for square, that should be faster than using pow(x, 2.0) */
static
DECL_MYPOW(square)
{
   assert(exponent == 2.0);
   return base*base;
}

/** process variable event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_Bool* ispropagated;

   assert(scip  != NULL);
   assert(event != NULL);
   assert(SCIPeventGetType(event) & SCIP_EVENTTYPE_BOUNDTIGHTENED);

   ispropagated = (SCIP_Bool*)eventdata;
   assert(ispropagated != NULL);

   *ispropagated = FALSE;

   return SCIP_OKAY;
}  /*lint !e715*/

/** catch variable bound tightening events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for variables */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* if z is multiaggregated, then bound changes on x could not be propagated, so we do not need to catch them */
   if( SCIPvarGetStatus(consdata->z) != SCIP_VARSTATUS_MULTAGGR )
   {
      eventtype = SCIP_EVENTTYPE_DISABLED;
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
      if( !SCIPisInfinity(scip,  consdata->rhs) )
         eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;

      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->x, eventtype, eventhdlr, (SCIP_EVENTDATA*)&consdata->isxpropagated, &consdata->xeventfilterpos) );

      consdata->isxpropagated = FALSE;
   }
   else
      consdata->isxpropagated = TRUE;

   /* if x is multiaggregated, then bound changes on z could not be propagated, so we do not need to catch them */
   if( SCIPvarGetStatus(consdata->x) != SCIP_VARSTATUS_MULTAGGR )
   {
      eventtype = SCIP_EVENTTYPE_DISABLED;
      if( consdata->zcoef > 0.0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
      }
      else
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
      }

      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->z, eventtype, eventhdlr, (SCIP_EVENTDATA*)&consdata->iszpropagated, &consdata->zeventfilterpos) );
      consdata->iszpropagated = FALSE;
   }
   else
      consdata->iszpropagated = TRUE;

   return SCIP_OKAY;
}

/** drop variable bound tightening events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler for variables */
   SCIP_CONS*            cons                /**< constraint for which to drop bound change events */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EVENTTYPE eventtype;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(eventhdlr != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( SCIPvarGetStatus(consdata->z) != SCIP_VARSTATUS_MULTAGGR )
   {
      eventtype = SCIP_EVENTTYPE_DISABLED;
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
      if( !SCIPisInfinity(scip,  consdata->rhs) )
         eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;

      SCIP_CALL( SCIPdropVarEvent(scip, consdata->x, eventtype, eventhdlr, (SCIP_EVENTDATA*)&consdata->isxpropagated, consdata->xeventfilterpos) );
      consdata->xeventfilterpos = -1;
   }

   if( SCIPvarGetStatus(consdata->x) != SCIP_VARSTATUS_MULTAGGR )
   {
      eventtype = SCIP_EVENTTYPE_DISABLED;
      if( consdata->zcoef > 0.0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
      }
      else
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            eventtype |= SCIP_EVENTTYPE_LBTIGHTENED;
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            eventtype |= SCIP_EVENTTYPE_UBTIGHTENED;
      }

      SCIP_CALL( SCIPdropVarEvent(scip, consdata->z, eventtype, eventhdlr, (SCIP_EVENTDATA*)&consdata->iszpropagated, consdata->zeventfilterpos) );
      consdata->zeventfilterpos = -1;
   }

   assert(consdata->xeventfilterpos == -1);
   assert(consdata->zeventfilterpos == -1);

   return SCIP_OKAY;
}

/** get key of hash element */
static
SCIP_DECL_HASHGETKEY(presolveFindDuplicatesGetKey)
{
   return elem;
}  /*lint !e715*/

/** checks if two constraints have the same x variable, the same exponent, and either the same offset or the same linear variable and are both equality constraint */
static
SCIP_DECL_HASHKEYEQ(presolveFindDuplicatesKeyEQ)
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   assert(consdata1 != NULL);
   assert(consdata2 != NULL);

   if( consdata1->x != consdata2->x )
      return FALSE;

   if( consdata1->exponent != consdata2->exponent )  /*lint !e777*/
      return FALSE;

   if( consdata1->xoffset != consdata2->xoffset && consdata1->z != consdata2->z )  /*lint !e777*/
      return FALSE;

   return TRUE;
}  /*lint !e715*/

/** get value of hash element when comparing on x */
static
SCIP_DECL_HASHKEYVAL(presolveFindDuplicatesKeyVal)
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);

   return ((unsigned int)(size_t)consdata->x << 16) + (unsigned int)(consdata->exponent*0x80);
}  /*lint !e715*/

/** checks if two constraints have the same z variable and the same exponent */
static
SCIP_DECL_HASHKEYEQ(presolveFindDuplicatesKeyEQ2)
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;

   consdata1 = SCIPconsGetData((SCIP_CONS*)key1);
   consdata2 = SCIPconsGetData((SCIP_CONS*)key2);
   assert(consdata1 != NULL);
   assert(consdata2 != NULL);

   if( consdata1->z != consdata2->z )
      return FALSE;

   if( consdata1->exponent != consdata2->exponent )  /*lint !e777*/
      return FALSE;

   return TRUE;
}  /*lint !e715*/

/** get value of hash element when comparing on z */
static
SCIP_DECL_HASHKEYVAL(presolveFindDuplicatesKeyVal2)
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData((SCIP_CONS*)key);
   assert(consdata != NULL);

   return ((unsigned int)(size_t)consdata->z << 16) + (unsigned int)(consdata->exponent*0x80);
}  /*lint !e715*/

/** upgrades a signpower constraint to a linear constraint if a second signpower constraint with same nonlinear term is available */
static
SCIP_RETCODE presolveFindDuplicatesUpgradeCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons1,              /**< constraint to upgrade to a linear constraint */
   SCIP_CONS*            cons2,              /**< constraint which defines a relation for x|x|^{n-1} */
   SCIP_Bool*            infeas,             /**< buffer where to indicate if infeasibility has been detected */
   int*                  nupgdconss,         /**< buffer where to add number of upgraded conss */
   int*                  ndelconss,          /**< buffer where to add number of deleted conss */
   int*                  naggrvars           /**< buffer where to add number of aggregated variables */
   )
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;
   SCIP_CONS* lincons;
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   SCIP_VAR*  vars[2];
   SCIP_Real  coefs[2];

   assert(scip  != NULL);
   assert(cons1 != NULL);
   assert(cons2 != NULL);
   assert(infeas != NULL);
   assert(nupgdconss != NULL);
   assert(ndelconss != NULL);
   assert(naggrvars != NULL);

   consdata1 = SCIPconsGetData(cons1);
   consdata2 = SCIPconsGetData(cons2);
   assert(consdata1 != NULL);
   assert(consdata2 != NULL);

   assert(SCIPisEQ(scip, consdata2->lhs, consdata2->rhs));
   assert(!SCIPisInfinity(scip, consdata2->lhs));
   assert(consdata1->x        == consdata2->x);
   assert(consdata1->exponent == consdata2->exponent);  /*lint !e777*/
   assert(consdata1->xoffset  == consdata2->xoffset);   /*lint !e777*/

   lhs = consdata1->lhs;
   if( !SCIPisInfinity(scip, -lhs) )
      lhs -= consdata2->lhs;
   rhs = consdata1->rhs;
   if( !SCIPisInfinity(scip,  rhs) )
      rhs -= consdata2->lhs;

   vars[0] = consdata1->z;
   vars[1] = consdata2->z;

   coefs[0] =  consdata1->zcoef;
   coefs[1] = -consdata2->zcoef;

   if( SCIPisEQ(scip, lhs, rhs) )
   {
      SCIP_Bool redundant;
      SCIP_Bool aggregated;

      /* try aggregation */
      SCIP_CALL( SCIPaggregateVars(scip, consdata1->z, consdata2->z, consdata1->zcoef, -consdata2->zcoef, rhs, infeas, &redundant, &aggregated) );

      /* if infeasibility has been detected, stop here */
      if( *infeas )
         return SCIP_OKAY;

      if( redundant )
      {
         /* if redundant is TRUE, then either the aggregation has been done, or it was redundant */
         if( aggregated )
            ++*naggrvars;

         ++*ndelconss;

         SCIP_CALL( SCIPdelCons(scip, cons1) );
         return SCIP_OKAY;
      }
   }

   /* if aggregation did not succeed, then either because some variable is multi-aggregated or due to numerics or because lhs != rhs
    * we then add a linear constraint instead
    */
   vars[0] = consdata1->z;
   vars[1] = consdata2->z;
   coefs[0] =  consdata1->zcoef;
   coefs[1] = -consdata2->zcoef;

   SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons1), 2, vars, coefs, lhs, rhs,
         SCIPconsIsInitial(cons1), SCIPconsIsSeparated(cons1), SCIPconsIsEnforced(cons1),
         SCIPconsIsChecked(cons1), SCIPconsIsPropagated(cons1),  SCIPconsIsLocal(cons1),
         SCIPconsIsModifiable(cons1), SCIPconsIsDynamic(cons1), SCIPconsIsRemovable(cons1),
         SCIPconsIsStickingAtNode(cons1)) );
   SCIP_CALL( SCIPaddCons(scip, lincons) );
   SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

   SCIP_CALL( SCIPdelCons(scip, cons1) );
   ++*nupgdconss;

   return SCIP_OKAY;
}

/** solves a system of two absolute power equations
 * Given:   (x+xoffset1)|x+xoffset1|^{exponent-1} + zcoef1 * z == rhs1
 *      and (x+xoffset2)|x+xoffset2|^{exponent-1} + zcoef2 * z == rhs2
 * with xoffset1 != xoffset2 and zcoef1 * rhs2 == zcoef2 * rhs1 and exponent == 2,
 * finds values for x and z that satisfy these equations, or reports infeasibility if no solution exists.
 *
 * Multiplying the second equation by -zcoef1/zcoef2 and adding it to the first one gives
 *   (x+xoffset1)|x+xoffset1| - zcoef1/zcoef2 (x+offset2)|x+offset2| == 0
 *
 * If zcoef1 == zcoef2, then there exists, due to monotonicity of x|x|, no x such that
 *   (x+xoffset1)|x+xoffset1| == (x+xoffset2)|x+xoffset2|.
 *
 * In general, for zcoef1 / zcoef2 > 0.0, we get
 *   x = (xoffset2 - xoffset1) / (sqrt(zcoef2 / zcoef1) - 1.0) - xoffset1,
 * and for zcoef1 / zcoef2 < 0.0, we get
 *   x = (xoffset2 - xoffset1) / (-sqrt(-zcoef2 / zcoef1) - 1.0) - xoffset1.
 *
 * This then yields z = (rhs1 - (x+xoffset1)|x+xoffset1|) / zcoef1.
 */
static
void presolveFindDuplicatesSolveEquations(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            infeas,             /**< buffer to indicate if the system of equations has no solution */
   SCIP_Real*            xval,               /**< buffer to store value of x in the solution, if any */
   SCIP_Real*            zval,               /**< buffer to store value of z in the solution, if any */
   SCIP_Real             exponent,           /**< exponent in absolute power equations */
   SCIP_Real             xoffset1,           /**< offset for x in first absolute power equation */
   SCIP_Real             zcoef1,             /**< coefficient of z in first absolute power equation */
   SCIP_Real             rhs1,               /**< right-hand-side in first absolute power equation */
   SCIP_Real             xoffset2,           /**< offset for x in second absolute power equation */
   SCIP_Real             zcoef2,             /**< coefficient of z in second absolute power equation */
   SCIP_Real             rhs2                /**< right-hand-side in second absolute power equation */
   )
{
   assert(scip != NULL);
   assert(infeas != NULL);
   assert(xval != NULL);
   assert(zval != NULL);
   assert(exponent == 2.0);
   assert(!SCIPisEQ(scip, xoffset1, xoffset2));
   assert(SCIPisEQ(scip, zcoef1 * rhs2, zcoef2 * rhs1));
   assert(zcoef1 != 0.0);
   assert(zcoef2 != 0.0);

   if( xoffset2 < xoffset1 )
   {
      presolveFindDuplicatesSolveEquations(scip, infeas, xval, zval, exponent, xoffset2, zcoef2, rhs2, xoffset1, zcoef1, rhs1);
      return;
   }

   if( SCIPisEQ(scip, zcoef1, zcoef2) )
   {
      *infeas = TRUE;
      return;
   }

   *infeas = FALSE;

   if( SCIPisEQ(scip, zcoef1, -zcoef2) )
   {
      *xval = - (xoffset1 + xoffset2) / 2.0;
   }
   else if( zcoef2 * zcoef1 > 0.0 )
   {
      *xval = (xoffset2 - xoffset1) / (sqrt(zcoef2 / zcoef1) - 1.0) - xoffset1;
   }
   else
   {
      assert(zcoef2 * zcoef1 < 0.0);
      *xval = (xoffset2 - xoffset1) / (-sqrt(-zcoef2 / zcoef1) - 1.0) - xoffset1;
   }

   *zval = rhs1 - (*xval + xoffset1) * REALABS(*xval + xoffset1);
   *zval /= zcoef1;

   assert(SCIPisFeasEQ(scip, (*xval + xoffset1) * REALABS(*xval + xoffset1) + zcoef1 * *zval, rhs1));
   assert(SCIPisFeasEQ(scip, (*xval + xoffset2) * REALABS(*xval + xoffset2) + zcoef2 * *zval, rhs2));
}

/** finds and removes duplicates in a set of absolute power constraints */
static
SCIP_RETCODE presolveFindDuplicates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for absolute power constraints */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int*                  nupgdconss,         /**< pointer where to add number of upgraded constraints */
   int*                  ndelconss,          /**< pointer where to add number of deleted constraints */
   int*                  naddconss,          /**< pointer where to add number of added constraints */
   int*                  nfixedvars,         /**< pointer where to add number of fixed variables */
   int*                  naggrvars,          /**< pointer where to add number of aggregated variables */
   SCIP_Bool*            success,            /**< pointer to store whether a duplicate was found (and removed) */
   SCIP_Bool*            infeas              /**< pointer to store whether infeasibility was detected */
   )
{
   SCIP_HASHTABLE*     hashtable;
   SCIP_HASHTABLELIST* hashtablelist;
   SCIP_CONSHDLRDATA*  conshdlrdata;
   int c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(nupgdconss != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(success  != NULL);
   assert(infeas != NULL);

   *success = FALSE;
   *infeas = FALSE;

   if( nconss <= 1 )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check all constraints in the given set for duplicates, dominance, or possible simplifications w.r.t. the x variable */

   SCIP_CALL( SCIPhashtableCreate(&hashtable, SCIPblkmem(scip), SCIPcalcHashtableSize(nconss),
         presolveFindDuplicatesGetKey, presolveFindDuplicatesKeyEQ, presolveFindDuplicatesKeyVal, (void*)scip) );

   for( c = 0; c < nconss && !*infeas; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;

      cons0 = conss[c];  /*lint !e613*/

      assert(!SCIPconsIsModifiable(cons0));  /* absolute power constraints aren't modifiable */
      assert(!SCIPconsIsLocal(cons0));       /* shouldn't have local constraints in presolve */
      assert(SCIPconsIsActive(cons0));       /* shouldn't get inactive constraints here */

      hashtablelist = NULL;

      do
      {
         SCIP_CONSDATA* consdata0;
         SCIP_CONSDATA* consdata1;

         /* get constraint from current hash table with same x variable as cons0 and same exponent */
         cons1 = (SCIP_CONS*)(SCIPhashtableRetrieveNext(hashtable, &hashtablelist, (void*)cons0));
         if( cons1 == NULL )
         {
            /* processed all constraints like cons0 from hash table, so insert cons0 and go to conss[c+1] */
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
            break;
         }

         assert(cons0 != cons1);

         consdata0 = SCIPconsGetData(cons0);
         consdata1 = SCIPconsGetData(cons1);
         assert(consdata0 != NULL);
         assert(consdata1 != NULL);

         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );

         assert(consdata0->x        == consdata1->x);
         assert(consdata0->exponent == consdata1->exponent);  /*lint !e777*/

         if( SCIPisEQ(scip, consdata0->xoffset, consdata1->xoffset) )
         {
            /* we have two constraints with the same (x+offset)|x+offset|^n term */

            /* if both constraints have the same functions; strengthen sides of cons1 and throw cons0 away */
            if( consdata0->z == consdata1->z && SCIPisEQ(scip, consdata0->zcoef, consdata1->zcoef) )
            {
               /* check if side strenghtening would result in inconsistency */
               if( SCIPisGT(scip, consdata0->lhs, consdata1->rhs) || SCIPisGT(scip, consdata1->lhs, consdata0->rhs) )
               {
                  SCIPdebugMessage("<%s> and <%s> are contradictory; declare infeasibility\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
                  *infeas = TRUE;
                  break;
               }

               SCIPdebugMessage("<%s> and <%s> are equivalent; dropping the first\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));

               /* if a side of cons1 gets finite via merging with cons0, then this changes locks and events */
               if( (SCIPisInfinity(scip, -consdata1->lhs) && !SCIPisInfinity(scip, -consdata0->lhs)) ||
                  ( SCIPisInfinity(scip,  consdata1->rhs) && !SCIPisInfinity(scip,  consdata0->rhs)) )
               {
                  SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons1) );
                  SCIP_CALL( SCIPunlockVarCons(scip, consdata1->x, cons1, !SCIPisInfinity(scip, -consdata1->lhs), !SCIPisInfinity(scip, consdata1->rhs)) );
                  if( consdata1->zcoef > 0.0 )
                     SCIP_CALL( SCIPunlockVarCons(scip, consdata1->z, cons1, !SCIPisInfinity(scip, -consdata1->lhs), !SCIPisInfinity(scip,  consdata1->rhs)) );
                  else
                     SCIP_CALL( SCIPunlockVarCons(scip, consdata1->z, cons1, !SCIPisInfinity(scip,  consdata1->rhs), !SCIPisInfinity(scip, -consdata1->lhs)) );

                  consdata1->lhs = MAX(consdata0->lhs, consdata1->lhs);
                  consdata1->rhs = MIN(consdata0->rhs, consdata1->rhs);

                  SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons1) );
                  SCIP_CALL( SCIPlockVarCons(scip, consdata1->x, cons1, !SCIPisInfinity(scip, -consdata1->lhs), !SCIPisInfinity(scip, consdata1->rhs)) );
                  if( consdata1->zcoef > 0.0 )
                     SCIP_CALL( SCIPunlockVarCons(scip, consdata1->z, cons1, !SCIPisInfinity(scip, -consdata1->lhs), !SCIPisInfinity(scip,  consdata1->rhs)) );
                  else
                     SCIP_CALL( SCIPunlockVarCons(scip, consdata1->z, cons1, !SCIPisInfinity(scip,  consdata1->rhs), !SCIPisInfinity(scip, -consdata1->lhs)) );
               }
               else
               {
                  consdata1->lhs = MAX(consdata0->lhs, consdata1->lhs);
                  consdata1->rhs = MIN(consdata0->rhs, consdata1->rhs);
               }

               SCIP_CALL( SCIPdelCons(scip, cons0) );
               ++*ndelconss;
               *success = TRUE;

               break;
            }

            /* if cons1 defines a linear expression for sign(x+offset)|x+offset|^n, use it to replace cons0 by a linear constraint */
            if( SCIPisEQ(scip, consdata1->lhs, consdata1->rhs) )
            {
               SCIPdebugMessage("substitute <%s> in <%s> to make linear constraint\n", SCIPconsGetName(cons1), SCIPconsGetName(cons0));
               SCIP_CALL( presolveFindDuplicatesUpgradeCons(scip, cons0, cons1, infeas, nupgdconss, ndelconss, naggrvars) );

               *success = TRUE;
               break;
            }

            /* if cons0 defines a linear expression for sign(x+offset)|x+offset|^n, use it to replace cons1 by a linear constraint */
            if( SCIPisEQ(scip, consdata0->lhs, consdata0->rhs) )
            {
               SCIPdebugMessage("substitute <%s> in <%s> to make linear constraint\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
               SCIP_CALL( presolveFindDuplicatesUpgradeCons(scip, cons1, cons0, infeas, nupgdconss, ndelconss, naggrvars) );

               SCIP_CALL( SCIPhashtableRemove(hashtable, cons1) );
               *success = TRUE;

               if( *infeas )
                  break;
            }
            else
            {
               /* introduce a new equality constraint for sign(x+offset)|x+offset|^n and use it to replace cons0 and cons1 */
               /* @todo maybe we could be more clever by looking which constraint sides are finite */
               SCIP_VAR*  auxvar;
               SCIP_CONS* auxcons;
               char       name[SCIP_MAXSTRLEN];
               SCIP_VAR*  vars[2];
               SCIP_Real  coefs[2];

               SCIPdebugMessage("introduce new auxvar for signpower(%s+%g, %g) to make <%s> and <%s> linear constraint\n", SCIPvarGetName(consdata0->x), consdata0->exponent, consdata0->xoffset, SCIPconsGetName(cons0), SCIPconsGetName(cons1));

               /* create auxiliary variable to represent sign(x+offset)|x+offset|^n */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "auxvar_abspower%s_%g_%g", SCIPvarGetName(consdata0->x), consdata0->exponent, consdata0->xoffset);
               SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
                     TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
               SCIP_CALL( SCIPaddVar(scip, auxvar) );

               /* create auxiliary constraint auxvar = sign(x+offset)|x+offset|^n */
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "auxcons_abspower%s_%g_%g", SCIPvarGetName(consdata0->x), consdata0->exponent, consdata0->xoffset);
               SCIP_CALL( SCIPcreateConsAbspower(scip, &auxcons, name, consdata0->x, auxvar, consdata0->exponent, consdata0->xoffset, -1.0, 0.0, 0.0,
                     SCIPconsIsInitial(cons0) || SCIPconsIsInitial(cons1),
                     SCIPconsIsSeparated(cons0) || SCIPconsIsSeparated(cons1),
                     SCIPconsIsEnforced(cons0) || SCIPconsIsEnforced(cons1),
                     SCIPconsIsChecked(cons0) || SCIPconsIsChecked(cons1),
                     SCIPconsIsPropagated(cons0) || SCIPconsIsPropagated(cons1),
                     FALSE,
                     FALSE,
                     SCIPconsIsDynamic(cons0) || SCIPconsIsDynamic(cons1),
                     SCIPconsIsRemovable(cons0) || SCIPconsIsRemovable(cons1),
                     SCIPconsIsStickingAtNode(cons0) || SCIPconsIsStickingAtNode(cons1)
                     ) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
               ++*naddconss;

#ifdef SCIP_DEBUG_SOLUTION
               {
                  SCIP_Real xval;

                  SCIP_CALL( SCIPdebugGetSolVal(scip, consdata0->x, &xval) );
                  SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SIGN(xval + consdata0->xoffset) * pow(REALABS(xval + consdata0->xoffset), consdata0->exponent)) );
               }
#endif

               /* create linear constraint equivalent for cons0 */
               vars[0]  = auxvar;
               vars[1]  = consdata0->z;
               coefs[0] = 1.0;
               coefs[1] = consdata0->zcoef;
               SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, SCIPconsGetName(cons0), 2, vars, coefs, consdata0->lhs, consdata0->rhs,
                     SCIPconsIsInitial(cons0), SCIPconsIsSeparated(cons0), SCIPconsIsEnforced(cons0),
                     SCIPconsIsChecked(cons0), SCIPconsIsPropagated(cons0), SCIPconsIsLocal(cons0),
                     SCIPconsIsModifiable(cons0), SCIPconsIsDynamic(cons0), SCIPconsIsRemovable(cons0),
                     SCIPconsIsStickingAtNode(cons0)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
               ++*nupgdconss;

               /* create linear constraint equivalent for cons1 */
               vars[1]  = consdata1->z;
               coefs[1] = consdata1->zcoef;
               SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, SCIPconsGetName(cons1), 2, vars, coefs, consdata1->lhs, consdata1->rhs,
                     SCIPconsIsInitial(cons1), SCIPconsIsSeparated(cons1), SCIPconsIsEnforced(cons1),
                     SCIPconsIsChecked(cons1), SCIPconsIsPropagated(cons1), SCIPconsIsLocal(cons1),
                     SCIPconsIsModifiable(cons1), SCIPconsIsDynamic(cons1), SCIPconsIsRemovable(cons1),
                     SCIPconsIsStickingAtNode(cons1)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );
               ++*nupgdconss;

               SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );

               SCIP_CALL( SCIPdelCons(scip, cons0) );
               SCIP_CALL( SCIPdelCons(scip, cons1) );
               SCIP_CALL( SCIPhashtableRemove(hashtable, cons1) );
               *success = TRUE;

               break;
            }
         }
         else if( consdata0->z == consdata1->z &&
            consdata0->exponent == 2.0 &&
            !SCIPisZero(scip, consdata0->zcoef) &&
            !SCIPisZero(scip, consdata1->zcoef) &&
            SCIPisEQ(scip, consdata0->lhs, consdata0->rhs) &&
            SCIPisEQ(scip, consdata1->lhs, consdata1->rhs) &&
            SCIPisEQ(scip, consdata0->lhs * consdata1->zcoef, consdata1->lhs * consdata0->zcoef) )
         {
            /* If we have two equality constraints with the same variables and the same exponent and compatible constants,
             * then this system of equations should have either no or a single solution.
             * Thus, we can report cutoff or fix the variables to this solution, and forget about the constraints.
             * @todo think about inequalities, differing exponents, and exponents != 2
             */

            SCIP_Real xval;
            SCIP_Real zval;

            assert(consdata0->x == consdata1->x);
            assert(consdata0->exponent == consdata1->exponent);  /*lint !e777*/
            assert(!SCIPisEQ(scip, consdata0->xoffset, consdata1->xoffset));

            presolveFindDuplicatesSolveEquations(scip, infeas, &xval, &zval,
               consdata0->exponent,
               consdata0->xoffset, consdata0->zcoef, consdata0->lhs,
               consdata1->xoffset, consdata1->zcoef, consdata1->lhs);

            if( *infeas )
            {
               SCIPdebugMessage("infeasibility detected while solving the equations, no solution exists\n");
               SCIPdebug( SCIPprintCons(scip, cons0, NULL) );
               SCIPdebug( SCIPprintCons(scip, cons1, NULL) );
               break;
            }

            SCIPdebugMessage("fixing variables <%s>[%g, %g] to %g and <%s>[%g, %g] to %g due to equations\n",
               SCIPvarGetName(consdata0->x), SCIPvarGetLbLocal(consdata0->x), SCIPvarGetUbLocal(consdata0->x), xval,
               SCIPvarGetName(consdata0->z), SCIPvarGetLbLocal(consdata0->z), SCIPvarGetUbLocal(consdata0->z), zval);
            SCIPdebug( SCIPprintCons(scip, cons0, NULL) );
            SCIPdebug( SCIPprintCons(scip, cons1, NULL) );

            if( SCIPvarGetStatus(SCIPvarGetProbvar(consdata0->x)) != SCIP_VARSTATUS_MULTAGGR )
            {
               SCIP_Bool fixed;

               SCIP_CALL( SCIPfixVar(scip, consdata0->x, xval, infeas, &fixed) );
               ++*ndelconss;

               if( fixed )
                  ++*nfixedvars;

               if( *infeas )
               {
                  SCIPdebugMessage("infeasibility detected after fixing <%s>\n", SCIPvarGetName(consdata0->x));
                  break;
               }
            }
            else
            {
               SCIP_CONS* lincons;
               SCIP_Real  one;

               one = 1.0;
               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons0), 1, &consdata0->x, &one, xval, xval,
                     TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
               ++*nupgdconss;
            }

            if( SCIPvarGetStatus(SCIPvarGetProbvar(consdata0->z)) != SCIP_VARSTATUS_MULTAGGR )
            {
               SCIP_Bool fixed;

               SCIP_CALL( SCIPfixVar(scip, consdata0->z, zval, infeas, &fixed) );
               ++*ndelconss;

               if( fixed )
                  ++*nfixedvars;

               if( *infeas )
               {
                  SCIPdebugMessage("infeasibility detected after fixing <%s>\n", SCIPvarGetName(consdata0->z));
                  break;
               }
            }
            else
            {
               SCIP_CONS* lincons;
               SCIP_Real  one;

               one = 1.0;
               SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(cons1), 1, &consdata0->z, &one, zval, zval,
                     TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPaddCons(scip, lincons) );
               SCIP_CALL( SCIPreleaseCons(scip, &lincons) );
               ++*nupgdconss;
            }

            SCIP_CALL( SCIPdelCons(scip, cons0) );
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            SCIP_CALL( SCIPhashtableRemove(hashtable, cons1) );
            *success = TRUE;

            break;
         }

         if( hashtablelist == NULL )
         {
            /* processed all constraints like cons0 from hash table, but cons0 could not be removed, so insert cons0 into hashmap and go to conss[c+1] */
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
            break;
         }
      }
      while( TRUE );  /*lint !e506*/
   }

   /* free hash table */
   SCIPhashtableFree(&hashtable);

   if( *infeas )
      return SCIP_OKAY;


   /* check all constraints in the given set for duplicates, dominance, or possible simplifications w.r.t. the z variable */

   SCIP_CALL( SCIPhashtableCreate(&hashtable, SCIPblkmem(scip), SCIPcalcHashtableSize(nconss),
         presolveFindDuplicatesGetKey, presolveFindDuplicatesKeyEQ2, presolveFindDuplicatesKeyVal2, (void*) scip) );

   for( c = 0; c < nconss && !*infeas; ++c )
   {
      SCIP_CONS* cons0;
      SCIP_CONS* cons1;
      SCIP_CONSDATA* consdata0;

      cons0 = conss[c];  /*lint !e613*/

      assert(!SCIPconsIsModifiable(cons0));  /* absolute power constraints aren't modifiable */
      assert(!SCIPconsIsLocal(cons0));       /* shouldn't have local constraints in presolve */

      /* do not consider constraints that we have deleted in the above loop */
      if( SCIPconsIsDeleted(cons0) )
         continue;
      assert(SCIPconsIsActive(cons0));       /* shouldn't get inactive constraints here */

      consdata0 = SCIPconsGetData(cons0);
      assert(consdata0 != NULL);

      /* consider only equality constraints so far
       * @todo do also something with inequalities
       */
      if( !SCIPisEQ(scip, consdata0->lhs, consdata0->rhs) )
         continue;

      hashtablelist = NULL;

      do
      {
         SCIP_CONSDATA* consdata1;

         /* get constraint from current hash table with same z variable as cons0 and same exponent */
         cons1 = (SCIP_CONS*)(SCIPhashtableRetrieveNext(hashtable, &hashtablelist, (void*)cons0));
         if( cons1 == NULL )
         {
            /* processed all constraints like cons0 from hash table, so insert cons0 and go to conss[c+1] */
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
            break;
         }

         assert(cons0 != cons1);
         assert(!SCIPconsIsDeleted(cons1));

         consdata1 = SCIPconsGetData(cons1);
         assert(consdata1 != NULL);

         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );

         assert(consdata0->z        == consdata1->z);
         assert(consdata0->exponent == consdata1->exponent);  /*lint !e777*/
         assert(SCIPisEQ(scip, consdata1->lhs, consdata1->rhs));
         assert(!SCIPisZero(scip, consdata1->zcoef));

         if( SCIPisEQ(scip, consdata0->lhs*consdata1->zcoef, consdata1->lhs*consdata0->zcoef) )
         {
            /* have two absolute power equations with same z and compatible constants
             * we can then reduce this to one absolute power and one linear equation
             * -> x0 + xoffset0 = signpower(zcoef0/zcoef1, 1/exponent) (x1 + xoffset1)
             * -> keep cons1
             * the latter can be realized as an aggregation (if x0 and x1 are not multiaggregated) or linear constraint
             */
            SCIP_Bool redundant;
            SCIP_Bool aggregated;
            SCIP_Real coef;
            SCIP_Real rhs;

            SCIPdebugMessage("<%s> and <%s> can be reformulated to one abspower and one aggregation\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );

            if( consdata0->exponent == 2.0 )
               coef = SIGN(consdata0->zcoef / consdata1->zcoef) * sqrt(REALABS(consdata0->zcoef / consdata1->zcoef));
            else
               coef = SIGN(consdata0->zcoef / consdata1->zcoef) * pow(REALABS(consdata0->zcoef / consdata1->zcoef), 1.0/consdata0->exponent);
            rhs = coef * consdata1->xoffset - consdata0->xoffset;

            /* try aggregation */
            SCIP_CALL( SCIPaggregateVars(scip, consdata0->x, consdata1->x, 1.0, -coef, rhs, infeas, &redundant, &aggregated) );
            if( *infeas )
            {
               /* if infeasibility has been detected, stop here */
               break;
            }
            else if( redundant )
            {
               /* if redundant is TRUE, then either the aggregation has been done, or it was redundant */
               if( aggregated )
                  ++*naggrvars;

               ++*ndelconss;
            }
            else
            {
               /* if aggregation did not succeed, then either because some variable is multi-aggregated or due to numerics
                * we then add a linear constraint instead
                */
               SCIP_CONS* auxcons;
               SCIP_VAR* vars[2];
               SCIP_Real coefs[2];

               vars[0] = consdata0->x;
               vars[1] = consdata1->x;
               coefs[0] = 1.0;
               coefs[1] = -coef;

               /* create linear constraint equivalent for cons0 */
               SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, SCIPconsGetName(cons0), 2, vars, coefs, rhs, rhs,
                     SCIPconsIsInitial(cons0), SCIPconsIsSeparated(cons0), SCIPconsIsEnforced(cons0),
                     SCIPconsIsChecked(cons0), SCIPconsIsPropagated(cons0), SCIPconsIsLocal(cons0),
                     SCIPconsIsModifiable(cons0), SCIPconsIsDynamic(cons0), SCIPconsIsRemovable(cons0),
                     SCIPconsIsStickingAtNode(cons0)) );
               SCIP_CALL( SCIPaddCons(scip, auxcons) );
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, auxcons, NULL) ) );
               SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );

               ++*nupgdconss;
            }
            SCIP_CALL( SCIPdelCons(scip, cons0) );

            *success = TRUE;
            break;
         }

         if( hashtablelist == NULL )
         {
            /* processed all constraints like cons0 from hash table, but cons0 could not be removed, so insert cons0 into hashmap and go to conss[c+1] */
            SCIP_CALL( SCIPhashtableInsert(hashtable, (void*) cons0) );
            break;
         }
      }
      while( TRUE );  /*lint !e506*/
   }

   /* free hash table */
   SCIPhashtableFree(&hashtable);

   return SCIP_OKAY;
}

/** given a variable and an interval, tightens the local bounds of this variable to the given interval */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_VAR*             var,               /**< variable which bounds to tighten */
   SCIP_INTERVAL         bounds,            /**< new bounds */
   SCIP_Bool             force,             /**< force tightening even if below bound strengthening tolerance */
   SCIP_CONS*            cons,              /**< constraint that is propagated */
   SCIP_RESULT*          result,            /**< pointer to store the result of the propagation call */
   int*                  nchgbds,           /**< buffer where to add the number of changed bounds */
   int*                  nfixedvars,        /**< buffer where to add the number of fixed variables, can be equal to nchgbds */
   int*                  naddconss          /**< buffer where to add the number of added constraints, can be NULL if force is FALSE */
   )
{
   SCIP_Bool infeas;
   SCIP_Bool tightened;

   assert(scip    != NULL);
   assert(var     != NULL);
   assert(cons    != NULL);
   assert(result  != NULL);
   assert(nchgbds != NULL);
   assert(nfixedvars != NULL);

   *result = SCIP_DIDNOTFIND;

   if( SCIPisInfinity(scip, SCIPintervalGetInf(bounds)) || SCIPisInfinity(scip, -SCIPintervalGetSup(bounds)) )
   {
      /* domain outside [-infty, +infty] -> declare as infeasible */
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* if variable is not multiaggregated (or aggregated to a multiaggregated), then try SCIPfixVar or SCIPtightenVarLb/Ub
    * otherwise, if bound tightening is forced, add a linear constraint
    * otherwise, forget about the bound tightening
    */
   if( SCIPvarIsActive(SCIPvarGetProbvar(var)) )
   {
      /* check if variable can be fixed */
      if( SCIPisEQ(scip, bounds.inf, bounds.sup) )
      {
         if( !SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            /* if variable not fixed yet, then do so now */
            SCIP_Real fixval;

            fixval = SCIPselectSimpleValue(bounds.inf - SCIPepsilon(scip), bounds.sup + SCIPepsilon(scip), MAXDNOM);
            SCIP_CALL( SCIPfixVar(scip, var, fixval, &infeas, &tightened) );

            if( infeas )
            {
               SCIPdebugMessage("found <%s> infeasible due to fixing variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(var));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            if( tightened )
            {
               SCIPdebugMessage("fixed variable <%s> in constraint <%s> to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetLbLocal(var));
               ++*nfixedvars;
               *result = SCIP_REDUCEDDOM;
            }
         }
         else
         {
            /* only check if new fixing value is consistent with variable bounds, otherwise cutoff */
            if( SCIPisLT(scip, bounds.sup, SCIPvarGetUbLocal(var)) || SCIPisGT(scip, bounds.inf, SCIPvarGetLbLocal(var)) )
            {
               SCIPdebugMessage("found <%s> infeasible due to fixing fixed variable <%s>[%.20g,%.20g] to [%.20g,%.20g]\n",
                  SCIPconsGetName(cons), SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), bounds.inf, bounds.sup);
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
         }

         return SCIP_OKAY;
      }

      /* check if lower bound can be tightened */
      if( SCIPintervalGetInf(bounds) > SCIPvarGetLbLocal(var) )
      {
         assert(!SCIPisInfinity(scip, -SCIPintervalGetInf(bounds)));
         SCIP_CALL( SCIPtightenVarLb(scip, var, SCIPintervalGetInf(bounds), force, &infeas, &tightened) );
         if( infeas )
         {
            SCIPdebugMessage("found %s infeasible due to domain propagation for variable %s in constraint %s\n", SCIPconsGetName(cons), SCIPvarGetName(var), SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         if( tightened )
         {
            SCIPdebugMessage("tightened lower bound of variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetLbLocal(var));
            ++*nchgbds;
            *result = SCIP_REDUCEDDOM;
         }
      }

      /* check if upper bound can be tightened */
      if( SCIPintervalGetSup(bounds) < SCIPvarGetUbLocal(var) )
      {
         assert(!SCIPisInfinity(scip, SCIPintervalGetSup(bounds)));
         SCIP_CALL( SCIPtightenVarUb(scip, var, SCIPintervalGetSup(bounds), force, &infeas, &tightened) );
         if( infeas )
         {
            SCIPdebugMessage("found %s infeasible due to domain propagation for linear variable %s in constraint %s\n", SCIPconsGetName(cons), SCIPvarGetName(var), SCIPconsGetName(cons));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         if( tightened )
         {
            SCIPdebugMessage("tightened upper bound of variable %s in constraint %s to %g\n", SCIPvarGetName(var), SCIPconsGetName(cons), SCIPvarGetUbLocal(var));
            ++*nchgbds;
            *result = SCIP_REDUCEDDOM;
         }
      }
   }
   else if( force && (SCIPisLT(scip, SCIPvarGetLbLocal(var), bounds.inf) || SCIPisGT(scip, SCIPvarGetUbLocal(var), bounds.sup)) )
   {
      /* add a linear constraint bounds.inf <= x <= bounds.sup */
      SCIP_CONS* auxcons;
      SCIP_Bool local;
      SCIP_Real one;

      assert(naddconss != NULL);

      /* we add constraint as local constraint if we are during probing or if we are during solve and not at the root node */
      local = SCIPinProbing(scip) || (SCIPgetStage(scip) == SCIP_STAGE_SOLVING && (SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0));

      one = 1.0;
      SCIP_CALL( SCIPcreateConsLinear(scip, &auxcons, SCIPconsGetName(cons), 1, &var, &one, bounds.inf, bounds.sup,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), local,
            FALSE, FALSE, TRUE, FALSE) );

      if( local )
      {
         SCIP_CALL( SCIPaddConsLocal(scip, auxcons, NULL) );
      }
      else
      {
         SCIP_CALL( SCIPaddCons(scip, auxcons) );
      }
      SCIP_CALL( SCIPreleaseCons(scip, &auxcons) );

      ++*naddconss;
      *result = SCIP_CONSADDED;
   }

   return SCIP_OKAY;
}

/** computes bounds on z in a absolute power constraints for given bounds on x */
static
void computeBoundsZ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_INTERVAL         xbnds,              /**< bounds on x that are to be propagated */
   SCIP_INTERVAL*        zbnds               /**< buffer to store corresponding bounds on z */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real bnd;
   SCIP_Real x;

   assert(scip  != NULL);
   assert(cons  != NULL);
   assert(zbnds != NULL);
   assert(!SCIPintervalIsEmpty(xbnds));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPintervalSetEntire(SCIPinfinity(scip), zbnds);

   /* apply zcoef*z <= rhs - signedpow(xbnds.inf + offset, n) */
   if( !SCIPisInfinity(scip,  consdata->rhs) && !SCIPisInfinity(scip, -xbnds.inf) )
   {
      x = xbnds.inf - PROPVARTOL + consdata->xoffset;
      bnd = consdata->rhs + PROPSIDETOL - SIGN(x) * consdata->pow(REALABS(x), consdata->exponent);

      if( consdata->zcoef > 0.0 )
         zbnds->sup = bnd / consdata->zcoef;
      else
         zbnds->inf = bnd / consdata->zcoef;
   }

   /* apply zcoef*z >= lhs - signedpow(xbnds.sup + offset, n) */
   if( !SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip,  xbnds.sup) )
   {
      x = xbnds.sup + PROPVARTOL + consdata->xoffset;
      bnd = consdata->lhs - PROPSIDETOL - SIGN(x) * consdata->pow(REALABS(x), consdata->exponent);

      if( consdata->zcoef > 0.0 )
         zbnds->inf = bnd / consdata->zcoef;
      else
         zbnds->sup = bnd / consdata->zcoef;
   }

   SCIPdebugMessage("given x = [%.20g, %.20g], computed z = [%.20g, %.20g] via", xbnds.inf, xbnds.sup, zbnds->inf, zbnds->sup);
   SCIPdebug( SCIPprintCons(scip, cons, NULL) );

   assert(!SCIPintervalIsEmpty(*zbnds));
}

/** computes bounds on x in a absolute power constraints for given bounds on z */
static
void computeBoundsX(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_INTERVAL         zbnds,              /**< bounds on x that are to be propagated */
   SCIP_INTERVAL*        xbnds               /**< buffer to store corresponding bounds on z */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real bnd;
   SCIP_Real z;

   assert(scip  != NULL);
   assert(cons  != NULL);
   assert(xbnds != NULL);
   assert(!SCIPintervalIsEmpty(zbnds));

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPintervalSetEntire(SCIPinfinity(scip), xbnds);

   /* apply signedpow(x+offset, n) <= rhs - (zcoef * zbnds).inf */
   z = (consdata->zcoef > 0.0 ? zbnds.inf : zbnds.sup);
   if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, REALABS(z)) )
   {
      bnd = consdata->rhs + PROPSIDETOL - consdata->zcoef * z + REALABS(consdata->zcoef) * PROPVARTOL;
      if( consdata->exponent == 2.0 )
         bnd = SIGN(bnd) * sqrt(REALABS(bnd));
      else
         bnd = SIGN(bnd) * pow(REALABS(bnd), 1.0/consdata->exponent);
      xbnds->sup = bnd - consdata->xoffset;
   }

   /* apply signedpow(x+offset, n) >= lhs - (zcoef * zbnds).sup */
   z = (consdata->zcoef > 0.0 ? zbnds.sup : zbnds.inf);
   if( !SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, REALABS(z)) )
   {
      bnd = consdata->lhs - PROPSIDETOL - consdata->zcoef * z - REALABS(consdata->zcoef) * PROPVARTOL;
      if( consdata->exponent == 2.0 )
         bnd = SIGN(bnd) * sqrt(REALABS(bnd));
      else
         bnd = SIGN(bnd) * pow(REALABS(bnd), 1.0/consdata->exponent);
      xbnds->inf = bnd - consdata->xoffset;
   }

   SCIPdebugMessage("given z = [%.20g, %.20g], computed x = [%.20g, %.20g] via", zbnds.inf, zbnds.sup, xbnds->inf, xbnds->sup);
   SCIPdebug( SCIPprintCons(scip, cons, NULL) );

   assert(!SCIPintervalIsEmpty(*xbnds));
}

/** checks if x or z is fixed and replaces them or deletes constraint */
static
SCIP_RETCODE checkFixedVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for absolute power constraints */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  ndelconss,          /**< counter for number of deleted constraints */
   int*                  nupgdconss,         /**< counter for number of upgraded constraints */
   int*                  nchgbds,            /**< counter for number of variable bound changes */
   int*                  nfixedvars,         /**< counter for number of variable fixations */
   SCIP_RESULT*          result              /**< to store result if we detect infeasibility or remove constraint */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real scalar;
   SCIP_Real constant;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(ndelconss  != NULL);
   assert(nupgdconss != NULL);
   assert(nchgbds    != NULL);
   assert(nfixedvars != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *result = SCIP_DIDNOTFIND;

   if( !SCIPvarIsActive(consdata->x) && SCIPvarGetStatus(consdata->x) != SCIP_VARSTATUS_MULTAGGR )
   {
      /* replace x variable */

      /* get relation x = scalar * var + constant */
      var = consdata->x;
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, &scalar, &constant) );

      if( scalar == 0.0 )
      {
         SCIP_INTERVAL xbnds;
         SCIP_INTERVAL zbnds;
         int naddconss;

         naddconss = 0;

         /* x has been fixed to constant */
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->x), constant));
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(consdata->x), constant));

         /* compute corresponding bounds on z */
         SCIPintervalSet(&xbnds, constant);
         computeBoundsZ(scip, cons, xbnds, &zbnds);

         SCIPdebugMessage("in cons <%s>: x = <%s> fixed to %g -> tighten <%s> to [%g, %g]\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->x), constant, SCIPvarGetName(consdata->z), zbnds.inf, zbnds.sup);

         if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
         {
            /* if sides are equal, then we should either fix z, or declare infeasibility */
            if( SCIPisFeasLT(scip, SCIPvarGetUbGlobal(consdata->z), zbnds.inf) || SCIPisFeasGT(scip, SCIPvarGetLbGlobal(consdata->z), zbnds.sup) )
            {
               SCIPdebugMessage("bounds inconsistent -> cutoff\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            else
            {
               /* compute fixing value for z as value corresponding to fixing of x, projected onto bounds of z */
               SCIP_Real zfix;

               zfix = consdata->rhs - SIGN(constant + consdata->xoffset) * consdata->pow(REALABS(constant + consdata->xoffset), consdata->exponent);
               zfix /= consdata->zcoef;
               assert(SCIPisLE(scip, zbnds.inf, zfix));
               assert(SCIPisGE(scip, zbnds.sup, zfix));
               zfix = MIN(SCIPvarGetUbGlobal(consdata->z), MAX(SCIPvarGetLbGlobal(consdata->z), zfix));  /*lint !e666*/

               zbnds.inf = zfix;
               zbnds.sup = zfix;
               SCIP_CALL( tightenBounds(scip, consdata->z, zbnds, TRUE, cons, result, nchgbds, nfixedvars, &naddconss) );
            }
         }
         else
         {
            /* tighten bounds on z accordingly */
            SCIP_CALL( tightenBounds(scip, consdata->z, zbnds, TRUE, cons, result, nchgbds, nfixedvars, &naddconss) );
         }

         /* delete constraint */
         SCIP_CALL( SCIPdelCons(scip, cons) );

         /* if tightenBounds added a constraint (because z was multiaggregated), then count this as constraint upgrade, otherwise as constraint deletion */
         if( naddconss > 0 )
            ++*nupgdconss;
         else
            ++*ndelconss;

         return SCIP_OKAY;
      }

      SCIPdebugMessage("in cons <%s>: x = <%s> replaced by %g*<%s> + %g\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->x), scalar, SCIPvarGetName(var), constant);

      /* we drop here the events for both variables, because if x is replaced by a multiaggregated variable here, then we do not need to catch bound tightenings on z anymore */
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      SCIP_CALL( SCIPunlockVarCons(scip, consdata->x, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );

      consdata->x = var;
      if( SCIPvarIsActive(consdata->x) )
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, consdata->x) );
      }

      /* add constant to offset */
      consdata->xoffset += constant;

      /* divide constraint by scalar*pow(|scalar|,exponent-1), if not 1.0 */
      if( scalar == 1.0 ) ;
      else if( scalar > 0.0 )
      {
         SCIP_Real factor;

         factor = consdata->pow(scalar, consdata->exponent);
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhs /= factor;
         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->rhs /= factor;
         consdata->zcoef /= factor;
         consdata->xoffset /= scalar;
      }
      else
      {
         SCIP_Real factor;
         SCIP_Real oldlhs;

         assert(scalar < 0.0);

         factor = -consdata->pow(-scalar, consdata->exponent);
         assert(factor < 0.0);

         oldlhs = consdata->lhs;

         if( !SCIPisInfinity(scip,  consdata->rhs) )
            consdata->lhs = consdata->rhs / factor;
         else
            consdata->lhs = -SCIPinfinity(scip);
         if( !SCIPisInfinity(scip, -oldlhs) )
            consdata->rhs = oldlhs / factor;
         else
            consdata->rhs = SCIPinfinity(scip);
         consdata->zcoef /= factor;
         consdata->xoffset /= scalar;
         /* since we flip both constraint sides and the sign of zcoef, the events catched for z remain the same, so update necessary there */
      }

      SCIP_CALL( SCIPlockVarCons(scip, consdata->x, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip, consdata->rhs)) );
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );

      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

      /* rerun constraint comparison */
      conshdlrdata->comparedpairwise = FALSE;
   }

   if( !SCIPvarIsActive(consdata->z) && SCIPvarGetStatus(consdata->z) != SCIP_VARSTATUS_MULTAGGR )
   {
      /* replace z variable */

      /* get relation z = scalar * var + constant */
      var = consdata->z;
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&var, &scalar, &constant) );

      if( scalar == 0.0 )
      {
         SCIP_INTERVAL xbnds;
         SCIP_INTERVAL zbnds;
         int naddconss;

         naddconss = 0;

         /* z has been fixed to constant */
         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(consdata->z), constant));
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(consdata->z), constant));

         /* compute corresponding bounds on x */
         SCIPintervalSet(&zbnds, constant);
         computeBoundsX(scip, cons, zbnds, &xbnds);

         SCIPdebugMessage("in cons <%s>: z = <%s> fixed to %g -> tighten <%s> to [%g, %g]\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->z), constant, SCIPvarGetName(consdata->x), xbnds.inf, xbnds.sup);

         if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
         {
            /* if sides are equal, then we should either fix x, or declare infeasibility */
            if( SCIPisFeasLT(scip, SCIPvarGetUbGlobal(consdata->x), xbnds.inf) || SCIPisFeasGT(scip, SCIPvarGetLbGlobal(consdata->x), xbnds.sup) )
            {
               SCIPdebugMessage("bounds inconsistent -> cutoff\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            else
            {
               /* compute fixing value for x as value corresponding to fixing of z, projected onto bounds of x */
               SCIP_Real xfix;

               xfix = consdata->rhs - consdata->zcoef * constant;
               if( consdata->exponent == 2.0 )
                  xfix = SIGN(xfix) * sqrt(REALABS(xfix)) - consdata->xoffset;
               else
                  xfix = SIGN(xfix) * pow(REALABS(xfix), 1.0/consdata->exponent) - consdata->xoffset;
               assert(SCIPisLE(scip, xbnds.inf, xfix));
               assert(SCIPisGE(scip, xbnds.sup, xfix));
               xfix = MIN(SCIPvarGetUbGlobal(consdata->x), MAX(SCIPvarGetLbGlobal(consdata->x), xfix));  /*lint !e666*/

               xbnds.inf = xfix;
               xbnds.sup = xfix;
               SCIP_CALL( tightenBounds(scip, consdata->x, xbnds, TRUE, cons, result, nchgbds, nfixedvars, &naddconss) );
            }
         }
         else
         {
            /* tighten bounds on x accordingly */
            SCIP_CALL( tightenBounds(scip, consdata->x, xbnds, TRUE, cons, result, nchgbds, nfixedvars, &naddconss) );
         }

         /* delete constraint */
         SCIP_CALL( SCIPdelCons(scip, cons) );

         /* if tightenBounds added a constraint (because x was multiaggregated), then count this as constraint upgrade, otherwise as constraint deletion */
         if( naddconss > 0 )
            ++*nupgdconss;
         else
            ++*ndelconss;

         return SCIP_OKAY;
      }

      SCIPdebugMessage("in cons <%s>: z = <%s> replaced by %g*<%s> + %g\n", SCIPconsGetName(cons), SCIPvarGetName(consdata->z), scalar, SCIPvarGetName(var), constant);

      /* we drop here the events for both variables, because if z is replaced by a multiaggregated variable here, then we do not need to catch bound tightenings on x anymore */
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      if( consdata->zcoef > 0.0 )
         SCIP_CALL( SCIPunlockVarCons(scip, consdata->z, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
      else
         SCIP_CALL( SCIPunlockVarCons(scip, consdata->z, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );

      consdata->z = var;
      if( SCIPvarIsActive(consdata->z) )
      {
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, consdata->z) );
      }

      /* substract constant from constraint sides */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
         consdata->lhs -= consdata->zcoef * constant;
      if( !SCIPisInfinity(scip,  consdata->rhs) )
         consdata->rhs -= consdata->zcoef * constant;

      /* multiply zcoef by scalar */
      consdata->zcoef *= scalar;

      if( consdata->zcoef > 0.0 )
         SCIP_CALL( SCIPlockVarCons(scip, consdata->z, cons, !SCIPisInfinity(scip, -consdata->lhs), !SCIPisInfinity(scip,  consdata->rhs)) );
      else
         SCIP_CALL( SCIPlockVarCons(scip, consdata->z, cons, !SCIPisInfinity(scip,  consdata->rhs), !SCIPisInfinity(scip, -consdata->lhs)) );
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );

      /* rerun constraint comparison */
      conshdlrdata->comparedpairwise = FALSE;
   }

   assert(
      (SCIPvarIsActive(consdata->x) || SCIPvarGetStatus(consdata->x) == SCIP_VARSTATUS_MULTAGGR) &&
      (SCIPvarIsActive(consdata->z) || SCIPvarGetStatus(consdata->z) == SCIP_VARSTATUS_MULTAGGR));

   return SCIP_OKAY;
}

/** gets maximal absolute value in gradient of quadratic function
 * thus, gives \f$max(n |x+offset|^{n-1}, |zcoef|)\f$.
 */
static
SCIP_Real getGradientMaxElement(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol                 /**< solution or NULL if LP solution should be used */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xval;
   SCIP_Real val;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   xval = SCIPgetSolVal(scip, sol, consdata->x);
   assert(!SCIPisInfinity(scip, REALABS(xval)));

   if( consdata->exponent == 2.0 )
      val = consdata->exponent * REALABS(xval + consdata->xoffset);
   else
      val = consdata->exponent * pow(REALABS(xval + consdata->xoffset), consdata->exponent - 1.0);

   return MAX(val, REALABS(consdata->zcoef));  /*lint !e666*/
}

/** computes violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   SCIP_Real*            viol                /**< buffer to store absolute (unscaled) constraint violation */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real val;
   SCIP_Real xval;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(viol != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   xval = SCIPgetSolVal(scip, sol, consdata->x);

   if( SCIPisInfinity(scip, REALABS(xval)) )
   {
      consdata->lhsviol = (SCIPisInfinity(scip, -consdata->lhs) ? 0.0 : SCIPinfinity(scip));
      consdata->rhsviol = (SCIPisInfinity(scip,  consdata->rhs) ? 0.0 : SCIPinfinity(scip));

      return SCIP_OKAY;
   }

   xval += consdata->xoffset;

   val  = SIGN(xval) * consdata->pow(REALABS(xval), consdata->exponent);
   val += consdata->zcoef * SCIPgetSolVal(scip, sol, consdata->z);

   if( val < consdata->lhs && !SCIPisInfinity(scip, -consdata->lhs) )
      consdata->lhsviol = *viol = consdata->lhs - val;
   else
      consdata->lhsviol = 0.0;

   if( val > consdata->rhs && !SCIPisInfinity(scip,  consdata->rhs) )
      consdata->rhsviol = *viol = val - consdata->rhs;
   else
      consdata->rhsviol = 0.0;

   if( consdata->lhsviol > 0.0 || consdata->rhsviol > 0.0 )
   {
      SCIP_Real norm;

      norm = getGradientMaxElement(scip, cons, sol);
      if( norm > 1.0 )
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

   *maxviolcon = NULL;

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      SCIP_CALL( computeViolation(scip, conss[c], sol, &viol) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      viol = MAX(consdata->lhsviol, consdata->rhsviol);
      if( viol > maxviol && SCIPisGT(scip, viol, SCIPfeastol(scip)) )
      {
         maxviol = viol;
         *maxviolcon = conss[c];
      }
   }

   return SCIP_OKAY;
}

/** proposes branching point for constraint */
static
SCIP_Real proposeBranchingPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint which variable to get branching point for */
   int                   preferzero,         /**< how much we prefer branching on -xoffset (0, 1, or 2) if sign is not fixed */
   SCIP_Bool             branchminconverror  /**< whether to minimize convexification error if sign is fixed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*      x;
   SCIP_Real      xref;
   SCIP_Real      xlb;
   SCIP_Real      xub;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   x = consdata->x;
   xlb = SCIPvarGetLbLocal(x);
   xub = SCIPvarGetUbLocal(x);

   /* check if sign of x is not fixed yet */
   if( SCIPisLT(scip, xlb, -consdata->xoffset) && SCIPisGT(scip, xub, -consdata->xoffset) )
   {
      /* if preferzero is 0, just return SCIP_INVALID
       * if preferzero is 1, then propose -xoffset if branching on -xoffset would cut off solution in both child nodes, otherwise return SCIP_INVALID
       * if preferzero is >1, then always propose -xoffset
       */
      assert(preferzero >= 0);

      if( preferzero == 0 )
         return SCIP_INVALID;

      if( preferzero > 1 || SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
         return -consdata->xoffset;

      xlb += consdata->xoffset;
      xub += consdata->xoffset;

      xref = SCIPgetVarSol(scip, x) + consdata->xoffset;
      if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         /* signpow(x,n,offset) + c*z <= 0 is violated
          *  if we are close to or right of -offset, then branching on -offset gives a convex function on the right branch, this is good
          *  otherwise if branching on -offset yields a violated secant cut in left branch, then current solution would be cutoff there, this is also still good
          */
         if( !SCIPisFeasNegative(scip, xref) || SCIPisFeasPositive(scip, -consdata->pow(-xlb, consdata->exponent)*xref/xlb + consdata->zcoef * SCIPgetVarSol(scip, consdata->z)) )
            return -consdata->xoffset;
         return SCIP_INVALID;
      }

      assert(SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) );
      /* signpow(x,n) + c*z >= 0 is violated
       *  if we are close to or left of zero, then branching on 0.0 gives a concave function on the left branch, this is good
       *  otherwise if branching on 0.0 yields a violated secant cut in right branch, then current solution would be cutoff there, this is also still good
       */
      if( !SCIPisFeasPositive(scip, xref) || SCIPisFeasNegative(scip, -consdata->pow(xub, consdata->exponent)*xref/xub + consdata->zcoef * SCIPgetVarSol(scip, consdata->z)) )
         return -consdata->xoffset;
      return SCIP_INVALID;
   }

   if( branchminconverror )
   {
      /* given x^n with xlb <= x <= xub, then the sum of the integrals between the function and its secant on the left and right branches are minimized
       * for branching on ( (ub^n - lb^n) / (n*(ub - lb)) ) ^ (1/(n-1))
       */
      if( SCIPisGE(scip, xlb, -consdata->xoffset) )
      {
         SCIP_Real ref;
         xlb = MAX(0.0, xlb + consdata->xoffset);
         xub = MAX(0.0, xub + consdata->xoffset);

         ref = (consdata->pow(xub, consdata->exponent) - consdata->pow(xlb, consdata->exponent)) / (consdata->exponent * (xub - xlb));
         ref = pow(ref, 1.0/(consdata->exponent-1.0));
         ref -= consdata->xoffset;
         assert(SCIPisGE(scip, ref, SCIPvarGetLbLocal(x)));
         assert(SCIPisLE(scip, ref, SCIPvarGetUbLocal(x)));

         return ref;
      }
      else
      {
         SCIP_Real ref;

         assert(SCIPisLE(scip, xub, -consdata->xoffset));

         xlb = MIN(0.0, xlb + consdata->xoffset);
         xub = MIN(0.0, xub + consdata->xoffset);

         ref = (consdata->pow(-xlb, consdata->exponent) - consdata->pow(-xub, consdata->exponent)) / (consdata->exponent * (-xlb + xub));
         ref = -pow(ref, 1.0/(consdata->exponent-1.0));
         ref -= consdata->xoffset;
         assert(SCIPisGE(scip, ref, SCIPvarGetLbLocal(x)));
         assert(SCIPisLE(scip, ref, SCIPvarGetUbLocal(x)));

         return ref;
      }
   }

   return SCIP_INVALID;
}

/** registers branching variable candidates
 * registers x for all violated absolute power constraints where x is not in convex region
 */
static
SCIP_RETCODE registerBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          onlynonfixedsign;
   int                c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *nnotify = 0;

   onlynonfixedsign = conshdlrdata->preferzerobranch == 3;

   do
   {
      for( c = 0; c < nconss; ++c )
      {
         assert(conss[c] != NULL);  /*lint !e613*/

         consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
         assert(consdata != NULL);

         SCIPdebugMessage("cons <%s> violation: %g %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);  /*lint !e613*/

         /* domain propagation should have removed constraints with fixed x */
         assert(!SCIPisRelEQ(scip, SCIPvarGetLbLocal(consdata->x), SCIPvarGetUbLocal(consdata->x)));

         /* skip variables which sign is already fixed, if we are only interested in variables with unfixed sign here */
         if( onlynonfixedsign &&
            (  !SCIPisLT(scip, SCIPvarGetLbLocal(consdata->x), -consdata->xoffset) ||
               !SCIPisGT(scip, SCIPvarGetUbLocal(consdata->x),  consdata->xoffset)) )
            continue;

         /* if the value of x lies in a concave range (i.e., where a secant approximation is used), then register x as branching variable */
         if( (SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) && (SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->x)) || SCIPgetSolVal(scip, NULL, consdata->x) + consdata->xoffset <= -consdata->root * (SCIPvarGetLbLocal(consdata->x) + consdata->xoffset))) ||
            ( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && (SCIPisInfinity(scip,  SCIPvarGetUbLocal(consdata->x)) || SCIPgetSolVal(scip, NULL, consdata->x) + consdata->xoffset >= -consdata->root * (SCIPvarGetUbLocal(consdata->x) + consdata->xoffset))) )
         {
            SCIPdebugMessage("register var <%s> in cons <%s> with violation %g %g\n", SCIPvarGetName(consdata->x), SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);  /*lint !e613*/
            SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->x, MAX(consdata->lhsviol, consdata->rhsviol), proposeBranchingPoint(scip, conss[c], conshdlrdata->preferzerobranch, conshdlrdata->branchminconverror)) );  /*lint !e613*/
            ++*nnotify;
         }
      }

      if( onlynonfixedsign && *nnotify == 0 )
      {
         /* if we could not a variable in a violated constraint which sign is not already fixed, do another round where we consider all variables again */
         onlynonfixedsign = FALSE;
         continue;
      }

   } while( FALSE );

   return SCIP_OKAY;
}

/** registers a variable from a violated constraint as branching candidate that has a large absolute value in the LP relaxation */
static
SCIP_RETCODE registerLargeLPValueVariableForBranching(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_VAR**            brvar               /**< buffer to store branching variable */
   )
{
   SCIP_CONSDATA*      consdata;
   SCIP_Real           val;
   SCIP_Real           brvarval;
   int                 c;

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

      val = SCIPgetSolVal(scip, NULL, consdata->x) + consdata->xoffset;
      if( REALABS(val) > brvarval )
      {
         brvarval = ABS(val);
         *brvar   = consdata->x;
      }
   }

   if( *brvar != NULL )
   {
      SCIP_CALL( SCIPaddExternBranchCand(scip, *brvar, brvarval, SCIP_INVALID) );
   }

   return SCIP_OKAY;
}

/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *  see cons_varbound
 */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that inferred the bound change */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE              proprule,           /**< propagation rule that deduced the bound change */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the changed bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(infervar != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, consdata->zcoef));

   switch( proprule )
   {
   case PROPRULE_1:
      /* lhs <= sign(x+offset)|x+offset|^n + c*z: left hand side and bounds on z -> lower bound on x */
      assert(infervar == consdata->x);
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      if( consdata->zcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->z, bdchgidx) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->z, bdchgidx) );
      }
      break;

   case PROPRULE_2:
      /* lhs <= sign(x+offset)|x+offset|^n + c*z: left hand side and upper bound on x -> bound on z */
      assert(infervar == consdata->z);
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      SCIP_CALL( SCIPaddConflictUb(scip, consdata->x, bdchgidx) );
      break;

   case PROPRULE_3:
      /* sign(x+offset)|x+offset|^n + c*z <= rhs: right hand side and bounds on z -> upper bound on x */
      assert(infervar == consdata->x);
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      if( consdata->zcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->z, bdchgidx) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->z, bdchgidx) );
      }
      break;

   case PROPRULE_4:
      /* sign(x+offset)|x+offset|^n + c*z <= rhs: right hand side and lower bound on x -> bound on z */
      assert(infervar == consdata->z);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      SCIP_CALL( SCIPaddConflictLb(scip, consdata->x, bdchgidx) );
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in absolute power constraint <%s>\n", proprule, SCIPconsGetName(cons));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** analyze infeasibility */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_VAR*             infervar,           /**< variable that was deduced */
   PROPRULE              proprule,           /**< propagation rule that deduced the bound change */
   SCIP_BOUNDTYPE        boundtype           /**< the type of the changed bound (lower or upper bound) */
   )
{
   /* conflict analysis can only be applied in solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip) );

   /* add the bound which got violated */
   if( boundtype == SCIP_BOUNDTYPE_LOWER )
   {
      SCIP_CALL( SCIPaddConflictUb(scip, infervar, NULL) );
   }
   else
   {
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      SCIP_CALL( SCIPaddConflictLb(scip, infervar, NULL) );
   }

   /* add the reason for the violated of the bound */
   SCIP_CALL( resolvePropagation(scip, cons, infervar, proprule, boundtype, NULL) );

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, NULL) );

   return SCIP_OKAY;
}

/** propagation method for absolute power constraint
 * SCIPinferVarXbCons to allow for repropagation
 */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool             canaddcons,         /**< are we allowed to add a linear constraint when enforcing bounds for a multiaggregated variable? */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real zlb;
   SCIP_Real zub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Bool tightened;
   SCIP_Bool tightenedround;
   SCIP_Real minact;
   SCIP_Real maxact;

   assert(conshdlr != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("propagating absolute power constraint <%s>\n", SCIPconsGetName(cons));

   *cutoff = FALSE;

   /* get current bounds of variables */
   xlb = SCIPvarGetLbLocal(consdata->x);
   xub = SCIPvarGetUbLocal(consdata->x);
   zlb = SCIPvarGetLbLocal(consdata->z);
   zub = SCIPvarGetUbLocal(consdata->z);

   /* if some bound is not tightened, tighten bounds of variables as long as possible */
   tightenedround = !consdata->isxpropagated || !consdata->iszpropagated;
   while( tightenedround )
   {
      tightenedround = FALSE;

      /* propagate left hand side inequality: lhs <= (x+offset)*|x+offset|^n + c*z */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         assert(!*cutoff);

         /* propagate bounds on x (if not multiaggregated):
          *  (1) left hand side and bounds on z -> lower bound on x
          */
         if( SCIPvarIsActive(SCIPvarGetProbvar(consdata->x)) )
         {
            /* if z is fixed, first compute new lower bound on x without tolerances
             * if that is feasible, project new lower bound onto current bounds
             *   otherwise, recompute with tolerances and continue as usual
             */
            if( SCIPisFeasEQ(scip, zlb, zub) )
            {
               assert(!SCIPisInfinity(scip, -zlb));
               assert(!SCIPisInfinity(scip,  zub));

               newlb = consdata->lhs - consdata->zcoef * (consdata->zcoef > 0.0 ? zub : zlb);

               /* invert sign(x+offset)|x+offset|^(n-1) = y -> x = sign(y)|y|^(1/n) - offset */
               if( consdata->exponent == 2.0 )
                  newlb = SIGN(newlb) * sqrt(ABS(newlb));
               else
                  newlb = SIGN(newlb) * pow(ABS(newlb), 1.0/consdata->exponent);
               newlb -= consdata->xoffset;

               if( SCIPisFeasGT(scip, newlb, xub) )
               {
                  /* if new lower bound for x would yield cutoff, recompute with tolerances */
                  newlb = consdata->lhs - PROPSIDETOL - consdata->zcoef * (consdata->zcoef > 0.0 ? (zub + PROPVARTOL) : (zlb - PROPVARTOL));

                  /* invert sign(x+offset)|x+offset|^(n-1) = y -> x = sign(y)|y|^(1/n) - offset */
                  if( consdata->exponent == 2.0 )
                     newlb = SIGN(newlb) * sqrt(ABS(newlb));
                  else
                     newlb = SIGN(newlb) * pow(ABS(newlb), 1.0/consdata->exponent);
                  newlb -= consdata->xoffset;
               }
               else
               {
                  /* project new lower bound onto current bounds */
                  newlb = MIN(newlb, xub);
               }
            }
            else
            {
               if( consdata->zcoef > 0.0 )
               {
                  if( !SCIPisInfinity(scip, zub) )
                     newlb = consdata->lhs - PROPSIDETOL - consdata->zcoef * (zub + PROPVARTOL);
                  else
                     newlb = -SCIPinfinity(scip);
               }
               else
               {
                  if( !SCIPisInfinity(scip, -zlb) )
                     newlb = consdata->lhs - PROPSIDETOL - consdata->zcoef * (zlb - PROPVARTOL);
                  else
                     newlb = -SCIPinfinity(scip);
               }

               if( !SCIPisInfinity(scip, -newlb) )
               {
                  /* invert sign(x+offset)|x+offset|^(n-1) = y -> x = sign(y)|y|^(1/n) - offset */
                  if( consdata->exponent == 2.0 )
                     newlb = SIGN(newlb) * sqrt(ABS(newlb));
                  else
                     newlb = SIGN(newlb) * pow(ABS(newlb), 1.0/consdata->exponent);
                  newlb -= consdata->xoffset;
               }
            }

            if( !SCIPisInfinity(scip, -newlb) )
            {
               if( SCIPisLbBetter(scip, newlb, xlb, xub) )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->x), xlb, xub, newlb, xub);
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->x, newlb, cons, (int)PROPRULE_1, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisInfinity(scip, newlb) || SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->x)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->x, PROPRULE_1, SCIP_BOUNDTYPE_LOWER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  xlb = SCIPvarGetLbLocal(consdata->x);
               }
            }
         }

         assert(!*cutoff);

         /* propagate bounds on z:
          *  (2) left hand side and upper bound on x -> bound on z
          */
         if( SCIPvarGetStatus(consdata->z) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, xub) ) /* cannot change bounds of multaggr vars */
         {
            SCIP_Real newbd;

            /* if x is fixed, first compute new bound on z without tolerances
             * if that is feasible, project new bound onto current bounds
             *   otherwise, recompute with tolerances and continue as usual
             */
            if( SCIPisFeasEQ(scip, xlb, xub) )
            {
               newbd  = xub + consdata->xoffset;
               newbd  = consdata->lhs - SIGN(newbd) * consdata->pow(REALABS(newbd), consdata->exponent);
               newbd /= consdata->zcoef;

               if( (consdata->zcoef > 0.0 && SCIPisFeasGT(scip, newbd, zub)) || (consdata->zcoef < 0.0 && SCIPisFeasLT(scip, newbd, zlb)) )
               {
                  /* if infeasible, recompute with tolerances */
                  newbd  = xub + PROPVARTOL + consdata->xoffset;
                  newbd  = consdata->lhs - PROPSIDETOL - SIGN(newbd) * consdata->pow(REALABS(newbd), consdata->exponent);
                  newbd /= consdata->zcoef;
               }
               else
               {
                  /* project onto current bounds of z */
                  newbd = MIN(zub, MAX(zlb, newbd) );
               }
            }
            else
            {
               newbd  = xub + PROPVARTOL + consdata->xoffset;
               newbd  = consdata->lhs - PROPSIDETOL - SIGN(newbd) * consdata->pow(REALABS(newbd), consdata->exponent);
               newbd /= consdata->zcoef;
            }

            if( consdata->zcoef > 0.0 )
            {
               newlb = newbd;
               if( SCIPisLbBetter(scip, newlb, zlb, zub) )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->z), zlb, zub, newlb, zub);
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->z, newlb, cons, (int)PROPRULE_2, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisInfinity(scip, newlb) || SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->z)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->z, PROPRULE_2, SCIP_BOUNDTYPE_LOWER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  zlb = SCIPvarGetLbLocal(consdata->z);
               }
            }
            else
            {
               newub = newbd;
               if( SCIPisUbBetter(scip, newub, zlb, zub) )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->z), zlb, zub, zlb, newub);
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->z, newub, cons, (int)PROPRULE_2, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisInfinity(scip, -newub) || SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->z)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->z, PROPRULE_2, SCIP_BOUNDTYPE_UPPER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  zub = SCIPvarGetUbLocal(consdata->z);
               }
            }
         }
      }

      assert(!*cutoff);

      /* propagate right hand side inequality: sign(x+offset)|x+offset|^n + c*z <= rhs */
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         /* propagate bounds on x:
          *  (3) right hand side and bounds on z -> upper bound on x
          */
         if( SCIPvarIsActive(SCIPvarGetProbvar(consdata->x)) ) /* cannot change bounds of multaggr or fixed vars */
         {
            /* if z is fixed, first compute new upper bound on x without tolerances
             * if that is feasible, project new upper bound onto current bounds
             *   otherwise, recompute with tolerances and continue as usual
             */
            if( SCIPisFeasEQ(scip, zlb, zub) )
            {
               assert(!SCIPisInfinity(scip, -zlb));
               assert(!SCIPisInfinity(scip,  zub));

               newub = consdata->rhs - consdata->zcoef * (consdata->zcoef > 0.0 ? zlb : zub);

               /* invert sign(x+offset)|x+offset|^(n-1) = y -> x = sign(y)|y|^(1/n) - offset */
               if( consdata->exponent == 2.0 )
                  newub = SIGN(newub) * sqrt(ABS(newub));
               else
                  newub = SIGN(newub) * pow(ABS(newub), 1.0/consdata->exponent);
               newub -= consdata->xoffset;

               if( SCIPisFeasLT(scip, newub, xlb) )
               {
                  /* if new lower bound for x would yield cutoff, recompute with tolerances */
                  newub = consdata->rhs + PROPSIDETOL - consdata->zcoef * (consdata->zcoef > 0.0 ? (zlb - PROPVARTOL) : (zub + PROPVARTOL));

                  /* invert sign(x+offset)|x+offset|^(n-1) = y -> x = sign(y)|y|^(1/n) - offset */
                  if( consdata->exponent == 2.0 )
                     newub = SIGN(newub) * sqrt(ABS(newub));
                  else
                     newub = SIGN(newub) * pow(ABS(newub), 1.0/consdata->exponent);
                  newub -= consdata->xoffset;
               }
               else
               {
                  /* project new upper bound onto current bounds */
                  newub = MAX(newub, xlb);
               }
            }
            else
            {
               if( consdata->zcoef > 0.0 )
               {
                  if( !SCIPisInfinity(scip, -zlb) )
                     newub = consdata->rhs + PROPSIDETOL - consdata->zcoef * (zlb - PROPVARTOL);
                  else
                     newub = SCIPinfinity(scip);
               }
               else
               {
                  if( !SCIPisInfinity(scip, zub) )
                     newub = consdata->rhs + PROPSIDETOL - consdata->zcoef * (zub + PROPVARTOL);
                  else
                     newub = SCIPinfinity(scip);
               }
               if( !SCIPisInfinity(scip, -newub) )
               {
                  /* invert sign(x+offset)|x+offset|^(n-1) = y -> x = sign(y)|y|^(1/n) - offset */
                  if( consdata->exponent == 2.0 )
                     newub = SIGN(newub) * sqrt(ABS(newub));
                  else
                     newub = SIGN(newub) * pow(ABS(newub), 1.0/consdata->exponent);
                  newub -= consdata->xoffset;
               }
            }

            if( !SCIPisInfinity(scip, newub) )
            {
               if( SCIPisUbBetter(scip, newub, xlb, xub) )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->x), xlb, xub, xlb, newub);
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->x, newub, cons, (int)PROPRULE_3, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisInfinity(scip, -newub) || SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->x)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->x, PROPRULE_3, SCIP_BOUNDTYPE_UPPER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  xub = SCIPvarGetUbLocal(consdata->x);
               }
            }
         }

         assert(!*cutoff);

         /* propagate bounds on z:
          *  (4) right hand side and lower bound on x -> bound on z
          */
         if( SCIPvarGetStatus(consdata->z) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, -xlb) ) /* cannot change bounds of multaggr vars */
         {
            SCIP_Real newbd;

            /* if x is fixed, first compute new bound on z without tolerances
             * if that is feasible, project new bound onto current bounds
             *   otherwise, recompute with tolerances and continue as usual
             */
            if( SCIPisFeasEQ(scip, xlb, xub) )
            {
               newbd  = xlb + consdata->xoffset;
               newbd  = consdata->rhs - SIGN(newbd) * consdata->pow(REALABS(newbd), consdata->exponent);
               newbd /= consdata->zcoef;

               if( (consdata->zcoef > 0.0 && SCIPisFeasLT(scip, newbd, zlb)) || (consdata->zcoef < 0.0 && SCIPisFeasGT(scip, newbd, zub)) )
               {
                  /* if infeasible, recompute with tolerances */
                  newbd  = xlb - PROPVARTOL + consdata->xoffset;
                  newbd  = consdata->rhs + PROPSIDETOL - SIGN(newbd) * consdata->pow(REALABS(newbd), consdata->exponent);
                  newbd /= consdata->zcoef;
               }
               else
               {
                  /* project onto current bounds of z */
                  newbd = MIN(zub, MAX(zlb, newbd) );
               }
            }
            else
            {
               newbd  = xlb - PROPVARTOL + consdata->xoffset;
               newbd  = consdata->rhs + PROPSIDETOL - SIGN(newbd) * consdata->pow(REALABS(newbd), consdata->exponent);
               newbd /= consdata->zcoef;
            }

            if( consdata->zcoef > 0.0 )
            {
               newub = newbd;
               if( SCIPisUbBetter(scip, newub, zlb, zub) )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->z), zlb, zub, zlb, newub);
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->z, newub, cons, (int)PROPRULE_4, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisInfinity(scip, -newub) || SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->z)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->z, PROPRULE_4, SCIP_BOUNDTYPE_UPPER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  zub = SCIPvarGetUbLocal(consdata->z);
               }
            }
            else
            {
               newlb = newbd;
               if( SCIPisLbBetter(scip, newlb, zlb, zub) )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->z), zlb, zub, newlb, zub);
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->z, newlb, cons, (int)PROPRULE_4, FALSE, cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisInfinity(scip, newlb) || SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->z)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->z, PROPRULE_4, SCIP_BOUNDTYPE_LOWER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  zlb = SCIPvarGetLbLocal(consdata->z);
               }
            }
         }
      }

      assert(!*cutoff);
   }

   /* mark the constraint propagated */
   consdata->isxpropagated = TRUE;
   consdata->iszpropagated = TRUE;

   if( *cutoff )
      return SCIP_OKAY;

   /* check for redundancy */
   if( !SCIPisInfinity(scip, -xlb) && !SCIPisInfinity(scip, consdata->zcoef > 0.0 ? -zlb :  zub) )
      minact = SIGN(xlb + consdata->xoffset) * consdata->pow(REALABS(xlb + consdata->xoffset), consdata->exponent) + consdata->zcoef * (consdata->zcoef > 0.0 ? zlb : zub);
   else
      minact = -SCIPinfinity(scip);

   if( !SCIPisInfinity(scip,  xub) && !SCIPisInfinity(scip, consdata->zcoef > 0.0 ?  zub : -zlb) )
      maxact = SIGN(xub + consdata->xoffset) * consdata->pow(REALABS(xub + consdata->xoffset), consdata->exponent) + consdata->zcoef * (consdata->zcoef > 0.0 ? zub : zlb);
   else
      maxact = SCIPinfinity(scip);

   if( SCIPisFeasGE(scip, minact, consdata->lhs) && SCIPisFeasLE(scip, maxact, consdata->rhs) )
   {
      SCIPdebugMessage("absolute power constraint <%s> is redundant: <%s>[%.15g,%.15g], <%s>[%.15g,%.15g]\n",
         SCIPconsGetName(cons),
         SCIPvarGetName(consdata->x), SCIPvarGetLbLocal(consdata->x), SCIPvarGetUbLocal(consdata->x),
         SCIPvarGetName(consdata->z), SCIPvarGetLbLocal(consdata->z), SCIPvarGetUbLocal(consdata->z));

      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      ++*ndelconss;

      return SCIP_OKAY;
   }

   /* delete constraint if x has been fixed */
   if( SCIPisRelEQ(scip, xlb, xub) && (SCIPvarIsActive(consdata->z) || canaddcons) )
   {
      SCIP_RESULT tightenresult;
      SCIP_INTERVAL xbnds;
      SCIP_INTERVAL zbnds;

      SCIPdebugMessage("x-variable in constraint <%s> is fixed: x = <%s>[%.15g,%.15g], z = <%s>[%.15g,%.15g]\n",
         SCIPconsGetName(cons), SCIPvarGetName(consdata->x), xlb, xub, SCIPvarGetName(consdata->z), zlb, zub);

      SCIPintervalSetBounds(&xbnds, MIN(xlb, xub), MAX(xlb, xub));
      computeBoundsZ(scip, cons, xbnds, &zbnds);

      /* in difference to the loop above, here we enforce a possible bound tightening on z, and may add a linear cons if z is multiaggregated */
      SCIP_CALL( tightenBounds(scip, consdata->z, zbnds, TRUE, cons, &tightenresult, nchgbds, nchgbds, naddconss) );
      if( tightenresult == SCIP_CUTOFF )
         *cutoff = TRUE;

      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      ++*ndelconss;

      return SCIP_OKAY;
   }

   /* delete constraint if z has been fixed */
   if( SCIPisRelEQ(scip, zlb, zub) && (SCIPvarIsActive(consdata->x) || canaddcons) )
   {
      SCIP_RESULT tightenresult;
      SCIP_INTERVAL xbnds;
      SCIP_INTERVAL zbnds;

      SCIPdebugMessage("z-variable in constraint <%s> is fixed: x = <%s>[%.15g,%.15g], z = <%s>[%.15g,%.15g]\n",
         SCIPconsGetName(cons), SCIPvarGetName(consdata->x), xlb, xub, SCIPvarGetName(consdata->z), zlb, zub);

      SCIPintervalSetBounds(&zbnds, MIN(zlb, zub), MAX(zlb, zub));
      computeBoundsX(scip, cons, zbnds, &xbnds);

      /* in difference to the loop above, here we enforce a possible bound tightening on x, and may add a linear cons if x is multiaggregated */
      SCIP_CALL( tightenBounds(scip, consdata->x, xbnds, TRUE, cons, &tightenresult, nchgbds, nchgbds, naddconss) );
      if( tightenresult == SCIP_CUTOFF )
         *cutoff = TRUE;

      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      ++*ndelconss;

      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** notifies SCIP about a variable bound lhs <= x + c*y <= rhs */
static
SCIP_RETCODE addVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< absolute power constraint this variable bound is derived form */
   SCIP_VAR*             var,                /**< variable x for which we want to add a variable bound */
   SCIP_VAR*             vbdvar,             /**< variable y which makes the bound a variable bound */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable vbdvar */
   SCIP_Real             lhs,                /**< left  hand side of varbound constraint */
   SCIP_Real             rhs,                /**< right hand side of varbound constraint */
   SCIP_Bool*            infeas,             /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer where to add number of performed bound changes, or NULL */
   )
{
   int nbdchgs_local;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(var  != NULL);
   assert(vbdvar != NULL);
   assert(!SCIPisZero(scip, vbdcoef));
   assert(!SCIPisInfinity(scip, ABS(vbdcoef)));
   assert(infeas != NULL);

   *infeas = FALSE;

   /* make sure vbdvar is active, so we can search for it in SCIPvarGetVxbdVars() */
   if( !SCIPvarIsActive(vbdvar) )
   {
      SCIP_Real constant;

      constant = 0.0;
      SCIP_CALL( SCIPvarGetProbvarSum(&vbdvar, &vbdcoef, &constant) );
      if( !SCIPvarIsActive(vbdvar) || (vbdcoef == 0.0) )
         return SCIP_OKAY;

      if( !SCIPisInfinity(scip, -lhs) )
         lhs -= constant;
      if( !SCIPisInfinity(scip,  rhs) )
         rhs -= constant;
   }

   /* vbdvar should be a non-fixed binary variable */
   assert(SCIPvarIsIntegral(vbdvar));
   assert(SCIPisZero(scip, SCIPvarGetLbGlobal(vbdvar)));
   assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(vbdvar), 1.0));

   SCIPdebugMessage("-> %g <= <%s> + %g*<%s> <= %g\n", lhs, SCIPvarGetName(var), vbdcoef, SCIPvarGetName(vbdvar), rhs);

   if( !SCIPisInfinity(scip, -lhs) )
   {
      SCIP_CALL( SCIPaddVarVlb(scip, var, vbdvar, -vbdcoef, lhs, infeas, &nbdchgs_local) );
      if( *infeas )
         return SCIP_OKAY;
      if( nbdchgs )
         *nbdchgs += nbdchgs_local;
   }

   if( !SCIPisInfinity(scip,  rhs) )
   {
      SCIP_CALL( SCIPaddVarVub(scip, var, vbdvar, -vbdcoef, rhs, infeas, &nbdchgs_local) );
      if( *infeas )
         return SCIP_OKAY;
      if( nbdchgs )
         *nbdchgs += nbdchgs_local;
   }

   return SCIP_OKAY;
}

/** propagates varbounds of variables
 * Let f(x) = sign(x+offset)|x+offset|^n,  f^{-1}(y) = sign(y)|y|^(1/n) - offset.
 * Thus, constraint is lhs <= f(x) + c*z <= rhs.
 *
 * Given a variable bound constraint x <= a*y + b with y a binary variable, one obtains
 * y = 0 -> f(x) <= f(b)   -> lhs <= f(b)   + c*z
 * y = 1 -> f(x) <= f(a+b) -> lhs <= f(a+b) + c*z
 * => lhs <= f(b) + y * (f(a+b)-f(b)) + c*z
 *
 * Given a variable bound constraint x >= a*y + b with y a binary variable, one obtains analogously
 * f(b) + y * (f(a+b)-f(b)) + c*z <= rhs
 *
 * Given a variable bound constraint c*z <= a*y + b with y a binary variable, one obtains
 * y = 0 -> lhs <= f(x) + b   -> x >= f^{-1}(lhs - b)
 * y = 1 -> lhs <= f(x) + a+b -> x >= f^{-1}(lhs - (a+b))
 * => x >= f^{-1}(lhs - b) + y * (f^{-1}(lhs - (a+b)) - f^{-1}(lhs - b))
 *
 * Given a variable bound constraint c*z >= a*y + b with y a binary variable, one obtains analogously
 *    x <= f^{-1}(rhs - b) + y * (f^{-1}(rhs - (a+b)) - f^{-1}(rhs - b))
 */
static
SCIP_RETCODE propagateVarbounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS*            cons,               /**< absolute power constraint */
   SCIP_Bool*            infeas,             /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer where to add number of performed bound changes, or NULL */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*      y;
   SCIP_Real      a;
   SCIP_Real      b;
   SCIP_Real      fb;
   SCIP_Real      fab;
   SCIP_Real      vbcoef;
   SCIP_Real      vbconst;
   int            i;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(cons     != NULL);
   assert(infeas   != NULL);

   *infeas =  FALSE;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->z != NULL);
   assert(!SCIPisZero(scip, consdata->zcoef));

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      /* propagate varbounds x <= a*y+b onto z
       *    lhs <= f(b) + y * (f(a+b)-f(b)) + c*z
       * -> c*z >= lhs-f(b) + y * (f(b)-f(a+b))
       */
      for( i = 0; i < SCIPvarGetNVubs(consdata->x); ++i )
      {
         y = SCIPvarGetVubVars(consdata->x)[i];
         a = SCIPvarGetVubCoefs(consdata->x)[i];
         b = SCIPvarGetVubConstants(consdata->x)[i];

         /* skip variable bound if y is not integer or its valid values are not {0,1}
          * @todo extend to arbitrary integer variables
          */
         if( !SCIPvarIsBinary(y) || SCIPvarGetLbGlobal(y) > 0.5 || SCIPvarGetUbGlobal(y) < 0.5 )
            continue;

         /* skip variable bound if coefficient is very small */
         if( SCIPisFeasZero(scip, consdata->pow(a, consdata->exponent)) )
            continue;

         SCIPdebugMessage("propagate variable bound <%s> <= %g*<%s> + %g\n", SCIPvarGetName(consdata->x), a, SCIPvarGetName(y), b);

         fb  = SIGN(  b + consdata->xoffset) * consdata->pow(  b + consdata->xoffset, consdata->exponent);  /* f(  b) = sign(  b) |  b|^n */
         fab = SIGN(a+b + consdata->xoffset) * consdata->pow(a+b + consdata->xoffset, consdata->exponent);  /* f(a+b) = sign(a+b) |a+b|^n */

         vbcoef  = (fb - fab) / consdata->zcoef;
         vbconst = (consdata->lhs - fb) / consdata->zcoef;

         if( consdata->zcoef > 0.0 )
         {
            /* add varbound z >= (lhs-f(b))/c + y * (f(b)-f(a+b))/c */
            SCIP_CALL( addVarbound(scip, cons, consdata->z, y, -vbcoef, vbconst,  SCIPinfinity(scip), infeas, nbdchgs) );
         }
         else
         {
            /* add varbound z <= (lhs-f(b))/c + y * (f(b)-f(a+b))/c */
            SCIP_CALL( addVarbound(scip, cons, consdata->z, y, -vbcoef, -SCIPinfinity(scip), vbconst, infeas, nbdchgs) );
         }
         if( *infeas )
            return SCIP_OKAY;
      }
   }

   /* propagate varbounds x >= a*y+b onto z
    *    f(b) + y * (f(a+b)-f(b)) + c*z <= rhs
    * -> c*z <= rhs-f(b) + y * (f(b)-f(a+b))
    */
   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      for( i = 0; i < SCIPvarGetNVlbs(consdata->x); ++i )
      {
         y = SCIPvarGetVlbVars(consdata->x)[i];
         a = SCIPvarGetVlbCoefs(consdata->x)[i];
         b = SCIPvarGetVlbConstants(consdata->x)[i];

         /* skip variable bound if y is not integer or its valid values are not {0,1}
          * @todo extend to arbitrary integer variables
          */
         if( !SCIPvarIsBinary(y) || SCIPvarGetLbGlobal(y) > 0.5 || SCIPvarGetUbGlobal(y) < 0.5 )
            continue;

         /* skip variable bound if coefficient is very small */
         if( SCIPisFeasZero(scip, consdata->pow(a, consdata->exponent)) )
            continue;

         SCIPdebugMessage("propagate variable bound <%s> >= %g*<%s> + %g\n", SCIPvarGetName(consdata->x), a, SCIPvarGetName(y), b);

         fb  = SIGN(  b + consdata->xoffset) * consdata->pow(  b + consdata->xoffset, consdata->exponent);  /* f(  b) = sign(  b) |  b|^n */
         fab = SIGN(a+b + consdata->xoffset) * consdata->pow(a+b + consdata->xoffset, consdata->exponent);  /* f(a+b) = sign(a+b) |a+b|^n */

         vbcoef  = (fb - fab) / consdata->zcoef;
         vbconst = (consdata->rhs - fb) / consdata->zcoef;

         if( consdata->zcoef > 0.0 )
         {
            /* add varbound z <= (rhs-f(b))/c + y * (f(b)-f(a+b))/c */
            SCIP_CALL( addVarbound(scip, cons, consdata->z, y, -vbcoef, -SCIPinfinity(scip), vbconst, infeas, nbdchgs) );
         }
         else
         {
            /* add varbound z >= (rhs-f(b))/c + y * (f(b)-f(a+b))/c */
            SCIP_CALL( addVarbound(scip, cons, consdata->z, y, -vbcoef, vbconst,  SCIPinfinity(scip), infeas, nbdchgs) );
         }
         if( *infeas )
            return SCIP_OKAY;
      }
   }

   /* propagate variable upper bounds on z onto x
    * c*z <= a*y+b -> x >= f^{-1}(lhs - b) + y * (f^{-1}(lhs - (a+b)) - f^{-1}(lhs - b))
    * c*z >= a*y+b -> x <= f^{-1}(rhs - b) + y * (f^{-1}(rhs - (a+b)) - f^{-1}(rhs - b))
    */
   if( (consdata->zcoef > 0.0 && !SCIPisInfinity(scip, -consdata->lhs)) ||
      ( consdata->zcoef < 0.0 && !SCIPisInfinity(scip,  consdata->rhs)) )
      for( i = 0; i < SCIPvarGetNVubs(consdata->z); ++i )
      {
         y = SCIPvarGetVubVars(consdata->z)[i];
         a = SCIPvarGetVubCoefs(consdata->z)[i] * consdata->zcoef;
         b = SCIPvarGetVubConstants(consdata->z)[i] * consdata->zcoef;

         SCIPdebugMessage("propagate variable bound %g*<%s> %c= %g*<%s> + %g\n", consdata->zcoef, SCIPvarGetName(consdata->z), consdata->zcoef > 0 ? '<' : '>', a, SCIPvarGetName(y), b);

         /* skip variable bound if y is not integer or its valid values are not {0,1}
          * @todo extend to arbitrary integer variables
          */
         if( !SCIPvarIsBinary(y) || SCIPvarGetLbGlobal(y) > 0.5 || SCIPvarGetUbGlobal(y) < 0.5 )
            continue;

         if( consdata->zcoef > 0.0 )
         {
            fb = consdata->lhs - b;
            fb = SIGN(fb) * pow(ABS(fb), 1.0/consdata->exponent);
            fab = consdata->lhs - (a+b);
            fab = SIGN(fab) * pow(ABS(fab), 1.0/consdata->exponent);
            SCIP_CALL( addVarbound(scip, cons, consdata->x, y, fb - fab, fb - consdata->xoffset, SCIPinfinity(scip),  infeas, nbdchgs) );
         }
         else
         {
            fb = consdata->rhs - b;
            fb = SIGN(fb) * pow(ABS(fb), 1.0/consdata->exponent);
            fab = consdata->rhs - (a+b);
            fab = SIGN(fab) * pow(ABS(fab), 1.0/consdata->exponent);
            SCIP_CALL( addVarbound(scip, cons, consdata->x, y, fb - fab, -SCIPinfinity(scip), fb - consdata->xoffset, infeas, nbdchgs) );
         }
         if( *infeas )
            return SCIP_OKAY;
      }

   /* propagate variable lower bounds on z onto x
    * c*z <= a*y+b -> x >= f^{-1}(lhs - b) + y * (f^{-1}(lhs - (a+b)) - f^{-1}(lhs - b))
    * c*z >= a*y+b -> x <= f^{-1}(rhs - b) + y * (f^{-1}(rhs - (a+b)) - f^{-1}(rhs - b))
    */
   if( (consdata->zcoef < 0.0 && !SCIPisInfinity(scip, -consdata->lhs)) ||
      ( consdata->zcoef > 0.0 && !SCIPisInfinity(scip,  consdata->rhs)) )
      for( i = 0; i < SCIPvarGetNVlbs(consdata->z); ++i )
      {
         y = SCIPvarGetVlbVars(consdata->z)[i];
         a = SCIPvarGetVlbCoefs(consdata->z)[i] * consdata->zcoef;
         b = SCIPvarGetVlbConstants(consdata->z)[i] * consdata->zcoef;

         SCIPdebugMessage("propagate variable bound %g*<%s> %c= %g*<%s> + %g\n", consdata->zcoef, SCIPvarGetName(consdata->z), consdata->zcoef > 0 ? '>' : '<', a, SCIPvarGetName(y), b);

         /* skip variable bound if y is not integer or its valid values are not {0,1}
          * @todo extend to arbitrary integer variables
          */
         if( !SCIPvarIsBinary(y) || SCIPvarGetLbGlobal(y) > 0.5 || SCIPvarGetUbGlobal(y) < 0.5 )
            continue;

         if( consdata->zcoef > 0.0 )
         {
            fb = consdata->rhs - b;
            fb = SIGN(fb) * pow(ABS(fb), 1.0/consdata->exponent);
            fab = consdata->rhs - (a+b);
            fab = SIGN(fab) * pow(ABS(fab), 1.0/consdata->exponent);
            SCIP_CALL( addVarbound(scip, cons, consdata->x, y, fb - fab, -SCIPinfinity(scip), fb - consdata->xoffset, infeas, nbdchgs) );
         }
         else
         {
            fb = consdata->lhs - b;
            fb = SIGN(fb) * pow(ABS(fb), 1.0/consdata->exponent);
            fab = consdata->lhs - (a+b);
            fab = SIGN(fab) * pow(ABS(fab), 1.0/consdata->exponent);
            SCIP_CALL( addVarbound(scip, cons, consdata->x, y, fb - fab, fb - consdata->xoffset,  SCIPinfinity(scip), infeas, nbdchgs) );
         }
         if( *infeas )
            return SCIP_OKAY;
      }

   return SCIP_OKAY;
}

/** adds variable bounds on a variable as LP rows into the cutpool */
static
SCIP_RETCODE addVarboundsToCutPool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable */
   )
{
   SCIP_Real constant;
   SCIP_Real scalar;
   SCIP_ROW* row;
   SCIP_VAR* rowvars[2];
   SCIP_Real rowcoef[2];
   char      name[SCIP_MAXSTRLEN];
   int i;

   assert(scip != NULL);
   assert(var  != NULL);

   constant = 0.0;
   scalar   = 1.0;

   /* make sure var is active */
   SCIP_CALL( SCIPvarGetProbvarSum(&var, &scalar, &constant) );
   if( !SCIPvarIsActive(var) || scalar == 0.0 )
      return SCIP_OKAY;

   rowvars[0] = var;
   rowcoef[0] = scalar;

   for( i = 0; i < SCIPvarGetNVlbs(var); ++i )
   {
      /* variable lower bound scalar*var+constant >= b_iz_i + d_i */
      rowvars[1] = SCIPvarGetVlbVars(var)[i];
      rowcoef[1] = -SCIPvarGetVlbCoefs(var)[i];
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s_varlb", SCIPvarGetName(rowvars[0]), SCIPvarGetName(rowvars[1]));

      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, SCIPvarGetVlbConstants(var)[i] - constant, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPaddVarsToRow(scip, row, 2, rowvars, rowcoef) );

      SCIP_CALL( SCIPaddPoolCut(scip, row) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   for( i = 0; i < SCIPvarGetNVubs(var); ++i )
   {
      /* variable upper bound scalar*var+constant <= b_iz_i + d_i */
      rowvars[1] = SCIPvarGetVubVars(var)[i];
      rowcoef[1] = -SCIPvarGetVubCoefs(var)[i];
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s_varub", SCIPvarGetName(rowvars[0]), SCIPvarGetName(rowvars[1]));

      SCIP_CALL( SCIPcreateEmptyRow(scip, &row, name, -SCIPinfinity(scip), SCIPvarGetVubConstants(var)[i] - constant, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPaddVarsToRow(scip, row, 2, rowvars, rowcoef) );

      SCIP_CALL( SCIPaddPoolCut(scip, row) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** computes linear underestimator for (x+offset)^n + c*z <= rhs by linearization in x
 *
 * the generated cut is xmul * n * (refpoint+offset)^(n-1) * x + c*z <= rhs + ((n-1)*refpoint-offset) * (refpoint+offset)^(n-1)
 */
static
SCIP_RETCODE generateLinearizationCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store row */
   SCIP_Real             refpoint,           /**< base point for linearization */
   SCIP_Real             exponent,           /**< exponent n in sign(x)abs(x)^n */
   SCIP_Real             xoffset,            /**< offset of x */
   SCIP_Real             xmult,              /**< multiplier for coefficient of x */
   SCIP_Real             zcoef,              /**< coefficient of z */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_VAR*             x,                  /**< variable x */
   SCIP_VAR*             z,                  /**< variable z */
   SCIP_Bool             islocal             /**< whether the cut is valid only locally */
   )
{
   SCIP_Real tmp;

   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, refpoint+xoffset));
   assert(!SCIPisInfinity(scip, refpoint));

   if( refpoint < -xoffset )
      refpoint = -xoffset;

   tmp = exponent == 2.0 ? refpoint+xoffset : pow(refpoint+xoffset, exponent-1);
   if( SCIPisInfinity(scip, tmp) )
   {
      SCIPdebugMessage("skip linearization cut because (refpoint+offset)^(exponent-1) > infinity\n");
      *row = NULL;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, "signpowlinearizecut", -SCIPinfinity(scip), SCIPinfinity(scip), islocal, FALSE /* modifiable */, TRUE /* removable */ ) );

   SCIP_CALL( SCIPaddVarToRow(scip, *row, x, xmult*exponent*tmp) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, z, zcoef) );
   SCIP_CALL( SCIPchgRowRhs(scip, *row, rhs + ((exponent-1)*refpoint-xoffset)*tmp) );

   return SCIP_OKAY;
}

/** computes linear underestimator for (x+xoffset)^n + c*z <= rhs by linearization in x
 *
 * the generated cut is xmul * n * (refpoint+offset)^(n-1) * x + c*z <= rhs + ((n-1)*refpoint-offset) * (refpoint+offset)^(n-1)
 * where refpoint is computed by projecting (xref, zref) onto the graph of (x+offset)^n w.r.t. euclidean norm
 *
 * Thus, the projection is computed by minimizing 1/2(x-xref)^2 + 1/2(((x+offset)^n-rhs)/(-c) - zref)^2.
 * I.e., we aim to find a root of
 *   g(x) = x - xref + n/c (x+offset)^(n-1) (zref - rhs/c) + n/c^2 (x+offset)^(2n-1)
 * We do this numerically by executing up to five newton iterations. It is
 *  g'(x) = 1 + n(n-1)/c (x+offset)^(n-2) (zref - rhs/c) + n(2n-1)/c^2 (x+offset)^(2n-2)
 */
static
SCIP_RETCODE generateLinearizationCutProject(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store row */
   SCIP_Real             xref,               /**< reference point for x */
   SCIP_Real             zref,               /**< reference point for z */
   SCIP_Real             xmin,               /**< minimal value x is allowed to take */
   SCIP_Real             exponent,           /**< exponent n in sign(x+offset)abs(x+offset)^n */
   SCIP_Real             xoffset,            /**< offset of x */
   SCIP_Real             xmult,              /**< multiplier for coefficient of x */
   SCIP_Real             zcoef,              /**< coefficient of z */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_VAR*             x,                  /**< variable x */
   SCIP_VAR*             z,                  /**< variable z */
   SCIP_Bool             islocal             /**< whether the cut is valid only locally */
   )
{
   SCIP_Real tmp;
   SCIP_Real xproj;
   SCIP_Real gval;
   SCIP_Real gderiv;
   int       iter;

   assert(scip != NULL);
   assert(!SCIPisFeasNegative(scip, xref+xoffset));
   assert(!SCIPisInfinity(scip, xref));

   if( xref < xmin )
      xref = xmin;

   xproj = xref;
   iter = 0;
   if( exponent == 2.0 )
      do
      {
         tmp = (xproj+xoffset) * (xproj+xoffset);
         gval = xproj - xref + 2*(xproj+xoffset) / zcoef * ((tmp-rhs)/zcoef + zref);
         if( !SCIPisFeasPositive(scip, ABS(gval)) )
            break;

         gderiv = 1 + 6 * tmp / (zcoef*zcoef) + 2 / zcoef * (zref - rhs/zcoef);
         xproj -= gval / gderiv;

      }
      while( ++iter <= 5 );
   else
      do
      {
         tmp = pow(xproj + xoffset, exponent-1);
         gval = xproj - xref + exponent / zcoef * (pow(xproj+xoffset, 2*exponent-1)/zcoef + tmp * (zref-rhs/zcoef));
         if( !SCIPisFeasPositive(scip, ABS(gval)) )
            break;

         gderiv = 1 + exponent / zcoef * ( (2*exponent-1)*tmp*tmp/zcoef + (exponent-1)*pow(xproj+xoffset, exponent-2) * (zref-rhs/zcoef) );
         xproj -= gval / gderiv;

      }
      while( ++iter <= 5 );

   if( xproj < xmin )
      xproj = xmin;

   SCIP_CALL( generateLinearizationCut(scip, row, xproj, exponent, xoffset, xmult, zcoef, rhs, x, z, islocal) );

   return SCIP_OKAY;
}

/** computes secant underestimator for sign(x+offset)abs(x+offset)^n + c*z <= rhs
 *
 * the generated cut is slope*xmult*x + c*z <= rhs + (-xlb-offset)^n + slope*xlb,
 * where slope = (sign(xub+offset)*abs(xub+offset)^n + (-xlb-offset)^n) / (xub - xlb).
 *
 * the cut is not generated if the given solution (or the LP solution) would not be cutoff
 */
static
SCIP_RETCODE generateSecantCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store row */
   SCIP_SOL*             sol,                /**< point we want to cut off, or NULL for LP solution */
   SCIP_Real             xlb,                /**< lower bound of x */
   SCIP_Real             xub,                /**< upper bound of x */
   SCIP_Real             exponent,           /**< exponent n in sign(x+offset)abs(x+offset)^n */
   SCIP_Real             xoffset,            /**< offset of x */
   DECL_MYPOW            ((*mypow)),         /**< function to use for computing power */
   SCIP_Real             xmult,              /**< multiplier for coefficient of x */
   SCIP_Real             zcoef,              /**< coefficient of z */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_VAR*             x,                  /**< variable x */
   SCIP_VAR*             z                   /**< variable z */
   )
{
   SCIP_Real slope, tmp, val;

   assert(scip != NULL);
   assert(!SCIPisRelEQ(scip, xlb, xub));   /* constraints with fixed x should have been removed before separation is called */

   assert(xlb < xub);
   assert(!SCIPisPositive(scip, xlb+xoffset));

   if( xlb > -xoffset )
      xlb = -xoffset;

   tmp = mypow(-xlb-xoffset, exponent);
   slope  = SIGN(xub+xoffset) * mypow(ABS(xub+xoffset), exponent) + tmp;
   slope /= xub - xlb;

   /* check if cut would violated solution, check that slope is not above value of infinity */
   val = -tmp + slope * (xmult * SCIPgetSolVal(scip, sol, x) - xlb) + zcoef * SCIPgetSolVal(scip, sol, z) - rhs;
   if( !SCIPisFeasPositive(scip, val) || SCIPisInfinity(scip, REALABS(slope)) )
   {
      *row = NULL;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, "signpowsecantcut", -SCIPinfinity(scip), SCIPinfinity(scip),
         SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0 /* local */, FALSE /* modifiable */, TRUE /* removable */ ) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, x, xmult*slope) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, z, zcoef) );
   SCIP_CALL( SCIPchgRowRhs(scip, *row, rhs + tmp + slope*xlb) );

   return SCIP_OKAY;
}

/** computes secant underestimator for sign(x+xoffset)abs(x+xoffset)^n + c*z <= rhs
 *
 * the generated cut is slope*xmult*x + c*z <= rhs + (-xlb-xoffset)^n + slope*xlb,
 * where slope = (sign(xub+xoffset)*abs(xub+xoffset)^n + (-xlb-xoffset)^n) / (xub - xlb).
 */
static
SCIP_RETCODE generateSecantCutNoCheck(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< buffer to store row */
   SCIP_Real             xlb,                /**< lower bound of x */
   SCIP_Real             xub,                /**< upper bound of x */
   SCIP_Real             exponent,           /**< exponent n in sign(x)abs(x)^n */
   SCIP_Real             xoffset,            /**< offset of x */
   DECL_MYPOW            ((*mypow)),         /**< function to use for computing power */
   SCIP_Real             xmult,              /**< multiplier for coefficient of x */
   SCIP_Real             zcoef,              /**< coefficient of z */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_VAR*             x,                  /**< variable x */
   SCIP_VAR*             z                   /**< variable z */
   )
{
   SCIP_Real slope, tmp;

   assert(scip != NULL);
   assert(!SCIPisRelEQ(scip, xlb, xub));   /* constraints with fixed x should have been removed before separation is called */
   assert(xlb < xub);
   assert(!SCIPisPositive(scip, xlb + xoffset));

   if( xlb > -xoffset )
      xlb = -xoffset;

   tmp = mypow(-xlb-xoffset, exponent);
   slope  = SIGN(xub+xoffset) * mypow(ABS(xub+xoffset), exponent) + tmp;
   slope /= xub - xlb;

   SCIP_CALL( SCIPcreateEmptyRow(scip, row, "signpowcut", -SCIPinfinity(scip), SCIPinfinity(scip),
         SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0 /* local */, FALSE /* modifiable */, TRUE /* removable */ ) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, x, xmult*slope) );
   SCIP_CALL( SCIPaddVarToRow(scip, *row, z, zcoef) );
   SCIP_CALL( SCIPchgRowRhs(scip, *row, rhs + tmp + slope*xlb) );

   return SCIP_OKAY;
}

/** generates a cut
 * based on Liberti and Pantelides, Convex Envelopes of Monomials of Odd Degree, J. Global Optimization 25, 157-168, 2003, and previous publications
 */
static
SCIP_RETCODE generateCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_ROW**            row,                /**< storage for cut */
   SCIP_Bool             onlyinbounds        /**< whether linearization is allowed only in variable bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_SIDETYPE  violside;
   SCIP_Real      c;
   SCIP_Real      xlb;
   SCIP_Real      xglb;
   SCIP_Real      xub;
   SCIP_Real      xval;
   SCIP_Real      xoffset;
   SCIP_Real      xmult;
   SCIP_Real      zcoef;
   SCIP_Real      rhs;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(row  != NULL);

   conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   assert(SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)));

   violside = SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT;
   *row = NULL;

   SCIPdebugMessage("generate cut for constraint <%s> with violated side %d\n", SCIPconsGetName(cons), violside);
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );
   SCIPdebugMessage("xlb = %g  xub = %g  xval = %g\n", SCIPvarGetLbLocal(consdata->x), SCIPvarGetUbLocal(consdata->x), SCIPgetSolVal(scip, sol, consdata->x));

   if( violside == SCIP_SIDETYPE_RIGHT )
   {
      xglb  = SCIPvarGetLbGlobal(consdata->x);
      xlb   = SCIPvarGetLbLocal(consdata->x);
      xub   = SCIPvarGetUbLocal(consdata->x);
      xval  = SCIPgetSolVal(scip, sol, consdata->x);
      xoffset = consdata->xoffset;
      xmult = 1.0;
      zcoef = consdata->zcoef;
      rhs   = consdata->rhs;
   }
   else
   {
      xglb  = -SCIPvarGetUbGlobal(consdata->x);
      xlb   = -SCIPvarGetUbLocal(consdata->x);
      xub   = -SCIPvarGetLbLocal(consdata->x);
      xval  = -SCIPgetSolVal(scip, sol, consdata->x);
      xoffset = -consdata->xoffset;
      xmult = -1.0;
      zcoef = -consdata->zcoef;
      rhs   = -consdata->lhs;
   }

   if( SCIPisInfinity(scip, REALABS(xval)) )
   {
      SCIPdebugMessage("skip separation since x is at infinity\n");
      return SCIP_OKAY;
   }

   if( !SCIPisNegative(scip, xlb+xoffset) )
   {
      /* [xlb, xub] completely in positive orthant -> function is convex on whole domain */
      SCIP_Bool islocal;

      islocal = (!SCIPconsIsGlobal(cons) || SCIPisNegative(scip, xglb+xoffset)) && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0;
      if( conshdlrdata->projectrefpoint && !onlyinbounds )
         SCIP_CALL( generateLinearizationCutProject(scip, row, xval, SCIPgetSolVal(scip, sol, consdata->z), -xoffset, consdata->exponent, xoffset, xmult, zcoef, rhs, consdata->x, consdata->z, islocal) );
      else if( !onlyinbounds )
         SCIP_CALL( generateLinearizationCut(scip, row, xval, consdata->exponent, xoffset, xmult, zcoef, rhs, consdata->x, consdata->z, islocal) );
      else
         SCIP_CALL( generateLinearizationCut(scip, row, 2.0*xval > xlb + xub ? xub : xlb, consdata->exponent, xoffset, xmult, zcoef, rhs, consdata->x, consdata->z, islocal) );
   }
   else if( !SCIPisPositive(scip, xub+xoffset) )
   {
      /* [xlb, xub] completely in negative orthant -> function is concave on whole domain */
      if( SCIPisInfinity(scip, -xlb) )
         return SCIP_OKAY;
      SCIP_CALL( generateSecantCut(scip, row, sol, xlb, xub, consdata->exponent, xoffset, consdata->pow, xmult, zcoef, rhs, consdata->x, consdata->z) );
   }
   else if( (c = - consdata->root * (xlb+xoffset) - xoffset) > xub )
   {
      /* c is right of xub -> use secant */
      if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, xub) )
         return SCIP_OKAY;
      SCIP_CALL( generateSecantCut(scip, row, sol, xlb, xub, consdata->exponent, xoffset, consdata->pow, xmult, zcoef, rhs, consdata->x, consdata->z) );
   }
   else if( xval >= c )
   {
      /* xval is right of c -> use linearization */
      if( conshdlrdata->projectrefpoint && !onlyinbounds )
         SCIP_CALL( generateLinearizationCutProject(scip, row, xval, SCIPgetSolVal(scip, sol, consdata->z), c, consdata->exponent, xoffset, xmult, zcoef, rhs, consdata->x, consdata->z, SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0) );
      else if( !onlyinbounds )
         SCIP_CALL( generateLinearizationCut(scip, row, xval, consdata->exponent, xoffset, xmult, zcoef, rhs, consdata->x, consdata->z, xval+xoffset < - consdata->root * (xglb+xoffset) && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0) );
      else
         SCIP_CALL( generateLinearizationCut(scip, row, xub, consdata->exponent, xoffset, xmult, zcoef, rhs, consdata->x, consdata->z, xval+xoffset < - consdata->root * (xglb+xoffset) && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) > 0) );
   }
   else
   {
      /* xval between xlb and c -> use secant */
      if( SCIPisInfinity(scip, -xlb) || SCIPisInfinity(scip, c) )
         return SCIP_OKAY;
      SCIP_CALL( generateSecantCut(scip, row, sol, xlb, c, consdata->exponent, xoffset, consdata->pow, xmult, zcoef, rhs, consdata->x, consdata->z) );
   }

   /* check numerics */
   if( *row != NULL )
   {
      SCIPdebug( SCIPprintRow(scip, *row, NULL) );

      /* check range of coefficients */
      SCIPdebugMessage(" -> found cut rhs=%f, min=%f, max=%f range=%g\n",
         SCIProwGetRhs(*row),
         SCIPgetRowMinCoef(scip, *row), SCIPgetRowMaxCoef(scip, *row),
         SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));

      if( SCIPisInfinity(scip, REALABS(SCIProwGetRhs(*row))) )
      {
         SCIPdebugMessage("skip cut for constraint <%s> because of very large right hand side: %g\n", SCIPconsGetName(cons), SCIProwGetRhs(*row));
         SCIP_CALL( SCIPreleaseRow(scip, row) );
         return SCIP_OKAY;
      }

      if( SCIPgetRowMaxCoef(scip, *row) / SCIPgetRowMinCoef(scip, *row) >= conshdlrdata->cutmaxrange )
      {
         SCIPdebugMessage("skip cut for constraint <%s> because of very large range: %g\n", SCIPconsGetName(cons), SCIPgetRowMaxCoef(scip, *row)/SCIPgetRowMinCoef(scip, *row));
         SCIP_CALL( SCIPreleaseRow(scip, row) );
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 *  assumes that constraint violations have been computed
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< quadratic constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Real             minefficacy,        /**< minimal efficacy of a cut if it should be added to the LP */
   SCIP_Bool             convexalways,       /**< whether to ignore minefficacy criteria for a convex constraint (and use feastol instead) */
   SCIP_Bool             onlyinbounds,       /**< whether linearization is allowed only in variable bounds */
   SCIP_Bool*            success,            /**< result of separation: separated point (TRUE) or not (FALSE) */
   SCIP_Real*            bestefficacy        /**< buffer to store best efficacy of a cut that was added to the LP, if found; or NULL if not of interest */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          efficacy;
   SCIP_Real          feasibility;
   SCIP_Real          norm;
   SCIP_Bool          convex;
   int                c;
   SCIP_ROW*          row;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(success  != NULL);

   *success = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( bestefficacy != NULL )
      *bestefficacy = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         /* try to generate a cut */
         SCIP_CALL( generateCut(scip, conss[c], sol, &row, onlyinbounds) );  /*lint !e613*/
         if( row == NULL ) /* failed to generate cut */
            continue;

         /* check if we separate in convex area */
         if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         {
            convex = !SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->x))
               && (!SCIPisNegative(scip, SCIPvarGetLbLocal(consdata->x)+consdata->xoffset)
                  || SCIPgetSolVal(scip, NULL, consdata->x)+consdata->xoffset >= -consdata->root*(SCIPvarGetLbLocal(consdata->x)+consdata->xoffset));
         }
         else
         {
            convex = !SCIPisInfinity(scip,  SCIPvarGetUbLocal(consdata->x))
               && (!SCIPisPositive(scip, SCIPvarGetUbLocal(consdata->x)+consdata->xoffset)
                  || SCIPgetSolVal(scip, NULL, consdata->x)+consdata->xoffset <= -consdata->root*(SCIPvarGetUbLocal(consdata->x)+consdata->xoffset));
         }

         feasibility = SCIPgetRowSolFeasibility(scip, row, sol);

         /* in difference to SCIPgetCutEfficacy, we scale by norm only if the norm is > 1.0
          * this avoid finding cuts efficiant which are only very slightly violated
          * CPLEX does not seem to scale row coefficients up too
          */
         norm = SCIPgetRowMaxCoef(scip, row);
         if( norm > 1.0 )
            efficacy = -feasibility / norm;
         else
            efficacy = -feasibility;

         /* if cut is strong or it's weak but we are convex and desperate, then add */
         if( efficacy > minefficacy || (convexalways && convex && SCIPisFeasPositive(scip, efficacy)) )
         {
            SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE) );
            *success = TRUE;
            if( bestefficacy != NULL && efficacy > *bestefficacy )
               *bestefficacy = efficacy;

            /* notify indicator constraint handler about this cut */
            if( conshdlrdata->conshdlrindicator != NULL && !SCIProwIsLocal(row) )
            {
               SCIP_CALL( SCIPaddRowIndicator(scip, conshdlrdata->conshdlrindicator, row) );
            }
         }

         SCIP_CALL( SCIPreleaseRow (scip, &row) );
      }

      /* enforce only useful constraints
       * others are only checked and enforced if we are still feasible or have not found a separating cut yet
       */
      if( c >= nusefulconss && *success )
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
   SCIP_CONSDATA* consdata;
   SCIP_ROW* row;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   if( separatedlpsol != NULL )
      *separatedlpsol = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      if( SCIPconsIsLocal(conss[c]) )  /*lint !e613*/
         continue;

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( !SCIPisGT(scip, SCIPvarGetUbGlobal(consdata->x), -consdata->xoffset) && !SCIPisInfinity(scip, -consdata->lhs) )
      {
         /* constraint function is concave for x+offset <= 0.0, so can linearize w.r.t. lhs */
         consdata->lhsviol = 1.0;
         consdata->rhsviol = 0.0;
         SCIP_CALL( generateCut(scip, conss[c], ref, &row, FALSE) );  /*lint !e613*/
      }
      else if( !SCIPisLT(scip, SCIPvarGetLbGlobal(consdata->x), -consdata->xoffset) && !SCIPisInfinity(scip, -consdata->rhs) )
      {
         /* constraint function is convex for x+offset >= 0.0, so can linearize w.r.t. rhs */
         consdata->lhsviol = 0.0;
         consdata->rhsviol = 1.0;
         SCIP_CALL( generateCut(scip, conss[c], ref, &row, FALSE) );  /*lint !e613*/
      }
      else
      {
         /* sign not fixed or nothing to linearize */
         continue;
      }

      if( row == NULL )
         continue;

      assert(!SCIProwIsLocal(row));

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

      SCIP_CALL( SCIPaddPoolCut(scip, row) );

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

/** given a solution, try to make absolute power constraints feasible by shifting the linear variable z and pass this solution to the trysol heuristic */
static
SCIP_RETCODE proposeFeasibleSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol                 /**< solution to process */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_SOL* newsol;
   SCIP_Real xtermval;
   SCIP_Real zval;
   SCIP_Real viol;
   int c;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

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
      assert(consdata->z != NULL);
      assert(!SCIPisZero(scip, consdata->zcoef));

      /* recompute violation w.r.t. current solution */
      SCIP_CALL( computeViolation(scip, conss[c], newsol, &viol) );  /*lint !e613*/

      /* do nothing if constraint is satisfied */
      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      /* @todo could also adjust x while keeping z fixed */

      /* if variable is multiaggregated, then cannot set its solution value, so give up */
      if( SCIPvarGetStatus(consdata->z) == SCIP_VARSTATUS_MULTAGGR )
         break;

      /* compute value of x-term */
      xtermval  = SCIPgetSolVal(scip, newsol, consdata->x);
      xtermval += consdata->xoffset;
      xtermval  = SIGN(xtermval) * consdata->pow(ABS(xtermval), consdata->exponent);

      /* if left hand side is violated, try to set z such that lhs is active */
      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
      {
         assert(!SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip))); /* should only have one side violated (otherwise some variable is at infinity) */

         zval = (consdata->lhs - xtermval)/consdata->zcoef;
         /* bad luck: z would get value outside of its domain */
         if( SCIPisFeasLT(scip, zval, SCIPvarGetLbGlobal(consdata->z)) || SCIPisFeasGT(scip, zval, SCIPvarGetUbGlobal(consdata->z)) )
            break;
         SCIP_CALL( SCIPsetSolVal(scip, newsol, consdata->z, zval) );
      }

      /* if right hand side is violated, try to set z such that rhs is active */
      if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         zval = (consdata->rhs - xtermval)/consdata->zcoef;
         /* bad luck: z would get value outside of its domain */
         if( SCIPisFeasLT(scip, zval, SCIPvarGetLbGlobal(consdata->z)) || SCIPisFeasGT(scip, zval, SCIPvarGetUbGlobal(consdata->z)) )
            break;
         SCIP_CALL( SCIPsetSolVal(scip, newsol, consdata->z, zval) );
      }
   }

   /* if we have a solution that should satisfy all absolute power constraints and has a better objective than the current upper bound, then pass it to the trysol heuristic */
   if( c == nconss )
   {
      SCIP_CONSHDLRDATA* conshdlrdata;

      SCIPdebugMessage("pass solution with objective val %g to trysol heuristic\n", SCIPgetSolTransObj(scip, newsol));

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);
      assert(conshdlrdata->trysolheur != NULL);

      SCIP_CALL( SCIPheurPassSolTrySol(scip, conshdlrdata->trysolheur, newsol) );
   }

   SCIP_CALL( SCIPfreeSol(scip, &newsol) );

   return SCIP_OKAY;
}

/** create a nonlinear row representation of the constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_EXPRTREE* exprtree;
   SCIP_QUADELEM quadelem;
   SCIP_VAR* linvars[2];
   SCIP_Real lincoefs[2];
   SCIP_VAR* quadvar;
   SCIP_Real constant;
   SCIP_Bool expisint;
   int sign;
   int nlinvars;
   int nquadvars;
   int nquadelems;
   int n;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   nlinvars = 0;
   nquadvars = 0;
   nquadelems = 0;
   exprtree = NULL;
   constant = 0.0;

   /* check if sign of x is fixed */
   if( !SCIPisNegative(scip, SCIPvarGetLbGlobal(consdata->x)+consdata->xoffset) )
      sign =  1;
   else if( !SCIPisPositive(scip, SCIPvarGetUbGlobal(consdata->x)+consdata->xoffset) )
      sign = -1;
   else
      sign =  0;

   /* check if exponent is integral */
   expisint = SCIPisIntegral(scip, consdata->exponent);
   n = (int)SCIPround(scip, consdata->exponent);

   /* create quadelem or expression tree for nonlinear part sign(x+offset)abs(x+offset)^n */
   if( sign != 0 || (expisint && (n % 2 == 1)) )
   {
      /* sign is fixes or exponent is odd integer */
      if( expisint && n == 2 )
      {
         /* sign of x is clear and exponent is 2.0 -> generate quadratic, linear, and constant term for +/- (x+offset)^n */
         assert(sign == -1 || sign == 1);
         nquadelems = 1;
         quadelem.idx1 = 0;
         quadelem.idx2 = 0;
         quadelem.coef = (SCIP_Real)sign;
         nquadvars = 1;
         quadvar = consdata->x;

         if( consdata->xoffset != 0.0 )
         {
            linvars[0] = consdata->x;
            lincoefs[0] = sign * 2.0 * consdata->xoffset;
            constant = sign * consdata->xoffset * consdata->xoffset;
         }
      }
      else
      {
         /* exponent is odd or sign of x is clear, generate expression tree for +/- (+/-(x+offset))^exponent */
         SCIP_EXPR* expr;
         SCIP_EXPR* expr2;

         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) ); /* x */
         if( consdata->xoffset != 0.0 )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->xoffset) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PLUS, expr, expr2) ); /* x + offset */
         }
         if( sign == -1 && !expisint )
         {
            /* if exponent is not integer and x is negative, then negate */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, -1.0) );         /* -1 */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL,   expr, expr2) );  /* -(x+offset) */
         }
         /* use intpower for integer exponent and realpower for fractional exponent */
         if( expisint )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_INTPOWER, expr, n) );  /* (x+offset)^n */
         }
         else
         {
            assert(sign == 1 || sign == -1);
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_REALPOWER, expr, consdata->exponent) );  /* abs(x+offset)^exponent */
         }
         /* if exponent is even integer, then negate result; if it's an odd integer, then intpower already takes care of correct sign */
         if( sign == -1 && !(expisint && n % 2 == 1) )
         {
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, -1.0) );         /* -1 */
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr,  SCIP_EXPR_MUL,   expr, expr2) );  /* -abs(x+offset)^exponent */
         }
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 1, 0, NULL) );
      }
   }
   else
   {
      /* exponent is not odd integer and sign of x is not fixed -> generate expression tree for signpower(x+offset, n) */
      SCIP_EXPR* expr;
      SCIP_EXPR* expr2;

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_VARIDX, 0) ); /* x */
      if( consdata->xoffset != 0.0 )
      {
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr2, SCIP_EXPR_CONST, consdata->xoffset) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_PLUS, expr, expr2) );   /* x + offset */
      }
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_SIGNPOWER, expr, (SCIP_Real)consdata->exponent) ); /* signpower(x+offset, n) */

      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 1, 0, NULL) );
   }
   assert(exprtree != NULL || nquadelems > 0);

   /* tell expression tree, which is its variable */
   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, 1, &consdata->x) );
   }

   assert(nlinvars < 2);
   linvars[nlinvars] = consdata->z;
   lincoefs[nlinvars] = consdata->zcoef;
   ++nlinvars;

   /* create nlrow */
   SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), constant,
         nlinvars, linvars, lincoefs,
         nquadvars, &quadvar, nquadelems, &quadelem,
         exprtree, consdata->lhs, consdata->rhs
         ) );

   if( exprtree != NULL )
   {
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** upgrades a quadratic constraint where the quadratic part is only a single square term and the quadratic variable sign is fixed to a signpower constraint */
static
SCIP_DECL_QUADCONSUPGD(quadconsUpgdAbspower)
{  /*lint --e{715}*/
   SCIP_QUADVARTERM quadvarterm;
   SCIP_VAR* x;
   SCIP_VAR* z;
   SCIP_Real xoffset;
   SCIP_Real zcoef;
   SCIP_Real signpowcoef;
   SCIP_Real lhs;
   SCIP_Real rhs;

   *nupgdconss = 0;

   /* need at least one linear variable */
   if( SCIPgetNLinearVarsQuadratic(scip, cons) == 0 )
      return SCIP_OKAY;

   /* consider only quadratic constraints with a single square term */
   if( SCIPgetNQuadVarTermsQuadratic(scip, cons) != 1 )
      return SCIP_OKAY;
   assert(SCIPgetNBilinTermsQuadratic(scip, cons) == 0);

   quadvarterm = SCIPgetQuadVarTermsQuadratic(scip, cons)[0];
   if( SCIPisZero(scip, quadvarterm.sqrcoef) )
      return SCIP_OKAY;

   x = quadvarterm.var;
   xoffset = quadvarterm.lincoef / (2.0 * quadvarterm.sqrcoef);

   /* check that x has fixed sign */
   if( SCIPisNegative(scip, SCIPvarGetLbGlobal(x) + xoffset) && SCIPisPositive(scip, SCIPvarGetUbGlobal(x) + xoffset) )
      return SCIP_OKAY;

   /* check whether upgdconss array has enough space to store 1 or 2 constraints */
   if( SCIPgetNLinearVarsQuadratic(scip, cons) > 1 )
      *nupgdconss = -2;
   else
      *nupgdconss = -1;
   if( -*nupgdconss > upgdconsssize )
      return SCIP_OKAY;

   *nupgdconss = 0;

   SCIPdebugMessage("upgrade quadratic constraint <%s> to absolute power, x = [%g,%g], offset = %g\n", SCIPconsGetName(cons), SCIPvarGetLbGlobal(x), SCIPvarGetUbGlobal(x), xoffset);
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons, NULL) ) );

   lhs = SCIPgetLhsQuadratic(scip, cons);
   rhs = SCIPgetRhsQuadratic(scip, cons);

   /* get z and its coefficient */
   if( SCIPgetNLinearVarsQuadratic(scip, cons) > 1 )
   {
      /* create auxiliary variable and constraint for linear part, since we can handle only at most one variable in cons_signpower */
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR* auxvar;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_linpart", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
            SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );

      SCIP_CALL( SCIPcreateConsLinear(scip, &upgdconss[0], name, SCIPgetNLinearVarsQuadratic(scip, cons),
            SCIPgetLinearVarsQuadratic(scip, cons), SCIPgetCoefsLinearVarsQuadratic(scip, cons),
            SCIPisInfinity(scip, -lhs) ? -SCIPinfinity(scip) : 0.0,
            SCIPisInfinity(scip,  rhs) ?  SCIPinfinity(scip) : 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCoefLinear(scip, upgdconss[*nupgdconss], auxvar, -1.0) );

      z = auxvar;
      zcoef = 1.0;

      ++*nupgdconss;

      /* compute and set value of auxvar in debug solution */
#ifdef SCIP_DEBUG_SOLUTION
      {
         SCIP_Real debugval;
         SCIP_Real debugvarval;
         int i;

         debugval = 0.0;
         for( i = 0; i < SCIPgetNLinearVarsQuadratic(scip, cons); ++i )
         {
            SCIP_CALL( SCIPdebugGetSolVal(scip, SCIPgetLinearVarsQuadratic(scip, cons)[i], &debugvarval) );
            debugval += SCIPgetCoefsLinearVarsQuadratic(scip, cons)[i] * debugvarval;
         }

         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
      }
#endif

      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
   }
   else
   {
      assert(SCIPgetNLinearVarsQuadratic(scip, cons) == 1);
      z = SCIPgetLinearVarsQuadratic(scip, cons)[0];
      zcoef = SCIPgetCoefsLinearVarsQuadratic(scip, cons)[0];
   }

   /* we now have lhs <= sqrcoef * (x + offset)^2 - sqrcoef * offset^2 + zcoef * z <= rhs */

   /* move sqrcoef * offset^2 into lhs and rhs */
   if( !SCIPisInfinity(scip, -lhs) )
      lhs += quadvarterm.sqrcoef * xoffset * xoffset;
   if( !SCIPisInfinity(scip,  rhs) )
      rhs += quadvarterm.sqrcoef * xoffset * xoffset;

   /* divide by sqrcoef if x+offset > 0 and by -sqrcoef if < 0 */
   signpowcoef = quadvarterm.sqrcoef;
   if( SCIPisNegative(scip, SCIPvarGetLbGlobal(x) + xoffset) )
      signpowcoef = -signpowcoef;
   if( signpowcoef > 0.0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
         lhs /= signpowcoef;
      if( !SCIPisInfinity(scip,  rhs) )
         rhs /= signpowcoef;
   }
   else
   {
      SCIP_Real newrhs;

      if( !SCIPisInfinity(scip, -lhs) )
         newrhs = lhs / signpowcoef;
      else
         newrhs = SCIPinfinity(scip);
      if( !SCIPisInfinity(scip,  rhs) )
         lhs = rhs / signpowcoef;
      else
         lhs = -SCIPinfinity(scip);
      rhs = newrhs;
   }
   zcoef /= signpowcoef;

   /* create the absolute power constraint */
   SCIP_CALL( SCIPcreateConsAbspower(scip, &upgdconss[*nupgdconss], SCIPconsGetName(cons), x, z, 2.0,
         xoffset, zcoef, lhs, rhs,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );
   SCIPdebug( SCIP_CALL( SCIPprintCons(scip, upgdconss[*nupgdconss], NULL) ) );
   ++*nupgdconss;

   return SCIP_OKAY;
}

/** tries to upgrade a nonlinear constraint into a absolute power constraint */
static
SCIP_DECL_NONLINCONSUPGD(nonlinconsUpgdAbspower)
{
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   SCIP_EXPRGRAPHNODE* child;
   SCIP_Real exponent;
   SCIP_VAR* x;
   SCIP_VAR* z;
   SCIP_Real signpowcoef;
   SCIP_Real zcoef;
   SCIP_Real xoffset;
   SCIP_Real constant;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   /* absolute power needs at least one linear variable (constraint is trivial, otherwise) */
   if( SCIPgetNLinearVarsNonlinear(scip, cons) == 0 )
      return SCIP_OKAY;

   node = SCIPgetExprgraphNodeNonlinear(scip, cons);

   /* no interest in linear constraints */
   if( node == NULL )
      return SCIP_OKAY;

   /* need exactly one argument */
   if( SCIPexprgraphGetNodeNChildren(node) != 1 )
      return SCIP_OKAY;

   constant = 0.0;
   signpowcoef = 1.0; /* coefficient of sign(x)abs(x)^n term, to be reformulated away... */

   child = SCIPexprgraphGetNodeChildren(node)[0];

   /* check if node expression fits to absolute power constraint */
   switch( SCIPexprgraphGetNodeOperator(node) )
   {
   case SCIP_EXPR_REALPOWER:
   {
      /* realpower with exponent > 1.0 can always be signpower, since it assumes that argument is >= 0.0 */
      exponent = SCIPexprgraphGetNodeRealPowerExponent(node);
      if( exponent <= 1.0 )
         return SCIP_OKAY;

      assert(SCIPexprgraphGetNodeBounds(child).inf >= 0.0);
      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      /* check if exponent > 1.0 and either odd or even with child having fixed sign */
      SCIP_INTERVAL childbounds;

      exponent = (SCIP_Real)SCIPexprgraphGetNodeIntPowerExponent(node);
      if( exponent <= 1.0 )
         return SCIP_OKAY;

      childbounds = SCIPexprgraphGetNodeBounds(child);
      if( (int)exponent % 2 == 0 && childbounds.inf < 0.0 && childbounds.sup > 0.0 )
         return SCIP_OKAY;

      /* use x^exponent = -sign(x) |x|^exponent if exponent is even and x always negative */
      if( (int)exponent % 2 == 0 && childbounds.inf < 0.0 )
         signpowcoef = -1.0;

      break;
   }

   case SCIP_EXPR_SQUARE:
   {
      /* check if child has fixed sign */
      SCIP_INTERVAL childbounds;

      childbounds = SCIPexprgraphGetNodeBounds(child);
      if( childbounds.inf < 0.0 && childbounds.sup > 0.0 )
         return SCIP_OKAY;

      /* use x^2 = -sign(x) |x|^2 if x is always negative */
      if( childbounds.inf < 0.0 )
         signpowcoef = -1.0;

      exponent = 2.0;
      break;
   }

   case SCIP_EXPR_SIGNPOWER:
   {
      /* check if exponent > 1.0 */
      exponent = SCIPexprgraphGetNodeSignPowerExponent(node);
      if( exponent <= 1.0 )
         return SCIP_OKAY;
      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_INTERVAL childbounds;

      /* check if only one univariate monomial with exponent > 1.0 */

      /* if sum of univariate monomials, then this should have been taken care of by exprgraphnodeReformSignpower */
      if( SCIPexprgraphGetNodePolynomialNMonomials(node) > 1 )
         return SCIP_OKAY;
      assert(SCIPexprgraphGetNodePolynomialNMonomials(node) == 1); /* assume simplified, i.e., no constant polynomial */

      monomial = SCIPexprgraphGetNodePolynomialMonomials(node)[0];
      assert(SCIPexprGetMonomialNFactors(monomial) == 1); /* since we have only one children and assume simplified */

      exponent = SCIPexprGetMonomialExponents(monomial)[0];
      if( exponent <= 1.0 )
         return SCIP_OKAY;

      /* if exponent is even integer and child has mixed sign, then cannot do
       * if exponent is even integer and child is always negative, then can do via multiplication by -1.0 */
      childbounds = SCIPexprgraphGetNodeBounds(child);
      if( SCIPisIntegral(scip, exponent) && ((int)SCIPround(scip, exponent) % 2 == 0) && childbounds.inf < 0.0 )
      {
         if( childbounds.sup > 0.0 )
            return SCIP_OKAY;
         signpowcoef = -1.0;
      }

      constant = SCIPexprgraphGetNodePolynomialConstant(node);
      signpowcoef *= SCIPexprGetMonomialCoef(monomial);

      break;
   }

   default:
      return SCIP_OKAY;
   }  /*lint !e788*/
   assert(SCIPexprgraphGetNodeNChildren(node) == 1);

   /* check magnitue of coefficient of z in signpower constraint */
   zcoef = 1.0;
   if( SCIPgetNLinearVarsNonlinear(scip, cons) == 1 )
      zcoef = SCIPgetLinearCoefsNonlinear(scip, cons)[0];
   zcoef /= signpowcoef;
   if( SCIPexprgraphGetNodeOperator(child) == SCIP_EXPR_LINEAR && SCIPexprgraphGetNodeNChildren(child) == 1 )
   {
      zcoef /= pow(REALABS(SCIPexprgraphGetNodeLinearCoefs(child)[0]), exponent);
   }
   if( SCIPisZero(scip, zcoef) )
   {
      SCIPdebugMessage("skip upgrade to signpower since |zcoef| = %g would be zero\n", zcoef);
      return SCIP_OKAY;
   }

   /* count how many constraints we need to add (use negative numbers, for convenience):
    * one constraint for absolute power,
    * plus one if we need to replace the linear part by single variable,
    * plus one if we need to replace the argument of absolute power by a single variable
    */
   *nupgdconss = -1;

   if( SCIPexprgraphGetNodeOperator(child) != SCIP_EXPR_VARIDX && (SCIPexprgraphGetNodeOperator(child) != SCIP_EXPR_LINEAR || SCIPexprgraphGetNodeNChildren(child) > 1) )
   {
      /* if node has known curvature and we would add auxiliary var for child, then don't upgrade
       * it's not really necessary, but may introduce more numerical troubles
       * @todo maybe still do if child is linear?
       */
      if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
      {
         *nupgdconss = 0;
         return SCIP_OKAY;
      }

      --*nupgdconss;
   }

   if( SCIPgetNLinearVarsNonlinear(scip, cons) > 1 )
      --*nupgdconss;

   /* request larger upgdconss array */
   if( upgdconsssize < -*nupgdconss )
      return SCIP_OKAY;

   SCIPdebugMessage("upgrading constraint <%s>\n", SCIPconsGetName(cons));

   /* start counting at zero again */
   *nupgdconss = 0;

   exprgraph = SCIPgetExprgraphNonlinear(scip, SCIPconsGetHdlr(cons));

   lhs = SCIPgetLhsNonlinear(scip, cons);
   rhs = SCIPgetRhsNonlinear(scip, cons);

   /* get x and it's offset */
   if( SCIPexprgraphGetNodeOperator(child) == SCIP_EXPR_VARIDX )
   {
      x = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, child);
      xoffset = 0.0;
   }
   else if( SCIPexprgraphGetNodeOperator(child) == SCIP_EXPR_LINEAR && SCIPexprgraphGetNodeNChildren(child) == 1 )
   {
      SCIP_Real xcoef;

      x = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(child)[0]);
      xcoef = SCIPexprgraphGetNodeLinearCoefs(child)[0];
      assert(!SCIPisZero(scip, xcoef));

      signpowcoef *= (xcoef < 0.0 ? -1.0 : 1.0) * pow(REALABS(xcoef), exponent);
      xoffset = SCIPexprgraphGetNodeLinearConstant(child) / xcoef;
   }
   else
   {
      /* reformulate by adding auxiliary variable and constraint for child */
      char name[SCIP_MAXSTRLEN];
      SCIP_INTERVAL bounds;
      SCIP_VAR* auxvar;
      SCIP_Real minusone;

      bounds = SCIPexprgraphGetNodeBounds(child);
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_powerarg", SCIPconsGetName(cons));

      SCIPdebugMessage("add auxiliary variable and constraint %s for node %p(%d,%d)\n", name, (void*)child, SCIPexprgraphGetNodeDepth(child), SCIPexprgraphGetNodePosition(child));

      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, SCIPintervalGetInf(bounds), SCIPintervalGetSup(bounds), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );

      /* create new constraint child == auxvar
       * since signpower is monotonic, we need only child <= auxvar or child >= auxvar, if not both sides are finite, and depending on signpowcoef
       * i.e., we need child - auxvar <= 0.0 if rhs is finite and signpowcoef > 0.0 or lhs is finite and signpowcoef < 0.0
       *   and we need 0.0 <= child - auxvar if lhs is finite and signpowcoef > 0.0 or rhs is finite and signpowcoef < 0.0
       */
      minusone = -1.0;
      assert(upgdconsssize > *nupgdconss);
      SCIP_CALL( SCIPcreateConsNonlinear2(scip, &upgdconss[*nupgdconss], name, 1, &auxvar, &minusone, child,
            ((signpowcoef > 0.0 && !SCIPisInfinity(scip, -lhs)) || (signpowcoef < 0.0 && !SCIPisInfinity(scip,  rhs))) ? 0.0 : -SCIPinfinity(scip),
            ((signpowcoef > 0.0 && !SCIPisInfinity(scip,  rhs)) || (signpowcoef < 0.0 && !SCIPisInfinity(scip, -lhs))) ? 0.0 :  SCIPinfinity(scip),
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      ++*nupgdconss;

      /* use auxvar to setup absolute power constraint */
      x = auxvar;
      xoffset = 0.0;

      /* compute and set value of auxvar in debug solution, if debugging is enabled */
      SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(child)) );  /*lint !e506 !e774*/

      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
   }

   /* get z and its coefficient */
   if( SCIPgetNLinearVarsNonlinear(scip, cons) > 1 )
   {
      /* create auxiliary variable and constraint for linear part, since we can handle only at most one variable in cons_signpower */
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR* auxvar;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_linpart", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
            SCIPconsIsInitial(cons), SCIPconsIsRemovable(cons), NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );

      assert(upgdconsssize > *nupgdconss);
      SCIP_CALL( SCIPcreateConsLinear(scip, &upgdconss[*nupgdconss], name, SCIPgetNLinearVarsNonlinear(scip, cons),
            SCIPgetLinearVarsNonlinear(scip, cons), SCIPgetLinearCoefsNonlinear(scip, cons),
            SCIPisInfinity(scip, -lhs) ? -SCIPinfinity(scip) : 0.0,
            SCIPisInfinity(scip,  rhs) ?  SCIPinfinity(scip) : 0.0,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
            SCIPconsIsStickingAtNode(cons)) );
      SCIP_CALL( SCIPaddCoefLinear(scip, upgdconss[*nupgdconss], auxvar, -1.0) );

      z = auxvar;
      zcoef = 1.0;

      ++*nupgdconss;

      /* compute and set value of auxvar in debug solution */
#ifdef SCIP_DEBUG_SOLUTION
      {
         SCIP_Real debugval;
         SCIP_Real debugvarval;
         int i;

         debugval = 0.0;
         for( i = 0; i < SCIPgetNLinearVarsNonlinear(scip, cons); ++i )
         {
            SCIP_CALL( SCIPdebugGetSolVal(scip, SCIPgetLinearVarsNonlinear(scip, cons)[i], &debugvarval) );
            debugval += SCIPgetLinearCoefsNonlinear(scip, cons)[i] * debugvarval;
         }

         SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, debugval) );
      }
#endif

      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
   }
   else
   {
      assert(SCIPgetNLinearVarsNonlinear(scip, cons) == 1);
      z = SCIPgetLinearVarsNonlinear(scip, cons)[0];
      zcoef = SCIPgetLinearCoefsNonlinear(scip, cons)[0];
   }

   if( constant != 0.0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
         lhs -= constant;
      if( !SCIPisInfinity(scip,  rhs) )
         rhs -= constant;
   }

   /* divide absolute power constraint by signpowcoef */
   if( signpowcoef != 1.0 )
   {
      zcoef /= signpowcoef;
      if( signpowcoef < 0.0 )
      {
         SCIP_Real newrhs;

         newrhs = SCIPisInfinity(scip, -lhs) ? SCIPinfinity(scip) : lhs/signpowcoef;
         lhs = SCIPisInfinity(scip, rhs) ? -SCIPinfinity(scip) : rhs/signpowcoef;
         rhs = newrhs;
      }
      else
      {
         if( !SCIPisInfinity(scip, -lhs) )
            lhs /= signpowcoef;
         if( !SCIPisInfinity(scip,  rhs) )
            rhs /= signpowcoef;
      }
   }

   /* finally setup a absolute power constraint */

   assert(*nupgdconss < upgdconsssize);
   SCIP_CALL( SCIPcreateConsAbspower(scip, &upgdconss[*nupgdconss], SCIPconsGetName(cons),
         x, z, exponent, xoffset, zcoef, lhs, rhs,
         SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
         SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
         SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
         SCIPconsIsStickingAtNode(cons)) );
   ++*nupgdconss;

   return SCIP_OKAY;
}

/** tries to reformulate a expression graph node via introducing a absolute power constraint
 * if node fits to absolute power and has indefinte curvature and has no nonlinear parents and has siblings, then replace by auxvar and absolute power constraint
 * if it still has nonlinear parents, then we wait to see if reformulation code move node into auxiliary constraint,
 *   so we do not add unnessary auxiliary variables for something like an x^2 in an exp(x^2)
 * if it has no siblings, then we let the upgrading for nonlinear constraints take care of it,
 *   since it may be able to upgrade the constraint as a whole and can take the constraint sides into account too (may need only <=/>= auxcons)
 */
static
SCIP_DECL_EXPRGRAPHNODEREFORM(exprgraphnodeReformAbspower)
{
   SCIP_EXPRGRAPHNODE* child;
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   SCIP_Real exponent;
   SCIP_VAR* x;
   SCIP_VAR* z;
   SCIP_Real signpowcoef;
   SCIP_Real xoffset;
   SCIP_Real constant;

   assert(scip != NULL);
   assert(exprgraph != NULL);
   assert(node != NULL);
   assert(naddcons != NULL);
   assert(reformnode != NULL);

   *reformnode = NULL;

   if( SCIPexprgraphGetNodeCurvature(node) != SCIP_EXPRCURV_UNKNOWN )
      return SCIP_OKAY;

   constant = 0.0;
   signpowcoef = 1.0; /* coefficient of sign(x)abs(x)^n term, to be move in from of z... */

   /* check if node expression fits to absolute power constraint */
   switch( SCIPexprgraphGetNodeOperator(node) )
   {
   case SCIP_EXPR_REALPOWER:
   {
      /* realpower with exponent > 1.0 can always be absolute power, since it assumes that argument is >= 0.0 */
      exponent = SCIPexprgraphGetNodeRealPowerExponent(node);
      if( exponent <= 1.0 )
         return SCIP_OKAY;

      assert(SCIPexprgraphGetNodeBounds(SCIPexprgraphGetNodeChildren(node)[0]).inf >= 0.0);
      break;
   }

   case SCIP_EXPR_INTPOWER:
   {
      /* check if exponent > 1.0 and either odd or even with child having fixed sign */
      SCIP_INTERVAL childbounds;

      exponent = (SCIP_Real)SCIPexprgraphGetNodeIntPowerExponent(node);
      if( exponent <= 1.0 )
         return SCIP_OKAY;

      childbounds = SCIPexprgraphGetNodeBounds(SCIPexprgraphGetNodeChildren(node)[0]);
      if( (int)exponent % 2 == 0 && childbounds.inf < 0.0 && childbounds.sup > 0.0 )
         return SCIP_OKAY;

      /* use x^exponent = -sign(x) |x|^exponent if exponent is even and x always negative */
      if( (int)exponent % 2 == 0 && childbounds.inf < 0.0 )
         signpowcoef = -1.0;

      break;
   }

   case SCIP_EXPR_SQUARE:
   {
      /* check if child has fixed sign */
      SCIP_INTERVAL childbounds;

      childbounds = SCIPexprgraphGetNodeBounds(SCIPexprgraphGetNodeChildren(node)[0]);
      if( childbounds.inf < 0.0 && childbounds.sup > 0.0 )
         return SCIP_OKAY;

      /* use x^2 = -sign(x) |x|^2 if x is always negative */
      if( childbounds.inf < 0.0 )
         signpowcoef = -1.0;

      exponent = 2.0;
      break;
   }

   case SCIP_EXPR_SIGNPOWER:
   {
      /* check if exponent > 1.0 */
      exponent = SCIPexprgraphGetNodeSignPowerExponent(node);
      if( exponent <= 1.0 )
         return SCIP_OKAY;
      break;
   }

   case SCIP_EXPR_POLYNOMIAL:
   {
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_INTERVAL childbounds;

      /* check if only one univariate monomial with exponent > 1.0 */
      if( SCIPexprgraphGetNodeNChildren(node) > 1 )
         return SCIP_OKAY;
      assert(SCIPexprgraphGetNodeNChildren(node) == 1);

      if( SCIPexprgraphGetNodePolynomialNMonomials(node) > 1 )
         return SCIP_OKAY;
      assert(SCIPexprgraphGetNodePolynomialNMonomials(node) == 1); /* assume simplified, i.e., no constant polynomial */

      monomial = SCIPexprgraphGetNodePolynomialMonomials(node)[0];
      assert(SCIPexprGetMonomialNFactors(monomial) == 1); /* since we have only one children and assume simplified */

      exponent = SCIPexprGetMonomialExponents(monomial)[0];
      if( exponent <= 1.0 )
         return SCIP_OKAY;

      /* if exponent is even integer and child has mixed sign, then cannot do
       * if exponent is even integer and child is always negative, then can do via multiplication by -1.0 */
      childbounds = SCIPexprgraphGetNodeBounds(SCIPexprgraphGetNodeChildren(node)[0]);
      if( SCIPisIntegral(scip, exponent) && ((int)SCIPround(scip, exponent) % 2 == 0) && childbounds.inf < 0.0 )
      {
         if( childbounds.sup > 0.0 )
            return SCIP_OKAY;
         signpowcoef = -1.0;
      }

      constant = SCIPexprgraphGetNodePolynomialConstant(node);
      signpowcoef *= SCIPexprGetMonomialCoef(monomial);

      break;
   }

   default:
      return SCIP_OKAY;
   }  /*lint !e788*/
   assert(SCIPexprgraphGetNodeNChildren(node) == 1);

   if( SCIPexprgraphHasNodeNonlinearAncestor(node) )
      return SCIP_OKAY;
   if( !SCIPexprgraphHasNodeSibling(node) )
      return SCIP_OKAY;

   SCIPdebugMessage("reformulate node %p via signpower\n", (void*)node);

   /* get x and its offset */
   child = SCIPexprgraphGetNodeChildren(node)[0];
   if( SCIPexprgraphGetNodeOperator(child) == SCIP_EXPR_VARIDX )
   {
      x = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, child);
      xoffset = 0.0;
   }
   else if( SCIPexprgraphGetNodeOperator(child) == SCIP_EXPR_LINEAR && SCIPexprgraphGetNodeNChildren(child) == 1 )
   {
      SCIP_Real xcoef;

      x = (SCIP_VAR*)SCIPexprgraphGetNodeVar(exprgraph, SCIPexprgraphGetNodeChildren(child)[0]);
      xcoef = SCIPexprgraphGetNodeLinearCoefs(child)[0];
      assert(!SCIPisZero(scip, xcoef));

      signpowcoef *= (xcoef < 0.0 ? -1.0 : 1.0) * pow(REALABS(xcoef), exponent);
      xoffset = SCIPexprgraphGetNodeLinearConstant(child) / xcoef;
   }
   else
   {
      /* reformulate by adding auxiliary variable and constraint for child */
      SCIP_INTERVAL bounds;
      SCIP_VAR* auxvar;
      SCIP_Real minusone;

      bounds = SCIPexprgraphGetNodeBounds(child);
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%dsp", *naddcons);

      SCIPdebugMessage("add auxiliary variable and constraint %s for node %p(%d,%d)\n", name, (void*)child, SCIPexprgraphGetNodeDepth(child), SCIPexprgraphGetNodePosition(child));

      SCIP_CALL( SCIPcreateVar(scip, &auxvar, name, SCIPintervalGetInf(bounds), SCIPintervalGetSup(bounds), 0.0, SCIP_VARTYPE_CONTINUOUS,
            TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, auxvar) );

      /* create new constraint child == auxvar */
      minusone = -1.0;
      SCIP_CALL( SCIPcreateConsNonlinear2(scip, &cons, name, 1, &auxvar, &minusone, child, 0.0, 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      ++*naddcons;

      /* use auxvar to setup signpower constraint */
      x = auxvar;
      xoffset = 0.0;

      SCIP_CALL( SCIPdebugAddSolVal(scip, auxvar, SCIPexprgraphGetNodeVal(child)) );  /*lint !e506 !e774*/

      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      SCIP_CALL( SCIPreleaseVar(scip, &auxvar) );
   }

   /* create auxiliary variable z and add to expression graph */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nlreform%dsp", *naddcons);
   SCIP_CALL( SCIPcreateVar(scip, &z, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS,
         TRUE, TRUE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, z) );
   SCIP_CALL( SCIPexprgraphAddVars(exprgraph, 1, (void**)&z, reformnode) );

   /* setup a absolute power constraint */
   if( REALABS(signpowcoef) * SCIPfeastol(scip) < 1.0 )
   {
      /* if signpowcoef is not huge (<10^6), then put it into absolute power constraint */
      SCIP_CALL( SCIPcreateConsAbspower(scip, &cons, name,
            x, z, exponent, xoffset, -1.0/signpowcoef, -constant/signpowcoef, -constant/signpowcoef,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIPprintCons(scip, cons, NULL) );
      ++*naddcons;

      /* compute value of z and reformnode and set in debug solution and expression graph, resp. */
#ifdef SCIP_DEBUG_SOLUTION
      {
         SCIP_Real xval;
         SCIP_Real zval;

         SCIP_CALL( SCIPdebugGetSolVal(scip, x, &xval) );
         zval = signpowcoef * SIGN(xval + xoffset) * pow(REALABS(xval + xoffset), exponent) + constant;

         SCIP_CALL( SCIPdebugAddSolVal(scip, z, zval) );
         SCIPexprgraphSetVarNodeValue(*reformnode, zval);
      }
#endif
   }
   else
   {
      /* if signpowcoef is huge, then avoid very small coefficient of z
       * instead create additional node on top of current reformnode */
      SCIP_EXPRGRAPHNODE* linnode;

      SCIP_CALL( SCIPcreateConsAbspower(scip, &cons, name,
            x, z, exponent, xoffset, -1.0, 0.0, 0.0,
            TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIPdebug( SCIPprintCons(scip, cons, NULL) );
      ++*naddcons;

      /* compute value of z and reformnode and set in debug solution and expression graph, resp. */
#ifdef SCIP_DEBUG_SOLUTION
      {
         SCIP_Real xval;
         SCIP_Real zval;

         SCIP_CALL( SCIPdebugGetSolVal(scip, x, &xval) );
         zval = SIGN(xval + xoffset) * pow(REALABS(xval + xoffset), exponent);

         SCIP_CALL( SCIPdebugAddSolVal(scip, z, zval) );
         SCIPexprgraphSetVarNodeValue(*reformnode, zval);
      }
#endif

      SCIP_CALL( SCIPexprgraphCreateNodeLinear(SCIPblkmem(scip), &linnode, 1, &signpowcoef, constant) );
      SCIP_CALL( SCIPexprgraphAddNode(exprgraph, linnode, -1, 1, reformnode) );

      *reformnode = linnode;
   }

   SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );

   return SCIP_OKAY;
}

/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 1
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyAbspower)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrAbspower(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}
#else
#define conshdlrCopyAbspower NULL
#endif

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_CONSFREE(consFreeAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip     != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}
#else
#define consFreeAbspower NULL
#endif


/** initialization method of constraint handler (called after problem was transformed) */
#if 1
static
SCIP_DECL_CONSINIT(consInitAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = SCIPfindHeur(scip, "subnlp");
   conshdlrdata->trysolheur = SCIPfindHeur(scip, "trysol");
   conshdlrdata->conshdlrindicator = SCIPfindConshdlr(scip, "indicator");

   return SCIP_OKAY;
}
#else
#define consInitAbspower NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 1
static
SCIP_DECL_CONSEXIT(consExitAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->subnlpheur = NULL;
   conshdlrdata->trysolheur = NULL;
   conshdlrdata->conshdlrindicator = NULL;

   return SCIP_OKAY;
}
#else
#define consExitAbspower NULL
#endif


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 1
static
SCIP_DECL_CONSINITPRE(consInitpreAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* initialize comparedpairwise flag to TRUE, if at most one constraint, otherwise 0 */
   conshdlrdata->comparedpairwise = (nconss <= 1);

   return SCIP_OKAY;
}
#else
#define consInitpreAbspower NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#if 1
static
SCIP_DECL_CONSEXITPRE(consExitpreAbspower)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_RESULT replaceresult;
   int dummy;
   int c;

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* ensure that there are no fixed variables in the constraint (except for multiaggregated) */
      replaceresult = SCIP_DIDNOTFIND;
      SCIP_CALL( checkFixedVariables(scip, conshdlr, conss[c], &dummy, &dummy, &dummy, &dummy, &replaceresult) );  /*lint !e613*/
      if( replaceresult == SCIP_CUTOFF )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if( SCIPconsIsDeleted(conss[c]) )  /*lint !e613*/
         continue;

      /* tell SCIP that we have something nonlinear, and whether we are nonlinear in a continuous variable */
      if( SCIPvarGetType(consdata->x) >= SCIP_VARTYPE_CONTINUOUS )
         SCIPmarkContinuousNonlinearitiesPresent(scip);
      else
         SCIPmarkNonlinearitiesPresent(scip);
   }

   return SCIP_OKAY;
}
#else
#define consExitpreAbspower NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 1
static
SCIP_DECL_CONSINITSOL(consInitsolAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int            c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert(scip  != NULL);
   assert(conss != NULL || nconss == 0);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      assert(consdata->exponent > 1.0);

      /* setup root that corresponds to exponent */
      if( SCIPisIntegral(scip, consdata->exponent) && consdata->exponent-0.5 < ROOTS_KNOWN )
      {
         consdata->root = roots[(int)SCIPfloor(scip, consdata->exponent+0.5)];
      }
      else if( SCIPisEQ(scip, consdata->exponent, 1.852) )
      {
         consdata->root = 0.398217;
      }
      else
      {
         SCIP_Real root;
         SCIP_Real polyval;
         SCIP_Real gradval;
         int iter;

         /* search for a positive root of (n-1) y^n + n y^(n-1) - 1
          * use the closest precomputed root as starting value */
         if( consdata->exponent >= ROOTS_KNOWN )
            root = roots[ROOTS_KNOWN];
         else if( consdata->exponent <= 2.0 )
            root = roots[2];
         else
            root = roots[(int)SCIPfloor(scip, consdata->exponent)];
         iter = 0;
         do
         {
            polyval = (consdata->exponent - 1.0) * consdata->pow(root, consdata->exponent) + consdata->exponent * pow(root, consdata->exponent-1.0) - 1.0;
            if( SCIPisZero(scip, polyval) )
               break;

            /* gradient of (n-1) y^n + n y^(n-1) - 1 is n(n-1)y^(n-1) + n(n-1)y^(n-2) */
            gradval = (consdata->exponent - 1.0) * consdata->exponent * (pow(root, consdata->exponent - 1.0) + pow(root, consdata->exponent - 2.0));
            if( SCIPisZero(scip, gradval) )
               break;

            /* update root by adding -polyval/gradval (Newton's method) */
            root -= polyval / gradval;
            if( root < 0.0 )
               root = 0.0;
         }
         while( ++iter < 1000 );

         if( !SCIPisZero(scip, polyval) )
         {
            SCIPerrorMessage("failed to compute root for exponent %g\n", consdata->exponent);
            return SCIP_ERROR;
         }
         SCIPdebugMessage("root for %g is %.20g, certainty = %g\n", consdata->exponent, root, polyval);
         /* @todo cache root value?? (they are actually really fast to compute...) */

         consdata->root = root;
      }

      if( conshdlrdata->addvarbounds )
      {
         SCIP_CALL( addVarboundsToCutPool(scip, consdata->x) );
         SCIP_CALL( addVarboundsToCutPool(scip, consdata->z) );
      }

      /* add nlrow respresentation to NLP, if NLP had been constructed */
      if( SCIPisNLPConstructed(scip) )
      {
         if( consdata->nlrow == NULL )
         {
            SCIP_CALL( createNlRow(scip, conss[c]) );  /*lint !e613*/
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
#else
#define consInitsolAbspower NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 1
static
SCIP_DECL_CONSEXITSOL(consExitsolAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip  != NULL);
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
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }
   }

   return SCIP_OKAY;
}
#else
#define consExitsolAbspower NULL
#endif


/** frees specific constraint data */
#if 1
static
SCIP_DECL_CONSDELETE(consDeleteAbspower)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   assert(consdata != NULL);
   assert((*consdata)->x != NULL);
   assert((*consdata)->z != NULL);
   assert((*consdata)->xeventfilterpos == -1);
   assert((*consdata)->zeventfilterpos == -1);

   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeMemory(scip, consdata);

   return SCIP_OKAY;
}
#else
#define consDeleteAbspower NULL
#endif


/** transforms constraint data into data belonging to the transformed problem */
#if 1
static
SCIP_DECL_CONSTRANS(consTransAbspower)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   SCIP_CALL( SCIPduplicateMemory(scip, &targetdata, sourcedata) );
   assert(targetdata->xeventfilterpos == -1);
   assert(targetdata->zeventfilterpos == -1);

   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->x, &targetdata->x) );
   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->z, &targetdata->z) );

   /* branching on multiaggregated variables does not seem to work well, so avoid multiagg. x */
   assert( SCIPvarIsActive(targetdata->x) );
   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, targetdata->x) );

   /* cannot propagate on multiaggregated vars, so avoid multiagg. z */
   assert( SCIPvarIsActive(targetdata->z) );
   SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, targetdata->z) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}
#else
#define consTransAbspower NULL
#endif


/** LP initialization method of constraint handler
 *
 * we add secant underestimators
 */
#if 1
static
SCIP_DECL_CONSINITLP(consInitlpAbspower)
{  /*lint --e{715}*/
   SCIP_CONSDATA*     consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_ROW*          row;
   int                c;
   SCIP_Real          xlb;
   SCIP_Real          xub;

   assert(scip  != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      xlb = SCIPvarGetLbGlobal(consdata->x);
      xub = SCIPvarGetUbGlobal(consdata->x);

      if( SCIPisRelEQ(scip, xlb, xub) )
         continue;

      if( !SCIPisInfinity(scip,  consdata->rhs) )
      {
         if( !SCIPisInfinity(scip, -xlb) )
         {
            if( SCIPisNegative(scip, xlb + consdata->xoffset) )
            {
               /* generate secant between xlb and right changepoint */
               SCIP_CALL( generateSecantCutNoCheck(scip, &row, xlb, MIN(-consdata->root * (xlb+consdata->xoffset) - consdata->xoffset, xub), consdata->exponent, consdata->xoffset, consdata->pow, 1.0, consdata->zcoef, consdata->rhs, consdata->x, consdata->z) );
               assert(row != NULL);
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row) < conshdlrdata->cutmaxrange )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );

                  if( conshdlrdata->conshdlrindicator != NULL && !SCIProwIsLocal(row) )
                  {
                     SCIP_CALL( SCIPaddRowIndicator(scip, conshdlrdata->conshdlrindicator, row) );
                  }
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
            else if( xlb < INITLPMAXVARVAL )
            {
               /* generate tangent in lower bound */
               SCIP_CALL( generateLinearizationCut(scip, &row, xlb, consdata->exponent, consdata->xoffset, 1.0, consdata->zcoef, consdata->rhs, consdata->x, consdata->z, FALSE) );
               assert(row != NULL);
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row) < conshdlrdata->cutmaxrange )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );

                  if( conshdlrdata->conshdlrindicator != NULL )
                  {
                     SCIP_CALL( SCIPaddRowIndicator(scip, conshdlrdata->conshdlrindicator, row) );
                  }
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
         }

         if( !SCIPisInfinity(scip, xub) )
         {
            /* generate tangent in upper bound */
            if( -consdata->root * (xlb+consdata->xoffset) - consdata->xoffset < xub && xub <= INITLPMAXVARVAL )
            {
               SCIP_CALL( generateLinearizationCut(scip, &row, xub, consdata->exponent, consdata->xoffset, 1.0, consdata->zcoef, consdata->rhs, consdata->x, consdata->z, FALSE) );
               assert(row != NULL);
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row) < conshdlrdata->cutmaxrange )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );

                  if( conshdlrdata->conshdlrindicator != NULL )
                  {
                     SCIP_CALL( SCIPaddRowIndicator(scip, conshdlrdata->conshdlrindicator, row) );
                  }
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
         }
      }

      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         if( !SCIPisInfinity(scip, xub) )
         {
            if( SCIPisPositive(scip, xub + consdata->xoffset) )
            {
               /* generate secant between left change point and upper bound */
               SCIP_CALL( generateSecantCutNoCheck(scip, &row, -xub, MIN(consdata->root * (xub+consdata->xoffset) + consdata->xoffset, -xlb), consdata->exponent, -consdata->xoffset, consdata->pow, -1.0, -consdata->zcoef, -consdata->lhs, consdata->x, consdata->z) );
               assert(row != NULL);
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row) < conshdlrdata->cutmaxrange )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );

                  if( conshdlrdata->conshdlrindicator != NULL && !SCIProwIsLocal(row) )
                  {
                     SCIP_CALL( SCIPaddRowIndicator(scip, conshdlrdata->conshdlrindicator, row) );
                  }
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
            else if( xub >= -INITLPMAXVARVAL )
            {
               /* generate tangent in upper bound */
               SCIP_CALL( generateLinearizationCut(scip, &row, -xub, consdata->exponent, -consdata->xoffset, -1.0, -consdata->zcoef, -consdata->lhs, consdata->x, consdata->z, FALSE) );
               assert(row != NULL);
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row) < conshdlrdata->cutmaxrange )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );

                  if( conshdlrdata->conshdlrindicator != NULL )
                  {
                     SCIP_CALL( SCIPaddRowIndicator(scip, conshdlrdata->conshdlrindicator, row) );
                  }
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
         }

         if( !SCIPisInfinity(scip, -xlb) )
         {
            /* generate tangent in lower bound */
            if( -consdata->root * (xub+consdata->xoffset) - consdata->xoffset > xlb && xlb >= -INITLPMAXVARVAL )
            {
               SCIP_CALL( generateLinearizationCut(scip, &row, -xlb, consdata->exponent, -consdata->xoffset, -1.0, -consdata->zcoef, -consdata->lhs, consdata->x, consdata->z, FALSE) );
               assert(row != NULL);
               if( !SCIPisInfinity(scip, SCIProwGetRhs(row)) && SCIPgetRowMaxCoef(scip, row)/SCIPgetRowMinCoef(scip, row) < conshdlrdata->cutmaxrange )
               {
                  SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE /* forcecut */) );

                  if( conshdlrdata->conshdlrindicator != NULL )
                  {
                     SCIP_CALL( SCIPaddRowIndicator(scip, conshdlrdata->conshdlrindicator, row) );
                  }
               }
               SCIP_CALL( SCIPreleaseRow(scip, &row) );
            }
         }
      }
   }

   return SCIP_OKAY;
}
#else
#define consInitlpAbspower NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 1
static
SCIP_DECL_CONSSEPALP(consSepalpAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          success;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   /* at root, check if we want to solve the NLP relaxation and use its solutions as reference point
    * if there is something convex, then linearizing in the solution of the NLP relaxation can be very useful
    */
   if( SCIPgetDepth(scip) == 0 && !conshdlrdata->sepanlp &&
      (SCIPgetNContVars(scip) >= conshdlrdata->sepanlpmincont * SCIPgetNVars(scip) || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY) &&
      SCIPisNLPConstructed(scip) && SCIPgetNNlpis(scip) > 0 )
   {
      SCIP_CONSDATA*  consdata;
      SCIP_NLPSOLSTAT solstat;
      SCIP_Bool       solvednlp;
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

            if( (!SCIPisGT(scip, SCIPvarGetUbGlobal(consdata->x), -consdata->xoffset) && !SCIPisInfinity(scip, -consdata->lhs)) ||
               ( !SCIPisLT(scip, SCIPvarGetLbGlobal(consdata->x), -consdata->xoffset) && !SCIPisInfinity(scip, -consdata->rhs)) )
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

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, conshdlrdata->mincutefficacysepa, FALSE, conshdlrdata->sepainboundsonly, &success, NULL) );
   if( success )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}
#else
#define consSepalpAbspower NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 1
static
SCIP_DECL_CONSSEPASOL(consSepasolAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_Bool          success;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(sol      != NULL);
   assert(result   != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeViolations(scip, conss, nconss, sol, &maxviolcon) );
   if( maxviolcon == NULL )
      return SCIP_OKAY;

   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, conshdlrdata->mincutefficacysepa, FALSE, FALSE, &success, NULL) );
   if( success )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}
#else
#define consSepasolAbspower NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcons;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          success;
   SCIP_Real          minefficacy;
   SCIP_Real          sepaefficacy;
   SCIP_Real          maxviol;
   int                nnotify;
   int                c;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, &maxviolcons) );

   if( maxviolcons == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

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

   /* run domain propagation for violated constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_Bool cutoff;
      int       nchgbds;
      int       ndelconss;
      int       naddconss;

      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      nchgbds = 0;
      ndelconss = 0;
      naddconss = 0;
      SCIP_CALL( propagateCons(scip, conshdlr, conss[c], TRUE, &cutoff, &nchgbds, &ndelconss, &naddconss) );  /*lint !e613*/
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( nchgbds )
         *result = SCIP_REDUCEDDOM;
      if( naddconss )
         *result = SCIP_CONSADDED;
   }
   if( *result == SCIP_REDUCEDDOM || *result == SCIP_CONSADDED )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(maxviolcons);
   assert(consdata != NULL);
   maxviol = consdata->lhsviol + consdata->rhsviol;
   assert(SCIPisGT(scip, maxviol, SCIPfeastol(scip)));

   /* we would like a cut that is efficient enough that it is not redundant in the LP (>feastol)
    * however, if the maximal violation is very small, also the best cut efficacy cannot be large
    * thus, in the latter case, we are also happy if the efficacy is at least, say, 75% of the maximal violation
    * but in any case we need an efficacy that is at least feastol
    */
   minefficacy = MIN(0.75*maxviol, conshdlrdata->mincutefficacyenfofac * SCIPfeastol(scip));  /*lint !e666*/
   minefficacy = MAX(minefficacy, SCIPfeastol(scip));  /*lint !e666*/
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, minefficacy, TRUE, FALSE, &success, &sepaefficacy) );
   if( success )
   {
      SCIPdebugMessage("separation succeeded (bestefficacy = %g, minefficacy = %g)\n", sepaefficacy, minefficacy);
      *result = SCIP_SEPARATED;
      return SCIP_OKAY;
   }
   SCIPdebugMessage("separation failed (bestefficacy = %g < %g = minefficacy ); max viol: %g\n", sepaefficacy, minefficacy, maxviol);

   /* we are not feasible, the whole node is not infeasible, and we cannot find a reasonable cut
    * -> collect variables for branching
    */
   SCIP_CALL( registerBranchingCandidates(scip, conshdlr, conss, nconss, &nnotify) );

   if( nnotify == 0 && !solinfeasible && minefficacy > SCIPfeastol(scip) )
   {
      /* fallback 1: we also have no branching candidates, so try to find a weak cut */
      SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, SCIPfeastol(scip), TRUE, FALSE, &success, &sepaefficacy) );
      if( success )
      {
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   }

   if( nnotify == 0 && !solinfeasible )
   {
      /* fallback 2: separation probably failed because of numerical difficulties with a convex constraint;
         if noone declared solution infeasible yet and we had not even found a weak cut, try to resolve by branching */
      SCIP_VAR* brvar = NULL;
      SCIP_CALL( registerLargeLPValueVariableForBranching(scip, conss, nconss, &brvar) );
      if( brvar == NULL )
      {
         SCIPwarningMessage("Could not find any branching variable candidate. Cutting off node. Max viol = %g.\n", SCIPconsGetData(maxviolcons)->lhsviol+SCIPconsGetData(maxviolcons)->rhsviol);
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage("Could not find any usual branching variable candidate. Proposed variable %s with LP value %g for branching. Max. viol. cons. <%s>: %g+%g\n", SCIPvarGetName(brvar), SCIPgetSolVal(scip, NULL, brvar), SCIPconsGetName(maxviolcons), SCIPconsGetData(maxviolcons)->lhsviol, SCIPconsGetData(maxviolcons)->rhsviol);
         nnotify = 1;
      }
   }

   assert(*result == SCIP_INFEASIBLE && (solinfeasible || nnotify > 0));
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS*         maxviolcon;
   SCIP_CONSDATA*     consdata;
   int                c;
   int                nnotify;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( computeViolations(scip, conss, nconss, NULL, &maxviolcon) );
   if( maxviolcon == NULL )
   {
      *result = SCIP_FEASIBLE;
      return SCIP_OKAY;
   }

   *result = SCIP_INFEASIBLE;

   /* run domain propagation for violated constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_Bool cutoff;
      int       nchgbds;
      int       ndelconss;
      int       naddconss;

      assert(conss[c] != NULL);  /*lint !e613*/

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      nchgbds = 0;
      ndelconss = 0;
      naddconss = 0;
      SCIP_CALL( propagateCons(scip, conshdlr, conss[c], TRUE, &cutoff, &nchgbds, &ndelconss, &naddconss) );  /*lint !e613*/
      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( nchgbds )
         *result = SCIP_REDUCEDDOM;
      if( naddconss )
         *result = SCIP_CONSADDED;
   }
   if( *result == SCIP_REDUCEDDOM || *result == SCIP_CONSADDED )
      return SCIP_OKAY;

   /* we are not feasible and we cannot proof that the whole node is infeasible
    * -> branch on all unfixed variables in violated constraints
    */
   nnotify = 0;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);
      SCIPdebugMessage("cons <%s> violation: %g %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);

      if( !SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) && !SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
         continue;

      SCIPdebugMessage("cons <%s> violation: %g %g\n", SCIPconsGetName(conss[c]), consdata->lhsviol, consdata->rhsviol);

      /* domain propagation should have removed cons when x is fixed */
      assert(!SCIPisRelEQ(scip, SCIPvarGetLbLocal(consdata->x), SCIPvarGetUbLocal(consdata->x)));

      SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->x, consdata->lhsviol + consdata->rhsviol, proposeBranchingPoint(scip, conss[c], conshdlrdata->preferzerobranch, conshdlrdata->branchminconverror)) );
      ++nnotify;
   }

   if( nnotify == 0 )
   {
      SCIPdebugMessage("All variables in violated constraints fixed (up to epsilon). Cannot find branching candidate. Forcing solution of LP.\n");
      *result = SCIP_SOLVELP;
   }

   assert(*result == SCIP_SOLVELP || (*result == SCIP_INFEASIBLE && nnotify > 0));
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 1
static
SCIP_DECL_CONSPROP(consPropAbspower)
{  /*lint --e{715}*/
   int         c;
   int         nchgbds;
   int         ndelconss;
   int         naddconss;
   SCIP_Bool   cutoff = FALSE;

   assert(scip     != NULL);
   assert(conshdlr != NULL);
   assert(conss    != NULL || nconss == 0);
   assert(result   != NULL);

   *result = SCIP_DIDNOTFIND;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);

      /* propagate constraint, but do not allow to add a constraint for tightening a multiaggregated variable (not allowed in CONSPROP) */
      nchgbds = 0;
      ndelconss = 0;
      naddconss = 0;
      SCIP_CALL( propagateCons(scip, conshdlr, conss[c], FALSE, &cutoff, &nchgbds, &ndelconss, &naddconss) );
      assert(naddconss == 0);

      if( cutoff )
      {
         *result = SCIP_CUTOFF;
         break;
      }

      if( nchgbds )
         *result = SCIP_REDUCEDDOM;

      if( c >= nusefulconss && *result != SCIP_DIDNOTFIND )
         break;
   }

   return SCIP_OKAY;
}
#else
#define consPropAbspower NULL
#endif

/** presolving method of constraint handler */
#if 1
static
SCIP_DECL_CONSPRESOL(consPresolAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_RESULT    replaceresult;
   SCIP_Bool      success;
   SCIP_Bool      infeas;
   int            localnchgbds;
   int            localndelconss;
   int            localnaddconss;
   int            c;

   assert(scip   != NULL);
   assert(conss  != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;

   /* check for duplicates, if not done yet or if absolute power constraints were modified (variable fixings) or new absolute power constraints had been added */
   if( !conshdlrdata->comparedpairwise )
   {
      SCIP_CALL( presolveFindDuplicates(scip, conshdlr, conss, nconss, nupgdconss, ndelconss, naddconss, nfixedvars, naggrvars, &success, &infeas) );
      if( infeas )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( success )
      {
         *result = SCIP_SUCCESS;
      }
   }

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);  /*lint !e613*/

      if( SCIPconsIsDeleted(conss[c]) )  /*lint !e613*/
         continue;

      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      SCIPdebugMessage("presolving constraint <%s>\n", SCIPconsGetName(conss[c]));  /*lint !e613*/
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );  /*lint !e613*/

      /* check if we can upgrade to a linear constraint */
      if( consdata->exponent == 1.0 )
      {
         SCIP_VAR*  vars[2];
         SCIP_Real  coefs[2];
         SCIP_CONS* lincons;
         SCIP_Real  lhs;
         SCIP_Real  rhs;

         vars[0] = consdata->x;
         vars[1] = consdata->z;
         coefs[0] = 1.0;
         coefs[1] = consdata->zcoef;
         lhs = consdata->lhs;
         rhs = consdata->rhs;
         if( !SCIPisInfinity(scip, -lhs) )
            lhs -= consdata->xoffset;
         if( !SCIPisInfinity(scip,  rhs) )
            rhs -= consdata->xoffset;

         SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, SCIPconsGetName(conss[c]), 2, vars, coefs, lhs, rhs,
               SCIPconsIsInitial(conss[c]), SCIPconsIsSeparated(conss[c]), SCIPconsIsEnforced(conss[c]),
               SCIPconsIsChecked(conss[c]), SCIPconsIsPropagated(conss[c]),  SCIPconsIsLocal(conss[c]),
               SCIPconsIsModifiable(conss[c]), SCIPconsIsDynamic(conss[c]), SCIPconsIsRemovable(conss[c]),
               SCIPconsIsStickingAtNode(conss[c])) );  /*lint !e613*/
         SCIP_CALL( SCIPaddCons(scip, lincons) );
         SCIP_CALL( SCIPreleaseCons(scip, &lincons) );

         SCIP_CALL( SCIPdelCons(scip, conss[c]) );  /*lint !e613*/
         ++*nupgdconss;
         continue;
      }

      /* check for fixed variables */
      replaceresult = SCIP_DIDNOTFIND;
      SCIP_CALL( checkFixedVariables(scip, conshdlr, conss[c], ndelconss, nupgdconss, nchgbds, nfixedvars, &replaceresult) );  /*lint !e613*/
      switch( replaceresult )
      {
      case SCIP_DIDNOTFIND:
         break;

      case SCIP_CUTOFF:
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;

      case SCIP_REDUCEDDOM:
      case SCIP_CONSADDED:
         *result = SCIP_SUCCESS;
         break;

      default:
         SCIPerrorMessage("invalid result from checkFixedVariables\n");
         SCIPABORT();
      }  /*lint !e788*/
      if( SCIPconsIsDeleted(conss[c]) )  /*lint !e613*/
      {
         *result = SCIP_SUCCESS;
         continue;
      }

      /* run domain propagation, also checks for redundancy */
      localnchgbds = 0;
      localndelconss = 0;
      localnaddconss = 0;
      SCIP_CALL( propagateCons(scip, conshdlr, conss[c], TRUE, &infeas, &localnchgbds, &localndelconss, &localnaddconss) );  /*lint !e613*/
      if( infeas )
      {
         SCIPdebugMessage("propagation on constraint <%s> says problem is infeasible in presolve\n", SCIPconsGetName(conss[c]));  /*lint !e613*/
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if( localnchgbds > 0 || localndelconss > 0 || localnaddconss > 0 )
      {
         *nchgbds   += localnchgbds;
         *ndelconss += localndelconss;
         *naddconss += localnaddconss;
         *result = SCIP_SUCCESS;
      }

      /* propagate variable bound constraints */
      if( !consdata->propvarbounds )
      {
         SCIP_CALL( propagateVarbounds(scip, conshdlr, conss[c], &infeas, nchgbds) );  /*lint !e613*/

         if( infeas )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         consdata->propvarbounds = TRUE;
      }

      /* check if we can make z implicit integer
       * if constraint is signpow(x,n) + c*z = rhs with x integer, |c| = 1, rhs and n integral, then z is implicit integral
       */
      if( SCIPvarGetType(consdata->z) == SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(consdata->x) != SCIP_VARTYPE_CONTINUOUS &&
         SCIPisEQ(scip, consdata->lhs, consdata->rhs) && SCIPisIntegral(scip, consdata->rhs) && SCIPisEQ(scip, REALABS(consdata->zcoef), 1.0) && SCIPisIntegral(scip, consdata->exponent)
         )
      {
         SCIPdebugMessage("make z = <%s> implicit integer in cons <%s>\n", SCIPvarGetName(consdata->z), SCIPconsGetName(conss[c]));  /*lint !e613*/
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );
         SCIP_CALL( SCIPchgVarType(scip, consdata->z, SCIP_VARTYPE_IMPLINT, &infeas) );
         if( infeas )
         {
            SCIPdebugMessage("problem found infeasible in presolve when making <%s> implicit integer\n", SCIPvarGetName(consdata->z));
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else
         {
            ++*nchgvartypes;
         }
      }
   }

   return SCIP_OKAY;
}
#else
#define consPresolAbspower NULL
#endif


/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) left hand side and bounds on z -> lower bound on x
 *   (2) left hand side and upper bound on x -> bound on z
 *   (3) right hand side and bounds on z -> upper bound on x
 *   (4) right hand side and lower bound on x -> bound on z
 */
#if 1
static
SCIP_DECL_CONSRESPROP(consRespropAbspower)
{
   assert(result != NULL);

   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, boundtype, bdchgidx) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}  /*lint !e715*/
#else
#define consRespropAbspower NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockAbspower)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Bool      haslb;
   SCIP_Bool      hasub;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   haslb = !SCIPisInfinity(scip, -consdata->lhs);
   hasub = !SCIPisInfinity(scip,  consdata->rhs);

   if( consdata->x != NULL )
   {
      if( haslb )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->x, nlockspos, nlocksneg) );
      }
      if( hasub )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->x, nlocksneg, nlockspos) );
      }
   }

   if( consdata->z != NULL )
   {
      if( consdata->zcoef > 0 )
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlockspos, nlocksneg) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlocksneg, nlockspos) );
         }
      }
      else
      {
         if( haslb )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlocksneg, nlockspos) );
         }
         if( hasub )
         {
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->z, nlockspos, nlocksneg) );
         }
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 1
static
SCIP_DECL_CONSACTIVE(consActiveAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* (re)run constraint comparison, since new constraint is added */
   conshdlrdata->comparedpairwise = FALSE;

   return SCIP_OKAY;
}
#else
#define consActiveAbspower NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveAbspower)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of absolute power constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveAbspower NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 1
static
SCIP_DECL_CONSENABLE(consEnableAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );

   return SCIP_OKAY;
}
#else
#define consEnableAbspower NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 1
static
SCIP_DECL_CONSDISABLE(consDisableAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->eventhdlr != NULL);

   SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );

   return SCIP_OKAY;
}
#else
#define consDisableAbspower NULL
#endif


/** variable deletion method of constraint handler */
#define consDelvarsAbspower NULL

/** constraint display method of constraint handler */
#if 1
static
SCIP_DECL_CONSPRINT(consPrintAbspower)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

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
   SCIPinfoMessage(scip, file, "signpower(");
   SCIP_CALL( SCIPwriteVarName(scip, file, consdata->x, TRUE) );
   SCIPinfoMessage(scip, file, " %+.15g, %.15g) ", consdata->xoffset, consdata->exponent);

   SCIPinfoMessage(scip, file, "%+.15g", consdata->zcoef);
   SCIP_CALL( SCIPwriteVarName(scip, file, consdata->z, TRUE) );

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
#else
#define consPrintAbspower NULL
#endif

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckAbspower)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Bool          dolinfeasshift;
   SCIP_Real          maxviol;
   SCIP_Real          viol;
   int                c;

   assert(scip   != NULL);
   assert(conss  != NULL || nconss == 0);
   assert(sol    != NULL);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;

   maxviol = 0.0;

   dolinfeasshift = conshdlrdata->linfeasshift && (conshdlrdata->trysolheur != NULL) && SCIPgetStage(scip) > SCIP_STAGE_PROBLEM && SCIPgetStage(scip) < SCIP_STAGE_SOLVED;
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, &viol) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIPinfoMessage(scip, NULL, "absolute power constraint <%s> violated by %g (scaled = %g)\n\t",
               SCIPconsGetName(conss[c]), viol, MAX(consdata->lhsviol, consdata->rhsviol));
            SCIP_CALL( consPrintAbspower(scip, conshdlr, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
         }

         if( conshdlrdata->subnlpheur == NULL && !dolinfeasshift )
            return SCIP_OKAY;
         if( consdata->lhsviol > maxviol || consdata->rhsviol > maxviol )
            maxviol = MAX(consdata->lhsviol, consdata->rhsviol);
      }
   }

   if( *result == SCIP_INFEASIBLE && dolinfeasshift )
   {
      SCIP_CALL( proposeFeasibleSolution(scip, conshdlr, conss, nconss, sol) );
   }

   if( *result == SCIP_INFEASIBLE && conshdlrdata->subnlpheur != NULL )
   {
      SCIP_CALL( SCIPupdateStartpointHeurSubNlp(scip, conshdlrdata->subnlpheur, sol, maxviol) );
   }

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
#if 1
static
SCIP_DECL_CONSCOPY(consCopyAbspower)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*      x;
   SCIP_VAR*      z;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(sourcescip != NULL);
   assert(sourcecons != NULL);
   assert(varmap != NULL);
   assert(valid != NULL);

   consdata = SCIPconsGetData(sourcecons);
   assert(consdata != NULL);

   *valid = TRUE;
   *cons  = NULL;

   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->x, &x, varmap, consmap, global, valid) );

   if( *valid )
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, consdata->z, &z, varmap, consmap, global, valid) );
   }

   if( *valid )
   {
      SCIP_CALL( SCIPcreateConsAbspower(scip, cons, name != NULL ? name : SCIPconsGetName(sourcecons),
            x, z, consdata->exponent, consdata->xoffset, consdata->zcoef, consdata->lhs, consdata->rhs,
            initial, separate, enforce, check, propagate, local, FALSE, dynamic, removable, stickingatnode) );  /*lint !e644*/
   }

   return SCIP_OKAY;
}
#else
#define consCopyAbspower NULL
#endif


/** constraint parsing method of constraint handler */
#if 1
static
SCIP_DECL_CONSPARSE(consParseAbspower)
{
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real xoffset;
   SCIP_Real exponent;
   SCIP_Real zcoef;
   SCIP_Real value;
   char*     endptr;
   char      sense;
   SCIP_VAR* x;
   SCIP_VAR* z;

   *success = TRUE;

   /* set right hand and left side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   SCIPdebugMessage("start parsing absolute power constraint expression %s\n", str);

   if( strncmp(str, "signpower(", 10) != 0 )
   {
      /* str does not start with signpower string, so may be left-hand-side of ranged constraint */
      if( !SCIPstrToRealValue(str, &lhs, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Syntax error: left-hand-side or 'signpower(' expected at begin on '%s'\n", str);
         *success = FALSE;
         return SCIP_OKAY;
      }
      str = endptr;
   }
   else
   {
      str += 10;
   }

   /* parse (x +offset, exponent) +coef z */

   /* parse variable name */
   SCIP_CALL( SCIPparseVarName(scip, str, &x, &endptr) );
   if( x == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str = endptr;

   /* skip whitespace */
   while( isspace((int)*str) )
      ++str;

   /* parse offset */
   if( !SCIPstrToRealValue(str, &xoffset, &endptr) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected coefficient at begin of '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str = endptr;

   if( *str != ',' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected ',' at begin of '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   ++str;

   /* skip whitespace */
   while( isspace((int)*str) )
      ++str;

   /* parse exponent */
   if( !SCIPstrToRealValue(str, &exponent, &endptr) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected coefficient at begin of '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str = endptr;

   if( *str != ')' )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected ')' at begin of '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   ++str;

   /* skip whitespace */
   while( isspace((int)*str) )
      ++str;

   /* parse coefficient */
   if( !SCIPstrToRealValue(str, &zcoef, &endptr) )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected coefficient at begin of '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str = endptr;

   /* parse variable name */
   SCIP_CALL( SCIPparseVarName(scip, str, &z, &endptr) );
   if( z == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "unknown variable name at '%s'\n", str);
      *success = FALSE;
      return SCIP_OKAY;
   }
   str = endptr;

   /* skip whitespace */
   while( isspace((int)*str) )
      ++str;

   if( strncmp(str, "[free]", 6) != 0 )
   {
      /* parse sense */
      if( (*str != '<' && *str != '>' && *str != '=') || str[1] != '=' )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected sense at begin of '%s'\n", str);
         *success = FALSE;
         return SCIP_OKAY;
      }
      sense = *str;
      str += 2;

      /* parse value at rhs */
      if( !SCIPstrToRealValue(str, &value, &endptr) )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "expected rhs value at begin of '%s'\n", str);
         *success = FALSE;
         return SCIP_OKAY;
      }

      switch( sense )
      {
      case '<' :
         rhs = value;
         break;
      case '>' :
         lhs = value;
         break;
      case '=' :
         lhs = rhs = value;
         break;
      default:
         SCIPABORT(); /* checked above that this cannot happen */
      }
   }

   SCIP_CALL( SCIPcreateConsAbspower(scip, cons, name, x, z, exponent, xoffset, zcoef, lhs, rhs,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}  /*lint !e715*/
#else
#define consParseAbspower NULL
#endif

/*
 * constraint specific interface methods
 */

/** creates the handler for absolute power constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrAbspower(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create absolute power constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   BMSclearMemory(conshdlrdata);

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyAbspower,
         consFreeAbspower, consInitAbspower, consExitAbspower,
         consInitpreAbspower, consExitpreAbspower, consInitsolAbspower, consExitsolAbspower,
         consDeleteAbspower, consTransAbspower, consInitlpAbspower,
         consSepalpAbspower, consSepasolAbspower, consEnfolpAbspower, consEnfopsAbspower, consCheckAbspower,
         consPropAbspower, consPresolAbspower, consRespropAbspower, consLockAbspower,
         consActiveAbspower, consDeactiveAbspower,
         consEnableAbspower, consDisableAbspower, consDelvarsAbspower,
         consPrintAbspower, consCopyAbspower, consParseAbspower,
         conshdlrdata) );

   /* include the quadratic constraint upgrade in the quadratic constraint handler */
   SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, quadconsUpgdAbspower, QUADCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );

   /* include the absolute power constraint upgrade and node reform in the nonlinear constraint handler
    * we give it higher priority as quadratic, so it also takes care of x^2 constraints, if possible
    */
   SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, nonlinconsUpgdAbspower, exprgraphnodeReformAbspower, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );

   /* add absolute power constraint handler parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacysepa",
         "minimal efficacy for a cut to be added to the LP during separation; overwrites separating/efficacy",
         &conshdlrdata->mincutefficacysepa, FALSE, 0.0001, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/minefficacyenfofac",
         "minimal target efficacy of a cut in order to add it to relaxation during enforcement as factor of feasibility tolerance (may be ignored)",
         &conshdlrdata->mincutefficacyenfofac, FALSE, 2.0, 1.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/cutmaxrange",
         "maximal coef range of a cut (maximal coefficient divided by minimal coefficient) in order to be added to LP relaxation",
         &conshdlrdata->cutmaxrange, FALSE, 1e+7, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/projectrefpoint",
         "whether to project the reference point when linearizing an absolute power constraint in a convex region",
         &conshdlrdata->projectrefpoint, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/"CONSHDLR_NAME"/preferzerobranch",
         "how much to prefer branching on 0.0 when sign of variable is not fixed yet: 0 no preference, 1 prefer if LP solution will be cutoff in both child nodes, 2 prefer always, 3 ensure always",
         &conshdlrdata->preferzerobranch, FALSE, 1, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/branchminconverror",
         "whether to compute branching point such that the convexification error is minimized (after branching on 0.0)",
         &conshdlrdata->branchminconverror, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/addvarbounds",
         "should variable bounds be added to the cutpool?",
         &conshdlrdata->addvarbounds, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/linfeasshift",
         "whether to try to make solutions in check function feasible by shifting the linear variable z",
         &conshdlrdata->linfeasshift, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/"CONSHDLR_NAME"/sepainboundsonly",
         "whether to separate linearization cuts only in the variable bounds (does not affect enforcement)",
         &conshdlrdata->sepainboundsonly, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/"CONSHDLR_NAME"/sepanlpmincont",
         "minimal required fraction of continuous variables in problem to use solution of NLP relaxation in root for separation",
         &conshdlrdata->sepanlpmincont, FALSE, 1.0, 0.0, 2.0, NULL, NULL) );

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME, "signals a bound change on a variable to an absolute power constraint",
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, processVarEvent, NULL) );
   conshdlrdata->eventhdlr = SCIPfindEventhdlr(scip, CONSHDLR_NAME);

   SCIP_CALL( SCIPincludeEventhdlr(scip, CONSHDLR_NAME"_newsolution", "handles the event that a new primal solution has been found",
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, processNewSolutionEvent, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a absolute power constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             x,                  /**< nonlinear variable x in constraint */
   SCIP_VAR*             z,                  /**< linear variable z in constraint */
   SCIP_Real             exponent,           /**< exponent n of |x+offset|^n term in constraint */
   SCIP_Real             xoffset,            /**< offset in |x+offset|^n term in constraint */
   SCIP_Real             zcoef,              /**< coefficient of z in constraint */
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

   assert(x != NULL);
   assert(z != NULL);
   assert(exponent > 1.0);
   assert(!SCIPisZero(scip, zcoef));
   assert(!SCIPisInfinity(scip, REALABS(zcoef)));
   assert(!modifiable); /* we do not support column generation */

   /* find the absolute power constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("absolute power constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory( scip, &consdata) );
   BMSclearMemory(consdata);
   consdata->xeventfilterpos = -1;
   consdata->zeventfilterpos = -1;

   consdata->x        = x;
   consdata->z        = z;
   consdata->xoffset  = xoffset;
   consdata->zcoef    = zcoef;
   consdata->lhs      = lhs;
   consdata->rhs      = rhs;

   if( SCIPisEQ(scip, exponent, 2.0) )
   {
      consdata->exponent = 2.0;
      consdata->pow      = square;
   }
   else
   {
      consdata->exponent = exponent;
      consdata->pow      = pow;
   }

   /* branching on multiaggregated variables does not seem to work well, so try to avoid multiagg. x */
   if( SCIPvarIsActive(x) )
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, x) );

   /* cannot propagate on multiaggregated vars, so avoid multiagg. z */
   if( SCIPvarIsActive(z) )
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, z) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** gets the absolute power constraint as a nonlinear row representation */
SCIP_RETCODE SCIPgetNlRowAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_NLROW**          nlrow               /**< a buffer where to store pointer to nonlinear row */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons  != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);
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

/** gets nonlinear variable x in absolute power constraint */
SCIP_VAR* SCIPgetNonlinearVarAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->x;
}

/** gets linear variable z in absolute power constraint */
SCIP_VAR* SCIPgetLinearVarAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->z;
}

/** gets exponent in power term in absolute power constraint */
SCIP_Real SCIPgetExponentAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->exponent;
}

/** gets offset in power term in absolute power constraint */
SCIP_Real SCIPgetOffsetAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->xoffset;
}

/** gets coefficient of linear variable in absolute power constraint */
SCIP_Real SCIPgetCoefLinearAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->zcoef;
}

/** gets left hand side in absolute power constraint */
SCIP_Real SCIPgetLhsAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side in absolute power constraint */
SCIP_Real SCIPgetRhsAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< absolute power constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets the absolute violation of a absolute power constraint by a solution */
SCIP_Real SCIPgetViolationAbspower(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< absolute power constraint */
   SCIP_SOL*             sol                 /**< LP solution */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real z_val;
   SCIP_Real x_val;
   SCIP_Real rhs;
   SCIP_Real proj_val;

   assert(cons != NULL);
   assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->lhs == 0.0);
   assert(consdata->rhs == 0.0);

   z_val = SCIPgetSolVal(scip, sol, consdata->z);
   x_val = SCIPgetSolVal(scip, sol, consdata->x);

   rhs = -1.0 * consdata->zcoef * z_val;
   proj_val = SIGN(rhs) * pow(REALABS(rhs), 1.0 / consdata->exponent) - consdata->xoffset;

   SCIPdebugMessage("computing slack: linear: %f, power: %f, projected: %f\n", z_val, x_val, proj_val);

   return x_val - proj_val;
}
