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

/**@file   cons_varbound.c
 * @brief  Constraint handler for variable bound constraints \f$lhs \le x + c y \le rhs\f$.
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>

#include "scip/cons_varbound.h"
#include "scip/cons_linear.h"


/* constraint handler properties */
/* @note: despite the description due to aggregations x can be upgraded to a binary variable */
#define CONSHDLR_NAME          "varbound"
#define CONSHDLR_DESC          "variable bounds  lhs <= x + c*y <= rhs, x non-binary, y non-continuous"
#define CONSHDLR_SEPAPRIORITY   +900000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY   -500000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  -500000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             0 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING             SCIP_PROPTIMING_BEFORELP

#define EVENTHDLR_NAME         "varbound"
#define EVENTHDLR_DESC         "bound change event handler for variable bound constraints"

#define LINCONSUPGD_PRIORITY     +50000 /**< priority of the constraint handler for upgrading of linear constraints */

#define DEFAULT_PRESOLPAIRWISE     TRUE /**< should pairwise constraint comparison be performed in presolving? */

#define MAXSCALEDCOEF            1000LL /**< maximal coefficient value after scaling */


/** variable bound constraint data */
struct SCIP_ConsData
{
   SCIP_Real             vbdcoef;            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs;                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs;                /**< right hand side of variable bound inequality */
   SCIP_VAR*             var;                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar;             /**< binary, integer or implicit integer bounding variable y */
   SCIP_ROW*             row;                /**< LP row, if constraint is already stored in LP row format */
   unsigned int          propagated:1;       /**< is the variable bound constraint already propagated? */
   unsigned int          presolved:1;        /**< is the variable bound constraint already presolved? */
   unsigned int          addvarbounds:1;     /**< are the globally valid variable bound are added? */
   unsigned int          changed:1;          /**< was constraint changed since last aggregation round in preprocessing? */
   unsigned int          tightened:1;        /**< were the vbdcoef and all sides already tightened? */
};


/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             presolpairwise;     /**< should pairwise constraint comparison be performed in presolving? */
};

/*
 * Propagation rules
 */

enum Proprule
{
   PROPRULE_1,                          /**< left hand side and bounds on y -> lower bound on x */
   PROPRULE_2,                          /**< left hand side and upper bound on x -> bound on y */
   PROPRULE_3,                          /**< right hand side and bounds on y -> upper bound on x */
   PROPRULE_4,                          /**< right hand side and lower bound on x -> bound on y */
   PROPRULE_INVALID                     /**< propagation was applied without a specific propagation rule */
};
typedef enum Proprule PROPRULE;




/*
 * Local methods
 */

/** compares the index of both variables in a constraint,
 *
 *  returns -1 if index of first consdata->var is smaller than the index of the second consdata->var or both indices are
 *  equal and the index of the first consdata->vbdvar is smaller than the index of the second consdata->vbdvar,
 *
 *  returns 0 if the indices are pairwise equal,
 *
 *  and returns +1 otherwise
 */
static
SCIP_DECL_SORTPTRCOMP(consVarboundComp)
{
   SCIP_CONSDATA* consdata1;
   SCIP_CONSDATA* consdata2;

   assert(elem1 != NULL);
   assert(elem2 != NULL);

   consdata1 = SCIPconsGetData((SCIP_CONS*) elem1);
   consdata2 = SCIPconsGetData((SCIP_CONS*) elem2);

   assert(consdata1 != NULL);
   assert(consdata2 != NULL);

   /* comparison is done over 3 ordered criteria:
    *  (i) variable index of variable 1
    *  (ii) variable index of variable 2.
    *  (iii) changed status
    */
   if( SCIPvarGetIndex(consdata1->var) < SCIPvarGetIndex(consdata2->var)
      || (SCIPvarGetIndex(consdata1->var) == SCIPvarGetIndex(consdata2->var)
         && SCIPvarGetIndex(consdata1->vbdvar) < SCIPvarGetIndex(consdata2->vbdvar))
      || (SCIPvarGetIndex(consdata1->var) == SCIPvarGetIndex(consdata2->var)
         && SCIPvarGetIndex(consdata1->vbdvar) == SCIPvarGetIndex(consdata2->vbdvar)
         && !consdata1->changed && consdata2->changed) )
      return -1;
   else if( SCIPvarGetIndex(consdata1->var) == SCIPvarGetIndex(consdata2->var)
      && SCIPvarGetIndex(consdata1->vbdvar) == SCIPvarGetIndex(consdata2->vbdvar)
      && (consdata1->changed == consdata2->changed) )
      return 0;
   else
      return +1;
}

/** creates constraint handler data for varbound constraint handler */
static
SCIP_RETCODE conshdlrdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to store the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, conshdlrdata) );

   return SCIP_OKAY;
}

/** frees constraint handler data for varbound constraint handler */
static
void conshdlrdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA**   conshdlrdata        /**< pointer to the constraint handler data */
   )
{
   assert(scip != NULL);
   assert(conshdlrdata != NULL);
   assert(*conshdlrdata != NULL);

   SCIPfreeMemory(scip, conshdlrdata);
}

/** catches events for variables */
static
SCIP_RETCODE catchEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< variable bound constraint data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* catch bound change events on variables */
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr,
         (SCIP_EVENTDATA*)consdata, NULL) );
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr,
         (SCIP_EVENTDATA*)consdata, NULL) );

   return SCIP_OKAY;
}

/** drops events for variables */
static
SCIP_RETCODE dropEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata            /**< variable bound constraint data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(consdata != NULL);

   eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   assert(eventhdlr != NULL);

   /* drop events on variables */
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->var, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr,
         (SCIP_EVENTDATA*)consdata, -1) );
   SCIP_CALL( SCIPdropVarEvent(scip, consdata->vbdvar, SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED, eventhdlr,
         (SCIP_EVENTDATA*)consdata, -1) );

   return SCIP_OKAY;
}

/** creates a variable bound constraint data object */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the variable bound constraint data */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs                 /**< right hand side of variable bound inequality */
   )
{
   assert(consdata != NULL);
   assert(SCIPvarGetType(vbdvar) != SCIP_VARTYPE_CONTINUOUS);

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   if( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -rhs) )
      rhs = -SCIPinfinity(scip);

   if( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, lhs) )
      lhs = SCIPinfinity(scip);

   if( SCIPisGT(scip, lhs, rhs) )
   {
      SCIPerrorMessage("left hand side of varbound constraint greater than right hand side\n");
      SCIPerrorMessage(" -> lhs=%g, rhs=%g\n", lhs, rhs);
      return SCIP_INVALIDDATA;
   }

   if( SCIPisInfinity(scip, vbdcoef) )
      vbdcoef = SCIPinfinity(scip);
   else if( SCIPisInfinity(scip, -vbdcoef) )
      vbdcoef = -SCIPinfinity(scip);

   (*consdata)->var = var;
   (*consdata)->vbdvar = vbdvar;
   (*consdata)->vbdcoef = vbdcoef;
   (*consdata)->lhs = lhs;
   (*consdata)->rhs = rhs;
   (*consdata)->row = NULL;
   (*consdata)->propagated = FALSE;
   (*consdata)->presolved = FALSE;
   (*consdata)->addvarbounds = FALSE;
   (*consdata)->changed = TRUE;
   (*consdata)->tightened = FALSE;

   /* if we are in the transformed problem, get transformed variables, add variable bound information, and catch events */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->var, &(*consdata)->var) );
      SCIP_CALL( SCIPgetTransformedVar(scip, (*consdata)->vbdvar, &(*consdata)->vbdvar) );

      /* catch events for variables */
      SCIP_CALL( catchEvents(scip, *consdata) );
   }

   /* capture variables */
   SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->var) );
   SCIP_CALL( SCIPcaptureVar(scip, (*consdata)->vbdvar) );

   return SCIP_OKAY;
}

/** frees a variable bound constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the variable bound constraint */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* release the row */
   if( (*consdata)->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(*consdata)->row) );
   }

   /* drop events */
   if( SCIPisTransformed(scip) )
   {
      SCIP_CALL( dropEvents(scip, *consdata) );
   }

   /* release variables */
   SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->var) );
   SCIP_CALL( SCIPreleaseVar(scip, &(*consdata)->vbdvar) );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** creates LP row corresponding to variable bound constraint */
static
SCIP_RETCODE createRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< variable bound constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->row == NULL);

   SCIP_CALL( SCIPcreateEmptyRow(scip, &consdata->row, SCIPconsGetName(cons), consdata->lhs, consdata->rhs,
         SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), SCIPconsIsRemovable(cons)) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->var, 1.0) );
   SCIP_CALL( SCIPaddVarToRow(scip, consdata->row, consdata->vbdvar, consdata->vbdcoef) );

   return SCIP_OKAY;
}

/** adds linear relaxation of variable bound constraint to the LP */
static
SCIP_RETCODE addRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< variable bound constraint */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   if( !SCIProwIsInLP(consdata->row) )
   {
      SCIPdebugMessage("adding relaxation of variable bound constraint <%s>: ", SCIPconsGetName(cons));
      SCIPdebug( SCIProwPrint(consdata->row, NULL) );
      SCIP_CALL( SCIPaddCut(scip, NULL, consdata->row, FALSE) );
   }

   return SCIP_OKAY;
}

/** separates the given variable bound constraint */
static
SCIP_RETCODE separateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   SCIP_Bool*            separated           /**< pointer to store whether a cut was found */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real feasibility;

   assert(separated != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("separating variable bound constraint <%s>\n", SCIPconsGetName(cons));

   *separated = FALSE;

   /* create LP relaxation if not yet existing */
   if( consdata->row == NULL )
   {
      SCIP_CALL( createRelaxation(scip, cons) );
   }
   assert(consdata->row != NULL);

   /* check non-LP rows for feasibility and add them as cut, if violated */
   if( sol != NULL && !SCIProwIsInLP(consdata->row) )
   {
      feasibility = SCIPgetRowSolFeasibility(scip, consdata->row, sol);
      if( SCIPisFeasNegative(scip, feasibility) )
      {
         SCIP_CALL( SCIPaddCut(scip, sol, consdata->row, FALSE) );
         *separated = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** returns whether the given solution is feasible for the given variable bound constraint */
static
SCIP_Bool checkCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_SOL*             sol,                /**< solution to check, NULL for current solution */
   SCIP_Bool             checklprows         /**< should LP rows be checked? */
   )
{
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("checking variable bound constraint <%s> for feasibility of solution %p (lprows=%u)\n",
      SCIPconsGetName(cons), (void*)sol, checklprows);

   if( checklprows || consdata->row == NULL || !SCIProwIsInLP(consdata->row) )
   {
      SCIP_Real sum;

      sum = SCIPgetSolVal(scip, sol, consdata->var);
      sum += consdata->vbdcoef * SCIPgetSolVal(scip, sol, consdata->vbdvar);

      return SCIPisFeasGE(scip, sum, consdata->lhs) && SCIPisFeasLE(scip, sum, consdata->rhs);
   }
   else
      return TRUE;
}


/** resolves a propagation on the given variable by supplying the variables needed for applying the corresponding
 *  propagation rule (see propagateCons()):
 *   (1) left hand side and bounds on y -> lower bound on x
 *   (2) left hand side and upper bound on x -> bound on y
 *   (3) right hand side and bounds on y -> upper bound on x
 *   (4) right hand side and lower bound on x -> bound on y
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

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(!SCIPisZero(scip, consdata->vbdcoef));

   switch( proprule )
   {
   case PROPRULE_1:
      /* lhs <= x + c*y: left hand side and bounds on y -> lower bound on x */
      assert(infervar == consdata->var);
      assert(boundtype == SCIP_BOUNDTYPE_LOWER);
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->vbdvar, bdchgidx) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->vbdvar, bdchgidx) );
      }
      break;

   case PROPRULE_2:
      /* lhs <= x + c*y: left hand side and upper bound on x -> bound on y */
      assert(infervar == consdata->vbdvar);
      assert(!SCIPisInfinity(scip, -consdata->lhs));
      SCIP_CALL( SCIPaddConflictUb(scip, consdata->var, bdchgidx) );
      break;

   case PROPRULE_3:
      /* x + c*y <= rhs: right hand side and bounds on y -> upper bound on x */
      assert(infervar == consdata->var);
      assert(boundtype == SCIP_BOUNDTYPE_UPPER);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->vbdvar, bdchgidx) );
      }
      else
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->vbdvar, bdchgidx) );
      }
      break;

   case PROPRULE_4:
      /* x + c*y <= rhs: right hand side and lower bound on x -> bound on y */
      assert(infervar == consdata->vbdvar);
      assert(!SCIPisInfinity(scip, consdata->rhs));
      SCIP_CALL( SCIPaddConflictLb(scip, consdata->var, bdchgidx) );
      break;

   case PROPRULE_INVALID:
   default:
      SCIPerrorMessage("invalid inference information %d in variable bound constraint <%s>\n", proprule, SCIPconsGetName(cons));
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


/** sets left hand side of varbound constraint */
static
SCIP_RETCODE chgLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, lhs));

   /* adjust value to not be smaller than -inf */
   if ( SCIPisInfinity(scip, -lhs) )
      lhs = -SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->var != NULL && consdata->vbdvar != NULL);
   assert(!SCIPisInfinity(scip, consdata->lhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->lhs, lhs) )
      return SCIP_OKAY;

   assert(consdata->row == NULL);

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */
   if( SCIPisEQ(scip, lhs, consdata->rhs) )
      consdata->rhs = lhs;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      /* the left hand side switched from -infinity to a non-infinite value -> install rounding locks */
      if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, -lhs) )
      {
	 SCIP_CALL( SCIPlockVarCons(scip, consdata->var, cons, TRUE, FALSE) );

	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
      }
      /* the left hand side switched from a non-infinite value to -infinity -> remove rounding locks */
      else if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, -lhs) )
      {
	 SCIP_CALL( SCIPunlockVarCons(scip, consdata->var, cons, TRUE, FALSE) );
	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
      }
   }

   /* if left hand side got tighter, we want to do additional presolving on this constraint */
   if( SCIPisLT(scip, consdata->lhs, lhs) )
   {
      consdata->propagated = FALSE;
      consdata->addvarbounds = FALSE;
      consdata->tightened = FALSE;
   }

   consdata->presolved = FALSE;
   consdata->lhs = lhs;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** sets right hand side of linear constraint */
static
SCIP_RETCODE chgRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linear constraint */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(!SCIPisInfinity(scip, -rhs));

   /* adjust value to not be larger than inf */
   if ( SCIPisInfinity(scip, rhs) )
      rhs = SCIPinfinity(scip);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->var != NULL && consdata->vbdvar != NULL);
   assert(!SCIPisInfinity(scip, -consdata->rhs));

   /* check whether the side is not changed */
   if( SCIPisEQ(scip, consdata->rhs, rhs) )
      return SCIP_OKAY;

   assert(consdata->row == NULL);

   /* ensure that rhs >= lhs is satisfied without numerical tolerance */
   if( SCIPisEQ(scip, rhs, consdata->lhs) )
      consdata->lhs = rhs;

   /* if necessary, update the rounding locks of variables */
   if( SCIPconsIsLocked(cons) )
   {
      assert(SCIPconsIsTransformed(cons));

      /* the right hand side switched from infinity to a non-infinite value -> install rounding locks */
      if( SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, rhs) )
      {
	 SCIP_CALL( SCIPlockVarCons(scip, consdata->var, cons, FALSE, TRUE) );

	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
      }
      /* the right hand side switched from a non-infinite value to infinity -> remove rounding locks */
      else if( !SCIPisInfinity(scip, consdata->rhs) && SCIPisInfinity(scip, rhs) )
      {
	 SCIP_CALL( SCIPunlockVarCons(scip, consdata->var, cons, FALSE, TRUE) );
	 if( SCIPisPositive(scip, consdata->vbdcoef) )
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, FALSE, TRUE) );
	 }
	 else
	 {
	    SCIP_CALL( SCIPunlockVarCons(scip, consdata->vbdvar, cons, TRUE, FALSE) );
	 }
      }
   }

   /* if right hand side got tighter, we want to do additional presolving on this constraint */
   if( SCIPisGT(scip, consdata->rhs, rhs) )
   {
      consdata->propagated = FALSE;
      consdata->addvarbounds = FALSE;
      consdata->tightened = FALSE;
   }

   consdata->presolved = FALSE;
   consdata->rhs = rhs;
   consdata->changed = TRUE;

   return SCIP_OKAY;
}

/** propagation method for variable bound constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the node can be cut off */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss           /**< pointer to count number of deleted constraints, or NULL */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xlb;
   SCIP_Real xub;
   SCIP_Real ylb;
   SCIP_Real yub;
   SCIP_Real newlb;
   SCIP_Real newub;
   SCIP_Bool tightened;
   SCIP_Bool tightenedround;

   assert(cutoff != NULL);
   assert(nchgbds != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIPdebugMessage("propagating variable bound constraint <%s>: %.15g <= <%s> + %.15g<%s> <= %.15g\n",
      SCIPconsGetName(cons), consdata->lhs, SCIPvarGetName(consdata->var), consdata->vbdcoef,
      SCIPvarGetName(consdata->vbdvar), consdata->rhs);

   *cutoff = FALSE;

   /* check, if constraint is already propagated */
   if( consdata->propagated )
      return SCIP_OKAY;

   /* get current bounds of variables */
   xlb = SCIPvarGetLbLocal(consdata->var);
   xub = SCIPvarGetUbLocal(consdata->var);
   ylb = SCIPvarGetLbLocal(consdata->vbdvar);
   yub = SCIPvarGetUbLocal(consdata->vbdvar);

   /* tighten bounds of variables as long as possible */
   do
   {
      tightenedround = FALSE;

      /* propagate left hand side inequality: lhs <= x + c*y */
      if( !SCIPisInfinity(scip, -consdata->lhs) )
      {
         assert(!(*cutoff));

         /* propagate bounds on x:
          *  (1) left hand side and bounds on y -> lower bound on x
          */
         if( SCIPvarGetStatus(consdata->var) != SCIP_VARSTATUS_MULTAGGR ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               if( !SCIPisInfinity(scip, yub) )
                  newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * yub);
               else
                  newlb = -SCIPinfinity(scip);
            }
            else
            {
               if( !SCIPisInfinity(scip, -ylb) )
                  newlb = SCIPadjustedVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * ylb);
               else
                  newlb = -SCIPinfinity(scip);
            }

            SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
               SCIPvarGetName(consdata->var), xlb, xub, newlb, xub);
            SCIP_CALL( SCIPinferVarLbCons(scip, consdata->var, newlb, cons, (int)PROPRULE_1, FALSE,
                  cutoff, &tightened) );

            if( *cutoff )
            {
               assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->var)));

               /* analyze infeasibility */
               SCIP_CALL( analyzeConflict(scip, cons, consdata->var, PROPRULE_1, SCIP_BOUNDTYPE_LOWER) );
               break;
            }

            if( tightened )
            {
               tightenedround = TRUE;
               (*nchgbds)++;
            }
            xlb = SCIPvarGetLbLocal(consdata->var);
         }

         assert(!*cutoff);

         /* propagate bounds on y:
          *  (2) left hand side and upper bound on x -> bound on y
          */
         if( SCIPvarGetStatus(consdata->vbdvar) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, xub) ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
               if( newlb > ylb + 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, (int)PROPRULE_2, FALSE,
                        cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->vbdvar)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_2, SCIP_BOUNDTYPE_LOWER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  ylb = SCIPvarGetLbLocal(consdata->vbdvar);
               }
            }
            else
            {
               newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->lhs - xub)/consdata->vbdcoef);
               if( newub < yub - 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, (int)PROPRULE_2, FALSE,
                        cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->vbdvar)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_2, SCIP_BOUNDTYPE_UPPER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  yub = SCIPvarGetUbLocal(consdata->vbdvar);
               }
            }
         }
      }

      assert(!*cutoff);

      /* propagate right hand side inequality: x + c*y <= rhs */
      if( !SCIPisInfinity(scip, consdata->rhs) )
      {
         /* propagate bounds on x:
          *  (3) right hand side and bounds on y -> upper bound on x
          */
         if( SCIPvarGetStatus(consdata->var) != SCIP_VARSTATUS_MULTAGGR ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               if( !SCIPisInfinity(scip, -ylb) )
                  newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * ylb);
               else
                  newub = SCIPinfinity(scip);
            }
            else
            {
               if( !SCIPisInfinity(scip, yub) )
                  newub = SCIPadjustedVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * yub);
               else
                  newub = SCIPinfinity(scip);
            }

            SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
               SCIPvarGetName(consdata->var), xlb, xub, xlb, newub);
            SCIP_CALL( SCIPinferVarUbCons(scip, consdata->var, newub, cons, (int)PROPRULE_3, FALSE,
                  cutoff, &tightened) );

            if( *cutoff )
            {
               assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->var)));

               /* analyze infeasibility */
               SCIP_CALL( analyzeConflict(scip, cons, consdata->var, PROPRULE_3, SCIP_BOUNDTYPE_UPPER) );
               break;
            }

            if( tightened )
            {
               tightenedround = TRUE;
               (*nchgbds)++;
            }
            xub = SCIPvarGetUbLocal(consdata->var);
         }

         assert(!*cutoff);

         /* propagate bounds on y:
          *  (4) right hand side and lower bound on x -> bound on y
          */
         if( SCIPvarGetStatus(consdata->vbdvar) != SCIP_VARSTATUS_MULTAGGR && !SCIPisInfinity(scip, -xlb) ) /* cannot change bounds of multaggr vars */
         {
            if( consdata->vbdcoef > 0.0 )
            {
               newub = SCIPadjustedVarUb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
               if( newub < yub - 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, ylb, newub);
                  SCIP_CALL( SCIPinferVarUbCons(scip, consdata->vbdvar, newub, cons, (int)PROPRULE_4, FALSE,
                        cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisLT(scip, newub, SCIPvarGetLbLocal(consdata->vbdvar)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_4, SCIP_BOUNDTYPE_UPPER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  yub = SCIPvarGetUbLocal(consdata->vbdvar);
               }
            }
            else
            {
               newlb = SCIPadjustedVarLb(scip, consdata->vbdvar, (consdata->rhs - xlb)/consdata->vbdcoef);
               if( newlb > ylb + 0.5 )
               {
                  SCIPdebugMessage(" -> tighten <%s>[%.15g,%.15g] -> [%.15g,%.15g]\n",
                     SCIPvarGetName(consdata->vbdvar), ylb, yub, newlb, yub);
                  SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vbdvar, newlb, cons, (int)PROPRULE_4, FALSE,
                        cutoff, &tightened) );

                  if( *cutoff )
                  {
                     assert(SCIPisGT(scip, newlb, SCIPvarGetUbLocal(consdata->vbdvar)));

                     /* analyze infeasibility */
                     SCIP_CALL( analyzeConflict(scip, cons, consdata->vbdvar, PROPRULE_4, SCIP_BOUNDTYPE_LOWER) );
                     break;
                  }

                  if( tightened )
                  {
                     tightenedround = TRUE;
                     (*nchgbds)++;
                  }
                  ylb = SCIPvarGetLbLocal(consdata->vbdvar);
               }
            }
         }
      }
      assert(!(*cutoff));
   }
   while( tightenedround );

   /* check for redundancy */
   if( !(*cutoff) && (SCIPisInfinity(scip, -consdata->lhs)
         || (consdata->vbdcoef > 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * ylb, consdata->lhs))
         || (consdata->vbdcoef < 0.0 && SCIPisFeasGE(scip, xlb + consdata->vbdcoef * yub, consdata->lhs)))
      && (SCIPisInfinity(scip, consdata->rhs)
         || (consdata->vbdcoef > 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * yub, consdata->rhs))
         || (consdata->vbdcoef < 0.0 && SCIPisFeasLE(scip, xub + consdata->vbdcoef * ylb, consdata->rhs))) )
   {
      SCIPdebugMessage("variable bound constraint <%s> is redundant: <%s>[%.15g,%.15g], <%s>[%.15g,%.15g]\n",
         SCIPconsGetName(cons),
         SCIPvarGetName(consdata->var), SCIPvarGetLbLocal(consdata->var), SCIPvarGetUbLocal(consdata->var),
         SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbLocal(consdata->vbdvar), SCIPvarGetUbLocal(consdata->vbdvar));
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      if( ndelconss != NULL )
         (*ndelconss)++;
   }

   /* mark the constraint propagated */
   consdata->propagated = TRUE;

   return SCIP_OKAY;
}

/** updates the flags of the first constraint according to the ones of the second constraint */
static
SCIP_RETCODE updateFlags(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons0,              /**< constraint that should stay */
   SCIP_CONS*            cons1               /**< constraint that should be deleted */
   )
{
   if( SCIPconsIsInitial(cons1) )
   {
      SCIP_CALL( SCIPsetConsInitial(scip, cons0, TRUE) );
   }
   if( SCIPconsIsSeparated(cons1) )
   {
      SCIP_CALL( SCIPsetConsSeparated(scip, cons0, TRUE) );
   }
   if( SCIPconsIsEnforced(cons1) )
   {
      SCIP_CALL( SCIPsetConsEnforced(scip, cons0, TRUE) );
   }
   if( SCIPconsIsChecked(cons1) )
   {
      SCIP_CALL( SCIPsetConsChecked(scip, cons0, TRUE) );
   }
   if( SCIPconsIsPropagated(cons1) )
   {
      SCIP_CALL( SCIPsetConsPropagated(scip, cons0, TRUE) );
   }
   if( !SCIPconsIsDynamic(cons1) )
   {
      SCIP_CALL( SCIPsetConsDynamic(scip, cons0, FALSE) );
   }
   if( !SCIPconsIsRemovable(cons1) )
   {
      SCIP_CALL( SCIPsetConsRemovable(scip, cons0, FALSE) );
   }
   if( SCIPconsIsStickingAtNode(cons1) )
   {
      SCIP_CALL( SCIPsetConsStickingAtNode(scip, cons0, TRUE) );
   }

   return SCIP_OKAY;
}

/* check whether one constraints side is redundant to another constraints side by calculating extreme values for
 * variables
 */
static
void checkRedundancySide(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             coef0,              /**< coefficient c0 of bounding variable y for constraint 0 */
   SCIP_Real             coef1,              /**< coefficient c1 of bounding variable y for constraint 1 */
   SCIP_Real             side0,              /**< one side of variable bound inequality for constraint 0 */
   SCIP_Real             side1,              /**< one side of variable bound inequality for constraint 1 */
   SCIP_Bool*            sideequal,          /**< pointer to store if both constraints have the same redundancy on the
                                              *   given side */
   SCIP_Bool*            cons0sidered,       /**< pointer to store if side of constraint 0 is redundant */
   SCIP_Bool*            cons1sidered,       /**< pointer to store if side of constraint 1 is redundant */
   SCIP_Bool             islhs               /**< do we check the left or the right hand side */
   )
{
   SCIP_Real lbvar;
   SCIP_Real ubvar;
   SCIP_Real lbvbdvar;
   SCIP_Real ubvbdvar;
   SCIP_Real boundxlb1;
   SCIP_Real boundxlb2;
   SCIP_Real boundylb1;
   SCIP_Real boundylb2;
   SCIP_Real boundxub1;
   SCIP_Real boundxub2;
   SCIP_Real boundyub1;
   SCIP_Real boundyub2;
   SCIP_Real boundvaluex1;
   SCIP_Real boundvaluex2;
   SCIP_Real boundvaluey1;
   SCIP_Real boundvaluey2;
   SCIP_Real valuex1;
   SCIP_Real valuex2;
   SCIP_Real valuey1;
   SCIP_Real valuey2;
   SCIP_Bool* redundant0;
   SCIP_Bool* redundant1;

   assert(scip != NULL);
   assert(var != NULL);
   assert(vbdvar != NULL);
   assert(sideequal != NULL);
   assert(cons0sidered != NULL);
   assert(cons1sidered != NULL);

   *cons0sidered = SCIPisInfinity(scip, REALABS(side0));
   *cons1sidered = SCIPisInfinity(scip, REALABS(side1));
   *sideequal = FALSE;

   if( islhs )
   {
      redundant0 = cons1sidered;
      redundant1 = cons0sidered;
   }
   else
   {
      redundant0 = cons0sidered;
      redundant1 = cons1sidered;
   }

   lbvar = SCIPvarGetLbGlobal(var);
   ubvar = SCIPvarGetUbGlobal(var);
   lbvbdvar = SCIPvarGetLbGlobal(vbdvar);
   ubvbdvar = SCIPvarGetUbGlobal(vbdvar);

   /* if both constraint have this side */
   if( !*redundant0 && !*redundant1 )
   {
      /* calculate extreme values, which are reached by setting the other variable to their lower/upper bound */
      boundxlb1 = side0 - lbvbdvar*coef0;
      boundxlb2 = side1 - lbvbdvar*coef1;
      boundylb1 = (side0 - lbvar)/coef0;
      boundylb2 = (side1 - lbvar)/coef1;

      boundxub1 = side0 - ubvbdvar*coef0;
      boundxub2 = side1 - ubvbdvar*coef1;
      boundyub1 = (side0 - ubvar)/coef0;
      boundyub2 = (side1 - ubvar)/coef1;

      if( islhs )
      {
	 boundvaluex1 = MAX(boundxlb1, boundxlb2);
	 boundvaluex2 = MAX(boundxub1, boundxub2);
      }
      else
      {
	 boundvaluex1 = MIN(boundxlb1, boundxlb2);
	 boundvaluex2 = MIN(boundxub1, boundxub2);
      }

      /* calculate important values for variables */
      if( SCIPisPositive(scip, coef0) )
      {
         valuex1 = MIN(boundvaluex1, ubvar);
         valuex1 = MAX(valuex1, lbvar);
         valuex2 = MAX(boundvaluex2, lbvar);
         valuex2 = MIN(valuex2, ubvar);

         /* if variable is of integral type make values integral too */
         if( SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS )
         {
            if( !SCIPisIntegral(scip, valuex1) )
               valuex1 = SCIPfloor(scip, valuex1);
            if( !SCIPisIntegral(scip, valuex2) )
               valuex2 = SCIPceil(scip, valuex2);
         }
      }
      else
      {
         valuex1 = MAX(boundvaluex1, lbvar);
         valuex1 = MIN(valuex1, ubvar);
         valuex2 = MIN(boundvaluex2, ubvar);
         valuex2 = MAX(valuex2, lbvar);

         /* if variable is of integral type make values integral too */
         if( SCIPvarGetType(var) < SCIP_VARTYPE_CONTINUOUS )
         {
            if( !SCIPisIntegral(scip, valuex1) )
               valuex1 = SCIPceil(scip, valuex1);
            if( !SCIPisIntegral(scip, valuex2) )
               valuex2 = SCIPfloor(scip, valuex2);
         }
      }

      /* calculate resulting values of variable y by setting x to valuex1 */
      valuey1 = (side0 - valuex1)/coef0;
      valuey2 = (side1 - valuex1)/coef1;

      /* determine redundancy of one constraints side */
      if( SCIPisPositive(scip, coef0) )
      {
         if( SCIPisEQ(scip, valuey1, valuey2) )
            *sideequal = TRUE;
         else if( SCIPisLT(scip, valuey1, valuey2) )
            *redundant1 = TRUE;
         else
            *redundant0 = TRUE;
      }
      else
      {
         if( SCIPisEQ(scip, valuey1, valuey2) )
            *sideequal = TRUE;
         else if( SCIPisLT(scip, valuey1, valuey2) )
            *redundant0 = TRUE;
         else
            *redundant1 = TRUE;
      }

      /* calculate resulting values of variable y by setting x to valuex2 */
      valuey1 = (side0 - valuex2)/coef0;
      valuey2 = (side1 - valuex2)/coef1;

      /* determine redundancy of one constraints side by checking for the first valuex2 */
      if( SCIPisPositive(scip, coef0) )
      {
         /* if both constraints are weaker than the other on one value, we have no redundancy */
         if( (*redundant1 && SCIPisGT(scip, valuey1, valuey2)) || (*redundant0 && SCIPisLT(scip, valuey1, valuey2)) )
         {
            *sideequal = FALSE;
            *redundant0 = FALSE;
            *redundant1 = FALSE;
            return;
         }
         else if( *sideequal )
         {
            if( SCIPisLT(scip, valuey1, valuey2) )
            {
               *sideequal = FALSE;
               *redundant1 = TRUE;
            }
            else if( SCIPisGT(scip, valuey1, valuey2) )
            {
               *sideequal = FALSE;
               *redundant0 = TRUE;
            }
         }
      }
      else
      {
         /* if both constraints are weaker than the other on one value, we have no redundancy */
         if( (*redundant1 && SCIPisLT(scip, valuey1, valuey2)) || (*redundant0 && SCIPisGT(scip, valuey1, valuey2)) )
         {
            *sideequal = FALSE;
            *redundant0 = FALSE;
            *redundant1 = FALSE;
            return;
         }
         else if( *sideequal )
         {
            if( SCIPisLT(scip, valuey1, valuey2) )
            {
               *sideequal = FALSE;
               *redundant0 = TRUE;
            }
            else if( SCIPisGT(scip, valuey1, valuey2) )
            {
               *sideequal = FALSE;
               *redundant1 = TRUE;
            }
         }
      }
      assert(*sideequal || *redundant0 || *redundant1);

      /* calculate feasability domain values for variable y concerning these both constraints */
      if( SCIPisPositive(scip, coef0) )
      {
	 if( islhs )
	 {
	    boundvaluey1 = MAX(boundylb1, boundylb2);
	    boundvaluey2 = MAX(boundyub1, boundyub2);
	 }
	 else
	 {
	    boundvaluey1 = MIN(boundylb1, boundylb2);
	    boundvaluey2 = MIN(boundyub1, boundyub2);
	 }

         valuey1 = MIN(boundvaluey1, ubvbdvar);
         valuey1 = MAX(valuey1, lbvbdvar);
         valuey2 = MAX(boundvaluey2, lbvbdvar);
         valuey2 = MIN(valuey2, ubvbdvar);

         if( !SCIPisIntegral(scip, valuey1) )
            valuey1 = SCIPfloor(scip, valuey1);
         if( !SCIPisIntegral(scip, valuey2) )
            valuey2 = SCIPceil(scip, valuey2);
      }
      else
      {
	 if( islhs )
	 {
	    boundvaluey1 = MIN(boundylb1, boundylb2);
	    boundvaluey2 = MIN(boundyub1, boundyub2);
	 }
	 else
	 {
	    boundvaluey1 = MAX(boundylb1, boundylb2);
	    boundvaluey2 = MAX(boundyub1, boundyub2);
	 }

         valuey1 = MAX(boundvaluey1, lbvbdvar);
         valuey1 = MIN(valuey1, ubvbdvar);
         valuey2 = MIN(boundvaluey2, ubvbdvar);
         valuey2 = MAX(valuey2, lbvbdvar);

         /* if variable is of integral type make values integral too */
         if( !SCIPisIntegral(scip, valuey1) )
            valuey1 = SCIPceil(scip, valuey1);
         if( !SCIPisIntegral(scip, valuey2) )
            valuey2 = SCIPfloor(scip, valuey2);
      }

      /* calculate resulting values of variable x by setting y to valuey1 */
      valuex1 = side0 - valuey1*coef0;
      valuex2 = side1 - valuey1*coef1;

      /* determine redundancy of one constraints side by checking for the first valuey1 */
      if( (*redundant1 && SCIPisGT(scip, valuex1, valuex2)) || (*redundant0 && SCIPisLT(scip, valuex1, valuex2)) )
      {
         *sideequal = FALSE;
         *redundant0 = FALSE;
         *redundant1 = FALSE;
         return;
      }
      if( *sideequal )
      {
         if( SCIPisLT(scip, valuex1, valuex2) )
         {
            *sideequal = FALSE;
            *redundant1 = TRUE;
         }
         else if( SCIPisGT(scip, valuex1, valuex2) )
         {
            *sideequal = FALSE;
            *redundant0 = TRUE;
         }
      }

      /* calculate resulting values of variable x by setting y to valuey2 */
      valuex1 = side0 - valuey2*coef0;
      valuex2 = side1 - valuey2*coef1;

      /* determine redundancy of one constraints side by checking for the first valuey1 */
      if( (*redundant1 && SCIPisGT(scip, valuex1, valuex2)) || (*redundant0 && SCIPisLT(scip, valuex1, valuex2)) )
      {
         *sideequal = FALSE;
         *redundant0 = FALSE;
         *redundant1 = FALSE;
         return;
      }
      if( *sideequal )
      {
         if( SCIPisLT(scip, valuex1, valuex2) )
         {
            *sideequal = FALSE;
            *redundant1 = TRUE;
         }
         else if( SCIPisGT(scip, valuex1, valuex2) )
         {
            *sideequal = FALSE;
            *redundant0 = TRUE;
         }
      }
      assert(*redundant0 || *redundant1 || *sideequal);
   }
}

/** compares each constraint with all other constraints for possible redundancy and removes or changes constraint
 *
 *  we will order all constraint to have constraints with same variables next to each other and fasten presolving
 *
 *  consider two constraints like lhs1 <= x + b1*y <= rhs1 and lhs2 <= x + b2*y <= rhs2
 *  we are doing the following presolving steps:
 *
 *  if( b1 == b2 )
 *      newlhs = MAX(lhs1, lhs2)
 *      newrhs = MIN(rhs1, rhs2)
 *      updateSides
 *      delete one constraint
 *  else if( ((b1 > 0) == (b2 > 0)) && (lhs1 != -inf && lhs2 != -inf) || (rhs1 != inf && rhs2 != inf) )
 *
 *       (i.e. both constraint have either a valid lhs or a valid rhs and infinity is on the same side and the
 *             coeffcients have the same size )
 *
 *      if( y is binary )
 *          if( lhs1 != -inf )
 *              newlhs = MAX(lhs1, lhs2)
 *              newb = newlhs - MAX(lhs1 - b1, lhs2 - b2)
 *          else
 *              newrhs = MIN(lhs1, lhs2)
 *              newb = newrhs - MIN(rhs1 - b1, rhs2 - b2)
 *          updateSidesAndCoef
 *          delete one constraint
 *      else
 *          we calculate possible values for both variables and check which constraint is tighter
 *  else
 *      nothing possible
 */
static
SCIP_RETCODE preprocessConstraintPairs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   SCIP_Bool*            cutoff,             /**< pointer to store TRUE, if a cutoff was found */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_CONS** sortedconss;
   int c;
   int s;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   /* create our temporary working array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &sortedconss, conss, nconss) );

   /* sort all constraints, so that all constraints with same variables stand next to each other */
   SCIPsortPtr((void**)sortedconss, consVarboundComp, nconss);

   /* check all constraints for redundancy */
   for( c = nconss - 1; c > 0 && !(*cutoff); --c )
   {
      SCIP_CONS* cons0;
      SCIP_CONSDATA* consdata0;

      cons0 = sortedconss[c];

      if( !SCIPconsIsActive(cons0) || SCIPconsIsModifiable(cons0) )
         continue;

      consdata0 = SCIPconsGetData(cons0);
      assert(consdata0 != NULL);
      assert(consdata0->var != NULL);
      assert(consdata0->vbdvar != NULL);

      /* do not check for already redundant constraints */
      assert(!SCIPisZero(scip, consdata0->vbdcoef));
      assert(!SCIPisInfinity(scip, -consdata0->lhs) || !SCIPisInfinity(scip, consdata0->rhs));

      if( !consdata0->changed )
         continue;

      consdata0->changed = FALSE;

      for( s = c - 1; s >= 0; --s )
      {
         SCIP_CONS* cons1;
         SCIP_CONSDATA* consdata1;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real coef;
         SCIP_Bool deletecons1;

         cons1 = sortedconss[s];

         if( !SCIPconsIsActive(cons1) || SCIPconsIsModifiable(cons1) )
            continue;

         consdata1 = SCIPconsGetData(cons1);
         assert(consdata1 != NULL);
         assert(consdata1->var != NULL);
         assert(consdata1->vbdvar != NULL);

         /* do not check for already redundant constraints */
         assert(!SCIPisZero(scip, consdata1->vbdcoef));
         assert(!SCIPisInfinity(scip, -consdata1->lhs) || !SCIPisInfinity(scip, consdata1->rhs));

         /* check for equal variables */
         if( consdata0->var != consdata1->var || consdata0->vbdvar != consdata1->vbdvar )
            break;

         /* mark constraint1 for deletion if possible */
         deletecons1 = TRUE;

         lhs = consdata0->lhs;
         rhs = consdata0->rhs;
         coef = consdata0->vbdcoef;

         /* the coefficients of both constraints are equal */
         if( SCIPisEQ(scip, coef, consdata1->vbdcoef) )
         {
            lhs = MAX(consdata1->lhs, lhs);
            rhs = MIN(consdata1->rhs, rhs);
         }
         /* now only one side and in both constraints the same side should be infinity and the vbdvar should be binary
          * then we neither do not need to have the same side nor the same coefficient
          */
         else if( SCIPvarIsBinary(consdata0->vbdvar)
            && (SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs))
            && (SCIPisInfinity(scip, -consdata1->lhs) || SCIPisInfinity(scip, consdata1->rhs))
            && (SCIPisInfinity(scip, -lhs) == SCIPisInfinity(scip, -consdata1->lhs)) )
         {
            /* lhs <= x + b*y <= +inf */
            if( !SCIPisInfinity(scip, -lhs) )
            {
               lhs = MAX(consdata1->lhs, lhs);
               coef = lhs - MAX(consdata1->lhs - consdata1->vbdcoef, consdata0->lhs - coef);
            }
            /* -inf <= x + b*y <= rhs */
            else
            {
               rhs = MIN(consdata1->rhs, rhs);
               coef = rhs - MIN(consdata1->rhs - consdata1->vbdcoef, consdata0->rhs - coef);
            }
	    consdata0->propagated = FALSE;
         }
         else if( SCIPisPositive(scip, coef) == SCIPisPositive(scip, consdata1->vbdcoef)
            && ((!SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, -consdata1->lhs))
               || (!SCIPisInfinity(scip, rhs) && !SCIPisInfinity(scip, consdata1->rhs))) )
         {
            SCIP_Bool cons0lhsred;
            SCIP_Bool cons0rhsred;
            SCIP_Bool cons1lhsred;
            SCIP_Bool cons1rhsred;
            SCIP_Bool lhsequal;
            SCIP_Bool rhsequal;

            assert(!SCIPisInfinity(scip, lhs));
            assert(!SCIPisInfinity(scip, consdata1->lhs));
            assert(!SCIPisInfinity(scip, -rhs));
            assert(!SCIPisInfinity(scip, -consdata1->rhs));

            /* check if a left hand side of one constraints is redundant */
            checkRedundancySide(scip, consdata0->var, consdata0->vbdvar, coef, consdata1->vbdcoef, lhs, consdata1->lhs, &lhsequal, &cons0lhsred, &cons1lhsred, TRUE);

            /* check if a right hand side of one constraints is redundant */
            checkRedundancySide(scip, consdata0->var, consdata0->vbdvar, coef, consdata1->vbdcoef, rhs, consdata1->rhs, &rhsequal, &cons0rhsred, &cons1rhsred, FALSE);

            /* if cons0 is redundant, update cons1 and delete cons0 */
            if( (lhsequal || cons0lhsred) && (rhsequal || cons0rhsred) )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( updateFlags(scip, cons1, cons0) );

               SCIPdebugMessage("constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
               SCIPdebugMessage("is redundant to constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );

               SCIP_CALL( SCIPdelCons(scip, cons0) );
               ++(*ndelconss);

               /* get next cons0 */
               break;
            }
            /* if cons1 is redundant, update cons0 and delete cons1 */
            else if( cons1lhsred && cons1rhsred )
            {
               /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
               SCIP_CALL( updateFlags(scip, cons0, cons1) );

               SCIPdebugMessage("constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
               SCIPdebugMessage("is redundant to constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );

               SCIP_CALL( SCIPdelCons(scip, cons1) );
               ++(*ndelconss);

               /* get next cons1 */
               continue;
            }
            /* if left hand side of cons0 is redundant set it to -infinity */
            else if( (lhsequal || cons0lhsred) && !SCIPisInfinity(scip, -lhs) )
            {
               lhs = -SCIPinfinity(scip);

	       /* if right hand side of cons1 is redundant too, set it to infinity */
	       if( cons1rhsred && !SCIPisInfinity(scip, consdata1->rhs) )
	       {
		  SCIP_CALL( chgRhs(scip, cons1, SCIPinfinity(scip)) );
		  ++(*nchgsides);

		  SCIPdebugMessage("deleted rhs of constraint: ");
		  SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
		  SCIPdebugMessage("due to constraint: ");
		  SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
	       }

               /* later on we cannot not want to delete cons1 */
               deletecons1 = FALSE;
            }
            /* if right hand side of cons0 is redundant set it to infinity */
            else if( (rhsequal || cons0rhsred) && !SCIPisInfinity(scip, rhs) )
            {
               rhs = SCIPinfinity(scip);

	       /* if left hand side of cons1 is redundant too, set it to -infinity */
	       if( cons1lhsred && !SCIPisInfinity(scip, -consdata1->lhs) )
	       {
		  SCIP_CALL( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
		  ++(*nchgsides);

		  SCIPdebugMessage("deleted lhs of constraint: ");
		  SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
		  SCIPdebugMessage("due to constraint: ");
		  SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
	       }

               /* later on we cannot not want to delete cons1 */
               deletecons1 = FALSE;
            }
            /* if left hand side of cons1 is redundant set it to -infinity */
            else if( cons1lhsred && !SCIPisInfinity(scip, -consdata1->lhs) )
            {
	       SCIP_CALL( chgLhs(scip, cons1, -SCIPinfinity(scip)) );
	       ++(*nchgsides);

               SCIPdebugMessage("deleted lhs of constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
               SCIPdebugMessage("due to constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );

               continue;
            }
            /* if right hand side of cons1 is redundant set it to infinity */
            else if( cons1rhsred && !SCIPisInfinity(scip, consdata1->rhs) )
            {
	       SCIP_CALL( chgRhs(scip, cons1, SCIPinfinity(scip)) );
	       ++(*nchgsides);

               SCIPdebugMessage("deleted rhs of constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
               SCIPdebugMessage("due to constraint: ");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );

               continue;
            }
            else /* nothing was redundant */
               continue;
         }
         else
         {
            /* there is no redundancy in both constraints with same variables */
            continue;
         }

         if( SCIPisFeasLT(scip, rhs, lhs) )
         {
            SCIPdebugMessage("constraint <%s> and <%s> lead to infeasibility due to their sides\n", SCIPconsGetName(cons0), SCIPconsGetName(cons1));
            *cutoff = TRUE;
            break;
         }

         /* ensure that lhs <= rhs holds without tolerances as we only allow such rows to enter the LP */
         if( lhs > rhs )
         {
            rhs = (lhs + rhs)/2;
            lhs = rhs;
         }

         /* we decide to let constraint cons0 stay, so update data structure consdata0 */

         /* update coefficient of cons0 */

         /* special case if new coefficient becomes zero, both constraints are redundant but we may tighten the bounds */
         if( SCIPisZero(scip, coef) )
         {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;

            SCIPdebugMessage("constraint: ");
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
            SCIPdebugMessage("and constraint: ");
            SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );
            SCIPdebugMessage("are both redundant and lead to bounding of <%s> in [%g, %g]\n", SCIPvarGetName(consdata0->var), lhs, rhs);

            /* delete cons1 */
            SCIP_CALL( SCIPdelCons(scip, cons1) );
            ++(*ndelconss);

            /* update upper bound if possible */
            SCIP_CALL( SCIPtightenVarUb(scip, consdata0->var, rhs, FALSE, &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               ++(*nchgbds);

            /* update lower bound if possible */
            SCIP_CALL( SCIPtightenVarLb(scip, consdata0->var, lhs, FALSE, &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened )
               ++(*nchgbds);

            /* delete cons0 */
            SCIP_CALL( SCIPdelCons(scip, cons0) );
            ++(*ndelconss);

            /* get next cons0 */
            break;
         }

         SCIPdebugMessage("constraint: ");
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons1, NULL) ) );
         SCIPdebugMessage("and constraint: ");
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );

         /* if sign of coefficient switches, update the rounding locks of the variable */
         if( SCIPconsIsLocked(cons0) && consdata0->vbdcoef * coef < 0.0 )
         {
            assert(SCIPconsIsTransformed(cons0));

            /* remove rounding locks for variable with old coefficient and install rounding locks for variable with new
             * coefficient
             */
            if( SCIPisPositive(scip, consdata0->vbdcoef) )
            {
               SCIP_CALL( SCIPunlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, -consdata0->lhs), !SCIPisInfinity(scip, consdata0->rhs)) );
               SCIP_CALL( SCIPlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, consdata0->rhs), !SCIPisInfinity(scip, -consdata0->lhs)) );
            }
            else
            {
               SCIP_CALL( SCIPunlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, consdata0->rhs), !SCIPisInfinity(scip, -consdata0->lhs)) );
               SCIP_CALL( SCIPlockVarCons(scip, consdata0->vbdvar, cons0,
                     !SCIPisInfinity(scip, -consdata0->lhs), !SCIPisInfinity(scip, consdata0->rhs)) );
            }
         }

         /* now change the coefficient */
         if( !SCIPisEQ(scip, consdata0->vbdcoef, coef) )
         {
            ++(*nchgcoefs);

            /* mark to add new varbound information */
            consdata0->addvarbounds = FALSE;
	    consdata0->tightened = FALSE;
	    consdata0->propagated = FALSE;
	    consdata0->presolved = FALSE;
	    consdata0->changed = FALSE;

	    consdata0->vbdcoef = coef;
         }

         /* update lhs and rhs of cons0 */
         if( !SCIPisEQ(scip, consdata0->lhs, lhs) )
         {
	    SCIP_CALL( chgLhs(scip, cons0, lhs) );
	    ++(*nchgsides);
	 }
         if( !SCIPisEQ(scip, consdata0->rhs, rhs) )
         {
	    SCIP_CALL( chgRhs(scip, cons0, rhs) );
	    ++(*nchgsides);
	 }

         /* update flags of constraint which caused the redundancy s.t. nonredundant information doesn't get lost */
         SCIP_CALL( updateFlags(scip, cons0, cons1) );

         SCIPdebugMessage("lead to new constraint: ");
         SCIPdebug( SCIP_CALL( SCIPprintCons(scip, cons0, NULL) ) );

	 /* if cons1 is still marked for deletion, delete it */
         if( deletecons1 )
         {
	    /* delete cons1 */
	    SCIP_CALL( SCIPdelCons(scip, cons1) );
	    ++(*ndelconss);
	 }

         assert(SCIPconsIsActive(cons0));
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &sortedconss);

   return SCIP_OKAY;
}

/* for varbound constraint with two integer variables make coefficients integral */
static
SCIP_RETCODE prettifyConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraint set */
   int                   nconss,             /**< number of constraints in constraint set */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides           /**< pointer to count number of changed left/right hand sides */
   )
{
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);

   /* if we cannot find any constraint for prettifying, stop */
   if( SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip) < 2 )
      return SCIP_OKAY;

   for( c = nconss - 1; c >= 0; --c )
   {
      assert(conss != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* check for integer variables and one coefficient with an absolute value smaller than 1 */
      /* @note: we allow that the variable type of the bounded variable can be smaller than the variable type of the
       *        bounding variable
       */
      if( (SCIPvarGetType(consdata->var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(consdata->var) == SCIP_VARTYPE_INTEGER
	    || SCIPvarGetType(consdata->var) == SCIP_VARTYPE_IMPLINT)
	 && (SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_IMPLINT)
	 && SCIPisLT(scip, REALABS(consdata->vbdcoef), 1.0) )
      {
         SCIP_Real epsilon;
         SCIP_Longint nominator;
         SCIP_Longint denominator;
         SCIP_Longint maxmult;
         SCIP_Bool success;

         epsilon = SCIPepsilon(scip) * 0.9;  /* slightly decrease epsilon to be safe in rational conversion below */
         maxmult = (SCIP_Longint)(SCIPfeastol(scip)/epsilon + SCIPfeastol(scip));
         maxmult = MIN(maxmult, MAXSCALEDCOEF);

         success = SCIPrealToRational(consdata->vbdcoef, -epsilon, epsilon , maxmult, &nominator, &denominator);

         if( success )
         {
            /* it is possible that the dominator is a multiple of the nominator */
            if( SCIPisIntegral(scip, denominator / (SCIP_Real)nominator) )
            {
               denominator /= nominator;
               nominator = 1;
            }

            success = success && (denominator <= maxmult);

            /* scale the constraint denominator/nominator */
            if( success && ABS(denominator) > 1 && nominator == 1)
            {
               SCIP_VAR* swapvar;

               /* print constraint before scaling */
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );

               assert(SCIPisEQ(scip, consdata->vbdcoef * denominator, 1.0));

               /* need to switch sides if coefficient is smaller then 0 */
               if( consdata->vbdcoef < 0 )
               {
                  assert(denominator < 0);

                  /* compute new sides */

                  /* only right hand side exists */
                  if( SCIPisInfinity(scip, -consdata->lhs) )
                  {
                     consdata->lhs = consdata->rhs * denominator;
                     assert(!SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->lhs));

                     consdata->rhs = SCIPinfinity(scip);
                  }
                  /* only left hand side exists */
                  else if( SCIPisInfinity(scip, consdata->rhs) )
                  {
                     consdata->rhs = consdata->lhs * denominator;
                     assert(!SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->rhs));

                     consdata->lhs = -SCIPinfinity(scip);
                  }
                  /* both sides exist */
                  else
                  {
                     SCIP_Real tmp;

                     tmp = consdata->lhs;
                     consdata->lhs = consdata->rhs * denominator;
                     consdata->rhs = tmp * denominator;
		     consdata->tightened = FALSE;

                     assert(!SCIPisInfinity(scip, consdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs));
                     assert(SCIPisGE(scip, consdata->rhs, consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs));
                  }
                  *nchgsides += 2;
               }
               /* coefficient > 0 */
               else
               {
                  assert(denominator > 0);

                  /* compute new left hand side */
                  if( !SCIPisInfinity(scip, -consdata->lhs) )
                  {
                     consdata->lhs *= denominator;
                     assert(!SCIPisInfinity(scip, consdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs));
                     ++(*nchgsides);
                  }

                  /* compute new right hand side */
                  if( !SCIPisInfinity(scip, consdata->rhs) )
                  {
                     consdata->rhs *= denominator;
                     assert(!SCIPisInfinity(scip, consdata->rhs) && !SCIPisInfinity(scip, -consdata->rhs));
                     ++(*nchgsides);
                  }

                  assert(SCIPisGE(scip, consdata->rhs, consdata->lhs));
               }

               /* swap both variables */
               swapvar = consdata->var;
               consdata->var = consdata->vbdvar;
               consdata->vbdvar = swapvar;

               /* swap coefficient */
               consdata->vbdcoef = (SCIP_Real)denominator;
               ++(*nchgcoefs);

               /* mark to add new varbound information */
               consdata->addvarbounds = FALSE;
	       consdata->tightened = FALSE;

               /* print constraint after scaling */
               SCIPdebugMessage("transformed into:");
               SCIPdebug( SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) ) );
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** replaces fixed and aggregated variables in variable bound constraint by active problem variables */
static
SCIP_RETCODE applyFixings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   SCIP_Bool*            cutoff,             /**< pointer to store whether an infeasibility was detected */
   int*                  nchgbds,            /**< pointer to count number of bound changes */
   int*                  ndelconss,          /**< pointer to count number of deleted constraints */
   int*                  naddconss           /**< pointer to count number of added constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   SCIP_Real varscalar;
   SCIP_Real varconstant;
   SCIP_VAR* vbdvar;
   SCIP_Real vbdvarscalar;
   SCIP_Real vbdvarconstant;
   SCIP_Bool varschanged;
   SCIP_Bool redundant;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(cutoff != NULL);
   assert(nchgbds != NULL);
   assert(ndelconss != NULL);
   assert(naddconss != NULL);

   *cutoff = FALSE;
   redundant = FALSE;

   /* the variable bound constraint is: lhs <= x + c*y <= rhs */
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* get active problem variables of x and y */
   var = consdata->var;
   varscalar = 1.0;
   varconstant = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&var, &varscalar, &varconstant) );
   vbdvar = consdata->vbdvar;
   vbdvarscalar = 1.0;
   vbdvarconstant = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&vbdvar, &vbdvarscalar, &vbdvarconstant) );
   varschanged = (var != consdata->var || vbdvar != consdata->vbdvar);

   /* if the variables are equal, the variable bound constraint reduces to standard bounds on the single variable */
   if( var == vbdvar && SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      SCIPdebugMessage("variable bound constraint <%s> has equal variable and vbd variable <%s>\n",
         SCIPconsGetName(cons), SCIPvarGetName(var));

      /*      lhs <= a1*z + b1 + c(a2*z + b2) <= rhs
       * <=>  lhs <= (a1 + c*a2)z + (b1 + c*b2) <= rhs
       */
      scalar = varscalar + consdata->vbdcoef * vbdvarscalar;
      constant = varconstant + consdata->vbdcoef * vbdvarconstant;
      if( SCIPisZero(scip, scalar) )
      {
         /* no variable is left: the constraint is redundant or infeasible */
         if( SCIPisFeasLT(scip, constant, consdata->lhs) || SCIPisFeasGT(scip, constant, consdata->rhs) )
            *cutoff = TRUE;
      }
      else if( scalar > 0.0 )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarLb(scip, var, (consdata->lhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
               (*nchgbds)++;
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarUb(scip, var, (consdata->rhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetUbGlobal(var));
               (*nchgbds)++;
            }
         }
      }
      else
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarUb(scip, var, (consdata->lhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetUbGlobal(var));
               (*nchgbds)++;
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
         {
            SCIP_Bool tightened;
            
            SCIP_CALL( SCIPtightenVarLb(scip, var, (consdata->rhs - constant)/scalar, TRUE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                  SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
               (*nchgbds)++;
            }
         }
      }
      redundant = TRUE;
   }
   else
   {
      /* if the variables should be replaced, drop the events and catch the events on the new variables afterwards */
      if( varschanged )
      {
         SCIP_CALL( dropEvents(scip, consdata) );
      }

      /* apply aggregation on x */
      if( SCIPisZero(scip, varscalar) )
      {
         SCIPdebugMessage("variable bound constraint <%s>: variable <%s> is fixed to %.15g\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->var), varconstant);

         /* cannot change bounds on multi-aggregated variables */
         if( SCIPvarGetStatus(vbdvar) != SCIP_VARSTATUS_MULTAGGR )
         {
            /* x is fixed to varconstant: update bounds of y and delete the variable bound constraint */
            if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarLb(scip, consdata->vbdvar, (consdata->lhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
               else
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vbdvar, (consdata->lhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetUbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarUb(scip, consdata->vbdvar, (consdata->rhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetUbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
               else
               {
                  SCIP_Bool tightened;

                  SCIP_CALL( SCIPtightenVarLb(scip, consdata->vbdvar, (consdata->rhs - varconstant)/consdata->vbdcoef,
                        TRUE, cutoff, &tightened) );
                  if( tightened )
                  {
                     SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                        SCIPvarGetName(consdata->vbdvar), SCIPvarGetLbGlobal(consdata->vbdvar));
                     (*nchgbds)++;
                  }
               }
            }
            redundant = TRUE;
         }
      }
      else if( var != consdata->var )
      {
         /* replace aggregated variable x in the constraint by its aggregation */
         if( varscalar > 0.0 )
         {
            /* lhs := (lhs - varconstant) / varscalar
             * rhs := (rhs - varconstant) / varscalar
             * c   := c / varscalar
             */
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs = (consdata->lhs - varconstant)/varscalar;
            if( !SCIPisInfinity(scip, consdata->rhs) )
               consdata->rhs = (consdata->rhs - varconstant)/varscalar;
            consdata->vbdcoef /= varscalar;

	    consdata->tightened = FALSE;
         }
         else
         {
            SCIP_Real lhs;
            
            assert(varscalar != 0.0);

            /* lhs := (rhs - varconstant) / varscalar
             * rhs := (lhs - varconstant) / varscalar
             * c   := c / varscalar
             */
            lhs = consdata->lhs;
            consdata->lhs = -consdata->rhs;
            consdata->rhs = -lhs;
            if( !SCIPisInfinity(scip, -consdata->lhs) )
               consdata->lhs = (consdata->lhs + varconstant)/(-varscalar);
            if( !SCIPisInfinity(scip, consdata->rhs) )
               consdata->rhs = (consdata->rhs + varconstant)/(-varscalar);
            consdata->vbdcoef /= varscalar;

	    consdata->tightened = FALSE;
         }
         /* release old variable */
         SCIP_CALL( SCIPreleaseVar(scip, &(consdata->var)) );
         consdata->var = var;
         /* capture new variable */
         SCIP_CALL( SCIPcaptureVar(scip, consdata->var) );
      }

      /* apply aggregation on y */
      if( SCIPisZero(scip, vbdvarscalar) )
      {
         SCIPdebugMessage("variable bound constraint <%s>: vbd variable <%s> is fixed to %.15g\n",
            SCIPconsGetName(cons), SCIPvarGetName(consdata->vbdvar), vbdvarconstant);

         /* cannot change bounds on multi-aggregated variables */
         if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR )
         {
            /* y is fixed to vbdvarconstant: update bounds of x and delete the variable bound constraint */
            if( !SCIPisInfinity(scip, -consdata->lhs) && !(*cutoff) )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarLb(scip, consdata->var, consdata->lhs - consdata->vbdcoef * vbdvarconstant,
                     TRUE, cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMessage(" -> tightened lower bound: <%s> >= %.15g\n", 
                     SCIPvarGetName(consdata->var), SCIPvarGetLbGlobal(consdata->var));
                  (*nchgbds)++;
               }
            }
            if( !SCIPisInfinity(scip, consdata->rhs) && !(*cutoff) )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarUb(scip, consdata->var, consdata->rhs - consdata->vbdcoef * vbdvarconstant,
                     TRUE, cutoff, &tightened) );
               if( tightened )
               {
                  SCIPdebugMessage(" -> tightened upper bound: <%s> <= %.15g\n", 
                     SCIPvarGetName(consdata->var), SCIPvarGetUbGlobal(consdata->var));
                  (*nchgbds)++;
               }
            }
            redundant = TRUE;
         }
      }
      else if( vbdvar != consdata->vbdvar )
      {
         /* replace aggregated variable y in the constraint by its aggregation:
          * lhs := lhs - c * vbdvarconstant
          * rhs := rhs - c * vbdvarconstant
          * c   := c * vbdvarscalar
          */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
            consdata->lhs -= consdata->vbdcoef * vbdvarconstant;
         if( !SCIPisInfinity(scip, consdata->rhs) )
            consdata->rhs -= consdata->vbdcoef * vbdvarconstant;
         consdata->vbdcoef *= vbdvarscalar;

         /* release old variable */
         SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vbdvar)) );
         consdata->vbdvar = vbdvar;
         /* capture new variable */
         SCIP_CALL( SCIPcaptureVar(scip, consdata->vbdvar) );
      }

      /* catch the events again on the new variables */
      if( varschanged )
      {
         SCIP_CALL( catchEvents(scip, consdata) );
      }
   }

   /* mark constraint changed, if a variable was exchanged */
   if( varschanged )
   {
      consdata->changed = TRUE;
   }

   /* active multi aggregations are now resolved by creating a new linear constraint */
   if( !(*cutoff) && !redundant && (SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(vbdvar) == SCIP_VARSTATUS_MULTAGGR) )
   {
      SCIP_CONS* newcons;
      SCIP_Real lhs;
      SCIP_Real rhs;

      lhs = consdata->lhs;
      rhs = consdata->rhs;

      assert(var == consdata->var);
      assert(vbdvar == consdata->vbdvar);

      /* create upgraded linear constraint */
      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, SCIPconsGetName(cons), 0, NULL, NULL, lhs, rhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons),
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );

      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, consdata->var, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, newcons, consdata->vbdvar, consdata->vbdcoef) );

      SCIP_CALL( SCIPaddCons(scip, newcons) );

      SCIPdebugMessage("resolved multi aggregation in varbound constraint <%s> by creating a new linear constraint\n", SCIPconsGetName(cons));
      SCIPdebug( SCIP_CALL( SCIPprintCons(scip, newcons) ) );

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );

      redundant = TRUE;
      ++(*naddconss);
   }

   /* delete a redundant constraint */
   if( !(*cutoff) && redundant )
   {
      SCIPdebugMessage(" -> variable bound constraint <%s> is redundant\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      (*ndelconss)++;
   }

   return SCIP_OKAY;
}

/** tightens variable bound coefficient by inspecting the global bounds of the involved variables; note: this is also
 *  performed by the linear constraint handler - only necessary if the user directly creates variable bound constraints
 */
static
SCIP_RETCODE tightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< variable bound constraint */
   int*                  nchgcoefs,          /**< pointer to count the number of changed coefficients */
   int*                  nchgsides,          /**< pointer to count the number of left and right hand sides */
   int*                  ndelconss           /**< pointer to count number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real xlb;
   SCIP_Real xub;
   int oldnchgcoefs;
   int oldnchgsides;

   assert(nchgcoefs != NULL);
   assert(nchgsides != NULL);
   assert(ndelconss != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* tightening already done */
   if( consdata->tightened )
      return SCIP_OKAY;

   consdata->tightened = TRUE;

   /* if values and variable are integral the sides should it be too */
   if( SCIPvarGetType(consdata->var) <= SCIP_VARTYPE_IMPLINT
      && SCIPvarGetType(consdata->vbdvar) <= SCIP_VARTYPE_IMPLINT
      && SCIPisIntegral(scip, consdata->vbdcoef) )
   {
      if( !SCIPisIntegral(scip, consdata->lhs) )
      {
         consdata->lhs = SCIPceil(scip, consdata->lhs);
         ++(*nchgsides);
         consdata->changed = TRUE;
      }
      if( !SCIPisIntegral(scip, consdata->rhs) )
      {
         consdata->rhs = SCIPfloor(scip, consdata->rhs);
         ++(*nchgsides);
         consdata->changed = TRUE;
      }
   }

   /* coefficient tightening only works for binary bound variable */
   if( !SCIPvarIsBinary(consdata->vbdvar) )
      return SCIP_OKAY;

   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   /* coefficients tightening when all variables are integer */
   /* we consider the following varbound constraint: lhs <= x + b*y <= rhs (sides are possibly infinity)
    * y should always be binary and x of integral type and b not integral, we also need at least one side with infinity
    * or not integral value.
    *
    * 1. if( (lhs is integral and not -infinity) and ((rhs is infinity) or (b - floor(b) <= rhs - floor(rhs))) ):
    *
    *        lhs <= x + b*y <= rhs =>   lhs <= x + floor(b)*y <= floor(rhs)
    *
    * 2. if( (rhs is integral and not infinity) and ((lhs is -infinity) or (b - floor(b) >= lhs - floor(lhs))) ):
    *
    *        lhs <= x + b*y <= rhs   =>  ceil(lhs) <= x + ceil(b)*y <= rhs
    *
    * 3. if( ((lhs is -infinity) or (b - floor(b) >= lhs - floor(lhs)))
    *       and ((rhs is infinity) or (b - floor(b) > rhs - floor(rhs))) ):
    *
    *        lhs <= x + b*y <= rhs  =>   ceil(lhs) <= x + ceil(b)*y <= floor(rhs)
    *
    * 4. if( ((lhs is -infinity) or (b - floor(b) < lhs - floor(lhs)))
    *       and ((rhs is infinity) or (b - floor(b) <= rhs - floor(rhs))) ):
    *
    *        lhs <= x + b*y <= rhs  =>   ceil(lhs) <= x + floor(b)*y <= floor(rhs)
    *
    * 5. if( (lhs is not integral) or (rhs is not integral) )
    *
    *       if (lhs is not -infinity)
    *          if (b - floor(b) < lhs - floor(lhs)):
    *
    *             lhs <= x + b*y  =>   ceil(lhs) <= x + b*y
    *
    *          else if (b - floor(b) > lhs - floor(lhs)):
    *
    *             lhs <= x + b*y  =>   floor(lhs) + b - floor(b) <= x + b*y
    *
    *       if (rhs is not infinity)
    *          if (b - floor(b) < rhs - floor(rhs)):
    *
    *             x + b*y <= rhs  =>   x + b*y <= floor(rhs) + b - floor(b)
    *
    *          else if (b - floor(b) > rhs - floor(rhs)):
    *
    *             x + b*y <= rhs  =>   x + b*y <= floor(rhs)
    */
   if( (SCIPvarGetType(consdata->var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(consdata->var) == SCIP_VARTYPE_IMPLINT)
      && !SCIPisIntegral(scip, consdata->vbdcoef)
      && (!SCIPisIntegral(scip, consdata->lhs) || SCIPisInfinity(scip, -consdata->lhs)
	 || !SCIPisIntegral(scip, consdata->rhs) || SCIPisInfinity(scip, consdata->rhs)) )
   {
      /* infinity should be an integral value */
      assert(!SCIPisInfinity(scip, -consdata->lhs) || SCIPisIntegral(scip, consdata->lhs));
      assert(!SCIPisInfinity(scip, consdata->rhs) || SCIPisIntegral(scip, consdata->rhs));

      /* should not be a redundant constraint */
      assert(!SCIPisInfinity(scip, consdata->rhs) || !SCIPisInfinity(scip, -consdata->lhs));

      /* case 1 */
      if( SCIPisIntegral(scip, consdata->lhs) && !SCIPisInfinity(scip, -consdata->lhs) &&
         (SCIPisInfinity(scip, consdata->rhs) || SCIPisLE(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfloor(scip, consdata->rhs))) )
      {
         consdata->vbdcoef = SCIPfloor(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisIntegral(scip, consdata->rhs) )
         {
            consdata->rhs = SCIPfloor(scip, consdata->rhs);
            ++(*nchgsides);
         }
      }
      /* case 2 */
      else if( SCIPisIntegral(scip, consdata->rhs) && !SCIPisInfinity(scip, consdata->rhs) &&
         (SCIPisInfinity(scip, -consdata->lhs) || SCIPisGE(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfloor(scip, consdata->lhs))) )

      {
         consdata->vbdcoef = SCIPceil(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisIntegral(scip, consdata->lhs) )
         {
            consdata->lhs = SCIPceil(scip, consdata->lhs);
            ++(*nchgsides);
         }
      }
      /* case 3 */
      else if( (SCIPisInfinity(scip, -consdata->lhs) || SCIPisGE(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfloor(scip, consdata->lhs))) && (SCIPisInfinity(scip, consdata->rhs) || SCIPisGT(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfloor(scip, consdata->rhs))) )
      {
         consdata->vbdcoef = SCIPceil(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisIntegral(scip, consdata->lhs) )
         {
            consdata->lhs = SCIPceil(scip, consdata->lhs);
            ++(*nchgsides);
         }
         if( !SCIPisIntegral(scip, consdata->rhs) )
         {
            consdata->rhs = SCIPfloor(scip, consdata->rhs);
            ++(*nchgsides);
         }
      }
      /* case 4 */
      else if( (SCIPisInfinity(scip, -consdata->lhs) || SCIPisLT(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfloor(scip, consdata->lhs))) && (SCIPisInfinity(scip, consdata->rhs) || SCIPisLE(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfloor(scip, consdata->rhs))) )
      {
         consdata->vbdcoef = SCIPfloor(scip, consdata->vbdcoef);
         ++(*nchgcoefs);

         if( !SCIPisIntegral(scip, consdata->lhs) )
         {
            consdata->lhs = SCIPceil(scip, consdata->lhs);
            ++(*nchgsides);
         }
         if( !SCIPisIntegral(scip, consdata->rhs) )
         {
            consdata->rhs = SCIPfloor(scip, consdata->rhs);
            ++(*nchgsides);
         }
      }
      /* case 5 */
      if( !SCIPisIntegral(scip, consdata->lhs) || !SCIPisIntegral(scip, consdata->rhs) )
      {
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            if( SCIPisLT(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfloor(scip, consdata->lhs)) )
            {
               consdata->lhs = SCIPceil(scip, consdata->lhs);
               ++(*nchgsides);
            }
            else if( SCIPisGT(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->lhs - SCIPfloor(scip, consdata->lhs)) )
            {
               consdata->lhs = SCIPfloor(scip, consdata->lhs) + (consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef));
               ++(*nchgsides);
            }
         }
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            if( SCIPisLT(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfloor(scip, consdata->rhs)) )
            {
               consdata->rhs = SCIPfloor(scip, consdata->rhs) + (consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef));
               ++(*nchgsides);
            }
            else if( SCIPisGT(scip, consdata->vbdcoef - SCIPfloor(scip, consdata->vbdcoef), consdata->rhs - SCIPfloor(scip, consdata->rhs)) )
            {
               consdata->rhs = SCIPfloor(scip, consdata->rhs);
               ++(*nchgsides);
            }
         }
      }
   }

   /* check if due to tightening the constraint got redundant */
   if( SCIPisZero(scip, consdata->vbdcoef) )
   {
      assert(SCIPisGE(scip, SCIPvarGetLbGlobal(consdata->var), consdata->lhs));
      assert(SCIPisLE(scip, SCIPvarGetUbGlobal(consdata->var), consdata->rhs));

      SCIPdebugMessage(" -> variable bound constraint <%s> is redundant\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );
      ++(*ndelconss);

      return SCIP_OKAY;
   }

   /* get bounds of variable x */
   xlb = SCIPvarGetLbGlobal(consdata->var);
   xub = SCIPvarGetUbGlobal(consdata->var);

   /* modification of coefficient can only be applied if only one side is finite */
   if( !SCIPisInfinity(scip, -consdata->lhs) && SCIPisInfinity(scip, consdata->rhs) )
   {
      /* lhs <= x + c*y  =>  x >= lhs - c*y */
      if( consdata->vbdcoef > 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs - consdata->vbdcoef) )
      {
         /* constraint has positive slack for the non-restricting case y = 1
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 1 and equivalent in the restricting case y = 0
          * -> c' = lhs - xlb
          */
         SCIPdebugMessage("tighten binary VLB <%s>[%.15g,%.15g] %+.15g<%s> >= %.15g to <%s> %+.15g<%s> >= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs,
            SCIPvarGetName(consdata->var), consdata->lhs - xlb, SCIPvarGetName(consdata->vbdvar), consdata->lhs);

	 /* we cnnot allow that the coefficient changes the sign because of the rounding locks */
	 assert(consdata->vbdcoef * (consdata->lhs - xlb) > 0);

         consdata->vbdcoef = consdata->lhs - xlb;
         (*nchgcoefs)++;
      }
      else if( consdata->vbdcoef < 0.0 && SCIPisFeasGT(scip, xlb, consdata->lhs) )
      {
         /* constraint has positive slack for the non-restricting case y = 0
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 0 and equivalent in the restricting case y = 1
          * -> c' = c - lhs + xlb, lhs' = xlb
          */
         SCIPdebugMessage("tighten binary VLB <%s>[%.15g,%.15g] %+.15g<%s> >= %.15g to <%s> %+.15g<%s> >= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs,
            SCIPvarGetName(consdata->var), consdata->vbdcoef - consdata->lhs + xlb, SCIPvarGetName(consdata->vbdvar), xlb);

	 /* we cnnot allow that the coefficient changes the sign because of the rounding locks */
	 assert(consdata->vbdcoef * (consdata->vbdcoef - consdata->lhs + xlb) > 0);

         consdata->vbdcoef = consdata->vbdcoef - consdata->lhs + xlb;
         consdata->lhs = xlb;
         (*nchgcoefs)++;
         (*nchgsides)++;
      }
   }
   else if( SCIPisInfinity(scip, -consdata->lhs) && !SCIPisInfinity(scip, consdata->rhs) )
   {
      /* x + c*y <= rhs  =>  x <= rhs - c*y */
      if( consdata->vbdcoef < 0.0 && SCIPisFeasLT(scip, xub, consdata->rhs - consdata->vbdcoef) )
      {
         /* constraint has positive slack for the non-restricting case y = 1
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 1 and equivalent in the restricting case y = 0
          * -> c' = rhs - xub
          */
         SCIPdebugMessage("tighten binary VUB <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to <%s> %+.15g<%s> <= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            SCIPvarGetName(consdata->var), consdata->rhs - xub, SCIPvarGetName(consdata->vbdvar), consdata->rhs);

	 /* we cnnot allow that the coefficient changes the sign because of the rounding locks */
	 assert(consdata->vbdcoef * (consdata->rhs - xub) > 0);

         consdata->vbdcoef = consdata->rhs - xub;
         (*nchgcoefs)++;
      }
      else if( consdata->vbdcoef > 0.0 && SCIPisFeasLT(scip, xub, consdata->rhs) )
      {
         /* constraint has positive slack for the non-restricting case y = 0
          * -> modify coefficients such that constraint is tight in the non-restricting case y = 0 and equivalent in the restricting case y = 1
          * -> c' = c - rhs + xub, rhs' = xub
          */
         SCIPdebugMessage("tighten binary VUB <%s>[%.15g,%.15g] %+.15g<%s> <= %.15g to <%s> %+.15g<%s> <= %.15g\n",
            SCIPvarGetName(consdata->var), xlb, xub, consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs,
            SCIPvarGetName(consdata->var), consdata->vbdcoef - consdata->rhs + xub, SCIPvarGetName(consdata->vbdvar), xub);

	 /* we cnnot allow that the coefficient changes the sign because of the rounding locks */
	 assert(consdata->vbdcoef * (consdata->vbdcoef - consdata->rhs + xub) > 0);

         consdata->vbdcoef = consdata->vbdcoef - consdata->rhs + xub;
         consdata->rhs = xub;
         (*nchgcoefs)++;
         (*nchgsides)++;
      }
   }

   /* if something a coefficient or side of the varbound constraint was changed, ensure that the variable lower or
    * upper bounds of the variables are informed
    */
   if( *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
   {
      consdata->addvarbounds = FALSE;
      consdata->changed = TRUE;
   }

   return SCIP_OKAY;
}

/*
 * Linear constraint upgrading
 */

/** tries to upgrade a linear constraint into a variable bound constraint */
static
SCIP_DECL_LINCONSUPGD(linconsUpgdVarbound)
{  /*lint --e{715}*/
   SCIP_Bool upgrade;

   assert(upgdcons != NULL);
   
   /* check, if linear constraint can be upgraded to a variable bound constraint  lhs <= x + a*y <= rhs
    * - there are exactly two variables
    * - one of the variables is non-binary (called the bounded variable x)
    * - one of the variables is non-continuous (called the bounding variable y)
    */
   upgrade = (nvars == 2) && (nposbin + nnegbin <= 1) && (nposcont + nnegcont <= 1);

   if( upgrade )
   {
      SCIP_VAR* var;
      SCIP_VAR* vbdvar;
      SCIP_Real vbdcoef;
      SCIP_Real vbdlhs;
      SCIP_Real vbdrhs;
      int vbdind;

      SCIPdebugMessage("upgrading constraint <%s> to variable bound constraint\n", SCIPconsGetName(cons));

      /* decide which variable we want to use as bounding variable y */
      if( SCIPvarGetType(vars[0]) < SCIPvarGetType(vars[1]) )
         vbdind = 0;
      else if( SCIPvarGetType(vars[0]) > SCIPvarGetType(vars[1]) )
         vbdind = 1;
      else if( SCIPisIntegral(scip, vals[0]) && !SCIPisIntegral(scip, vals[1]) )
         vbdind = 0;
      else if( !SCIPisIntegral(scip, vals[0]) && SCIPisIntegral(scip, vals[1]) )
         vbdind = 1;
      else if( REALABS(REALABS(vals[0]) - 1.0) < REALABS(REALABS(vals[1]) - 1.0) )
         vbdind = 0;
      else
         vbdind = 1;

      var = vars[1-vbdind];
      vbdvar = vars[vbdind];

      assert(!SCIPisZero(scip, vals[1-vbdind]));
      vbdcoef = vals[vbdind]/vals[1-vbdind];

      if( vals[1-vbdind] > 0.0 )
      {
         vbdlhs = SCIPisInfinity(scip, -lhs) ? -SCIPinfinity(scip) : lhs/vals[1-vbdind];
         vbdrhs = SCIPisInfinity(scip, rhs) ? SCIPinfinity(scip) : rhs/vals[1-vbdind];
      }
      else
      {
         vbdlhs = SCIPisInfinity(scip, rhs) ? -SCIPinfinity(scip) : rhs/vals[1-vbdind];
         vbdrhs = SCIPisInfinity(scip, -lhs) ? SCIPinfinity(scip) : lhs/vals[1-vbdind];
      }

      /* create the bin variable bound constraint (an automatically upgraded constraint is always unmodifiable) */
      assert(!SCIPconsIsModifiable(cons));
      SCIP_CALL( SCIPcreateConsVarbound(scip, upgdcons, SCIPconsGetName(cons), var, vbdvar, vbdcoef, vbdlhs, vbdrhs,
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), 
            SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), 
            SCIPconsIsLocal(cons), SCIPconsIsModifiable(cons), 
            SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), SCIPconsIsStickingAtNode(cons)) );
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyVarbound)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrVarbound(scip) );
 
   *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* free constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdataFree(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
#define consInitVarbound NULL


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#define consExitVarbound NULL


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#define consInitpreVarbound NULL


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
#define consExitpreVarbound NULL


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#define consInitsolVarbound NULL


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int c;

   /* release the rows of all constraints */
   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->row != NULL )
      {
         SCIP_CALL( SCIPreleaseRow(scip, &consdata->row) );
      }
   }

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteVarbound)
{  /*lint --e{715}*/
   SCIP_CALL( consdataFree(scip, consdata) );
   
   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */ 
static
SCIP_DECL_CONSTRANS(consTransVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create target constraint data */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->var, sourcedata->vbdvar, sourcedata->vbdcoef, 
         sourcedata->lhs, sourcedata->rhs) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons), 
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler */
static
SCIP_DECL_CONSINITLP(consInitlpVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      assert(SCIPconsIsInitial(conss[i]));
      SCIP_CALL( addRelaxation(scip, conss[i]) );
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpVarbound)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int i;

   *result = SCIP_DIDNOTFIND;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* separate remaining constraints */
   for( i = nusefulconss; i < nconss && *result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolVarbound)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int i;

   *result = SCIP_DIDNOTFIND;

   /* separate useful constraints */
   for( i = 0; i < nusefulconss; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], sol, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   /* separate remaining constraints */
   for( i = nusefulconss; i < nconss && *result == SCIP_DIDNOTFIND; ++i )
   {
      SCIP_CALL( separateCons(scip, conss[i], sol, &separated) );
      if( separated )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpVarbound)
{  /*lint --e{715}*/
   SCIP_Bool separated;
   int i;

   *result = SCIP_FEASIBLE;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, FALSE) )
      {
         SCIP_CALL( separateCons(scip, conss[i], NULL, &separated) );
         if( separated )
         {
            *result = SCIP_SEPARATED;
            break;
         }
         else
            *result = SCIP_INFEASIBLE;
      }
   } 

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], NULL, TRUE) )
      {
         *result = SCIP_INFEASIBLE;
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;  
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckVarbound)
{  /*lint --e{715}*/
   int i;

   for( i = 0; i < nconss; i++ )
   {
      if( !checkCons(scip, conss[i], sol, checklprows) )
      {
         *result = SCIP_INFEASIBLE;

         if( printreason )
         {
            SCIP_Real sum;
            SCIP_CONSDATA* consdata;

            consdata = SCIPconsGetData(conss[i]);
            assert( consdata != NULL );
            
            sum = SCIPgetSolVal(scip, sol, consdata->var);
            sum += consdata->vbdcoef * SCIPgetSolVal(scip, sol, consdata->vbdvar);   
            
            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
            if( !SCIPisFeasGE(scip, sum, consdata->lhs) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhs - sum);
            }
            if( !SCIPisFeasLE(scip, sum, consdata->rhs) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", sum - consdata->rhs);
            }
         }
         return SCIP_OKAY;
      }
   } 
   *result = SCIP_FEASIBLE;
   
   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropVarbound)
{  /*lint --e{715}*/
   SCIP_Bool cutoff;
   int nchgbds;
   int i;

   cutoff = FALSE;
   nchgbds = 0;

   for( i = 0; i < nusefulconss && !cutoff; i++ )
   {
      SCIP_CALL( propagateCons(scip, conss[i], &cutoff, &nchgbds, NULL) );
   } 

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolVarbound)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   int oldnchgbds;
   int oldndelconss;
   int oldnaddconss;
   int oldnchgcoefs;
   int oldnchgsides;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   /* get constraint handler data */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   cutoff = FALSE;
   oldnchgbds = *nchgbds;
   oldndelconss = *ndelconss;
   oldnaddconss = *naddconss;
   oldnchgcoefs = *nchgcoefs;
   oldnchgsides = *nchgsides;

   for( i = 0; i < nconss && !cutoff && !SCIPisStopped(scip); i++ )
   {
      assert(!SCIPconsIsModifiable(conss[i]));

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* force presolving the constraint in the initial round */
      if( nrounds == 0 )
         consdata->presolved = FALSE;

      if( consdata->presolved )
         continue;
      consdata->presolved = TRUE;

      /* make sure that the constraint is propagated */
      consdata->propagated = FALSE;

      /* incorporate fixings and aggregations in constraint */
      SCIP_CALL( applyFixings(scip, conss[i], &cutoff, nchgbds, ndelconss, naddconss) );
      if( cutoff || !SCIPconsIsActive(conss[i]) )
         continue;

      /* propagate constraint */
      SCIP_CALL( propagateCons(scip, conss[i], &cutoff, nchgbds, ndelconss) );
      if( cutoff || !SCIPconsIsActive(conss[i]) )
         continue;

      /* tighten variable bound coefficient */
      SCIP_CALL( tightenCoefs(scip, conss[i], nchgcoefs, nchgsides, ndelconss) );
      if( !SCIPconsIsActive(conss[i]) )
         continue;

      /** informs once variable x about a globally valid variable lower or upper bound */
      if( !consdata->addvarbounds )
      {
         SCIP_Bool infeasible;
         int nlocalchgbds;
         
         nlocalchgbds = 0;

         /* if lhs is finite, we have a variable lower bound: lhs <= x + c*y  =>  x >= -c*y + lhs */
         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIPdebugMessage("adding variable lower bound <%s> >= %g<%s> + %g\n", 
               SCIPvarGetName(consdata->var), -consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->lhs);

            SCIP_CALL( SCIPaddVarVlb(scip, consdata->var, consdata->vbdvar, -consdata->vbdcoef, consdata->lhs,
                  &infeasible, &nlocalchgbds) );
            assert(!infeasible);
            
            *nchgbds += nlocalchgbds;

            /* if lhs is finite, and x is not continuous we can add more variable bounds */
            if( SCIPvarGetType(consdata->var) != SCIP_VARTYPE_CONTINUOUS )
            {
               if( consdata->vbdcoef >= 0.0 )
               {
                  assert(consdata->vbdcoef != 0.0);

                  SCIPdebugMessage("adding variable lower bound <%s> >= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->lhs/consdata->vbdcoef);

                  /* if c > 0, we have a variable lower bound: lhs <= x + c*y  =>  y >= (lhs-x)/c */
                  SCIP_CALL( SCIPaddVarVlb(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->lhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
               else
               {
                  SCIPdebugMessage("adding variable upper bound <%s> <= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->lhs/consdata->vbdcoef);

                  /* if c < 0, we have a variable upper bound: lhs <= x + c*y  =>  y <= (lhs-x)/c */
                  SCIP_CALL( SCIPaddVarVub(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->lhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
            }
         }
         
         /* if rhs is finite, we have a variable upper bound: x + c*y <= rhs  =>  x <= -c*y + rhs */
         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIPdebugMessage("adding variable upper bound <%s> <= %g<%s> + %g\n", 
               SCIPvarGetName(consdata->var), -consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar), consdata->rhs);

            SCIP_CALL( SCIPaddVarVub(scip, consdata->var, consdata->vbdvar, -consdata->vbdcoef, consdata->rhs,
                  &infeasible, &nlocalchgbds) );
            assert(!infeasible);

            *nchgbds += nlocalchgbds;

            /* if rhs is finite, and x is not continuous we can add more variable bounds */
            if( SCIPvarGetType(consdata->var) != SCIP_VARTYPE_CONTINUOUS )
            {
               if( consdata->vbdcoef > 0.0 )
               {
                  assert(consdata->vbdcoef != 0.0);

                  SCIPdebugMessage("adding variable upper bound <%s> <= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->rhs/consdata->vbdcoef);

                  /* if c > 0 we have a variable upper bound: x + c*y <= rhs  =>  y <= (rhs-x)/c */
                  SCIP_CALL( SCIPaddVarVub(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->rhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
               else
               {
                  SCIPdebugMessage("adding variable lower bound <%s> >= %g<%s> + %g\n", 
                     SCIPvarGetName(consdata->vbdvar), -1.0/consdata->vbdcoef, SCIPvarGetName(consdata->var), 
                     consdata->rhs/consdata->vbdcoef);
                  
                  /* if c < 0 we have a variable lower bound: x + c*y <= rhs  =>  y >= (rhs-x)/c */
                  SCIP_CALL( SCIPaddVarVlb(scip, consdata->vbdvar, consdata->var, 
                        -1.0/consdata->vbdcoef, consdata->rhs/consdata->vbdcoef, &infeasible, &nlocalchgbds) );
                  assert(!infeasible);

                  *nchgbds += nlocalchgbds;
               }
            }
         }
         consdata->addvarbounds = TRUE;
      }
   }

   if( !cutoff )
   {
      /* for varbound constraint with two integer variables make coefficients integral */
      SCIP_CALL( prettifyConss(scip, conss, nconss, nchgcoefs, nchgsides) );

      if( conshdlrdata->presolpairwise )
      {
	 /* preprocess pairs of variable bound constraints */
	 SCIP_CALL( preprocessConstraintPairs(scip, conss, nconss, &cutoff, nchgbds, ndelconss, nchgcoefs, nchgsides) );
      }
   }

   /* return the correct result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( *nchgbds > oldnchgbds || *ndelconss > oldndelconss || *naddconss > oldnaddconss
      || *nchgcoefs > oldnchgcoefs || *nchgsides > oldnchgsides )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropVarbound)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, cons, infervar, (PROPRULE)inferinfo, boundtype, bdchgidx) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( !SCIPisInfinity(scip, -consdata->lhs) )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->var, nlockspos, nlocksneg) );
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
   }

   if( !SCIPisInfinity(scip, consdata->rhs) )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->var, nlocksneg, nlockspos) );
      if( consdata->vbdcoef > 0.0 )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlocksneg, nlockspos) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vbdvar, nlockspos, nlocksneg) );
      }
   }

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#define consActiveVarbound NULL


/** constraint deactivation notification method of constraint handler */
#define consDeactiveVarbound NULL


/** constraint enabling notification method of constraint handler */
#define consEnableVarbound NULL


/** constraint disabling notification method of constraint handler */
#define consDisableVarbound NULL


/** variable deletion method of constraint handler:
 *  varbound constraints are not modifiable and must have exactly two variables,
 *  so is is also not allowed to delete variables from them
 */
#define consDelvarsVarbound NULL


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(cons != NULL);
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged rows */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print coefficients and variables */
   SCIPinfoMessage(scip, file, "<%s>[%c] %+.15g<%s>[%c]", SCIPvarGetName(consdata->var), 
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_BINARY ? SCIP_VARTYPE_BINARY_CHAR :
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_INTEGER ? SCIP_VARTYPE_INTEGER_CHAR :
      SCIPvarGetType(consdata->var) == SCIP_VARTYPE_IMPLINT ? SCIP_VARTYPE_IMPLINT_CHAR : SCIP_VARTYPE_CONTINUOUS_CHAR,
      consdata->vbdcoef, SCIPvarGetName(consdata->vbdvar),
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_BINARY ? SCIP_VARTYPE_BINARY_CHAR :
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_INTEGER ? SCIP_VARTYPE_INTEGER_CHAR :
      SCIPvarGetType(consdata->vbdvar) == SCIP_VARTYPE_IMPLINT ? SCIP_VARTYPE_IMPLINT_CHAR : SCIP_VARTYPE_CONTINUOUS_CHAR);
   
   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}

/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyVarbound)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   const char* consname;
   
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, 2) );

   vars[0] = SCIPgetVarVarbound(sourcescip, sourcecons);
   vars[1] = SCIPgetVbdvarVarbound(sourcescip, sourcecons);

   coefs[0] = 1.0;
   coefs[1] = SCIPgetVbdcoefVarbound(sourcescip, sourcecons);

   if( name != NULL )
      consname = name;
   else
      consname = SCIPconsGetName(sourcecons);

   /* copy the varbound using the linear constraint copy method */
   SCIP_CALL( SCIPcopyConsLinear(scip, cons, sourcescip, consname, 2, vars, coefs,
         SCIPgetLhsVarbound(sourcescip, sourcecons), SCIPgetRhsVarbound(sourcescip, sourcecons), varmap, consmap, 
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, valid) );
   
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseVarbound)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real lhs;
   SCIP_Real rhs;
   char* endstr;
   int requiredsize;
   int nvars;

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   (*success) = FALSE;

   /* return of string empty */
   if( !*str )
      return SCIP_OKAY;

   /* ignore whitespace */
   while( isspace(*str) )
      ++str;

   if( isdigit(str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit(str[1])) )
   {
      if( !SCIPstrToRealValue(str, &lhs, &endstr) )
      {
         SCIPerrorMessage("error parsing left hand side\n");
         return SCIP_OKAY;
      }

      /* ignore whitespace */
      while( isspace(*endstr) )
         ++endstr;

      if( endstr[0] != '<' || endstr[1] != '=' )
      {
         SCIPerrorMessage("missing \"<=\" after left hand side(, found %c%c)\n", endstr[0], endstr[1]);
         return SCIP_OKAY;
      }

      SCIPdebugMessage("found left hand side <%g>\n", lhs);

      /* it was indeed a left-hand-side, so continue parsing after it */
      str = endstr + 2;
   }

   /* pares x + c*y as linear sum */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, 2) );

   /* parse linear sum to get variables and coefficients */
   SCIP_CALL( SCIPparseVarsLinearsum(scip, str, vars, coefs, &nvars, 2, &requiredsize, &endstr, success) );

   if( requiredsize == 2 && *success )
   {
      SCIP_Bool foundvalue;
      SCIP_Real value;

      assert(nvars == 2);
      assert(SCIPisEQ(scip, coefs[0], 1.0));

      SCIPdebugMessage("found linear sum <%s> + %g <%s>\n", SCIPvarGetName(vars[0]), coefs[1], SCIPvarGetName(vars[1]));

      /* ignore whitespace */
      while( isspace(*endstr) )
         ++endstr;

      str = endstr;

      foundvalue = SCIPstrToRealValue(str+2, &value, &endstr);

      if( foundvalue )
      {
         /* search for end of linear sum: either '<=', '>=', '==', or '[free]' */
         switch( *str )
         {
         case '<':
            assert(str[1] == '=');
            rhs = value;
            break;
         case '=':
            assert(str[1] == '=');
            assert(SCIPisInfinity(scip, -lhs));
            lhs = value;
            rhs = value;
            break;
         case '>':
            assert(str[1] == '=');
            assert(SCIPisInfinity(scip, -lhs));
            lhs = value;
            break;
         default:
            SCIPerrorMessage("missing relation symbol after linear sum\n");
            *success = FALSE;
         }
      }
      else if( !strncmp(str, "[free]", 6) == 0 )
         *success = FALSE;
   }

   if( *success )
   {
      SCIP_CALL( SCIPcreateConsVarbound(scip, cons, name, vars[0], vars[1], coefs[1], lhs, rhs,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}





/*
 * Event Handler
 */

/** execution method of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecVarbound)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   consdata = (SCIP_CONSDATA*)eventdata;
   assert(consdata != NULL);

   consdata->propagated = FALSE;
   consdata->presolved = FALSE;

   return SCIP_OKAY;
}




/*
 * constraint specific interface methods
 */

/** creates the handler for variable bound constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrVarbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create variable bound constraint handler data */
   SCIP_CALL( conshdlrdataCreate(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS, 
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_DELAYPRESOL, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING,
         conshdlrCopyVarbound,
         consFreeVarbound, consInitVarbound, consExitVarbound, 
         consInitpreVarbound, consExitpreVarbound, consInitsolVarbound, consExitsolVarbound,
         consDeleteVarbound, consTransVarbound, consInitlpVarbound,
         consSepalpVarbound, consSepasolVarbound, consEnfolpVarbound, consEnfopsVarbound, consCheckVarbound, 
         consPropVarbound, consPresolVarbound, consRespropVarbound, consLockVarbound,
         consActiveVarbound, consDeactiveVarbound, 
         consEnableVarbound, consDisableVarbound,
         consDelvarsVarbound, consPrintVarbound, consCopyVarbound, consParseVarbound,
         conshdlrdata) );

   if( SCIPfindConshdlr(scip,"linear") != NULL )
   {
      /* include the linear constraint to varbound constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeLinconsUpgrade(scip, linconsUpgdVarbound, LINCONSUPGD_PRIORITY, CONSHDLR_NAME) );
   }

   /* include event handler for bound change events */
   eventhdlrdata = NULL;
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecVarbound,
         eventhdlrdata) );

   /* add varbound constraint handler parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "constraints/"CONSHDLR_NAME"/presolpairwise",
         "should pairwise constraint comparison be performed in presolving?",
         &conshdlrdata->presolpairwise, TRUE, DEFAULT_PRESOLPAIRWISE, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a variable bound constraint: lhs <= x + c*y <= rhs
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             var,                /**< variable x that has variable bound */
   SCIP_VAR*             vbdvar,             /**< binary, integer or implicit integer bounding variable y */
   SCIP_Real             vbdcoef,            /**< coefficient c of bounding variable y */
   SCIP_Real             lhs,                /**< left hand side of variable bound inequality */
   SCIP_Real             rhs,                /**< right hand side of variable bound inequality */
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
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the variable bound constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("variable bound constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( consdataCreate(scip, &consdata, var, vbdvar, vbdcoef, lhs, rhs) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** gets left hand side of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetLhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->lhs;
}

/** gets right hand side of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetRhsVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->rhs;
}

/** gets bounded variable x of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_VAR* SCIPgetVarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->var;
}

/** gets bounding variable y of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_VAR* SCIPgetVbdvarVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vbdvar;
}

/** gets bound coefficient c of variable bound constraint lhs <= x + c*y <= rhs */
SCIP_Real SCIPgetVbdcoefVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->vbdcoef;
}

/** gets the dual solution of the variable bound constraint in the current LP */
SCIP_Real SCIPgetDualsolVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualsol(consdata->row);
   else
      return 0.0;
}

/** gets the dual Farkas value of the variable bound constraint in the current infeasible LP */
SCIP_Real SCIPgetDualfarkasVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }
   
   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->row != NULL )
      return SCIProwGetDualfarkas(consdata->row);
   else
      return 0.0;
}

/** returns the linear relaxation of the given variable bound constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
SCIP_ROW* SCIPgetRowVarbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not a variable bound constraint\n");
      SCIPABORT();
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->row;
}

