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

/**@file   sepa_intobj.c
 * @brief  integer objective value separator
 * @author Tobias Achterberg
 * @author Timo Berthold
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_intobj.h"


#define SEPA_NAME              "intobj"
#define SEPA_DESC              "integer objective value separator"
#define SEPA_PRIORITY              -100
#define SEPA_FREQ                    -1
#define SEPA_MAXBOUNDDIST           0.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define EVENTHDLR_NAME         "intobj"
#define EVENTHDLR_DESC         "objective change event handler for integer objective value separator"




/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_ROW*             objrow;             /**< objective value inequality */
   SCIP_VAR*             objvar;             /**< objective value variable */
   SCIP_Real             setoff;             /**< setoff of the inequality */
};




/*
 * Local methods
 */

/** creates separator data */
static
SCIP_RETCODE sepadataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< pointer to store separator data */
   )
{
   assert(sepadata != NULL);

   SCIP_CALL( SCIPallocMemory(scip, sepadata) );
   (*sepadata)->objrow = NULL;
   (*sepadata)->objvar = NULL;
   (*sepadata)->setoff = 0.0;

   return SCIP_OKAY;
}

/** frees separator data */
static
SCIP_RETCODE sepadataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< pointer to separator data */
   )
{
   assert(sepadata != NULL);
   assert(*sepadata != NULL);
   assert((*sepadata)->objrow == NULL);
   assert((*sepadata)->objvar == NULL);

   SCIPfreeMemory(scip, sepadata);

   return SCIP_OKAY;
}

/** creates the objective value inequality and the objective value variable, if not yet existing */
static
SCIP_RETCODE createObjRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data */
    )
{
   assert(sepadata != NULL);

   if( sepadata->objrow == NULL )
   {
      SCIP_VAR** vars;
      SCIP_Real obj;
      SCIP_Real intobjval;
      int nvars;
      int v;
      SCIP_Bool attendobjvarbound;
      
      attendobjvarbound = FALSE;
      /* create and add objective value variable */
      if( sepadata->objvar == NULL )
      {
         SCIP_CALL( SCIPcreateVar(scip, &sepadata->objvar, "objvar", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_IMPLINT, FALSE, TRUE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, sepadata->objvar) );
         SCIP_CALL( SCIPaddVarLocks(scip, sepadata->objvar, +1, +1) );
      }
      else
         attendobjvarbound = TRUE;
      
      /* get problem variables */
      vars = SCIPgetOrigVars(scip);
      nvars = SCIPgetNOrigVars(scip);

      /* create objective value inequality */
      if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
      {
         if( attendobjvarbound )
            intobjval = SCIPceil(scip, SCIPgetDualbound(scip)) - SCIPvarGetLbGlobal(sepadata->objvar);
         else
            intobjval = SCIPceil(scip, SCIPgetDualbound(scip));
         SCIP_CALL( SCIPcreateEmptyRow(scip, &sepadata->objrow, "objrow", intobjval, SCIPinfinity(scip),
               FALSE, !SCIPallVarsInProb(scip), TRUE) );
         sepadata->setoff = intobjval;
      }
      else
      {
         if( attendobjvarbound )
            intobjval = SCIPceil(scip, SCIPgetDualbound(scip)) - SCIPvarGetUbGlobal(sepadata->objvar);
         else
            intobjval = SCIPfloor(scip, SCIPgetDualbound(scip));
         SCIP_CALL( SCIPcreateEmptyRow(scip, &sepadata->objrow, "objrow", -SCIPinfinity(scip), intobjval,
               FALSE, !SCIPallVarsInProb(scip), TRUE) );
         sepadata->setoff = intobjval;
      }         

      SCIP_CALL( SCIPcacheRowExtensions(scip, sepadata->objrow) );
      for( v = 0; v < nvars; ++v )
      {
         obj = SCIPvarGetObj(vars[v]);
         if( !SCIPisZero(scip, obj) )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, sepadata->objrow, vars[v], obj) );
         }
      }
      SCIP_CALL( SCIPaddVarToRow(scip, sepadata->objrow, sepadata->objvar, -1.0) );
      SCIP_CALL( SCIPflushRowExtensions(scip, sepadata->objrow) );

      SCIPdebugMessage("created objective value row: ");
      SCIPdebug(SCIPprintRow(scip, sepadata->objrow, NULL));
   }

   return SCIP_OKAY;
}

/** searches and adds integral objective cuts that separate the given primal solution */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_SEPA*            sepa,               /**< the intobj separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_Real objval;
   SCIP_Real intbound;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert(result != NULL);
   assert(*result == SCIP_DIDNOTRUN);

   /* if the objective value may be fractional, we cannot do anything */
   if( !SCIPisObjIntegral(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* if the current objective value is integral, there is no integral objective value cut */
   if( sol == NULL )
      objval = SCIPretransformObj(scip, SCIPgetLPObjval(scip));
   else
      objval = SCIPgetSolOrigObj(scip, sol);
   if( SCIPisFeasIntegral(scip, objval) )
      return SCIP_OKAY;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* the objective value is fractional: create the objective value inequality, if not yet existing */
   SCIP_CALL( createObjRow(scip, sepadata) );

   /* adjust the bounds of the objective value variable */
   if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
   {
      intbound = SCIPceil(scip, objval) - sepadata->setoff;
      SCIP_CALL( SCIPtightenVarLb(scip, sepadata->objvar, intbound, FALSE, &infeasible, &tightened) );
      SCIPdebugMessage("new objective variable lower bound: <%s>[%g,%g]\n", 
         SCIPvarGetName(sepadata->objvar), SCIPvarGetLbLocal(sepadata->objvar), SCIPvarGetUbLocal(sepadata->objvar));
   }
   else
   {
      intbound = SCIPfloor(scip, objval) - sepadata->setoff;
      SCIP_CALL( SCIPtightenVarUb(scip, sepadata->objvar, intbound, FALSE, &infeasible, &tightened) );
      SCIPdebugMessage("new objective variable upper bound: <%s>[%g,%g]\n", 
         SCIPvarGetName(sepadata->objvar), SCIPvarGetLbLocal(sepadata->objvar), SCIPvarGetUbLocal(sepadata->objvar));
   }

   /* add the objective value inequality as a cut to the LP */
   if( infeasible )
      *result = SCIP_CUTOFF;
   else
   {
      if( !SCIProwIsInLP(sepadata->objrow) )
      {
         SCIP_CALL( SCIPaddCut(scip, sol, sepadata->objrow, FALSE) );
      }
      if( tightened )
         *result = SCIP_REDUCEDDOM;
      else
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}




/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyIntobj)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaIntobj(scip) );
 
   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeIntobj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( sepadataFree(scip, &sepadata) );

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#define sepaInitIntobj NULL


/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitIntobj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* release objective variable */
   if( sepadata->objvar != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &sepadata->objvar) );
   }

   return SCIP_OKAY;
}


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolIntobj NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolIntobj)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* release objective row */
   if( sepadata->objrow != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &sepadata->objrow) );
   }

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpIntobj)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   /* only call separator, if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only call separator, if there are fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( separateCuts(scip, sepa, NULL, result) );

   return SCIP_OKAY;
}


/** arbitrary primal solution separation method of separator */
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolIntobj)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( separateCuts(scip, sepa, sol, result) );

   return SCIP_OKAY;
}




/*
 * event handler for objective changes
 */

/** destructor of event handler to free user data (called when SCIP is exiting) */
#define eventFreeIntobj NULL

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitIntobj)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_OBJCHANGED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitIntobj)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED | SCIP_EVENTTYPE_OBJCHANGED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** solving process initialization method of event handler (called when branch and bound process is about to begin) */
#define eventInitsolIntobj NULL

/** solving process deinitialization method of event handler (called before branch and bound process data is freed) */
#define eventExitsolIntobj NULL

/** frees specific event data */
#define eventDeleteIntobj NULL

/** execution method of objective change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecIntobj)
{  /*lint --e{715}*/
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_SEPADATA* sepadata;
   SCIP_VAR* var;
   SCIP_Real objdelta;

   eventhdlrdata = SCIPeventhdlrGetData(eventhdlr);
   sepadata = (SCIP_SEPADATA*)eventhdlrdata;
   assert(sepadata != NULL);

   /* we don't have anything to do, if the objective value inequality doesn't yet exist */
   if( sepadata->objrow == NULL )
      return SCIP_OKAY;

   var = SCIPeventGetVar(event);

   switch( SCIPeventGetType(event) )
   {
   case SCIP_EVENTTYPE_VARADDED:
      SCIPdebugMessage("variable <%s> with obj=%g was added to the problem\n", SCIPvarGetName(var), SCIPvarGetObj(var));
      objdelta = SCIPvarGetObj(var);
      if( !SCIPisZero(scip, objdelta) )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, sepadata->objrow, var, SCIPvarGetObj(var)) );
      }
      break;

   case SCIP_EVENTTYPE_OBJCHANGED:
      SCIPdebugMessage("variable <%s> changed objective value from %g to %g\n", 
         SCIPvarGetName(var), SCIPeventGetOldobj(event), SCIPeventGetNewobj(event));
      objdelta = SCIPeventGetNewobj(event) - SCIPeventGetOldobj(event);
      SCIP_CALL( SCIPaddVarToRow(scip, sepadata->objrow, var, objdelta) );
      break;

   default:
      SCIPerrorMessage("invalid event type %x\n", SCIPeventGetType(event));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}




/*
 * separator specific interface methods
 */

/** creates the integer objective value separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaIntobj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_EVENTHDLRDATA* eventhdlrdata;

   /* create intobj separator data */
   SCIP_CALL( sepadataCreate(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyIntobj, sepaFreeIntobj, sepaInitIntobj, sepaExitIntobj, 
         sepaInitsolIntobj, sepaExitsolIntobj, 
         sepaExeclpIntobj, sepaExecsolIntobj,
         sepadata) );

   /* include event handler for objective change events */
   eventhdlrdata = (SCIP_EVENTHDLRDATA*)sepadata;
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC, 
         NULL,
         eventFreeIntobj, eventInitIntobj, eventExitIntobj, 
         eventInitsolIntobj, eventExitsolIntobj, eventDeleteIntobj, eventExecIntobj,
         eventhdlrdata) );

   return SCIP_OKAY;
}
