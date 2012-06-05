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

/**@file   heur_trysol.c
 * @brief  primal heuristic that tries a given solution
 * @author Marc Pfetsch
 *
 * This heuristic takes a solution from somewhere else via the function SCIPheurPassSolTrySol(). It
 * then tries to commit this solution. It is mainly used by cons_indicator, which tries to correct a
 * given solution, but cannot directly submit this solution, because it is a constraint handler and
 * not a heuristic.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_trysol.h"


#define HEUR_NAME             "trysol"
#define HEUR_DESC             "try solution heuristic"
#define HEUR_DISPCHAR         'y'
#define HEUR_PRIORITY         -3000000     /* should process after all other heuristics */
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_DURINGLPLOOP | SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */




/*
 * Data structures
 */


/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*      sol;                /**< storing solution passed to heuristic (NULL if none) */
   SCIP_Bool      rec;                /**< whether we are within our own call */
};




/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTrySol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTrySol(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeTrySol)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   SCIPdebugMessage("free method of trysol primal heuristic.\n");

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeMemory(scip, &heurdata);

   return SCIP_OKAY;
}


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitTrySol)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   SCIPdebugMessage("exit method of trysol primal heuristic.\n");

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* free solution if one is still present */
   if( heurdata->sol != NULL )
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );
   assert( heurdata->sol == NULL );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTrySol)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Bool stored;
#ifdef SCIP_DEBUG
   SCIP_Real obj;
#endif

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only run if solution present */
   if( heurdata->sol == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage("exec method of trysol primal heuristic.\n");
   *result = SCIP_DIDNOTFIND;
   heurdata->rec = TRUE;

   /* try solution and free it - check everything, because we are not sure */
#ifdef SCIP_DEBUG
   obj = SCIPgetSolOrigObj(scip, heurdata->sol);
#endif
   SCIP_CALL( SCIPtrySolFree(scip, &heurdata->sol, FALSE, TRUE, TRUE, TRUE, &stored) );
   assert( heurdata->sol == NULL );

   if( stored )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMessage("Found feasible solution of value %g.\n", obj);
#endif
      *result = SCIP_FOUNDSOL;
   }
   heurdata->rec = FALSE;

   return SCIP_OKAY;
}


#define heurInitTrySol NULL
#define heurExitsolTrySol NULL
#define heurInitsolTrySol NULL


/*
 * primal heuristic specific interface methods
 */

/** creates the trysol primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTrySol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   heurdata->sol = NULL;
   heurdata->rec = FALSE;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyTrySol,heurFreeTrySol, heurInitTrySol, heurExitTrySol,
         heurInitsolTrySol, heurExitsolTrySol, heurExecTrySol,
         heurdata) );

   return SCIP_OKAY;
}


/** pass solution to trysol heuristic */
SCIP_RETCODE SCIPheurPassSolTrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< trysol heuristic */
   SCIP_SOL*             sol                 /**< solution to be passed */
   )
{
   SCIP_HEURDATA* heurdata;

   assert( scip != NULL );
   assert( heur != NULL );
   assert( sol != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* only store solution if we are not within our own SCIPtrySol() call */
   if( ! heurdata->rec )
   {
      if( heurdata->sol == NULL || SCIPisLT(scip, SCIPgetSolOrigObj(scip, sol), SCIPgetSolOrigObj(scip, heurdata->sol)) )
      {
         if( heurdata->sol != NULL )
            SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

         SCIPdebugMessage("Received solution of value %g.\n", SCIPgetSolOrigObj(scip, sol));
         SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->sol, sol) );
         SCIP_CALL( SCIPunlinkSol(scip, heurdata->sol) );
         SCIPsolSetHeur(heurdata->sol, heur);
      }
   }

   return SCIP_OKAY;
}
