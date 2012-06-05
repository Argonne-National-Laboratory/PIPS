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

/**@file   heur_rens.c
 * @brief  LNS heuristic that finds the optimal rounding to a given point
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "scip/scip.h"
#include "scip/heur_rens.h"
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "scip/cons_linear.h"          /* needed if the LP relaxation gets copied into linear constraints */
#include "scip/pub_misc.h"

/* default values for standard parameters that every primal heuristic has in SCIP */
#define HEUR_NAME             "rens"
#define HEUR_DESC             "LNS exploring fractional neighborhood of relaxation's optimum"
#define HEUR_DISPCHAR         'E'
#define HEUR_PRIORITY         -1100000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

/* default values for RENS-specific plugins */
#define DEFAULT_BINARYBOUNDS  TRUE      /* should general integers get binary bounds [floor(.),ceil(.)] ?      */
#define DEFAULT_MAXNODES      5000LL    /* maximum number of nodes to regard in the subproblem                 */
#define DEFAULT_MINFIXINGRATE 0.5       /* minimum percentage of integer variables that have to be fixed       */
#define DEFAULT_MINIMPROVE    0.01      /* factor by which RENS should at least improve the incumbent          */
#define DEFAULT_MINNODES      500LL     /* minimum number of nodes to regard in the subproblem                 */
#define DEFAULT_NODESOFS      500LL     /* number of nodes added to the contingent of the total nodes          */
#define DEFAULT_NODESQUOT     0.1       /* subproblem nodes in relation to nodes of the original problem       */
#define DEFAULT_STARTSOL      'l'       /* solution that is used for fixing values                             */
#define STARTSOL_CHOICES      "nl"      /* possible values for startsol ('l'p relaxation, 'n'lp relaxation)    */
#define DEFAULT_USELPROWS    FALSE      /* should subproblem be created out of the rows in the LP rows,
                                         * otherwise, the copy constructors of the constraints handlers are used
                                         */
#define DEFAULT_COPYCUTS      TRUE      /* if DEFAULT_USELPROWS is FALSE, then should all active cuts from the cutpool
                                         * of the original scip be copied to constraints of the subscip
                                         */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          maxnodes;           /**< maximum number of nodes to regard in the subproblem                 */
   SCIP_Longint          minnodes;           /**< minimum number of nodes to regard in the subproblem                 */
   SCIP_Longint          nodesofs;           /**< number of nodes added to the contingent of the total nodes          */
   SCIP_Longint          usednodes;          /**< nodes already used by RENS in earlier calls                         */
   SCIP_Real             minfixingrate;      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove;         /**< factor by which RENS should at least improve the incumbent          */
   SCIP_Real             nodesquot;          /**< subproblem nodes in relation to nodes of the original problem       */
   char                  startsol;           /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Bool             binarybounds;       /**< should general integers get binary bounds [floor(.),ceil(.)] ?      */
   SCIP_Bool             uselprows;          /**< should subproblem be created out of the rows in the LP rows?        */
   SCIP_Bool             copycuts;           /**< if uselprows == FALSE, should all active cuts from cutpool be copied
                                              *   to constraints in subproblem?
                                              */
};


/*
 * Local methods
 */

/** creates a subproblem for subscip by fixing a number of variables */
static
SCIP_RETCODE createSubproblem(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP*                 subscip,            /**< SCIP data structure for the subproblem                              */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                                     */
   SCIP_Real             minfixingrate,      /**< percentage of integer variables that have to be fixed               */
   char                  startsol,           /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)] ?      */
   SCIP_Bool             uselprows,          /**< should subproblem be created out of the rows in the LP rows?        */
   SCIP_Bool*            success             /**< pointer to store whether the problem was created successfully       */
   )
{
   SCIP_VAR** vars;                          /* original SCIP variables */

   SCIP_Real fixingrate;

   int nvars;
   int nbinvars;
   int nintvars;
   int i;
   int fixingcounter;

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);

   assert(0.0 <= minfixingrate && minfixingrate <= 1.0);
   assert(startsol == 'l' || startsol == 'n');

   /* if there is no NLP relaxation available (e.g., because the presolved problem is linear), use LP relaxation */
   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPdebugMessage("no NLP present, use LP relaxation instead\n");
      startsol = 'l';
   }

   /* get required variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, NULL, NULL) );

   fixingcounter = 0;

   /* try to solve nlp relaxation */
   if( startsol == 'n')
   {
      SCIP_NLPSOLSTAT stat;
      SCIPdebug( int nlpverblevel; )

      /* only call this function if nlp relaxation is available */
      assert(SCIPisNLPConstructed(scip));

      /* activate nlp solver output if we are in SCIP's debug mode */
      SCIPdebug( SCIP_CALL( SCIPgetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, &nlpverblevel) ) );
      SCIPdebug( SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, MAX(1,nlpverblevel)) ) );

      SCIPdebugMessage("try to solve nlp relaxation to obtain fixing values\n");

      /* set starting point to lp solution */
      SCIP_CALL( SCIPsetNLPInitialGuessSol(scip, NULL) );

      /* solve nlp relaxation */
      SCIP_CALL( SCIPsolveNLP(scip) );

      /* get solution status of NLP solver */
      stat = SCIPgetNLPSolstat(scip);
      *success = (stat == SCIP_NLPSOLSTAT_GLOBOPT) || (stat == SCIP_NLPSOLSTAT_LOCOPT) || stat == (SCIP_NLPSOLSTAT_FEASIBLE);
      SCIPdebugMessage("solving nlp relaxation was %s successful (stat=%d)\n", *success ? "" : "not", stat);

      /* reset nlp verblevel to the value it had before */
      SCIPdebug( SCIP_CALL( SCIPsetNLPIntPar(scip, SCIP_NLPPAR_VERBLEVEL, nlpverblevel) ) );
   }

   /* change bounds of variables of the subproblem */
   for( i = 0; i < nbinvars + nintvars; i++ )
   {
      SCIP_Real solval;
      SCIP_Real lb;
      SCIP_Real ub;

      /* get the current LP solution for each variable */
      if( startsol == 'l')
         solval = SCIPvarGetLPSol(vars[i]);
      else
         solval = SCIPvarGetNLPSol(vars[i]);

      if( SCIPisFeasIntegral(scip, solval) )
      {
         /* fix variables to current LP solution if it is integral,
          * use exact integral value, if the variable is only integral within numerical tolerances
          */
         lb = SCIPfloor(scip, solval+0.5);
         ub = lb;
         fixingcounter++;
      }
      else if( binarybounds )
      {
         /* if the sub problem should be a binary problem, change the bounds to nearest integers */
         lb = SCIPfeasFloor(scip,solval);
         ub = SCIPfeasCeil(scip,solval);
      }
      else
      {
         /* otherwise just copy bounds */
         lb =  SCIPvarGetLbGlobal(vars[i]);
         ub =  SCIPvarGetUbGlobal(vars[i]);
      }

      /* perform the bound change */
      SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], lb) );
      SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], ub) );
   }

   fixingrate = 0.0;

   /* abort, if all integer variables were fixed (which should not happen for MIP) */
   if( fixingcounter == nbinvars + nintvars )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   else
      fixingrate = fixingcounter / (SCIP_Real)(MAX(nbinvars + nintvars, 1));
   SCIPdebugMessage("fixing rate: %g = %d of %d\n", fixingrate, fixingcounter, nbinvars + nintvars);

   /* abort, if the amount of fixed variables is insufficient */
   if( fixingrate < minfixingrate )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   if( uselprows )
   {
      SCIP_ROW** rows;                          /* original scip rows                         */
      int nrows;

      /* get the rows and their number */
      SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

      /* copy all rows to linear constraints */
      for( i = 0; i < nrows; i++ )
      {
         SCIP_CONS* cons;
         SCIP_VAR** consvars;
         SCIP_COL** cols;
         SCIP_Real constant;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real* vals;
         int nnonz;
         int j;

         /* ignore rows that are only locally valid */
         if( SCIProwIsLocal(rows[i]) )
            continue;

         /* get the row's data */
         constant = SCIProwGetConstant(rows[i]);
         lhs = SCIProwGetLhs(rows[i]) - constant;
         rhs = SCIProwGetRhs(rows[i]) - constant;
         vals = SCIProwGetVals(rows[i]);
         nnonz = SCIProwGetNNonz(rows[i]);
         cols = SCIProwGetCols(rows[i]);

         assert( lhs <= rhs );

         /* allocate memory array to be filled with the corresponding subproblem variables */
         SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, nnonz) );
         for( j = 0; j < nnonz; j++ )
            consvars[j] = subvars[SCIPvarGetProbindex(SCIPcolGetVar(cols[j]))];

         /* create a new linear constraint and add it to the subproblem */
         SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, SCIProwGetName(rows[i]), nnonz, consvars, vals, lhs, rhs,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, cons) );
         SCIP_CALL( SCIPreleaseCons(subscip, &cons) );

         /* free temporary memory */
         SCIPfreeBufferArray(subscip, &consvars);
      }
   }

   *success = TRUE;
   return SCIP_OKAY;
}


/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< RENS heuristic structure                            */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;                         /* the original problem's number of variables      */
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */

   assert(scip != NULL);
   assert(subscip != NULL);
   assert(subvars != NULL);
   assert(subsol != NULL);

   /* get variables' data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */
   assert(nvars <= SCIPgetNOrigVars(subscip));

   SCIP_CALL( SCIPallocBufferArray(scip, &subsolvals, nvars) );

   /* copy the solution */
   SCIP_CALL( SCIPgetSolVals(subscip, subsol, nvars, subvars, subsolvals) );

   /* create new solution for the original problem */
   SCIP_CALL( SCIPcreateSol(scip, &newsol, heur) );
   SCIP_CALL( SCIPsetSolVals(scip, newsol, nvars, vars, subsolvals) );

   /* try to add new solution to scip and free it immediately */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );

   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/** main procedure of the RENS heuristic, creates and solves a sub-SCIP */
SCIP_RETCODE SCIPapplyRens(
   SCIP*                 scip,               /**< original SCIP data structure                                        */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                            */
   SCIP_RESULT*          result,             /**< result data structure                                               */
   SCIP_Real             minfixingrate,      /**< minimum percentage of integer variables that have to be fixed       */
   SCIP_Real             minimprove,         /**< factor by which RENS should at least improve the incumbent          */
   SCIP_Longint          maxnodes,           /**< maximum number of  nodes for the subproblem                         */
   SCIP_Longint          nstallnodes,        /**< number of stalling nodes for the subproblem                         */
   char                  startsol,           /**< solution used for fixing values ('l'p relaxation, 'n'lp relaxation) */
   SCIP_Bool             binarybounds,       /**< should general integers get binary bounds [floor(.),ceil(.)]?       */
   SCIP_Bool             uselprows           /**< should subproblem be created out of the rows in the LP rows?        */
   )
{
   SCIP* subscip;                            /* the subproblem created by RENS                  */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */
   SCIP_VAR** vars;                          /* original problem's variables                    */
   SCIP_VAR** subvars;                       /* subproblem's variables                          */

   SCIP_Real cutoff;                         /* objective cutoff for the subproblem             */
   SCIP_Real timelimit;
   SCIP_Real memorylimit;

   int nvars;
   int i;

   SCIP_Bool success;
   SCIP_RETCODE retcode;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   assert(maxnodes >= 0);
   assert(nstallnodes >= 0);

   assert(0.0 <= minfixingrate && minfixingrate <= 1.0);
   assert(0.0 <= minimprove && minimprove <= 1.0);
   assert(startsol == 'l' || startsol == 'n');

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );

   if( uselprows )
   {
      char probname[SCIP_MAXSTRLEN];

      /* copy all plugins */
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

      /* get name of the original problem and add the string "_renssub" */
      (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_renssub", SCIPgetProbName(scip));

      /* create the subproblem */
      SCIP_CALL( SCIPcreateProb(subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* copy all variables */
      SCIP_CALL( SCIPcopyVars(scip, subscip, varmapfw, NULL, TRUE) );
   }
   else
   {
      SCIP_Bool valid;
      SCIP_HEURDATA* heurdata;

      valid = FALSE;

      SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "rens", TRUE, FALSE, &valid) );

      /* get heuristic's data */
      heurdata = SCIPheurGetData(heur);
      assert( heurdata != NULL );

      if( heurdata->copycuts )
      {
         /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
         SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, TRUE) );
      }

      SCIPdebugMessage("Copying the SCIP instance was %s complete.\n", valid ? "" : "not ");
   }

   for( i = 0; i < nvars; i++ )
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* create a new problem, which fixes variables with same value in bestsol and LP relaxation */
   SCIP_CALL( createSubproblem(scip, subscip, subvars, minfixingrate, startsol, binarybounds, uselprows, &success) );
   SCIPdebugMessage("RENS subproblem: %d vars, %d cons, success=%u\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip), success);

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* check whether there is enough time and memory left */
   timelimit = 0.0;
   memorylimit = 0.0;
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit <= 0.0 || memorylimit <= 0.0 )
      goto TERMINATE;

   /* set limits for the subproblem */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", nstallnodes) );
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", maxnodes) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );

   /* forbid recursive call of heuristics and separators solving sub-SCIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* use best estimate node selection */
   if( SCIPfindNodesel(scip, "estimate") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/estimate/stdpriority", INT_MAX/4) );
   }

   /* use inference branching */
   if( SCIPfindBranchrule(scip, "inference") != NULL )
   {
      SCIP_CALL( SCIPsetIntParam(subscip, "branching/inference/priority", INT_MAX/4) );
   }

   /* disable conflict analysis */
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useprop", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useinflp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/useboundlp", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usesb", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "conflict/usepseudo", FALSE) );


#ifdef SCIP_DEBUG
   /* for debugging RENS, enable MIP output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#endif


   /* if the subproblem could not be created, free memory and return */
   if( !success )
   {
      *result = SCIP_DIDNOTRUN;
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY;
   }

   /* if there is already a solution, add an objective cutoff */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;
      cutoff = SCIPinfinity(scip);
      assert( !SCIPisInfinity(scip,SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip,-1.0*SCIPgetLowerbound(scip)) )
      {
         cutoff = (1-minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound ( scip ) >= 0 )
            cutoff = ( 1 - minimprove ) * SCIPgetUpperbound ( scip );
         else
            cutoff = ( 1 + minimprove ) * SCIPgetUpperbound ( scip );
      }
      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL( SCIPsetObjlimit(subscip, cutoff) );
   }

   /* presolve the subproblem */
   retcode = SCIPpresolve(subscip);

   /* Errors in solving the subproblem should not kill the overall solving process
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   {
#ifndef NDEBUG
      SCIP_CALL( retcode );
#endif
      SCIPwarningMessage("Error while presolving subproblem in RENS heuristic; sub-SCIP terminated with code <%d>\n",retcode);
   }

   SCIPdebugMessage("RENS presolved subproblem: %d vars, %d cons, success=%u\n", SCIPgetNVars(subscip), SCIPgetNConss(subscip), success);

   /* after presolving, we should have at least reached a certain fixing rate over ALL variables (including continuous)
    * to ensure that not only the MIP but also the LP relaxation is easy enough
    */
   if( ( nvars - SCIPgetNVars(subscip) ) / (SCIP_Real)nvars >= minfixingrate / 2.0 )
   {
      SCIP_SOL** subsols;
      int nsubsols;

      /* solve the subproblem */
      SCIPdebugMessage("solving subproblem: nstallnodes=%"SCIP_LONGINT_FORMAT", maxnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, maxnodes);
      retcode = SCIPsolve(subscip);

      /* Errors in solving the subproblem should not kill the overall solving process
       * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      if( retcode != SCIP_OKAY )
      {
#ifndef NDEBUG
         SCIP_CALL( retcode );
#endif
         SCIPwarningMessage("Error while solving subproblem in RENS heuristic; sub-SCIP terminated with code <%d>\n",retcode);
      }

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until one was accepted
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      success = FALSE;
      for( i = 0; i < nsubsols && !success; ++i )
      {
         SCIP_CALL( createNewSol(scip, subscip, subvars, heur, subsols[i], &success) );
      }
      if( success )
         *result = SCIP_FOUNDSOL;
   }

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyRens)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurRens(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeRens)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitRens)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* initialize data */
   heurdata->usednodes = 0;

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitRens NULL

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#define heurInitsolRens NULL

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#define heurExitsolRens NULL

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecRens)
{  /*lint --e{715}*/

   SCIP_HEURDATA* heurdata;                  /* heuristic's data                    */
   SCIP_Longint nstallnodes;                 /* number of stalling nodes for the subproblem */

   assert( heur != NULL );
   assert( scip != NULL );
   assert( result != NULL );
   assert( SCIPhasCurrentNodeLP(scip) );

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   *result = SCIP_DELAYED;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* only continue with some fractional variables */
   if( SCIPgetNLPBranchCands(scip) == 0 )
      return SCIP_OKAY;

   /* do not proceed, when we should use the NLP relaxation, but there is no NLP solver included in SCIP */
   if( heurdata->startsol == 'n' && SCIPgetNNlpis(scip) == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of branching nodes until heuristic is aborted */
   nstallnodes = (SCIP_Longint)(heurdata->nodesquot * SCIPgetNNodes(scip));

   /* reward RENS if it succeeded often */
   nstallnodes = (SCIP_Longint)(nstallnodes * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
   nstallnodes -= 100 * SCIPheurGetNCalls(heur);  /* count the setup costs for the sub-SCIP as 100 nodes */
   nstallnodes += heurdata->nodesofs;

   /* determine the node limit for the current process */
   nstallnodes -= heurdata->usednodes;
   nstallnodes = MIN(nstallnodes, heurdata->maxnodes);

   /* check whether we have enough nodes left to call subproblem solving */
   if( nstallnodes < heurdata->minnodes )
   {
      SCIPdebugMessage("skipping RENS: nstallnodes=%"SCIP_LONGINT_FORMAT", minnodes=%"SCIP_LONGINT_FORMAT"\n", nstallnodes, heurdata->minnodes);
      return SCIP_OKAY;
   }

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPapplyRens(scip, heur, result, heurdata->minfixingrate, heurdata->minimprove,
         heurdata->maxnodes, nstallnodes, heurdata->startsol, heurdata->binarybounds, heurdata->uselprows) );

   return SCIP_OKAY;
}



/*
 * primal heuristic specific interface methods
 */

/** creates the rens primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurRens(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyRens,
         heurFreeRens, heurInitRens, heurExitRens,
         heurInitsolRens, heurExitsolRens, heurExecRens,
         heurdata) );

   /* add rens primal heuristic parameters */

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minfixingrate",
         "minimum percentage of integer variables that have to be fixable ",
         &heurdata->minfixingrate, FALSE, DEFAULT_MINFIXINGRATE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/maxnodes",
         "maximum number of nodes to regard in the subproblem",
         &heurdata->maxnodes,  TRUE,DEFAULT_MAXNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/nodesofs",
         "number of nodes added to the contingent of the total nodes",
         &heurdata->nodesofs, FALSE, DEFAULT_NODESOFS, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/"HEUR_NAME"/minnodes",
         "minimum number of nodes required to start the subproblem",
         &heurdata->minnodes, TRUE, DEFAULT_MINNODES, 0LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nodesquot",
         "contingent of sub problem nodes in relation to the number of nodes of the original problem",
         &heurdata->nodesquot, FALSE, DEFAULT_NODESQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which RENS should at least improve the incumbent  ",
         &heurdata->minimprove, TRUE, DEFAULT_MINIMPROVE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip, "heuristics/"HEUR_NAME"/startsol",
         "solution that is used for fixing values ('l'p relaxation, 'n'lp relaxation)",
         &heurdata->startsol, FALSE, DEFAULT_STARTSOL, STARTSOL_CHOICES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/binarybounds",
         "should general integers get binary bounds [floor(.),ceil(.)] ?",
         &heurdata->binarybounds, TRUE, DEFAULT_BINARYBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/uselprows",
         "should subproblem be created out of the rows in the LP rows?",
         &heurdata->uselprows, TRUE, DEFAULT_USELPROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/copycuts",
         "if uselprows == FALSE, should all active cuts from cutpool be copied to constraints in subproblem?",
         &heurdata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );

   return SCIP_OKAY;
}
