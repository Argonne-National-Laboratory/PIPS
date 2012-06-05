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

/**@file   sepa_rapidlearning.c
 * @brief  rapidlearning separator
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#ifndef NDEBUG
#include <string.h>
#endif

#include "scip/sepa_rapidlearning.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_var.h"

#define SEPA_NAME              "rapidlearning"
#define SEPA_DESC               "rapid learning heuristic and separator"
#define SEPA_PRIORITY          -1200000
#define SEPA_FREQ                    -1 
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP           TRUE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_APPLYCONFLICTS     TRUE /**< should the found conflicts be applied in the original SCIP?                 */
#define DEFAULT_APPLYBDCHGS        TRUE /**< should the found global bound deductions be applied in the original SCIP?   
                                         *   apply only if conflicts and incumbent solution will be copied too
                                         */
#define DEFAULT_APPLYINFERVALS     TRUE /**< should the inference values be used as initialization in the original SCIP? */
#define DEFAULT_REDUCEDINFER      FALSE /**< should the inference values only be used when rapid learning found other reductions? */
#define DEFAULT_APPLYPRIMALSOL     TRUE /**< should the incumbent solution be copied to the original SCIP?               */
#define DEFAULT_APPLYSOLVED        TRUE /**< should a solved status be copied to the original SCIP?                      */

#define DEFAULT_MAXNVARS          10000 /**< maximum problem size (variables) for which rapid learning will be called */
#define DEFAULT_MAXNCONSS         10000 /**< maximum problem size (constraints) for which rapid learning will be called */

#define DEFAULT_MINNODES            500 /**< minimum number of nodes considered in rapid learning run */
#define DEFAULT_MAXNODES           5000 /**< maximum number of nodes considered in rapid learning run */

#define DEFAULT_CONTVARS          FALSE /**< should rapid learning be applied when there are continuous variables? */
#define DEFAULT_CONTVARSQUOT        0.3 /**< maximal portion of continuous variables to apply rapid learning       */
#define DEFAULT_LPITERQUOT          0.2 /**< maximal fraction of LP iterations compared to node LP iterations      */
#define DEFAULT_COPYCUTS           TRUE /**< should all active cuts from the cutpool of the
                                         *   original scip be copied to constraints of the subscip
                                         */


/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   SCIP_Bool             applyconflicts;     /**< should the found conflicts be applied in the original SCIP?                 */
   SCIP_Bool             applybdchgs;        /**< should the found global bound deductions be applied in the original SCIP?   */
   SCIP_Bool             applyinfervals;     /**< should the inference values be used as initialization in the original SCIP? */
   SCIP_Bool             reducedinfer;       /**< should the inference values only be used when rapid learning found other reductions? */
   SCIP_Bool             applyprimalsol;     /**< should the incumbent solution be copied to the original SCIP?               */
   SCIP_Bool             applysolved;        /**< should a solved status ba copied to the original SCIP?                      */
   int                   maxnvars;           /**< maximum problem size (variables) for which rapid learning will be called   */
   int                   maxnconss;          /**< maximum problem size (constraints) for which rapid learning will be called */
   int                   minnodes;           /**< minimum number of nodes considered in rapid learning run */
   int                   maxnodes;           /**< maximum number of nodes considered in rapid learning run */
   SCIP_Bool             contvars;           /**< should rapid learning be applied when there are continuous variables? */
   SCIP_Real             contvarsquot;       /**< maximal portion of continuous variables to apply rapid learning       */
   SCIP_Real             lpiterquot;         /**< maximal fraction of LP iterations compared to node LP iterations      */
   SCIP_Bool             copycuts;           /**< should all active cuts from cutpool be copied to constraints in
                                              *   subproblem?
                                              */
};

/** creates a new solution for the original problem by copying the solution of the subproblem */
static
SCIP_RETCODE createNewSol(
   SCIP*                 scip,               /**< original SCIP data structure                        */
   SCIP*                 subscip,            /**< SCIP structure of the subproblem                    */
   SCIP_VAR**            subvars,            /**< the variables of the subproblem                     */
   SCIP_HEUR*            heur,               /**< trysol heuristic structure                          */
   SCIP_SOL*             subsol,             /**< solution of the subproblem                          */
   SCIP_Bool*            success             /**< used to store whether new solution was found or not */
   )
{
   SCIP_VAR** vars;                          /* the original problem's variables                */
   int        nvars;
   SCIP_Real* subsolvals;                    /* solution values of the subproblem               */
   SCIP_SOL*  newsol;                        /* solution to be created for the original problem */
        
   assert( scip != NULL );
   assert( subscip != NULL );
   assert( subvars != NULL );
   assert( heur != NULL );
   assert( subsol != NULL );
   assert( success != NULL );

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

   /* check feasibility of new solution and pass it to trysol heuristic */
   SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, TRUE, TRUE, TRUE, success) );
   
   SCIPfreeBufferArray(scip, &subsolvals);

   return SCIP_OKAY;
}

/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyRapidlearning)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaRapidlearning(scip) );
 
   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeRapidlearning)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   assert(scip != NULL);

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);
   SCIPfreeMemory(scip, &sepadata);
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#define sepaInitRapidlearning NULL

/** deinitialization method of separator (called before transformed problem is freed) */
#define sepaExitRapidlearning NULL

/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolRapidlearning NULL

/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolRapidlearning NULL

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpRapidlearning)
{/*lint --e{715}*/
   SCIP* subscip;                            /* the subproblem created by rapid learning       */
   SCIP_SEPADATA* sepadata;                  /* separator's private data                       */

   SCIP_VAR** vars;                          /* original problem's variables                   */
   SCIP_VAR** subvars;                       /* subproblem's variables                         */
   SCIP_HASHMAP* varmapfw;                   /* mapping of SCIP variables to sub-SCIP variables */    
   SCIP_HASHMAP* varmapbw;                   /* mapping of sub-SCIP variables to SCIP variables */

   SCIP_CONSHDLR** conshdlrs;                /* array of constraint handler's that might that might obtain conflicts */
   int* oldnconss;                           /* number of constraints without rapid learning conflicts               */

   SCIP_Longint nodelimit;                   /* node limit for the subproblem                  */
   SCIP_Real timelimit;                      /* time limit for the subproblem                  */
   SCIP_Real memorylimit;                    /* memory limit for the subproblem                */

   int nconshdlrs;                           /* size of conshdlr and oldnconss array                      */
   int nfixedvars;                           /* number of variables that could be fixed by rapid learning */
   int nvars;                                /* number of variables                                       */           
   int restartnum;                           /* maximal number of conflicts that should be created        */
   int i;                                    /* counter                                                   */

   SCIP_Bool success;                        /* was problem creation / copying constraint successful? */
   SCIP_RETCODE retcode;                     /* used for catching sub-SCIP errors in debug mode */

   int nconflicts;                          /* statistic: number of conflicts applied         */
   int nbdchgs;                             /* statistic: number of bound changes applied     */
   int n1startinfers;                       /* statistic: number of one side infer values     */
   int n2startinfers;                       /* statistic: number of both side infer values    */

   SCIP_Bool soladded;                      /* statistic: was a new incumbent found?          */
   SCIP_Bool dualboundchg;                  /* statistic: was a new dual bound found?         */
   SCIP_Bool disabledualreductions;         /* TRUE, if dual reductions in sub-SCIP are not valid for original SCIP,
                                             * e.g., because a constraint could not be copied or a primal solution
                                             * could not be copied back 
                                             */

   int ndiscvars;

   soladded = FALSE;

   assert(sepa != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   
   ndiscvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip)+SCIPgetNImplVars(scip);

   /* only run when still not fixed binary variables exists */
   if( ndiscvars == 0 )
      return SCIP_OKAY;

   /* get separator's data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* only run for integer programs */
   if( !sepadata->contvars && ndiscvars != SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* only run if there are few enough continuous variables */
   if( sepadata->contvars && SCIPgetNContVars(scip) > sepadata->contvarsquot * SCIPgetNVars(scip) )
      return SCIP_OKAY;

   /* do not run if pricers are present */
   if( SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   /* if the separator should be exclusive to the root node, this prevents multiple calls due to restarts */
   if(  SCIPsepaGetFreq(sepa) == 0 && SCIPsepaGetNCalls(sepa) > 0)
      return SCIP_OKAY;

   /* call separator at most once per node */
   if( SCIPsepaGetNCallsAtNode(sepa) > 0 )
      return SCIP_OKAY;

   /* do not call rapid learning, if the problem is too big */
   if( SCIPgetNVars(scip) > sepadata->maxnvars || SCIPgetNConss(scip) > sepadata->maxnconss )
      return SCIP_OKAY; 

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;
   
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* initializing the subproblem */  
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) ); 
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPcalcHashtableSize(5 * nvars)) );
   success = FALSE;

   /* copy the subproblem */
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, "rapid", FALSE, FALSE, &success) );
   
   if( sepadata->copycuts )
   {
      /** copies all active cuts from cutpool of sourcescip to linear constraints in targetscip */
      SCIP_CALL( SCIPcopyCuts(scip, subscip, varmapfw, NULL, FALSE) );
   }

   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) (size_t) SCIPhashmapGetImage(varmapfw, vars[i]);
   
   SCIPhashmapFree(&varmapfw);
   
   /* this avoids dual presolving */
   if( !success )
   {
      for( i = 0; i < nvars; i++ )
      {     
         SCIP_CALL( SCIPaddVarLocks(subscip, subvars[i], 1, 1 ) );
      }
   }

   SCIPdebugMessage("Copying SCIP was%s successful.\n", success ? "" : " not");
   
   /* mimic an FD solver: DFS, no LP solving, 1-FUIP instead of all-FUIP */
   SCIP_CALL( SCIPsetIntParam(subscip, "lp/solvefreq", -1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "conflict/fuiplevels", 1) );
   SCIP_CALL( SCIPsetIntParam(subscip, "nodeselection/dfs/stdpriority", INT_MAX/4) ); 
   SCIP_CALL( SCIPsetBoolParam(subscip, "constraints/disableenfops", TRUE) );
   SCIP_CALL( SCIPsetIntParam(subscip, "propagating/pseudoobj/freq", -1) );

   /* use inference branching */
   SCIP_CALL( SCIPsetBoolParam(subscip, "branching/inference/useweightedsum", FALSE) );

   /* only create short conflicts */
   SCIP_CALL( SCIPsetRealParam(subscip, "conflict/maxvarsfac", 0.05) );
  
   /* set limits for the subproblem */
   nodelimit = SCIPgetNLPIterations(scip);
   nodelimit = MAX(sepadata->minnodes, nodelimit);
   nodelimit = MIN(sepadata->maxnodes, nodelimit);

   restartnum = 1000;
   
   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   SCIP_CALL( SCIPgetRealParam(scip, "limits/memory", &memorylimit) );
   if( !SCIPisInfinity(scip, memorylimit) )   
      memorylimit -= SCIPgetMemUsed(scip)/1048576.0;
   if( timelimit <= 0.0 || memorylimit <= 0.0 )
      goto TERMINATE;

   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit/5) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   SCIP_CALL( SCIPsetRealParam(subscip, "limits/memory", memorylimit) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/restarts", 0) );
   SCIP_CALL( SCIPsetIntParam(subscip, "conflict/restartnum", restartnum) );

   /* forbid recursive call of heuristics and separators solving subMIPs */
   SCIP_CALL( SCIPsetSubscipsOff(subscip, TRUE) );

   /* disable cutting plane separation */
   SCIP_CALL( SCIPsetSeparating(subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* disable expensive presolving */
   SCIP_CALL( SCIPsetPresolving(subscip, SCIP_PARAMSETTING_FAST, TRUE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

#ifndef SCIP_DEBUG
   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

   /* add an objective cutoff */
   SCIP_CALL( SCIPsetObjlimit(subscip, SCIPgetUpperbound(scip)) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapbw, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * nvars)) );

   /* store reversing mapping of variables */
   SCIP_CALL( SCIPtransformProb(subscip) );
   for( i = 0; i < nvars; ++i)
   {  
      SCIP_CALL( SCIPhashmapInsert(varmapbw, SCIPvarGetTransVar(subvars[i]), vars[i]) );
   }

   /** allocate memory for constraints storage. Each constraint that will be created from now on will be a conflict.
    *  Therefore, we need to remember oldnconss to get the conflicts from the FD search. 
    */
   nconshdlrs = 4;
   SCIP_CALL( SCIPallocBufferArray(scip, &conshdlrs, nconshdlrs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldnconss, nconshdlrs) );

   /* store number of constraints before rapid learning search */
   conshdlrs[0] = SCIPfindConshdlr(subscip, "bounddisjunction");
   conshdlrs[1] = SCIPfindConshdlr(subscip, "setppc");
   conshdlrs[2] = SCIPfindConshdlr(subscip, "linear");
   conshdlrs[3] = SCIPfindConshdlr(subscip, "logicor");

   /* redundant constraints might be eliminated in presolving */
   SCIP_CALL( SCIPpresolve(subscip));

   for( i = 0; i < nconshdlrs; ++i)
   {
      if( conshdlrs[i] != NULL )
         oldnconss[i] = SCIPconshdlrGetNConss(conshdlrs[i]);
   }

   nfixedvars = SCIPgetNFixedVars(scip);
   
   /* solve the subproblem */
   retcode = SCIPsolve(subscip);
   
   /* Errors in solving the subproblem should not kill the overall solving process 
    * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
    */
   if( retcode != SCIP_OKAY )
   { 
#ifndef NDEBUG
      SCIP_CALL( retcode );     
#endif
      SCIPwarningMessage("Error while solving subproblem in rapid learning separator; sub-SCIP terminated with code <%d>\n",retcode);
   }
 
   /* abort solving, if limit of applied conflicts is reached */
   if( SCIPgetNConflictConssApplied(subscip) >= restartnum )
   {
      SCIPdebugMessage("finish after %lld successful conflict calls.\n", SCIPgetNConflictConssApplied(subscip)); 
   }
   /* if the first 20% of the solution process were successful, proceed */
   else if( (sepadata->applyprimalsol && SCIPgetNSols(subscip) > 0 && SCIPisFeasLT(scip, SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip) ) )
      || (sepadata->applybdchgs && SCIPgetNFixedVars(subscip) > nfixedvars)
      || (sepadata->applyconflicts && SCIPgetNConflictConssApplied(subscip) > 0) ) 
   {
      SCIPdebugMessage("proceed solving after the first 20%% of the solution process, since:\n");

      if( SCIPgetNSols(subscip) > 0 && SCIPisFeasLE(scip, SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip) ) )
      {
         SCIPdebugMessage("   - there was a better solution (%f < %f)\n",SCIPgetUpperbound(subscip), SCIPgetUpperbound(scip));
      }
      if( SCIPgetNFixedVars(subscip) > nfixedvars )
      {
         SCIPdebugMessage("   - there were %d variables fixed\n", SCIPgetNFixedVars(scip)-nfixedvars );
      }
      if( SCIPgetNConflictConssFound(subscip) > 0 )
      {
         SCIPdebugMessage("   - there were %lld conflict constraints created\n", SCIPgetNConflictConssApplied(subscip));
      }

      /* set node limit to 100% */
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );

      /* solve the subproblem */
      retcode = SCIPsolve(subscip);
   
      /* Errors in solving the subproblem should not kill the overall solving process 
       * Hence, the return code is caught and a warning is printed, only in debug mode, SCIP will stop.
       */
      if( retcode != SCIP_OKAY )
      { 
#ifndef NDEBUG
         SCIP_CALL( retcode );     
#endif
         SCIPwarningMessage("Error while solving subproblem in rapid learning separator; sub-SCIP terminated with code <%d>\n",retcode);
      }
   }
   else
   {
      SCIPdebugMessage("do not proceed solving after the first 20%% of the solution process.\n");
   }

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

   disabledualreductions = FALSE;

   /* check, whether a solution was found */
   if( sepadata->applyprimalsol && SCIPgetNSols(subscip) > 0 && SCIPfindHeur(scip, "trysol") != NULL )
   {
      SCIP_HEUR* heurtrysol;
      SCIP_SOL** subsols;
      int nsubsols;

      /* check, whether a solution was found;
       * due to numerics, it might happen that not all solutions are feasible -> try all solutions until was declared to be feasible 
       */
      nsubsols = SCIPgetNSols(subscip);
      subsols = SCIPgetSols(subscip);
      soladded = FALSE;
      heurtrysol = SCIPfindHeur(scip, "trysol");

      /* sequentially add solutions to trysol heuristic */
      for( i = 0; i < nsubsols && !soladded; ++i )
      {
         SCIPdebugMessage("Try to create new solution by copying subscip solution.\n");
         SCIP_CALL( createNewSol(scip, subscip, subvars, heurtrysol, subsols[i], &soladded) );
      }
      if( !soladded || !SCIPisEQ(scip, SCIPgetSolOrigObj(subscip, subsols[i-1]), SCIPgetSolOrigObj(subscip, subsols[0])) )
         disabledualreductions = TRUE;
   }

   /* if the sub problem was solved completely, we update the dual bound */
   dualboundchg = FALSE;
   if( sepadata->applysolved && !disabledualreductions 
      && (SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL || SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE) )
   {
      /* we need to multiply the dualbound with the scaling factor and add the offset, 
       * because this information has been disregarded in the sub-SCIP */
      SCIPdebugMessage("Update old dualbound %g to new dualbound %g.\n", SCIPgetDualbound(scip), SCIPgetTransObjscale(scip) * SCIPgetDualbound(subscip) + SCIPgetTransObjoffset(scip));

      SCIP_CALL( SCIPupdateLocalDualbound(scip, SCIPgetDualbound(subscip) * SCIPgetTransObjscale(scip) + SCIPgetTransObjoffset(scip)) );
      dualboundchg = TRUE;
   }

   /* check, whether conflicts were created */
   nconflicts = 0;
   if( sepadata->applyconflicts && !disabledualreductions && SCIPgetNConflictConssApplied(subscip) > 0 )
   {
      SCIP_HASHMAP* consmap;
      int hashtablesize;

      assert(SCIPgetNConflictConssApplied(subscip) < (SCIP_Longint) INT_MAX);
      hashtablesize = (int) SCIPgetNConflictConssApplied(subscip);
      assert(hashtablesize < INT_MAX/5);
      hashtablesize *= 5;

      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&consmap, SCIPblkmem(scip), SCIPcalcHashtableSize(hashtablesize)) );

      /* loop over all constraint handlers that might contain conflict constraints */
      for( i = 0; i < nconshdlrs; ++i)
      {
         /* copy constraints that have been created in FD run */
         if( conshdlrs[i] != NULL && SCIPconshdlrGetNConss(conshdlrs[i]) > oldnconss[i] )
         {
            SCIP_CONS** conss;
            int c;
            int nconss;
            
            nconss = SCIPconshdlrGetNConss(conshdlrs[i]);
            conss = SCIPconshdlrGetConss(conshdlrs[i]);

            /* loop over all constraints that have been added in sub-SCIP run, these are the conflicts */            
            for( c = oldnconss[i]; c < nconss; ++c)
            {
               SCIP_CONS* cons;
               SCIP_CONS* conscopy;
               
               cons = conss[c];
               assert(cons != NULL);        

               success = FALSE;

               SCIP_CALL( SCIPgetConsCopy(subscip, scip, cons, &conscopy, conshdlrs[i], varmapbw, consmap, NULL,
                     SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons),
                     SCIPconsIsPropagated(cons), TRUE, FALSE, SCIPconsIsDynamic(cons), 
                     SCIPconsIsRemovable(cons), FALSE, TRUE, &success) );

               if( success )
               {
                  nconflicts++;
                  SCIP_CALL( SCIPaddCons(scip, conscopy) );
                  SCIP_CALL( SCIPreleaseCons(scip, &conscopy) );
               }
               else
               {
                  SCIPdebugMessage("failed to copy conflict constraint %s back to original SCIP\n", SCIPconsGetName(cons));
               }
            }
         }
      }   
      SCIPhashmapFree(&consmap);
   }

   /* check, whether tighter global bounds were detected */
   nbdchgs = 0;
   if( sepadata->applybdchgs && !disabledualreductions )
      for( i = 0; i < nvars; ++i )
      {
         SCIP_Bool infeasible;
         SCIP_Bool tightened;
         
         assert(SCIPisLE(scip, SCIPvarGetLbGlobal(vars[i]), SCIPvarGetLbGlobal(subvars[i]))); 
         assert(SCIPisLE(scip, SCIPvarGetLbGlobal(subvars[i]), SCIPvarGetUbGlobal(subvars[i])));
         assert(SCIPisLE(scip, SCIPvarGetUbGlobal(subvars[i]), SCIPvarGetUbGlobal(vars[i])));  
         
         /* update the bounds of the original SCIP, if a better bound was proven in the sub-SCIP */
         SCIP_CALL( SCIPtightenVarUb(scip, vars[i], SCIPvarGetUbGlobal(subvars[i]), FALSE, &infeasible, &tightened) );
         if( tightened ) 
            nbdchgs++;
         
         SCIP_CALL( SCIPtightenVarLb(scip, vars[i], SCIPvarGetLbGlobal(subvars[i]), FALSE, &infeasible, &tightened) );
         if( tightened )
            nbdchgs++;   
      }

   n1startinfers = 0;
   n2startinfers = 0;

   /* install start values for inference branching */
   if( sepadata->applyinfervals && (!sepadata->reducedinfer || soladded || nbdchgs+nconflicts > 0) )
   {
      for( i = 0; i < nvars; ++i )
      {
         SCIP_Real downinfer;
         SCIP_Real upinfer;
         SCIP_Real downvsids;
         SCIP_Real upvsids;
         SCIP_Real downconflen;
         SCIP_Real upconflen;
        
         /* copy downwards branching statistics */
         downvsids = SCIPgetVarVSIDS(subscip, subvars[i], SCIP_BRANCHDIR_DOWNWARDS);            
         downconflen = SCIPgetVarAvgConflictlength(subscip, subvars[i], SCIP_BRANCHDIR_DOWNWARDS);
         downinfer = SCIPgetVarAvgInferences(subscip, subvars[i], SCIP_BRANCHDIR_DOWNWARDS);            
         
         /* copy upwards branching statistics */
         upvsids = SCIPgetVarVSIDS(subscip, subvars[i], SCIP_BRANCHDIR_UPWARDS);                     
         upconflen = SCIPgetVarAvgConflictlength(subscip, subvars[i], SCIP_BRANCHDIR_UPWARDS);
         upinfer = SCIPgetVarAvgInferences(subscip, subvars[i], SCIP_BRANCHDIR_UPWARDS);            
        
         /* memorize statistics */
         if( downinfer+downconflen+downvsids > 0.0 || upinfer+upconflen+upvsids != 0 )
            n1startinfers++;
         
         if( downinfer+downconflen+downvsids > 0.0 && upinfer+upconflen+upvsids != 0 )
            n2startinfers++;
         
         SCIP_CALL( SCIPinitVarBranchStats(scip, vars[i], 0.0, 0.0, downvsids, upvsids, downconflen, upconflen, downinfer, upinfer, 0.0, 0.0) );
      }   
   }
   
   SCIPdebugPrintf("XXX Rapidlearning added %d conflicts, changed %d bounds, %s primal solution, %s dual bound improvement.\n", nconflicts, nbdchgs, soladded ? "found" : "no", 
      dualboundchg ? "found" : "no");

   SCIPdebugPrintf("YYY Infervalues initialized on one side: %5.2f %% of variables, %5.2f %% on both sides\n", 
      100.0 * n1startinfers/(SCIP_Real)nvars, 100.0 * n2startinfers/(SCIP_Real)nvars);

   /* change result pointer */
   if( nconflicts > 0 || dualboundchg )
      *result = SCIP_CONSADDED;
   else if( nbdchgs > 0 )
      *result = SCIP_REDUCEDDOM;
  
   /* free local data */
   SCIPfreeBufferArray(scip, &oldnconss);
   SCIPfreeBufferArray(scip, &conshdlrs);

   SCIPhashmapFree(&varmapbw);

 TERMINATE:
   /* free subproblem */
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );
  
   return SCIP_OKAY;
}

/** arbitrary primal solution separation method of separator */
#define sepaExecsolRapidlearning NULL

/*
 * separator specific interface methods
 */

/** creates the rapidlearning separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaRapidlearning(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create rapidlearning separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyRapidlearning, sepaFreeRapidlearning, sepaInitRapidlearning, sepaExitRapidlearning, 
         sepaInitsolRapidlearning, sepaExitsolRapidlearning,
         sepaExeclpRapidlearning, sepaExecsolRapidlearning,
         sepadata) );

   /* add rapidlearning separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/applyconflicts",
         "should the found conflicts be applied in the original SCIP?",
         &sepadata->applyconflicts, TRUE, DEFAULT_APPLYCONFLICTS, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/applybdchgs",
         "should the found global bound deductions be applied in the original SCIP?",
         &sepadata->applybdchgs, TRUE, DEFAULT_APPLYBDCHGS, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/applyinfervals",
         "should the inference values be used as initialization in the original SCIP?",
         &sepadata->applyinfervals, TRUE, DEFAULT_APPLYINFERVALS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/reducedinfer",
         "should the inference values only be used when "SEPA_NAME" found other reductions?",
         &sepadata->reducedinfer, TRUE, DEFAULT_REDUCEDINFER, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/applyprimalsol",
         "should the incumbent solution be copied to the original SCIP?",
         &sepadata->applyprimalsol, TRUE, DEFAULT_APPLYPRIMALSOL, NULL, NULL) );
  
   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/applysolved",
         "should a solved status be copied to the original SCIP?",
         &sepadata->applysolved, TRUE, DEFAULT_APPLYSOLVED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/contvars",
         "should rapid learning be applied when there are continuous variables?",
         &sepadata->contvars, TRUE, DEFAULT_CONTVARS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/"SEPA_NAME"/contvarsquot",
         "maximal portion of continuous variables to apply rapid learning",
         &sepadata->contvarsquot, TRUE, DEFAULT_CONTVARSQUOT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/"SEPA_NAME"/lpiterquot",
         "maximal fraction of LP iterations compared to node LP iterations",
         &sepadata->lpiterquot, TRUE, DEFAULT_LPITERQUOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );
 
   SCIP_CALL( SCIPaddIntParam(scip, "separating/"SEPA_NAME"/maxnvars",
         "maximum problem size (variables) for which rapid learning will be called",
         &sepadata->maxnvars, TRUE, DEFAULT_MAXNVARS, 0, INT_MAX, NULL, NULL) );
   
   SCIP_CALL( SCIPaddIntParam(scip, "separating/"SEPA_NAME"/maxnconss",
         "maximum problem size (constraints) for which rapid learning will be called",
         &sepadata->maxnconss, TRUE, DEFAULT_MAXNCONSS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/"SEPA_NAME"/maxnodes",
         "maximum number of nodes considered in rapid learning run",
         &sepadata->maxnodes, TRUE, DEFAULT_MAXNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/"SEPA_NAME"/minnodes",
         "minimum number of nodes considered in rapid learning run",
         &sepadata->minnodes, TRUE, DEFAULT_MINNODES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/"SEPA_NAME"/copycuts",
         "should all active cuts from cutpool be copied to constraints in subproblem?",
         &sepadata->copycuts, TRUE, DEFAULT_COPYCUTS, NULL, NULL) );
   
   return SCIP_OKAY;
}
