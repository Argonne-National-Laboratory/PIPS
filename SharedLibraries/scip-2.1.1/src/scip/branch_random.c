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

/**@file   branch_random.c
 * @brief  random variable branching rule
 * @author Tobias Achterberg
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/branch_random.h"


#define BRANCHRULE_NAME          "random"
#define BRANCHRULE_DESC          "random variable branching"
#define BRANCHRULE_PRIORITY      -100000
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

#define DEFAULT_INITSEED                0 /**< initial random seed */

/** branching rule data */
struct SCIP_BranchruleData
{
   int                   initseed;           /**< initial random seed value */
   unsigned int          seed;               /**< seed value for random number generator */
};

/*
 * Local methods
 */

/** selects a random active variable from a given list of variables */
static
void getRandomVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            cands,              /**< array of branching candidates */
   SCIP_Real*            candssol,           /**< relaxation solution values of branching candidates, or NULL */
   int                   ncands,             /**< number of branching candidates */
   SCIP_VAR**            bestcand,           /**< buffer to store pointer to best candidate */
   SCIP_Real*            bestcandsol,        /**< buffer to store solution value of best candidate */
   unsigned int*         seed                /**< seed for random number generator */
   )
{
   int idx;
   int firstidx;
   
   assert(scip != NULL);
   assert(cands != NULL);
   assert(ncands > 0);
   assert(bestcand != NULL);
   assert(bestcandsol != NULL);
   assert(seed != NULL);
 
   idx = SCIPgetRandomInt(0, ncands-1, seed);
   assert(idx >= 0);
   
   /* handle case where cands[idx] is fixed by selecting next idx with unfixed var
    * this may happen if we are inside a multi-aggregation */
   firstidx = idx;
   while( SCIPisEQ(scip, SCIPvarGetLbLocal(cands[idx]), SCIPvarGetUbLocal(cands[idx])) )
   {
      ++idx;
      if( idx == ncands )
         idx = 0;
      if( idx == firstidx )
      {
         /* odd: all variables seem to be fixed */
         SCIPdebugMessage("Warning: all branching candidates seem to be fixed\n");
         return;
      }
   }
   
   /* a branching variable candidate should either be an active problem variable or a multi-aggregated variable */
   assert(SCIPvarIsActive(SCIPvarGetProbvar(cands[idx])) ||
      SCIPvarGetStatus(SCIPvarGetProbvar(cands[idx])) == SCIP_VARSTATUS_MULTAGGR);
   
   if( SCIPvarGetStatus(SCIPvarGetProbvar(cands[idx])) == SCIP_VARSTATUS_MULTAGGR )
   {
      /* for a multi-aggregated variable, we call the getRandomVariable function recursively with all variables in the multi-aggregation */
      SCIP_VAR* cand;
      
      cand = SCIPvarGetProbvar(cands[idx]);
      
      getRandomVariable(scip, SCIPvarGetMultaggrVars(cand), NULL, SCIPvarGetMultaggrNVars(cand), bestcand, bestcandsol, seed);
      return;
   }
   
   assert(idx >= 0 && idx < ncands);
   
   *bestcand = cands[idx];
   assert(*bestcand != NULL);
   
   if( candssol != NULL )
      *bestcandsol = candssol[idx];
}

/*
 * Callback methods
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
static
SCIP_DECL_BRANCHCOPY(branchCopyRandom)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);

   /* call inclusion method of branchrule */
   SCIP_CALL( SCIPincludeBranchruleRandom(scip) );

   return SCIP_OKAY;
}

/** destructor of branching rule to free user data (called when SCIP is exiting) */
static
SCIP_DECL_BRANCHFREE(branchFreeRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   /* free branching rule data */
   branchruledata = SCIPbranchruleGetData(branchrule);
   SCIPfreeMemory(scip, &branchruledata);
   SCIPbranchruleSetData(branchrule, NULL);

   return SCIP_OKAY;
}


/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* set the seed value to the initial random seed value */
   branchruledata->seed = (unsigned int) branchruledata->initseed;
   
   return SCIP_OKAY;
}


/** deinitialization method of branching rule (called before transformed problem is freed) */
#define branchExitRandom NULL


/** solving process initialization method of branching rule (called when branch and bound process is about to begin) */
#define branchInitsolRandom NULL


/** solving process deinitialization method of branching rule (called before branch and bound process data is freed) */
#define branchExitsolRandom NULL


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** lpcands;
   int nlpcands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execlp method of random branching in depth %d\n", SCIPgetDepth(scip));

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, NULL, NULL, &nlpcands) );
   assert(nlpcands > 0);

   /* get random branching candidate */
   bestcand = SCIPgetRandomInt(0, nlpcands-1, &branchruledata->seed);
   assert(bestcand >= 0);

   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s>\n",
      nlpcands, bestcand, SCIPvarGetName(lpcands[bestcand]));

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, lpcands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}


/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** externcands;
   SCIP_Real* externcandssol;
   int nprioexterncands;
   SCIP_VAR* bestcand;
   SCIP_Real bestcandsol;
   SCIP_Real brpoint;
   SCIP_NODE* downchild;
   SCIP_NODE* eqchild;
   SCIP_NODE* upchild;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execrel method of random branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   bestcand = NULL;
   bestcandsol = 0.0;

   /* get branching candidates */
   SCIP_CALL( SCIPgetExternBranchCands(scip, &externcands, &externcandssol, NULL, NULL, &nprioexterncands, NULL, NULL, NULL) );
   assert(nprioexterncands > 0);

   /* get random branching candidate
    *
    * since variables can occur several times in the list of candidates, variables that have been added more often have
    * a higher probability to be chosen for branching
    */
   getRandomVariable(scip, externcands, externcandssol, nprioexterncands, &bestcand, &bestcandsol, &branchruledata->seed);

   if( bestcand == NULL )
   {
      SCIPerrorMessage("branchExecrelRandom failed to select a branching variable from %d candidates\n", nprioexterncands);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   brpoint = SCIPgetBranchingPoint(scip, bestcand, bestcandsol);

   SCIPdebugMessage(" -> %d candidates, selected variable <%s> with solution value %g, branching point=%g\n",
      nprioexterncands, SCIPvarGetName(bestcand), bestcandsol, brpoint);

   SCIP_CALL( SCIPbranchVarVal(scip, bestcand, brpoint, &downchild, &eqchild, &upchild) );

   if( downchild != NULL || eqchild != NULL || upchild != NULL )
   {
      *result = SCIP_BRANCHED;
   }
   else
   {
      /* if there are no children, then variable should have been fixed by SCIPbranchVarVal */
      assert(SCIPisEQ(scip, SCIPvarGetLbLocal(bestcand), SCIPvarGetUbLocal(bestcand)));
      *result = SCIP_REDUCEDDOM;
   }
   
   return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsRandom)
{  /*lint --e{715}*/
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_VAR** pseudocands;
   int npseudocands;
   int bestcand;

   assert(branchrule != NULL);
   assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
   assert(scip != NULL);
   assert(result != NULL);

   SCIPdebugMessage("Execps method of random branching\n");

   branchruledata = SCIPbranchruleGetData(branchrule);
   assert(branchruledata != NULL);

   /* get branching candidates */
   SCIP_CALL( SCIPgetPseudoBranchCands(scip, &pseudocands, NULL, &npseudocands) );
   assert(npseudocands > 0);

   /* get random branching candidate */
   bestcand = SCIPgetRandomInt(0, npseudocands-1, &branchruledata->seed);
   assert(bestcand >= 0);

   SCIPdebugMessage(" -> %d candidates, selected candidate %d: variable <%s>\n",
      npseudocands, bestcand, SCIPvarGetName(pseudocands[bestcand]));

   /* perform the branching */
   SCIP_CALL( SCIPbranchVar(scip, pseudocands[bestcand], NULL, NULL, NULL) );
   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}





/*
 * branching specific interface methods
 */

/** creates the random branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRandom(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;

   /* create random branching rule data */
   SCIP_CALL( SCIPallocMemory(scip, &branchruledata) );
   branchruledata->seed = 0;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST,
         branchCopyRandom,
         branchFreeRandom, branchInitRandom, branchExitRandom, branchInitsolRandom, branchExitsolRandom, 
         branchExeclpRandom, branchExecextRandom, branchExecpsRandom,
         branchruledata) );

   SCIP_CALL( SCIPaddIntParam(scip, "branching/"BRANCHRULE_NAME"/seed", "initial random seed value",
         &branchruledata->initseed, FALSE, DEFAULT_INITSEED, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
