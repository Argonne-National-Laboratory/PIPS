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

/**@file   prop_probing.c
 * @brief  probing propagator
 * @author Tobias Achterberg
 * @author Matthias Miltenberger
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_probing.h"


#define PROP_NAME               "probing"
#define PROP_DESC               "probing propagator on binary variables"
#define PROP_TIMING             SCIP_PROPTIMING_AFTERLPLOOP
#define PROP_PRIORITY           -100000 /**< propagation priority */
#define PROP_FREQ                    -1 /**< propagation frequency */
#define PROP_DELAY                 TRUE /**< should propagation method be delayed, if other propagators found
                                         *   reductions? */
#define PROP_PRESOL_PRIORITY    -100000 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOL_DELAY          TRUE /**< should presolving be delay, if other presolvers found reductions?  */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */

#define MAXDNOM                 10000LL /**< maximal denominator for simple rational fixed values */



/*
 * Default parameter settings
 */

#define DEFAULT_MAXRUNS              1  /**< maximal number of runs, probing participates in (-1: no limit) */
#define DEFAULT_PROPROUNDS          -1  /**< maximal number of propagation rounds in probing subproblems */
#define DEFAULT_MAXFIXINGS          25  /**< maximal number of fixings found, until probing is interrupted
                                         *   (0: don't interrupt) */
#define DEFAULT_MAXUSELESS        1000  /**< maximal number of successive probings without fixings,
                                         *   until probing is aborted (0: don't abort) */
#define DEFAULT_MAXTOTALUSELESS     50  /**< maximal number of successive probings without fixings, bound changes,
                                         *   and implications, until probing is aborted (0: don't abort) */
#define DEFAULT_MAXSUMUSELESS        0  /**< maximal number of probings without fixings, until probing is aborted
                                         *   (0: don't abort) */



/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_VAR**            sortedvars;         /**< problem variables sorted by number of rounding locks, used in presolving */
   int                   nsortedvars;        /**< number of problem variables, used in presolving */
   int                   nsortedbinvars;     /**< number of binary problem variables, used in presolving */
   int                   maxruns;            /**< maximal number of runs, probing participates in (-1: no limit) */
   int                   proprounds;         /**< maximal number of propagation rounds in probing subproblems */
   int                   maxfixings;         /**< maximal number of fixings found, until probing is interrupted
                                              *   (0: don't interrupt) */
   int                   maxuseless;         /**< maximal number of successive probings without fixings,
                                              *   until probing is aborted (0: don't abort) */
   int                   maxtotaluseless;    /**< maximal number of successive probings without fixings, bound changes,
                                              *   and implications, until probing is aborted (0: don't abort) */
   int                   maxsumuseless;      /**< maximal number of probings without fixings, until probing is aborted
                                              *   (0: don't abort) */
   int                   startidx;           /**< starting variable index of next call, used in presolving */
   int                   lastsortstartidx;   /**< last starting variable index where the variables have been sorted, used in presolving */
   int                   nfixings;           /**< total number of fixings found in probing */
   int                   naggregations;      /**< total number of aggregations found in probing */
   int                   nimplications;      /**< total number of implications found in probing */
   int                   nbdchgs;            /**< total number of bound changes found in probing */
   int                   nuseless;           /**< current number of successive useless probings */
   int                   ntotaluseless;      /**< current number of successive totally useless probings */
   int                   nsumuseless;        /**< current number of useless probings */
   SCIP_Longint          lastnode;           /**< last node where probing was applied, or -1 for presolving, and -2 for not applied yet */
};




/*
 * Local methods
 */
/** initializes the propagator data */
static
void initPropdata(
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   propdata->sortedvars = NULL;
   propdata->nsortedvars = 0;
   propdata->nsortedbinvars = 0;
   propdata->startidx = 0;
   propdata->lastsortstartidx = -1;
   propdata->nfixings = 0;
   propdata->naggregations = 0;
   propdata->nimplications = 0;
   propdata->nbdchgs = 0;
   propdata->nuseless = 0;
   propdata->ntotaluseless = 0;
   propdata->nsumuseless = 0;
   propdata->lastnode = -2;
}

/** frees the sorted vars array */
static
SCIP_RETCODE freeSortedvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   assert(propdata != NULL);

   if( propdata->sortedvars != NULL )
   {
      int i;

      /* release variables */
      for( i = 0; i < propdata->nsortedvars; ++i )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &propdata->sortedvars[i]) );
      }
      SCIPfreeMemoryArray(scip, &propdata->sortedvars);
      propdata->nsortedvars = 0;
      propdata->nsortedbinvars = 0;
   }

   return SCIP_OKAY;
}

/** sorts the binary variables starting with the given index by rounding locks and implications */
static
SCIP_RETCODE sortVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< problem variables to be sorted */
   int                   nvars,              /**< number of problem variables to be sorted */
   int                   firstidx            /**< first index that should be subject to sorting */
   )
{
   SCIP_VAR** sortedvars;
   int nsortedvars;
   int* scores;
   int i;

   assert(vars != NULL || nvars == 0);

   nsortedvars = nvars - firstidx;
   if( nsortedvars <= 0 )
      return SCIP_OKAY;

   assert(vars != NULL);

   sortedvars = &(vars[firstidx]);

   SCIPdebugMessage("resorting probing variables %d to %d\n", firstidx, nvars-1);

   /* sort the variables by number of rounding locks and implications */
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, nsortedvars) );

   for( i = 0; i < nsortedvars; ++i )
   {
      SCIP_VAR* var;
      var = sortedvars[i];

      assert(SCIPvarIsBinary(var));
      if( SCIPvarIsActive(var) )
         scores[i] = SCIPvarGetNLocksDown(var) + SCIPvarGetNLocksUp(var)
            + SCIPvarGetNImpls(var, FALSE) + SCIPvarGetNImpls(var, TRUE)
            + 5*(SCIPvarGetNCliques(var, FALSE) + SCIPvarGetNCliques(var, TRUE));
      else
         scores[i] = -1;
   }

   SCIPsortDownIntPtr(scores, (void**) sortedvars, nsortedvars);

   SCIPfreeBufferArray(scip, &scores);

   return SCIP_OKAY;
}

/** applies and evaluates probing of a single variable in the given direction */
static
SCIP_RETCODE applyProbingVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   probingpos,         /**< variable number to apply probing on */
   SCIP_Bool             probingdir,         /**< value to fix probing variable to */
   SCIP_Real*            impllbs,            /**< array to store lower bounds after applying implications and cliques */
   SCIP_Real*            implubs,            /**< array to store upper bounds after applying implications and cliques */
   SCIP_Real*            proplbs,            /**< array to store lower bounds after full propagation */
   SCIP_Real*            propubs,            /**< array to store upper bounds after full propagation */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing direction is infeasible */
   )
{
   assert(propdata != NULL);
   assert(impllbs != NULL);
   assert(implubs != NULL);
   assert(proplbs != NULL);
   assert(propubs != NULL);
   assert(cutoff != NULL);
   assert(0 <= probingpos && probingpos < nvars);
   assert(SCIPvarIsBinary(vars[probingpos]));
   assert(SCIPvarGetLbLocal(vars[probingpos]) < 0.5);
   assert(SCIPvarGetUbLocal(vars[probingpos]) > 0.5);

   SCIPdebugMessage("applying probing on variable <%s> == %u (nlocks=%d/%d, impls=%d/%d, clqs=%d/%d)\n",
      SCIPvarGetName(vars[probingpos]), probingdir,
      SCIPvarGetNLocksDown(vars[probingpos]), SCIPvarGetNLocksUp(vars[probingpos]),
      SCIPvarGetNImpls(vars[probingpos], FALSE), SCIPvarGetNImpls(vars[probingpos], TRUE),
      SCIPvarGetNCliques(vars[probingpos], FALSE), SCIPvarGetNCliques(vars[probingpos], TRUE));

   /* start probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* fix variable */
   if( probingdir == FALSE )
   {
      SCIP_CALL( SCIPchgVarUbProbing(scip, vars[probingpos], 0.0) );
   }
   else
   {
      SCIP_CALL( SCIPchgVarLbProbing(scip, vars[probingpos], 1.0) );
   }

   /* apply propagation of implication graph and clique table */
   SCIP_CALL( SCIPpropagateProbingImplications(scip, cutoff) );
   if( !(*cutoff) )
   {
      int i;

      for( i = 0; i < nvars; ++i )
      {
         impllbs[i] = SCIPvarGetLbLocal(vars[i]);
         implubs[i] = SCIPvarGetUbLocal(vars[i]);
      }
   }

   /* apply propagation */
   if( !(*cutoff) )
   {
      SCIP_CALL( SCIPpropagateProbing(scip, propdata->proprounds, cutoff, NULL) );
   }

   /* evaluate propagation */
   if( !(*cutoff) )
   {
      int i;

      for( i = 0; i < nvars; ++i )
      {
         proplbs[i] = SCIPvarGetLbLocal(vars[i]);
         propubs[i] = SCIPvarGetUbLocal(vars[i]);
#if 0
#ifdef SCIP_DEBUG
         if( SCIPisGT(scip, proplbs[i], SCIPvarGetLbGlobal(vars[i])) )
         {
            SCIPdebugMessage(" -> <%s>[%g,%g] >= %g\n", SCIPvarGetName(vars[i]),
               SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]), proplbs[i]);
         }
         if( SCIPisLT(scip, propubs[i], SCIPvarGetUbGlobal(vars[i])) )
         {
            SCIPdebugMessage(" -> <%s>[%g,%g] <= %g\n", SCIPvarGetName(vars[i]),
               SCIPvarGetLbGlobal(vars[i]), SCIPvarGetUbGlobal(vars[i]), propubs[i]);
         }
#endif
#endif
      }
   }

   /* exit probing mode */
   SCIP_CALL( SCIPendProbing(scip) );

   return SCIP_OKAY;
}


/** the main probing loop */
static
SCIP_RETCODE applyProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR**            vars,               /**< problem variables */
   int                   nvars,              /**< number of problem variables */
   int                   nbinvars,           /**< number of binary variables */
   int*                  startidx,           /**< pointer to store starting variable index of next call */
   int*                  nfixedvars,         /**< pointer to store number of fixed variables */
   int*                  naggrvars,          /**< pointer to store number of aggregated variables */
   int*                  nchgbds,            /**< pointer to store number of changed bounds */
   int                   oldnfixedvars,      /**< number of previously fixed variables */
   int                   oldnaggrvars,       /**< number of previously aggregated variables */
   SCIP_Bool*            delay,              /**< pointer to store whether propagator should be delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether cutoff occured */
   )
{
   SCIP_Real* zeroimpllbs;
   SCIP_Real* zeroimplubs;
   SCIP_Real* zeroproplbs;
   SCIP_Real* zeropropubs;
   SCIP_Real* oneimpllbs;
   SCIP_Real* oneimplubs;
   SCIP_Real* oneproplbs;
   SCIP_Real* onepropubs;
   int localnfixedvars;
   int localnaggrvars;
   int localnchgbds;
   int localnimplications;
   int maxfixings;
   int maxuseless;
   int maxtotaluseless;
   int maxsumuseless;
   int i;

   assert(vars != NULL);
   assert(nbinvars > 0);

   maxfixings = (propdata->maxfixings > 0 ? propdata->maxfixings : INT_MAX);
   maxuseless = (propdata->maxuseless > 0 ? propdata->maxuseless : INT_MAX);
   maxtotaluseless = (propdata->maxtotaluseless > 0 ? propdata->maxtotaluseless : INT_MAX);
   maxsumuseless = (propdata->maxsumuseless > 0 ? propdata->maxsumuseless : INT_MAX);

   /* get temporary memory for storing probing results */
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimpllbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroimplubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeroproplbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &zeropropubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimpllbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneimplubs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oneproplbs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &onepropubs, nvars) );

   /* for each binary variable, probe fixing the variable to zero and one */
   *delay = FALSE;
   *cutoff = FALSE;
   for( i = *startidx; i < nbinvars && !(*cutoff); ++i )
   {
      SCIP_Bool localcutoff;
      SCIP_Bool probingzero;
      SCIP_Bool probingone;

      /* check whether probing should be aborted */
      if( SCIPisStopped(scip) || propdata->nuseless >= maxuseless
         || propdata->ntotaluseless >= maxtotaluseless || propdata->nsumuseless >= maxsumuseless )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) probing: %d/%d (%.1f%%) - %d fixings, %d aggregations, %d implications, %d bound changes\n",
            SCIPgetSolvingTime(scip), i+1, nbinvars, 100.0*(SCIP_Real)(i+1)/(SCIP_Real)nbinvars,
            propdata->nfixings, propdata->naggregations, propdata->nimplications, propdata->nbdchgs);

         if( SCIPisStopped(scip) )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: solving stopped\n", SCIPgetSolvingTime(scip));
         }
         else if( propdata->nuseless >= maxuseless )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: %d/%d successive useless probings\n", SCIPgetSolvingTime(scip),
               propdata->nuseless, maxuseless);
         }
         else if( propdata->ntotaluseless >= maxtotaluseless )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: %d/%d successive totally useless probings\n", SCIPgetSolvingTime(scip),
               propdata->ntotaluseless, maxtotaluseless);
         }
         else if( propdata->nsumuseless >= maxsumuseless )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
               "   (%.1fs) probing aborted: %d/%d useless probings in total\n", SCIPgetSolvingTime(scip),
               propdata->nsumuseless, maxsumuseless);
         }
         break;
      }

      /* check if we already fixed enough variables for this round */
      if( *nfixedvars - oldnfixedvars + *naggrvars - oldnaggrvars >= maxfixings )
      {
         *delay = TRUE;
         break;
      }

      /* display probing status */
      if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING && (i+1) % 100 == 0 )
      {
         SCIP_VERBLEVEL verblevel;

         verblevel = ((i+1) % 1000 == 0 ? SCIP_VERBLEVEL_HIGH : SCIP_VERBLEVEL_FULL);
         SCIPverbMessage(scip, verblevel, NULL,
            "   (%.1fs) probing: %d/%d (%.1f%%) - %d fixings, %d aggregations, %d implications, %d bound changes\n",
            SCIPgetSolvingTime(scip), i+1, nbinvars, 100.0*(SCIP_Real)(i+1)/(SCIP_Real)nbinvars,
            propdata->nfixings, propdata->naggregations, propdata->nimplications, propdata->nbdchgs);
      }

      /* ignore variables, that were fixed, aggregated, or deleted in prior probings */
      if( !SCIPvarIsActive(vars[i]) || SCIPvarIsDeleted(vars[i])
         || SCIPvarGetLbLocal(vars[i]) > 0.5 || SCIPvarGetUbLocal(vars[i]) < 0.5 )
         continue;

      if( propdata->nuseless > 0 )
         propdata->nsumuseless++;
      else
         propdata->nsumuseless = MAX(propdata->nsumuseless-1, 0);
      propdata->nuseless++;
      propdata->ntotaluseless++;

      /* determine whether zero probing should happen */
      probingzero = TRUE;
      if( SCIPvarGetNLocksDown(vars[i]) == 0 )
         probingzero = FALSE;

      if( probingzero )
      {
         /* apply probing for fixing the variable to zero */
         SCIP_CALL( applyProbingVar(scip, propdata, vars, nvars, i, FALSE, zeroimpllbs, zeroimplubs, zeroproplbs, zeropropubs,
               &localcutoff) );

         if( localcutoff )
         {
            SCIP_Bool fixed;

            if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
            {
               /* the variable can be fixed to TRUE */
               SCIP_CALL( SCIPfixVar(scip, vars[i], 1.0, cutoff, &fixed) );
            }
            else
            {
               SCIP_CALL( SCIPtightenVarLb(scip, vars[i], 1.0, TRUE, cutoff, &fixed) );
            }

            if( fixed )
            {
               SCIPdebugMessage("fixed probing variable <%s> to 1.0, nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
               (*nfixedvars)++;
               propdata->nfixings++;
               propdata->nuseless = 0;
               propdata->ntotaluseless = 0;
            }
            continue; /* don't try upwards direction, because the variable is already fixed */
         }

         /* ignore variables, that were fixed, aggregated, or deleted in prior probings
          * (propagators in zero-probe might have found global fixings but did not trigger the localcutoff)
          */
         if( !SCIPvarIsActive(vars[i]) || SCIPvarIsDeleted(vars[i])
            || SCIPvarGetLbLocal(vars[i]) > 0.5 || SCIPvarGetUbLocal(vars[i]) < 0.5 )
            continue;
      }

      /* determine whether one probing should happen */
      probingone = TRUE;
      if( SCIPvarGetNLocksUp(vars[i]) == 0 )
         probingone = FALSE;

      if( probingone )
      {
         /* apply probing for fixing the variable to one */
         SCIP_CALL( applyProbingVar(scip, propdata, vars, nvars, i, TRUE, oneimpllbs, oneimplubs, oneproplbs, onepropubs,
               &localcutoff) );

         if( localcutoff )
         {
            SCIP_Bool fixed;

            if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
            {
               /* the variable can be fixed to FALSE */
               SCIP_CALL( SCIPfixVar(scip, vars[i], 0.0, cutoff, &fixed) );
               assert(fixed);
            }
            else
            {
               SCIP_CALL( SCIPtightenVarUb(scip, vars[i], 0.0, TRUE, cutoff, &fixed) );
            }

            if( fixed )
            {
               SCIPdebugMessage("fixed probing variable <%s> to 0.0, nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[i]), SCIPvarGetNLocksDown(vars[i]), SCIPvarGetNLocksUp(vars[i]));
               (*nfixedvars)++;
               propdata->nfixings++;
               propdata->nuseless = 0;
               propdata->ntotaluseless = 0;
            }
            continue; /* don't analyze probing deductions, because the variable is already fixed */
         }
      }

      /* not have to check deductions if only one probing direction has been checked */
      if( !probingzero || !probingone )
         continue;

      /* analyze probing deductions */
      localnfixedvars    = 0;
      localnaggrvars     = 0;
      localnimplications = 0;
      localnchgbds       = 0;
      SCIP_CALL( SCIPanalyzeDeductionsProbing(scip, vars[i], 0.0, 1.0,
         nvars, vars, zeroimpllbs, zeroimplubs, zeroproplbs, zeropropubs, oneimpllbs, oneimplubs, oneproplbs, onepropubs,
         &localnfixedvars, &localnaggrvars, &localnimplications, &localnchgbds, cutoff) );

      *nfixedvars += localnfixedvars;
      *naggrvars  += localnaggrvars;
      *nchgbds    += localnchgbds;
      propdata->nfixings      += localnfixedvars;
      propdata->naggregations += localnaggrvars;
      propdata->nbdchgs       += localnchgbds;
      propdata->nimplications += localnimplications;

      if( localnfixedvars > 0 || localnaggrvars > 0 )
      {
         propdata->nuseless = 0;
         propdata->ntotaluseless = 0;
      }
      if( localnimplications > 0 || localnchgbds > 0 )
         propdata->ntotaluseless = 0;
   }

   *startidx = i;

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &onepropubs);
   SCIPfreeBufferArray(scip, &oneproplbs);
   SCIPfreeBufferArray(scip, &oneimplubs);
   SCIPfreeBufferArray(scip, &oneimpllbs);
   SCIPfreeBufferArray(scip, &zeropropubs);
   SCIPfreeBufferArray(scip, &zeroproplbs);
   SCIPfreeBufferArray(scip, &zeroimplubs);
   SCIPfreeBufferArray(scip, &zeroimpllbs);

   return SCIP_OKAY;
}



/*
 * Callback methods of propagator
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyProbing)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method for propagator */
   SCIP_CALL( SCIPincludePropProbing(scip) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);
   assert(propdata->sortedvars == NULL);
   assert(propdata->nsortedvars == 0);
   assert(propdata->nsortedbinvars == 0);

   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   initPropdata(propdata);

   return SCIP_OKAY;
}


/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( freeSortedvars(scip, propdata) );
   assert(propdata->sortedvars == NULL);
   assert(propdata->nsortedvars == 0);
   assert(propdata->nsortedbinvars == 0);


   return SCIP_OKAY;
}

/** presolving initialization method of propagator (called when presolving is about to begin) */
static
SCIP_DECL_PROPINITPRE(propInitpreProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->lastnode = -2;

   return SCIP_OKAY;
}


/** presolving deinitialization method of propagator (called after presolving has been finished) */
static
SCIP_DECL_PROPEXITPRE(propExitpreProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* delete the vars array, if the maximal number of runs are exceeded */
   if( propdata->maxruns >= 0 && SCIPgetNRuns(scip) >= propdata->maxruns )
   {
      SCIP_CALL( freeSortedvars(scip, propdata) );
      assert(propdata->sortedvars == NULL);
      assert(propdata->nsortedvars == 0);
      assert(propdata->nsortedbinvars == 0);
   }

   return SCIP_OKAY;
}


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolProbing)
{
  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* reset all propdata elements for stopping propagation earlier */
   propdata->nuseless = 0;
   propdata->ntotaluseless = 0;
   propdata->nsumuseless = 0;
#if 0
   propdata->nfixings = 0;
   propdata->naggregations = 0;
   propdata->nimplications = 0;
   propdata->nbdchgs = 0;
#endif

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
#define propExitsolProbing NULL


/** presolve method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int nvars;
   int nbinvars;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldnimplications;
   SCIP_Bool delay;
   SCIP_Bool cutoff;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* check, if probing should be applied in the current run */
   if( propdata->maxruns >= 0 && SCIPgetNRuns(scip) > propdata->maxruns )
      return SCIP_OKAY;

   /* if no domains changed since the last call, we don't need to probe */
   if( propdata->lastnode == -1 && nnewfixedvars == 0 && nnewaggrvars == 0 && nnewchgbds == 0 && nnewholes == 0 )
      return SCIP_OKAY;

   /**@todo currently we only perform probing on variables which are of binary type; there are, however, more
    *       implicit binary variables which can be found with the method SCIPvarIsBinary(); for these variable we could
    *       also perform probing */

   SCIPdebugMessage("executing probing (used %.1f sec)\n", SCIPpropGetTime(prop));

   *result = SCIP_DIDNOTFIND;

   /* allow some additional probings */
   propdata->nuseless -= propdata->nuseless/10;
   propdata->ntotaluseless -= propdata->ntotaluseless/10;

   /* get variable data */
   if( propdata->sortedvars == NULL )
   {
      SCIP_VAR** probvars;
      int i;

      assert(propdata->startidx == 0);

      SCIP_CALL( SCIPgetVarsData(scip, &probvars, &nvars, &nbinvars, NULL, NULL, NULL) );
      if( nbinvars == 0 )
         return SCIP_OKAY;

      SCIP_CALL( SCIPduplicateMemoryArray(scip, &propdata->sortedvars, probvars, nvars) );
      propdata->nsortedvars = nvars;
      propdata->nsortedbinvars = nbinvars;

      /* capture variables to make sure, the variables are not deleted */
      for( i = 0; i < propdata->nsortedvars; ++i )
      {
         SCIP_CALL( SCIPcaptureVar(scip, propdata->sortedvars[i]) );
      }
   }
   else
      nbinvars = SCIPgetNBinVars(scip);

   /* if we probed all binary variables in previous runs, start again with the first one */
   if( propdata->lastnode != -1 && propdata->startidx >= nbinvars )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "   (%.1fs) probing cycle finished: starting next cycle\n", SCIPgetSolvingTime(scip));
      propdata->startidx = 0;
      propdata->lastsortstartidx = -1;
      propdata->nuseless = 0;
      propdata->ntotaluseless = 0;
   }
   propdata->lastnode = -1;

   /* sort the binary variables by number of rounding locks, if at least 100 variables were probed since last sort */
   if( propdata->lastsortstartidx < 0 || propdata->startidx - propdata->lastsortstartidx >= 100 )
   {
      SCIP_CALL( sortVariables(scip, propdata->sortedvars, propdata->nsortedbinvars, propdata->startidx) );
      propdata->lastsortstartidx = propdata->startidx;
   }

   oldnfixedvars = *nfixedvars;
   oldnaggrvars = *naggrvars;
   oldnchgbds = *nchgbds;
   oldnimplications = propdata->nimplications;

   /* start probing on variables */
   SCIP_CALL( applyProbing(scip, propdata, propdata->sortedvars, propdata->nsortedvars, propdata->nsortedbinvars, &(propdata->startidx), nfixedvars, naggrvars, nchgbds, oldnfixedvars, oldnaggrvars, &delay, &cutoff) );

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( delay )
   {
      *result = SCIP_DELAYED;
      /* probing was interrupted because it reached the maximal fixings parameter, so we want to rerun it at the next call */
      propdata->lastnode = -2;
   }
   else if( *nfixedvars > oldnfixedvars || *naggrvars > oldnaggrvars || *nchgbds > oldnchgbds
      || propdata->nimplications > oldnimplications )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecProbing)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR** binvars;
   int nvars;
   int nbinvars;
   int i;
   int nfixedvars;
   int naggrvars;
   int nchgbds;
   int oldnfixedvars;
   int oldnaggrvars;
   int oldnchgbds;
   int oldnimplications;
   int startidx;
   SCIP_Bool delay;
   SCIP_Bool cutoff;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* avoid recursive infinity loop */
   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* only call heuristic, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* if already called stop */
   if( propdata->lastnode == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) )
      return SCIP_OKAY;

   propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   /* get number of fractional variables that should be integral */
   nvars = SCIPgetNLPBranchCands(scip);
   nbinvars = 0;

   /* alloc array for fractional binary variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nvars) );

   /* get fractional variables that should be integral */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &vars, NULL, NULL, &nvars, NULL) );

   /* copy binary variables to array binvars */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      var = vars[i];

      assert(var != NULL);
      if( SCIPvarIsBinary(var) )
      {
         assert(SCIPvarGetLbLocal(var) < 0.5);
         assert(SCIPvarGetUbLocal(var) > 0.5);

         binvars[nbinvars] = var;
         ++nbinvars;
      }
   }
   SCIPdebugMessage("problem <%s> node %"SCIP_LONGINT_FORMAT" probing propagation found %d of %d possible probing candidates\n", SCIPgetProbName(scip), SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), nbinvars, nvars);

   if( nbinvars == 0 )
   {
      *result = SCIP_DIDNOTFIND;
      goto TERMINATE;
   }

   /* sort binary variables */
   SCIP_CALL( sortVariables(scip, binvars, nbinvars, 0) );

   oldnfixedvars = 0;
   oldnaggrvars = 0;
   oldnchgbds = 0;
   nfixedvars = 0;
   naggrvars = 0;
   nchgbds = 0;
   startidx = 0;
   oldnimplications = propdata->nimplications;
   
   /* start probing on found variables */
   SCIP_CALL( applyProbing(scip, propdata, binvars, nbinvars, nbinvars, &startidx, &nfixedvars, &naggrvars, &nchgbds, oldnfixedvars, oldnaggrvars, &delay, &cutoff) );
   SCIPdebugMessage("probing propagation found %d fixings, %d aggregation, %d nchgbds, and %d implications\n", nfixedvars, naggrvars, nchgbds, (propdata->nimplications) - oldnimplications);

   if( delay )
   {
      /* probing was interrupted because it reached the maximal fixings parameter, so we want to rerun it at the next call */
      propdata->lastnode = -2;
   }

   /* adjust result code */
   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > oldnfixedvars || naggrvars > oldnaggrvars || nchgbds > oldnchgbds )
      *result = SCIP_REDUCEDDOM;

 TERMINATE:
   SCIPfreeBufferArray(scip, &binvars);

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropProbing)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}




/*
 * propagator specific interface methods
 */

/** creates the probing propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;

   /* create probing propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   initPropdata(propdata);

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOL_DELAY,
         propCopyProbing,
         propFreeProbing, propInitProbing, propExitProbing, propInitpreProbing, propExitpreProbing, propInitsolProbing, propExitsolProbing, propPresolProbing, propExecProbing, propRespropProbing, propdata) );

   /* add probing propagator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/probing/maxruns",
         "maximal number of runs, probing participates in (-1: no limit)",
         &propdata->maxruns, FALSE, DEFAULT_MAXRUNS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/probing/proprounds",
         "maximal number of propagation rounds in probing subproblems (-1: no limit, 0: auto)",
         &propdata->proprounds, TRUE, DEFAULT_PROPROUNDS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/probing/maxfixings",
         "maximal number of fixings found, until probing is interrupted (0: don't iterrupt)",
         &propdata->maxfixings, TRUE, DEFAULT_MAXFIXINGS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/probing/maxuseless",
         "maximal number of successive probings without fixings, until probing is aborted (0: don't abort)",
         &propdata->maxuseless, TRUE, DEFAULT_MAXUSELESS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/probing/maxtotaluseless",
         "maximal number of successive probings without fixings, bound changes, and implications, until probing is aborted (0: don't abort)",
         &propdata->maxtotaluseless, TRUE, DEFAULT_MAXTOTALUSELESS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/probing/maxsumuseless",
         "maximal number of probings without fixings, until probing is aborted (0: don't abort)",
         &propdata->maxsumuseless, TRUE, DEFAULT_MAXSUMUSELESS, 0, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}


/** analyses boundchanges resulting from probing on a variable and performs deduced fixations, aggregations, and domain tightenings
 * Given a variable probingvar with domain [l,u] and bound tightening results from reducing the domain
 * once to [l,leftub] and once to [rightlb,u], the method computes and applies resulting variable fixations, aggregations,
 * implications, and bound changes. Variable probingvar does not need to be binary.
 * The whole domain of probingvar need to be covered by the left and right branches, i.e.,
 * we assume leftub >= rightlb for continuous variables or floor(leftub) >= ceil(rightlb)-1 for discrete variables.
 * Bounds after applying implications and cliques do not need to be provided, but if they are omitted and probingvar is a binary variable,
 * then already existing implications may be added.
 */
SCIP_RETCODE SCIPanalyzeDeductionsProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             probingvar,         /**< the probing variable */
   SCIP_Real             leftub,             /**< upper bound of probing variable in left branch */
   SCIP_Real             rightlb,            /**< lower bound of probing variable in right branch */
   int                   nvars,              /**< number of variables which bound changes should be analyzed */
   SCIP_VAR**            vars,               /**< variables which bound changes should be analyzed */
   SCIP_Real*            leftimpllbs,        /**< lower bounds after applying implications and cliques in left branch, or NULL */
   SCIP_Real*            leftimplubs,        /**< upper bounds after applying implications and cliques in left branch, or NULL */
   SCIP_Real*            leftproplbs,        /**< lower bounds after applying domain propagation in left branch */
   SCIP_Real*            leftpropubs,        /**< upper bounds after applying domain propagation in left branch */
   SCIP_Real*            rightimpllbs,       /**< lower bounds after applying implications and cliques in right branch, or NULL */
   SCIP_Real*            rightimplubs,       /**< upper bounds after applying implications and cliques in right branch, or NULL */
   SCIP_Real*            rightproplbs,       /**< lower bounds after applying domain propagation in right branch */
   SCIP_Real*            rightpropubs,       /**< upper bounds after applying domain propagation in right branch */
   int*                  nfixedvars,         /**< pointer to counter which is increased by the number of deduced variable fixations */
   int*                  naggrvars,          /**< pointer to counter which is increased by the number of deduced variable aggregations */
   int*                  nimplications,      /**< pointer to counter which is increased by the number of deduced implications */
   int*                  nchgbds,            /**< pointer to counter which is increased by the number of deduced bound tightenings */
   SCIP_Bool*            cutoff              /**< buffer to store whether a cutoff is detected */
   )
{
   SCIP_Bool fixedleft;
   SCIP_Bool fixedright;
   int j;

   assert(scip != NULL);
   assert(probingvar != NULL);
   assert(SCIPisGE(scip, leftub,  SCIPvarGetLbLocal(probingvar))); /* left  branch should not be empty by default */
   assert(SCIPisLE(scip, rightlb, SCIPvarGetUbLocal(probingvar))); /* right branch should not be empty by default */
   assert(vars != NULL || nvars == 0);
   assert(leftproplbs != NULL);
   assert(leftpropubs != NULL);
   assert(rightproplbs != NULL);
   assert(rightpropubs != NULL);
   assert(nfixedvars != NULL);
   assert(naggrvars != NULL);
   assert(nimplications != NULL);
   assert(nchgbds != NULL);
   assert(cutoff != NULL);

   /* @todo the asserts below could be relaxed by taking domain holes into account */
   if( SCIPvarGetType(probingvar) != SCIP_VARTYPE_CONTINUOUS )
   {
      /* adjust bounds to actually used ones */
      leftub  = SCIPfloor(scip, leftub);
      rightlb = SCIPceil(scip, rightlb);
      /* assert dichotomy in case of discrete var: leftub >= rightlb - 1.0 */
      assert(SCIPisGE(scip, leftub, rightlb - 1.0));
   }
   else
   {
      /* assert dichotomy in case of continuous var: leftub >= rightlb */
      assert(SCIPisGE(scip, leftub, rightlb));
   }

   /* check if probing variable was fixed in the branches */
   fixedleft  = SCIPisEQ(scip, SCIPvarGetLbLocal(probingvar), leftub);
   fixedright = SCIPisEQ(scip, SCIPvarGetUbLocal(probingvar), rightlb);

   *cutoff = FALSE;

   for( j = 0; j < nvars && !*cutoff; ++j )
   {
      SCIP_Real newlb;
      SCIP_Real newub;

      assert(vars != NULL); /* for flexelint */
      assert(vars[j] != NULL);

      /* @todo: add holes, and even add holes if x was the probing variable and it followed a better bound on x itself */
      /* @todo: check if we probed on an integer variable, that this maybe led to aggregation on two other variables, i.e
       *        probing on x <= 1 and x >= 2 led to y = 1, z = 1 and y = 0, z = 0 resp., which means y = Z
       */

      /* if probing variable is binary, then there is nothing we could deduce here (variable should be fixed in both branches)
       * if it is not binary, we want to see if we found bound tightenings, even though it seems quite unlikely */
      if( vars[j] == probingvar && SCIPvarIsBinary(probingvar) )
         continue;

      /* new bounds of the variable is the union of the propagated bounds of the left and right case */
      newlb = MIN(leftproplbs[j], rightproplbs[j]);
      newub = MAX(leftpropubs[j], rightpropubs[j]);

      /* check for fixed variables */
      if( SCIPisEQ(scip, newlb, newub) )
      {
         SCIP_Real fixval;
         SCIP_Bool fixed;

         /* in both probings, variable j is deduced to the same value: fix variable to this value */
         fixval = SCIPselectSimpleValue(newlb - SCIPepsilon(scip), newub + SCIPepsilon(scip), MAXDNOM);
         if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 )
         {
            SCIP_CALL( SCIPfixVar(scip, vars[j], fixval, cutoff, &fixed) );
         }
         else
         {
            SCIP_CALL( SCIPtightenVarLb(scip, vars[j], fixval, TRUE, cutoff, &fixed) );
            if( !*cutoff )
            {
               SCIP_Bool tightened;

               SCIP_CALL( SCIPtightenVarUb(scip, vars[j], fixval, TRUE, cutoff, &tightened) );
               fixed &= tightened;
            }
         }
         if( fixed )
         {
            SCIPdebugMessage("fixed variable <%s> to %g due to probing on <%s> with nlocks=(%d/%d)\n",
               SCIPvarGetName(vars[j]), fixval,
               SCIPvarGetName(probingvar), SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));
            (*nfixedvars)++;
         }
         continue;
      }
      else
      {
         /* check for bound tightenings */
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Bool tightened;

         oldlb = SCIPvarGetLbLocal(vars[j]);
         oldub = SCIPvarGetUbLocal(vars[j]);
         if( SCIPisLbBetter(scip, newlb, oldlb, oldub) )
         {
            /* in both probings, variable j is deduced to be at least newlb: tighten lower bound */
            SCIP_CALL( SCIPtightenVarLb(scip, vars[j], newlb, FALSE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage("tightened lower bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[j]), oldlb, oldub, newlb,
                  SCIPvarGetName(probingvar), SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));
               (*nchgbds)++;
            }
         }
         if( SCIPisUbBetter(scip, newub, oldlb, oldub) && !*cutoff )
         {
            /* in both probings, variable j is deduced to be at most newub: tighten upper bound */
            SCIP_CALL( SCIPtightenVarUb(scip, vars[j], newub, FALSE, cutoff, &tightened) );
            if( tightened )
            {
               SCIPdebugMessage("tightened upper bound of variable <%s>[%g,%g] to %g due to probing on <%s> with nlocks=(%d/%d)\n",
                  SCIPvarGetName(vars[j]), oldlb, oldub, newub,
                  SCIPvarGetName(probingvar), SCIPvarGetNLocksDown(probingvar), SCIPvarGetNLocksUp(probingvar));
               (*nchgbds)++;
            }
         }
         if( *cutoff )
            break;
      }

      /* below we add aggregations and implications between probingvar and vars[j],
       * we don't want this if both variables are the same
       */
      if( vars[j] == probingvar )
         continue;

      /* check for aggregations and implications */
      if( fixedleft && fixedright &&
          SCIPisEQ(scip, leftproplbs[j],  leftpropubs[j]) && SCIPisEQ(scip, rightproplbs[j], rightpropubs[j]) )
      {
         /* vars[j] is fixed whenever probingvar is fixed, i.e.,
          *   vars[j] = leftproplbs[j] + (rightproplbs[j] - leftproplbs[j]) / (rightlb - leftub) * (probingvar - leftub)
          * -> both variables can be aggregated:
          *    (rightlb - leftub) * (vars[j] - leftproplbs[j]) = (rightproplbs[j] - leftproplbs[j]) * (probingvar - leftub)
          * -> (rightlb - leftub) * vars[j] - (rightproplbs[j] - leftproplbs[j]) * probingvar = leftproplbs[j] * rightlb - rightproplbs[j] * leftub
          *
          * check for case where both variables are binary: leftub = 1, rightlb = 0
          * case leftproplbs[j] = 0, rightproplbs[j] = 1, i.e., vars[j] and probingvar are fixed to same value
          *    -> aggregation is 1 * vars[j] - 1 * probingvar = 0 * 1 - 1 * 0 = 0 -> correct
          * case leftproplbs[j] = 1, rightproblbs[j] = 0, i.e., vars[j] and probingvar are fixed to opposite values
          *    -> aggregation is 1 * vars[j] + 1 * probingvar = 1 * 1 - 0 * 0 = 0 -> correct
          */
         if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
         {
            SCIP_Bool aggregated;
            SCIP_Bool redundant;

            SCIP_CALL( SCIPaggregateVars(scip, vars[j], probingvar,
               rightlb - leftub, -(rightproplbs[j] - leftproplbs[j]), leftproplbs[j] * rightlb - rightproplbs[j] * leftub,
               cutoff, &redundant, &aggregated) );

            if( aggregated )
            {
               SCIPdebugMessage("aggregated variables %g<%s> - %g<%s> == %g, nlocks=(%d/%d)\n",
                  rightlb - leftub, SCIPvarGetName(vars[j]),
                  rightproplbs[j] - leftproplbs[j], SCIPvarGetName(probingvar),
                  leftproplbs[j] * rightlb - rightproplbs[j] * leftub,
                  SCIPvarGetNLocksDown(vars[j]), SCIPvarGetNLocksUp(probingvar));
               (*naggrvars)++;
            }
         }
         else if( SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0 && (SCIPvarGetType(probingvar) == SCIP_VARTYPE_BINARY || SCIPvarGetType(probingvar) == SCIP_VARTYPE_INTEGER) )
         {
            /* if we are not in presolving, then we cannot do aggregations
             * but we can use variable bounds to code the same equality
             * vars[j] == ((leftproplbs[j] * rightlb - rightproplbs[j] * leftub) + (rightproplbs[j] - leftproplbs[j]) * probingvar) / (rightlb - leftub)
             */
            int nboundchanges;

            assert(!SCIPisEQ(scip, leftub, rightlb));

            SCIP_CALL( SCIPaddVarVlb(scip, vars[j], probingvar, (rightproplbs[j] - leftproplbs[j]) / (rightlb - leftub), (leftproplbs[j] * rightlb - rightproplbs[j] * leftub) / (rightlb - leftub), cutoff, &nboundchanges) );
            (*nchgbds) += nboundchanges;
            if( !*cutoff )
            {
               SCIP_CALL( SCIPaddVarVub(scip, vars[j], probingvar, (rightproplbs[j] - leftproplbs[j]) / (rightlb - leftub), (leftproplbs[j] * rightlb - rightproplbs[j] * leftub) / (rightlb - leftub), cutoff, &nboundchanges) );
               (*nchgbds) += nboundchanges;
            }
            (*nimplications)++;
         }
         /* if probingvar is continuous and we are in solving stage, then we do nothing, but it's unlikely that we get here (fixedleft && fixedright) with a continuous variable */
      }
      /* @todo: check if we can add variable lowerbounds/upperbounds on integer variables */
      /* can only add implications on binary variables which are globally valid */
      else if( SCIPvarGetType(probingvar) == SCIP_VARTYPE_BINARY &&  (SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0) )
      {
         /* implications can be added only for binary variables */
         int nboundchanges;

         /* since probing var is binary variable, probing should have fixed variable in both branches,
          * which is to 0.0 in the left branch and to 1.0 in the right branch */
         assert(fixedleft);
         assert(fixedright);
         assert(SCIPisZero(scip, leftub));
         assert(SCIPisEQ(scip, rightlb, 1.0));

         if( SCIPisEQ(scip, newlb, leftpropubs[j]) && (leftimplubs == NULL || leftimplubs[j] > leftpropubs[j]) )
         {
            /* vars[j] is fixed to lower bound whenever probingvar is fixed to 0.0
             * and implication is not already known
             * -> insert implication: probingvar == 0  =>  vars[j] <= leftpropubs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), leftpropubs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, leftpropubs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         else if( SCIPisEQ(scip, newub, leftproplbs[j]) && (leftimpllbs == NULL || leftimpllbs[j] < leftproplbs[j]) )
         {
            /* vars[j] is fixed to upper bound whenever probingvar is fixed to 0.0
             * and implication is not already known
             * -> insert implication: probingvar == 0  =>  vars[j] >= leftproplbs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), leftproplbs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, leftproplbs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         /* we can do an else here, since the case where vars[j] is fixed for both fixings of probingvar had been handled as aggregation */
         else if( SCIPisEQ(scip, newlb, rightpropubs[j]) && (rightimplubs == NULL || rightimplubs[j] > rightpropubs[j]) )
         {
            /* vars[j] is fixed to lower bound whenever probingvar is fixed to 1.0
             * and implication is not already known
             * -> insert implication: probingvar == 1  =>  vars[j] <= rightpropubs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), rightpropubs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, rightpropubs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         else if( SCIPisEQ(scip, newub, rightproplbs[j]) && (rightimpllbs == NULL || rightimpllbs[j] < rightproplbs[j]) )
         {
            /* vars[j] is fixed to upper bound whenever probingvar is fixed to 1.0
             * and implication is not already known
             * -> insert implication: probingvar == 1  =>  vars[j] >= leftproplbs[j]
             */
            /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s> == %g\n",
              SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), rightproplbs[j]);*/
            SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, rightproplbs[j],
               cutoff, &nboundchanges) );
            (*nimplications)++;
            (*nchgbds) += nboundchanges;
         }
         else if( SCIPvarGetType(vars[j]) != SCIP_VARTYPE_BINARY )
         {
            /* check for implications for lower or upper bounds (only store implications with bounds tightened at least by 0.5)
             * in case of binary variables, this should have been handled in the previous cases, since every boundchange also fixes the variable
             */
            if( leftpropubs[j] < newub - 0.5 && (leftimplubs == NULL || leftpropubs[j] < leftimplubs[j]) )
            {
               /* insert implication: probingvar == 0  =>  vars[j] <= leftpropubs[j] */
               /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] <= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, leftpropubs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_UPPER, leftpropubs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
            if( leftproplbs[j] > newlb + 0.5 && (leftimpllbs == NULL || leftproplbs[j] > leftimpllbs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 0  =>  vars[j] >= leftproplbs[j] */
               /*SCIPdebugMessage("found implication <%s> == 0  =>  <%s>[%g,%g] >= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, leftproplbs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, FALSE, vars[j], SCIP_BOUNDTYPE_LOWER, leftproplbs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
            if( rightpropubs[j] < newub - 0.5 && (rightimplubs == NULL || rightpropubs[j] < rightimplubs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 1  =>  vars[j] <= rightpropubs[j] */
               /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] <= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, rightpropubs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_UPPER, rightpropubs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
            if( rightproplbs[j] > newlb + 0.5 && (rightimpllbs == NULL || rightproplbs[j] > rightimpllbs[j]) && !*cutoff )
            {
               /* insert implication: probingvar == 1  =>  vars[j] >= rightproplbs[j] */
               /*SCIPdebugMessage("found implication <%s> == 1  =>  <%s>[%g,%g] >= %g\n",
                 SCIPvarGetName(probingvar), SCIPvarGetName(vars[j]), newlb, newub, rightproplbs[j]);*/
               SCIP_CALL( SCIPaddVarImplication(scip, probingvar, TRUE, vars[j], SCIP_BOUNDTYPE_LOWER, rightproplbs[j],
                     cutoff, &nboundchanges) );
               (*nimplications)++;
               (*nchgbds) += nboundchanges;
            }
         }
      }
   }

   return SCIP_OKAY;
}
