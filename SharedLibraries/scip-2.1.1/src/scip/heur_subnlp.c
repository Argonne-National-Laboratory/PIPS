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

/**@file    heur_subnlp.c
 * @brief   NLP local search primal heuristic using sub-SCIPs
 * @author  Stefan Vigerske
 * 
 * @todo set cutoff or similar in NLP
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_subnlp.h"
#include "nlpi/nlpi.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_bounddisjunction.h"

#define HEUR_NAME        "subnlp"
#define HEUR_DESC        "primal heuristic that performs a local search in an NLP after fixing integer variables and presolving"
#define HEUR_DISPCHAR    'q'
#define HEUR_PRIORITY    -2000000
#define HEUR_FREQ        1
#define HEUR_FREQOFS     0
#define HEUR_MAXDEPTH    -1
#define HEUR_TIMING      SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP FALSE               /**< does the heuristic use a secondary SCIP instance? we set this to FALSE because we want this heuristic to also run within other heuristics */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP*                 subscip;            /**< copy of CIP where presolving and NLP solving is done */
   SCIP_Bool             triedsetupsubscip;  /**< whether we have tried to setup a sub-SCIP */
   SCIP_Bool             subscipisvalid;     /**< whether all constraints have been copied */
   int                   nseriousnlpierror;  /**< number of consecutive serious NLP solver failures (memout, ...) */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */

   int                   nvars;              /**< number of active transformed variables in SCIP */
   int                   nsubvars;           /**< number of original variables in sub-SCIP */
   SCIP_VAR**            var_subscip2scip;   /**< mapping variables in sub-SCIP to SCIP variables */
   SCIP_VAR**            var_scip2subscip;   /**< mapping variables in SCIP to sub-SCIP variables */

   SCIP_SOL*             startcand;          /**< candidate for start point for heuristic */
   SCIP_Real             startcandviol;      /**< violation of start point candidate w.r.t. constraint that reported this candidate */

   SCIP_NLPSTATISTICS*   nlpstatistics;      /**< statistics from NLP solver */
   SCIP_Bool             comblinearconsadded;/**< whether the linear constraint adding method has been called for combinatorial constraints already */
   SCIP_Bool             contlinearconsadded;/**< whether the linear constraint adding method has been called for continuous constraints already */

   int                   nlpverblevel;       /**< verbosity level of NLP solver */
   int                   nlpiterlimit;       /**< iteration limit of NLP solver; 0 for off */
   SCIP_Real             nlptimelimit;       /**< time limit of NLP solver; 0 for off */
   SCIP_Real             resolvetolfactor;   /**< factor for feasibility tolerance when resolving NLP due to disagreement of feasibility */
   SCIP_Bool             resolvefromscratch; /**< whether a resolve of an NLP due to disagreement of feasibility should be from the original starting point or the infeasible solution */
   char*                 nlpoptfile;         /**< name of NLP solver specific option file */
   SCIP_Real             minimprove;         /**< desired minimal improvement in objective function value when running heuristic */
   int                   maxpresolverounds;  /**< limit on number of presolve rounds in sub-SCIP */
   SCIP_Bool             forbidfixings;      /**< whether to add constraints that forbid specific fixations that turned out to be infeasible */
   SCIP_Bool             keepcopy;           /**< whether to keep SCIP copy or to create new copy each time heuristic is applied */

   SCIP_Longint          iterused;           /**< number of iterations used so far */
   int                   iteroffset;         /**< number of iterations added to the contingent of the total number of iterations */
   SCIP_Real             iterquot;           /**< contingent of NLP iterations in relation to the number of nodes in SCIP */
   int                   itermin;            /**< minimal number of iterations required to start local search */
   SCIP_Bool             runalways;          /**< whether to run NLP heuristic always (independent of iteroffset,iterquot,itermin) */
};


/*
 * Local methods
 */

/** creates copy of CIP from problem in SCIP */
static
SCIP_RETCODE createSubSCIP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   int nvars;
   SCIP_VAR** vars;
   SCIP_VAR** subvars;
   SCIP_VAR*  var;
   SCIP_VAR*  subvar;
   SCIP_Bool success;
   char probname[SCIP_MAXSTRLEN];
   int i;
   SCIP_HASHMAP* varsmap;
   SCIP_HASHMAP* conssmap;
   SCIP_HASHMAPLIST* list;
#ifdef SCIP_DEBUG
   static const SCIP_Bool copydisplays = FALSE;
   static const SCIP_Bool copyreader = FALSE;
#else
   static const SCIP_Bool copydisplays = TRUE;
   static const SCIP_Bool copyreader = TRUE;
#endif

   assert(heurdata != NULL);
   assert(heurdata->subscip == NULL);

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   heurdata->triedsetupsubscip = TRUE;

   /* initializing the subproblem */
   SCIP_CALL( SCIPcreate(&heurdata->subscip) );

   /* create variable hash mapping scip -> subscip */
   SCIP_CALL( SCIPhashmapCreate(&varsmap, SCIPblkmem(scip), MAX(nvars, 5)) );

   /* create sub-SCIP copy of CIP */

   /* copy interesting plugins */
   success = TRUE;
   SCIP_CALL( SCIPcopyPlugins(scip, heurdata->subscip,
         copyreader, /* readers */
         FALSE, /* pricers */
         TRUE,  /* conshdlrs */
         FALSE, /* conflicthdlrs */
         TRUE,  /* presolvers */
         FALSE, /* relaxators */
         FALSE, /* separators */
         TRUE,  /* propagators */
         FALSE, /* heuristics */
         TRUE,  /* eventhandler */
         TRUE,  /* nodeselectors (SCIP gives an error if there is none) */
         FALSE,  /* branchrules */
         copydisplays, /* displays */
         FALSE, /* dialogs */
         TRUE,  /* nlpis */
         &success) );
   if( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "In heur_subnlp: failed to copy some plugins to sub-SCIP, continue anyway\n");
   }

   /* check if we still have NLPI's in subscip */
   if( SCIPgetNNlpis(heurdata->subscip) <= 0 )
   {
      SCIPdebugMessage("some NLPIs from main SCIP did not copy into sub-SCIP, give up heuristic.\n");
      SCIP_CALL( SCIPfree(&heurdata->subscip) );
      SCIPhashmapFree(&varsmap);

      return SCIP_OKAY;
   }

   /* copy parameter settings */
   SCIP_CALL( SCIPcopyParamSettings(scip, heurdata->subscip) );

   /* create problem in sub-SCIP */
   /* get name of the original problem and add "subnlp" */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s_subnlp", SCIPgetProbName(scip));
   SCIP_CALL( SCIPcreateProb(heurdata->subscip, probname, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   /* copy all variables */
   SCIP_CALL( SCIPcopyVars(scip, heurdata->subscip, varsmap, NULL, TRUE) );

   /* copy as many constraints as possible */
   SCIP_CALL( SCIPhashmapCreate(&conssmap, SCIPblkmem(scip), SCIPcalcHashtableSize(2 * SCIPgetNConss(scip))) );
   SCIP_CALL( SCIPcopyConss(scip, heurdata->subscip, varsmap, conssmap, TRUE, FALSE, &heurdata->subscipisvalid) );
   SCIPhashmapFree(&conssmap);
   if( !heurdata->subscipisvalid )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "In heur_subnlp: failed to copy some constraints to sub-SCIP, continue anyway\n");
      SCIPdebugMessage("In heur_subnlp: failed to copy some constraints to sub-SCIP, continue anyway\n");
   }

   /* create arrays translating scip transformed vars to subscip original vars, and vice versa
    * capture variables in SCIP and sub-SCIP
    * catch global bound change events
    */

   SCIP_CALL( SCIPgetVarsData(heurdata->subscip, &subvars, &heurdata->nsubvars, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->var_subscip2scip, heurdata->nsubvars) );
#ifndef NDEBUG
   BMSclearMemoryArray(heurdata->var_subscip2scip, heurdata->nsubvars);
#endif

   heurdata->nvars = nvars;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &heurdata->var_scip2subscip, heurdata->nvars) );
#ifndef NDEBUG
   BMSclearMemoryArray(heurdata->var_scip2subscip, heurdata->nvars);
#endif

   /* we need to get all subscip variables, also those which are copies of fixed variables from the main scip
    * therefore we iterate over the hashmap
    */
   for( i = 0; i < SCIPhashmapGetNLists(varsmap); ++i )
   {
      for( list = SCIPhashmapGetList(varsmap, i); list != NULL; list = SCIPhashmapListGetNext(list) )
      {
         var    = (SCIP_VAR*)SCIPhashmapListGetOrigin(list);
         subvar = (SCIP_VAR*)SCIPhashmapListGetImage(list);

         assert(SCIPvarGetProbindex(subvar) >= 0);
         assert(SCIPvarGetProbindex(subvar) <= heurdata->nsubvars);

         if( SCIPvarIsActive(var) )
         {
            assert(SCIPvarGetProbindex(var) <= heurdata->nvars);
            assert(heurdata->var_scip2subscip[SCIPvarGetProbindex(var)] == NULL);  /* assert that we have no mapping for this var yet */
            heurdata->var_scip2subscip[SCIPvarGetProbindex(var)] = subvar;
         }

         assert(heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)] == NULL);  /* assert that we have no mapping for this subvar yet */
         heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)] = var;

         SCIP_CALL( SCIPcaptureVar(scip, var) );
         SCIP_CALL( SCIPcaptureVar(heurdata->subscip, subvar) );

         assert(SCIPisFeasEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetLbGlobal(subvar)));
         assert(SCIPisFeasEQ(scip, SCIPvarGetUbGlobal(var), SCIPvarGetUbGlobal(subvar)));

         SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_GBDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, NULL) );
      }
   }

#ifndef NDEBUG
   for( i = 0; i < heurdata->nvars; ++i )
   {
      assert(heurdata->var_scip2subscip[i] != NULL);
      assert((SCIP_VAR*)SCIPhashmapGetImage(varsmap, (void*)vars[i]) == heurdata->var_scip2subscip[i]);
   }
   for( i = 0; i < heurdata->nsubvars; ++i )
   {
      assert(heurdata->var_subscip2scip[i] != NULL);
      assert((SCIP_VAR*)SCIPhashmapGetImage(varsmap, (void*)heurdata->var_subscip2scip[i]) == subvars[i]);
   }
#endif

   /* do not need hashmap anymore */
   SCIPhashmapFree(&varsmap);

   /* initialize data structure for NLP solve statistics */
   SCIP_CALL( SCIPnlpStatisticsCreate(&heurdata->nlpstatistics) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(heurdata->subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "display/verblevel", 0) );

   /* disable conflict analysis and separation 
    * keep normal presolving, but disable probing and restarts
    * disable LP solve
    * set nodelimit to 0
    * heuristics and separators were not copied into subscip, so should not need to switch off
    */
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "presolving/maxrounds", heurdata->maxpresolverounds) );
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "propagating/probing/maxprerounds", 0) );
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "presolving/maxrestarts", 0) );

#ifdef SCIP_DEBUG
   /* for debugging, enable SCIP output */
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "display/verblevel", 5) );
#endif

   return SCIP_OKAY;
}

/** free sub-SCIP data structure */
static
SCIP_RETCODE freeSubSCIP(
   SCIP*                 scip,               /**< SCIP data structure */ 
   SCIP_HEURDATA*        heurdata            /**< heuristic data structure */
   )
{
   SCIP_VAR** subvars;
   int        nsubvars;
   int        i;
   SCIP_VAR*  var;
   SCIP_VAR*  subvar;

   assert(scip != NULL);
   assert(heurdata != NULL);

   assert(heurdata->subscip != NULL);

   /* free NLP statistics */
   if( heurdata->nlpstatistics != NULL )
      SCIPnlpStatisticsFree(&heurdata->nlpstatistics);
   assert(heurdata->nlpstatistics == NULL);

   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, NULL, NULL, NULL, NULL) );
   assert(nsubvars == heurdata->nsubvars);

   /* drop global bound change events 
    * release variables in SCIP and sub-SCIP
    */
   for( i = 0; i < heurdata->nsubvars; ++i )
   {
      subvar = subvars[i];
      assert(subvar != NULL);
      assert(SCIPvarGetProbindex(subvar) == i);

      var = heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) <= heurdata->nvars);
      assert(!SCIPvarIsActive(var) || heurdata->var_scip2subscip[SCIPvarGetProbindex(var)] == subvar);

      SCIP_CALL( SCIPdropVarEvent(scip, var, SCIP_EVENTTYPE_GBDCHANGED, heurdata->eventhdlr, (SCIP_EVENTDATA*)heurdata, -1) );

      SCIP_CALL( SCIPreleaseVar(heurdata->subscip, &subvar) );
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* free variable mappings subscip -> scip and scip -> subscip */
   SCIPfreeBlockMemoryArray(scip, &heurdata->var_subscip2scip, heurdata->nsubvars);
   SCIPfreeBlockMemoryArray(scip, &heurdata->var_scip2subscip, heurdata->nvars);
   heurdata->nsubvars = 0;
   heurdata->nvars = 0;

   /* free sub-SCIP */
   SCIP_CALL( SCIPfree(&heurdata->subscip) );

   return SCIP_OKAY;
}

/** process variable global bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            idx;

   assert(scip      != NULL);
   assert(event     != NULL);
   assert(eventdata != NULL);
   assert(eventhdlr != NULL);

   heurdata = (SCIP_HEURDATA*)eventdata;
   assert(heurdata  != NULL);

   var = SCIPeventGetVar(event);
   assert(var != NULL);

   idx = SCIPvarGetProbindex(var);
   /* if event corresponds to an active variable, we can easily look up the corresponding subvar
    * if it is an inactive variable that has been copied to the subproblem,
    * then we need to check the subscip2scip mapping
    * @todo we could do this faster if we keep the variables mapping from SCIPcopy around
    */
   if( idx >= 0 )
   {
      assert(idx < heurdata->nvars);

      subvar = heurdata->var_scip2subscip[idx];
   }
   else
   {
      for( idx = 0; idx < heurdata->nsubvars; ++idx )
      {
         if( heurdata->var_subscip2scip[idx] == var )
            break;
      }
      assert(idx < heurdata->nsubvars);
      subvar = SCIPgetVars(heurdata->subscip)[idx];
   }
   assert(subvar != NULL);

   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_GLBCHANGED )
   {
      SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, SCIPvarGetLbGlobal(var)) );
   }

   if( SCIPeventGetType(event) & SCIP_EVENTTYPE_GUBCHANGED )
   {
      SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, SCIPvarGetUbGlobal(var)) );
   }

   return SCIP_OKAY;
}

/** adds linear constraints from a SCIP instance to its NLP */
static
SCIP_RETCODE addLinearConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for linear constraints */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints to NLP */
   SCIP_Bool             addcontconss        /**< whether to add continuous    linear constraints to NLP */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   int           j;
   SCIP_Bool     iscombinatorial;
   int           nvars;
   SCIP_VAR**    vars;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNConss(conshdlr);
   conss  = SCIPconshdlrGetConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      /* under some circumstances, this method may be called even though the problem has been shown to be infeasible in presolve already
       * this infeasibility may come from a linear constraint with lhs > rhs
       * the NLP does not allow such constraints, so we skip them here
       */
      if( !SCIPisRelLE(scip, SCIPgetLhsLinear(scip, conss[i]), SCIPgetRhsLinear(scip, conss[i])) )
         continue;

      nvars = SCIPgetNVarsLinear(scip, conss[i]);
      vars  = SCIPgetVarsLinear(scip, conss[i]);

      /* check if constraint should be added, only need this check if we do not wanna any constraint anyway */
      if( !addcombconss || !addcontconss )
      {
         iscombinatorial = TRUE;

         for( j = 0; j < nvars; ++j )
            if( SCIPvarGetType(vars[j]) >= SCIP_VARTYPE_CONTINUOUS )
            {
               iscombinatorial = FALSE;
               break;
            }

         /* skip constraint, if not of interest */
         if( (iscombinatorial && !addcombconss) || (!iscombinatorial && !addcontconss) )
            continue;
      }

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            SCIPgetNVarsLinear(scip, conss[i]), SCIPgetVarsLinear(scip, conss[i]), SCIPgetValsLinear(scip, conss[i]),
            0, NULL, 0, NULL, NULL,
            SCIPgetLhsLinear(scip, conss[i]), SCIPgetRhsLinear(scip, conss[i])) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }

   return SCIP_OKAY;
}

/** adds variable bound constraints from a SCIP instance to its NLP */
static
SCIP_RETCODE addVarboundConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for linear constraints */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints to NLP */
   SCIP_Bool             addcontconss        /**< whether to add continuous    linear constraints to NLP */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   SCIP_VAR*     vars[2];
   SCIP_Real     coefs[2];
   SCIP_Bool     iscombinatorial;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNConss(conshdlr);
   conss  = SCIPconshdlrGetConss(conshdlr);

   if( nconss == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      vars[0] = SCIPgetVarVarbound(scip, conss[i]);
      vars[1] = SCIPgetVbdvarVarbound(scip, conss[i]);

      iscombinatorial = SCIPvarGetType(vars[0]) < SCIP_VARTYPE_CONTINUOUS && SCIPvarGetType(vars[1]) < SCIP_VARTYPE_CONTINUOUS;

      /* skip constraint, if not of interest */
      if( (iscombinatorial && !addcombconss) || (!iscombinatorial && !addcontconss) )
         continue;

      coefs[0] = 1.0;
      coefs[1] = SCIPgetVbdcoefVarbound(scip, conss[i]);

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            2, vars, coefs,
            0, NULL, 0, NULL, NULL,
            SCIPgetLhsVarbound(scip, conss[i]), SCIPgetRhsVarbound(scip, conss[i])) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }

   return SCIP_OKAY;
}

/** adds logic-or constraints to NLP */
static
SCIP_RETCODE addLogicOrConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler for linear constraints */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   int           j;
   SCIP_Real*    coefs;
   int           coefssize;
   int           nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNConss(conshdlr);
   if( !nconss )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);

   coefs = NULL;
   coefssize = 0;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      nvars = SCIPgetNVarsLogicor(scip, conss[i]);

      if( coefssize < nvars )
      {
         if( coefs == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, nvars) );
         }
         for( j = coefssize; j < nvars; ++j )
            coefs[j] = 1.0;
         coefssize = nvars;
      }

      /* logic or constraints: 1 == sum_j x_j */

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            nvars, SCIPgetVarsLogicor(scip, conss[i]), coefs,
            0, NULL, 0, NULL, NULL,
            1.0, 1.0) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }

   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** adds setppc constraints to NLP */
static
SCIP_RETCODE addSetppcConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler for linear constraints */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   int           j;
   SCIP_Real*    coefs;
   int           coefssize;
   int           nvars;
   SCIP_Real     lhs;
   SCIP_Real     rhs;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNConss(conshdlr);
   if( !nconss )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);

   coefs = NULL;
   coefssize = 0;

   for( i = 0; i < nconss; ++i )
   {
      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      nvars = SCIPgetNVarsSetppc(scip, conss[i]);

      if( coefssize < nvars )
      {
         if( coefs == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, nvars) );
         }
         for( j = coefssize; j < nvars; ++j )
            coefs[j] = 1.0;
         coefssize = nvars;
      }

      /* setppc constraint: 1 ~ sum_j x_j */

      switch( SCIPgetTypeSetppc(scip, conss[i]) )
      {
      case SCIP_SETPPCTYPE_PARTITIONING:
         lhs = 1.0;
         rhs = 1.0;
         break;

      case SCIP_SETPPCTYPE_PACKING:
         lhs = -SCIPinfinity(scip);
         rhs = 1.0;
         break;

      case SCIP_SETPPCTYPE_COVERING:
         lhs = 1.0;
         rhs = SCIPinfinity(scip);
         break;

      default:
         SCIPerrorMessage("unexpected setppc type\n");
         return SCIP_ERROR;
      }

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            nvars, SCIPgetVarsSetppc(scip, conss[i]), coefs,
            0, NULL, 0, NULL, NULL,
            lhs, rhs) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }

   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** adds knapsack constraints to NLP */
static
SCIP_RETCODE addKnapsackConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler for linear constraints */
   )
{
   SCIP_CONS**   conss;
   int           nconss;
   SCIP_NLROW*   nlrow;
   int           i;
   int           j;
   SCIP_Real*    coefs;
   int           coefssize;
   int           nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);

   nconss = SCIPconshdlrGetNConss(conshdlr);
   if( !nconss )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert(conss != NULL);

   coefs = NULL;
   coefssize = 0;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_Longint* weights;

      /* skip local and redundant constraints */
      if( !SCIPconsIsEnabled(conss[i]) || !SCIPconsIsChecked(conss[i]) )
         continue;

      nvars = SCIPgetNVarsKnapsack(scip, conss[i]);

      if( coefssize < nvars )
      {
         if( coefs == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nvars) );
         }
         else
         {
            SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, nvars) );
         }
         coefssize = nvars;
      }

      weights = SCIPgetWeightsKnapsack(scip, conss[i]);
      for( j = 0; j < nvars; ++j )
         coefs[j] = (SCIP_Real)weights[j];  /*lint !e613*/

      SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0,
            nvars, SCIPgetVarsKnapsack(scip, conss[i]), coefs,
            0, NULL, 0, NULL, NULL,
            -SCIPinfinity(scip), (SCIP_Real)SCIPgetCapacityKnapsack(scip, conss[i])) );

      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }

   SCIPfreeBufferArrayNull(scip, &coefs);

   return SCIP_OKAY;
}

/** adds combinatorial and/or continuous variants of linear constraints from a SCIP instance to its NLP */
static
SCIP_RETCODE addLinearConstraintsToNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints to NLP */
   SCIP_Bool             addcontconss        /**< whether to add continuous    linear constraints to NLP */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* add linear constraints */
   conshdlr = SCIPfindConshdlr(scip, "linear");
   if( conshdlr != NULL )
   {
      SCIP_CALL( addLinearConstraints(scip, conshdlr, addcombconss, addcontconss) );
   }

   /* add varbound constraints */
   conshdlr = SCIPfindConshdlr(scip, "varbound");
   if( conshdlr != NULL )
   {
      SCIP_CALL( addVarboundConstraints(scip, conshdlr, addcombconss, addcontconss) );
   }

   if( addcombconss )
   {
      /* add logic-or constraints */
      conshdlr = SCIPfindConshdlr(scip, "logicor");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addLogicOrConstraints(scip, conshdlr) );
      }

      /* add setppc constraints */
      conshdlr = SCIPfindConshdlr(scip, "setppc");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addSetppcConstraints(scip, conshdlr) );
      }

      /* add knapsack constraints */
      conshdlr = SCIPfindConshdlr(scip, "knapsack");
      if( conshdlr != NULL )
      {
         SCIP_CALL( addKnapsackConstraints(scip, conshdlr) );
      }
   }

   return SCIP_OKAY;
}

/* creates a SCIP_SOL in our SCIP space out of the solution from NLP solver in sub-SCIP */
static
SCIP_RETCODE createSolFromNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL**            sol                 /**< buffer to store solution value; if pointing to NULL, then a new solution is created, otherwise values in the given one are overwritten */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol  != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( *sol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );
   }

   /* sub-SCIP may have more variables than the number of active (transformed) variables in the main SCIP
    * since constraint copying may have required the copy of variables that are fixed in the main SCIP
    */ 
   assert(heurdata->nsubvars <= SCIPgetNOrigVars(heurdata->subscip));

   for( i = 0; i < heurdata->nsubvars; ++i )
   {
      var = heurdata->var_subscip2scip[i];
      if( var == NULL || !SCIPvarIsActive(var) )
         continue;

      subvar = SCIPgetOrigVars(heurdata->subscip)[i];
      assert(subvar != NULL);

      assert(SCIPvarGetNLPSol(subvar) != SCIP_INVALID);  /*lint !e777*/
      SCIP_CALL( SCIPsetSolVal(scip, *sol, var, SCIPvarGetNLPSol(subvar)) );
   }

   return SCIP_OKAY;
}

/* creates a SCIP_SOL in our SCIP space out of the SCIP_SOL from a sub-SCIP */
static
SCIP_RETCODE createSolFromSubScipSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL**            sol,                /**< buffer to store solution value; if pointing to NULL, then a new solution is created, otherwise values in the given one are overwritten */
   SCIP_SOL*             subsol              /**< solution of sub-SCIP */
   )
{
   SCIP_HEURDATA* heurdata;
   int            i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol  != NULL);
   assert(subsol != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( *sol == NULL )
   {
      SCIP_CALL( SCIPcreateSol(scip, sol, heur) );
   }

   assert(heurdata->nsubvars == SCIPgetNOrigVars(heurdata->subscip));
   for( i = 0; i < heurdata->nsubvars; ++i )
   {
      if( heurdata->var_subscip2scip[i] == NULL || !SCIPvarIsActive(heurdata->var_subscip2scip[i]) )
         continue;
      SCIP_CALL( SCIPsetSolVal(scip, *sol, heurdata->var_subscip2scip[i],
            SCIPgetSolVal(heurdata->subscip, subsol, SCIPgetOrigVars(heurdata->subscip)[i])) );
   }

   return SCIP_OKAY;
}

/* solves the subNLP specified in subscip */
static
SCIP_RETCODE solveSubNLP(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< buffer to store result, DIDNOTFIND, FOUNDSOL, or CUTOFF        */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver, or -1 for default of NLP heuristic */
   SCIP_Real             timelimit,          /**< time limit for NLP solver                                      */
   SCIP_Longint*         iterused,           /**< buffer to store number of iterations used by NLP solver, or NULL if not of interest */
   SCIP_Bool             tighttolerances,    /**< whether to use tight feasibility tolerances and reduce presolve */
   SCIP_SOL*             resultsol           /**< a solution where to store found solution values, if any, or NULL if to try adding to SCIP */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_Real*     startpoint;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            i;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(result != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( tighttolerances )
   {
      /* reduce feasibility tolerance of sub-SCIP and do less aggressive presolve */
      SCIP_CALL( SCIPsetRealParam(heurdata->subscip, "numerics/feastol", heurdata->resolvetolfactor*SCIPfeastol(scip)) );
      SCIP_CALL( SCIPsetRealParam(heurdata->subscip, "numerics/epsilon", heurdata->resolvetolfactor*SCIPepsilon(scip)) );
      SCIP_CALL( SCIPsetPresolving(heurdata->subscip, SCIP_PARAMSETTING_FAST, TRUE) );
      SCIP_CALL( SCIPsetBoolParam(heurdata->subscip, "constraints/linear/aggregatevariables", FALSE) );
   }

   /* transform sub-SCIP */
   SCIP_CALL( SCIPtransformProb(heurdata->subscip) );

   /* presolve sub-SCIP
    *  set node limit to 1 so that presolve can go
    *  reset maxpresolverounds, in case user changed
    */
   SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 1LL) );
   SCIP_CALL( SCIPsetIntParam(heurdata->subscip, "presolving/maxrounds", heurdata->maxpresolverounds) );
   SCIP_CALL( SCIPpresolve(heurdata->subscip) );
   if( SCIPpressedCtrlC(heurdata->subscip) )
   {
      SCIPdebugMessage("SCIP presolve interrupted by user\n");
      goto CLEANUP;
   }
   if( SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED )
   {
      /* presolve probably found the subproblem infeasible */
      SCIPdebugMessage("SCIP returned from presolve in stage solved with status %d\n", SCIPgetStatus(heurdata->subscip));
      /* if presolve found subproblem infeasible, report this to caller by setting *result to cutoff */
      if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_INFEASIBLE )
         *result = SCIP_CUTOFF;
      goto CLEANUP;
   }
   assert(SCIPgetStage(heurdata->subscip) == SCIP_STAGE_PRESOLVED);

   if( SCIPgetNVars(heurdata->subscip) > 0 )
   {
      /* do initial solve, i.e., "solve" root node with node limit 0 (should do scip.c::initSolve and then stop immediately in solve.c::SCIPsolveCIP) */
      SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 0LL) );
      SCIP_CALL( SCIPsolve(heurdata->subscip) );

      /* If no NLP was constructed, then there were no nonlinearities after presolve.
       * So we increase the nodelimit to 1 and hope that SCIP will find some solution to this probably linear subproblem.
       */
      if( !SCIPisNLPConstructed(heurdata->subscip) )
      {
         SCIP_CALL( SCIPsetLongintParam(heurdata->subscip, "limits/nodes", 1LL) );
         SCIP_CALL( SCIPsolve(heurdata->subscip) );
      }
   }
   else
   {
      /* If all variables were removed by presolve, but presolve did not end with status SOLVED,
       * then we run solve, still with nodelimit=1, and hope to find some (maybe trivial) solution.
       */
      SCIP_CALL( SCIPsolve(heurdata->subscip) );
   }

   /* if sub-SCIP found solutions already, then pass them to main scip */
   for( i = 0; i < SCIPgetNSols(heurdata->subscip); ++i )
   {
      SCIP_Bool stored;

      if( resultsol == NULL )
      {
         SCIP_SOL* sol;

         sol = NULL;
         SCIP_CALL( createSolFromSubScipSol(scip, heur, &sol, SCIPgetSols(heurdata->subscip)[i]) );

         SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, FALSE, TRUE, &stored) );
         if( stored )
         {
            SCIPdebugMessage("SCIP stored solution from sub-SCIP root node\n");
            *result = SCIP_FOUNDSOL;
            break;
         }
         else
         {
            SCIPdebugMessage("SCIP did not store sub-SCIP optimal solution\n");
         }
      }
      else
      {
         SCIP_CALL( createSolFromSubScipSol(scip, heur, &resultsol, SCIPgetSols(heurdata->subscip)[i]) );

         SCIP_CALL( SCIPcheckSol(scip, resultsol, FALSE, TRUE, FALSE, TRUE, &stored) );
         if( stored )
         {
            SCIPdebugMessage("SCIP solution from sub-SCIP root node is feasible\n");
            *result = SCIP_FOUNDSOL;
            break;
         }
         else
         {
            SCIPdebugMessage("SCIP solution form sub-SCIP root node is not feasible\n");
         }
      }
   }

   /* we should either have variables, or the problem was trivial, in which case it should have been solved */
   assert(SCIPgetNVars(heurdata->subscip) > 0 || SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED);

   /* if subscip is infeasible here, we signal this to the caller */
   if( SCIPgetStatus(heurdata->subscip) == SCIP_STATUS_INFEASIBLE )
   {
      assert(SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED);
      *result = SCIP_CUTOFF;
      goto CLEANUP;
   }

   /* if we stopped for some other reason, or there is no NLP, we also stop */
   if( SCIPgetStage(heurdata->subscip) == SCIP_STAGE_SOLVED || !SCIPisNLPConstructed(heurdata->subscip) )
      goto CLEANUP;

   /* in most cases, the status should be nodelimit
    * in some cases, if the sub-SCIP is very easy, it may report optimal, so we do not need invoke an NLP solver
    * if the presolve found the problem infeasible, then there is no use in solving an NLP 
    * if the user interrupted or a timelimit was reached, then we should also stop here
    * unbounded is very unlikely to happen, in most cases, it should have been concluded in the main scip already
    */
   switch( SCIPgetStatus(heurdata->subscip) )
   {
   case SCIP_STATUS_NODELIMIT:
      break; /* this is the status that is most likely happening */
   case SCIP_STATUS_OPTIMAL: 
   case SCIP_STATUS_INFEASIBLE: 
   case SCIP_STATUS_USERINTERRUPT:
   case SCIP_STATUS_TIMELIMIT:
   case SCIP_STATUS_MEMLIMIT:
   case SCIP_STATUS_UNBOUNDED:
   case SCIP_STATUS_INFORUNBD:
      goto CLEANUP;
   default:
      SCIPerrorMessage("unexpected status of sub-SCIP: <%d>\n", SCIPgetStatus(heurdata->subscip));
      return SCIP_ERROR;
   } /*lint !e788*/

   /* if NLP timelimit is set to 0.0, then caller just wanted to see if the instance is still feasible after presolve */
   if( timelimit == 0.0 )
      goto CLEANUP;

   /* add non-combinatorial linear constraints from subscip into subNLP (shall be replaced by catching row events in NLP) */
   SCIP_CALL( addLinearConstraintsToNlp(heurdata->subscip, FALSE, TRUE) );

   /* set starting values (=refpoint, if not NULL; otherwise LP solution (or pseudo solution)) */
   SCIP_CALL( SCIPallocBufferArray(scip, &startpoint, SCIPgetNNLPVars(heurdata->subscip)) );
   for( i = 0; i < SCIPgetNNLPVars(heurdata->subscip); ++i )
   {
      SCIP_Real scalar;
      SCIP_Real constant;

      subvar = SCIPgetNLPVars(heurdata->subscip)[i];

      /* gets corresponding original variable */
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&subvar, &scalar, &constant) );
      if( subvar == NULL )
      {
         startpoint[i] = constant;
         continue;
      }

      assert(SCIPvarGetProbindex(subvar) >= 0);
      assert(SCIPvarGetProbindex(subvar) <  heurdata->nsubvars);
      var = heurdata->var_subscip2scip[SCIPvarGetProbindex(subvar)];
      if( var == NULL || REALABS(SCIPgetSolVal(scip, refpoint, var)) > 1.0e+12 )
         startpoint[i] = MIN(MAX(0.0, SCIPvarGetLbGlobal(subvar)), SCIPvarGetUbGlobal(subvar));  /*lint !e666*/
      else
         /* scalar*subvar+constant corresponds to nlpvar[i], so nlpvar[i] gets value scalar*varval+constant */
         startpoint[i] = scalar * SCIPgetSolVal(scip, refpoint, var) + constant;
   }
   SCIP_CALL( SCIPsetNLPInitialGuess(heurdata->subscip, startpoint) );

   SCIPfreeBufferArray(scip, &startpoint);

   *result = SCIP_DIDNOTFIND;

   /* setup NLP parameters */

   if( tighttolerances )
   {
      /* set feasibility tolerance, if tighttolerances is set */
      SCIP_CALL( SCIPsetNLPRealPar(heurdata->subscip, SCIP_NLPPAR_FEASTOL, heurdata->resolvetolfactor*SCIPfeastol(scip)) );
   }

   /* set option file to use by NLP solver */
   if( heurdata->nlpoptfile != NULL && *heurdata->nlpoptfile != '\0' )
   {
      SCIP_CALL( SCIPsetNLPStringPar(heurdata->subscip, SCIP_NLPPAR_OPTFILE, heurdata->nlpoptfile) );
   }

   /* set iteration limit for NLP solver */
   if( itercontingent == -1 && heurdata->nlpiterlimit > 0 )
      itercontingent = heurdata->nlpiterlimit;
   if( itercontingent > 0 )
   {
      SCIP_CALL( SCIPsetNLPIntPar(heurdata->subscip, SCIP_NLPPAR_ITLIM, (int)MIN(INT_MAX, itercontingent)) );
   }

   /* set time limit for NLP solver */
   SCIP_CALL( SCIPsetNLPRealPar(heurdata->subscip, SCIP_NLPPAR_TILIM, timelimit) );

   /* set verbosity of NLP solver */
   SCIP_CALL( SCIPsetNLPIntPar(heurdata->subscip, SCIP_NLPPAR_VERBLEVEL, heurdata->nlpverblevel) );


   /* let the NLP solver do its magic */
   SCIPdebugMessage("start NLP solve with iteration limit %"SCIP_LONGINT_FORMAT" and timelimit %g\n", itercontingent, timelimit);
   SCIP_CALL( SCIPsolveNLP(heurdata->subscip) );

   SCIPdebugMessage("NLP solver returned with termination status %d and solution status %d, objective value is %g\n",
      SCIPgetNLPTermstat(heurdata->subscip), SCIPgetNLPSolstat(heurdata->subscip), SCIPgetNLPObjval(heurdata->subscip));

   if( SCIPgetNLPTermstat(heurdata->subscip) >= SCIP_NLPTERMSTAT_MEMERR )
   {
      /* oops, something did not go well at all */
      ++heurdata->nseriousnlpierror;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, 
         "NLP solver in subNLP heuristic for problem <%s> returned with bad termination status %d. This was the %d%s successive time.\n",
         SCIPgetProbName(scip), SCIPgetNLPTermstat(heurdata->subscip), heurdata->nseriousnlpierror,
         heurdata->nseriousnlpierror == 1 ? "st" : heurdata->nseriousnlpierror == 2 ? "nd" : heurdata->nseriousnlpierror == 3 ? "rd" : "th");
      if( heurdata->nseriousnlpierror >= 5 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_MINIMAL, NULL, "Will not run NLP heuristic again for this run.\n");
         SCIP_CALL( freeSubSCIP(scip, heurdata) );
      }
      goto CLEANUP;
   }
   heurdata->nseriousnlpierror = 0;

   SCIP_CALL( SCIPgetNLPStatistics(heurdata->subscip, heurdata->nlpstatistics) );

   if( iterused != NULL )
      *iterused += SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics);
   SCIPdebugMessage("NLP solver used %d iterations and %g seconds\n",
      SCIPnlpStatisticsGetNIterations(heurdata->nlpstatistics), SCIPnlpStatisticsGetTotalTime(heurdata->nlpstatistics));

   /* NLP solver claims it found a feasible (maybe even optimal) solution
    * if the objective value is better than our cutoff, then try to add it
    * if we do not plan to add the solution (resultsol != NULL), then also check it if objective value is not better than objlimit
    */
   if( SCIPgetNLPSolstat(heurdata->subscip) <= SCIP_NLPSOLSTAT_FEASIBLE && (resultsol != NULL || SCIPisLE(scip, SCIPgetNLPObjval(heurdata->subscip), SCIPgetObjlimit(heurdata->subscip))) )
   {
      if( resultsol == NULL )
      {
         SCIP_SOL*  sol;
         SCIP_Bool  stored;

         sol = NULL;
         SCIP_CALL( createSolFromNLP(scip, heur, &sol) );

         if( heurdata->resolvefromscratch )
         {
#ifdef SCIP_DEBUG
            /* print the infeasibilities to stdout */
            SCIP_CALL( SCIPtrySolFree(scip, &sol, TRUE, TRUE, FALSE, TRUE, &stored) );
#else
            SCIP_CALL( SCIPtrySolFree(scip, &sol, FALSE, TRUE, FALSE, TRUE, &stored) );
#endif
         }
         else
         {
#ifdef SCIP_DEBUG
            /* print the infeasibilities to stdout */
            SCIP_CALL( SCIPtrySol(scip, sol, TRUE, TRUE, FALSE, TRUE, &stored) );
#else
            SCIP_CALL( SCIPtrySol(scip, sol, FALSE, TRUE, FALSE, TRUE, &stored) );
#endif
         }

         if( stored )
         {  /* SCIP stored solution (yippi!), so we are done */
            SCIPdebugMessage("SCIP stored solution\n");
            *result = SCIP_FOUNDSOL;
         }
         else if( !tighttolerances && heurdata->resolvetolfactor < 1.0 )
         {
            /* if SCIP does not like solution, we try again with tighter tolerances recreate subproblem and resolve with tighter tolerances */
            SCIPdebugMessage("solution reported by NLP solver not feasible for SCIP, resolve with feasibility tolerance %g and epsilon %g\n", heurdata->resolvetolfactor*SCIPfeastol(scip), heurdata->resolvetolfactor*SCIPepsilon(scip));

            /* free transformed problem */
            SCIP_CALL( SCIPfreeTransform(heurdata->subscip) );

            SCIP_CALL( solveSubNLP(scip, heur, result, heurdata->resolvefromscratch ? refpoint : sol, itercontingent, timelimit, iterused, TRUE, resultsol) );
         }
         else
         {
            SCIPdebugMessage("solution reported by NLP solver not stored by SCIP\n");
         }

         if( sol != NULL )
         {
            SCIP_CALL( SCIPfreeSol(scip, &sol) );
         }
      }
      else
      {
         SCIP_Bool feasible;

         SCIP_CALL( createSolFromNLP(scip, heur, &resultsol) );

#ifdef SCIP_DEBUG
         /* print the infeasibilities to stdout */
         SCIP_CALL( SCIPcheckSol(scip, resultsol, TRUE, TRUE, FALSE, TRUE, &feasible) );
#else
         SCIP_CALL( SCIPcheckSol(scip, resultsol, FALSE, TRUE, FALSE, TRUE, &feasible) );
#endif
         if( feasible )
         {
            /* SCIP find solution feasible, so we are done */
            SCIPdebugMessage("solution reported by NLP solver feasible for SCIP\n");
            *result = SCIP_FOUNDSOL;
         }
         else if( !tighttolerances && heurdata->resolvetolfactor < 1.0 )
         {
            /* if SCIP does not like solution, we try again with tighter tolerances
             * recreate subproblem and resolve with tighter tolerances 
             */
            SCIPdebugMessage("solution reported by NLP solver not feasible for SCIP, resolve with feasibility tolerance %g and epsilon %g\n", heurdata->resolvetolfactor*SCIPfeastol(scip), heurdata->resolvetolfactor*SCIPepsilon(scip));

            /* free transformed problem */
            SCIP_CALL( SCIPfreeTransform(heurdata->subscip) );

            SCIP_CALL( solveSubNLP(scip, heur, result, heurdata->resolvefromscratch ? refpoint : resultsol, itercontingent, timelimit, iterused, TRUE, resultsol) );
         }
         else
         {
            SCIPdebugMessage("solution reported by NLP solver not feasible for SCIP\n");
         }
      }
   }

 CLEANUP:
   if( heurdata->subscip != NULL )
   {
      SCIP_CALL( SCIPfreeTransform(heurdata->subscip) );
      if( tighttolerances )
      {
         /* reset feasibility tolerance of sub-SCIP and reset to normal presolve */
         SCIP_CALL( SCIPsetRealParam(heurdata->subscip, "numerics/feastol", SCIPfeastol(scip)) );
         SCIP_CALL( SCIPsetRealParam(heurdata->subscip, "numerics/epsilon", SCIPepsilon(scip)) );
         SCIP_CALL( SCIPsetPresolving(heurdata->subscip, SCIP_PARAMSETTING_DEFAULT, TRUE) );
         SCIP_CALL( SCIPresetParam(heurdata->subscip, "constraints/linear/aggregatevariables") );
      }
   }

   return SCIP_OKAY;
}


/** adds a set covering or bound disjunction constraint to the original problem */
static
SCIP_RETCODE forbidFixation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEURDATA*        heurdata            /**< heuristic data */
   )
{
   SCIP_VAR** subvars;
   int nsubvars;
   int nsubbinvars;
   int nsubintvars;
   SCIP_VAR* var;
   SCIP_VAR* subvar;
   SCIP_CONS* cons;
   SCIP_VAR** consvars;
   int nconsvars;
   char name[SCIP_MAXSTRLEN];
   int i;
   SCIP_Real fixval;

   assert(scip != NULL);

   SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
   assert(nsubvars == heurdata->nsubvars);

   if( nsubbinvars == 0 && nsubintvars == 0 )
   {
      /* If we did not fix any discrete variables but found the "sub"CIP infeasible, then also the CIP is infeasible. */
      SCIPwarningMessage("heur_subnlp found subCIP infeasible after fixing no variables, something is strange here...\n");
      return SCIP_OKAY;
   }

   /* initialize */
   cons = NULL;
   consvars = NULL;

   /* create constraint name */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "subnlp_cutoff");

   /* if all discrete variables in the CIP are binary, then we create a set covering constraint
    *  sum_{x_i fixed at 0} x_i + sum_{x_i fixed at 1} ~x_i >= 1
    */
   if( nsubintvars == 0 )
   {
      /* allocate memory for constraint variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nsubbinvars) );

      /* get fixations of discrete variables 
       * to be sure, we take the values that were put into the subCIP before
       */
      nconsvars = 0;
      for( i = nsubbinvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL || SCIPisEQ(scip, SCIPvarGetLbGlobal(subvar), SCIPvarGetUbGlobal(subvar))); /* otherwise we should have exited in the variable fixation loop */
         if( var == NULL )
            continue;

         fixval = SCIPvarGetLbGlobal(subvar);
         assert(fixval == SCIPvarGetUbGlobal(subvar)); /* variable should be fixed in sub-SCIP */  /*lint !e777*/
         assert(fixval == 0.0 || fixval == 1.0); /* we have rounded values before fixing */

         if( fixval == 0.0 )
         {
            /* variable fixed at lower bound */
            consvars[nconsvars] = var;
         }
         else
         {
            SCIP_CALL( SCIPgetNegatedVar(scip, var, &consvars[nconsvars]) );
         }

         ++nconsvars;
      }

      /* create conflict constraint
       * In undercover, ConsLogicor is used, since then the inequality is not added to the LP.
       * However, I may want to use Setcover to avoid that the same fixing is computed by some LP based heuristic again.
       */
      SCIP_CALL( SCIPcreateConsSetcover(scip, &cons, name, nconsvars, consvars,
            FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
   }
   else
   {
      /* if there are also integer variable, then create a bound disjunction constraint
       * x_1 >= fixval_1 + 1 || x_1 <= fixval_1 - 1 || x_2 >= fixval_2 + 1 || x_2 <= fixval_2 - 1 || ...
       */
      SCIP_BOUNDTYPE* boundtypes;
      SCIP_Real* bounds;

      /* allocate memory for constraint variables, boundtypes, and bounds
       * (there should be at most two literals for each integer variable)
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars,   nsubbinvars + 2*nsubintvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, nsubbinvars + 2*nsubintvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &bounds,     nsubbinvars + 2*nsubintvars) );

      /* get fixations of discrete variables 
       * to be sure, we take the values that were put into the subCIP before
       */
      nconsvars = 0;
      for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL || SCIPisEQ(scip, SCIPvarGetLbGlobal(subvar), SCIPvarGetUbGlobal(subvar))); /* otherwise we should have exited in the variable fixation loop */

         if( var == NULL )
            continue;

         fixval = SCIPvarGetLbGlobal(subvar);
         assert(fixval == SCIPvarGetUbGlobal(subvar)); /* variable should be fixed in sub-SCIP */   /*lint !e777*/
         assert((int)fixval == fixval); /* we have rounded values before fixing */
         assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY || SCIPvarGetLbGlobal(var) == fixval || SCIPvarGetUbGlobal(var) == fixval); /* for binaries, the fixval should be either 0.0 or 1.0 */  /*lint !e777*/ 

         if( SCIPvarGetLbGlobal(var) < fixval )
         {
            assert(nconsvars < nsubbinvars + 2*nsubintvars);

            /* literal x_i <= fixval-1 */
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_UPPER;
            bounds[nconsvars]     = fixval - 1.0;
            consvars[nconsvars]   = var;
            ++nconsvars;
         }

         if( SCIPvarGetUbGlobal(var) > fixval )
         {
            assert(nconsvars < nsubbinvars + 2*nsubintvars);

            /* literal x_i >= fixval+1 */
            boundtypes[nconsvars] = SCIP_BOUNDTYPE_LOWER;
            bounds[nconsvars]     = fixval + 1.0;
            consvars[nconsvars]   = var;
            ++nconsvars;
         }
      }

      /* create conflict constraint */
      SCIP_CALL( SCIPcreateConsBounddisjunction(scip, &cons, name, nconsvars, consvars, boundtypes, bounds,
            FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

      SCIPfreeBufferArray(scip, &consvars);
      SCIPfreeBufferArray(scip, &boundtypes);
      SCIPfreeBufferArray(scip, &bounds);
   }

   /* add and release constraint if created successfully */
   if( cons != NULL )
   {
      SCIPdebugMessage("adding constraint to forbid fixation in main problem\n");
      /* SCIPdebug( SCIPprintCons(scip, cons, NULL) ); */
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* free memory */
   SCIPfreeBufferArrayNull(scip, &consvars);

   return SCIP_OKAY;
}


/** main procedure of the subNLP heuristic */
SCIP_RETCODE SCIPapplyHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_RESULT*          result,             /**< pointer to store result of: did not run, solution found, no solution found, or fixing is infeasible (cutoff) */
   SCIP_SOL*             refpoint,           /**< point to take fixation of discrete variables from, and startpoint for NLP solver; if NULL, then LP solution is used */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver, or -1 for default of NLP heuristic */
   SCIP_Real             timelimit,          /**< time limit for NLP solver                                      */
   SCIP_Real             minimprove,         /**< desired minimal relative improvement in objective function value */
   SCIP_Longint*         iterused            /**< buffer to store number of iterations used by NLP solver, or NULL if not of interest */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            i;
   SCIP_Real      cutoff;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* try to setup NLP if not tried before */
   if( heurdata->subscip == NULL && !heurdata->triedsetupsubscip )
   {
      SCIP_CALL( createSubSCIP(scip, heurdata) );
   }

   *result = SCIP_DIDNOTRUN;

   /* not initialized */
   if( heurdata->subscip == NULL )
      return SCIP_OKAY;

   assert(heurdata->nsubvars > 0);
   assert(heurdata->var_subscip2scip != NULL);
   assert(!SCIPisTransformed(heurdata->subscip));

   if( iterused != NULL )
      *iterused = 0;

   /* fix discrete variables in sub-SCIP */
   if( SCIPgetNBinVars(heurdata->subscip) || SCIPgetNIntVars(heurdata->subscip) )
   {
      SCIP_Real  fixval;
      SCIP_VAR** subvars;
      int        nsubvars;
      int        nsubbinvars;
      int        nsubintvars;

      SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
      assert(nsubvars == heurdata->nsubvars);

      /* fix discrete variables to values in startpoint */
      for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL);

         /* at this point, variables in subscip and in our scip should have same bounds */
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(subvar), SCIPvarGetLbGlobal(var)));
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(subvar), SCIPvarGetUbGlobal(var)));

         fixval = SCIPgetSolVal(scip, refpoint, var);

         /* only run heuristic on integer feasible points */
         if( !SCIPisFeasIntegral(scip, fixval) )
         {
            if( refpoint || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL )
            {
               SCIPdebugMessage("skip NLP heuristic because start candidate not integer feasible: var <%s> has value %g\n", SCIPvarGetName(var), fixval);
               goto CLEANUP;
            }
            /* otherwise we desperately wanna run the NLP heur, so we continue and round what we have */
         }
         /* if we do not really have a startpoint, then we should take care that we do not fix variables to very large values
          *  thus, we set to 0.0 here and project on bounds below
          */
         if( ABS(fixval) > 1E+10 && !refpoint && SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
            fixval = 0.0;

         /* round fractional variables to the nearest integer,
          *  use exact integral value, if the variable is only integral within numerical tolerances
          */
         fixval = SCIPfloor(scip, fixval+0.5);

         /* adjust value to the global bounds of the corresponding SCIP variable */
         fixval = MAX(fixval, SCIPvarGetLbGlobal(var));  /*lint !e666*/
         fixval = MIN(fixval, SCIPvarGetUbGlobal(var));  /*lint !e666*/

         /* SCIPdebugMessage("fix variable <%s> to %g\n", SCIPvarGetName(var), fixval); */
         SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, fixval) );
         SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, fixval) );
      }
   }

   /* if there is already a solution, add an objective cutoff in sub-SCIP */
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIP_Real upperbound;

      cutoff = SCIPinfinity(scip);
      assert( !SCIPisInfinity(scip, SCIPgetUpperbound(scip)) );

      upperbound = SCIPgetUpperbound(scip) - SCIPsumepsilon(scip);

      if( !SCIPisInfinity(scip, -SCIPgetLowerbound(scip)) )
      {
         cutoff = (1-minimprove)*SCIPgetUpperbound(scip) + minimprove*SCIPgetLowerbound(scip);
      }
      else
      {
         if( SCIPgetUpperbound(scip) >= 0 )
            cutoff = ( 1.0 - minimprove ) * SCIPgetUpperbound(scip);
         else
            cutoff = ( 1.0 + minimprove ) * SCIPgetUpperbound(scip);
      }
      cutoff = MIN(upperbound, cutoff);
      SCIP_CALL( SCIPsetObjlimit(heurdata->subscip, cutoff) );
      SCIPdebugMessage("set objective limit %g\n", cutoff);
   }
   else
      cutoff = SCIPinfinity(scip);

   /* solve the subNLP and try to add solution to SCIP */
   SCIP_CALL( solveSubNLP(scip, heur, result, refpoint, itercontingent, timelimit, iterused, FALSE, NULL) );

   if( heurdata->subscip == NULL )
   {
      /* something horrible must have happened that we decided to give up completely on this heuristic */
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }
   assert(!SCIPisTransformed(heurdata->subscip));

   if( *result == SCIP_CUTOFF )
   {
      if( heurdata->subscipisvalid && SCIPgetNActivePricers(scip) == 0 )
      {
         /* if the subNLP is valid and turned out to be globally infeasible (i.e., proven by SCIP), then we forbid this fixation in the main problem */
         if( SCIPisInfinity(scip, cutoff) && heurdata->forbidfixings )
         {
            SCIP_CALL( forbidFixation(scip, heurdata) );
         }
      }
      else
      {
         /* if the subNLP turned out to be globally infeasible but we are not sure that we have a valid copy, we change to DIDNOTFIND */
         *result = SCIP_DIDNOTFIND;
      }
   }

 CLEANUP:
   /* if the heuristic was applied before solving has started, then destroy subSCIP, since EXITSOL may not be called
    * also if keepcopy is disabled, then destroy subSCIP
    */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING || !heurdata->keepcopy )
   {
      SCIP_CALL( freeSubSCIP(scip, heurdata) );
      heurdata->triedsetupsubscip = FALSE;
   }
   else if( SCIPgetNBinVars(heurdata->subscip) || SCIPgetNIntVars(heurdata->subscip) )
   {
      /* undo fixing of discrete variables in sub-SCIP */
      SCIP_VAR** subvars;
      int        nsubvars;
      int        nsubbinvars;
      int        nsubintvars;

      SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
      assert(nsubvars == heurdata->nsubvars);

      /* set bounds of discrete variables to original values */
      for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL); 

         SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, SCIPvarGetLbGlobal(var)) );
         SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, SCIPvarGetUbGlobal(var)) );
      }
   }

   return SCIP_OKAY;
}

/** for a given solution, resolves the corresponding subNLP and updates solution values for continuous variables, if NLP solution is feasible in original problem */
SCIP_RETCODE SCIPresolveSolHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic data structure */
   SCIP_SOL*             sol,                /**< solution for which to solve NLP, and where to store resolved solution values */
   SCIP_Bool*            success,            /**< buffer where to store whether a feasible solution was found */
   SCIP_Longint          itercontingent,     /**< iteration limit for NLP solver, or -1 for default of NLP heuristic */
   SCIP_Real             timelimit           /**< time limit for NLP solver */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_VAR*      var;
   SCIP_VAR*      subvar;
   int            i;
   SCIP_Real      cutoff;
   SCIP_RESULT    result;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(sol  != NULL);
   assert(success != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* try to setup NLP if not tried before */
   if( heurdata->subscip == NULL && !heurdata->triedsetupsubscip )
   {
      SCIP_CALL( createSubSCIP(scip, heurdata) );
   }

   *success = FALSE;

   /* not initialized */
   if( heurdata->subscip == NULL )
      return SCIP_OKAY;

   assert(heurdata->nsubvars > 0);
   assert(heurdata->var_subscip2scip != NULL);
   assert(!SCIPisTransformed(heurdata->subscip));

   result = SCIP_DIDNOTRUN;

   /* fix discrete variables in subSCIP */
   if( SCIPgetNBinVars(heurdata->subscip) || SCIPgetNIntVars(heurdata->subscip) )
   {
      SCIP_Real  fixval;
      SCIP_VAR** subvars;
      int        nsubvars;
      int        nsubbinvars;
      int        nsubintvars;

      SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
      assert(nsubvars == heurdata->nsubvars);

      /* fix discrete variables to values in startpoint */
      for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL);

         /* at this point, variables in subscip and in our scip should have same bounds */
         assert(SCIPisEQ(scip, SCIPvarGetLbGlobal(subvar), SCIPvarGetLbGlobal(var)));
         assert(SCIPisEQ(scip, SCIPvarGetUbGlobal(subvar), SCIPvarGetUbGlobal(var)));

         fixval = SCIPgetSolVal(scip, sol, var);

         /* only run heuristic on integer feasible points */
         if( !SCIPisFeasIntegral(scip, fixval) )
         {
            SCIPdebugMessage("skip NLP heuristic because start candidate not integer feasible: var <%s> has value %g\n", SCIPvarGetName(var), fixval);
            goto CLEANUP;
            /* otherwise we desperately want to run the NLP heur, so we continue and round what we have */
         }

         /* round fractional variables to the nearest integer,
          *  use exact integral value, if the variable is only integral within numerical tolerances
          */
         fixval = SCIPround(scip, fixval);

         /* adjust value to the global bounds of the corresponding SCIP variable */
         fixval = MAX(fixval, SCIPvarGetLbGlobal(var));  /*lint !e666*/
         fixval = MIN(fixval, SCIPvarGetUbGlobal(var));  /*lint !e666*/

         /* SCIPdebugMessage("fix variable <%s> to %g\n", SCIPvarGetName(var), fixval); */
         SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, fixval) );
         SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, fixval) );
      }
   }

   /* if there is already a solution, add an objective cutoff in subSCIP */
   cutoff = SCIPgetSolOrigObj(scip, sol);
   if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE )
   {
      cutoff += 0.01 * REALABS(cutoff);
   }
   else
   {
      cutoff -= 0.01 * REALABS(cutoff);
   }
   cutoff = SCIPtransformObj(scip, cutoff);
   SCIPdebugMessage("set objective limit %g\n", cutoff);
   SCIP_CALL( SCIPsetObjlimit(heurdata->subscip, cutoff) );

   /* solve the subNLP and try to add solution to SCIP */
   SCIP_CALL( solveSubNLP(scip, heur, &result, sol, itercontingent, timelimit, NULL, FALSE, sol) );

   if( heurdata->subscip == NULL )
   {
      /* something horrible must have happened that we decided to give up completely on this heuristic */
      return SCIP_OKAY;
   }
   assert(!SCIPisTransformed(heurdata->subscip));

   if( result == SCIP_FOUNDSOL )
      *success = TRUE;

 CLEANUP:
   /* if the heuristic was applied before solving has started, then destroy subSCIP, since EXITSOL may not be called
    * also if keepcopy is not set, then destroy subSCIP
    */
   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING || !heurdata->keepcopy )
   {
      SCIP_CALL( freeSubSCIP(scip, heurdata) );
      heurdata->triedsetupsubscip = FALSE;
   }
   else if( SCIPgetNBinVars(heurdata->subscip) || SCIPgetNIntVars(heurdata->subscip) )
   {
      /* undo fixing of discrete variables in subSCIP */
      SCIP_VAR** subvars;
      int        nsubvars;
      int        nsubbinvars;
      int        nsubintvars;

      SCIP_CALL( SCIPgetOrigVarsData(heurdata->subscip, &subvars, &nsubvars, &nsubbinvars, &nsubintvars, NULL, NULL) );
      assert(nsubvars == heurdata->nsubvars);

      /* set bounds of discrete variables to original values */
      for( i = nsubbinvars + nsubintvars - 1; i >= 0; --i )
      {
         subvar = subvars[i];
         assert(SCIPvarGetProbindex(subvar) == i);

         var = heurdata->var_subscip2scip[i];
         assert(var != NULL);

         SCIP_CALL( SCIPchgVarLbGlobal(heurdata->subscip, subvar, SCIPvarGetLbGlobal(var)) );
         SCIP_CALL( SCIPchgVarUbGlobal(heurdata->subscip, subvar, SCIPvarGetUbGlobal(var)) );
      }
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySubNlp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSubNlp(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSubNlp)
{
   SCIP_HEURDATA* heurdata;
   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->subscip == NULL);
   assert(heurdata->var_subscip2scip == NULL);
   assert(heurdata->var_scip2subscip == NULL);
   assert(heurdata->startcand == NULL);

   SCIPfreeMemory(scip, &heurdata);

   return SCIP_OKAY;
}


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitSubNlp)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return SCIP_OKAY;
}
#else
#define heurInitSubNlp NULL
#endif

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitSubNlp)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return SCIP_OKAY;
}
#else
#define heurExitSubNlp NULL
#endif

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
static
SCIP_DECL_HEURINITSOL(heurInitsolSubNlp)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);

   /* do not setup sub-SCIP if heuristic is never called by SCIP */
   if( SCIPheurGetFreq(heur) < 0 )
      return SCIP_OKAY;

   /* do not setup sub-SCIP if no NLP solver is available */
   if( SCIPgetNNlpis(scip) <= 0 )
      return SCIP_OKAY;

   /* do not setup sub-SCIP if no continuous nonlinear variables are present */
   if( !SCIPhasContinuousNonlinearitiesPresent(scip) )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);
   assert(heurdata->subscip == NULL);

   if( heurdata->keepcopy )
   {
      /* create sub-SCIP for later use */
      SCIP_CALL( createSubSCIP(scip, heurdata) );

      /* creating sub-SCIP may fail if the NLP solver interfaces did not copy into subscip */
      if( heurdata->subscip == NULL )
         return SCIP_OKAY;
   }

   /* if the heuristic is called at the root node, we want to be called directly after the initial root LP solve */
   if( SCIPheurGetFreqofs(heur) == 0 )
      SCIPheurSetTimingmask(heur, SCIP_HEURTIMING_DURINGLPLOOP | HEUR_TIMING);

   return SCIP_OKAY;
}

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
static
SCIP_DECL_HEUREXITSOL(heurExitsolSubNlp)
{
   SCIP_HEURDATA* heurdata;
   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic's data */  
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   if( heurdata->subscip != NULL )
   {
      SCIP_CALL( freeSubSCIP(scip, heurdata) );
   }

   /* free start candidate */
   if( heurdata->startcand != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
   }

   SCIPheurSetTimingmask(heur, HEUR_TIMING);

   /* reset some flags and counters */
   heurdata->triedsetupsubscip = FALSE;
   heurdata->comblinearconsadded = FALSE;
   heurdata->contlinearconsadded = FALSE;
   heurdata->nseriousnlpierror = 0;
   heurdata->iterused = 0;

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSubNlp)
{  /*lint --e{666,715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_Longint   itercontingent;
   SCIP_Real      timelimit;
   SCIP_Longint   iterused;

   assert(scip != NULL);
   assert(heur != NULL);

   /* get heuristic's data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* obviously, we did not do anything yet */
   *result = SCIP_DIDNOTRUN;

   /* if keepcopy and subscip == NULL, then InitsolNlp decided that we do not need an NLP solver,
    *   probably because we do not have nonlinear continuous or implicit integer variables
    * if triedsetupsubscip and subscip == NULL, then we run the heuristic already, but gave up due to some serious error
    * in both cases, we do not want to run
    *
    * otherwise, we continue and let SCIPapplyHeurSubNlp try to create subscip
    */
   if( heurdata->subscip == NULL && (heurdata->keepcopy || heurdata->triedsetupsubscip) )
      return SCIP_OKAY;

   if( heurdata->startcand == NULL )
   {
      /* if no start candidate is given, we consider the LP solution of the current node */

      /* at least if we are not called the first time, we call the heuristic only if an optimal LP solution is available 
       * if we are called the first time and the LP is unbounded, then we are quite desperate and still give the NLP a try
       */
      if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      {
         if( SCIPgetNNodes(scip) > 1 || SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
         {
            *result = SCIP_DELAYED;
            SCIPdebugMessage("NLP heuristic delayed because no start candidate given and no LP solution available; LP status = %d\n", SCIPgetLPSolstat(scip));
            return SCIP_OKAY;
         }
         else
         {
            SCIPdebugMessage("LP is unbounded in root node, so we are quite desperate; run NLP heuristic and pray\n");
         }
      }
      else if( SCIPgetNLPBranchCands(scip) > 0 )
      {
         /* only call heuristic, if there are no fractional variables */
         *result = SCIP_DELAYED;
         SCIPdebugMessage("NLP heuristic delayed because no start candidate given and current LP solution is fractional\n");
         return SCIP_OKAY;
      }
      else if( !SCIPisInfinity(scip, SCIPgetPrimalbound(scip)) && SCIPisEQ(scip, SCIPgetLocalDualbound(scip), SCIPgetPrimalbound(scip)) )
      {
         /* only call heuristic, if there is still room for improvement in the current node */
         SCIPdebugMessage("NLP heuristic delayed because lower and upper bound coincide in current node\n");
         return SCIP_OKAY;
      }
   }
   else
   {
      SCIPdebugMessage("have startcand from heur %s\n", SCIPsolGetHeur(heurdata->startcand) ? SCIPheurGetName(SCIPsolGetHeur(heurdata->startcand)) : "NULL");
   }

   if( !heurdata->runalways )
   {
      /* check if enough nodes have been processed so that we want to run the heuristic again */

      /* compute the contingent on number of iterations that the NLP solver is allowed to use
       * we make it depending on the current number of processed nodes
       */
      itercontingent = (SCIP_Longint)(heurdata->iterquot * SCIPgetNNodes(scip));

      /* weight by previous success of heuristic */
      itercontingent = (SCIP_Longint)(itercontingent * 3.0 * (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
      /* add the fixed offset */
      itercontingent += heurdata->iteroffset;
      /* subtract the number of iterations used so far */
      itercontingent -= heurdata->iterused;

      if( itercontingent < heurdata->itermin )
      {
         /* not enough iterations left to start NLP solver */
         SCIPdebugMessage("skip NLP heuristic; contingent=%"SCIP_LONGINT_FORMAT"; minimal number of iterations=%d; success ratio=%g\n",
            itercontingent, heurdata->itermin, (SCIPheurGetNBestSolsFound(heur)+1.0)/(SCIPheurGetNCalls(heur) + 1.0));
         return SCIP_OKAY;
      }

      /* enforce user given iteration limit, if given */
      if( heurdata->nlpiterlimit > 0 )
         itercontingent = MIN(itercontingent, heurdata->nlpiterlimit);
   }
   else
   {
      itercontingent = -1;
   }

   /* check whether there is enough time left */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 0.0 )
      {
         SCIPdebugMessage("skip NLP heuristic; no time left\n");
         return SCIP_OKAY;
      }
   }
   /* enforce user given time limit, if given */
   if( heurdata->nlptimelimit > 0 )
      timelimit = MIN(heurdata->nlptimelimit, timelimit);

   /* so far we have not found any solution, but now we are willing to search for one */
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPapplyHeurSubNlp(scip, heur, result, heurdata->startcand, itercontingent, timelimit, heurdata->minimprove, &iterused) );
   heurdata->iterused += iterused;

   /* SCIP does not like cutoff as return, so we say didnotfind, since we did not find a solution */
   if( *result == SCIP_CUTOFF )
      *result = SCIP_DIDNOTFIND;

   /* forget startcand */
   if( heurdata->startcand != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
   }

   /* reset timing, if it was changed temporary (at the root node) */
   if( heurtiming != HEUR_TIMING )
      SCIPheurSetTimingmask(heur, HEUR_TIMING);

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the NLP local search primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSubNlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create Nlp primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );
   BMSclearMemory(heurdata);

   /* include variable event handler */
   SCIP_CALL( SCIPincludeEventhdlr(scip, HEUR_NAME, "propagates a global bound change to the sub-SCIP",
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, processVarEvent, NULL) );
   heurdata->eventhdlr = SCIPfindEventhdlr(scip, HEUR_NAME);

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopySubNlp,
         heurFreeSubNlp, heurInitSubNlp, heurExitSubNlp, heurInitsolSubNlp, heurExitsolSubNlp, heurExecSubNlp,
         heurdata) );

   /* add Nlp primal heuristic parameters */
   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpverblevel",
         "verbosity level of NLP solver",
         &heurdata->nlpverblevel, FALSE, 0, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/nlpiterlimit",
         "iteration limit of NLP solver; 0 to use solver default",
         &heurdata->nlpiterlimit, FALSE, 0, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/nlptimelimit",
         "time limit of NLP solver; 0 to use solver default",
         &heurdata->nlptimelimit, FALSE, 0.0, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip, "heuristics/"HEUR_NAME"/nlpoptfile",
         "name of an NLP solver specific options file",
         &heurdata->nlpoptfile, TRUE, "", NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/resolvetolfactor",
         "if SCIP does not accept a NLP feasible solution, resolve NLP with feas. tolerance reduced by this factor (set to 1.0 to turn off resolve)",
         &heurdata->resolvetolfactor, TRUE, 0.001, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/"HEUR_NAME"/resolvefromscratch",
         "should the NLP resolve be started from the original starting point or the infeasible solution?",
         &heurdata->resolvefromscratch, TRUE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/iteroffset",
         "number of iterations added to the contingent of the total number of iterations",
         &heurdata->iteroffset, FALSE, 500, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/iterquotient",
         "contingent of NLP iterations in relation to the number of nodes in SCIP",
         &heurdata->iterquot,  FALSE, 0.1, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam (scip, "heuristics/"HEUR_NAME"/itermin",
         "contingent of NLP iterations in relation to the number of nodes in SCIP",
         &heurdata->itermin,   FALSE, 300, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam (scip, "heuristics/"HEUR_NAME"/runalways",
         "whether to run NLP heuristic always if starting point available (does not use iteroffset,iterquot,itermin)",
         &heurdata->runalways, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "heuristics/"HEUR_NAME"/minimprove",
         "factor by which NLP heuristic should at least improve the incumbent",
         &heurdata->minimprove, TRUE, 0.01, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/"HEUR_NAME"/maxpresolverounds",
         "limit on number of presolve rounds in sub-SCIP (-1 for unlimited, 0 for no presolve)",
         &heurdata->maxpresolverounds, TRUE, -1, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam (scip, "heuristics/"HEUR_NAME"/forbidfixings",
         "whether to add constraints that forbid specific fixings that turned out to be infeasible",
         &heurdata->forbidfixings, FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam (scip, "heuristics/"HEUR_NAME"/keepcopy",
         "whether to keep SCIP copy or to create new copy each time heuristic is applied",
         &heurdata->keepcopy, TRUE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}

/** adds all known linear constraint to the NLP, if initialized and not done already
 * This function is temporary and will hopefully become obsolete in the near future.
 */ 
SCIP_RETCODE SCIPaddLinearConsToNlpHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur,               /**< heuristic data structure                                       */
   SCIP_Bool             addcombconss,       /**< whether to add combinatorial linear constraints, i.e., linear constraints that involve only discrete variables */
   SCIP_Bool             addcontconss        /**< whether to add continuous    linear constraints, i.e., linear constraints that involve not only discrete variables */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* return, if nothing to do */
   if( (!addcombconss || heurdata->comblinearconsadded) && (!addcontconss || heurdata->contlinearconsadded) )
      return SCIP_OKAY;

   SCIP_CALL( addLinearConstraintsToNlp(scip,
         addcombconss && !heurdata->comblinearconsadded,
         addcontconss && !heurdata->contlinearconsadded) );

   return SCIP_OKAY;
}

/** updates the starting point for the NLP heuristic
 * 
 * Is called by a constraint handler that handles nonlinear constraints when a check on feasibility of a solution fails.
 */
SCIP_RETCODE SCIPupdateStartpointHeurSubNlp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< NLP heuristic */
   SCIP_SOL*             solcand,            /**< solution candidate */
   SCIP_Real             violation           /**< constraint violation of solution candidate */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(solcand != NULL);
   assert(SCIPisPositive(scip, violation));

   /* too early or the game is over already: no more interest in starting points */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   /* we do not have a sub-SCIP, so we also do not need a starting point */
   if( heurdata->subscip == NULL )
      return SCIP_OKAY;

   /* if the solution is from our heuristic, then it is useless to use it as starting point again */
   if( SCIPsolGetHeur(solcand) == heur )
      return SCIP_OKAY;

   SCIPdebugMessage("consider solution candidate with violation %g and objective %g from %s\n",
      violation, SCIPgetSolTransObj(scip, solcand), SCIPsolGetHeur(solcand) ? SCIPheurGetName(SCIPsolGetHeur(solcand)) : "tree");

   /* if we have no point yet, or the new point has a lower constraint violation, or it has a better objective function value, then take the new point */
   if( heurdata->startcand == NULL || SCIPisGT(scip, heurdata->startcandviol, violation) ||
      SCIPisRelGT(scip, SCIPgetSolTransObj(scip, heurdata->startcand), SCIPgetSolTransObj(scip, solcand)) )
   {
      if( heurdata->startcand != NULL )
      {
         SCIP_CALL( SCIPfreeSol(scip, &heurdata->startcand) );
      }
      SCIP_CALL( SCIPcreateSolCopy(scip, &heurdata->startcand, solcand) );
      SCIP_CALL( SCIPunlinkSol(scip, heurdata->startcand) );
      heurdata->startcandviol = violation;
   }

   return SCIP_OKAY;
}

/** gets sub-SCIP used by NLP heuristic, or NULL if none */
SCIP* SCIPgetSubScipHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->subscip;
}

/** gets mapping of SCIP variables to sub-SCIP variables */
SCIP_VAR** SCIPgetVarMappingScip2SubScipHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->var_scip2subscip;
}

/** gets mapping of sub-SCIP variables to SCIP variables */
SCIP_VAR** SCIPgetVarMappingSubScip2ScipHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->var_subscip2scip;
}

/** gets startpoint candidate to be used in next call to NLP heuristic, or NULL if none */
SCIP_SOL* SCIPgetStartCandidateHeurSubNlp(
   SCIP*                 scip,               /**< original SCIP data structure                                   */
   SCIP_HEUR*            heur                /**< heuristic data structure                                       */
   )
{
   SCIP_HEURDATA* heurdata;

   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   return heurdata->startcand;
}
