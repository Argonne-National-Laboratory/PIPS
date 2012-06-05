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

/**@file   scip.c
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 *
 * @todo check all checkStage() calls, use bit flags instead of the SCIP_Bool parameters
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE 
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <ctype.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>

#ifndef NPARASCIP
#include <pthread.h>
#include <time.h>
#endif

#ifdef WITH_ZLIB
#include <zlib.h>
#endif

#include "scip/def.h"
#include "scip/retcode.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/vbc.h"
#include "scip/interrupt.h"
#include "scip/lpi.h"
#include "scip/mem.h"
#include "scip/misc.h"
#include "scip/history.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/nlp.h"
#include "scip/var.h"
#include "scip/implics.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/pricestore.h"
#include "scip/sepastore.h"
#include "scip/cutpool.h"
#include "scip/solve.h"
#include "scip/scipgithash.h"
#include "scip/scip.h"

#include "scip/branch.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/dialog.h"
#include "scip/disp.h"
#include "scip/heur.h"
#include "scip/nodesel.h"
#include "scip/reader.h"
#include "scip/presol.h"
#include "scip/pricer.h"
#include "scip/relax.h"
#include "scip/sepa.h"
#include "scip/prop.h"
#include "nlpi/nlpi.h"
#include "nlpi/exprinterpret.h"
#include "scip/debug.h"
#include "scip/dialog_default.h"

/* We include the linear constraint handler to be able to copy a (multi)aggregation of variables (to a linear constraint).
 * The better way would be to handle the distinction between original and transformed variables via a flag 'isoriginal' 
 * in the variable data structure. This would allow to have (multi)aggregated variables in the original problem.
 *
 * A second reason for including the linear constraint handler is for copying cuts to linear constraints.
 */ 
#include "scip/cons_linear.h"

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons)
 */
#ifndef NDEBUG
#include "scip/struct_scip.h"
#endif


/* 
 * Local methods
 */


/** checks, if SCIP is in one of the feasible stages */
#ifndef NDEBUG
static
SCIP_RETCODE checkStage(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           method,             /**< method that was called */
   SCIP_Bool             init,               /**< may method be called in the INIT stage? */
   SCIP_Bool             problem,            /**< may method be called in the PROBLEM stage? */
   SCIP_Bool             transforming,       /**< may method be called in the TRANSFORMING stage? */
   SCIP_Bool             transformed,        /**< may method be called in the TRANSFORMED stage? */
   SCIP_Bool             presolving,         /**< may method be called in the PRESOLVING stage? */
   SCIP_Bool             presolved,          /**< may method be called in the PRESOLVED stage? */
   SCIP_Bool             initsolve,          /**< may method be called in the INITSOLVE stage? */
   SCIP_Bool             solving,            /**< may method be called in the SOLVING stage? */
   SCIP_Bool             solved,             /**< may method be called in the SOLVED stage? */
   SCIP_Bool             freesolve,          /**< may method be called in the FREESOLVE stage? */
   SCIP_Bool             freetrans           /**< may method be called in the FREETRANS stage? */
   )
{
   assert(scip != NULL);
   assert(method != NULL);
   
   /*SCIPdebugMessage("called method <%s> at stage %d ------------------------------------------------\n",
     method, scip->set->stage);*/
   
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->interrupt != NULL);
   assert(scip->dialoghdlr != NULL);
   assert(scip->totaltime != NULL);
   
   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
      assert(scip->stat == NULL);
      assert(scip->origprob == NULL);
      assert(scip->eventfilter == NULL);
      assert(scip->eventqueue == NULL);
      assert(scip->branchcand == NULL);
      assert(scip->lp == NULL);
      assert(scip->nlp == NULL);
      assert(scip->primal == NULL);
      assert(scip->tree == NULL);
      assert(scip->conflict == NULL);
      assert(scip->transprob == NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);
      
      if( !init )
      {
         SCIPerrorMessage("cannot call method <%s> in initialization stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;
      
   case SCIP_STAGE_PROBLEM:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter == NULL);
      assert(scip->eventqueue == NULL);
      assert(scip->branchcand == NULL);
      assert(scip->lp == NULL);
      assert(scip->nlp == NULL);
      assert(scip->primal == NULL);
      assert(scip->tree == NULL);
      assert(scip->conflict == NULL);
      assert(scip->transprob == NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);

      if( !problem )
      {
         SCIPerrorMessage("cannot call method <%s> in problem creation stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->conflict != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);

      if( !transforming )
      {
         SCIPerrorMessage("cannot call method <%s> in problem transformation stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->conflict != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);

      if( !transformed )
      {
         SCIPerrorMessage("cannot call method <%s> in problem transformed stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->conflict != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);

      if( !presolving )
      {
         SCIPerrorMessage("cannot call method <%s> in presolving stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->conflict != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);

      if( !presolved )
      {
         SCIPerrorMessage("cannot call method <%s> in problem presolved stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_INITSOLVE:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->transprob != NULL);

      if( !initsolve )
      {
         SCIPerrorMessage("cannot call method <%s> in init solve stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->conflict != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore != NULL);
      assert(scip->sepastore != NULL);
      assert(scip->cutpool != NULL);

      if( !solving )
      {
         SCIPerrorMessage("cannot call method <%s> in solving stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVED:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->conflict != NULL);
      assert(scip->transprob != NULL);
      assert(scip->pricestore != NULL);
      assert(scip->sepastore != NULL);
      assert(scip->cutpool != NULL);

      if( !solved )
      {
         SCIPerrorMessage("cannot call method <%s> in problem solved stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREESOLVE:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->eventfilter != NULL);
      assert(scip->eventqueue != NULL);
      assert(scip->branchcand != NULL);
      assert(scip->lp != NULL);
      assert(scip->primal != NULL);
      assert(scip->tree != NULL);
      assert(scip->transprob != NULL);

      if( !freesolve )
      {
         SCIPerrorMessage("cannot call method <%s> in solve deinitialization stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREETRANS:
      assert(scip->stat != NULL);
      assert(scip->origprob != NULL);
      assert(scip->pricestore == NULL);
      assert(scip->sepastore == NULL);
      assert(scip->cutpool == NULL);

      if( !freetrans )
      {
         SCIPerrorMessage("cannot call method <%s> in free transformed problem stage\n", method);
         return SCIP_INVALIDCALL;
      }
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }
}
#else
#define checkStage(scip,method,init,problem,transforming,transformed,presolving,presolved,initsolve,solving,solved, \
   freesolve,freetrans) SCIP_OKAY
#endif


/** gets global primal bound (objective value of best solution or user objective limit) */
static
SCIP_Real getPrimalbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIPprobExternObjval(scip->transprob, scip->set, scip->primal->upperbound);
}

/** gets global dual bound */
static
SCIP_Real getDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real lowerbound;

   lowerbound = SCIPtreeGetLowerbound(scip->tree, scip->set);

   if( SCIPsetIsInfinity(scip->set, lowerbound) )
   {
      /* in case we could not prove whether the problem is unbounded or infeasible, we want to terminate with 
       * dual bound = -inf instead of dual bound = primal bound = +inf
       * also in case we prove that the problem is unbounded, it seems to make sense to return with dual bound = -inf,
       * since -infinity is the only valid lower bound
       */
      if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD || SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED )
         return SCIPprobExternObjval(scip->transprob, scip->set, -SCIPinfinity(scip));
      else
         return getPrimalbound(scip);
   }
   else
      return SCIPprobExternObjval(scip->transprob, scip->set, lowerbound);
}

/** gets global lower (dual) bound in transformed problem */
static
SCIP_Real getLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return SCIPtreeGetLowerbound(scip->tree, scip->set);
}

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit) */
static
SCIP_Real getUpperbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   return scip->primal->upperbound;
}


/*
 * miscellaneous methods
 */

/** returns SCIP version number */
SCIP_Real SCIPversion(
   void
   )
{
   return (SCIP_Real)(SCIP_VERSION)/100.0;
}

/** returns SCIP major version */
int SCIPmajorVersion(
   void
   )
{
   return SCIP_VERSION/100;
}

/** returns SCIP minor version */
int SCIPminorVersion(
   void
   )
{
   return (SCIP_VERSION/10) % 10;
}

/** returns SCIP technical version */
int SCIPtechVersion(
   void
   )
{
   return SCIP_VERSION % 10;
}

/** returns SCIP sub version number */
int SCIPsubversion(
   void
   )
{
   return SCIP_SUBVERSION;
}

/** prints a version information line to a file stream */
void SCIPprintVersion(
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIPmessageFPrintInfo(file, "SCIP version %d.%d.%d",
      SCIPmajorVersion(), SCIPminorVersion(), SCIPtechVersion());
#if SCIP_SUBVERSION > 0
   SCIPmessageFPrintInfo(file, ".%d", SCIPsubversion());
#endif
   
   SCIPmessageFPrintInfo(file, " [precision: %d byte]", (int)sizeof(SCIP_Real));

#ifndef BMS_NOBLOCKMEM
   SCIPmessageFPrintInfo(file, " [memory: block]");
#else
   SCIPmessageFPrintInfo(file, " [memory: standard]");
#endif
#ifndef NDEBUG
   SCIPmessageFPrintInfo(file, " [mode: debug]");
#else
   SCIPmessageFPrintInfo(file, " [mode: optimized]");
#endif
   SCIPmessageFPrintInfo(file, " [LP solver: %s]", SCIPlpiGetSolverName());
   SCIPmessageFPrintInfo(file, " [GitHash: %s]", SCIPgetGitHash());
   SCIPmessageFPrintInfo(file, "\n");
   SCIPmessageFPrintInfo(file, "%s\n", SCIP_COPYRIGHT);
}

/** prints error message for the given SCIP return code */
void SCIPprintError(
   SCIP_RETCODE          retcode,            /**< SCIP return code causing the error */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIPmessageFPrintInfo(file, "SCIP Error (%d): ", retcode);
   SCIPretcodePrint(file, retcode);
   SCIPmessageFPrintInfo(file, "\n");
}


/*
 * general SCIP methods
 */

/** creates and initializes SCIP data structures */
SCIP_RETCODE SCIPcreate(
   SCIP**                scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_ALLOC( BMSallocMemory(scip) );

   SCIP_CALL( SCIPmemCreate(&(*scip)->mem) );
   SCIP_CALL( SCIPsetCreate(&(*scip)->set, (*scip)->mem->setmem, *scip) );
   SCIP_CALL( SCIPinterruptCreate(&(*scip)->interrupt) );
   SCIP_CALL( SCIPdialoghdlrCreate((*scip)->set, &(*scip)->dialoghdlr) );
   SCIP_CALL( SCIPclockCreate(&(*scip)->totaltime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIPclockStart((*scip)->totaltime, (*scip)->set);
   (*scip)->stat = NULL;
   (*scip)->origprob = NULL;
   (*scip)->origprimal = NULL;
   (*scip)->eventfilter = NULL;
   (*scip)->eventqueue = NULL;
   (*scip)->branchcand = NULL;
   (*scip)->lp = NULL;
   (*scip)->nlp = NULL;
   (*scip)->primal = NULL;
   (*scip)->tree = NULL;
   (*scip)->conflict = NULL;
   (*scip)->transprob = NULL;
   (*scip)->pricestore = NULL;
   (*scip)->sepastore = NULL;
   (*scip)->cutpool = NULL;

   SCIP_CALL( SCIPnlpInclude((*scip)->set, SCIPblkmem(*scip)) );

   if( strcmp(SCIPlpiGetSolverName(), "NONE") != 0 )
   {
      SCIP_CALL( SCIPsetIncludeExternalCode((*scip)->set, SCIPlpiGetSolverName(), SCIPlpiGetSolverDesc()) );
   }
   if( strcmp(SCIPexprintGetName(), "NONE") != 0 )
   {
      SCIP_CALL( SCIPsetIncludeExternalCode((*scip)->set, SCIPexprintGetName(), SCIPexprintGetDesc()) );
   }

#ifdef WITH_ZLIB
   SCIP_CALL( SCIPsetIncludeExternalCode((*scip)->set, "ZLIB " ZLIB_VERSION, "General purpose compression library by J. Gailly and M. Adler (zlib.net)") );
#endif

   return SCIP_OKAY;
}

/** frees SCIP data structures */
SCIP_RETCODE SCIPfree(
   SCIP**                scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( checkStage(*scip, "SCIPfree", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPfreeProb(*scip) );
   assert((*scip)->set->stage == SCIP_STAGE_INIT);

   SCIP_CALL( SCIPsetFree(&(*scip)->set, (*scip)->mem->setmem) );
   SCIP_CALL( SCIPdialoghdlrFree(*scip, &(*scip)->dialoghdlr) );
   SCIPclockFree(&(*scip)->totaltime);
   SCIPinterruptFree(&(*scip)->interrupt);
   SCIP_CALL( SCIPmemFree(&(*scip)->mem) );

   BMSfreeMemory(scip);

   return SCIP_OKAY;
}

/** returns current stage of SCIP */
SCIP_STAGE SCIPgetStage(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->stage;
}

/** outputs SCIP stage and solution status if applicable */
SCIP_RETCODE SCIPprintStage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPprintStage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
      SCIPmessageFPrintInfo(file, "initialization");
      break;
   case SCIP_STAGE_PROBLEM:
      SCIPmessageFPrintInfo(file, "problem creation / modification");
      break;
   case SCIP_STAGE_TRANSFORMING:
      SCIPmessageFPrintInfo(file, "problem transformation");
      break;
   case SCIP_STAGE_TRANSFORMED:
      SCIPmessageFPrintInfo(file, "problem transformed");
      break;
   case SCIP_STAGE_PRESOLVING:
      if( SCIPsolveIsStopped(scip->set, scip->stat, TRUE) )
      {
         SCIPmessageFPrintInfo(file, "solving was interrupted [");
         SCIP_CALL( SCIPprintStatus(scip, file) );
         SCIPmessageFPrintInfo(file, "]");
      }
      else
         SCIPmessageFPrintInfo(file, "presolving process is running");
      break;
   case SCIP_STAGE_PRESOLVED:
      SCIPmessageFPrintInfo(file, "problem is presolved");
      break;
   case SCIP_STAGE_INITSOLVE:
      SCIPmessageFPrintInfo(file, "solving process initialization");
      break;
   case SCIP_STAGE_SOLVING:
      if( SCIPsolveIsStopped(scip->set, scip->stat, TRUE) )
      {
         SCIPmessageFPrintInfo(file, "solving was interrupted [");
         SCIP_CALL( SCIPprintStatus(scip, file) );
         SCIPmessageFPrintInfo(file, "]");
      }
      else
         SCIPmessageFPrintInfo(file, "solving process is running");
      break;
   case SCIP_STAGE_SOLVED:
      SCIPmessageFPrintInfo(file, "problem is solved [");
      SCIP_CALL( SCIPprintStatus(scip, file) );
      SCIPmessageFPrintInfo(file, "]");
      break;
   case SCIP_STAGE_FREESOLVE:
      SCIPmessageFPrintInfo(file, "solving process deinitialization");
      break;
   case SCIP_STAGE_FREETRANS:
      SCIPmessageFPrintInfo(file, "freeing transformed problem");
      break;
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** gets solution status */
SCIP_STATUS SCIPgetStatus(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetStatus", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( scip->set->stage == SCIP_STAGE_INIT )
      return SCIP_STATUS_UNKNOWN;
   else
   {
      assert(scip->stat != NULL);

      return scip->stat->status;
   }
}

/** outputs solution status */
SCIP_RETCODE SCIPprintStatus(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPprintStatus", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( SCIPgetStatus(scip) )
   {
   case SCIP_STATUS_UNKNOWN:
      SCIPmessageFPrintInfo(file, "unknown");
      break;
   case SCIP_STATUS_USERINTERRUPT:
      SCIPmessageFPrintInfo(file, "user interrupt");
      break;
   case SCIP_STATUS_NODELIMIT:
      SCIPmessageFPrintInfo(file, "node limit reached");
      break;
   case SCIP_STATUS_STALLNODELIMIT:
      SCIPmessageFPrintInfo(file, "stall node limit reached");
      break;
   case SCIP_STATUS_TIMELIMIT:
      SCIPmessageFPrintInfo(file, "time limit reached");
      break;
   case SCIP_STATUS_MEMLIMIT:
      SCIPmessageFPrintInfo(file, "memory limit reached");
      break;
   case SCIP_STATUS_GAPLIMIT:
      SCIPmessageFPrintInfo(file, "gap limit reached");
      break;
   case SCIP_STATUS_SOLLIMIT:
      SCIPmessageFPrintInfo(file, "solution limit reached");
      break;
   case SCIP_STATUS_BESTSOLLIMIT:
      SCIPmessageFPrintInfo(file, "solution improvement limit reached");
      break;
   case SCIP_STATUS_OPTIMAL:
      SCIPmessageFPrintInfo(file, "optimal solution found");
      break;
   case SCIP_STATUS_INFEASIBLE:
      SCIPmessageFPrintInfo(file, "infeasible");
      break;
   case SCIP_STATUS_UNBOUNDED:
      SCIPmessageFPrintInfo(file, "unbounded");
      break;
   case SCIP_STATUS_INFORUNBD:
      SCIPmessageFPrintInfo(file, "infeasible or unbounded");
      break;
   default:
      SCIPerrorMessage("invalid status code <%d>\n", SCIPgetStatus(scip));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** returns whether the current stage belongs to the transformed problem space */
SCIP_Bool SCIPisTransformed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return ((int)scip->set->stage >= (int)SCIP_STAGE_TRANSFORMING);
}

/** returns whether the solution process should be provably correct */
SCIP_Bool SCIPisExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return (scip->set->misc_exactsolve);
}

/** returns whether the user pressed CTRL-C to interrupt the solving process */
SCIP_Bool SCIPpressedCtrlC(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPpressedCtrlC", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPinterrupted();
}

/** returns whether the solving process should be / was stopped before proving optimality;
 *  if the solving process should be / was stopped, the status returned by SCIPgetStatus() yields
 *  the reason for the premature abort
 */
SCIP_Bool SCIPisStopped(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisStopped", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolveIsStopped(scip->set, scip->stat, FALSE);
}




/*
 * message output methods
 */

/** creates a message handler; this method can already be called before SCIPcreate() */
SCIP_RETCODE SCIPcreateMessagehdlr(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   SCIP_DECL_MESSAGEERROR((*messageerror)),  /**< error message print method of message handler */
   SCIP_DECL_MESSAGEWARNING((*messagewarning)),/**< warning message print method of message handler */
   SCIP_DECL_MESSAGEDIALOG((*messagedialog)),/**< dialog message print method of message handler */
   SCIP_DECL_MESSAGEINFO ((*messageinfo)),   /**< info message print method of message handler */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< message handler data */
   )
{
   SCIP_CALL( SCIPmessagehdlrCreate(messagehdlr, bufferedoutput,
         messageerror, messagewarning, messagedialog, messageinfo, messagehdlrdata) );

   return SCIP_OKAY;
}

/** frees message handler; this method can be called after SCIPfree() */
SCIP_RETCODE SCIPfreeMessagehdlr(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   SCIPmessagehdlrFree(messagehdlr);

   return SCIP_OKAY;
}

/** installs the given message handler, such that all messages are passed to this handler;
 *  this method can already be called before SCIPcreate()
 */
SCIP_RETCODE SCIPsetMessagehdlr(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   )
{
   SCIPmessageSetHandler(messagehdlr);

   return SCIP_OKAY;
}

/** installs the default message handler, such that all messages are printed to stdout and stderr;
 *  this method can already be called before SCIPcreate()
 */
SCIP_RETCODE SCIPsetDefaultMessagehdlr(
   void
   )
{
   SCIPmessageSetDefaultHandler();

   return SCIP_OKAY;
}

/** returns the currently installed message handler, or NULL if messages are currently suppressed;
 *  this method can already be called before SCIPcreate()
 */
SCIP_MESSAGEHDLR* SCIPgetMessagehdlr(
   void
   )
{
   return SCIPmessageGetHandler();
}

/** prints a dialog message that requests user interaction or is a direct response to a user interactive command */
void SCIPdialogMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPdialogMessage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintDialog(file, formatstr, ap);
   va_end(ap);
}

/** prints a message */
void SCIPinfoMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPinfoMessage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(file, formatstr, ap);
   va_end(ap);
}

/** prints a message depending on the verbosity level */
void SCIPverbMessage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPverbMessage", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintVerbInfo(scip->set->disp_verblevel, msgverblevel, file, formatstr, ap);
   va_end(ap);
}

/** returns the current message verbosity level */
SCIP_VERBLEVEL SCIPgetVerbLevel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVerbLevel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->disp_verblevel;
}

#ifndef NPARASCIP
/** allocates memory for all message handlers for number of given threads 
 *
 * @note it is necessary to call SCIPmessageSetDefaultHandler or SCIPmessageSetHandler for each thread
 */
SCIP_RETCODE SCIPcreateMesshdlrPThreads(
   int                   nthreads            /**< number of threads to allocate memory for */
   )
{
   assert(nthreads > 0);

   return SCIPmesshdlrCreatePThreads(nthreads);
}

/** frees memory for all message handlers */
void SCIPfreeMesshdlrPThreads(
   void
   )
{
   SCIPmesshdlrFreePThreads();
}
#endif


/*
 * SCIP copy methods
 */


#define HASHTABLESIZE_FACTOR 5

/** copies plugins from sourcescip to targetscip; in case that a constraint handler which does not need constraints
 *  cannot be copied, valid will return FALSE. All plugins can declare that, if their copy process failed, the 
 *  copied SCIP instance might not represent the same problem semantics as the original. 
 *  Note that in this case dual reductions might be invalid. 
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPcopyPlugins(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_Bool             copyreaders,        /**< should the file readers be copied */
   SCIP_Bool             copypricers,        /**< should the variable pricers be copied */
   SCIP_Bool             copyconshdlrs,      /**< should the constraint handlers be copied */
   SCIP_Bool             copyconflicthdlrs,  /**< should the conflict handlers be copied */
   SCIP_Bool             copypresolvers,     /**< should the presolvers be copied */
   SCIP_Bool             copyrelaxators,     /**< should the relaxation handlers  be copied */
   SCIP_Bool             copyseparators,     /**< should the separators be copied */
   SCIP_Bool             copypropagators,    /**< should the propagators be copied */
   SCIP_Bool             copyheuristics,     /**< should the heuristics be copied */
   SCIP_Bool             copyeventhdlrs,     /**< should the event handlers be copied */
   SCIP_Bool             copynodeselectors,  /**< should the node selectors be copied */
   SCIP_Bool             copybranchrules,    /**< should the branchrules be copied */
   SCIP_Bool             copydisplays,       /**< should the display columns be copied */
   SCIP_Bool             copydialogs,        /**< should the dialogs be copied */
   SCIP_Bool             copynlpis,          /**< should the NLPIs be copied */
   SCIP_Bool*            valid               /**< pointer to store whether plugins, in particular all constraint
                                              *   handlers which do not need constraints were validly copied */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourcescip->set != NULL);
   assert(targetscip->set != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPcopyPlugins", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPcopyPlugins", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsetCopyPlugins(sourcescip->set, targetscip->set, 
         copyreaders, copypricers, copyconshdlrs, copyconflicthdlrs, copypresolvers, copyrelaxators, copyseparators, copypropagators,
         copyheuristics, copyeventhdlrs, copynodeselectors, copybranchrules, copydisplays, copydialogs, copynlpis, valid) );

   return SCIP_OKAY;
}

/** create a problem by copying the problem data of the source SCIP 
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPcopyProb(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL  */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   const char*           name                /**< problem name of target */
   )
{
   SCIP_PROB* sourceprob;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(targetscip, "SCIPcopyProb", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   SCIP_CALL( checkStage(sourcescip, "SCIPcopyProb", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   /* free old problem */
   SCIP_CALL( SCIPfreeProb(targetscip) );
   assert(targetscip->set->stage == SCIP_STAGE_INIT);

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNVars(sourcescip))) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNConss(sourcescip))) );
   }
   else
      localconsmap = consmap;

   /* switch stage to PROBLEM */
   targetscip->set->stage = SCIP_STAGE_PROBLEM;

   if( SCIPisTransformed(sourcescip) )
      sourceprob = sourcescip->transprob;
   else
      sourceprob = sourcescip->origprob;

   /* create the statistics data structure */
   SCIP_CALL( SCIPstatCreate(&targetscip->stat, targetscip->mem->probmem, targetscip->set) );
   targetscip->stat->subscipdepth = sourcescip->stat->subscipdepth + 1;

   /* create the problem by copying the source problem */
   SCIP_CALL( SCIPprobCopy(&targetscip->origprob, targetscip->mem->probmem, targetscip->set, name, sourcescip, sourceprob, localvarmap, localconsmap, global) );

   /* creating the solution candidates storage */
   /**@todo copy solution of source SCIP as candidates for the target SCIP */
   SCIP_CALL( SCIPprimalCreate(&targetscip->origprimal) );

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}

/** returns copy of the source variable; if there already is a copy of the source variable in the variable hash map,
 *  it is just returned as target variable; elsewise a new variable will be created and added to the target SCIP; this
 *  created variable is added to the variable hash map and returned as target variable
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 *  @note if a new variable was created, this variable will be added to the target-SCIP, but it is not captured
 */
SCIP_RETCODE SCIPgetVarCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_VAR*             sourcevar,          /**< source variable */
   SCIP_VAR**            targetvar,          /**< pointer to store the target variable */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< should global or local bounds be used? */
   SCIP_Bool*            success             /**< pointer to store whether the copying was successful or not */
   )
{
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_VAR* var;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPgetVarCopy", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPgetVarCopy", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(targetscip != NULL);
   assert(sourcevar != NULL);
   assert(targetvar != NULL);

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);
   *success = TRUE;

   /* try to retrieve copied variable from hashmap */
   if( !uselocalvarmap )
   {
      *targetvar = (SCIP_VAR*) SCIPhashmapGetImage(varmap, sourcevar);
      if( *targetvar != NULL )
         return SCIP_OKAY;
   }

   /* if the target SCIP is already in solving stage we currently are not copying the variable!
    * this has to be done because we cannot simply add variables to SCIP during solving and thereby enlarge the search
    * space. 
    * unlike column generation we cannot assume here that the variable could be implicitly set to zero in all prior
    * computations
    */
   if( SCIPgetStage(targetscip) > SCIP_STAGE_PROBLEM )
   {
      *success = FALSE;
      *targetvar = NULL;

      return SCIP_OKAY;
   }

   /* create the variable mapping hash map */
   if( uselocalvarmap )
   {
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNVars(sourcescip))) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNConss(sourcescip))) );
   }
   else
      localconsmap = consmap;
   
   /* if variable does not exists yet in target SCIP, create it */
   switch( SCIPvarGetStatus(sourcevar) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_COLUMN:
   case SCIP_VARSTATUS_LOOSE:
   case SCIP_VARSTATUS_FIXED:
      SCIP_CALL( SCIPvarCopy(&var, targetscip->mem->probmem, targetscip->set, targetscip->stat, 
            sourcescip, sourcevar, localvarmap, localconsmap, global) );
      break;
      
   case SCIP_VARSTATUS_AGGREGATED:
   {
      SCIP_CONS* cons;
      char name[SCIP_MAXSTRLEN];
         
      SCIP_VAR* sourceaggrvar;
      SCIP_VAR* targetaggrvar;
      SCIP_Real aggrcoef;
      SCIP_Real constant;

      /* get aggregation data */
      sourceaggrvar = SCIPvarGetAggrVar(sourcevar);
      aggrcoef = SCIPvarGetAggrScalar(sourcevar);
      constant = SCIPvarGetAggrConstant(sourcevar);

      /* get copy of the aggregation variable */
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourceaggrvar, &targetaggrvar, localvarmap, localconsmap, global, success) );
      assert(*success);

      /* create copy of the aggregated variable */
      SCIP_CALL( SCIPvarCopy(&var, targetscip->mem->probmem, targetscip->set, targetscip->stat, 
            sourcescip, sourcevar, localvarmap, localconsmap, global) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_aggr", SCIPvarGetName(sourcevar));
         
      /* add aggregation x = a*y + c as linear constraint x - a*y = c */
      SCIP_CALL( SCIPcreateConsLinear(targetscip, &cons, name, 0, NULL, NULL, constant,
            constant, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCoefLinear(targetscip, cons, var, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(targetscip, cons, targetaggrvar, -aggrcoef) );

      SCIP_CALL( SCIPaddCons(targetscip, cons) );
      SCIP_CALL( SCIPreleaseCons(targetscip, &cons) );               
         
      break;
   }
   case SCIP_VARSTATUS_MULTAGGR:
   {
      SCIP_CONS* cons;
      char name[SCIP_MAXSTRLEN];
         
      SCIP_VAR** sourceaggrvars;
      SCIP_VAR** targetaggrvars;
      SCIP_Real* aggrcoefs;
      SCIP_Real constant;
         
      int naggrvars;
      int i;
        
      /* get the active representation */
      SCIP_CALL( SCIPflattenVarAggregationGraph(sourcescip, sourcevar) );
    
      /* get multi-aggregation data */
      naggrvars = SCIPvarGetMultaggrNVars(sourcevar);
      sourceaggrvars = SCIPvarGetMultaggrVars(sourcevar);
      aggrcoefs = SCIPvarGetMultaggrScalars(sourcevar);
      constant = SCIPvarGetMultaggrConstant(sourcevar);
   
      SCIP_CALL( SCIPallocBufferArray(targetscip, &targetaggrvars, naggrvars) );
            
      /* get copies of the active variables of the multi-aggregation */
      for( i = 0; i < naggrvars; ++i )
      {
         SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourceaggrvars[i], &targetaggrvars[i], localvarmap, localconsmap, global, success) );
         assert(*success);
      }

      /* create copy of the multi-aggregated variable */
      SCIP_CALL( SCIPvarCopy(&var, targetscip->mem->probmem, targetscip->set, targetscip->stat, 
            sourcescip, sourcevar, localvarmap, localconsmap, global) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_multaggr", SCIPvarGetName(sourcevar));
         
      /* add multi-aggregation x = a^T y + c as linear constraint a^T y - x = -c */
      SCIP_CALL( SCIPcreateConsLinear(targetscip, &cons, name, naggrvars, targetaggrvars, aggrcoefs, -constant,
            -constant, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCoefLinear(targetscip, cons, var, -1.0) );
      SCIP_CALL( SCIPaddCons(targetscip, cons) );
      SCIP_CALL( SCIPreleaseCons(targetscip, &cons) );               

      SCIPfreeBufferArray(targetscip, &targetaggrvars);
            
      break;
   }
   case SCIP_VARSTATUS_NEGATED:
   {
      SCIP_VAR* sourcenegatedvar;
      SCIP_VAR* targetnegatedvar;
         
      /* get negated source variable */
      sourcenegatedvar = SCIPvarGetNegationVar(sourcevar);
      assert(sourcenegatedvar != NULL);
      assert(SCIPvarGetStatus(sourcenegatedvar) != SCIP_VARSTATUS_NEGATED);

      /* get copy of negated source variable */         
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcenegatedvar, &targetnegatedvar, localvarmap, localconsmap, global, success) );
      assert(*success);
      assert(SCIPvarGetStatus(targetnegatedvar) != SCIP_VARSTATUS_NEGATED);

      /* get negation of copied negated source variable, this is the target variable */
      SCIP_CALL( SCIPgetNegatedVar(targetscip, targetnegatedvar, targetvar) );
      assert(SCIPvarGetStatus(*targetvar) == SCIP_VARSTATUS_NEGATED);

      /* free local hash maps if necessary */
      if( uselocalvarmap )
         SCIPhashmapFree(&localvarmap);

      if( uselocalconsmap )
         SCIPhashmapFree(&localconsmap);

      /* we have to return right away, to avoid adding the negated variable to the problem since the "not negated"
       * variable was already added */
      return SCIP_OKAY;
   }
   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      break;
   }

   /* add the (new) target variable to the target problem */
   SCIP_CALL( SCIPaddVar(targetscip, var) );

   *targetvar = var;

   /* remove the variable capture which was done due to the creation of the variable */
   SCIP_CALL( SCIPreleaseVar(targetscip, &var) );

   /* free local hash maps if necessary */
   if( uselocalvarmap )
      SCIPhashmapFree(&localvarmap);

   if( uselocalconsmap )
      SCIPhashmapFree(&localconsmap);

   return SCIP_OKAY;
}

/** copies all active variables from source-SCIP and adds these variable to the target-SCIP; the mapping between these
 *  variables are stored in the variable hashmap, target-SCIP has to be in problem creation stage, fixed and aggregated
 *  variables do not get copied
 *
 *  @note the variables are added to the target-SCIP but not captured
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPcopyVars(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global              /**< should global or local bounds be used? */
   )
{
   SCIP_VAR** sourcevars;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   int nsourcevars;   
   int i;
  
   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPcopyVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPcopyVars", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* get active variables of the source SCIP */
   SCIP_CALL( SCIPgetVarsData(sourcescip, &sourcevars, &nsourcevars, NULL, NULL, NULL, NULL) );

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNVars(sourcescip))) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNConss(sourcescip))) );
   }
   else
      localconsmap = consmap;

   /* create the variables of the target SCIP */
   for( i = 0; i < nsourcevars; ++i )
   {          
      SCIP_Bool success;
      SCIP_VAR* targetvar;

      /* copy variable and add this copy to the target SCIP if the copying was valid */
      SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcevars[i], &targetvar, localvarmap, localconsmap, global, &success) );
      assert(success);
      assert(targetvar != NULL);
   }

   /* integer variables that are fixed to zero or one or have bounds [0,1] will be converted to binaries */
#ifndef NDEBUG
   {
      SCIP_VAR** fixedvars;
      int nfixedvars;
      int nfixedbinvars;
      int nfixedintvars;
      int nfixedimplvars;
      int nfixedcontvars;

      fixedvars = SCIPgetFixedVars(sourcescip);
      nfixedvars = SCIPgetNFixedVars(sourcescip);
      nfixedbinvars = 0;
      nfixedintvars = 0;
      nfixedimplvars = 0;
      nfixedcontvars = 0;
      
      /* count number of fixed variables for all variable types */
      for( i = 0; i < nfixedvars; ++i )
      {
         switch( SCIPvarGetType(fixedvars[i]) )
         {
         case SCIP_VARTYPE_BINARY:
            nfixedbinvars++;
            break;
         case SCIP_VARTYPE_INTEGER:
            nfixedintvars++;
            break;
         case SCIP_VARTYPE_IMPLINT:
            nfixedimplvars++;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            nfixedcontvars++;
            break;
         default:
            SCIPerrorMessage("unknown variable type\n");
            return SCIP_INVALIDDATA;
         }
      }
      assert(nfixedvars == nfixedbinvars + nfixedintvars + nfixedimplvars + nfixedcontvars);
      assert(SCIPgetNBinVars(sourcescip) <= SCIPgetNBinVars(targetscip));
      assert(SCIPgetNIntVars(sourcescip) + SCIPgetNBinVars(sourcescip) <= SCIPgetNIntVars(targetscip) + SCIPgetNBinVars(targetscip)
         && SCIPgetNIntVars(targetscip) + SCIPgetNBinVars(targetscip) <= SCIPgetNIntVars(sourcescip) + SCIPgetNBinVars(sourcescip) + nfixedbinvars + nfixedintvars );
      assert(SCIPgetNImplVars(sourcescip) <= SCIPgetNImplVars(targetscip) 
         && SCIPgetNImplVars(targetscip) <= SCIPgetNImplVars(sourcescip) + nfixedimplvars);     
      assert(SCIPgetNContVars(sourcescip) <= SCIPgetNContVars(targetscip)
         && SCIPgetNContVars(targetscip) <= SCIPgetNContVars(targetscip) + nfixedcontvars);
   }
#endif

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}

/** returns copy of the source constraint; if there already is a copy of the source constraint in the constraint hash
 *  map, it is just returned as target constraint; elsewise a new constraint will be created and captured in the target
 *  SCIP; this created constraint is added to the constraint hash map and returned as target constraint; the variable
 *  map is used to map the variables of the source SCIP to the variables of the target SCIP
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 *
 *  @note The constraint is not added to the target SCIP
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPgetConsCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_CONS*            sourcecons,         /**< source constraint of the source SCIP */
   SCIP_CONS**           targetcons,         /**< pointer to store the created target constraint */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint handler for this constraint */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           name,               /**< name of constraint, or NULL if the name of the source constraint should be used */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing? */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing? */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility? */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing? */
   SCIP_Bool             local,              /**< is constraint only valid locally? */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging? */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            success             /**< pointer to store whether the copying was successful or not */
   )
{
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;

   assert(targetcons != NULL);
   assert(sourceconshdlr != NULL);
   
   SCIP_CALL( checkStage(targetscip, "SCIPgetConsCopy", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE) );
   
   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   /* a variables map and a constraint map is needed to avoid infinite recursion */
   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNVars(sourcescip))) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNConss(sourcescip))) );
   }
   else
   { 
      localconsmap = consmap;

      /* try to retrieve copied constraint from hash map */
      *targetcons = (SCIP_CONS*) SCIPhashmapGetImage(localconsmap, sourcecons);
      if( *targetcons != NULL )
      {
         *success =TRUE;
         return SCIP_OKAY;
      }
   }

   /*  copy the constraint */
   SCIP_CALL( SCIPconsCopy(targetcons, targetscip->set, name, sourcescip, sourceconshdlr, sourcecons, localvarmap, localconsmap,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, global, success) );

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}

/** copies constraints from the source-SCIP and adds these to the target-SCIP; for mapping the
 *  variables between the source and the target SCIP a hash map can be given; if the variable hash
 *  map is NULL or necessary variable mapping is missing, the required variables are created in the
 *  target-SCIP and added to the hash map, if not NULL; all variables which are created are added to
 *  the target-SCIP but not (user) captured; if the constraint hash map is not NULL the mapping
 *  between the constraints of the source and target-SCIP is stored
 *
 *  @note the constraints are added to the target-SCIP but are not (user) captured in the target SCIP
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPcopyConss(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? 
                                              *   If TRUE, the modifiable flag of constraints will be copied. */
   SCIP_Bool*            valid               /**< pointer to store whether all constraints were validly copied */
   )
{
   SCIP_CONSHDLR** sourceconshdlrs;
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   int nsourceconshdlrs;   
   int i;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(valid      != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPcopyConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPcopyConss", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check if we locally need to create a variable or constraint hash map */
   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);
   
   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNVars(sourcescip))) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNConss(sourcescip))) );
   }
   else
      localconsmap = consmap;

   nsourceconshdlrs = SCIPgetNConshdlrs(sourcescip);
   sourceconshdlrs = SCIPgetConshdlrs(sourcescip);
   assert(nsourceconshdlrs == 0 || sourceconshdlrs != NULL);
   assert(SCIPisTransformed(sourcescip));

   *valid = TRUE;

   /* copy constraints: loop through all (source) constraint handlers */  
   for( i = 0; i < nsourceconshdlrs; ++i )
   {
      SCIP_CONS** sourceconss;
      SCIP_CONS* targetcons;
      SCIP_Bool succeed;      
      int nsourceconss;
      int c;         
   
      assert(sourceconshdlrs[i] != NULL);

      /* constraint handlers have to explicitly set the succeed pointer to TRUE */
      succeed = FALSE; 
      
      /* Get all active constraints for copying; this array contains all active constraints;
       * constraints are active if they are globally valid and not deleted after presolving OR they
       * were locally added during the search and we are currently in a node which belongs to the
       * corresponding subtree.
       */
      nsourceconss = SCIPconshdlrGetNActiveConss(sourceconshdlrs[i]);
      sourceconss = SCIPconshdlrGetConss(sourceconshdlrs[i]);

#if 0
      /* @todo using the following might reduce the number of copied constraints - check whether this is better */
      /* Get all checked constraints for copying; this included local constraints */
      if( !global )
      {
         nsourceconss = SCIPconshdlrGetNCheckConss(sourceconshdlrs[i]);
         sourceconss = SCIPconshdlrGetCheckConss(sourceconshdlrs[i]);
      }
#endif

      assert(nsourceconss == 0 || sourceconss != NULL);

      if( nsourceconss > 0 )
      {
         SCIPdebugMessage("Attempting to copy %d %s constraints\n", nsourceconss, SCIPconshdlrGetName(sourceconshdlrs[i]));
      }

      /* copy all constraints of one constraint handler */  
      for( c = 0; c < nsourceconss; ++c )
      {
         /* all constraints have to be active */
         assert(sourceconss[c] != NULL);
         assert(SCIPconsIsActive(sourceconss[c]));
         assert(!SCIPconsIsDeleted(sourceconss[c]));

         /* in case of copying the global problem we have to ignore the local constraints which are active */
         if( global && SCIPconsIsLocal(sourceconss[c]) )
         {
            SCIPdebugMessage("did not copy active local constraint <%s> when creating global copy\n", SCIPconsGetName(sourceconss[c]));
            continue;
         }

         /* use the copy constructor of the constraint handler and creates and captures the constraint if possible */
         targetcons = NULL;
         SCIP_CALL( SCIPgetConsCopy(sourcescip, targetscip, sourceconss[c], &targetcons, sourceconshdlrs[i], localvarmap, localconsmap, NULL,
               SCIPconsIsInitial(sourceconss[c]), SCIPconsIsSeparated(sourceconss[c]),
               SCIPconsIsEnforced(sourceconss[c]), SCIPconsIsChecked(sourceconss[c]),
               SCIPconsIsPropagated(sourceconss[c]), FALSE, SCIPconsIsModifiable(sourceconss[c]), 
               SCIPconsIsDynamic(sourceconss[c]), SCIPconsIsRemovable(sourceconss[c]), FALSE, global, &succeed) );

         /* add the copied constraint to target SCIP if the copying process was valid */
         if( succeed )
         {
            assert(targetcons != NULL);

            if( !enablepricing )
               SCIPconsSetModifiable(targetcons, FALSE);

            /* add constraint to target SCIP */
            SCIP_CALL( SCIPaddCons(targetscip, targetcons) );

            /* insert constraint into mapping between source SCIP and the target SCIP */
            SCIP_CALL( SCIPhashmapInsert(localconsmap, sourceconss[c], targetcons) );

            /* release constraint once for the creation capture */
            SCIP_CALL( SCIPreleaseCons(targetscip, &targetcons) );
         }
         else
         {
            *valid = FALSE;
            SCIPdebugMessage("failed to copy constraint %s\n", SCIPconsGetName(sourceconss[c]));
         }
      }
   }

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }
   
   return SCIP_OKAY;
}


/** convert all active cuts from cutpool of sourcescip to linear constraints in targetscip, sourcescip and targetscip
 *  could be the same
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPconvertCutsToConss(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   int*                  ncutsadded          /**< pointer to store number of added cuts */
   )
{
   SCIP_CUT** cuts;
   int ncuts;
   int c;

   assert(sourcescip != NULL);
   assert(sourcescip->set != NULL);
   assert(targetscip != NULL);
   assert(ncutsadded != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPconvertCutsToConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPconvertCutsToConss", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE) );

   /* if we do not have any cuts, nothing can be converted */
   if( sourcescip->set->stage < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   if( SCIPfindConshdlr(targetscip, "linear") == NULL )
   {
      SCIPdebugMessage("No linear constraint handler available. Cannot convert cuts.\n");
      return SCIP_OKAY;
   }

   cuts = SCIPgetPoolCuts(sourcescip);
   ncuts = SCIPgetNPoolCuts(sourcescip);

   for( c = 0; c < ncuts; ++c )
   {
      SCIP_ROW* row;

      row = SCIPcutGetRow(cuts[c]);
      assert(!SCIProwIsLocal(row));
      assert(!SCIProwIsModifiable(row));

      /* create a linear constraint out of the cut */
      if( SCIPcutGetAge(cuts[c]) == 0 && SCIProwIsInLP(row) )
      {
         char name[SCIP_MAXSTRLEN];
         SCIP_CONS* cons;
         SCIP_COL** cols;
         SCIP_VAR** vars;
         int ncols;
         int i;
         
         cols = SCIProwGetCols(row);
         ncols = SCIProwGetNNonz(row);
         
         /* get all variables of the row */
         SCIP_CALL( SCIPallocBufferArray(targetscip, &vars, ncols) );
         for( i = 0; i < ncols; ++i )
            vars[i] = SCIPcolGetVar(cols[i]);

         /* get corresponding variables in targetscip if necessary */
         if( sourcescip != targetscip )
         {
            SCIP_Bool success;

            for( i = 0; i < ncols; ++i )
            {
               SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, vars[i], &vars[i], varmap, consmap, global, &success) );

               if( !success )
               {
                  SCIPdebugMessage("Converting cuts to constraints failed.\n");
                  
                  /* free temporary memory */
                  SCIPfreeBufferArray(targetscip, &vars);
                  return SCIP_OKAY;
               }
            }
         }
         
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%d", SCIProwGetName(row), SCIPgetNRuns(sourcescip));
         SCIP_CALL( SCIPcreateConsLinear(targetscip, &cons, name, ncols, vars, SCIProwGetVals(row),
               SCIProwGetLhs(row) - SCIProwGetConstant(row), SCIProwGetRhs(row) - SCIProwGetConstant(row),
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE) );
         SCIP_CALL( SCIPaddCons(targetscip, cons) );

         SCIPdebugMessage("Converted cut <%s> to constraint <%s>.\n", SCIProwGetName(row), SCIPconsGetName(cons));
         SCIPdebug( SCIP_CALL( SCIPprintCons(targetscip, cons, NULL) ) );
         SCIP_CALL( SCIPreleaseCons(targetscip, &cons) );

         /* free temporary memory */
         SCIPfreeBufferArray(targetscip, &vars);

         ++(*ncutsadded);
      }
   }

   return SCIP_OKAY;
}

/** copies all active cuts from cutpool of sourcescip to constraints in targetscip 
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPcopyCuts(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global              /**< create a global or a local copy? */
   )
{
   int ncutsadded;

   assert(sourcescip != NULL);
   assert(targetscip != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPcopyCuts", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPcopyCuts", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   ncutsadded = 0;

   /* create out of all active cuts in cutpool linear constraints in targetscip */
   SCIP_CALL( SCIPconvertCutsToConss(sourcescip, targetscip, varmap, consmap, global, &ncutsadded) );

   SCIPdebugMessage("Converted %d active cuts to constraints.\n", ncutsadded);

   return SCIP_OKAY;
}

/** copies parameter settings from sourcescip to targetscip 
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPcopyParamSettings(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip          /**< target SCIP data structure */
   )
{
   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(sourcescip->set != NULL);
   assert(targetscip->set != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPcopyParamSettings", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPcopyParamSettings", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsetCopyParams(sourcescip->set, targetscip->set) );

   return SCIP_OKAY;
}

/** gets depth of current scip instance (increased by each copy call) */
int SCIPgetSubscipDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );
   assert( scip->stat != NULL );

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSubscipDepath", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->subscipdepth;
}

/** copies source SCIP to target SCIP; the copying process is done in the following order:
 *  1) copy the plugins
 *  2) copy the settings
 *  3) create problem data in target-SCIP and copy the problem data of the source-SCIP
 *  4) copy all active variables
 *  5) copy all constraints
 *
 *  @note all variables and constraints which are created in the target-SCIP are not (user) captured 
 *
 *  @note In a multi thread case, you need to lock the copying procedure from outside with a mutex.
 *  @note Do not change the source SCIP environment during the copying process
 */
SCIP_RETCODE SCIPcopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables corresponding
                                              *   target variables, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   const char*           suffix,             /**< suffix which will be added to the names of the target SCIP, might be empty */          
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool             enablepricing,      /**< should pricing be enabled in copied SCIP instance? If TRUE, pricer
                                              *   plugins will be copied and activated, and the modifiable flag of
                                              *   constraints will be respected. If FALSE, valid will be set to FALSE, when
                                              *   there are pricers present */
   SCIP_Bool*            valid               /**< pointer to store whether the copying was valid or not */
   )
{
   SCIP_HASHMAP* localvarmap;
   SCIP_HASHMAP* localconsmap;
   SCIP_Bool uselocalvarmap;
   SCIP_Bool uselocalconsmap;
   SCIP_Bool consscopyvalid;
   char name[SCIP_MAXSTRLEN];

   assert(sourcescip != NULL);
   assert(targetscip != NULL);
   assert(suffix != NULL);
   assert(valid != NULL);

   /* check stages for both, the source and the target SCIP data structure */
   SCIP_CALL( checkStage(sourcescip, "SCIPcopy", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
   SCIP_CALL( checkStage(targetscip, "SCIPcopy", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* in case there are active pricers and pricing is disabled the target SCIP will not be a valid copy of the source
    * SCIP 
    */
   (*valid) = enablepricing || (SCIPgetNActivePricers(sourcescip) == 0);

   /* copy all plugins */
   SCIP_CALL( SCIPcopyPlugins(sourcescip, targetscip, TRUE, enablepricing, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
         TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, valid) );

   SCIPdebugMessage("Copying plugins was%s valid.\n", *valid ? "" : " not");

   uselocalvarmap = (varmap == NULL);
   uselocalconsmap = (consmap == NULL);

   if( uselocalvarmap )
   {
      /* create the variable mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localvarmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNVars(sourcescip))) );
   }
   else
      localvarmap = varmap;

   if( uselocalconsmap )
   {
      /* create the constraint mapping hash map */
      SCIP_CALL( SCIPhashmapCreate(&localconsmap, SCIPblkmem(targetscip), SCIPcalcHashtableSize(HASHTABLESIZE_FACTOR * SCIPgetNConss(sourcescip))) );
   }
   else
      localconsmap = consmap;

   /* construct name for the target SCIP using the source problem name and the given suffix string */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%s", SCIPgetProbName(sourcescip), suffix);

   /* copy all settings */
   SCIP_CALL( SCIPcopyParamSettings(sourcescip, targetscip) );

   /* create problem in the target SCIP and copying the source problem data */
   SCIP_CALL( SCIPcopyProb(sourcescip, targetscip, localvarmap, localconsmap, global, name) );

   /* copy all active variables */
   SCIP_CALL( SCIPcopyVars(sourcescip, targetscip, localvarmap, localconsmap, global) );

   /* copy all constraints */
   SCIP_CALL( SCIPcopyConss(sourcescip, targetscip, localvarmap, localconsmap, global, enablepricing, &consscopyvalid) );
   (*valid) = *valid && consscopyvalid;

   SCIPdebugMessage("Copying constraints was%s valid.\n", consscopyvalid ? "" : " not");

   if( uselocalvarmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localvarmap);
   }

   if( uselocalconsmap )
   {
      /* free hash map */
      SCIPhashmapFree(&localconsmap);
   }

   return SCIP_OKAY;
}



/*
 * parameter settings
 */

/** creates a SCIP_Bool parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPaddBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Bool*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Bool             defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetAddBoolParam(scip->set, scip->mem->setmem, name, desc, valueptr, isadvanced, defaultvalue,
         paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a int parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPaddIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   int*                  valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   int                   defaultvalue,       /**< default value of the parameter */
   int                   minvalue,           /**< minimum value for parameter */
   int                   maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetAddIntParam(scip->set, scip->mem->setmem, name, desc, valueptr, isadvanced, defaultvalue, minvalue,
         maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a SCIP_Longint parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPaddLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Longint*         valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Longint          defaultvalue,       /**< default value of the parameter */
   SCIP_Longint          minvalue,           /**< minimum value for parameter */
   SCIP_Longint          maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetAddLongintParam(scip->set, scip->mem->setmem, name, desc,
         valueptr, isadvanced, defaultvalue, minvalue, maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a SCIP_Real parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPaddRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   SCIP_Real*            valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   SCIP_Real             defaultvalue,       /**< default value of the parameter */
   SCIP_Real             minvalue,           /**< minimum value for parameter */
   SCIP_Real             maxvalue,           /**< maximum value for parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetAddRealParam(scip->set, scip->mem->setmem, name, desc, valueptr, isadvanced, defaultvalue, minvalue,
         maxvalue, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a char parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPaddCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char*                 valueptr,           /**< pointer to store the current parameter value, or NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   char                  defaultvalue,       /**< default value of the parameter */
   const char*           allowedvalues,      /**< array with possible parameter values, or NULL if not restricted */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetAddCharParam(scip->set, scip->mem->setmem, name, desc, valueptr, isadvanced, defaultvalue,
         allowedvalues, paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** creates a string parameter, sets it to its default value, and adds it to the parameter set */
SCIP_RETCODE SCIPaddStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           desc,               /**< description of the parameter */
   char**                valueptr,           /**< pointer to store the current parameter value, or NULL; if not NULL then *valueptr should be NULL */
   SCIP_Bool             isadvanced,         /**< is this parameter an advanced parameter? */
   const char*           defaultvalue,       /**< default value of the parameter */
   SCIP_DECL_PARAMCHGD   ((*paramchgd)),     /**< change information method of parameter */
   SCIP_PARAMDATA*       paramdata           /**< locally defined parameter specific data */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetAddStringParam(scip->set, scip->mem->setmem, name, desc, valueptr, isadvanced, defaultvalue,
         paramchgd, paramdata) );

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPgetBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool*            value               /**< pointer to store the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetGetBoolParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Int parameter */
SCIP_RETCODE SCIPgetIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   int*                  value               /**< pointer to store the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetGetIntParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPgetLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint*         value               /**< pointer to store the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetGetLongintParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPgetRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Real*            value               /**< pointer to store the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetGetRealParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing Char parameter */
SCIP_RETCODE SCIPgetCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char*                 value               /**< pointer to store the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetGetCharParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** gets the value of an existing String parameter */
SCIP_RETCODE SCIPgetStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char**                value               /**< pointer to store the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetGetStringParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing parameter */
SCIP_RETCODE SCIPsetParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   void*                 value               /**< new value of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Bool parameter */
SCIP_RETCODE SCIPsetBoolParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Bool             value               /**< new value of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetBoolParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetBoolParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Int parameter */
SCIP_RETCODE SCIPsetIntParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   int                   value               /**< new value of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetIntParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetIntParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Longint parameter */
SCIP_RETCODE SCIPsetLongintParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Longint          value               /**< new value of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetLongintParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetLongintParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing SCIP_Real parameter */
SCIP_RETCODE SCIPsetRealParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   SCIP_Real             value               /**< new value of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetRealParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetRealParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing Char parameter */
SCIP_RETCODE SCIPsetCharParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   char                  value               /**< new value of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetCharParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetCharParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** changes the value of an existing String parameter */
SCIP_RETCODE SCIPsetStringParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of the parameter */
   const char*           value               /**< new value of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetStringParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetStringParam(scip->set, name, value) );

   return SCIP_OKAY;
}

/** reads parameters from a file */
SCIP_RETCODE SCIPreadParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPreadParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetReadParams(scip->set, filename) );

   return SCIP_OKAY;
}

/** writes all parameters in the parameter set to a file */
SCIP_RETCODE SCIPwriteParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< file name, or NULL for stdout */
   SCIP_Bool             comments,           /**< should parameter descriptions be written as comments? */
   SCIP_Bool             onlychanged         /**< should only the parameters been written, that are changed from default? */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPwriteParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetWriteParams(scip->set, filename, comments, onlychanged) );

   return SCIP_OKAY;
}

/** resets a single parameter to its default value */
SCIP_RETCODE SCIPresetParam(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the parameter */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPresetParam", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetResetParam(scip->set, name) );

   return SCIP_OKAY;
}

/** resets all parameters to their default values */
SCIP_RETCODE SCIPresetParams(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPresetParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetResetParams(scip->set) );

   return SCIP_OKAY;
}

/** sets parameters to 
 *  - SCIP_PARAMEMPHASIS_DEFAULT to use default values (see also SCIPresetParams())
 *  - SCIP_PARAMEMPHASIS_COUNTER to get feasible and "fast" counting process
 *  - SCIP_PARAMEMPHASIS_CPSOLVER to get CP like search (e.g. no LP relaxation)
 *  - SCIP_PARAMEMPHASIS_EASYCIP to solve easy problems fast
 *  - SCIP_PARAMEMPHASIS_FEASIBILITY to detect feasibility fast 
 *  - SCIP_PARAMEMPHASIS_HARDLP to be capable to handle hard LPs
 *  - SCIP_PARAMEMPHASIS_OPTIMALITY to prove optimality fast
 */
SCIP_RETCODE SCIPsetEmphasis(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMEMPHASIS    paramemphasis,      /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetEmphasis", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetSetEmphasis(scip->set, paramemphasis, quiet) );
   
   return SCIP_OKAY;
}

/** sets parameters to deactivate separators and heuristics that use auxiliary SCIP instances; should be called for
 *  auxiliary SCIP instances to avoid recursion
 */
SCIP_RETCODE SCIPsetSubscipsOff(
   SCIP*                 scip,               /**< (auxiliary) SCIP data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetSubscipsOff", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   SCIP_CALL( SCIPsetSetSubscipsOff(scip->set, quiet) );
   
   return SCIP_OKAY;
}

/** sets heuristic parameters values to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all heuristic parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for heuristic is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the heuristic are called more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all heuristics
 */
SCIP_RETCODE SCIPsetHeuristics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetHeuristics", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);

   SCIP_CALL( SCIPsetSetHeuristics(scip->set, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** sets presolving parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all presolving parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for presolving is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the presolving is more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all presolving
 */
SCIP_RETCODE SCIPsetPresolving(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetPresolving", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);
   
   SCIP_CALL( SCIPsetSetPresolving(scip->set, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** sets separating parameters to 
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all separating parameters 
 *  - SCIP_PARAMSETTING_FAST such that the time spend for separating is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the separating is done more aggregative
 *  - SCIP_PARAMSETTING_OFF which turn off all separating
 */
SCIP_RETCODE SCIPsetSeparating(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetSeparating", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   assert(paramsetting == SCIP_PARAMSETTING_DEFAULT || paramsetting == SCIP_PARAMSETTING_FAST
      || paramsetting == SCIP_PARAMSETTING_AGGRESSIVE || paramsetting == SCIP_PARAMSETTING_OFF);
   
   SCIP_CALL( SCIPsetSetSeparating(scip->set, paramsetting, quiet) );

   return SCIP_OKAY;
}

/** returns the array of all available SCIP parameters */
SCIP_PARAM** SCIPgetParams(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetGetParams(scip->set);
}

/** returns the total number of all available SCIP parameters */
int SCIPgetNParams(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNParams", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetGetNParams(scip->set);
}




/*
 * SCIP user functionality methods: managing plugins
 */

/** creates a reader and includes it in SCIP */
SCIP_RETCODE SCIPincludeReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of reader */
   const char*           desc,               /**< description of reader */
   const char*           extension,          /**< file extension that reader processes */
   SCIP_DECL_READERCOPY  ((*readercopy)),    /**< copy method of reader or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_READERFREE  ((*readerfree)),    /**< destructor of reader */
   SCIP_DECL_READERREAD  ((*readerread)),    /**< read method */
   SCIP_DECL_READERWRITE ((*readerwrite)),   /**< write method */
   SCIP_READERDATA*      readerdata          /**< reader data */
   )
{
   SCIP_READER* reader;

   SCIP_CALL( checkStage(scip, "SCIPincludeReader", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether reader is already present */
   if( SCIPfindReader(scip, name) != NULL )
   {
      SCIPerrorMessage("reader <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPreaderCreate(&reader, name, desc, extension, readercopy, readerfree, readerread, readerwrite, readerdata) );
   SCIP_CALL( SCIPsetIncludeReader(scip->set, reader) );

   return SCIP_OKAY;
}

/** returns the reader of the given name, or NULL if not existing */
SCIP_READER* SCIPfindReader(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint handler */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindReader", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindReader(scip->set, name);
}

/** returns the array of currently available readers */
SCIP_READER** SCIPgetReaders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetReaders", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->readers;
}

/** returns the number of currently available readers */
int SCIPgetNReaders(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNReaders", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nreaders;
}

/** creates a variable pricer and includes it in SCIP
 *  To use the variable pricer for solving a problem, it first has to be activated with a call to SCIPactivatePricer().
 *  This should be done during the problem creation stage.
 */
SCIP_RETCODE SCIPincludePricer(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of variable pricer */
   const char*           desc,               /**< description of variable pricer */
   int                   priority,           /**< priority of the variable pricer */
   SCIP_Bool             delay,              /**< should the pricer be delayed until no other pricers or already existing
                                              *   problem variables with negative reduced costs are found?
                                              *   if this is set to FALSE it may happen that the pricer produces columns
                                              *   that already exist in the problem (which are also priced in by the
                                              *   default problem variable pricing in the same round)
                                              */
   SCIP_DECL_PRICERCOPY  ((*pricercopy)),    /**< copy method of variable pricer or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRICERFREE  ((*pricerfree)),    /**< destructor of variable pricer */
   SCIP_DECL_PRICERINIT  ((*pricerinit)),    /**< initialize variable pricer */
   SCIP_DECL_PRICEREXIT  ((*pricerexit)),    /**< deinitialize variable pricer */
   SCIP_DECL_PRICERINITSOL((*pricerinitsol)),/**< solving process initialization method of variable pricer */
   SCIP_DECL_PRICEREXITSOL((*pricerexitsol)),/**< solving process deinitialization method of variable pricer */
   SCIP_DECL_PRICERREDCOST((*pricerredcost)),/**< reduced cost pricing method of variable pricer for feasible LPs */
   SCIP_DECL_PRICERFARKAS((*pricerfarkas)),  /**< Farkas pricing method of variable pricer for infeasible LPs */
   SCIP_PRICERDATA*      pricerdata          /**< variable pricer data */
   )
{
   SCIP_PRICER* pricer;

   SCIP_CALL( checkStage(scip, "SCIPincludePricer", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether pricer is already present */
   if( SCIPfindPricer(scip, name) != NULL )
   {
      SCIPerrorMessage("pricer <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpricerCreate(&pricer, scip->set, scip->mem->setmem,
         name, desc, priority, delay,
         pricercopy,
         pricerfree, pricerinit, pricerexit, pricerinitsol, pricerexitsol, pricerredcost, pricerfarkas, pricerdata) );
   SCIP_CALL( SCIPsetIncludePricer(scip->set, pricer) );

   return SCIP_OKAY;
}

/** returns the variable pricer of the given name, or NULL if not existing */
SCIP_PRICER* SCIPfindPricer(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable pricer */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindPricer", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindPricer(scip->set, name);
}

/** returns the array of currently available variable pricers; active pricers are in the first slots of the array */
SCIP_PRICER** SCIPgetPricers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortPricers(scip->set);

   return scip->set->pricers;
}

/** returns the number of currently available variable pricers */
int SCIPgetNPricers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->npricers;
}

/** returns the number of currently active variable pricers, that are used in the LP solving loop */
int SCIPgetNActivePricers(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNAcvitePricers", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nactivepricers;
}

/** sets the priority of a variable pricer */
SCIP_RETCODE SCIPsetPricerPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer,             /**< variable pricer */
   int                   priority            /**< new priority of the variable pricer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetPricerPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpricerSetPriority(pricer, scip->set, priority);

   return SCIP_OKAY;
}

/** activates pricer to be used for the current problem
 *  This method should be called during the problem creation stage for all pricers that are necessary to solve
 *  the problem model.
 *  The pricers are automatically deactivated when the problem is freed.
 */
SCIP_RETCODE SCIPactivatePricer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPactivatePricer", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPpricerActivate(pricer, scip->set) );

   return SCIP_OKAY;
}

/** deactivates pricer */
SCIP_RETCODE SCIPdeactivatePricer(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRICER*          pricer              /**< variable pricer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdeactivatePricer", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE) );

   SCIP_CALL( SCIPpricerDeactivate(pricer, scip->set) );

   return SCIP_OKAY;
}

/** creates a constraint handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of constraint handler */
   const char*           desc,               /**< description of constraint handler */
   int                   sepapriority,       /**< priority of the constraint handler for separation */
   int                   enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int                   chckpriority,       /**< priority of the constraint handler for checking feasibility (and propagation) */
   int                   sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int                   propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   int                   eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
   int                   maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   SCIP_Bool             delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
   SCIP_Bool             delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
   SCIP_Bool             delaypresol,        /**< should presolving method be delayed, if other presolvers found reductions? */
   SCIP_Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagators should be executed */
   SCIP_DECL_CONSHDLRCOPY((*conshdlrcopy)),  /**< copy method of constraint handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   SCIP_DECL_CONSINIT    ((*consinit)),      /**< initialize constraint handler */
   SCIP_DECL_CONSEXIT    ((*consexit)),      /**< deinitialize constraint handler */
   SCIP_DECL_CONSINITPRE ((*consinitpre)),   /**< presolving initialization method of constraint handler */
   SCIP_DECL_CONSEXITPRE ((*consexitpre)),   /**< presolving deinitialization method of constraint handler */
   SCIP_DECL_CONSINITSOL ((*consinitsol)),   /**< solving process initialization method of constraint handler */
   SCIP_DECL_CONSEXITSOL ((*consexitsol)),   /**< solving process deinitialization method of constraint handler */
   SCIP_DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   SCIP_DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   SCIP_DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   SCIP_DECL_CONSSEPALP  ((*conssepalp)),    /**< separate cutting planes for LP solution */
   SCIP_DECL_CONSSEPASOL ((*conssepasol)),   /**< separate cutting planes for arbitrary primal solution */
   SCIP_DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   SCIP_DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   SCIP_DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   SCIP_DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   SCIP_DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   SCIP_DECL_CONSRESPROP ((*consresprop)),   /**< propagation conflict resolving method */
   SCIP_DECL_CONSLOCK    ((*conslock)),      /**< variable rounding lock method */
   SCIP_DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   SCIP_DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
   SCIP_DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   SCIP_DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
   SCIP_DECL_CONSDELVARS ((*consdelvars)),   /**< variable deletion method */
   SCIP_DECL_CONSPRINT   ((*consprint)),     /**< constraint display method */
   SCIP_DECL_CONSCOPY    ((*conscopy)),      /**< constraint copying method */
   SCIP_DECL_CONSPARSE   ((*consparse)),     /**< constraint parsing method */
   SCIP_CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   )
{
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( checkStage(scip, "SCIPincludeConshdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether constraint handler is already present */
   if( SCIPfindConshdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("constraint handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconshdlrCreate(&conshdlr, scip->set, scip->mem->setmem,
         name, desc, sepapriority, enfopriority, chckpriority, sepafreq, propfreq, eagerfreq, maxprerounds,
         delaysepa, delayprop, delaypresol, needscons,
         timingmask,
         conshdlrcopy,
         consfree, consinit, consexit, consinitpre, consexitpre, consinitsol, consexitsol,
         consdelete, constrans, consinitlp, conssepalp, conssepasol, consenfolp, consenfops, conscheck, consprop,
         conspresol, consresprop, conslock, consactive, consdeactive, consenable, consdisable, consdelvars, consprint,
         conscopy, consparse, conshdlrdata) );
   SCIP_CALL( SCIPsetIncludeConshdlr(scip->set, conshdlr) );

   return SCIP_OKAY;
}

/** returns the constraint handler of the given name, or NULL if not existing */
SCIP_CONSHDLR* SCIPfindConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint handler */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindConshdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindConshdlr(scip->set, name);
}

/** returns the array of currently available constraint handlers */
SCIP_CONSHDLR** SCIPgetConshdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetConshdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->conshdlrs;
}

/** returns the number of currently available constraint handlers */
int SCIPgetNConshdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNConshdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nconshdlrs;
}

/** creates a conflict handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeConflicthdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTCOPY((*conflictcopy)),  /**< copy method of conflict handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   )
{
   SCIP_CONFLICTHDLR* conflicthdlr;

   SCIP_CALL( checkStage(scip, "SCIPincludeConflicthdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether conflict handler is already present */
   if( SCIPfindConflicthdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("conflict handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPconflicthdlrCreate(&conflicthdlr, scip->set, scip->mem->setmem, name, desc, priority,
         conflictcopy,
         conflictfree, conflictinit, conflictexit, conflictinitsol, conflictexitsol, conflictexec,
         conflicthdlrdata) );
   SCIP_CALL( SCIPsetIncludeConflicthdlr(scip->set, conflicthdlr) );

   return SCIP_OKAY;
}

/** returns the conflict handler of the given name, or NULL if not existing */
SCIP_CONFLICTHDLR* SCIPfindConflicthdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of conflict handler */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindConflicthdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindConflicthdlr(scip->set, name);
}

/** returns the array of currently available conflict handlers */
SCIP_CONFLICTHDLR** SCIPgetConflicthdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetConflicthdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortConflicthdlrs(scip->set);

   return scip->set->conflicthdlrs;
}

/** returns the number of currently available conflict handlers */
int SCIPgetNConflicthdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNConflicthdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nconflicthdlrs;
}

/** sets the priority of a conflict handler */
SCIP_RETCODE SCIPsetConflicthdlrPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   int                   priority            /**< new priority of the conflict handler */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConflicthdlrPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPconflicthdlrSetPriority(conflicthdlr, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of presolver */
   const char*           desc,               /**< description of presolver */
   int                   priority,           /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
   int                   maxrounds,          /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
   SCIP_Bool             delay,              /**< should presolver be delayed, if other presolvers found reductions? */
   SCIP_DECL_PRESOLCOPY  ((*presolcopy)),    /**< copy method of presolver or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PRESOLFREE  ((*presolfree)),    /**< destructor of presolver to free user data (called when SCIP is exiting) */
   SCIP_DECL_PRESOLINIT  ((*presolinit)),    /**< initialization method of presolver (called after problem was transformed) */
   SCIP_DECL_PRESOLEXIT  ((*presolexit)),    /**< deinitialization method of presolver (called before transformed problem is freed) */
   SCIP_DECL_PRESOLINITPRE((*presolinitpre)),/**< presolving initialization method of presolver (called when presolving is about to begin) */
   SCIP_DECL_PRESOLEXITPRE((*presolexitpre)),/**< presolving deinitialization method of presolver (called after presolving has been finished) */
   SCIP_DECL_PRESOLEXEC  ((*presolexec)),    /**< execution method of presolver */
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   SCIP_PRESOL* presol;

   SCIP_CALL( checkStage(scip, "SCIPincludePresol", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether presolver is already present */
   if( SCIPfindPresol(scip, name) != NULL )
   {
      SCIPerrorMessage("presolver <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpresolCreate(&presol, scip->set, scip->mem->setmem, name, desc, priority, maxrounds, delay,
         presolcopy,
         presolfree, presolinit, presolexit, presolinitpre, presolexitpre, presolexec, presoldata) );
   SCIP_CALL( SCIPsetIncludePresol(scip->set, presol) );

   return SCIP_OKAY;
}

/** returns the presolver of the given name, or NULL if not existing */
SCIP_PRESOL* SCIPfindPresol(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of presolver */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindPresol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindPresol(scip->set, name);
}

/** returns the array of currently available presolvers */
SCIP_PRESOL** SCIPgetPresols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPresols", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortPresols(scip->set);

   return scip->set->presols;
}

/** returns the number of currently available presolvers */
int SCIPgetNPresols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPresols", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->npresols;
}

/** sets the priority of a presolver */
SCIP_RETCODE SCIPsetPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOL*          presol,             /**< presolver */
   int                   priority            /**< new priority of the presolver */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetPresolPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpresolSetPriority(presol, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a relaxation handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of relaxation handler */
   const char*           desc,               /**< description of relaxation handler */
   int                   priority,           /**< priority of the relaxation handler (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxation handler */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy)),     /**< copy method of relaxation handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxation handler */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxation handler */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxation handler */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxation handler */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxation handler */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxation handler */
   SCIP_RELAXDATA*       relaxdata           /**< relaxation handler data */
   )
{
   SCIP_RELAX* relax;

   SCIP_CALL( checkStage(scip, "SCIPincludeRelax", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether relaxation handler is already present */
   if( SCIPfindRelax(scip, name) != NULL )
   {
      SCIPerrorMessage("relaxation handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPrelaxCreate(&relax, scip->set, scip->mem->setmem,
         name, desc, priority, freq,
         relaxcopy,
         relaxfree, relaxinit, relaxexit, relaxinitsol, relaxexitsol, relaxexec, relaxdata) );
   SCIP_CALL( SCIPsetIncludeRelax(scip->set, relax) );

   return SCIP_OKAY;
}

/** returns the relaxation handler of the given name, or NULL if not existing */
SCIP_RELAX* SCIPfindRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of relaxation handler */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindRelax", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindRelax(scip->set, name);
}

/** returns the array of currently available relaxation handlers  */
SCIP_RELAX** SCIPgetRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRelaxs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortRelaxs(scip->set);

   return scip->set->relaxs;
}

/** returns the number of currently available relaxation handlers  */
int SCIPgetNRelaxs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNRelaxs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nrelaxs;
}

/** sets the priority of a relaxation handler */
SCIP_RETCODE SCIPsetRelaxPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax,              /**< relaxation handler */
   int                   priority            /**< new priority of the relaxation handler */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetRelaxPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPrelaxSetPriority(relax, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of separator */
   const char*           desc,               /**< description of separator */
   int                   priority,           /**< priority of separator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling separator */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound compared
                                              *   to best node's dual bound for applying separation */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_Bool             delay,              /**< should separator be delayed, if other separators found cuts? */
   SCIP_DECL_SEPACOPY    ((*sepacopy)),      /**< copy method of separator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_SEPAFREE    ((*sepafree)),      /**< destructor of separator */
   SCIP_DECL_SEPAINIT    ((*sepainit)),      /**< initialize separator */
   SCIP_DECL_SEPAEXIT    ((*sepaexit)),      /**< deinitialize separator */
   SCIP_DECL_SEPAINITSOL ((*sepainitsol)),   /**< solving process initialization method of separator */
   SCIP_DECL_SEPAEXITSOL ((*sepaexitsol)),   /**< solving process deinitialization method of separator */
   SCIP_DECL_SEPAEXECLP  ((*sepaexeclp)),    /**< LP solution separation method of separator */
   SCIP_DECL_SEPAEXECSOL ((*sepaexecsol)),   /**< arbitrary primal solution separation method of separator */
   SCIP_SEPADATA*        sepadata            /**< separator data */
   )
{
   SCIP_SEPA* sepa;

   SCIP_CALL( checkStage(scip, "SCIPincludeSepa", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether separator is already present */
   if( SCIPfindSepa(scip, name) != NULL )
   {
      SCIPerrorMessage("separator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPsepaCreate(&sepa, scip->set, scip->mem->setmem,
         name, desc, priority, freq, maxbounddist, usessubscip, delay,
         sepacopy,
         sepafree, sepainit, sepaexit, sepainitsol, sepaexitsol, sepaexeclp, sepaexecsol, sepadata) );
   SCIP_CALL( SCIPsetIncludeSepa(scip->set, sepa) );

   return SCIP_OKAY;
}

/** returns the separator of the given name, or NULL if not existing */
SCIP_SEPA* SCIPfindSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of separator */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindSepa", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindSepa(scip->set, name);
}

/** returns the array of currently available separators */
SCIP_SEPA** SCIPgetSepas(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSepas", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortSepas(scip->set);

   return scip->set->sepas;
}

/** returns the number of currently available separators */
int SCIPgetNSepas(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNSepas", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nsepas;
}

/** sets the priority of a separator */
SCIP_RETCODE SCIPsetSepaPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   int                   priority            /**< new priority of the separator */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetSepaPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsepaSetPriority(sepa, scip->set, priority);

   return SCIP_OKAY;
}

/** creates a propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludeProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of propagator */
   const char*           desc,               /**< description of propagator */
   int                   priority,           /**< priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   freq,               /**< frequency for calling propagator */
   SCIP_Bool             delay,              /**< should propagator be delayed, if other propagators found reductions? */
   SCIP_PROPTIMING       timingmask,         /**< positions in the node solving loop where propagators should be executed */
   int                   presolpriority,     /**< presolving priority of the propagator (>= 0: before, < 0: after constraint handlers) */
   int                   presolmaxrounds,    /**< maximal number of presolving rounds the propagator participates in (-1: no limit) */
   SCIP_Bool             presoldelay,        /**< should presolving be delayed, if other presolvers found reductions? */
   SCIP_DECL_PROPCOPY    ((*propcopy)),      /**< copy method of propagator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_PROPFREE    ((*propfree)),      /**< destructor of propagator */
   SCIP_DECL_PROPINIT    ((*propinit)),      /**< initialize propagator */
   SCIP_DECL_PROPEXIT    ((*propexit)),      /**< deinitialize propagator */
   SCIP_DECL_PROPINITPRE ((*propinitpre)),   /**< presolving initialization method of propagator */
   SCIP_DECL_PROPEXITPRE ((*propexitpre)),   /**< presolving deinitialization method of propagator */
   SCIP_DECL_PROPINITSOL ((*propinitsol)),   /**< solving process initialization method of propagator */
   SCIP_DECL_PROPEXITSOL ((*propexitsol)),   /**< solving process deinitialization method of propagator */
   SCIP_DECL_PROPPRESOL  ((*proppresol)),    /**< presolving method */
   SCIP_DECL_PROPEXEC    ((*propexec)),      /**< execution method of propagator */
   SCIP_DECL_PROPRESPROP ((*propresprop)),   /**< propagation conflict resolving method */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_PROP* prop;

   SCIP_CALL( checkStage(scip, "SCIPincludeProp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether propagator is already present */
   if( SCIPfindProp(scip, name) != NULL )
   {
      SCIPerrorMessage("propagator <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPpropCreate(&prop, scip->set, scip->mem->setmem,
         name, desc, priority, freq, delay, timingmask, presolpriority, presolmaxrounds, presoldelay,
         propcopy,
         propfree, propinit, propexit, propinitpre, propexitpre, propinitsol, propexitsol,
         proppresol, propexec, propresprop, propdata) );
   SCIP_CALL( SCIPsetIncludeProp(scip->set, prop) );

   return SCIP_OKAY;
}

/** returns the propagator of the given name, or NULL if not existing */
SCIP_PROP* SCIPfindProp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of propagator */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindProp", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindProp(scip->set, name);
}

/** returns the array of currently available propagators */
SCIP_PROP** SCIPgetProps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetProps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortProps(scip->set);

   return scip->set->props;
}

/** returns the number of currently available propagators */
int SCIPgetNProps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNProps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nprops;
}

/** sets the priority of a propagator */
SCIP_RETCODE SCIPsetPropPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   priority            /**< new priority of the propagator */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetPropPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpropSetPriority(prop, scip->set, priority);

   return SCIP_OKAY;
}

/** sets the presolving priority of a propagator */
SCIP_RETCODE SCIPsetPropPresolPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   int                   presolpriority      /**< new presol priority of the propagator */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetPropPresolPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPpropSetPresolPriority(prop, scip->set, presolpriority);

   return SCIP_OKAY;
}

/** creates a primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   unsigned int          timingmask,         /**< positions in the node solving loop where heuristic should be executed;
                                              *   see definition of SCIP_HeurTiming for possible values */
   SCIP_Bool             usessubscip,        /**< does the heuristic use a secondary SCIP instance? */
   SCIP_DECL_HEURCOPY    ((*heurcopy)),      /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol)),   /**< solving process initialization method of primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol)),   /**< solving process deinitialization method of primal heuristic */
   SCIP_DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   )
{
   SCIP_HEUR* heur;

   SCIP_CALL( checkStage(scip, "SCIPincludeHeur", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether heuristic is already present */
   if( SCIPfindHeur(scip, name) != NULL )
   {
      SCIPerrorMessage("heuristic <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPheurCreate(&heur, scip->set, scip->mem->setmem,
         name, desc, dispchar, priority, freq, freqofs, maxdepth, timingmask, usessubscip, 
         heurcopy, heurfree, heurinit, heurexit, heurinitsol, heurexitsol, heurexec, heurdata) );
   SCIP_CALL( SCIPsetIncludeHeur(scip->set, heur) );

   return SCIP_OKAY;
}

/** returns the primal heuristic of the given name, or NULL if not existing */
SCIP_HEUR* SCIPfindHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of primal heuristic */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindHeur", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindHeur(scip->set, name);
}

/** returns the array of currently available primal heuristics */
SCIP_HEUR** SCIPgetHeurs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetHeurs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortHeurs(scip->set);

   return scip->set->heurs;
}

/** returns the number of currently available primal heuristics */
int SCIPgetNHeurs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNHeurs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nheurs;
}

/** sets the priority of a primal heuristic */
SCIP_RETCODE SCIPsetHeurPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   priority            /**< new priority of the primal heuristic */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetHeurPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPheurSetPriority(heur, scip->set, priority);

   return SCIP_OKAY;
}

/** creates an event handler and includes it in SCIP */
SCIP_RETCODE SCIPincludeEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of event handler */
   const char*           desc,               /**< description of event handler */
   SCIP_DECL_EVENTCOPY   ((*eventcopy)),     /**< copy method of event handler or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_EVENTFREE   ((*eventfree)),     /**< destructor of event handler */
   SCIP_DECL_EVENTINIT   ((*eventinit)),     /**< initialize event handler */
   SCIP_DECL_EVENTEXIT   ((*eventexit)),     /**< deinitialize event handler */
   SCIP_DECL_EVENTINITSOL((*eventinitsol)),  /**< solving process initialization method of event handler */
   SCIP_DECL_EVENTEXITSOL((*eventexitsol)),  /**< solving process deinitialization method of event handler */
   SCIP_DECL_EVENTDELETE ((*eventdelete)),   /**< free specific event data */
   SCIP_DECL_EVENTEXEC   ((*eventexec)),     /**< execute event handler */
   SCIP_EVENTHDLRDATA*   eventhdlrdata       /**< event handler data */
   )
{
   SCIP_EVENTHDLR* eventhdlr;

   SCIP_CALL( checkStage(scip, "SCIPincludeEventhdlr", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether event handler is already present */
   if( SCIPfindEventhdlr(scip, name) != NULL )
   {
      SCIPerrorMessage("event handler <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPeventhdlrCreate(&eventhdlr, name, desc,
         eventcopy,
         eventfree, eventinit, eventexit, eventinitsol, eventexitsol, eventdelete, eventexec,
         eventhdlrdata) );
   SCIP_CALL( SCIPsetIncludeEventhdlr(scip->set, eventhdlr) );

   return SCIP_OKAY;
}

/** returns the event handler of the given name, or NULL if not existing */
SCIP_EVENTHDLR* SCIPfindEventhdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of event handler */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindEventhdlr", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindEventhdlr(scip->set, name);
}

/** returns the array of currently available event handlers */
SCIP_EVENTHDLR** SCIPgetEventhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetEventhdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->eventhdlrs;
}

/** returns the number of currently available event handlers */
int SCIPgetNEventhdlrs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNEventhdlrs", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->neventhdlrs;
}

/** creates a node selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of node selector */
   const char*           desc,               /**< description of node selector */
   int                   stdpriority,        /**< priority of the node selector in standard mode */
   int                   memsavepriority,    /**< priority of the node selector in memory saving mode */
   SCIP_DECL_NODESELCOPY ((*nodeselcopy)),   /**< copy method of node selector or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_NODESELFREE ((*nodeselfree)),   /**< destructor of node selector */
   SCIP_DECL_NODESELINIT ((*nodeselinit)),   /**< initialize node selector */
   SCIP_DECL_NODESELEXIT ((*nodeselexit)),   /**< deinitialize node selector */
   SCIP_DECL_NODESELINITSOL((*nodeselinitsol)),/**< solving process initialization method of node selector */
   SCIP_DECL_NODESELEXITSOL((*nodeselexitsol)),/**< solving process deinitialization method of node selector */
   SCIP_DECL_NODESELSELECT((*nodeselselect)),/**< node selection method */
   SCIP_DECL_NODESELCOMP ((*nodeselcomp)),   /**< node comparison method */
   SCIP_NODESELDATA*     nodeseldata         /**< node selector data */
   )
{
   SCIP_NODESEL* nodesel;

   SCIP_CALL( checkStage(scip, "SCIPincludeNodesel", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether node selector is already present */
   if( SCIPfindNodesel(scip, name) != NULL )
   {
      SCIPerrorMessage("node selector <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPnodeselCreate(&nodesel, scip->set, scip->mem->setmem, name, desc, stdpriority, memsavepriority,
         nodeselcopy,
         nodeselfree, nodeselinit, nodeselexit, nodeselinitsol, nodeselexitsol,
         nodeselselect, nodeselcomp, nodeseldata) );
   SCIP_CALL( SCIPsetIncludeNodesel(scip->set, nodesel) );

   return SCIP_OKAY;
}

/** returns the node selector of the given name, or NULL if not existing */
SCIP_NODESEL* SCIPfindNodesel(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of node selector */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindNodesel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindNodesel(scip->set, name);
}

/** returns the array of currently available node selectors */
SCIP_NODESEL** SCIPgetNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNodesels", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nodesels;
}

/** returns the number of currently available node selectors */
int SCIPgetNNodesels(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNodesels", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nnodesels;
}

/** sets the priority of a node selector in standard mode */
SCIP_RETCODE SCIPsetNodeselStdPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new standard priority of the node selector */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNodeselStdPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPnodeselSetStdPriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** sets the priority of a node selector in memory saving mode */
SCIP_RETCODE SCIPsetNodeselMemsavePriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODESEL*         nodesel,            /**< node selector */
   int                   priority            /**< new memory saving priority of the node selector */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNodeselMemsavePriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPnodeselSetMemsavePriority(nodesel, scip->set, priority);

   return SCIP_OKAY;
}

/** returns the currently used node selector */
SCIP_NODESEL* SCIPgetNodesel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNodesel", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetGetNodesel(scip->set, scip->stat);
}

/** creates a branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of branching rule */
   const char*           desc,               /**< description of branching rule */
   int                   priority,           /**< priority of the branching rule */
   int                   maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   SCIP_Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                              *   compared to best node's dual bound for applying branching rule
                                              *   (0.0: only on current best node, 1.0: on all nodes) */
   SCIP_DECL_BRANCHCOPY  ((*branchcopy)),    /**< copy method of branching rule or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   SCIP_DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   SCIP_DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   SCIP_DECL_BRANCHINITSOL((*branchinitsol)),/**< solving process initialization method of branching rule */
   SCIP_DECL_BRANCHEXITSOL((*branchexitsol)),/**< solving process deinitialization method of branching rule */
   SCIP_DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   SCIP_DECL_BRANCHEXECEXT((*branchexecext)),/**< branching execution method for external candidates */
   SCIP_DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   SCIP_BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   SCIP_BRANCHRULE* branchrule;

   SCIP_CALL( checkStage(scip, "SCIPincludeBranchrule", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether branching rule is already present */
   if( SCIPfindBranchrule(scip, name) != NULL )
   {
      SCIPerrorMessage("branching rule <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbranchruleCreate(&branchrule, scip->mem->setmem, scip->set, name, desc, priority, maxdepth,
         maxbounddist, branchcopy, branchfree, branchinit, branchexit, branchinitsol, branchexitsol,
         branchexeclp, branchexecext, branchexecps, branchruledata) );
   SCIP_CALL( SCIPsetIncludeBranchrule(scip->set, branchrule) );

   return SCIP_OKAY;
}

/** returns the branching rule of the given name, or NULL if not existing */
SCIP_BRANCHRULE* SCIPfindBranchrule(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of branching rule */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindBranchrule", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortBranchrules(scip->set);

   return SCIPsetFindBranchrule(scip->set, name);
}

/** returns the array of currently available branching rules */
SCIP_BRANCHRULE** SCIPgetBranchrules(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBranchrules", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->branchrules;
}

/** returns the number of currently available branching rules */
int SCIPgetNBranchrules(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNBranchrules", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nbranchrules;
}

/** sets the priority of a branching rule */
SCIP_RETCODE SCIPsetBranchrulePriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   priority            /**< new priority of the branching rule */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetBranchrulePriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPbranchruleSetPriority(branchrule, scip->set, priority);

   return SCIP_OKAY;
}

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
SCIP_RETCODE SCIPsetBranchruleMaxdepth(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   int                   maxdepth            /**< new maxdepth of the branching rule */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetBranchruleMaxdepth", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPbranchruleSetMaxdepth(branchrule, maxdepth);

   return SCIP_OKAY;
}

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
SCIP_RETCODE SCIPsetBranchruleMaxbounddist(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP_Real             maxbounddist        /**< new maxbounddist of the branching rule */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetBranchruleMaxbounddist", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPbranchruleSetMaxbounddist(branchrule, maxbounddist);

   return SCIP_OKAY;
}

/** creates a display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of display column */
   const char*           desc,               /**< description of display column */
   const char*           header,             /**< head line of display column */
   SCIP_DISPSTATUS       dispstatus,         /**< display activation status of display column */
   SCIP_DECL_DISPCOPY    ((*dispcopy)),      /**< copy method of display column or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   SCIP_DECL_DISPINIT    ((*dispinit)),      /**< initialize display column */
   SCIP_DECL_DISPEXIT    ((*dispexit)),      /**< deinitialize display column */
   SCIP_DECL_DISPINITSOL ((*dispinitsol)),   /**< solving process initialization method of display column */
   SCIP_DECL_DISPEXITSOL ((*dispexitsol)),   /**< solving process deinitialization method of display column */
   SCIP_DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   SCIP_DISPDATA*        dispdata,           /**< display column data */
   int                   width,              /**< width of display column (no. of chars used) */
   int                   priority,           /**< priority of display column */
   int                   position,           /**< relative position of display column */
   SCIP_Bool             stripline           /**< should the column be separated with a line from its right neighbor? */
   )
{
   SCIP_DISP* disp;

   SCIP_CALL( checkStage(scip, "SCIPincludeDisp", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether display column is already present */
   if( SCIPfindDisp(scip, name) != NULL )
   {
      SCIPerrorMessage("display column <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPdispCreate(&disp, scip->set, scip->mem->setmem,
         name, desc, header, dispstatus, 
         dispcopy,
         dispfree, dispinit, dispexit, dispinitsol, dispexitsol, dispoutput, dispdata,
         width, priority, position, stripline) );
   SCIP_CALL( SCIPsetIncludeDisp(scip->set, disp) );

   return SCIP_OKAY;
}

/** returns the display column of the given name, or NULL if not existing */
SCIP_DISP* SCIPfindDisp(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of display column */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindDisp", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindDisp(scip->set, name);
}

/** returns the array of currently available display columns */
SCIP_DISP** SCIPgetDisps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->disps;
}

/** returns the number of currently available display columns */
int SCIPgetNDisps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->ndisps;
}

/** automatically selects display columns for being shown w.r.t. the display width parameter */
SCIP_RETCODE SCIPautoselectDisps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPautoselectDisps", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPdispAutoActivate(scip->set) );

   return SCIP_OKAY;
}

/** method to call, when the priority of an NLPI was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdNlpiPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetSetPriorityNlpi() to mark the nlpis unsorted */
   SCIP_CALL( SCIPsetNlpiPriority(scip, (SCIP_NLPI*)paramdata, SCIPparamGetInt(param)) );

   return SCIP_OKAY;
}

/** includes an NLPI in SCIP */
SCIP_RETCODE SCIPincludeNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi                /**< NLPI data structure */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(nlpi != NULL);

   SCIP_CALL( checkStage(scip, "SCIPincludeNlpi", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* check whether NLPI is already present */
   if( SCIPfindNlpi(scip, SCIPnlpiGetName(nlpi)) != NULL )
   {
      SCIPerrorMessage("NLPI <%s> already included.\n", SCIPnlpiGetName(nlpi));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPsetIncludeNlpi(scip->set, nlpi) );

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "nlpi/%s/priority", SCIPnlpiGetName(nlpi));
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of NLPI <%s>", SCIPnlpiGetName(nlpi));
   SCIP_CALL( SCIPaddIntParam(scip, paramname, paramdesc,
         NULL, FALSE, SCIPnlpiGetPriority(nlpi), INT_MIN/4, INT_MAX/4,
         paramChgdNlpiPriority, (SCIP_PARAMDATA*)nlpi) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** returns the NLPI of the given name, or NULL if not existing */
SCIP_NLPI* SCIPfindNlpi(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of NLPI */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindNlpi", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetFindNlpi(scip->set, name);
}

/** returns the array of currently available NLPIs (sorted by priority) */
SCIP_NLPI** SCIPgetNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNlpis", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSortNlpis(scip->set);

   return scip->set->nlpis;
}

/** returns the number of currently available NLPIs */
int SCIPgetNNlpis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNlpis", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nnlpis;
}

/** sets the priority of an NLPI */
SCIP_RETCODE SCIPsetNlpiPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI*            nlpi,               /**< NLPI */
   int                   priority            /**< new priority of the NLPI */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNlpiPriority", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetSetPriorityNlpi(scip->set, nlpi, priority);

   return SCIP_OKAY;
}

/** includes information about an external code linked into the SCIP library */
SCIP_RETCODE SCIPincludeExternalCodeInformation(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< name of external code */
   const char*           description         /**< description of external code, or NULL */
   )
{
   assert(scip != NULL);
   assert(name != NULL);

   SCIP_CALL( checkStage(scip, "SCIPincludeExternalCodeInformation", TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsetIncludeExternalCode(scip->set, name, description) );

   return SCIP_OKAY;
}

/** returns an array of names of currently included external codes */
char** SCIPgetExternalCodeNames(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetExternalCodeNames", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->extcodenames;
}

/** returns an array of the descriptions of currently included external codes
 *
 *  @note some descriptions may be NULL
 */
char** SCIPgetExternalCodeDescriptions(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetExternalCodeDescriptions", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->extcodedescs;
}

/** returns the number of currently included information on external codes */
int SCIPgetNExternalCodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNExternalCodes", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->set->nextcodes;
}

/** prints information on external codes to a file stream */
void SCIPprintExternalCodes(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   int i;

   SCIPmessageFPrintInfo(file, "External codes: ");
   if( scip->set->nextcodes == 0 )
   {
      SCIPinfoMessage(scip, file, "none\n");
      return;
   }
   SCIPinfoMessage(scip, file, "\n");

   for( i = 0; i < scip->set->nextcodes; ++i )
   {
      SCIPinfoMessage(scip, file, "  %-20s %s\n", scip->set->extcodenames[i], scip->set->extcodedescs[i] != NULL ? scip->set->extcodedescs[i] : "");
   }
}


/*
 * user interactive dialog methods
 */

/** creates and includes dialog */
SCIP_RETCODE SCIPincludeDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog,             /**< pointer to store the dialog */
   SCIP_DECL_DIALOGCOPY  ((*dialogcopy)),    /**< copy method of dialog or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_DIALOGEXEC  ((*dialogexec)),    /**< execution method of dialog */
   SCIP_DECL_DIALOGDESC  ((*dialogdesc)),    /**< description output method of dialog, or NULL */
   SCIP_DECL_DIALOGFREE  ((*dialogfree)),    /**< destructor of dialog to free user data, or NULL */
   const char*           name,               /**< name of dialog: command name appearing in parent's dialog menu */
   const char*           desc,               /**< description of dialog used if description output method is NULL */
   SCIP_Bool             issubmenu,          /**< is the dialog a submenu? */
   SCIP_DIALOGDATA*      dialogdata          /**< user defined dialog data */
   )
{
   assert(scip != NULL);
   assert(dialog != NULL);

   SCIP_CALL( checkStage(scip, "SCIPincludeDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* check whether display column is already present */
   if( dialogcopy != NULL && SCIPexistsDialog(scip, *dialog) )
   {
      SCIPerrorMessage("dialog <%s> already included.\n", name);
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPdialogCreate(dialog, dialogcopy, dialogexec, dialogdesc, dialogfree, name, desc, issubmenu, dialogdata) );
   SCIP_CALL( SCIPsetIncludeDialog(scip->set, *dialog) );

   return SCIP_OKAY;
}

/** returns if the dialog already exists */
SCIP_Bool SCIPexistsDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPexistsDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetExistsDialog(scip->set, dialog);
}

/** captures a dialog */
SCIP_RETCODE SCIPcaptureDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcaptureDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPdialogCapture(dialog);

   return SCIP_OKAY;
}

/** releases a dialog */
SCIP_RETCODE SCIPreleaseDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG**         dialog              /**< pointer to the dialog */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPreleaseDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPdialogRelease(scip, dialog) );

   return SCIP_OKAY;
}

/** makes given dialog the root dialog of SCIP's interactive user shell; captures dialog and releases former root dialog */
SCIP_RETCODE SCIPsetRootDialog(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog              /**< dialog to be the root */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetRootDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPdialoghdlrSetRoot(scip, scip->dialoghdlr, dialog) );

   return SCIP_OKAY;
}

/** returns the root dialog of SCIP's interactive user shell */
SCIP_DIALOG* SCIPgetRootDialog(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRootDialog", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPdialoghdlrGetRoot(scip->dialoghdlr);
}

/** adds a sub dialog to the given dialog as menu entry and captures it */
SCIP_RETCODE SCIPaddDialogEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DIALOG*          dialog,             /**< dialog to extend, or NULL for root dialog */
   SCIP_DIALOG*          subdialog           /**< subdialog to add as menu entry in dialog */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddDialogEntry", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( dialog == NULL )
      dialog = SCIPdialoghdlrGetRoot(scip->dialoghdlr);

   SCIP_CALL( SCIPdialogAddEntry(dialog, scip->set, subdialog) );

   return SCIP_OKAY;
}

/** adds a single line of input which is treated as if the user entered the command line */
SCIP_RETCODE SCIPaddDialogInputLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddDialogInputLine", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPdialoghdlrAddInputLine(scip->dialoghdlr, inputline) );

   return SCIP_OKAY;
}

/** adds a single line of input to the command history which can be accessed with the cursor keys */
SCIP_RETCODE SCIPaddDialogHistoryLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           inputline           /**< input line to add */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddDialogHistoryLine", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPdialoghdlrAddHistory(scip->dialoghdlr, NULL, inputline, FALSE) );

   return SCIP_OKAY;
}

/** starts interactive mode of SCIP by executing the root dialog */
SCIP_RETCODE SCIPstartInteraction(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstartInteraction", TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /** includes or updates the default dialog menus in SCIP */
   SCIP_CALL( SCIPincludeDialogDefault(scip) );

   SCIP_CALL( SCIPdialoghdlrExec(scip->dialoghdlr, scip->set) );

   return SCIP_OKAY;
}




/*
 * global problem methods
 */

/** creates empty problem and initializes all solving data structures (the objective sense is set to MINIMIZE)
 *  If the problem type requires the use of variable pricers, these pricers should be added to the problem with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 */
SCIP_RETCODE SCIPcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */
   SCIP_DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   SCIP_DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   SCIP_DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   SCIP_DECL_PROBINITSOL ((*probinitsol)),   /**< solving process initialization method of transformed data */
   SCIP_DECL_PROBEXITSOL ((*probexitsol)),   /**< solving process deinitialization method of transformed data */
   SCIP_DECL_PROBCOPY    ((*probcopy)),      /**< copies user data if you want to copy it to a subscip, or NULL */
   SCIP_PROBDATA*        probdata            /**< user problem data set by the reader */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateProb", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   /* free old problem */
   SCIP_CALL( SCIPfreeProb(scip) );
   assert(scip->set->stage == SCIP_STAGE_INIT);

   /* switch stage to PROBLEM */
   scip->set->stage = SCIP_STAGE_PROBLEM;

   SCIP_CALL( SCIPstatCreate(&scip->stat, scip->mem->probmem, scip->set) );

   SCIP_CALL( SCIPprobCreate(&scip->origprob, scip->mem->probmem, scip->set, name,
         probdelorig, probtrans, probdeltrans, probinitsol, probexitsol, probcopy, probdata, FALSE) );

   /* create solution pool for original solution candidates */
   SCIP_CALL( SCIPprimalCreate(&scip->origprimal) );
   
   return SCIP_OKAY;
}

/** reads problem from file and initializes all solving data structures */
SCIP_RETCODE SCIPreadProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< problem file name */
   const char*           extension           /**< extension of the desired file reader, 
                                              *   or NULL if file extension should be used */
   )
{
   SCIP_RETCODE retcode;
   SCIP_RESULT result;
   SCIP_Bool usevartable;
   SCIP_Bool useconstable;
   int i;
   char* tmpfilename;
   char* fileextension;

   assert(scip != NULL);  
   assert(filename != NULL);
   
   SCIP_CALL( checkStage(scip, "SCIPreadProb", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPgetBoolParam(scip, "misc/usevartable", &usevartable) );
   SCIP_CALL( SCIPgetBoolParam(scip, "misc/useconstable", &useconstable) );

   if( !usevartable || !useconstable )
   {
      SCIPerrorMessage("Cannot read problem if vartable or constable is disabled. Make sure parameters 'misc/usevartable' and 'misc/useconstable' are set to TRUE.\n");
      return SCIP_READERROR;
   }

   /* try all readers until one could read the file */
   result = SCIP_DIDNOTRUN;

   /* copy filename */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &tmpfilename, filename, (int)strlen(filename)+1) );
   
   fileextension = NULL;
   if( extension == NULL )
   {
      /* get extension from filename */
      SCIPsplitFilename(tmpfilename, NULL, NULL, &fileextension, NULL);
   }
   
   for( i = 0; i < scip->set->nreaders && result == SCIP_DIDNOTRUN; ++i )
   {
      retcode = SCIPreaderRead(scip->set->readers[i], scip->set, filename, 
         extension != NULL ? extension : fileextension, &result);

      /* check for reader errors */
      if( retcode == SCIP_NOFILE || retcode == SCIP_READERROR )
         goto TERMINATE;
      SCIP_CALL( retcode );
   }

   switch( result )
   {
   case SCIP_DIDNOTRUN:
      retcode = SCIP_PLUGINNOTFOUND;
      break;
   case SCIP_SUCCESS:
      if( scip->origprob != NULL )
      {
         SCIP_Real readingtime;

         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
            "original problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
            scip->origprob->nvars, scip->origprob->nbinvars, scip->origprob->nintvars,
            scip->origprob->nimplvars, scip->origprob->ncontvars,
            scip->origprob->nconss);

         /* get reading time */
         readingtime = SCIPgetReadingTime(scip);
         
         /* display timing statistics */
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "Reading Time: %.2f\n", readingtime);
      }
      retcode = SCIP_OKAY;
      break;
   default:
      assert(i < scip->set->nreaders);
      SCIPerrorMessage("invalid result code <%d> from reader <%s> reading file <%s>\n",
         result, SCIPreaderGetName(scip->set->readers[i]), filename);
      retcode = SCIP_READERROR;
   }  /*lint !e788*/

 TERMINATE:
   /* free buffer array */
   SCIPfreeBufferArray(scip, &tmpfilename);

   /* check if reading time should belong to solving time */
   if( scip->set->time_reading )
   {
      SCIP_Real readingtime;

      /* get reading time */
      readingtime = SCIPgetReadingTime(scip);
    
      /* add reading time to solving time */
      SCIPclockSetTime(scip->stat->solvingtime, readingtime);
   }
   
   return retcode;
}

/* write original or transformed problem */
static
SCIP_RETCODE writeProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader, 
                                              *   or NULL if file extension should be used */
   SCIP_Bool             transformed,        /**< output the transformed problem? */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;
   char* tmpfilename;
   char* fileextension;
   char* compression;
   FILE* file;
   
   assert(scip != NULL );

   fileextension = NULL;
   compression = NULL;
   file = NULL;
   tmpfilename = NULL;
   retcode = SCIP_OKAY;

   if( filename != NULL &&  filename[0] != '\0' )
   {
      int success;

      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPerrorMessage("cannot create file <%s> for writing\n", filename);
         SCIPprintSysError(filename);
         return SCIP_FILECREATEERROR;
      }

      /* get extension from filename,
       * if an error occurred, close the file before returning */
      if( BMSduplicateMemoryArray(&tmpfilename, filename, strlen(filename)+1) == NULL )
      {
         (void) fclose(file);
         SCIPerrorMessage("Error <%d> in function call\n", SCIP_NOMEMORY);
         return SCIP_NOMEMORY;
      }

      SCIPsplitFilename(tmpfilename, NULL, NULL, &fileextension, &compression);
      
      if( compression != NULL )
      {
         SCIPwarningMessage("currently it is not possible to write files with any compression\n");
         BMSfreeMemoryArray(&tmpfilename);
         (void) fclose(file);
         return SCIP_FILECREATEERROR;
      }
      
      if( extension == NULL && fileextension == NULL )
      {
         SCIPwarningMessage("filename <%s> has no file extension, select default <cip> format for writing\n", filename);
      }
   
      if( transformed )
         retcode = SCIPprintTransProblem(scip, file, extension != NULL ? extension : fileextension, genericnames);
      else
         retcode = SCIPprintOrigProblem(scip, file, extension != NULL ? extension : fileextension, genericnames);
      
      BMSfreeMemoryArray(&tmpfilename);
      
      success = fclose(file);
      if( success != 0 )
      {
         SCIPerrorMessage("An error occurred while closing file <%s>\n", filename);
         return SCIP_FILECREATEERROR;
      }         
   }
   else
   {
      /* print to stdout */      
      if( transformed )
         retcode = SCIPprintTransProblem(scip, NULL, extension, genericnames);
      else
         retcode = SCIPprintOrigProblem(scip, NULL, extension, genericnames);
   }

   /* check for write errors */
   if( retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** writes original problem to file  */
SCIP_RETCODE SCIPwriteOrigProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader, 
                                              *   or NULL if file extension should be used */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( checkStage(scip, "SCIPwriteOrigProblem", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   assert( scip != NULL );
   assert( scip->origprob != NULL );

   retcode = writeProblem(scip, filename, extension, FALSE, genericnames);

   /* check for write errors */
   if( retcode == SCIP_FILECREATEERROR || retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** writes transformed problem which are valid in the current node to file */
SCIP_RETCODE SCIPwriteTransProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< output file (or NULL for standard output) */
   const char*           extension,          /**< extension of the desired file reader, 
                                              *   or NULL if file extension should be used */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( checkStage(scip, "SCIPwriteTransProblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   assert( scip != NULL );
   assert( scip->transprob != NULL );

   retcode = writeProblem(scip, filename, extension, TRUE, genericnames);

   /* check for write errors */
   if( retcode == SCIP_FILECREATEERROR || retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** frees problem and solution process data */
SCIP_RETCODE SCIPfreeProb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeProb", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPfreeTransform(scip) );
   assert(scip->set->stage == SCIP_STAGE_INIT || scip->set->stage == SCIP_STAGE_PROBLEM);
   
   /* free all debug data */ 
   SCIP_CALL( SCIPdebugFreeDebugData(scip->set) );

   if( scip->set->stage == SCIP_STAGE_PROBLEM )
   {
      int i;

      /* deactivate all pricers */
      for( i = scip->set->nactivepricers-1; i >= 0; --i )
      {
         SCIP_CALL( SCIPpricerDeactivate(scip->set->pricers[i], scip->set) );
      }
      assert(scip->set->nactivepricers == 0);

      /* free original primal solution candidate pool, original problem and problem statistics data structures */
      SCIP_CALL( SCIPprimalFree(&scip->origprimal, scip->mem->probmem) );
      SCIP_CALL( SCIPprobFree(&scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
      SCIP_CALL( SCIPstatFree(&scip->stat, scip->mem->probmem) );

      /* readers */
      for( i = 0; i < scip->set->nreaders; ++i )
      {
         SCIP_CALL( SCIPreaderResetReadingTime(scip->set->readers[i]) );
      }

      /* switch stage to INIT */
      scip->set->stage = SCIP_STAGE_INIT;
   }
   assert(scip->set->stage == SCIP_STAGE_INIT);

   return SCIP_OKAY;
}


/** permutes parts of the problem data structure */
SCIP_RETCODE SCIPpermuteProb(
   SCIP*                 scip,              /**< SCIP data structure */
   unsigned int          randseed,          /**< seed value for random generator */
   SCIP_Bool             permuteconss,      /**< should the list of constraints in each constraint handler be permuted? */
   SCIP_Bool             permutebinvars,    /**< should the list of binary variables be permuted? */
   SCIP_Bool             permuteintvars,    /**< should the list of integer variables be permuted? */
   SCIP_Bool             permuteimplvars,   /**< should the list of implicit integer variables be permuted? */
   SCIP_Bool             permutecontvars    /**< should the list of continuous integer variables be permuted? */
   )
{
   SCIP_VAR** vars;
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;   
   int nbinvars;
   int nintvars;
   int nimplvars;
   int nvars;
   int j;

   assert(scip != NULL);
   SCIP_CALL( checkStage(scip, "SCIPpermuteProb", FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, &nintvars, &nimplvars, NULL) );
      
   assert(nvars == 0 || vars != NULL);
   assert(nvars == nbinvars+nintvars+nimplvars+SCIPgetNContVars(scip));

   conshdlrs = SCIPgetConshdlrs(scip);
   nconshdlrs = SCIPgetNConshdlrs(scip);
   assert(nconshdlrs == 0 || conshdlrs != NULL);
   
   /*@note the constraint handler should not be permuted since they are called w.r.t. to certain properties; besides
    *      that the "conshdlrs" array should stay in the order as it is since this array is used to copy the plugins for
    *      sub-SCIPs and contains the dependencies between the constraint handlers; for example the linear constraint
    *      handler stays in front of all constraint handler which can upgrade a linear constraint (such as logicor,
    *      setppc, and knapsack)
    */

   /* for each constraint handler, permute its constraints */
   if( permuteconss )
   {
      int i;
         
      /* loop over all constraint handlers */
      for( i = 0; i < nconshdlrs; ++i )      
      {
         SCIP_CONS** conss;
         int nconss;

         conss = SCIPconshdlrGetConss(conshdlrs[i]);
         nconss = SCIPconshdlrGetNConss(conshdlrs[i]);
         assert(nconss == 0 || conss != NULL);

         SCIPpermuteArray((void**)conss, 0, nconss, &randseed);

         /* readjust the mapping of constraints to array positions */
         for( j = 0; j < nconss; ++j )      
            conss[j]->consspos = j;
      }
   }

   /* permute binary variables */   
   if( permutebinvars )
   {
      SCIPpermuteArray((void**)vars, 0, nbinvars, &randseed);
      
      /* readjust the mapping of variables to array positions */
      for( j = 0; j < nbinvars; ++j )      
         vars[j]->probindex = j;
   }

   /* permute general integer variables */
   if( permuteintvars )
   {
      SCIPpermuteArray((void**)vars, nbinvars, nbinvars+nintvars, &randseed);

      /* readjust the mapping of variables to array positions */
      for( j = nbinvars; j < nbinvars+nintvars; ++j )      
         vars[j]->probindex = j;
   }

   /* permute general integer variables */
   if( permuteimplvars )
   {
      SCIPpermuteArray((void**)vars, nbinvars+nintvars, nbinvars+nintvars+nimplvars, &randseed);

      /* readjust the mapping of variables to array positions */
      for( j = nbinvars+nintvars; j < nbinvars+nintvars+nimplvars; ++j )      
         vars[j]->probindex = j;
   }

   /* permute general integer variables */
   if( permutecontvars )
   {
      SCIPpermuteArray((void**)vars, nbinvars+nintvars+nimplvars, nvars, &randseed);

      /* readjust the mapping of variables to array positions */
      for( j = nbinvars+nintvars+nimplvars; j < nvars; ++j )      
         vars[j]->probindex = j;
   }

   return SCIP_OKAY;
}

/** gets user problem data */
SCIP_PROBDATA* SCIPgetProbData(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobGetData(scip->origprob);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      return SCIPprobGetData(scip->transprob);

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** sets user problem data */
SCIP_RETCODE SCIPsetProbData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< user problem data to use */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetProbData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetData(scip->origprob, probdata);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      SCIPprobSetData(scip->transprob, probdata);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets name of the current problem instance */
const char* SCIPgetProbName(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetProbName", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPprobGetName(scip->origprob);
}

/** sets name of the current problem instance */
SCIP_RETCODE SCIPsetProbName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name to be set */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetProbName", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPprobSetName(scip->origprob, name);
}

/** gets objective sense of original problem */
SCIP_OBJSENSE SCIPgetObjsense(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetObjsense", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->objsense;
}

/** sets objective sense of problem */
SCIP_RETCODE SCIPsetObjsense(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_OBJSENSE         objsense            /**< new objective sense */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetObjsense", FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   if( objsense != SCIP_OBJSENSE_MAXIMIZE && objsense != SCIP_OBJSENSE_MINIMIZE )
   {
      SCIPerrorMessage("invalid objective sense\n");
      return SCIP_INVALIDDATA;
   }

   SCIPprobSetObjsense(scip->origprob, objsense);

   return SCIP_OKAY;
}

/** adds offset of objective function */
SCIP_RETCODE SCIPaddObjoffset(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             addval              /**< value to add to objective offset */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddObjoffset", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   
   scip->transprob->objoffset += addval;
   
   return SCIP_OKAY;
}

/** returns the objective offset of the original problem */
SCIP_Real SCIPgetOrigObjoffset(
   SCIP*                 scip                /**< SCIP data structure */
   )   
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetOrigObjoffset", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   
   return scip->origprob->objoffset;
}

/** returns the objective scale of the original problem */
SCIP_Real SCIPgetOrigObjscale(
   SCIP*                 scip                /**< SCIP data structure */
   )   
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetOrigObjscale", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   
   return scip->origprob->objscale;
}

/** returns the objective offset of the transformed problem */
SCIP_Real SCIPgetTransObjoffset(
   SCIP*                 scip                /**< SCIP data structure */
   )   
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetTransObjoffset", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   
   return scip->transprob->objoffset;
}

/** returns the objective scale of the transformed problem */
SCIP_Real SCIPgetTransObjscale(
   SCIP*                 scip                /**< SCIP data structure */
   )   
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetTransObjscale", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   
   return scip->transprob->objscale;
}

/** sets limit on objective function, such that only solutions better than this limit are accepted */
SCIP_RETCODE SCIPsetObjlimit(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             objlimit            /**< new primal objective limit */
   )
{
   SCIP_Real oldobjlimit;

   SCIP_CALL( checkStage(scip, "SCIPsetObjlimit", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetObjlim(scip->origprob, objlimit);
      break;
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      oldobjlimit = SCIPprobGetObjlim(scip->origprob, scip->set);
      assert(oldobjlimit == SCIPprobGetObjlim(scip->transprob, scip->set)); /*lint !e777*/
      if( SCIPtransformObj(scip, objlimit) > SCIPprobInternObjval(scip->transprob, scip->set, oldobjlimit) )
      {
         SCIPerrorMessage("cannot relax objective limit from %.15g to %.15g after problem was transformed\n", oldobjlimit, objlimit);
         return SCIP_INVALIDDATA;
      }
      SCIPprobSetObjlim(scip->origprob, objlimit);
      SCIPprobSetObjlim(scip->transprob, objlimit);
      SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob,
            scip->tree, scip->lp) );
      break;
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets current limit on objective function */
SCIP_Real SCIPgetObjlimit(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetObjlimit", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobGetObjlim(scip->origprob, scip->set);
}

/** informs SCIP, that the objective value is always integral in every feasible solution */
SCIP_RETCODE SCIPsetObjIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetObjIntegral", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIPprobSetObjIntegral(scip->origprob);
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIPprobSetObjIntegral(scip->transprob);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** returns whether the objective value is known to be integral in every feasible solution */
SCIP_Bool SCIPisObjIntegral(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisObjIntegral", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobIsObjIntegral(scip->origprob);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      return SCIPprobIsObjIntegral(scip->transprob);

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }  /*lint !e788*/
}

/** returns the Euclidean norm of the objective function vector (available only for transformed problem) */
SCIP_Real SCIPgetObjNorm(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetObjNorm", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( scip->lp->objsqrnormunreliable )                                                          
      SCIPlpRecalculateObjSqrNorm(scip->set, scip->lp);                                          
   assert(!scip->lp->objsqrnormunreliable); 

   return SCIPlpGetObjNorm(scip->lp);
}

/** adds variable to the problem */
SCIP_RETCODE SCIPaddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to add */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* avoid inserting the same variable twice */
   if( SCIPvarGetProbindex(var) != -1 )
      return SCIP_OKAY;

   /* insert the negation variable x instead of the negated variable x' in x' = offset - x */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetNegationVar(var) != NULL);
      SCIP_CALL( SCIPaddVar(scip, SCIPvarGetNegationVar(var)) );
      return SCIP_OKAY;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot add transformed variables to original problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobAddVar(scip->origprob, scip->mem->probmem, scip->set, scip->lp, scip->branchcand,
            scip->eventfilter, scip->eventqueue, var) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      /* check variable's status */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot add original variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      else if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot add fixed or aggregated variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobAddVar(scip->transprob, scip->mem->probmem, scip->set, scip->lp,
            scip->branchcand, scip->eventfilter, scip->eventqueue, var) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** adds variable to the problem and uses it as pricing candidate to enter the LP */
SCIP_RETCODE SCIPaddPricedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to add */
   SCIP_Real             score               /**< pricing score of variable (the larger, the better the variable) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddPricedVar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* insert the negation variable x instead of the negated variable x' in x' = offset - x */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetNegationVar(var) != NULL);
      SCIP_CALL( SCIPaddPricedVar(scip, SCIPvarGetNegationVar(var), score) );
      return SCIP_OKAY;
   }

   /* add variable to problem if not yet inserted */
   if( SCIPvarGetProbindex(var) == -1 )
   {
      /* check variable's status */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot add original variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      else if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot add fixed or aggregated variables to transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobAddVar(scip->transprob, scip->mem->probmem, scip->set, scip->lp,
            scip->branchcand, scip->eventfilter, scip->eventqueue, var) );
   }

   /* add variable to pricing storage */
   SCIP_CALL( SCIPpricestoreAddVar(scip->pricestore, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, var, score,
         (SCIPtreeGetCurrentDepth(scip->tree) == 0)) );

   return SCIP_OKAY;
}

/** removes variable from the problem; however, the variable is NOT removed from the constraints */
SCIP_RETCODE SCIPdelVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to delete */
   SCIP_Bool*            deleted             /**< pointer to store whether marking variable to be deleted was successful */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(deleted != NULL);

   SCIP_CALL( checkStage(scip, "SCIPdelVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot remove transformed variables from original problem\n");
         return SCIP_INVALIDDATA;
      }
      SCIP_CALL( SCIPprobDelVar(scip->origprob, scip->mem->probmem, scip->set, scip->eventqueue, var, deleted) );

      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_FREESOLVE:
      /* check variable's status */
      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL )
      {
         SCIPerrorMessage("cannot remove original variables from transformed problem\n");
         return SCIP_INVALIDDATA;
      }
      else if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_LOOSE && SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot remove fixed or aggregated variables from transformed problem\n");
         return SCIP_INVALIDDATA;
      }

      /* fix the variable to 0, first */
      assert(!SCIPisFeasPositive(scip, SCIPvarGetLbGlobal(var)));
      assert(!SCIPisFeasNegative(scip, SCIPvarGetUbGlobal(var)));

      if ( !SCIPisFeasZero(scip, SCIPvarGetLbGlobal(var)) )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(scip, var, 0.0));
      }
      if ( !SCIPisFeasZero(scip, SCIPvarGetUbGlobal(var)) )
      {
         SCIP_CALL( SCIPchgVarUbGlobal(scip, var, 0.0));
      }

      SCIP_CALL( SCIPprobDelVar(scip->transprob, scip->mem->probmem, scip->set, scip->eventqueue, var, deleted) );

      return SCIP_OKAY;
   case SCIP_STAGE_FREETRANS:
      /* in FREETRANS stage, we don't need to remove the variable, because the transformed problem is freed anyways */
      *deleted = FALSE;

      return SCIP_OKAY;
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets variables of the problem along with the numbers of different variable types; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
SCIP_RETCODE SCIPgetVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetVarsData", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      if( vars != NULL )
         *vars = scip->origprob->vars;
      if( nvars != NULL )
         *nvars = scip->origprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->origprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->origprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->origprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->origprob->ncontvars;
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      if( vars != NULL )
         *vars = scip->transprob->vars;
      if( nvars != NULL )
         *nvars = scip->transprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->transprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->transprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->transprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->transprob->ncontvars;
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets array with active problem variables; data may become invalid after
 *  calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
SCIP_VAR** SCIPgetVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->vars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->vars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of active problem variables */
int SCIPgetNVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of binary active problem variables */
int SCIPgetNBinVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNBinVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nbinvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nbinvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of integer active problem variables */
int SCIPgetNIntVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNIntVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nintvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nintvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of implicit integer active problem variables */
int SCIPgetNImplVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNImplVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nimplvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nimplvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of continuous active problem variables */
int SCIPgetNContVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNContVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->ncontvars;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->ncontvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets array with fixed and aggregated problem variables; data may become invalid after
 *  calls to SCIPfixVar(), SCIPaggregateVars(), and SCIPmultiaggregateVar()
 */
SCIP_VAR** SCIPgetFixedVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetFixedVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return NULL;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->fixedvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of fixed or aggregated problem variables */
int SCIPgetNFixedVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNFixedVars", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return 0;
      
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nfixedvars;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets variables of the original problem along with the numbers of different variable types; data may become invalid
 *  after a call to SCIPchgVarType()
 */
SCIP_RETCODE SCIPgetOrigVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetOrigVarsData", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( vars != NULL )
      *vars = scip->origprob->vars;
   if( nvars != NULL )
      *nvars = scip->origprob->nvars;
   if( nbinvars != NULL )
      *nbinvars = scip->origprob->nbinvars;
   if( nintvars != NULL )
      *nintvars = scip->origprob->nintvars;
   if( nimplvars != NULL )
      *nimplvars = scip->origprob->nimplvars;
   if( ncontvars != NULL )
      *ncontvars = scip->origprob->ncontvars;

   return SCIP_OKAY;
}

/** gets array with original problem variables; data may become invalid after
 *  a call to SCIPchgVarType()
 */
SCIP_VAR** SCIPgetOrigVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->vars;
}

/** gets number of original problem variables */
int SCIPgetNOrigVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNOrigVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nvars;
}

/** gets number of binary original problem variables */
int SCIPgetNOrigBinVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNOrigBinVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nbinvars;
}

/** gets number of integer original problem variables */
int SCIPgetNOrigIntVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNOrigIntVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nintvars;
}

/** gets number of implicit integer original problem variables */
int SCIPgetNOrigImplVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNOrigImplVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nimplvars;
}

/** gets number of continuous original problem variables */
int SCIPgetNOrigContVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNOrigContVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->ncontvars;
}

/** gets number of all problem variables created during creation and solving of problem;
 *  this includes also variables that were deleted in the meantime
 */
int SCIPgetNTotalVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNTotalVars", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   assert(scip->stat != NULL);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      return scip->stat->nvaridx;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/

}


/** gets variables of the original or transformed problem along with the numbers of different variable types;
 *  the returned problem space (original or transformed) corresponds to the given solution;
 *  data may become invalid after calls to SCIPchgVarType(), SCIPfixVar(), SCIPaggregateVars(), and
 *  SCIPmultiaggregateVar()
 */
SCIP_RETCODE SCIPgetSolVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that selects the problem space, NULL for current solution */
   SCIP_VAR***           vars,               /**< pointer to store variables array or NULL if not needed */
   int*                  nvars,              /**< pointer to store number of variables or NULL if not needed */
   int*                  nbinvars,           /**< pointer to store number of binary variables or NULL if not needed */
   int*                  nintvars,           /**< pointer to store number of integer variables or NULL if not needed */
   int*                  nimplvars,          /**< pointer to store number of implicit integral vars or NULL if not needed */
   int*                  ncontvars           /**< pointer to store number of continuous variables or NULL if not needed */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetSolVarsData", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( scip->set->stage == SCIP_STAGE_PROBLEM || (sol != NULL && SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL) )
   {
      if( vars != NULL )
         *vars = scip->origprob->vars;
      if( nvars != NULL )
         *nvars = scip->origprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->origprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->origprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->origprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->origprob->ncontvars;
   }
   else
   {
      if( vars != NULL )
         *vars = scip->transprob->vars;
      if( nvars != NULL )
         *nvars = scip->transprob->nvars;
      if( nbinvars != NULL )
         *nbinvars = scip->transprob->nbinvars;
      if( nintvars != NULL )
         *nintvars = scip->transprob->nintvars;
      if( nimplvars != NULL )
         *nimplvars = scip->transprob->nimplvars;
      if( ncontvars != NULL )
         *ncontvars = scip->transprob->ncontvars;
   }

   return SCIP_OKAY;
}

/** returns variable of given name in the problem, or NULL if not existing */
SCIP_VAR* SCIPfindVar(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of variable to find */
   )
{
   SCIP_VAR* var;

   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindVar(scip->origprob, name);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      var = SCIPprobFindVar(scip->transprob, name);
      if( var == NULL )
         return SCIPprobFindVar(scip->origprob, name);
      else
         return var;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** returns TRUE iff all potential variables exist in the problem, and FALSE, if there may be additional variables,
 *  that will be added in pricing and improve the objective value
 */
SCIP_Bool SCIPallVarsInProb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPallVarsInProb", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return (scip->set->nactivepricers == 0);
}

/** adds constraint to the problem; if constraint is only valid locally, it is added to the local subproblem of the
 *  current node (and all of its subnodes); otherwise it is added to the global problem;
 *  if a local constraint is added at the root node, it is automatically upgraded into a global constraint
 */
SCIP_RETCODE SCIPaddCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to add */
   )
{
   assert(cons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddCons", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPprobAddCons(scip->origprob, scip->set, scip->stat, cons) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      assert( SCIPtreeGetCurrentDepth(scip->tree) >= 0 ||  scip->set->stage == SCIP_STAGE_PRESOLVED );
      if( SCIPtreeGetCurrentDepth(scip->tree) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
         SCIPconsSetLocal(cons, FALSE);
      if( SCIPconsIsGlobal(cons) )
      {
         SCIP_CALL( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
      }
      else
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) > SCIPtreeGetEffectiveRootDepth(scip->tree));
         SCIP_CALL( SCIPnodeAddCons(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->tree, cons) );
      }
      return SCIP_OKAY;

   case SCIP_STAGE_FREESOLVE:
      SCIP_CALL( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was added, or from the problem, if it was a problem constraint
 */
SCIP_RETCODE SCIPdelCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to delete */
   )
{
   assert(cons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPdelCons", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->origprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** returns original constraint of given name in the problem, or NULL if not existing */
SCIP_CONS* SCIPfindOrigCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint to find */
   )
{
   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindOrigCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      return SCIPprobFindCons(scip->origprob, name);
   
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** returns constraint of given name in the problem, or NULL if not existing */
SCIP_CONS* SCIPfindCons(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint to find */
   )
{
   SCIP_CONS* cons;

   assert(name != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfindCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return SCIPprobFindCons(scip->origprob, name);

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      cons = SCIPprobFindCons(scip->transprob, name);
      if( cons == NULL )
         return SCIPprobFindCons(scip->origprob, name);
      else
         return cons;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets number of upgraded constraints */
int SCIPgetNUpgrConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNUpgrConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return 0;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->stat->npresolupgdconss;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets total number of globally valid constraints currently in the problem */
int SCIPgetNConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->nconss;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->nconss;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return 0; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets array of globally valid constraints currently in the problem */
SCIP_CONS** SCIPgetConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetConss", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprob->conss;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      return scip->transprob->conss;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      SCIPABORT();
      return NULL; /*lint !e527*/
   }  /*lint !e788*/
}

/** gets total number of constraints in the original problem */
int SCIPgetNOrigConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNOrigConss", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->nconss;
}

/** gets array of constraints in the original problem */
SCIP_CONS** SCIPgetOrigConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetOrigConss", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->origprob->conss;
}




/*
 * local subproblem methods
 */

/** adds constraint to the given node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *  In this case, one should pass the more global node where the constraint is valid as "validnode".
 *  Note that the same constraint cannot be added twice to the branching tree with different "validnode" parameters.
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 */
SCIP_RETCODE SCIPaddConsNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to add constraint to */
   SCIP_CONS*            cons,               /**< constraint to add */
   SCIP_NODE*            validnode           /**< node at which the constraint is valid, or NULL */
   )
{
   assert(cons != NULL);
   assert(node != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( validnode != NULL )
   {
      int validdepth;

      validdepth = SCIPnodeGetDepth(validnode);
      if( validdepth > SCIPnodeGetDepth(node) )
      {
         SCIPerrorMessage("cannot add constraint <%s> valid in depth %d to a node of depth %d\n",
            SCIPconsGetName(cons), validdepth, SCIPnodeGetDepth(node));
         return SCIP_INVALIDDATA;
      }
      if( cons->validdepth != -1 && cons->validdepth != validdepth )
      {
         SCIPerrorMessage("constraint <%s> is already marked to be valid in depth %d - cannot mark it to be valid in depth %d\n",
            SCIPconsGetName(cons), cons->validdepth, validdepth);
         return SCIP_INVALIDDATA;
      }
      if( validdepth <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
         SCIPconsSetLocal(cons, FALSE);
      else
         cons->validdepth = validdepth;
   }

   if( SCIPnodeGetDepth(node) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
   {
      SCIPconsSetLocal(cons, FALSE);
      SCIP_CALL( SCIPprobAddCons(scip->transprob, scip->set, scip->stat, cons) );
   }
   else
   {
      SCIP_CALL( SCIPnodeAddCons(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, cons) );
   }

   return SCIP_OKAY;
}

/** adds constraint locally to the current node (and all of its subnodes), even if it is a global constraint;
 *  It is sometimes desirable to add the constraint to a more local node (i.e., a node of larger depth) even if
 *  the constraint is also valid higher in the tree, for example, if one wants to produce a constraint which is
 *  only active in a small part of the tree although it is valid in a larger part.
 *  In this case, one should pass the more global node where the constraint is valid as "validnode".
 *  Note that the same constraint cannot be added twice to the branching tree with different "validnode" parameters.
 *  If the constraint is valid at the same node as it is inserted (the usual case), one should pass NULL as "validnode".
 *  If the "validnode" is the root node, it is automatically upgraded into a global constraint, but still only added to
 *  the given node. If a local constraint is added to the root node, it is added to the global problem instead.
 */
SCIP_RETCODE SCIPaddConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to add */
   SCIP_NODE*            validnode           /**< node at which the constraint is valid, or NULL */
   )
{
   assert(cons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddConsLocal", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPaddConsNode(scip, SCIPtreeGetCurrentNode(scip->tree), cons, validnode) );

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the given node (and all subnodes);
 *  if the method is called at the root node, the constraint is globally deleted from the problem;
 *  the constraint deletion is being remembered at the given node, s.t. after leaving the node's subtree, the constraint
 *  is automatically enabled again, and after entering the node's subtree, it is automatically disabled;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 */
SCIP_RETCODE SCIPdelConsNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to disable constraint in */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   )
{
   assert(cons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPdelConsNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPnodeGetDepth(node) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
   {
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob) );
   }
   else
   {
      SCIP_CALL( SCIPnodeDelCons(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, cons) );
   }

   return SCIP_OKAY;
}

/** disables constraint's separation, enforcing, and propagation capabilities at the current node (and all subnodes);
 *  if the method is called during problem modification or at the root node, the constraint is globally deleted from
 *  the problem;
 *  the constraint deletion is being remembered at the current node, s.t. after leaving the current subtree, the
 *  constraint is automatically enabled again, and after reentering the current node's subtree, it is automatically
 *  disabled again;
 *  this may improve performance because redundant checks on this constraint are avoided, but it consumes memory;
 *  alternatively, use SCIPdisableCons()
 */
SCIP_RETCODE SCIPdelConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to locally delete */
   )
{
   SCIP_NODE* node;

   assert(cons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPdelConsLocal", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(cons->addconssetchg == NULL);
      SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->origprob) );
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_SOLVING:
      node = SCIPtreeGetCurrentNode(scip->tree);
      if( SCIPnodeGetDepth(node) <= SCIPtreeGetEffectiveRootDepth(scip->tree) )
      {
         SCIP_CALL( SCIPconsDelete(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob) );
      }
      else
      {
         SCIP_CALL( SCIPnodeDelCons(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, cons) );
      }
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets estimate of best primal solution w.r.t. original problem contained in current subtree */
SCIP_Real SCIPgetLocalOrigEstimate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLocalOrigEstimate", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);
   return node != NULL ? SCIPprobExternObjval(scip->transprob, scip->set, SCIPnodeGetEstimate(node)) : SCIP_INVALID;
}

/** gets estimate of best primal solution w.r.t. transformed problem contained in current subtree */
SCIP_Real SCIPgetLocalTransEstimate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLocalTransEstimate", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);

   return node != NULL ? SCIPnodeGetEstimate(node) : SCIP_INVALID;
}

/** gets dual bound of current node */
SCIP_Real SCIPgetLocalDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLocalDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);
   return node != NULL ? SCIPprobExternObjval(scip->transprob, scip->set, SCIPnodeGetLowerbound(node)) : SCIP_INVALID;
}

/** gets lower bound of current node in transformed problem */
SCIP_Real SCIPgetLocalLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODE* node;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   node = SCIPtreeGetCurrentNode(scip->tree);

   return node != NULL ? SCIPnodeGetLowerbound(node) : SCIP_INVALID;
}

/** gets dual bound of given node */
SCIP_Real SCIPgetNodeDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node to get dual bound for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->set, SCIPnodeGetLowerbound(node));
}

/** gets lower bound of given node in transformed problem */
SCIP_Real SCIPgetNodeLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node to get dual bound for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPnodeGetLowerbound(node);
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the current node's dual bound,
 *  sets the current node's dual bound to the new value
 */
SCIP_RETCODE SCIPupdateLocalDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPupdateLocalDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(SCIPtreeGetCurrentNode(scip->tree), scip->stat,
      SCIPprobInternObjval(scip->transprob, scip->set, newbound));

   return SCIP_OKAY;
}

/** if given value is larger than the current node's lower bound (in transformed problem), sets the current node's
 *  lower bound to the new value
 */
SCIP_RETCODE SCIPupdateLocalLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPupdateLocalLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(SCIPtreeGetCurrentNode(scip->tree), scip->stat, newbound);

   return SCIP_OKAY;
}

/** if given value is tighter (larger for minimization, smaller for maximization) than the node's dual bound,
 *  sets the node's dual bound to the new value
 */
SCIP_RETCODE SCIPupdateNodeDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to update dual bound for */
   SCIP_Real             newbound            /**< new dual bound for the node (if it's tighter than the old one) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPupdateNodeDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(node, scip->stat, SCIPprobInternObjval(scip->transprob, scip->set, newbound));

   return SCIP_OKAY;
}

/** if given value is larger than the node's lower bound (in transformed problem), sets the node's lower bound
 *  to the new value
 */
SCIP_RETCODE SCIPupdateNodeLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to update lower bound for */
   SCIP_Real             newbound            /**< new lower bound for the node (if it's larger than the old one) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPupdateNodeLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeUpdateLowerbound(node, scip->stat, newbound);

   return SCIP_OKAY;
}

/** change the node selection priority of the given child */
SCIP_RETCODE SCIPchgChildPrio(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            child,              /**< child to update the node selection priority */
   SCIP_Real             priority            /**< node selection priority value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgChildPrio", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPnodeGetType(child) != SCIP_NODETYPE_CHILD )
      return SCIP_INVALIDDATA;
   
   SCIPchildChgNodeselPrio(scip->tree, child, priority);
   
   return SCIP_OKAY;
}




/*
 * solve methods
 */

/** checks solution for feasibility in original problem without adding it to the solution store */
static
SCIP_RETCODE checkSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool             completely,         /**< should all violation be checked? */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             checkmodifiable     /**< have modifiable constraint to be checked? */
   )
{
   SCIP_RESULT result;
   int v;
   int c;
   int h;

   assert(scip != NULL);
   assert(sol != NULL);
   assert(feasible != NULL);

   SCIP_CALL( checkStage(scip, "checkSolOrig", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   *feasible = TRUE;

   /* check bounds */
   if( checkbounds )
   {
      for( v = 0; v < scip->origprob->nvars && (*feasible || printreason); ++v )
      {
         SCIP_VAR* var;
         SCIP_Real solval;
         SCIP_Real lb;
         SCIP_Real ub;
         
         var = scip->origprob->vars[v];
         solval = SCIPsolGetVal(sol, scip->set, scip->stat, var);
         lb = SCIPvarGetLbOriginal(var);
         ub = SCIPvarGetUbOriginal(var);
         if( SCIPsetIsFeasLT(scip->set, solval, lb) || SCIPsetIsFeasGT(scip->set, solval, ub) )
         {
            *feasible = FALSE;
            
            if( printreason )
            {
               SCIPmessagePrintInfo("solution violates original bounds of variable <%s> [%g,%g] solution value <%g>\n",
                  SCIPvarGetName(var), lb, ub, solval);
            }
            
            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   /* check original constraints
    *
    * in general modifiable constraints can not be checked, because the variables to fulfill them might be missing in
    * the original problem; however, if the solution comes from a heuristic during presolving modifiable constraints
    * have to be checked;
    */
   for( c = 0; c < scip->origprob->nconss; ++c )
   {
      if( SCIPconsIsChecked(scip->origprob->conss[c]) && (checkmodifiable || !SCIPconsIsModifiable(scip->origprob->conss[c])) )
      {
         /* check solution */
         SCIP_CALL( SCIPconsCheck(scip->origprob->conss[c], scip->set, sol, 
               checkintegrality, checklprows, printreason, &result) );

         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;
            if( !completely )
               return SCIP_OKAY;
         }
      }
   }

   /* call constraint handlers that don't need constraints */
   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      if( !SCIPconshdlrNeedsCons(scip->set->conshdlrs[h]) )
      {
         SCIP_CALL( SCIPconshdlrCheck(scip->set->conshdlrs[h], scip->mem->probmem, scip->set, scip->stat, sol,
               checkintegrality, checklprows, printreason, &result) );
         if( result != SCIP_FEASIBLE )
         {
            *feasible = FALSE;
            if( !completely) 
               return SCIP_OKAY;
         }
      }
   }
   
   return SCIP_OKAY;
}

/** initializes solving data structures and transforms problem */
SCIP_RETCODE SCIPtransformProb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int nfeassols;
   int ncandsols;
   int h;
   int s;

   SCIP_CALL( checkStage(scip, "SCIPtransformProb", FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* check, if the problem was already transformed */
   if( scip->set->stage >= SCIP_STAGE_TRANSFORMED )
      return SCIP_OKAY;

   assert(scip->stat->status == SCIP_STATUS_UNKNOWN);

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      SCIPerrorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* call garbage collector on original problem and parameter settings memory spaces */
   BMSgarbagecollectBlockMemory(scip->mem->setmem);
   BMSgarbagecollectBlockMemory(scip->mem->probmem);

   /* remember number of constraints */
   SCIPprobMarkNConss(scip->origprob);

   /* switch stage to TRANSFORMING */
   scip->set->stage = SCIP_STAGE_TRANSFORMING;

   /* mark statistics before solving */
   SCIPstatMark(scip->stat);

   /* init solve data structures */
   SCIP_CALL( SCIPeventfilterCreate(&scip->eventfilter, scip->mem->probmem) );
   SCIP_CALL( SCIPeventqueueCreate(&scip->eventqueue) );
   SCIP_CALL( SCIPbranchcandCreate(&scip->branchcand) );
   SCIP_CALL( SCIPlpCreate(&scip->lp, scip->set, scip->stat, SCIPprobGetName(scip->origprob)) );
   SCIP_CALL( SCIPrelaxationCreate(&scip->relaxation) );
   SCIP_CALL( SCIPprimalCreate(&scip->primal) );
   SCIP_CALL( SCIPtreeCreate(&scip->tree, scip->set, SCIPsetGetNodesel(scip->set, scip->stat)) );
   SCIP_CALL( SCIPconflictCreate(&scip->conflict, scip->mem->probmem, scip->set) );
   SCIP_CALL( SCIPcliquetableCreate(&scip->cliquetable) );

   /* copy problem in solve memory */
   SCIP_CALL( SCIPprobTransform(scip->origprob, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->lp,
         scip->branchcand, scip->eventfilter, scip->eventqueue, &scip->transprob) );

   /* switch stage to TRANSFORMED */
   scip->set->stage = SCIP_STAGE_TRANSFORMED;
  
   /* check solution of solution candidate storage */
   nfeassols = 0;
   ncandsols = scip->origprimal->nsols;
   for( s = scip->origprimal->nsols - 1; s >= 0; --s )
   {
      SCIP_Bool feasible;
      SCIP_Bool stored;
      SCIP_SOL* sol ;

      sol =  scip->origprimal->sols[s];

      /* SCIPprimalTrySol() can only be called on transformed solutions; therefore check solutions in original problem
       * including modifiable constraints
       */
      SCIP_CALL( checkSolOrig(scip, sol, &feasible, scip->set->misc_printreason, FALSE, TRUE, TRUE, TRUE, TRUE) );
      
      if( feasible )
      {
         SCIPsolRecomputeObj(sol, scip->set, scip->stat, scip->origprob);

         /* add primal solution to solution storage by copying it */
         SCIP_CALL( SCIPprimalAddSol(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob,
               scip->tree, scip->lp, scip->eventqueue, scip->eventfilter, sol, &stored) );

         if( stored )
            nfeassols++;
      }

      SCIP_CALL( SCIPsolFree(&sol, scip->mem->probmem, scip->origprimal) );
      scip->origprimal->nsols--;
   }
   
   assert(scip->origprimal->nsols == 0);
   assert(scip->origprimal->nexistingsols == 0);
   
   if( nfeassols > 0 )
   {
      SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "%d feasible solutions given by solution candidate storage\n", nfeassols);                    
   } 
   else if( ncandsols > 0 )
   {
      SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "all solutions given by solution candidate storage are infeasible\n");                    
   } 
      

   /* update upper bound and cutoff bound due to objective limit in primal data */
   SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
         scip->transprob, scip->tree, scip->lp) );

   /* print transformed problem statistics */
   SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
      "transformed problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
      scip->transprob->nvars, scip->transprob->nbinvars, scip->transprob->nintvars, scip->transprob->nimplvars,
      scip->transprob->ncontvars, scip->transprob->nconss);

   for( h = 0; h < scip->set->nconshdlrs; ++h )
   {
      int nactiveconss;

      nactiveconss = SCIPconshdlrGetNActiveConss(scip->set->conshdlrs[h]);
      if( nactiveconss > 0 )
      {
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "%7d constraints of type <%s>\n", nactiveconss, SCIPconshdlrGetName(scip->set->conshdlrs[h]));
      }
   }
   SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL, "\n");

   /* call initialization methods of plugins */
   SCIP_CALL( SCIPsetInitPlugins(scip->set, scip->mem->probmem, scip->stat) );

   /* in case the permutation seed is different to -1, permute the transformed problem */
   if( scip->set->misc_permutationseed > -1 )
   {
      SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "permute problem using random seed %d\n", scip->set->misc_permutationseed);
      
      SCIP_CALL( SCIPpermuteProb(scip, (unsigned int)scip->set->misc_permutationseed, TRUE, TRUE, TRUE, TRUE, TRUE) );
   }

   return SCIP_OKAY;
}

/** initializes presolving */
static
SCIP_RETCODE initPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            unbounded,          /**< pointer to store whether presolving detected unboundedness */
   SCIP_Bool*            infeasible          /**< pointer to store whether presolving detected infeasibility */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);
   assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);
   assert(unbounded != NULL);
   assert(infeasible != NULL);

   *unbounded = FALSE;
   *infeasible = FALSE;

   /* retransform all existing solutions to original problem space, because the transformed problem space may
    * get modified in presolving and the solutions may become invalid for the transformed problem
    */
   SCIP_CALL( SCIPprimalRetransformSolutions(scip->primal, scip->set, scip->stat, scip->origprob) );

   /* reset statistics for presolving and current branch and bound run */
   SCIPstatResetPresolving(scip->stat);

   /* increase number of branch and bound runs */
   scip->stat->nruns++;

   /* remember problem size of previous run */
   scip->stat->prevrunnvars = scip->transprob->nvars;

   /* switch stage to PRESOLVING */
   scip->set->stage = SCIP_STAGE_PRESOLVING;

   /* create temporary presolving root node */
   SCIP_CALL( SCIPtreeCreatePresolvingRoot(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
         scip->primal, scip->lp, scip->branchcand, scip->conflict, scip->eventfilter, scip->eventqueue) );

   /* inform plugins that the presolving is abound to begin */
   SCIP_CALL( SCIPsetInitprePlugins(scip->set, scip->mem->probmem, scip->stat, unbounded, infeasible) );
   assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);

   /* remove empty and single variable cliques from the clique table, and convert all two variable cliques
    * into implications
    * delete the variables from the problems that were marked to be deleted
    */
   if( !(*unbounded) && !(*infeasible) )
   {
      SCIP_Bool infeas;

      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp, scip->branchcand) );

      SCIP_CALL( SCIPcliquetableCleanup(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, &infeas) );
      if( infeas )
      {
         *infeasible = TRUE;
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "clique table cleanup detected infeasibility\n");
      }
   }

   return SCIP_OKAY;
}

/** deinitializes presolving */
static
SCIP_RETCODE exitPresolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            unbounded,          /**< pointer to store whether presolving detected unboundedness */
   SCIP_Bool*            infeasible,         /**< pointer to store whether presolving detected infeasibility */
   SCIP_Bool             isunbounded,        /**< was unboundedness already detected */
   SCIP_Bool             isinfeasible        /**< was infeasibility already detected */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);
   assert(scip->set->stage == SCIP_STAGE_PRESOLVING);

   if( !isunbounded && !isinfeasible )
   {
      /* flatten all variables */
      vars = SCIPgetFixedVars(scip);
      nvars = SCIPgetNFixedVars(scip);
      assert(nvars == 0 || vars != NULL);

      for( v = nvars - 1; v >= 0; --v )
      {
	 SCIP_VAR* var;
#ifndef NDEBUG
	 SCIP_VAR** multvars;
	 int i;
#endif
	 var = vars[v];
	 assert(var != NULL);

	 if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
	 {
	    /** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later-on */
	    SCIP_CALL( SCIPvarFlattenAggregationGraph(var, scip->mem->probmem, scip->set) );

#ifndef NDEBUG
	    multvars = SCIPvarGetMultaggrVars(var);
	    for( i = SCIPvarGetMultaggrNVars(var) - 1; i >= 0; --i)
	       assert(SCIPvarGetStatus(multvars[i]) != SCIP_VARSTATUS_MULTAGGR);
#endif
	 }
      }
   }

   *unbounded = isunbounded;
   *infeasible = isinfeasible;

   /* inform plugins that the presolving is finished, and perform final modifications */
   SCIP_CALL( SCIPsetExitprePlugins(scip->set, scip->mem->probmem, scip->stat, unbounded, infeasible) );
   assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);

   /* remove empty and single variable cliques from the clique table, and convert all two variable cliques
    * into implications
    * delete the variables from the problems that were marked to be deleted
    */
   if( !(*unbounded) && !(*infeasible) )
   {
      SCIP_Bool infeas;

      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp, scip->branchcand) );

      SCIP_CALL( SCIPcliquetableCleanup(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, &infeas) );
      if( infeas )
      {
         *infeasible = TRUE;
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "clique table cleanup detected infeasibility\n");
      }
   }

   /* exit presolving */
   SCIP_CALL( SCIPprobExitPresolve(scip->transprob,  scip->set) );
   assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);
   
   /* check, whether objective value is always integral by inspecting the problem, if it is the case adjust the
    * cutoff bound if primal solution is already known 
    */
   SCIP_CALL( SCIPprobCheckObjIntegral(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->primal,
         scip->tree, scip->lp, scip->eventqueue) );

   /* if possible, scale objective function such that it becomes integral with gcd 1 */
   SCIP_CALL( SCIPprobScaleObj(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->primal,
         scip->tree, scip->lp, scip->eventqueue) );

   /* free temporary presolving root node */
   SCIP_CALL( SCIPtreeFreePresolvingRoot(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
         scip->primal, scip->lp, scip->branchcand, scip->conflict, scip->eventfilter, scip->eventqueue) );

   /* switch stage to PRESOLVED */
   scip->set->stage = SCIP_STAGE_PRESOLVED;

   return SCIP_OKAY;
}

/** returns whether the presolving should be aborted */
static
SCIP_Bool isPresolveFinished(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             abortfac,           /**< presolving is finished, if changed portion is smaller than this factor */
   int                   maxnrounds,         /**< maximal number of presolving rounds to perform */
   int                   lastnfixedvars,     /**< number of fixed variables in last presolving round */
   int                   lastnaggrvars,      /**< number of aggregated variables in last presolving round */
   int                   lastnchgvartypes,   /**< number of changed variable types in last presolving round */
   int                   lastnchgbds,        /**< number of changed bounds in last presolving round */
   int                   lastnaddholes,      /**< number of added holes in last presolving round */
   int                   lastndelconss,      /**< number of deleted constraints in last presolving round */
   int                   lastnaddconss,      /**< number of added constraints in last presolving round */
   int                   lastnupgdconss,     /**< number of upgraded constraints in last presolving round */
   int                   lastnchgcoefs,      /**< number of changed coefficients in last presolving round */
   int                   lastnchgsides,      /**< number of changed sides in last presolving round */
   /*int                   lastnimplications,*/  /**< number of implications in last presolving round */
   /*int                   lastncliques,*/       /**< number of cliques in last presolving round */
   SCIP_Bool             unbounded,          /**< has presolving detected unboundedness? */
   SCIP_Bool             infeasible          /**< has presolving detected infeasibility? */
   )
{
   SCIP_Bool finished;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);

   /* don't abort, if enough changes were applied to the variables */
   finished = (scip->transprob->nvars == 0
      || (scip->stat->npresolfixedvars - lastnfixedvars
         + scip->stat->npresolaggrvars - lastnaggrvars
         + scip->stat->npresolchgvartypes - lastnchgvartypes
         + (scip->stat->npresolchgbds - lastnchgbds)/10.0
         + (scip->stat->npresoladdholes - lastnaddholes)/10.0 <= abortfac * scip->transprob->nvars)); /*lint !e653*/

   /* don't abort, if enough changes were applied to the constraints */
   finished = finished
      && (scip->transprob->nconss == 0
         || (scip->stat->npresoldelconss - lastndelconss
            + scip->stat->npresoladdconss - lastnaddconss
            + scip->stat->npresolupgdconss - lastnupgdconss
            + scip->stat->npresolchgsides - lastnchgsides
            <= abortfac * scip->transprob->nconss));

   /* don't abort, if enough changes were applied to the coefficients (assume a 1% density of non-zero elements) */
   finished = finished
      && (scip->transprob->nvars == 0 || scip->transprob->nconss == 0
         || (scip->stat->npresolchgcoefs - lastnchgcoefs
            <= abortfac * 0.01 * scip->transprob->nvars * scip->transprob->nconss));

#if 0
   /* don't abort, if enough new implications or cliques were found (assume 100 implications per variable) */
   finished = finished
      && (scip->stat->nimplications - lastnimplications <= abortfac * 100 * scip->transprob->nbinvars)
      && (SCIPcliquetableGetNCliques(scip->cliquetable) - lastncliques <= abortfac * scip->transprob->nbinvars);
#endif

   /* abort if problem is infeasible or unbounded */
   finished = finished || unbounded || infeasible;

   /* abort if maximal number of presolving rounds is reached */
   finished = finished || (scip->stat->npresolrounds >= maxnrounds);

   return finished;
}

/** applies one round of presolving */
static
SCIP_RETCODE presolveRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             onlydelayed,        /**< should only delayed presolvers be called? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a presolver was delayed */
   SCIP_Bool*            unbounded,          /**< pointer to store whether presolving detected unboundedness */
   SCIP_Bool*            infeasible          /**< pointer to store whether presolving detected infeasibility */
   )
{
   SCIP_RESULT result;
   SCIP_EVENT event;
   SCIP_Bool aborted;
   int i;
   int j;
   int priopresol;
   int prioprop;
   SCIP_Bool lastranpresol;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(delayed != NULL);
   assert(unbounded != NULL);
   assert(infeasible != NULL);

   *delayed = FALSE;
   *unbounded = FALSE;
   *infeasible = FALSE;
   aborted = FALSE;
   lastranpresol = FALSE;

   /* call included presolvers with nonnegative priority */
   for( i = 0, j = 0; !(*unbounded) && !(*infeasible) && !aborted && (i < scip->set->npresols || j < scip->set->nprops);  )
   {
      if( i < scip->set->npresols )
         priopresol = SCIPpresolGetPriority(scip->set->presols[i]);
      else
         priopresol = -1;

      if( j < scip->set->nprops )
         prioprop = SCIPpropGetPresolPriority(scip->set->props[j]);
      else
         prioprop = -1;

      /* choose presolving */
      if( prioprop >= priopresol )
      {
         /* only presolving methods which have non-negative priority will be called before constraint handlers */
         if( prioprop < 0 )
            break;
         
         if( onlydelayed && !SCIPpropWasPresolDelayed(scip->set->props[j]) )
         {
            ++j;
            continue;
         }

         SCIPdebugMessage("executing presolving of propagator <%s>\n", SCIPpropGetName(scip->set->props[j]));
         SCIP_CALL( SCIPpropPresol(scip->set->props[j], scip->set, onlydelayed, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs, 
               &scip->stat->npresolchgsides, &result) );
         assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);

         lastranpresol = FALSE;
         ++j;
      }
      else
      {
         /* only presolving methods which have non-negative priority will be called before constraint handlers */
         if( priopresol < 0 )
            break;

         if( onlydelayed && !SCIPpresolWasDelayed(scip->set->presols[i]) )
         {
            ++i;
            continue;
         }

         SCIPdebugMessage("executing presolver <%s>\n", SCIPpresolGetName(scip->set->presols[i]));
         SCIP_CALL( SCIPpresolExec(scip->set->presols[i], scip->set, onlydelayed, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs, 
               &scip->stat->npresolchgsides, &result) );
         assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);

         lastranpresol = TRUE;
         ++i;
      }

      if( result == SCIP_CUTOFF )
      {
         *infeasible = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected infeasibility\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected infeasibility\n", SCIPpropGetName(scip->set->props[j-1]));
      }
      else if( result == SCIP_UNBOUNDED )
      {
         *unbounded = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected unboundedness (or infeasibility)\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected  unboundedness (or infeasibility)\n", SCIPpropGetName(scip->set->props[j-1]));
      }
      *delayed = *delayed || (result == SCIP_DELAYED);

      /* delete the variables from the problems that were marked to be deleted */
      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp,
            scip->branchcand) );

      /* if we work off the delayed presolvers, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_SUCCESS )
      {
         *delayed = TRUE;
         aborted = TRUE;
      }
   }

   /* call presolve methods of constraint handlers */
   for( i = 0; i < scip->set->nconshdlrs && !(*unbounded) && !(*infeasible) && !aborted; ++i )
   {
      if( onlydelayed && !SCIPconshdlrWasPresolvingDelayed(scip->set->conshdlrs[i]) )
         continue;

      SCIPdebugMessage("executing presolve method of constraint handler <%s>\n",
         SCIPconshdlrGetName(scip->set->conshdlrs[i]));
      SCIP_CALL( SCIPconshdlrPresolve(scip->set->conshdlrs[i], scip->mem->probmem, scip->set, scip->stat,
            onlydelayed, scip->stat->npresolrounds,
            &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
            &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
            &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs, 
            &scip->stat->npresolchgsides, &result) );
      assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);
      if( result == SCIP_CUTOFF )
      {
         *infeasible = TRUE;
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "constraint handler <%s> detected infeasibility\n", SCIPconshdlrGetName(scip->set->conshdlrs[i]));
      }
      else if( result == SCIP_UNBOUNDED )
      {
         *unbounded = TRUE;
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "constraint handler <%s> detected unboundedness (or infeasibility)\n",
            SCIPconshdlrGetName(scip->set->conshdlrs[i]));
      }
      *delayed = *delayed || (result == SCIP_DELAYED);

      /* delete the variables from the problems that were marked to be deleted */
      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp,
            scip->branchcand) );

      /* if we work off the delayed presolvers, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_SUCCESS )
      {
         *delayed = TRUE;
         aborted = TRUE;
      }
   }

   /* call included presolvers with negative priority */
   for( i = 0, j = 0; !(*unbounded) && !(*infeasible) && !aborted && (i < scip->set->npresols || j < scip->set->nprops);  )
   {
      if( i < scip->set->npresols )
         priopresol = SCIPpresolGetPriority(scip->set->presols[i]);
      else
         priopresol = -INT_MAX;

      if( j < scip->set->nprops )
         prioprop = SCIPpropGetPresolPriority(scip->set->props[j]);
      else
         prioprop = -INT_MAX;

      /* choose presolving */
      if( prioprop >= priopresol )
      {
         /* only presolving methods which have negative priority will be called after constraint handlers */
         if( prioprop >= 0 || (onlydelayed && !SCIPpropWasPresolDelayed(scip->set->props[j])) )
         {
            ++j;
            continue;
         }

         SCIPdebugMessage("executing presolving of propagator <%s>\n", SCIPpropGetName(scip->set->props[j]));
         SCIP_CALL( SCIPpropPresol(scip->set->props[j], scip->set, onlydelayed, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs, 
               &scip->stat->npresolchgsides, &result) );
         assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);

         lastranpresol = FALSE;
         ++j;
      }
      else
      {
         /* only presolving methods which have negative priority will be called after constraint handlers */
         if( priopresol >= 0 || (onlydelayed && !SCIPpresolWasDelayed(scip->set->presols[i])) )
         {
            ++i;
            continue;
         }

         SCIPdebugMessage("executing presolver <%s>\n", SCIPpresolGetName(scip->set->presols[i]));
         SCIP_CALL( SCIPpresolExec(scip->set->presols[i], scip->set, onlydelayed, scip->stat->npresolrounds,
               &scip->stat->npresolfixedvars, &scip->stat->npresolaggrvars, &scip->stat->npresolchgvartypes,
               &scip->stat->npresolchgbds, &scip->stat->npresoladdholes, &scip->stat->npresoldelconss,
               &scip->stat->npresoladdconss, &scip->stat->npresolupgdconss, &scip->stat->npresolchgcoefs, 
               &scip->stat->npresolchgsides, &result) );
         assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);

         lastranpresol = TRUE;
         ++i;
      }

      if( result == SCIP_CUTOFF )
      {
         *infeasible = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected infeasibility\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected infeasibility\n", SCIPpropGetName(scip->set->props[j-1]));
      }
      else if( result == SCIP_UNBOUNDED )
      {
         *unbounded = TRUE;

         if( lastranpresol )
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "presolver <%s> detected unboundedness (or infeasibility)\n", SCIPpresolGetName(scip->set->presols[i-1]));
         else
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
               "propagator <%s> detected  unboundedness (or infeasibility)\n", SCIPpropGetName(scip->set->props[j-1]));
      }
      *delayed = *delayed || (result == SCIP_DELAYED);

      /* delete the variables from the problems that were marked to be deleted */
      SCIP_CALL( SCIPprobPerformVarDeletions(scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp,
            scip->branchcand) );

      /* if we work off the delayed presolvers, we stop immediately if a reduction was found */
      if( onlydelayed && result == SCIP_SUCCESS )
      {
         *delayed = TRUE;
         aborted = TRUE;
      }
   }

   /* remove empty and single variable cliques from the clique table, and convert all two variable cliques
    * into implications
    */
   if( !(*unbounded) && !(*infeasible) )
   {
      SCIP_Bool infeas;

      SCIP_CALL( SCIPcliquetableCleanup(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, &infeas) );
      if( infeas )
      {
         *infeasible = TRUE;
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "clique table cleanup detected infeasibility\n");
      }
      else if( scip->set->nheurs > 0 )
      {
         /* call primal heuristics that are applicable during presolving */
         SCIP_Bool foundsol;
         
         SCIPdebugMessage("calling primal heuristics during presolving\n");
         
         /* call primal heuristics */
         SCIP_CALL( SCIPprimalHeuristics(scip->set, scip->stat, scip->primal, NULL, NULL, NULL, 
               SCIP_HEURTIMING_DURINGPRESOLLOOP, &foundsol) );
         
         /* output a message, if a solution was found */
         if( foundsol )
         {
            SCIP_SOL* sol;
            
            assert(SCIPgetNSols(scip) > 0);         
            sol = SCIPgetBestSol(scip);
            assert(sol != NULL);           
            assert(SCIPgetSolOrigObj(scip,sol) != SCIP_INVALID); /*lint !e777*/
            
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
               "feasible solution found by %s heuristic, objective value %13.6e\n",
               SCIPheurGetName(SCIPsolGetHeur(sol)), SCIPgetSolOrigObj(scip,sol));                    
         }
      }
   }
   
   /* issue PRESOLVEROUND event */
   SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_PRESOLVEROUND) );
   SCIP_CALL( SCIPeventProcess(&event, scip->set, NULL, NULL, NULL, scip->eventfilter) );

   return SCIP_OKAY;
}

/** loops through the included presolvers and constraint's presolve methods, until changes are too few */
static
SCIP_RETCODE presolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            unbounded,          /**< pointer to store whether presolving detected unboundedness */
   SCIP_Bool*            infeasible          /**< pointer to store whether presolving detected infeasibility */
   )
{
   SCIP_Bool delayed;
   SCIP_Bool finished;
   SCIP_Bool stopped;
   SCIP_Real abortfac;
   int maxnrounds;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->primal != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->transprob != NULL);
   assert(scip->set->stage == SCIP_STAGE_TRANSFORMED || scip->set->stage == SCIP_STAGE_PRESOLVING);
   assert(unbounded != NULL);
   assert(infeasible != NULL);

   *unbounded = FALSE;
   *infeasible = FALSE;

   /* switch status to unknown */
   scip->stat->status = SCIP_STATUS_UNKNOWN;

   /* update upper bound and cutoff bound due to objective limit in primal data */
   SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
         scip->transprob, scip->tree, scip->lp) );

   /* start presolving timer */
   SCIPclockStart(scip->stat->presolvingtime, scip->set);

   /* initialize presolving */
   if( scip->set->stage == SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( initPresolve(scip, unbounded, infeasible) );
      if( *infeasible )
      {
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "presolve initialization detected infeasibility\n");
      }
      else if( *unbounded )
      {
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "presolve initialization detected unboundedness\n");
      }
   }
   assert(scip->set->stage == SCIP_STAGE_PRESOLVING);
   
   /* call primal heuristics that are applicable before presolving */
   if( !(*infeasible) && !(*unbounded) && scip->set->nheurs > 0 )
   {
      SCIP_Bool foundsol;

      SCIPdebugMessage("calling primal heuristics before presolving\n");

      /* call primal heuristics */
      SCIP_CALL( SCIPprimalHeuristics(scip->set, scip->stat, scip->primal, NULL, NULL, NULL, SCIP_HEURTIMING_BEFOREPRESOL, &foundsol) );

      /* output a message, if a solution was found */
      if( foundsol )
      {
         SCIP_SOL* sol;

         assert(SCIPgetNSols(scip) > 0);         
         sol = SCIPgetBestSol(scip);
         assert(sol != NULL);           
         assert(SCIPgetSolOrigObj(scip,sol) != SCIP_INVALID);  /*lint !e777*/
         
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "feasible solution found by %s heuristic, objective value %13.6e\n",
            SCIPheurGetName(SCIPsolGetHeur(sol)), SCIPgetSolOrigObj(scip,sol));                    
      }
   }

   maxnrounds = scip->set->presol_maxrounds;
   if( maxnrounds == -1 )
      maxnrounds = INT_MAX;

   abortfac = scip->set->presol_abortfac;

   SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "presolving:\n");

   finished = (*unbounded || *infeasible || scip->stat->npresolrounds >= maxnrounds);
   stopped = SCIPsolveIsStopped(scip->set, scip->stat, TRUE);

   /* perform presolving rounds */
   while( !finished && !stopped )
   {
      int lastnfixedvars;
      int lastnaggrvars;
      int lastnchgvartypes;
      int lastnchgbds;
      int lastnaddholes;
      int lastndelconss;
      int lastnaddconss;
      int lastnupgdconss;
      int lastnchgcoefs;
      int lastnchgsides;
      /*int lastnimplications;*/
      /*int lastncliques;*/


      lastnfixedvars = scip->stat->npresolfixedvars;
      lastnaggrvars = scip->stat->npresolaggrvars;
      lastnchgvartypes = scip->stat->npresolchgvartypes;
      lastnchgbds = scip->stat->npresolchgbds;
      lastnaddholes = scip->stat->npresoladdholes;
      lastndelconss = scip->stat->npresoldelconss;
      lastnaddconss = scip->stat->npresoladdconss;
      lastnupgdconss = scip->stat->npresolupgdconss;
      lastnchgcoefs = scip->stat->npresolchgcoefs;
      lastnchgsides = scip->stat->npresolchgsides;
      /*lastnimplications = scip->stat->nimplications;*/
      /*lastncliques = SCIPcliquetableGetNCliques(scip->cliquetable);*/

      /* sort propagators */
      SCIPsetSortPropsPresol(scip->set);

      /* sort presolvers by priority */
      SCIPsetSortPresols(scip->set);

      /* perform the presolving round by calling the presolvers and constraint handlers */
      assert(!(*unbounded));
      assert(!(*infeasible));
      SCIP_CALL( presolveRound(scip, FALSE, &delayed, unbounded, infeasible) );

      /* check, if we should abort presolving due to not enough changes in the last round */
      finished = isPresolveFinished(scip, abortfac, maxnrounds, lastnfixedvars, lastnaggrvars, lastnchgvartypes,
         lastnchgbds, lastnaddholes, lastndelconss, lastnaddconss, lastnupgdconss, lastnchgcoefs, lastnchgsides,
         /*lastnimplications, lastncliques,*/ *unbounded, *infeasible);

      /* if the presolving will be terminated, call the delayed presolvers */
      while( delayed && finished && !(*unbounded) && !(*infeasible) )
      {
         /* call the delayed presolvers and constraint handlers */
         SCIP_CALL( presolveRound(scip, TRUE, &delayed, unbounded, infeasible) );

         /* check again, if we should abort presolving due to not enough changes in the last round */
         finished = isPresolveFinished(scip, abortfac, maxnrounds, lastnfixedvars, lastnaggrvars, lastnchgvartypes,
            lastnchgbds, lastnaddholes, lastndelconss, lastnaddconss, lastnupgdconss, lastnchgcoefs, lastnchgsides,
            /*lastnimplications, lastncliques,*/ *unbounded, *infeasible);
      }

      /* increase round number */
      scip->stat->npresolrounds++;

      if( !finished )
      {
         /* print presolving statistics */
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
            "(round %d) %d del vars, %d del conss, %d add conss, %d chg bounds, %d chg sides, %d chg coeffs, %d upgd conss, %d impls, %d clqs\n",
            scip->stat->npresolrounds, scip->stat->npresolfixedvars + scip->stat->npresolaggrvars,
            scip->stat->npresoldelconss, scip->stat->npresoladdconss, 
            scip->stat->npresolchgbds, scip->stat->npresolchgsides,
            scip->stat->npresolchgcoefs, scip->stat->npresolupgdconss,
            scip->stat->nimplications, SCIPcliquetableGetNCliques(scip->cliquetable));
      }

      /* abort if time limit was reached or user interrupted */
      stopped = SCIPsolveIsStopped(scip->set, scip->stat, TRUE);
   }

   /* deinitialize presolving */
   if( finished )
   {
      SCIP_Bool unbd;
      SCIP_Bool infeas;

#if 0
      /* @todo the idea is to resort the variables w.r.t. to their index; this means to sort the variables in the same
       * order within their categories/region as they were generated in the beginning; the hope was to get a more
       * deterministic behavior; computational experiments showed that this resorting slows down the overall solution
       * process by 10%; therefore, this code is not used and more tests are needed to overcome this draw back of 10% */
      
      SCIP_VAR** vars;
      int nbinvars;
      int nintvars;
      int nimplvars;
      int ncontvars;
      int nvars;
      int i;

      /* sort variables, which appear in four categories (binary, integer, implicit, continuous) after presolve, again
       * with respect to their original index and in their categories. Adjust the problem index afterwards which is
       * supposed to reflect the position in the variable array. This additional sorting is supposed to reobtain a
       * possible block structure induced by the user model */
      vars = scip->transprob->vars;
      nvars = scip->transprob->nvars;
      nbinvars = scip->transprob->nbinvars;
      nintvars = scip->transprob->nintvars;
      nimplvars = scip->transprob->nimplvars;
      ncontvars = scip->transprob->ncontvars;

      assert(vars != NULL);
      assert(nbinvars + nintvars + nimplvars + ncontvars == nvars);

      SCIPdebugMessage("entering sorting with respect to original block structure! \n");
     
      /* check for every variable type whether it appears in the presolved problem. if yes, sort the specific part of
       * the variables array */

      /* sort binaries */
      if( nbinvars > 0 )
         SCIPsortPtr((void**)vars, SCIPvarComp, nbinvars);
      
      /* sort integers */
      if( nintvars > 0 )
         SCIPsortPtr((void**)&vars[nbinvars], SCIPvarComp, nintvars);
      
      /* sort implicit variables  */
      if( nimplvars > 0 )
         SCIPsortPtr((void**)&vars[nbinvars + nintvars], SCIPvarComp, nimplvars);
      
      /* sort continuous variables*/
      if( ncontvars > 0 )
         SCIPsortPtr((void**)&vars[nbinvars + nintvars + nimplvars], SCIPvarComp, ncontvars);

      /* after sorting, the problem index of each variable has to be adjusted */
      for( i = 0; i < nvars; i++ )
      {
         vars[i]->probindex = i;
         SCIPdebugMessage("Variable: Problem index <%d>, original index <%d> \n", vars[i]->probindex, vars[i]->index);
      }
#endif
      
      SCIP_CALL( exitPresolve(scip, &unbd, &infeas, *unbounded, *infeasible) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);
      if( infeas && !(*infeasible) )
      {
         *infeasible = TRUE;
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "presolve deinitialization detected infeasibility\n");
      }
      else if( unbd && !(*infeasible) && !(*unbounded) )
      {
         *unbounded = TRUE;
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_FULL,
            "presolve deinitialization detected unboundedness\n");
      }
   }
   assert(SCIPbufferGetNUsed(scip->set->buffer) == 0);

   /* stop presolving time */
   SCIPclockStop(scip->stat->presolvingtime, scip->set);

   /* print presolving statistics */
   SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      "presolving (%d rounds):\n", scip->stat->npresolrounds);
   SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      " %d deleted vars, %d deleted constraints, %d added constraints, %d tightened bounds, %d added holes, %d changed sides, %d changed coefficients\n",
      scip->stat->npresolfixedvars + scip->stat->npresolaggrvars, scip->stat->npresoldelconss, scip->stat->npresoladdconss,
      scip->stat->npresolchgbds, scip->stat->npresoladdholes, scip->stat->npresolchgsides, scip->stat->npresolchgcoefs);
   SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
      " %d implications, %d cliques\n", scip->stat->nimplications, SCIPcliquetableGetNCliques(scip->cliquetable));

   /* remember number of constraints */
   SCIPprobMarkNConss(scip->transprob);

   return SCIP_OKAY;
}

/** initializes solution process data structures */
static
SCIP_RETCODE initSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->nlp == NULL);
   assert(scip->set->stage == SCIP_STAGE_PRESOLVED);

   /* reset statistics for current branch and bound run */
   SCIPstatResetCurrentRun(scip->stat);
   SCIPstatEnforceLPUpdates(scip->stat);

   /* LP is empty anyway; mark empty LP to be solved and update validsollp counter */
   SCIP_CALL( SCIPlpReset(scip->lp, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter) );

   /* update upper bound and cutoff bound due to objective limit in primal data */
   SCIP_CALL( SCIPprimalUpdateObjlimit(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
         scip->transprob, scip->tree, scip->lp) );

   /* switch stage to INITSOLVE */
   scip->set->stage = SCIP_STAGE_INITSOLVE;

   /* initialize NLP if there are nonlinearities and there is someone who can make use of it and the user did not switch it off */
   assert(!scip->set->continnonlinpresent || scip->set->nonlinearitypresent); /* if there is continuous nonlinearity, then there should be nonlinearity in general */
   if( scip->set->nonlinearitypresent && !scip->set->nlp_disable )
   {
      SCIPdebugMessage("constructing empty NLP\n");

      SCIP_CALL( SCIPnlpCreate(&scip->nlp, scip->mem->probmem, scip->set, scip->stat, SCIPprobGetName(scip->transprob), scip->transprob->nvars) );
      assert(scip->nlp != NULL);

      SCIP_CALL( SCIPnlpAddVars(scip->nlp, scip->mem->probmem, scip->set, scip->transprob->nvars, scip->transprob->vars) );
   }

   /* create VBC output file */
   SCIP_CALL( SCIPvbcInit(scip->stat->vbc, scip->mem->probmem, scip->set) );

   /* initialize solution process data structures */
   SCIP_CALL( SCIPpricestoreCreate(&scip->pricestore) );
   SCIP_CALL( SCIPsepastoreCreate(&scip->sepastore) );
   SCIP_CALL( SCIPcutpoolCreate(&scip->cutpool, scip->mem->probmem, scip->set, scip->set->sepa_cutagelimit, TRUE) );
   SCIP_CALL( SCIPtreeCreateRoot(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );

   /* switch stage to SOLVING */
   scip->set->stage = SCIP_STAGE_SOLVING;

   /* inform the transformed problem that the branch and bound process starts now */
   SCIP_CALL( SCIPprobInitSolve(scip->transprob, scip->set) );

   /* inform plugins that the branch and bound process starts now */
   SCIP_CALL( SCIPsetInitsolPlugins(scip->set, scip->mem->probmem, scip->stat) );

   /* remember number of constraints */
   SCIPprobMarkNConss(scip->transprob);

   /* if all variables are known, calculate a trivial primal bound by setting all variables to their worst bound */
   if( scip->set->nactivepricers == 0 )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_Real objbound;
      SCIP_Real bd;
      int v;

      objbound = 0.0;
      for( v = 0; v < scip->transprob->nvars && !SCIPsetIsInfinity(scip->set, objbound); ++v )
      {
         var = scip->transprob->vars[v];
         obj = SCIPvarGetObj(var);
         if( !SCIPsetIsZero(scip->set, obj) )
         {
            bd = SCIPvarGetWorstBound(var);
            if( SCIPsetIsInfinity(scip->set, REALABS(bd)) )
               objbound = SCIPsetInfinity(scip->set);
            else
               objbound += obj * bd;
         }
      }

      /* update primal bound (add 1.0 to primal bound, such that solution with worst bound may be found) */
      if( !SCIPsetIsInfinity(scip->set, objbound) && SCIPsetIsLT(scip->set, objbound + 1.0, scip->primal->cutoffbound) )
      {
         SCIP_CALL( SCIPprimalSetCutoffbound(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
               scip->tree, scip->lp, objbound + 1.0) );
      }
   }

   return SCIP_OKAY;
}

/** frees solution process data structures */
static
SCIP_RETCODE freeSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             restart             /**< was this free solve call triggered by a restart? */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);
   assert(scip->stat != NULL);
   assert(scip->set->stage == SCIP_STAGE_SOLVING || scip->set->stage == SCIP_STAGE_SOLVED);

   /* mark that we are currently restarting */
   if( restart )
      scip->stat->inrestart = TRUE;

   /* remove focus from the current focus node */
   if( SCIPtreeGetFocusNode(scip->tree) != NULL )
   {
      SCIP_NODE* node = NULL;
      SCIP_Bool cutoff;

      SCIP_CALL( SCIPnodeFocus(&node, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
            scip->tree, scip->lp, scip->branchcand, scip->conflict, scip->eventfilter, scip->eventqueue, &cutoff) );
      assert(!cutoff);
   }

   /* switch stage to FREESOLVE */
   scip->set->stage = SCIP_STAGE_FREESOLVE;

   /* inform plugins that the branch and bound process is finished */
   SCIP_CALL( SCIPsetExitsolPlugins(scip->set, scip->mem->probmem, scip->stat, restart) );

   /* free the NLP, if there is one, and reset the flags indicating nonlinearity */
   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpFree(&scip->nlp, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp) );
   }
   scip->set->continnonlinpresent = FALSE;
   scip->set->nonlinearitypresent = FALSE;
   
   /* clear the LP, and flush the changes to clear the LP of the solver */
   SCIP_CALL( SCIPlpReset(scip->lp, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter) );
   SCIPlpInvalidateRootObjval(scip->lp);

   /* clear all row references in internal data structures */
   SCIP_CALL( SCIPcutpoolClear(scip->cutpool, scip->mem->probmem, scip->set, scip->lp) );

   /* we have to clear the tree prior to the problem deinitialization, because the rows stored in the forks and
    * subroots have to be released
    */
   SCIP_CALL( SCIPtreeClear(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );

   /* deinitialize transformed problem */
   SCIP_CALL( SCIPprobExitSolve(scip->transprob, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, restart) );

   /* free solution process data structures */
   SCIP_CALL( SCIPcutpoolFree(&scip->cutpool, scip->mem->probmem, scip->set, scip->lp) );
   SCIP_CALL( SCIPsepastoreFree(&scip->sepastore) );
   SCIP_CALL( SCIPpricestoreFree(&scip->pricestore) );

   /* close VBC output file */
   SCIPvbcExit(scip->stat->vbc, scip->set);

   /* reset statistics for current branch and bound run */
   SCIPstatResetCurrentRun(scip->stat);

   /* switch stage to TRANSFORMED */
   scip->set->stage = SCIP_STAGE_TRANSFORMED;

   /* restart finished */
   assert( ! restart || scip->stat->inrestart );
   scip->stat->inrestart = FALSE;

   return SCIP_OKAY;
}

/** free transformed problem */
static
SCIP_RETCODE freeTransform(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool stored;
   int nsols;
   int s;
   
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->stat != NULL);
   assert(scip->set->stage == SCIP_STAGE_TRANSFORMED || scip->set->stage == SCIP_STAGE_PRESOLVING);

   /* call exit methods of plugins */
   SCIP_CALL( SCIPsetExitPlugins(scip->set, scip->mem->probmem, scip->stat) );

   /* switch stage to FREETRANS */
   scip->set->stage = SCIP_STAGE_FREETRANS;
   
   assert(scip->origprimal->nsols == 0);

   nsols = MIN(scip->set->limit_maxorigsol, scip->primal->nsols);
   stored = TRUE;
   
   /* copy best primal solution to original solution candidate list */
   for( s = 0; s < nsols && stored; ++s )
   {
      SCIP_SOL* sol;

      sol = scip->primal->sols[s];
      assert(sol != NULL);
      
      if( SCIPsolGetOrigin(sol) != SCIP_SOLORIGIN_ORIGINAL )
      {
         /* retransform solution into the original problem space */
         SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob) );
      }
   
      /* add solution to original candidate solution storage */
      SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
   }
   
   /* free transformed problem data structures */
   SCIP_CALL( SCIPprobFree(&scip->transprob, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
   SCIP_CALL( SCIPcliquetableFree(&scip->cliquetable, scip->mem->probmem) );
   SCIP_CALL( SCIPconflictFree(&scip->conflict, scip->mem->probmem) );
   SCIP_CALL( SCIPtreeFree(&scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
   SCIP_CALL( SCIPprimalFree(&scip->primal, scip->mem->probmem) );
   SCIP_CALL( SCIPrelaxationFree(&scip->relaxation) );
   SCIP_CALL( SCIPlpFree(&scip->lp, scip->mem->probmem, scip->set, scip->eventqueue, scip->eventfilter) );
   SCIP_CALL( SCIPbranchcandFree(&scip->branchcand) );
   SCIP_CALL( SCIPeventfilterFree(&scip->eventfilter, scip->mem->probmem, scip->set) );
   SCIP_CALL( SCIPeventqueueFree(&scip->eventqueue) );

   /* free all debug data */ 
   SCIP_CALL( SCIPdebugFreeDebugData(scip->set) );

   if( scip->set->misc_resetstat )
   {
      /* reset statistics to the point before the problem was transformed */
      SCIPstatReset(scip->stat);
   }

   /* switch stage to PROBLEM */
   scip->set->stage = SCIP_STAGE_PROBLEM;

   /* reset original variable's local and global bounds to their original values */
   SCIP_CALL( SCIPprobResetBounds(scip->origprob, scip->mem->probmem, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** transforms and presolves problem */
SCIP_RETCODE SCIPpresolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool unbounded;
   SCIP_Bool infeasible;

   SCIP_CALL( checkStage(scip, "SCIPpresolve", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* start solving timer */
   SCIPclockStart(scip->stat->solvingtime, scip->set);

   /* capture the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptCapture(scip->interrupt);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* initialize solving data structures and transform problem */
      SCIP_CALL( SCIPtransformProb(scip) );
      assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);

      /*lint -fallthrough*/

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* presolve problem */
      SCIP_CALL( presolve(scip, &unbounded, &infeasible) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED || scip->set->stage == SCIP_STAGE_PRESOLVING);

      if( scip->set->stage == SCIP_STAGE_PRESOLVED )
      {
         if( infeasible || unbounded )
         {
	    /* first change status of scip, so that all plugins in their initsol callbacks can ask SCIP for the correct
	     * status
	     */
            if( infeasible )
            {
	       /* switch status to OPTIMAL */
               if( scip->primal->nsols > 0
                  && SCIPsetIsLT(scip->set, SCIPsolGetObj(scip->primal->sols[0], scip->set, scip->transprob),
                     SCIPprobInternObjval(scip->transprob, scip->set, SCIPprobGetObjlim(scip->transprob, scip->set))) )
                  scip->stat->status = SCIP_STATUS_OPTIMAL;
               else /* switch status to INFEASIBLE */
                  scip->stat->status = SCIP_STATUS_INFEASIBLE;
            }
            else if( scip->primal->nsols >= 1 ) /* switch status to UNBOUNDED */
               scip->stat->status = SCIP_STATUS_UNBOUNDED;
            else /* switch status to INFORUNBD */
               scip->stat->status = SCIP_STATUS_INFORUNBD;

            /* initialize solving process data structures to be able to switch to SOLVED stage */
            SCIP_CALL( initSolve(scip) );

            /* switch stage to SOLVED */
            scip->set->stage = SCIP_STAGE_SOLVED;

            /* print solution message */
            if( infeasible )
            {
               /* infeasibility in this round means, that the current best solution is optimal (if existing and not
                * worse than user objective limit)
                */
               if( scip->primal->nsols > 0
                  && SCIPsetIsLT(scip->set, SCIPsolGetObj(scip->primal->sols[0], scip->set, scip->transprob),
                     SCIPprobInternObjval(scip->transprob, scip->set, SCIPprobGetObjlim(scip->transprob, scip->set))) )
               {
                  /* switch status to OPTIMAL */
                  scip->stat->status = SCIP_STATUS_OPTIMAL;

                  /* remove the root node from the tree, s.t. the lower bound is set to +infinity */
                  SCIP_CALL( SCIPtreeClear(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->lp) );
               }
               else
               {
                  SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
                     "presolving detected infeasibility\n");
                  
                  /* switch status to INFEASIBLE */
                  scip->stat->status = SCIP_STATUS_INFEASIBLE;
               }
            }
            else if( scip->primal->nsols >= 1 )
            {
               SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
                  "presolving detected unboundedness\n");

               /* switch status to UNBOUNDED */
               scip->stat->status = SCIP_STATUS_UNBOUNDED;
            }
            else
            {
               SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
                  "presolving detected unboundedness (or infeasibility)\n");

               /* switch status to INFORUNBD */
               scip->stat->status = SCIP_STATUS_INFORUNBD;
            }
         }
         else
         {
            int h;

            /* print presolved problem statistics */
            SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL,
               "presolved problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints\n",
               scip->transprob->nvars, scip->transprob->nbinvars, scip->transprob->nintvars, scip->transprob->nimplvars,
               scip->transprob->ncontvars, scip->transprob->nconss);

            for( h = 0; h < scip->set->nconshdlrs; ++h )
            {
               int nactiveconss;

               nactiveconss = SCIPconshdlrGetNActiveConss(scip->set->conshdlrs[h]);
               if( nactiveconss > 0 )
               {
                  SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
                     "%7d constraints of type <%s>\n", nactiveconss, SCIPconshdlrGetName(scip->set->conshdlrs[h]));
               }
            }

            if( SCIPprobIsObjIntegral(scip->transprob) )
            {
               SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
                  "transformed objective value is always integral (scale: %.15g)\n", scip->transprob->objscale);
            }
         }
      }
      else
      {
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH, "presolving was interrupted.\n");
      }

      /* display timing statistics */
      SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_HIGH,
         "Presolving Time: %.2f\n", SCIPclockGetTime(scip->stat->presolvingtime));
      break;

   case SCIP_STAGE_PRESOLVED:
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   /* release the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptRelease(scip->interrupt);

   /* stop solving timer */
   SCIPclockStop(scip->stat->solvingtime, scip->set);

   return SCIP_OKAY;
}


/** transforms, presolves, and solves problem */
SCIP_RETCODE SCIPsolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool restart;

   SCIP_CALL( checkStage(scip, "SCIPsolve", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   /* check, if a node selector exists */
   if( SCIPsetGetNodesel(scip->set, scip->stat) == NULL )
   {
      SCIPerrorMessage("no node selector available\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* start solving timer */
   SCIPclockStart(scip->stat->solvingtime, scip->set);

   /* capture the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptCapture(scip->interrupt);

   /* automatic restarting loop */
   restart = scip->stat->userrestart;
   do
   {
      if( restart )
      {
         /* free the solving process data in order to restart */
         assert(scip->set->stage == SCIP_STAGE_SOLVING);
         if( scip->stat->userrestart )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "(run %d, node %lld) performing user restart\n",
               scip->stat->nruns, scip->stat->nnodes, scip->stat->nrootintfixingsrun);
         else
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
               "(run %d, node %lld) restarting after %d global fixings of integer variables\n",
               scip->stat->nruns, scip->stat->nnodes, scip->stat->nrootintfixingsrun);
         /* an extra blank line should be printed separately since the buffer message handler only handles up to one line
          * correctly */
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\n");
         SCIP_CALL( freeSolve(scip, TRUE) );
         assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);
      }
      restart = FALSE;
      scip->stat->userrestart = FALSE;
      
      switch( scip->set->stage )
      {
      case SCIP_STAGE_PROBLEM:
      case SCIP_STAGE_TRANSFORMED:
      case SCIP_STAGE_PRESOLVING:
         /* initialize solving data structures, transform and problem */
         SCIP_CALL( SCIPpresolve(scip) );
         if( scip->set->stage == SCIP_STAGE_SOLVED || scip->set->stage == SCIP_STAGE_PRESOLVING )
            break;
         assert(scip->set->stage == SCIP_STAGE_PRESOLVED);

         /*lint -fallthrough*/

      case SCIP_STAGE_PRESOLVED:
         /* initialize solving process data structures */
         SCIP_CALL( initSolve(scip) );
         assert(scip->set->stage == SCIP_STAGE_SOLVING);
         SCIPmessagePrintVerbInfo(scip->set->disp_verblevel, SCIP_VERBLEVEL_NORMAL, "\n");

         /*lint -fallthrough*/

      case SCIP_STAGE_SOLVING:
         /* reset display */
         SCIPstatResetDisplay(scip->stat);

         /* continue solution process */
         SCIP_CALL( SCIPsolveCIP(scip->mem->probmem, scip->set, scip->stat, scip->mem, scip->origprob, scip->transprob,
               scip->primal, scip->tree, scip->lp, scip->relaxation, scip->pricestore, scip->sepastore, 
               scip->cutpool, scip->branchcand, scip->conflict, scip->eventfilter, scip->eventqueue, &restart) );

         /* detect, whether problem is solved */
         if( SCIPtreeGetNNodes(scip->tree) == 0 && SCIPtreeGetCurrentNode(scip->tree) == NULL )
         {
            assert(scip->stat->status == SCIP_STATUS_OPTIMAL
               || scip->stat->status == SCIP_STATUS_INFEASIBLE
               || scip->stat->status == SCIP_STATUS_UNBOUNDED
               || scip->stat->status == SCIP_STATUS_INFORUNBD);
            assert(!restart);

            /* tree is empty, and no current node exists -> problem is solved */
            scip->set->stage = SCIP_STAGE_SOLVED;
         }
         break;

      case SCIP_STAGE_SOLVED:
         assert(scip->stat->status == SCIP_STATUS_OPTIMAL
            || scip->stat->status == SCIP_STATUS_INFEASIBLE
            || scip->stat->status == SCIP_STATUS_UNBOUNDED
            || scip->stat->status == SCIP_STATUS_INFORUNBD);
         break;

      default:
         SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
         return SCIP_ERROR;
      }  /*lint !e788*/
   }
   while( restart && !SCIPsolveIsStopped(scip->set, scip->stat, TRUE) 
      && (scip->set->limit_restarts == -1 || scip->stat->nruns <= scip->set->limit_restarts ) );
      
   /* release the CTRL-C interrupt */
   if( scip->set->misc_catchctrlc )
      SCIPinterruptRelease(scip->interrupt);

   /* stop solving timer */
   SCIPclockStop(scip->stat->solvingtime, scip->set);

   /* display most relevant statistics */
   if( scip->set->disp_verblevel >= SCIP_VERBLEVEL_NORMAL )
   {
      SCIPmessagePrintInfo("\n");
      SCIPmessagePrintInfo("SCIP Status        : ");
      SCIP_CALL( SCIPprintStage(scip, NULL) );
      SCIPmessagePrintInfo("\n");
      SCIPmessagePrintInfo("Solving Time (sec) : %.2f\n", SCIPclockGetTime(scip->stat->solvingtime));
      if( scip->stat->nruns > 1 )
         SCIPmessagePrintInfo("Solving Nodes      : %"SCIP_LONGINT_FORMAT" (total of %"SCIP_LONGINT_FORMAT" nodes in %d runs)\n",
            scip->stat->nnodes, scip->stat->ntotalnodes, scip->stat->nruns);
      else
         SCIPmessagePrintInfo("Solving Nodes      : %"SCIP_LONGINT_FORMAT"\n", scip->stat->nnodes);
      if( scip->set->stage >= SCIP_STAGE_TRANSFORMED && scip->set->stage <= SCIP_STAGE_FREESOLVE )
         SCIPmessagePrintInfo("Primal Bound       : %+.14e (%"SCIP_LONGINT_FORMAT" solutions)\n",
            getPrimalbound(scip), scip->primal->nsolsfound);
      if( scip->set->stage >= SCIP_STAGE_SOLVING && scip->set->stage <= SCIP_STAGE_SOLVED )
      {
         SCIPmessagePrintInfo("Dual Bound         : %+.14e\n", getDualbound(scip));
         SCIPmessagePrintInfo("Gap                : ");
         if( SCIPsetIsInfinity(scip->set, SCIPgetGap(scip)) )
            SCIPmessagePrintInfo("infinite\n");
         else
            SCIPmessagePrintInfo("%.2f %%\n", 100.0*SCIPgetGap(scip));
      }

      /* check solution for feasibility in original problem */
      if( scip->set->stage >= SCIP_STAGE_TRANSFORMED )
      {
         SCIP_SOL* sol;

         sol = SCIPgetBestSol(scip);
         if( sol != NULL )
         {
            SCIP_Bool feasible;
            SCIP_CALL( SCIPcheckSolOrig(scip, sol, &feasible, TRUE, FALSE) );
            
            if( !feasible )
            {
               SCIPmessagePrintInfo("best solution is not feasible in original problem\n");
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** frees branch and bound tree and all solution process data; statistics, presolving data and transformed problem is
 *  preserved
 */
SCIP_RETCODE SCIPfreeSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             restart             /**< should certain data be preserved for improved restarting? */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeSolve", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMED:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   {
      SCIP_Bool unbounded;
      SCIP_Bool infeasible;
      SCIP_Bool isunbounded;
      SCIP_Bool isinfeasible;

      if( scip->stat->status == SCIP_STATUS_INFEASIBLE )
      {
	 isinfeasible = TRUE;
	 isunbounded = FALSE;
      }
      else if( scip->stat->status == SCIP_STATUS_INFORUNBD || scip->stat->status == SCIP_STATUS_UNBOUNDED  )
      {
	 isinfeasible = FALSE;
	 isunbounded = TRUE;
      }
      else
      {
	 isinfeasible = FALSE;
	 isunbounded = FALSE;
      }

      /* exit presolving */
      SCIP_CALL( exitPresolve(scip, &unbounded, &infeasible, isunbounded, isinfeasible) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);
   }
   /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
      /* switch stage to TRANSFORMED */
      scip->set->stage = SCIP_STAGE_TRANSFORMED;
      return SCIP_OKAY;

   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* free solution process data structures */
      SCIP_CALL( freeSolve(scip, restart) );
      assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** frees all solution process data including presolving and transformed problem, only original problem is kept */
SCIP_RETCODE SCIPfreeTransform(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeTransform", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
   {
      SCIP_Bool unbounded;
      SCIP_Bool infeasible;
      SCIP_Bool isunbounded;
      SCIP_Bool isinfeasible;

      if( scip->stat->status == SCIP_STATUS_INFEASIBLE )
      {
	 isinfeasible = TRUE;
	 isunbounded = FALSE;
      }
      else if( scip->stat->status == SCIP_STATUS_INFORUNBD || scip->stat->status == SCIP_STATUS_UNBOUNDED  )
      {
	 isinfeasible = FALSE;
	 isunbounded = TRUE;
      }
      else
      {
	 isinfeasible = FALSE;
	 isunbounded = FALSE;
      }

      /* exit presolving */
      SCIP_CALL( exitPresolve(scip, &unbounded, &infeasible, isunbounded, isinfeasible) );
      assert(scip->set->stage == SCIP_STAGE_PRESOLVED);
   }
   /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      /* free solution process data */
      SCIP_CALL( SCIPfreeSolve(scip, FALSE) );
      assert(scip->set->stage == SCIP_STAGE_TRANSFORMED);

      /*lint -fallthrough*/

   case SCIP_STAGE_TRANSFORMED:
      /* free transformed problem data structures */
      SCIP_CALL( freeTransform(scip) );
      assert(scip->set->stage == SCIP_STAGE_PROBLEM);
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** interrupts solving process as soon as possible (e.g., after the current node has been solved) */
SCIP_RETCODE SCIPinterruptSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPinterruptSolve", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE) );

   /* set the userinterrupt flag */
   scip->stat->userinterrupt = TRUE;

   return SCIP_OKAY;
}

/** restarts solving process as soon as possible (e.g., after the current node has been solved) */
SCIP_RETCODE SCIPrestartSolve(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrestartSolve", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* set the userrestart flag */
   scip->stat->userrestart = TRUE;

   return SCIP_OKAY;
}

/** whether we are in the restarting phase */
SCIP_Bool SCIPisInRestart(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisInRestart", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   /* return the restart status */
   return scip->stat->inrestart;
}


/*
 * variable methods
 */

/** creates and captures problem variable; if variable is of integral type, fractional bounds are automatically rounded;
 *  an integer variable with bounds zero and one is automatically converted into a binary variable;
 *
 *  @warning When doing column generation and the original problem is a maximization problem, notice that SCIP will
 *           transform the problem into a minimization problem by multiplying the objective function by -1.  Thus, the
 *           original objective function value of variables created during the solving process has to be multiplied by
 *           -1, too.
 *
 *  @note the variable gets captured, hence at one point you have to release it using the method SCIPreleaseVar()
 */
SCIP_RETCODE SCIPcreateVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARTYPE          vartype,            /**< type of variable */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable, or NULL */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data, or NULL */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable, or NULL */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(var != NULL);
   assert(lb <= ub);

   SCIP_CALL( checkStage(scip, "SCIPcreateVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPvarCreateOriginal(var, scip->mem->probmem, scip->set, scip->stat,
            name, lb, ub, obj, vartype, initial, removable, vardelorig, vartrans, vardeltrans, varcopy, vardata) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPvarCreateTransformed(var, scip->mem->probmem, scip->set, scip->stat,
            name, lb, ub, obj, vartype, initial, removable, vardelorig, vartrans, vardeltrans, varcopy, vardata) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** outputs the variable name to the file stream */
SCIP_RETCODE SCIPwriteVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR*             var,                /**< variable array to output */
   SCIP_Bool             type                /**< should the variable type be also posted */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   
   SCIP_CALL( checkStage(scip, "SCIPwriteVarName", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   /* print variable name */
   if( SCIPvarIsNegated(var) )
   {
      SCIP_VAR* negatedvar;
      
      SCIP_CALL( SCIPgetNegatedVar(scip, var, &negatedvar) );
      SCIPinfoMessage(scip, file, "<~%s>", SCIPvarGetName(negatedvar));
   }
   else
   {
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(var));
   }

   if( type )
   {
      /* print variable type */
      SCIPinfoMessage(scip, file, "[%c]", 
         SCIPvarGetType(var) == SCIP_VARTYPE_BINARY ? SCIP_VARTYPE_BINARY_CHAR :
         SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER ? SCIP_VARTYPE_INTEGER_CHAR :
         SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT ? SCIP_VARTYPE_IMPLINT_CHAR : SCIP_VARTYPE_CONTINUOUS_CHAR);
   }
   
   return SCIP_OKAY;
}

/** print the given list of variables to output stream separated by the given delimiter character;
 *
 *  i. e. the variables x1, x2, ..., xn with given delimiter ',' are written as: \<x1\>, \<x2\>, ..., \<xn\>;
 *
 *  the method SCIPparseVarsList() can parse such a string
 */
SCIP_RETCODE SCIPwriteVarsList(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR**            vars,               /**< variable array to output */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             type,               /**< should the variable type be also posted */
   char                  delimiter           /**< character which is used for delimitation */
   )
{
   int v;
   
   SCIP_CALL( checkStage(scip, "SCIPwriteVarsList", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   
   for( v = 0; v < nvars; ++v )
   {
      if( v > 0 )
      {
         SCIPinfoMessage(scip, file, "%c ", delimiter);
      }
      
      /* print variable name */
      SCIP_CALL( SCIPwriteVarName(scip, file, vars[v], type) );
   }

   return SCIP_OKAY;
}

/** print the given variables and coefficients as linear sum in the following form 
 *  c1 \<x1\> + c2 \<x2\>   ... + cn \<xn\>
 *  
 *  This string can be parsed by the method SCIPparseVarsLinearsum().
 */
SCIP_RETCODE SCIPwriteVarsLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR**            vars,               /**< variable array to output */
   SCIP_Real*            vals,               /**< array of coefficients or NULL if all coefficients are 1.0 */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             type                /**< should the variable type be also posted */
   )
{
   int v;
   
   SCIP_CALL( checkStage(scip, "SCIPwriteVarsLinearsum", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   for( v = 0; v < nvars; ++v )
   {
      if( vals != NULL )
      {
         if( vals[v] == 1.0 )
         {
            if( v > 0 )
               SCIPinfoMessage(scip, file, " +");
         }
         else if( vals[v] == -1.0 )
            SCIPinfoMessage(scip, file, " -");
         else
            SCIPinfoMessage(scip, file, " %+.15g", vals[v]);
      }
      else if( nvars > 0 )
         SCIPinfoMessage(scip, file, " +");
      
      /* print variable name */
      SCIP_CALL( SCIPwriteVarName(scip, file, vars[v], type) );
   }
   
   return SCIP_OKAY;
}

/** print the given monomials as polynomial in the following form
 *  c1 \<x11\>^e11 \<x12\>^e12 ... \<x1n\>^e1n + c2 \<x21\>^e21 \<x22\>^e22 ... + ... + cn \<xn1\>^en1 ...
 *
 *  This string can be parsed by the method SCIPparseVarsPolynomial().
 */
SCIP_RETCODE SCIPwriteVarsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR***           monomialvars,       /**< arrays with variables for each monomial */
   SCIP_Real**           monomialexps,       /**< arrays with variable exponents, or NULL if always 1.0 */
   SCIP_Real*            monomialcoefs,      /**< array with monomial coefficients */
   int*                  monomialnvars,      /**< array with number of variables for each monomial */
   int                   nmonomials,         /**< number of monomials */
   SCIP_Bool             type                /**< should the variable type be also posted */
   )
{
   int i;
   int v;

   assert(scip != NULL);
   assert(monomialvars  != NULL || nmonomials == 0);
   assert(monomialcoefs != NULL || nmonomials == 0);
   assert(monomialnvars != NULL || nmonomials == 0);

   SCIP_CALL( checkStage(scip, "SCIPwriteVarsPolynomial", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( nmonomials == 0 )
   {
      SCIPinfoMessage(scip, file, " 0 ");
      return SCIP_OKAY;
   }

   for( i = 0; i < nmonomials; ++i )
   {
      if( monomialcoefs[i] == 1.0 )
      {
         if( i > 0 )
            SCIPinfoMessage(scip, file, " +");
      }
      else if( monomialcoefs[i] == -1.0 )
         SCIPinfoMessage(scip, file, " -");
      else
         SCIPinfoMessage(scip, file, " %+.15g", monomialcoefs[i]);

      assert(monomialvars[i] != NULL || monomialnvars[i] == 0);

      for( v = 0; v < monomialnvars[i]; ++v )
      {
         SCIP_CALL( SCIPwriteVarName(scip, file, monomialvars[i][v], type) );
         if( monomialexps != NULL && monomialexps[i] != NULL && monomialexps[i][v] != 1.0 )
         {
            SCIPinfoMessage(scip, file, "^%.15g", monomialexps[i][v]);
         }
      }
   }

   return SCIP_OKAY;
}

/** parses variable information (in cip format) out of a string; if the parsing process was successful a variable is
 *  created and captured; if variable is of integral type, fractional bounds are automatically rounded; an integer
 *  variable with bounds zero and one is automatically converted into a binary variable
 */
SCIP_RETCODE SCIPparseVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to store the problem variable */
   const char*           str,                /**< string to parse */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_DECL_VARCOPY     ((*varcopy)),       /**< copies variable data if wanted to subscip, or NULL */
   SCIP_DECL_VARDELORIG  ((*vardelorig)),    /**< frees user data of original variable */
   SCIP_DECL_VARTRANS    ((*vartrans)),      /**< creates transformed user data by transforming original user data */
   SCIP_DECL_VARDELTRANS ((*vardeltrans)),   /**< frees user data of transformed variable */
   SCIP_VARDATA*         vardata,            /**< user data for this specific variable */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   )
{
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPparseVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPvarParseOriginal(var, scip->mem->probmem, scip->set, scip->stat,
            str, initial, removable, varcopy, vardelorig, vartrans, vardeltrans, vardata, success) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPvarParseTransformed(var, scip->mem->probmem, scip->set, scip->stat,
            str, initial, removable, varcopy, vardelorig, vartrans, vardeltrans, vardata, success) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** parses the given string for a variable name and stores the variable in the corresponding pointer if such a variable
 *  exits and returns the position where the parsing stopped 
 */
SCIP_RETCODE SCIPparseVarName(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            var,                /**< pointer to store the problem variable, or NULL if it does not exit */
   char**                endptr              /**< pointer to store the final string position if successfully */
   )
{
   char varname[SCIP_MAXSTRLEN];
   
   assert(str != NULL);
   assert(var != NULL);
   assert(endptr != NULL);

   SCIP_CALL( checkStage(scip, "SCIPparseVarName", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   SCIPstrCopySection(str, '<', '>', varname, SCIP_MAXSTRLEN, endptr); 
   assert(*endptr != NULL);

   if( *varname == '\0' )
   {
      SCIPerrorMessage("invalid variable name string given: could not find '<'\n");
      return SCIP_INVALIDDATA;
   }
   
   /* check if we have a negated variable */
   if( *varname == '~' )
   {
      SCIPdebugMessage("parsed negated variable name <%s>\n", &varname[1]);

      /* search for the variable and ignore '~' */
      (*var) = SCIPfindVar(scip, &varname[1]);
      
      if( *var  != NULL )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, *var, var) );
      }
   }
   else
   {
      SCIPdebugMessage("parsed variable name <%s>\n", varname);
      
      /* search for the variable */
      (*var) = SCIPfindVar(scip, varname);
   }

   str = *endptr;
   
   /* skip additional variable type marker */
   if( *str == '[' &&
      (str[1] == SCIP_VARTYPE_BINARY_CHAR || str[1] == SCIP_VARTYPE_INTEGER_CHAR || 
         str[1] == SCIP_VARTYPE_IMPLINT_CHAR || str[1] == SCIP_VARTYPE_CONTINUOUS_CHAR )  && str[2] == ']' )
      (*endptr) += 3;
      
   return SCIP_OKAY;
}

/** parse the given string as variable list (here ',' is the delimiter)) (\<x1\>, \<x2\>, ..., \<xn\>) (see
 *  SCIPwriteVarsList() ); if it was successful, the pointer success is set to TRUE
 *
 *  @note the pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
SCIP_RETCODE SCIPparseVarsList(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variable */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   char                  delimiter,          /**< character which is used for delimitation */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successfully or not */
   )
{
   SCIP_VAR** tmpvars;
   SCIP_VAR* var;
   int ntmpvars;
   int v;
   
   SCIP_CALL( checkStage(scip, "SCIPparseVarsList", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* allocate buffer memory for temporary storing the parsed variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvars, varssize) );

   ntmpvars = 0;
   (*success) = TRUE;
   
   do
   {
      *endptr = (char*)str;
      
      /* parse variable name */ 
      SCIP_CALL( SCIPparseVarName(scip, str, &var, endptr) );

      if( var == NULL )
      {
         SCIPdebugMessage("variable with name <%s> does not exist\n", SCIPvarGetName(var));
         (*success) = FALSE;
         break;
      }
      
      /* store the variable in the tmp array */
      if( ntmpvars < varssize )
         tmpvars[ntmpvars] = var;
      
      ntmpvars++;

      str = *endptr;
      
      while( isspace((unsigned char)*str) )
         str++;
   }
   while( *str == delimiter );
   
   *endptr = (char*)str;

   /* if all variable name searches were successfully and the variable array has enough slots copy the collected
    * variables 
    */
   if( (*success) && ntmpvars <= varssize )
   {
      for( v = 0; v < ntmpvars; ++v )
         vars[v] = tmpvars[v];
      
      (*nvars) = ntmpvars;
   }
   else
      (*nvars) = 0;
   
   (*requiredsize) = ntmpvars;
   
   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &tmpvars);

   return SCIP_OKAY;
}

/** parse the given string as linear sum of variables and coefficients (c1 \<x1\> + c2 \<x2\> + ... + cn \<xn\>)
 *  (see SCIPwriteVarsLinearsum() ); if it was successful, the pointer success is set to TRUE
 *
 *  @note the pointer success in only set to FALSE in the case that a variable with a parsed variable name does not exist
 *
 *  @note If the number of (parsed) variables is greater than the available slots in the variable array, nothing happens
 *        except that the required size is stored in the corresponding integer; the reason for this approach is that we
 *        cannot reallocate memory, since we do not know how the memory has been allocated (e.g., by a C++ 'new' or SCIP
 *        memory functions).
 */
SCIP_RETCODE SCIPparseVarsLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR**            vars,               /**< array to store the parsed variables */
   SCIP_Real*            vals,               /**< array to store the parsed coefficients */
   int*                  nvars,              /**< pointer to store number of parsed variables */
   int                   varssize,           /**< size of the variable array */
   int*                  requiredsize,       /**< pointer to store the required array size for the active variables */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successfully or not */
   )
{
   SCIP_VAR*** monomialvars;
   SCIP_Real** monomialexps;
   SCIP_Real*  monomialcoefs;
   int*        monomialnvars;
   int         nmonomials;

   SCIP_CALL( checkStage(scip, "SCIPparseVarsLinearsum", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(scip != NULL);
   assert(str != NULL);
   assert(vars != NULL || varssize == 0);
   assert(vals != NULL || varssize == 0);
   assert(nvars != NULL);
   assert(requiredsize != NULL);
   assert(endptr != NULL);
   assert(success != NULL);

   *requiredsize = 0;

   SCIP_CALL( SCIPparseVarsPolynomial(scip, str, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, &nmonomials, endptr, success) );

   if( !*success )
   {
      assert(nmonomials == 0); /* SCIPparseVarsPolynomial should have freed all buffers, so no need to call free here */
      return SCIP_OKAY;
   }

   /* check if linear sum is just "0" */
   if( nmonomials == 1 && monomialnvars[0] == 0 && monomialcoefs[0] == 0.0 )
   {
      *nvars = 0;
      *requiredsize = 0;

      SCIPfreeParseVarsPolynomialData(scip, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, nmonomials);

      return SCIP_OKAY;
   }

   *nvars = nmonomials;
   *requiredsize = nmonomials;

   /* if we have enough slots in the variables array, copy variables over */
   if( varssize >= nmonomials )
   {
      int v;

      for( v = 0; v < nmonomials; ++v )
      {
         if( monomialnvars[v] == 0 )
         {
            SCIPerrorMessage("constant in linear sum\n");
            *success = FALSE;
            break;
         }
         if( monomialnvars[v] > 1 || monomialexps[v][0] != 1.0 )
         {
            SCIPerrorMessage("nonlinear monomial in linear sum\n");
            *success = FALSE;
            break;
         }
         assert(monomialnvars[v]   == 1);
         assert(monomialvars[v][0] != NULL);
         assert(monomialexps[v][0] == 1.0);

         vars[v] = monomialvars[v][0];
         vals[v] = monomialcoefs[v];
      }
   }

   SCIPfreeParseVarsPolynomialData(scip, &monomialvars, &monomialexps, &monomialcoefs, &monomialnvars, nmonomials);

   return SCIP_OKAY;
}

/** parse the given string as polynomial of variables and coefficients
 *  (c1 \<x11\>^e11 \<x12\>^e12 ... \<x1n\>^e1n + c2 \<x21\>^e21 \<x22\>^e22 ... + ... + cn \<xn1\>^en1 ...)
 *  (see SCIPwriteVarsPolynomial()); if it was successful, the pointer success is set to TRUE
 *
 *  The user has to call SCIPfreeParseVarsPolynomialData(scip, monomialvars, monomialexps,
 *  monomialcoefs, monomialnvars, *nmonomials) short after SCIPparseVarsPolynomial to free all the
 *  allocated memory again.  Do not keep the arrays created by SCIPparseVarsPolynomial around, since
 *  they use buffer memory that is intended for short term use only.
 *
 *  Parsing is stopped at the end of string (indicated by the \\0-character) or when no more monomials
 *  are recognized.
 */
SCIP_RETCODE SCIPparseVarsPolynomial(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           str,                /**< string to parse */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int*                  nmonomials,         /**< pointer to store number of parsed monomials */
   char**                endptr,             /**< pointer to store the final string position if successfully */
   SCIP_Bool*            success             /**< pointer to store the whether the parsing was successfully or not */
   )
{
   typedef enum
   {
      SCIPPARSEPOLYNOMIAL_STATE_BEGIN,       /* we are at the beginning of a monomial */
      SCIPPARSEPOLYNOMIAL_STATE_INTERMED,    /* we are in between the factors of a monomial */
      SCIPPARSEPOLYNOMIAL_STATE_COEF,        /* we parse the coefficient of a monomial */
      SCIPPARSEPOLYNOMIAL_STATE_VARS,        /* we parse monomial variables */
      SCIPPARSEPOLYNOMIAL_STATE_EXPONENT,    /* we parse the exponent of a variable */
      SCIPPARSEPOLYNOMIAL_STATE_END,         /* we are at the end the polynomial */
      SCIPPARSEPOLYNOMIAL_STATE_ERROR        /* a parsing error occured */
   } SCIPPARSEPOLYNOMIAL_STATES;

   SCIPPARSEPOLYNOMIAL_STATES state;
   int monomialssize;

   /* data of currently parsed monomial */
   int varssize;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Real* exponents;
   SCIP_Real coef;

   assert(scip != NULL);
   assert(str != NULL);
   assert(monomialvars != NULL);
   assert(monomialexps != NULL);
   assert(monomialnvars != NULL);
   assert(monomialcoefs != NULL);
   assert(nmonomials != NULL);
   assert(endptr != NULL);
   assert(success != NULL);

   SCIP_CALL( checkStage(scip, "SCIPparseVarsPolynomial", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *success = FALSE;
   *nmonomials = 0;
   monomialssize = 0;
   *monomialvars = NULL;
   *monomialexps = NULL;
   *monomialcoefs = NULL;
   *monomialnvars = NULL;

   /* initialize state machine */
   state = SCIPPARSEPOLYNOMIAL_STATE_BEGIN;
   varssize = 0;
   nvars = 0;
   vars = NULL;
   exponents = NULL;
   coef = SCIP_INVALID;

   SCIPdebugMessage("parsing polynomial from '%s'\n", str);

   while( *str && state != SCIPPARSEPOLYNOMIAL_STATE_END && state != SCIPPARSEPOLYNOMIAL_STATE_ERROR )
   {
      /* skip white space */
      while( isspace((unsigned char)*str) )
         str++;

      assert(state != SCIPPARSEPOLYNOMIAL_STATE_END);

      switch( state )
      {
      case SCIPPARSEPOLYNOMIAL_STATE_BEGIN:
      {
         if( coef != SCIP_INVALID  )
         {
            SCIPdebugMessage("push monomial with coefficient <%g> and <%d> vars\n", coef, nvars);
            /* push previous monomial */
            if( monomialssize <= *nmonomials )
            {
               monomialssize = SCIPcalcMemGrowSize(scip, *nmonomials+1);

               SCIP_CALL( SCIPreallocBufferArray(scip, monomialvars,  monomialssize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, monomialexps,  monomialssize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, monomialnvars, monomialssize) );
               SCIP_CALL( SCIPreallocBufferArray(scip, monomialcoefs, monomialssize) );
            }

            if( nvars > 0 )
            {
               SCIP_CALL( SCIPduplicateBufferArray(scip, &(*monomialvars)[*nmonomials], vars, nvars) );
               SCIP_CALL( SCIPduplicateBufferArray(scip, &(*monomialexps)[*nmonomials], exponents, nvars) );
            }
            else
            {
               (*monomialvars)[*nmonomials] = NULL;
               (*monomialexps)[*nmonomials] = NULL;
            }
            (*monomialcoefs)[*nmonomials] = coef;
            (*monomialnvars)[*nmonomials] = nvars;
            ++*nmonomials;

            nvars = 0;
            coef = SCIP_INVALID;
         }

         if( *str == '<' )
         {
            /* there seem to come a variable at the beginning of a monomial
             * so assume the coefficient is 1.0
             */
            state = SCIPPARSEPOLYNOMIAL_STATE_VARS;
            coef = 1.0;
            break;
         }
         if( *str == '-' || *str == '+' || isdigit(*str) )
         {
            state = SCIPPARSEPOLYNOMIAL_STATE_COEF;
            break;
         }

         state = SCIPPARSEPOLYNOMIAL_STATE_END;

         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_INTERMED:
      {
         if( *str == '<' )
         {
            /* there seem to come another variable */
            state = SCIPPARSEPOLYNOMIAL_STATE_VARS;
            break;
         }

         if( *str == '-' || *str == '+' || isdigit(*str) )
         {
            /* there seem to come a coefficient, which means the next monomial */
            state = SCIPPARSEPOLYNOMIAL_STATE_BEGIN;
            break;
         }

         /* since we cannot detect the symbols we stop parsing the polynomial */
         state = SCIPPARSEPOLYNOMIAL_STATE_END;
         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_COEF:
      {
         if( *str == '+' && !isdigit(str[1]) )
         {
            /* only a plus sign, without number */
            coef =  1.0;
            ++str;
         }
         else if( *str == '-' && !isdigit(str[1]) )
         {
            /* only a minus sign, without number */
            coef = -1.0;
            ++str;
         }
         else if( SCIPstrToRealValue(str, &coef, endptr) )
         {
            str = *endptr;
         }
         else
         {
            SCIPerrorMessage("could not parse number in the beginning of '%s'\n", str);
            state = SCIPPARSEPOLYNOMIAL_STATE_ERROR;
            break;
         }

         /* after the coefficient we go into the intermediate state, i.e., expecting next variables */
         state = SCIPPARSEPOLYNOMIAL_STATE_INTERMED;

         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_VARS:
      {
         SCIP_VAR* var;

         assert(*str == '<');

         /* parse variable name */
         SCIP_CALL( SCIPparseVarName(scip, str, &var, endptr) );

         /* check if variable name was parsed */
         if( *endptr == str )
         {
            state = SCIPPARSEPOLYNOMIAL_STATE_END;
            break;
         }

         if( var == NULL )
         {
            SCIPerrorMessage("did not find variable in the beginning of %s\n", str);
            state = SCIPPARSEPOLYNOMIAL_STATE_ERROR;
            break;
         }

         str = *endptr;

         /* add variable to vars array */
         if( nvars + 1 > varssize )
         {
            varssize = SCIPcalcMemGrowSize(scip, nvars+1);
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars,      varssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &exponents, varssize) );
         }
         vars[nvars] = var;
         exponents[nvars] = 1.0;
         ++nvars;

         str = *endptr;

         if( *str == '^' )
            state = SCIPPARSEPOLYNOMIAL_STATE_EXPONENT;
         else
            state = SCIPPARSEPOLYNOMIAL_STATE_INTERMED;

         break;
      }

      case SCIPPARSEPOLYNOMIAL_STATE_EXPONENT:
      {
         assert(*str == '^');
         assert(nvars > 0); /* we should be in a monomial that has already a variable */
         ++str;

         if( !SCIPstrToRealValue(str, &exponents[nvars-1], endptr) )
         {
            SCIPerrorMessage("could not parse number in the beginning of '%s'\n", str);
            state = SCIPPARSEPOLYNOMIAL_STATE_ERROR;
            break;
         }
         str = *endptr;

         /* after the exponent we go into the intermediate state, i.e., expecting next variables */
         state = SCIPPARSEPOLYNOMIAL_STATE_INTERMED;
         break;
      }

      default:
         SCIPerrorMessage("unexpected state\n");
         return SCIP_ERROR;
      }
   }

   /* set end pointer */
   *endptr = (char*)str;

   /* check state at end of string */
   switch( state )
   {
   case SCIPPARSEPOLYNOMIAL_STATE_BEGIN:
   case SCIPPARSEPOLYNOMIAL_STATE_END:
   case SCIPPARSEPOLYNOMIAL_STATE_INTERMED:
   {
      if( coef != SCIP_INVALID )
      {
         /* push last monomial */
         SCIPdebugMessage("push monomial with coefficient <%g> and <%d> vars\n", coef, nvars);
         if( monomialssize <= *nmonomials )
         {
            monomialssize = *nmonomials+1;
            SCIP_CALL( SCIPreallocBufferArray(scip, monomialvars,  monomialssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, monomialexps,  monomialssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, monomialnvars, monomialssize) );
            SCIP_CALL( SCIPreallocBufferArray(scip, monomialcoefs, monomialssize) );
         }

         if( nvars > 0 )
         {
            /* shrink vars and exponents array to needed size and take over ownership */
            SCIP_CALL( SCIPreallocBufferArray(scip, &vars,      nvars) );
            SCIP_CALL( SCIPreallocBufferArray(scip, &exponents, nvars) );
            (*monomialvars)[*nmonomials] = vars;
            (*monomialexps)[*nmonomials] = exponents;
            vars = NULL;
            exponents = NULL;
            varssize = 0;
         }
         else
         {
            (*monomialvars)[*nmonomials] = NULL;
            (*monomialexps)[*nmonomials] = NULL;
         }
         (*monomialcoefs)[*nmonomials] = coef;
         (*monomialnvars)[*nmonomials] = nvars;
         ++*nmonomials;
      }

      *success = TRUE;
      break;
   }

   case SCIPPARSEPOLYNOMIAL_STATE_COEF:
   case SCIPPARSEPOLYNOMIAL_STATE_VARS:
   case SCIPPARSEPOLYNOMIAL_STATE_EXPONENT:
   {
      SCIPerrorMessage("unexpected parsing state at end of polynomial string\n");
   }

   case SCIPPARSEPOLYNOMIAL_STATE_ERROR:
      break;
   }

   /* free memory to store current monomial, if still existing */
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &exponents);
   varssize = 0;

   if( *success && *nmonomials > 0 )
   {
      /* shrink arrays to required size, so we do not need to keep monomialssize around */
      assert(*nmonomials <= monomialssize);
      SCIP_CALL( SCIPreallocBufferArray(scip, monomialvars,  *nmonomials) );
      SCIP_CALL( SCIPreallocBufferArray(scip, monomialexps,  *nmonomials) );
      SCIP_CALL( SCIPreallocBufferArray(scip, monomialnvars, *nmonomials) );
      SCIP_CALL( SCIPreallocBufferArray(scip, monomialcoefs, *nmonomials) );

      /* SCIPwriteVarsPolynomial(scip, NULL, *monomialvars, *monomialexps, *monomialcoefs, *monomialnvars, *nmonomials, FALSE); */
   }
   else
   {
      /* in case of error, cleanup all data here */
      SCIPfreeParseVarsPolynomialData(scip, monomialvars, monomialexps, monomialcoefs, monomialnvars, *nmonomials);
      *nmonomials = 0;
   }

   return SCIP_OKAY;
}

/** frees memory allocated when parsing a polynomial from a string */
void SCIPfreeParseVarsPolynomialData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR****          monomialvars,       /**< pointer to store arrays with variables for each monomial */
   SCIP_Real***          monomialexps,       /**< pointer to store arrays with variable exponents */
   SCIP_Real**           monomialcoefs,      /**< pointer to store array with monomial coefficients */
   int**                 monomialnvars,      /**< pointer to store array with number of variables for each monomial */
   int                   nmonomials          /**< pointer to store number of parsed monomials */
   )
{
   int i;

   assert(scip != NULL);
   assert(monomialvars  != NULL);
   assert(monomialexps  != NULL);
   assert(monomialcoefs != NULL);
   assert(monomialnvars != NULL);
   assert((*monomialvars  != NULL) == (nmonomials > 0));
   assert((*monomialexps  != NULL) == (nmonomials > 0));
   assert((*monomialcoefs != NULL) == (nmonomials > 0));
   assert((*monomialnvars != NULL) == (nmonomials > 0));

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfreeParseVarsPolynomialData", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( nmonomials == 0 )
      return;

   for( i = 0; i < nmonomials; ++i )
   {
      SCIPfreeBufferArrayNull(scip, &(*monomialvars)[i]);
      SCIPfreeBufferArrayNull(scip, &(*monomialexps)[i]);
   }

   SCIPfreeBufferArray(scip, monomialvars);
   SCIPfreeBufferArray(scip, monomialexps);
   SCIPfreeBufferArray(scip, monomialcoefs);
   SCIPfreeBufferArray(scip, monomialnvars);
}

/** increases usage counter of variable */
SCIP_RETCODE SCIPcaptureVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to capture */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcaptureVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   SCIPvarCapture(var);

   return SCIP_OKAY;
}

/** decreases usage counter of variable, if the usage pointer reaches zero the variable gets freed
 *
 *  @note the pointer of the variable will be NULLed
 */
SCIP_RETCODE SCIPreleaseVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var                 /**< pointer to variable */
   )
{
   assert(var != NULL);
   assert(*var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPreleaseVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      if( !SCIPvarIsTransformed(*var) && (*var)->nuses == 1 )
      {
         SCIPerrorMessage("cannot release last use of original variable while the transformed problem exists\n");
         return SCIP_INVALIDCALL;
      }
      SCIP_CALL( SCIPvarRelease(var, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets and captures transformed variable of a given variable; if the variable is not yet transformed,
 *  a new transformed variable for this variable is created
 */
SCIP_RETCODE SCIPtransformVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get/create transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   assert(transvar != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtransformVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPvarIsTransformed(var) )
   {
      *transvar = var;
      SCIPvarCapture(*transvar);
   }
   else
   {
      SCIP_CALL( SCIPvarTransform(var, scip->mem->probmem, scip->set, scip->stat, scip->origprob->objsense, transvar) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed variables for an array of variables;
 *  if a variable of the array is not yet transformed, a new transformed variable for this variable is created;
 *  it is possible to call this method with vars == transvars
 */
SCIP_RETCODE SCIPtransformVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get/create transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get/create transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || transvars != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtransformVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarIsTransformed(vars[v]) )
      {
         transvars[v] = vars[v];
         SCIPvarCapture(transvars[v]);
      }
      else
      {
         SCIP_CALL( SCIPvarTransform(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->origprob->objsense,
               &transvars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed variable of a given variable;
 *  returns NULL as transvar, if transformed variable is not yet existing
 */
SCIP_RETCODE SCIPgetTransformedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get transformed variable for */
   SCIP_VAR**            transvar            /**< pointer to store the transformed variable */
   )
{
   assert(transvar != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetTransformedVar", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPvarIsTransformed(var) )
      *transvar = var;
   else
   {
      SCIP_CALL( SCIPvarGetTransformed(var, scip->mem->probmem, scip->set, scip->stat, transvar) );
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed variables for an array of variables;
 *  stores NULL in a transvars slot, if the transformed variable is not yet existing;
 *  it is possible to call this method with vars == transvars, but remember that variables that are not
 *  yet transformed will be replaced with NULL
 */
SCIP_RETCODE SCIPgetTransformedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get transformed variables for */
   SCIP_VAR**            vars,               /**< array with variables to get transformed variables for */
   SCIP_VAR**            transvars           /**< array to store the transformed variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || transvars != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetTransformedVars", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarIsTransformed(vars[v]) )
         transvars[v] = vars[v];
      else
      {
         SCIP_CALL( SCIPvarGetTransformed(vars[v], scip->mem->probmem, scip->set, scip->stat, &transvars[v]) );
      }
   }

   return SCIP_OKAY;
}

/** gets negated variable x' = lb + ub - x of variable x; negated variable is created, if not yet existing */
SCIP_RETCODE SCIPgetNegatedVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get negated variable for */
   SCIP_VAR**            negvar              /**< pointer to store the negated variable */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNegatedVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPvarNegate(var, scip->mem->probmem, scip->set, scip->stat, negvar) );

   return SCIP_OKAY;
}

/** gets negated variables x' = lb + ub - x of variables x; negated variables are created, if not yet existing */
SCIP_RETCODE SCIPgetNegatedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get negated variables for */
   SCIP_VAR**            vars,               /**< array of variables to get negated variables for */
   SCIP_VAR**            negvars             /**< array to store the negated variables */
   )
{
   int v;

   SCIP_CALL( checkStage(scip, "SCIPgetNegatedVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarNegate(vars[v], scip->mem->probmem, scip->set, scip->stat, &(negvars[v])) );
   }

   return SCIP_OKAY;
}

/** gets a binary variable that is equal to the given binary variable, and that is either active, fixed, or
 *  multi-aggregated, or the negated variable of an active, fixed, or multi-aggregated variable
 */
SCIP_RETCODE SCIPgetBinvarRepresentative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to get binary representative for */
   SCIP_VAR**            repvar,             /**< pointer to store the binary representative */
   SCIP_Bool*            negated             /**< pointer to store whether the negation of an active variable was returned */
   )
{
   assert(scip != NULL);
   assert(var != NULL);
   assert(repvar != NULL);
   assert(negated != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetBinvarRepresentative", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* get the active representative of the given variable */
   *repvar = var;
   *negated = FALSE;
   SCIP_CALL( SCIPvarGetProbvarBinary(repvar, negated) );

   /* negate the representative, if it corresponds to the negation of the given variable */
   if( *negated )
   {
      SCIP_CALL( SCIPgetNegatedVar(scip, *repvar, repvar) );
   }

   return SCIP_OKAY;
}

/** gets binary variables that are equal to the given binary variables, and which are either active, fixed, or 
 *  multi-aggregated, or the negated variables of active, fixed, or multi-aggregated variables
 */
extern
SCIP_RETCODE SCIPgetBinvarRepresentatives(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of binary variables to get representatives for */
   SCIP_VAR**            vars,               /**< binary variables to get binary representatives for */
   SCIP_VAR**            repvars,            /**< array to store the binary representatives */
   SCIP_Bool*            negated             /**< array to store whether the negation of an active variable was returned */
   )
{
   int v;

   assert(scip != NULL);
   assert(vars != NULL || nvars == 0);
   assert(repvars != NULL || nvars == 0);
   assert(negated != NULL || nvars == 0);

   SCIP_CALL( checkStage(scip, "SCIPgetBinvarRepresentatives", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( nvars == 0 )
      return SCIP_OKAY;

   /* get the active representative of the given variable */
   BMScopyMemoryArray(repvars, vars, nvars);
   BMSclearMemoryArray(negated, nvars);
   SCIP_CALL( SCIPvarsGetProbvarBinary(&repvars, &negated, nvars) );

   /* negate the representatives, if they correspond to the negation of the given variables */
   for( v = nvars - 1; v >= 0; --v )
      if( negated[v] )
      {
         SCIP_CALL( SCIPgetNegatedVar(scip, repvars[v], &(repvars[v])) );
      }
   
   return SCIP_OKAY;
}

/** flattens aggregation graph of multi-aggregated variable in order to avoid exponential recursion later on */
SCIP_RETCODE SCIPflattenVarAggregationGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   assert( scip != NULL );
   assert( var != NULL );
   SCIP_CALL( checkStage(scip, "SCIPflattenVarAggregationGraph", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarFlattenAggregationGraph(var, scip->mem->probmem, scip->set) );
   return SCIP_OKAY;
}

/** Transforms a given linear sum of variables, that is a_1*x_1 + ... + a_n*x_n + c into a corresponding linear sum of
 *  active variables, that is b_1*y_1 + ... + b_m*y_m + d.
 *
 *  If the number of needed active variables is greater than the available slots in the variable array, nothing happens
 *  except that the required size is stored in the corresponding variable (requiredsize). Otherwise, the active variable
 *  representation is stored in the variable array, scalar array and constant.
 *
 *  The reason for this approach is that we cannot reallocate memory, since we do not know how the memory has been
 *  allocated (e.g., by a C++ 'new' or SCIP functions).
 *
 *  @note The resulting linear sum is stored into the given variable array, scalar array, and constant. That means the
 *        given entries are overwritten.
 *
 *  @note That method can be used to convert a single variables into variable space of active variables. Therefore call
 *        the method with the linear sum 1.0*x + 0.0.
 */
SCIP_RETCODE SCIPgetProbvarLinearSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array x_1, ..., x_n in the linear sum which will be
                                              *   overwritten by the variable array y_1, ..., y_m in the linear sum
                                              *   w.r.t. active variables */
   SCIP_Real*            scalars,            /**< scalars a_1, ..., a_n in linear sum which will be overwritten to the
                                              *   scalars b_1, ..., b_m in the linear sum of the active variables  */
   int*                  nvars,              /**< pointer to number of variables in the linear sum which will be
                                              *   overwritten by the number of variables in the linear sum corresponding
                                              *   to the active variables */
   int                   varssize,           /**< available slots in vars and scalars array which is needed to check if
                                              *   the array are large enough for the linear sum w.r.t. active
                                              *   variables */
   SCIP_Real*            constant,           /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c which
                                              *   will chnage to constant d in the linear sum b_1*y_1 + ... + b_m*y_m +
                                              *   d w.r.t. the active variables */
   int*                  requiredsize,       /**< pointer to store the required array size for the linear sum w.r.t. the
                                              *   active variables */
   SCIP_Bool             mergemultiples      /**< should multiple occurrences of a var be replaced by a single coeff? */
   )
{
   assert( scip != NULL );
   assert( nvars != NULL );
   assert( vars != NULL || *nvars == 0 );
   assert( scalars != NULL || *nvars == 0 );
   assert( constant != NULL );
   assert( requiredsize != NULL );
   assert( *nvars <= varssize );

   SCIP_CALL( checkStage(scip, "SCIPgetProbvarLinearSum", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );
   SCIP_CALL( SCIPvarGetActiveRepresentatives(scip->set, vars, scalars, nvars, varssize, constant, requiredsize, mergemultiples) );

   return SCIP_OKAY;
}


/** returns the reduced costs of the variable in the current node's LP relaxation;
 *  the current node has to have a feasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 */
SCIP_Real SCIPgetVarRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPgetVarRedcost(scip,var->data.original.transvar);

   case SCIP_VARSTATUS_COLUMN:
      return SCIPgetColRedcost(scip,SCIPvarGetCol(var));

   case SCIP_VARSTATUS_LOOSE:
      return SCIP_INVALID;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      return 0;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** returns the Farkas coefficient of the variable in the current node's LP relaxation;
 *  the current node has to have an infeasible LP.
 *
 *  returns SCIP_INVALID if the variable is active but not in the current LP;
 *  returns 0 if the variable has been aggregated out or fixed in presolving.
 */
SCIP_Real SCIPgetVarFarkasCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get reduced costs, should be a column in current node LP */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
      if( var->data.original.transvar == NULL )
         return SCIP_INVALID;
      return SCIPgetVarFarkasCoef(scip,var->data.original.transvar);

   case SCIP_VARSTATUS_COLUMN:
      return SCIPgetColFarkasCoef(scip,SCIPvarGetCol(var));

   case SCIP_VARSTATUS_LOOSE:
      return SCIP_INVALID;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      return 0;

   default:
      SCIPerrorMessage("unknown variable status\n");
      SCIPABORT();
      return SCIP_INVALID; /*lint !e527*/
   }
}

/** gets solution value for variable in current node */
SCIP_Real SCIPgetVarSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get solution value for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPvarGetSol(var, SCIPtreeHasCurrentNodeLP(scip->tree));
}

/** gets solution values of multiple variables in current node */
SCIP_RETCODE SCIPgetVarSols(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetVarSols", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPvarGetLPSol(vars[v]);
   }
   else
   {
      for( v = 0; v < nvars; ++v )
         vals[v] = SCIPvarGetPseudoSol(vars[v]);
   }

   return SCIP_OKAY;
}

/** sets the solution value of all variables in the global relaxation solution to zero */
SCIP_RETCODE SCIPclearRelaxSolVals(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(scip != NULL);

   SCIP_CALL( checkStage(scip, "SCIPclearRelaxSolVals", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* the relaxation solution is already cleared */
   if( SCIPrelaxationIsSolZero(scip->relaxation) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_CALL( SCIPvarSetRelaxSol(vars[v], scip->set, scip->relaxation, 0.0, FALSE) );
   }
   
   SCIPrelaxationSetSolObj(scip->relaxation, 0.0);
   SCIPrelaxationSetSolZero(scip->relaxation, TRUE);

   return SCIP_OKAY;
}

/** sets the value of the given variable in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  You can use SCIPclearRelaxSolVals() to set all values to zero, initially;
 *  after setting all solution values, you have to call SCIPmarkRelaxSolValid() 
 *  to inform SCIP that the stored solution is valid
 */
SCIP_RETCODE SCIPsetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to set value for */
   SCIP_Real             val                 /**< solution value of variable */
   )
{
   assert(scip != NULL);

   SCIP_CALL( checkStage(scip, "SCIPsetRelaxSolVal", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarSetRelaxSol(var, scip->set, scip->relaxation, val, TRUE) );

   if( val != 0.0 )
      SCIPrelaxationSetSolZero(scip->relaxation, FALSE);
   SCIPrelaxationSetSolValid(scip->relaxation, FALSE);

   return SCIP_OKAY;
}

/** sets the values of the given variables in the global relaxation solution;
 *  this solution can be filled by the relaxation handlers  and can be used by heuristics and for separation;
 *  the solution is automatically cleared, s.t. all other variables get value 0.0
 */
SCIP_RETCODE SCIPsetRelaxSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables to set relaxation solution value for */
   SCIP_VAR**            vars,               /**< array with variables to set value for */
   SCIP_Real*            vals                /**< array with solution values of variables */
   )
{
   int v;

   assert(scip != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( checkStage(scip, "SCIPsetRelaxSolVals", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPclearRelaxSolVals(scip) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_CALL( SCIPvarSetRelaxSol(vars[v], scip->set, scip->relaxation, vals[v], TRUE) );
   }

   SCIPrelaxationSetSolZero(scip->relaxation, FALSE);
   SCIPrelaxationSetSolValid(scip->relaxation, TRUE); 

   return SCIP_OKAY;
}

/** sets the values of the variables in the global relaxation solution to the values 
 *  in the given primal solution; the relaxation solution can be filled by the relaxation hanlders
 *  and might be used by heuristics and for separation
 */
SCIP_RETCODE SCIPsetRelaxSolValsSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal relaxation solution */ 
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int nvars;
   int v;

   assert(scip != NULL);

   SCIP_CALL( checkStage(scip, "SCIPsetRelaxSolValsSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* alloc buffer array for solution values of the variables and get the values */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, vars, vals) );

   SCIP_CALL( SCIPclearRelaxSolVals(scip) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_CALL( SCIPvarSetRelaxSol(vars[v], scip->set, scip->relaxation, vals[v], FALSE) );
   }       

   SCIPrelaxationSetSolObj(scip->relaxation, SCIPsolGetObj(sol, scip->set, scip->transprob));

   SCIPrelaxationSetSolZero(scip->relaxation, FALSE);
   SCIPrelaxationSetSolValid(scip->relaxation, TRUE);

   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

/** returns whether the relaxation solution is valid */
SCIP_Bool SCIPisRelaxSolValid(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPisRelaxSolValid", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPrelaxationIsSolValid(scip->relaxation);
}

/** informs SCIP, that the relaxation solution is valid */
SCIP_RETCODE SCIPmarkRelaxSolValid(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( checkStage(scip, "SCIPmarkRelaxSolValid", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPrelaxationSetSolValid(scip->relaxation, TRUE);

   return SCIP_OKAY;
}

/** informs SCIP, that the relaxation solution is invalid */
SCIP_RETCODE SCIPmarkRelaxSolInvalid(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL( checkStage(scip, "SCIPmarkRelaxSolInvalid", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPrelaxationSetSolValid(scip->relaxation, FALSE);

   return SCIP_OKAY;
}

/** gets the relaxation solution value of the given variable */
SCIP_Real SCIPgetRelaxSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRelaxSolVal", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("Relaxation Solution is not valid!\n");
      SCIPABORT();
   }

   return SCIPvarGetRelaxSol(var, scip->set);
}

/** gets the relaxation solution objective value */
SCIP_Real SCIPgetRelaxSolObj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRelaxSolObj", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("Relaxation Solution is not valid!\n");
      SCIPABORT();
   }

   return SCIPrelaxationGetSolObj(scip->relaxation);
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPstartStrongbranch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );
   SCIP_CALL( checkStage(scip, "SCIPstartStrongbranch", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpStartStrongbranch(scip->lp) );

   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPendStrongbranch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   SCIP_CALL( checkStage(scip, "SCIPendStrongbranch", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpEndStrongbranch(scip->lp) );

   return SCIP_OKAY;
}

/** gets strong branching information on column variable with fractional value */
SCIP_RETCODE SCIPgetVarStrongbranchFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL* col;

   assert(lperror != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetVarStrongbranchFrac", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( downvalid != NULL )
      *downvalid = FALSE;
   if( upvalid != NULL )
      *upvalid = FALSE;
   if( downinf != NULL )
      *downinf = FALSE;
   if( upinf != NULL )
      *upinf = FALSE;
   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
   {
      SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
      return SCIP_OKAY;
   }

   /* call strong branching for column */
   SCIP_CALL( SCIPcolGetStrongbranchFrac(col, scip->set, scip->stat, scip->lp, itlim,
         down, up, downvalid, upvalid, lperror) );

   /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
    * declare the sub nodes infeasible
    */
   if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
   {
      SCIP_Bool downcutoff;
      SCIP_Bool upcutoff;

      downcutoff = col->sbdownvalid && SCIPsetIsGE(scip->set, col->sbdown, scip->lp->cutoffbound);
      upcutoff = col->sbupvalid && SCIPsetIsGE(scip->set, col->sbup, scip->lp->cutoffbound);
      if( downinf != NULL )
         *downinf = downcutoff;
      if( upinf != NULL )
         *upinf = upcutoff;

      /* analyze infeasible strong branching sub problems:
       * because the strong branching's bound change is necessary for infeasibility, it cannot be undone;
       * therefore, infeasible strong branchings on non-binary variables will not produce a valid conflict constraint
       */
      if( scip->set->conf_enable && scip->set->conf_usesb && scip->set->nconflicthdlrs > 0
         && SCIPvarIsBinary(var) && SCIPtreeGetCurrentDepth(scip->tree) > 0 )
      {
         if( (downcutoff && SCIPsetFeasCeil(scip->set, col->primsol-1.0) >= col->lb - 0.5)
            || (upcutoff && SCIPsetFeasFloor(scip->set, col->primsol+1.0) <= col->ub + 0.5) )
         {
            SCIP_CALL( SCIPconflictAnalyzeStrongbranch(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
                  scip->transprob, scip->tree, scip->lp, col, downconflict, upconflict) );
         }
      }
   }

   return SCIP_OKAY;
}

/** gets strong branching information on column variable with integral value */
SCIP_RETCODE SCIPgetVarStrongbranchInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get strong branching values for */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< pointer to store whether the downwards branch is infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< pointer to store whether the upwards branch is infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict,         /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL* col;

   assert(lperror != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetVarStrongbranchInt", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( downvalid != NULL )
      *downvalid = FALSE;
   if( upvalid != NULL )
      *upvalid = FALSE;
   if( downinf != NULL )
      *downinf = FALSE;
   if( upinf != NULL )
      *upinf = FALSE;
   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   col = SCIPvarGetCol(var);
   assert(col != NULL);

   if( !SCIPcolIsInLP(col) )
   {
      SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
      return SCIP_OKAY;
   }

   /* call strong branching for column */
   SCIP_CALL( SCIPcolGetStrongbranchInt(col, scip->set, scip->stat, scip->lp, itlim,
         down, up, downvalid, upvalid, lperror) );

   /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
    * declare the sub nodes infeasible
    */
   if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
   {
      SCIP_Bool downcutoff;
      SCIP_Bool upcutoff;

      downcutoff = col->sbdownvalid && SCIPsetIsGE(scip->set, col->sbdown, scip->lp->cutoffbound);
      upcutoff = col->sbupvalid && SCIPsetIsGE(scip->set, col->sbup, scip->lp->cutoffbound);
      if( downinf != NULL )
         *downinf = downcutoff;
      if( upinf != NULL )
         *upinf = upcutoff;

      /* analyze infeasible strong branching sub problems:
       * because the strong branching's bound change is necessary for infeasibility, it cannot be undone;
       * therefore, infeasible strong branchings on non-binary variables will not produce a valid conflict constraint
       */
      if( scip->set->conf_enable && scip->set->conf_usesb && scip->set->nconflicthdlrs > 0
         && SCIPvarIsBinary(var) && SCIPtreeGetCurrentDepth(scip->tree) > 0 )
      {
         if( (downcutoff && SCIPsetFeasCeil(scip->set, col->primsol-1.0) >= col->lb - 0.5)
            || (upcutoff && SCIPsetFeasFloor(scip->set, col->primsol+1.0) <= col->ub + 0.5) )
         {
            SCIP_CALL( SCIPconflictAnalyzeStrongbranch(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
                  scip->transprob, scip->tree, scip->lp, col, downconflict, upconflict) );
         }
      }
   }

   return SCIP_OKAY;
}

/** gets strong branching information on column variables with fractional values */
SCIP_RETCODE SCIPgetVarsStrongbranchesFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables to get strong branching values for */
   int                   nvars,              /**< number of variables */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching variables down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching variables up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< array to store whether the downward branches are infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< array to store whether the upward branches are infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< array to store whether conflict constraints were created for
                                              *   infeasible downward branches, or NULL */
   SCIP_Bool*            upconflict,         /**< array to store whether conflict constraints were created for
                                              *   infeasible upward branches, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL** cols;
   int j;

   assert(lperror != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetVarsStrongbranchesFrac", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert( vars != NULL );

   /* set up data */
   cols = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &cols, nvars) );
   assert(cols != NULL);
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      SCIP_COL* col;

      if( downvalid != NULL )
         downvalid[j] = FALSE;
      if( upvalid != NULL )
         upvalid[j] = FALSE;
      if( downinf != NULL )
         downinf[j] = FALSE;
      if( upinf != NULL )
         upinf[j] = FALSE;
      if( downconflict != NULL )
         downconflict[j] = FALSE;
      if( upconflict != NULL )
         upconflict[j] = FALSE;
      
      var = vars[j];
      assert( var != NULL );
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      cols[j] = col;

      if( !SCIPcolIsInLP(col) )
      {
         SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
   }
   else
   {
      /* call strong branching for columns */
      SCIP_CALL( SCIPcolGetStrongbranchesFrac(cols, nvars, scip->set, scip->stat, scip->lp, itlim,
            down, up, downvalid, upvalid, lperror) );
      
      /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
       * declare the sub nodes infeasible
       */
      if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
      {
         for( j = 0; j < nvars; ++j )
         {
            SCIP_Bool downcutoff;
            SCIP_Bool upcutoff;
            SCIP_COL* col;
            
            col = cols[j];
            downcutoff = col->sbdownvalid && SCIPsetIsGE(scip->set, col->sbdown, scip->lp->cutoffbound);
            upcutoff = col->sbupvalid && SCIPsetIsGE(scip->set, col->sbup, scip->lp->cutoffbound);
            if( downinf != NULL )
               downinf[j] = downcutoff;
            if( upinf != NULL )
               upinf[j] = upcutoff;
            
            /* analyze infeasible strong branching sub problems:
             * because the strong branching's bound change is necessary for infeasibility, it cannot be undone;
             * therefore, infeasible strong branchings on non-binary variables will not produce a valid conflict constraint
             */
            if( scip->set->conf_enable && scip->set->conf_usesb && scip->set->nconflicthdlrs > 0
               && SCIPvarIsBinary(vars[j]) && SCIPtreeGetCurrentDepth(scip->tree) > 0 )
            {
               if( (downcutoff && SCIPsetFeasCeil(scip->set, col->primsol-1.0) >= col->lb - 0.5)
                  || (upcutoff && SCIPsetFeasFloor(scip->set, col->primsol+1.0) <= col->ub + 0.5) )
               {
                  assert(downconflict != NULL);
                  assert(upconflict   != NULL);
                  SCIP_CALL( SCIPconflictAnalyzeStrongbranch(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
                        scip->transprob, scip->tree, scip->lp, col, &(downconflict[j]), &(upconflict[j])) );
               }
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &cols);

   return SCIP_OKAY;
}

/** gets strong branching information on column variables with integral values */
SCIP_RETCODE SCIPgetVarsStrongbranchesInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables to get strong branching values for */
   int                   nvars,              /**< number of variables */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching variables down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching variables up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are valid dual bounds, or NULL;
                                              *   otherwise, they can only be used as an estimate value */
   SCIP_Bool*            downinf,            /**< array to store whether the downward branches are infeasible, or NULL */
   SCIP_Bool*            upinf,              /**< array to store whether the upward branches are infeasible, or NULL */
   SCIP_Bool*            downconflict,       /**< array to store whether conflict constraints were created for
                                              *   infeasible downward branches, or NULL */
   SCIP_Bool*            upconflict,         /**< array to store whether conflict constraints were created for
                                              *   infeasible upward branches, or NULL */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred or the
                                              *   solving process should be stopped (e.g., due to a time limit) */
   )
{
   SCIP_COL** cols;
   int j;

   assert(lperror != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetVarsStrongbranchesInt", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert( vars != NULL );

   /* set up data */
   cols = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &cols, nvars) );
   assert(cols != NULL);
   for( j = 0; j < nvars; ++j )
   {
      SCIP_VAR* var;
      SCIP_COL* col;

      if( downvalid != NULL )
         downvalid[j] = FALSE;
      if( upvalid != NULL )
         upvalid[j] = FALSE;
      if( downinf != NULL )
         downinf[j] = FALSE;
      if( upinf != NULL )
         upinf[j] = FALSE;
      if( downconflict != NULL )
         downconflict[j] = FALSE;
      if( upconflict != NULL )
         upconflict[j] = FALSE;
      
      var = vars[j];
      assert( var != NULL );
      if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      {
         SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable <%s>\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }

      col = SCIPvarGetCol(var);
      assert(col != NULL);
      cols[j] = col;

      if( !SCIPcolIsInLP(col) )
      {
         SCIPerrorMessage("cannot get strong branching information on variable <%s> not in current LP\n", SCIPvarGetName(var));
         SCIPfreeBufferArray(scip, &cols);
         return SCIP_INVALIDDATA;
      }
   }

   /* check if the solving process should be aborted */
   if( SCIPsolveIsStopped(scip->set, scip->stat, FALSE) )
   {
      /* mark this as if the LP failed */
      *lperror = TRUE;
   }
   else
   {
      /* call strong branching for columns */
      SCIP_CALL( SCIPcolGetStrongbranchesInt(cols, nvars, scip->set, scip->stat, scip->lp, itlim,
            down, up, downvalid, upvalid, lperror) );
      
      /* check, if the branchings are infeasible; in exact solving mode, we cannot trust the strong branching enough to
       * declare the sub nodes infeasible
       */
      if( !(*lperror) && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) && !scip->set->misc_exactsolve )
      {
         for( j = 0; j < nvars; ++j )
         {
            SCIP_Bool downcutoff;
            SCIP_Bool upcutoff;
            SCIP_COL* col;
            
            col = cols[j];
            downcutoff = col->sbdownvalid && SCIPsetIsGE(scip->set, col->sbdown, scip->lp->cutoffbound);
            upcutoff = col->sbupvalid && SCIPsetIsGE(scip->set, col->sbup, scip->lp->cutoffbound);
            if( downinf != NULL )
               downinf[j] = downcutoff;
            if( upinf != NULL )
               upinf[j] = upcutoff;
            
            /* analyze infeasible strong branching sub problems:
             * because the strong branching's bound change is necessary for infeasibility, it cannot be undone;
             * therefore, infeasible strong branchings on non-binary variables will not produce a valid conflict constraint
             */
            if( scip->set->conf_enable && scip->set->conf_usesb && scip->set->nconflicthdlrs > 0
               && SCIPvarIsBinary(vars[j]) && SCIPtreeGetCurrentDepth(scip->tree) > 0 )
            {
               if( (downcutoff && SCIPsetFeasCeil(scip->set, col->primsol-1.0) >= col->lb - 0.5)
                  || (upcutoff && SCIPsetFeasFloor(scip->set, col->primsol+1.0) <= col->ub + 0.5) )
               {
                  assert(downconflict != NULL);
                  assert(upconflict   != NULL);
                  SCIP_CALL( SCIPconflictAnalyzeStrongbranch(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
                        scip->transprob, scip->tree, scip->lp, col, &(downconflict[j]), &(upconflict[j])) );
               }
            }
         }
      }
   }
   SCIPfreeBufferArray(scip, &cols);

   return SCIP_OKAY;
}

/** gets strong branching information on COLUMN variable of the last SCIPgetVarStrongbranch() call;
 *  returns values of SCIP_INVALID, if strong branching was not yet called on the given variable;
 *  keep in mind, that the returned old values may have nothing to do with the current LP solution
 */
SCIP_RETCODE SCIPgetVarStrongbranchLast(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to get last strong branching values for */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound, or NULL;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Real*            solval,             /**< stores LP solution value of variable at the last strong branching call */
   SCIP_Real*            lpobjval            /**< stores LP objective value at last strong branching call, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetVarsStrongbranchLast", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
   {
      SCIPerrorMessage("cannot get strong branching information on non-COLUMN variable\n");
      return SCIP_INVALIDDATA;
   }

   SCIPcolGetStrongbranchLast(SCIPvarGetCol(var), down, up, downvalid, upvalid, solval, lpobjval);

   return SCIP_OKAY;
}

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given variable, or -1 if strong branching was never applied to the variable in current run
 */
SCIP_Longint SCIPgetVarStrongbranchNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarStrongbranchNode", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return -1;

   return SCIPcolGetStrongbranchNode(SCIPvarGetCol(var));
}

/** if strong branching was already applied on the variable at the current node, returns the number of LPs solved after
 *  the LP where the strong branching on this variable was applied;
 *  if strong branching was not yet applied on the variable at the current node, returns INT_MAX
 */
int SCIPgetVarStrongbranchLPAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get strong branching LP age for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarStrongbranchLPAge", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return INT_MAX;

   return SCIPcolGetStrongbranchLPAge(SCIPvarGetCol(var), scip->stat);
}

/** gets number of times, strong branching was applied in current run on the given variable */
int SCIPgetVarNStrongbranchs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get last strong branching node for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarNStrongbranchs", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_COLUMN )
      return 0;

   return SCIPcolGetNStrongbranchs(SCIPvarGetCol(var));
}

/** adds given values to lock numbers of variable for rounding */
SCIP_RETCODE SCIPaddVarLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   nlocksdown,         /**< modification in number of rounding down locks */
   int                   nlocksup            /**< modification in number of rounding up locks */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarLocks", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      SCIP_CALL( SCIPvarAddLocks(var, scip->mem->probmem, scip->set, scip->eventqueue, nlocksdown, nlocksup) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** locks rounding of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 */
SCIP_RETCODE SCIPlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be locked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be locked in upwards direction? */
   )
{
   int nlocksdown;
   int nlocksup;

   SCIP_CALL( checkStage(scip, "SCIPlockVarCons", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE) );

   nlocksdown = 0;
   nlocksup = 0;
   if( SCIPconsIsLockedPos(cons) )
   {
      if( lockdown )
         nlocksdown++;
      if( lockup )
         nlocksup++;
   }
   if( SCIPconsIsLockedNeg(cons) )
   {
      if( lockdown )
         nlocksup++;
      if( lockup )
         nlocksdown++;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      SCIP_CALL( SCIPvarAddLocks(var, scip->mem->probmem, scip->set, scip->eventqueue, nlocksdown, nlocksup) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** unlocks rounding of variable with respect to the lock status of the constraint and its negation;
 *  this method should be called whenever the lock status of a variable in a constraint changes, for example if
 *  the coefficient of the variable changed its sign or if the left or right hand sides of the constraint were
 *  added or removed
 */
SCIP_RETCODE SCIPunlockVarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             lockdown,           /**< should the rounding be unlocked in downwards direction? */
   SCIP_Bool             lockup              /**< should the rounding be unlocked in upwards direction? */
   )
{
   int nlocksdown;
   int nlocksup;

   SCIP_CALL( checkStage(scip, "SCIPunlockVarCons", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE) );

   nlocksdown = 0;
   nlocksup = 0;
   if( SCIPconsIsLockedPos(cons) )
   {
      if( lockdown )
         nlocksdown++;
      if( lockup )
         nlocksup++;
   }
   if( SCIPconsIsLockedNeg(cons) )
   {
      if( lockdown )
         nlocksup++;
      if( lockup )
         nlocksdown++;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      /*lint -fallthrough*/
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      SCIP_CALL( SCIPvarAddLocks(var, scip->mem->probmem, scip->set, scip->eventqueue, -nlocksdown, -nlocksup) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** changes variable's objective value */
SCIP_RETCODE SCIPchgVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarObj", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgObj(var, scip->mem->probmem, scip->set, scip->primal, scip->lp, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
      SCIP_CALL( SCIPvarChgObj(var, scip->mem->probmem, scip->set, scip->primal, scip->lp, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** adds value to variable's objective value */
SCIP_RETCODE SCIPaddVarObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             addobj              /**< additional objective value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarObj", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarAddObj(var, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->primal,
            scip->tree, scip->lp, scip->eventqueue, addobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
      SCIP_CALL( SCIPvarAddObj(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
            scip->tree, scip->lp, scip->eventqueue, addobj) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) lower bound value;
 *  does not change the bounds of the variable
 */
SCIP_Real SCIPadjustedVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             lb                  /**< lower bound value to adjust */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPadjustedVarLb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPvarAdjustLb(var, scip->set, &lb);

   return lb;
}

/** returns the adjusted (i.e. rounded, if the given variable is of integral type) upper bound value;
 *  does not change the bounds of the variable
 */
SCIP_Real SCIPadjustedVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to adjust the bound for */
   SCIP_Real             ub                  /**< upper bound value to adjust */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPadjustedVarUb", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPvarAdjustUb(var, scip->set, &ub);

   return ub;
}

/** depending on SCIP's stage, changes lower bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarLb", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &newbound);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** depending on SCIP's stage, changes upper bound of variable in the problem, in preprocessing, or in current node;
 *  if possible, adjusts bound to integral value; doesn't store any inference information in the bound change, such
 *  that in conflict analysis, this change is treated like a branching decision
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarUb", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustUb(var, scip->set, &newbound);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** changes lower bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
SCIP_RETCODE SCIPchgVarLbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarLbNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( node == NULL )
   {
      SCIP_CALL( SCIPchgVarLb(scip, var, newbound) );
   }  
   else
   {
      SCIPvarAdjustLb(var, scip->set, &newbound);
      SCIP_CALL( SCIPnodeAddBoundchg(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
            scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
   }

   return SCIP_OKAY;
}

/** changes upper bound of variable in the given node; if possible, adjust bound to integral value; doesn't store any
 *  inference information in the bound change, such that in conflict analysis, this change is treated like a branching
 *  decision
 */
SCIP_RETCODE SCIPchgVarUbNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node to change bound at, or NULL for current node */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarUbNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( node == NULL )
   {
      SCIP_CALL( SCIPchgVarUb(scip, var, newbound) );
   }
   else
   {
      SCIPvarAdjustUb(var, scip->set, &newbound);
      SCIP_CALL( SCIPnodeAddBoundchg(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
            scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
   }
   
   return SCIP_OKAY;
}

/** changes global lower bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarLbGlobal", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarAdjustLb(var, scip->set, &newbound);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
            scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** changes global upper bound of variable; if possible, adjust bound to integral value; also tightens the local bound,
 *  if the global bound is better than the local bound
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPchgVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarUbGlobal", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   SCIPvarAdjustUb(var, scip->set, &newbound);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_Bool infeasible;

            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, &infeasible) );
            assert(!infeasible);
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
            scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** changes lazy lower bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  lazy bounds are bounds, that are enforced by constraints and the objective function; hence, these bounds do not need
 *  to be put into the LP explicitly.
 *
 *  @note lazy bounds are useful for branch-and-price since the the corresponding variable bounds are not part of the LP
 */
SCIP_RETCODE SCIPchgVarLbLazy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lazylb              /**< the lazy lower bound to be set */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPchgVarLbLazy", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarChgLbLazy(var, scip->set, lazylb) );

   return SCIP_OKAY;
}

/** changes lazy upper bound of the variable, this is only possible if the variable is not in the LP yet
 *
 *  lazy bounds are bounds, that are enforced by constraints and the objective function; hence, these bounds do not need
 *  to be put into the LP explicitly.
 *
 *  @note lazy bounds are useful for branch-and-price since the the corresponding variable bounds are not part of the LP
 */
SCIP_RETCODE SCIPchgVarUbLazy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             lazyub              /**< the lazy lower bound to be set */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPchgVarUbLazy", FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarChgUbLazy(var, scip->set, lazyub) );

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtightenVarLb", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPcomputeVarLbLocal(scip, var);
   ub = SCIPcomputeVarUbLocal(scip, var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( (force && SCIPsetIsLE(scip->set, newbound, lb)) || (!force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  doesn't store any inference information in the bound change, such that in conflict analysis, this change
 *  is treated like a branching decision
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtightenVarUb", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPcomputeVarLbLocal(scip, var);
   ub = SCIPcomputeVarUbLocal(scip, var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( (force && SCIPsetIsGE(scip->set, newbound, ub)) || (!force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarLbCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPinferVarLbCons", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( (force && SCIPsetIsLE(scip->set, newbound, lb)) || (!force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;
   
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER,
            infercons, NULL, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarUbCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_CONS*            infercons,          /**< constraint that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPinferVarUbCons", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( (force && SCIPsetIsGE(scip->set, newbound, ub)) || (!force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub)) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER,
            infercons, NULL, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference constraint is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 */
SCIP_RETCODE SCIPinferBinvarCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to fix */
   SCIP_Bool             fixedval,           /**< value to fix binary variable to */
   SCIP_CONS*            infercons,          /**< constraint that deduced the fixing */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(SCIPvarIsBinary(var));
   assert(fixedval == TRUE || fixedval == FALSE);
   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPinferBinvarCons", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsEQ(scip->set, lb, 0.0) || SCIPsetIsEQ(scip->set, lb, 1.0));
   assert(SCIPsetIsEQ(scip->set, ub, 0.0) || SCIPsetIsEQ(scip->set, ub, 1.0));
   assert(SCIPsetIsLE(scip->set, lb, ub));

   /* check, if variable is already fixed */
   if( (lb > 0.5) || (ub < 0.5) )
   {
      *infeasible = (fixedval == (lb < 0.5));

      return SCIP_OKAY;
   }

   /* apply the fixing */
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      {
         SCIP_Bool fixed;

         SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal, scip->tree,
               scip->lp, scip->branchcand, scip->eventqueue, (SCIP_Real)fixedval, infeasible, &fixed) );
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 1.0, SCIP_BOUNDTYPE_LOWER,
               infercons, NULL, inferinfo, FALSE) );
      }
      else
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 0.0, SCIP_BOUNDTYPE_UPPER,
               infercons, NULL, inferinfo, FALSE) );
      }
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarLbProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPinferVarLbProp", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( !force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub) )
      return SCIP_OKAY;
   
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER,
            NULL, inferprop, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
   
   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current bound; if possible, adjusts bound to integral value;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason
 *  for the deduction of the bound change
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPinferVarUbProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the bound change */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the bound change is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPinferVarUbProp", FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( !force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub) )
      return SCIP_OKAY;
   
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER,
            NULL, inferprop, inferinfo, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** depending on SCIP's stage, fixes binary variable in the problem, in preprocessing, or in current node;
 *  the given inference propagator is stored, such that the conflict analysis is able to find out the reason for the
 *  deduction of the fixing
 */
SCIP_RETCODE SCIPinferBinvarProp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< binary variable to fix */
   SCIP_Bool             fixedval,           /**< value to fix binary variable to */
   SCIP_PROP*            inferprop,          /**< propagator that deduced the fixing */
   int                   inferinfo,          /**< user information for inference to help resolving the conflict */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            tightened           /**< pointer to store whether the fixing tightened the local bounds, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(SCIPvarIsBinary(var));
   assert(fixedval == TRUE || fixedval == FALSE);
   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPinferBinvarProp", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   /* get current bounds */
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   assert(SCIPsetIsEQ(scip->set, lb, 0.0) || SCIPsetIsEQ(scip->set, lb, 1.0));
   assert(SCIPsetIsEQ(scip->set, ub, 0.0) || SCIPsetIsEQ(scip->set, ub, 1.0));
   assert(SCIPsetIsLE(scip->set, lb, ub));

   /* check, if variable is already fixed */
   if( (lb > 0.5) || (ub < 0.5) )
   {
      *infeasible = (fixedval == (lb < 0.5));

      return SCIP_OKAY;
   }

   /* apply the fixing */
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      {
         SCIP_Bool fixed;

         SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal, scip->tree,
               scip->lp, scip->branchcand, scip->eventqueue, (SCIP_Real)fixedval, infeasible, &fixed) );
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      if( fixedval == TRUE )
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 1.0, SCIP_BOUNDTYPE_LOWER,
               NULL, inferprop, inferinfo, FALSE) );
      }
      else
      {
         SCIP_CALL( SCIPnodeAddBoundinfer(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, 0.0, SCIP_BOUNDTYPE_UPPER,
               NULL, inferprop, inferinfo, FALSE) );
      }
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes global lower bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current global bound; if possible, adjusts bound to integral value;
 *  also tightens the local bound, if the global bound is better than the local bound
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtightenVarLbGlobal", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustLb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasGT(scip->set, newbound, ub) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MIN(newbound, ub);

   if( !force && !SCIPsetIsLbBetter(scip->set, newbound, lb, ub) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgLbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgLbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
            scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   /* coverity: unreachable code */
   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/** changes global upper bound of variable in preprocessing or in the current node, if the new bound is tighter
 *  (w.r.t. bound strengthening epsilon) than the current global bound; if possible, adjusts bound to integral value;
 *  also tightens the local bound, if the global bound is better than the local bound
 *
 *  @note If SCIP is in presolving stage, it can happen that the internal variable array (which get be accessed via
 *        SCIPgetVars()) gets resorted.
 *
 *  @note During presolving, an integer variable which bound changes to {0,1} is upgraded to a binary variable.
 */
SCIP_RETCODE SCIPtightenVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_Bool             force,              /**< force tightening even if below bound strengthening tolerance */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the new domain is empty */
   SCIP_Bool*            tightened           /**< pointer to store whether the bound was tightened, or NULL */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(infeasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtightenVarUbGlobal", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( tightened != NULL )
      *tightened = FALSE;

   SCIPvarAdjustUb(var, scip->set, &newbound);

   /* get current bounds */
   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);
   assert(SCIPsetIsLE(scip->set, lb, ub));

   if( SCIPsetIsFeasLT(scip->set, newbound, lb) )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }
   newbound = MAX(newbound, lb);

   if( !force && !SCIPsetIsUbBetter(scip->set, newbound, lb, ub) )
      return SCIP_OKAY;

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbLocal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      SCIP_CALL( SCIPvarChgUbOriginal(var, scip->set, newbound) );
      break;

   case SCIP_STAGE_TRANSFORMING:
      SCIP_CALL( SCIPvarChgUbGlobal(var, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, newbound) );
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPinProbing(scip) )
      {
         assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);
         assert(scip->tree->root == SCIPtreeGetCurrentNode(scip->tree));

         SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );

         if( var->vartype == SCIP_VARTYPE_INTEGER && SCIPvarIsBinary(var) )
         {
            SCIP_CALL( SCIPchgVarType(scip, var, SCIP_VARTYPE_BINARY, infeasible) );
            assert(!(*infeasible));
         }
         break;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPnodeAddBoundchg(scip->tree->root, scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
            scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, FALSE) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   /* coverity: unreachable code */
   if( tightened != NULL )
      *tightened = TRUE;

   return SCIP_OKAY;
}

/* some simple variable functions implemented as defines */
#undef SCIPcomputeVarLbGlobal
#undef SCIPcomputeVarUbGlobal
#undef SCIPcomputeVarLbLocal
#undef SCIPcomputeVarUbLocal

/** for a multi-aggregated variable, returns the global lower bound computed by adding the global bounds from all aggregation variables
 * this global bound may be tighter than the one given by SCIPvarGetLbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbGlobal
 */
SCIP_Real SCIPcomputeVarLbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcomputeVarLbGlobal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrLbGlobal(var, scip->set);
   else
      return SCIPvarGetLbGlobal(var);
}

/** for a multi-aggregated variable, returns the global upper bound computed by adding the global bounds from all aggregation variables
 * this global bound may be tighter than the one given by SCIPvarGetUbGlobal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbGlobal
 */
SCIP_Real SCIPcomputeVarUbGlobal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcomputeVarUbGlobal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrUbGlobal(var, scip->set);
   else
      return SCIPvarGetUbGlobal(var);
}

/** for a multi-aggregated variable, returns the local lower bound computed by adding the local bounds from all aggregation variables
 * this local bound may be tighter than the one given by SCIPvarGetLbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetLbLocal
 */
SCIP_Real SCIPcomputeVarLbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcomputeVarLbLocal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrLbLocal(var, scip->set);
   else
      return SCIPvarGetLbLocal(var);
}

/** for a multi-aggregated variable, returns the local upper bound computed by adding the local bounds from all aggregation variables
 * this local bound may be tighter than the one given by SCIPvarGetUbLocal, since the latter is not updated if bounds of aggregation variables are changing
 * calling this function for a non-multi-aggregated variable results in a call to SCIPvarGetUbLocal
 */
SCIP_Real SCIPcomputeVarUbLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to compute the bound for */
   )
{
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcomputeVarUbLocal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
      return SCIPvarGetMultaggrUbLocal(var, scip->set);
   else
      return SCIPvarGetUbLocal(var);
}

/** returns solution value and index of variable lower bound that is closest to the variable's value in the given primal
 *  solution or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is
 *  available
 */
SCIP_RETCODE SCIPgetVarClosestVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvlb,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest variable lower bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetVarClosestVlb", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarGetClosestVlb(var, sol, scip->set, scip->stat, closestvlb, closestvlbidx);

   return SCIP_OKAY;
}

/** returns solution value and index of variable upper bound that is closest to the variable's value in the given primal solution;
 *  or current LP solution if no primal solution is given; returns an index of -1 if no variable upper bound is available
 */
SCIP_RETCODE SCIPgetVarClosestVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< active problem variable */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for LP solution */
   SCIP_Real*            closestvub,         /**< pointer to store the value of the closest variable lower bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest variable lower bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetVarClosestVub", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarGetClosestVub(var, sol, scip->set, scip->stat, closestvub, closestvubidx);

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable lower bound x >= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  improves the global bounds of the variable and the vlb variable if possible
 */
SCIP_RETCODE SCIPaddVarVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VAR*             vlbvar,             /**< variable z    in x >= b*z + d */
   SCIP_Real             vlbcoef,            /**< coefficient b in x >= b*z + d */
   SCIP_Real             vlbconstant,        /**< constant d    in x >= b*z + d */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarVlb", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarAddVlb(var, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
         scip->branchcand, scip->eventqueue, vlbvar, vlbcoef, vlbconstant, TRUE, infeasible, nbdchgs) );

   return SCIP_OKAY;
}

/** informs variable x about a globally valid variable upper bound x <= b*z + d with integer variable z;
 *  if z is binary, the corresponding valid implication for z is also added;
 *  improves the global bounds of the variable and the vlb variable if possible
 */
SCIP_RETCODE SCIPaddVarVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_VAR*             vubvar,             /**< variable z    in x <= b*z + d */
   SCIP_Real             vubcoef,            /**< coefficient b in x <= b*z + d */
   SCIP_Real             vubconstant,        /**< constant d    in x <= b*z + d */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarVub", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarAddVub(var, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
         scip->branchcand, scip->eventqueue, vubvar, vubcoef, vubconstant, TRUE, infeasible, nbdchgs) );

   return SCIP_OKAY;
}

/** informs binary variable x about a globally valid implication:  x == 0 or x == 1  ==>  y <= b  or  y >= b;
 *  also adds the corresponding implication or variable bound to the implied variable;
 *  if the implication is conflicting, the variable is fixed to the opposite value;
 *  if the variable is already fixed to the given value, the implication is performed immediately;
 *  if the implication is redundant with respect to the variables' global bounds, it is ignored
 */
SCIP_RETCODE SCIPaddVarImplication(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Bool             varfixing,          /**< FALSE if y should be added in implications for x <= 0, TRUE for x >= 1 */
   SCIP_VAR*             implvar,            /**< variable y in implication y <= b or y >= b */
   SCIP_BOUNDTYPE        impltype,           /**< type       of implication y <= b (SCIP_BOUNDTYPE_UPPER)
                                              *                          or y >= b (SCIP_BOUNDTYPE_LOWER) */
   SCIP_Real             implbound,          /**< bound b    in implication y <= b or y >= b */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarImplication", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY )
   {
      SCIPerrorMessage("can't add implication for nonbinary variable\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPvarAddImplic(var, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
         scip->branchcand, scip->eventqueue, varfixing, implvar, impltype, implbound, TRUE, infeasible, nbdchgs) );

   return SCIP_OKAY;
}

/** adds a clique information to SCIP, stating that at most one of the given binary variables can be set to 1;
 *  if a variable appears twice in the same clique, the corresponding implications are performed
 */
SCIP_RETCODE SCIPaddClique(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   SCIP_Bool*            values,             /**< values of the variables in the clique; NULL to use TRUE for all vars */
   int                   nvars,              /**< number of variables in the clique */
   SCIP_Bool*            infeasible,         /**< pointer to store whether an infeasibility was detected */
   int*                  nbdchgs             /**< pointer to store the number of performed bound changes, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddClique", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   if( nbdchgs != NULL )
      *nbdchgs = 0;

   if( nvars == 2 )
   {
      SCIP_Bool val0;
      SCIP_Bool val1;

      assert(vars != NULL);
      if( values == NULL )
      {
         val0 = TRUE;
         val1 = TRUE;
      }
      else
      {
         val0 = values[0];
         val1 = values[1];
      }
      
      /* add the implications instead of the clique */
      if( SCIPvarGetType(vars[0]) == SCIP_VARTYPE_BINARY )
      {
         /* this function call adds the implication form vars[0] to vars[1] as well as the implication from vars[1] to
          * vars[0] if vars[1] in of binary type
          */
         SCIP_CALL( SCIPvarAddImplic(vars[0], scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
               scip->branchcand, scip->eventqueue, val0, vars[1], val1 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER,
               val1 ? 0.0 : 1.0, TRUE, infeasible, nbdchgs) );
      }
      else if( SCIPvarGetType(vars[1]) == SCIP_VARTYPE_BINARY )
      {
         /* this function call adds the implication form vars[1] to vars[0] */
         SCIP_CALL( SCIPvarAddImplic(vars[1], scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
               scip->branchcand, scip->eventqueue, val1, vars[0], val0 ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER,
               val0 ? 0.0 : 1.0, TRUE, infeasible, nbdchgs) );
      }
      /** in case one both variables are not of binary type we have to add the implication as variable bounds */
      else
      {
         /* both variables are not of binary type but are implicit binary; in that case we can only add the this
          * implication as variable bounds
          */
         assert(SCIPvarGetType(vars[0]) != SCIP_VARTYPE_BINARY && SCIPvarIsBinary(vars[0]));
         assert(SCIPvarGetType(vars[1]) != SCIP_VARTYPE_BINARY && SCIPvarIsBinary(vars[1]));

         /* add variable upper or rather variable lower bound on vars[0] */
         if( val0 )
         {
            SCIP_CALL( SCIPvarAddVub(vars[0], scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
                  scip->branchcand, scip->eventqueue, vars[1], val1 ? -1.0 : 1.0, val1 ? 1.0 : 0.0, TRUE, infeasible, nbdchgs) );
         }
         else
         {
            SCIP_CALL( SCIPvarAddVlb(vars[0], scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
                  scip->branchcand, scip->eventqueue, vars[1], val1 ? 1.0 : -1.0, val1 ? 0.0 : 1.0, TRUE, infeasible, nbdchgs) );
         }
         
         /* add variable upper or rather variable lower bound on vars[1] */
         if( val1 )
         {
            SCIP_CALL( SCIPvarAddVub(vars[1], scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
                  scip->branchcand, scip->eventqueue, vars[0], val0 ? -1.0 : 1.0, val0 ? 1.0 : 0.0, TRUE, infeasible, nbdchgs) );
         }
         else
         {
            SCIP_CALL( SCIPvarAddVlb(vars[1], scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->cliquetable,
                  scip->branchcand, scip->eventqueue, vars[0], val0 ? 1.0 : -1.0, val0 ? 0.0 : 1.0, TRUE, infeasible, nbdchgs) );
         }
      }
   }
   else if( nvars >= 3 )
   {
      /* add the clique to the clique table */
      SCIP_CALL( SCIPcliquetableAdd(scip->cliquetable, scip->mem->probmem, scip->set, scip->stat, scip->lp,
            scip->branchcand, scip->eventqueue, vars, values, nvars, infeasible, nbdchgs) );
   }

   return SCIP_OKAY;
}

/* calculate clique partition for a maximal amount of comparisons on variables due to expensive algorithm
 * @todo: check for a good value, maybe it's better to check parts of variables
 */
#define MAXNCLIQUEVARSCOMP 1000000

/** calculates a partition of the given set of binary variables into cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same clique;
 *  the first variable is always assigned to clique 0, and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique at most 1 variables can be set to TRUE in a feasible solution;
 */
SCIP_RETCODE SCIPcalcCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   )
{
   SCIP_VAR** tmpvars;
   SCIP_VAR** cliquevars;
   SCIP_Bool* cliquevalues;
   SCIP_Bool* tmpvalues;
   int ncliquevars;
   int i;
   int maxncliquevarscomp;

   assert(scip != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || cliquepartition != NULL);
   assert(ncliques != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcalcCliquePartition", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( nvars == 0 )
   {
      *ncliques = 0;
      return SCIP_OKAY;
   }

   /* allocate temporary memory for storing the variables of the current clique */
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &cliquevars, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &cliquevalues, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &tmpvalues, nvars) );
   SCIP_CALL( SCIPsetDuplicateBufferArray(scip->set, &tmpvars, vars, nvars) );
   ncliquevars = 0;

   /* initialize the cliquepartition array with -1 */
   /* initialize the tmpvalues array */
   for( i = nvars - 1; i >= 0; --i )
   {  
      tmpvalues[i] = TRUE;
      cliquepartition[i] = -1;
   }

   /* get corresponding active problem variables */
   SCIP_CALL( SCIPvarsGetProbvarBinary(&tmpvars, &tmpvalues, nvars) );

   maxncliquevarscomp = MIN(nvars*nvars, MAXNCLIQUEVARSCOMP);

   /* calculate the clique partition */
   *ncliques = 0;
   for( i = 0; i < nvars; ++i )
   {
      if( cliquepartition[i] == -1 )
      {
         int j;

         /* variable starts a new clique */
         cliquepartition[i] = *ncliques;
         cliquevars[0] = tmpvars[i];
         cliquevalues[0] = tmpvalues[i];
         ncliquevars = 1;

         /* if variable is not active (multi-aggregated or fixed), it cannot be in any clique */
         if( SCIPvarIsActive(tmpvars[i]) )
         {
            /* greedily fill up the clique */
            for( j = i+1; j < nvars; ++j )
            {
               /* if variable is not active (multi-aggregated or fixed), it cannot be in any clique */
               if( cliquepartition[j] == -1 && SCIPvarIsActive(tmpvars[j]) )
               {
                  int k;

                  /* check if every variable in the actual clique is in clique with the new variable */ 
                  for( k = ncliquevars - 1; k >= 0; --k )
                  {
                     if( !SCIPvarsHaveCommonClique(tmpvars[j], tmpvalues[j], cliquevars[k], cliquevalues[k], TRUE) )
                        break;
                  }
                                    
                  if( k == -1 )
                  {
                     /* put the variable into the same clique */
                     cliquepartition[j] = cliquepartition[i];
                     cliquevars[ncliquevars] = tmpvars[j];
                     cliquevalues[ncliquevars] = tmpvalues[j];
                     ++ncliquevars;
                  }
               }
            }
         }

         /* this clique is finished */
         ++(*ncliques);
      }
      assert(cliquepartition[i] >= 0 && cliquepartition[i] < i+1);

      /* break if we reached the maximal number of comparisons */
      if( i * nvars > maxncliquevarscomp )
         break;
   }
   /* if we had to much variables fill up the cliquepartition and put each variable in a separate clique */
   for( ; i < nvars; ++i )
   {
      if( cliquepartition[i] == -1 )
      {
         cliquepartition[i] = *ncliques;
         ++(*ncliques);
      }
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(scip->set, &tmpvars);
   SCIPsetFreeBufferArray(scip->set, &tmpvalues);
   SCIPsetFreeBufferArray(scip->set, &cliquevalues);
   SCIPsetFreeBufferArray(scip->set, &cliquevars);

   return SCIP_OKAY;
}

/** calculates a partition of the given set of binary variables into negated cliques;
 *  afterwards the output array contains one value for each variable, such that two variables got the same value iff they
 *  were assigned to the same negated clique;
 *  the first variable is always assigned to clique 0 and a variable can only be assigned to clique i if at least one of
 *  the preceding variables was assigned to clique i-1;
 *  for each clique with n_c variables at least n_c-1 variables can be set to TRUE in a feasible solution;
 */
SCIP_RETCODE SCIPcalcNegatedCliquePartition(
   SCIP*const            scip,               /**< SCIP data structure */
   SCIP_VAR**const       vars,               /**< binary variables in the clique from which at most one can be set to 1 */
   int const             nvars,              /**< number of variables in the clique */
   int*const             cliquepartition,    /**< array of length nvars to store the clique partition */
   int*const             ncliques            /**< pointer to store the number of cliques actually contained in the partition */
   )
{
   SCIP_VAR** negvars;
   int v;

   assert(scip != NULL);
   assert(cliquepartition != NULL || nvars == 0);
   assert(ncliques != NULL);

   if( nvars == 0 )
   {
      *ncliques = 0;
      return SCIP_OKAY;
   }
   assert(vars != NULL);
   
   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(scip->set, &negvars, nvars) );
   
   /* get all negated variables */
   for( v = nvars - 1; v >= 0; --v )
   {
      SCIP_CALL( SCIPgetNegatedVar(scip, vars[v], &(negvars[v])) );
   }

   /* calculate cliques on negated variables, which are "negated" cliques on normal variables array */
   SCIP_CALL( SCIPcalcCliquePartition( scip, negvars, nvars, cliquepartition, ncliques) );

   /* free temporary memory */
   SCIPsetFreeBufferArray(scip->set, &negvars);
   
   return SCIP_OKAY;
}

/** gets the number of cliques in the clique table */
int SCIPgetNCliques(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNCliques", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPcliquetableGetNCliques(scip->cliquetable);
}

/** gets the array of cliques in the clique table */
SCIP_CLIQUE** SCIPgetCliques(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetCliques", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPcliquetableGetCliques(scip->cliquetable);
}

/** sets the branch factor of the variable; this value can be used in the branching methods to scale the score
 *  values of the variables; higher factor leads to a higher probability that this variable is chosen for branching
 */
SCIP_RETCODE SCIPchgVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             branchfactor        /**< factor to weigh variable's branching score with */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchFactor(var, scip->set, branchfactor);

   return SCIP_OKAY;
}

/** scales the branch factor of the variable with the given value */
SCIP_RETCODE SCIPscaleVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             scale               /**< factor to scale variable's branching factor with */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPscaleVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchFactor(var, scip->set, scale * SCIPvarGetBranchFactor(var));

   return SCIP_OKAY;
}

/** adds the given value to the branch factor of the variable */
SCIP_RETCODE SCIPaddVarBranchFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             addfactor           /**< value to add to the branch factor of the variable */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarBranchFactor", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchFactor(var, scip->set, addfactor + SCIPvarGetBranchFactor(var));

   return SCIP_OKAY;
}

/** sets the branch priority of the variable; variables with higher branch priority are always preferred to variables
 *  with lower priority in selection of branching variable
 *
 * @note the default branching priority is 0
 */
SCIP_RETCODE SCIPchgVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< branch priority of the variable */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchPriority(var, branchpriority);

   return SCIP_OKAY;
}

/** changes the branch priority of the variable to the given value, if it is larger than the current priority */
SCIP_RETCODE SCIPupdateVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   branchpriority      /**< new branch priority of the variable, if it is larger than current priority */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPupdateVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( branchpriority > SCIPvarGetBranchPriority(var) )
      SCIPvarChgBranchPriority(var, branchpriority);

   return SCIP_OKAY;
}

/** adds the given value to the branch priority of the variable */
SCIP_RETCODE SCIPaddVarBranchPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   int                   addpriority         /**< value to add to the branch priority of the variable */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarBranchPriority", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchPriority(var, addpriority + SCIPvarGetBranchPriority(var));

   return SCIP_OKAY;
}

/** sets the branch direction of the variable (-1: prefer downwards branch, 0: automatic selection, +1: prefer upwards
 *  branch)
 */
SCIP_RETCODE SCIPchgVarBranchDirection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        branchdirection     /**< preferred branch direction of the variable (downwards, upwards, auto) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarBranchDirection", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPvarChgBranchDirection(var, branchdirection);

   return SCIP_OKAY;
}

/** tightens the variable bounds due a new variable type */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_VARTYPE          vartype,            /**< new type of variable */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected (, due to
                                              *   integrality condition of the new variable type) */
   )
{
   assert(scip != NULL);
   assert(SCIPgetStage(scip) == SCIP_STAGE_PROBLEM || SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING);
   assert(scip->set->stage == SCIP_STAGE_PROBLEM || SCIPvarIsTransformed(var));

   *infeasible = FALSE;

   /* adjusts bounds if the variable type changed form continuous to non-continuous (integral) */
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS && vartype != SCIP_VARTYPE_CONTINUOUS )
   {
      SCIP_Bool tightened;

      /* we adjust variable bounds to integers first, since otherwise a later bound tightening with a fractional old
       * bound may give an assert because SCIP expects non-continuous variables to have non-fractional bounds
       */
      if( !SCIPisFeasIntegral(scip, SCIPvarGetLbGlobal(var)) )
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, SCIPfeasCeil(scip, SCIPvarGetLbGlobal(var)), TRUE, infeasible, &tightened) );
         if( *infeasible )
            return SCIP_OKAY;

         assert(tightened);
      }
      if( !SCIPisFeasIntegral(scip, SCIPvarGetUbGlobal(var)) )
      {
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, SCIPfeasFloor(scip, SCIPvarGetUbGlobal(var)), TRUE, infeasible, &tightened) );
         if( *infeasible )
            return SCIP_OKAY;

         assert(tightened);
      }
   }

   return SCIP_OKAY;
}

/** changes type of variable in the problem;
 *
 * @note this type change might change the variable array returned from SCIPgetVars() and SCIPgetVarsData();
 *
 * @note if SCIP is already beyond the SCIP_STAGE_PROBLEM and a original variable is passed, the variable type of the
 *       corresponding transformed variable is changed; the type of the original variable does not change
 *
 * @note if the type changes from a continuous variable to a non-continuous variable the bounds of the variable get
 *       adjusted w.r.t. to integrality information
 */
SCIP_RETCODE SCIPchgVarType(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_VARTYPE          vartype,            /**< new type of variable */
   SCIP_Bool*            infeasible          /**< pointer to store whether an infeasibility was detected (, due to
                                              *   integrality condition of the new variable type) */
   )
{
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPchgVarType", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* change variable type */
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));

      /* first adjust the variable due new integrality information */
      SCIP_CALL( tightenBounds(scip, var, vartype, infeasible) );

      if( *infeasible )
         return SCIP_OKAY;

      /* second change variable type */
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         SCIP_CALL( SCIPprobChgVarType(scip->origprob, scip->mem->probmem, scip->set, scip->branchcand, var, vartype) );
      }
      else
      {
         SCIP_CALL( SCIPvarChgType(var, vartype) );
      }
      break;

   case SCIP_STAGE_PRESOLVING:
      if( !SCIPvarIsTransformed(var) )
      {
         SCIP_VAR* transvar;

         SCIP_CALL( SCIPgetTransformedVar(scip, var, &transvar) );
         assert(transvar != NULL);

         /* recall method with transformed variable */
         SCIP_CALL( SCIPchgVarType(scip, transvar, vartype, infeasible) );
         return SCIP_OKAY;
      }

      /* first adjust the variable due new integrality information */
      SCIP_CALL( tightenBounds(scip, var, vartype, infeasible) );

      if( *infeasible )
         return SCIP_OKAY;

      /* second change variable type */
      if( SCIPvarGetProbindex(var) >= 0 )
      {
         SCIP_CALL( SCIPprobChgVarType(scip->transprob, scip->mem->probmem, scip->set, scip->branchcand, var, vartype) );
      }
      else
      {
         SCIP_CALL( SCIPvarChgType(var, vartype) );
      }
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** in problem creation and solving stage, both bounds of the variable are set to the given value;
 *  in presolving stage, the variable is converted into a fixed variable, and bounds are changed respectively;
 *  conversion into a fixed variable changes the vars array returned from SCIPgetVars() and SCIPgetVarsData(),
 *  and also renders arrays returned from the SCIPvarGetImpl...() methods invalid
 */
SCIP_RETCODE SCIPfixVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to fix */
   SCIP_Real             fixedval,           /**< value to fix variable to */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the fixing is infeasible */
   SCIP_Bool*            fixed               /**< pointer to store whether the fixing was performed (variable was unfixed) */
   )
{
   assert(var != NULL);
   assert(infeasible != NULL);
   assert(fixed != NULL);

   SCIP_CALL( checkStage(scip, "SCIPfixVar", FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   *fixed = FALSE;

   /* in the problem creation stage, modify the bounds as requested, independently from the current bounds */
   if( scip->set->stage != SCIP_STAGE_PROBLEM )
   {
      if( (SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS && !SCIPsetIsFeasIntegral(scip->set, fixedval))
         || SCIPsetIsFeasLT(scip->set, fixedval, SCIPvarGetLbLocal(var))
         || SCIPsetIsFeasGT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
      {
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
      else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
      {
         *infeasible = !SCIPsetIsFeasEQ(scip->set, fixedval, SCIPvarGetLbLocal(var));
         return SCIP_OKAY;
      }
   }
   else
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      /* in the problem creation stage, modify the bounds as requested, independently from the current bounds;
       * we have to make sure, that the order of the bound changes does not intermediately produce an invalid
       * interval lb > ub
       */
      if( fixedval <= SCIPvarGetLbLocal(var) )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, fixedval) );
         SCIP_CALL( SCIPchgVarUb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, fixedval) );
         SCIP_CALL( SCIPchgVarLb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      return SCIP_OKAY;

   case SCIP_STAGE_PRESOLVING:
      if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      {
         SCIP_CALL( SCIPvarFix(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal, scip->tree,
               scip->lp, scip->branchcand, scip->eventqueue, fixedval, infeasible, fixed) );
         return SCIP_OKAY;
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      if( SCIPsetIsFeasGT(scip->set, fixedval, SCIPvarGetLbLocal(var)) )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      if( SCIPsetIsFeasLT(scip->set, fixedval, SCIPvarGetUbLocal(var)) )
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, fixedval) );
         *fixed = TRUE;
      }
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** From a given equality a*x + b*y == c, aggregates one of the variables and removes it from the set of
 *  active problem variables. This changes the vars array returned from SCIPgetVars() and SCIPgetVarsData(),
 *  and also renders the arrays returned from the SCIPvarGetImpl...() methods for the two variables invalid.
 *  In the first step, the equality is transformed into an equality with active problem variables
 *  a'*x' + b'*y' == c'. If x' == y', this leads to the detection of redundancy if a' == -b' and c' == 0,
 *  of infeasibility, if a' == -b' and c' != 0, or to a variable fixing x' == c'/(a'+b') (and possible
 *  infeasibility) otherwise.
 *  In the second step, the variable to be aggregated is chosen among x' and y', prefering a less strict variable
 *  type as aggregation variable (i.e. continuous variables are preferred over implicit integers, implicit integers
 *  over integers, and integers over binaries). If none of the variables is continuous, it is tried to find an integer
 *  aggregation (i.e. integral coefficients a'' and b'', such that a''*x' + b''*y' == c''). This can lead to
 *  the detection of infeasibility (e.g. if c'' is fractional), or to a rejection of the aggregation (denoted by
 *  aggregated == FALSE), if the resulting integer coefficients are too large and thus numerically instable.
 *
 *  The output flags have the following meaning:
 *  - infeasible: the problem is infeasible
 *  - redundant:  the equality can be deleted from the constraint set
 *  - aggregated: the aggregation was successfully performed (the variables were not aggregated before)
 */
SCIP_RETCODE SCIPaggregateVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             varx,               /**< variable x in equality a*x + b*y == c */
   SCIP_VAR*             vary,               /**< variable y in equality a*x + b*y == c */
   SCIP_Real             scalarx,            /**< multiplier a in equality a*x + b*y == c */
   SCIP_Real             scalary,            /**< multiplier b in equality a*x + b*y == c */
   SCIP_Real             rhs,                /**< right hand side c in equality a*x + b*y == c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            redundant,          /**< pointer to store whether the equality is (now) redundant */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_Real constantx;
   SCIP_Real constanty;

   assert(infeasible != NULL);
   assert(redundant != NULL);
   assert(aggregated != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaggregateVars", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   *infeasible = FALSE;
   *redundant = FALSE;
   *aggregated = FALSE;

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot aggregate variables during probing\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);

   /* do not perform aggregation if it is globally deactivated */
   if( scip->set->presol_donotaggr )
      return SCIP_OKAY;

   /* get the corresponding equality in active problem variable space:
    * transform both expressions "a*x + 0" and "b*y + 0" into problem variable space
    */
   constantx = 0.0;
   constanty = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&varx, &scalarx, &constantx) );
   SCIP_CALL( SCIPvarGetProbvarSum(&vary, &scalary, &constanty) );

   /* we cannot aggregate multi-aggregated variables */
   if( SCIPvarGetStatus(varx) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(vary) == SCIP_VARSTATUS_MULTAGGR )
      return SCIP_OKAY;

   /* move the constant to the right hand side to acquire the form "a'*x' + b'*y' == c'" */
   rhs -= (constantx + constanty);

   /* if a scalar is zero, treat the variable as fixed-to-zero variable */
   if( SCIPsetIsZero(scip->set, scalarx) )
      varx = NULL;
   if( SCIPsetIsZero(scip->set, scalary) )
      vary = NULL;

   /* capture the special cases that less than two variables are left, due to resolutions to a fixed variable or
    * to the same active variable
    */
   if( varx == NULL && vary == NULL )
   {
      /* both variables were resolved to fixed variables */
      *infeasible = !SCIPsetIsZero(scip->set, rhs);
      *redundant = TRUE;
   }
   else if( varx == NULL )
   {
      assert(SCIPsetIsZero(scip->set, scalarx));
      assert(!SCIPsetIsZero(scip->set, scalary));

      /* variable x was resolved to fixed variable: variable y can be fixed to c'/b' */
      SCIP_CALL( SCIPvarFix(vary, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal, scip->tree,
            scip->lp, scip->branchcand, scip->eventqueue, rhs/scalary, infeasible, aggregated) );
      *redundant = TRUE;
   }
   else if( vary == NULL )
   {
      assert(SCIPsetIsZero(scip->set, scalary));
      assert(!SCIPsetIsZero(scip->set, scalarx));
      
      /* variable y was resolved to fixed variable: variable x can be fixed to c'/a' */
      SCIP_CALL( SCIPvarFix(varx, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal, scip->tree,
            scip->lp, scip->branchcand, scip->eventqueue, rhs/scalarx, infeasible, aggregated) );
      *redundant = TRUE;
   }
   else if( varx == vary )
   {
      /* both variables were resolved to the same active problem variable: this variable can be fixed */
      scalarx += scalary;
      if( SCIPsetIsZero(scip->set, scalarx) )
      {
         /* left hand side of equality is zero: equality is potentially infeasible */
         *infeasible = !SCIPsetIsZero(scip->set, rhs);
      }
      else
      {
         /* sum of scalars is not zero: fix variable x' == y' to c'/(a'+b') */
         SCIP_CALL( SCIPvarFix(varx, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
               scip->tree, scip->lp, scip->branchcand, scip->eventqueue, rhs/scalarx, infeasible, aggregated) );
      }
      *redundant = TRUE;
   }
   else
   {
      /* both variables are different active problem variables, and both scalars are non-zero: try to aggregate them */
      SCIP_CALL( SCIPvarTryAggregateVars(scip->set, scip->mem->probmem, scip->stat, scip->transprob, scip->primal,
	    scip->tree, scip->lp, scip->cliquetable, scip->branchcand, scip->eventfilter, scip->eventqueue,
	    varx, vary, scalarx, scalary, rhs, infeasible, aggregated) );
      *redundant = *aggregated;
   }

   return SCIP_OKAY;
}

/** converts variable into multi-aggregated variable; this changes the variable array returned from
 *  SCIPgetVars() and SCIPgetVarsData();
 *
 *  @warning The integrality condition is not checked anymore on the multi-aggregated variable. You must not
 *           multi-aggregate an integer variable without being sure, that integrality on the aggregation variables
 *           implies integrality on the aggregated variable.
 *
 *  The output flags have the following meaning:
 *  - infeasible: the problem is infeasible
 *  - aggregated: the aggregation was successfully performed (the variables were not aggregated before)
 */
SCIP_RETCODE SCIPmultiaggregateVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable x to aggregate */
   int                   naggvars,           /**< number n of variables in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_VAR**            aggvars,            /**< variables y_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real*            scalars,            /**< multipliers a_i in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Real             constant,           /**< constant shift c in aggregation x = a_1*y_1 + ... + a_n*y_n + c */
   SCIP_Bool*            infeasible,         /**< pointer to store whether the aggregation is infeasible */
   SCIP_Bool*            aggregated          /**< pointer to store whether the aggregation was successful */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPmultiaggregateVar", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot multi-aggregate variables during probing\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeGetCurrentDepth(scip->tree) == 0);

   SCIP_CALL( SCIPvarMultiaggregate(var, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal,
         scip->tree, scip->lp, scip->cliquetable, scip->branchcand, scip->eventfilter, scip->eventqueue,
         naggvars, aggvars, scalars, constant, infeasible, aggregated) );

   return SCIP_OKAY;
}

/** returns whether aggregation of variables is not allowed */
SCIP_Bool SCIPdoNotAggr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->set->presol_donotaggr;
}

/** returns whether variable is not allowed to be multi-aggregated */
SCIP_Bool SCIPdoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable x to aggregate */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   return scip->set->presol_donotmultaggr || SCIPvarDoNotMultaggr(var);
}

/** marks the variable to not to be multi-aggregated */
SCIP_RETCODE SCIPmarkDoNotMultaggrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to delete */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPmarkDoNotMultiaggrVar", TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarMarkDoNotMultaggr(var) );

   return SCIP_OKAY;
}

/** updates the pseudo costs of the given variable and the global pseudo costs after a change of "solvaldelta" in the
 *  variable's solution value and resulting change of "objdelta" in the in the LP's objective value;
 *  the update is ignored, if the objective value difference is infinite
 */
SCIP_RETCODE SCIPupdateVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta,        /**< difference of variable's new LP value - old LP value */
   SCIP_Real             objdelta,           /**< difference of new LP's objective value - old LP's objective value */
   SCIP_Real             weight              /**< weight in (0,1] of this update in pseudo cost sum */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPupdateVarPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( !SCIPsetIsInfinity(scip->set, 2*objdelta) ) /* differences  infinity - eps  should also be treated as infinity */
   {
      SCIP_CALL( SCIPvarUpdatePseudocost(var, scip->set, scip->stat, solvaldelta, objdelta, weight) );
   }

   return SCIP_OKAY;
}

/** gets the variable's pseudo cost value for the given change of the variable's LP value */
SCIP_Real SCIPgetVarPseudocostVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocostVal", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetPseudocost(var, scip->stat, solvaldelta);
}

/** gets the variable's pseudo cost value for the given change of the variable's LP value,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetVarPseudocostValCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocostValCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
 
   return SCIPvarGetPseudocostCurrentRun(var, scip->stat, solvaldelta);
}

/** gets the variable's pseudo cost value for the given direction */
SCIP_Real SCIPgetVarPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocost", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   return SCIPvarGetPseudocost(var, scip->stat, dir == SCIP_BRANCHDIR_DOWNWARDS ? -1.0 : 1.0);
}

/** gets the variable's pseudo cost value for the given direction,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetVarPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocostCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);
   
   return SCIPvarGetPseudocostCurrentRun(var, scip->stat, dir == SCIP_BRANCHDIR_DOWNWARDS ? -1.0 : 1.0);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction */
SCIP_Real SCIPgetVarPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocostCount", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   return SCIPvarGetPseudocostCount(var, dir);
}

/** gets the variable's (possible fractional) number of pseudo cost updates for the given direction,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetVarPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocostCountCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   assert(dir == SCIP_BRANCHDIR_DOWNWARDS || dir == SCIP_BRANCHDIR_UPWARDS);

   return SCIPvarGetPseudocostCountCurrentRun(var, dir);
}

/** gets the variable's pseudo cost score value for the given LP solution value */
SCIP_Real SCIPgetVarPseudocostScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   )
{
   SCIP_Real downsol;
   SCIP_Real upsol;
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocostScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   downsol = SCIPsetFeasCeil(scip->set, solval-1.0);
   upsol = SCIPsetFeasFloor(scip->set, solval+1.0);
   pscostdown = SCIPvarGetPseudocost(var, scip->stat, downsol-solval);
   pscostup = SCIPvarGetPseudocost(var, scip->stat, upsol-solval);

   return SCIPbranchGetScore(scip->set, var, pscostdown, pscostup);
}

/** gets the variable's pseudo cost score value for the given LP solution value,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetVarPseudocostScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             solval              /**< variable's LP solution value */
   )
{
   SCIP_Real downsol;
   SCIP_Real upsol;
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarPseudocostScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   downsol = SCIPsetFeasCeil(scip->set, solval-1.0);
   upsol = SCIPsetFeasFloor(scip->set, solval+1.0);
   pscostdown = SCIPvarGetPseudocostCurrentRun(var, scip->stat, downsol-solval);
   pscostup = SCIPvarGetPseudocostCurrentRun(var, scip->stat, upsol-solval);

   return SCIPbranchGetScore(scip->set, var, pscostdown, pscostup);
}

/** returns the variable's conflict score value */
SCIP_Real SCIPgetVarVSIDS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarVSIDS", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );
   
   return SCIPvarGetVSIDS(var, scip->stat, dir);
}

/** returns the variable's conflict score value only using conflicts of the current run */
SCIP_Real SCIPgetVarVSIDSCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarVSIDSCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetVSIDSCurrentRun(var, scip->stat, dir);
}

/** returns the variable's conflict score value */
SCIP_Real SCIPgetVarConflictScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarConflictScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   downscore = SCIPvarGetVSIDS(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetVSIDS(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's conflict score value only using conflicts of the current run */
SCIP_Real SCIPgetVarConflictScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarConflictScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   downscore = SCIPvarGetVSIDSCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetVSIDSCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's conflict length score */
SCIP_Real SCIPgetVarConflictlengthScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarConflictlengthScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   downscore = SCIPvarGetAvgConflictlength(var, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetAvgConflictlength(var, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's conflict length score only using conflicts of the current run */
SCIP_Real SCIPgetVarConflictlengthScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real downscore;
   SCIP_Real upscore;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarConflictlengthScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   downscore = SCIPvarGetAvgConflictlengthCurrentRun(var, SCIP_BRANCHDIR_DOWNWARDS);
   upscore = SCIPvarGetAvgConflictlengthCurrentRun(var, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, downscore, upscore);
}

/** returns the variable's average conflict length */
SCIP_Real SCIPgetVarAvgConflictlength(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgConflictlength", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgConflictlength(var, dir);
}

/** returns the variable's average  conflict length only using conflicts of the current run */
SCIP_Real SCIPgetVarAvgConflictlengthCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgConflictlengthCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgConflictlengthCurrentRun(var, dir);
}

/** returns the average number of inferences found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 */
SCIP_Real SCIPgetVarAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgInferences", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgInferences(var, scip->stat, dir);
}

/** returns the average number of inferences found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of inferences
 *  over all variables for branching in the given direction is returned
 */
SCIP_Real SCIPgetVarAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgInferencesCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, dir);
}

/** returns the variable's average inference score value */
SCIP_Real SCIPgetVarAvgInferenceScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real inferdown;
   SCIP_Real inferup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   inferdown = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, inferdown, inferup);
}

/** returns the variable's average inference score value only using inferences of the current run */
SCIP_Real SCIPgetVarAvgInferenceScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real inferdown;
   SCIP_Real inferup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   inferdown = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, inferdown, inferup);
}

/** initializes the upwards and downwards pseudocosts, conflict scores, conflict lengths, inference scores, cutoff scores
 *  of a variable to the given values 
 */
SCIP_RETCODE SCIPinitVarBranchStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which should be initialized */
   SCIP_Real             downpscost,         /**< value to which pseudocosts for downwards branching should be initialized */
   SCIP_Real             uppscost,           /**< value to which pseudocosts for upwards branching should be initialized */
   SCIP_Real             downvsids,          /**< value to which VSIDS score for downwards branching should be initialized */
   SCIP_Real             upvsids,            /**< value to which VSIDS score for upwards branching should be initialized */
   SCIP_Real             downconflen,        /**< value to which conflict length score for downwards branching should be initialized */
   SCIP_Real             upconflen,          /**< value to which conflict length score for upwards branching should be initialized */
   SCIP_Real             downinfer,          /**< value to which inference counter for downwards branching should be initialized */
   SCIP_Real             upinfer,            /**< value to which inference counter for upwards branching should be initialized */
   SCIP_Real             downcutoff,         /**< value to which cutoff counter for downwards branching should be initialized */
   SCIP_Real             upcutoff            /**< value to which cutoff counter for upwards branching should be initialized */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPinitVarBranchStats", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   assert(downpscost >= 0.0 && uppscost >= 0.0);
   assert(downvsids >= 0.0 && upvsids >= 0.0);
   assert(downconflen >= 0.0 && upconflen >= 0.0);
   assert(downinfer >= 0.0 && upinfer >= 0.0);
   assert(downcutoff >= 0.0 && upcutoff >= 0.0);
   
   if( !SCIPisFeasZero(scip, downpscost) || !SCIPisFeasZero(scip, downvsids)
      || !SCIPisFeasZero(scip, downinfer) || !SCIPisFeasZero(scip, downcutoff) )
   {
      SCIP_CALL( SCIPvarIncNBranchings(var, scip->stat, 1, SCIP_BRANCHDIR_DOWNWARDS) );
      SCIP_CALL( SCIPvarUpdatePseudocost(var, scip->set, scip->stat, -1.0, downpscost, 1.0) );  
      SCIP_CALL( SCIPvarIncInferenceSum(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, downinfer) );  
      SCIP_CALL( SCIPvarIncVSIDS(var, SCIP_BRANCHDIR_DOWNWARDS, downvsids) ); 
      SCIP_CALL( SCIPvarIncCutoffSum(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS, downcutoff) ); 
   }
   
   if( !SCIPisFeasZero(scip, downconflen) )
   {
      SCIP_CALL( SCIPvarIncNActiveConflicts(var, SCIP_BRANCHDIR_DOWNWARDS, downconflen) );
   }

   if( !SCIPisFeasZero(scip, uppscost) || !SCIPisFeasZero(scip, upvsids)
      || !SCIPisFeasZero(scip, upinfer) || !SCIPisFeasZero(scip, upcutoff) )
   {
      SCIP_CALL( SCIPvarIncNBranchings(var, scip->stat, 1, SCIP_BRANCHDIR_UPWARDS) );
      SCIP_CALL( SCIPvarUpdatePseudocost(var, scip->set, scip->stat, 1.0, uppscost, 1.0) );  
      SCIP_CALL( SCIPvarIncInferenceSum(var, scip->stat, SCIP_BRANCHDIR_UPWARDS, upinfer) );
      SCIP_CALL( SCIPvarIncVSIDS(var, SCIP_BRANCHDIR_UPWARDS, upvsids) );     
      SCIP_CALL( SCIPvarIncCutoffSum(var, scip->stat, SCIP_BRANCHDIR_UPWARDS, upcutoff) );   
   }
   
   if( !SCIPisFeasZero(scip, upconflen) )
   {
      SCIP_CALL( SCIPvarIncNActiveConflicts(var, SCIP_BRANCHDIR_UPWARDS, upconflen) );
   }

   return SCIP_OKAY;
}

/** returns the average number of cutoffs found after branching on the variable in given direction;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 */
SCIP_Real SCIPgetVarAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgCutoffs(var, scip->stat, dir);
}

/** returns the average number of cutoffs found after branching on the variable in given direction in the current run;
 *  if branching on the variable in the given direction was yet evaluated, the average number of cutoffs
 *  over all variables for branching in the given direction is returned
 */
SCIP_Real SCIPgetVarAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffsCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, dir);
}

/** returns the variable's average cutoff score value */
SCIP_Real SCIPgetVarAvgCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   cutoffdown = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, cutoffdown, cutoffup);
}

/** returns the variable's average cutoff score value, only using cutoffs of the current run */
SCIP_Real SCIPgetVarAvgCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   cutoffdown = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var, cutoffdown, cutoffup);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor
 */
SCIP_Real SCIPgetVarAvgInferenceCutoffScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   SCIP_Real avginferdown;
   SCIP_Real avginferup;
   SCIP_Real avginfer;
   SCIP_Real inferdown;
   SCIP_Real inferup;
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceCutoffScore", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   avginferdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   avginferup = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);
   avginfer = (avginferdown + avginferup)/2.0;
   inferdown = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferences(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);
   cutoffdown = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffs(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var,
      inferdown + cutoffweight * avginfer * cutoffdown, inferup + cutoffweight * avginfer * cutoffup);
}

/** returns the variable's average inference/cutoff score value, weighting the cutoffs of the variable with the given
 *  factor, only using inferences and cutoffs of the current run
 */
SCIP_Real SCIPgetVarAvgInferenceCutoffScoreCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             cutoffweight        /**< factor to weigh average number of cutoffs in branching score */
   )
{
   SCIP_Real avginferdown;
   SCIP_Real avginferup;
   SCIP_Real avginfer;
   SCIP_Real inferdown;
   SCIP_Real inferup;
   SCIP_Real cutoffdown;
   SCIP_Real cutoffup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarAvgInferenceCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   avginferdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   avginferup = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);
   avginfer = (avginferdown + avginferup)/2.0;
   inferdown = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   inferup = SCIPvarGetAvgInferencesCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);
   cutoffdown = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffup = SCIPvarGetAvgCutoffsCurrentRun(var, scip->stat, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, var,
      inferdown + cutoffweight * avginfer * cutoffdown, inferup + cutoffweight * avginfer * cutoffup);
}

/** outputs variable information to file stream */
SCIP_RETCODE SCIPprintVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPprintVar", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPvarPrint(var, scip->set, file);

   return SCIP_OKAY;
}

/*
 * conflict analysis methods
 */

/** initializes the conflict analysis by clearing the conflict candidate queue; this method must be called before
 *  you enter the conflict variables by calling SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(),
 *  or SCIPaddConflictBinvar();
 */
SCIP_RETCODE SCIPinitConflictAnalysis(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPinitConflictAnalysis", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictInit(scip->conflict, scip->set, scip->stat, scip->transprob) );

   return SCIP_OKAY;
}

/** adds lower bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictLb() should be called for each lower bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictLb() should be called
 *      for each lower bound, whose current assignment lead to the deduction of the given conflict bound.
 */
SCIP_RETCODE SCIPaddConflictLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose lower bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddConflictLb", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->set, scip->stat, var, SCIP_BOUNDTYPE_LOWER, bdchgidx) );

   return SCIP_OKAY;
}

/** adds upper bound of variable at the time of the given bound change index to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictUb() should be called for each upper bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictUb() should be called
 *      for each upper bound, whose current assignment lead to the deduction of the given conflict bound.
 */
SCIP_RETCODE SCIPaddConflictUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddConflictUb", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->set, scip->stat, var, SCIP_BOUNDTYPE_UPPER, bdchgidx) );

   return SCIP_OKAY;
}

/** adds lower or upper bound of variable at the time of the given bound change index to the conflict analysis' candidate
 *  storage; this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBd() should be called for each bound
 *      that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBd() should be called
 *      for each bound, whose current assignment lead to the deduction of the given conflict bound.
 */
SCIP_RETCODE SCIPaddConflictBd(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable whose upper bound should be added to conflict candidate queue */
   SCIP_BOUNDTYPE        boundtype,          /**< the type of the conflicting bound (lower or upper bound) */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index representing time on path to current node, when the
                                              *   conflicting bound was valid, NULL for current local bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddConflictBd", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->set, scip->stat, var, boundtype, bdchgidx) );

   return SCIP_OKAY;
}

/** adds changed bound of fixed binary variable to the conflict analysis' candidate storage;
 *  this method should be called in one of the following two cases:
 *   1. Before calling the SCIPanalyzeConflict() method, SCIPaddConflictBinvar() should be called for each fixed binary
 *      variable that lead to the conflict (e.g. the infeasibility of globally or locally valid constraint).
 *   2. In the propagation conflict resolving method of a constraint handler, SCIPaddConflictBinvar() should be called
 *      for each binary variable, whose current fixing lead to the deduction of the given conflict bound.
 */
SCIP_RETCODE SCIPaddConflictBinvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< binary variable whose changed bound should be added to conflict queue */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddConflictBinvar", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(SCIPvarIsBinary(var));
   if( SCIPvarGetLbLocal(var) > 0.5 )
   {
      SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->set, scip->stat, var, SCIP_BOUNDTYPE_LOWER, NULL) );
   }
   else if( SCIPvarGetUbLocal(var) < 0.5 )
   {
      SCIP_CALL( SCIPconflictAddBound(scip->conflict, scip->set, scip->stat, var, SCIP_BOUNDTYPE_UPPER, NULL) );
   }

   return SCIP_OKAY;
}

/** analyzes conflict bounds that were added after a call to SCIPinitConflictAnalysis() with calls to
 *  SCIPconflictAddLb(), SCIPconflictAddUb(), SCIPconflictAddBd(), or SCIPaddConflictBinvar();
 *  on success, calls the conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  the given valid depth must be a depth level, at which the conflict set defined by calls to SCIPaddConflictLb(),
 *  SCIPaddConflictUb(), SCIPconflictAddBd(), and SCIPaddConflictBinvar() is valid for the whole subtree;
 *  if the conflict was found by a violated constraint, use SCIPanalyzeConflictCons() instead of SCIPanalyzeConflict()
 *  to make sure, that the correct valid depth is used
 */
SCIP_RETCODE SCIPanalyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPanalyzeConflict", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconflictAnalyze(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
         scip->transprob, scip->tree, validdepth, success) );

   return SCIP_OKAY;
}

/** analyzes conflict bounds that were added with calls to SCIPconflictAddLb(), SCIPconflictAddUb(), SCIPconflictAddBd(),
 *  or SCIPaddConflictBinvar(); on success, calls the conflict handlers to create a conflict constraint out of the
 *  resulting conflict set;
 *  the given constraint must be the constraint that detected the conflict, i.e. the constraint that is infeasible
 *  in the local bounds of the initial conflict set (defined by calls to SCIPaddConflictLb(), SCIPaddConflictUb(),
 *  SCIPconflictAddBd(), and SCIPaddConflictBinvar())
 */
SCIP_RETCODE SCIPanalyzeConflictCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that detected the conflict */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPanalyzeConflictCons", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPconsIsGlobal(cons) )
   {
      SCIP_CALL( SCIPconflictAnalyze(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->tree, 0, success) );
   }
   else if( SCIPconsIsActive(cons) )
   {
      SCIP_CALL( SCIPconflictAnalyze(scip->conflict, scip->mem->probmem, scip->set, scip->stat,
            scip->transprob, scip->tree, SCIPconsGetValidDepth(cons), success) );
   }

   return SCIP_OKAY;
}




/*
 * constraint methods
 */

/** creates and captures a constraint of the given constraint handler
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP_CONSDATA*        consdata,           /**< data for this specific constraint */
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
   assert(cons != NULL);
   assert(name != NULL);
   assert(conshdlr != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcreateCons", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPconsCreate(cons, scip->mem->probmem, scip->set, name, conshdlr, consdata,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, TRUE, TRUE) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_FREESOLVE:
      SCIP_CALL( SCIPconsCreate(cons, scip->mem->probmem, scip->set, name, conshdlr, consdata,
            initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, FALSE, TRUE) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** parses constraint information (in cip format) out of a string; if the parsing process was successful a constraint is
 *  creates and captures;
 *
 *  @warning If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution may
 *           be declared feasible even if it violates this particular constraint.  This constellation should only be
 *           used, if no LP or pseudo solution can violate the constraint -- e.g. if a local constraint is redundant due
 *           to the variable's local bounds.
 */
SCIP_RETCODE SCIPparseCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to constraint */
   const char*           str,                /**< name of constraint */
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
   SCIP_Bool             stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   SCIP_Bool*            success             /**< pointer store if the paring process was successful */
   )
{
   assert(cons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPparseCons", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE) );
   
   SCIP_CALL( SCIPconsParse(cons, scip->set, str,
         initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode, success) );

   return SCIP_OKAY;
}

/** increases usage counter of constraint */
SCIP_RETCODE SCIPcaptureCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to capture */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcaptureCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIPconsCapture(cons);

   return SCIP_OKAY;
}

/** decreases usage counter of constraint, if the usage pointer reaches zero the constraint gets freed
 *
 *  @note the pointer of the constraint will be NULLed
 */
SCIP_RETCODE SCIPreleaseCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons                /**< pointer to constraint */
   )
{
   assert(cons != NULL);
   assert(*cons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPreleaseCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPconsRelease(cons, scip->mem->probmem, scip->set) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      if( SCIPconsIsOriginal(*cons) && (*cons)->nuses == 1 )
      {
         SCIPerrorMessage("cannot release last use of original constraint while the transformed problem exists\n");
         return SCIP_INVALIDCALL;
      }
      SCIP_CALL( SCIPconsRelease(cons, scip->mem->probmem, scip->set) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** sets the initial flag of the given constraint */
SCIP_RETCODE SCIPsetConsInitial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             initial             /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsInitial", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetInitial(cons, scip->set, initial) );

   return SCIP_OKAY;
}

/** sets the separate flag of the given constraint */
SCIP_RETCODE SCIPsetConsSeparated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             separate            /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsSeparated", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetSeparated(cons, scip->set, separate) );

   return SCIP_OKAY;
}

/** sets the enforce flag of the given constraint */
SCIP_RETCODE SCIPsetConsEnforced(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             enforce             /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsEnforced", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetEnforced(cons, scip->set, enforce) );

   return SCIP_OKAY;
}

/** sets the check flag of the given constraint */
SCIP_RETCODE SCIPsetConsChecked(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             check               /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsChecked", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetChecked(cons, scip->set, check) );

   return SCIP_OKAY;
}

/** sets the propagate flag of the given constraint */
SCIP_RETCODE SCIPsetConsPropagated(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             propagate           /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsPropagated", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsSetPropagated(cons, scip->set, propagate) );

   return SCIP_OKAY;
}

/** sets the local flag of the given constraint */
SCIP_RETCODE SCIPsetConsLocal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             local               /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsLocal", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPconsSetLocal(cons, local);

   return SCIP_OKAY;
}

/** sets the modifiable flag of the given constraint */
SCIP_RETCODE SCIPsetConsModifiable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             modifiable          /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsModifiable", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE) );

   SCIPconsSetModifiable(cons, modifiable);

   return SCIP_OKAY;
}

/** sets the dynamic flag of the given constraint */
SCIP_RETCODE SCIPsetConsDynamic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             dynamic             /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsDynamic", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPconsSetDynamic(cons, dynamic);

   return SCIP_OKAY;
}

/** sets the removable flag of the given constraint */
SCIP_RETCODE SCIPsetConsRemovable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             removable           /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsRemovable", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPconsSetRemovable(cons, removable);

   return SCIP_OKAY;
}

/** sets the stickingatnode flag of the given constraint */
SCIP_RETCODE SCIPsetConsStickingAtNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Bool             stickingatnode      /**< new value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetConsStickingAtNode", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPconsSetStickingAtNode(cons, stickingatnode);

   return SCIP_OKAY;
}

/** gets and captures transformed constraint of a given constraint; if the constraint is not yet transformed,
 *  a new transformed constraint for this constraint is created
 */
SCIP_RETCODE SCIPtransformCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get/create transformed constraint for */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(transcons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtransformCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPconsIsTransformed(cons) )
   {
      *transcons = cons;
      SCIPconsCapture(*transcons);
   }
   else
   {
      SCIP_CALL( SCIPconsTransform(cons, scip->mem->probmem, scip->set, transcons) );
   }

   return SCIP_OKAY;
}

/** gets and captures transformed constraints for an array of constraints;
 *  if a constraint in the array is not yet transformed, a new transformed constraint for this constraint is created;
 *  it is possible to call this method with conss == transconss
 */
SCIP_RETCODE SCIPtransformConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nconss,             /**< number of constraints to get/create transformed constraints for */
   SCIP_CONS**           conss,              /**< array with constraints to get/create transformed constraints for */
   SCIP_CONS**           transconss          /**< array to store the transformed constraints */
   )
{
   int c;

   assert(nconss == 0 || conss != NULL);
   assert(nconss == 0 || transconss != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtransformConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsTransformed(conss[c]) )
      {
         transconss[c] = conss[c];
         SCIPconsCapture(transconss[c]);
      }
      else
      {
         SCIP_CALL( SCIPconsTransform(conss[c], scip->mem->probmem, scip->set, &transconss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** gets corresponding transformed constraint of a given constraint;
 *  returns NULL as transcons, if transformed constraint is not yet existing
 */
SCIP_RETCODE SCIPgetTransformedCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get the transformed constraint for */
   SCIP_CONS**           transcons           /**< pointer to store the transformed constraint */
   )
{
   assert(transcons != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetTransformedCons", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPconsIsTransformed(cons) )
      *transcons = cons;
   else
      *transcons = SCIPconsGetTransformed(cons);

   return SCIP_OKAY;
}

/** gets corresponding transformed constraints for an array of constraints;
 *  stores NULL in a transconss slot, if the transformed constraint is not yet existing;
 *  it is possible to call this method with conss == transconss, but remember that constraints that are not
 *  yet transformed will be replaced with NULL
 */
SCIP_RETCODE SCIPgetTransformedConss(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nconss,             /**< number of constraints to get the transformed constraints for */
   SCIP_CONS**           conss,              /**< constraints to get the transformed constraints for */
   SCIP_CONS**           transconss          /**< array to store the transformed constraints */
   )
{
   int c;

   assert(nconss == 0 || conss != NULL);
   assert(nconss == 0 || transconss != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetTransformedConss", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   for( c = 0; c < nconss; ++c )
   {
      if( SCIPconsIsTransformed(conss[c]) )
         transconss[c] = conss[c];
      else
         transconss[c] = SCIPconsGetTransformed(conss[c]);
   }

   return SCIP_OKAY;
}

/** adds given value to age of constraint, but age can never become negative;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
SCIP_RETCODE SCIPaddConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real             deltaage            /**< value to add to the constraint's age */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddConsAge", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsAddAge(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob, deltaage) );

   return SCIP_OKAY;
}

/** increases age of constraint by 1.0;
 *  should be called
 *   - in constraint separation, if no cut was found for this constraint,
 *   - in constraint enforcing, if constraint was feasible, and
 *   - in constraint propagation, if no domain reduction was deduced;
 */
SCIP_RETCODE SCIPincConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPincConsAge", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsIncAge(cons, scip->mem->probmem, scip->set, scip->stat, scip->transprob) );

   return SCIP_OKAY;
}

/** resets age of constraint to zero;
 *  should be called
 *   - in constraint separation, if a cut was found for this constraint,
 *   - in constraint enforcing, if the constraint was violated, and
 *   - in constraint propagation, if a domain reduction was deduced;
 */
SCIP_RETCODE SCIPresetConsAge(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPresetConsAge", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsResetAge(cons, scip->set) );

   return SCIP_OKAY;
}

/** enables constraint's separation, propagation, and enforcing capabilities */
SCIP_RETCODE SCIPenableCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPenableCons", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnable(cons, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** disables constraint's separation, propagation, and enforcing capabilities, s.t. the constraint is not propagated,
 *  separated, and enforced anymore until it is enabled again with a call to SCIPenableCons();
 *  in contrast to SCIPdelConsLocal() and SCIPdelConsNode(), the disabling is not associated to a node in the tree and
 *  does not consume memory; therefore, the constraint is neither automatically enabled on leaving the node nor
 *  automatically disabled again on entering the node again;
 *  note that the constraints enforcing capabilities are necessary for the solution's feasibility, if the constraint
 *  is a model constraint; that means, you must be sure that the constraint cannot be violated in the current subtree,
 *  and you have to enable it again manually by calling SCIPenableCons(), if this subtree is left (e.g. by using
 *  an appropriate event handler that watches the corresponding variables' domain changes)
 */
SCIP_RETCODE SCIPdisableCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdisableCons", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsDisable(cons, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** enables constraint's separation capabilities */
SCIP_RETCODE SCIPenableConsSeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPenableConsSeparation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnableSeparation(cons, scip->set) );

   return SCIP_OKAY;
}

/** disables constraint's separation capabilities s.t. the constraint is not propagated anymore until the separation
 *  is enabled again with a call to SCIPenableConsSeparation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 */
SCIP_RETCODE SCIPdisableConsSeparation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdisableConsSeparation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsDisableSeparation(cons, scip->set) );

   return SCIP_OKAY;
}

/** enables constraint's propagation capabilities */
SCIP_RETCODE SCIPenableConsPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPenableConsPropagation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsEnablePropagation(cons, scip->set) );

   return SCIP_OKAY;
}

/** disables constraint's propagation capabilities s.t. the constraint is not propagated anymore until the propagation
 *  is enabled again with a call to SCIPenableConsPropagation(); in contrast to SCIPdelConsLocal() and SCIPdelConsNode(),
 *  the disabling is not associated to a node in the tree and does not consume memory; therefore, the constraint
 *  is neither automatically enabled on leaving the node nor automatically disabled again on entering the node again
 */
SCIP_RETCODE SCIPdisableConsPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdisableConsPropagation", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsDisablePropagation(cons, scip->set) );

   return SCIP_OKAY;
}

/** adds given values to lock status of the constraint and updates the rounding locks of the involved variables */
SCIP_RETCODE SCIPaddConsLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   int                   nlockspos,          /**< increase in number of rounding locks for constraint */
   int                   nlocksneg           /**< increase in number of rounding locks for constraint's negation */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddConsLocks", FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE) );

   SCIP_CALL( SCIPconsAddLocks(cons, scip->set, nlockspos, nlocksneg) );

   return SCIP_OKAY;
}

/** checks single constraint for feasibility of the given solution */
SCIP_RETCODE SCIPcheckCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to check */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcheckCons", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconsCheck(cons, scip->set, sol, checkintegrality, checklprows, printreason, result) );

   return SCIP_OKAY;
}

/** outputs constraint information to file stream */
SCIP_RETCODE SCIPprintCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPprintCons", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPconsPrint(cons, scip->set, file) );

   return SCIP_OKAY;
}

/*
 * LP methods
 */

/** returns, whether the LP was or is to be solved in the current node */
SCIP_Bool SCIPhasCurrentNodeLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPhasCurrentNodeLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeHasCurrentNodeLP(scip->tree);
}

/** returns, whether the LP of the current node is already constructed */
SCIP_Bool SCIPisLPConstructed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisLPConstructed", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeIsFocusNodeLPConstructed(scip->tree);
}

/** makes sure that the LP of the current node is loaded and may be accessed through the LP information methods */
SCIP_RETCODE SCIPconstructLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPconstructLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPconstructCurrentLP(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->tree, scip->lp,
         scip->pricestore, scip->sepastore, scip->branchcand, scip->eventqueue, scip->eventfilter, cutoff) );

   return SCIP_OKAY;
}

/** makes sure that the LP of the current node is flushed */
SCIP_RETCODE SCIPflushLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPflushLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpFlush(scip->lp, scip->mem->probmem, scip->set, scip->eventqueue) );

   return SCIP_OKAY;
}

/** gets solution status of current LP */
SCIP_LPSOLSTAT SCIPgetLPSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPSolstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetSolstat(scip->lp);
   else
      return SCIP_LPSOLSTAT_NOTSOLVED;
}

/** returns whether the current lp is a relaxation of the current problem and its optimal objective value is a local lower bound */
SCIP_Bool SCIPisLPRelax(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisLPRelax", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpIsRelax(scip->lp);
}

/** gets objective value of current LP (which is the sum of column and loose objective value) */
SCIP_Real SCIPgetLPObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetObjval(scip->lp, scip->set);
}

/** gets part of objective value of current LP that results from COLUMN variables only */
SCIP_Real SCIPgetLPColumnObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPColumnObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetColumnObjval(scip->lp);
}

/** gets part of objective value of current LP that results from LOOSE variables only */
SCIP_Real SCIPgetLPLooseObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPLooseObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetLooseObjval(scip->lp, scip->set);
}

/** gets pseudo objective value of the current LP */
SCIP_Real SCIPgetPseudoObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPseudoObjval", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetPseudoObjval(scip->lp, scip->set);
}

/** returns whether the root lp is a relaxation of the problem and its optimal objective value is a global lower bound */
SCIP_Bool SCIPisRootLPRelax(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisRootLPRelax", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpIsRootLPRelax(scip->lp);
}

/** gets the objective value of the root node LP; returns SCIP_INVALID if the root node LP was not (yet) solved */
SCIP_Real SCIPgetLPRootObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPRootObjval", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetRootObjval(scip->lp);
}

/** gets part of the objective value of the root node LP that results from COLUMN variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
SCIP_Real SCIPgetLPRootColumnObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPRootColumnObjval", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetRootColumnObjval(scip->lp);
}

/** gets part of the objective value of the root node LP that results from LOOSE variables only;
 *  returns SCIP_INVALID if the root node LP was not (yet) solved
 */
SCIP_Real SCIPgetLPRootLooseObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPRootLooseObjval", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpGetRootLooseObjval(scip->lp);
}

/** gets current LP columns along with the current number of LP columns */
SCIP_RETCODE SCIPgetLPColsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL***           cols,               /**< pointer to store the array of LP columns, or NULL */
   int*                  ncols               /**< pointer to store the number of LP columns, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPColsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      if( cols != NULL )
         *cols = SCIPlpGetCols(scip->lp);
      if( ncols != NULL )
         *ncols = SCIPlpGetNCols(scip->lp);
   }
   else
   {
      if( cols != NULL )
         *cols = NULL;
      if( ncols != NULL )
         *ncols = 0;
   }

   return SCIP_OKAY;
}

/** gets current LP columns */
SCIP_COL** SCIPgetLPCols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetCols(scip->lp);
   else
      return NULL;
}

/** gets current number of LP columns */
int SCIPgetNLPCols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPCols", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetNCols(scip->lp);
   else
      return 0;
}

/** gets current LP rows along with the current number of LP rows */
SCIP_RETCODE SCIPgetLPRowsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW***           rows,               /**< pointer to store the array of LP rows, or NULL */
   int*                  nrows               /**< pointer to store the number of LP rows, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPRowsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      if( rows != NULL )
         *rows = SCIPlpGetRows(scip->lp);
      if( nrows != NULL )
         *nrows = SCIPlpGetNRows(scip->lp);
   }
   else
   {
      if( rows != NULL )
         *rows = NULL;
      if( nrows != NULL )
         *nrows = 0;
   }

   return SCIP_OKAY;
}

/** gets current LP rows */
SCIP_ROW** SCIPgetLPRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetRows(scip->lp);
   else
      return NULL;
}

/** gets current number of LP rows */
int SCIPgetNLPRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
      return SCIPlpGetNRows(scip->lp);
   else
      return 0;
}

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 */
SCIP_Bool SCIPallColsInLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPallColsInLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp);
}

/** returns whether the current LP solution is basic, i.e. is defined by a valid simplex basis */
SCIP_Bool SCIPisLPSolBasic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisLPSolBasic", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPlpIsSolBasic(scip->lp);
}


/** gets all indices of basic columns and rows: index i >= 0 corresponds to column i, index i < 0 to row -i-1 */
SCIP_RETCODE SCIPgetLPBasisInd(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  basisind            /**< pointer to store the basis indices */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPBasisInd", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBasisInd(scip->lp, basisind) );
   
   return SCIP_OKAY;
}

/** gets a row from the inverse basis matrix B^-1 */
SCIP_RETCODE SCIPgetLPBInvRow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPBInvRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvRow(scip->lp, r, coef) );

   /* debug check if the coef is the r-th line of the inverse matrix B^-1 */
   SCIP_CALL( SCIPdebugCheckBInvRow(scip, r, coef) );

   return SCIP_OKAY;
}

/** gets a column from the inverse basis matrix B^-1 */
SCIP_RETCODE SCIPgetLPBInvCol(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the LP
                                              *   returned by SCIPcolGetLPPos(); you have to call SCIPgetBasisInd()
                                              *   to get the array which links the B^-1 column numbers to the row and
                                              *   column numbers of the LP! c must be between 0 and nrows-1, since the
                                              *   basis has the size nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPBInvCol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvCol(scip->lp, c, coef) );

   return SCIP_OKAY;
}

/** gets a row from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A) */
SCIP_RETCODE SCIPgetLPBInvARow(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   r,                  /**< row number */
   SCIP_Real*            binvrow,            /**< row in B^-1 from prior call to SCIPgetLPBInvRow(), or NULL */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPBInvARow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvARow(scip->lp, r, binvrow, coef) );

   return SCIP_OKAY;
}

/** gets a column from the product of inverse basis matrix B^-1 and coefficient matrix A (i.e. from B^-1 * A),
 *  i.e., it computes B^-1 * A_c with A_c being the c'th column of A
 */
SCIP_RETCODE SCIPgetLPBInvACol(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   c,                  /**< column number which can be accessed by SCIPcolGetLPPos() */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPBInvACol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolBasic(scip->lp) )
   {
      SCIPerrorMessage("current LP solution is not basic\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPlpGetBInvACol(scip->lp, c, coef) );

   return SCIP_OKAY;
}

/** calculates a weighted sum of all LP rows; for negative weights, the left and right hand side of the corresponding
 *  LP row are swapped in the summation
 */
SCIP_RETCODE SCIPsumLPRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   SCIP_REALARRAY*       sumcoef,            /**< array to store sum coefficients indexed by variables' probindex */
   SCIP_Real*            sumlhs,             /**< pointer to store the left hand side of the row summation */
   SCIP_Real*            sumrhs              /**< pointer to store the right hand side of the row summation */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsumLPRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpSumRows(scip->lp, scip->set, scip->transprob, weights, sumcoef, sumlhs, sumrhs) );

   return SCIP_OKAY;
}

/* calculates a MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
SCIP_RETCODE SCIPcalcMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcalcMIR", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpCalcMIR(scip->lp, scip->set, scip->stat, scip->transprob, sol,
         boundswitch, usevbds, allowlocal, fixintegralrhs, boundsfortrans, boundtypesfortrans, maxmksetcoefs, 
         maxweightrange, minfrac, maxfrac, weights, scale, mksetcoefs, mksetcoefsvalid, mircoef, mirrhs, cutactivity, 
         success, cutislocal) );

   return SCIP_OKAY;
}

/* calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because these
 * rows cannot participate in a MIR cut.
 */
SCIP_RETCODE SCIPcalcStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mircoef,            /**< array to store strong CG coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the strong CG row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   SCIP_Bool*            cutislocal          /**< pointer to store whether the returned cut is only valid locally */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcalcStrongCG", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPlpCalcStrongCG(scip->lp, scip->set, scip->stat, scip->transprob,
         boundswitch, usevbds, allowlocal, maxmksetcoefs, maxweightrange, minfrac, maxfrac, weights, scale,
         mircoef, mirrhs, cutactivity, success, cutislocal) );

   return SCIP_OKAY;
}

/** writes current LP to a file */
SCIP_RETCODE SCIPwriteLP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname               /**< file name */
   )
{

   SCIP_Bool cutoff;
   
   SCIP_CALL( checkStage(scip, "SCIPwriteLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   if( !SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      SCIP_CALL( SCIPconstructCurrentLP(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->tree, scip->lp,
            scip->pricestore, scip->sepastore, scip->branchcand, scip->eventqueue, scip->eventfilter, &cutoff) );
   }

   /* we need a flushed lp to write the current lp */
   SCIP_CALL( SCIPlpFlush(scip->lp, scip->mem->probmem, scip->set, scip->eventqueue) );

   SCIP_CALL( SCIPlpWrite(scip->lp, fname) );
   
   return SCIP_OKAY;
}

/** writes MIP relaxation of the current branch-and-bound node to a file */
SCIP_RETCODE SCIPwriteMIP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname,              /**< file name */
   SCIP_Bool             genericnames,       /**< should generic names like x_i and row_j be used in order to avoid
                                              *   troubles with reserved symbols? */
   SCIP_Bool             origobj             /**< should the original objective function be used? */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPwriteMIP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* we need a flushed lp to write the current mip */
   SCIP_CALL( SCIPlpFlush(scip->lp, scip->mem->probmem, scip->set, scip->eventqueue) );
   
   SCIP_CALL( SCIPlpWriteMip(scip->lp, scip->set, fname, genericnames,
         origobj, scip->origprob->objsense, scip->transprob->objscale, scip->transprob->objoffset) );

   return SCIP_OKAY;
}

/** gets the LP interface of SCIP;
 *  with the LPI you can use all of the methods defined in scip/lpi.h;
 *
 *  @warning You have to make sure, that the full internal state of the LPI does not change or is recovered completely
 *           after the end of the method that uses the LPI. In particular, if you manipulate the LP or its solution
 *           (e.g. by calling one of the SCIPlpiAdd...() or one of the SCIPlpiSolve...() methods), you have to check in
 *           advance with SCIPlpiWasSolved() whether the LP is currently solved. If this is the case, you have to make
 *           sure, the internal solution status is recovered completely at the end of your method. This can be achieved
 *           by getting the LPI state before applying any LPI manipulations with SCIPlpiGetState() and restoring it
 *           afterwards with SCIPlpiSetState() and SCIPlpiFreeState(). Additionally you have to resolve the LP with the
 *           appropriate SCIPlpiSolve...() call in order to reinstall the internal solution status.
 *
 *  @warning Make also sure, that all parameter values that you have changed are set back to their original values.
 */
SCIP_RETCODE SCIPgetLPI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_LPI**            lpi                 /**< pointer to store the LP interface */
   )
{
   assert(lpi != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetLPI", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   *lpi = SCIPlpGetLPI(scip->lp);

   return SCIP_OKAY;
}

/** displays quality information about the current LP solution
 * an LP solution need to be available
 * information printed is subject to what the LP solver supports
 */
SCIP_RETCODE SCIPprintLPSolutionQuality(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_LPI* lpi;
   SCIP_Real quality;

   SCIP_CALL( checkStage(scip, "SCIPprintLPSolutionQuality", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
      case SCIP_STAGE_INIT:
      case SCIP_STAGE_PROBLEM:
      case SCIP_STAGE_TRANSFORMED:
      case SCIP_STAGE_PRESOLVING:
      case SCIP_STAGE_PRESOLVED:
         SCIPmessageFPrintInfo(file, "Problem not solving yet, no LP available.\n");
         return SCIP_OKAY;

      case SCIP_STAGE_SOLVING:
      case SCIP_STAGE_SOLVED:
         break;

      default:
         SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
         return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   lpi = SCIPlpGetLPI(scip->lp);
   assert(lpi != NULL);

   SCIP_CALL( SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_ESTIMCONDITION, &quality) );
   SCIPmessageFPrintInfo(file, "Basis matrix condition (estimated): ");
   if( quality != SCIP_INVALID )
      SCIPmessageFPrintInfo(file, "%.6e\n", quality);
   else
      SCIPmessageFPrintInfo(file, "not available\n", quality);

   SCIP_CALL( SCIPlpiGetRealSolQuality(lpi, SCIP_LPSOLQUALITY_EXACTCONDITION, &quality) );
   SCIPmessageFPrintInfo(file, "Basis matrix condition (exact):     ");
   if( quality != SCIP_INVALID )
      SCIPmessageFPrintInfo(file, "%.6e\n", quality);
   else
      SCIPmessageFPrintInfo(file, "not available\n", quality);

   return SCIP_OKAY;
}

/** Compute relative interior point to current LP w.r.t. one-norm
 *
 *  We use the approach of@par
 *  R. Freund, R. Roundy, M. J. Todd@par
 *  "Identifying the Set of Always-Active Constraints in a System of Linear Inequalities by a Single Linear Program"@par
 *  Tech. Rep, No. 1674-85, Sloan School, M.I.T., 1985 
 *
 *  to compute a relative interior point for the current LP.
 *
 *  Assume the orginal LP looks as follows:
 *  \f[
 *     \begin{array}{rrl}
 *        \min & c^T x &\\
 *             & A x & \geq a\\
 *             & B x & \leq b\\
 *             & D x & = d.
 *     \end{array}
 *  \f]
 *  Note that bounds should be included in the system. The following artificial LP does the job:
 *  \f[
 *     \begin{array}{rrl}
 *        \max & 1^T y &\\
 *             & A x - y - \alpha a & \geq 0\\
 *             & B x + y - \alpha b & \leq 0\\
 *             & D x - \alpha d & = 0\\
 *             & 0 \leq y & \leq 1\\
 *             & alpha & \geq 1.
 *     \end{array}
 *  \f]
 *  If the original LP is feasible, this LP is feasible as well. Any optimal solution yields the
 *  relative interior point \f$x^*_j/\alpha^*\f$.
 */
static
SCIP_RETCODE computeLPRelIntPointOneNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             relaxrows,          /**< should the rows be relaxed */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   SCIP_SOL**            point               /**< relative interior point on exit */
   )
{
   SCIP_LP* lp;
   SCIP_LPI* lpi;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* primal;
   SCIP_Real* colvals;
   SCIP_Real zero;
   SCIP_Real minusinf;
   SCIP_Real plusinf;
   SCIP_Real objval;
   SCIP_Real alpha;
   int* rowinds;
   int* colinds;
   int nnewcols;
#ifndef NDEBUG
   int nslacks;
#endif
   int beg;
   int cnt;
   int i;
   int j;

   assert( scip != NULL );
   assert( point != NULL );

   SCIP_CALL( checkStage(scip, "computeLPRelIntPoint", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   /* initialize point */
   *point = NULL;

   /* exit if there are no columns */
   lp = scip->lp;
   assert( lp->nrows >= 0 && lp->ncols >= 0 );
   if( lp->ncols == 0 )
      return SCIP_OKAY;

   /* disable objective cutoff if we have none */
   if( inclobjcutoff && (SCIPgetCutoffbound(scip) >= SCIPinfinity(scip) || SCIPlpGetLooseObjval(lp, scip->set) <= -SCIPinfinity(scip) || SCIPlpGetLooseObjval(lp, scip->set) == SCIP_INVALID) )
      inclobjcutoff = FALSE;

   SCIPdebugMessage("Computing relative interior point to current LP.\n");

   /* if there are no rows, we return the zero point */
   if( lp->nrows == 0 && !inclobjcutoff )
   {
      /* create zero point */
      SCIP_CALL( SCIPcreateSol(scip, point, NULL) );
      return SCIP_OKAY;
   }

   /* create auxiliary LP */
   SCIP_CALL( SCIPlpiCreate(&lpi, "relativeInterior", SCIP_OBJSEN_MAXIMIZE) );

   /* get storage */
   nnewcols = 3*lp->ncols + 2*lp->nrows + (inclobjcutoff ? 1 : 0) + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nnewcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nnewcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nnewcols) );

   /* create original columns (bounds are relaxed below, unless the variable is fixed) */
   for( j = 0; j < lp->ncols; ++j )
   {
      /* note: if the variable is fixed we cannot simply fix the variables (because alpha scales the problem) */
      obj[j] = 0.0;
      lb[j] = -SCIPlpiInfinity(lpi);
      ub[j] = SCIPlpiInfinity(lpi);
   }

   /* add artificial alpha variable */
   nnewcols = lp->ncols;
   obj[nnewcols] = 0.0;
   lb[nnewcols] = 1.0;
   ub[nnewcols] = SCIPlpiInfinity(lpi);
   ++nnewcols;
   
   /* create slacks for rows */
   for( i = 0; i < lp->nrows; ++i )
   {
      SCIP_ROW* row;

      row = lp->rows[i];
      assert( row != NULL );

      if( SCIProwIsModifiable(row) )
         continue;

      /* check whether we have an equation */
      if( SCIPisEQ(scip, row->lhs, row->rhs) )
      {
         assert( !SCIPisInfinity(scip, ABS(row->lhs)) );
         assert( !SCIPisInfinity(scip, ABS(row->rhs)) );
      }
      else
      {
         if( relaxrows )
         {
            /* otherwise add slacks for each side if necessary */
            if( !SCIPisInfinity(scip, ABS(row->lhs)) )
            {
               obj[nnewcols] = 1.0;
               lb[nnewcols] = 0.0;
               ub[nnewcols] = 1.0;
               ++nnewcols;
            }
            if( !SCIPisInfinity(scip, ABS(row->rhs)) )
            {
               obj[nnewcols] = 1.0;
               lb[nnewcols] = 0.0;
               ub[nnewcols] = 1.0;
               ++nnewcols;
            }
         }
      }
   }

   /* create slacks for objective cutoff row */
   if( inclobjcutoff )
   {
      if( relaxrows )
      {
         /* add slacks for right hand side */
         obj[nnewcols] = 1.0;
         lb[nnewcols] = 0.0;
         ub[nnewcols] = 1.0;
         ++nnewcols;
      }
   }

   /* create slacks for bounds */
   for( j = 0; j < lp->ncols; ++j )
   {
      SCIP_COL* col;

      col = lp->cols[j];
      assert( col != NULL );

      /* add slacks for each bound if necessary */
      if( !SCIPisInfinity(scip, ABS(col->lb)) )
      {
         obj[nnewcols] = 1.0;
         lb[nnewcols] = 0.0;
         ub[nnewcols] = 1.0;
         ++nnewcols;
      }
      if( !SCIPisInfinity(scip, ABS(col->ub)) )
      {
         obj[nnewcols] = 1.0;
         lb[nnewcols] = 0.0;
         ub[nnewcols] = 1.0;
         ++nnewcols;
      }
   }
#ifndef NDEBUG
   nslacks = nnewcols - lp->ncols - 1;
   assert( nslacks >= 0 );
   assert( nnewcols <= 3*lp->ncols + 2*lp->nrows + (inclobjcutoff ? 1 : 0) + 1 );
#endif

   /* add columns */
   SCIP_CALL( SCIPlpiAddCols(lpi, nnewcols, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

   /* free storage */
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   /* prepare storage for rows */
   SCIP_CALL( SCIPallocBufferArray(scip, &rowinds, lp->nrows + (inclobjcutoff ? 1 : 0) + 2*lp->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colinds, lp->ncols+2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colvals, lp->ncols+2) );

   /* create rows arising from original rows */
   cnt = 0;
   beg = 0;
   zero = 0.0;
   minusinf = -SCIPlpiInfinity(lpi);
   plusinf = SCIPlpiInfinity(lpi);
   for( i = 0; i < lp->nrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nnonz;

      row = lp->rows[i];
      assert( row != NULL );

      if( SCIProwIsModifiable(row) )
         continue;

      /* get row data */
      lhs = row->lhs - (SCIPisInfinity(scip, -row->lhs) ? 0.0 : row->constant);
      rhs = row->rhs - (SCIPisInfinity(scip,  row->rhs) ? 0.0 : row->constant);
      nnonz = row->nlpcols;
      assert( nnonz <= lp->ncols );
      rowcols = row->cols;
      rowvals = row->vals;

      /* set up indices */
      for( j = 0; j < nnonz; ++j )
      {
         assert( rowcols[j] != NULL );
         assert( 0 <= rowcols[j]->lppos && rowcols[j]->lppos < lp->ncols );
         assert( lp->cols[rowcols[j]->lppos] == rowcols[j] );
         colinds[j] = rowcols[j]->lppos;
         colvals[j] = rowvals[j];
      }

      /* if we have an equation */
      if( SCIPisEQ(scip, lhs, rhs) )
      {
         /* add artificial variable */
         colinds[nnonz] = lp->ncols;
         colvals[nnonz] = -rhs;

         /* add row */
         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &zero, &zero, NULL, nnonz+1, &beg, colinds, colvals) );
      }
      else
      {
         /* treat lhs */
         if( !SCIPisInfinity(scip, ABS(lhs)) )
         {
            assert( ! SCIPisEQ(scip, lhs, rhs) );

            /* add artificial variable */
            colinds[nnonz] = lp->ncols;
            colvals[nnonz] = -lhs;

            if( relaxrows )
            {
               /* add slack variable */
               colinds[nnonz+1] = lp->ncols + 1 + cnt;
               colvals[nnonz+1] = -1.0;
               ++cnt;
               SCIP_CALL( SCIPlpiAddRows(lpi, 1, &zero, &plusinf, NULL, nnonz+2, &beg, colinds, colvals) );
            }
            else
            {
               SCIP_CALL( SCIPlpiAddRows(lpi, 1, &zero, &plusinf, NULL, nnonz+1, &beg, colinds, colvals) );
            }
         }

         /* treat rhs */
         if( !SCIPisInfinity(scip, ABS(rhs)) )
         {
            assert( ! SCIPisEQ(scip, lhs, rhs) );

            /* add artificial variable */
            colinds[nnonz] = lp->ncols;
            colvals[nnonz] = -rhs;

            if( relaxrows )
            {
               /* add slack variable */
               colinds[nnonz+1] = lp->ncols + 1 + cnt;
               colvals[nnonz+1] = 1.0;
               ++cnt;
               SCIP_CALL( SCIPlpiAddRows(lpi, 1, &minusinf, &zero, NULL, nnonz+2, &beg, colinds, colvals) );
            }
            else
            {
               SCIP_CALL( SCIPlpiAddRows(lpi, 1, &minusinf, &zero, NULL, nnonz+1, &beg, colinds, colvals) );
            }
         }
      }
   }

   /* create row arising from objective cutoff */
   if( inclobjcutoff )
   {
      SCIP_Real rhs;
      int nnonz;

      /* get row data */
      rhs = SCIPgetCutoffbound(scip) - SCIPlpGetLooseObjval(lp, scip->set);

      /* set up indices and coefficients */
      nnonz = 0;
      for( j = 0; j < lp->ncols; ++j )
      {
         assert( lp->cols[j] != NULL );
         assert( 0 <= lp->cols[j]->lppos && lp->cols[j]->lppos < lp->ncols );
         assert( lp->cols[lp->cols[j]->lppos] == lp->cols[j] );
         if( lp->cols[j]->obj == 0.0 )
            continue;
         colinds[nnonz] = lp->cols[j]->lppos;
         colvals[nnonz] = lp->cols[j]->obj;
         ++nnonz;
      }

      /* treat rhs */
      /* add artificial variable */
      colinds[nnonz] = lp->ncols;
      colvals[nnonz] = -rhs;
      ++nnonz;

      if( relaxrows )
      {
         /* add slack variable */
         colinds[nnonz] = lp->ncols + 1 + cnt;
         colvals[nnonz] = 1.0;
         ++nnonz;
         ++cnt;
      }
      SCIP_CALL( SCIPlpiAddRows(lpi, 1, &minusinf, &zero, NULL, nnonz, &beg, colinds, colvals) );
   }

   /* create rows arising from bounds */
   for( j = 0; j < lp->ncols; ++j )
   {
      SCIP_COL* col;

      col = lp->cols[j];
      assert( col != NULL );
      assert( col->lppos == j );

      /* set up index of column */
      colinds[0] = j;
      colvals[0] = 1.0;
      
      /* lower bound */
      if( !SCIPisInfinity(scip, ABS(col->lb)) )
      {
         /* add artificial variable */
         colinds[1] = lp->ncols;
         colvals[1] = -col->lb;
         
         /* add slack variable */
         colinds[2] = lp->ncols + 1 + cnt;
         colvals[2] = -1.0;
         ++cnt;
         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &zero, &plusinf, NULL, 3, &beg, colinds, colvals) );
      }

      /* upper bound */
      if( !SCIPisInfinity(scip, ABS(col->ub)) )
      {
         /* add artificial variable */
         colinds[1] = lp->ncols;
         colvals[1] = -col->ub;
         
         /* add slack variable */
         colinds[2] = lp->ncols + 1 + cnt;
         colvals[2] = 1.0;
         ++cnt;
         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &minusinf, &zero, NULL, 3, &beg, colinds, colvals) );
      }
   }
   assert( cnt == nslacks );

   SCIPfreeBufferArray(scip, &colvals);
   SCIPfreeBufferArray(scip, &colinds);
   SCIPfreeBufferArray(scip, &rowinds);

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPlpiWriteLP(lpi, "relativeInteriorOneNorm.lp") );
#endif

   /* solve and store point */
   /* SCIP_CALL( SCIPlpiSolvePrimal(lpi) ); */
   SCIP_CALL( SCIPlpiSolveDual(lpi) );  /* dual is usually faster */

   if( SCIPlpiIsOptimal(lpi) )
   {   
      /* get primal solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &primal, nnewcols) );
      SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primal, NULL, NULL, NULL) );
      alpha = primal[lp->ncols];
      assert( SCIPisFeasGE(scip, alpha, 1.0) );

      SCIPdebugMessage("Solved relative interior lp with objective %g.\n", objval);

      /* construct relative interior point */
      SCIP_CALL( SCIPcreateSol(scip, point, NULL) );
      for( j = 0; j < lp->ncols; ++j )
      {
         if( !SCIPisFeasZero(scip, primal[j]) )
         {
            SCIP_VAR* var;
            var = lp->cols[j]->var;
            SCIP_CALL( SCIPsetSolVal(scip, *point, var, primal[j]/alpha) );
         }
      }

#ifdef SCIP_DEBUG
      /* check whether the point is a relative interior point */
      cnt = 0;
      if( relaxrows )
      {
         for( i = 0; i < lp->nrows; ++i )
         {
            SCIP_ROW* row;
            SCIP_COL** rowcols;
            SCIP_Real* rowvals;
            SCIP_Real lhs;
            SCIP_Real rhs;
            SCIP_Real sum;
            int nnonz;

            row = lp->rows[i];
            assert( row != NULL );

            /* get row data */
            lhs = row->lhs - (SCIPisInfinity(scip, -row->lhs) ? 0.0 : row->constant);
            rhs = row->rhs - (SCIPisInfinity(scip,  row->rhs) ? 0.0 : row->constant);
            nnonz = row->nlpcols;
            assert( nnonz <= lp->ncols );
            rowcols = row->cols;
            rowvals = row->vals;

            sum = 0.0;
            for( j = 0; j < nnonz; ++j )
               sum += rowvals[j] * primal[rowcols[j]->lppos];
            sum /= alpha;

            /* if we have an equation */
            if( SCIPisEQ(scip, lhs, rhs) )
            {
               assert( SCIPisFeasEQ(scip, sum, lhs) );
            }
            else
            {
               /* treat lhs */
               if( !SCIPisInfinity(scip, ABS(lhs)) )
               {
                  assert( SCIPisFeasZero(scip, primal[lp->ncols+1+cnt]) || SCIPisFeasGT(scip, sum, lhs) );
                  ++cnt;
               }
               /* treat rhs */
               if( !SCIPisInfinity(scip, ABS(rhs)) )
               {
                  assert( SCIPisFeasZero(scip, primal[lp->ncols+1+cnt]) || SCIPisFeasLT(scip, sum, rhs) );
                  ++cnt;
               }
            }
         }
         if( inclobjcutoff )
         {
            SCIP_Real sum;
            SCIP_Real rhs;

            sum = 0.0;
            for( j = 0; j < lp->ncols; ++j )
               sum += lp->cols[j]->obj * primal[lp->cols[j]->lppos];
            sum /= alpha;

            rhs = SCIPgetCutoffbound(scip) - SCIPlpGetLooseObjval(lp, scip->set);
            assert( SCIPisFeasZero(scip, primal[lp->ncols+1+cnt]) || SCIPisFeasLT(scip, sum, rhs) );
            ++cnt;
         }
      }
      /* check bounds */
      for( j = 0; j < lp->ncols; ++j )
      {
         SCIP_COL* col;
         SCIP_Real val;

         col = lp->cols[j];
         assert( col != NULL );
         val = primal[col->lppos] / alpha;

         /* if the variable is not fixed */
         if( !SCIPisEQ(scip, col->lb, col->ub) )
         {
            /* treat lb */
            if( !SCIPisInfinity(scip, ABS(col->lb)) )
            {
               assert( SCIPisFeasZero(scip, primal[lp->ncols+1+cnt]) || SCIPisFeasGT(scip, val, col->lb) );
               ++cnt;
            }
            /* treat rhs */
            if( !SCIPisInfinity(scip, ABS(col->ub)) )
            {
               assert( SCIPisFeasZero(scip, primal[lp->ncols+1+cnt]) || SCIPisFeasLT(scip, val, col->ub) );
               ++cnt;
            }
         }
      }
#endif

      /* free */
      SCIPfreeBufferArray(scip, &primal);
   }
   SCIP_CALL( SCIPlpiFree(&lpi) );

   return SCIP_OKAY;
}

/** Compute relative interior point to current LP w.r.t. supremum norm
 *
 *  Assume the orginal LP looks as follows:
 *  \f[
 *     \begin{array}{rrl}
 *        \min & c^T x &\\
 *             & A x & \geq a\\
 *             & B x & \leq b\\
 *             & D x & = d.
 *     \end{array}
 *  \f]
 *  Note that bounds should be included in the system. The following artificial LP does the job:
 *  \f[
 *     \begin{array}{rrl}
 *        \max & sigma &\\
 *             & <A_i, x> - sigma ||A_i|| & \geq a_i\\
 *             & <B_i, x> + sigma ||B_i|| & \leq b_i\\
 *             & D x & = d\\
 *             & sigma & \geq 0.
 *     \end{array}
 *  \f]
 *  If the original LP is feasible, this LP is feasible as well. Any optimal solution yields a
 *  relative interior point.
 */
static
SCIP_RETCODE computeLPRelIntPointSupNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   SCIP_SOL**            point               /**< relative interior point on exit */
   )
{
   SCIP_LP* lp;
   SCIP_LPI* lpi;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* primal;
   SCIP_Real* colvals;
   SCIP_Real minusinf;
   SCIP_Real plusinf;
   SCIP_Real objval;
   int* colinds;
   int ncols;
   int beg;
   int i;
   int j;
#ifndef NDEBUG
   SCIP_Real sigma;
#endif

   assert( scip != NULL );
   assert( point != NULL );

   SCIP_CALL( checkStage(scip, "computeLPRelIntPointSup", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   /* initialize point */
   *point = NULL;

   /* exit if there are no columns */
   lp = scip->lp;
   assert( lp->nrows >= 0 && lp->ncols >= 0 );
   if( lp->ncols == 0 )
      return SCIP_OKAY;

   /* disable objective cutoff if we have none */
   if( inclobjcutoff && (SCIPgetCutoffbound(scip) >= SCIPinfinity(scip) || SCIPlpGetLooseObjval(lp, scip->set) <= -SCIPinfinity(scip) || SCIPlpGetLooseObjval(lp, scip->set) == SCIP_INVALID) )
      inclobjcutoff = FALSE;

   SCIPdebugMessage("Computing relative interior point to current LP.\n");

   /* if there are no rows, we return the zero point */
   if( lp->nrows == 0 && !inclobjcutoff )
   {
      /* create zero point */
      SCIP_CALL( SCIPcreateSol(scip, point, NULL) );
      return SCIP_OKAY;
   }

   /* create auxiliary LP */
   SCIP_CALL( SCIPlpiCreate(&lpi, "relativeInteriorSup", SCIP_OBJSEN_MAXIMIZE) );

   /* get storage */
   ncols = lp->ncols + 1;
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, ncols) );

   /* create original columns (bounds are relaxed below, unless the variable is fixed) */
   for( j = 0; j < lp->ncols; ++j )
   {
      obj[j] = 0.0;
      lb[j]  = lp->cols[j]->lb;
      ub[j]  = lp->cols[j]->ub;
   }

   /* add artificial sigma variable */
   obj[lp->ncols] = 1.0;
   lb[lp->ncols]  = 0.0;
   ub[lp->ncols]  = SCIPlpiInfinity(lpi);

   /* add columns */
   SCIP_CALL( SCIPlpiAddCols(lpi, ncols, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

   /* free storage */
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   /* prepare storage for rows */
   SCIP_CALL( SCIPallocBufferArray(scip, &colinds, lp->ncols+1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colvals, lp->ncols+1) );

   /* create rows arising from original rows */
   beg = 0;
   minusinf = -SCIPlpiInfinity(lpi);
   plusinf  =  SCIPlpiInfinity(lpi);
   for( i = 0; i < lp->nrows; ++i )
   {
      SCIP_ROW* row;
      SCIP_COL** rowcols;
      SCIP_Real* rowvals;
      SCIP_Real lhs;
      SCIP_Real rhs;
      int nnonz;

      row = lp->rows[i];
      assert( row != NULL );

      if( SCIProwIsModifiable(row) )
         continue;

      /* get row data */
      lhs = row->lhs - (SCIPisInfinity(scip, -row->lhs) ? 0.0 : row->constant);
      rhs = row->rhs - (SCIPisInfinity(scip,  row->rhs) ? 0.0 : row->constant);
      nnonz = row->nlpcols;
      assert( nnonz <= lp->ncols );
      rowcols = row->cols;
      rowvals = row->vals;

      /* set up indices */
      for( j = 0; j < nnonz; ++j )
      {
         assert( rowcols[j] != NULL );
         assert( 0 <= rowcols[j]->lppos && rowcols[j]->lppos < lp->ncols );
         assert( lp->cols[rowcols[j]->lppos] == rowcols[j] );
         colinds[j] = rowcols[j]->lppos;
         colvals[j] = rowvals[j];
      }

      if( SCIPisEQ(scip, lhs, rhs) )
      {
         /* add original row to LP -> this is for the 'relative' in the relative interior */
         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &rhs, NULL, nnonz, &beg, colinds, colvals) );
      }
      else
      {
         /* add row <A_i,x> + sigma ||A_i|| <= rhs, if rhs is finite */
         if( !SCIPisInfinity(scip,  rhs) )
         {
            colinds[nnonz] = lp->ncols;
            colvals[nnonz] =  SCIProwGetNorm(row);

            SCIP_CALL( SCIPlpiAddRows(lpi, 1, &minusinf, &rhs, NULL, nnonz+1, &beg, colinds, colvals) );
         }

         /* add row <A_i,x> - sigma ||A_i|| >= lhs, if lhs is finite */
         if( !SCIPisInfinity(scip, -lhs) )
         {
            colinds[nnonz] = lp->ncols;
            colvals[nnonz] = -SCIProwGetNorm(row);

            SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs, &plusinf, NULL, nnonz+1, &beg, colinds, colvals) );
         }
      }
   }

   /* create row arising from objective cutoff */
   if( inclobjcutoff )
   {
      SCIP_Real rhs;
      SCIP_Real norm;
      int nnonz;

      /* get row data */
      rhs = SCIPgetCutoffbound(scip) - SCIPlpGetLooseObjval(lp, scip->set);

      /* set up indices and coefficients */
      nnonz = 0;
      norm = 0.0;
      for( j = 0; j < lp->ncols; ++j )
      {
         assert( lp->cols[j] != NULL );
         assert( 0 <= lp->cols[j]->lppos && lp->cols[j]->lppos < lp->ncols );
         assert( lp->cols[lp->cols[j]->lppos] == lp->cols[j] );
         if( lp->cols[j]->obj == 0.0 )
            continue;
         colinds[nnonz] = lp->cols[j]->lppos;
         colvals[nnonz] = lp->cols[j]->obj;
         norm += lp->cols[j]->obj * lp->cols[j]->obj;
         ++nnonz;
      }

      /* add artificial variable */
      colinds[nnonz] = lp->ncols;
      colvals[nnonz] = sqrt(norm);
      ++nnonz;

      SCIP_CALL( SCIPlpiAddRows(lpi, 1, &minusinf, &rhs, NULL, nnonz, &beg, colinds, colvals) );
   }

   /* create rows arising from bounds */
   for( j = 0; j < lp->ncols; ++j )
   {
      SCIP_COL* col;

      col = lp->cols[j];
      assert( col != NULL );
      assert( col->lppos == j );

      /* skip fixed variables, column bounds have been set above already */
      if( SCIPisEQ(scip, col->lb, col->ub) )
         continue;

      /* set up index of column */
      colinds[0] = j;
      colvals[0] = 1.0;

      if( !SCIPisInfinity(scip,  col->ub) )
      {
         /* add artificial variable with coefficient +1.0 */
         colinds[1] = lp->ncols;
         colvals[1] = 1.0;

         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &minusinf, &col->ub, NULL, 2, &beg, colinds, colvals) );
      }

      if( !SCIPisInfinity(scip, -col->lb) )
      {
         /* add artificial variable with coefficient -1.0 */
         colinds[1] = lp->ncols;
         colvals[1] = -1.0;

         SCIP_CALL( SCIPlpiAddRows(lpi, 1, &col->lb, &plusinf, NULL, 2, &beg, colinds, colvals) );
      }
   }

   SCIPfreeBufferArray(scip, &colvals);
   SCIPfreeBufferArray(scip, &colinds);

#ifdef SCIP_OUTPUT
   SCIP_CALL( SCIPlpiWriteLP(lpi, "relativeInteriorSupNorm.lp") );
#endif

   /* solve and store point */
   SCIP_CALL( SCIPlpiSolveDual(lpi) );  /* dual is usually faster */

   if( SCIPlpiIsOptimal(lpi) )
   {
      /* get primal solution */
      SCIP_CALL( SCIPallocBufferArray(scip, &primal, ncols) );
      SCIP_CALL( SCIPlpiGetSol(lpi, &objval, primal, NULL, NULL, NULL) );
#ifndef NDEBUG
      sigma = primal[lp->ncols];
      assert( !SCIPisNegative(scip, sigma) );
#endif

      SCIPdebugMessage("Solved relative interior lp with objective %g.\n", objval);

      /* construct relative interior point */
      SCIP_CALL( SCIPcreateSol(scip, point, NULL) );
      for( j = 0; j < lp->ncols; ++j )
      {
         if( !SCIPisFeasZero(scip, primal[j]) )
         {
            SCIP_VAR* var;
            var = lp->cols[j]->var;
            SCIP_CALL( SCIPsetSolVal(scip, *point, var, primal[j]) );
         }
      }

#ifdef SCIP_DEBUG
      /* check whether the point is a relative interior point */
      for( i = 0; i < lp->nrows; ++i )
      {
         SCIP_ROW* row;
         SCIP_COL** rowcols;
         SCIP_Real* rowvals;
         SCIP_Real lhs;
         SCIP_Real rhs;
         SCIP_Real sum;
         int nnonz;

         row = lp->rows[i];
         assert( row != NULL );

         /* get row data */
         lhs = row->lhs - (SCIPisInfinity(scip, -row->lhs) ? 0.0 : row->constant);
         rhs = row->rhs - (SCIPisInfinity(scip,  row->rhs) ? 0.0 : row->constant);
         nnonz = row->nlpcols;
         assert( nnonz <= lp->ncols );
         rowcols = row->cols;
         rowvals = row->vals;

         sum = 0.0;
         for( j = 0; j < nnonz; ++j )
            sum += rowvals[j] * primal[rowcols[j]->lppos];

         /* if we have an equation */
         if( SCIPisEQ(scip, lhs, rhs) )
         {
            assert( SCIPisFeasEQ(scip, sum, lhs) );
         }
         else
         {
            assert( SCIPisInfinity(scip, lhs) || SCIPisFeasZero(scip, sigma) || SCIPisFeasGT(scip, sum, lhs) );
            assert( SCIPisInfinity(scip, rhs) || SCIPisFeasZero(scip, sigma) || SCIPisFeasLT(scip, sum, rhs) );
         }
      }
      if( inclobjcutoff )
      {
         SCIP_Real sum;
         SCIP_Real rhs;

         sum = 0.0;
         for( j = 0; j < lp->ncols; ++j )
            sum += lp->cols[j]->obj * primal[lp->cols[j]->lppos];

         rhs = SCIPgetCutoffbound(scip) - SCIPlpGetLooseObjval(lp, scip->set);
         assert( SCIPisFeasZero(scip, sigma) || SCIPisFeasLT(scip, sum, rhs) );
      }

      /* check bounds */
      for( j = 0; j < lp->ncols; ++j )
      {
         SCIP_COL* col;
         SCIP_Real val;
         
         col = lp->cols[j];
         assert( col != NULL );

         /* if the variable is not fixed */
         if( !SCIPisEQ(scip, col->lb, col->ub) )
         {
            val = primal[col->lppos];
            assert( SCIPisInfinity(scip, -col->lb) || SCIPisFeasZero(scip, sigma) || SCIPisFeasGT(scip, val, col->lb) );
            assert( SCIPisInfinity(scip,  col->ub) || SCIPisFeasZero(scip, sigma) || SCIPisFeasLT(scip, val, col->ub) );
         }
      }
#endif
  
      /* free */
      SCIPfreeBufferArray(scip, &primal);
   }
   SCIP_CALL( SCIPlpiFree(&lpi) );

   return SCIP_OKAY;
}

/** compute relative interior point to current LP
 * @see computeLPRelIntPointOneNorm
 * @see computeLPRelIntPointSupNorm
 */
SCIP_RETCODE SCIPcomputeLPRelIntPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             relaxrows,          /**< should the rows be relaxed, can be FALSE only if normtype is 's' */
   SCIP_Bool             inclobjcutoff,      /**< should a row for the objective cutoff be included */
   char                  normtype,           /**< which norm to use: 'o'ne-norm or 's'upremum-norm */
   SCIP_SOL**            point               /**< relative interior point on exit */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcomputeLPRelIntPoint", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( normtype == 'o' )
   {
      SCIP_CALL( computeLPRelIntPointOneNorm(scip, relaxrows, inclobjcutoff, point) );
   }
   else
   {
      assert(normtype == 's');
      assert(relaxrows);

      SCIP_CALL( computeLPRelIntPointSupNorm(scip, inclobjcutoff, point) );
   }

   return SCIP_OKAY;
}

/*
 * LP column methods
 */

/** returns the reduced costs of a column in the last (feasible) LP */
SCIP_Real SCIPgetColRedcost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetColRedcost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("cannot get reduced costs, because node LP is not processed\n");
      SCIPABORT();
   }

   return SCIPcolGetRedcost(col, scip->stat, scip->lp);
}


/** returns the Farkas coefficient of a column in the last (infeasible) LP */
SCIP_Real SCIPgetColFarkasCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_COL*             col                 /**< LP column */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetColFarkasCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("cannot get Farkas coeff, because node LP is not processed\n");
      SCIPABORT();
   }

   return SCIPcolGetFarkasCoef(col, scip->stat, scip->lp);
}


/*
 * LP row methods
 */

/** creates and captures an LP row */
SCIP_RETCODE SCIPcreateRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   int                   len,                /**< number of nonzeros in the row */
   SCIP_COL**            cols,               /**< array with columns of row entries */
   SCIP_Real*            vals,               /**< array with coefficients of row entries */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, len, cols, vals, lhs, rhs, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** creates and captures an LP row without any coefficients */
SCIP_RETCODE SCIPcreateEmptyRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row,                /**< pointer to row */
   const char*           name,               /**< name of row */
   SCIP_Real             lhs,                /**< left hand side of row */
   SCIP_Real             rhs,                /**< right hand side of row */
   SCIP_Bool             local,              /**< is row only valid locally? */
   SCIP_Bool             modifiable,         /**< is row modifiable during node processing (subject to column generation)? */
   SCIP_Bool             removable           /**< should the row be removed from the LP due to aging or cleanup? */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateEmptyRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCreate(row, scip->mem->probmem, scip->set, scip->stat,
         name, 0, NULL, NULL, lhs, rhs, local, modifiable, removable) );

   return SCIP_OKAY;
}

/** increases usage counter of LP row */
SCIP_RETCODE SCIPcaptureRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row to capture */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcaptureRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIProwCapture(row);

   return SCIP_OKAY;
}

/** decreases usage counter of LP row, and frees memory if necessary */
SCIP_RETCODE SCIPreleaseRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW**            row                 /**< pointer to LP row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPreleaseRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE) );

   SCIP_CALL( SCIProwRelease(row, scip->mem->probmem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** changes left hand side of LP row */
SCIP_RETCODE SCIPchgRowLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgRowLhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwChgLhs(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of LP row */
SCIP_RETCODE SCIPchgRowRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgRowRhs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwChgRhs(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp, rhs) );

   return SCIP_OKAY;
}

/** informs row, that all subsequent additions of variables to the row should be cached and not directly applied;
 *  after all additions were applied, SCIPflushRowExtensions() must be called;
 *  while the caching of row extensions is activated, information methods of the row give invalid results;
 *  caching should be used, if a row is build with SCIPaddVarToRow() calls variable by variable to increase
 *  the performance
 */
SCIP_RETCODE SCIPcacheRowExtensions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcacheRowExtension", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   return SCIP_OKAY;
}

/** flushes all cached row extensions after a call of SCIPcacheRowExtensions() and merges coefficients with
 *  equal columns into a single coefficient
 */
SCIP_RETCODE SCIPflushRowExtensions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPflushRowExtension", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* force the row sorting, and merge equal column entries */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** resolves variable to columns and adds them with the coefficient to the row */
SCIP_RETCODE SCIPaddVarToRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value of coefficient */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddVarToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPvarAddToRow(var, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lp, row, val) );

   return SCIP_OKAY;
}

/** resolves variables to columns and adds them with the coefficients to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
SCIP_RETCODE SCIPaddVarsToRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real*            vals                /**< values of coefficients */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddVarsToRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   SCIP_CALL( SCIProwEnsureSize(row, scip->mem->probmem, scip->set, SCIProwGetNNonz(row) + nvars) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarAddToRow(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lp,
            row, vals[v]) );
   }

   /* force the row sorting */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** resolves variables to columns and adds them with the same single coefficient to the row;
 *  this method caches the row extensions and flushes them afterwards to gain better performance
 */
SCIP_RETCODE SCIPaddVarsToRowSameCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real             val                 /**< unique value of all coefficients */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddVarsToRowSameCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* resize the row to be able to store all variables (at least, if they are COLUMN variables) */
   SCIP_CALL( SCIProwEnsureSize(row, scip->mem->probmem, scip->set, SCIProwGetNNonz(row) + nvars) );

   /* delay the row sorting */
   SCIProwDelaySort(row);

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPvarAddToRow(vars[v], scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->transprob, scip->lp,
            row, val) );
   }

   /* force the row sorting */
   SCIProwForceSort(row, scip->set);

   return SCIP_OKAY;
}

/** tries to find a value, such that all row coefficients, if scaled with this value become integral */
SCIP_RETCODE SCIPcalcRowIntegralScalar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal allowed scalar */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Real*            intscalar,          /**< pointer to store scalar that would make the coefficients integral, or NULL */
   SCIP_Bool*            success             /**< stores whether returned value is valid */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcalcRowIntegralScalar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwCalcIntegralScalar(row, scip->set, mindelta, maxdelta, maxdnom, maxscale,
         usecontvars, intscalar, success) );

   return SCIP_OKAY;
}

/** tries to scale row, s.t. all coefficients (of integer variables) become integral */
SCIP_RETCODE SCIPmakeRowIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta,           /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Longint          maxdnom,            /**< maximal denominator allowed in rational numbers */
   SCIP_Real             maxscale,           /**< maximal value to scale row with */
   SCIP_Bool             usecontvars,        /**< should the coefficients of the continuous variables also be made integral? */
   SCIP_Bool*            success             /**< stores whether row could be made rational */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPmakeRowIntegral", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIProwMakeIntegral(row, scip->mem->probmem, scip->set, scip->eventqueue, scip->stat, scip->lp, mindelta, maxdelta, maxdnom, maxscale,
         usecontvars, success) );

   return SCIP_OKAY;
}

/** returns minimal absolute value of row vector's non-zero coefficients */
SCIP_Real SCIPgetRowMinCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowMinCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetMinval(row, scip->set);
}

/** returns maximal absolute value of row vector's non-zero coefficients */
SCIP_Real SCIPgetRowMaxCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowMaxCoef", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetMaxval(row, scip->set);
}

/** returns the minimal activity of a row w.r.t. the column's bounds */
SCIP_Real SCIPgetRowMinActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowMinActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetMinActivity(row, scip->set, scip->stat);
}

/** returns the maximal activity of a row w.r.t. the column's bounds */
SCIP_Real SCIPgetRowMaxActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowMaxActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetMaxActivity(row, scip->set, scip->stat);
}

/** recalculates the activity of a row in the last LP solution */
SCIP_RETCODE SCIPrecalcRowLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrecalcRowLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIProwRecalcLPActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row in the last LP solution */
SCIP_Real SCIPgetRowLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetLPActivity(row, scip->set, scip->stat, scip->lp);
}

/** returns the feasibility of a row in the last LP solution: negative value means infeasibility */
SCIP_Real SCIPgetRowLPFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowLPFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetLPFeasibility(row, scip->set, scip->stat, scip->lp);
}

/** recalculates the activity of a row for the current pseudo solution */
SCIP_RETCODE SCIPrecalcRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrecalcRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIProwRecalcPseudoActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row for the current pseudo solution */
SCIP_Real SCIPgetRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetPseudoActivity(row, scip->set, scip->stat);
}

/** returns the feasibility of a row for the current pseudo solution: negative value means infeasibility */
SCIP_Real SCIPgetRowPseudoFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowPseudoFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIProwGetPseudoFeasibility(row, scip->set, scip->stat);
}

/** recalculates the activity of a row in the last LP or pseudo solution */
SCIP_RETCODE SCIPrecalcRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrecalcRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      SCIProwRecalcLPActivity(row, scip->stat);
   else
      SCIProwRecalcPseudoActivity(row, scip->stat);

   return SCIP_OKAY;
}

/** returns the activity of a row in the last LP or pseudo solution */
SCIP_Real SCIPgetRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPActivity(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoActivity(row, scip->set, scip->stat);
}

/** returns the feasibility of a row in the last LP or pseudo solution */
SCIP_Real SCIPgetRowFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< LP row */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPFeasibility(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoFeasibility(row, scip->set, scip->stat);
}

/** returns the activity of a row for the given primal solution */
SCIP_Real SCIPgetRowSolActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowSolActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( sol != NULL )
      return SCIProwGetSolActivity(row, scip->set, scip->stat, sol);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPActivity(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoActivity(row, scip->set, scip->stat);
}

/** returns the feasibility of a row for the given primal solution */
SCIP_Real SCIPgetRowSolFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
      return SCIProwGetSolFeasibility(row, scip->set, scip->stat, sol);
   else if( SCIPtreeHasCurrentNodeLP(scip->tree) )
      return SCIProwGetLPFeasibility(row, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetPseudoFeasibility(row, scip->set, scip->stat);
}

/** output row to file stream */
SCIP_RETCODE SCIPprintRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(row != NULL);

   SCIP_CALL( checkStage(scip, "SCIPprintRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   SCIProwPrint(row, file);

   return SCIP_OKAY;
}


/*
 * NLP methods
 */

/** marks that there are constraints that are representable by nonlinear constraints involving continuous variables
 * This method should be called by a constraint handler if it has constraints that have a representation as nonlinear rows
 * and that are still nonlinear after fixing all discrete variables in the CIP.
 * 
 * The function should be called before the branch-and-bound process is initialized, e.g., when presolve is exiting.
 * 
 */ 
void SCIPmarkContinuousNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPmarkContinuousNonlinearitiesPresent", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   
   scip->set->continnonlinpresent = TRUE;
   scip->set->nonlinearitypresent = TRUE;
}

/** marks that there are constraints that are representable by nonlinear constraints
 * This method should be called by a constraint handler if it has constraints that have a representation as nonlinear rows.
 * 
 * The function should be called before the branch-and-bound process is initialized, e.g., when presolve is exiting.
 * 
 * Calling SCIPmarkContinuousNonlinearitiesPresent makes a call to SCIPmarkNonlinearitiesPresent dispensable.
 */ 
void SCIPmarkNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPmarkNonlinearitiesPresent", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   
   scip->set->nonlinearitypresent = TRUE;
}

/** returns whether constraints representable as nonlinear rows are present that involve continuous nonlinear variables
 * @see SCIPmarkContinuousNonlinearitiesPresent
 */
SCIP_Bool SCIPhasContinuousNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPhasContinuousNonlinearitiesPresent", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   
   return scip->set->continnonlinpresent;
}

/** returns whether constraints representable as nonlinear rows are present
 * @see SCIPmarkNonlinearitiesPresent
 */
SCIP_Bool SCIPhasNonlinearitiesPresent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPhasNonlinearitiesPresent", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   
   return scip->set->nonlinearitypresent;
}

/** returns, whether an NLP has been constructed */
SCIP_Bool SCIPisNLPConstructed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisNLPConstructed", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   return (scip->nlp != NULL);
}

/** gets current NLP variables along with the current number of NLP variables */
SCIP_RETCODE SCIPgetNLPVarsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to store the array of NLP variables, or NULL */
   int*                  nvars               /**< pointer to store the number of NLP variables, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPVarsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      if( vars != NULL )
         *vars = SCIPnlpGetVars(scip->nlp);
      if( nvars != NULL )
         *nvars = SCIPnlpGetNVars(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets array with variables of the NLP */
SCIP_VAR** SCIPgetNLPVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPVars", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetVars(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return NULL;
}

/** gets current number of variables in NLP */
int SCIPgetNNLPVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNLPVars", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetNVars(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return 0;
}

/** computes for each variables the number of NLP rows in which the variable appears in a nonlinear var */
SCIP_RETCODE SCIPgetNLPVarsNonlinearity(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  nlcount             /**< an array of length at least SCIPnlpGetNVars() to store nonlinearity counts of variables */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPVarsNonlinearity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpGetVarsNonlinearity(scip->nlp, nlcount) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gives dual solution values associated with lower bounds of NLP variables */
SCIP_Real* SCIPgetNLPVarsLbDualsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPVarsLbDualsol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetVarsLbDualsol(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return NULL;
}

/** gives dual solution values associated with upper bounds of NLP variables */
SCIP_Real* SCIPgetNLPVarsUbDualsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPVarsUbDualsol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetVarsUbDualsol(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return NULL;
}

/** gets current NLP nonlinear rows along with the current number of NLP nonlinear rows */
SCIP_RETCODE SCIPgetNLPNlRowsData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW***         nlrows,             /**< pointer to store the array of NLP nonlinear rows, or NULL */
   int*                  nnlrows             /**< pointer to store the number of NLP nonlinear rows, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPNlRowsData", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      if( nlrows != NULL )
         *nlrows = SCIPnlpGetNlRows(scip->nlp);
      if( nnlrows != NULL )
         *nnlrows = SCIPnlpGetNNlRows(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets array with nonlinear rows of the NLP */
SCIP_NLROW** SCIPgetNLPNlRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPNlRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetNlRows(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return NULL;
}

/** gets current number of nonlinear rows in NLP */
int SCIPgetNNLPNlRows(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNLPNlRows", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetNNlRows(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return 0;
}

/** adds a nonlinear row to the NLP
 * row is captured by NLP */
SCIP_RETCODE SCIPaddNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< nonlinear row to add to NLP */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddNlRow", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpAddNlRow(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, nlrow) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** makes sure that the NLP of the current node is flushed */
SCIP_RETCODE SCIPflushNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPflushNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpFlush(scip->nlp, scip->mem->probmem, scip->set) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** sets or clears initial primal guess for NLP solution (start point for NLP solver) */
SCIP_RETCODE SCIPsetNLPInitialGuess(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            initialguess        /**< values of initial guess (corresponding to variables from SCIPgetNLPVarsData), or NULL to use no start point */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNLPInitialGuess", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpSetInitialGuess(scip->nlp, SCIPblkmem(scip), initialguess) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** sets initial primal guess for NLP solution (start point for NLP solver) */
SCIP_RETCODE SCIPsetNLPInitialGuessSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution which values should be taken as initial guess, or NULL for LP solution */
   )
{
   SCIP_Real* vals;

   SCIP_CALL( checkStage(scip, "SCIPsetNLPInitialGuessSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, SCIPnlpGetNVars(scip->nlp)) );
      SCIP_CALL( SCIPgetSolVals(scip, sol, SCIPnlpGetNVars(scip->nlp), SCIPnlpGetVars(scip->nlp), vals) );
      SCIP_CALL( SCIPnlpSetInitialGuess(scip->nlp, SCIPblkmem(scip), vals) );
      SCIPfreeBufferArray(scip, &vals);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** solves the current NLP */
SCIP_RETCODE SCIPsolveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsolveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpSolve(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets solution status of current NLP */
SCIP_NLPSOLSTAT SCIPgetNLPSolstat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPSolstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetSolstat(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return SCIP_NLPSOLSTAT_UNKNOWN;
}

/** gets termination status of last NLP solve */
SCIP_NLPTERMSTAT SCIPgetNLPTermstat(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPTermstat", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetTermstat(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return SCIP_NLPTERMSTAT_OTHER;
}

/** gives statistics (number of iterations, solving time, ...) of last NLP solve */
SCIP_RETCODE SCIPgetNLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPSTATISTICS*   statistics          /**< pointer to store statistics */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPStatistics", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetStatistics(scip->nlp, statistics);
   }

   SCIPerrorMessage("NLP has not been not constructed.\n");
   return SCIP_ERROR;
}

/** gets objective value of current NLP */
SCIP_Real SCIPgetNLPObjval(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPObjval", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpGetObjval(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_INVALID;
   }
}

/** indicates whether a feasible solution for the current NLP is available
 * thus, returns whether the solution status <= feasible  */
SCIP_Bool SCIPhasNLPSolution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPhasNLPSolution", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      return SCIPnlpHasSolution(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      SCIPABORT();
   }

   return FALSE;
}

/** gets fractional variables of last NLP solution along with solution values and fractionalities */
SCIP_RETCODE SCIPgetNLPFracVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           fracvars,           /**< pointer to store the array of NLP fractional variables, or NULL */
   SCIP_Real**           fracvarssol,        /**< pointer to store the array of NLP fractional variables solution values, or NULL */
   SCIP_Real**           fracvarsfrac,       /**< pointer to store the array of NLP fractional variables fractionalities, or NULL */
   int*                  nfracvars,          /**< pointer to store the number of NLP fractional variables , or NULL */
   int*                  npriofracvars       /**< pointer to store the number of NLP fractional variables with maximal branching priority, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPFracVars", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpGetFracVars(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, fracvars, fracvarssol, fracvarsfrac, nfracvars, npriofracvars) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets integer parameter of NLP */
SCIP_RETCODE SCIPgetNLPIntPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPIntPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpGetIntPar(scip->nlp, type, ival) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** sets integer parameter of NLP */
SCIP_RETCODE SCIPsetNLPIntPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNLPIntPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpSetIntPar(scip->nlp, type, ival) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets floating point parameter of NLP */
SCIP_RETCODE SCIPgetNLPRealPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPRealPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpGetRealPar(scip->nlp, type, dval) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of NLP */
SCIP_RETCODE SCIPsetNLPRealPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNLPRealPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpSetRealPar(scip->nlp, type, dval) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets string parameter of NLP */
SCIP_RETCODE SCIPgetNLPStringPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char**          sval                /**< pointer to store the parameter value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNLPStringPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpGetStringPar(scip->nlp, type, sval) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** sets string parameter of NLP */
SCIP_RETCODE SCIPsetNLPStringPar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPPARAM         type,               /**< parameter number */
   const char*           sval                /**< parameter value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNLPStringPar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpSetStringPar(scip->nlp, type, sval) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** writes current NLP to a file */
SCIP_RETCODE SCIPwriteNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           fname               /**< file name */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPwriteNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpWrite(scip->nlp, scip->set, fname) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** gets the NLP interface and problem used by the SCIP NLP;
 *  with the NLPI and its problem you can use all of the methods defined in nlpi/nlpi.h;
 *
 *  @warning You have to make sure, that the full internal state of the NLPI does not change or is recovered completely
 *           after the end of the method that uses the NLPI. In particular, if you manipulate the NLP or its solution
 *           (e.g. by calling one of the SCIPnlpiAdd...() or the SCIPnlpiSolve() method), you have to check in advance
 *           whether the NLP is currently solved.  If this is the case, you have to make sure, the internal solution
 *           status is recovered completely at the end of your method. Additionally you have to resolve the NLP with
 *           SCIPnlpiSolve() in order to reinstall the internal solution status.
 *
 *  @warning Make also sure, that all parameter values that you have changed are set back to their original values.
 */
SCIP_RETCODE SCIPgetNLPI(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLPI**           nlpi,               /**< pointer to store the NLP solver interface */
   SCIP_NLPIPROBLEM**    nlpiproblem         /**< pointer to store the NLP solver interface problem */
   )
{
   assert(nlpi != NULL);
   assert(nlpiproblem != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetNLPI", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      *nlpi = SCIPnlpGetNLPI(scip->nlp);
      *nlpiproblem = SCIPnlpGetNLPIProblem(scip->nlp);
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}


/*
 * NLP diving methods
 */

/**@name NLP Diving Methods */
/**@{ */

/** initiates NLP diving
 * making methods SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), SCIPchgVarsBoundsDiveNLP(), and SCIPsolveDiveNLP() available */
SCIP_RETCODE SCIPstartDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstartDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpStartDive(scip->nlp, SCIPblkmem(scip), scip->set) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** ends NLP diving
 * resets changes made by SCIPchgVarObjDiveNLP(), SCIPchgVarBoundsDiveNLP(), and SCIPchgVarsBoundsDiveNLP() */
SCIP_RETCODE SCIPendDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPendDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpEndDive(scip->nlp, SCIPblkmem(scip), scip->set) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** changes linear objective coefficient of a variable in diving NLP */
SCIP_RETCODE SCIPchgVarObjDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which coefficient to change */
   SCIP_Real             coef                /**< new value for coefficient */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarObjDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpChgVarObjDive(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat, var, coef) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** changes bounds of a variable in diving NLP */
SCIP_RETCODE SCIPchgVarBoundsDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable which bounds to change */
   SCIP_Real             lb,                 /**< new lower bound */
   SCIP_Real             ub                  /**< new upper bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarBoundsDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpChgVarBoundsDive(scip->nlp, var, lb, ub) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** changes bounds of a set of variables in diving NLP */
SCIP_RETCODE SCIPchgVarsBoundsDiveNLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables which bounds to changes */
   SCIP_VAR**            vars,               /**< variables which bounds to change */
   SCIP_Real*            lbs,                /**< new lower bounds */
   SCIP_Real*            ubs                 /**< new upper bounds */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarsBoundsDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpChgVarsBoundsDive(scip->nlp, scip->set, nvars, vars, lbs, ubs) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** solves diving NLP */
SCIP_RETCODE SCIPsolveDiveNLP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsolveDiveNLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL )
   {
      SCIP_CALL( SCIPnlpSolveDive(scip->nlp, SCIPblkmem(scip), scip->set, scip->stat) );
   }
   else
   {
      SCIPerrorMessage("NLP has not been not constructed.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/**@} */


/*
 * NLP nonlinear row methods
 */

/** creates and captures an NLP row */
SCIP_RETCODE SCIPcreateNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             constant,           /**< constant */
   int                   nlinvars,           /**< number of linear variables */
   SCIP_VAR**            linvars,            /**< linear variables, or NULL if nlinvars == 0 */
   SCIP_Real*            lincoefs,           /**< linear coefficients, or NULL if nlinvars == 0 */
   int                   nquadvars,          /**< number variables in quadratic terms */
   SCIP_VAR**            quadvars,           /**< variables in quadratic terms, or NULL if nquadvars == 0 */
   int                   nquadelems,         /**< number of elements in quadratic term */
   SCIP_QUADELEM*        quadelems,          /**< elements (i.e., monomials) in quadratic term, or NULL if nquadelems == 0 */
   SCIP_EXPRTREE*        expression,         /**< nonlinear expression, or NULL */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowCreate(nlrow, scip->mem->probmem, scip->set,
         name, constant, nlinvars, linvars, lincoefs, nquadvars, quadvars, nquadelems, quadelems, expression, lhs, rhs) );

   return SCIP_OKAY;
}

/** creates and captures an NLP nonlinear row without any coefficients */
SCIP_RETCODE SCIPcreateEmptyNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   const char*           name,               /**< name of nonlinear row */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateEmptyNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowCreate(nlrow, scip->mem->probmem, scip->set,
         name, 0.0, 0, NULL, NULL, 0, NULL, 0, NULL, NULL, lhs, rhs) );

   return SCIP_OKAY;
}

/** creates and captures an NLP row from a linear row */
SCIP_RETCODE SCIPcreateNlRowFromRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow,              /**< buffer to store pointer to nonlinear row */
   SCIP_ROW*             row                 /**< the linear row to copy */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateNlRowFromRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowCreateFromRow(nlrow, scip->mem->probmem, scip->set, row) );

   return SCIP_OKAY;
}

/** increases usage counter of NLP nonlinear row */
SCIP_RETCODE SCIPcaptureNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< nonlinear row to capture */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcaptureNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnlrowCapture(nlrow);

   return SCIP_OKAY;
}

/** decreases usage counter of NLP nonlinear row, and frees memory if necessary */
SCIP_RETCODE SCIPreleaseNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW**          nlrow               /**< nonlinear row to release */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPreleaseNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE) );

   SCIP_CALL( SCIPnlrowRelease(nlrow, scip->mem->probmem, scip->set) );

   return SCIP_OKAY;
}

/** changes left hand side of NLP nonlinear row */
SCIP_RETCODE SCIPchgNlRowLhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real             lhs                 /**< new left hand side */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgNlRowLhs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgLhs(nlrow, scip->set, scip->stat, scip->nlp, lhs) );

   return SCIP_OKAY;
}

/** changes right hand side of NLP nonlinear row */
SCIP_RETCODE SCIPchgNlRowRhs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real             rhs                 /**< new right hand side */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgNlRowRhs", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgRhs(nlrow, scip->set, scip->stat, scip->nlp, rhs) );

   return SCIP_OKAY;
}

/** changes constant of NLP nonlinear row */
SCIP_RETCODE SCIPchgNlRowConstant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real             constant            /**< new value for constant */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgNlRowConstant", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgConstant(nlrow, scip->set, scip->stat, scip->nlp, constant) );

   return SCIP_OKAY;
}

/** adds variable with a linear coefficient to the nonlinear row */
SCIP_RETCODE SCIPaddLinearCoefToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_Real             val                 /**< value of coefficient in linear part of row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddLinearCoefToNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowAddLinearCoef(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, var, val) );

   return SCIP_OKAY;
}

/** adds variables with linear coefficients to the row */
SCIP_RETCODE SCIPaddLinearCoefsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nvars,              /**< number of variables to add to the row */
   SCIP_VAR**            vars,               /**< problem variables to add */
   SCIP_Real*            vals                /**< values of coefficients in linear part of row */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddLinearCoefsToNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* add the variables to the row */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPnlrowAddLinearCoef(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, vars[v], vals[v]) );
   }

   return SCIP_OKAY;
}

/** changes linear coefficient of a variables in a row
 * setting the coefficient to 0.0 means that it is removed from the row
 * the variable does not need to exists before */
SCIP_RETCODE SCIPchgNlRowLinearCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var,                /**< variable */
   SCIP_Real             coef                /**< new value of coefficient */
   )
{
   assert(var != NULL);

   SCIP_CALL( checkStage(scip, "SCIPchgNlRowLinearCoef", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgLinearCoef(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, var, coef) );

   return SCIP_OKAY;
}

/** adds quadratic variable to the nonlinear row
 * after adding a quadratic variable, it can be used to add quadratic elements */
SCIP_RETCODE SCIPaddQuadVarToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_VAR*             var                 /**< problem variable */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddQuadVarToNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, scip->mem->probmem, scip->set, var) );

   return SCIP_OKAY;
}

/** adds quadratic variables to the nonlinear row
 * after adding quadratic variables, they can be used to add quadratic elements */
SCIP_RETCODE SCIPaddQuadVarsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nvars,              /**< number of problem variables */
   SCIP_VAR**            vars                /**< problem variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddQuadVarsToNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowEnsureQuadVarsSize(nlrow, scip->mem->probmem, scip->set, SCIPnlrowGetNQuadVars(nlrow) + nvars) );
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPnlrowAddQuadVar(nlrow, scip->mem->probmem, scip->set, vars[v]) );
   }

   return SCIP_OKAY;
}

/** add a quadratic element to the nonlinear row
 * variable indices of the quadratic element need to be relative to quadratic variables array of row */
SCIP_RETCODE SCIPaddQuadElementToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_QUADELEM         quadelem            /**< quadratic element */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddQuadElementToNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, quadelem) );

   return SCIP_OKAY;
}

/** adds quadratic elements to the nonlinear row
 * variable indices of the quadratic elements need to be relative to quadratic variables array of row */
SCIP_RETCODE SCIPaddQuadElementsToNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   nquadelems,         /**< number of quadratic elements */
   SCIP_QUADELEM*        quadelems           /**< quadratic elements */
   )
{
   int v;

   assert(nquadelems == 0 || quadelems != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddQuadElementsToNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowEnsureQuadElementsSize(nlrow, scip->mem->probmem, scip->set, SCIPnlrowGetNQuadElems(nlrow) + nquadelems) );
   for( v = 0; v < nquadelems; ++v )
   {
      SCIP_CALL( SCIPnlrowAddQuadElement(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, quadelems[v]) );
   }

   return SCIP_OKAY;
}

/** changes coefficient in quadratic part of a row
 * setting the coefficient in the quadelement to 0.0 means that it is removed from the row
 * the element does not need to exists before */
SCIP_RETCODE SCIPchgNlRowQuadElement(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_QUADELEM         quadelement         /**< new quadratic element, or update for existing one */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgNlRowQuadElement", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgQuadElem(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, quadelement) );

   return SCIP_OKAY;
}

/** sets or deletes expression tree in the nonlinear row */
SCIP_RETCODE SCIPsetNlRowExprtree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_EXPRTREE*        exprtree            /**< expression tree, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNlRowExprtree", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgExprtree(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, exprtree) );

   return SCIP_OKAY;
}

/** sets a parameter of expression tree in the nonlinear row */
SCIP_RETCODE SCIPsetNlRowExprtreeParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   int                   paramidx,           /**< index of parameter in expression tree */
   SCIP_Real             paramval            /**< new value of parameter in expression tree */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNlRowExprtreeParam", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgExprtreeParam(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, paramidx, paramval) );

   return SCIP_OKAY;
}

/** sets parameters of expression tree in the nonlinear row */
SCIP_RETCODE SCIPsetNlRowExprtreeParams(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real*            paramvals           /**< new values of parameter in expression tree */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetNlRowExprtreeParams", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowChgExprtreeParams(nlrow, scip->mem->probmem, scip->set, scip->stat, scip->nlp, paramvals) );

   return SCIP_OKAY;
}

/** recalculates the activity of a nonlinear row in the last NLP solution */
SCIP_RETCODE SCIPrecalcNlRowNLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrecalcNlRowNLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("do not have NLP for computing NLP activity\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlrowRecalcNLPActivity(nlrow, scip->set, scip->stat, scip->nlp) );

   return SCIP_OKAY;
}

/** returns the activity of a nonlinear row in the last NLP solution */
SCIP_RETCODE SCIPgetNlRowNLPActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetMlRowNLPActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("do not have NLP for computing NLP activity\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, scip->set, scip->stat, scip->nlp, activity) );

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row in the last NLP solution: negative value means infeasibility */
SCIP_RETCODE SCIPgetNlRowNLPFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowNLPFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp == NULL )
   {
      SCIPerrorMessage("do not have NLP for computing NLP feasibility\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, scip->set, scip->stat, scip->nlp, feasibility) );

   return SCIP_OKAY;
}

/** recalculates the activity of a nonlinear row for the current pseudo solution */
SCIP_RETCODE SCIPrecalcNlRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrecalcNlRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowRecalcPseudoActivity(nlrow, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** gives the activity of a nonlinear row for the current pseudo solution */
SCIP_RETCODE SCIPgetNlRowPseudoActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            pseudoactivity      /**< buffer to store pseudo activity value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowPseudoActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, scip->set, scip->stat, pseudoactivity) );

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row for the current pseudo solution: negative value means infeasibility */
SCIP_RETCODE SCIPgetNlRowPseudoFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            pseudofeasibility   /**< buffer to store pseudo feasibility value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowPseudoFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowGetPseudoFeasibility(nlrow, scip->set, scip->stat, pseudofeasibility) );

   return SCIP_OKAY;
}

/** recalculates the activity of a nonlinear row in the last NLP or pseudo solution */
SCIP_RETCODE SCIPrecalcNlRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow               /**< NLP nonlinear row */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrecalcNlRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowRecalcNLPActivity(nlrow, scip->set, scip->stat, scip->nlp) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowRecalcPseudoActivity(nlrow, scip->set, scip->stat) );
   }

   return SCIP_OKAY;
}

/** gives the activity of a nonlinear row in the last NLP or pseudo solution */
SCIP_RETCODE SCIPgetNlRowActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, scip->set, scip->stat, scip->nlp, activity) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, scip->set, scip->stat, activity) );
   }

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row in the last NLP or pseudo solution */
SCIP_RETCODE SCIPgetNlRowFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, scip->set, scip->stat, scip->nlp, feasibility) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoFeasibility(nlrow, scip->set, scip->stat, feasibility) );
   }

   return SCIP_OKAY;
}

/** gives the activity of a nonlinear row for the given primal solution or NLP solution or pseudo solution */
SCIP_RETCODE SCIPgetNlRowSolActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for NLP solution of pseudo solution */
   SCIP_Real*            activity            /**< buffer to store activity value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowSolActivity", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
   {
      SCIP_CALL( SCIPnlrowGetSolActivity(nlrow, scip->set, scip->stat, sol, activity) );
   }
   else if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPActivity(nlrow, scip->set, scip->stat, scip->nlp, activity) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoActivity(nlrow, scip->set, scip->stat, activity) );
   }

   return SCIP_OKAY;
}

/** gives the feasibility of a nonlinear row for the given primal solution */
SCIP_RETCODE SCIPgetNlRowSolFeasibility(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP nonlinear row */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Real*            feasibility         /**< buffer to store feasibility value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowSolFeasibility", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   if( sol != NULL )
   {
      SCIP_CALL( SCIPnlrowGetSolFeasibility(nlrow, scip->set, scip->stat, sol, feasibility) );
   }
   else if( scip->nlp != NULL && SCIPnlpHasCurrentNodeNLP(scip->nlp) && SCIPnlpHasSolution(scip->nlp) )
   {
      SCIP_CALL( SCIPnlrowGetNLPFeasibility(nlrow, scip->set, scip->stat, scip->nlp, feasibility) );
   }
   else
   {
      SCIP_CALL( SCIPnlrowGetPseudoFeasibility(nlrow, scip->set, scip->stat, feasibility) );
   }

   return SCIP_OKAY;
}

/** gives the minimal and maximal activity of a nonlinear row w.r.t. the variable's bounds */
SCIP_RETCODE SCIPgetNlRowActivityBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   SCIP_Real*            minactivity,        /**< buffer to store minimal activity, or NULL */
   SCIP_Real*            maxactivity         /**< buffer to store maximal activity, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetNlRowActivityBounds", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowGetActivityBounds(nlrow, scip->set, scip->stat, minactivity, maxactivity) );

   return SCIP_OKAY;
}

/** output nonlinear row to file stream */
SCIP_RETCODE SCIPprintNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLROW*           nlrow,              /**< NLP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(nlrow != NULL);

   SCIP_CALL( checkStage(scip, "SCIPprintNlRow", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnlrowPrint(nlrow, file) );

   return SCIP_OKAY;
}

/**@name Expression tree methods */
/**@{ */

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)

/** replaced array of variables in expression tree by corresponding transformed variables */
SCIP_RETCODE SCIPgetExprtreeTransformedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{
   assert(scip != NULL);
   assert(tree != NULL);

   if( SCIPexprtreeGetNVars(tree) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetTransformedVars(scip, SCIPexprtreeGetNVars(tree), SCIPexprtreeGetVars(tree), SCIPexprtreeGetVars(tree)) );

   return SCIP_OKAY;
}

/** evaluates an expression tree for a primal solution or LP solution */
SCIP_RETCODE SCIPevalExprtreeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SOL*             sol,                /**< a solution, or NULL for current LP solution */
   SCIP_Real*            val                 /**< buffer to store value */
   )
{
   SCIP_Real* varvals;
   int nvars;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   nvars = SCIPexprtreeGetNVars(tree);
   
   if( nvars == 0 )
   {
      SCIP_CALL( SCIPexprtreeEval(tree, NULL, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, nvars, SCIPexprtreeGetVars(tree), varvals) );

   SCIP_CALL( SCIPexprtreeEval(tree, varvals, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current global bounds */
SCIP_RETCODE SCIPevalExprtreeGlobalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
   )
{
   SCIP_INTERVAL* varvals;
   SCIP_VAR**     vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   nvars = SCIPexprtreeGetNVars(tree);

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, NULL, val) );
      return SCIP_OKAY;
   }
   
   vars = SCIPexprtreeGetVars(tree);
   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -SCIPvarGetLbGlobal(vars[i])),  /*lint !e666*/
          infty2infty(SCIPinfinity(scip), infinity,  SCIPvarGetUbGlobal(vars[i]))); /*lint !e666*/
   }

   SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, varvals, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current local bounds */
SCIP_RETCODE SCIPevalExprtreeLocalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
   )
{
   SCIP_INTERVAL* varvals;
   SCIP_VAR**     vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   nvars = SCIPexprtreeGetNVars(tree);

   if( nvars == 0 )
   {
      SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, NULL, val) );
      return SCIP_OKAY;
   }

   vars = SCIPexprtreeGetVars(tree);
   assert(vars != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      /* due to numerics, the lower bound on a variable in SCIP can be slightly higher than the upper bound
       * in this case, we take the most conservative way and switch the bounds
       * further, we translate SCIP's value for infinity to the users value for infinity
       */
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -MIN(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i]))),  /*lint !e666*/
          infty2infty(SCIPinfinity(scip), infinity,  MAX(SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])))); /*lint !e666*/
   }

   SCIP_CALL( SCIPexprtreeEvalInt(tree, infinity, varvals, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

#undef infty2infty

/**@} */

/*
 * cutting plane methods
 */

/** returns efficacy of the cut with respect to the given primal solution or the current LP solution:
 *  e = -feasibility/norm
 */
SCIP_Real SCIPgetCutEfficacy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetCutEfficacy", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( sol == NULL )
      return SCIProwGetLPEfficacy(cut, scip->set, scip->stat, scip->lp);
   else
      return SCIProwGetSolEfficacy(cut, scip->set, scip->stat, sol);
}

/** returns whether the cut's efficacy with respect to the given primal solution or the current LP solution is greater
 *  than the minimal cut efficacy
 */
SCIP_Bool SCIPisCutEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution, or NULL for current LP solution */
   SCIP_ROW*             cut                 /**< separated cut */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisCutEfficacious", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( sol == NULL )
      return SCIProwIsLPEfficacious(cut, scip->set, scip->stat, scip->lp, (SCIPtreeGetCurrentDepth(scip->tree) == 0));
   else
      return SCIProwIsSolEfficacious(cut, scip->set, scip->stat, sol, (SCIPtreeGetCurrentDepth(scip->tree) == 0));
}

/** checks, if the given cut's efficacy is larger than the minimal cut efficacy */
SCIP_Bool SCIPisEfficacious(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             efficacy            /**< efficacy of the cut */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisEfficacious", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsetIsEfficacious(scip->set, (SCIPtreeGetCurrentDepth(scip->tree) == 0), efficacy);
}

/** calculates the efficacy norm of the given vector, which depends on the "separating/efficacynorm" parameter */
SCIP_Real SCIPgetVectorEfficacyNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            vals,               /**< array of values */
   int                   nvals               /**< number of values */
   )
{
   SCIP_Real norm;
   int i;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVectorEfficacyNorm", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   norm = 0.0;
   switch( scip->set->sepa_efficacynorm )
   {
   case 'e':
      for( i = 0; i < nvals; ++i )
         norm += SQR(vals[i]);
      norm = SQRT(norm);
      break;
   case 'm':
      for( i = 0; i < nvals; ++i )
      {
         SCIP_Real absval;

         absval = REALABS(vals[i]);
         norm = MAX(norm, absval);
      }
      break;
   case 's':
      for( i = 0; i < nvals; ++i )
         norm += REALABS(vals[i]);
      break;
   case 'd':
      for( i = 0; i < nvals; ++i )
      {
         if( !SCIPisZero(scip, vals[i]) )
         {
            norm = 1.0;
            break;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", scip->set->sepa_efficacynorm);
      assert(FALSE);
   }

   return norm;
}

/** adds cut to separation storage */
SCIP_RETCODE SCIPaddCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that was separated, or NULL for LP solution */
   SCIP_ROW*             cut,                /**< separated cut */
   SCIP_Bool             forcecut            /**< should the cut be forced to enter the LP? */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(SCIPtreeGetCurrentNode(scip->tree) != NULL);

   SCIP_CALL( SCIPsepastoreAddCut(scip->sepastore, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
         scip->eventfilter, scip->lp, sol, cut, forcecut, (SCIPtreeGetCurrentDepth(scip->tree) == 0)) );

   return SCIP_OKAY;
}

/** if not already existing, adds row to global cut pool */
SCIP_RETCODE SCIPaddPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolAddRow(scip->cutpool, scip->mem->probmem, scip->set, row) );

   return SCIP_OKAY;
}

/** removes the row from the global cut pool */
SCIP_RETCODE SCIPdelPoolCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdelPoolCut", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolDelRow(scip->cutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** gets current cuts in the global cut pool */
SCIP_CUT** SCIPgetPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPoolCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPcutpoolGetCuts(scip->cutpool);
}

/** gets current number of rows in the global cut pool */
int SCIPgetNPoolCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPoolCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPcutpoolGetNCuts(scip->cutpool);
}

/** gets the global cut pool used by SCIP */
SCIP_CUTPOOL* SCIPgetGlobalCutpool(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetGlobalCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE) );

   return scip->cutpool;
}

/** creates a cut pool */
SCIP_RETCODE SCIPcreateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool,            /**< pointer to store cut pool */
   int                   agelimit            /**< maximum age a cut can reach before it is deleted from the pool */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateCutpool", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolCreate(cutpool, scip->mem->probmem, scip->set, agelimit, FALSE) );

   return SCIP_OKAY;
}

/** frees a cut pool */
SCIP_RETCODE SCIPfreeCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL**        cutpool             /**< pointer to store cut pool */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeCutpool", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPcutpoolFree(cutpool, scip->mem->probmem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** if not already existing, adds row to a cut pool and captures it */
SCIP_RETCODE SCIPaddRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolAddRow(cutpool, scip->mem->probmem, scip->set, row) );

   return SCIP_OKAY;
}

/** adds row to a cut pool and captures it; doesn't check for multiple cuts */
SCIP_RETCODE SCIPaddNewRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< cutting plane to add */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddNewRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolAddNewRow(cutpool, scip->mem->probmem, scip->set, row) );

   return SCIP_OKAY;
}

/** removes the LP row from a cut pool */
SCIP_RETCODE SCIPdelRowCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_ROW*             row                 /**< row to remove */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdelRowCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   SCIP_CALL( SCIPcutpoolDelRow(cutpool, scip->mem->probmem, scip->set, scip->stat, scip->lp, row) );

   return SCIP_OKAY;
}

/** separates cuts from a cut pool */
SCIP_RETCODE SCIPseparateCutpool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CUTPOOL*         cutpool,            /**< cut pool */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPseparateCutpool", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   assert(SCIPtreeGetCurrentNode(scip->tree) != NULL);

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("cannot add cuts, because node LP is not processed\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPcutpoolSeparate(cutpool, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter,
         scip->lp, scip->sepastore, (SCIPtreeGetCurrentDepth(scip->tree) == 0), result) );
   return SCIP_OKAY;
}

/** separates the given primal solution or the current LP solution by calling the separators and constraint handlers'
 *  separation methods;
 *  the generated cuts are stored in the separation storage and can be accessed with the methods SCIPgetCuts() and
 *  SCIPgetNCuts();
 *  after evaluating the cuts, you have to call SCIPclearCuts() in order to remove the cuts from the
 *  separation storage;
 *  it is possible to call SCIPseparateSol() multiple times with different solutions and evaluate the found cuts
 *  afterwards
 */
SCIP_RETCODE SCIPseparateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution that should be separated, or NULL for LP solution */
   SCIP_Bool             pretendroot,        /**< should the cut separators be called as if we are at the root node? */
   SCIP_Bool             onlydelayed,        /**< should only separators be called that were delayed in the previous round? */
   SCIP_Bool*            delayed,            /**< pointer to store whether a separator was delayed */
   SCIP_Bool*            cutoff              /**< pointer to store whether the node can be cut off */
   )
{
   int actdepth;

   SCIP_CALL( checkStage(scip, "SCIPseparateSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* get current depth */
   actdepth = (pretendroot ? 0 : SCIPtreeGetCurrentDepth(scip->tree));

   /* apply separation round */
   SCIP_CALL( SCIPseparationRound(scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter, scip->transprob, scip->primal, scip->tree, scip->lp, scip->sepastore,
         sol, actdepth, onlydelayed, delayed, cutoff) );

   return SCIP_OKAY;
}

/** gets the array of cuts currently stored in the separation storage */
SCIP_ROW** SCIPgetCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetCuts", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetCuts(scip->sepastore);
}

/** get current number of cuts in the separation storage */
int SCIPgetNCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNCuts", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCuts(scip->sepastore);
}

/** clears the separation storage */
SCIP_RETCODE SCIPclearCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPclearCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsepastoreClearCuts(scip->sepastore, scip->mem->probmem, scip->set, scip->eventqueue, scip->eventfilter, scip->lp) );

   return SCIP_OKAY;
}

/** removes cuts that are inefficacious w.r.t. the current LP solution from separation storage without adding the cuts to the LP */
SCIP_RETCODE SCIPremoveInefficaciousCuts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Bool isroot;

   SCIP_CALL( checkStage(scip, "SCIPremoveInefficaciousCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   isroot = FALSE;
   if( SCIPtreeGetCurrentDepth(scip->tree) == 0 )
      isroot = TRUE;
   SCIP_CALL( SCIPsepastoreRemoveInefficaciousCuts(scip->sepastore, scip->mem->probmem, scip->set, scip->stat, 
         scip->eventqueue, scip->eventfilter, scip->lp, isroot) );

   return SCIP_OKAY;
}

#if 0
/**@todo make this method available; implement methods SCIPaddProbingRow() and SCIPaddProbingCol();
 *       -> the probing node needs a counter (like the fork) on the number of added rows and cols,
 *          and treeBacktrackProbing() needs to shrink the LP;
 *       this is useful, e.g., for the feaspump heuristic, s.t. it can temporarily add the auxiliary columns
 *       needed to represent the objective function for integer variables not on one of their bounds
 */
/** applies the cuts in the separation storage to the LP and clears the storage afterwards;
 *  this method can only be applied during probing; the user should resolve the probing LP afterwards
 *  in order to get a new solution
 */
SCIP_RETCODE SCIPapplyCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether an empty domain was created */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPapplyCuts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsepastoreApplyCuts(scip->sepastore, scip->mem->probmem, scip->set, scip->stat,
         scip->eventqueue, scip->eventfilter, scip->tree, scip->lp, scip->branchcand, scip->eventqueue, cutoff) );

   return SCIP_OKAY;
}
#endif




/*
 * LP diving methods
 */

/** initiates LP diving, making methods SCIPchgVarObjDive(), SCIPchgVarLbDive(), and SCIPchgVarUbDive() available */
SCIP_RETCODE SCIPstartDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstartDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   assert(SCIPnodeGetType(SCIPgetCurrentNode(scip)) == SCIP_NODETYPE_FOCUSNODE);

   if( SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("already in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("cannot start diving while being in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   if( !SCIPtreeIsFocusNodeLPConstructed(scip->tree) )
   {
      SCIPerrorMessage("cannot start diving if LP has not been constructed\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPtreeHasCurrentNodeLP(scip->tree));

   SCIP_CALL( SCIPlpStartDive(scip->lp, scip->mem->probmem, scip->set, scip->stat) );

   return SCIP_OKAY;
}

/** quits LP diving and resets bounds and objective values of columns to the current node's values */
SCIP_RETCODE SCIPendDive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPendDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* unmark the diving flag in the LP and reset all variables' objective and bound values */
   SCIP_CALL( SCIPlpEndDive(scip->lp, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter,
         scip->transprob, scip->transprob->vars, scip->transprob->nvars) );

   /* the lower bound may have changed slightly due to LP resolve in SCIPlpEndDive() */
   if( !scip->lp->resolvelperror && scip->tree->focusnode != NULL && SCIPlpIsRelax(scip->lp) )
   {
      assert(SCIPtreeIsFocusNodeLPConstructed(scip->tree));
      SCIP_CALL( SCIPnodeUpdateLowerboundLP(scip->tree->focusnode, scip->set, scip->stat, scip->lp) );
   }
   /* reset the probably changed LP's cutoff bound */
   SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, scip->primal->cutoffbound) );
   assert(scip->lp->cutoffbound == scip->primal->cutoffbound); /*lint !e777*/

   /* if a new best solution was created, the cutoff of the tree was delayed due to diving;
    * the cutoff has to be done now.
    */
   if( scip->tree->cutoffdelayed )
   {
      SCIP_CALL( SCIPtreeCutoff(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue,
            scip->lp, scip->primal->cutoffbound) );
   }

   return SCIP_OKAY;
}

/** changes variable's objective value in current dive */
SCIP_RETCODE SCIPchgVarObjDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Real             newobj              /**< new objective value */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* invalidate the LP's cutoff bound, since this has nothing to do with the current objective value anymore;
    * the cutoff bound is reset in SCIPendDive()
    */
   SCIP_CALL( SCIPlpSetCutoffbound(scip->lp, scip->set, SCIPsetInfinity(scip->set)) );

   /* mark the LP's objective function invalid */
   SCIPlpMarkDivingObjChanged(scip->lp);

   /* change the objective value of the variable in the diving LP */
   SCIP_CALL( SCIPvarChgObjDive(var, scip->set, scip->lp, newobj) );

   return SCIP_OKAY;
}

/** changes variable's lower bound in current dive */
SCIP_RETCODE SCIPchgVarLbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPvarChgLbDive(var, scip->set, scip->lp, newbound) );

   return SCIP_OKAY;
}

/** changes variable's upper bound in current dive */
SCIP_RETCODE SCIPchgVarUbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPvarChgUbDive(var, scip->set, scip->lp, newbound) );

   return SCIP_OKAY;
}

/** gets variable's objective value in current dive */
SCIP_Real SCIPgetVarObjDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarObjDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      SCIPABORT();
   }

   return SCIPvarGetObjLP(var);
}

/** gets variable's lower bound in current dive */
SCIP_Real SCIPgetVarLbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarLbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      SCIPABORT();
   }

   return SCIPvarGetLbLP(var, scip->set);
}

/** gets variable's upper bound in current dive */
SCIP_Real SCIPgetVarUbDive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get the bound for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetVarUbDive", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      SCIPABORT();
   }

   return SCIPvarGetUbLP(var, scip->set);
}

/** solves the LP of the current dive; no separation or pricing is applied */
SCIP_RETCODE SCIPsolveDiveLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsolveDiveLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("not in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   /* solve diving LP */
   SCIP_CALL( SCIPlpSolveAndEval(scip->lp, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter, scip->transprob,
         itlim, FALSE, FALSE, FALSE, lperror) );

   /* analyze an infeasible LP (not necessary in the root node)
    * the infeasibility in diving is only proven, if all columns are in the LP (and no external pricers exist)
    */
   if( !scip->set->misc_exactsolve && SCIPtreeGetCurrentDepth(scip->tree) > 0
      && (SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_INFEASIBLE
         || (SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_OBJLIMIT && !SCIPlpDivingObjChanged(scip->lp)))
      && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) )
   {
      SCIP_CALL( SCIPconflictAnalyzeLP(scip->conflict, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->tree, scip->lp, NULL) );
   }

   return SCIP_OKAY;
}

/** returns the number of the node in the current branch and bound run, where the last LP was solved in diving
 *  or probing mode
 */
SCIP_Longint SCIPgetLastDivenode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLastDivenode", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->lastdivenode;
}




/*
 * probing methods
 */

/** returns whether we are in probing mode */
SCIP_Bool SCIPinProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPinProbing", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   
   return SCIPtreeProbing(scip->tree);
}

/** initiates probing, making methods SCIPnewProbingNode(), SCIPbacktrackProbing(), SCIPchgVarLbProbing(),
 *  SCIPchgVarUbProbing(), SCIPfixVarProbing(), SCIPpropagateProbing(), and SCIPsolveProbingLP() available
 */
SCIP_RETCODE SCIPstartProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstartProbing", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("already in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   if( scip->lp != NULL && SCIPlpDiving(scip->lp) )
   {
      SCIPerrorMessage("cannot start probing while in diving mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPtreeStartProbing(scip->tree, scip->mem->probmem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** creates a new probing sub node, whose changes can be undone by backtracking to a higher node in the probing path
 *  with a call to SCIPbacktrackProbing();
 *  using a sub node for each set of probing bound changes can improve conflict analysis
 */
SCIP_RETCODE SCIPnewProbingNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPnewProbingNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPtreeCreateProbingNode(scip->tree, scip->mem->probmem, scip->set, scip->lp) );

   return SCIP_OKAY;
}

/** returns the current probing depth, i.e. the number of probing sub nodes existing in the probing path */
int SCIPgetProbingDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetProbingDepth", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      SCIPABORT();
   }

   return SCIPtreeGetProbingDepth(scip->tree);
}

/** undoes all changes to the problem applied in probing up to the given probing depth;
 *  the changes of the probing node of the given probing depth are the last ones that remain active;
 *  changes that were applied before calling SCIPnewProbingNode() cannot be undone
 */
SCIP_RETCODE SCIPbacktrackProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   probingdepth        /**< probing depth of the node in the probing path that should be reactivated */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPbacktrackProbing", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   if( probingdepth < 0 || probingdepth > SCIPtreeGetProbingDepth(scip->tree) )
   {
      SCIPerrorMessage("backtracking probing depth %d out of current probing range [0,%d]\n",
         probingdepth, SCIPtreeGetProbingDepth(scip->tree));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBacktrackProbing(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->lp,
         scip->branchcand, scip->eventqueue, scip->eventfilter, probingdepth) );

   return SCIP_OKAY;
}

/** quits probing and resets bounds and constraints to the focus node's environment */
SCIP_RETCODE SCIPendProbing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPendProbing", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* switch back from probing to normal operation mode and restore variables and constraints to focus node */
   SCIP_CALL( SCIPtreeEndProbing(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->lp,
         scip->branchcand, scip->eventqueue, scip->eventfilter) );

   return SCIP_OKAY;
}

/** injects a change of variable's lower bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarLb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 */
SCIP_RETCODE SCIPchgVarLbProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarLbProbing", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(scip->tree)) == SCIP_NODETYPE_PROBINGNODE);

   SCIPvarAdjustLb(var, scip->set, &newbound);

   SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
         scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_LOWER, TRUE) );

   return SCIP_OKAY;
}

/** injects a change of variable's upper bound into current probing node; the same can also be achieved with a call to
 *  SCIPchgVarUb(), but in this case, the bound change would be treated like a deduction instead of a branching decision
 */
SCIP_RETCODE SCIPchgVarUbProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             newbound            /**< new value for bound */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgVarUbProbing", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(scip->tree)) == SCIP_NODETYPE_PROBINGNODE);

   SCIPvarAdjustUb(var, scip->set, &newbound);

   SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
         scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, newbound, SCIP_BOUNDTYPE_UPPER, TRUE) );

   return SCIP_OKAY;
}

/** injects a change of variable's bounds into current probing node to fix the variable to the specified value;
 *  the same can also be achieved with a call to SCIPfixVar(), but in this case, the bound changes would be treated
 *  like deductions instead of branching decisions
 */
SCIP_RETCODE SCIPfixVarProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the bound for */
   SCIP_Real             fixedval            /**< value to fix variable to */
   )
{
   SCIP_Real fixlb;
   SCIP_Real fixub;

   SCIP_CALL( checkStage(scip, "SCIPfixVarProbing", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }
   assert(SCIPnodeGetType(SCIPtreeGetCurrentNode(scip->tree)) == SCIP_NODETYPE_PROBINGNODE);

   /** we adjust the fixing value here and compare the old bound with the adjusted values because otherwise, 
    *  it might happen that the unadjusted value is better and we add the boundchange,
    *  but within SCIPnodeAddBoundchg() the bounds are adjusted - using the feasibility epsilon for integer variables -
    *  and it is asserted, that the bound is still better than the old one which might then be incorrect.
    */
   fixlb = fixedval;
   fixub = fixedval;
   SCIPvarAdjustLb(var, scip->set, &fixlb);
   SCIPvarAdjustUb(var, scip->set, &fixub);
   assert(SCIPsetIsEQ(scip->set, fixlb, fixub));

   if( SCIPsetIsGT(scip->set, fixlb, SCIPvarGetLbLocal(var)) )
   {
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, fixlb, SCIP_BOUNDTYPE_LOWER, TRUE) );
   }
   if( SCIPsetIsLT(scip->set, fixub, SCIPvarGetUbLocal(var)) )
   {
      SCIP_CALL( SCIPnodeAddBoundchg(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
            scip->tree, scip->lp, scip->branchcand, scip->eventqueue, var, fixub, SCIP_BOUNDTYPE_UPPER, TRUE) );
   }

   return SCIP_OKAY;
}

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 */
SCIP_RETCODE SCIPpropagateProbing(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxproprounds,      /**< maximal number of propagation rounds (-1: no limit, 0: parameter settings) */
   SCIP_Bool*            cutoff,             /**< pointer to store whether the probing node can be cut off */
   SCIP_Longint*         ndomredsfound       /**< pointer to store the number of domain reductions found, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPpropagateProbing", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   if( ndomredsfound != NULL )
      *ndomredsfound = -(scip->stat->nprobboundchgs + scip->stat->nprobholechgs);

   SCIP_CALL( SCIPpropagateDomains(scip->mem->probmem, scip->set, scip->stat, scip->transprob, 
         scip->primal, scip->tree, scip->conflict, 0, maxproprounds, SCIP_PROPTIMING_ALWAYS, cutoff) );

   if( ndomredsfound != NULL )
      *ndomredsfound += scip->stat->nprobboundchgs + scip->stat->nprobholechgs;

   return SCIP_OKAY;
}

/** applies domain propagation on the probing sub problem, that was changed after SCIPstartProbing() was called;
 *  only propagations of the binary variables fixed at the current probing node that are triggered by the implication
 *  graph and the clique table are applied;
 *  the propagated domains of the variables can be accessed with the usual bound accessing calls SCIPvarGetLbLocal()
 *  and SCIPvarGetUbLocal(); the propagation is only valid locally, i.e. the local bounds as well as the changed
 *  bounds due to SCIPchgVarLbProbing(), SCIPchgVarUbProbing(), and SCIPfixVarProbing() are used for propagation
 */
SCIP_RETCODE SCIPpropagateProbingImplications(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool*            cutoff              /**< pointer to store whether the probing node can be cut off */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPpropagateProbingImplications", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPnodePropagateImplics(SCIPtreeGetCurrentNode(scip->tree), scip->mem->probmem, scip->set, scip->stat,
         scip->tree, scip->lp, scip->branchcand, scip->eventqueue, cutoff) );

   return SCIP_OKAY;
}

/** solves the LP at the current probing node (cannot be applied at preprocessing stage) with or without pricing */
static
SCIP_RETCODE solveProbingLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool             pricing,            /**< should pricing be applied? */
   SCIP_Bool             pretendroot,        /**< should the pricers be called as if we are at the root node? */
   SCIP_Bool             displayinfo,        /**< should info lines be displayed after each pricing round? */
   int                   maxpricerounds,     /**< maximal number of pricing rounds (-1: no limit) */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   assert(lperror != NULL);
   assert(SCIPtreeIsFocusNodeLPConstructed(scip->tree));

   if( !SCIPtreeProbing(scip->tree) )
   {
      SCIPerrorMessage("not in probing mode\n");
      return SCIP_INVALIDCALL;
   }

   /* load the LP state (if necessary) */
   SCIP_CALL( SCIPtreeLoadProbingLPState(scip->tree, scip->mem->probmem, scip->set, scip->eventqueue, scip->lp) );

   /* solve probing LP */
   SCIP_CALL( SCIPlpSolveAndEval(scip->lp, scip->mem->probmem, scip->set, scip->stat, scip->eventqueue, scip->eventfilter, scip->transprob,
         itlim, FALSE, FALSE, FALSE, lperror) );

   /* mark the probing node to have a solved LP */
   if( !(*lperror) )
   {
      SCIP_CALL( SCIPtreeMarkProbingNodeHasLP(scip->tree, scip->mem->probmem, scip->lp) );

      /* call pricing */
      if( pricing )
      {
         SCIP_Bool mustsepa;
         int npricedcolvars;
         SCIP_Real lowerbound;
         SCIP_Bool result;

         mustsepa = FALSE;
         SCIP_CALL( SCIPpriceLoop(scip->mem->probmem, scip->set, scip->stat, scip->transprob, scip->primal, scip->tree, scip->lp,
               scip->pricestore, scip->sepastore, scip->branchcand, scip->eventqueue, scip->eventfilter, pretendroot, displayinfo, maxpricerounds,
               &npricedcolvars, &mustsepa, &lowerbound, lperror, &result) );

         /* mark the probing node again to update the LP size in the node and the tree path */
         if( !(*lperror) )
         {
            SCIP_CALL( SCIPtreeMarkProbingNodeHasLP(scip->tree, scip->mem->probmem, scip->lp) );
         }
      }
   }

   /* remember that probing might have changed the LPi state; this holds even if solving returned with an LP error */
   scip->tree->probingsolvedlp = TRUE;

   /* analyze an infeasible LP (not necessary in the root node)
    * the infeasibility in probing is only proven, if all columns are in the LP (and no external pricers exist)
    */
   if( !(*lperror) && !scip->set->misc_exactsolve && SCIPtreeGetCurrentDepth(scip->tree) > 0 && SCIPlpIsRelax(scip->lp)
      && (SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_INFEASIBLE
         || SCIPlpGetSolstat(scip->lp) == SCIP_LPSOLSTAT_OBJLIMIT)
      && SCIPprobAllColsInLP(scip->transprob, scip->set, scip->lp) )
   {
      SCIP_CALL( SCIPconflictAnalyzeLP(scip->conflict, scip->mem->probmem, scip->set, scip->stat, scip->transprob,
            scip->tree, scip->lp, NULL) );
   }

   return SCIP_OKAY;
}

/** solves the LP at the current probing node (cannot be applied at preprocessing stage);
 *  no separation or pricing is applied
 */
SCIP_RETCODE SCIPsolveProbingLP(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   itlim,              /**< maximal number of LP iterations to perform, or -1 for no limit */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsolveProbingLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( solveProbingLP(scip, itlim, FALSE, FALSE, FALSE, -1, lperror) );

   return SCIP_OKAY;
}

/** solves the LP at the current probing node (cannot be applied at preprocessing stage) and applies pricing
 *  until the LP is solved to optimality; no separation is applied
 */
SCIP_RETCODE SCIPsolveProbingLPWithPricing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             pretendroot,        /**< should the pricers be called as if we are at the root node? */
   SCIP_Bool             displayinfo,        /**< should info lines be displayed after each pricing round? */
   int                   maxpricerounds,     /**< maximal number of pricing rounds (-1: no limit);
                                              *   a finite limit means that the LP might not be solved to optimality! */
   SCIP_Bool*            lperror             /**< pointer to store whether an unresolved LP error occurred */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsolveProbingLPWithPricing", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( solveProbingLP(scip, -1, TRUE, pretendroot, displayinfo, maxpricerounds, lperror) );

   return SCIP_OKAY;
}




/*
 * branching methods
 */

/** gets branching candidates for LP solution branching (fractional variables) along with solution values,
 *  fractionalities, and number of branching candidates;
 *  branching rules should always select the branching candidate among the first npriolpcands of the candidate
 *  list
 */
SCIP_RETCODE SCIPgetLPBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   SCIP_Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   SCIP_Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*                  nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*                  npriolpcands        /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPerrorMessage("LP not solved to optimality - solstat=%d\n", SCIPlpGetSolstat(scip->lp));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
         lpcands, lpcandssol, lpcandsfrac, nlpcands, npriolpcands) );

   return SCIP_OKAY;
}

/** gets number of branching candidates for LP solution branching (number of fractional variables) */
int SCIPgetNLPBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int nlpcands;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPerrorMessage("LP not solved to optimality\n");
      SCIPABORT();
   }

   SCIP_CALL_ABORT( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
         NULL, NULL, NULL, &nlpcands, NULL) );

   return nlpcands;
}

/** gets number of branching candidates with maximal priority for LP solution branching */
int SCIPgetNPrioLPBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int npriolpcands;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioLPBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_OPTIMAL && SCIPlpGetSolstat(scip->lp) != SCIP_LPSOLSTAT_UNBOUNDEDRAY )
   {
      SCIPerrorMessage("LP not solved to optimality\n");
      SCIPABORT();
   }

   SCIP_CALL_ABORT( SCIPbranchcandGetLPCands(scip->branchcand, scip->set, scip->stat, scip->lp,
         NULL, NULL, NULL, NULL, &npriolpcands) );

   return npriolpcands;
}

/** gets external branching candidates along with solution values, scores, and number of branching candidates;
 *  these branching candidates can be used by relaxations or nonlinear constraint handlers
 *  branching rules should always select the branching candidate among the first nprioexterncands of the candidate
 *  list
 */
SCIP_RETCODE SCIPgetExternBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           externcands,        /**< pointer to store the array of extern branching candidates, or NULL */
   SCIP_Real**           externcandssol,     /**< pointer to store the array of extern candidate solution values, or NULL */
   SCIP_Real**           externcandsscore,   /**< pointer to store the array of extern candidate scores, or NULL */
   int*                  nexterncands,       /**< pointer to store the number of extern branching candidates, or NULL */
   int*                  nprioexterncands,   /**< pointer to store the number of candidates with maximal priority, or NULL */
   int*                  nprioexternbins,    /**< pointer to store the number of binary candidates with maximal priority, or NULL */
   int*                  nprioexternints,    /**< pointer to store the number of integer candidates with maximal priority, or NULL */
   int*                  nprioexternimpls    /**< pointer to store the number of implicit integer candidates with maximal priority, 
                                              *   or NULL */
   )
{
   assert(scip != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchcandGetExternCands(scip->branchcand, externcands, externcandssol, externcandsscore, nexterncands, 
         nprioexterncands, nprioexternbins, nprioexternints, nprioexternimpls) );

   return SCIP_OKAY;
}

/** gets number of external branching candidates */
int SCIPgetNExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNExternCands(scip->branchcand);
}

/** gets number of external branching candidates with maximal branch priority */
int SCIPgetNPrioExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternCands(scip->branchcand);
}

/** gets number of binary external branching candidates with maximal branch priority */
int SCIPgetNPrioExternBranchBins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioExternBranchBins", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternBins(scip->branchcand);
}


/** gets number of integer external branching candidates with maximal branch priority */
int SCIPgetNPrioExternBranchInts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioExternBranchInts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternInts(scip->branchcand);
}

/** gets number of implicit integer external branching candidates with maximal branch priority */
int SCIPgetNPrioExternBranchImpls(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioExternBranchImpls", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternImpls(scip->branchcand);
}

/** gets number of continuous external branching candidates with maximal branch priority */
int SCIPgetNPrioExternBranchConts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioExternBranchConts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioExternConts(scip->branchcand);
}

/** insert variable, its score and its solution value into the external branching candidate storage
 * the relative difference of the current lower and upper bounds of a continuous variable must be at least epsilon
 */
SCIP_RETCODE SCIPaddExternBranchCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to insert */
   SCIP_Real             score,              /**< score of external candidate, e.g. infeasibility */
   SCIP_Real             solval              /**< value of the variable in the current solution */
   )
{
   assert(scip != NULL);

   SCIP_CALL( checkStage(scip, "SCIPaddExternBranchCand", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchcandAddExternCand(scip->branchcand, scip->set, var, score, solval) );
   
   return SCIP_OKAY;
}

/** removes all external candidates from the storage for external branching */
void SCIPclearExternBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPclearExternBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPbranchcandClearExternCands(scip->branchcand);
}

/** checks whether the given variable is contained in the candidate storage for external branching */
SCIP_Bool SCIPcontainsExternBranchCand(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to look for */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPcontainsExternBranchCand", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandContainsExternCand(scip->branchcand, var);
}



/** gets branching candidates for pseudo solution branching (non-fixed variables) along with the number of candidates */
SCIP_RETCODE SCIPgetPseudoBranchCands(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*                  npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*                  npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchcandGetPseudoCands(scip->branchcand, scip->set, scip->transprob,
         pseudocands, npseudocands, npriopseudocands) );

   return SCIP_OKAY;
}

/** gets branching candidates for pseudo solution branching (non-fixed variables) */
int SCIPgetNPseudoBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPseudoCands(scip->branchcand);
}

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchCands(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchCands", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoCands(scip->branchcand);
}

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchBins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchBins", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoBins(scip->branchcand);
}

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchInts(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchInts", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoInts(scip->branchcand);
}

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPgetNPrioPseudoBranchImpls(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrioPseudoBranchImpls", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchcandGetNPrioPseudoImpls(scip->branchcand);
}

/** calculates the branching score out of the gain predictions for a binary branching */
SCIP_Real SCIPgetBranchScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   SCIP_Real             downgain,           /**< prediction of objective gain for rounding downwards */
   SCIP_Real             upgain              /**< prediction of objective gain for rounding upwards */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBranchScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchGetScore(scip->set, var, downgain, upgain);
}

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
SCIP_Real SCIPgetBranchScoreMultiple(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int                   nchildren,          /**< number of children that the branching will create */
   SCIP_Real*            gains               /**< prediction of objective gain for each child */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBranchScoreMultiple", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchGetScoreMultiple(scip->set, var, nchildren, gains);
}

/** computes a branching point for a continuous or discrete variable
 * @see SCIPbranchGetBranchingPoint
 */
SCIP_Real SCIPgetBranchingPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable, of which the branching point should be computed */
   SCIP_Real             suggestion          /**< suggestion for branching point, or SCIP_INVALID if no suggestion */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBranchingPoint", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPbranchGetBranchingPoint(scip->set, scip->tree, var, suggestion);
}

/** calculates the node selection priority for moving the given variable's LP value to the given target value;
 *  this node selection priority can be given to the SCIPcreateChild() call
 */
SCIP_Real SCIPcalcNodeselPriority(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the branching is applied */
   SCIP_BRANCHDIR        branchdir,          /**< type of branching that was performed: upwards, downwards, or fixed;
                                              *   fixed should only be used, when both bounds changed
                                              */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPcalcNodeselPriority", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeCalcNodeselPriority(scip->tree, scip->set, scip->stat, var, branchdir, targetvalue);
}

/** calculates an estimate for the objective of the best feasible solution contained in the subtree after applying the given
 *  branching; this estimate can be given to the SCIPcreateChild() call
 */
SCIP_Real SCIPcalcChildEstimate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable on which the branching is applied */
   SCIP_Real             targetvalue         /**< new value of the variable in the child node */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPcalcChildEstimate", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeCalcChildEstimate(scip->tree, scip->set, scip->stat, var, targetvalue);
}

/** creates a child node of the focus node */
SCIP_RETCODE SCIPcreateChild(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE**           node,               /**< pointer to node data structure */
   SCIP_Real             nodeselprio,        /**< node selection priority of new node */
   SCIP_Real             estimate            /**< estimate for(transformed) objective value of best feasible solution in subtree */
   )
{
   assert(node != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcreateChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPnodeCreateChild(node, scip->mem->probmem, scip->set, scip->stat, scip->tree, nodeselprio, estimate) );

   return SCIP_OKAY;
}

/** branches on a non-continuous variable v using the current LP or pseudo solution;
 *  if solution value x' is fractional, two child nodes will be created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if solution value is integral, the x' is equal to lower or upper bound of the branching
 *  variable and the bounds of v are finite, then two child nodes will be created
 *  (x <= x", x >= x"+1 with x" = floor((lb + ub)/2)),
 *  otherwise (up to) three child nodes will be created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
SCIP_RETCODE SCIPbranchVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPbranchVar", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
   {
      SCIPerrorMessage("cannot branch on continuous variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   if( SCIPsetIsEQ(scip->set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      SCIPerrorMessage("cannot branch on variable <%s> with fixed domain [%.15g,%.15g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBranchVar(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->branchcand,
         scip->eventqueue, var, SCIP_INVALID, downchild, eqchild, upchild) );

   return SCIP_OKAY;
}

/** branches on a variable x using a given value x'; 
 *  for continuous variables with relative domain width larger epsilon, x' must not be one of the bounds;
 *  two child nodes (x <= x', x >= x') are created;
 *  for integer variables, if solution value x' is fractional, two child nodes are created
 *  (x <= floor(x'), x >= ceil(x')),
 *  if x' is integral, three child nodes are created
 *  (x <= x'-1, x == x', x >= x'+1)
 */
SCIP_RETCODE SCIPbranchVarVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on */
   SCIP_NODE**           downchild,          /**< pointer to return the left child with variable rounded down, or NULL */
   SCIP_NODE**           eqchild,            /**< pointer to return the middle child with variable fixed, or NULL */
   SCIP_NODE**           upchild             /**< pointer to return the right child with variable rounded up, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPbranchVarVal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
   
   /* for a continuous variable, their will either be variable fixing or a branching
    * fixing is done if RelEQ(lb,ub)
    * in the other case, the given branching value should be such that it does not sits on one of the bounds
    * we assert this by requiring that it is at least eps/2 away from each bound
    * the /2 is there, because ub-lb may be in (eps, 2eps], in which case there is no way to choose a branching value that is at least eps away from both bounds
    * however, if variable bounds are below/above -/+infinity/2.1, then SCIPisLT will give an assert, so we omit the check then
    */
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS ||
      SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) ||
      SCIPisInfinity(scip, -2.1*SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, 2.1*SCIPvarGetUbLocal(var)) ||
      (SCIPisLT(scip, 2.1*SCIPvarGetLbLocal(var), 2.1*val) && SCIPisLT(scip, 2.1*val, 2.1*SCIPvarGetUbLocal(var)) ) );

   if( SCIPsetIsEQ(scip->set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      SCIPerrorMessage("cannot branch on variable <%s> with fixed domain [%.15g,%.15g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBranchVar(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->branchcand,
         scip->eventqueue, var, val, downchild, eqchild, upchild) );

   return SCIP_OKAY;
}

/** n-ary branching on a variable x using a given value
 * Branches on variable x such that up to n/2 children are created on each side of the usual branching value.
 * The branching value is selected as in SCIPbranchVarVal().
 * The parameters minwidth and widthfactor determine the domain width of the branching variable in the child nodes.
 * If n is odd, one child with domain width 'width' and having the branching value in the middle is created.
 * Otherwise, two children with domain width 'width' and being left and right of the branching value are created.
 * Next further nodes to the left and right are created, where width is multiplied by widthfactor with increasing distance from the first nodes.
 * The initial width is calculated such that n/2 nodes are created to the left and to the right of the branching value.
 * If this value is below minwidth, the initial width is set to minwidth, which may result in creating less than n nodes.
 *
 * Giving a large value for widthfactor results in creating children with small domain when close to the branching value
 * and large domain when closer to the current variable bounds. That is, setting widthfactor to a very large value and n to 3
 * results in a ternary branching where the branching variable is mostly fixed in the middle child.
 * Setting widthfactor to 1.0 results in children where the branching variable always has the same domain width
 * (except for one child if the branching value is not in the middle).
 */
SCIP_RETCODE SCIPbranchVarValNary(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to branch on */
   SCIP_Real             val,                /**< value to branch on */
   int                   n,                  /**< attempted number of children to be created, must be >= 2 */
   SCIP_Real             minwidth,           /**< minimal domain width in children */
   SCIP_Real             widthfactor,        /**< multiplier for children domain width with increasing distance from val, must be >= 1.0 */
   int*                  nchildren           /**< buffer to store number of created children, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPbranchVarValNary", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* see comment in SCIPbranchVarVal */
   assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS ||
      SCIPisRelEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) ||
      SCIPisInfinity(scip, -2.1*SCIPvarGetLbLocal(var)) || SCIPisInfinity(scip, 2.1*SCIPvarGetUbLocal(var)) ||
      (SCIPisLT(scip, 2.1*SCIPvarGetLbLocal(var), 2.1*val) && SCIPisLT(scip, 2.1*val, 2.1*SCIPvarGetUbLocal(var)) ) );

   if( SCIPsetIsEQ(scip->set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      SCIPerrorMessage("cannot branch on variable <%s> with fixed domain [%.15g,%.15g]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPtreeBranchVarNary(scip->tree, scip->mem->probmem, scip->set, scip->stat, scip->lp, scip->branchcand,
         scip->eventqueue, var, val, n, minwidth, widthfactor, nchildren) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
SCIP_RETCODE SCIPbranchLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPbranchLP", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecLP(scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
         scip->sepastore, scip->branchcand, scip->eventqueue, scip->primal->cutoffbound, TRUE, result) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on a external candidates; if no such candidates exist, the result is SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPbranchExtern(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPbranchExtern", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecExtern(scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
         scip->sepastore, scip->branchcand, scip->eventqueue, scip->primal->cutoffbound, TRUE, result) );

   return SCIP_OKAY;
}

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
SCIP_RETCODE SCIPbranchPseudo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPbranchPseudo", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPbranchExecPseudo(scip->mem->probmem, scip->set, scip->stat, scip->tree, scip->lp,
         scip->branchcand, scip->eventqueue, scip->primal->cutoffbound, TRUE, result) );

   return SCIP_OKAY;
}




/*
 * primal solutions
 */

/** creates a primal solution, initialized to zero */
SCIP_RETCODE SCIPcreateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );
   
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPsolCreateOriginal(sol, scip->mem->probmem, scip->set, scip->stat, scip->origprimal, NULL, heur) );
      return SCIP_OKAY;
      
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPsolCreate(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, heur) );
      return SCIP_OKAY;
      
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** creates a primal solution, initialized to the current LP solution */
SCIP_RETCODE SCIPcreateLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPtreeHasCurrentNodeLP(scip->tree) )
   {
      SCIPerrorMessage("LP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolCreateLPSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->lp,
         heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current NLP solution */
SCIP_RETCODE SCIPcreateNLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateNLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPisNLPConstructed(scip) )
   {
      SCIPerrorMessage("NLP does not exist\n");
      return SCIP_INVALIDCALL;
   }
   assert(scip->nlp != NULL);

   if( !SCIPnlpHasSolution(scip->nlp) )
   {
      SCIPerrorMessage("NLP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolCreateNLPSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->nlp,
         heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current relaxation solution */
SCIP_RETCODE SCIPcreateRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateRelaxSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("relaxation solution is not valid\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolCreateRelaxSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->relaxation, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current pseudo solution */
SCIP_RETCODE SCIPcreatePseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreatePseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolCreatePseudoSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->lp,
         heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to the current LP or pseudo solution, depending on whether the LP was solved
 *  at the current node
 */
SCIP_RETCODE SCIPcreateCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolCreateCurrentSol(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, scip->lp,
         heur) );

   return SCIP_OKAY;
}

/** creates a primal solution, initialized to unknown values */
SCIP_RETCODE SCIPcreateUnknownSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateUnknownSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolCreateUnknown(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, heur) );

   return SCIP_OKAY;
}

/** creates a primal solution living in the original problem space, initialized to zero;
 *  a solution in original space allows to set original variables to values, that would be invalid in the
 *  transformed problem due to preprocessing fixings or aggregations
 */
SCIP_RETCODE SCIPcreateOrigSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateOrigSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPsolCreateOriginal(sol, scip->mem->probmem, scip->set, scip->stat, scip->origprimal, NULL, heur) );
      return SCIP_OKAY;
      
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPsolCreateOriginal(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree, heur) );
      return SCIP_OKAY;
      
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** creates a copy of a primal solution; note that a copy of a linked solution is also linked and needs to be unlinked
 *  if it should stay unaffected from changes in the LP or pseudo solution
 */
SCIP_RETCODE SCIPcreateSolCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to store the solution */
   SCIP_SOL*             sourcesol           /**< primal CIP solution to copy */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateSolCopy", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* check if we want to copy the current solution, which is the same as creating a current solution */
   if( sourcesol == NULL )
   {
      SCIP_CALL( SCIPcreateCurrentSol(scip, sol, NULL) );
   }
   else
   {
      SCIP_CALL( SCIPsolCopy(sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, sourcesol) );
   }

   return SCIP_OKAY;
}

/** frees primal CIP solution */
SCIP_RETCODE SCIPfreeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol                 /**< pointer to the solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPsolFree(sol, scip->mem->probmem, scip->origprimal) );
      break;
   case SCIP_STAGE_FREETRANS:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      SCIP_CALL( SCIPsolFree(sol, scip->mem->probmem, scip->primal) );
      break;
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** links a primal solution to the current LP solution */
SCIP_RETCODE SCIPlinkLPSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPlinkLPSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPlpIsSolved(scip->lp) )
   {
      SCIPerrorMessage("LP solution does not exist\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolLinkLPSol(sol, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current relaxation solution */
SCIP_RETCODE SCIPlinkRelaxSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPlinkRelaxSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( !SCIPrelaxationIsSolValid(scip->relaxation) )
   {
      SCIPerrorMessage("relaxation solution is not valid\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolLinkRelaxSol(sol, scip->set, scip->stat, scip->tree, scip->relaxation) );

   return SCIP_OKAY;
}

/** links a primal solution to the current pseudo solution */
SCIP_RETCODE SCIPlinkPseudoSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPlinkPseudoSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolLinkPseudoSol(sol, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** links a primal solution to the current LP or pseudo solution */
SCIP_RETCODE SCIPlinkCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPlinkCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPsolLinkCurrentSol(sol, scip->set, scip->stat, scip->tree, scip->lp) );

   return SCIP_OKAY;
}

/** clears a primal solution */
SCIP_RETCODE SCIPclearSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPclearSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsolClear(sol, scip->stat, scip->tree) );

   return SCIP_OKAY;
}

/** stores solution values of variables in solution's own array */
SCIP_RETCODE SCIPunlinkSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPunlinkSol", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsolUnlink(sol, scip->set, scip->transprob) );

   return SCIP_OKAY;
}

/** sets value of variable in primal CIP solution */
SCIP_RETCODE SCIPsetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to add to solution */
   SCIP_Real             val                 /**< solution value of variable */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetSolVal", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL && SCIPvarIsTransformed(var) )
   {
      SCIPerrorMessage("cannot set value of transformed variable <%s> in original space solution\n",
         SCIPvarGetName(var));
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolSetVal(sol, scip->set, scip->stat, scip->tree, var, val) );

   return SCIP_OKAY;
}

/** sets values of multiple variables in primal CIP solution */
SCIP_RETCODE SCIPsetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   int                   nvars,              /**< number of variables to set solution value for */
   SCIP_VAR**            vars,               /**< array with variables to add to solution */
   SCIP_Real*            vals                /**< array with solution values of variables */
   )
{
   int v;

   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( checkStage(scip, "SCIPsetSolVals", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      for( v = 0; v < nvars; ++v )
      {
         if( SCIPvarIsTransformed(vars[v]) )
         {
            SCIPerrorMessage("cannot set value of transformed variable <%s> in original space solution\n",
               SCIPvarGetName(vars[v]));
            return SCIP_INVALIDCALL;
         }
      }
   }

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPsolSetVal(sol, scip->set, scip->stat, scip->tree, vars[v], vals[v]) );
   }

   return SCIP_OKAY;
}

/** increases value of variable in primal CIP solution */
SCIP_RETCODE SCIPincSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_VAR*             var,                /**< variable to increase solution value for */
   SCIP_Real             incval              /**< increment for solution value of variable */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPincSolVal", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL && SCIPvarIsTransformed(var) )
   {
      SCIPerrorMessage("cannot increase value of transformed variable <%s> in original space solution\n",
         SCIPvarGetName(var));
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolIncVal(sol, scip->set, scip->stat, scip->tree, var, incval) );

   return SCIP_OKAY;
}

/** returns value of variable in primal CIP solution, or in current LP/pseudo solution */
SCIP_Real SCIPgetSolVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolVal", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL && SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL && SCIPvarIsTransformed(var) )
   {
      SCIP_VAR* origvar;
      SCIP_Real scalar;
      SCIP_Real constant;

      /* we cannot get the value of a transformed variable for a solution that lives in the original problem space
       * -> get the corresponding original variable first
       */
      origvar = var;
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL_ABORT( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );
      if( origvar == NULL )
      {
         /* the variable has no original counterpart: in the original solution, it has a value of zero */
         return 0.0;
      }
      assert(!SCIPvarIsTransformed(origvar));
      return scalar * SCIPgetSolVal(scip, sol, origvar) + constant;
   }
   else if( sol != NULL )
      return SCIPsolGetVal(sol, scip->set, scip->stat, var);
   else
   {
      SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolVal(sol==NULL)", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      return SCIPvarGetSol(var, SCIPtreeHasCurrentNodeLP(scip->tree));
   }
}

/** gets values of multiple variables in primal CIP solution */
SCIP_RETCODE SCIPgetSolVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   int                   nvars,              /**< number of variables to get solution value for */
   SCIP_VAR**            vars,               /**< array with variables to get value for */
   SCIP_Real*            vals                /**< array to store solution values of variables */
   )
{
   assert(nvars == 0 || vars != NULL);
   assert(nvars == 0 || vals != NULL);

   SCIP_CALL( checkStage(scip, "SCIPgetSolVals", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
   {
      int v;

      if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
      {
         SCIP_VAR* origvar;
         SCIP_Real scalar;
         SCIP_Real constant;

         for( v = 0; v < nvars; ++v )
         {
            /* we cannot get the value of a transformed variable for a solution that lives in the original problem space
             * -> get the corresponding original variable first
             */
            origvar = vars[v];
            scalar = 1.0;
            constant = 0.0;
            SCIP_CALL_ABORT( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );
            if( origvar == NULL )
            {
               /* the variable has no original counterpart: in the original solution, it has a value of zero */
               vals[v] = 0.0;
            }
            else
            {
               assert(!SCIPvarIsTransformed(origvar));
               vals[v] = scalar * SCIPgetSolVal(scip, sol, origvar) + constant;
            }
         }
      }
      else
      {
         for( v = 0; v < nvars; ++v )
            vals[v] = SCIPsolGetVal(sol, scip->set, scip->stat, vars[v]);
      }
   }
   else
   {
      SCIP_CALL( SCIPgetVarSols(scip, nvars, vars, vals) );
   }

   return SCIP_OKAY;
}

/** returns objective value of primal CIP solution w.r.t. original problem, or current LP/pseudo objective value */
SCIP_Real SCIPgetSolOrigObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   )
{
   /* for original solutions, an original objective value is already available in SCIP_STAGE_PROBLEM
    * for all other solutions, we should be at least in SCIP_STAGE_TRANSFORMING
    */
   if( sol != NULL && SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolOrigObj", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

      return SCIPsolGetOrigObj(sol);
   }

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolOrigObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      return SCIPprobExternObjval(scip->transprob, scip->set, SCIPsolGetObj(sol, scip->set, scip->transprob));
   else
   {
      SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolOrigObj(sol==NULL)",
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      if( SCIPtreeHasCurrentNodeLP(scip->tree) )
         return SCIPprobExternObjval(scip->transprob, scip->set, SCIPlpGetObjval(scip->lp, scip->set));
      else
         return SCIPprobExternObjval(scip->transprob, scip->set, SCIPlpGetPseudoObjval(scip->lp, scip->set));
   }
}

/** returns transformed objective value of primal CIP solution, or transformed current LP/pseudo objective value */
SCIP_Real SCIPgetSolTransObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution, or NULL for current LP/pseudo objective value */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolTransObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( sol != NULL )
      return SCIPsolGetObj(sol, scip->set, scip->transprob);
   else
   {
      SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolTransObj(sol==NULL)",
            FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );
      if( SCIPtreeHasCurrentNodeLP(scip->tree) )
         return SCIPlpGetObjval(scip->lp, scip->set);
      else
         return SCIPlpGetPseudoObjval(scip->lp, scip->set);
   }
}

/** maps original space objective value into transformed objective value */
SCIP_Real SCIPtransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< original space objective value to transform */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPtransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobInternObjval(scip->transprob, scip->set, obj);
}

/** maps transformed objective value into original space */
SCIP_Real SCIPretransformObj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             obj                 /**< transformed objective value to retransform in original space */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPretransformObj", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->set, obj);
}

/** gets clock time, when this solution was found */
SCIP_Real SCIPgetSolTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolTime", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetTime(sol);
}

/** gets branch and bound run number, where this solution was found */
int SCIPgetSolRunnum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolRunnum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetRunnum(sol);
}

/** gets node number of the specific branch and bound run, where this solution was found */
SCIP_Longint SCIPgetSolNodenum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolNodenum", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetNodenum(sol);
}

/** gets heuristic, that found this solution (or NULL if it's from the tree) */
SCIP_HEUR* SCIPgetSolHeur(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal solution */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolHeur", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolGetHeur(sol);
}

/** returns whether two given solutions are exactly equal */
SCIP_Bool SCIPareSolsEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol1,               /**< first primal CIP solution */
   SCIP_SOL*             sol2                /**< second primal CIP solution */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPareSolsEqual", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPsolsAreEqual(sol1, sol2, scip->set, scip->stat, scip->origprob, scip->transprob);
}

/** outputs non-zero variables of solution in original problem space to file stream */
SCIP_RETCODE SCIPprintSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Real objvalue;
   SCIP_Bool currentsol;

   assert(SCIPisTransformed(scip) || sol != NULL);

   SCIP_CALL( checkStage(scip, "SCIPprintSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   
   currentsol = (sol == NULL);
   if( currentsol )
   {
      SCIP_CALL( checkStage(scip, "SCIPprintSol(sol==NULL)",
            FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
      
      /* create a temporary solution that is linked to the current solution */
      SCIP_CALL( SCIPsolCreateCurrentSol(&sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree,
            scip->lp, NULL) );
   }

   SCIPmessageFPrintInfo(file, "objective value:                 ");
   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
      objvalue = SCIPsolGetOrigObj(sol);
   else
      objvalue = SCIPprobExternObjval(scip->transprob, scip->set, SCIPsolGetObj(sol, scip->set, scip->transprob));

   SCIPprintReal(scip, file, objvalue, 20, 15);
   SCIPmessageFPrintInfo(file, "\n");

   SCIP_CALL( SCIPsolPrint(sol, scip->set, scip->stat, scip->origprob, scip->transprob, file, printzeros) );

   if( currentsol )
   {
      /* free temporary solution */
      SCIP_CALL( SCIPsolFree(&sol, scip->mem->probmem, scip->primal) );
   }

   return SCIP_OKAY;
}

/** outputs non-zero variables of solution in transformed problem space to file stream */
SCIP_RETCODE SCIPprintTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution, or NULL for current LP/pseudo solution */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_Bool currentsol;

   SCIP_CALL( checkStage(scip, "SCIPprintTransSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   currentsol = (sol == NULL);
   if( currentsol )
   {
      /* create a temporary solution that is linked to the current solution */
      SCIP_CALL( SCIPsolCreateCurrentSol(&sol, scip->mem->probmem, scip->set, scip->stat, scip->primal, scip->tree,
            scip->lp, NULL) );
   }

   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIPerrorMessage("cannot print original space solution as transformed solution\n");
      return SCIP_INVALIDCALL;
   }

   SCIPmessageFPrintInfo(file, "objective value:                 ");
   SCIPprintReal(scip, file, SCIPsolGetObj(sol, scip->set, scip->transprob), 20, 9);
   SCIPmessageFPrintInfo(file, "\n");

   SCIP_CALL( SCIPsolPrint(sol, scip->set, scip->stat, scip->transprob, NULL, file, printzeros) );

   if( currentsol )
   {
      /* free temporary solution */
      SCIP_CALL( SCIPsolFree(&sol, scip->mem->probmem, scip->primal) );
   }

   return SCIP_OKAY;
}

/** gets number of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is SCIP_STAGE_PROBLEM, it returns the number solution in the original solution candidate
 *  storage
 */
int SCIPgetNSols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNSols", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSols", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprimal->nsols;
      
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      return scip->primal->nsols;
      
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** gets array of feasible primal solutions stored in the solution storage in case the problem is transformed; in case
 *  if the problem stage is in SCIP_STAGE_PROBLEM, it returns the number array of solution candidate stored
 */
SCIP_SOL** SCIPgetSols(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSols", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      return scip->origprimal->sols;
      
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      return scip->primal->sols;
      
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return NULL;
   }  /*lint !e788*/
}

/** gets best feasible primal solution found so far if the problem is transformed; in case the problem is in
 *  SCIP_STAGE_PROBLEM it returns the best solution candidate, or NULL if no solution has been found or the candidate
 *  store is empty;
 */
SCIP_SOL* SCIPgetBestSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBestSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );
   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(scip->origprimal != NULL);
      if(  scip->origprimal->nsols > 0 )
      {
         assert(scip->origprimal->sols != NULL);
         assert(scip->origprimal->sols[0] != NULL);
         return scip->origprimal->sols[0];
      }
      break;
      
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
      assert(scip->primal != NULL);
      if(  scip->primal->nsols > 0 )
      {
         assert(scip->primal->sols != NULL);
         assert(scip->primal->sols[0] != NULL);
         return scip->primal->sols[0];
      }
      break;
      
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_FREETRANS:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return NULL;
   }

   return NULL;
}

/** outputs best feasible primal solution found so far to file stream */
SCIP_RETCODE SCIPprintBestSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_SOL* sol;

   SCIP_CALL( checkStage(scip, "SCIPprintBestSol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   sol = SCIPgetBestSol(scip);

   if( sol == NULL )
      SCIPmessageFPrintInfo(file, "no solution available\n");
   else
   {
      SCIP_CALL( SCIPprintSol(scip, sol, file, printzeros) );
   }

   return SCIP_OKAY;
}

/** outputs best feasible primal solution found so far in transformed variables to file stream */
SCIP_RETCODE SCIPprintBestTransSol(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Bool             printzeros          /**< should variables set to zero be printed? */
   )
{
   SCIP_SOL* sol;

   SCIP_CALL( checkStage(scip, "SCIPprintBestTransSol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   sol = SCIPgetBestSol(scip);

   if( sol != NULL && SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIPerrorMessage("best solution is defined in original space - cannot print it as transformed solution\n");
      return SCIP_INVALIDCALL;
   }

   if( sol == NULL )
      SCIPmessageFPrintInfo(file, "no solution available\n");
   else
   {
      SCIP_CALL( SCIPprintTransSol(scip, sol, file, printzeros) );
   }

   return SCIP_OKAY;
}

/** try to round given solution */
SCIP_RETCODE SCIProundSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal solution */
   SCIP_Bool*            success             /**< pointer to store whether rounding was successful */
   )
{
   SCIP_CALL( checkStage(scip, "SCIProundSol", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIPerrorMessage("cannot round original space solution\n");
      return SCIP_INVALIDCALL;
   }

   SCIP_CALL( SCIPsolRound(sol, scip->set, scip->stat, scip->transprob, scip->tree, success) );

   return SCIP_OKAY;
}

/** retransforms solution to original problem space */
SCIP_RETCODE SCIPretransformSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< primal CIP solution */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPretransformSol", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   switch ( SCIPsolGetOrigin(sol) )
   {
   case SCIP_SOLORIGIN_ORIGINAL:
      /* nothing to do */
      return SCIP_OKAY;

   case SCIP_SOLORIGIN_LPSOL:
   case SCIP_SOLORIGIN_NLPSOL:
   case SCIP_SOLORIGIN_RELAXSOL:
   case SCIP_SOLORIGIN_PSEUDOSOL:

      /* first unlink solution */
      SCIP_CALL( SCIPunlinkSol(scip, sol) );
      
      /*lint -fallthrough*/
   case SCIP_SOLORIGIN_ZERO:

      SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob) );
      break;

   case SCIP_SOLORIGIN_UNKNOWN:
      SCIPerrorMessage("unkown solution origin.\n");
      return SCIP_INVALIDCALL;

   default:
      SCIPerrorMessage("invalid solution origin <%d>\n", SCIPsolGetOrigin(sol));
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** reads a given solution file, problem has to be transformed in advance */
SCIP_RETCODE SCIPreadSol(
   SCIP*                 scip,              /**< SCIP data structure */
   const char*           fname              /**< name of the input file */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPreadSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE) );

   /* we pass the reading of the solution file on to reader_sol via the following call */
   SCIP_CALL( SCIPreadProb(scip, fname, "sol") );

   return SCIP_OKAY;
}

/** adds feasible primal solution to solution storage by copying it */
SCIP_RETCODE SCIPaddSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddSol", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_FREETRANS:
      assert( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL );
      SCIP_CALL( SCIPprimalAddOrigSol(scip->origprimal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, sol, stored) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* if the solution is added during presolving and it is not defined on original variables, 
       * presolving operations will destroy its validity, so we retransform it to the original space
       */
      if( SCIPsolGetOrigin(sol) != SCIP_SOLORIGIN_ORIGINAL )
      {
         SCIP_CALL( SCIPsolUnlink(sol, scip->set, scip->transprob) );
         SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob) );      
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPprimalAddSol(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob, scip->tree,
            scip->lp, scip->eventqueue, scip->eventfilter, sol, stored) );
      return SCIP_OKAY;
      
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** adds primal solution to solution storage, frees the solution afterwards */
SCIP_RETCODE SCIPaddSolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddSolFree", FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_FREETRANS:
      assert( SCIPsolGetOrigin(*sol) == SCIP_SOLORIGIN_ORIGINAL );
      SCIP_CALL( SCIPprimalAddOrigSolFree(scip->origprimal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, sol, stored) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
      /* if the solution is added during presolving and it is not defined on original variables, 
       * presolving operations will destroy its validity, so we retransform it to the original space
       */
      if( SCIPsolGetOrigin(*sol) != SCIP_SOLORIGIN_ORIGINAL )
      {
         SCIP_CALL( SCIPsolUnlink(*sol, scip->set, scip->transprob) );
         SCIP_CALL( SCIPsolRetransform(*sol, scip->set, scip->stat, scip->origprob) );      
      }
      /*lint -fallthrough*/
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPprimalAddSolFree(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob, scip->tree,
            scip->lp, scip->eventqueue, scip->eventfilter, sol, stored) );
      return SCIP_OKAY;
      
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_ERROR;
   }  /*lint !e788*/
}

/** adds current LP/pseudo solution to solution storage */
SCIP_RETCODE SCIPaddCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPaddCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPprimalAddCurrentSol(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob,
         scip->tree, scip->lp, scip->eventqueue, scip->eventfilter, heur, stored) );

   return SCIP_OKAY;
}

/** checks solution for feasibility; if possible, adds it to storage by copying */
SCIP_RETCODE SCIPtrySol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< should all reasons of violation be printed? */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   assert(stored != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtrySol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* if the solution is added during presolving and it is not defined on original variables, 
    * presolving operations will destroy its validity, so we retransform it to the original space
    */ 
   if( scip->set->stage == SCIP_STAGE_PRESOLVING && SCIPsolGetOrigin(sol) != SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIP_CALL( SCIPsolUnlink(sol, scip->set, scip->transprob) );
      SCIP_CALL( SCIPsolRetransform(sol, scip->set, scip->stat, scip->origprob) );      
   }

   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIP_Bool feasible;

      /* SCIPprimalTrySol() can only be called on transformed solutions; therefore check solutions in original problem
       * including modifiable constraints */
      SCIP_CALL( checkSolOrig(scip, sol, &feasible, printreason, FALSE, checkbounds, checkintegrality, checklprows, TRUE) );
      if( feasible )
      {
         SCIP_CALL( SCIPprimalAddSol(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob,
               scip->tree, scip->lp, scip->eventqueue, scip->eventfilter, sol, stored) );
      }
      else
         *stored = FALSE;
   }
   else
   {
      SCIP_CALL( SCIPprimalTrySol(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob, scip->tree,
            scip->lp, scip->eventqueue, scip->eventfilter, sol, printreason, checkbounds, checkintegrality, checklprows, stored) );
   }

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
SCIP_RETCODE SCIPtrySolFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   )
{
   assert(stored != NULL);
   assert(sol != NULL);

   SCIP_CALL( checkStage(scip, "SCIPtrySolFree", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   /* if the solution is added during presolving and it is not defined on original variables, 
    * presolving operations will destroy its validity, so we retransform it to the original space
    */
   if( scip->set->stage == SCIP_STAGE_PRESOLVING && SCIPsolGetOrigin(*sol) != SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIP_CALL( SCIPsolUnlink(*sol, scip->set, scip->transprob) );
      SCIP_CALL( SCIPsolRetransform(*sol, scip->set, scip->stat, scip->origprob) );      
   }

   if( SCIPsolGetOrigin(*sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      SCIP_Bool feasible;

      /* SCIPprimalTrySol() can only be called on transformed solutions; therefore check solutions in original problem 
       * including modifiable constraints 
       */
      SCIP_CALL( checkSolOrig(scip, *sol, &feasible, printreason, FALSE, checkbounds, checkintegrality, checklprows, TRUE) );
      
      if( feasible )
      {
         SCIP_CALL( SCIPprimalAddSolFree(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob,
               scip->tree, scip->lp, scip->eventqueue, scip->eventfilter, sol, stored) );
      }
      else
      {
         SCIP_CALL( SCIPsolFree(sol, scip->mem->probmem, scip->primal) );
         *stored = FALSE;
      }
   }
   else
   {
      SCIP_CALL( SCIPprimalTrySolFree(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob,
            scip->tree, scip->lp, scip->eventqueue, scip->eventfilter, sol, printreason, checkbounds, checkintegrality, checklprows, stored) );
   }

   return SCIP_OKAY;
}

/** checks current LP/pseudo solution for feasibility; if possible, adds it to storage */
SCIP_RETCODE SCIPtryCurrentSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< heuristic that found the solution */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPtryCurrentSol", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPprimalTryCurrentSol(scip->primal, scip->mem->probmem, scip->set, scip->stat, scip->origprob, scip->transprob,
         scip->tree, scip->lp, scip->eventqueue, scip->eventfilter, heur, printreason, checkintegrality, checklprows, stored) );

   return SCIP_OKAY;
}

/** checks solution for feasibility without adding it to the solution store */
SCIP_RETCODE SCIPcheckSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool             printreason,        /**< should all reasons of violations be printed? */
   SCIP_Bool             checkbounds,        /**< should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool*            feasible            /**< stores whether given solution is feasible */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcheckSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || scip->set->misc_exactsolve;

   if( SCIPsolGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
   {
      /* SCIPsolCheck() can only be called on transformed solutions */
      SCIP_CALL( checkSolOrig(scip, sol, feasible, printreason, FALSE, checkbounds, checkintegrality, checklprows, FALSE) );
   }
   else
   {
      SCIP_CALL( SCIPsolCheck(sol, scip->mem->probmem, scip->set, scip->stat, scip->transprob, printreason,
            checkbounds, checkintegrality, checklprows, feasible) );
   }

   return SCIP_OKAY;
}

/** checks solution for feasibility in original problem without adding it to the solution store;
 *  this method is used to double check a solution in order to validate the presolving process
 */
SCIP_RETCODE SCIPcheckSolOrig(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< primal CIP solution */
   SCIP_Bool*            feasible,           /**< stores whether given solution is feasible */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool             completely          /**< should all violation be checked? */
   )
{
   assert(scip != NULL);
   assert(sol != NULL);
   assert(feasible != NULL);

   SCIP_CALL( checkStage(scip, "SCIPcheckSolOrig", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   /* check solution in original problem; that includes bounds, integrality, and non modifiable constraints */
   SCIP_CALL( checkSolOrig(scip, sol, feasible, printreason, completely, TRUE, TRUE, TRUE, FALSE) );
   
   return SCIP_OKAY;
}

/** return whether a primal ray is stored that proves unboundedness of the LP relaxation */
SCIP_Bool SCIPhasPrimalRay(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPhasPrimalRay", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->primal->primalray != NULL;
}

/** gets value of given variable in primal ray causing unboundedness of the LP relaxation; 
 *  should only be called if such a ray is stored (check with SCIPhasPrimalRay()) */
SCIP_Real SCIPgetPrimalRayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var                 /**< variable to get value for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPrimalRayVal", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   assert(var != NULL);
   assert(scip->primal->primalray != NULL);

   return SCIPsolGetRayVal(scip->primal->primalray, scip->set, scip->stat, var);
}



/*
 * event methods
 */

/** catches a global (not variable or row dependent) event */
SCIP_RETCODE SCIPcatchEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcatchEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPeventfilterAdd(scip->eventfilter, scip->mem->probmem, scip->set,
         eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops a global event (stops to track event) */
SCIP_RETCODE SCIPdropEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchEvent(), or -1 */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdropEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPeventfilterDel(scip->eventfilter, scip->mem->probmem, scip->set,
         eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** catches an objective value or domain change event on the given transformed variable */
SCIP_RETCODE SCIPcatchVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< transformed variable to catch event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcatchVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( (eventtype & SCIP_EVENTTYPE_VARCHANGED) == 0 )
   {
      SCIPerrorMessage("event does not operate on a single variable\n");
      return SCIP_INVALIDDATA;
   }

   if( SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("cannot catch events on original variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPvarCatchEvent(var, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops an objective value or domain change event (stops to track event) on the given transformed variable */
SCIP_RETCODE SCIPdropVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< transformed variable to drop event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdropVarEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( SCIPvarIsOriginal(var) )
   {
      SCIPerrorMessage("cannot drop events on original variable <%s>\n", SCIPvarGetName(var));
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIPvarDropEvent(var, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** catches a row coefficient, constant, or side change event on the given row */
SCIP_RETCODE SCIPcatchRowEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< linear row to catch event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask to select events to catch */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int*                  filterpos           /**< pointer to store position of event filter entry, or NULL */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcatchRowEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   if( (eventtype & SCIP_EVENTTYPE_ROWCHANGED) == 0 )
   {
      SCIPerrorMessage("event does not operate on a single row\n");
      return SCIP_INVALIDDATA;
   }

   SCIP_CALL( SCIProwCatchEvent(row, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );

   return SCIP_OKAY;
}

/** drops a row coefficient, constant, or side change event (stops to track event) on the given row */
SCIP_RETCODE SCIPdropRowEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row,                /**< linear row to drop event for */
   SCIP_EVENTTYPE        eventtype,          /**< event type mask of dropped event */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler to process events with */
   SCIP_EVENTDATA*       eventdata,          /**< event data to pass to the event handler when processing this event */
   int                   filterpos           /**< position of event filter entry returned by SCIPcatchVarEvent(), or -1 */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPdropRowEvent", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIProwDropEvent(row, scip->mem->probmem, scip->set, eventtype, eventhdlr, eventdata, filterpos) );
   
   return SCIP_OKAY;
}


/*
 * tree methods
 */

/** gets current node in the tree */
SCIP_NODE* SCIPgetCurrentNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetCurrentNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetCurrentNode(scip->tree);
}

/** gets the root node of the tree */
SCIP_NODE* SCIPgetRootNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRootNode", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetRootNode(scip->tree);
}

/** returns whether the current node is already solved and only propagated again */
SCIP_Bool SCIPinRepropagation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPinRepropagation", FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeInRepropagation(scip->tree);
}

/** gets children of focus node along with the number of children */
SCIP_RETCODE SCIPgetChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          children,           /**< pointer to store children array, or NULL if not needed */
   int*                  nchildren           /**< pointer to store number of children, or NULL if not needed */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetChildren", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( children != NULL )
      *children = scip->tree->children;
   if( nchildren != NULL )
      *nchildren = scip->tree->nchildren;

   return SCIP_OKAY;
}

/** gets number of children of focus node */
int SCIPgetNChildren(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNChildren", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->nchildren;
}

/** gets siblings of focus node along with the number of siblings */
SCIP_RETCODE SCIPgetSiblings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          siblings,           /**< pointer to store siblings array, or NULL if not needed */
   int*                  nsiblings           /**< pointer to store number of siblings, or NULL if not needed */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( siblings != NULL )
      *siblings = scip->tree->siblings;
   if( nsiblings != NULL )
      *nsiblings = scip->tree->nsiblings;

   return SCIP_OKAY;
}

/** gets number of siblings of focus node */
int SCIPgetNSiblings(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNSiblings", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->nsiblings;
}

/** gets leaves of the tree along with the number of leaves */
SCIP_RETCODE SCIPgetLeaves(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE***          leaves,             /**< pointer to store leaves array, or NULL if not needed */
   int*                  nleaves             /**< pointer to store number of leaves, or NULL if not needed */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPgetLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( leaves != NULL )
      *leaves = SCIPnodepqNodes(scip->tree->leaves);
   if( nleaves != NULL )
      *nleaves = SCIPnodepqLen(scip->tree->leaves);

   return SCIP_OKAY;
}

/** gets number of leaves in the tree */
int SCIPgetNLeaves(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLeaves", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPnodepqLen(scip->tree->leaves);
}

/** gets the best child of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPgetPrioChild(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPrioChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetPrioChild(scip->tree);
}

/** gets the best sibling of the focus node w.r.t. the node selection priority assigned by the branching rule */
SCIP_NODE* SCIPgetPrioSibling(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPrioSibling", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetPrioSibling(scip->tree);
}

/** gets the best child of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPgetBestChild(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBestChild", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestChild(scip->tree, scip->set);
}

/** gets the best sibling of the focus node w.r.t. the node selection strategy */
SCIP_NODE* SCIPgetBestSibling(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBestSibling", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestSibling(scip->tree, scip->set);
}

/** gets the best leaf from the node queue w.r.t. the node selection strategy */
SCIP_NODE* SCIPgetBestLeaf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBestLeaf", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestLeaf(scip->tree);
}

/** gets the best node from the tree (child, sibling, or leaf) w.r.t. the node selection strategy */
SCIP_NODE* SCIPgetBestNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBestNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetBestNode(scip->tree, scip->set);
}

/** gets the node with smallest lower bound from the tree (child, sibling, or leaf) */
SCIP_NODE* SCIPgetBestboundNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBestboundNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return SCIPtreeGetLowerboundNode(scip->tree, scip->set);
}

/** cuts off node and whole sub tree from branch and bound tree */
SCIP_RETCODE SCIPcutoffNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be cut off */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcutoffNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodeCutoff(node, scip->set, scip->stat, scip->tree);

   return SCIP_OKAY;
}

/** marks the given node to be propagated again the next time a node of its subtree is processed */
SCIP_RETCODE SCIPrepropagateNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node                /**< node that should be propagated again */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPrepropagateNode", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   SCIPnodePropagateAgain(node, scip->set, scip->stat, scip->tree);

   return SCIP_OKAY;
}

/** returns depth of first node in active path that is marked being cutoff */
int SCIPgetCutoffdepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetCutoffdepth", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->cutoffdepth;
}

/** returns depth of first node in active path that has to be propagated again */
int SCIPgetRepropdepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRepropdepth", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->tree->repropdepth;
}

/* prints all branching decisions on variables from the root to the given node */
SCIP_RETCODE SCIPprintNodeRootPath(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NODE*            node,               /**< node data */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR**            branchvars;         /**< array of variables on which the branchings has been performed in all ancestors */                        
   SCIP_Real*            branchbounds;       /**< array of bounds which the branchings in all ancestors set */                                             
   SCIP_BOUNDTYPE*       boundtypes;         /**< array of boundtypes which the branchings in all ancestors set */                                         
   int*                  nodeswitches;       /**< marks, where in the arrays the branching decisions of the next node on the path start                    
                                              * branchings performed at the parent of node always start at position 0. For single variable branching,      
                                              * nodeswitches[i] = i holds */                                                                               
   int                   nbranchvars;        /**< number of variables on which branchings have been performed in all ancestors                             
                                              *   if this is larger than the array size, arrays should be reallocated and method should be called again */ 
   int                   branchvarssize;     /**< available slots in arrays */                                                                             
   int                   nnodes;             /* number of nodes in the nodeswitch array */                                                                 
   int                   nodeswitchsize;     /**< available slots in node switch array */                                                                  
   
   branchvarssize = SCIPnodeGetDepth(node);
   nodeswitchsize = branchvarssize;
   
   /* memory allocation */
   SCIP_CALL( SCIPallocBufferArray(scip, &branchvars, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &branchbounds, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtypes, branchvarssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeswitches, nodeswitchsize) );

   SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize );
   
   /* if the arrays were to small, we have to reallocate them and recall SCIPnodeGetAncestorBranchingPath */
   if( nbranchvars > branchvarssize || nnodes > nodeswitchsize )
   {
      branchvarssize = nbranchvars;
      nodeswitchsize = nnodes;

      /* memory reallocation */
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchvars, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &branchbounds, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &boundtypes, branchvarssize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &nodeswitches, nodeswitchsize) );
      
      SCIPnodeGetAncestorBranchingPath(node, branchvars, branchbounds, boundtypes, &nbranchvars, branchvarssize, nodeswitches, &nnodes, nodeswitchsize);
      assert(nbranchvars == branchvarssize);
   }
   
   /* we only want to create output, if branchings were performed */
   if( nbranchvars >= 1 )
   {
      int i;
      int j;

      /* print all nodes, starting from the root, which is last in the arrays */
      for( j = nnodes-1; j >= 0; --j)
      {
         int end;
         if(j == nnodes-1)
            end =  nbranchvars;
         else 
            end =  nodeswitches[j+1];
         
         for( i = nodeswitches[j]; i < end; ++i )
         {
            if( i > nodeswitches[j] )
               SCIPmessageFPrintInfo(file, " AND ");
            SCIPmessageFPrintInfo(file, "<%s> %s %.1f",SCIPvarGetName(branchvars[i]), boundtypes[i] == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", branchbounds[i]);
         }
         SCIPmessageFPrintInfo(file, "\n");            
         if( j > 0 )
         {
            if(  nodeswitches[j]-nodeswitches[j-1] != 1 )
               SCIPmessageFPrintInfo(file, " |\n |\n");
            else if( boundtypes[i-1] == SCIP_BOUNDTYPE_LOWER )
               SCIPmessageFPrintInfo(file, "\\ \n \\\n");
            else
               SCIPmessageFPrintInfo(file, " /\n/ \n");
         }
      }
   }
   
   /* free all local memory */
   SCIPfreeBufferArray(scip, &nodeswitches);
   SCIPfreeBufferArray(scip, &boundtypes);
   SCIPfreeBufferArray(scip, &branchbounds);
   SCIPfreeBufferArray(scip, &branchvars);
   
   return SCIP_OKAY;
}


/*
 * statistic methods
 */

/** gets number of branch and bound runs performed, including the current run */
int SCIPgetNRuns(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNRuns", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->nruns;
}

/** gets number of processed nodes in current run, including the focus node */
SCIP_Longint SCIPgetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNodes", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->nnodes;
}

/** gets total number of processed nodes in all runs, including the focus node */
SCIP_Longint SCIPgetNTotalNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNTotalNodes", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return scip->stat->ntotalnodes;
}

/** gets number of nodes left in the tree (children + siblings + leaves) */
int SCIPgetNNodesLeft(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNodesLeft", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetNNodes(scip->tree);
}

/** gets total number of LPs solved so far */
int SCIPgetNLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nlps;
}

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm */
SCIP_Longint SCIPgetNLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nlpiterations;
}

/** gets total number of iterations used so far in primal and dual simplex and barrier algorithm for the root node */
SCIP_Longint SCIPgetNRootLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNRootLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nrootlpiterations;
}

/** gets total number of primal LPs solved so far */
int SCIPgetNPrimalLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrimalLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nprimallps;
}

/** gets total number of iterations used so far in primal simplex */
SCIP_Longint SCIPgetNPrimalLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrimalLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nprimallpiterations;
}

/** gets total number of dual LPs solved so far */
int SCIPgetNDualLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNDualLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nduallps;
}

/** gets total number of iterations used so far in dual simplex */
SCIP_Longint SCIPgetNDualLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNDualLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nduallpiterations;
}

/** gets total number of barrier LPs solved so far */
int SCIPgetNBarrierLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNBarrierLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nbarrierlps;
}

/** gets total number of iterations used so far in barrier algorithm */
SCIP_Longint SCIPgetNBarrierLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNBarrierLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nbarrierlpiterations;
}

/** gets total number of LPs solved so far that were resolved from an advanced start basis */
int SCIPgetNResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNResolveLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelps + scip->stat->ndualresolvelps;
}

/** gets total number of simplex iterations used so far in primal and dual simplex calls where an advanced start basis
 *  was available
 */
SCIP_Longint SCIPgetNResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNResolveLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelpiterations + scip->stat->ndualresolvelpiterations;
}

/** gets total number of primal LPs solved so far that were resolved from an advanced start basis */
int SCIPgetNPrimalResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrimalResolveLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelps;
}

/** gets total number of simplex iterations used so far in primal simplex calls where an advanced start basis
 *  was available
 */
SCIP_Longint SCIPgetNPrimalResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPrimalResolveLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nprimalresolvelpiterations;
}

/** gets total number of dual LPs solved so far that were resolved from an advanced start basis */
int SCIPgetNDualResolveLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNDualResolveLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ndualresolvelps;
}

/** gets total number of simplex iterations used so far in dual simplex calls where an advanced start basis
 *  was available
 */
SCIP_Longint SCIPgetNDualResolveLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNDualResolveLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ndualresolvelpiterations;
}

/** gets total number of LPs solved so far for node relaxations */
int SCIPgetNNodeLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNodeLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nnodelps;
}

/** gets total number of simplex iterations used so far for node relaxations */
SCIP_Longint SCIPgetNNodeLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNodeLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nnodelpiterations;
}

/** gets total number of LPs solved so far for initial LP in node relaxations */
int SCIPgetNNodeInitLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNInitNodeLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ninitlps;
}

/** gets total number of simplex iterations used so far for initial LP in node relaxations */
SCIP_Longint SCIPgetNNodeInitLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNNodeInitLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ninitlpiterations;
}

/** gets total number of LPs solved so far during diving and probing */
int SCIPgetNDivingLPs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNDivingLPs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ndivinglps;
}

/** gets total number of simplex iterations used so far during diving and probing */
SCIP_Longint SCIPgetNDivingLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNDivingLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->ndivinglpiterations;
}

/** gets total number of times, strong branching was called (each call represents solving two LPs) */
int SCIPgetNStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNStrongbranchs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nstrongbranchs;
}

/** gets total number of simplex iterations used so far in strong branching */
SCIP_Longint SCIPgetNStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNStrongbranchLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nsblpiterations;
}

/** gets total number of times, strong branching was called at the root node (each call represents solving two LPs) */
int SCIPgetNRootStrongbranchs(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNRootStrongbranchs", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nrootstrongbranchs;
}

/** gets total number of simplex iterations used so far in strong branching at the root node */
SCIP_Longint SCIPgetNRootStrongbranchLPIterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNRootStrongbranchLPIterations", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nrootsblpiterations;
}

/** gets number of pricing rounds performed so far at the current node */
int SCIPgetNPriceRounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPriceRounds", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->npricerounds;
}

/** get current number of variables in the pricing store */
int SCIPgetNPricevars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPricevars", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPpricestoreGetNVars(scip->pricestore);
}

/** get total number of pricing variables found so far */
int SCIPgetNPricevarsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPricevarsFound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPpricestoreGetNVarsFound(scip->pricestore);
}

/** get total number of pricing variables applied to the LPs */
int SCIPgetNPricevarsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNPricevarsApplied", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPpricestoreGetNVarsApplied(scip->pricestore);
}

/** gets number of separation rounds performed so far at the current node */
int SCIPgetNSepaRounds(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNSepaRounds", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nseparounds;
}

/** get total number of cuts found so far */
int SCIPgetNCutsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNCutsFound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCutsFound(scip->sepastore);
}

/** get number of cuts found so far in current separation round */
int SCIPgetNCutsFoundRound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNCutsFoundRound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCutsFoundRound(scip->sepastore);
}

/** get total number of cuts applied to the LPs */
int SCIPgetNCutsApplied(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNCutsApplied", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPsepastoreGetNCutsApplied(scip->sepastore);
}

/** get total number of constraints found in conflict analysis (conflict and reconvergence constraints) */
SCIP_Longint SCIPgetNConflictConssFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNConflictConssFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPconflictGetNPropConflictConss(scip->conflict)
      + SCIPconflictGetNPropReconvergenceConss(scip->conflict)
      + SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict)
      + SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict)
      + SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict)
      + SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict)
      + SCIPconflictGetNStrongbranchConflictConss(scip->conflict)
      + SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict)
      + SCIPconflictGetNPseudoConflictConss(scip->conflict)
      + SCIPconflictGetNPseudoReconvergenceConss(scip->conflict);
}

/** get number of conflict constraints found so far at the current node */
int SCIPgetNConflictConssFoundNode(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNConflictConssFoundNode", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPconflictGetNConflicts(scip->conflict);
}

/** get total number of conflict constraints added to the problem */
SCIP_Longint SCIPgetNConflictConssApplied(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNConflictConssApplied", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPconflictGetNAppliedConss(scip->conflict);
}

/** gets depth of current node, or -1 if no current node exists; in probing, the current node is the last probing node,
 *  such that the depth includes the probing path
 */
int SCIPgetDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPtreeGetCurrentDepth(scip->tree);
}

/** gets depth of the focus node, or -1 if no focus node exists; the focus node is the currently processed node in the
 *  branching tree, excluding the nodes of the probing path
 */
int SCIPgetFocusDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetFocusDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPtreeGetFocusDepth(scip->tree);
}

/** gets maximal depth of all processed nodes in current branch and bound run (excluding probing nodes) */
int SCIPgetMaxDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetMaxDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->maxdepth;
}

/** gets maximal depth of all processed nodes over all branch and bound runs */
int SCIPgetMaxTotalDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetMaxTotalDepth", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->maxtotaldepth;
}

/** gets total number of backtracks, i.e. number of times, the new node was selected from the leaves queue */
SCIP_Longint SCIPgetNBacktracks(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNBacktracks", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->stat->nbacktracks;
}

/** gets current plunging depth (successive times, a child was selected as next node) */
int SCIPgetPlungeDepth(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPlungeDepth", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->plungedepth;
}

/** gets total number of active constraints at the current node */
int SCIPgetNActiveConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNActiveConss", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nactiveconss;
}

/** gets total number of enabled constraints at the current node */
int SCIPgetNEnabledConss(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNEnabledConss", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   return scip->stat->nenabledconss;
}

/** gets average dual bound of all unprocessed nodes for original problem */
SCIP_Real SCIPgetAvgDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPprobExternObjval(scip->transprob, scip->set,
      SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound));
}

/** gets average lower (dual) bound of all unprocessed nodes in transformed problem */
SCIP_Real SCIPgetAvgLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPtreeGetAvgLowerbound(scip->tree, scip->primal->cutoffbound);
}

/** gets global dual bound */
SCIP_Real SCIPgetDualbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetDualbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return getDualbound(scip);
}

/** gets global lower (dual) bound in transformed problem */
SCIP_Real SCIPgetLowerbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLowerbound", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return getLowerbound(scip);
}

/** gets dual bound of the root node for the original problem */
SCIP_Real SCIPgetDualboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetDualboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPsetIsInfinity(scip->set, scip->stat->rootlowerbound) )
      return getPrimalbound(scip);
   else
      return SCIPprobExternObjval(scip->transprob, scip->set, scip->stat->rootlowerbound);
}

/** gets lower (dual) bound in transformed problem of the root node */
SCIP_Real SCIPgetLowerboundRoot(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetLowerboundRoot", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPnodeGetLowerbound(scip->tree->root);
}

/** gets global primal bound (objective value of best solution or user objective limit) for the original problem */
SCIP_Real SCIPgetPrimalbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPrimalbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return getPrimalbound(scip);
}

/** gets global upper (primal) bound in transformed problem (objective value of best solution or user objective limit) */
SCIP_Real SCIPgetUpperbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetUpperbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return getUpperbound(scip);
}

/** gets global cutoff bound in transformed problem: a sub problem with lower bound larger than the cutoff
 *  cannot contain a better feasible solution; usually, this bound is equal to the upper bound, but if the
 *  objective value is always integral, the cutoff bound is (nearly) one less than the upper bound;
 *  additionally, due to objective function domain propagation, the cutoff bound can be further reduced
 */
SCIP_Real SCIPgetCutoffbound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetCutoffbound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->cutoffbound;
}

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
SCIP_Bool SCIPisPrimalboundSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPisPrimalboundSol", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return SCIPprimalUpperboundIsSol(scip->primal, scip->set, scip->transprob);
}

/** gets current gap |(primalbound - dualbound)/min(|primalbound|,|dualbound|)| if both bounds have same sign,
 *  or infinity, if they have opposite sign
 */
SCIP_Real SCIPgetGap(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real primalbound;
   SCIP_Real dualbound;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetGap", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   if( SCIPsetIsInfinity(scip->set, getLowerbound(scip)) )
   {
      /* in case we could not prove whether the problem is unbounded or infeasible, we want to terminate with
       * gap = +inf instead of gap = 0
       */
      if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD )
         return SCIPsetInfinity(scip->set);
      else
         return 0.0;
   } 
   
   primalbound = getPrimalbound(scip);
   dualbound = getDualbound(scip);

   if( SCIPsetIsEQ(scip->set, primalbound, dualbound) )
      return 0.0;
   else if( SCIPsetIsZero(scip->set, dualbound)
      || SCIPsetIsZero(scip->set, primalbound)
      || SCIPsetIsInfinity(scip->set, REALABS(primalbound))
      || SCIPsetIsInfinity(scip->set, REALABS(dualbound))
      || primalbound * dualbound < 0.0 )
      return SCIPsetInfinity(scip->set);
   else
      return REALABS((primalbound - dualbound)/MIN(REALABS(dualbound),REALABS(primalbound)));
}

/** gets current gap |(upperbound - lowerbound)/min(|upperbound|,|lowerbound|)| in transformed problem if both bounds
 *  have same sign, or infinity, if they have opposite sign
 */
SCIP_Real SCIPgetTransGap(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real upperbound;
   SCIP_Real lowerbound;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetTransGap", FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   upperbound = getUpperbound(scip);
   lowerbound = getLowerbound(scip);

   if( SCIPsetIsInfinity(scip->set, lowerbound) )
      /* in case we could not prove whether the problem is unbounded or infeasible, we want to terminate with
       * gap = +inf instead of gap = 0
       */
      if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD )
         return SCIPsetInfinity(scip->set);
      else
         return 0.0;
   else if( SCIPsetIsEQ(scip->set, upperbound, lowerbound) )
      return 0.0;
   else if( SCIPsetIsZero(scip->set, lowerbound)
      || SCIPsetIsZero(scip->set, upperbound)
      || SCIPsetIsInfinity(scip->set, upperbound)
      || SCIPsetIsInfinity(scip->set, -lowerbound)
      || lowerbound * upperbound < 0.0 )
      return SCIPsetInfinity(scip->set);
   else
      return REALABS((upperbound - lowerbound)/MIN(REALABS(lowerbound),REALABS(upperbound)));
}

/** gets number of feasible primal solutions found so far */
SCIP_Longint SCIPgetNSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNSolsFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->nsolsfound;
}

/** gets number of feasible primal solutions found so far, that improved the primal bound at the time they were found */
SCIP_Longint SCIPgetNBestSolsFound(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNBestSolsFound", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   return scip->primal->nbestsolsfound;
}

/** gets the average pseudo cost value for the given direction over all variables */
SCIP_Real SCIPgetAvgPseudocost(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgPseudocost", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocost(scip->stat->glbhistory, solvaldelta);
}

/** gets the average pseudo cost value for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetAvgPseudocostCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             solvaldelta         /**< difference of variable's new LP value - old LP value */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgPseudocostCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocost(scip->stat->glbhistorycrun, solvaldelta);
}

/** gets the average number of pseudo cost updates for the given direction over all variables */
SCIP_Real SCIPgetAvgPseudocostCount(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgPseudocostCount", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocostCount(scip->stat->glbhistory, dir)
      / MAX(scip->transprob->nbinvars + scip->transprob->nintvars, 1);
}

/** gets the average number of pseudo cost updates for the given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetAvgPseudocostCountCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgPseudocostCountCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetPseudocostCount(scip->stat->glbhistorycrun, dir)
      / MAX(scip->transprob->nbinvars + scip->transprob->nintvars, 1);
}

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5 */
SCIP_Real SCIPgetAvgPseudocostScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgPseudocostScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   pscostdown = SCIPhistoryGetPseudocost(scip->stat->glbhistory, -0.5);
   pscostup = SCIPhistoryGetPseudocost(scip->stat->glbhistory, +0.5);

   return SCIPbranchGetScore(scip->set, NULL, pscostdown, pscostup);
}

/** gets the average pseudo cost score value over all variables, assuming a fractionality of 0.5,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetAvgPseudocostScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real pscostdown;
   SCIP_Real pscostup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgPseudocostScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   pscostdown = SCIPhistoryGetPseudocost(scip->stat->glbhistorycrun, -0.5);
   pscostup = SCIPhistoryGetPseudocost(scip->stat->glbhistorycrun, +0.5);

   return SCIPbranchGetScore(scip->set, NULL, pscostdown, pscostup);
}

/** gets the average conflict score value over all variables */
SCIP_Real SCIPgetAvgConflictScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictscoredown;
   SCIP_Real conflictscoreup;
   SCIP_Real scale;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgConflictScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   scale = scip->transprob->nvars * scip->stat->vsidsweight;
   conflictscoredown = SCIPhistoryGetVSIDS(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) / scale;
   conflictscoreup = SCIPhistoryGetVSIDS(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) / scale;

   return SCIPbranchGetScore(scip->set, NULL, conflictscoredown, conflictscoreup);
}

/** gets the average conflict score value over all variables, only using the pseudo cost information of the current run */
SCIP_Real SCIPgetAvgConflictScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictscoredown;
   SCIP_Real conflictscoreup;
   SCIP_Real scale;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgConflictScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   scale = scip->transprob->nvars * scip->stat->vsidsweight;
   conflictscoredown = SCIPhistoryGetVSIDS(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS) / scale;
   conflictscoreup = SCIPhistoryGetVSIDS(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS) / scale;

   return SCIPbranchGetScore(scip->set, NULL, conflictscoredown, conflictscoreup);
}

/** gets the average inference score value over all variables */
SCIP_Real SCIPgetAvgConflictlengthScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictlengthdown;
   SCIP_Real conflictlengthup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgConflictlengthScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   conflictlengthdown = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   conflictlengthup = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, conflictlengthdown, conflictlengthup);
}

/** gets the average conflictlength score value over all variables, only using the pseudo cost information of the current run */
SCIP_Real SCIPgetAvgConflictlengthScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real conflictlengthdown;
   SCIP_Real conflictlengthup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgConflictlengthScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   conflictlengthdown = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   conflictlengthup = SCIPhistoryGetAvgConflictlength(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, conflictlengthdown, conflictlengthup);
}

/** returns the average number of inferences found after branching in given direction over all variables */
SCIP_Real SCIPgetAvgInferences(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgInferences", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetAvgInferences(scip->stat->glbhistory, dir);
}

/** returns the average number of inferences found after branching in given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetAvgInferencesCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgInferencesCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, dir);
}

/** gets the average inference score value over all variables */
SCIP_Real SCIPgetAvgInferenceScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real inferencesdown;
   SCIP_Real inferencesup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgInferenceScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   inferencesdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   inferencesup = SCIPhistoryGetAvgInferences(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, inferencesdown, inferencesup);
}

/** gets the average inference score value over all variables, only using the inference information information of the
 *  current run
 */
SCIP_Real SCIPgetAvgInferenceScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real inferencesdown;
   SCIP_Real inferencesup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgInferenceScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   inferencesdown = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   inferencesup = SCIPhistoryGetAvgInferences(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, inferencesdown, inferencesup);
}

/** returns the average number of cutoffs found after branching in given direction over all variables */
SCIP_Real SCIPgetAvgCutoffs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgCutoffs", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetAvgCutoffs(scip->stat->glbhistory, dir);
}

/** returns the average number of cutoffs found after branching in given direction over all variables,
 *  only using the pseudo cost information of the current run
 */
SCIP_Real SCIPgetAvgCutoffsCurrentRun(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BRANCHDIR        dir                 /**< branching direction (downwards, or upwards) */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgCutoffsCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPhistoryGetAvgCutoffs(scip->stat->glbhistorycrun, dir);
}

/** gets the average cutoff score value over all variables */
SCIP_Real SCIPgetAvgCutoffScore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real cutoffsdown;
   SCIP_Real cutoffsup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgCutoffScore", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   cutoffsdown = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffsup = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, cutoffsdown, cutoffsup);
}

/** gets the average cutoff score value over all variables, only using the pseudo cost information of the current run */
SCIP_Real SCIPgetAvgCutoffScoreCurrentRun(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real cutoffsdown;
   SCIP_Real cutoffsup;

   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetAvgCutoffScoreCurrentRun", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   cutoffsdown = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_DOWNWARDS);
   cutoffsup = SCIPhistoryGetAvgCutoffs(scip->stat->glbhistorycrun, SCIP_BRANCHDIR_UPWARDS);

   return SCIPbranchGetScore(scip->set, NULL, cutoffsdown, cutoffsup);
}

/** outputs problem to file stream */
static
SCIP_RETCODE printProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROB*            prob,               /**< problem data */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format) */
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RESULT result;
   int i;
   assert(scip != NULL);
   assert(prob != NULL);

   /* try all readers until one could read the file */
   result = SCIP_DIDNOTRUN;
   for( i = 0; i < scip->set->nreaders && result == SCIP_DIDNOTRUN; ++i )
   {
      SCIP_RETCODE retcode;

      if( extension != NULL )
         retcode = SCIPreaderWrite(scip->set->readers[i], prob, scip->set, file, extension, genericnames, &result);
      else
         retcode = SCIPreaderWrite(scip->set->readers[i], prob, scip->set, file, "cip", genericnames, &result);

      /* check for reader errors */
      if( retcode == SCIP_WRITEERROR )
         return retcode;

      SCIP_CALL( retcode );
   }

   switch( result )
   {
   case SCIP_DIDNOTRUN:
      return SCIP_PLUGINNOTFOUND;

   case SCIP_SUCCESS:
      return SCIP_OKAY;

   default:
      assert(i < scip->set->nreaders);
      SCIPerrorMessage("invalid result code <%d> from reader <%s> writing <%s> format\n",
         result, SCIPreaderGetName(scip->set->readers[i]), extension);
      return SCIP_READERROR;
   }  /*lint !e788*/
}

/** outputs original problem to file stream */
SCIP_RETCODE SCIPprintOrigProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( checkStage(scip, "SCIPprintOrigProblem", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   assert(scip != NULL);
   assert( scip->origprob != NULL );

   retcode = printProblem(scip, scip->origprob, file, extension, genericnames);

   /* check for write errors */
   if( retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** outputs transformed problem of the current node to file stream */
SCIP_RETCODE SCIPprintTransProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   const char*           extension,          /**< file format (or NULL for default CIP format)*/
   SCIP_Bool             genericnames        /**< using generic variable and constraint names? */
   )
{
   SCIP_RETCODE retcode;

   SCIP_CALL( checkStage(scip, "SCIPprintTransProblem", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   assert(scip != NULL);
   assert(scip->transprob != NULL );

   retcode = printProblem(scip, scip->transprob, file, extension, genericnames);

   /* check for write errors */
   if( retcode == SCIP_WRITEERROR || retcode == SCIP_PLUGINNOTFOUND )
      return retcode;
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

static
void printPresolverStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPmessageFPrintInfo(file, "Presolvers         :       Time  FixedVars   AggrVars   ChgTypes  ChgBounds   AddHoles    DelCons    AddCons   ChgSides   ChgCoefs\n");

   /* presolver statistics */
   for( i = 0; i < scip->set->npresols; ++i )
   {
      SCIP_PRESOL* presol;
      presol = scip->set->presols[i];
      SCIPmessageFPrintInfo(file, "  %-17.17s:", SCIPpresolGetName(presol));
      SCIPmessageFPrintInfo(file, " %10.2f %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
         SCIPpresolGetTime(presol),
         SCIPpresolGetNFixedVars(presol),
         SCIPpresolGetNAggrVars(presol),
         SCIPpresolGetNChgVarTypes(presol),
         SCIPpresolGetNChgBds(presol),
         SCIPpresolGetNAddHoles(presol),
         SCIPpresolGetNDelConss(presol),
         SCIPpresolGetNAddConss(presol),
         SCIPpresolGetNChgSides(presol),
         SCIPpresolGetNChgCoefs(presol));
   }

   /* presolver statistics */
   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop;
      prop = scip->set->props[i];
      if( SCIPpropDoesPresolve(prop) 
#if 0
         && ( SCIPpropGetNFixedVars(prop) > 0
            || SCIPpropGetNAggrVars(prop) > 0
            || SCIPpropGetNChgVarTypes(prop) > 0
            || SCIPpropGetNChgBds(prop) > 0
            || SCIPpropGetNAddHoles(prop) > 0
            || SCIPpropGetNDelConss(prop) > 0
            || SCIPpropGetNAddConss(prop) > 0
            || SCIPpropGetNChgSides(prop) > 0
            || SCIPpropGetNChgCoefs(prop) > 0
            || SCIPpropGetNUpgdConss(prop) > 0) 
#endif
         )
      {
         SCIPmessageFPrintInfo(file, "  %-17.17s:", SCIPpropGetName(prop));
         SCIPmessageFPrintInfo(file, " %10.2f %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
            SCIPpropGetPresolTime(prop),
            SCIPpropGetNFixedVars(prop),
            SCIPpropGetNAggrVars(prop),
            SCIPpropGetNChgVarTypes(prop),
            SCIPpropGetNChgBds(prop),
            SCIPpropGetNAddHoles(prop),
            SCIPpropGetNDelConss(prop),
            SCIPpropGetNAddConss(prop),
            SCIPpropGetNChgSides(prop),
            SCIPpropGetNChgCoefs(prop));
      }
   }

   /* constraint handler presolving methods statistics */
   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      int maxnactiveconss;

      conshdlr = scip->set->conshdlrs[i];
      maxnactiveconss = SCIPconshdlrGetMaxNActiveConss(conshdlr);
      if( SCIPconshdlrDoesPresolve(conshdlr)
         && (maxnactiveconss > 0 || !SCIPconshdlrNeedsCons(conshdlr)
            || SCIPconshdlrGetNFixedVars(conshdlr) > 0
            || SCIPconshdlrGetNAggrVars(conshdlr) > 0
            || SCIPconshdlrGetNChgVarTypes(conshdlr) > 0
            || SCIPconshdlrGetNChgBds(conshdlr) > 0
            || SCIPconshdlrGetNAddHoles(conshdlr) > 0
            || SCIPconshdlrGetNDelConss(conshdlr) > 0
            || SCIPconshdlrGetNAddConss(conshdlr) > 0
            || SCIPconshdlrGetNChgSides(conshdlr) > 0
            || SCIPconshdlrGetNChgCoefs(conshdlr) > 0
            || SCIPconshdlrGetNUpgdConss(conshdlr) > 0) )
      {
         SCIPmessageFPrintInfo(file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         SCIPmessageFPrintInfo(file, " %10.2f %10d %10d %10d %10d %10d %10d %10d %10d %10d\n",
            SCIPconshdlrGetPresolTime(conshdlr),
            SCIPconshdlrGetNFixedVars(conshdlr),
            SCIPconshdlrGetNAggrVars(conshdlr),
            SCIPconshdlrGetNChgVarTypes(conshdlr),
            SCIPconshdlrGetNChgBds(conshdlr),
            SCIPconshdlrGetNAddHoles(conshdlr),
            SCIPconshdlrGetNDelConss(conshdlr),
            SCIPconshdlrGetNAddConss(conshdlr),
            SCIPconshdlrGetNChgSides(conshdlr),
            SCIPconshdlrGetNChgCoefs(conshdlr));
      }
   }

   /* root node bound changes */
   SCIPmessageFPrintInfo(file, "  root node        :          - %10d          -          - %10d          -          -          -          -          -\n",
      scip->stat->nrootintfixings, scip->stat->nrootboundchgs);
}

/** print constraint statistics to output file */
static
void printConstraintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   /** Add maximal number of constraints of the same type? So far this information is not added because of lack of space. */
   SCIPmessageFPrintInfo(file, "Constraints        :     Number  MaxNumber  #Separate #Propagate    #EnfoLP    #EnfoPS     #Check   #ResProp    Cutoffs    DomReds       Cuts      Conss   Children\n");

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      int startnactiveconss;
      int maxnactiveconss;

      conshdlr = scip->set->conshdlrs[i];
      startnactiveconss = SCIPconshdlrGetStartNActiveConss(conshdlr);
      maxnactiveconss = SCIPconshdlrGetMaxNActiveConss(conshdlr);
      if( maxnactiveconss > 0 || !SCIPconshdlrNeedsCons(conshdlr) )
      {
         SCIPmessageFPrintInfo(file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         SCIPmessageFPrintInfo(file, " %10d%c%10d %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT"\n",
            startnactiveconss,
            maxnactiveconss > startnactiveconss ? '+' : ' ',
            maxnactiveconss,
            SCIPconshdlrGetNSepaCalls(conshdlr),
            SCIPconshdlrGetNPropCalls(conshdlr),
            SCIPconshdlrGetNEnfoLPCalls(conshdlr),
            SCIPconshdlrGetNEnfoPSCalls(conshdlr),
            SCIPconshdlrGetNCheckCalls(conshdlr),
            SCIPconshdlrGetNRespropCalls(conshdlr),
            SCIPconshdlrGetNCutoffs(conshdlr),
            SCIPconshdlrGetNDomredsFound(conshdlr),
            SCIPconshdlrGetNCutsFound(conshdlr),
            SCIPconshdlrGetNConssFound(conshdlr),
            SCIPconshdlrGetNChildren(conshdlr));
      }
   }
}

/** print constraint timing statistics to output file */
static
void printConstraintTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPmessageFPrintInfo(file, "Constraint Timings :  TotalTime   Separate  Propagate     EnfoLP     EnfoPS      Check    ResProp\n");

   for( i = 0; i < scip->set->nconshdlrs; ++i )
   {
      SCIP_CONSHDLR* conshdlr;
      int maxnactiveconss;

      conshdlr = scip->set->conshdlrs[i];
      maxnactiveconss = SCIPconshdlrGetMaxNActiveConss(conshdlr);
      if( maxnactiveconss > 0 || !SCIPconshdlrNeedsCons(conshdlr) )
      {
         SCIP_Real totaltime;

         totaltime = SCIPconshdlrGetSepaTime(conshdlr) + SCIPconshdlrGetPropTime(conshdlr) 
            + SCIPconshdlrGetEnfoLPTime(conshdlr) 
            + SCIPconshdlrGetEnfoPSTime(conshdlr) 
            + SCIPconshdlrGetCheckTime(conshdlr)
            + SCIPconshdlrGetRespropTime(conshdlr);
         
         SCIPmessageFPrintInfo(file, "  %-17.17s:", SCIPconshdlrGetName(conshdlr));
         SCIPmessageFPrintInfo(file, " %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
            totaltime,
            SCIPconshdlrGetSepaTime(conshdlr),
            SCIPconshdlrGetPropTime(conshdlr),
            SCIPconshdlrGetEnfoLPTime(conshdlr),
            SCIPconshdlrGetEnfoPSTime(conshdlr),
            SCIPconshdlrGetCheckTime(conshdlr),
            SCIPconshdlrGetRespropTime(conshdlr));
      }
   }
}

/** print propagator statistics to output file */
static
void printPropagatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPmessageFPrintInfo(file, "Propagators        : #Propagate   #ResProp    Cutoffs    DomReds\n");

   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop;
      prop = scip->set->props[i];

      SCIPmessageFPrintInfo(file, "  %-17.17s: %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT"\n",
         SCIPpropGetName(prop),
         SCIPpropGetNCalls(prop),
         SCIPpropGetNRespropCalls(prop),
         SCIPpropGetNCutoffs(prop),
         SCIPpropGetNDomredsFound(prop));
   }
   
   SCIPmessageFPrintInfo(file, "Propagator Timings :  TotalTime   Presolve  Propagate    ResProp\n");

   for( i = 0; i < scip->set->nprops; ++i )
   {
      SCIP_PROP* prop;
      SCIP_Real totaltime;
         
      prop = scip->set->props[i];
      totaltime = SCIPpropGetPresolTime(prop) + SCIPpropGetTime(prop) + SCIPpropGetRespropTime(prop);
      
      SCIPmessageFPrintInfo(file, "  %-17.17s:", SCIPpropGetName(prop));
      SCIPmessageFPrintInfo(file, " %10.2f %10.2f %10.2f %10.2f\n",
         totaltime, SCIPpropGetPresolTime(prop), SCIPpropGetTime(prop), SCIPpropGetRespropTime(prop));
   }
}

/** print conflict statistic to given output stream */
static
void printConflictStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIPmessageFPrintInfo(file, "Conflict Analysis  :       Time      Calls    Success  Conflicts   Literals    Reconvs ReconvLits   LP Iters\n");
   SCIPmessageFPrintInfo(file, "  propagation      : %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT" %10.1f          -\n",
      SCIPconflictGetPropTime(scip->conflict),
      SCIPconflictGetNPropCalls(scip->conflict),
      SCIPconflictGetNPropSuccess(scip->conflict),
      SCIPconflictGetNPropConflictConss(scip->conflict),
      SCIPconflictGetNPropConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPropConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPropConflictConss(scip->conflict) : 0,
      SCIPconflictGetNPropReconvergenceConss(scip->conflict),
      SCIPconflictGetNPropReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPropReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPropReconvergenceConss(scip->conflict) : 0);
   SCIPmessageFPrintInfo(file, "  infeasible LP    : %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT"\n",
      SCIPconflictGetInfeasibleLPTime(scip->conflict),
      SCIPconflictGetNInfeasibleLPCalls(scip->conflict),
      SCIPconflictGetNInfeasibleLPSuccess(scip->conflict),
      SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict),
      SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNInfeasibleLPConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNInfeasibleLPConflictConss(scip->conflict) : 0,
      SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict),
      SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNInfeasibleLPReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNInfeasibleLPReconvergenceConss(scip->conflict) : 0,
      SCIPconflictGetNInfeasibleLPIterations(scip->conflict));
   SCIPmessageFPrintInfo(file, "  bound exceed. LP : %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT"\n",
      SCIPconflictGetBoundexceedingLPTime(scip->conflict),
      SCIPconflictGetNBoundexceedingLPCalls(scip->conflict),
      SCIPconflictGetNBoundexceedingLPSuccess(scip->conflict),
      SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict),
      SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNBoundexceedingLPConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNBoundexceedingLPConflictConss(scip->conflict) : 0,
      SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict),
      SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNBoundexceedingLPReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNBoundexceedingLPReconvergenceConss(scip->conflict) : 0,
      SCIPconflictGetNBoundexceedingLPIterations(scip->conflict));
   SCIPmessageFPrintInfo(file, "  strong branching : %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT"\n",
      SCIPconflictGetStrongbranchTime(scip->conflict),
      SCIPconflictGetNStrongbranchCalls(scip->conflict),
      SCIPconflictGetNStrongbranchSuccess(scip->conflict),
      SCIPconflictGetNStrongbranchConflictConss(scip->conflict),
      SCIPconflictGetNStrongbranchConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNStrongbranchConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNStrongbranchConflictConss(scip->conflict) : 0,
      SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict),
      SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNStrongbranchReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNStrongbranchReconvergenceConss(scip->conflict) : 0,
      SCIPconflictGetNStrongbranchIterations(scip->conflict));
   SCIPmessageFPrintInfo(file, "  pseudo solution  : %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10.1f %10"SCIP_LONGINT_FORMAT" %10.1f          -\n",
      SCIPconflictGetPseudoTime(scip->conflict),
      SCIPconflictGetNPseudoCalls(scip->conflict),
      SCIPconflictGetNPseudoSuccess(scip->conflict),
      SCIPconflictGetNPseudoConflictConss(scip->conflict),
      SCIPconflictGetNPseudoConflictConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPseudoConflictLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPseudoConflictConss(scip->conflict) : 0,
      SCIPconflictGetNPseudoReconvergenceConss(scip->conflict),
      SCIPconflictGetNPseudoReconvergenceConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNPseudoReconvergenceLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNPseudoReconvergenceConss(scip->conflict) : 0);
   SCIPmessageFPrintInfo(file, "  applied globally :          -          -          - %10"SCIP_LONGINT_FORMAT" %10.1f          -          -          -\n",
      SCIPconflictGetNAppliedGlobalConss(scip->conflict),
      SCIPconflictGetNAppliedGlobalConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNAppliedGlobalLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNAppliedGlobalConss(scip->conflict) : 0);
   SCIPmessageFPrintInfo(file, "  applied locally  :          -          -          - %10"SCIP_LONGINT_FORMAT" %10.1f          -          -          -\n",
      SCIPconflictGetNAppliedLocalConss(scip->conflict),
      SCIPconflictGetNAppliedLocalConss(scip->conflict) > 0
      ? (SCIP_Real)SCIPconflictGetNAppliedLocalLiterals(scip->conflict)
      / (SCIP_Real)SCIPconflictGetNAppliedLocalConss(scip->conflict) : 0);
}

static
void printSeparatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPmessageFPrintInfo(file, "Separators         :       Time      Calls    Cutoffs    DomReds       Cuts      Conss\n");
   SCIPmessageFPrintInfo(file, "  cut pool         : %10.2f %10"SCIP_LONGINT_FORMAT"          -          - %10"SCIP_LONGINT_FORMAT"          -    (maximal pool size: %d)\n",
      SCIPcutpoolGetTime(scip->cutpool),
      SCIPcutpoolGetNCalls(scip->cutpool),
      SCIPcutpoolGetNCutsFound(scip->cutpool),
      SCIPcutpoolGetMaxNCuts(scip->cutpool));

   for( i = 0; i < scip->set->nsepas; ++i )
      SCIPmessageFPrintInfo(file, "  %-17.17s: %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT"\n",
         SCIPsepaGetName(scip->set->sepas[i]),
         SCIPsepaGetTime(scip->set->sepas[i]),
         SCIPsepaGetNCalls(scip->set->sepas[i]),
         SCIPsepaGetNCutoffs(scip->set->sepas[i]),
         SCIPsepaGetNDomredsFound(scip->set->sepas[i]),
         SCIPsepaGetNCutsFound(scip->set->sepas[i]),
         SCIPsepaGetNConssFound(scip->set->sepas[i]));
}

static
void printPricerStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPmessageFPrintInfo(file, "Pricers            :       Time      Calls       Vars\n");
   SCIPmessageFPrintInfo(file, "  problem variables: %10.2f %10d %10d\n",
      SCIPpricestoreGetProbPricingTime(scip->pricestore),
      SCIPpricestoreGetNProbPricings(scip->pricestore),
      SCIPpricestoreGetNProbvarsFound(scip->pricestore));

   for( i = 0; i < scip->set->nactivepricers; ++i )
      SCIPmessageFPrintInfo(file, "  %-17.17s: %10.2f %10d %10d\n",
         SCIPpricerGetName(scip->set->pricers[i]),
         SCIPpricerGetTime(scip->set->pricers[i]),
         SCIPpricerGetNCalls(scip->set->pricers[i]),
         SCIPpricerGetNVarsFound(scip->set->pricers[i]));
}

static
void printBranchruleStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIPmessageFPrintInfo(file, "Branching Rules    :       Time      Calls    Cutoffs    DomReds       Cuts      Conss   Children\n");

   for( i = 0; i < scip->set->nbranchrules; ++i )
      SCIPmessageFPrintInfo(file, "  %-17.17s: %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT"\n",
         SCIPbranchruleGetName(scip->set->branchrules[i]),
         SCIPbranchruleGetTime(scip->set->branchrules[i]),
         SCIPbranchruleGetNLPCalls(scip->set->branchrules[i]) + SCIPbranchruleGetNPseudoCalls(scip->set->branchrules[i]) + SCIPbranchruleGetNExternCalls(scip->set->branchrules[i])
         ,
         SCIPbranchruleGetNCutoffs(scip->set->branchrules[i]),
         SCIPbranchruleGetNDomredsFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNCutsFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNConssFound(scip->set->branchrules[i]),
         SCIPbranchruleGetNChildren(scip->set->branchrules[i]));
}

static
void printHeuristicStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->tree != NULL);

   SCIPmessageFPrintInfo(file, "Primal Heuristics  :       Time      Calls      Found\n");
   SCIPmessageFPrintInfo(file, "  LP solutions     : %10.2f          - %10"SCIP_LONGINT_FORMAT"\n",
      SCIPclockGetTime(scip->stat->lpsoltime),
      scip->stat->nlpsolsfound);
   SCIPmessageFPrintInfo(file, "  pseudo solutions : %10.2f          - %10"SCIP_LONGINT_FORMAT"\n",
      SCIPclockGetTime(scip->stat->pseudosoltime),
      scip->stat->npssolsfound);

   SCIPsetSortHeurs(scip->set);

   for( i = 0; i < scip->set->nheurs; ++i )
      SCIPmessageFPrintInfo(file, "  %-17.17s: %10.2f %10"SCIP_LONGINT_FORMAT" %10"SCIP_LONGINT_FORMAT"\n",
         SCIPheurGetName(scip->set->heurs[i]),
         SCIPheurGetTime(scip->set->heurs[i]),
         SCIPheurGetNCalls(scip->set->heurs[i]),
         SCIPheurGetNSolsFound(scip->set->heurs[i]));
}

static
void printLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->lp != NULL);

   SCIPmessageFPrintInfo(file, "LP                 :       Time      Calls Iterations  Iter/call   Iter/sec  Time-0-It Calls-0-It\n");

   SCIPmessageFPrintInfo(file, "  primal LP        : %10.2f %10d %10"SCIP_LONGINT_FORMAT" %10.2f",
      SCIPclockGetTime(scip->stat->primallptime),
      scip->stat->nprimallps + scip->stat->nprimalzeroitlps,
      scip->stat->nprimallpiterations,
      scip->stat->nprimallps > 0 ? (SCIP_Real)scip->stat->nprimallpiterations/(SCIP_Real)scip->stat->nprimallps : 0.0);
   if( SCIPclockGetTime(scip->stat->primallptime) >= 0.01 )
      SCIPmessageFPrintInfo(file, " %10.2f", (SCIP_Real)scip->stat->nprimallpiterations/SCIPclockGetTime(scip->stat->primallptime));
   else
      SCIPmessageFPrintInfo(file, "          -");
   SCIPmessageFPrintInfo(file, " %10.2f %10d\n",
      scip->stat->primalzeroittime,
      scip->stat->nprimalzeroitlps);

   SCIPmessageFPrintInfo(file, "  dual LP          : %10.2f %10d %10"SCIP_LONGINT_FORMAT" %10.2f",
      SCIPclockGetTime(scip->stat->duallptime),
      scip->stat->nduallps + scip->stat->ndualzeroitlps,
      scip->stat->nduallpiterations,
      scip->stat->nduallps > 0 ? (SCIP_Real)scip->stat->nduallpiterations/(SCIP_Real)scip->stat->nduallps : 0.0);
   if( SCIPclockGetTime(scip->stat->duallptime) >= 0.01 )
      SCIPmessageFPrintInfo(file, " %10.2f", (SCIP_Real)scip->stat->nduallpiterations/SCIPclockGetTime(scip->stat->duallptime));
   else
      SCIPmessageFPrintInfo(file, "          -");
   SCIPmessageFPrintInfo(file, " %10.2f %10d\n",
      scip->stat->dualzeroittime,
      scip->stat->ndualzeroitlps);

   SCIPmessageFPrintInfo(file, "  lex dual LP      : %10.2f %10d %10"SCIP_LONGINT_FORMAT" %10.2f",
      SCIPclockGetTime(scip->stat->lexduallptime),
      scip->stat->nlexduallps,
      scip->stat->nlexduallpiterations,
      scip->stat->nlexduallps > 0 ? (SCIP_Real)scip->stat->nlexduallpiterations/(SCIP_Real)scip->stat->nlexduallps : 0.0);
   if( SCIPclockGetTime(scip->stat->lexduallptime) >= 0.01 )
      SCIPmessageFPrintInfo(file, " %10.2f\n", (SCIP_Real)scip->stat->nlexduallpiterations/SCIPclockGetTime(scip->stat->lexduallptime));
   else
      SCIPmessageFPrintInfo(file, "          -\n");

   SCIPmessageFPrintInfo(file, "  barrier LP       : %10.2f %10d %10"SCIP_LONGINT_FORMAT" %10.2f",
      SCIPclockGetTime(scip->stat->barrierlptime),
      scip->stat->nbarrierlps,
      scip->stat->nbarrierlpiterations,
      scip->stat->nbarrierlps > 0 ? (SCIP_Real)scip->stat->nbarrierlpiterations/(SCIP_Real)scip->stat->nbarrierlps : 0.0);
   if( SCIPclockGetTime(scip->stat->barrierlptime) >= 0.01 )
      SCIPmessageFPrintInfo(file, " %10.2f", (SCIP_Real)scip->stat->nbarrierlpiterations/SCIPclockGetTime(scip->stat->barrierlptime));
   else
      SCIPmessageFPrintInfo(file, "          -");
   SCIPmessageFPrintInfo(file, " %10.2f %10d\n",
      scip->stat->barrierzeroittime,
      scip->stat->nbarrierzeroitlps);

   SCIPmessageFPrintInfo(file, "  diving/probing LP: %10.2f %10d %10"SCIP_LONGINT_FORMAT" %10.2f",
      SCIPclockGetTime(scip->stat->divinglptime),
      scip->stat->ndivinglps,
      scip->stat->ndivinglpiterations,
      scip->stat->ndivinglps > 0 ? (SCIP_Real)scip->stat->ndivinglpiterations/(SCIP_Real)scip->stat->ndivinglps : 0.0);
   if( SCIPclockGetTime(scip->stat->divinglptime) >= 0.01 )
      SCIPmessageFPrintInfo(file, " %10.2f\n", (SCIP_Real)scip->stat->ndivinglpiterations/SCIPclockGetTime(scip->stat->divinglptime));
   else
      SCIPmessageFPrintInfo(file, "          -\n");

   SCIPmessageFPrintInfo(file, "  strong branching : %10.2f %10d %10"SCIP_LONGINT_FORMAT" %10.2f",
      SCIPclockGetTime(scip->stat->strongbranchtime),
      scip->stat->nstrongbranchs,
      scip->stat->nsblpiterations,
      scip->stat->nstrongbranchs > 0 ? (SCIP_Real)scip->stat->nsblpiterations/(SCIP_Real)scip->stat->nstrongbranchs : 0.0);
   if( SCIPclockGetTime(scip->stat->strongbranchtime) >= 0.01 )
      SCIPmessageFPrintInfo(file, " %10.2f\n", (SCIP_Real)scip->stat->nsblpiterations/SCIPclockGetTime(scip->stat->strongbranchtime));
   else
      SCIPmessageFPrintInfo(file, "          -\n");

   SCIPmessageFPrintInfo(file, "    (at root node) :          - %10d %10"SCIP_LONGINT_FORMAT" %10.2f          -\n",
      scip->stat->nrootstrongbranchs,
      scip->stat->nrootsblpiterations,
      scip->stat->nrootstrongbranchs > 0
      ? (SCIP_Real)scip->stat->nrootsblpiterations/(SCIP_Real)scip->stat->nrootstrongbranchs : 0.0);

   SCIPmessageFPrintInfo(file, "  conflict analysis: %10.2f %10d %10"SCIP_LONGINT_FORMAT" %10.2f",
      SCIPclockGetTime(scip->stat->conflictlptime),
      scip->stat->nconflictlps,
      scip->stat->nconflictlpiterations,
      scip->stat->nconflictlps > 0 ? (SCIP_Real)scip->stat->nconflictlpiterations/(SCIP_Real)scip->stat->nconflictlps : 0.0);
   if( SCIPclockGetTime(scip->stat->conflictlptime) >= 0.01 )
      SCIPmessageFPrintInfo(file, " %10.2f\n", (SCIP_Real)scip->stat->nconflictlpiterations/SCIPclockGetTime(scip->stat->conflictlptime));
   else
      SCIPmessageFPrintInfo(file, "          -\n");
   
}

static
void printNLPStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);

   if( scip->nlp == NULL )
      return;

   SCIPmessageFPrintInfo(file, "NLP                :       Time      Calls\n");

   SCIPmessageFPrintInfo(file, "  all NLPs         : %10.2f %10d\n",
      SCIPclockGetTime(scip->stat->nlpsoltime),
      scip->stat->nnlps);
}

static
void printRelaxatorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   if( scip->set->nrelaxs == 0 )
      return;

   SCIPmessageFPrintInfo(file, "Relaxators         :       Time      Calls\n");

   for( i = 0; i < scip->set->nrelaxs; ++i )
      SCIPmessageFPrintInfo(file, "  %-17.17s: %10.2f %10"SCIP_LONGINT_FORMAT"\n",
         SCIPrelaxGetName(scip->set->relaxs[i]),
         SCIPrelaxGetTime(scip->set->relaxs[i]),
         SCIPrelaxGetNCalls(scip->set->relaxs[i]));
}

static
void printTreeStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->tree != NULL);

   SCIPmessageFPrintInfo(file, "B&B Tree           :\n");
   SCIPmessageFPrintInfo(file, "  number of runs   : %10d\n", scip->stat->nruns);
   SCIPmessageFPrintInfo(file, "  nodes            : %10"SCIP_LONGINT_FORMAT"\n", scip->stat->nnodes);
   SCIPmessageFPrintInfo(file, "  nodes (total)    : %10"SCIP_LONGINT_FORMAT"\n", scip->stat->ntotalnodes);
   SCIPmessageFPrintInfo(file, "  nodes left       : %10d\n", SCIPtreeGetNNodes(scip->tree));
   SCIPmessageFPrintInfo(file, "  max depth        : %10d\n", scip->stat->maxdepth);
   SCIPmessageFPrintInfo(file, "  max depth (total): %10d\n", scip->stat->maxtotaldepth);
   SCIPmessageFPrintInfo(file, "  backtracks       : %10"SCIP_LONGINT_FORMAT" (%.1f%%)\n", scip->stat->nbacktracks,
      scip->stat->nnodes > 0 ? 100.0 * (SCIP_Real)scip->stat->nbacktracks / (SCIP_Real)scip->stat->nnodes : 0.0);
   SCIPmessageFPrintInfo(file, "  delayed cutoffs  : %10"SCIP_LONGINT_FORMAT"\n", scip->stat->ndelayedcutoffs);
   SCIPmessageFPrintInfo(file, "  repropagations   : %10"SCIP_LONGINT_FORMAT" (%"SCIP_LONGINT_FORMAT" domain reductions, %"SCIP_LONGINT_FORMAT" cutoffs)\n",
      scip->stat->nreprops, scip->stat->nrepropboundchgs, scip->stat->nrepropcutoffs);
   SCIPmessageFPrintInfo(file, "  avg switch length: %10.2f\n",
      scip->stat->nnodes > 0
      ? (SCIP_Real)(scip->stat->nactivatednodes + scip->stat->ndeactivatednodes) / (SCIP_Real)scip->stat->nnodes : 0.0);
   SCIPmessageFPrintInfo(file, "  switching time   : %10.2f\n", SCIPclockGetTime(scip->stat->nodeactivationtime));
}

/** display solution statistics */
static
void printSolutionStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_Real primalbound;
   SCIP_Real dualbound;
   SCIP_Real dualboundroot;
   SCIP_Real bestsol;
   SCIP_Real gap;
   SCIP_Real firstprimalbound;

   assert(scip != NULL);
   assert(scip->stat != NULL);
   assert(scip->primal != NULL);

   primalbound = getPrimalbound(scip);
   dualbound = getDualbound(scip);
   dualboundroot = SCIPgetDualboundRoot(scip);
   gap = SCIPgetGap(scip);

   SCIPmessageFPrintInfo(file, "Solution           :\n");
   SCIPmessageFPrintInfo(file, "  Solutions found  : %10"SCIP_LONGINT_FORMAT" (%d improvements)\n",
      scip->primal->nsolsfound, scip->primal->nbestsolsfound);

   if( SCIPsetIsInfinity(scip->set, REALABS(primalbound)) )
   {
      if( scip->set->stage == SCIP_STAGE_SOLVED )
      {
         if( scip->primal->nsols == 0 )
         {
            if( SCIPgetStatus(scip) == SCIP_STATUS_INFORUNBD )
               SCIPmessageFPrintInfo(file, "  Primal Bound     : infeasible or unbounded\n");
            else
            {
               assert(SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE);
               SCIPmessageFPrintInfo(file, "  Primal Bound     : infeasible\n");
            }
         }
         else
         {
            assert(SCIPgetStatus(scip) == SCIP_STATUS_UNBOUNDED);
            SCIPmessageFPrintInfo(file, "  Primal Bound     :  unbounded\n");
         }
      }
      else
         SCIPmessageFPrintInfo(file, "  Primal Bound     :          -\n");
   }
   else
   {
      if( scip->primal->nsols == 0 )
      {
         SCIPmessageFPrintInfo(file, "  Primal Bound     : %+21.14e", primalbound);
         SCIPmessageFPrintInfo(file, "   (user objective limit)\n");
      }
      else
      {
         /* display first primal bound line */
         firstprimalbound = scip->stat->firstprimalbound;
         SCIPmessageFPrintInfo(file, "  First Solution   : %+21.14e", firstprimalbound);
         
         SCIPmessageFPrintInfo(file, "   (in run %d, after %"SCIP_LONGINT_FORMAT" nodes, %.2f seconds, depth %d, found by <%s>)\n",
            scip->stat->nrunsbeforefirst,
            scip->stat->nnodesbeforefirst,
            scip->stat->firstprimaltime,
            scip->stat->firstprimaldepth,
            ( scip->stat->firstprimalheur != NULL )
            ? ( SCIPheurGetName(scip->stat->firstprimalheur) )
            : (( scip->stat->nrunsbeforefirst == 0 ) ? "initial" : "relaxation"));
         
         SCIPmessageFPrintInfo(file, "  Primal Bound     : %+21.14e", primalbound);
         
         /* display (best) primal bound */
         bestsol = SCIPsolGetObj(scip->primal->sols[0], scip->set, scip->transprob);
         bestsol = SCIPretransformObj(scip, bestsol);
         if( SCIPsetIsGT(scip->set, bestsol, primalbound) )
         {
            SCIPmessageFPrintInfo(file, "   (user objective limit)\n");
            SCIPmessageFPrintInfo(file, "  Best Solution    : %+21.14e", bestsol);
         }
         SCIPmessageFPrintInfo(file, "   (in run %d, after %"SCIP_LONGINT_FORMAT" nodes, %.2f seconds, depth %d, found by <%s>)\n",
            SCIPsolGetRunnum(scip->primal->sols[0]),
            SCIPsolGetNodenum(scip->primal->sols[0]),
            SCIPsolGetTime(scip->primal->sols[0]),
            SCIPsolGetDepth(scip->primal->sols[0]),
            SCIPsolGetHeur(scip->primal->sols[0]) != NULL
            ? SCIPheurGetName(SCIPsolGetHeur(scip->primal->sols[0]))
            : (SCIPsolGetRunnum(scip->primal->sols[0]) == 0 ? "initial" : "relaxation"));
      }
   }
   if( SCIPsetIsInfinity(scip->set, REALABS(dualbound)) )
      SCIPmessageFPrintInfo(file, "  Dual Bound       :          -\n");
   else
      SCIPmessageFPrintInfo(file, "  Dual Bound       : %+21.14e\n", dualbound);
   if( SCIPsetIsInfinity(scip->set, gap) )
      SCIPmessageFPrintInfo(file, "  Gap              :   infinite\n");
   else
      SCIPmessageFPrintInfo(file, "  Gap              : %10.2f %%\n", 100.0 * gap);
   if( SCIPsetIsInfinity(scip->set, REALABS(dualboundroot)) )
      SCIPmessageFPrintInfo(file, "  Root Dual Bound  :          -\n");
   else
      SCIPmessageFPrintInfo(file, "  Root Dual Bound  : %+21.14e\n", dualboundroot);

   SCIPmessageFPrintInfo(file, "  Root Iterations  : %10"SCIP_LONGINT_FORMAT"\n", scip->stat->nrootlpiterations);
}
      
/** display timing statistics */
static
void printTimingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file */
   )
{
   SCIP_Real readingtime;

   assert(SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM);
   
   readingtime = SCIPgetReadingTime(scip);
   
   if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
   {
      SCIPmessageFPrintInfo(file, "Total Time         : %10.2f\n", readingtime);
      SCIPmessageFPrintInfo(file, "  reading          : %10.2f\n", readingtime);
   }
   else
   {
      SCIP_Real totaltime;
      SCIP_Real solvingtime;
      
      solvingtime  = SCIPclockGetTime(scip->stat->solvingtime);

      if( scip->set->time_reading ) 
         totaltime = solvingtime;
      else
         totaltime = solvingtime + readingtime;
   
      SCIPmessageFPrintInfo(file, "Total Time         : %10.2f\n", totaltime);
      SCIPmessageFPrintInfo(file, "  solving          : %10.2f\n", solvingtime);
      SCIPmessageFPrintInfo(file, "  presolving       : %10.2f (included in solving)\n", SCIPclockGetTime(scip->stat->presolvingtime));
      SCIPmessageFPrintInfo(file, "  reading          : %10.2f%s\n", readingtime, scip->set->time_reading ? " (included in solving)" : "");
   }
}

/** outputs solving statistics */
SCIP_RETCODE SCIPprintStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPprintStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   SCIPmessageFPrintInfo(file, "SCIP Status        : ");
   SCIP_CALL( SCIPprintStage(scip, file) );
   SCIPmessageFPrintInfo(file, "\n");

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
      SCIPmessageFPrintInfo(file, "Original Problem   : no problem exists.\n");
      return SCIP_OKAY;

   case SCIP_STAGE_PROBLEM:
   {
      printTimingStatistics(scip, file);
      SCIPmessageFPrintInfo(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      return SCIP_OKAY;
   }
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   {
      printTimingStatistics(scip, file);
      SCIPmessageFPrintInfo(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      SCIPmessageFPrintInfo(file, "Presolved Problem  :\n");
      SCIPprobPrintStatistics(scip->transprob, file);
      printPresolverStatistics(scip, file);
      printConstraintStatistics(scip, file);
      printConstraintTimingStatistics(scip, file);
      printPropagatorStatistics(scip, file);
      printConflictStatistics(scip, file);
      return SCIP_OKAY;
   }
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   {
      printTimingStatistics(scip, file);
      SCIPmessageFPrintInfo(file, "Original Problem   :\n");
      SCIPprobPrintStatistics(scip->origprob, file);
      SCIPmessageFPrintInfo(file, "Presolved Problem  :\n");
      SCIPprobPrintStatistics(scip->transprob, file);
      printPresolverStatistics(scip, file);
      printConstraintStatistics(scip, file);
      printConstraintTimingStatistics(scip, file);
      printPropagatorStatistics(scip, file);
      printConflictStatistics(scip, file);
      printSeparatorStatistics(scip, file);
      printPricerStatistics(scip, file);
      printBranchruleStatistics(scip, file);
      printHeuristicStatistics(scip, file);
      printLPStatistics(scip, file);
      printNLPStatistics(scip, file);
      printRelaxatorStatistics(scip, file);
      printTreeStatistics(scip, file);
      printSolutionStatistics(scip, file);
      return SCIP_OKAY;
   }
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** outputs history statistics about branchings on variables */
SCIP_RETCODE SCIPprintBranchingStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** vars;
   int totalnstrongbranchs;
   int v;

   SCIP_CALL( checkStage(scip, "SCIPprintBranchingStatistics", TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
      SCIPmessageFPrintInfo(file, "problem not yet solved. branching statistics not available.\n");
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, scip->transprob->nvars) );
      for( v = 0; v < scip->transprob->nvars; ++v )
      {
         SCIP_VAR* var;
         int i;

         var = scip->transprob->vars[v];
         for( i = v; i > 0 && strcmp(SCIPvarGetName(var), SCIPvarGetName(vars[i-1])) < 0; i-- )
            vars[i] = vars[i-1];
         vars[i] = var;
      }

      SCIPmessageFPrintInfo(file, "                                      locks              branchings              inferences      cutoffs            LP gain       pscostcount\n");
      SCIPmessageFPrintInfo(file, "variable          prio   factor   down     up  depth    down      up    sb     down       up   down     up      down        up    down      up\n");

      totalnstrongbranchs = 0;
      for( v = 0; v < scip->transprob->nvars; ++v )
      {
         if( SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS) > 0
            || SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS) > 0
            || SCIPgetVarNStrongbranchs(scip, vars[v]) > 0 )
         {
            int nstrongbranchs;

            nstrongbranchs = SCIPgetVarNStrongbranchs(scip, vars[v]);
            totalnstrongbranchs += nstrongbranchs;
            SCIPmessageFPrintInfo(file, "%-16s %5d %8.1f %6d %6d %6.1f %7"SCIP_LONGINT_FORMAT" %7"SCIP_LONGINT_FORMAT" %5d %8.1f %8.1f %5.1f%% %5.1f%% %9.4f %9.4f %7.1f %7.1f\n",
               SCIPvarGetName(vars[v]),
               SCIPvarGetBranchPriority(vars[v]),
               SCIPvarGetBranchFactor(vars[v]),
               SCIPvarGetNLocksDown(vars[v]),
               SCIPvarGetNLocksUp(vars[v]),
               (SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_DOWNWARDS)
                  + SCIPvarGetAvgBranchdepth(vars[v], SCIP_BRANCHDIR_UPWARDS))/2.0 - 1.0,
               SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetNBranchings(vars[v], SCIP_BRANCHDIR_UPWARDS),
               nstrongbranchs,
               SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetAvgInferences(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS),
               100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_DOWNWARDS),
               100.0 * SCIPvarGetAvgCutoffs(vars[v], scip->stat, SCIP_BRANCHDIR_UPWARDS),
               SCIPvarGetPseudocost(vars[v], scip->stat, -1.0),
               SCIPvarGetPseudocost(vars[v], scip->stat, +1.0),
               SCIPvarGetPseudocostCount(vars[v], SCIP_BRANCHDIR_DOWNWARDS),
               SCIPvarGetPseudocostCount(vars[v], SCIP_BRANCHDIR_UPWARDS));
         }
      }
      SCIPmessageFPrintInfo(file, "total                                                %7"SCIP_LONGINT_FORMAT" %7"SCIP_LONGINT_FORMAT" %5d %8.1f %8.1f %5.1f%% %5.1f%% %9.4f %9.4f %7.1f %7.1f\n",
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS),
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS),
         totalnstrongbranchs,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) > 0
         ? SCIPhistoryGetInferenceSum(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) : 0.0,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) > 0
         ? SCIPhistoryGetInferenceSum(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) : 0.0,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) > 0
         ? SCIPhistoryGetCutoffSum(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS) : 0.0,
         SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) > 0
         ? SCIPhistoryGetCutoffSum(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS)
         / (SCIP_Real)SCIPhistoryGetNBranchings(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS) : 0.0,
         SCIPhistoryGetPseudocost(scip->stat->glbhistory, -1.0),
         SCIPhistoryGetPseudocost(scip->stat->glbhistory, +1.0),
         SCIPhistoryGetPseudocostCount(scip->stat->glbhistory, SCIP_BRANCHDIR_DOWNWARDS),
         SCIPhistoryGetPseudocostCount(scip->stat->glbhistory, SCIP_BRANCHDIR_UPWARDS));

      SCIPfreeBufferArray(scip, &vars);

      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** outputs node information display line */
SCIP_RETCODE SCIPprintDisplayLine(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_VERBLEVEL        verblevel           /**< minimal verbosity level to actually display the information line */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPprintDisplayLine", FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) );

   if( (SCIP_VERBLEVEL)scip->set->disp_verblevel >= verblevel )
   {
      SCIP_CALL( SCIPdispPrintLine(scip->set, scip->stat, file, TRUE) );
   }

   return SCIP_OKAY;
}

/** gets total number of implications between variables that are stored in the implication graph */
int SCIPgetNImplications(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetNImplications", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return scip->stat->nimplications;
}

/** stores conflict graph of binary variables' implications into a file, which can be used as input for the DOT tool */
SCIP_RETCODE SCIPwriteImplicationConflictGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< file name, or NULL for stdout */
   )
{
   FILE* file;
   SCIP_VAR** vars;
   int nvars;
   int v;

   SCIP_CALL( checkStage(scip, "SCIPwriteImplicationConflictGraph", FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   /* open file */
   if( filename == NULL )
      file = NULL;
   else
   {
      file = fopen(filename, "w");
      if( file == NULL )
      {
         SCIPerrorMessage("cannot create file <%s>\n", filename);
         SCIPprintSysError(filename);
         return SCIP_FILECREATEERROR;
      }
   }

   vars = scip->transprob->vars;
   nvars = scip->transprob->nbinvars;

   /* write header */
   SCIPmessageFPrintInfo(file, "digraph implconfgraph {\n");

   /* store nodes */
   for( v = 0; v < nvars; ++v )
   {
      if( SCIPvarGetNImpls(vars[v], TRUE) > 0 )
         SCIPmessageFPrintInfo(file, "pos%d [label=\"%s\"];\n", v, SCIPvarGetName(vars[v]));
      if( SCIPvarGetNImpls(vars[v], FALSE) > 0 )
         SCIPmessageFPrintInfo(file, "neg%d [style=filled,fillcolor=red,label=\"%s\"];\n", v, SCIPvarGetName(vars[v]));
      if( SCIPvarGetNImpls(vars[v], TRUE) > 0 && SCIPvarGetNImpls(vars[v], FALSE) > 0 )
         SCIPmessageFPrintInfo(file, "pos%d -> neg%d [dir=both];\n", v, v);
   }

   /* store edges */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Bool fix;

      fix = FALSE;
      do
      {
         SCIP_VAR** implvars;
         SCIP_BOUNDTYPE* impltypes;
         int nimpls;
         int i;

         nimpls = SCIPvarGetNBinImpls(vars[v], fix);
         implvars = SCIPvarGetImplVars(vars[v], fix);
         impltypes = SCIPvarGetImplTypes(vars[v], fix);
         for( i = 0; i < nimpls; ++i )
         {
            int implidx;

            implidx = SCIPvarGetProbindex(implvars[i]);
            if( implidx > v )
               SCIPmessageFPrintInfo(file, "%s%d -> %s%d [dir=none];\n", fix == TRUE ? "pos" : "neg", v,
                  impltypes[i] == SCIP_BOUNDTYPE_UPPER ? "pos" : "neg", implidx);
         }
         fix = !fix;
      }
      while( fix == TRUE );
   }

   /* write footer */
   SCIPmessageFPrintInfo(file, "}\n");

   /* close file */
   if( filename != NULL )
   {
      assert(file != NULL);
      fclose(file);
   }

   return SCIP_OKAY;
}




/*
 * timing methods
 */

/** gets current time of day in seconds (standard time zone) */
SCIP_Real SCIPgetTimeOfDay(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetTimeOfDay", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTimeOfDay();
}

/** creates a clock using the default clock type */
SCIP_RETCODE SCIPcreateClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPclockCreate(clck, SCIP_CLOCKTYPE_DEFAULT) );

   return SCIP_OKAY;
}

/** creates a clock counting the CPU user seconds */
SCIP_RETCODE SCIPcreateCPUClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateCPUClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPclockCreate(clck, SCIP_CLOCKTYPE_CPU) );

   return SCIP_OKAY;
}

/** creates a clock counting the wall clock seconds */
SCIP_RETCODE SCIPcreateWallClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateWallClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPclockCreate(clck, SCIP_CLOCKTYPE_WALL) );

   return SCIP_OKAY;
}

/** frees a clock */
SCIP_RETCODE SCIPfreeClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK**          clck                /**< pointer to clock timer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockFree(clck);

   return SCIP_OKAY;
}

/** resets the time measurement of a clock to zero and completely stops the clock */
SCIP_RETCODE SCIPresetClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPresetClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockReset(clck);

   return SCIP_OKAY;
}

/** starts the time measurement of a clock */
SCIP_RETCODE SCIPstartClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstartClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStart(clck, scip->set);

   return SCIP_OKAY;
}

/** stops the time measurement of a clock */
SCIP_RETCODE SCIPstopClock(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstopClock", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStop(clck, scip->set);

   return SCIP_OKAY;
}

/** starts the current solving time */
SCIP_RETCODE SCIPstartSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstartSolvingTime", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStart(scip->stat->solvingtime, scip->set);

   return SCIP_OKAY;
}

/** stops the current solving time in seconds */
SCIP_RETCODE SCIPstopSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPstopSolvingTime", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockStop(scip->stat->solvingtime, scip->set);

   return SCIP_OKAY;
}

/** gets the measured time of a clock in seconds */
SCIP_Real SCIPgetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck                /**< clock timer */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetClockTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTime(clck);
}

/** sets the measured time of a clock to the given value in seconds */
SCIP_RETCODE SCIPsetClockTime(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CLOCK*           clck,               /**< clock timer */
   SCIP_Real             sec                 /**< time in seconds to set the clock's timer to */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetClockTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPclockSetTime(clck, sec);

   return SCIP_OKAY;
}

/** gets the current total SCIP time in seconds */
SCIP_Real SCIPgetTotalTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetTotalTime", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPclockGetTime(scip->totaltime);
}

/** gets the current solving time in seconds */
SCIP_Real SCIPgetSolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetSolvingTime", FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPclockGetTime(scip->stat->solvingtime);
}

/** gets the current reading time in seconds */
SCIP_Real SCIPgetReadingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real readingtime;
   int r;
   
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetReadingTime", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   readingtime = 0.0;
   
   /* sum up the reading time of all readers */
   for( r = 0; r < scip->set->nreaders; ++r )
   {
      assert(scip->set->readers[r] != NULL);
      assert(!SCIPisNegative(scip, SCIPreaderGetReadingTime(scip->set->readers[r])));
      readingtime += SCIPreaderGetReadingTime(scip->set->readers[r]);
   }
   
   return readingtime;
}

/** gets the current presolving time in seconds */
SCIP_Real SCIPgetPresolvingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPresolvingTime", FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   return SCIPclockGetTime(scip->stat->presolvingtime);
}




/*
 * numeric values and comparisons
 */

/** returns value treated as infinity */
SCIP_Real SCIPinfinity(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetInfinity(scip->set);
}

/** returns value treated as zero */
SCIP_Real SCIPepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetEpsilon(scip->set);
}

/** returns value treated as zero for sums of floating point values */
SCIP_Real SCIPsumepsilon(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetSumepsilon(scip->set);
}

/** returns feasibility tolerance for constraints */
SCIP_Real SCIPfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeastol(scip->set);
}

/** returns feasibility tolerance for reduced costs */
SCIP_Real SCIPdualfeastol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetDualfeastol(scip->set);
}

/** returns convergence tolerance used in barrier algorithm */
SCIP_Real SCIPbarrierconvtol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetBarrierconvtol(scip->set);
}

/** sets the feasibility tolerance for constraints */
SCIP_RETCODE SCIPchgFeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             feastol             /**< new feasibility tolerance for constraints */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgFeastol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* mark the LP unsolved, if the feasibility tolerance was tightened */
   if( scip->lp != NULL && feastol < SCIPsetFeastol(scip->set) )
      scip->lp->solved = FALSE;

   /* change the settings */
   SCIP_CALL( SCIPsetSetFeastol(scip->set, feastol) );

   return SCIP_OKAY;
}

/** sets the feasibility tolerance for reduced costs */
SCIP_RETCODE SCIPchgDualfeastol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             dualfeastol         /**< new feasibility tolerance for reduced costs */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgDualfeastol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* mark the LP unsolved, if the dual feasibility tolerance was tightened */
   if( scip->lp != NULL && dualfeastol < SCIPsetDualfeastol(scip->set) )
      scip->lp->solved = FALSE;

   /* change the settings */
   SCIP_CALL( SCIPsetSetDualfeastol(scip->set, dualfeastol) );

   return SCIP_OKAY;
}

/** sets the convergence tolerance used in barrier algorithm */
SCIP_RETCODE SCIPchgBarrierconvtol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             barrierconvtol      /**< new convergence tolerance used in barrier algorithm */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPchgBarrierconvtol", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* mark the LP unsolved, if the convergence tolerance was tightened, and the LP was solved with the barrier algorithm */
   if( scip->lp != NULL && barrierconvtol < SCIPsetBarrierconvtol(scip->set)
      && (scip->lp->lastlpalgo == SCIP_LPALGO_BARRIER || scip->lp->lastlpalgo == SCIP_LPALGO_BARRIERCROSSOVER) )
      scip->lp->solved = FALSE;

   /* change the settings */
   SCIP_CALL( SCIPsetSetBarrierconvtol(scip->set, barrierconvtol) );

   return SCIP_OKAY;
}

/** marks that some limit parameter was changed */
void SCIPmarkLimitChanged(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPmarkLimitChanged", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   /* change the settings */
   SCIPsetSetLimitChanged(scip->set);
}

/** outputs a real number, or "+infinity", or "-infinity" to a file */
void SCIPprintReal(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file (or NULL for standard output) */
   SCIP_Real             val,                /**< value to print */
   int                   width,              /**< width of the field */
   int                   precision           /**< number of significant digits printed */
   )
{
   char s[SCIP_MAXSTRLEN];
   char strformat[SCIP_MAXSTRLEN];

   assert(scip != NULL);

   if( SCIPsetIsInfinity(scip->set, val) )
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "+infinity");
   else if( SCIPsetIsInfinity(scip->set, -val) )
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "-infinity");
   else
   {
      (void) SCIPsnprintf(strformat, SCIP_MAXSTRLEN, "%%.%dg", precision);
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, (const char*)strformat, val);
   }
   (void) SCIPsnprintf(strformat, SCIP_MAXSTRLEN, "%%%ds", width);
   SCIPmessageFPrintInfo(file, (const char*)strformat, s);
}




/*
 * memory management
 */

/** returns block memory to use at the current time */
BMS_BLKMEM* SCIPblkmem(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPblkmem", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   assert(scip != NULL);
   assert(scip->set != NULL);
   assert(scip->mem != NULL);

   switch( scip->set->stage )
   {
   case SCIP_STAGE_INIT:
   case SCIP_STAGE_PROBLEM:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_INITSOLVE:
   case SCIP_STAGE_SOLVING:
   case SCIP_STAGE_SOLVED:
   case SCIP_STAGE_FREESOLVE:
   case SCIP_STAGE_FREETRANS:
      return scip->mem->probmem;
      
   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return NULL;
   }  /*lint !e788*/
}

/** returns the total number of bytes used in block memory */
SCIP_Longint SCIPgetMemUsed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetMemUsed", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPmemGetUsed(scip->mem);
}

/** calculate memory size for dynamically allocated arrays */
int SCIPcalcMemGrowSize(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(scip != NULL);

   return SCIPsetCalcMemGrowSize(scip->set, num);
}

/** extends a dynamically allocated block memory array to be able to store at least the given number of elements;
 *  use SCIPensureBlockMemoryArray() define to call this method!
 */
SCIP_RETCODE SCIPensureBlockMemoryArray_call(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                arrayptr,           /**< pointer to dynamically sized array */
   size_t                elemsize,           /**< size in bytes of each element in array */
   int*                  arraysize,          /**< pointer to current array size */
   int                   minsize             /**< required minimal array size */
   )
{
   assert(scip != NULL);
   assert(arrayptr != NULL);
   assert(elemsize > 0);
   assert(arraysize != NULL);

   if( minsize > *arraysize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(scip->set, minsize);
      SCIP_ALLOC( BMSreallocBlockMemorySize(SCIPblkmem(scip), arrayptr, *arraysize * elemsize, newsize * elemsize) );
      *arraysize = newsize;
   }

   return SCIP_OKAY;
}

/** gets a memory buffer with at least the given size */
SCIP_RETCODE SCIPallocBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to store the buffer */
   int                   size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   SCIP_CALL( checkStage(scip, "SCIPallocBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetAllocBufferSize(scip->set, ptr, size) );

   return SCIP_OKAY;
}

/** allocates a memory buffer with at least the given size and copies the given memory into the buffer */
SCIP_RETCODE SCIPduplicateBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to store the buffer */
   const void*           source,             /**< memory block to copy into the buffer */
   int                   size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   SCIP_CALL( checkStage(scip, "SCIPduplicateBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetDuplicateBufferSize(scip->set, ptr, source, size) );

   return SCIP_OKAY;
}

/** reallocates a memory buffer to at least the given size */
SCIP_RETCODE SCIPreallocBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to the buffer */
   int                   size                /**< required size in bytes of buffer */
   )
{
   assert(ptr != NULL);

   SCIP_CALL( checkStage(scip, "SCIPreallocBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPsetReallocBufferSize(scip->set, ptr, size) );

   return SCIP_OKAY;
}

/** frees a memory buffer */
void SCIPfreeBufferSize(
   SCIP*                 scip,               /**< SCIP data structure */
   void**                ptr,                /**< pointer to the buffer */
   int                   dummysize           /**< used to get a safer define for SCIPfreeBuffer() and SCIPfreeBufferArray() */
   )
{  /*lint --e{715}*/
   assert(ptr != NULL);
   assert(dummysize == 0);

   SCIP_CALL_ABORT( checkStage(scip, "SCIPfreeBufferSize", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIPsetFreeBufferSize(scip->set, ptr);
}

/** prints output about used memory */
void SCIPprintMemoryDiagnostic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(scip->set != NULL);

   BMSdisplayMemory();

   SCIPmessagePrintInfo("\nParameter Block Memory (%p):\n", scip->mem->setmem);
   BMSdisplayBlockMemory(scip->mem->setmem);

   SCIPmessagePrintInfo("\nSolution Block Memory (%p):\n", scip->mem->probmem);
   BMSdisplayBlockMemory(scip->mem->probmem);

   SCIPmessagePrintInfo("\nMemory Buffers:\n");
   SCIPbufferPrint(scip->set->buffer);
}




/*
 * simple functions implemented as defines
 */

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 * In optimized mode, the methods are implemented as defines to improve performance.
 * However, we want to have them in the library anyways, so we have to undef the defines.
 */

#undef SCIPisInfinity
#undef SCIPisEQ
#undef SCIPisLT
#undef SCIPisLE
#undef SCIPisGT
#undef SCIPisGE
#undef SCIPisZero
#undef SCIPisPositive
#undef SCIPisNegative
#undef SCIPisIntegral
#undef SCIPisScalingIntegral
#undef SCIPisFracIntegral
#undef SCIPfloor
#undef SCIPceil
#undef SCIPround
#undef SCIPfrac
#undef SCIPisSumEQ
#undef SCIPisSumLT
#undef SCIPisSumLE
#undef SCIPisSumGT
#undef SCIPisSumGE
#undef SCIPisSumZero
#undef SCIPisSumPositive
#undef SCIPisSumNegative
#undef SCIPisFeasEQ
#undef SCIPisFeasLT
#undef SCIPisFeasLE
#undef SCIPisFeasGT
#undef SCIPisFeasGE
#undef SCIPisFeasZero
#undef SCIPisFeasPositive
#undef SCIPisFeasNegative
#undef SCIPisFeasIntegral
#undef SCIPisFeasFracIntegral
#undef SCIPfeasFloor
#undef SCIPfeasCeil
#undef SCIPfeasRound
#undef SCIPfeasFrac
#undef SCIPisLbBetter
#undef SCIPisUbBetter
#undef SCIPisRelEQ
#undef SCIPisRelLT
#undef SCIPisRelLE
#undef SCIPisRelGT
#undef SCIPisRelGE
#undef SCIPisSumRelEQ
#undef SCIPisSumRelLT
#undef SCIPisSumRelLE
#undef SCIPisSumRelGT
#undef SCIPisSumRelGE
#undef SCIPcreateRealarray
#undef SCIPfreeRealarray
#undef SCIPextendRealarray
#undef SCIPclearRealarray
#undef SCIPgetRealarrayVal
#undef SCIPsetRealarrayVal
#undef SCIPincRealarrayVal
#undef SCIPgetRealarrayMinIdx
#undef SCIPgetRealarrayMaxIdx
#undef SCIPcreateIntarray
#undef SCIPfreeIntarray
#undef SCIPextendIntarray
#undef SCIPclearIntarray
#undef SCIPgetIntarrayVal
#undef SCIPsetIntarrayVal
#undef SCIPincIntarrayVal
#undef SCIPgetIntarrayMinIdx
#undef SCIPgetIntarrayMaxIdx
#undef SCIPcreateBoolarray
#undef SCIPfreeBoolarray
#undef SCIPextendBoolarray
#undef SCIPclearBoolarray
#undef SCIPgetBoolarrayVal
#undef SCIPsetBoolarrayVal
#undef SCIPgetBoolarrayMinIdx
#undef SCIPgetBoolarrayMaxIdx
#undef SCIPcreatePtrarray
#undef SCIPfreePtrarray
#undef SCIPextendPtrarray
#undef SCIPclearPtrarray
#undef SCIPgetPtrarrayVal
#undef SCIPsetPtrarrayVal
#undef SCIPgetPtrarrayMinIdx
#undef SCIPgetPtrarrayMaxIdx

/** checks, if values are in range of epsilon */
SCIP_Bool SCIPisEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/
   
   return SCIPsetIsEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) lower than val2 */
SCIP_Bool SCIPisLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) greater than val2 */
SCIP_Bool SCIPisLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than epsilon) greater than val2 */
SCIP_Bool SCIPisGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than epsilon) lower than val2 */
SCIP_Bool SCIPisGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsGE(scip->set, val1, val2);
}

/** checks, if value is (positive) infinite */
SCIP_Bool SCIPisInfinity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to be compared against infinity */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsInfinity(scip->set, val);
}

/** checks, if value is in range epsilon of 0.0 */
SCIP_Bool SCIPisZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsZero(scip->set, val);
}

/** checks, if value is greater than epsilon */
SCIP_Bool SCIPisPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsPositive(scip->set, val);
}

/** checks, if value is lower than -epsilon */
SCIP_Bool SCIPisNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsNegative(scip->set, val);
}

/** checks, if value is integral within epsilon */
SCIP_Bool SCIPisIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsIntegral(scip->set, val);
}

/** checks whether the product val * scalar is integral in epsilon scaled by scalar */
SCIP_Bool SCIPisScalingIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val,                /**< unscaled value to check for scaled integrality */
   SCIP_Real             scalar              /**< value to scale val with for checking for integrality */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsScalingIntegral(scip->set, val, scalar);
}

/** checks, if given fractional part is smaller than epsilon */
SCIP_Bool SCIPisFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFracIntegral(scip->set, val);
}

/** rounds value + epsilon down to the next integer */
SCIP_Real SCIPfloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFloor(scip->set, val);
}

/** rounds value - epsilon up to the next integer */
SCIP_Real SCIPceil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetCeil(scip->set, val);
}

/** rounds value to the nearest integer with epsilon tolerance */
SCIP_Real SCIPround(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetRound(scip->set, val);
}

/** returns fractional part of value, i.e. x - floor(x) in epsilon tolerance */
SCIP_Real SCIPfrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to return fractional part for */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFrac(scip->set, val);
}

/** checks, if values are in range of sumepsilon */
SCIP_Bool SCIPisSumEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/
   
   return SCIPsetIsSumEQ(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPisSumLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumLT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPisSumLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);
   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumLE(scip->set, val1, val2);
}

/** checks, if val1 is (more than sumepsilon) greater than val2 */
SCIP_Bool SCIPisSumGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumGT(scip->set, val1, val2);
}

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
SCIP_Bool SCIPisSumGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumGE(scip->set, val1, val2);
}

/** checks, if value is in range sumepsilon of 0.0 */
SCIP_Bool SCIPisSumZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumZero(scip->set, val);
}

/** checks, if value is greater than sumepsilon */
SCIP_Bool SCIPisSumPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumPositive(scip->set, val);
}

/** checks, if value is lower than -sumepsilon */
SCIP_Bool SCIPisSumNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsSumNegative(scip->set, val);
}

/** checks, if relative difference of values is in range of feasibility tolerance */
SCIP_Bool SCIPisFeasEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsFeasEQ(scip->set, val1, val2);
}

/** checks, if relative difference val1 and val2 is lower than feasibility tolerance */
SCIP_Bool SCIPisFeasLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsFeasLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than feasibility tolerance */
SCIP_Bool SCIPisFeasLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsFeasLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than feastol */
SCIP_Bool SCIPisFeasGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsFeasGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -feastol */
SCIP_Bool SCIPisFeasGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsFeasGE(scip->set, val1, val2);
}

/** checks, if value is in range feasibility tolerance of 0.0 */
SCIP_Bool SCIPisFeasZero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasZero(scip->set, val);
}

/** checks, if value is greater than feasibility tolerance */
SCIP_Bool SCIPisFeasPositive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasPositive(scip->set, val);
}

/** checks, if value is lower than -feasibility tolerance */
SCIP_Bool SCIPisFeasNegative(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasNegative(scip->set, val);
}

/** checks, if value is integral within the LP feasibility bounds */
SCIP_Bool SCIPisFeasIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasIntegral(scip->set, val);
}

/** checks, if given fractional part is smaller than feastol */
SCIP_Bool SCIPisFeasFracIntegral(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetIsFeasFracIntegral(scip->set, val);
}

/** rounds value + feasibility tolerance down to the next integer */
SCIP_Real SCIPfeasFloor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasFloor(scip->set, val);
}

/** rounds value - feasibility tolerance up to the next integer */
SCIP_Real SCIPfeasCeil(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasCeil(scip->set, val);
}

/** rounds value - feasibility tolerance up to the next integer in feasibility tolerance */
SCIP_Real SCIPfeasRound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasRound(scip->set, val);
}

/** returns fractional part of value, i.e. x - floor(x) */
SCIP_Real SCIPfeasFrac(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val                 /**< value to process */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFeasFrac(scip->set, val);
}

/** checks, if the given new lower bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
SCIP_Bool SCIPisLbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   assert(scip != NULL);

   return SCIPsetIsLbBetter(scip->set, newlb, oldlb, oldub);
}

/** checks, if the given new upper bound is tighter (w.r.t. bound strengthening epsilon) than the old one */
SCIP_Bool SCIPisUbBetter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             newub,              /**< new upper bound */
   SCIP_Real             oldlb,              /**< old lower bound */
   SCIP_Real             oldub               /**< old upper bound */
   )
{
   assert(scip != NULL);

   return SCIPsetIsUbBetter(scip->set, newub, oldlb, oldub);
}

/** checks, if relative difference of values is in range of epsilon */
SCIP_Bool SCIPisRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsRelEQ(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is lower than epsilon */
SCIP_Bool SCIPisRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsRelLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
SCIP_Bool SCIPisRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsRelLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than epsilon */
SCIP_Bool SCIPisRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsRelGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
SCIP_Bool SCIPisRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsRelGE(scip->set, val1, val2);
}

/** checks, if relative difference of values is in range of sumepsilon */
SCIP_Bool SCIPisSumRelEQ(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumRelEQ(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
SCIP_Bool SCIPisSumRelLT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumRelLT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
SCIP_Bool SCIPisSumRelLE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumRelLE(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
SCIP_Bool SCIPisSumRelGT(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumRelGT(scip->set, val1, val2);
}

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
SCIP_Bool SCIPisSumRelGE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             val1,               /**< first value to be compared */
   SCIP_Real             val2                /**< second value to be compared */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   /* avoid to compare two different infinities; the reason for that is
    * that such a comparison can lead to unexpected results */
   assert( ((!SCIPisInfinity(scip, val1) || !SCIPisInfinity(scip, val2))
         && (!SCIPisInfinity(scip, -val1) || !SCIPisInfinity(scip, -val2)))
      || val1 == val2 );    /*lint !e777*/

   return SCIPsetIsSumRelGE(scip->set, val1, val2);
}

/** creates a dynamic array of real values */
SCIP_RETCODE SCIPcreateRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to store the real array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPrealarrayCreate(realarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of real values */
SCIP_RETCODE SCIPfreeRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY**      realarray           /**< pointer to the real array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPrealarrayFree(realarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPextendRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPextendRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPrealarrayExtend(realarray, scip->set, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic real array */
SCIP_RETCODE SCIPclearRealarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPclearRealarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPrealarrayClear(realarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Real SCIPgetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx                 /**< array index to get value for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRealarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPrealarrayGetVal(realarray, idx);
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPsetRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to set value for */
   SCIP_Real             val                 /**< value to set array index to */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetRealarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPrealarraySetVal(realarray, scip->set, idx, val) );

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPincRealarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray,          /**< dynamic real array */
   int                   idx,                /**< array index to increase value for */
   SCIP_Real             incval              /**< value to increase array index */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPincRealarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPrealarrayIncVal(realarray, scip->set, idx, incval) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPgetRealarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRealarrayMinIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPrealarrayGetMinIdx(realarray);
}

/** returns the maximal index of all stored non-zero elements */
int SCIPgetRealarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_REALARRAY*       realarray           /**< dynamic real array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetRealarrayMaxIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPrealarrayGetMaxIdx(realarray);
}

/** creates a dynamic array of int values */
SCIP_RETCODE SCIPcreateIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to store the int array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPintarrayCreate(intarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of int values */
SCIP_RETCODE SCIPfreeIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY**       intarray            /**< pointer to the int array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPintarrayFree(intarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPextendIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPextendIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPintarrayExtend(intarray, scip->set, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic int array */
SCIP_RETCODE SCIPclearIntarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPclearIntarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPintarrayClear(intarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
int SCIPgetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetIntarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPintarrayGetVal(intarray, idx);
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPsetIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   int                   val                 /**< value to set array index to */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetIntarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPintarraySetVal(intarray, scip->set, idx, val) );

   return SCIP_OKAY;
}

/** increases value of entry in dynamic array */
SCIP_RETCODE SCIPincIntarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray,           /**< dynamic int array */
   int                   idx,                /**< array index to increase value for */
   int                   incval              /**< value to increase array index */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPincIntarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPintarrayIncVal(intarray, scip->set, idx, incval) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPgetIntarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetIntarrayMinIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPintarrayGetMinIdx(intarray);
}

/** returns the maximal index of all stored non-zero elements */
int SCIPgetIntarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_INTARRAY*        intarray            /**< dynamic int array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetIntarrayMaxIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPintarrayGetMaxIdx(intarray);
}

/** creates a dynamic array of bool values */
SCIP_RETCODE SCIPcreateBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to store the bool array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreateBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPboolarrayCreate(boolarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of bool values */
SCIP_RETCODE SCIPfreeBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY**      boolarray           /**< pointer to the bool array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreeBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPboolarrayFree(boolarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPextendBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPextendBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPboolarrayExtend(boolarray, scip->set, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic bool array */
SCIP_RETCODE SCIPclearBoolarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPclearBoolarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPboolarrayClear(boolarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
SCIP_Bool SCIPgetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx                 /**< array index to get value for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBoolarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPboolarrayGetVal(boolarray, idx);
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPsetBoolarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray,          /**< dynamic bool array */
   int                   idx,                /**< array index to set value for */
   SCIP_Bool             val                 /**< value to set array index to */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetBoolarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPboolarraySetVal(boolarray, scip->set, idx, val) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPgetBoolarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBoolarrayMinIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPboolarrayGetMinIdx(boolarray);
}

/** returns the maximal index of all stored non-zero elements */
int SCIPgetBoolarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BOOLARRAY*       boolarray           /**< dynamic bool array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetBoolarrayMaxIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPboolarrayGetMaxIdx(boolarray);
}

/** creates a dynamic array of pointers */
SCIP_RETCODE SCIPcreatePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to store the int array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPcreatePtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPptrarrayCreate(ptrarray, SCIPblkmem(scip)) );

   return SCIP_OKAY;
}

/** frees a dynamic array of pointers */
SCIP_RETCODE SCIPfreePtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY**       ptrarray            /**< pointer to the int array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPfreePtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPptrarrayFree(ptrarray) );

   return SCIP_OKAY;
}

/** extends dynamic array to be able to store indices from minidx to maxidx */
SCIP_RETCODE SCIPextendPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   minidx,             /**< smallest index to allocate storage for */
   int                   maxidx              /**< largest index to allocate storage for */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPextendPtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPptrarrayExtend(ptrarray, scip->set, minidx, maxidx) );

   return SCIP_OKAY;
}

/** clears a dynamic pointer array */
SCIP_RETCODE SCIPclearPtrarray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic int array */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPclearPtrarray", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPptrarrayClear(ptrarray) );

   return SCIP_OKAY;
}

/** gets value of entry in dynamic array */
void* SCIPgetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx                 /**< array index to get value for */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPtrarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPptrarrayGetVal(ptrarray, idx);
}

/** sets value of entry in dynamic array */
SCIP_RETCODE SCIPsetPtrarrayVal(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray,           /**< dynamic int array */
   int                   idx,                /**< array index to set value for */
   void*                 val                 /**< value to set array index to */
   )
{
   SCIP_CALL( checkStage(scip, "SCIPsetPtrarrayVal", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   SCIP_CALL( SCIPptrarraySetVal(ptrarray, scip->set, idx, val) );

   return SCIP_OKAY;
}

/** returns the minimal index of all stored non-zero elements */
int SCIPgetPtrarrayMinIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPtrarrayMinIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPptrarrayGetMinIdx(ptrarray);
}

/** returns the maximal index of all stored non-zero elements */
int SCIPgetPtrarrayMaxIdx(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PTRARRAY*        ptrarray            /**< dynamic ptr array */
   )
{
   SCIP_CALL_ABORT( checkStage(scip, "SCIPgetPtrarrayMaxIdx", TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE) );

   return SCIPptrarrayGetMaxIdx(ptrarray);
}
