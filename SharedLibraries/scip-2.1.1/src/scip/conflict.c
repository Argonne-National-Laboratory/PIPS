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

/**@file   conflict.c
 * @brief  methods and datastructures for conflict analysis
 * @author Tobias Achterberg
 *
 * This file implements a conflict analysis method like the one used in modern
 * SAT solvers like zchaff. The algorithm works as follows:
 *
 * Given is a set of bound changes that are not allowed being applied simultaneously, because they
 * render the current node infeasible (e.g. because a single constraint is infeasible in the these
 * bounds, or because the LP relaxation is infeasible).  The goal is to deduce a clause on variables
 * -- a conflict clause -- representing the "reason" for this conflict, i.e., the branching decisions
 * or the deductions (applied e.g. in domain propagation) that lead to the conflict. This clause can
 * then be added to the constraint set to help cutting off similar parts of the branch and bound
 * tree, that would lead to the same conflict.  A conflict clause can also be generated, if the
 * conflict was detected by a locally valid constraint. In this case, the resulting conflict clause
 * is also locally valid in the same depth as the conflict detecting constraint. If all involved
 * variables are binary, a linear (set covering) constraint can be generated, otherwise a bound
 * disjunction constraint is generated. Details are given in
 *
 * Tobias Achterberg, Conflict Analysis in Mixed Integer Programming@n
 * Discrete Optimization, 4, 4-20 (2007)
 *
 * See also @ref CONF. Here is an outline of the algorithm:
 *
 * -#  Put all the given bound changes to a priority queue, which is ordered,
 *     such that the bound change that was applied last due to branching or deduction
 *     is at the top of the queue. The variables in the queue are always active
 *     problem variables. Because binary variables are prefered over general integer 
 *     variables, integer variables are put on the priority queue prior to the binary 
 *     variables. Create an empty conflict set.
 * -#  Remove the top bound change b from the priority queue.
 * -#  Perform the following case distinction:
 *     -#  If the remaining queue is non-empty, and bound change b' (the one that is now
 *         on the top of the queue) was applied at the same depth level as b, and if
 *         b was a deduction with known inference reason, and if the inference constraint's
 *         valid depth is smaller or equal to the conflict detecting constraint's valid
 *         depth:
 *          - Resolve bound change b by asking the constraint that inferred the
 *            bound change to put all the bound changes on the priority queue, that
 *            lead to the deduction of b.
 *            Note that these bound changes have at most the same inference depth
 *            level as b, and were deduced earlier than b.
 *     -#  Otherwise, the bound change b was a branching decision or a deduction with
 *         missing inference reason, or the inference constraint's validity is more local
 *         than the one of the conflict detecing constraint.
 *          - If a the bound changed corresponds to a binary variable, add it or its 
 *            negation to the conflict set, depending on which of them is currently fixed to
 *            FALSE (i.e., the conflict set consists of literals that cannot be FALSE
 *            altogether at the same time).
 *          - Otherwise put the bound change iinto the conflict set.
 *         Note that if the bound change was a branching, all deduced bound changes
 *         remaining in the priority queue have smaller inference depth level than b,
 *         since deductions are always applied after the branching decisions. However,
 *         there is the possibility, that b was a deduction, where the inference
 *         reason was not given or the inference constraint was too local.
 *         With this lack of information, we must treat the deduced bound change like
 *         a branching, and there may exist other deduced bound changes of the same
 *         inference depth level in the priority queue.
 * -#  If priority queue is non-empty, goto step 2.
 * -#  The conflict set represents the conflict clause saying that at least one
 *     of the conflict variables must take a different value. The conflict set is then passed
 *     to the conflict handlers, that may create a corresponding constraint (e.g. a logicor 
 *     constraint or bound disjunction constraint) out of these conflict variables and 
 *     add it to the problem.
 *
 * If all deduced bound changes come with (global) inference information, depending on
 * the conflict analyzing strategy, the resulting conflict set has the following property:
 *  - 1-FirstUIP: In the depth level where the conflict was found, at most one variable
 *    assigned at that level is member of the conflict set. This conflict variable is the
 *    first unique implication point of its depth level (FUIP).
 *  - All-FirstUIP: For each depth level, at most one variable assigned at that level is
 *    member of the conflict set. This conflict variable is the first unique implication
 *    point of its depth level (FUIP).
 *
 * The user has to do the following to get the conflict analysis running in its
 * current implementation:
 *  - A constraint handler or propagator supporting the conflict analysis must implement
 *    the CONSRESPROP/PROPRESPROP call, that processes a bound change inference b and puts all
 *    the reason bounds leading to the application of b with calls to
 *    SCIPaddConflictBound() on the conflict queue (algorithm step 3.(a)).
 *  - If the current bounds lead to a deduction of a bound change (e.g. in domain
 *    propagation), a constraint handler should call SCIPinferVarLbCons() or
 *    SCIPinferVarUbCons(), thus providing the constraint that infered the bound change.
 *    A propagator should call SCIPinferVarLbProp() or SCIPinferVarUbProp() instead,
 *    thus providing a pointer to itself.
 *  - If (in the current bounds) an infeasibility is detected, the constraint handler or
 *    propagator should
 *     1. call SCIPinitConflictAnalysis() to initialize the conflict queue,
 *     2. call SCIPaddConflictBound() for each bound that lead to the conflict,
 *     3. call SCIPanalyzeConflictCons() or SCIPanalyzeConflict() to analyze the conflict
 *        and add an appropriate conflict constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/vbc.h"
#include "scip/lpi.h"
#include "scip/pub_misc.h"
#include "scip/history.h"
#include "scip/paramset.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/prop.h"
#include "scip/debug.h"

#include "scip/struct_conflict.h"



#define CONFLICTSETSCORE(conflictset) (-(conflictset)->nbdchginfos - 100*(conflictset)->repropdepth \
      - 1000*(conflictset)->validdepth)




#ifdef SCIP_CONFGRAPH
/*
 * Output of Conflict Graph
 */

#include <stdio.h>

static FILE*             confgraphfile = NULL;              /**< output file for current conflict graph */
static SCIP_BDCHGINFO*   confgraphcurrentbdchginfo = NULL;  /**< currently resolved bound change */
static int               confgraphnconflictsets = 0;        /**< number of conflict sets marked in the graph */
 
/** writes a node section to the conflict graph file */
static
void confgraphWriteNode(
   void*                 idptr,              /**< id of the node */
   const char*           label,              /**< label of the node */
   const char*           nodetype,           /**< type of the node */
   const char*           fillcolor,          /**< color of the node's interior */
   const char*           bordercolor         /**< color of the node's border */
   )
{
   assert(confgraphfile != NULL);

   fprintf(confgraphfile, "  node\n");
   fprintf(confgraphfile, "  [\n");
   fprintf(confgraphfile, "    id      %d\n", (int)(size_t)idptr);
   fprintf(confgraphfile, "    label   \"%s\"\n", label);
   fprintf(confgraphfile, "    graphics\n");
   fprintf(confgraphfile, "    [\n");
   fprintf(confgraphfile, "      w       120.0\n");
   fprintf(confgraphfile, "      h       30.0\n");
   fprintf(confgraphfile, "      type    \"%s\"\n", nodetype);
   fprintf(confgraphfile, "      fill    \"%s\"\n", fillcolor);
   fprintf(confgraphfile, "      outline \"%s\"\n", bordercolor);
   fprintf(confgraphfile, "    ]\n");
   fprintf(confgraphfile, "    LabelGraphics\n");
   fprintf(confgraphfile, "    [\n");
   fprintf(confgraphfile, "      text      \"%s\"\n", label);
   fprintf(confgraphfile, "      fontSize  13\n");
   fprintf(confgraphfile, "      fontName  \"Dialog\"\n");
   fprintf(confgraphfile, "      anchor    \"c\"\n");
   fprintf(confgraphfile, "    ]\n");
   fprintf(confgraphfile, "  ]\n");
}

/** writes an edge section to the conflict graph file */
static
void confgraphWriteEdge(
   void*                 source,             /**< source node of the edge */
   void*                 target,             /**< target node of the edge */
   const char*           color               /**< color of the edge */
   )
{
   assert(confgraphfile != NULL);


   fprintf(confgraphfile, "  edge\n");
   fprintf(confgraphfile, "  [\n");
   fprintf(confgraphfile, "    source  %d\n", (int)(size_t)source);
   fprintf(confgraphfile, "    target  %d\n", (int)(size_t)target);
   fprintf(confgraphfile, "    graphics\n");
   fprintf(confgraphfile, "    [\n");
   fprintf(confgraphfile, "      fill    \"%s\"\n", color);
   fprintf(confgraphfile, "      targetArrow     \"standard\"\n");
   fprintf(confgraphfile, "    ]\n");
   fprintf(confgraphfile, "  ]\n");
}

/** creates a file to output the current conflict graph into; adds the conflict vertex to the graph */
static
void confgraphCreate(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   char fname[SCIP_MAXSTRLEN];

   assert(conflict != NULL);
   assert(confgraphfile == NULL);

   (void) SCIPsnprintf(fname, SCIP_MAXSTRLEN, "conf%d.gml", conflict->count);
   printf("storing conflict graph in file <%s>\n", fname);
   confgraphfile = fopen(fname, "w");
   if( confgraphfile == NULL )
   {
      SCIPerrorMessage("cannot open conflict graph file <%s>\n", fname);
      SCIPABORT();
   }

   fprintf(confgraphfile, "graph\n");
   fprintf(confgraphfile, "[\n");
   fprintf(confgraphfile, "  hierarchic      1\n");
   fprintf(confgraphfile, "  directed        1\n");
 
   confgraphWriteNode(NULL, "conflict", "ellipse", "#ff0000", "#000000");

   confgraphcurrentbdchginfo = NULL;
}

/** closes conflict graph file */
static
void confgraphFree(
   void
   )
{
   if( confgraphfile != NULL )
   {
      fprintf(confgraphfile, "]\n");
      
      fclose(confgraphfile);
      confgraphfile = NULL;
      confgraphnconflictsets = 0;
   }
}

/** adds a bound change node to the conflict graph and links it to the currently resolved bound change */
static
void confgraphAddBdchg(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict graph */
   )
{
   const char* colors[] = {
      "#8888ff", /**< blue for constraint resolving */
      "#ffff00", /**< yellow for propagator resolving */
      "#55ff55"  /**< green branching decision */
   };
   char label[SCIP_MAXSTRLEN];
   char depth[SCIP_MAXSTRLEN];
   int col;

   
   switch( SCIPbdchginfoGetChgtype(bdchginfo) )
   {
   case SCIP_BOUNDCHGTYPE_BRANCHING:
      col = 2;
      break;
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      col = 0;
      break;
   case SCIP_BOUNDCHGTYPE_PROPINFER:
      col = (SCIPbdchginfoGetInferProp(bdchginfo) == NULL ? 1 : 0);
      break;
   default:
      SCIPerrorMessage("invalid bound change type\n");
      col = 0;
      SCIPABORT();
      break;
   }

   if( SCIPbdchginfoGetDepth(bdchginfo) == INT_MAX )
      (void) SCIPsnprintf(depth, SCIP_MAXSTRLEN, "dive");
   else
      (void) SCIPsnprintf(depth, SCIP_MAXSTRLEN, "%d", SCIPbdchginfoGetDepth(bdchginfo));
   (void) SCIPsnprintf(label, SCIP_MAXSTRLEN, "%s %s %g\n[%s:%d]", SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
      SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
      SCIPbdchginfoGetNewbound(bdchginfo), depth, SCIPbdchginfoGetPos(bdchginfo));
   confgraphWriteNode(bdchginfo, label, "ellipse", colors[col], "#000000");
   confgraphWriteEdge(bdchginfo, confgraphcurrentbdchginfo, "#000000");
}

/** links the already existing bound change node to the currently resolved bound change */
static
void confgraphLinkBdchg(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict graph */
   )
{
   confgraphWriteEdge(bdchginfo, confgraphcurrentbdchginfo, "#000000");
}

/** marks the given bound change to be the currently resolved bound change */
static
void confgraphSetCurrentBdchg(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict graph */
   )
{
   confgraphcurrentbdchginfo = bdchginfo;
}

/** marks given conflict set in the conflict graph */
static
void confgraphMarkConflictset(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   char label[SCIP_MAXSTRLEN];
   int i;

   assert(conflictset != NULL);

   confgraphnconflictsets++;
   (void) SCIPsnprintf(label, SCIP_MAXSTRLEN, "conf %d (%d)", confgraphnconflictsets, conflictset->validdepth);
   confgraphWriteNode((void*)(size_t)confgraphnconflictsets, label, "rectangle", "#ff00ff", "#000000");
   for( i = 0; i < conflictset->nbdchginfos; ++i )
      confgraphWriteEdge((void*)(size_t)confgraphnconflictsets, conflictset->bdchginfos[i], "#ff00ff");
}

#endif



/*
 * Conflict Handler
 */

/** compares two conflict handlers w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrComp)
{  /*lint --e{715}*/
   return ((SCIP_CONFLICTHDLR*)elem2)->priority - ((SCIP_CONFLICTHDLR*)elem1)->priority;
}

/** method to call, when the priority of a conflict handler was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdConflicthdlrPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetConflicthdlrPriority() to mark the conflicthdlrs unsorted */
   SCIP_CALL( SCIPsetConflicthdlrPriority(scip, (SCIP_CONFLICTHDLR*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** copies the given conflict handler to a new scip */
SCIP_RETCODE SCIPconflicthdlrCopyInclude(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);
   assert(set->scip != NULL);

   if( conflicthdlr->conflictcopy != NULL )
   {
      SCIPdebugMessage("including conflict handler %s in subscip %p\n", SCIPconflicthdlrGetName(conflicthdlr), (void*)set->scip);
      SCIP_CALL( conflicthdlr->conflictcopy(set->scip, conflicthdlr) );
   }

   return SCIP_OKAY;
}

/** creates a conflict handler */
SCIP_RETCODE SCIPconflicthdlrCreate(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
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
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(conflicthdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocMemory(conflicthdlr) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conflicthdlr)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conflicthdlr)->desc, desc, strlen(desc)+1) );
   (*conflicthdlr)->priority = priority;
   (*conflicthdlr)->conflictcopy = conflictcopy;
   (*conflicthdlr)->conflictfree = conflictfree;
   (*conflicthdlr)->conflictinit = conflictinit;
   (*conflicthdlr)->conflictexit = conflictexit;
   (*conflicthdlr)->conflictinitsol = conflictinitsol;
   (*conflicthdlr)->conflictexitsol = conflictexitsol;
   (*conflicthdlr)->conflictexec = conflictexec;
   (*conflicthdlr)->conflicthdlrdata = conflicthdlrdata;
   (*conflicthdlr)->initialized = FALSE;

   /* add parameters */
   (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "conflict/%s/priority", name);
   (void) SCIPsnprintf(paramdesc, SCIP_MAXSTRLEN, "priority of conflict handler <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*conflicthdlr)->priority, TRUE, priority, INT_MIN, INT_MAX,
         paramChgdConflicthdlrPriority, (SCIP_PARAMDATA*)(*conflicthdlr)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of conflict handler */
SCIP_RETCODE SCIPconflicthdlrFree(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(*conflicthdlr != NULL);
   assert(!(*conflicthdlr)->initialized);
   assert(set != NULL);

   /* call destructor of conflict handler */
   if( (*conflicthdlr)->conflictfree != NULL )
   {
      SCIP_CALL( (*conflicthdlr)->conflictfree(set->scip, *conflicthdlr) );
   }

   BMSfreeMemoryArray(&(*conflicthdlr)->name);
   BMSfreeMemoryArray(&(*conflicthdlr)->desc);
   BMSfreeMemory(conflicthdlr);

   return SCIP_OKAY;
}

/** calls initialization method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   if( conflicthdlr->initialized )
   {
      SCIPerrorMessage("conflict handler <%s> already initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call initialization method of conflict handler */
   if( conflicthdlr->conflictinit != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictinit(set->scip, conflicthdlr) );
   }
   conflicthdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   if( !conflicthdlr->initialized )
   {
      SCIPerrorMessage("conflict handler <%s> not initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of conflict handler */
   if( conflicthdlr->conflictexit != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictexit(set->scip, conflicthdlr) );
   }
   conflicthdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs conflict handler that the branch and bound process is being started */
SCIP_RETCODE SCIPconflicthdlrInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   /* call solving process initialization method of conflict handler */
   if( conflicthdlr->conflictinitsol != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictinitsol(set->scip, conflicthdlr) );
   }

   return SCIP_OKAY;
}

/** informs conflict handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPconflicthdlrExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of conflict handler */
   if( conflicthdlr->conflictexitsol != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictexitsol(set->scip, conflicthdlr) );
   }

   return SCIP_OKAY;
}

/** calls execution method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExec(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node to add conflict constraint to */
   SCIP_NODE*            validnode,          /**< node at which the constraint is valid */
   SCIP_BDCHGINFO**      bdchginfos,         /**< bound change resembling the conflict set */
   int                   nbdchginfos,        /**< number of bound changes in the conflict set */
   SCIP_Bool             resolved,           /**< was the conflict set already used to create a constraint? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);
   assert(bdchginfos != NULL || nbdchginfos == 0);
   assert(result != NULL);

   /* call solution start method of conflict handler */
   *result = SCIP_DIDNOTRUN;
   if( conflicthdlr->conflictexec != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictexec(set->scip, conflicthdlr, node, validnode, bdchginfos, nbdchginfos,
            set->conf_seperate, (SCIPnodeGetDepth(validnode) > 0), set->conf_dynamic, set->conf_removable, resolved, result) );

      if( *result != SCIP_CONSADDED
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("execution method of conflict handler <%s> returned invalid result <%d>\n",
            conflicthdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** gets user data of conflict handler */
SCIP_CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->conflicthdlrdata;
}

/** sets user data of conflict handler; user has to free old data in advance! */
void SCIPconflicthdlrSetData(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflicthdlrdata = conflicthdlrdata;
}

/** gets name of conflict handler */
const char* SCIPconflicthdlrGetName(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->name;
}

/** gets description of conflict handler */
const char* SCIPconflicthdlrGetDesc(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->desc;
}

/** gets priority of conflict handler */
int SCIPconflicthdlrGetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->priority;
}

/** sets priority of conflict handler */
void SCIPconflicthdlrSetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the conflict handler */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   conflicthdlr->priority = priority;
   set->conflicthdlrssorted = FALSE;
}

/** is conflict handler initialized? */
SCIP_Bool SCIPconflicthdlrIsInitialized(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->initialized;
}




/*
 * Conflict Sets
 */

/** resizes the array of the temporary bound change informations to be able to store at least num bound change entries */
static
SCIP_RETCODE conflictEnsureTmpbdchginfosMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in arrays */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->tmpbdchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->tmpbdchginfos, newsize) );
      conflict->tmpbdchginfossize = newsize;
   }
   assert(num <= conflict->tmpbdchginfossize);

   return SCIP_OKAY;
}

/** creates a temporary bound change information object that is destroyed after the conflict sets are flushed */
static
SCIP_RETCODE conflictCreateTmpBdchginfo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< active variable that changed the bounds */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound for var: lower or upper bound */
   SCIP_Real             oldbound,           /**< old value for bound */
   SCIP_Real             newbound,           /**< new value for bound */
   SCIP_BDCHGINFO**      bdchginfo           /**< pointer to store bound change information */
   )
{
   assert(conflict != NULL);

   SCIP_CALL( conflictEnsureTmpbdchginfosMem(conflict, set, conflict->ntmpbdchginfos+1) );
   SCIP_CALL( SCIPbdchginfoCreate(&conflict->tmpbdchginfos[conflict->ntmpbdchginfos], blkmem, 
         var, boundtype, oldbound, newbound) );
   *bdchginfo = conflict->tmpbdchginfos[conflict->ntmpbdchginfos];
   conflict->ntmpbdchginfos++;

   return SCIP_OKAY;
}

/** frees all temporarily created bound change information data */
static
void conflictFreeTmpBdchginfos(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int i;

   assert(conflict != NULL);

   for( i = 0; i < conflict->ntmpbdchginfos; ++i )
      SCIPbdchginfoFree(&conflict->tmpbdchginfos[i], blkmem);
   conflict->ntmpbdchginfos = 0;
}

/** clears the given conflict set */
static
void conflictsetClear(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   assert(conflictset != NULL);

   conflictset->nbdchginfos = 0;
   conflictset->validdepth = 0;
   conflictset->insertdepth = 0;
   conflictset->conflictdepth = 0;
   conflictset->repropdepth = 0;
   conflictset->repropagate = TRUE;
}

/** creates an empty conflict set */
static
SCIP_RETCODE conflictsetCreate(
   SCIP_CONFLICTSET**    conflictset,        /**< pointer to store the conflict set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflictset != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, conflictset) );
   (*conflictset)->bdchginfos = NULL;
   (*conflictset)->sortvals = NULL;
   (*conflictset)->bdchginfossize = 0;

   conflictsetClear(*conflictset);

   return SCIP_OKAY;
}

/** creates a copy of the given conflict set, allocating an additional amount of memory */
static
SCIP_RETCODE conflictsetCopy(
   SCIP_CONFLICTSET**    targetconflictset,  /**< pointer to store the conflict set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_CONFLICTSET*     sourceconflictset,  /**< source conflict set */
   int                   nadditionalelems    /**< number of additional elements to allocate memory for */
   )
{
   int targetsize;

   assert(targetconflictset != NULL);
   assert(sourceconflictset != NULL);

   targetsize = sourceconflictset->nbdchginfos + nadditionalelems;
   SCIP_ALLOC( BMSallocBlockMemory(blkmem, targetconflictset) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetconflictset)->bdchginfos, targetsize) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*targetconflictset)->sortvals, targetsize) );
   (*targetconflictset)->bdchginfossize = targetsize;

   BMScopyMemoryArray((*targetconflictset)->bdchginfos, sourceconflictset->bdchginfos, sourceconflictset->nbdchginfos);
   BMScopyMemoryArray((*targetconflictset)->sortvals, sourceconflictset->sortvals, sourceconflictset->nbdchginfos);

   (*targetconflictset)->nbdchginfos = sourceconflictset->nbdchginfos;
   (*targetconflictset)->validdepth = sourceconflictset->validdepth;
   (*targetconflictset)->insertdepth = sourceconflictset->insertdepth;
   (*targetconflictset)->conflictdepth = sourceconflictset->conflictdepth;
   (*targetconflictset)->repropdepth = sourceconflictset->repropdepth;

   return SCIP_OKAY;
}

/** frees a conflict set */
static
void conflictsetFree(
   SCIP_CONFLICTSET**    conflictset,        /**< pointer to the conflict set */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflictset != NULL);
   assert(*conflictset != NULL);

   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictset)->bdchginfos, (*conflictset)->bdchginfossize);
   BMSfreeBlockMemoryArrayNull(blkmem, &(*conflictset)->sortvals, (*conflictset)->bdchginfossize);
   BMSfreeBlockMemory(blkmem, conflictset);
}

/** resizes the arrays of the conflict set to be able to store at least num bound change entries */
static
SCIP_RETCODE conflictsetEnsureBdchginfosMem(
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in arrays */
   )
{
   assert(conflictset != NULL);
   assert(set != NULL);

   if( num > conflictset->bdchginfossize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictset->bdchginfos, conflictset->bdchginfossize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(blkmem, &conflictset->sortvals, conflictset->bdchginfossize, newsize) );
      conflictset->bdchginfossize = newsize;
   }
   assert(num <= conflictset->bdchginfossize);

   return SCIP_OKAY;
}

/** updates the score of the conflict set */
static
SCIP_Real conflictsetCalcScore(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   assert(conflictset != NULL);

   return (SCIP_Real)CONFLICTSETSCORE(conflictset); /*lint !e790*/
}

/** adds a bound change to a conflict set */
static
SCIP_RETCODE conflictsetAddBound(
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict set */
   )
{
   SCIP_BDCHGINFO** bdchginfos;
   int* sortvals;
   SCIP_VAR* var;
   SCIP_BOUNDTYPE boundtype;
   int idx;
   int sortval;
   int i;

   assert(conflictset != NULL);
   assert(bdchginfo != NULL);

   /* allocate memory for additional element */
   SCIP_CALL( conflictsetEnsureBdchginfosMem(conflictset, blkmem, set, conflictset->nbdchginfos+1) );

   /* insert the new bound change in the arrays sorted by increasing variable index and by bound type */
   bdchginfos = conflictset->bdchginfos;
   sortvals = conflictset->sortvals;
   var = SCIPbdchginfoGetVar(bdchginfo);
   boundtype = SCIPbdchginfoGetBoundtype(bdchginfo);
   idx = SCIPvarGetIndex(var);
   assert(idx < INT_MAX/2);
   assert((int)boundtype == 0 || (int)boundtype == 1);
   sortval = 2*idx + (int)boundtype; /* first sorting criteria: variable index, second criteria: boundtype */
   for( i = conflictset->nbdchginfos; i > 0 && sortval < sortvals[i-1]; --i )
   {
      bdchginfos[i] = bdchginfos[i-1];
      sortvals[i] = sortvals[i-1];
   }
   bdchginfos[i] = bdchginfo;
   sortvals[i] = sortval;
   conflictset->nbdchginfos++;

   /* merge multiple bound changes */
   if( i > 0 && sortval == sortvals[i-1] )
   {
      /* this is a multiple bound change: only keep the one with the tighter bound */
      if( SCIPbdchginfoIsTighter(bdchginfo, bdchginfos[i-1]) )
         bdchginfos[i-1] = bdchginfo;

      /* remove the redundant bound change by moving the later ones one slot to the front */
      conflictset->nbdchginfos--;
      for( ; i < conflictset->nbdchginfos; ++i )
      {
         bdchginfos[i] = bdchginfos[i+1];
         sortvals[i] = sortvals[i+1];
      }
   }

   return SCIP_OKAY;
}

/** calculates the conflict and the repropagation depths of the conflict set */
static
void conflictsetCalcConflictDepth(
   SCIP_CONFLICTSET*     conflictset         /**< conflict set */
   )
{
   int maxdepth[2];
   int i;

   assert(conflictset != NULL);
   assert(conflictset->validdepth <= conflictset->insertdepth);
   
   /* get the depth of the last and last but one bound change */
   maxdepth[0] = conflictset->validdepth;
   maxdepth[1] = conflictset->validdepth;
   for( i = 0; i < conflictset->nbdchginfos; ++i )
   {
      int depth;

      depth = SCIPbdchginfoGetDepth(conflictset->bdchginfos[i]);
      assert(depth >= 0);
      if( depth > maxdepth[0] )
      {
         maxdepth[1] = maxdepth[0];
         maxdepth[0] = depth;
      }
      else if( depth > maxdepth[1] )
         maxdepth[1] = depth;
   }
   assert(maxdepth[0] >= maxdepth[1]);

   conflictset->conflictdepth = maxdepth[0];
   conflictset->repropdepth = maxdepth[1];
}

/** identifies the depth, at which the conflict set should be added:
 *  - if the branching rule operates on variables only, and if all branching variables up to a certain
 *    depth level are member of the conflict, the conflict constraint can only be violated in the subtree
 *    of the node at that depth, because in all other nodes, at least one of these branching variables
 *    violates its conflicting bound, such that the conflict constraint is feasible
 *  - if there is at least one branching variable in a node, we assume, that this branching was performed
 *    on variables, and that the siblings of this node are disjunct w.r.t. the branching variables' fixings
 *  - we have to add the conflict set at least in the valid depth of the initial conflict set,
 *    so we start searching at the first branching after this depth level, i.e. validdepth+1
 */
static
SCIP_RETCODE conflictsetCalcInsertDepth(
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_Bool* branchingincluded;
   int currentdepth;
   int i;

   assert(conflictset != NULL);
   assert(set != NULL);
   assert(tree != NULL);

   /* the conflict set must not be inserted prior to its valid depth */
   conflictset->insertdepth = conflictset->validdepth;
   assert(conflictset->insertdepth >= 0);

   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);

   /* mark the levels for which a branching variable is included in the conflict set */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &branchingincluded, currentdepth+2) );
   BMSclearMemoryArray(branchingincluded, currentdepth+2);
   for( i = 0; i < conflictset->nbdchginfos; ++i )
   {
      int depth;
         
      depth = SCIPbdchginfoGetDepth(conflictset->bdchginfos[i]);
      depth = MIN(depth, currentdepth+1); /* put diving/probing/strong branching changes in this depth level */
      branchingincluded[depth] = TRUE;
   }

   /* skip additional depth levels where branching on the conflict variables was applied */
   while( conflictset->insertdepth < currentdepth && branchingincluded[conflictset->insertdepth+1] )
      conflictset->insertdepth++;

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &branchingincluded);

   assert(conflictset->validdepth <= conflictset->insertdepth && conflictset->insertdepth <= currentdepth);

   return SCIP_OKAY;
}

/** checks whether the first conflict set is redundant to the second one */
static
SCIP_Bool conflictsetIsRedundant(
   SCIP_CONFLICTSET*     conflictset1,       /**< first conflict conflict set */
   SCIP_CONFLICTSET*     conflictset2        /**< second conflict conflict set */
   )
{
   int i1;
   int i2;

   assert(conflictset1 != NULL);
   assert(conflictset2 != NULL);

   /* if conflictset1 has smaller validdepth, it is definitely not redundant to conflictset2 */
   if( conflictset1->validdepth < conflictset2->validdepth )
      return FALSE;

   /* check, if all bound changes in conflictset2 are also present at least as tight in conflictset1;
    * we can stop immediately, if more bound changes are remaining in conflictset2 than in conflictset1
    */
   for( i1 = 0, i2 = 0; i2 < conflictset2->nbdchginfos && conflictset1->nbdchginfos - i1 >= conflictset2->nbdchginfos - i2;
        ++i1, ++i2 )
   {
      int sortval;

      assert(i2 == 0 || conflictset2->sortvals[i2-1] < conflictset2->sortvals[i2]);

      sortval = conflictset2->sortvals[i2];
      for( ; i1 < conflictset1->nbdchginfos && conflictset1->sortvals[i1] < sortval; ++i1 )
      {
         /* while scanning conflictset1, check consistency */
         assert(i1 == 0 || conflictset1->sortvals[i1-1] < conflictset1->sortvals[i1]);
      }
      if( i1 >= conflictset1->nbdchginfos || conflictset1->sortvals[i1] > sortval
         || SCIPbdchginfoIsTighter(conflictset2->bdchginfos[i2], conflictset1->bdchginfos[i1]) )
         return FALSE;
   }

   return (i2 == conflictset2->nbdchginfos);
}

#ifdef SCIP_DEBUG
/** prints a conflict set to the screen */
static
void conflictsetPrint(
   SCIP_CONFLICTSET*          conflictset              /**< conflict set */
   )
{
   int i;

   assert(conflictset != NULL);
   for( i = 0; i < conflictset->nbdchginfos; ++i )
      SCIPdebugPrintf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(conflictset->bdchginfos[i]),
         SCIPvarGetName(SCIPbdchginfoGetVar(conflictset->bdchginfos[i])),
         SCIPbdchginfoGetBoundtype(conflictset->bdchginfos[i]) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(conflictset->bdchginfos[i]));
   SCIPdebugPrintf("\n");
}
#endif

/** resizes conflictsets array to be able to store at least num entries */
static
SCIP_RETCODE conflictEnsureConflictsetsMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->conflictsetssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->conflictsets, newsize) );
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->conflictsetscores, newsize) );
      conflict->conflictsetssize = newsize;
   }
   assert(num <= conflict->conflictsetssize);

   return SCIP_OKAY;
}

/** inserts conflict set into sorted conflictsets array and deletes the conflict set pointer */
static
SCIP_RETCODE conflictInsertConflictset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_CONFLICTSET**    conflictset         /**< pointer to conflict set to insert */
   )
{
   SCIP_Real score;
   int pos;
   int i;
   int j;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(conflictset != NULL);
   assert(*conflictset != NULL);
   assert((*conflictset)->validdepth <= (*conflictset)->insertdepth);
   assert(set->conf_allowlocal || (*conflictset)->validdepth == 0);

   /* calculate conflict and repropagation depth */
   conflictsetCalcConflictDepth(*conflictset);

   /* if we apply repropagations, the conflict set should be inserted at most at its repropdepth */
   if( set->conf_repropagate )
      (*conflictset)->insertdepth = MIN((*conflictset)->insertdepth, (*conflictset)->repropdepth);
   else
      (*conflictset)->repropdepth = INT_MAX;
   assert((*conflictset)->insertdepth <= (*conflictset)->repropdepth);

   SCIPdebugMessage("inserting conflict set (valid: %d, insert: %d, conf: %d, reprop: %d):",
      (*conflictset)->validdepth, (*conflictset)->insertdepth, (*conflictset)->conflictdepth, (*conflictset)->repropdepth);
   SCIPdebug(conflictsetPrint(*conflictset));

   /* get the score of the conflict set */
   score = conflictsetCalcScore(*conflictset);

   /* check, if conflict set is redundant to a better conflict set */
   for( pos = 0; pos < conflict->nconflictsets && score < conflict->conflictsetscores[pos]; ++pos )
   {
      /* check if conflict set is redundant with respect to conflictsets[pos] */
      if( conflictsetIsRedundant(*conflictset, conflict->conflictsets[pos]) )
      {
         SCIPdebugMessage(" -> conflict set is redundant to: ");
         SCIPdebug(conflictsetPrint(conflict->conflictsets[pos]));
         conflictsetFree(conflictset, blkmem);
         return SCIP_OKAY;
      }

      /**@todo like in sepastore.c: calculate overlap between conflictsets -> large overlap reduces score */

   }

   /* insert conflictset into the sorted conflictsets array*/
   SCIP_CALL( conflictEnsureConflictsetsMem(conflict, set, conflict->nconflictsets + 1) );
   for( i = conflict->nconflictsets; i > pos; --i )
   {
      assert(score >= conflict->conflictsetscores[i-1]);
      conflict->conflictsets[i] = conflict->conflictsets[i-1];
      conflict->conflictsetscores[i] = conflict->conflictsetscores[i-1];
   }
   conflict->conflictsets[pos] = *conflictset;
   conflict->conflictsetscores[pos] = score;
   conflict->nconflictsets++;

   /* remove worse conflictsets that are redundant to the new conflictset */
   for( i = pos+1, j = pos+1; i < conflict->nconflictsets; ++i )
   {
      if( conflictsetIsRedundant(conflict->conflictsets[i], *conflictset) )
      {
         SCIPdebugMessage(" -> conflict set dominates: ");
         SCIPdebug(conflictsetPrint(conflict->conflictsets[i]));
         conflictsetFree(&conflict->conflictsets[i], blkmem);
      }
      else
      {
         assert(j <= i);
         conflict->conflictsets[j] = conflict->conflictsets[i];
         conflict->conflictsetscores[j] = conflict->conflictsetscores[i];
         j++;
      }
   }
   assert(j <= conflict->nconflictsets);
   conflict->nconflictsets = j;

#ifdef SCIP_CONFGRAPH
   confgraphMarkConflictset(*conflictset);
#endif

   *conflictset = NULL; /* ownership of pointer is now in the conflictsets array */

   return SCIP_OKAY;
}

/** calculates the maximal size of conflict sets to be used */
static
int conflictCalcMaxsize(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   int maxsize;

   assert(set != NULL);
   assert(prob != NULL);

   maxsize = (int)(set->conf_maxvarsfac * (prob->nvars - prob->ncontvars));
   maxsize = MAX(maxsize, set->conf_minmaxvars);
   
   return maxsize;
}

/** adds the given conflict set as conflict constraint to the problem */
static
SCIP_RETCODE conflictAddConflictCons(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CONFLICTSET*     conflictset,        /**< conflict set to add to the tree */
   int                   insertdepth,        /**< depth level at which the conflict set should be added */
   SCIP_Bool*            success             /**< pointer to store whether the addition was successful */
   )
{
   int h;

   assert(conflict != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);
   assert(conflictset != NULL);
   assert(conflictset->validdepth <= insertdepth);
   assert(success != NULL);

   /* sort conflict handlers by priority */
   SCIPsetSortConflicthdlrs(set);

   /* call conflict handlers to create a conflict constraint */
   *success = FALSE;
   for( h = 0; h < set->nconflicthdlrs; ++h )
   {
      SCIP_RESULT result;

      SCIP_CALL( SCIPconflicthdlrExec(set->conflicthdlrs[h], set, tree->path[insertdepth],
            tree->path[conflictset->validdepth], conflictset->bdchginfos, conflictset->nbdchginfos, *success, &result) );
      if( result == SCIP_CONSADDED )
      {
         
         *success = TRUE;
         if( insertdepth > 0 )
         {
            conflict->nappliedlocconss++;
            conflict->nappliedlocliterals += conflictset->nbdchginfos;
         }
         else
         {
            int i;
            int conflictlength;
            conflictlength = conflictset->nbdchginfos;

            for( i = 0; i < conflictlength; i++ )
            {
               SCIP_VAR* var;
               SCIP_BRANCHDIR branchdir;
               unsigned int boundtype;
               
               var = conflictset->bdchginfos[i]->var;
               boundtype =  conflictset->bdchginfos[i]->boundtype;
               assert(stat != NULL);               
               branchdir = (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS); /*lint !e641*/
               
               SCIP_CALL( SCIPvarIncNActiveConflicts(var, branchdir, (SCIP_Real)conflictlength) );
               SCIPhistoryIncNActiveConflicts(stat->glbhistory, branchdir, (SCIP_Real)conflictlength);
               SCIPhistoryIncNActiveConflicts(stat->glbhistorycrun, branchdir, (SCIP_Real)conflictlength);
            }
            conflict->nappliedglbconss++;
            conflict->nappliedglbliterals += conflictset->nbdchginfos;
         }
      }
      SCIPdebugMessage(" -> call conflict handler <%s> (prio=%d) to create conflict set with %d bounds returned result %d\n",
         SCIPconflicthdlrGetName(set->conflicthdlrs[h]), SCIPconflicthdlrGetPriority(set->conflicthdlrs[h]),
         conflictset->nbdchginfos, result);
   }

   return SCIP_OKAY;
}

/** adds the collected conflict constraints to the corresponding nodes; the best set->conf_maxconss conflict constraints
 *  are added to the node of their validdepth; additionally (if not yet added, and if repropagation is activated), the
 *  conflict constraint that triggers the earliest repropagation is added to the node of its validdepth
 */
SCIP_RETCODE SCIPconflictFlushConss(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);

   /* is there anything to do? */
   if( conflict->nconflictsets > 0 )
   {
      SCIP_CONFLICTSET* repropconflictset;
      int nconflictsetsused;
      int focusdepth;
#ifndef NDEBUG
      int currentdepth;
#endif
      int cutoffdepth;
      int repropdepth;
      int maxconflictsets;
      int maxsize;
      int i;

      /* calculate the maximal number of conflict sets to accept, and the maximal size of each accepted conflict set */
      maxconflictsets = (set->conf_maxconss == -1 ? INT_MAX : set->conf_maxconss);
      maxsize = conflictCalcMaxsize(set, prob);

      focusdepth = SCIPtreeGetFocusDepth(tree);
#ifndef NDEBUG
      currentdepth = SCIPtreeGetCurrentDepth(tree);
      assert(focusdepth <= currentdepth);
      assert(currentdepth == tree->pathlen-1);
#endif

      SCIPdebugMessage("flushing %d conflict sets at focus depth %d (maxconflictsets: %d, maxsize: %d)\n",
         conflict->nconflictsets, focusdepth, maxconflictsets, maxsize);

      /* mark the focus node to have produced conflict sets in the VBC tool output */
      SCIPvbcFoundConflict(stat->vbc, stat, tree->path[focusdepth]);

      /* insert the conflict sets at the corresponding nodes */
      nconflictsetsused = 0;
      cutoffdepth = INT_MAX;
      repropdepth = INT_MAX;
      repropconflictset = NULL;
      for( i = 0; i < conflict->nconflictsets && nconflictsetsused < maxconflictsets; ++i )
      {
         SCIP_CONFLICTSET* conflictset;

         conflictset = conflict->conflictsets[i];
         assert(conflictset != NULL);
         assert(0 <= conflictset->validdepth);
         assert(conflictset->validdepth <= conflictset->insertdepth);
         assert(conflictset->insertdepth <= focusdepth);
         assert(conflictset->insertdepth <= conflictset->repropdepth);
         assert(conflictset->repropdepth <= currentdepth || conflictset->repropdepth == INT_MAX); /* INT_MAX for dive/probing/strong */
         assert(conflictset->conflictdepth <= currentdepth || conflictset->conflictdepth == INT_MAX); /* INT_MAX for dive/probing/strong */

         /* ignore conflict sets that are only valid at a node that was already cut off */
         if( conflictset->insertdepth >= cutoffdepth )
         {
            SCIPdebugMessage(" -> ignoring conflict set with insertdepth %d >= cutoffdepth %d\n",
               conflictset->validdepth, cutoffdepth);
            continue;
         }

         /* if no conflict bounds exist, the node and its sub tree in the conflict set's valid depth can be
          * cut off completely
          */
         if( conflictset->nbdchginfos == 0 )
         {
            SCIPdebugMessage(" -> empty conflict set in depth %d cuts off sub tree at depth %d\n",
               focusdepth, conflictset->validdepth);

            SCIPnodeCutoff(tree->path[conflictset->validdepth], set, stat, tree);
            cutoffdepth = conflictset->validdepth;
            continue;
         }

         /* if the conflict set is too long, use the conflict set only if it decreases the repropagation depth */
         if( conflictset->nbdchginfos > maxsize )
         {
            SCIPdebugMessage(" -> conflict set is too long: %d > %d literals\n", conflictset->nbdchginfos, maxsize);
            if( set->conf_keepreprop && conflictset->repropagate && conflictset->repropdepth < repropdepth )
            {
               repropdepth = conflictset->repropdepth;
               repropconflictset = conflictset;
            }
         }
         else
         {
            SCIP_Bool success;

            /* call conflict handlers to create a conflict constraint */
            SCIP_CALL( conflictAddConflictCons(conflict, set, stat, tree, conflictset, conflictset->insertdepth, &success) );

            if( success )
            {
               SCIPdebugMessage(" -> conflict set %d/%d added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf:%d, reprop:%d, len:%d):\n",
                  nconflictsetsused+1, maxconflictsets, SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
                  conflictset->insertdepth, conflictset->validdepth, conflictset->conflictdepth, conflictset->repropdepth,
                  conflictset->nbdchginfos);
               SCIPdebug(conflictsetPrint(conflictset));

               if( conflictset->repropagate && conflictset->repropdepth <= repropdepth )
               {
                  repropdepth = conflictset->repropdepth;
                  repropconflictset = NULL;
               }
               nconflictsetsused++;
            }
         }
      }

      /* reactivate propagation on the first node where one of the new conflict sets trigger a deduction */
      if( set->conf_repropagate && repropdepth < cutoffdepth && repropdepth < tree->pathlen )
      {
         assert(0 <= repropdepth && repropdepth < tree->pathlen);
         assert(tree->path[repropdepth]->depth == repropdepth);

         /* if the conflict constraint of smallest repropagation depth was not yet added, insert it now */
         if( repropconflictset != NULL )
         {
            SCIP_Bool success;

            assert(repropconflictset->repropagate);
            assert(repropconflictset->repropdepth == repropdepth);

            SCIP_CALL( conflictAddConflictCons(conflict, set, stat, tree, repropconflictset, repropdepth, &success) );
#ifdef SCIP_DEBUG
            if( success )
            {
               SCIPdebugMessage(" -> additional reprop conflict set added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf:%d, reprop:%d, len:%d):\n",
                  SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
                  repropconflictset->insertdepth, repropconflictset->validdepth, repropconflictset->conflictdepth,
                  repropconflictset->repropdepth, repropconflictset->nbdchginfos);
               SCIPdebug(conflictsetPrint(repropconflictset));
            }
#endif
         }

         /* mark the node in the repropdepth to be propagated again */
         SCIPnodePropagateAgain(tree->path[repropdepth], set, stat, tree);

         SCIPdebugMessage("marked node %p in depth %d to be repropagated due to conflicts found in depth %d\n",
            (void*)tree->path[repropdepth], repropdepth, focusdepth);
      }

      /* free the conflict storage */
      for( i = 0; i < conflict->nconflictsets; ++i )
      {
         conflictsetFree(&conflict->conflictsets[i], blkmem);
      }
      conflict->nconflictsets = 0;
   }

   /* free all temporarily created bound change information data */
   conflictFreeTmpBdchginfos(conflict, blkmem);

   return SCIP_OKAY;
}

/** returns the current number of conflict sets in the conflict set storage */
int SCIPconflictGetNConflicts(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nconflictsets;
}

/** returns the total number of conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbconss + conflict->nappliedlocconss;
}

/** returns the total number of literals in conflict constraints that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals + conflict->nappliedlocliterals;
}

/** returns the total number of conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbconss;
}

/** returns the total number of literals in conflict constraints that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals;
}

/** returns the total number of conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocconss;
}

/** returns the total number of literals in conflict constraints that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocliterals;
}




/*
 * Propagation Conflict Analysis
 */

/** returns whether bound change has a valid reason that can be resolved in conflict analysis */
static
SCIP_Bool bdchginfoIsResolvable(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   return (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
      || (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_PROPINFER
         && SCIPbdchginfoGetInferProp(bdchginfo) != NULL));
}

/** compares two conflict set entries, such that bound changes infered later are
 *  ordered prior to ones that were infered earlier
 */
static
SCIP_DECL_SORTPTRCOMP(conflictBdchginfoComp)
{  /*lint --e{715}*/
   SCIP_BDCHGINFO* bdchginfo1;
   SCIP_BDCHGINFO* bdchginfo2;

   bdchginfo1 = (SCIP_BDCHGINFO*)elem1;
   bdchginfo2 = (SCIP_BDCHGINFO*)elem2;
   assert(bdchginfo1 != NULL);
   assert(bdchginfo2 != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo1));
   assert(!SCIPbdchginfoIsRedundant(bdchginfo2));

   if( !SCIPbdchgidxIsEarlierNonNull(SCIPbdchginfoGetIdx(bdchginfo1), SCIPbdchginfoGetIdx(bdchginfo2)) )
      return -1;
   else
      return +1;
}

/** creates conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictCreate(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflict != NULL);

   SCIP_ALLOC( BMSallocMemory(conflict) );

   SCIP_CALL( SCIPclockCreate(&(*conflict)->propanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->inflpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->boundlpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->sbanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->pseudoanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->bdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp) );
   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->forcedbdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp) );
   SCIP_CALL( conflictsetCreate(&(*conflict)->conflictset, blkmem) );
   (*conflict)->conflictsets = NULL;
   (*conflict)->conflictsetscores = NULL;
   (*conflict)->tmpbdchginfos = NULL;
   (*conflict)->conflictsetssize = 0;
   (*conflict)->nconflictsets = 0;
   (*conflict)->tmpbdchginfossize = 0;
   (*conflict)->ntmpbdchginfos = 0;
   (*conflict)->count = 0;
   (*conflict)->nappliedglbconss = 0;
   (*conflict)->nappliedglbliterals = 0;
   (*conflict)->nappliedlocconss = 0;
   (*conflict)->nappliedlocliterals = 0;
   (*conflict)->npropcalls = 0;
   (*conflict)->npropsuccess = 0;
   (*conflict)->npropconfconss = 0;
   (*conflict)->npropconfliterals = 0;
   (*conflict)->npropreconvconss = 0;
   (*conflict)->npropreconvliterals = 0;
   (*conflict)->ninflpcalls = 0;
   (*conflict)->ninflpsuccess = 0;
   (*conflict)->ninflpconfconss = 0;
   (*conflict)->ninflpconfliterals = 0;
   (*conflict)->ninflpreconvconss = 0;
   (*conflict)->ninflpreconvliterals = 0;
   (*conflict)->ninflpiterations = 0;
   (*conflict)->nboundlpcalls = 0;
   (*conflict)->nboundlpsuccess = 0;
   (*conflict)->nboundlpconfconss = 0;
   (*conflict)->nboundlpconfliterals = 0;
   (*conflict)->nboundlpreconvconss = 0;
   (*conflict)->nboundlpreconvliterals = 0;
   (*conflict)->nboundlpiterations = 0;
   (*conflict)->nsbcalls = 0;
   (*conflict)->nsbsuccess = 0;
   (*conflict)->nsbconfconss = 0;
   (*conflict)->nsbconfliterals = 0;
   (*conflict)->nsbreconvconss = 0;
   (*conflict)->nsbreconvliterals = 0;
   (*conflict)->nsbiterations = 0;
   (*conflict)->npseudocalls = 0;
   (*conflict)->npseudosuccess = 0;
   (*conflict)->npseudoconfconss = 0;
   (*conflict)->npseudoconfliterals = 0;
   (*conflict)->npseudoreconvconss = 0;
   (*conflict)->npseudoreconvliterals = 0;

   return SCIP_OKAY;
}

/** frees conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictFree(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(conflict != NULL);
   assert(*conflict != NULL);
   assert((*conflict)->nconflictsets == 0);
   assert((*conflict)->ntmpbdchginfos == 0);

#ifdef SCIP_CONFGRAPH
   confgraphFree();
#endif

   SCIPclockFree(&(*conflict)->propanalyzetime);
   SCIPclockFree(&(*conflict)->inflpanalyzetime);
   SCIPclockFree(&(*conflict)->boundlpanalyzetime);
   SCIPclockFree(&(*conflict)->sbanalyzetime);
   SCIPclockFree(&(*conflict)->pseudoanalyzetime);
   SCIPpqueueFree(&(*conflict)->bdchgqueue);
   SCIPpqueueFree(&(*conflict)->forcedbdchgqueue);
   conflictsetFree(&(*conflict)->conflictset, blkmem);
   BMSfreeMemoryArrayNull(&(*conflict)->conflictsets);
   BMSfreeMemoryArrayNull(&(*conflict)->conflictsetscores);
   BMSfreeMemoryArrayNull(&(*conflict)->tmpbdchginfos);
   BMSfreeMemory(conflict);

   return SCIP_OKAY;
}

/** clears the conflict queue and the current conflict set */
static
void conflictClear(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   SCIPpqueueClear(conflict->bdchgqueue);
   SCIPpqueueClear(conflict->forcedbdchgqueue);
   conflictsetClear(conflict->conflictset);
}

/** initializes the propagation conflict analysis by clearing the conflict candidate queue */
SCIP_RETCODE SCIPconflictInit(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);

   SCIPdebugMessage("initializing conflict analysis\n");

   /* clear the conflict candidate queue and the conflict set */
   conflictClear(conflict);

   /* increase the conflict counter, such that binary variables of new conflict set and new conflict queue are labeled
    * with this new counter
    */
   conflict->count++;
   if( conflict->count == 0 ) /* make sure, 0 is not a valid conflict counter (may happen due to integer overflow) */
      conflict->count = 1;

   /* increase the conflict score weight for history updates of future conflict reasons */
   if( stat->nnodes > stat->lastconflictnode )
   {
      assert(0.0 < set->conf_scorefac && set->conf_scorefac <= 1.0);
      stat->vsidsweight /= set->conf_scorefac;
      assert(stat->vsidsweight > 0.0);

      /* if the conflict score for the next conflict exceeds 1000.0, rescale all history conflict scores */
      if( stat->vsidsweight >= 1000.0 )
      {
         int v;

         for( v = 0; v < prob->nvars; ++v )
         {
            SCIP_CALL( SCIPvarScaleVSIDS(prob->vars[v], 1.0/stat->vsidsweight) );
         }
         SCIPhistoryScaleVSIDS(stat->glbhistory, 1.0/stat->vsidsweight);
         SCIPhistoryScaleVSIDS(stat->glbhistorycrun, 1.0/stat->vsidsweight);
         stat->vsidsweight = 1.0;
      }
      stat->lastconflictnode = stat->nnodes;
   }

#ifdef SCIP_CONFGRAPH
   confgraphFree();
   confgraphCreate(conflict);
#endif

   return SCIP_OKAY;
}

/** marks bound to be present in the current conflict and returns whether a bound which is at least as tight was already
 *  member of the current conflict (i.e., the given bound change does not need to be added)
 */
static
SCIP_Bool conflictMarkBoundCheckPresence(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict set */
   )
{
   SCIP_VAR* var;
   SCIP_Real newbound;

   assert(conflict != NULL);

   var = SCIPbdchginfoGetVar(bdchginfo);
   newbound = SCIPbdchginfoGetNewbound(bdchginfo);
   assert(var != NULL);

   switch( SCIPbdchginfoGetBoundtype(bdchginfo) )
   {
   case SCIP_BOUNDTYPE_LOWER:
      if( var->conflictlbcount == conflict->count && var->conflictlb >= newbound )
      {
         SCIPdebugMessage("ignoring redundant bound change <%s> >= %g\n", SCIPvarGetName(var), newbound);
         return TRUE;
      }

      var->conflictlbcount = conflict->count;
      var->conflictlb = newbound;

      return FALSE;

   case SCIP_BOUNDTYPE_UPPER:
      if( var->conflictubcount == conflict->count && var->conflictub <= newbound )
      {
         SCIPdebugMessage("ignoring redundant bound change <%s> <= %g\n", SCIPvarGetName(var), newbound);
         return TRUE;
      }

      var->conflictubcount = conflict->count;
      var->conflictub = newbound;

      return FALSE;

   default:
      SCIPerrorMessage("invalid bound type %d\n", SCIPbdchginfoGetBoundtype(bdchginfo));
      SCIPABORT();
      return FALSE; /*lint !e527*/
   }
}

/** puts bound change into the current conflict set */
static
SCIP_RETCODE conflictAddConflictBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change to add to the conflict set */
   )
{
   assert(conflict != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   SCIPdebugMessage("putting bound change <%s> %s %g at depth %d to current conflict set\n",
      SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)), 
      SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=", SCIPbdchginfoGetNewbound(bdchginfo),
      SCIPbdchginfoGetDepth(bdchginfo));

   /* mark the bound to be member of the conflict and check if a bound which is at least as tight is already member of
    * the conflict
    */
   if( !conflictMarkBoundCheckPresence(conflict, bdchginfo) )
   {
      /* add the bound change to the current conflict set */
      SCIP_CALL( conflictsetAddBound(conflict->conflictset, blkmem, set, bdchginfo) );

#ifdef SCIP_CONFGRAPH
      if( bdchginfo != confgraphcurrentbdchginfo )
         confgraphAddBdchg(bdchginfo);
#endif
   }
#ifdef SCIP_CONFGRAPH
   else
      confgraphLinkBdchg(bdchginfo);
#endif

   return SCIP_OKAY;
}

/** returns whether the negation of the given bound change would lead to a globally valid literal */
static
SCIP_Bool isBoundchgUseless(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   SCIP_VAR* var;
   SCIP_BOUNDTYPE boundtype;
   SCIP_Real bound;

   var = SCIPbdchginfoGetVar(bdchginfo);
   boundtype = SCIPbdchginfoGetBoundtype(bdchginfo);
   bound = SCIPbdchginfoGetNewbound(bdchginfo);

   return (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
      && ((boundtype == SCIP_BOUNDTYPE_LOWER && SCIPsetIsFeasGE(set, bound, SCIPvarGetUbGlobal(var)))
         || (boundtype == SCIP_BOUNDTYPE_UPPER && SCIPsetIsFeasLE(set, bound, SCIPvarGetLbGlobal(var)))));
}

/** adds given bound change information to the conflict candidate queue */
static
SCIP_RETCODE conflictQueueBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);
   assert(bdchginfo != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   /* mark the bound to be member of the conflict and check if a bound which is at least as tight is already member of
    * the conflict
    */
   if( !conflictMarkBoundCheckPresence(conflict, bdchginfo) )
   {
      /* insert the bound change into the conflict queue */
      if( (!set->conf_preferbinary || SCIPvarIsBinary(SCIPbdchginfoGetVar(bdchginfo)))
         && !isBoundchgUseless(set, bdchginfo) )
      {
         SCIP_CALL( SCIPpqueueInsert(conflict->bdchgqueue, (void*)bdchginfo) );
      }
      else
      {
         SCIP_CALL( SCIPpqueueInsert(conflict->forcedbdchgqueue, (void*)bdchginfo) );
      }

#ifdef SCIP_CONFGRAPH
      confgraphAddBdchg(bdchginfo);
#endif
   }
#ifdef SCIP_CONFGRAPH
   else
      confgraphLinkBdchg(bdchginfo);
#endif

   return SCIP_OKAY;
}

/** increases the conflict score of the variable in the given direction */
static
SCIP_RETCODE incVSIDS(
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype           /**< type of bound for which the score should be increased */
   )
{
   SCIP_BRANCHDIR branchdir;

   assert(stat != NULL);

   branchdir = (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS);
   SCIP_CALL( SCIPvarIncVSIDS(var, branchdir, stat->vsidsweight) );
   SCIPhistoryIncVSIDS(stat->glbhistory, branchdir, stat->vsidsweight);
   SCIPhistoryIncVSIDS(stat->glbhistorycrun, branchdir, stat->vsidsweight);

   return SCIP_OKAY;
}

/** adds variable's bound to conflict candidate queue */
SCIP_RETCODE SCIPconflictAddBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real scalar;
   SCIP_Real constant;

   assert(conflict != NULL);
   assert(stat != NULL);
   assert(var != NULL);

   /* get active problem variable */
   scalar = 1.0;
   constant = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&var, &scalar, &constant) );

   /* we can ignore fixed variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
      return SCIP_OKAY;

   assert(!SCIPsetIsZero(set, scalar));

   /* if the scalar of the aggregation is negative, we have to switch the bound type */
   if( scalar < 0.0 )
      boundtype = SCIPboundtypeOpposite(boundtype);

   /* if the variable is multi-aggregated, add the bounds of all aggregation variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_VAR** vars;
      SCIP_Real* scalars;
      int nvars;
      int i;

      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      nvars = SCIPvarGetMultaggrNVars(var);
      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPconflictAddBound(conflict, set, stat, vars[i],
               (scalars[i] < 0.0 ? SCIPboundtypeOpposite(boundtype) : boundtype), bdchgidx) );
      }

      return SCIP_OKAY;
   }
   assert(SCIPvarIsActive(var));

   /* get bound change information */
   bdchginfo = SCIPvarGetBdchgInfo(var, boundtype, bdchgidx, FALSE);

   /* if bound of variable was not changed, we can ignore the conflicting bound */
   if( bdchginfo == NULL )
      return SCIP_OKAY;
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   SCIPdebugMessage(" -> adding bound <%s> %s %.15g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] to candidates\n",
      SCIPvarGetName(var),
      boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
      boundtype == SCIP_BOUNDTYPE_LOWER ?
      SCIPvarGetLbAtIndex(var, bdchgidx, FALSE) : SCIPvarGetUbAtIndex(var, bdchgidx, FALSE),
      SCIPvarGetStatus(var), SCIPvarGetType(var),
      SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
      SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
      : (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
         ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
         : (SCIPbdchginfoGetInferProp(bdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo))
            : "none")),
      SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_BRANCHING ? SCIPbdchginfoGetInferInfo(bdchginfo) : -1);

   /* the local bound change may be resolved and has to be put on the candidate queue;
    * we even put bound changes without inference information on the queue in order to automatically
    * eliminate multiple insertions of the same bound change
    */
   assert(SCIPbdchginfoGetVar(bdchginfo) == var);
   assert(SCIPbdchginfoGetBoundtype(bdchginfo) == boundtype);
   assert(SCIPbdchginfoGetDepth(bdchginfo) >= 0);
   assert(SCIPbdchginfoGetPos(bdchginfo) >= 0);
   assert(SCIPbdchgidxIsEarlier(SCIPbdchginfoGetIdx(bdchginfo), bdchgidx));

   /* put bound change information into priority queue */
   SCIP_CALL( conflictQueueBound(conflict, set, bdchginfo) );
   SCIP_CALL( incVSIDS(stat, var, boundtype) );

   return SCIP_OKAY;
}

/** check if the bound change info (which is the potential next candidate which is queued) is valid for the current
 *  conflict analysis; a bound change info can get invalid if after this one was added to the queue, a weaker bound
 *  change was added to the queue (due the bound widening idea) which immediately makes this bound change redundant; due
 *  to the priority we did not removed that bound change info since that cost O(log(n)); hence we have to skip/ignore it
 *  now
 *
 *  The following situations can occur before for example the bound change info (x >= 3) is potentially popped from the
 *  queue.
 *
 *  Postcondition: the reason why (x >= 3) was queued is that at this time point no lower bound of x was involved yet in
 *                 the current conflict or the lower bound which was involved until then was stronger, e.g., (x >= 2).
 *
 *  1) during the time until (x >= 3) gets potentially popped no weaker lower bound was added to the queue, in that case
 *     the conflictlbcount is valid and conflictlb is 3; that is (var->conflictlbcount == conflict->count &&
 *     var->conflictlb == 3)
 *
 *  2) a weaker bound change info gets queued (e.g., x >= 4); this bound change is popped before (x >= 3) since it has
 *     higher priority (which is the time stamp of the bound change info and (x >= 4) has to be done after (x >= 3)
 *     during propagation or branching)
 *
 *    a) if (x >= 4) is popped and added to the conflict set the conflictlbcount is still valid and conflictlb is at
 *      most 4; that is (var->conflictlbcount == conflict->count && var->conflictlb >= 4); it follows that any bound
 *      change info which is stronger than (x >= 4) gets ignored (for example x >= 2)
 *
 *    b) if (x >= 4) is popped and resolved without introducing a new lower bound on x until (x >= 3) is a potentially
 *       candidate the conflictlbcount indicates that bound change is currently not present; that is
 *       (var->conflictlbcount != conflict->count)
 *
 *    c) if (x >= 4) is popped and resolved and a new lower bound on x (e.g., x >= 2) is introduced until (x >= 3) is
 *       pooped, the conflictlbcount indicates that bound change is currently present; that is (var->conflictlbcount ==
 *       conflict->count); however the (x >= 3) only has be explained if conflictlb matches that one; that is
 *       (var->conflictlb == bdchginfo->newbound); otherwise it redundant/invalid.
 */
static
SCIP_Bool bdchginfoIsInvalid(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   SCIP_VAR* var;

   assert(bdchginfo != NULL);

   var = SCIPbdchginfoGetVar(bdchginfo);
   assert(var != NULL);

   /* the bound change info of a binary (domained) variable can never be invalid since the concepts of relaxed bounds
    * and bound widening do not make sense for these type of variables
    */
   if( SCIPvarIsBinary(var) )
      return FALSE;

   /* check if the bdchginfo is invaild since a tight/weaker bound change was already explained */
   if( SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER )
   {
      if( var->conflictlbcount != conflict->count || var->conflictlb != SCIPbdchginfoGetNewbound(bdchginfo) ) /*lint !e777*/
      {
         assert(!SCIPvarIsBinary(var));
         return TRUE;
      }
   }
   else
   {
      assert(SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_UPPER);

      if( var->conflictubcount != conflict->count || var->conflictub != SCIPbdchginfoGetNewbound(bdchginfo) ) /*lint !e777*/
      {
         assert(!SCIPvarIsBinary(var));
         return TRUE;
      }
   }

   return FALSE;
}

/** removes and returns next conflict analysis candidate from the candidate queue */
static
SCIP_BDCHGINFO* conflictRemoveCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_VAR* var;

   assert(conflict != NULL);

   if( SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0 )
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->forcedbdchgqueue));
   else
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->bdchgqueue));

   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   /* if we have a candidate this one should be valid for the current conflict analysis */
   assert(!bdchginfoIsInvalid(conflict, bdchginfo));

   /* mark the bound change to be no longer in the conflict (it will be either added again to the conflict set or
    * replaced by resolving, which might add a weaker change on the same bound to the queue)
    */
   var = SCIPbdchginfoGetVar(bdchginfo);
   if( SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER )
      var->conflictlbcount = 0;
   else
      var->conflictubcount = 0;

#ifdef SCIP_CONFGRAPH
   confgraphSetCurrentBdchg(bdchginfo);
#endif

   return bdchginfo;
}

/** returns next conflict analysis candidate from the candidate queue without removing it */
static
SCIP_BDCHGINFO* conflictFirstCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   if( SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0 )
   {
      /* get next potetioal candidate */
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->forcedbdchgqueue));

      /* check if this candidate is valid */
      if( bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         SCIPdebugMessage("bound change info [%d:<%s> %s %g] is invaild -> pop it from the force queue\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));

         /* pop the invalid bound change info from the queue */
         (void)(SCIPpqueueRemove(conflict->forcedbdchgqueue));

         /* call method recursively to get next conflict analysis candidate */
         bdchginfo = conflictFirstCand(conflict);
      }
   }
   else
   {
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->bdchgqueue));

      /* check if this candidate is valid */
      if( bdchginfo != NULL && bdchginfoIsInvalid(conflict, bdchginfo) )
      {
         SCIPdebugMessage("bound change info [%d:<%s> %s %g] is invaild -> pop it from the queue\n", SCIPbdchginfoGetDepth(bdchginfo),
            SCIPvarGetName(SCIPbdchginfoGetVar(bdchginfo)),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo));

         /* pop the invalid bound change info from the queue */
         (void)(SCIPpqueueRemove(conflict->bdchgqueue));

         /* call method recursively to get next conflict analysis candidate */
         bdchginfo = conflictFirstCand(conflict);
      }
   }
   assert(bdchginfo == NULL || !SCIPbdchginfoIsRedundant(bdchginfo));

   return bdchginfo;
}

/** adds the current conflict set (extended by all remaining bound changes in the queue) to the pool of conflict sets */
static
SCIP_RETCODE conflictAddConflictset(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the conflict set is valid */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   SCIP_Bool             repropagate,        /**< should the constraint trigger a repropagation? */
   SCIP_Bool*            success,            /**< pointer to store whether the conflict set is valid */
   int*                  nliterals           /**< pointer to store the number of literals in the generated conflictset */
   )
{
   SCIP_CONFLICTSET* conflictset;
   SCIP_BDCHGINFO** bdchginfos;
   int nbdchginfos;
   int currentdepth;
   int focusdepth;
   int i;

   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(success != NULL);
   assert(nliterals != NULL);
   assert(SCIPpqueueNElems(conflict->forcedbdchgqueue) == 0);

   *success = FALSE;
   *nliterals = 0;

   /* check, whether local conflicts are allowed */
   validdepth = MAX(validdepth, conflict->conflictset->validdepth);
   if( !set->conf_allowlocal && validdepth > 0 )
      return SCIP_OKAY;

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);
   assert(0 <= conflict->conflictset->validdepth && conflict->conflictset->validdepth <= currentdepth);
   assert(0 <= validdepth && validdepth <= currentdepth);

   /* get the elements of the bound change queue */
   bdchginfos = (SCIP_BDCHGINFO**)SCIPpqueueElems(conflict->bdchgqueue);
   nbdchginfos = SCIPpqueueNElems(conflict->bdchgqueue);

   /* create a copy of the current conflict set, allocating memory for the additional elements of the queue */
   SCIP_CALL( conflictsetCopy(&conflictset, blkmem, conflict->conflictset, nbdchginfos) );
   conflictset->validdepth = validdepth;
   conflictset->repropagate = repropagate;

   /* add the valid queue elements to the conflict set */
   SCIPdebugMessage("adding %d variables from the queue as temporary conflict variables\n", nbdchginfos);
   for( i = 0; i < nbdchginfos; ++i )
   {
      assert(!SCIPbdchginfoIsRedundant(bdchginfos[i]));

      if( !bdchginfoIsInvalid(conflict, bdchginfos[i]) )
      {
         SCIP_CALL( conflictsetAddBound(conflictset, blkmem, set, bdchginfos[i]) );
      }
   }

   /* calculate the depth, at which the conflictset should be inserted */
   SCIP_CALL( conflictsetCalcInsertDepth(conflictset, set, tree) );
   assert(conflictset->validdepth <= conflictset->insertdepth && conflictset->insertdepth <= currentdepth);
   SCIPdebugMessage(" -> conflict with %d literals found at depth %d is active in depth %d and valid in depth %d\n",
      conflictset->nbdchginfos, currentdepth, conflictset->insertdepth, conflictset->validdepth);

   /* if all branching variables are in the conflict set, the conflict set is of no use;
    * don't use conflict sets that are only valid in the probing path but not in the problem tree
    */
   if( (diving || conflictset->insertdepth < currentdepth) && conflictset->insertdepth <= focusdepth )
   {
      /* if the conflict should not be located only in the subtree where it is useful, put it to its valid depth level */
      if( !set->conf_settlelocal )
         conflictset->insertdepth = conflictset->validdepth;

      *nliterals = conflictset->nbdchginfos;
      SCIPdebugMessage(" -> final conflict set has %d literals\n", *nliterals);

      /* check conflict set on debugging solution */
      SCIP_CALL( SCIPdebugCheckConflict(blkmem, set, tree->path[validdepth],
            conflictset->bdchginfos, conflictset->nbdchginfos) ); /*lint !e506 !e774*/

      /* move conflictset to the conflictset storage */
      SCIP_CALL( conflictInsertConflictset(conflict, blkmem, set, &conflictset) );
      *success = TRUE;
   }
   else
   {
      /* free the temporary conflict set */
      conflictsetFree(&conflictset, blkmem);
   }

   return SCIP_OKAY;
}

/** tries to resolve given bound change
 *   - resolutions on local constraints are only applied, if the constraint is valid at the
 *     current minimal valid depth level, because this depth level is the topmost level to add the conflict
 *     constraint to anyways
 */
static
SCIP_RETCODE conflictResolveBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to resolve */
   int                   validdepth,         /**< minimal depth level at which the conflict is valid */
   SCIP_Bool*            resolved            /**< pointer to store whether the bound change was resolved */
   )
{
   SCIP_VAR* actvar;
   SCIP_CONS* infercons;
   SCIP_PROP* inferprop;
   SCIP_RESULT result;

   assert(conflict != NULL);
   assert(resolved != NULL);
   assert(!SCIPbdchginfoIsRedundant(bdchginfo));

   *resolved = FALSE;

   actvar = SCIPbdchginfoGetVar(bdchginfo);
   assert(actvar != NULL);
   assert(SCIPvarIsActive(actvar));

#ifdef SCIP_DEBUG
   {
      int i;
      int nforcedbdchgqueue;
      int nbdchgqueue;

      SCIPdebugMessage("processing next conflicting bound (depth: %d, valid depth: %d, bdchgtype: %s [%s], vartype: %d): [<%s> %s %g]\n",
         SCIPbdchginfoGetDepth(bdchginfo), validdepth,
         SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
         : SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER ? "cons" : "prop",
         SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "-"
         : SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
         ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
         : SCIPbdchginfoGetInferProp(bdchginfo) == NULL ? "-"
         : SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo)),
         SCIPvarGetType(actvar), SCIPvarGetName(actvar),
         SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(bdchginfo));
      SCIPdebugMessage(" - conflict set       :");
      for( i = 0; i < conflict->conflictset->nbdchginfos; ++i )
         SCIPdebugPrintf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(conflict->conflictset->bdchginfos[i]), 
            SCIPvarGetName(SCIPbdchginfoGetVar(conflict->conflictset->bdchginfos[i])),
            SCIPbdchginfoGetBoundtype(conflict->conflictset->bdchginfos[i]) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(conflict->conflictset->bdchginfos[i]));
      SCIPdebugPrintf("\n");
      SCIPdebugMessage(" - forced candidates  :");
      nforcedbdchgqueue = SCIPpqueueNElems(conflict->forcedbdchgqueue);
      for( i = 0; i < nforcedbdchgqueue; ++i )
      {
         SCIP_BDCHGINFO* info = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->forcedbdchgqueue)[i]);
         SCIPdebugPrintf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
            SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(info));
      }
      SCIPdebugPrintf("\n");
      SCIPdebugMessage(" - optional candidates:");
      nbdchgqueue = SCIPpqueueNElems(conflict->bdchgqueue);
      for( i = 0; i < nbdchgqueue; ++i )
      {
         SCIP_BDCHGINFO* info = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->bdchgqueue)[i]);
         SCIPdebugPrintf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
            bdchginfoIsInvalid(conflict, info) ?  SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=" : "<!>",
            SCIPbdchginfoGetNewbound(info));
      }
      SCIPdebugPrintf("\n");
   }
#endif

   /* check, if the bound change can and should be resolved:
    *  - resolutions on local constraints should only be applied, if the constraint is valid at the
    *    current minimal valid depth level (which is initialized with the valid depth level of the initial
    *    conflict set), because this depth level is the topmost level to add the conflict constraint to anyways
    */
   switch( SCIPbdchginfoGetChgtype(bdchginfo) )
   {
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      infercons = SCIPbdchginfoGetInferCons(bdchginfo);
      assert(infercons != NULL);

      if( SCIPconsIsGlobal(infercons) || SCIPconsGetValidDepth(infercons) <= validdepth )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the constraint that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);

         SCIPdebugMessage("resolving bound <%s> %s %g [status:%d, type:%d, depth:%d, pos:%d]: <%s> %s %g [cons:<%s>(%s), info:%d]\n",
            SCIPvarGetName(actvar),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPvarGetStatus(actvar), SCIPvarGetType(actvar),
            SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar),
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPvarGetBdAtIndex(infervar, inferboundtype, bdchgidx, TRUE),
            SCIPconsGetName(infercons),
            SCIPconsIsGlobal(infercons) ? "global" : "local",
            inferinfo);

         SCIP_CALL( SCIPconsResolvePropagation(infercons, set, infervar, inferinfo, inferboundtype, bdchgidx, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_PROPINFER:
      inferprop = SCIPbdchginfoGetInferProp(bdchginfo);
      if( inferprop != NULL )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the propagator that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);

         SCIPdebugMessage("resolving bound <%s> %s %g [status:%d, depth:%d, pos:%d]: <%s> %s %g [prop:<%s>, info:%d]\n",
            SCIPvarGetName(actvar),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPvarGetStatus(actvar), SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar),
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPvarGetBdAtIndex(infervar, inferboundtype, bdchgidx, TRUE),
            SCIPpropGetName(inferprop), inferinfo);

         SCIP_CALL( SCIPpropResolvePropagation(inferprop, set, infervar, inferinfo, inferboundtype, bdchgidx, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_BRANCHING:
      assert(!(*resolved));
      break;

   default:
      SCIPerrorMessage("invalid bound change type <%d>\n", SCIPbdchginfoGetChgtype(bdchginfo));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** if only one conflicting bound change of the last depth level was used, and if this can be resolved,
 *  creates GRASP-like reconvergence conflict constraints in the conflict graph up to the branching variable of this
 *  depth level
 */
static
SCIP_RETCODE conflictCreateReconvergenceConss(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_BDCHGINFO*       firstuip,           /**< first UIP of conflict graph */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   )
{
   SCIP_BDCHGINFO* uip;
   int firstuipdepth;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;

   assert(conflict != NULL);
   assert(firstuip != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);
   assert(!SCIPbdchginfoIsRedundant(firstuip));

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   /* check, whether local constraints are allowed; however, don't generate reconvergence constraints that are only valid
    * in the probing path and not in the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_allowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   firstuipdepth = SCIPbdchginfoGetDepth(firstuip);

   /* for each succeeding UIP pair of the last depth level, create one reconvergence constraint */
   uip = firstuip;
   while( uip != NULL && SCIPbdchginfoGetDepth(uip) == SCIPbdchginfoGetDepth(firstuip) && bdchginfoIsResolvable(uip) )
   {
      SCIP_BDCHGINFO* oppositeuip;
      SCIP_BDCHGINFO* bdchginfo;
      SCIP_BDCHGINFO* nextuip;
      SCIP_VAR* uipvar;
      SCIP_Real oppositeuipbound;
      SCIP_BOUNDTYPE oppositeuipboundtype;
      int nresolutions;

      assert(!SCIPbdchginfoIsRedundant(uip));

      SCIPdebugMessage("creating reconvergence constraint for UIP <%s> %s %g in depth %d pos %d\n",
         SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPbdchginfoGetBoundtype(uip) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(uip), SCIPbdchginfoGetDepth(uip), SCIPbdchginfoGetPos(uip));

      /* initialize conflict data */
      SCIP_CALL( SCIPconflictInit(conflict, set, stat, prob) );

      /* create a temporary bound change information for the negation of the UIP's bound change;
       * this bound change information is freed in the SCIPconflictFlushConss() call;
       * for reconvergence constraints for continuous variables we can only use the "negation" !(x <= u) == (x >= u);
       * during conflict analysis, we treat a continuous bound "x >= u" in the conflict set as "x > u", and in the
       * generated constraint this is negated again to "x <= u" which is correct.
       */
      uipvar = SCIPbdchginfoGetVar(uip);
      oppositeuipboundtype = SCIPboundtypeOpposite(SCIPbdchginfoGetBoundtype(uip));
      oppositeuipbound = SCIPbdchginfoGetNewbound(uip);
      if( SCIPvarIsIntegral(uipvar) )
      {
         assert(SCIPsetIsIntegral(set, oppositeuipbound));
         oppositeuipbound += (oppositeuipboundtype == SCIP_BOUNDTYPE_LOWER ? +1.0 : -1.0);
      }
      SCIP_CALL( conflictCreateTmpBdchginfo(conflict, blkmem, set, uipvar,
            oppositeuipboundtype, oppositeuipboundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_REAL_MIN : SCIP_REAL_MAX,
            oppositeuipbound, &oppositeuip) );

      /* put the negated UIP into the conflict set */
      SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, oppositeuip) );

      /* put positive UIP into priority queue */
      SCIP_CALL( conflictQueueBound(conflict, set, uip) );

      /* resolve the queue until the next UIP is reached */
      bdchginfo = conflictFirstCand(conflict);
      nextuip = NULL;
      nresolutions = 0;
      while( bdchginfo != NULL && validdepth <= maxvaliddepth )
      {
         SCIP_BDCHGINFO* nextbdchginfo;
         SCIP_Bool forceresolve;
         int bdchgdepth;

         /* check if the next bound change must be resolved in every case */
         forceresolve = (SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0);

         /* remove currently processed candidate and get next conflicting bound from the conflict candidate queue */
         assert(bdchginfo == conflictFirstCand(conflict));
         bdchginfo = conflictRemoveCand(conflict);
         nextbdchginfo = conflictFirstCand(conflict);
         bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
         assert(bdchginfo != NULL);
         assert(!SCIPbdchginfoIsRedundant(bdchginfo));
         assert(nextbdchginfo == NULL || SCIPbdchginfoGetDepth(bdchginfo) >= SCIPbdchginfoGetDepth(nextbdchginfo)
            || forceresolve);
         assert(bdchgdepth <= firstuipdepth);

         /* bound changes that are higher in the tree than the valid depth of the conflict can be ignored;
          * multiple insertions of the same bound change can be ignored
          */
         if( bdchgdepth > validdepth && bdchginfo != nextbdchginfo )
         {
            SCIP_VAR* actvar;
            SCIP_Bool resolved;

            actvar = SCIPbdchginfoGetVar(bdchginfo);
            assert(actvar != NULL);
            assert(SCIPvarIsActive(actvar));

            /* check if we have to resolve the bound change in this depth level
             *  - the starting uip has to be resolved
             *  - a bound change should be resolved, if it is in the fuip's depth level and not the
             *    next uip (i.e., if it is not the last bound change in the fuip's depth level)
             *  - a forced bound change must be resolved in any case
             */
            resolved = FALSE;
            if( bdchginfo == uip
               || (bdchgdepth == firstuipdepth
                  && nextbdchginfo != NULL
                  && SCIPbdchginfoGetDepth(nextbdchginfo) == bdchgdepth)
               || forceresolve )
            {
               SCIP_CALL( conflictResolveBound(conflict, set, bdchginfo, validdepth, &resolved) );
            }

            if( resolved )
               nresolutions++;
            else if( forceresolve )
            {
               /* variable cannot enter the conflict clause: we have to make the conflict clause local, s.t.
                * the unresolved bound change is active in the whole sub tree of the conflict clause
                */
               assert(bdchgdepth >= validdepth);
               validdepth = bdchgdepth;

               SCIPdebugMessage("couldn't resolve forced bound change on <%s> -> new valid depth: %d\n",
                  SCIPvarGetName(actvar), validdepth);
            }
            else if( bdchginfo != uip )
            {
               assert(conflict->conflictset != NULL);
               assert(conflict->conflictset->nbdchginfos >= 1); /* starting UIP is already member of the conflict set */

               /* if this is the first variable of the conflict set besides the current starting UIP, it is the next
                * UIP (or the first unresolvable bound change)
                */
               if( bdchgdepth == firstuipdepth && conflict->conflictset->nbdchginfos == 1 )
               {
                  assert(nextuip == NULL);
                  nextuip = bdchginfo;
               }

               /* put bound change into the conflict set */
               SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, bdchginfo) );
               assert(conflict->conflictset->nbdchginfos >= 2);
            }
            else
               assert(conflictFirstCand(conflict) == NULL); /* the starting UIP was not resolved */
         }

         /* get next conflicting bound from the conflict candidate queue (this does not need to be nextbdchginfo, because
          * due to resolving the bound changes, a variable could be added to the queue which must be
          * resolved before nextbdchginfo)
          */
         bdchginfo = conflictFirstCand(conflict);
      }
      assert(nextuip != uip);

      /* if only one propagation was resolved, the reconvergence constraint is already member of the constraint set
       * (it is exactly the constraint that produced the propagation)
       */
      if( nextuip != NULL && nresolutions >= 2 && bdchginfo == NULL && validdepth <= maxvaliddepth )
      {
         int nlits;
         SCIP_Bool success;

         assert(SCIPbdchginfoGetDepth(nextuip) == SCIPbdchginfoGetDepth(uip));

         SCIPdebugMessage("creating reconvergence constraint from UIP <%s> to UIP <%s> in depth %d with %d literals after %d resolutions\n",
            SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPvarGetName(SCIPbdchginfoGetVar(nextuip)),
            SCIPbdchginfoGetDepth(uip), conflict->conflictset->nbdchginfos, nresolutions);

         /* call the conflict handlers to create a conflict set */
         SCIP_CALL( conflictAddConflictset(conflict, blkmem, set, stat, tree, validdepth, diving, FALSE,
               &success, &nlits) );
         if( success )
         {
            (*nreconvconss)++;
            (*nreconvliterals) += nlits;
         }
      }

      /* clear the conflict candidate queue and the conflict set (to make sure, oppositeuip is not referenced anymore) */
      conflictClear(conflict);

      uip = nextuip;
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  afterwards the conflict queue and the conflict set is cleared
 */
static
SCIP_RETCODE conflictAnalyze(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool             mustresolve,        /**< should the conflict set only be used, if a resolution was applied? */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGINFO** firstuips;
   int nfirstuips;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;
   int resolvedepth;
   int nresolutions;
   int lastconsnresolutions;
   int lastconsresoldepth;

   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(conflict->conflictset->nbdchginfos >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(0 <= validdepth && validdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(nconss != NULL);
   assert(nliterals != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   resolvedepth = ((set->conf_fuiplevels >= 0 && set->conf_fuiplevels <= currentdepth)
      ? currentdepth - set->conf_fuiplevels + 1 : 0);
   assert(0 <= resolvedepth && resolvedepth <= currentdepth + 1);

   /* if we must resolve at least one bound change, find the first UIP at least in the last depth level */
   if( mustresolve )
      resolvedepth = MIN(resolvedepth, currentdepth);

   SCIPdebugMessage("analyzing conflict with %d+%d conflict candidates and starting conflict set of size %d in depth %d (resolvedepth=%d)\n",
      SCIPpqueueNElems(conflict->forcedbdchgqueue), SCIPpqueueNElems(conflict->bdchgqueue),
      conflict->conflictset->nbdchginfos, currentdepth, resolvedepth);

   *nconss = 0;
   *nliterals = 0;
   *nreconvconss = 0;
   *nreconvliterals = 0;

   /* check, whether local conflicts are allowed; however, don't generate conflict constraints that are only valid in the
    * probing path and not in the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_allowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   /* allocate temporary memory for storing first UIPs (in each depth level, at most two bound changes can be flagged
    * as UIP, namely a binary and a non-binary bound change)
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &firstuips, 2*(currentdepth+1)) );

   /* process all bound changes in the conflict candidate queue */
   nresolutions = 0;
   lastconsnresolutions = (mustresolve ? 0 : -1);
   lastconsresoldepth = (mustresolve ? currentdepth : INT_MAX);
   bdchginfo = conflictFirstCand(conflict);
   nfirstuips = 0;
   while( bdchginfo != NULL && validdepth <= maxvaliddepth )
   {
      SCIP_BDCHGINFO* nextbdchginfo;
      SCIP_Bool forceresolve;
      int bdchgdepth;

      assert(!SCIPbdchginfoIsRedundant(bdchginfo));

      /* check if the next bound change must be resolved in every case */
      forceresolve = (SCIPpqueueNElems(conflict->forcedbdchgqueue) > 0);

      /* resolve next bound change in queue */
      bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
      assert(0 <= bdchgdepth && bdchgdepth <= currentdepth);
      assert(SCIPvarIsActive(SCIPbdchginfoGetVar(bdchginfo)));
      assert(bdchgdepth < tree->pathlen);
      assert(tree->path[bdchgdepth] != NULL);
      assert(tree->path[bdchgdepth]->domchg != NULL);
      assert(SCIPbdchginfoGetPos(bdchginfo) < (int)tree->path[bdchgdepth]->domchg->domchgbound.nboundchgs);
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].var
         == SCIPbdchginfoGetVar(bdchginfo));
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].newbound
         == SCIPbdchginfoGetNewbound(bdchginfo)
         || (SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER
            ? SCIPvarGetLbGlobal(SCIPbdchginfoGetVar(bdchginfo)) : SCIPvarGetUbGlobal(SCIPbdchginfoGetVar(bdchginfo)))
         == SCIPbdchginfoGetNewbound(bdchginfo)); /*lint !e777*/
      assert((SCIP_BOUNDTYPE)tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].boundtype
         == SCIPbdchginfoGetBoundtype(bdchginfo));

      /* create intermediate conflict constraint */
      assert(nresolutions >= lastconsnresolutions);
      if( !forceresolve )
      {
         if( nresolutions == lastconsnresolutions )
            lastconsresoldepth = bdchgdepth; /* all intermediate depth levels consisted of only unresolved bound changes */
         else if( bdchgdepth < lastconsresoldepth && (set->conf_interconss == -1 || *nconss < set->conf_interconss) )
         {
            int nlits;
            SCIP_Bool success;

            /* call the conflict handlers to create a conflict set */
            SCIPdebugMessage("creating intermediate conflictset after %d resolutions up to depth %d (valid at depth %d): %d conflict bounds, %d bounds in queue\n",
               nresolutions, bdchgdepth, validdepth, conflict->conflictset->nbdchginfos,
               SCIPpqueueNElems(conflict->bdchgqueue));
            SCIP_CALL( conflictAddConflictset(conflict, blkmem, set, stat, tree, validdepth, diving, TRUE,
                  &success, &nlits) );
            lastconsnresolutions = nresolutions;
            lastconsresoldepth = bdchgdepth;
            if( success )
            {
               (*nconss)++;
               (*nliterals) += nlits;
            }
         }
      }

      /* remove currently processed candidate and get next conflicting bound from the conflict candidate queue */
      assert(bdchginfo == conflictFirstCand(conflict));
      bdchginfo = conflictRemoveCand(conflict);
      nextbdchginfo = conflictFirstCand(conflict);
      assert(bdchginfo != NULL);
      assert(!SCIPbdchginfoIsRedundant(bdchginfo));
      assert(nextbdchginfo == NULL || SCIPbdchginfoGetDepth(bdchginfo) >= SCIPbdchginfoGetDepth(nextbdchginfo)
         || forceresolve);

      /* we don't need to resolve bound changes that are already active in the valid depth of the current conflict set,
       * because the conflict set can only be added locally at the valid depth, and all bound changes applied in this
       * depth or earlier can be removed from the conflict constraint, since they are already applied in the constraint's
       * subtree;
       * if the next bound change on the remaining queue is equal to the current bound change,
       * this is a multiple insertion in the conflict candidate queue and we can ignore the current
       * bound change
       */
      if( bdchgdepth > validdepth && bdchginfo != nextbdchginfo )
      {
         SCIP_VAR* actvar;
         SCIP_Bool resolved;

         actvar = SCIPbdchginfoGetVar(bdchginfo);
         assert(actvar != NULL);
         assert(SCIPvarIsActive(actvar));

         /* check if we want to resolve the bound change in this depth level
          *  - bound changes should be resolved, if
          *     (i)   we must apply at least one resolution and didn't resolve a bound change yet, or
          *     (ii)  their depth level is at least equal to the minimal resolving depth, and
          *           they are not the last remaining conflicting bound change in their depth level
          *     (iii) the bound change resolving is forced (i.e., the forced queue was non-empty)
          */
         resolved = FALSE;
         if( (mustresolve && nresolutions == 0)
            || (bdchgdepth >= resolvedepth
               && nextbdchginfo != NULL
               && SCIPbdchginfoGetDepth(nextbdchginfo) == bdchgdepth)
            || forceresolve )
         {
            SCIP_CALL( conflictResolveBound(conflict, set, bdchginfo, validdepth, &resolved) );
         }

         if( resolved )
            nresolutions++;
         else if( forceresolve )
         {
            /* variable cannot enter the conflict clause: we have to make the conflict clause local, s.t.
             * the unresolved bound change is active in the whole sub tree of the conflict clause
             */
            assert(bdchgdepth >= validdepth);
            validdepth = bdchgdepth;

            SCIPdebugMessage("couldn't resolve forced bound change on <%s> -> new valid depth: %d\n",
               SCIPvarGetName(actvar), validdepth);
         }
         else
         {
            /* if this is a UIP (the last bound change in its depth level), it can be used to generate a
             * UIP reconvergence constraint
             */
            if( nextbdchginfo == NULL || SCIPbdchginfoGetDepth(nextbdchginfo) != bdchgdepth )
            {
               assert(nfirstuips < 2*(currentdepth+1));
               firstuips[nfirstuips] = bdchginfo;
               nfirstuips++;
            }

            /* put variable into the conflict set, using the literal that is currently fixed to FALSE */
            SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, bdchginfo) );
         }
      }

      /* get next conflicting bound from the conflict candidate queue (this needs not to be nextbdchginfo, because
       * due to resolving the bound changes, a bound change could be added to the queue which must be
       * resolved before nextbdchginfo)
       */
      bdchginfo = conflictFirstCand(conflict);
   }

   /* check, if a valid conflict set was found */
   if( bdchginfo == NULL
      && nresolutions > lastconsnresolutions
      && validdepth <= maxvaliddepth
      && (!mustresolve || nresolutions > 0 || conflict->conflictset->nbdchginfos == 0)
      && SCIPpqueueNElems(conflict->forcedbdchgqueue) == 0 )
   {
      int nlits;
      SCIP_Bool success;

      /* call the conflict handlers to create a conflict set */
      SCIP_CALL( conflictAddConflictset(conflict, blkmem, set, stat, tree, validdepth, diving, TRUE, &success, &nlits) );
      if( success )
      {
         (*nconss)++;
         (*nliterals) += nlits;
      }
   }

   /* produce reconvergence constraints defined by succeeding UIP's of the last depth level */
   if( set->conf_reconvlevels != 0 && validdepth <= maxvaliddepth )
   {
      int reconvlevels;
      int i;

      reconvlevels = (set->conf_reconvlevels == -1 ? INT_MAX : set->conf_reconvlevels);
      for( i = 0; i < nfirstuips; ++i )
      {
         if( SCIPbdchginfoHasInferenceReason(firstuips[i])
            && currentdepth - SCIPbdchginfoGetDepth(firstuips[i]) < reconvlevels )
         {
            SCIP_CALL( conflictCreateReconvergenceConss(conflict, blkmem, set, stat, prob, tree, diving,
                  validdepth, firstuips[i], nreconvconss, nreconvliterals) );
         }
      }
   }

   /* free the temporary memory */
   SCIPsetFreeBufferArray(set, &firstuips);

   /* clear the conflict candidate queue and the conflict set */
   conflictClear(conflict);

   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  updates statistics for propagation conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyze(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(conflict->conflictset != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if propagation conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_useprop )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* check, if the conflict set will get too large with high probability */
   if( conflict->conflictset->nbdchginfos + SCIPpqueueNElems(conflict->bdchgqueue)
      + SCIPpqueueNElems(conflict->forcedbdchgqueue) >= 2*conflictCalcMaxsize(set, prob) )
      return SCIP_OKAY;

   SCIPdebugMessage("analyzing conflict after infeasible propagation in depth %d\n", SCIPtreeGetCurrentDepth(tree));

   /* start timing */
   SCIPclockStart(conflict->propanalyzetime, set);

   conflict->npropcalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   SCIP_CALL( conflictAnalyze(conflict, blkmem, set, stat, prob, tree, FALSE, validdepth, TRUE,
         &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
   conflict->npropsuccess += (nconss > 0 ? 1 : 0);
   conflict->npropconfconss += nconss;
   conflict->npropconfliterals += nliterals;
   conflict->npropreconvconss += nreconvconss;
   conflict->npropreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0);

   /* stop timing */
   SCIPclockStop(conflict->propanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing propagation conflicts */
SCIP_Real SCIPconflictGetPropTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->propanalyzetime);
}

/** gets number of calls to propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropcalls;
}

/** gets number of calls to propagation conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPropSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropsuccess;
}

/** gets number of conflict constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfconss;
}

/** gets total number of literals in conflict constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfliterals;
}

/** gets number of reconvergence constraints detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvconss;
}

/** gets total number of literals in reconvergence constraints created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvliterals;
}




/*
 * Infeasible LP Conflict Analysis
 */

/** ensures, that side change arrays can store at least num entries */
static
SCIP_RETCODE ensureSidechgsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int                   num                 /**< minimal number of entries to be able to store in side change arrays */
   )
{
   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);

   if( num > *sidechgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechginds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldrhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewrhss, newsize) );
      *sidechgssize = newsize;
   }
   assert(num <= *sidechgssize);

   return SCIP_OKAY;
}

/** adds removal of row's side to side change arrays; finite sides are only replaced by near infinite sides, such
 *  that the row's sense in the LP solver is not changed
 */
static
SCIP_RETCODE addSideRemoval(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row,                /**< LP row to change the sides for */
   SCIP_Real             lpiinfinity,        /**< value treated as infinity in LP solver */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int*                  nsidechgs           /**< pointer to number of used slots in side change arrays */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real constant;

   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);
   assert(nsidechgs != NULL);

   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);
   constant = SCIProwGetConstant(row);
   assert(!SCIPsetIsInfinity(set, -lhs) || !SCIPsetIsInfinity(set, rhs));

   /* get memory to store additional side change */
   SCIP_CALL( ensureSidechgsSize(set, sidechginds, sidechgoldlhss, sidechgoldrhss, sidechgnewlhss, sidechgnewrhss,
         sidechgssize, (*nsidechgs)+1) );
   assert(*nsidechgs < *sidechgssize);
   assert(*sidechginds != NULL);
   assert(*sidechgoldlhss != NULL);
   assert(*sidechgoldrhss != NULL);
   assert(*sidechgnewlhss != NULL);
   assert(*sidechgnewrhss != NULL);

   /* store side change */
   (*sidechginds)[*nsidechgs] = SCIProwGetLPPos(row);
   if( SCIPsetIsInfinity(set, -lhs) )
   {
      (*sidechgoldlhss)[*nsidechgs] = -lpiinfinity;
      (*sidechgnewlhss)[*nsidechgs] = -lpiinfinity;
   }
   else
   {
      (*sidechgoldlhss)[*nsidechgs] = lhs - constant;
      (*sidechgnewlhss)[*nsidechgs] = -lpiinfinity/2;
   }
   if( SCIPsetIsInfinity(set, rhs) )
   {
      (*sidechgoldrhss)[*nsidechgs] = lpiinfinity;
      (*sidechgnewrhss)[*nsidechgs] = lpiinfinity;
   }
   else
   {
      (*sidechgoldrhss)[*nsidechgs] = rhs - constant;
      (*sidechgnewrhss)[*nsidechgs] = lpiinfinity/2;
   }
   (*nsidechgs)++;

   return SCIP_OKAY;
}

/** ensures, that bound change arrays can store at least num entries */
static
SCIP_RETCODE ensureBdchgsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int                   num                 /**< minimal number of entries to be able to store in side change arrays */
   )
{
   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);

   if( num > *sidechgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechginds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldrhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewrhss, newsize) );
      *sidechgssize = newsize;
   }
   assert(num <= *sidechgssize);

   return SCIP_OKAY;
}

/** inserts variable's new bounds into bound change arrays */
static
SCIP_RETCODE addBdchg(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to change the LP bounds for */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             newub,              /**< new upper bound */
   int**                 bdchginds,          /**< pointer to bound change index array */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays */
   int*                  nbdchgs             /**< pointer to number of used slots in bound change arrays */
   )
{
   assert(newlb <= newub);
   assert(bdchginds != NULL);
   assert(bdchgoldlbs != NULL);
   assert(bdchgoldubs != NULL);
   assert(bdchgnewlbs != NULL);
   assert(bdchgnewubs != NULL);
   assert(bdchgssize != NULL);
   assert(nbdchgs != NULL);
   assert(*nbdchgs <= *bdchgssize);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
   {
      SCIP_COL* col;
      int c;

      col = SCIPvarGetCol(var);
      c = SCIPcolGetLPPos(col);
      if( c >= 0 )
      {
         SCIP_CALL( ensureBdchgsSize(set, bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs,
               bdchgssize, (*nbdchgs)+1) );
         assert(*bdchginds != NULL);
         assert(*bdchgoldlbs != NULL);
         assert(*bdchgoldubs != NULL);
         assert(*bdchgnewlbs != NULL);
         assert(*bdchgnewubs != NULL);

         (*bdchginds)[*nbdchgs] = c;
         (*bdchgoldlbs)[*nbdchgs] = SCIPvarGetLbLP(var, set);
         (*bdchgoldubs)[*nbdchgs] = SCIPvarGetUbLP(var, set);
         (*bdchgnewlbs)[*nbdchgs] = newlb;
         (*bdchgnewubs)[*nbdchgs] = newub;
         (*nbdchgs)++;
      }
   }

   return SCIP_OKAY;
}

/** ensures, that candidate array can store at least num entries */
static
SCIP_RETCODE ensureCandsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR***           cands,              /**< pointer to candidate array */
   SCIP_Real**           candscores,         /**< pointer to candidate score array */
   SCIP_Real**           newbounds,          /**< pointer to candidate new bounds array */
   SCIP_Real**           proofactdeltas,     /**< pointer to candidate proof delta array */
   int*                  candssize,          /**< pointer to size of array */
   int                   num                 /**< minimal number of candidates to store in array */
   )
{
   assert(cands != NULL);
   assert(candssize != NULL);

   if( num > *candssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_CALL( SCIPsetReallocBufferArray(set, cands, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, candscores, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, newbounds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, proofactdeltas, newsize) );
      *candssize = newsize;
   }
   assert(num <= *candssize);

   return SCIP_OKAY;
}

/** adds variable to candidate list, if the current best bound corresponding to the proof coefficient is local;
 *  returns the array position in the candidate list, where the new candidate was inserted, or -1 if the
 *  variable can relaxed to global bounds immediately without increasing the proof's activity;
 *  the candidates are sorted with respect to the following two criteria:
 *  - prefer bound changes that have been applied deeper in the tree, to get a more global conflict
 *  - prefer variables with small Farkas coefficient to get rid of as many bound changes as possible
 */
static
SCIP_RETCODE addCand(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_VAR*             var,                /**< variable to add to candidate array */
   int                   lbchginfopos,       /**< positions of currently active lower bound change information in variable's array */
   int                   ubchginfopos,       /**< positions of currently active upper bound change information in variable's array */
   SCIP_Real             proofcoef,          /**< coefficient of variable in infeasibility/bound proof */
   SCIP_Real             prooflhs,           /**< left hand side of infeasibility/bound proof */
   SCIP_Real             proofact,           /**< activity of infeasibility/bound proof row */
   SCIP_VAR***           cands,              /**< pointer to candidate array for undoing bound changes */
   SCIP_Real**           candscores,         /**< pointer to candidate score array for undoing bound changes */
   SCIP_Real**           newbounds,          /**< pointer to candidate new bounds array for undoing bound changes */
   SCIP_Real**           proofactdeltas,     /**< pointer to proof activity increase array for undoing bound changes */
   int*                  candssize,          /**< pointer to size of cands arrays */
   int*                  ncands,             /**< pointer to count number of candidates in bound change list */
   int                   firstcand           /**< position of first unprocessed bound change candidate */
   )
{
   SCIP_Real oldbound;
   SCIP_Real newbound;
   SCIP_Real proofactdelta;
   SCIP_Real score;
   int depth;
   int i;
   SCIP_Bool resolvable;

   assert(set != NULL);
   assert(var != NULL);
   assert(-1 <= lbchginfopos && lbchginfopos <= var->nlbchginfos);
   assert(-1 <= ubchginfopos && ubchginfopos <= var->nubchginfos);
   assert(!SCIPsetIsZero(set, proofcoef));
   assert(SCIPsetIsGT(set, prooflhs, proofact));
   assert(cands != NULL);
   assert(candscores != NULL);
   assert(newbounds != NULL);
   assert(proofactdeltas != NULL);
   assert(candssize != NULL);
   assert(ncands != NULL);
   assert(*ncands <= *candssize);
   assert(0 <= firstcand && firstcand <= *ncands);

   /* in the infeasibility or dual bound proof, the variable's bound is chosen to maximize the proof's activity */
   if( proofcoef > 0.0 )
   {
      assert(ubchginfopos >= 0); /* otherwise, undoBdchgsProof() should already have relaxed the local bound */

      /* calculate the difference of current bound to the previous bound the variable was set to */
      if( ubchginfopos == var->nubchginfos )
      {
         /* current bound is the strong branching or diving bound */
         oldbound = SCIPvarGetUbLP(var, set);
         newbound = SCIPvarGetUbLocal(var);
         depth = currentdepth+1;
         resolvable = FALSE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         resolvable = bdchginfoIsResolvable(&var->ubchginfos[ubchginfopos]);
         depth = var->ubchginfos[ubchginfopos].bdchgidx.depth;
         oldbound = var->ubchginfos[ubchginfopos].newbound;
         newbound = var->ubchginfos[ubchginfopos].oldbound;
      }
   }
   else
   {
      assert(lbchginfopos >= 0); /* otherwise, undoBdchgsProof() should already have relaxed the local bound */

      /* calculate the difference of current bound to the previous bound the variable was set to */
      if( lbchginfopos == var->nlbchginfos )
      {
         /* current bound is the strong branching or diving bound */
         oldbound = SCIPvarGetLbLP(var, set);
         newbound = SCIPvarGetLbLocal(var);
         depth = currentdepth+1;
         resolvable = FALSE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         resolvable = bdchginfoIsResolvable(&var->lbchginfos[lbchginfopos]);
         depth = var->lbchginfos[lbchginfopos].bdchgidx.depth;
         oldbound = var->lbchginfos[lbchginfopos].newbound;
         newbound = var->lbchginfos[lbchginfopos].oldbound;
      }
   }

   /* calculate the increase in the proof's activity */
   proofactdelta = (newbound - oldbound)*proofcoef;
   assert(proofactdelta > 0.0);

   /* calculate score for undoing the bound change */
   score = 1.0 - proofactdelta/(prooflhs - proofact);
   score = MAX(score, 0.0);
   score += set->conf_depthscorefac * (SCIP_Real)(depth+1)/(SCIP_Real)(currentdepth+1);
   if( !resolvable )
   {
      score += 10.0;
      if( !SCIPvarIsBinary(var) )
         score += 10.0;
   }

   /* get enough memory to store new candidate */
   SCIP_CALL( ensureCandsSize(set, cands, candscores, newbounds, proofactdeltas, candssize, (*ncands)+1) );
   assert(*cands != NULL);
   assert(*candscores != NULL);
   assert(*newbounds != NULL);
   assert(*proofactdeltas != NULL);

   SCIPdebugMessage(" -> local <%s> %s %g, relax <%s> %s %g, proofcoef=%g, dpt=%d, resolve=%u, delta=%g, score=%g\n",
      SCIPvarGetName(var), proofcoef > 0.0 ? "<=" : ">=", oldbound,
      SCIPvarGetName(var), proofcoef > 0.0 ? "<=" : ">=", newbound,
      proofcoef, depth, resolvable, proofactdelta, score);

   /* insert variable in candidate list without touching the already processed candidates */
   for( i = *ncands; i > firstcand && score > (*candscores)[i-1]; --i )
   {
      (*cands)[i] = (*cands)[i-1];
      (*candscores)[i] = (*candscores)[i-1];
      (*newbounds)[i] = (*newbounds)[i-1];
      (*proofactdeltas)[i] = (*proofactdeltas)[i-1];
   }
   (*cands)[i] = var;
   (*candscores)[i] = score;
   (*newbounds)[i] = newbound;
   (*proofactdeltas)[i] = proofactdelta;
   (*ncands)++;

   return SCIP_OKAY;
}

/** after changing the global bound of a variable, the bdchginfos that are now redundant are replaced with
 *  oldbound = newbound = global bound; if the current bdchginfo is of such kind, the bound is equal to the
 *  global bound and we can ignore it by installing a -1 as the corresponding bound change info position
 */
static
void skipRedundantBdchginfos(
   SCIP_VAR*             var,                /**< problem variable */
   int*                  lbchginfopos,       /**< pointer to lower bound change information position */
   int*                  ubchginfopos        /**< pointer to upper bound change information position */
   )
{
   assert(var != NULL);
   assert(lbchginfopos != NULL);
   assert(ubchginfopos != NULL);
   assert(-1 <= *lbchginfopos && *lbchginfopos <= var->nlbchginfos);
   assert(-1 <= *ubchginfopos && *ubchginfopos <= var->nubchginfos);
   assert(*lbchginfopos == -1 || *lbchginfopos == var->nlbchginfos
      || var->lbchginfos[*lbchginfopos].redundant
      == (var->lbchginfos[*lbchginfopos].oldbound == var->lbchginfos[*lbchginfopos].newbound)); /*lint !e777*/
   assert(*ubchginfopos == -1 || *ubchginfopos == var->nubchginfos
      || var->ubchginfos[*ubchginfopos].redundant
      == (var->ubchginfos[*ubchginfopos].oldbound == var->ubchginfos[*ubchginfopos].newbound)); /*lint !e777*/

   if( *lbchginfopos >= 0 && *lbchginfopos < var->nlbchginfos && var->lbchginfos[*lbchginfopos].redundant )
   {
      assert(SCIPvarGetLbGlobal(var) == var->lbchginfos[*lbchginfopos].oldbound); /*lint !e777*/
      *lbchginfopos = -1;
   }
   if( *ubchginfopos >= 0 && *ubchginfopos < var->nubchginfos && var->ubchginfos[*ubchginfopos].redundant )
   {
      assert(SCIPvarGetUbGlobal(var) == var->ubchginfos[*ubchginfopos].oldbound); /*lint !e777*/
      *ubchginfopos = -1;
   }
}

/** undoes bound changes on variables, still leaving the given infeasibility proof valid */
static
SCIP_RETCODE undoBdchgsProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            proofcoefs,         /**< coefficients in infeasibility proof */
   SCIP_Real             prooflhs,           /**< left hand side of proof */
   SCIP_Real             proofact,           /**< current activity of proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int**                 bdchginds,          /**< pointer to bound change index array, or NULL */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*                  nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   SCIP_Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** cands;
   SCIP_Real* candscores;
   SCIP_Real* newbounds;
   SCIP_Real* proofactdeltas;
   int nvars;
   int ncands;
   int candssize;
   int v;
   int i;

   assert(prob != NULL);
   assert(proofcoefs != NULL);
   assert(SCIPsetIsFeasGT(set, prooflhs, proofact));
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);

   if( resolve != NULL )
      *resolve = FALSE;

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* calculate the order in which the bound changes are tried to be undone, and relax all bounds if this doesn't
    * increase the proof's activity
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cands, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &candscores, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &newbounds, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &proofactdeltas, nvars) );
   ncands = 0;
   candssize = nvars;
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Bool relaxed;

      var = vars[v];

      /* after changing the global bound of a variable, the bdchginfos that are now redundant are replaced with
       * oldbound = newbound = global bound; if the current bdchginfo is of such kind, the bound is equal to the
       * global bound and we can ignore it
       */
      skipRedundantBdchginfos(var, &lbchginfoposs[v], &ubchginfoposs[v]);

      /* ignore variables already relaxed to global bounds */
      if( lbchginfoposs[v] == -1 && ubchginfoposs[v] == -1 )
         continue;

      /* relax bounds that are not used in the proof to the global bounds */
      relaxed = FALSE;
      if( !SCIPsetIsNegative(set, proofcoefs[v]) )
      {
         /* the lower bound is not used */
         if( lbchginfoposs[v] >= 0 )
         {
            SCIPdebugMessage(" -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g\n",
               SCIPvarGetName(var), curvarlbs[v], curvarubs[v], SCIPvarGetLbGlobal(var), curvarubs[v],
               proofcoefs[v], prooflhs, proofact);
            curvarlbs[v] = SCIPvarGetLbGlobal(var);
            lbchginfoposs[v] = -1;
            relaxed = TRUE;
         }
      }
      if( !SCIPsetIsPositive(set, proofcoefs[v]) )
      {
         /* the upper bound is not used */
         if( ubchginfoposs[v] >= 0 )
         {
            SCIPdebugMessage(" -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g\n",
               SCIPvarGetName(var), curvarlbs[v], curvarubs[v], curvarlbs[v], SCIPvarGetUbGlobal(var),
               proofcoefs[v], prooflhs, proofact);
            curvarubs[v] = SCIPvarGetUbGlobal(var);
            ubchginfoposs[v] = -1;
            relaxed = TRUE;
         }
      }
      if( relaxed && nbdchgs != NULL )
      {
         SCIP_CALL( addBdchg(set, var, curvarlbs[v], curvarubs[v],
               bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs) );
      }

      /* add bound to candidate list */
      if( lbchginfoposs[v] >= 0 || ubchginfoposs[v] >= 0 )
      {
         SCIP_CALL( addCand(set, currentdepth, var, lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v],
               prooflhs, proofact, &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, 0) );
      }
   }

   /* try to undo remaining local bound changes while still keeping the proof row violated:
    * bound changes can be undone, if prooflhs > proofact + proofactdelta;
    * afterwards, the current proof activity has to be updated
    */
   for( i = 0; i < ncands; ++i )
   {
      assert(proofactdeltas[i] > 0.0);
      assert((lbchginfoposs[SCIPvarGetProbindex(cands[i])] >= 0) != (ubchginfoposs[SCIPvarGetProbindex(cands[i])] >= 0));

      if( SCIPsetIsGT(set, prooflhs, proofact + proofactdeltas[i]) )
      {
         v = SCIPvarGetProbindex(cands[i]);
         assert(0 <= v && v < nvars);
         assert((lbchginfoposs[v] >= 0) != (ubchginfoposs[v] >= 0));

         SCIPdebugMessage(" -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g + %g\n",
            SCIPvarGetName(cands[i]), curvarlbs[v], curvarubs[v],
            proofcoefs[v] > 0.0 ? curvarlbs[v] : newbounds[i],
            proofcoefs[v] > 0.0 ? newbounds[i] : curvarubs[v],
            proofcoefs[v], prooflhs, proofact, proofactdeltas[i]);

         assert((SCIPsetIsPositive(set, proofcoefs[v]) && SCIPsetIsGT(set, newbounds[i], curvarubs[v]))
            || (SCIPsetIsNegative(set, proofcoefs[v]) && SCIPsetIsLT(set, newbounds[i], curvarlbs[v])));
         assert((SCIPsetIsPositive(set, proofcoefs[v])
               && SCIPsetIsEQ(set, proofactdeltas[i], (newbounds[i] - curvarubs[v])*proofcoefs[v]))
            || (SCIPsetIsNegative(set, proofcoefs[v])
               && SCIPsetIsEQ(set, proofactdeltas[i], (newbounds[i] - curvarlbs[v])*proofcoefs[v])));
         assert(!SCIPsetIsZero(set, proofcoefs[v]));

         if( proofcoefs[v] > 0.0 )
         {
            assert(ubchginfoposs[v] >= 0);
            assert(lbchginfoposs[v] == -1);
            curvarubs[v] = newbounds[i];
            ubchginfoposs[v]--;
         }
         else
         {
            assert(lbchginfoposs[v] >= 0);
            assert(ubchginfoposs[v] == -1);
            curvarlbs[v] = newbounds[i];
            lbchginfoposs[v]--;
         }
         if( nbdchgs != NULL )
         {
            SCIP_CALL( addBdchg(set, cands[i], curvarlbs[v], curvarubs[v],
                  bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs) );
         }
         proofact += proofactdeltas[i];
         if( resolve != NULL && SCIPvarIsInLP(cands[i]) )
            *resolve = TRUE;

         /* after changing the global bound of a variable, the bdchginfos that are now redundant are replaced with
          * oldbound = newbound = global bound; if the current bdchginfo is of such kind, the bound is equal to the
          * global bound and we can ignore it
          */
         skipRedundantBdchginfos(cands[i], &lbchginfoposs[v], &ubchginfoposs[v]);

         /* insert the new local bound of the variable into the candidate list */
         if( lbchginfoposs[v] >= 0 || ubchginfoposs[v] >= 0 )
         {
            SCIP_CALL( addCand(set, currentdepth, cands[i], lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v],
                  prooflhs, proofact, &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, i+1) );
         }
      }
   }

   /* free the buffer for the sorted bound change candidates */
   SCIPsetFreeBufferArray(set, &proofactdeltas);
   SCIPsetFreeBufferArray(set, &newbounds);
   SCIPsetFreeBufferArray(set, &candscores);
   SCIPsetFreeBufferArray(set, &cands);

   return SCIP_OKAY;
}

/* because calculations might cancel out some values, we stop the infeasibility analysis if a value is bigger than
 * 2^53 = 9007199254740992
 */
#define NUMSTOP 9007199254740992.0

/** analyzes an infeasible LP and undoes additional bound changes while staying infeasible */
static
SCIP_RETCODE undoBdchgsDualfarkas(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int**                 bdchginds,          /**< pointer to bound change index array, or NULL */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*                  nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   SCIP_Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   SCIP_Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again */
   )
{
   SCIP_RETCODE retcode;
   SCIP_LPI* lpi;
   SCIP_ROW** rows;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   SCIP_VAR* var;
   SCIP_Real* dualfarkas;
   SCIP_Real* farkascoefs;
   SCIP_Real farkaslhs;
   SCIP_Real farkasact;
   int nrows;
   int nvars;
   int r;
   int v;
   int i;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(valid != NULL);
   assert(resolve != NULL);

   SCIPdebugMessage("undoing bound changes in infeasible LP: cutoff=%g\n", lp->cutoffbound);

   *valid = FALSE;
   *resolve = FALSE;

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   vars = prob->vars;
   nvars = prob->nvars;
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* allocate temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualfarkas, nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &farkascoefs, nvars) );

   /* if solve for some reason did not produce a dual ray, e.g. because of numerical instabilities, abort conflict analysis */
   if( ! SCIPlpiHasDualRay(lpi) )
      goto TERMINATE;

   /* get dual Farkas values of rows */
   retcode = SCIPlpiGetDualfarkas(lpi, dualfarkas);
   if( retcode == SCIP_LPERROR ) /* on an error in the LP solver, just abort the conflict analysis */
      goto TERMINATE;
   SCIP_CALL( retcode );

   /* calculate the Farkas row */
   BMSclearMemoryArray(farkascoefs, nvars);
   farkaslhs = 0.0;
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore local rows and rows with Farkas value 0.0 */
      if( !row->local && !SCIPsetIsZero(set, dualfarkas[r]) )
      {
#ifndef NDEBUG
         {
            SCIP_Real lpilhs;
            SCIP_Real lpirhs;

            SCIP_CALL( SCIPlpiGetSides(lpi, r, r, &lpilhs, &lpirhs) );
            assert((SCIPsetIsInfinity(set, -lpilhs) && SCIPsetIsInfinity(set, -row->lhs))
               || SCIPsetIsRelEQ(set, lpilhs, row->lhs - row->constant));
            assert((SCIPsetIsInfinity(set, lpirhs) && SCIPsetIsInfinity(set, row->rhs))
               || SCIPsetIsRelEQ(set, lpirhs, row->rhs - row->constant));
         }
#endif

         /* add row side to Farkas row lhs: dualfarkas > 0 -> lhs, dualfarkas < 0 -> rhs */
         if( dualfarkas[r] > 0.0 )
         {
            /* check if sign of dual Farkas value is valid */
            if( SCIPsetIsInfinity(set, -row->lhs) )
               continue;

            /* due to numerical reasons we want to stop */
            if( REALABS(dualfarkas[r] * (row->lhs - row->constant)) > NUMSTOP )
               goto TERMINATE;

            farkaslhs += dualfarkas[r] * (row->lhs - row->constant);
         }
         else
         {
            /* check if sign of dual Farkas value is valid */
            if( SCIPsetIsInfinity(set, row->rhs) )
               continue;

            /* due to numerical reasons we want to stop */
            if( REALABS(dualfarkas[r] * (row->rhs - row->constant)) > NUMSTOP )
               goto TERMINATE;

            farkaslhs += dualfarkas[r] * (row->rhs - row->constant);
         }
         SCIPdebugMessage(" -> farkaslhs: %g<%s>[%g,%g] -> %g\n", dualfarkas[r], SCIProwGetName(row),
            row->lhs - row->constant, row->rhs - row->constant, farkaslhs);

         /* due to numerical reasons we want to stop */
         if( REALABS(farkaslhs) > NUMSTOP )
            goto TERMINATE;

         /* add row coefficients to Farkas row */
         for( i = 0; i < row->len; ++i )
         {
            v = SCIPvarGetProbindex(SCIPcolGetVar(row->cols[i]));
            assert(0 <= v && v < nvars);
            farkascoefs[v] += dualfarkas[r] * row->vals[i];
         }
      }
#ifdef SCIP_DEBUG
      else if( !SCIPsetIsZero(set, dualfarkas[r]) )
      {
         SCIPdebugMessage(" -> ignoring %s row <%s> with dual Farkas value %.10f (lhs=%g, rhs=%g)\n",
            row->local ? "local" : "global", SCIProwGetName(row), dualfarkas[r],
            row->lhs - row->constant, row->rhs - row->constant);
      }
#endif
   }

   /* calculate the current Farkas activity, always using the best bound w.r.t. the Farkas coefficient */
   farkasact = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(SCIPvarGetProbindex(var) == v);

      /* ignore coefficients close to 0.0 */
      if( SCIPsetIsZero(set, farkascoefs[v]) )
         farkascoefs[v] = 0.0;
      else if( farkascoefs[v] > 0.0 )
      {
         assert((SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && SCIPcolGetLPPos(SCIPvarGetCol(var)) >= 0)
            || !SCIPsetIsPositive(set, SCIPvarGetUbLP(var, set)));
         if( SCIPsetIsInfinity(set, curvarubs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarubs[v];
         SCIPdebugMessage(" -> farkasact: %g<%s>[%g,%g] -> %g\n", farkascoefs[v], SCIPvarGetName(var),
            curvarlbs[v], curvarubs[v], farkasact);
      }
      else
      {
         assert((SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && SCIPcolGetLPPos(SCIPvarGetCol(var)) >= 0)
            || !SCIPsetIsNegative(set, SCIPvarGetLbLP(var, set)));
         if( SCIPsetIsInfinity(set, -curvarlbs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarlbs[v];
         SCIPdebugMessage(" -> farkasact: %g<%s>[%g,%g] -> %g\n", farkascoefs[v], SCIPvarGetName(var),
            curvarlbs[v], curvarubs[v], farkasact);
      }
   }
   SCIPdebugMessage(" -> farkaslhs=%g, farkasact=%g\n", farkaslhs, farkasact);

   /* check, if the Farkas row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, farkaslhs, farkasact) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      SCIP_CALL( undoBdchgsProof(set, prob, currentdepth, farkascoefs, farkaslhs, farkasact,
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
            bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs, resolve) );

      *valid = TRUE;

      /* resolving does not make sense: the old dual ray is still valid -> resolving will not change the solution */
      *resolve = FALSE;
   }

 TERMINATE:

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &farkascoefs);
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** analyzes an LP exceeding the objective limit and undoes additional bound changes while staying beyond the
 *  objective limit
 */
static
SCIP_RETCODE undoBdchgsDualsol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int**                 bdchginds,          /**< pointer to bound change index array, or NULL */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*                  nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   SCIP_Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   SCIP_Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again */
   )
{
   SCIP_RETCODE retcode;
   SCIP_LPI* lpi;
   SCIP_ROW** rows;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   SCIP_VAR* var;
   SCIP_Real* primsols;
   SCIP_Real* dualsols;
   SCIP_Real* redcosts;
   SCIP_Real* dualcoefs;
   SCIP_Real* varredcosts;
   SCIP_Real duallhs;
   SCIP_Real dualact;
   int nrows;
   int ncols;
   int nvars;
   int r;
   int v;
   int i;

   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(valid != NULL);
   assert(resolve != NULL);

   *valid = FALSE;
   *resolve = FALSE;

   SCIPdebugMessage("undoing bound changes in LP exceeding cutoff: cutoff=%g\n", lp->cutoffbound);

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   ncols = SCIPlpGetNCols(lp);
   vars = prob->vars;
   nvars = prob->nvars;
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* get temporary memory */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsols, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualsols, nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redcosts, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualcoefs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varredcosts, nvars) );

   /* get solution from LPI */
   retcode = SCIPlpiGetSol(lpi, NULL, primsols, dualsols, NULL, redcosts);
   if( retcode == SCIP_LPERROR ) /* on an error in the LP solver, just abort the conflict analysis */
      goto TERMINATE;
   SCIP_CALL( retcode );
#ifdef SCIP_DEBUG
   {
      SCIP_Real objval;
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      SCIPdebugMessage(" -> LP objval: %g\n", objval);
   }
#endif

   /* Let y be the dual solution and r be the reduced cost vector. Let z be defined as
    *    z_i := y_i if i is a global row,
    *    z_i := 0   if i is a local row.
    * Define the set X := {x | lhs <= Ax <= rhs, lb <= x <= ub, c^Tx <= c*}, with c* being the current primal bound.
    * Then the following inequalities are valid for all x \in X:
    *                                 - c* <= -c^Tx
    *   <=>                     z^TAx - c* <= (z^TA - c^T) x
    *   <=>                     z^TAx - c* <= (y^TA - c^T - (y-z)^TA) x
    *   <=>                     z^TAx - c* <= (-r^T - (y-z)^TA) x         (dual feasibility of (y,r): y^TA + r^T == c^T)
    * Because lhs <= Ax <= rhs and lb <= x <= ub, the inequality can be relaxed to give
    *     min{z^Tq | lhs <= q <= rhs} - c* <= max{(-r^T - (y-z)^TA) x | lb <= x <= ub}, or X = {}.
    *
    * The resulting dual row is:  z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub},
    * where lhs, rhs, lb, and ub are selected in order to maximize the feasibility of the row.
    */

   BMSclearMemoryArray(dualcoefs, nvars);

   /* use a slightly tighter cutoff bound, because solutions with equal objective value should also be declared
    * infeasible
    */
   duallhs = -(lp->cutoffbound - SCIPsetSumepsilon(set));
   dualact = 0.0;

   /* dual row: z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub}
    * process rows: add z^T{lhs,rhs} to the dual row's left hand side, and -(y-z)^TA to the dual row's coefficients
    */
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore dual solution values of 0.0 (in this case: y_i == z_i == 0) */
      if( SCIPsetIsZero(set, dualsols[r]) )
         continue;

      /* check dual feasibility */
      if( (SCIPsetIsInfinity(set, -row->lhs) && dualsols[r] > 0.0) || (SCIPsetIsInfinity(set, row->rhs) && dualsols[r] < 0.0) )
      {
         SCIPdebugMessage(" -> infeasible dual solution %g in row <%s>: lhs=%g, rhs=%g\n",
            dualsols[r], SCIProwGetName(row), row->lhs, row->rhs);
         goto TERMINATE;
      }

      /* local rows add up to the dual row's coefficients (because z_i == 0 => -(y_i - z_i) == -y_i),
       * global rows add up to the dual row's left hand side (because z_i == y_i != 0)
       */
      if( row->local )
      {
         /* add -y_i A_i to coefficients of dual row */
         for( i = 0; i < row->len; ++i )
         {
            v = SCIPvarGetProbindex(SCIPcolGetVar(row->cols[i]));
            assert(0 <= v && v < nvars);
            dualcoefs[v] -= dualsols[r] * row->vals[i];
         }
         SCIPdebugMessage(" -> local row <%s>: dual=%g\n", SCIProwGetName(row), dualsols[r]);
      }
      else
      {
         /* add minimal value to dual row's left hand side: z_i == y_i > 0 -> lhs, z_i == y_i < 0 -> rhs */
         if( dualsols[r] > 0.0 )
         {
            assert(!SCIPsetIsInfinity(set, -row->lhs));
            duallhs += dualsols[r] * (row->lhs - row->constant);
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, row->rhs));
            duallhs += dualsols[r] * (row->rhs - row->constant);
         }
         SCIPdebugMessage(" -> global row <%s>[%g,%g]: dual=%g -> duallhs=%g\n",
            SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, dualsols[r], duallhs);
      }
   }

   /* dual row: z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub}
    * process variables: subtract reduced costs from dual row's coefficients, and calculate current maximal dual
    *                    activity by multiplying the resultant coefficient with lb or ub
    */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(SCIPvarGetProbindex(var) == v);

      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         /* reduced costs for loose variables are equal to the objective value */
         varredcosts[v] = SCIPvarGetObj(var);
      }
      else
      {
         SCIP_COL* col;
         int c;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         col = SCIPvarGetCol(var);
         c = SCIPcolGetLPPos(col);
         assert(c == -1 || col == lp->cols[c]);
         assert(c == -1 || col == lp->lpicols[c]);

         /* get reduced costs from LPI, or calculate it manually if the column is not in current LP */
         varredcosts[v] = (c >= 0 ? redcosts[c] : SCIPcolCalcRedcost(col, dualsols));

         /* check dual feasibility */
         if( (SCIPsetIsGT(set, primsols[c], curvarlbs[v]) && SCIPsetIsFeasPositive(set, varredcosts[v]))
            || (SCIPsetIsLT(set, primsols[c], curvarubs[v]) && SCIPsetIsFeasNegative(set, varredcosts[v])) )
         {
            SCIPdebugMessage(" -> infeasible reduced costs %g in var <%s>: lb=%g, ub=%g\n",
               varredcosts[v], SCIPvarGetName(var), curvarlbs[v], curvarubs[v]);
            goto TERMINATE;
         }
      }

      /* subtract reduced costs from dual row's coefficients */
      dualcoefs[v] -= varredcosts[v];

      /* add maximal value to dual row's activity: dualcoef > 0 -> ub, dualcoef < 0 -> lb */
      if( dualcoefs[v] > 0.0 )
      {
         if( SCIPsetIsInfinity(set, curvarubs[v]) )
            goto TERMINATE;
         dualact += dualcoefs[v] * curvarubs[v];
      }
      else
      {
         if( SCIPsetIsInfinity(set, -curvarlbs[v]) )
            goto TERMINATE;
         dualact += dualcoefs[v] * curvarlbs[v];
      }
   }
   SCIPdebugMessage(" -> final dual values: lhs=%g, act=%g\n", duallhs, dualact);

   /* check, if the dual row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, duallhs, dualact) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      SCIP_CALL( undoBdchgsProof(set, prob, currentdepth, dualcoefs, duallhs, dualact,
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
            bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs, resolve) );

      *valid = TRUE;
   }

 TERMINATE:

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &varredcosts);
   SCIPsetFreeBufferArray(set, &dualcoefs);
   SCIPsetFreeBufferArray(set, &redcosts);
   SCIPsetFreeBufferArray(set, &dualsols);
   SCIPsetFreeBufferArray(set, &primsols);

   return SCIP_OKAY;
}

/** applies conflict analysis starting with given bound changes, that could not be undone during previous
 *  infeasibility analysis
 */
static
SCIP_RETCODE conflictAnalyzeRemainingBdchgs(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change information in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change information in variables' arrays */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int nvars;
   int v;
   int nbdchgs;
   int maxsize;

   assert(prob != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(nconss != NULL);
   assert(nliterals != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);

   *nconss = 0;
   *nliterals = 0;
   *nreconvconss = 0;
   *nreconvliterals = 0;

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   maxsize = 2*conflictCalcMaxsize(set, prob);

   /* initialize conflict data */
   SCIP_CALL( SCIPconflictInit(conflict, set, stat, prob) );

   /* add remaining bound changes to conflict queue */
   SCIPdebugMessage("initial conflict set after undoing bound changes:\n");
   nbdchgs = 0;
   for( v = 0; v < nvars && nbdchgs < maxsize; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(var->nlbchginfos >= 0);
      assert(var->nubchginfos >= 0);
      assert(-1 <= lbchginfoposs[v] && lbchginfoposs[v] <= var->nlbchginfos);
      assert(-1 <= ubchginfoposs[v] && ubchginfoposs[v] <= var->nubchginfos);

      if( lbchginfoposs[v] == var->nlbchginfos || ubchginfoposs[v] == var->nubchginfos )
      {
         SCIP_BDCHGINFO* bdchginfo;

         /* the strong branching or diving bound stored in the column is responsible for the conflict:
          * it cannot be resolved and therefore has to be directly put into the conflict set
          */
         assert((lbchginfoposs[v] == var->nlbchginfos) != (ubchginfoposs[v] == var->nubchginfos)); /* only one can be tight in the dual! */
         assert(lbchginfoposs[v] < var->nlbchginfos || SCIPvarGetLbLP(var, set) > SCIPvarGetLbLocal(var));
         assert(ubchginfoposs[v] < var->nubchginfos || SCIPvarGetUbLP(var, set) < SCIPvarGetUbLocal(var));

         /* create an artificial bound change information for the diving/strong branching bound change;
          * they are freed in the SCIPconflictFlushConss() call
          */
         if( lbchginfoposs[v] == var->nlbchginfos )
         {
            SCIP_CALL( conflictCreateTmpBdchginfo(conflict, blkmem, set, var, SCIP_BOUNDTYPE_LOWER,
                  SCIPvarGetLbLocal(var), SCIPvarGetLbLP(var, set), &bdchginfo) );
         }
         else
         {
            SCIP_CALL( conflictCreateTmpBdchginfo(conflict, blkmem, set, var, SCIP_BOUNDTYPE_UPPER,
                  SCIPvarGetUbLocal(var), SCIPvarGetUbLP(var, set), &bdchginfo) );
         }

         /* put variable into the conflict set */
         SCIPdebugMessage("   force: <%s> %s %g [status: %d, type: %d, dive/strong]\n",
            SCIPvarGetName(var), lbchginfoposs[v] == var->nlbchginfos ? ">=" : "<=",
            lbchginfoposs[v] == var->nlbchginfos ? SCIPvarGetLbLP(var, set) : SCIPvarGetUbLP(var, set),
            SCIPvarGetStatus(var), SCIPvarGetType(var));
         SCIP_CALL( conflictAddConflictBound(conflict, blkmem, set, bdchginfo) );
         SCIP_CALL( incVSIDS(stat, var, SCIPbdchginfoGetBoundtype(bdchginfo)) );
         nbdchgs++;
      }
      else
      {
         /* put remaining bound changes into conflict candidate queue */
         if( lbchginfoposs[v] >= 0 )
         {
            assert(!SCIPbdchginfoIsRedundant(&var->lbchginfos[lbchginfoposs[v]]));
            SCIPdebugMessage("   queue: <%s> >= %g [status: %d, type: %d, depth: %d, pos:%d, chgtype: %d]\n",
               SCIPvarGetName(var), SCIPbdchginfoGetNewbound(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPvarGetStatus(var), SCIPvarGetType(var),
               SCIPbdchginfoGetDepth(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPbdchginfoGetPos(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPbdchginfoGetChgtype(&var->lbchginfos[lbchginfoposs[v]]));
            SCIP_CALL( conflictQueueBound(conflict, set, &var->lbchginfos[lbchginfoposs[v]]) );
            SCIP_CALL( incVSIDS(stat, var, SCIP_BOUNDTYPE_LOWER) );
            nbdchgs++;
         }
         if( ubchginfoposs[v] >= 0 )
         {
            assert(!SCIPbdchginfoIsRedundant(&var->ubchginfos[ubchginfoposs[v]]));
            SCIPdebugMessage("   queue: <%s> <= %g [status: %d, type: %d, depth: %d, pos:%d, chgtype: %d]\n",
               SCIPvarGetName(var), SCIPbdchginfoGetNewbound(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPvarGetStatus(var), SCIPvarGetType(var),
               SCIPbdchginfoGetDepth(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPbdchginfoGetPos(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPbdchginfoGetChgtype(&var->ubchginfos[ubchginfoposs[v]]));
            SCIP_CALL( conflictQueueBound(conflict, set, &var->ubchginfos[ubchginfoposs[v]]) );
            SCIP_CALL( incVSIDS(stat, var, SCIP_BOUNDTYPE_UPPER) );
            nbdchgs++;
         }
      }
   }

   if( v == nvars )
   {
      /* analyze the conflict set, and create conflict constraints on success */
      SCIP_CALL( conflictAnalyze(conflict, blkmem, set, stat, prob, tree, diving, 0, FALSE,
            nconss, nliterals, nreconvconss, nreconvliterals) );
   }

   return SCIP_OKAY;
}

/** actually performs analyzation of infeasible LP */
static
SCIP_RETCODE conflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int*                  iterations,         /**< pointer to store the total number of LP iterations used */
   int*                  nconss,             /**< pointer to store the number of generated conflict constraints */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict constraints */
   int*                  nreconvconss,       /**< pointer to store the number of generated reconvergence constraints */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence constraints */
   )
{
   SCIP_LPI* lpi;
   SCIP_VAR** vars;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   int nvars;
   int v;
   int* bdchginds;
   SCIP_Real* bdchgoldlbs;
   SCIP_Real* bdchgoldubs;
   SCIP_Real* bdchgnewlbs;
   SCIP_Real* bdchgnewubs;
   SCIP_Bool resolve;
   SCIP_Bool solvelp;
   int bdchgssize;
   int nbdchgs;
   int lastnbdchgs;
   int ncols;
   SCIP_Bool valid;

   assert(conflict != NULL);
   assert(conflict->nconflictsets == 0);
   assert(set != NULL);
   assert(SCIPprobAllColsInLP(prob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(iterations != NULL);
   assert(nconss != NULL);
   assert(nliterals != NULL);
   assert(nreconvconss != NULL);
   assert(nreconvliterals != NULL);

   *iterations = 0;
   *nconss = 0;
   *nliterals = 0;
   *nreconvconss = 0;
   *nreconvliterals = 0;

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));
   assert(SCIPlpiIsPrimalInfeasible(lpi) || !SCIPlpDivingObjChanged(lp));

   if( !SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_Real objval;
         
      assert(!SCIPlpDivingObjChanged(lp));

      /* make sure, a dual feasible solution exists, that exceeds the objective limit;
       * With FASTMIP setting, CPLEX does not apply the final pivot to reach the dual solution exceeding the objective
       * limit. Therefore, we have to either turn off FASTMIP and resolve the problem or continue solving it without
       * objective limit for at least one iteration. It seems that the strategy to continue with FASTMIP for one
       * additional simplex iteration yields better results.
       */
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiuobjlim )
      {
         SCIP_RETCODE retcode;

         /* temporarily disable objective limit and install an iteration limit */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, SCIPlpiInfinity(lpi)) );
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, 1) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve LP */
         retcode = SCIPlpiSolveDual(lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         valid = (retcode != SCIP_LPERROR);
         if( valid )
         {
            int iter;

            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lpi, &iter) );
            (*iterations) += iter;
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            SCIPdebugMessage(" -> resolved objlim exceeding LP in %d iterations (total: %"SCIP_LONGINT_FORMAT") (infeasible:%u, objlim: %u, optimal:%u)\n",
               iter, stat->nconflictlpiterations, SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi),
               SCIPlpiIsOptimal(lpi));
            valid = (SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsDualFeasible(lpi));
         }

         /* reinstall old objective and iteration limits in LP solver */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

         /* abort, if the LP produced an error */
         if( !valid )
            return SCIP_OKAY;
      }
   }
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsDualFeasible(lpi));

   if( !SCIPlpiIsPrimalInfeasible(lpi) )
   {
      SCIP_Real objval;

      assert(!SCIPlpDivingObjChanged(lp));

      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiuobjlim )
      {
         SCIPdebugMessage(" -> LP does not exceed the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiuobjlim);
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage(" -> LP exceeds the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiuobjlim);
      }
   }
   
   SCIPdebugMessage("analyzing conflict on infeasible LP (infeasible: %u, objlimexc: %u, optimal:%u) in depth %d (diving: %u)\n",
      SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsOptimal(lpi), SCIPtreeGetCurrentDepth(tree), diving);
#ifdef SCIP_DEBUG
   {
      SCIP_Real uobjlim;

      SCIP_CALL( SCIPlpiGetRealpar(lpi, SCIP_LPPAR_UOBJLIM, &uobjlim) );
      SCIPdebugMessage(" -> objective limit in LP solver: %g (in LP: %g)\n", uobjlim, lp->lpiuobjlim);
   }
#endif

   /* get active problem variables */
   vars = prob->vars;
   nvars = prob->nvars;

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information
    * positions in variable's bound change information arrays
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* the following algorithm is used to find a subset of changed bounds leading to an infeasible LP:
    * 1. call undoBdchgsDualfarkas() or undoBdchgsDualsol()
    *    -> update lb/ubchginfoposs arrays
    *    -> store additional changes in bdchg and curvarlbs/ubs arrays
    *    -> apply additional changes to the LPI
    * 2. (optional) if additional bound changes were undone:
    *    -> resolve LP
    *    -> goto 1.
    * 3. redo all bound changes in the LPI to restore the LPI to its original state
    * 4. analyze conflict
    *    -> put remaining changed bounds (see lb/ubchginfoposs arrays) into starting conflict set
    */

   /* get current bounds and current positions in lb/ubchginfos arrays of variables */
   valid = TRUE;
   for( v = 0; v < nvars && valid; ++v )
   {
      SCIP_VAR* var;

      var = vars[v];

      curvarlbs[v] = SCIPvarGetLbLP(var, set);
      curvarubs[v] = SCIPvarGetUbLP(var, set);
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;
      assert(diving || SCIPsetIsEQ(set, curvarlbs[v], SCIPvarGetLbLocal(var)));
      assert(diving || SCIPsetIsEQ(set, curvarubs[v], SCIPvarGetUbLocal(var)));

      /* check, if last bound changes were due to strong branching or diving */
      if( diving )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         lb = SCIPvarGetLbLocal(var);
         ub = SCIPvarGetUbLocal(var);
         if( SCIPsetIsGT(set, curvarlbs[v], lb) )
            lbchginfoposs[v] = var->nlbchginfos;
         else if( SCIPsetIsLT(set, curvarlbs[v], lb) )
         {
            /* the bound in the diving LP was relaxed -> the LP is not a subproblem of the current node -> abort! */
            /**@todo we could still analyze such a conflict, but we would have to take care with our data structures */
            valid = FALSE;
         }
         if( SCIPsetIsLT(set, curvarubs[v], ub) )
            ubchginfoposs[v] = var->nubchginfos;
         else if( SCIPsetIsGT(set, curvarubs[v], ub) )
         {
            /* the bound in the diving LP was relaxed -> the LP is not a subproblem of the current node -> abort! */
            /**@todo we could still analyze such a conflict, but we would have to take care with our data structures */
            valid = FALSE;
         }
      }
   }

   /* get number of columns in the LP */
   ncols = SCIPlpGetNCols(lp);

   /* get temporary memory for remembering bound changes on LPI columns */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchginds, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgoldlbs, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgoldubs, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgnewlbs, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgnewubs, ncols) );
   bdchgssize = ncols;
   nbdchgs = 0;
   lastnbdchgs = 0;

   /* undo as many bound changes as possible with the current LP solution */
   resolve = FALSE;
   if( valid )
   {
      int currentdepth;
      currentdepth = SCIPtreeGetCurrentDepth(tree);
      if( SCIPlpiIsPrimalInfeasible(lpi) )
      {
         SCIP_CALL( undoBdchgsDualfarkas(set, prob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
               &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
               &valid, &resolve) );
      }
      else
      {
         assert(SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi));
         SCIP_CALL( undoBdchgsDualsol(set, prob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
               &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
               &valid, &resolve) );
      }
   }

   /* check if we want to solve the LP */
   assert(SCIPprobAllColsInLP(prob, set, lp));
   solvelp = (set->conf_maxlploops != 0 && set->conf_lpiterations != 0);

   if( valid && resolve && solvelp )
   {
      SCIP_RETCODE retcode;
      SCIP_ROW** rows;
      int* sidechginds;
      SCIP_Real* sidechgoldlhss;
      SCIP_Real* sidechgoldrhss;
      SCIP_Real* sidechgnewlhss;
      SCIP_Real* sidechgnewrhss;
      SCIP_Real lpiinfinity;
      int maxlploops;
      int lpiterations;
      int sidechgssize;
      int nsidechgs;
      int nrows;
      int nloops;
      int r;

      /* get infinity value of LP solver */
      lpiinfinity = SCIPlpiInfinity(lpi);

      /* temporarily disable objective limit and install an iteration limit */
      maxlploops = (set->conf_maxlploops >= 0 ? set->conf_maxlploops : INT_MAX);
      lpiterations = (set->conf_lpiterations >= 0 ? set->conf_lpiterations : INT_MAX);
      SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lpiinfinity) );
      SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, lpiterations) );

      /* get LP rows */
      rows = SCIPlpGetRows(lp);
      nrows = SCIPlpGetNRows(lp);
      assert(nrows == 0 || rows != NULL);

      /* get temporary memory for remembering side changes on LPI rows */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechginds, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgoldlhss, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgoldrhss, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgnewlhss, nrows) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgnewrhss, nrows) );
      sidechgssize = nrows;
      nsidechgs = 0;

      /* remove all local rows by setting their sides to infinity;
       * finite sides are only changed to near infinity, such that the row's sense in the LP solver
       * is not affected (e.g. CPLEX cannot handle free rows)
       */
      for( r = 0 ; r < nrows; ++r )
      {
         assert(SCIProwGetLPPos(rows[r]) == r);

         if( SCIProwIsLocal(rows[r]) )
         {
            SCIPdebugMessage(" -> removing local row <%s> [%g,%g]\n",
               SCIProwGetName(rows[r]), SCIProwGetLhs(rows[r]), SCIProwGetRhs(rows[r]));
            SCIP_CALL( addSideRemoval(set, rows[r], lpiinfinity, &sidechginds, &sidechgoldlhss, &sidechgoldrhss,
                  &sidechgnewlhss, &sidechgnewrhss, &sidechgssize, &nsidechgs) );
         }
      }

      /* apply changes of local rows to the LP solver */
      if( nsidechgs > 0 )
      {
         SCIP_CALL( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgnewlhss, sidechgnewrhss) );
      }

      /* undo as many additional bound changes as possible by resolving the LP */
      assert(valid);
      assert(resolve);
      nloops = 0;
      while( valid && resolve && nloops < maxlploops )
      {
         int iter;

         nloops++;
         resolve = FALSE;

         SCIPdebugMessage("infeasible LP conflict analysis loop %d (changed col bounds: %d)\n", nloops, nbdchgs);

         /* apply bound changes to the LP solver */
         assert(nbdchgs >= lastnbdchgs);
         if( nbdchgs > lastnbdchgs )
         {
            SCIPdebugMessage(" -> applying %d bound changes to the LP solver (total: %d)\n", nbdchgs - lastnbdchgs, nbdchgs);
            SCIP_CALL( SCIPlpiChgBounds(lpi, nbdchgs - lastnbdchgs,
                  &bdchginds[lastnbdchgs], &bdchgnewlbs[lastnbdchgs], &bdchgnewubs[lastnbdchgs]) );
            lastnbdchgs = nbdchgs;
         }

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve LP */
         retcode = SCIPlpiSolveDual(lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode == SCIP_LPERROR )
         {
            valid = FALSE;
            break;
         }
         SCIP_CALL( retcode );

         /* count number of LP iterations */
         SCIP_CALL( SCIPlpiGetIterations(lpi, &iter) );
         (*iterations) += iter;
         stat->nconflictlps++;
         stat->nconflictlpiterations += iter;
         SCIPdebugMessage(" -> resolved LP in %d iterations (total: %"SCIP_LONGINT_FORMAT") (infeasible:%u)\n",
            iter, stat->nconflictlpiterations, SCIPlpiIsPrimalInfeasible(lpi));

         /* evaluate result */
         if( SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi) )
         {
            SCIP_Real objval;

            SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
            valid = (objval >= lp->lpiuobjlim && !SCIPlpDivingObjChanged(lp));
         }
         else
            valid = SCIPlpiIsPrimalInfeasible(lpi);

         if( valid )
         {
            int currentdepth;
            currentdepth = SCIPtreeGetCurrentDepth(tree);            
            /* undo additional bound changes */
            if( SCIPlpiIsPrimalInfeasible(lpi) )
            {
               SCIP_CALL( undoBdchgsDualfarkas(set, prob, lp, currentdepth, curvarlbs, curvarubs,
                     lbchginfoposs, ubchginfoposs,
                     &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
                     &valid, &resolve) );
            }
            else
            {
               assert(SCIPlpiIsDualFeasible(lpi) || SCIPlpiIsObjlimExc(lpi));
               SCIP_CALL( undoBdchgsDualsol(set, prob, lp, currentdepth, curvarlbs, curvarubs,
                     lbchginfoposs, ubchginfoposs,
                     &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
                     &valid, &resolve) );
            }
         }
         assert(!resolve || valid);
         assert(!resolve || nbdchgs > lastnbdchgs);
         SCIPdebugMessage(" -> finished infeasible LP conflict analysis loop %d (iter: %d, nbdchgs: %d)\n",
            nloops, iter, nbdchgs - lastnbdchgs);
      }

      SCIPdebugMessage("finished undoing bound changes after %d loops (valid=%u, nbdchgs: %d)\n",
         nloops, valid, nbdchgs);

      /* reset variables to local bounds */
      if( nbdchgs > 0 )
      {
         SCIP_CALL( SCIPlpiChgBounds(lpi, nbdchgs, bdchginds, bdchgoldlbs, bdchgoldubs) );
      }

      /* reset changes of local rows */
      if( nsidechgs > 0 )
      {
         SCIP_CALL( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgoldlhss, sidechgoldrhss) );
      }

      /* mark the LP unsolved */
      if( nbdchgs > 0 || nsidechgs > 0 )
      {
         lp->solved = FALSE;
         lp->primalfeasible = FALSE;
         lp->dualfeasible = FALSE;
         lp->lpobjval = SCIP_INVALID;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      }

      /* reinstall old objective and iteration limits in LP solver */
      SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );
      SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_LPITLIM, lp->lpiitlim) );

      /* free temporary memory */
      SCIPsetFreeBufferArray(set, &sidechgnewrhss);
      SCIPsetFreeBufferArray(set, &sidechgnewlhss);
      SCIPsetFreeBufferArray(set, &sidechgoldrhss);
      SCIPsetFreeBufferArray(set, &sidechgoldlhss);
      SCIPsetFreeBufferArray(set, &sidechginds);
   }

   /* analyze the conflict starting with remaining bound changes */
   if( valid )
   {
      SCIP_CALL( conflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, prob, tree, diving, 
            lbchginfoposs, ubchginfoposs, nconss, nliterals, nreconvconss, nreconvliterals) );
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &bdchgnewubs);
   SCIPsetFreeBufferArray(set, &bdchgnewlbs);
   SCIPsetFreeBufferArray(set, &bdchgoldubs);
   SCIPsetFreeBufferArray(set, &bdchgoldlbs);
   SCIPsetFreeBufferArray(set, &bdchginds);
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);

   /* flush conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, prob, tree) );

   return SCIP_OKAY;
}

/** analyzes an infeasible LP to find out the bound changes on variables that were responsible for the infeasibility;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible LP conflict analysis
 */
static
SCIP_RETCODE conflictAnalyzeInfeasibleLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int iterations;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(SCIPprobAllColsInLP(prob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */

   if( success != NULL )
      *success = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_useinflp )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("analyzing conflict on infeasible LP in depth %d (solstat: %d, objchanged: %u)\n",
      SCIPtreeGetCurrentDepth(tree), SCIPlpGetSolstat(lp), SCIPlpDivingObjChanged(lp));

   /* start timing */
   SCIPclockStart(conflict->inflpanalyzetime, set);
   conflict->ninflpcalls++;

   /* perform conflict analysis */
   SCIP_CALL( conflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, SCIPlpDiving(lp),
         &iterations, &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
   conflict->ninflpsuccess += (nconss > 0 ? 1 : 0);
   conflict->ninflpiterations += iterations;
   conflict->ninflpconfconss += nconss;
   conflict->ninflpconfliterals += nliterals;
   conflict->ninflpreconvconss += nreconvconss;
   conflict->ninflpreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0);

   /* stop timing */
   SCIPclockStop(conflict->inflpanalyzetime, set);

   return SCIP_OKAY;
}

/** analyzes a bound exceeding LP to find out the bound changes on variables that were responsible for exceeding the
 *  primal bound;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for bound exceeding LP conflict analysis
 */
static
SCIP_RETCODE conflictAnalyzeBoundexceedingLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int iterations;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);
   assert(!SCIPlpDivingObjChanged(lp));
   assert(SCIPprobAllColsInLP(prob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */

   if( success != NULL )
      *success = FALSE;

   /* check, if bound exceeding LP conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_useboundlp )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("analyzing conflict on bound exceeding LP in depth %d (solstat: %d)\n",
      SCIPtreeGetCurrentDepth(tree), SCIPlpGetSolstat(lp));

   /* start timing */
   SCIPclockStart(conflict->boundlpanalyzetime, set);
   conflict->nboundlpcalls++;

   /* perform conflict analysis */
   SCIP_CALL( conflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, SCIPlpDiving(lp),
         &iterations, &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
   conflict->nboundlpsuccess += (nconss > 0 ? 1 : 0);
   conflict->nboundlpiterations += iterations;
   conflict->nboundlpconfconss += nconss;
   conflict->nboundlpconfliterals += nliterals;
   conflict->nboundlpreconvconss += nreconvconss;
   conflict->nboundlpreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nconss > 0);

   /* stop timing */
   SCIPclockStop(conflict->boundlpanalyzetime, set);

   return SCIP_OKAY;
}

/** analyzes an infeasible or bound exceeding LP to find out the bound changes on variables that were responsible for the
 *  infeasibility or for exceeding the primal bound;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible or bound exceeding LP conflict analysis;
 *  may only be called if SCIPprobAllColsInLP()
 */
SCIP_RETCODE SCIPconflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   /* LP conflict analysis is only valid, if all variables are known */
   assert(SCIPprobAllColsInLP(prob, set, lp));

   /* check, if the LP was infeasible or bound exceeding */
   if( SCIPlpiIsPrimalInfeasible(SCIPlpGetLPI(lp)) )
   {
      SCIP_CALL( conflictAnalyzeInfeasibleLP(conflict, blkmem, set, stat, prob, tree, lp, success) );
   }
   else
   {
      SCIP_CALL( conflictAnalyzeBoundexceedingLP(conflict, blkmem, set, stat, prob, tree, lp, success) );
   }

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing infeasible LP conflicts */
SCIP_Real SCIPconflictGetInfeasibleLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->inflpanalyzetime);
}

/** gets number of calls to infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpcalls;
}

/** gets number of calls to infeasible LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNInfeasibleLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpsuccess;
}

/** gets number of conflict constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpconfconss;
}

/** gets total number of literals in conflict constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpconfliterals;
}

/** gets number of reconvergence constraints detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpreconvconss;
}

/** gets total number of literals in reconvergence constraints created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpreconvliterals;
}

/** gets number of LP iterations in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNInfeasibleLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ninflpiterations;
}

/** gets time in seconds used for analyzing bound exceeding LP conflicts */
SCIP_Real SCIPconflictGetBoundexceedingLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->boundlpanalyzetime);
}

/** gets number of calls to bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpcalls;
}

/** gets number of calls to bound exceeding LP conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNBoundexceedingLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpsuccess;
}

/** gets number of conflict constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpconfconss;
}

/** gets total number of literals in conflict constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpconfliterals;
}

/** gets number of reconvergence constraints detected in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpreconvconss;
}

/** gets total number of literals in reconvergence constraints created in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpreconvliterals;
}

/** gets number of LP iterations in bound exceeding LP conflict analysis */
SCIP_Longint SCIPconflictGetNBoundexceedingLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nboundlpiterations;
}




/*
 * infeasible strong branching conflict analysis
 */

/** analyses infeasible strong branching sub problems for conflicts */
SCIP_RETCODE SCIPconflictAnalyzeStrongbranch(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible downwards branch, or NULL */
   SCIP_Bool*            upconflict          /**< pointer to store whether a conflict constraint was created for an
                                              *   infeasible upwards branch, or NULL */
   )
{
   int* cstat;
   int* rstat;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   int iter;
   int nconss;
   int nliterals;
   int nreconvconss;
   int nreconvliterals;

   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPprobAllColsInLP(prob, set, lp)); /* LP conflict analysis is only valid, if all variables are known */
   assert(col != NULL);
   assert((col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound)
         && SCIPsetFeasCeil(set, col->primsol-1.0) >= col->lb - 0.5)
      || (col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound)
         && SCIPsetFeasFloor(set, col->primsol+1.0) <= col->ub + 0.5));

   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_usesb )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->sbanalyzetime, set);

   /* get temporary memory for storing current LP basis */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, lp->nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, lp->nlpirows) );

   /* get current LP basis */
   SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   /* remember old bounds */
   oldlb = col->lb;
   oldub = col->ub;

   /* is down branch infeasible? */
   if( col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound) )
   {
      newub = SCIPsetFeasCeil(set, col->primsol-1.0);
      if( newub >= col->lb - 0.5 )
      {
         SCIP_RETCODE retcode;

         SCIPdebugMessage("analyzing conflict on infeasible downwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(SCIPcolGetVar(col)), SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)), 
            SCIPtreeGetCurrentDepth(tree));

         conflict->nsbcalls++;

         /* change the upper bound */
         col->ub = newub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPdebugMessage(" -> resolved downwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP */
            SCIP_CALL( conflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, TRUE,
                  &iter, &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
            conflict->nsbsuccess += (nconss > 0 ? 1 : 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfconss += nconss;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvconss += nreconvconss;
            conflict->nsbreconvliterals += nreconvliterals;
            if( downconflict != NULL )
               *downconflict = (nconss > 0);
         }

         /* reset the upper bound */
         col->ub = oldub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );

         /* mark the LP unsolved */
         lp->solved = FALSE;
         lp->primalfeasible = FALSE;
         lp->dualfeasible = FALSE;
         lp->lpobjval = SCIP_INVALID;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      }
   }

   /* is up branch infeasible? */
   if( col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound) )
   {
      newlb = SCIPsetFeasFloor(set, col->primsol+1.0);
      if( newlb <= col->ub + 0.5 )
      {
         SCIP_RETCODE retcode;

         SCIPdebugMessage("analyzing conflict on infeasible upwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(SCIPcolGetVar(col)), SCIPvarGetLbLocal(SCIPcolGetVar(col)), SCIPvarGetUbLocal(SCIPcolGetVar(col)), 
            SCIPtreeGetCurrentDepth(tree));

         conflict->nsbcalls++;

         /* change the lower bound */
         col->lb = newlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPdebugMessage(" -> resolved upwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP */
            SCIP_CALL( conflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, TRUE,
                  &iter, &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
            conflict->nsbsuccess += (nconss > 0 ? 1 : 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfconss += nconss;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvconss += nreconvconss;
            conflict->nsbreconvliterals += nreconvliterals;
            if( upconflict != NULL )
               *upconflict = (nconss > 0);
         }

         /* reset the lower bound */
         col->lb = oldlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );

         /* mark the LP unsolved */
         lp->solved = FALSE;
         lp->primalfeasible = FALSE;
         lp->dualfeasible = FALSE;
         lp->lpobjval = SCIP_INVALID;
         lp->lpsolstat = SCIP_LPSOLSTAT_NOTSOLVED;
      }
   }

   /* free temporary memory for storing current LP basis */
   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);

   /* stop timing */
   SCIPclockStop(conflict->sbanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing infeasible strong branching conflicts */
SCIP_Real SCIPconflictGetStrongbranchTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->sbanalyzetime);
}

/** gets number of calls to infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbcalls;
}

/** gets number of calls to infeasible strong branching conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNStrongbranchSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbsuccess;
}

/** gets number of conflict constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfconss;
}

/** gets total number of literals in conflict constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfliterals;
}

/** gets number of reconvergence constraints detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvconss;
}

/** gets total number of literals in reconvergence constraints created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvliterals;
}

/** gets number of LP iterations in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbiterations;
}




/*
 * pseudo solution conflict analysis
 */

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on
 *  variables that were responsible for the objective value degradation;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for pseudo solution conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyzePseudo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   SCIP_Real* pseudocoefs;
   SCIP_Real pseudolhs;
   SCIP_Real pseudoact;
   int nvars;
   int v;

   assert(conflict != NULL);
   assert(conflict->nconflictsets == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(!SCIPsetIsInfinity(set, -SCIPlpGetPseudoObjval(lp, set)));
   assert(!SCIPsetIsInfinity(set, lp->cutoffbound));

   if( success != NULL )
      *success = FALSE;

   /* check, if pseudo solution conflict analysis is enabled */
   if( !set->conf_enable || !set->conf_usepseudo )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("analyzing pseudo solution (obj: %g) that exceeds objective limit (%g)\n",
      SCIPlpGetPseudoObjval(lp, set), lp->cutoffbound);

   /* start timing */
   SCIPclockStart(conflict->pseudoanalyzetime, set);
   conflict->npseudocalls++;

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* The current primal bound c* gives an upper bound for the current pseudo objective value:
    *   min{c^T x | lb <= x <= ub} <= c*.
    * We have to transform this row into a >= inequality in order to use methods above:
    *                          -c* <= max{-c^T x | lb <= x <= ub}.
    * In the local subproblem, this row is violated. We want to undo bound changes while still keeping the
    * row violated.
    */

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information
    * positions in variable's bound change information arrays
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* get temporary memory for infeasibility proof coefficients */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &pseudocoefs, nvars) );

   /* use a slightly tighter cutoff bound, because solutions with equal objective value should also be declared
    * infeasible
    */
   pseudolhs = -(lp->cutoffbound - SCIPsetSumepsilon(set));

   /* store the objective values as infeasibility proof coefficients, and recalculate the pseudo activity */
   pseudoact = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      pseudocoefs[v] = -SCIPvarGetObj(var);
      curvarlbs[v] = SCIPvarGetLbLocal(var);
      curvarubs[v] = SCIPvarGetUbLocal(var);
      if( pseudocoefs[v] > 0.0 )
         pseudoact += pseudocoefs[v] * curvarubs[v];
      else
         pseudoact += pseudocoefs[v] * curvarlbs[v];
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;
   }
   assert(SCIPsetIsFeasEQ(set, pseudoact, -SCIPlpGetPseudoObjval(lp, set)));
   SCIPdebugMessage("  -> recalculated pseudo infeasibility proof:  %g <= %g\n", pseudolhs, pseudoact);

   /* check, if the pseudo row is still violated (after recalculation of pseudo activity) */
   if( SCIPsetIsFeasGT(set, pseudolhs, pseudoact) )
   {
      int nconss;
      int nliterals;
      int nreconvconss;
      int nreconvliterals;

      /* undo bound changes without destroying the infeasibility proof */
      SCIP_CALL( undoBdchgsProof(set, prob, SCIPtreeGetCurrentDepth(tree), pseudocoefs, pseudolhs, pseudoact,
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* analyze conflict on remaining bound changes */
      SCIP_CALL( conflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, prob, tree, FALSE,
            lbchginfoposs, ubchginfoposs, &nconss, &nliterals, &nreconvconss, &nreconvliterals) );
      conflict->npseudosuccess += (nconss > 0 ? 1 : 0);
      conflict->npseudoconfconss += nconss;
      conflict->npseudoconfliterals += nliterals;
      conflict->npseudoreconvconss += nreconvconss;
      conflict->npseudoreconvliterals += nreconvliterals;
      if( success != NULL )
         *success = (nconss > 0);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &pseudocoefs);
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);

   /* flush conflict set storage */
   SCIP_CALL( SCIPconflictFlushConss(conflict, blkmem, set, stat, prob, tree) );

   /* stop timing */
   SCIPclockStop(conflict->pseudoanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing pseudo solution conflicts */
SCIP_Real SCIPconflictGetPseudoTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->pseudoanalyzetime);
}

/** gets number of calls to pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudocalls;
}

/** gets number of calls to pseudo solution conflict analysis that yield at least one conflict constraint */
SCIP_Longint SCIPconflictGetNPseudoSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudosuccess;
}

/** gets number of conflict constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfconss;
}

/** gets total number of literals in conflict constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfliterals;
}

/** gets number of reconvergence constraints detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceConss(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvconss;
}

/** gets total number of literals in reconvergence constraints created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvliterals;
}
