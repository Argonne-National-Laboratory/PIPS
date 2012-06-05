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

/**@file   struct_scip.h
 * @brief  SCIP main data structure
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SCIP_H__
#define __SCIP_STRUCT_SCIP_H__


#include "scip/def.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_clock.h"
#include "scip/type_interrupt.h"
#include "scip/type_mem.h"
#include "scip/type_event.h"
#include "scip/type_lp.h"
#include "scip/type_nlp.h"
#include "scip/type_implics.h"
#include "scip/type_prob.h"
#include "scip/type_primal.h"
#include "scip/type_tree.h"
#include "scip/type_pricestore.h"
#include "scip/type_sepastore.h"
#include "scip/type_cutpool.h"
#include "scip/type_branch.h"
#include "scip/type_conflict.h"
#include "scip/type_dialog.h"

#ifdef __cplusplus
extern "C" {
#endif

/** SCIP main data structure */
struct Scip
{
   /* INIT */
   SCIP_MEM*             mem;                /**< block memory buffers */
   SCIP_SET*             set;                /**< global SCIP settings */
   SCIP_INTERRUPT*       interrupt;          /**< CTRL-C interrupt data */
   SCIP_DIALOGHDLR*      dialoghdlr;         /**< dialog handler for user interface */
   SCIP_CLOCK*           totaltime;          /**< total SCIP running time */

   /* PROBLEM */
   SCIP_STAT*            stat;               /**< dynamic problem statistics */
   SCIP_PROB*            origprob;           /**< original problem data */
   SCIP_PRIMAL*          origprimal;         /**< primal data and solution storage for solution candidates */

   /* TRANSFORMED */
   SCIP_EVENTFILTER*     eventfilter;        /**< event filter for global (not variable dependent) events */
   SCIP_EVENTQUEUE*      eventqueue;         /**< event queue to cache events and process them later (bound change events) */
   SCIP_BRANCHCAND*      branchcand;         /**< storage for branching candidates */
   SCIP_LP*              lp;                 /**< LP data */
   SCIP_NLP*             nlp;                /**< NLP data */
   SCIP_RELAXATION*      relaxation;         /**< global relaxation data */
   SCIP_PRIMAL*          primal;             /**< primal data and solution storage */
   SCIP_TREE*            tree;               /**< branch and bound tree */
   SCIP_CONFLICT*        conflict;           /**< conflict analysis data */
   SCIP_CLIQUETABLE*     cliquetable;        /**< collection of cliques */
   SCIP_PROB*            transprob;          /**< transformed problem after presolve */

   /* SOLVING */
   SCIP_PRICESTORE*      pricestore;         /**< storage for priced variables */
   SCIP_SEPASTORE*       sepastore;          /**< storage for separated cuts */
   SCIP_CUTPOOL*         cutpool;            /**< global cut pool */
};

#ifdef __cplusplus
}
#endif

#endif
