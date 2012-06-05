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

/**@file   presol.h
 * @brief  internal methods for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_H__
#define __SCIP_PRESOL_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_presol.h"
#include "scip/pub_presol.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given presolver to a new scip */
extern
SCIP_RETCODE SCIPpresolCopyInclude(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a presolver */
extern
SCIP_RETCODE SCIPpresolCreate(
   SCIP_PRESOL**         presol,             /**< pointer to store presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
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
   );

/** frees memory of presolver */
extern
SCIP_RETCODE SCIPpresolFree(
   SCIP_PRESOL**         presol,             /**< pointer to presolver data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes presolver */
extern
SCIP_RETCODE SCIPpresolInit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** deinitializes presolver */
extern
SCIP_RETCODE SCIPpresolExit(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs presolver that the presolving process is being started */
extern
SCIP_RETCODE SCIPpresolInitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             isunbounded,        /**< was unboundedness already detected */
   SCIP_Bool             isinfeasible,       /**< was infeasibility already detected */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** informs presolver that the presolving process is finished */
extern
SCIP_RETCODE SCIPpresolExitpre(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             isunbounded,        /**< was unboundedness already detected */
   SCIP_Bool             isinfeasible,       /**< was infeasibility already detected */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** executes presolver */
extern
SCIP_RETCODE SCIPpresolExec(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_Bool             execdelayed,        /**< execute presolver even if it is marked to be delayed */
   int                   nrounds,            /**< number of presolving rounds already done */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*                  nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*                  nchgbds,            /**< pointer to total number of variable bounds tightened of all presolvers */
   int*                  naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*                  naddconss,          /**< pointer to total number of added constraints of all presolvers */
   int*                  nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*                  nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*                  nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of presolver */
extern
void SCIPpresolSetPriority(
   SCIP_PRESOL*          presol,             /**< presolver */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the presolver */
   );

#ifdef __cplusplus
}
#endif

#endif
