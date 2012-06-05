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

/**@file   relax.h
 * @brief  internal methods for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAX_H__
#define __SCIP_RELAX_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_relax.h"
#include "scip/pub_relax.h"

#ifdef __cplusplus
extern "C" {
#endif

/** copies the given relaxator to a new scip */
extern
SCIP_RETCODE SCIPrelaxCopyInclude(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a relaxator */
extern
SCIP_RETCODE SCIPrelaxCreate(
   SCIP_RELAX**          relax,              /**< pointer to relaxator data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of relaxator */
   const char*           desc,               /**< description of relaxator */
   int                   priority,           /**< priority of the relaxator (negative: after LP, non-negative: before LP) */
   int                   freq,               /**< frequency for calling relaxator */
   SCIP_DECL_RELAXCOPY   ((*relaxcopy)),     /**< copy method of relaxator or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_RELAXFREE   ((*relaxfree)),     /**< destructor of relaxator */
   SCIP_DECL_RELAXINIT   ((*relaxinit)),     /**< initialize relaxator */
   SCIP_DECL_RELAXEXIT   ((*relaxexit)),     /**< deinitialize relaxator */
   SCIP_DECL_RELAXINITSOL((*relaxinitsol)),  /**< solving process initialization method of relaxator */
   SCIP_DECL_RELAXEXITSOL((*relaxexitsol)),  /**< solving process deinitialization method of relaxator */
   SCIP_DECL_RELAXEXEC   ((*relaxexec)),     /**< execution method of relaxator */
   SCIP_RELAXDATA*       relaxdata           /**< relaxator data */
   );

/** calls destructor and frees memory of relaxator */
extern
SCIP_RETCODE SCIPrelaxFree(
   SCIP_RELAX**          relax,              /**< pointer to relaxator data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes relaxator */
extern
SCIP_RETCODE SCIPrelaxInit(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of relaxator */
extern
SCIP_RETCODE SCIPrelaxExit(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs relaxator that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPrelaxInitsol(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs relaxator that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPrelaxExitsol(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls execution method of relaxator */
extern
SCIP_RETCODE SCIPrelaxExec(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   int                   depth,              /**< depth of current node */
   SCIP_Real*            lowerbound,         /**< pointer to lower bound computed by the relaxator */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of relaxator */
extern
void SCIPrelaxSetPriority(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the relaxator */
   );

/** returns whether the relaxation was completely solved at the current node */
extern
SCIP_Bool SCIPrelaxIsSolved(
   SCIP_RELAX*           relax,              /**< relaxator */
   SCIP_STAT*            stat                /**< dynamic problem statistics */
   );

/* 
 *  methods for the global relaxation data 
 */

/** creates global relaxation data */
extern
SCIP_RETCODE SCIPrelaxationCreate(
   SCIP_RELAXATION**     relaxation          /**< global relaxation data */
   );

/** frees global relaxation data */
extern
SCIP_RETCODE SCIPrelaxationFree(
   SCIP_RELAXATION**     relaxation          /**< global relaxation data */
   );

/** sets the relaxsolzero flag in the relaxation data to the given value */
extern
void SCIPrelaxationSetSolZero(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Bool             iszero              /**< are all values of the relaxation solution set to zero? */
   );

/** returns whether the global relaxation solution is cleared and all values are set to zero */
extern
SCIP_Bool SCIPrelaxationIsSolZero(
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   );

/** sets the relaxsolvalid flag in the relaxation data to the given value */
extern
void SCIPrelaxationSetSolValid(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Bool             isvalid             /**< is the stored solution valid? */
   );

/** returns whether the global relaxation solution is valid */
extern
SCIP_Bool SCIPrelaxationIsSolValid(
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   );

/** sets the objective value of the global relaxation solution */
extern
void SCIPrelaxationSetSolObj(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Real             obj                 /**< objective value */
   );

/** returns the objective value of the global relaxation solution w.r.t. the transformed problem */
extern
SCIP_Real SCIPrelaxationGetSolObj(
   SCIP_RELAXATION*      relaxation          /**< global relaxation data */
   );

/** adds the given value to the global relaxation solution's objective value */
extern
void SCIPrelaxationSolObjAdd(
   SCIP_RELAXATION*      relaxation,         /**< global relaxation data */
   SCIP_Real             val                 /**< value to add to the objective value */
   );

#ifdef __cplusplus
}
#endif

#endif
