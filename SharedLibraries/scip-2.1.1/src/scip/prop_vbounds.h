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

/**@file   prop_vbounds.h
 * @ingroup PROPAGATORS
 * @brief  variable upper and lower bound propagator
 * @author Stefan Heinz
 * @author Jens Schulz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_VBOUNDS_H__
#define __SCIP_PROP_VBOUNDS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the vbounds propagator and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePropVbounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** create a topological sorted variable array of the given variables and stores if (needed) the involved variables into
 *  the corresponding variable array and hash map
 *
 * @note: for all arrays and the hash map (if needed) you need to allocate enough memory before calling this method 
 */
extern
SCIP_RETCODE SCIPcreateTopoSortedVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variable which we want sort */
   int                   nvars,              /**< number of variables */
   SCIP_HASHMAP*         varHashmap,         /**< mapping a variable to its position in the (used) variable array, or NULL */    
   SCIP_VAR**            usedvars,           /**< array of variables which are involved in the propagation, or NULL */
   int*                  nusedvars,          /**< number of variables which are involved in the propagation, or NULL */
   SCIP_VAR**            topovars,           /**< array where the topological sorted variables are stored */
   int*                  ntopovars,          /**< pointer to store the number of topological sorted variables */
   SCIP_Bool             lowerbound          /**< topological sorted with respect to the variable lower bounds, otherwise variable upper bound */
   );

/** returns TRUE if the propagator has the status that all variable lower and upper bounds are propagated */
extern
SCIP_Bool SCIPisPropagatedVbounds(
   SCIP*                 scip                 /**< SCIP data structure */
   );

/** performs propagation of variables lower and upper bounds */
extern
SCIP_RETCODE SCIPexecPropVbounds(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_Bool             force,               /**< should domain changes be forced */
   SCIP_RESULT*          result               /**< pointer to store result */
   );

#ifdef __cplusplus
}
#endif

#endif
