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

/**@file   cons_linking.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for linking binary variables to an integer variable
 * @author Stefan Heinz
 * @author Jens Schulz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_LINKING_H__
#define __SCIP_CONS_LINKING_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for linking constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrLinking(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a linking constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
extern
SCIP_RETCODE SCIPcreateConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             intvar,             /**< integer variable which should be linked */
   SCIP_VAR**            binvars,            /**< binary variables */
   int                   nbinvars,           /**< number of binary starting variables */
   int                   offset,             /**< offset of the binary variable representation */
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
   );


/** checks if for the given integer variable a linking constraint exists */
extern
SCIP_Bool SCIPexistsConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             intvar              /**< integer variable which should be linked */
   );

/** returns the linking constraint belonging the given integer variable or NULL if it does not exist yet */
extern
SCIP_CONS* SCIPgetConsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             intvar              /**< integer variable which should be linked */
   );

/** returns the integer variable of the linking constraint */
extern
SCIP_VAR* SCIPgetIntvarLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/** returns the binary variables of the linking constraint */
extern
SCIP_RETCODE SCIPgetBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< linking constraint */
   SCIP_VAR***           binvars,            /**< pointer to store the binary variables array pointer */
   int*                  nbinvars            /**< pointer to store the number of returned binary variables */
   );

/** returns the number of binary variables of the linking constraint */
extern
int SCIPgetNBinvarsLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

/** returns the offset of the linking constraint */
extern
int SCIPgetOffsetLinking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< linking constraint */
   );

#ifdef __cplusplus
}
#endif

#endif
