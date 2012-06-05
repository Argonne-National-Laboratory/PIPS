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

/**@file   cons_orbitope.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for (partitioning/packing) orbitope constraints w.r.t. the full symmetric group
 * @author Timo Berthold
 * @author Marc Pfetsch
 *
 * This constraint handler can be used to handle symmetries in certain 0/1-programs. The principle
 * structure is that some variables can be ordered in matrix form, such that permuting columns does
 * not change the validity and objective function value of a solution. That is, the symmetry group
 * of the program contains the full symmetric group obtained by permuting the columns of this
 * matrix. The variables in each row have to be contained in set packing or partitioning
 * constraints.
 *
 * In more mathematical terms the structure has to be as follows: There are 0/1-variables
 * \f$x_{ij}\f$, \f$i \in \{1, \dots, p\}\f$, \f$j \in \{1, \dots, q\}\f$. The variables are coupled
 * through set packing or partitioning constraints:
 * \f[
 *    \sum_{j = 1}^q x_{ij} \leq 1  \quad \mbox{or} \quad \sum_{j = 1}^q x_{ij} = 1 \quad \mbox{for all }i = 1, \ldots, p.
 * \f]
 * Permuting columns of \f$x\f$ does not change the validity and objective function value of any feasible solution.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_ORBITOPE_H__
#define __SCIP_CONS_ORBITOPE_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the handler for orbitope constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrOrbitope(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a orbitope constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
extern
SCIP_RETCODE SCIPcreateConsOrbitope(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR***           vars,               /**< matrix of variables on which the symmetry acts */
   SCIP_Bool             ispart,             /**< whether we deal with the partitioning case (packing otherwise) */
   int                   nspcons,            /**< number of set partitioning/packing constraints  <=> p */
   int                   nblocks,            /**< number of symmetric variable blocks             <=> q */
   SCIP_Bool             resolveprop,        /**< should propagation be resolved? */
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

#ifdef __cplusplus
}
#endif

#endif
