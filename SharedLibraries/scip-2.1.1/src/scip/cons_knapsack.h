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

/**@file   cons_knapsack.h
 * @ingroup CONSHDLRS
 * @brief  Constraint handler for knapsack constraints of the form  \f$a^T x \le b\f$, x binary and \f$a \ge 0\f$.
 * @author Tobias Achterberg
 * @author Kati Wolter
 * @author Michael Winkler
 *
 * This constraint handler handles a special type of linear constraints, namely knapsack constraints.
 * A knapsack constraint has the form
 * \f[
 *   \sum_{i=1}^n a_i x_i \leq b
 * \f]
 * with non-negative integer coefficients \f$a_i\f$, integer right-hand side \f$b\f$, and binary variables \f$x_i\f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_KNAPSACK_H__
#define __SCIP_CONS_KNAPSACK_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for knapsack constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrKnapsack(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a knapsack constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
extern
SCIP_RETCODE SCIPcreateConsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of items in the knapsack */
   SCIP_VAR**            vars,               /**< array with item variables */
   SCIP_Longint*         weights,            /**< array with item weights */
   SCIP_Longint          capacity,           /**< capacity of knapsack (right hand side of inequality) */
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

/** adds new item to knapsack constraint */
extern
SCIP_RETCODE SCIPaddCoefKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint data */
   SCIP_VAR*             var,                /**< item variable */
   SCIP_Longint          weight              /**< item weight */
   );

/** gets the capacity of the knapsack constraint */
extern
SCIP_Longint SCIPgetCapacityKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the number of items in the knapsack constraint */
extern
int SCIPgetNVarsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the array of variables in the knapsack constraint; the user must not modify this array! */
extern
SCIP_VAR** SCIPgetVarsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the array of weights in the knapsack constraint; the user must not modify this array! */
extern
SCIP_Longint* SCIPgetWeightsKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the dual solution of the knapsack constraint in the current LP */
extern
SCIP_Real SCIPgetDualsolKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** gets the dual Farkas value of the knapsack constraint in the current infeasible LP */
extern
SCIP_Real SCIPgetDualfarkasKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** returns the linear relaxation of the given knapsack constraint; may return NULL if no LP row was yet created;
 *  the user must not modify the row!
 */
extern
SCIP_ROW* SCIPgetRowKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   );

/** solves knapsack problem in maximization form exactly using dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 */
extern
SCIP_RETCODE SCIPsolveKnapsackExactly(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval,             /**< pointer to store optimal solution value, or NULL */
   SCIP_Bool*            success             /**< pointer to store if an error occured during solving (normally a memory problem) */
   );

/** solves knapsack problem in maximization form approximately by solving the LP-relaxation of the problem using Dantzig's
 *  method and rounding down the solution; if needed, one can provide arrays to store all selected items and all not 
 *  selected items
 */
SCIP_RETCODE SCIPsolveKnapsackApproximately(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   );

/** separates lifted valid inequalities for given knapsack problem */
SCIP_RETCODE SCIPseparateKnapsackCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_SOL*             sol,                /**< primal CIP solution to separate, NULL for current LP solution */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   );

/** solves knapsack problem with dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 */
extern
SCIP_RETCODE SCIPsolveKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Longint*         weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers, or NULL */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   );

/** lifts given cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
 *  polytope by using uplifting for all variables not in the cover and downlifting for all variables in the cover that 
 *  are fixed to one (C2)
 */
extern
SCIP_RETCODE SCIPliftKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_Real*            solvals,            /**< solution values of all problem variables */
   int*                  covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*                  noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int                   ncovervars,         /**< number of cover variables */
   int                   ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int                   ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int                   nnoncovervars,      /**< number of noncover variables */
   int*                  liftcoefs,          /**< pointer to store lifting coefficient of variables in knapsack constraint */
   int*                  liftrhs,            /**< pointer to store right hand side of the lifted cover inequality */
   SCIP_Real*            liftlpval           /**< pointer to store LP solution value of lifted variables */  
   );

/** separates lifted cover inequalities for given knapsack problem */
extern
SCIP_RETCODE SCIPseparateKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the knapsack problem, or NULL */
   SCIP_VAR**            vars,               /**< variables in knapsack constraint */
   int                   nvars,              /**< number of variables in knapsack constraint */
   SCIP_Longint*         weights,            /**< weights of variables in knapsack constraint */
   SCIP_Longint          capacity,           /**< capacity of knapsack */
   SCIP_SOL*             sol,                /**< primal CIP solution to separate, NULL for current LP solution */
   int                   maxnumcardlift,     /**< maximal number of cardinality inequalities lifted per sepa round (-1: unlimited) */
   int*                  ncuts               /**< pointer to add up the number of found cuts */
   );

/* relaxes given general linear constraint into a knapsack constraint and separates lifted knapsack cover inequalities */
extern
SCIP_RETCODE SCIPseparateRelaxedKnapsack(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint that originates the linear constraint, or NULL */
   int                   nknapvars,          /**< number of variables in the continuous knapsack constraint */
   SCIP_VAR**            knapvars,           /**< variables in the continuous knapsack constraint */
   SCIP_Real*            knapvals,           /**< coefficients of the variables in the continuous knapsack constraint */
   SCIP_Real             valscale,           /**< -1.0 if lhs of row is used as rhs of c. k. constraint, +1.0 otherwise */
   SCIP_Real             rhs,                /**< right hand side of the continuous knapsack constraint */
   SCIP_SOL*             sol,                /**< primal CIP solution, NULL for current LP solution */
   int*                  ncuts,              /**< pointer to add up the number of found cuts */
   SCIP_Bool*            cutoff              /**< pointer to store whether a cutoff was found */
   );

#ifdef __cplusplus
}
#endif

#endif
