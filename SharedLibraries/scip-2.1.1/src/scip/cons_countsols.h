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

/**@file   cons_countsols.h
 * @ingroup CONSHDLRS
 * @brief  Constraint handler for counting feasible solutions
 * @author Stefan Heinz
 * @author Michael Winkler
 *
 * If this constraint handler is activated than it counts or collects all feasible solutions. We refer to \ref COUNTER for
 * more details about using SCIP for counting feasible solutions.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_COUNTSOLS_H__
#define __SCIP_CONS_COUNTSOLS_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for sparse solutions */ 
struct SparseSolution
{
   SCIP_Longint*         lbvalues;           /**< array of lower bounds */
   SCIP_Longint*         ubvalues;           /**< array of upper bounds */
};
typedef struct SparseSolution SPARSESOLUTION;


/** dialog execution method for the count command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCountPresolve);

/** dialog execution method for the count command */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecCount);
   
/** execution method of dialog for writing all solutions */
extern
SCIP_DECL_DIALOGEXEC(SCIPdialogExecWriteAllsolutions);

/** creates the handler for countsol constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrCountsols(
   SCIP*                 scip                /**< SCIP data structure */
   );

/* execute counting */
extern
SCIP_RETCODE SCIPcount(
   SCIP*                 scip                /**< SCIP data structure */
   );

#if 0
/* returns TRUE if the counting process was correct; otherwise FALSE */
extern
SCIP_Bool SCIPisCountValid(
   SCIP*                 scip                /**< SCIP data structure */
   );
#endif

/** returns number of feasible solutions found as SCIP_Longint; if the number does not fit into
 *  a SCIP_Longint the valid flag is set to FALSE
 */
extern
SCIP_Longint SCIPgetNCountedSols(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_Bool*            valid              /**< pointer to store if the return value is valid */
   );

/** returns number of counted solutions as string */
extern
void SCIPgetNCountedSolsstr(
   SCIP*                 scip,              /**< SCIP data structure */
   char**                buffer,             /**< buffer to store the number for counted solutions */
   int                   buffersize,         /**< buffer size */
   int*                  requiredsize        /**< pointer to store the required size */
   );

/** returns number of counted feasible subtrees */
extern
SCIP_Longint SCIPgetNCountedFeasSubtrees(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** Method to get the sparse solution.
 *
 *  @note You get the pointer to the sparse solutions stored in the constraint handler (not a copy).
 *
 *  @note The sparse solutions are stored w.r.t. the active variables. This are the variables which got not removed
 *        during presolving. For none active variables the value has to be computed depending on their aggregation
 *        type. See for more details about that \ref COLLECTALLFEASEBLES.
 */
extern
void SCIPgetCountedSparseSolutions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< pointer to variable array defining to variable order */
   int*                  nvars,              /**< number of variables */
   SPARSESOLUTION***     sols,               /**< pointer to the solutions */
   int*                  nsols               /**< pointer to number of solutions */
   );

/** setting SCIP parameters for such that a valid counting process is possible */
extern
SCIP_RETCODE SCIPsetParamsCountsols(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
