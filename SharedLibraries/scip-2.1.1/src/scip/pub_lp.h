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

/**@file   pub_lp.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for LP management
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_LP_H__
#define __SCIP_PUB_LP_H__


#include <stdio.h>

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_lpi.h"
#include "scip/type_lp.h"
#include "scip/type_var.h"
#include "scip/type_sol.h"

#ifdef NDEBUG
#include "scip/struct_lp.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Column methods
 */

/** output column to file stream */
extern
void SCIPcolPrint(
   SCIP_COL*             col,                /**< LP column */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** sorts column entries such that LP rows precede non-LP rows and inside both parts lower row indices precede higher ones
 */
extern
void SCIPcolSort(
   SCIP_COL*             col                 /**< column to be sorted */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** gets objective value of column */
extern
SCIP_Real SCIPcolGetObj(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets lower bound of column */
extern
SCIP_Real SCIPcolGetLb(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets upper bound of column */
extern
SCIP_Real SCIPcolGetUb(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets best bound of column with respect to the objective function */
extern
SCIP_Real SCIPcolGetBestBound(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the primal LP solution of a column */
extern
SCIP_Real SCIPcolGetPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the minimal LP solution value, this column ever assumed */
extern
SCIP_Real SCIPcolGetMinPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the maximal LP solution value, this column ever assumed */
extern
SCIP_Real SCIPcolGetMaxPrimsol(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets the basis status of a column in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_ZERO for columns not in the current SCIP_LP
 */
extern
SCIP_BASESTAT SCIPcolGetBasisStatus(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets variable this column represents */
extern
SCIP_VAR* SCIPcolGetVar(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets unique index of col */
extern
int SCIPcolGetIndex(
   SCIP_COL*             col                 /**< LP col */
   );

/** returns whether the associated variable is of integral type (binary, integer, implicit integer) */
extern
SCIP_Bool SCIPcolIsIntegral(
   SCIP_COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is removable from the LP (due to aging or cleanup) */
extern
SCIP_Bool SCIPcolIsRemovable(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets position of column in current LP, or -1 if it is not in LP */
extern
int SCIPcolGetLPPos(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets depth in the tree where the column entered the LP, or -1 if it is not in LP */
extern
int SCIPcolGetLPDepth(
   SCIP_COL*             col                 /**< LP column */
   );

/** returns TRUE iff column is member of current LP */
extern
SCIP_Bool SCIPcolIsInLP(
   SCIP_COL*             col                 /**< LP column */
   );

/** get number of nonzero entries in column vector */
extern
int SCIPcolGetNNonz(
   SCIP_COL*             col                 /**< LP column */
   );

/** get number of nonzero entries in column vector, that correspond to rows currently in the SCIP_LP;
 *  Warning! This method is only applicable on columns, that are completely linked to their rows (e.g. a column
 *  that is in the current LP and the LP was solved, or a column that was in a solved LP and didn't change afterwards
 */
extern
int SCIPcolGetNLPNonz(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets array with rows of nonzero entries */
extern
SCIP_ROW** SCIPcolGetRows(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets array with coefficients of nonzero entries */
extern
SCIP_Real* SCIPcolGetVals(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets node number of the last node in current branch and bound run, where strong branching was used on the
 *  given column, or -1 if strong branching was never applied to the column in current run
 */
extern
SCIP_Longint SCIPcolGetStrongbranchNode(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets number of times, strong branching was applied in current run on the given column */
extern
int SCIPcolGetNStrongbranchs(
   SCIP_COL*             col                 /**< LP column */
   );

/** gets opposite bound type of given bound type */
extern
SCIP_BOUNDTYPE SCIPboundtypeOpposite(
   SCIP_BOUNDTYPE        boundtype           /**< type of bound (lower or upper) */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPcolGetObj(col)              (col)->obj
#define SCIPcolGetLb(col)               (col)->lb
#define SCIPcolGetUb(col)               (col)->ub
#define SCIPcolGetBestBound(col)        ((col)->obj >= 0.0 ? (col)->lb : (col)->ub)
#define SCIPcolGetPrimsol(col)          ((col)->lppos >= 0 ? (col)->primsol : 0.0)
#define SCIPcolGetMinPrimsol(col)       ((col)->minprimsol)
#define SCIPcolGetMaxPrimsol(col)       ((col)->maxprimsol)
#define SCIPcolGetBasisStatus(col)      ((col)->basisstatus)
#define SCIPcolGetVar(col)              (col)->var
#define SCIPcolGetIndex(col)            (col)->index
#define SCIPcolIsIntegral(col)          (col)->integral
#define SCIPcolIsRemovable(col)         (col)->removable
#define SCIPcolGetLPPos(col)            (col)->lppos
#define SCIPcolGetLPDepth(col)          (col)->lpdepth
#define SCIPcolIsInLP(col)              ((col)->lppos >= 0)
#define SCIPcolGetNNonz(col)            (col)->len
#define SCIPcolGetNLPNonz(col)          (col)->nlprows
#define SCIPcolGetRows(col)             (col)->rows
#define SCIPcolGetVals(col)             (col)->vals
#define SCIPcolGetStrongbranchNode(col) (col)->sbnode
#define SCIPcolGetNStrongbranchs(col)   (col)->nsbcalls
#define SCIPboundtypeOpposite(boundtype) \
   ((boundtype) == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER)

#endif




/*
 * Row methods
 */

/** comparison method for sorting rows by non-decreasing index */
extern
SCIP_DECL_SORTPTRCOMP(SCIProwComp);

/** locks an unmodifiable row, which forbids further changes; has no effect on modifiable rows */
extern
void SCIProwLock(
   SCIP_ROW*             row                 /**< LP row */
   );

/** unlocks a lock of an unmodifiable row; a row with no sealed lock may be modified; has no effect on modifiable rows */
extern
void SCIProwUnlock(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the scalar product of the coefficient vectors of the two given rows */
extern
SCIP_Real SCIProwGetScalarProduct(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2                /**< second LP row */
   );

/** returns the degree of parallelism between the hyperplanes defined by the two row vectors v, w:
 *  p = |v*w|/(|v|*|w|);
 *  the hyperplanes are parallel, iff p = 1, they are orthogonal, iff p = 0
 */
extern
SCIP_Real SCIProwGetParallelism(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   );

/** returns the degree of orthogonality between the hyperplanes defined by the two row vectors v, w:
 *  o = 1 - |v*w|/(|v|*|w|);
 *  the hyperplanes are orthogonal, iff p = 1, they are parallel, iff p = 0
 */
extern
SCIP_Real SCIProwGetOrthogonality(
   SCIP_ROW*             row1,               /**< first LP row */
   SCIP_ROW*             row2,               /**< second LP row */
   char                  orthofunc           /**< function used for calc. scalar prod. ('e'uclidean, 'd'iscrete) */
   );

/** output row to file stream */
extern
void SCIProwPrint(
   SCIP_ROW*             row,                /**< LP row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** sorts row entries such that LP columns precede non-LP columns and inside both parts lower column indices precede
 *  higher ones
 */
extern
void SCIProwSort(
   SCIP_ROW*             row                 /**< row to be sorted */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** get number of nonzero entries in row vector */
extern
int SCIProwGetNNonz(
   SCIP_ROW*             row                 /**< LP row */
   );

/** get number of nonzero entries in row vector, that correspond to columns currently in the SCIP_LP;
 *  Warning! This method is only applicable on rows, that are completely linked to their columns (e.g. a row
 *  that is in the current LP and the LP was solved, or a row that was in a solved LP and didn't change afterwards
 */
extern
int SCIProwGetNLPNonz(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets array with columns of nonzero entries */
extern
SCIP_COL** SCIProwGetCols(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets array with coefficients of nonzero entries */
extern
SCIP_Real* SCIProwGetVals(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets constant shift of row */
extern
SCIP_Real SCIProwGetConstant(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets Euclidean norm of row vector */
extern
SCIP_Real SCIProwGetNorm(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets sum norm of row vector (sum of absolute values of coefficients) */
extern
SCIP_Real SCIProwGetSumNorm(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the left hand side of the row */
extern
SCIP_Real SCIProwGetLhs(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the right hand side of the row */
extern
SCIP_Real SCIProwGetRhs(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the dual LP solution of a row */
extern
SCIP_Real SCIProwGetDualsol(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the dual Farkas coefficient of a row in an infeasible LP */
extern
SCIP_Real SCIProwGetDualfarkas(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets the basis status of a row in the LP solution; only valid for LPs with status SCIP_LPSOLSTAT_OPTIMAL
 *  and with SCIPisLPSolBasic(scip) == TRUE; returns SCIP_BASESTAT_BASIC for rows not in the current SCIP_LP
 */
extern
SCIP_BASESTAT SCIProwGetBasisStatus(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns the name of the row */
extern
const char* SCIProwGetName(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets unique index of row */
extern
int SCIProwGetIndex(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets age of row */
extern
int SCIProwGetAge(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff the activity of the row (without the row's constant) is always integral in a feasible solution */
extern
SCIP_Bool SCIProwIsIntegral(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is only valid locally */
extern
SCIP_Bool SCIProwIsLocal(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is modifiable during node processing (subject to column generation) */
extern
SCIP_Bool SCIProwIsModifiable(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is removable from the LP (due to aging or cleanup) */
extern
SCIP_Bool SCIProwIsRemovable(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is member of the global cut pool */
extern
SCIP_Bool SCIProwIsInGlobalCutpool(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets position of row in current LP, or -1 if it is not in LP */
extern
int SCIProwGetLPPos(
   SCIP_ROW*             row                 /**< LP row */
   );

/** gets depth in the tree where the row entered the LP, or -1 if it is not in LP */
extern
int SCIProwGetLPDepth(
   SCIP_ROW*             row                 /**< LP row */
   );

/** returns TRUE iff row is member of current LP */
extern
SCIP_Bool SCIProwIsInLP(
   SCIP_ROW*             row                 /**< LP row */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIProwGetNNonz(row)            (row)->len
#define SCIProwGetNLPNonz(row)          (row)->nlpcols
#define SCIProwGetCols(row)             (row)->cols
#define SCIProwGetVals(row)             (row)->vals
#define SCIProwGetConstant(row)         (row)->constant
#define SCIProwGetNorm(row)             sqrt((row)->sqrnorm)
#define SCIProwGetSumNorm(row)          (row)->sumnorm
#define SCIProwGetLhs(row)              (row)->lhs
#define SCIProwGetRhs(row)              (row)->rhs
#define SCIProwGetDualsol(row)          ((row)->lppos >= 0 ? (row)->dualsol : 0.0)
#define SCIProwGetDualfarkas(row)       ((row)->lppos >= 0 ? (row)->dualfarkas : 0.0)
#define SCIProwGetBasisStatus(row)      (row)->basisstatus
#define SCIProwGetName(row)             (row)->name
#define SCIProwGetIndex(row)            (row)->index
#define SCIProwGetAge(row)            (row)->age
#define SCIProwIsIntegral(row)          (row)->integral
#define SCIProwIsLocal(row)             (row)->local
#define SCIProwIsModifiable(row)        (row)->modifiable
#define SCIProwIsRemovable(row)         (row)->removable
#define SCIProwIsInGlobalCutpool(row)   (row)->inglobalcutpool
#define SCIProwGetLPPos(row)            (row)->lppos
#define SCIProwGetLPDepth(row)          (row)->lpdepth
#define SCIProwIsInLP(row)              ((row)->lppos >= 0)

#endif

#ifdef __cplusplus
}
#endif

#endif
