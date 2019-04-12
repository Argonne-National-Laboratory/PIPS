/*
 * pipschecks.h
 *
 *  Created on: 23.02.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_
#define PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_

#include "StochVector.h"
#include "StochGenMatrix.h"
#include "SimpleVector.h"
#include <vector>

class sData;

// is the permuation vector valid?
bool permutationIsValid(const std::vector<unsigned int>& perm);

// are the columns of the given sub-matrix ordered?
bool subMatrixIsOrdered(const int* rowptr, const int* colidx,
      int rowstart, int rowend);

// compute residual norms for Ax=rhs with A in CSR with 1-indexing (Fortran)
void computeFortranCSRMatResidualNorms(const int* rowptr, const int* colidx, const double* vals, /*const*/ SimpleVector& rhs,
      /*const*/ SimpleVector& x, double& res_norm2, double& res_nrmInf, double& sol_inf, double& mat_max);

// is root node data of sData object same on all procs?
bool rootNodeInSyncSData(const sData& s_data);

// is root node data of StochVector same on all procs?
bool rootNodeInSyncStochVector(const StochVector& stoch_vec);

// is root node data of StochMatrix same on all procs?
// not checking dynamic storage !
bool rootNodeInSyncStochGenMatrix(const StochGenMatrix& stoch_mat);

#endif /* PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_ */
