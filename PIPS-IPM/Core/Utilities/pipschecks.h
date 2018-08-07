/*
 * pipschecks.h
 *
 *  Created on: 23.02.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_
#define PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_

#include "SimpleVector.h"

bool subMatrixIsOrdered(const int* rowptr, const int* colidx,
      int rowstart, int rowend);

void computeFortranCSRMatResidualNorms(const int* rowptr, const int* colidx, const double* vals, /*const*/ SimpleVector& rhs,
      /*const*/ SimpleVector& x, double& res_norm2, double& res_nrmInf, double& sol_inf, double& mat_max);


#endif /* PIPS_IPM_CORE_UTILITIES_PIPSCHECKS_H_ */
