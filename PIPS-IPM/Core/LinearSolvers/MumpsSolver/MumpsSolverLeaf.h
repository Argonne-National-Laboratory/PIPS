/*
 * MumpsSolverLeaf.h
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MumpsSolverLeaf_MumpsSolverLeafLEAF_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MumpsSolverLeaf_MumpsSolverLeafLEAF_H_

#include "dmumps_c.h"
#include "mpi.h"
#include "MumpsSolverBase.h"
#include "SparseSymMatrix.h"
#include "pipsport.h"


/** implements linear solver class for leaf nodes that uses the MUMPS solver
 */

class MumpsSolverLeaf : public MumpsSolverBase {

 public:
  MumpsSolverLeaf( SparseSymMatrix * sgm );

  ~MumpsSolverLeaf();

  void matrixChanged() override;

  // rhs need to be in CSC Fortran format
  void solve(/* const */ GenMatrix& rhs_f,  int startRow, int range, double* sol);
};



#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MumpsSolverLeaf_MumpsSolverLeafLEAF_H_ */
