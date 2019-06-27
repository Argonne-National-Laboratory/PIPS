/*
 * MumpsSolverRoot.h
 *
 *  Created on: 25.06.2019
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_

#include "dmumps_c.h"

#include "mpi.h"
#include "MumpsSolverBase.h"
#include "SparseSymMatrix.h"
#include "OoqpVector.h"
#include "pipsport.h"


/** implements linear solver class for root nodes that uses the MUMPS solver
 */

class MumpsSolverRoot : public MumpsSolverBase {

 public:
  MumpsSolverRoot( SparseSymMatrix * sgm );
  MumpsSolverRoot( MPI_Comm mpiComm, SparseSymMatrix * sgm );

  ~MumpsSolverRoot();

  void matrixRebuild( DoubleMatrix& matrixNew ) override;
  void matrixChanged() override;
  void solve( OoqpVector& rhs ) override;

 private:
  void factorize();
};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_ */
