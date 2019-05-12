/*
 * MumpsSolver.h
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_

//#include "dmumps_c.h"

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "OoqpVectorHandle.h"


/** implements linear solver class that uses the MUMPS solver
 */


class MumpsSolver : public DoubleLinearSolver {
 private:
  MumpsSolver() {};

 public:
  /** sets mStorage to refer to the argument sgm */
  MumpsSolver( SparseSymMatrix * sgm );
  ~MumpsSolver();

  void diagonalChanged( int idiag, int extent );
  void matrixChanged();
  void solve( OoqpVector& rhs );
  void solve( int nrhss, double* rhss, int* colSparsity );

 private:
#if 0
  static MUMPS_INT getFortranMPIComm(MPI_Comm mpiComm_c)
  {
     return MUMPS_INT(MPI_Comm_c2f(mpiComm_c));
  };

  DMUMPS_STRUC_C* mumps;
#endif
  SparseSymMatrix* Msys;


};



#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_ */
