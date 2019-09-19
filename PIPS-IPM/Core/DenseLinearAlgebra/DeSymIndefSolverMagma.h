
#ifndef DESYMINDEFSOLVERMAGMA_H
#define DESYMINDEFSOLVERMAGMA_H

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "SparseSymMatrix.h"
#include "DenseStorageHandle.h"
#include "pipsport.h"

/** A linear solver for dense, symmetric indefinite systems
 * @ingroup DenseLinearAlgebra
 * @ingroup LinearSolvers
 */
class DeSymIndefSolverMagma : public DoubleLinearSolver {
public:
  DenseStorageHandle mStorage;
protected:
  int *ipiv;
  SparseSymMatrix *sparseMat;
#ifdef GPUCODE
  double* mFact_gpu;
  double* mRhs_gpu;
#endif
public:
  DeSymIndefSolverMagma( DenseSymMatrix * storage );
  DeSymIndefSolverMagma( SparseSymMatrix * storage );
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();
  using DoubleLinearSolver::solve;
  void solve ( OoqpVector& vec ) override;
  void solve ( GenMatrix& vec ) override;
  virtual ~DeSymIndefSolverMagma();
};

#endif
