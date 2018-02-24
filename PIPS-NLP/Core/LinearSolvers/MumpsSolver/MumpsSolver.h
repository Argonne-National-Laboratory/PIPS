/* PIPS file - Written by Cosmin Petra, 2018 */

#ifndef MUMPSSOLVER_H
#define MUMPSSOLVER_H

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "SparseSymMatrix.h"
#include "DenseStorageHandle.h"

#include "dmumps_c.h"

/** A linear solver for symmetric indefinite systems using MUMPS
 *
 */
class MumpsSolver : public DoubleLinearSolver {
public:
  DenseStorageHandle mStorage;
protected:
  double* work; int lwork;
  int *ipiv;
  SparseSymMatrix *sparseMat;
public:
  MumpsSolver( DenseSymMatrix * storage );
  MumpsSolver( SparseSymMatrix * storage );
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve ( OoqpVector& vec );
  virtual void solve ( GenMatrix& vec );
  virtual ~MumpsSolver();
};

#endif
 
