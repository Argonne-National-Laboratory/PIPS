/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DESYMINDEFSOLVER_H
#define DESYMINDEFSOLVER_H

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "SparseSymMatrix.h"
#include "DenseStorageHandle.h"

/** A linear solver for dense, symmetric indefinite systems 
 * @ingroup DenseLinearAlgebra 
 * @ingroup LinearSolvers
 */
class DeSymIndefSolver : public DoubleLinearSolver {
public:
  DenseStorageHandle mStorage;
protected:
  double* work; int lwork;
  int *ipiv;
  SparseSymMatrix *sparseMat;
public:
  DeSymIndefSolver( DenseSymMatrix * storage );
  DeSymIndefSolver( SparseSymMatrix * storage );
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();
  virtual void solve ( OoqpVector& vec );
  virtual void solve ( GenMatrix& vec );
  virtual ~DeSymIndefSolver();
};

#endif
