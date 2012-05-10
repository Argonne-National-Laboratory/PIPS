/* PIPS                                                               *
 * Authors: Miles Lubin                                               *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DESYMINDEFSOLVER2_H
#define DESYMINDEFSOLVER2_H

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "DenseStorageHandle.h"

/** Specialized LDL^T solver for saddle point systems 
 * @ingroup DenseLinearAlgebra 
 * @ingroup LinearSolvers
 */
class DeSymIndefSolver2 : public DoubleLinearSolver {
public:
  DenseStorageHandle mStorage;
protected:
  int nx, ny, n;
public:
  DeSymIndefSolver2( DenseSymMatrix * storage, int nx );
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();
  virtual void solve ( OoqpVector& vec );
  virtual ~DeSymIndefSolver2();
};

#endif
