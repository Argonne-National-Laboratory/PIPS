/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP
 *
 * C. Petra Added multiple Lsolve with multiple right hand sides      */

#ifndef DESYMPSDSOLVER_H
#define DESYMPSDSOLVER_H

#include "DoubleLinearSolver.h"
#include "DenseStorageHandle.h"
#include "DenseSymMatrixHandle.h"
#include "OoqpVectorHandle.h"
#include "DenseGenMatrix.h"
#include "pipsport.h"

/** A linear solver for dense, symmetric positive-definite systems.
 *  @ingroup DenseLinearAlgebra
 *  @ingroup LinearSolvers
 */
class DeSymPSDSolver : public DoubleLinearSolver {
protected:
  DenseStorageHandle mStorage;
public:
  DeSymPSDSolver( DenseSymMatrix * dsm );
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();

  using DoubleLinearSolver::solve;
  void solve ( OoqpVector& x ) override;

  //specialized method that uses BLAS-3 function DTRSM for the triagular solve.
  using DoubleLinearSolver::Lsolve;
  void Lsolve ( GenMatrix& mat ) override;

  virtual ~DeSymPSDSolver();
};

#endif
