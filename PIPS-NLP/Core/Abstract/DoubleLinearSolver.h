/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef DOUBLELINEARSOLVER_H
#define DOUBLELINEARSOLVER_H

#include "OoqpVectorHandle.h"
#include "DoubleMatrix.h"


/**
 * @defgroup LinearSolvers
 */

/**
 * Implements the main solver for linear systems that arise in
 * primal-dual interior-point methods for QP.
 * @ingroup LinearSolvers 
 * @ingroup AbstractLinearAlgebra
 */
class DoubleLinearSolver {
public:
  DoubleLinearSolver(){negEigVal=-1; KryIter=0;};
  
  /** called if the diagonal elements of the matrix have
   *  changed. Triggers a refactorization of the matrix, if necessary.
   *
   *  @param idiag index of the first diagonal element that
   *               changed
   *  @param extent the number of diagonal element that changed.  */
  virtual void diagonalChanged( int idiag, int extent ) = 0;

  /** called if some elements of the matrix have changed.  Triggers a
   *  refactorization of the matrix, if necessary.  */
  virtual int matrixChanged() = 0;

  /** solves a linear system.
   *
   * @param x on entry the right hand side of the system to be solved.
   *           On exit, the solution.  */
  virtual void solve ( OoqpVector& x ) = 0;
  // solve trans sys, i.e: A'x=b, default is symmetric 
  virtual void solveTrans ( OoqpVector& x ) {this->solve(x);};

	// solve with multiple RHS
	virtual void solve ( GenMatrix& rhs ) { assert(0 && "Not implemented"); }

  virtual void Lsolve  ( OoqpVector& x ) {}
  virtual void Dsolve  ( OoqpVector& x ) { solve(x);}
  virtual void Ltsolve ( OoqpVector& x ) {}

  /** Destructor  */
  virtual ~DoubleLinearSolver() {};

  int negEigVal;
  int KryIter;  
};

/**
 * The abstract  interface to a mat-vec operation required by 
 * the iterative solvers.
 * @ingroup LinearSolvers 
 * @ingroup AbstractLinearAlgebra
 */

class MatTimesVec {
 public:
  /** y = beta * y + alpha * A * x */
  virtual void doIt(double beta, OoqpVector& y, double alpha, OoqpVector& x) = 0;
};

/**
 * The abstract interface for a linear iterative solver for linear systems 
 * that arise in primal-dual interior-point methods for QP.
 * @ingroup LinearSolvers 
 * @ingroup AbstractLinearAlgebra
 */

class DoubleIterativeLinearSolver : public DoubleLinearSolver {
 public:
  DoubleIterativeLinearSolver( MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2=NULL );

  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();

  virtual void solve( OoqpVector& rhs ) = 0;

  virtual ~DoubleIterativeLinearSolver();
  
 protected:
  DoubleIterativeLinearSolver();
  /** MatVec operation involving system matrix*/
  MatTimesVec *A;

  /** MatVec ops for left and right preconditioner */
  MatTimesVec *ML, *MR;

  /** Actual mat-vec operations */
  void applyA (double beta, OoqpVector& res, double alpha, OoqpVector& x);
  void applyM1(double beta, OoqpVector& res, double alpha, OoqpVector& x);
  void applyM2(double beta, OoqpVector& res, double alpha, OoqpVector& x);
};

/** 
 * An implementation of the abstract class @MatTimesVec that performs a mat-vec
 * with both the matrix and vector being on the same processor.
 * 
 * It can use OOQP matrix and the implementation is based on SimpleVector class.
 */
class StoredMatTimesVec : public MatTimesVec {
 public:
  StoredMatTimesVec(DoubleMatrix* mat);
  virtual ~StoredMatTimesVec() {};

  void doIt(double beta, OoqpVector& y, double alpha, OoqpVector& x);
 protected:
  DoubleMatrix* mMat;
};
/** 
 * An implementation of the abstract class @MatTimesVec that performs a 
 * mat transpose-vec with both the matrix and vector being on the same processor.
 * 
 * It can use OOQP matrix and the implementation is based on SimpleVector class.
 */
class StoredMatTransTimesVec : public MatTimesVec {
 public:
  StoredMatTransTimesVec(DoubleMatrix* mat);
  virtual ~StoredMatTransTimesVec() {};

  void doIt(double beta, OoqpVector& y, double alpha, OoqpVector& x);
 protected:
  DoubleMatrix* mMat;
};
#endif




