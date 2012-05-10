#ifndef BICGSTABLINSYS_H
#define BICGSTABLINSYS_H

#include "DoubleLinearSolver.h"
#include "DoubleMatrix.h"
#include "SimpleVector.h"

class MatTimesVec;

/** implements the linear solver class using an implementation of the stabilized BiCG
 *
 * @ingroup LinearSolvers 
 */
class BiCGStabSolver : public  DoubleIterativeLinearSolver {
public:
  /** build the class for the linear system 
   *  for Ax=b preconditioned with M1 (or M1 and M2)
   */
  BiCGStabSolver( MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2=NULL );

  virtual ~BiCGStabSolver() {};

  /** version of the main solve routine that takes argument as an
   * OoqpVector
   *
   * @param drhs on input contains the right-hand side; on output
   * contains the solution
   */
  void solve( OoqpVector& rhs );

protected:
  BiCGStabSolver() {};

  double tol; double iter; int maxit; int flag;
};


#endif



