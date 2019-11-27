#ifndef CGLINSYS_H
#define CGLINSYS_H

#include "DoubleLinearSolver.h"
#include "DoubleMatrix.h"
#include "SimpleVector.h"
#include "pipsport.h"

class MatTimesVec;

/** implements the linear solver class using an implementation of the stabilized BiCG
 *
 * @ingroup LinearSolvers 
 */
class CGSolver : public  DoubleIterativeLinearSolver {
public:
  /** build the class for the linear system 
   *  for Ax=b preconditioned with M1 (or M1 and M2)
   */
  CGSolver( MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2=nullptr );

  virtual ~CGSolver();

  /** version of the main solve routine that takes argument as an
   * OoqpVector
   *
   * @param drhs on input contains the right-hand side; on output
   * contains the solution
   */
  void solve( OoqpVector& rhs );

protected:
  CGSolver() {};

  double tol; double iter; int maxit; int flag;

  double *tmpVec1,*tmpVec2,*tmpVec3,*tmpVec4,*tmpVec5,*tmpVec6;
  //int firstSolve;
};


#endif



