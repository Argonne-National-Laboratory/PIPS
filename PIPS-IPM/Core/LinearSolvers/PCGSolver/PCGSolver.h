#ifndef ProjCGLINSYS_H
#define ProjCGLINSYS_H

#include "DoubleLinearSolver.h"
#include "DoubleMatrix.h"
#include "SimpleVector.h"

class MatTimesVec;

/** a linear solver class implementing the
 * projected CG for KKT (saddle) linear systems of the form
 *  [ H  A'] [x]  =  [b1], H is n x n
 *  [ A  0 ] [y]  =  [b2], A is m x n
 *
 * The algorithm is a modification of the projected CG in the sense
 * that it operates on the full x-space and does not need an actual
 * representation of the Ker(A).  For more details see 
 * N.I.M Gould, M.R. Hribar, and J.Nocedal, On the solution of
 * equality constrained quadratic programming problems arising in
 * optimization, SIAM. J Sci.Comput., Vol 23, No 4, pp.1376-1395.
 *
 * A preconditioner G of H should be available and the above saddle
 * point liner system is preconditioned by P given by
 *        [ G  A']
 *        [ A  0 ]
 *
 * @ingroup LinearSolvers 
 */
class PCGSolver : public  DoubleIterativeLinearSolver {
public:
  /** initialization constructor */
  PCGSolver( MatTimesVec* H, MatTimesVec* P, MatTimesVec* At, int n, int m);

  virtual ~PCGSolver();

  /** version of the main solve routine that takes argument as an
   * OoqpVector
   *
   * @param drhs on input contains the right-hand side; on output
   * contains the solution
   */
  void solve( OoqpVector& rhs );

  /** the standard projected CG */
  void solvefull(OoqpVector& rhs);

  /** an economical version of projected CG */
  void solveecon(OoqpVector& rhs);

  /** "outer" iterative refinement, used if the algorithm was able to reduce the
   *  residual but not to the desired level. */
  void refine(OoqpVector& rhs);
protected:
  PCGSolver() {};

  double tol; double iter; int maxit; int flag;

  double *tmpVec1,*tmpVec2,*tmpVec3,*tmpVec4,*tmpVec5,*tmpVec6;

  MatTimesVec* At;

  int n,m;
};


#endif



