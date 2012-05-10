/* PIPS                                                               *
 * Author: Miles Lubin                                                */

#ifndef WSMPLINSYS_H
#define WSMPLINSYS_H

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"

/** implements the linear solver class using the HSL MA57 solver
 *
 * @ingroup LinearSolvers 
 */
class WSMPSolver : public DoubleLinearSolver {
private:
  WSMPSolver() {};
protected:
	int    iparm[64];
	double dparm[64]; 

	bool first;
	int nthreads, instance;
  static int instances;
  /** the Threshold Pivoting parameter, stored as U in the ma27dd
   *  common block. Takes values in the range [0,1]. Larger values
   *  enforce greater stability in the factorization as they insist on
   *  larger pivots. Smaller values preserve sparsity at the cost of
   *  using smaller pivots.  */
  double   kThresholdPivoting;

  /** the Threshold Pivoting parameter may need to be increased during
   * the algorithm if poor precision is obtained from the linear
   * solves.  kThresholdPivoting indicates the largest value we are
   * willing to tolerate.  */
  double   kThresholdPivotingMax;

  /** the factor in the range (1,inf) by which kThresholdPivoting is
   * increased when it is found to be inadequate.  */
  double   kThresholdPivotingFactor;

  /** the "Treat As Zero" parameter, stored as pivtol in the common
   * block ma27td. The factorization will not accept a pivot whose
   * absolute value is less than this parameter as a 1x1 pivot or as
   * the off-diagonal in a 2x2 pivot.  */
  double   kTreatAsZero;

  /** precision we demand from the linear system solver. If it isn't
   * attained on the first solve, we use iterative refinement and
   * possibly refactorization with a higher value of
   * kThresholdPivoting. */
  double  kPrecision;

  /** storage for the transpose */
  int     *krowMt,    *jcolMt;
  double  *Mt;

  /** dimension of the whole matrix */
  int      n;

  /** number of nonzeros in the matrix */
  int      nnz;

  /** temporary storage for the factorization process */
  int     *perm, *invp, *diagmap;


  /** used to indicate when we need a fresh factorization (when
   * iterative refinement has failed to improve the precision of the
   * computed solution satisfactorily */
  int     freshFactor;

  /** store as a sparse symmetric matrix */
  SparseStorageHandle mStorage;

  /** called the very first time a matrix is factored. Allocates space
   * for the factorization and performs ordering */
  virtual void firstCall();
public:
  /** sets mStorage to refer to the argument sgm */
  WSMPSolver( SparseSymMatrix * sgm );
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();
  virtual void solve( OoqpVector& rhs );
	virtual void solve( GenMatrix& rhs);

  virtual void Lsolve  ( OoqpVector& x );
  virtual void Dsolve  ( OoqpVector& x );
  virtual void Ltsolve ( OoqpVector& x );
  //virtual void Refine  ( OoqpVector& x );
 private:
  void solve(int solveType, OoqpVector& rhs);


 public:
  

  /** destructor */
  virtual ~WSMPSolver();
};

#endif
