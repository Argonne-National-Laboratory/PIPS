/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#ifndef PARDISO_SCHUR_SOLVER
#define PARDISO_SCHUR_SOLVER

#include "DoubleLinearSolver.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "DenseSymMatrix.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"



#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif


/** implements the linear solver class using the Pardiso SC solver
 */
 
class PardisoSchurSolver : public DoubleLinearSolver {
private:
  PardisoSchurSolver() {};
  
 public:
  virtual void firstCall(); //first factorization call
  void firstSolveCall(); //first solve call
  
  /** sets mStorage to refer to the argument sgm */
  PardisoSchurSolver( SparseSymMatrix * sgm ); 
  
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();
  virtual void solve( OoqpVector& rhs );
  virtual void solve( GenMatrix& rhs);
 
  /** Functions specific to the Schur approach. The last argument is the Schur first
   * stage matrix that will be updated.
   * 1. schur_solver( rhs, SC)
   *  - this is the generic function
   *
   * 2. schur_solve(R,A,B, SC) 
   *  - avoids forming the matrix rhs [R' A' B']'
   *  - assumes rhs does not change
   */
  virtual void schur_solve(/*const*/ SparseGenMatrix& R, 
			   /*const*/ SparseGenMatrix& A,
			   /*const*/ SparseGenMatrix& C,
			   DenseSymMatrix& SC);
  
 private:
  SparseSymMatrix* Msys; // this is the (1,1) block in the augmented system
  bool first;
  bool firstSolve;
  void  *pt[64]; 
  int mtype;
  int solver;
  int iparm[64];
  int num_threads;
  
  double b[8], x[8];
  double dparm[64];
  int error;
  int nrhs;  //  Number of right-hand sides 
  int maxfct;
  int mnum;
  int phase;
  int n;

  /** storage for the upper triangular (in row-major format) */
  int     *krowM,    *jcolM;
  double  *M;
  
  
  /** number of nonzeros in the matrix */
  int      nnz;
  
  /** temporary storage for the factorization process */
  //int     *perm, *invp, *diagmap;
  double* nvec; //temporary vec
  
  virtual ~PardisoSchurSolver();
};

#endif
