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
#include "pipsport.h"
#include <map>

using namespace std;

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

 constexpr static int symbFactorIntervalDefault = 3;
 constexpr static int pivotPerturbationExpDefault = 8;

protected:
  PardisoSchurSolver() {};
  
 public:
  virtual void firstCall(); //first factorization call
  void firstSolveCall(SparseGenMatrix& R, 
		      SparseGenMatrix& A,
		      SparseGenMatrix& C,
		      SparseGenMatrix& F,
		      SparseGenMatrix& G,
		      int nSC0); //first solve call
  
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
			   /*const*/ SparseGenMatrix& F,
			   /*const*/ SparseGenMatrix& G,
			   DenseSymMatrix& SC);
  
  virtual void schur_solve_sparse(/*const*/ SparseGenMatrix& R,
            /*const*/ SparseGenMatrix& A,
            /*const*/ SparseGenMatrix& C,
            /*const*/ SparseGenMatrix& F,
            /*const*/ SparseGenMatrix& G,
            SparseSymMatrix& SC);

 protected:
  SparseSymMatrix* Msys; // this is the (1,1) block in the augmented system
  bool first;
  bool firstSolve;
  void  *pt[64];
  int iparm[64];
  bool useSparseRhs;
  int symbFactorInterval;


#ifndef WITH_MKL_PARDISO
  int num_threads;
  double dparm[64];
#endif

  int* shrinked2orgSC;

  int pivotPerturbationExp; // 10^-exp

  /* pardiso params */
  int maxfct, mnum, phase, msglvl, solver, mtype, nrhs;

  /** dimension of the PARDISO augmented system */
  int n; 
  /** dimension of the Schur complement (# of rhs) */
  int nSC; 
  /** number of nonzeros in the PARDISO augmented matrix */
  int nnz;
  /** storage for the upper triangular (in row-major format) */
  int *rowptrAug, *colidxAug;
  double *eltsAug;
  /** mapping from from the diagonals of the PIPS linear systems to 
      the diagonal elements of the (1,1) block  in the augmented system */
  map<int,int> diagMap;

  //temporary vector of size n
  double* nvec;
  double* nvec2;
  int nvec_size; // to be save
  

  void setIparm(int* iparm);
  bool iparmUnchanged();

  virtual void computeSC(
            int nSCO,
            /*const*/ SparseGenMatrix& R,
            /*const*/ SparseGenMatrix& A,
            /*const*/ SparseGenMatrix& C,
            /*const*/ SparseGenMatrix& F,
            /*const*/ SparseGenMatrix& G,
            int*& rowptrSC,
            int*& colidxSC,
            double*& eltsSC
  );

  virtual ~PardisoSchurSolver();
};

class PardisoSchur32Solver : public PardisoSchurSolver
{
 public:
    PardisoSchur32Solver( SparseSymMatrix * sgm ); 
 private: 
    PardisoSchur32Solver () {};
 public:
    virtual void firstCall(); //first factorization call
    virtual void solve (OoqpVector& rhs );
    
};

#endif
