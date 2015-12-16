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

#include <map>


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
protected:
  PardisoSchurSolver() {};
  void SetUpPardisoSchurSolver(SparseSymMatrix * sgm);
  
 public:
  virtual void firstCall(); //first factorization call
  void firstSolveCall(SparseGenMatrix& R, 
		      SparseGenMatrix& A,
		      SparseGenMatrix& C); //first solve call
  
  /** sets mStorage to refer to the argument sgm */
  PardisoSchurSolver( SparseSymMatrix * sgm ); 
  PardisoSchurSolver( SparseSymMatrix * sgm, const int numOfNegEigVal_in ); 
  
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve( OoqpVector& rhs );
  virtual void solve( GenMatrix& rhs);
 
  /** Functions specific to the Schur approach. The last argument is the Schur first
   * stage matrix that will be updated.
   * 1. schur_solver( rhs, SC)
   *  - this is the generic function
   *
   * 2. schur_solve(R,A,B, SC) 
   *  - avoids forming the matrix rhs [R' A' B']'
   *  - assumes rhs block matrices does not change
   */
  virtual void schur_solve(/*const*/ SparseGenMatrix& R, 
			   /*const*/ SparseGenMatrix& A,
			   /*const*/ SparseGenMatrix& C,
			   DenseSymMatrix& SC);
  
 protected:
  SparseSymMatrix* Msys; // this is the (1,1) block in the augmented system

//  SparseSymMatrix* MAug_LowTri; // this is the augmented system in lower triangular form
  
  bool first;
  bool firstSolve;
  void  *pt[64]; 
  int iparm[64];
  double dparm[64];
  int num_threads;

  /** dimension of the PARDISO augmented system */
  int n; 
  /** dimension of the Schur complement (# of rhs) */
  int nSC; 
  /** number of nonzeros in the PARDISO augmented matrix */
  int      nnz;
  /** storage for the upper triangular (in row-major format) */
  int     *rowptrAug, *colidxAug;
  double  *eltsAug;

  // index in Lower triangular form
  int	  *index_LowTri;
  /** mapping from from the diagonals of the PIPS linear systems to 
      the diagonal elements of the (1,1) block  in the augmented system */
  map<int,int> diagMap;


  /** mapping from the element of the PIPS linear systems to 
      the elements of the (1,1) block  in the augmented system, hence we can update Jac */
//  map<int,int> elementMap_Rt;  
//  map<int,int> elementMap_At; 
//  map<int,int> elementMap_Ct;   
//  map<int,int> elementMap_BD; 

  int* elementMap_R;
  int* elementMap_A; 
  int* elementMap_C;

  int numOfNegEigVal;

  //temporary vector of size n
  double* nvec;
  
  virtual ~PardisoSchurSolver();
};

class PardisoSchur32Solver : public PardisoSchurSolver
{
 public:
    PardisoSchur32Solver( SparseSymMatrix * sgm ); 
	PardisoSchur32Solver( SparseSymMatrix * sgm, const int numOfNegEigVal_in ); 
	
 private: 
    PardisoSchur32Solver () {};
 public:
    virtual void firstCall(); //first factorization call
    virtual void solve (OoqpVector& rhs );
    
};

#endif
