/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin, Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
 
/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef PARDISOLINSYS_H
#define PARDISOLINSYS_H

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"
#include "DenseSymMatrix.h"

#include <map>


#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

using namespace std;
/** implements the linear solver class using the Pardiso solver
 */
 
class PardisoSolver : public DoubleLinearSolver {
private:
  PardisoSolver() {};
  void SetUpPardisoSolver(SparseSymMatrix * sgm);
  
 public:
  virtual void firstCall();
  
  /** sets mStorage to refer to the argument sgm */
//  PardisoSolver( SparseSymMatrix * sgm ); 
  PardisoSolver( SparseSymMatrix * sgm, const int numOfNegEigVal_in =-1); 
  PardisoSolver( DenseSymMatrix* m);
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve( OoqpVector& rhs );
  virtual void solve( GenMatrix& rhs);
  virtual void solve( GenMatrix& rhs, int *colSparsity);  
 // virtual void Lsolve( OoqpVector& x );
 // virtual void Dsolve( OoqpVector& x );
 // virtual void Ltsolve( OoqpVector& x );

 private:
  //Helper functions for when the input is a dense matrix
  int getNumberOfNonZeros(DenseSymMatrix& m);
  
 private:
  SparseSymMatrix* Msys;
  DenseSymMatrix* Mdsys;
  bool first;
  void  *pt[64]; 
  int iparm[64];
  int num_threads;
  double dparm[64];
  int n;
  int numOfNegEigVal;
  
  /** storage for the upper triangular (in row-major format) */
  int     *krowM,    *jcolM;
  double  *M;
  
  /** number of nonzeros in the matrix */
  int      nnz;

  /** mapping from from the diagonals of the PIPS linear systems to 
      the diagonal elements of the (1,1) block  in the augmented system */
  std::map<int,int> diagMap;

  /** mapping from the element of the PIPS linear systems to 
    the elements of the (1,1) block  in the augmented system, hence we can update Jac */
  map<int,int> elementMap_A;  
  map<int,int> elementMap_BD; 
  map<int,int> elementMap_C; 

  /** temporary storage for the factorization process */
  double* nvec; //temporary vec
  double* sol; //solution
  int sz_sol; //allocated size
  virtual ~PardisoSolver();
};

#endif
