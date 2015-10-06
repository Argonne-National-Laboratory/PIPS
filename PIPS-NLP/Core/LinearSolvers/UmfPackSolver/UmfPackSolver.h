/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef UMFPACKLINSYS_H
#define UMFPACKLINSYS_H

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseGenMatrixHandle.h"

#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"

#include "umfpack.h"



/** implements the linear solver class using the HSL UmfPack solver
 *
 * @ingroup LinearSolvers 
 */
class UmfPackSolver : public DoubleLinearSolver {

public:
  bool firstCallFlag;

  double umf_Control [UMFPACK_CONTROL];
  double umf_Info [UMFPACK_INFO];
  void *Symbolic, *Numeric ;

  int *irowM,*kcolbegM, *eleMap;

  /** storage for the original matrix */
  double  *eleM;


  int matrixSingular;
  

  /** dimension of the whole matrix */
  int      n;
  /** number of nonzeros in the matrix */
  int      nnz;
  /** store as a sparse symmetric matrix */
  SparseStorageHandle mStorage;

//************************** methods *********************************
private:
  UmfPackSolver() {};

public:
  virtual ~UmfPackSolver();
  /** sets mStorage to refer to the argument sgm */
  UmfPackSolver( SparseGenMatrix * sgm );
  UmfPackSolver( SparseSymMatrix * sgm );
  UmfPackSolver( SparseSymMatrix * sgm, const int numOfNegEigVal_in ); 

  
  virtual void diagonalChanged( int idiag, int extent ){matrixChanged();};
  virtual int matrixChanged();
  virtual void solve( OoqpVector& rhs );
  virtual void solve( GenMatrix& rhs);

  virtual void solveTrans( OoqpVector& rhs_in );
  virtual void solveTrans( GenMatrix& rhs_in);

  virtual void Lsolve  ( OoqpVector& x );
  virtual void Dsolve  ( OoqpVector& x ){this->solve(x);};
  virtual void Ltsolve ( OoqpVector& x );

  void freeSymFactInfo(){umfpack_di_free_symbolic (&Symbolic) ;Symbolic=NULL;};
  void freeNumFactInfo(){umfpack_di_free_numeric (&Numeric) ; Numeric=NULL;};
private:
  void solve(int solveType, OoqpVector& rhs);




  
protected:
  virtual void firstCall();	
  virtual void firstCall_Sym();	
  int isSymm;
};

#endif
