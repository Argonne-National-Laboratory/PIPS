/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PARDISOLINSYS_H
#define PARDISOLINSYS_H

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"



#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif




/** implements the linear solver class using the HSL Pardiso solver
 *
 * @ingroup LinearSolvers 
 */
 
class PardisoSolver : public DoubleLinearSolver {
private:
  PardisoSolver() {};
  
public:
  virtual void firstCall();

  /** sets mStorage to refer to the argument sgm */
  PardisoSolver( SparseSymMatrix * sgm );
  
  void pardisoInitMurat(SparseSymMatrix * sgm);
  void pardisoMurat();
 
  
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual void matrixChanged();
  virtual void solve( OoqpVector& rhs );
  virtual void solve( GenMatrix& rhs);
  
 // virtual void Lsolve( OoqpVector& x );
 // virtual void Dsolve( OoqpVector& x );
 // virtual void Ltsolve( OoqpVector& x );
  
 private:
 	SparseSymMatrix * sgm_;
    bool first;
	void  *pt[64]; 
	int mtype;
	int solver;
	int iparm[64];
	
	double b[8], x[8];
	double dparm[64];
	int error;
	int nrhs;  //  Number of right-hand sides 
	int maxfct;
	int mnum;
	int phase;
	int n;
	



SparseStorageHandle mStorage;  // Murat
	
	/** storage for the transpose */
  int     *krowMt,    *jcolMt;
  double  *Mt;

  
  /** number of nonzeros in the matrix */
  int      nnz;

  /** temporary storage for the factorization process */
  int     *perm, *invp, *diagmap;

  virtual ~PardisoSolver();
};

#endif
