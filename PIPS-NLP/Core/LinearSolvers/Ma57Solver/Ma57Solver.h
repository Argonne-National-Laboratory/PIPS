/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Rewritten by Nai-Yuan Chiang for NLP*/

#ifndef MA57LINSYS_H
#define MA57LINSYS_H

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"

#ifdef WITH_PETSC
#include "PetscIterativeSolver_Sparse.h"
#endif

#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif


extern "C" {
  void FNAME(ma57id)( double cntl[],  int icntl[] );

  void FNAME(ma57ad)( int * n,        int * ne,       int irn[],     
		int jcn[],      int * lkeep,    int keep[],
		int iwork[],    int icntl[],    int info[],
		double rinfo[] );

  void FNAME(ma57bd)( int * n,        int * ne,       double a[],
		double fact[],  int * lfact,    int ifact[],
		int * lifact,   int * lkeep,    int keep[],
		int ppos[],     int * icntl,    double cntl[],
		int info[],     double rinfo[] );
  void FNAME(ma57cd)( int * job,      int * n,        double fact[],
		int * lfact,    int ifact[],    int * lifact,
		int * nrhs,     double rhs[],   int * lrhs,
		double w[],     int * lw,       int iw1[],
		int icntl[],    int info[]);
  void FNAME(ma57dd)( int * job,      int * n,        int * ne,
		double a[],     int irn[],      int jcn[],
		double fact[],  int * lfact,    int ifact[],
		int * lifact,   double rhs[],   double x[],
		double resid[], double w[],     int iw[],
		int icntl[],    double cntl[],  int info[],
		double rinfo[] );
  void FNAME(ma57ed)( int * n,        int * ic,       int keep[],
		double fact[],  int * lfact,    double * newfac,
		int * lnew,     int  ifact[],   int * lifact,
		int newifc[],   int * linew,    int * info );
}

/** implements the linear solver class using the HSL MA57 solver
 *
 * @ingroup LinearSolvers 
 */
class Ma57Solver : public DoubleLinearSolver {
private:
  Ma57Solver() {};
  void SetUpMa57Solver( SparseSymMatrix * sgm);
protected:
  int     icntl[20];
  int     info[40];
  double  cntl[5];
  double  rinfo[20];

  /** the Threshold Pivoting parameter, stored as U in the ma27dd
   *  common block. Takes values in the range [0,1]. Larger values
   *  enforce greater stability in the factorization as they insist on
   *  larger pivots. Smaller values preserve sparsity at the cost of
   *  using smaller pivots.  */
  double   kThresholdPivoting;
;

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

  /** index array for the factorization */
  int     *irowM,    *jcolM;

  /** storage for the original matrix */
  double  *M;

  /** dimension of the whole matrix */
  int      n;

  /** number of nonzeros in the matrix */
  int      nnz;

  /** temporary storage */
  int     lkeep, *keep;

  /** temporary storage for the factorization process */
  int     lifact, *ifact, lfact;

  /* storage for the factors */
  double *fact;

  /** amounts by which to increase allocated factorization space when
   * inadequate space is detected. ipessimism is for array "iw",
   * rpessimism is for the array "fact". */
  double  ipessimism, rpessimism;


  /** store as a sparse symmetric matrix */
  SparseStorageHandle mStorage;

  /** called the very first time a matrix is factored. Allocates space
   * for the factorization and performs ordering */
  virtual void firstCall();
public:

  Ma57Solver( SparseSymMatrix * sgm, const int numOfNegEigVal_in = -1 ); 
  virtual ~Ma57Solver();
  
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve( OoqpVector& rhs );
  virtual void solve( GenMatrix& rhs);

  //virtual void Lsolve  ( OoqpVector& x );
  //virtual void Dsolve  ( OoqpVector& x );
  //virtual void Ltsolve ( OoqpVector& x );
  //virtual void Refine  ( OoqpVector& x );
 private:
  void solve(int solveType, OoqpVector& rhs);

  int* iworkn, niworkn;  
  int* new_iworkn(int dim);

  double* dworkn; int ndworkn;
  double* new_dworkn(int dim);
 public:
  /** if pivot < this val, treat it as Zero in MA57  */
  void setValTreatAsZero(const double tolTreatAsZero) { cntl[1] = tolTreatAsZero; }

  /** if pivot tolerence in MA57  */
  void setPivotTol(const double tolPivot) {  cntl[0] = tolPivot; }  

#ifdef WITH_PETSC
  friend class PetscIterativeSolver_Sparse;
#endif
  
};

#endif
