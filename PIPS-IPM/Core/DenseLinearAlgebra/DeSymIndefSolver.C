/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "DeSymIndefSolver.h"
#include "SimpleVector.h"
#include <cassert>

#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"


#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

// declarations for LAPACK functions used to factor/solve:

// dsytrf_() factors a symmetric indefinite matrix A, see LAPACK 
// documentation for more details.
extern "C" void FNAME(dsytrf)(char *uplo, 
			int *n, 
			double A[], 
			int *lda, 
			int ipiv[], 
			double work[],
			int *lwork, 
			int *info);

// dsytrs_() solves the system Ax = b using the factor obtained by dsytrf_().
extern "C" void FNAME(dsytrs)(char *uplo, 
			int *n, 
			int *nrhs, 
			double A[], 
			int *lda, 
			int ipiv[], 
			double b[], 
			int *ldb,
			int *info);

#ifdef TIMING_FLOPS
extern "C" {
    void HPM_Init(void);
    void HPM_Start(char *);
    void HPM_Stop(char *);
    void HPM_Print(void);
    void HPM_Print_Flops(void);
    void HPM_Print_Flops_Agg(void);
    void HPM_Terminate(char*);
}
#endif

  
DeSymIndefSolver::DeSymIndefSolver( DenseSymMatrix * dm )
{
  mStorage = DenseStorageHandle( dm->getStorage() );

  int size = mStorage->n;
  ipiv = new int[size];
  lwork = -1;
  work = NULL;
  sparseMat = 0;
  
}

DeSymIndefSolver::DeSymIndefSolver( SparseSymMatrix * sm )
{
  int size = sm->size();
  mStorage = DenseStorageHandle( new DenseStorage(size,size) );

  ipiv = new int[size];
  lwork = -1;
  work = NULL;
  sparseMat = sm;

  

}


//#include "mpi.h"
void DeSymIndefSolver::matrixChanged()
{
  char fortranUplo = 'U';
  int info;

  int n = mStorage->n;

  if (sparseMat) {
    std::fill(mStorage->M[0],mStorage->M[0]+n*n,0.);

    const double *sM = sparseMat->M();
    const int *jcolM = sparseMat->jcolM();
    const int *krowM = sparseMat->krowM();
    for (int i = 0; i < n; i++) {
      for (int k = krowM[i]; k < krowM[i+1]; k++) {
        int col = jcolM[k];
        mStorage->M[i][col] = sM[k];
      }
    }
  }

#ifdef DENSE_USE_HALF
#ifndef NDEBUG
  for( int i = 0; i < n; i++ )
     for( int j = 0; j < n; j++ )
        assert(j <= i || mStorage->M[i][j] == 0.0);
#endif
#endif


  //query the size of workspace
  lwork=-1;
  double lworkNew;
  FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
	   ipiv, &lworkNew, &lwork, &info );
  
  lwork = (int)lworkNew; 
  if(work) delete[] work;
  work = new double[lwork];

#ifdef TIMING_FLOPS
  HPM_Start("DSYTRFFact");
#endif

  //factorize
  FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
	   ipiv, work, &lwork, &info );

#ifdef TIMING_FLOPS
  HPM_Stop("DSYTRFFact");
#endif
  if(info!=0)
      printf("DeSymIndefSolver::matrixChanged : error - dsytrf returned info=%d\n", info);
  //assert(info==0);

  //int piv2x2=0;
  //for(int i=0; i<n; i++)
  //  if(ipiv[i]<0) piv2x2++;
  //printf("%d 2x2 pivots were used\n", piv2x2);

}

void DeSymIndefSolver::solve ( OoqpVector& v )
{
  char fortranUplo = 'U';
  int info;
  int one = 1;

  int n = mStorage->n; SimpleVector &  sv = dynamic_cast<SimpleVector &>(v);
#ifdef TIMING_FLOPS
  HPM_Start("DSYTRSSolve");
#endif

  FNAME(dsytrs)( &fortranUplo, &n, &one,	&mStorage->M[0][0],	&n,
	   ipiv, &sv[0],	&n,	&info);

#ifdef TIMING_FLOPS
  HPM_Stop("DSYTRSSolve");
#endif
  assert(info==0);
}

void DeSymIndefSolver::solve ( GenMatrix& rhs_in )
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  char fortranUplo = 'U';
  int info;
  int nrows,ncols; rhs.getSize(ncols,nrows);

  int n = mStorage->n;

  FNAME(dsytrs)( &fortranUplo, &n, &ncols,	&mStorage->M[0][0],	&n,
	   ipiv, &rhs[0][0],	&n,	&info);

  assert(info==0);
}

void DeSymIndefSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

DeSymIndefSolver::~DeSymIndefSolver()
{
  delete[] ipiv;
  if(work) delete[] work;
}
