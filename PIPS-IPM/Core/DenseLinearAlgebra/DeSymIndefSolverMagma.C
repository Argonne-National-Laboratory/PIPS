/* PIPS-IPM                                                           *
 * Authors: Cosmin Petra, Olaf Schenk (USI), Dmitry Mikushin (USI)    *
 * (C) 2013 Argonne National Laboratory. See Copyright Notification   */

#include "pipsport.h"
#include "DeSymIndefSolverMagma.h"
#include "SimpleVector.h"
#include <cassert>

#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"

#ifdef GPUCODE
#include <cuda_runtime.h>
#include <cstdio>
#define CUDA_SAFE_CALL(x) \
  do { cudaError_t err = x; if (err != cudaSuccess) { \
  fprintf(stderr, "Error \"%s\" at %s:%d\n", cudaGetErrorString(err), \
	  __FILE__, __LINE__); return -1; \
    }} while (0);
#endif

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

#ifdef GPUCODE
extern "C" void FNAME(magmaf_dgetrf_gpu)(int* m,
				   int *n,
				   double** A,
				   int *lda,
				   int* ipiv,
				   int *info);

extern "C" void FNAME(magmaf_dgetrs_gpu)(char* trans,
				   int *n,
				   int *nrhs,
				   double** A,
				   int *lda,
				   int* ipiv,
				   double** b,
				   int *ldb,
				   int *info);
#else
extern "C" void FNAME(dgetrf)(int* m,
                        int *n,
                        double* A,
                        int *lda,
                        int* ipiv,
                        int *info);

extern "C" void FNAME(dgetrs)(char* trans,
			int *n,
			int *nrhs,
			double* A,
			int *lda,
			int* ipiv,
			double* b,
			int *ldb,
			int *info);
#endif

DeSymIndefSolverMagma::DeSymIndefSolverMagma( DenseSymMatrix * dm )
{

  cout << "creating GPUCODE solver" << endl;
  mStorage = dm->getStorageHandle();

  int size = mStorage->n;
  ipiv = new int[size];
  //lwork = -1;
  //work = nullptr;
  sparseMat = 0;


#ifdef GPUCODE
  mFact_gpu = nullptr;
  //CUDA_SAFE_CALL(cudaMalloc(&mFact_gpu, sizeof(double) * size * size));
  cudaMalloc(&mFact_gpu, sizeof(double) * size * size);
  assert(mFact_gpu!=nullptr);

  mRhs_gpu = nullptr;
  cudaMalloc(&mRhs_gpu, sizeof(double) * size);
  assert(mRhs_gpu!=nullptr);
#else
  cout << "DeSymIndefSolverMagma used but GPU code is not enabled" << endl;
#endif

  cout << "MAGMA: size of matrix is " << size << "done with MAGMA constructor" << endl;
}

DeSymIndefSolverMagma::DeSymIndefSolverMagma( SparseSymMatrix * sm )
{
  assert(false && "Not available! Can be implemented though :)");
  /*
    cout << "Magma solver : sparse matrix -------------" << endl;
  int size = sm->size();
  mStorage = DenseStorageHandle( new DenseStorage(size,size) );

  ipiv = new int[size];
  sparseMat = sm;
  */
}


//#include "mpi.h"
void DeSymIndefSolverMagma::matrixChanged()
{
  cout << "Magma::matrixChanged()" << endl;
  int info;

  int n = mStorage->n;

  //factorize
#ifdef GPUCODE
  cudaMemcpy(mFact_gpu,
	     &mStorage->M[0][0], sizeof(double) * n * n,
	     cudaMemcpyHostToDevice);
  FNAME(magmaf_dgetrf_gpu)( &n, &n, &mFact_gpu, &n, ipiv, &info );
#else
  FNAME(dgetrf)( &n, &n, &mStorage->M[0][0], &n,
	   ipiv, &info );
#endif
  if(info!=0)
      printf("DeSymIndefSolverMagma::matrixChanged : error - dsytrf returned info=%d\n", info);
  assert(info==0);
  cout << "MAGMA factorization done." << endl;
}

void DeSymIndefSolverMagma::solve ( OoqpVector& v )
{
  cout << "magma::solve ( OoqpVector& v )" << endl;
  char trans = 'N';
  int info;
  int one = 1;

  int n = mStorage->n; SimpleVector &  sv = dynamic_cast<SimpleVector &>(v);

#ifdef GPUCODE
  cudaMemcpy(mRhs_gpu, &sv[0], sizeof(double) * n, cudaMemcpyHostToDevice);
  FNAME(magmaf_dgetrs_gpu)(&trans, &n, &one, &mFact_gpu, &n, ipiv, &mRhs_gpu, &n, &info);
  cudaMemcpy(&sv[0], mRhs_gpu, sizeof(double) * n, cudaMemcpyDeviceToHost);
#else
  FNAME(dgetrs)( &trans, &n, &one,	&mStorage->M[0][0],	&n,
			ipiv, &sv[0],	&n,	&info);
#endif
  if(info!=0) cout << "dgetrs returned info=" << info << endl;
  assert(info==0);
  cout << "MAGMA : solve done" << endl;
}

void DeSymIndefSolverMagma::solve ( GenMatrix& rhs_in )
{
  assert(false && "Not implemented");
  /*
  assert(false);
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  char trans = 'U';
  int info;
  //RANS   (input) CHARACTER*1
  //          Specifies the form of the system of equations:
  //          = 'N':  A * X = B  (No transpose)
  //          = 'T':  A'* X = B  (Transpose)
  //          = 'C':  A'* X = B  (Conjugate transpose = Transpose)

  int nrows,ncols; rhs.getSize(ncols,nrows);

  int n = mStorage->n;

  FNAME(magmaf_dgetrs_gpu)( &trans, &n, &ncols,	&mStorage->M[0][0],	&n,
	   ipiv, &rhs[0][0],	&n,	&info);

  if(info!=0) cout << "dgetrs (nrhs) returned info=" << info << endl;
  assert(info==0);
  */
}

void DeSymIndefSolverMagma::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

DeSymIndefSolverMagma::~DeSymIndefSolverMagma()
{
  delete[] ipiv;
#ifdef GPUCODE
  cudaFree(mFact_gpu);
  cudaFree(mRhs_gpu);
#endif
  //if(work) delete[] work;
}
