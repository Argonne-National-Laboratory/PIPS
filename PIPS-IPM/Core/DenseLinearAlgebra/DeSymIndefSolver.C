/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "DeSymIndefSolver.h"
#include "SimpleVector.h"
#include <cassert>

#include "DenseSymMatrix.h"


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
  
DeSymIndefSolver::DeSymIndefSolver( DenseSymMatrix * dm )
{
  mStorage = DenseStorageHandle( dm->getStorage() );

  int size = mStorage->n;
  ipiv = new int[size];
  lwork = -1;
  work = NULL;
  
}
//#include "mpi.h"
void DeSymIndefSolver::matrixChanged()
{

  char fortranUplo = 'U';
  int info;

  int n = mStorage->n;

  //!log 
  /*
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0) {
    //!log 
    printf("DeSymIndefSolver::matrixChanged - matrix is:\n");
    for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) 
	printf("%24.16f ", mStorage->M[i][j]);
      printf(";\n");
    }
    //sv.writeToStream(cout);
  } 
  */

  //query the size of workspace
  lwork=-1;
  double lworkNew;
  FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
	   ipiv, &lworkNew, &lwork, &info );
  
  lwork = (int)lworkNew; 
  if(work) delete[] work;
  work = new double[lwork];

  //printf("%d allocated, n being [%d]\n", lwork, n);

  //factorize
  FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
	   ipiv, work, &lwork, &info );

  assert(info==0);

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

  FNAME(dsytrs)( &fortranUplo, &n, &one,	&mStorage->M[0][0],	&n,
	   ipiv, &sv[0],	&n,	&info);

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
