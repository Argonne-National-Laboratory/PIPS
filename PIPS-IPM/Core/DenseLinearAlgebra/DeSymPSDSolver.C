/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "DeSymPSDSolver.h"
#include <cassert>

#include "DenseStorage.h"
#include "DenseSymMatrix.h"
#include "SimpleVector.h"

// declarations for LAPACK functions used to factor/solve:

// dsytrf_() factors a symmetric indefinite matrix A, see LAPACK 
// documentation for more details.
extern "C" void dpotrf_( char * uplo, int * n,
			 double * A, int * lda, int * info );
 
extern "C" void dpotrs_( char * uplo, int * n, int * nrhs,
			 double * A, int * lda,
			 double * b, int * ldb, int * info ) ;

extern "C" void dtrsm_(char* side, 
		       char* uplo, 
		       char* transa,
		       char* diag,
		       int* m, int* n,
		       double* alpha,
		       double* A, int* LDA,
		       double* B, int* LDB);

  
DeSymPSDSolver::DeSymPSDSolver( DenseSymMatrix * dsm )
{
  mStorage = DenseStorageHandle( dsm->getStorage() );
}

void DeSymPSDSolver::matrixChanged()
{
  char fortranUplo = 'U';
  int info;

  int n = mStorage->n;

  dpotrf_( &fortranUplo, &n, &mStorage->M[0][0], &n, &info );

  assert(info==0);
}

void DeSymPSDSolver::solve ( OoqpVector& x_in )
{
  char fortranUplo = 'U';
  int info;
  int one = 1;

  int n = mStorage->n;
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  
  dpotrs_( &fortranUplo, &n, &one,	&mStorage->M[0][0],	&n,
	   &x[0], &n, &info);
}

void DeSymPSDSolver::Lsolve( DenseGenMatrix& B)
{
  /* 

  double A[9];
  //A[0]=1;A[3]=0;A[6]=0;
  //A[1]=2;A[4]=1;A[7]=0;
  //A[2]=3;A[5]=1;A[8]=1;
  A[0]=2;A[1]=0;A[2]=0;
  A[3]=2;A[4]=2;A[5]=0;
  A[6]=3;A[7]=1;A[8]=2;
  
  DenseSymMatrix AA(A,3);

  double b[6];
  //b[0]=1;  b[3]=1;
  //b[1]=3;  b[4]=2;
  //b[2]=5;  b[5]=3;
  b[0]=2;  b[1]=2;
  b[2]=4;  b[3]=2;
  b[4]=6;  b[5]=3;

  DenseGenMatrix BB(b, 3,2);
  //for(int i=0; i<3; i++) { for(int j=0; j<2; j++) printf("%5.3f ", BB[i][j]); printf("\n"); }; printf("----------------\n");
 
  // side will be 'R'
  char sidet='R';    //solve AX=B (L) or XA=B (R)

  //lower triangular but we say upper so FORTRAN will access the right part.
  char uplot='U';    //Upper or Lower triangle solve

  //X^T L^T=B^T 
  char transAt='N';  //solve with L or L'

  char diagt='N';    //is or is  not unit triangular

  int nt=3; int colst=2; 
  
  double onet = 1.0;
  
  int col=1; 
  SimpleVector rest(nt); 
  for(int i=0; i<nt; i++)
    rest[i]=BB[i][col-1];
  //printf("rhs:"); for(int i=0; i<nt; i++) printf("%5.3f ", rest[i]); printf("\n----------------\n");

  dtrsm_(&sidet, &uplot, &transAt, &diagt,
	 &colst,
	 &nt,
 	 &onet, 
	 &AA.getStorage()->M[0][0], &nt, 
	 &BB.getStorage()->M[0][0], 
	 &colst);

  //for(int i=0; i<3; i++) { for(int j=0; j<2; j++) printf("%5.3f ", BB[i][j]); printf("\n"); }; printf("----------------\n");

  SimpleVector solt(nt); 
  for(int i=0; i<nt; i++)
    solt[i]=BB[i][col-1];

  //printf("sol:"); for(int i=0; i<nt; i++) printf("%5.3f ", solt[i]); printf("\n----------------\n");

  
  DenseGenMatrix AAA(A, 3,3);
  //for(int i=0; i<3; i++) { for(int j=0; j<3; j++) printf("%5.3f ", AAA[i][j]); printf("\n"); }; printf("----------------\n");
  AAA.mult(-1.0,rest, 1.0,solt);

  //for(int i=0; i<nt; i++) printf("%5.3f ", rest[i]); printf("\n----------------\n");

  double nrm = rest.twonorm();
  printf("nrm=%g\n", nrm);

  int nn=mStorage->n;
  double* buf = new double[mStorage->n*mStorage->n];
  for(int i=0; i<mStorage->n; i++)
    for(int j=0; j<mStorage->n;j++)
      buf[i+j*nn] = mStorage->M[i][j];

  
  int mm; B.getSize(mm,nn);
  double* buf2 = new double[mm*nn];
  for(int i=0; i<mm; i++)
    for(int j=0; j<nn;j++)
      buf2[i+j*mm] = B[i][j];
  */


  //R is in row-major format, FORTRAN works with column-major
  //Passing the correct params is a PITA
  //
  // We want to solve LX=B, BU=T we will instruct BLAS to solve X^T L^T=B^T
  // This is because the different storage schemes between FORTRAN and C/C++.
  
  // side will be 'R'
  char side='R';    //solve AX=B (L) or XA=B (R)

  //lower triangular but we say upper so FORTRAN will access the right part.
  char uplo='U';    //Upper or Lower triangle solve

  //X^T L^T=B^T 
  char transA='N';  //solve with L or L'

  char diag='N';    //is or is  not unit triangular

  int n,cols; B.getSize(n, cols); assert(n==mStorage->n);
  int ldb = cols;
  
  /*
  SimpleVector res(n);  col=10;
  for(int i=0; i<n; i++) res[i]=B[i][col-1];
  */

  double one = 1.0;
  
  dtrsm_(&side, &uplo, &transA, &diag,
	 &cols,  //number of rows of B, lie and give the columns
	 &n,     //number of cols of B, lie again 
 	 &one, 
	 &mStorage->M[0][0], &n, 
	 &B.getStorageRef().M[0][0], 
	 &ldb);

  /*SimpleVector sol(n); 
  for(int i=0; i<n; i++) sol[i]=B[i][col-1];
  
  DenseGenMatrix THIS(&mStorage->M[0][0], n, n);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      if(i<j) THIS[i][j]=0.0;
  THIS.mult(-1.0, res, 1.0, sol);
  nrm=res.twonorm(); printf("NRM=%g\n", nrm);*/

  //assert(false);

}

void DeSymPSDSolver::diagonalChanged( int /* idiag */, int /*extent*/ )
{
  this->matrixChanged();
}

DeSymPSDSolver::~DeSymPSDSolver()
{
}
