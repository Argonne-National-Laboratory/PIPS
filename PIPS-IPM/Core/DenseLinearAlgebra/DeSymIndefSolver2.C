/* PIPS                                                               *
 * Authors: Miles Lubin                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "DeSymIndefSolver2.h"
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

extern "C" void FNAME(dpotrf)(char *uplo,
      int *n, double *A, int*lda, int*info);

extern "C" void FNAME(dtrsm)(char *side,
      char *uplo, char *transa, char *diag,
      int* m, int *n, double *alpha, double *A,
      int *lda, double *B, int *ldb);

extern "C" void FNAME(dsyrk)(char *uplo,
      char *trans, int* n, int* k, double *alpha,
      double *A, int *lda, double *beta,
      double *C, int *ldc);

extern "C" void FNAME(dtrsv)(char *uplo,
      char *trans, char* diag, int *n,
      double *A, int *lda, double *x,
      int *incx);

extern "C" void FNAME(dscal)(int *n,
      double *alpha, double *x, int *incx);


DeSymIndefSolver2::DeSymIndefSolver2( DenseSymMatrix * dm, int nx ) : nx(nx)
{
  mStorage = dm->getStorageHandle();

  n = mStorage->n;
  ny = n - nx;

}
//#include "mpi.h"
void DeSymIndefSolver2::matrixChanged()
{

  char fortranUplo = 'U';
  char fortranL = 'L';
  char fortranN = 'N';
  char fortranT = 'T';
  int info;
  double *mat = &mStorage->M[0][0];
  double one = 1, zero = 0;

  // treat matrix as column-major and assume upper
  // triangle is filled
  // this is different from the elemental version
  // which uses the lower triangle

  /*
  start with:
  [Q A^T
   *  0]

  cholesky on Q to get M^T
  trsm to get M^-1A^T
  syrk to form (M^-1A^T)^T(M^-1A^T)
  cholesky on bottom right
  */


  FNAME(dpotrf)(&fortranUplo,&nx,mat,&n,&info);
  if (info != 0) {
    cerr << "error factoring Q block: info = " << info << endl;
  }
  assert(info==0);
  if (ny == 0) return;

  printf("dtrsm\n");

  FNAME(dtrsm)(&fortranL,&fortranUplo,&fortranT,
    &fortranN,&nx,&ny,&one,mat,&n,
    mat+nx*n,&n);

  printf("dsyrk\n");

  FNAME(dsyrk)(&fortranUplo,&fortranT,&ny,&nx,&one,
    mat+nx*n,&n,&zero,mat+nx*n+nx,&n);

  printf("dpotrf2\n");

  FNAME(dpotrf)(&fortranUplo,&ny,mat+nx*n+nx,&n,&info);
  if (info != 0) {
    cerr << "error factoring AQ^-1A^T block: info = " << info << endl;
  }
  assert(info==0);

  printf("finished factorization\n");


}

void DeSymIndefSolver2::solve ( OoqpVector& v )
{
  char fortranUplo = 'U';
  char fortranT = 'T';
  char fortranN = 'N';

  double *mat = &mStorage->M[0][0];
  int one = 1;
  double minus1 = -1;

  SimpleVector &  sv = dynamic_cast<SimpleVector &>(v);
  double *rhs = &sv[0];

  FNAME(dtrsv)(&fortranUplo,&fortranT,&fortranN,
    &n,mat,&n,rhs,&one);

  if (ny > 0)
    FNAME(dscal)(&ny,&minus1,rhs+nx,&one);

  FNAME(dtrsv)(&fortranUplo,&fortranN,&fortranN,
    &n,mat,&n,rhs,&one);

  printf("finished rhs solve\n");

}

void DeSymIndefSolver2::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

DeSymIndefSolver2::~DeSymIndefSolver2()
{
}
