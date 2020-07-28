/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cmath>
#include <cassert>

#include "DenseSymMatrix.h"
#include "DeSymPSDSolver.h"
#include "DeSymIndefSolver.h"
#include "OoqpBlas.h"
#include "SimpleVector.h"

#include "DenseGenMatrix.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"

#include "DoubleMatrixTypes.h"

extern "C" void  dsyrk_(char* UPLO, char* TRANS,
			int* N, int* K,
			double* alpha, double* A, int* lda,
			double* beta,  double* C, int* ldc);

int DenseSymMatrix::isKindOf( int matrixType ) const
{
  return matrixType == kDenseSymMatrix || matrixType == kSymMatrix;
}


DenseSymMatrix::DenseSymMatrix( int size )
{ 
  mStorage = DenseStorageHandle( new DenseStorage( size, size ) );
}


DenseSymMatrix::DenseSymMatrix( double Q[], int size )
{
  mStorage = DenseStorageHandle( new DenseStorage( Q, size, size ) );
}


void DenseSymMatrix::symAtPutZeros( int row, int col,
				     int rowExtent, int colExtent )
{
  mStorage->atPutZeros( row, col, rowExtent, colExtent );
}


void DenseSymMatrix::putSparseTriple( int irow[], int len,
					  int jcol[], double A[], 
					  int& info )
{
  mStorage->putSparseTriple( irow, len, jcol, A, info );
}


void DenseSymMatrix::atAddOuterProductOf( int row, int col, double alpha,
					      double * x, int incx, int nx )
{
  mStorage->atAddOuterProductOf( row, col, alpha, x, incx, nx );
}


void DenseSymMatrix::getDiagonal( OoqpVector& vec )
{
  mStorage->getDiagonal( vec );
}


void DenseSymMatrix::setToDiagonal( OoqpVector& vec )
{
  mStorage->setToDiagonal( vec );
}


void DenseSymMatrix::fromGetSpRow( int row, int col,
				       double A[], int lenA,
				       int jcolA[], int& nnz,
				       int colExtent, int& info )
{
  if( col + colExtent < row + 1 ) {
    mStorage->fromGetSpRow( row, col, A, lenA, jcolA, nnz, colExtent, info );
  } else {
    if( col <= row ) {
      mStorage->fromGetSpRow( row, col, A, lenA, jcolA,
			      nnz, row - col + 1, info );
    }
  }
}


void DenseSymMatrix::getSize( long long& m, long long& n ) const
{
  m = mStorage->m;
  n = mStorage->n;
}

void DenseSymMatrix::getSize( int& m, int& n ) const
{
  m = mStorage->m;
  n = mStorage->n;
}


long long DenseSymMatrix::size()
{ 
   return mStorage->m;
}


void DenseSymMatrix::symAtPutSubmatrix( int destRow, int destCol,
					    DoubleMatrix& Mat,
					    int srcRow, int srcCol,
					    int rowExtent, int colExtent )
{
  int m = mStorage->m, n = mStorage->n;
  double ** M = mStorage->M;

  assert( destRow >= 0 && destRow + rowExtent <= m );
  assert( destCol >= 0 && destCol + colExtent <= n );

  // If assertions are turned off, clip to the actual size of this matrix
  destRow = ( destRow >= 0 ) ? destRow : 0;
  destCol = ( destCol >= 0 ) ? destCol : 0;
  rowExtent = ( destRow + rowExtent <= m ) ?  rowExtent : m - destRow;
  colExtent = ( destCol + colExtent <= n ) ?  colExtent : n - destCol;

  Mat.fromGetDense( srcRow, srcCol, &M[destRow][destCol], n,
		     rowExtent, colExtent );

  
}


void DenseSymMatrix::mult ( double beta,  double y[], int incy,
				double alpha, const double x[], int incx ) const
{
  char fortranUplo = 'U';
  int n = mStorage->n;
  
  dsymv_( &fortranUplo, &n, &alpha, &mStorage->M[0][0], &n,
	  x, &incx, &beta, y, &incy );
}


void DenseSymMatrix::mult ( double beta,  OoqpVector& y_in,
				double alpha, const OoqpVector& x_in ) const
{
  char fortranUplo = 'U';
  int n = mStorage->n;
  SimpleVector & y = (SimpleVector &) y_in;
  SimpleVector & x = (SimpleVector &) x_in;
  int incx = 1, incy = 1;
  
  if( n != 0 ) {
    dsymv_( &fortranUplo, &n, &alpha, &mStorage->M[0][0], &n,
	    &x[0], &incx, &beta, &y[0], &incy );
  } 
}


void DenseSymMatrix::transMult ( double beta,  OoqpVector& y,
				     double alpha, const OoqpVector& x ) const
{
  // We're symmetric silly
  this->mult( beta, y, alpha, x );
}


void DenseSymMatrix::transMult ( double beta,  double y[], int incy,
				     double alpha, const double x[], int incx ) const
{
  // We're symmetric silly
  this->mult( beta, y, incy, alpha, x, incx );
}


double DenseSymMatrix::abmaxnorm()
{
  double norm = 0;
  
  int i, j;
  double ** M = mStorage->M;
  int m = mStorage->m;
  double eltNorm;

  for ( i = 0; i < m; i++ ) {
    for ( j = 0; j <= i; j++ ) {
      eltNorm = fabs( M[i][j] );
      if ( eltNorm > norm ) norm = eltNorm;
    }
  }
  return norm;
}


void DenseSymMatrix::writeToStream(ostream& out) const
{
  int i, j;
  int n = mStorage->n;
  double ** M = mStorage->M;

  for( i = 0; i < n; i++ ) {
    for( j = 0; j <= i && j < n - 1; j++ ) {
      out << M[i][j] << "   ";
    }
    for(      ; j < n - 1; j++ ) {
      out << M[j][i] << "   ";
    }
    if ( j < n )     out << M[j][i];
    if ( i < n - 1 ) out << endl;
  }
}

void DenseSymMatrix::writeToStreamDense(std::ostream& out) const
{
   writeToStream(out);
}


void DenseSymMatrix:: randomizePSD(double * seed)
{  
  int n = mStorage->n;
  double ** M = mStorage->M;
  double drand(double*);
  int i, j, k;

  mStorage->atPutZeros( 0, 0, n, n );
  for(i=0;i<n;i++) {
    for(j=0;j<=i;j++) {
      M[i][j] = drand(seed);
    }
  }

  for( i = n-1; i >= 0; i-- ) {
    for( j = i; j >= 0; j--) {
      M[i][j] = M[i][j] * M[j][j];
      for( k = j - 1; k >= 0; k-- ) {
	M[i][j] += M[i][k]*M[j][k];
      }
    }
  }
}


void DenseSymMatrix::fromGetDense( int row, int col, double * A,
				       int lda,
				       int rowExtent, int colExtent )
{
  int m = mStorage->m, n = mStorage->n;
  double ** M = mStorage->M;

  int i, j;
  assert( row >= 0 && row + rowExtent <= m );
  assert( col >= 0 && col + colExtent <= n );

  // If assertions are turned off, clip to the actual size of this matrix
  row = ( row >= 0 ) ? row : 0;
  col = ( col >= 0 ) ? col : 0;
  rowExtent = ( row + rowExtent <= m ) ?  rowExtent : m - row;
  colExtent = ( col + colExtent <= n ) ?  colExtent : n - col;

  for ( i = row; i < row + rowExtent; i++ ) {
    for ( j = col; j <= i && j < col + colExtent; j++ ) {
      A[(i - row)*lda + j - col] = M[i][j];
    }
    for( ; j < col + colExtent; j++ ) {
      A[(i - row)*lda + j - col] = M[j][i];
    }
  }
}


void DenseSymMatrix::atPutDiagonal( int idiag,
					OoqpVector& v )
{
  mStorage->atPutDiagonal( idiag, v );
}


void DenseSymMatrix::fromGetDiagonal( int idiag,
					  OoqpVector& v )
{
  mStorage->fromGetDiagonal( idiag, v );
}


void DenseSymMatrix::symAtPutSpRow( int row, double A[], int lenA,
					int jcolA[], int& info )
{
  // Lower triangular put
  int lA = lenA;
  while( lA > 0 && jcolA[lA - 1] > row ) lA--;
  if( lA > 0 ) {
    mStorage->atPutSpRow( row, A, lA, jcolA, info );
  } else {
    info = 0;
  }
}

void DenseSymMatrix::symAtPutDense( int row, int col, double * A, int lda,
				     int rowExtent, int colExtent )
{
  mStorage->atPutDense( row, col, A, lda, rowExtent, colExtent );
}

void DenseSymMatrix::SymmetricScale( OoqpVector& vec )
{
  mStorage->SymmetricScale( vec );
}

void DenseSymMatrix::ColumnScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );
}

void DenseSymMatrix::RowScale( OoqpVector& vec )
{
  mStorage->RowScale( vec );
}

void DenseSymMatrix::scalarMult( double num )
{
  mStorage->scalarMult( num );
}

/* updates the upper left block only  if sizes does not matches */
void DenseSymMatrix::matMult(double alpha, 
			     GenMatrix& A_, int transA, 
			     GenMatrix& B_, int transB,
			     double beta)
{
  
  DenseGenMatrix& A = dynamic_cast<DenseGenMatrix&>(A_);
  DenseGenMatrix& B = dynamic_cast<DenseGenMatrix&>(B_);

  // the other way around since fortran stores column-wise and we store row-wise
  char forTransA = (transA==0?'T':'N');
  char forTransB = (transB==0?'T':'N');

  DenseSymMatrix& C = *this;

  int m,n,k,kB; int ldc = mStorage->m;

  if(!transA) {A.getSize(m,k);}
  else        {A.getSize(k,m);}

  if(!transB) {B.getSize(kB,n);}
  else        {B.getSize(n,kB);}

  assert(k==kB);

  assert(mStorage->m >= m);
  assert(mStorage->n >= n);

  double** AA = A.mStorage->M;
  double** BB = B.mStorage->M;
  double** CC = C.mStorage->M;

 
  //Ap[0] = 0.0; Ap[1] = 0.0; Ap[2] = 0.0; Ap[3] = 0.0; Ap[4] = 0.0; 
  //Ap[5] = 0.0; Ap[6] = 0.0; Ap[7] = 0.0; Ap[8] = 0.0; Ap[9] = 0.0; 
  //Ap[10]= 0.0; Ap[11]= 0.0; Ap[12]= 1.0; Ap[13]= 1.0; Ap[14]= 0.0; 
  //Ap[15]= 0.0; Ap[16]= 0.0; Ap[17]= 0.0; Ap[18]= 0.0; Ap[19]= 0.0; 

 
//   printf("\nA is \n");
//   for(int i=0; i<k; i++) {
//     for(int j=0; j<m; j++) 
//       printf("%7.4f ", AA[i][j]);
//     printf("\n");
//   }
 
//   printf("\nB is\n");
//   for(int i=0; i<k; i++) {
//     for(int j=0; j<n; j++) 
//       printf("%g ", BB[i][j]);
//     printf("\n");
//   }
  

  dgemm_(&forTransA, &forTransB, &m, &n, &k,
	 &alpha, 
	 &AA[0][0], &m,
	 &BB[0][0], &n,
	 &beta,
	 &CC[0][0], &ldc);


//   //!log
//   printf("\nUtV is\n");
//   for(int i=0; i<m; i++) {
//     for(int j=0; j<n; j++)
//       printf("%g ", CC[i][j]);
//     printf("\n");

//   }
}



void DenseSymMatrix::symAtPutSubmatrix( int destRow, int destCol,
					DoubleMatrix& Mat,
					int srcRow, int srcCol,
					int rowExtent, int colExtent,
					int forceSymUpdate)
{

  if(forceSymUpdate==0) {
    symAtPutSubmatrix(destRow, destCol, 
		      Mat, srcRow, srcCol, rowExtent, colExtent);

    return;
  }

  int m = mStorage->m, n = mStorage->n;
  double ** M = mStorage->M;

  assert( destRow >= 0 && destRow + rowExtent <= m );
  assert( destCol >= 0 && destCol + colExtent <= n );

  // If assertions are turned off, clip to the actual size of this matrix
  destRow = ( destRow >= 0 ) ? destRow : 0;
  destCol = ( destCol >= 0 ) ? destCol : 0;
  rowExtent = ( destRow + rowExtent <= m ) ?  rowExtent : m - destRow;
  colExtent = ( destCol + colExtent <= n ) ?  colExtent : n - destCol;

  Mat.fromGetDense( srcRow, srcCol, &M[destRow][destCol], n,
		    rowExtent, colExtent );

  //!exec
  for(int i=destRow; i<destRow + rowExtent; i++) {
    for(int j=destCol; j<destCol + colExtent; j++) {
      M[j][i] = M[i][j]; 
    }
  }

  
}

void DenseSymMatrix::atRankkUpdate( double alpha, double beta, DenseGenMatrix& U, int trans)
{
  int n, k; int ldu, lda;
  //-----------------------------------------------
  // setup if the U is stored in column-major form
  // (FORTRAN Style)
  //   char UPLO  = 'U'; char TRANS = trans==0?'N':'T';
  
  //   U.getSize(n,k); ldu=n;
  //   if(trans) k=n;
  
  //   n = mStorage->n; 
  //   lda=n;
  //----------------------------------------------

  // U and 'this' are stored in row-major form -> a little change in passing params to FORTRAN is needed
  char UPLO  = 'U'; //update LOWER triagular part for symmetric matrix 'this'
  //trans=1 -> this += U'*U -> tell BLAS to do U*U'
  //trans=0 -> this += U*U' -> tell BLAS to do U'*U
  char TRANS = trans==0?'T':'N';  

  int m;
  U.getSize(m,k);
  ldu=k; // change leading dim so that U in row-major  in col-major
  if(trans) k=m;

  n = mStorage->n; 
  lda=n;

#ifdef DEBUG
  //TRANS = 'N', k specifies the number of columns of the matrix U
  //we pass U' instead of U, so k should be the number of rows of U
  int r,c; U.getSize(r,c);
  if(TRANS=='N') assert(k==r);
  else if(TRANS=='T') assert(k==c);
  else assert(false);
#endif


  dsyrk_(&UPLO, &TRANS,
	 &n, &k,
	 &beta,   &U.getStorageRef().M[0][0], &ldu,
	 &alpha,  &mStorage->M[0][0],       &lda);
}

int DenseSymMatrix::getNumberOfNonZeros() const
{
   assert( mStorage->m == mStorage->n );
   int nnz = 0;
   for( int i = 0; i < mStorage->m; i++ )
   {
      for( int j = i + 1; j < mStorage->n; j++ )
         if( mStorage->M[i][j] != 0.0 )
            nnz++;
      nnz++; //always have diags
   }
   return nnz;
}
