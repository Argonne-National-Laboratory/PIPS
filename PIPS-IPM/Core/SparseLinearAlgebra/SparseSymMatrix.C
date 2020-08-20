/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseSymMatrix.h"
#include "SparseStorage.h"
#include <cassert>
#include <cmath>
#include "SimpleVector.h"
#include "pipsport.h"

#include "DoubleMatrixTypes.h"

int SparseSymMatrix::isKindOf( int type ) const
{
  return type == kSparseSymMatrix || type == kSymMatrix;
}

SparseSymMatrix::SparseSymMatrix()
 : isLower(true)
{
  mStorage = nullptr;
}

SparseSymMatrix::SparseSymMatrix( int size, int nnz, bool isLower )
: isLower(isLower)
{
  mStorage = SparseStorageHandle( new SparseStorage(size, size, nnz) );
}


/*SparseSymMatrix::SparseSymMatrix(const std::vector<SparseSymMatrix*> &blocks)
{
  vector<SparseStorage*> v(blocks.size());
  for (size_t i = 0; i < blocks.size(); i++) v[i] = blocks[i]->mStorage;
  mStorage = SparseStorageHandle(new SparseStorage(v,true)); // must be diagonal
}
*/

SparseSymMatrix::SparseSymMatrix( int size, int nnz,
				  int krowM[], int jcolM[], double M[],
				  int deleteElts, bool isLower )
: isLower(isLower)
{
  mStorage = SparseStorageHandle( new SparseStorage(size, size,
						    nnz, krowM, jcolM, M,
						    deleteElts) );
}


void SparseSymMatrix::putSparseTriple( int irow[], int len,
					   int jcol[], double A[], 
					   int& info )
{
  mStorage->putSparseTriple( irow, len, jcol, A, info );
}


void SparseSymMatrix::fromGetDense( int row, int col, double * A, int lda,
					int rowExtent, int colExtent )
{
  mStorage->fromGetDense( row, col, A, lda, rowExtent, colExtent );
}


void SparseSymMatrix::getDiagonal( OoqpVector& vec )
{
  mStorage->getDiagonal( vec );
}

void SparseSymMatrix::setToDiagonal( OoqpVector& vec )
{
  mStorage->setToDiagonal( vec );
}

void SparseSymMatrix::symAtPutSpRow( int row, double A[],
					 int lenA, int jcolA[],
					 int& info )
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

void SparseSymMatrix::symPutZeroes()
{
   assert(mStorage);
   mStorage->clear();
}

void SparseSymMatrix::fromGetSpRow( int row, int col,
					double A[], int lenA,
					int jcolA[], int& nnz,
					int colExtent, int& info )
{
  mStorage->fromGetSpRow( row, col, A, lenA, jcolA, nnz, colExtent, info );
}


void SparseSymMatrix::randomizePSD(double * seed)
{
  int k, NN, chosen, icurrent;
  int nnz;

  double drand( double * );
  int n       = mStorage->n;
  double * M  = mStorage->M;
  int * krowM = mStorage->krowM;
  int * jcolM = mStorage->jcolM;

  // We will always have non-zeros on the diagonal, so there
  // is no randomness there. In fact, choose the (0,0) element now
  krowM[0] = 0;
  jcolM[0] = 0;
  M[0]     = 1e-8 + drand( seed );
  krowM[1] = 1;
  nnz      = 1;

  // Knuth's algorithm for choosing len elements out of NN elts.
  // NN here is the number of elements in the strict lower triangle.
  NN        = n * ( n - 1 )/ 2;
  // len is the number of elements that can be stored, minus the number
  // of elements in the diagonal, which will always be in the matrix.
  int len = mStorage->len - n; 
  // but never more than NN elts
  len = (len <= NN) ? len : NN;

  // chosen is the number of elements that have already been chosen (now 0)
  chosen    = 0;
  // nnz is the number of non-zeros in the matrix (now 1, because the
  // (0,0) element is already in the matrix.
  nnz       = 1;
  // icurrent is the index of the last row whose start has been stored in 
  // krowM;
  icurrent  = 1;
  for ( k = 0; k < NN; k++ ) {
    double r = drand( seed );
	
    if( (NN - k) * r < len - chosen ) {
      // Element k is chosen. What row is it in?
      // In a lower triangular matrix (including a diagonal), it will be in
      // the largest row such that row ( row + 1 ) / 2 < k. In other words
      int row = (int) floor( ( -1 + sqrt( 1.0 + 8.0 * k ) ) / 2 );
      // and its column will be the remainder
      int col = k - row * (row + 1)/2;
      // but since we are only filling in the *strict* lower triangle of 
      // the matrix, we shift the row by 1
      row++;

      if ( row > icurrent ) {
	// We have chosen a row beyond the current row. 
	// Choose a diagonal elt for each intermediate row and fix the
	// data structure.
	for ( ; icurrent < row; icurrent++ ) {
	  // Choose the diagonal
	  M[nnz] = 0.0;
	  int ll;
	  for( ll = krowM[icurrent]; ll < nnz; ll++ ) {
	    M[nnz] += fabs( M[ll] );
	  }
	  M[nnz] +=  1e-8 + drand( seed );
	  jcolM[nnz] = icurrent;


	  nnz++;
	  krowM[icurrent + 1] = nnz;
	}
      } // end if we have chosen a row beyond the current row;
      M[nnz]     = drand(seed);
      jcolM[nnz] = col;
      // add the value of this element (which occurs symmetrically in the 
      // upper triangle) to the appropriate diagonal element
      M[ krowM[col+1] - 1 ] += fabs( M[nnz] );

      nnz++; // We have added another element to the matrix
      chosen++; // And finished choosing another element.
    }  	
  }
  // and of course, we must choose all remaining diagonal elts.
  for ( ; icurrent < n; icurrent++ ) {
    // Choose the diagonal
    M[nnz] = 0.0;
    int ll;
    for( ll = krowM[icurrent]; ll < nnz; ll++ ) {
      M[nnz] += fabs( M[ll] );
    }
    M[nnz] +=  1e-8 + drand( seed );
    jcolM[nnz] = icurrent;

    nnz++;
    krowM[icurrent + 1] = nnz;
  }

}


void SparseSymMatrix::symAtPutSubmatrix( int destRow, int destCol,
					     DoubleMatrix& M,
					     int srcRow, int srcCol,
					     int rowExtent, int colExtent )
{
  int i, k;
  int info, nnz;

  int *    ja = new int[colExtent];
  double * a = new double[colExtent];

  nnz = 0;
  for ( i = 0; i < rowExtent; i++ ) {
    M.fromGetSpRow( srcRow + i, srcCol, a, colExtent, ja,
		     nnz, colExtent, info );
    for( k = 0; k < nnz; k++ ) {
      ja[k] += (destCol - srcCol);
    }
    this->symAtPutSpRow( destRow + i, a, nnz, ja, info );
    assert( info == 0 );
  }

  delete [] a;
  delete [] ja;
}

// Pass these to storage
void SparseSymMatrix::getSize( long long& m, long long& n ) const
{
  int mint, nint;
  mStorage->getSize( mint, nint );
  m=mint; n=nint;
}

void SparseSymMatrix::getSize( int& m, int& n ) const
{
  mStorage->getSize( m, n);
}


long long SparseSymMatrix::size()
{
  return mStorage->rows();
}
void SparseSymMatrix::mult ( double beta,  OoqpVector& y_in,
				 double alpha, const OoqpVector& x_in ) const
{
  const SimpleVector & x = dynamic_cast<const SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert( x.n == mStorage->n &&  y.n == mStorage->m );
  
  const double *xv = nullptr;
  double* yv = nullptr;
  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  this->mult( beta, yv, 1, alpha, xv, 1 );
}

void SparseSymMatrix::transMult ( double beta,   OoqpVector& y_in,
				      double alpha,  const OoqpVector& x_in ) const
{
  const SimpleVector & x = dynamic_cast<const SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert( x.n == mStorage->n &&  y.n == mStorage->m );
  
  const double *xv = nullptr;
  double *yv = nullptr;
  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  this->mult( beta, yv, 1, alpha, xv, 1 );
}
  
void SparseSymMatrix::transMult ( double beta,  double y[], int incy,
				      double alpha, const double x[], int incx ) const
{
  this->mult( beta, y, incy, alpha, x, incx );
}

double SparseSymMatrix::abmaxnorm() const
{
  return mStorage->abmaxnorm();
}

void SparseSymMatrix::writeToStream( std::ostream& out ) const
{
  mStorage->writeToStream( out );
}

void SparseSymMatrix::writeToStreamDense( std::ostream& out ) const
{
  mStorage->writeToStreamDense( out );
}


void SparseSymMatrix::writeToStreamDenseRow( std::stringstream& out, int row ) const
{
   if( mStorage->n > 0 )
   {
      assert(row < mStorage->m);
      mStorage->writeToStreamDenseRow(out, row);
   }
}

std::string SparseSymMatrix::writeToStreamDenseRow( int row ) const
{
   std::stringstream out;

   if( mStorage->n > 0 )
   {
      assert(row < mStorage->m);
      mStorage->writeToStreamDenseRow(out, row);
   }

   return out.str();
}

void SparseSymMatrix::mult ( double beta,  double y[], int incy,
				 double alpha, const double x[], int incx ) const
{
  int m, n, i, j, k;
  this->getSize(m, n); 

  int * jcolM = mStorage->jcolM;
  int * krowM = mStorage->krowM;
  double * M  = mStorage->M;

  for ( i = 0; i < m; i++ ) {
    y[i * incy] *= beta;
  }
  for ( i = 0; i < n; i++ ) {
    for( k = krowM[i]; k < krowM[i+1]; k++ ) {
      j = jcolM[k];
	  
      y[i * incy] += alpha * M[k] * x[j * incx];
      // todo fixme won't work with Q or any other "really" symmetric matrix!
      // Necessary because CtDC from sLinsysRootAug is stored as a general matrix
      if ( i != j && 0 ) {
	y[j * incy] += alpha * M[k] * x[i * incx];
      }
    }
  }
}


void SparseSymMatrix::atPutDiagonal( int idiag, OoqpVector& v )
{
  mStorage->atPutDiagonal( idiag, v );
}

void SparseSymMatrix::fromGetDiagonal( int idiag, OoqpVector& v )
{
  mStorage->fromGetDiagonal( idiag, v );
}

void SparseSymMatrix::SymmetricScale( OoqpVector& vec )
{
  mStorage->SymmetricScale( vec );
}

void SparseSymMatrix::ColumnScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );
}

void SparseSymMatrix::RowScale( OoqpVector& vec )
{
  mStorage->RowScale( vec );
}

void SparseSymMatrix::scalarMult( double num )
{
  mStorage->scalarMult( num );
}

void SparseSymMatrix::reduceToLower()
{
  mStorage->reduceToLower();
}


void SparseSymMatrix::deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec)
{
   const SimpleVectorBase<int>& vec = dynamic_cast<const SimpleVectorBase<int>&>(nnzVec);
   mStorage->deleteEmptyRowsCols(vec.elements(), vec.elements());
}

void SparseSymMatrix::getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const
{
   mStorage->getSparseTriplet_c2fortran(irn, jcn, val);
}

void SparseSymMatrix::getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const
{
   mStorage->getSparseTriplet_fortran2fortran(irn, jcn, val);
}


void SparseSymMatrix::deleteZeroRowsCols(int*& new2orgIdx)
{
   mStorage->deleteZeroRowsColsSym(new2orgIdx);
}

