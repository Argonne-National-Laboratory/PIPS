/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseGenMatrix.h"
#include "SparseStorage.h"
#include <cassert>
#include "SimpleVector.h"

#include "DoubleMatrixTypes.h"
#include <limits>
#include "SparseSymMatrix.h"

int SparseGenMatrix::isKindOf( int type )
{
  return type == kSparseGenMatrix || type == kGenMatrix;
}


SparseGenMatrix::SparseGenMatrix( int rows, int cols, int nnz )
  : m_Mt(NULL)
{
  mStorage = SparseStorageHandle( new SparseStorage( rows, cols, nnz ) );
}


SparseGenMatrix::SparseGenMatrix( int rows, int cols, int nnz,
				  int krowM[], int jcolM[],  double M[],
				  int deleteElts)
  : m_Mt(NULL)
{
  //cout << "SparseGenMatrix1  " << rows << " " << cols << " " << nnz << endl;
  mStorage = SparseStorageHandle( new SparseStorage( rows, cols,nnz, 
						     krowM, jcolM, M,
						     deleteElts) );
}

/*SparseGenMatrix::SparseGenMatrix(const std::vector<SparseGenMatrix*> &blocks, bool diagonal)
  : m_Mt(NULL)
{
  vector<SparseStorage*> v(blocks.size());
  for (size_t i = 0; i < blocks.size(); i++) v[i] = blocks[i]->mStorage;
  mStorage = SparseStorageHandle(new SparseStorage(v,diagonal));
}
*/

SparseGenMatrix::~SparseGenMatrix()
{
  //cout << "~~~~~~~~~SparseGenMatrix " << mStorage->mRefs  << endl;
  if(m_Mt) delete m_Mt;
  
}


void SparseGenMatrix::atPutDense( int row, int col, double * A, int lda,
				      int rowExtent, int colExtent )
{
  mStorage->atPutDense( row, col, A, lda, rowExtent, colExtent );
  assert(m_Mt == NULL);
}


void SparseGenMatrix::fromGetDense( int row, int col, double * A, int lda,
					int rowExtent, int colExtent )
{
  mStorage->fromGetDense( row, col, A, lda, rowExtent, colExtent );
}
  

void SparseGenMatrix::fromGetSpRow( int row, int col,
				    double A[], int lenA,
				    int jcolA[], int& nnz,
				    int colExtent, int& info )
{
  mStorage->fromGetSpRow( row, col, A, lenA, jcolA, nnz,
			  colExtent, info );
}


void SparseGenMatrix::putSparseTriple( int irow[], int len,
					   int jcol[], double A[], 
					   int& info )
{
  mStorage->putSparseTriple( irow, len, jcol, A, info );
  assert(m_Mt == NULL);
}


void SparseGenMatrix::writeToStream(ostream& out) const
{
  mStorage->writeToStream( out );
}

void SparseGenMatrix::writeToStreamDense(ostream& out) const
{
  mStorage->writeToStreamDense( out );
}

void SparseGenMatrix::writeToStreamDenseRow( stringstream& out, int rowidx) const
{
   // check if mStorage might be empty
   if(mStorage->n > 0){
      assert( rowidx < mStorage->m);
      mStorage->writeToStreamDenseRow( out, rowidx );
   }
}

std::string SparseGenMatrix::writeToStreamDenseRow( int rowidx) const
{

   stringstream out;
   // check if mStorage might be empty
   if(mStorage->n > 0){
      assert( rowidx < mStorage->m);
      mStorage->writeToStreamDenseRow( out, rowidx );
   }
   return out.str();
}

void SparseGenMatrix::randomize( double alpha, double beta, double * seed )
{
  mStorage->randomize( alpha, beta, seed );
  assert(m_Mt == NULL);
}


void SparseGenMatrix::getDiagonal( OoqpVector& vec )
{
  mStorage->getDiagonal( vec );
}


void SparseGenMatrix::setToDiagonal( OoqpVector& vec )
{
  mStorage->setToDiagonal( vec );
  assert(m_Mt == NULL);
}


void SparseGenMatrix::atPutSpRow( int row, double A[],
				      int lenA, int jcolA[], int& info )
{
  mStorage->atPutSpRow( row, A, lenA, jcolA, info );
  assert(m_Mt == NULL);
}


int SparseGenMatrix::numberOfNonZeros()
{
  return mStorage->numberOfNonZeros();
}


void SparseGenMatrix::symmetrize( int& info ) 
{
  mStorage->symmetrize( info );
  assert(m_Mt == NULL);
}


void SparseGenMatrix::getSize( long long& m, long long& n )
{
  m = mStorage->m;
  n = mStorage->n;
}
void SparseGenMatrix::getSize( int& m, int& n )
{
  m = mStorage->m;
  n = mStorage->n;
}


void SparseGenMatrix::atPutSubmatrix( int destRow, int destCol,
					  DoubleMatrix& M,
					  int srcRow, int srcCol,
					  int rowExtent, int colExtent )
{
  int i, k;
  int info, nnz;

  int *    ja = new int[colExtent];
  double * a  = new double[colExtent];

  nnz = 0;
  for ( i = 0; i < rowExtent; i++ ) {
    M.fromGetSpRow( srcRow + i, srcCol, a, colExtent, ja,
		     nnz, colExtent, info );
    for( k = 0; k < nnz; k++ ) {
      ja[k] += (destCol - srcCol);
    }
    mStorage->atPutSpRow( destRow + i, a, nnz, ja, info );
    assert(m_Mt == NULL);
  }

  delete [] ja;
  delete [] a;
}


void SparseGenMatrix::mult ( double beta,  OoqpVector& y_in,
				 double alpha, OoqpVector& x_in )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert( x.n == mStorage->n && y.n == mStorage->m );

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->mult( beta, yv, 1, alpha, xv, 1 );
}

void SparseGenMatrix::mult ( double beta,  double y[], int incy,
			     double alpha, double x[], int incx )
{
  mStorage->mult( beta, y, incy, alpha, x, incx);
}


void SparseGenMatrix::transMult ( double beta,   OoqpVector& y_in,
				  double alpha,  OoqpVector& x_in )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);

  assert( x.n == mStorage->m && y.n == mStorage->n );

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->transMult( beta, yv, 1, alpha, xv, 1 );
}

// wrapper added by cpetra
void SparseGenMatrix::transMult( double beta,  OoqpVector& y_in, int incy,
			         double alpha, OoqpVector& x_in, int incx )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert(x.n>0 && y.n>0);
  assert(x.n>=incx*mStorage->m);
  assert(y.n>=incy*mStorage->n);

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->transMult( beta, yv, incy, alpha, xv, incx );
}

// wrapper added by cpetra
void SparseGenMatrix::transMult( double beta,  double yv[], int incy,
				 double alpha, double xv[], int incx )
{
  mStorage->transMult( beta, yv, incy, alpha, xv, incx );
}


double SparseGenMatrix::abmaxnorm()
{
  return mStorage->abmaxnorm();
}


void SparseGenMatrix::atPutDiagonal( int idiag, OoqpVector& vvec )
{
  SimpleVector & v = dynamic_cast<SimpleVector &>(vvec);

  mStorage->atPutDiagonal( idiag, &v[0], 1, v.n );

  assert(m_Mt == NULL);
}


void SparseGenMatrix::fromGetDiagonal( int idiag, OoqpVector& vvec )
{
  mStorage->fromGetDiagonal( idiag, vvec );
}

void SparseGenMatrix::ColumnScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );

  if( m_Mt != NULL )
     m_Mt->RowScale(vec);
}

void SparseGenMatrix::SymmetricScale( OoqpVector& vec )
{
  mStorage->SymmetricScale( vec );

  if( m_Mt != NULL )
     m_Mt->SymmetricScale(vec);
}

void SparseGenMatrix::RowScale( OoqpVector& vec )
{
  mStorage->RowScale( vec );

  if( m_Mt != NULL )
     m_Mt->ColumnScale(vec);
}

void SparseGenMatrix::scalarMult( double num )
{
  mStorage->scalarMult( num );

  if( m_Mt != NULL )
     m_Mt->scalarMult(num);
}

void SparseGenMatrix::matTransDMultMat(OoqpVector& d_, SymMatrix** res)
{
  SimpleVector& d = dynamic_cast<SimpleVector &>(d_);

  int m=mStorage->m; int n=mStorage->n; int nnz=mStorage->numberOfNonZeros();

  if(*res==NULL) {
    assert(m_Mt==NULL);
    //we need to form the transpose
    m_Mt=new SparseGenMatrix(n,m,nnz);
    mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());

    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=NULL; int* jcolMtM=NULL; double* dMtM=NULL;
    mStorage->matTransDSymbMultMat(&d[0], 
				   m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
				   &krowMtM, &jcolMtM, &dMtM);

    *res = new SparseSymMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
  }

  assert(res); 
  assert(m_Mt);

  SparseSymMatrix* MtDM = dynamic_cast<SparseSymMatrix*>(*res);

  mStorage->matTransDMultMat(&d[0], 
			     m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
			     MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}

void SparseGenMatrix::matTransDinvMultMat(OoqpVector& d_, SymMatrix** res)
{
  SimpleVector& d = dynamic_cast<SimpleVector &>(d_);
  int m=mStorage->m; int n=mStorage->n; int nnz=mStorage->numberOfNonZeros();

  if(*res==NULL) {

    // we need to form the transpose
    if(!m_Mt) {
      m_Mt=new SparseGenMatrix(n,m,nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
    }

    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=NULL; int* jcolMtM=NULL; double* dMtM=NULL;
    mStorage->matTransDSymbMultMat(&d[0], 
				   m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
				   &krowMtM, &jcolMtM, &dMtM);

    *res = new SparseSymMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
  }

  assert(res); 
  assert(m_Mt);

  SparseSymMatrix* MtDM = dynamic_cast<SparseSymMatrix*>(*res);

  mStorage->matTransDinvMultMat(&d[0], 
			     m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
			     MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}
void SparseGenMatrix::matMultTrans(SymMatrix** res)
{
  int m=mStorage->m; int n=mStorage->n; 
  int nnz=mStorage->numberOfNonZeros();

  SimpleVector d(n); d.setToConstant(1.0);

  if(*res==NULL) {
    //we need to form the transpose
    if(!m_Mt) {
      m_Mt = new SparseGenMatrix(n,m,nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
    }
    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=NULL; int* jcolMtM=NULL; double* dMtM=NULL;

    m_Mt->mStorage->matTransDSymbMultMat(&d[0], 
					 krowM(), jcolM(), M(),
					 &krowMtM, &jcolMtM, &dMtM);
    *res = new SparseSymMatrix(m, krowMtM[m], krowMtM, jcolMtM, dMtM, 1);
  }
  SparseSymMatrix* MMt = dynamic_cast<SparseSymMatrix*>(*res);
  m_Mt->mStorage->matTransDMultMat(&d[0],
				   krowM(), jcolM(), M(),
				   MMt->krowM(), MMt->jcolM(), MMt->M());
}

void
SparseGenMatrix::getMinMaxVec( bool getMin, bool initializeVec,
      const SparseStorage* storage, const OoqpVector* coScaleVec, OoqpVector& minmaxVec )
{
   SimpleVector& mvec = dynamic_cast<SimpleVector&>(minmaxVec);

   assert(mvec.length() == storage->m);

   if( initializeVec )
   {
      if( getMin )
         mvec.setToConstant(std::numeric_limits<double>::max());
      else
         mvec.setToConstant(0.0);
   }

   if( coScaleVec )
   {
      const SimpleVector* covec = dynamic_cast<const SimpleVector*>(coScaleVec);

      storage->getRowMinMaxVec(getMin, covec->elements(), mvec.elements());
   }
   else
   {
      storage->getRowMinMaxVec(getMin, NULL, mvec.elements());
   }
}

void SparseGenMatrix::getRowMinMaxVec(bool getMin, bool initializeVec,
      const OoqpVector* colScaleVec, OoqpVector& minmaxVec)
{
   getMinMaxVec(getMin, initializeVec, mStorage, colScaleVec, minmaxVec);
}

void SparseGenMatrix::getColMinMaxVec(bool getMin, bool initializeVec,
      const OoqpVector* rowScaleVec, OoqpVector& minmaxVec)
{
   // we may need to form the transpose
   if( !m_Mt )
   {
      const int nnz = mStorage->numberOfNonZeros();
      const int m = mStorage->m;
      const int n = mStorage->n;
      m_Mt = new SparseGenMatrix(n, m, nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
   }

   getMinMaxVec(getMin, initializeVec, m_Mt->mStorage, rowScaleVec, minmaxVec);
}

void SparseGenMatrix::updateTransposed()
{
   if( m_Mt )
   {
      delete m_Mt;
      m_Mt = NULL;
   }
   const int nnz = mStorage->numberOfNonZeros();

   const int m = mStorage->m;
   const int n = mStorage->n;
   m_Mt = new SparseGenMatrix(n, m, nnz);

   mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
}

void SparseGenMatrix::deleteTransposed()
{
   if( m_Mt )
   {
      delete m_Mt;
      m_Mt = NULL;
   }
}


void SparseGenMatrix::updateNonEmptyRowsCount(std::vector<int>& rowcount) const
{
   const int m = mStorage->m;
   const int* const rowStart = mStorage->krowM;

   assert(unsigned(m) == rowcount.size());

   for( int i = 0; i < m; i++ )
      if( rowStart[i] != rowStart[i + 1] )
         rowcount[i]++;
}

void SparseGenMatrix::updateNonEmptyRowsCount(int blockPosition, std::vector<int>& rowcount, std::vector<int>& linkBlockPos1,
      std::vector<int>& linkBlockPos2) const
{
   const int m = mStorage->m;
   const int* const rowStart = mStorage->krowM;

   assert(blockPosition >= 0 && m >= 0);
   assert(unsigned(m) == rowcount.size());
   assert(rowcount.size() == linkBlockPos1.size() && rowcount.size() == linkBlockPos2.size());

   for( int i = 0; i < m; i++ )
      if( rowStart[i] != rowStart[i + 1] )
      {
         if( linkBlockPos1[i] < 0 )
            linkBlockPos1[i] = blockPosition;
         else
            linkBlockPos2[i] = blockPosition;

         rowcount[i]++;
      }
}

SparseGenMatrix& SparseGenMatrix::getTranspose()
{
   if( m_Mt )
     return *m_Mt;

   updateTransposed();

   assert(m_Mt);

   return *m_Mt;
}

void SparseGenMatrix::permuteRows(const std::vector<unsigned int>& permvec)
{
   mStorage->permuteRows(permvec);

   // todo implement column permutation
   assert(!m_Mt);
}


