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
#include "pipsport.h"

int SparseGenMatrix::isKindOf( int type ) const
{
  return type == kSparseGenMatrix || type == kGenMatrix;
}


SparseGenMatrix::SparseGenMatrix( )
  : mStorageDynamic(nullptr), m_Mt(nullptr)
{
   SpNil(mStorage);
}

SparseGenMatrix::SparseGenMatrix( int rows, int cols, int nnz )
  : mStorageDynamic(nullptr), m_Mt(nullptr)
{
  mStorage = SparseStorageHandle( new SparseStorage( rows, cols, nnz ) );
}


SparseGenMatrix::SparseGenMatrix( int rows, int cols, int nnz,
				  int krowM[], int jcolM[],  double M[],
				  int deleteElts)
  : mStorageDynamic(nullptr), m_Mt(nullptr)
{
  //cout << "SparseGenMatrix1  " << rows << " " << cols << " " << nnz << endl;
  mStorage = SparseStorageHandle( new SparseStorage( rows, cols,nnz,
						     krowM, jcolM, M,
						     deleteElts) );
}

/*SparseGenMatrix::SparseGenMatrix(const std::vector<SparseGenMatrix*> &blocks, bool diagonal)
  : m_Mt(nullptr)
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

  delete mStorageDynamic;
}

/* create a matrix with the same amount of columns but no rows in it */
SparseGenMatrix* SparseGenMatrix::cloneEmptyRows(bool switchToDynamicStorage) const
{
  SparseGenMatrix* clone;

  if( switchToDynamicStorage )
  {
     clone = new SparseGenMatrix();
     clone->mStorageDynamic = new SparseStorageDynamic(0, mStorage->n, 0);
     assert(!m_Mt);
  }
  else
     clone = new SparseGenMatrix(0, mStorage->n, 0);

  if( m_Mt )
  {
     assert(clone->m_Mt == NULL);
     clone->m_Mt = new SparseGenMatrix(0, m_Mt->getStorageRef().n, 0);
  }

  return clone;
}

/* same as clone empty rows but transposes the matrix first */
SparseGenMatrix* SparseGenMatrix::cloneEmptyRowsTransposed(bool switchToDynamicStorage) const
{
  SparseGenMatrix* clone;
  if( switchToDynamicStorage )
  {
     clone = new SparseGenMatrix();
     clone->mStorageDynamic = new SparseStorageDynamic(0, mStorage->m, 0);
     assert(!m_Mt);
  }
  else
     clone = new SparseGenMatrix(0, mStorage->m, 0);

  if( m_Mt )
  {
     assert(clone->m_Mt == nullptr);
     clone->m_Mt = new SparseGenMatrix(0, m_Mt->getStorageRef().m, 0);
  }

  return clone;
}

SparseGenMatrix* SparseGenMatrix::cloneFull(bool switchToDynamicStorage) const
{
   SparseGenMatrix* clone;

   if( switchToDynamicStorage )
   {
      clone = new SparseGenMatrix();
      clone->mStorageDynamic = new SparseStorageDynamic(*mStorage);
      assert(!m_Mt);
   }
   else
   {
      clone = new SparseGenMatrix(mStorage->m, mStorage->n, mStorage->len);
      mStorage->copyFrom(clone->krowM(), clone->jcolM(), clone->M());
   }

   if( m_Mt )
   {
      assert(clone->m_Mt == nullptr);

      SparseStorage& storage_t = m_Mt->getStorageRef();
      clone->m_Mt = new SparseGenMatrix(storage_t.m, storage_t.n, storage_t.len);

      SparseGenMatrix* clone_t = clone->m_Mt;

      storage_t.copyFrom(clone_t->krowM(), clone_t->jcolM(), clone_t->M());
   }

   return clone;
}


void SparseGenMatrix::atPutDense( int row, int col, double * A, int lda,
				      int rowExtent, int colExtent )
{
  mStorage->atPutDense( row, col, A, lda, rowExtent, colExtent );
  assert(m_Mt == nullptr);
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
  assert(m_Mt == nullptr);
}


void SparseGenMatrix::writeToStream(ostream& out) const
{
  mStorage->writeToStream( out );
}

void SparseGenMatrix::writeToStreamDense(ostream& out) const
{
   if( mStorageDynamic != nullptr )
      mStorageDynamic->writeToStreamDense( out );

   else
      mStorage->writeToStreamDense( out );

}

void
SparseGenMatrix::writeToStreamDenseRow(stringstream& out, int rowidx) const
{
   if( mStorageDynamic != nullptr )
   {
      if( mStorageDynamic->getN() > 0 )
      {
         assert(rowidx < mStorageDynamic->getM());
         mStorageDynamic->writeToStreamDenseRow(out, rowidx);
      }
   }
   else
   {
      if( mStorage->n > 0 )
      {
         assert(rowidx < mStorage->m);
         mStorage->writeToStreamDenseRow(out, rowidx);
      }
   }
}

std::string SparseGenMatrix::writeToStreamDenseRow(int rowidx) const
{
   stringstream out;
   if( mStorageDynamic != nullptr )
   {
      if( mStorageDynamic->getN() > 0 )
      {
         assert(rowidx < mStorageDynamic->getM());
         mStorageDynamic->writeToStreamDenseRow(out, rowidx);
      }
   }
   else
   {
      if( mStorage->n > 0 )
      {
         assert(rowidx < mStorage->m);
         mStorage->writeToStreamDenseRow(out, rowidx);
      }

   }
   return out.str();
}

void SparseGenMatrix::randomize( double alpha, double beta, double * seed )
{
  mStorage->randomize( alpha, beta, seed );
  assert(m_Mt == nullptr);
}


void SparseGenMatrix::getDiagonal( OoqpVector& vec )
{
  mStorage->getDiagonal( vec );
}


void SparseGenMatrix::setToDiagonal( OoqpVector& vec )
{
  mStorage->setToDiagonal( vec );
  assert(m_Mt == nullptr);
}


void SparseGenMatrix::atPutSpRow( int row, double A[],
				      int lenA, int jcolA[], int& info )
{
  mStorage->atPutSpRow( row, A, lenA, jcolA, info );
  assert(m_Mt == nullptr);
}


int SparseGenMatrix::numberOfNonZeros()
{
  return mStorage->numberOfNonZeros();
}


void SparseGenMatrix::symmetrize( int& info )
{
  mStorage->symmetrize( info );
  assert(m_Mt == nullptr);
}


void SparseGenMatrix::getSize( long long& m, long long& n ) const
{
   if( mStorageDynamic != nullptr )
   {
      m = mStorageDynamic->getM();
      n = mStorageDynamic->getN();
   }
   else
   {
      m = mStorage->m;
      n = mStorage->n;
   }
}
void SparseGenMatrix::getSize( int& m, int& n ) const
{
  if( mStorageDynamic != nullptr )
  {
     m = mStorageDynamic->getM();
     n = mStorageDynamic->getN();
  }
  else
  {
     m = mStorage->m;
     n = mStorage->n;
  }
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
    assert(m_Mt == nullptr);
  }

  delete [] ja;
  delete [] a;
}


void SparseGenMatrix::mult ( double beta,  OoqpVector& y_in,
				 double alpha, const OoqpVector& x_in ) const
{
  const SimpleVector & x = dynamic_cast<const SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);

  assert( x.n == mStorage->n && y.n == mStorage->m );

  const double *xv = 0;
  double *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->mult( beta, yv, 1, alpha, xv, 1 );
}

void SparseGenMatrix::mult ( double beta,  double y[], int incy,
			     double alpha, double x[], int incx )
{
  mStorage->mult( beta, y, incy, alpha, x, incx);
}

void SparseGenMatrix::multMatSymUpper( double beta, SymMatrix& y,
      double alpha, double x[], int yrowstart, int ycolstart ) const
{
  SparseSymMatrix& y_sparse = dynamic_cast<SparseSymMatrix &>(y);
  assert(!y_sparse.isLower);

  mStorage->multMatSymUpper( beta, y_sparse.getStorageRef(), alpha, x, yrowstart, ycolstart );
}

void SparseGenMatrix::transmultMatSymUpper( double beta, SymMatrix& y,
      double alpha, double x[], int yrowstart, int ycolstart ) const
{
  assert(m_Mt);

  SparseSymMatrix& y_sparse = dynamic_cast<SparseSymMatrix &>(y);
  assert(!y_sparse.isLower);

  m_Mt->getStorageRef().multMatSymUpper( beta, y_sparse.getStorageRef(), alpha, x, yrowstart, ycolstart );
}

void SparseGenMatrix::transMult ( double beta,   OoqpVector& y_in,
				  double alpha,  const OoqpVector& x_in ) const
{
  const SimpleVector & x = dynamic_cast<const SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);

  assert( x.n == mStorage->m && y.n == mStorage->n );

  const double* xv = 0;
  double* yv = 0;

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

double SparseGenMatrix::abmaxnorm() const
{
  return mStorage->abmaxnorm();
}


void SparseGenMatrix::atPutDiagonal( int idiag, OoqpVector& vvec )
{
  SimpleVector & v = dynamic_cast<SimpleVector &>(vvec);

  mStorage->atPutDiagonal( idiag, &v[0], 1, v.n );

  assert(m_Mt == nullptr);
}


void SparseGenMatrix::fromGetDiagonal( int idiag, OoqpVector& vvec )
{
  mStorage->fromGetDiagonal( idiag, vvec );
}

void SparseGenMatrix::ColumnScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );

  if( m_Mt != nullptr )
     m_Mt->RowScale(vec);
}

void SparseGenMatrix::SymmetricScale( OoqpVector& vec )
{
  mStorage->SymmetricScale( vec );

  if( m_Mt != nullptr )
     m_Mt->SymmetricScale(vec);
}

void SparseGenMatrix::RowScale( OoqpVector& vec )
{
  mStorage->RowScale( vec );

  if( m_Mt != nullptr )
     m_Mt->ColumnScale(vec);
}

void SparseGenMatrix::scalarMult( double num )
{
  mStorage->scalarMult( num );

  if( m_Mt != nullptr )
     m_Mt->scalarMult(num);
}

void SparseGenMatrix::matTransDMultMat(OoqpVector& d_, SymMatrix** res)
{
  SimpleVector& d = dynamic_cast<SimpleVector &>(d_);

  int m=mStorage->m; int n=mStorage->n; int nnz=mStorage->numberOfNonZeros();

  if(*res==nullptr) {
    assert(m_Mt==nullptr);
    //we need to form the transpose
    m_Mt=new SparseGenMatrix(n,m,nnz);
    mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());

    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=nullptr; int* jcolMtM=nullptr; double* dMtM=nullptr;
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

void SparseGenMatrix::initTransposed(bool dynamic)
{
   assert(m_Mt == nullptr);

   if( dynamic )
   {
      assert(mStorageDynamic != nullptr);
      m_Mt = new SparseGenMatrix();

      m_Mt->mStorageDynamic = mStorageDynamic->getTranspose();
   }
   else
   {
      const int m = mStorage->m;
      const int n = mStorage->n;
      const int nnz = mStorage->numberOfNonZeros();

      m_Mt = new SparseGenMatrix(n, m, nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
   }
}

void SparseGenMatrix::matTransDinvMultMat(OoqpVector& d_, SymMatrix** res)
{
  SimpleVector& d = dynamic_cast<SimpleVector &>(d_);
  int m=mStorage->m; int n=mStorage->n; int nnz=mStorage->numberOfNonZeros();

  if(*res==nullptr) {

    // we need to form the transpose
    if(!m_Mt) {
      m_Mt=new SparseGenMatrix(n,m,nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
    }

    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=nullptr; int* jcolMtM=nullptr; double* dMtM=nullptr;
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

  if(*res==nullptr) {
    //we need to form the transpose
    if(!m_Mt) {
      m_Mt = new SparseGenMatrix(n,m,nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
    }
    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=nullptr; int* jcolMtM=nullptr; double* dMtM=nullptr;

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
SparseGenMatrix::addNnzPerRow(OoqpVectorBase<int>& nnzVec)
{
   SimpleVectorBase<int>& vec = dynamic_cast<SimpleVectorBase<int>&>(nnzVec);

   if( mStorageDynamic != nullptr )
   {
      assert(vec.length() == mStorageDynamic->getM());
      mStorageDynamic->addNnzPerRow(vec.elements());
   }
   else
   {
      assert(vec.length() == mStorage->m);
      mStorage->addNnzPerRow(vec.elements());
   }
}

void
SparseGenMatrix::addNnzPerCol(OoqpVectorBase<int>& nnzVec)
{
   SimpleVectorBase<int>& vec = dynamic_cast<SimpleVectorBase<int>&>(nnzVec);

   if( !m_Mt)
      initTransposed();

   if( m_Mt->mStorageDynamic != nullptr  )
   {
      assert(vec.length() == m_Mt->mStorageDynamic->getM());
      m_Mt->mStorageDynamic->addNnzPerRow(vec.elements());
   }
   else
   {
      assert(vec.length() == m_Mt->mStorage->m);
      m_Mt->mStorage->addNnzPerRow(vec.elements());
   }
}

void
SparseGenMatrix::addRowSums(OoqpVector& sumVec)
{
   SimpleVector& vec = dynamic_cast<SimpleVector&>(sumVec);

   assert(mStorageDynamic == nullptr && mStorage != nullptr);
   assert(vec.length() == mStorage->m);

   mStorage->addRowSums(vec.elements());
}

void
SparseGenMatrix::addColSums(OoqpVector& sumVec)
{
   SimpleVector& vec = dynamic_cast<SimpleVector&>(sumVec);

   if( !m_Mt )
      initTransposed();

   assert(m_Mt->mStorageDynamic == nullptr && m_Mt->mStorage != nullptr);
   assert(vec.length() == m_Mt->mStorage->m);

   m_Mt->mStorage->addRowSums(vec.elements());
}

void
SparseGenMatrix::getMinMaxVec( bool getMin, bool initializeVec,
      const SparseStorageDynamic* storage_dynamic, const OoqpVector* coScaleVec, OoqpVector& minmaxVec )
{
   SimpleVector& mvec = dynamic_cast<SimpleVector&>(minmaxVec);

   assert(mvec.length() == storage_dynamic->getM());

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

      storage_dynamic->getRowMinMaxVec(getMin, covec->elements(), mvec.elements());
   }
   else
   {
      storage_dynamic->getRowMinMaxVec(getMin, nullptr, mvec.elements());
   }
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
      storage->getRowMinMaxVec(getMin, nullptr, mvec.elements());
   }
}

void SparseGenMatrix::getRowMinMaxVec(bool getMin, bool initializeVec,
      const OoqpVector* colScaleVec, OoqpVector& minmaxVec)
{
   if( hasDynamicStorage() )
      getMinMaxVec(getMin, initializeVec, mStorageDynamic, colScaleVec, minmaxVec);
   else
      getMinMaxVec(getMin, initializeVec, mStorage, colScaleVec, minmaxVec);
}

void SparseGenMatrix::getColMinMaxVec(bool getMin, bool initializeVec,
      const OoqpVector* rowScaleVec, OoqpVector& minmaxVec)
{
   if( hasDynamicStorage() )
   {
      assert( getStorageDynamic() );
      assert( getStorageDynamicTransposed() );

      getMinMaxVec(getMin, initializeVec, getStorageDynamicTransposed(), rowScaleVec, minmaxVec);
   }
   else
   {
      // we may need to form the transpose
      if( !m_Mt )
         initTransposed();

      getMinMaxVec(getMin, initializeVec, m_Mt->mStorage, rowScaleVec, minmaxVec);
   }
}

void SparseGenMatrix::initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>* colNnzVec)
{
   assert(mStorageDynamic != nullptr);

   const SimpleVectorBase<int>& rowNnzVecSimple = dynamic_cast<const SimpleVectorBase<int>&>(rowNnzVec);
   const SimpleVectorBase<int>* colNnzVecSimple = dynamic_cast<const SimpleVectorBase<int>*>(colNnzVec);

   mStorageDynamic->restoreOrder();
   SparseStorageHandle staticStorage(mStorageDynamic->getStaticStorage(rowNnzVecSimple.elements(),
         (colNnzVecSimple == nullptr) ? nullptr : colNnzVecSimple->elements()));

   mStorage = staticStorage;

   assert(mStorage->refs() == 2);
}

void SparseGenMatrix::deleteEmptyRowsCols(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec)
{
   const SimpleVectorBase<int>& rowNnzVecSimple = dynamic_cast<const SimpleVectorBase<int>&>(rowNnzVec);
   const SimpleVectorBase<int>& colNnzVecSimple = dynamic_cast<const SimpleVectorBase<int>&>(colNnzVec);

   mStorage->deleteEmptyRowsCols(rowNnzVecSimple.elements(), colNnzVecSimple.elements());
}

void SparseGenMatrix::fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset,
       double* rowsArrayDense, int* rowSparsity)
{

   mStorage->fromGetRowsBlock(rowIndices, nRows, arrayLineSize, arrayLineOffset, rowsArrayDense, rowSparsity);
}

void SparseGenMatrix::deleteEmptyRows(int*& orgIndex)
{
   mStorage->deleteEmptyRows(orgIndex);

   if( m_Mt )
      this->updateTransposed();
}

void SparseGenMatrix::fromGetColsBlock(const int* colIndices, int nCols, int arrayLineSize, int arrayLineOffset,
       double* colsArrayDense, int* rowSparsity)
{
   if( !m_Mt )
      updateTransposed();

   m_Mt->getStorageRef().fromGetRowsBlock(colIndices, nCols, arrayLineSize, arrayLineOffset,
         colsArrayDense, rowSparsity);
}

bool SparseGenMatrix::hasTransposed() const
{
   return (m_Mt != nullptr);
}

void SparseGenMatrix::freeDynamicStorage()
{
   delete mStorageDynamic;
   mStorageDynamic = nullptr;
}

void SparseGenMatrix::updateTransposed()
{
   if( m_Mt )
   {
      delete m_Mt;
      m_Mt = nullptr;
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
      m_Mt = nullptr;
   }
}

void SparseGenMatrix::getLinkVarsNnz(std::vector<int>& vec) const
{
   mStorage->getLinkVarsNnz(vec);
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

   if( m_Mt )
      m_Mt->mStorage->permuteCols(permvec);
}

void SparseGenMatrix::permuteCols(const std::vector<unsigned int>& permvec)
{
   mStorage->permuteCols(permvec);

   if( m_Mt )
      m_Mt->mStorage->permuteRows(permvec);
}

int SparseGenMatrix::appendRow(const SparseGenMatrix& matrix_row, int row)
{
   assert( hasDynamicStorage() );
   assert( !hasTransposed() );
   assert( matrix_row.hasDynamicStorage() );

   mStorageDynamic->appendRow( matrix_row.getStorageDynamicRef(), row );

   return mStorageDynamic->getM() - 1;
}

int SparseGenMatrix::appendCol(const SparseGenMatrix& matrix_col, int col)
{
   assert( matrix_col.hasTransposed() );
   assert( matrix_col.hasDynamicStorage() );
   assert( hasDynamicStorage() );
   assert( !hasTransposed() );

   mStorageDynamic->appendRow( matrix_col.getStorageDynamicTransposedRef(), col );

   return mStorageDynamic->getM() - 1;
}

void SparseGenMatrix::axpyWithRowAt(double alpha, SimpleVector& y, int row) const
{
   assert( hasDynamicStorage() );

   mStorageDynamic->axpyWithRowAt(alpha, y.elements(), y.length(), row);
}

void SparseGenMatrix::axpyWithRowAtPosNeg( double alpha, SimpleVector& y_pos, SimpleVector& y_neg, int row) const
{
   assert( hasDynamicStorage() );
   assert( y_pos.length() == y_neg.length() );

   mStorageDynamic->axpyWithRowAtPosNeg(alpha, y_pos.elements(), y_neg.elements(), y_pos.length(), row);
}

double SparseGenMatrix::localRowTimesVec( const SimpleVector& vec, int row ) const
{
  assert( hasDynamicStorage() );

  return mStorageDynamic->rowTimesVec( vec.elements(), vec.length(), row);
}

void SparseGenMatrix::removeRow( int row )
{
   assert( hasDynamicStorage() );
   if( hasTransposed() )
   {
      assert( m_Mt->getStorageDynamic() );
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );
      removeRowUsingTransposed( row, m_Mt->getStorageDynamicRef() );
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );
   }
   else
      mStorageDynamic->clearRow(row);
}

void SparseGenMatrix::removeRowUsingTransposed( int row, SparseStorageDynamic& mat_trans)
{
   assert( hasDynamicStorage() );

   SparseStorageDynamic& mat = *mStorageDynamic;

   const int start = mat.getRowPtr(row).start;
   const int end = mat.getRowPtr(row).end;

   for( int i = start; i < end; ++i )
   {
      const int col = mat.getJcolM(i);
      mat_trans.removeEntryAtRowCol(col, row);
   }

   mat.clearRow(row);
}

void SparseGenMatrix::removeCol( int col )
{
   assert( hasDynamicStorage() );
   if( hasTransposed() )
   {
      assert(m_Mt->getStorageDynamic());
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );
      m_Mt->removeRowUsingTransposed( col, *mStorageDynamic);

      m_Mt->removeRow( col );
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );
   }
   else
      mStorageDynamic->clearCol(col);
}

void SparseGenMatrix::removeEntryAtRowCol( int row, int col )
{
   assert( hasDynamicStorage() );

   if( hasTransposed() )
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );

   mStorageDynamic->removeEntryAtRowCol(row, col);

   if( hasTransposed() )
   {
      m_Mt->removeEntryAtRowCol(col, row);
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );
   }
}

void SparseGenMatrix::removeEntryAtRowColIndex( int row, int col_index )
{
   assert( hasDynamicStorage() );

   if( hasTransposed() )
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );

   const int col = mStorageDynamic->getJcolM(col_index);
   mStorageDynamic->removeEntryAtIndex(row, col_index);

   if( hasTransposed() )
   {
      m_Mt->removeEntryAtRowCol(col, row);
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );
   }
}

void SparseGenMatrix::addColToRow( double coeff, int col, int row )
{
   assert( hasDynamicStorage() );

   if( hasTransposed() )
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );

   mStorageDynamic->addColToRow(coeff, col, row);

   if( hasTransposed() )
   {
      m_Mt->addColToRow(coeff, row, col);
      assert( mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals() );
   }
}
