/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSEGENMATRIX_H
#define SPARSEGENMATRIX_H

#include "OoqpVectorHandle.h"
#include "DoubleMatrix.h"
#include "SparseStorage.h"
#include "SparseStorageDynamic.h"
#include "SparseGenMatrixHandle.h"
#include <vector>

/** Represents sparse non-symmetric, possibly non-square matrices stored in
 *  row-major Harwell-Boeing format.
 *  @ingroup SparseLinearAlgebra
 */
class SparseGenMatrix : public GenMatrix {
private:
  static
  void getMinMaxVec( bool getMin, bool initializeVec,
        const SparseStorage* storage, const OoqpVector* coScaleVec, OoqpVector& minmaxVec );

protected:
  SparseStorageHandle mStorage;
  SparseStorageDynamic* mStorageDynamic;

  // in the case of A'*A we internally form the transpose only once
  SparseGenMatrix* m_Mt;

public:

  SparseGenMatrix();

  void updateTransposed();
  void deleteTransposed();

  SparseGenMatrix( int rows, int cols, int nnz );
  SparseGenMatrix( int rows, int cols, int nnz,
		   int krowM[], int jcolM[], double M[],
		   int deleteElts=0);

  virtual SparseGenMatrix* cloneEmptyRows(bool switchToDynamicStorage = false) const;
  virtual SparseGenMatrix* cloneFull(bool switchToDynamicStorage = false) const;

  virtual void getSize( long long& m, long long& n );
  virtual void getSize( int& m, int& n );

  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */
  virtual int numberOfNonZeros();

  virtual int isKindOf( int matType ) const;

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );
  virtual void ColumnScale( OoqpVector& vec );
  virtual void RowScale( OoqpVector& vec );
  virtual void SymmetricScale( OoqpVector &vec);
  virtual void scalarMult( double num);
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info );
  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent );
  virtual void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info );

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[],
				int& info );

  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );

  virtual void mult ( double beta,  OoqpVector& y,
                      double alpha, OoqpVector& x );
  virtual void mult ( double beta,  double y[], int incy,
                      double alpha, double x[], int incx );

  virtual void multMatSymUpper( double beta, SymMatrix& y,
        double alpha, double x[], int yrowstart, int ycolstart ) const;

  virtual void transmultMatSymUpper( double beta, SymMatrix& y,
        double alpha, double x[], int yrowstart, int ycolstart ) const;

  virtual void transMult( double beta,   OoqpVector& y,
			  double alpha,  OoqpVector& x );
  virtual void transMult( double beta,  OoqpVector& y_in, int incy,
			  double alpha, OoqpVector& x_in, int incx );
  virtual void transMult( double beta,  double y_in[], int incy,
			  double alpha, double x_in[], int incx );

  /** C = this^T * D * this where D=diag(d) is a diagonal matrix. */
  virtual void matTransDMultMat(OoqpVector& d, SymMatrix** res);
  /** C = this^T * inv(D) * this where D=diag(d) is a diagonal matrix. */
  virtual void matTransDinvMultMat(OoqpVector& d, SymMatrix** res);

  /** initialize (dynamic) transposed matrix */
  virtual void initTransposed(bool dynamic = false);

  /** C = this * this^T */
  virtual void matMultTrans(SymMatrix** res);

  virtual double abmaxnorm();

  virtual void writeToStream(ostream& out) const;
  virtual void writeToStreamDense(ostream& out) const;
  virtual void writeToStreamDenseRow( stringstream& out, int rowidx) const;
  virtual std::string writeToStreamDenseRow( int rowidx) const;

  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  virtual void symmetrize( int& info );

  virtual void randomize( double alpha, double beta, double * seed );

  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& v );

  SparseStorageHandle getStorageHandle() { return mStorage; }
  SparseStorage& getStorageRef() { return *mStorage; }
  int * krowM() { return mStorage->krowM; }
  int * jcolM() { return mStorage->jcolM; }
  double * M() { return mStorage->M; }

  SparseStorageDynamic * getStorageDynamic() { assert(mStorageDynamic != NULL); return mStorageDynamic; }
  const SparseStorageDynamic * getStorageDynamic() const { assert(mStorageDynamic != NULL); return mStorageDynamic; }
  SparseStorageDynamic& getStorageDynamicRef() { assert(mStorageDynamic != NULL); return *mStorageDynamic; }
  const SparseStorageDynamic& getStorageDynamicRef() const { assert(mStorageDynamic != NULL); return *mStorageDynamic; }
  SparseStorageDynamic * getStorageDynamicTransposed() { assert(m_Mt != NULL && m_Mt->hasDynamicStorage() ); return m_Mt->getStorageDynamic(); }
  const SparseStorageDynamic * getStorageDynamicTransposed() const { assert(m_Mt != NULL && m_Mt->hasDynamicStorage() ); return m_Mt->getStorageDynamic(); }
  SparseStorageDynamic& getStorageDynamicTransposedRef() { assert(m_Mt != NULL && m_Mt->hasDynamicStorage()); return m_Mt->getStorageDynamicRef(); }
  const SparseStorageDynamic& getStorageDynamicTransposedRef() const { assert(m_Mt != NULL && m_Mt->hasDynamicStorage()); return m_Mt->getStorageDynamicRef(); }
  bool hasDynamicStorage() { return (mStorageDynamic != NULL); }

  virtual void addNnzPerRow(OoqpVector& nnzVec);
  virtual void addNnzPerCol(OoqpVector& nnzVec);

  /** fill vector with absolute minimum/maximum value of each row */
  virtual void getRowMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* colScaleVec, OoqpVector& minmaxVec );

  /** fill vector with absolute minimum/maximum value of each column */
  virtual void getColMinMaxVec( bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, OoqpVector& minmaxVec );

  virtual void addRowSums( OoqpVector& sumVec );
  virtual void addColSums( OoqpVector& sumVec );

  void initStaticStorageFromDynamic(const OoqpVector& rowNnzVec, const OoqpVector* colNnzVec);

  void permuteRows(const std::vector<unsigned int>& permvec);

  void permuteCols(const std::vector<unsigned int>& permvec);

  void getLinkVarsNnz(std::vector<int>& vec) const;

  void updateNonEmptyRowsCount(std::vector<int>& rowcount) const;

  void updateNonEmptyRowsCount(int blockPosition, std::vector<int>& rowcount, std::vector<int>& linkBlockPos1,
     std::vector<int>& linkBlockPos2) const;

  SparseGenMatrix& getTranspose();

  void deleteEmptyRowsCols(const OoqpVector& rowNnzVec, const OoqpVector& colNnzVec);

  void deleteEmptyRows(int*& orgIndex);

  void fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset,
        double* rowsArrayDense, int* rowSparsity = NULL);

  void fromGetColsBlock(const int* colIndices, int nCols, int arrayLineSize, int arrayLineOffset,
        double* colsArrayDense, int* rowSparsity = NULL);

  bool hasTransposed() const;

  void freeDynamicStorage();

  virtual int appendRow( const OoqpVector& row );

  virtual ~SparseGenMatrix();

};

#endif
