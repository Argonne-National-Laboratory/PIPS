/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSESYMMATRIX_H
#define SPARSESYMMATRIX_H

#include "DoubleMatrix.h"

#include "SparseStorage.h"
#include "OoqpVector.h"
#include "SparseSymMatrixHandle.h"

/** Represents sparse symmetric matrices stored in
 *  row-major Harwell-Boeing format.
 *  @ingroup SparseLinearAlgebra
 */
class SparseSymMatrix : public SymMatrix {
  SparseStorageHandle mStorage;
public:
  SparseSymMatrix();
  SparseSymMatrix( int size, int nnz, bool isLower = true );
  SparseSymMatrix( int size, int nnz,
		   int krowM[], int jcolM[], double M[], int deleteElts = 0, bool isLower = true);
  //SparseSymMatrix(const std::vector<SparseSymMatrix*> &blocks); not needed anymore; cpetra

  SparseStorage&  getStorageRef() { return *mStorage; }
  SparseStorageHandle  getStorageHandle() { return mStorage; }

  // is lower part of matrix stored? (otherwise upper part is stored)
  const bool isLower;

  int * krowM() { return mStorage->krowM; }
  int * jcolM() { return mStorage->jcolM; }
  double * M() { return mStorage->M; }

  virtual int isKindOf( int type ) const;

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[],
				int& info );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info );
  virtual void SymmetricScale ( OoqpVector& vec );
  virtual void ColumnScale ( OoqpVector& vec );
  virtual void RowScale ( OoqpVector& vec );
  virtual void scalarMult( double num);

  virtual void symAtPutSpRow( int col, double A[], int lenA, int jcolA[],
			      int& info );

  virtual void symPutZeroes();

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;
  virtual long long size();

  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );

  virtual void symAtPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent );

  virtual void mult ( double beta,  double y[], int incy,
		      double alpha, const double x[], int incx ) const;
  virtual void transMult ( double beta,  double y[], int incy,
			   double alpha, const double x[], int incx ) const;

  void mult ( double beta,  OoqpVector& y,
                      double alpha, const OoqpVector& x ) const override;

  void transMult ( double beta,   OoqpVector& y,
                           double alpha,  const OoqpVector& x ) const override;

  double abmaxnorm() const override;

  void writeToStream( std::ostream& out ) const override;
  void writeToStreamDense( std::ostream& out ) const override;
  virtual std::string writeToStreamDenseRow(int row) const;
  virtual void writeToStreamDenseRow( std::stringstream& out, int row ) const;

  virtual void randomizePSD(double *);

  virtual void atPutDiagonal( int idiag, OoqpVector& v );

  virtual void fromGetDiagonal( int idiag, OoqpVector& v );

  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */
  int numberOfNonZeros() { return mStorage->numberOfNonZeros(); }

  /** Reduce the matrix to lower triangular */
  void reduceToLower();

  void deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec);

  void deleteZeroRowsCols(int*& new2orgIdx);

  void getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const;

  void getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const;


  virtual ~SparseSymMatrix() {};
};
#endif
