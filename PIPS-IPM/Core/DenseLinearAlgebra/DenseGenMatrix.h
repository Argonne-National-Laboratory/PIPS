/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSEGENMATRIX_H
#define DENSEGENMATRIX_H

#include "DoubleMatrix.h"
#include "DenseStorage.h"
#include "DenseGenMatrixHandle.h"

class DoubleLinearSolver;

/** A class of dense, non-symmetric, possibly non-square, matrices.
 *  @ingroup DenseLinearAlgebra
 */
class DenseGenMatrix : public GenMatrix {
public:
  DenseStorageHandle mStorage;

  DenseGenMatrix( int size );
  DenseGenMatrix( int m, int n );
  DenseGenMatrix( double A[], int m, int n );

  virtual int isKindOf( int matType ) const;

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );

  /** Fill a region of this matrix with zeros.
   *
   *  The region starts at (row, col) and extends rowExtent rows
   *  and colExtent columns.
   */
  virtual void atPutZeros( int row, int col,
			   int rowExtent, int colExtent );


  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );

  virtual void atPutSubmatrix( int destRow, int destCol,
			       DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent );
  virtual void atPutSpRow( int row, double A[], int lenA, int jcolA[],
			   int& info );

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[],
				int& info );

  void mult ( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const override;
  virtual void mult ( double beta,  double y[], int incy,
		      double alpha, const double x[], int incx ) const;

  void transMult ( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const override;
  virtual void transMult ( double beta,  double y[], int incy,
			   double alpha, const double x[], int incx ) const;

  virtual void matTransDMultMat(OoqpVector& d, SymMatrix** res);
  virtual void matTransDinvMultMat(OoqpVector& d, SymMatrix** res);

  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );

  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int rowExtent, int& info );

  virtual void ColumnScale( OoqpVector& vec );
  virtual void RowScale( OoqpVector& vec );
  virtual void SymmetricScale( OoqpVector &vec);
  virtual void scalarMult( double num);

  double abmaxnorm() const override;
  void writeToStream( std::ostream& out ) const override;
  void writeToStreamDense( std::ostream& out ) const override;
  virtual void randomize( double alpha, double beta, double * seed );

  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& v );
  /** Get a row of this matrix. */
  virtual void getRow ( int rowIndex, OoqpVector& v_in);

  double * operator[]( int index ) { return mStorage->M[index]; }

  const double * operator[]( int index ) const
  { return mStorage->M[index]; }

  /** Return a pointer to the first element in the matrix */
  double * elements() { return mStorage->M[0]; };
  /** Return mMat, an    */
  double **Mat() { return mStorage->M; };

  DenseStorage& getStorageRef() { return *mStorage; }
  DenseStorageHandle getStorageHandle() { return mStorage; }

  /* the following functions added by C.Petra 09/09 */

  /** MatMat product
   *
   * this = alpha* op(A) * op(B) + beta*this
   *
   * op(...) specifies whether to use the matrix or its transpose
   */
  virtual void matMult(double alpha,
		       DenseGenMatrix& A, int transA,
		       DenseGenMatrix& B, int transB,
		       double beta);

};

#endif
