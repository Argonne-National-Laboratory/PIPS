/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

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
	
  SparseSymMatrix( int size, int nnz );
  SparseSymMatrix( int size, int nnz,
		   int krowM[], int jcolM[], double M[], int deleteElts=0);
  //SparseSymMatrix(const std::vector<SparseSymMatrix*> &blocks); not needed anymore; cpetra

  SparseStorage&  getStorageRef() { return *mStorage; }
  SparseStorage*  getStorage() { return mStorage.ptr(); }
  int * krowM() { return mStorage->krowM; }
  int * jcolM() { return mStorage->jcolM; }
  double * M() { return mStorage->M; }

  virtual int isKindOf( int type );

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

  virtual void getSize( long long& m, long long& n );
  virtual void getSize( int& m, int& n );
  virtual long long size();

  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );

  virtual void symAtPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent );

  virtual void mult ( double beta,  double y[], int incy,
		      double alpha, double x[], int incx );
  virtual void transMult ( double beta,  double y[], int incy,
			   double alpha, double x[], int incx );
  
  virtual void mult ( double beta,  OoqpVector& y,
                      double alpha, OoqpVector& x );

  virtual void transMult ( double beta,   OoqpVector& y,
                           double alpha,  OoqpVector& x );

  virtual double abmaxnorm();
  
  virtual void writeToStream(std::ostream& out) const;
  
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

  virtual void copyMtxFromDouble(int copyLength,double *values);

  virtual ~SparseSymMatrix() {};

  virtual void setAdditiveDiagonal(OoqpVector& diag);

  virtual void copyDiagonalVal_From( int idiag, OoqpVector& v, bool firstCall, std::map<int,int> &ValIdxMap );

  virtual void symAtSetSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent,bool firstCall, std::map<int,int> &ValIdxMap );  

  virtual void symAtPutSpRow_CorrectMap( int col, double A[], int lenA, int jcolA[],
				int& info, std::map<int,int> &ValIdxMap,int const constIDX);

  virtual void fromGetSpRow_WithRowStart( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info, int & rowStart);   

  virtual double * getMatVal(){return mStorage->M; }


  virtual void fromGetDense_withMap( int row, int col, double * A,
				       int lda,
				       int rowExtent, int colExtent, int const FirstCall, std::map<int,int> &ValIdxMap ); 

  virtual void printMatrixInMatlab(char *name);

  virtual void symmetrize_valonly( double *val_lower,int *goffIDX);

};
#endif
