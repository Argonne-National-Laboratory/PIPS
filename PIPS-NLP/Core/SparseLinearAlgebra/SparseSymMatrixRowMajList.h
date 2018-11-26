/* PIPS-NLP                                                              *
 * Authors: Cosmin G. Petra                                              *
 */

#ifndef SPARSESYMMATRIXRML_H
#define SPARSESYMMATRIXRML_H

#include "DoubleMatrix.h"

#include "OoqpVector.h"
#include <vector>
#include <list>
/** Represents sparse symmetric matrices stored in
 *  row-major Harwell-Boeing format, but using lists 
 *  instead of arrays for fast insertion.
 */

struct ColVal
{
  ColVal(const int& jcol_, const double& M_)
  : jcol(jcol_), M(M_) {};
  int jcol;
  double M;
};

class SparseSymMatrixRowMajList : public SymMatrix {
 public:
  SparseSymMatrixRowMajList( int size );
  //SparseSymMatrixRowMajList( int size, int nnz,
  //		   int krowM[], int jcolM[], double M[], int deleteElts=0);

  virtual int isKindOf( int type ) {assert(false && "deprecated"); return 0; }

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
  virtual void atAddDiagonal( int idiag, OoqpVector& v );

  virtual void fromGetDiagonal( int idiag, OoqpVector& v );

  /** The actual number of structural non-zero elements in this sparse
   *  matrix. This includes so-called "accidental" zeros, elements that
   *  are treated as non-zero even though their value happens to be zero.
   */  
  int numberOfNonZeros() { return nnz;}

  /** Reduce the matrix to lower triangular */
  void reduceToLower();

  virtual void copyMtxFromDouble(int copyLength,double *values);

  virtual ~SparseSymMatrixRowMajList() {};

  virtual void setAdditiveDiagonal(OoqpVector& diag);

  virtual void copyDiagonalVal_From( int idiag, OoqpVector& v, bool firstCall, std::map<int,int> &ValIdxMap );

  virtual void symAtSetSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent);  

  virtual void symAtPutSpRow_CorrectMap( int col, double A[], int lenA, int jcolA[],
				int& info, std::map<int,int> &ValIdxMap,int const constIDX);

  virtual void fromGetSpRow_WithRowStart( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info, int & rowStart);   

  virtual void printMatrixInMatlab(char *name);

  void setToConstant(const double& val);

  virtual void atAddSpRow(const int& row, std::list<ColVal>& colvalSrc);
  virtual void atAddSpRow(const int& row, int* jcolSrc, double* M, const int& nelems);

  void symAtAddSubmatrix( int destRow, int destCol,
			  DoubleMatrix& Mat,
			  int srcRow, int srcCol,
			  int rowExtent, int colExtent );

  void atGetSparseTriplet(int* i, int* j, double* M, bool fortran=true);

  //copies entries to 'M' while checking that the sparsity pattern of 'this' is identical to pattern of
  // the the sparse triplet matrix (i,j,M). Returns false if the pattern does not match
  bool fromGetSparseTriplet_w_patternMatch(const int* i, const int* j, const int& nnz_, double* M); 

  //copies (i,j,M) to 'this'; returns false if an entry of (i,j,M) not found in 'this', otherwise true.
  bool atPutSparseTriplet(const int*i, const int* j, const double* M, const int& nnz_);

  inline void forceSymUpdate(bool val) { bSymUpdate=val; }
 protected:
  std::vector<std::list<ColVal> > vlmat;
  int nnz;
  bool bSymUpdate;
 protected:
  /* Precondition: both colDest and jcolSrc are supposed to be ordered */
  void mergeSetColumn(const int& rowDestIdx, const int& colDestIdxOffset,
		     int* jcolSrc, double* Msrc, const int& nElemsSrc, const int& colSrcIdxOffset,
		     const int& colExtent);
  //do not use this sparingly...
  void putElem(const int& row, const int& col, const double& M);
  void addElem(const int& row, const int& col, const double& M);
};
#endif
