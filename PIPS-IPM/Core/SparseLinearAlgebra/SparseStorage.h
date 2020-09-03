/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSESTORAGE_H
#define SPARSESTORAGE_H

#include "DoubleMatrix.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "pipsport.h"

#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

/** A class for managing the matrix elements used by sparse matrices.
 *  @ingroup SparseLinearAlgebra
 */
class SparseStorage : public DoubleStorage {
private:
  /** store absolute non-zero minimum entry of row i and vec[i] in vec[i]; empty rows get value 0.0  */
  void getRowMinVec(const double* colScaleVec, double* vec) const;

  /** store absolute non-zero maximum entry of row i and vec[i] in vec[i]; empty rows get value 0.0  */
  void getRowMaxVec(const double* colScaleVec, double* vec) const;

  class index_sort
  {
     private:
       const int* indices;
       const int maxsize;

     public:
       index_sort(const int* indices, int maxsize) : indices(indices), maxsize(maxsize) {}

       bool operator()(int i, int j) const
       {
          assert(i < maxsize && j < maxsize);
          return (indices[i] < indices[j]);
       }
  };
protected:
  int neverDeleteElts;
  
public:
  static int instances;

  int m;
  int n;
  int len;
  int * jcolM;
  int * krowM;
  double * M;

  SparseStorage( int m_, int n_, int len_ );
  SparseStorage( int m_, int n_, int len_,
		 int * krowM_, int * jcolM_, double * M_,
		 int deleteElts=0);
  //SparseStorage(const vector<SparseStorage*> &blocks, bool diagonal); -- not needed anymore; cpetra

  void copyFrom(int * krowM_, int * jcolM_, double * M_) const;

  void shiftRows( int row, int shift, int& info );
  void getSize( int& m, int& n ) const override;
  int rows() { return m; }
  int cols() { return n; }

  bool isValid(bool verbose = false) const;
  bool isSorted() const;


  int length() { return len; };
  int numberOfNonZeros() const {	return krowM[m]; };
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );
  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info );

  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );

  virtual void ColumnScale( OoqpVector& vec );
  virtual void RowScale( OoqpVector& vec );
  virtual void SymmetricScale( OoqpVector& vec );
  virtual void scalarMult( double num);

  virtual void atPutSpRow( int col, double A[], int lenA, int irowA[],
			   int& info );

  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int irowA[], int& nnz,
			     int rowExtent, int& info );

  virtual void randomize( double alpha, double beta, double * seed );

  virtual void clear();

  virtual void getTransposePat( int row, int col, int rowExtent, int colExtent,
				int kpat[], int krowM[], int jcolM[] );
  virtual void getFromPat( double data[], int n, int kpat[] );
  virtual void mult( double beta,  double y[], int incy,
		     double alpha, const double x[], int incx ) const;

  virtual void transMult ( double beta,  double y[], int incy,
			   double alpha, const double x[], int incx ) const;

  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& v );

  virtual void atPutDiagonal( int idiag,
			      double x[], int incx, int extent );

  virtual void writeToStream(ostream& out) const;
  virtual void writeToStreamDense(ostream& out) const;
  virtual void writeToStreamDenseRow( stringstream& out, int rowidx) const;

  virtual void symmetrize( int& info);
  double abmaxnorm() const override;

  /** Computes the sparsity pattern of MtM = M^T * D * M 
   *  where D=diag(d) is a diagonal matrix and M=this.
   *  
   *  Find the nonzero pattern of the product matrix, allocate it and returns it.
   *
   *  Also allocates, builds and returns this^T since it is needed later for
   *   numerical multiplication.
   */
  void matTransDSymbMultMat(double* d, 
			    int* krowMt, int* jcolMt, double* dMt,
			    int** krowMtM, int** jcolMtM, double** dMtM); 
			    

  /** Numerical multiplication MtM = M^T * D * M  where 
   *  D=diag(d) is a diagonal matrix and M=this.
   *  M^T and MtM buffers should be allocated before calling this method by calling
   *  method matTransDSymbMultMat.
   */
  void matTransDMultMat(double* d, 
			int* krowMt, int* jcolMt, double* dMt,
			int* krowMtM, int* jcolMtM, double* dMtM);
  void matTransDinvMultMat(double* d, 
			int* krowMt, int* jcolMt, double* dMt,
			int* krowMtM, int* jcolMtM, double* dMtM);
  /** Builds the transpose: Mt = this^T */
  void transpose(int* krowMt, int* jcolMt, double* dMt) const;

  void reduceToLower();

  void multMatSymUpper( double beta, SparseStorage& y,
                double alpha, double x[], int yrow, int ycolstart ) const;

  void transMultLower( double beta,  double y[],
			       double alpha, double x[], int firstrow );
  void transMultMat( double beta,  double* Y, int ny, int ldy,
						 double alpha, double *X, int ldx);
  void transMultMatLower( double* Y, int ny, int firstrow,
						 double alpha, double *X, int ldx);

  /** Y <- alpha* M^T X + beta*Y, where M is this 
   * Special update function, computes only the elements in Y that are lower
   * triangular elements in a larger matrix that contains Y (see impl file for 
   * more details)
   */
  void transMultMatLower( double beta,  double* Y, int ny, int ldy,
			  double alpha, double *X, int ldx, int colStart);

  void fromGetColBlock(int col, double *A, int lda, int colExtent, bool &allzero);
  void fromGetColBlock(int col, double *A, int lda, int colExtent, int* colSparsity, bool &allzero);

  void fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset,
        double* rowsArrayDense, int* rowSparsity = nullptr);

  /** add nnz per row to given array (of size nRows) */
  void addNnzPerRow(int* vec) const;
  void getLinkVarsNnz(std::vector<int>& vec) const;

  /** add abs. sum per row to given array (of size nRows) */
  void addRowSums( double* vec ) const;

  /** store absolute non-zero minimum/maximum entry of row i and vec[i] in vec[i];
   *  empty rows get value 0.0 for maximization and <double>::max() for minimization  */
  void getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const;

  void permuteRows(const std::vector<unsigned int>& permvec);
  void permuteCols(const std::vector<unsigned int>& permvec);

  void sortCols();

  void dump(const string& filename);

  void deleteEmptyRowsCols(const int* nnzRowVec, const int* nnzColVec);

  void getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const;

  void getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const;

  void deleteEmptyRows(int*& orgIndex);

  // should be used with care! other methods might nor work correctly todo: add flag to check in other methods
  void c2fortran();

  void fortran2c();

  bool fortranIndexed() const;

  void set2FortranIndexed();

  void deleteZeroRowsColsSym(int*& new2orgIdx);


  /*
   * computes the full sparse matrix representation from a upper triangular symmetric sparse representation
   *
   * Must be square, the storage for the full representation will be allocated within the matrix and must be released later
   */
  void fullMatrixFromUpperTriangular(int*& rowPtrFull, int*& colIdxFull, double*& valuesFull) const;


  virtual ~SparseStorage();

private:
  bool isFortranIndexed;
};

#endif
