#ifndef STOCHSYMMATRIX_H
#define STOCHSYMMATRIX_H

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "pipsport.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "mpi.h"

class StochSymMatrix : public SymMatrix {

private:

  // note: also used for dummy class!
  virtual void deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec, const OoqpVectorBase<int>* linkParent);
  virtual void writeToStreamDenseChild(stringstream& out, int offset) const;

public:
  /** Constructs a matrix with local size 'local_n' having 'local_nnz' local nonzeros
      and set the global size and the id to to 'global_n' and 'id', respectively.
      The parameter 'id' is used for output/debug purposes only.
      The created matrix will have no children.*/
  StochSymMatrix( int id, long long global_n, int local_n, int local_nnz, MPI_Comm mpiComm );
  StochSymMatrix( int id, long long global_n, 
		  int diag_n, int diag_nnz, 
		  int border_n, int border_nnz,
		  MPI_Comm mpiComm_);
  //StochSymMatrix(const vector<StochSymMatrix*> &blocks); -- not needed anymore; petra
  virtual ~StochSymMatrix();

  std::vector<StochSymMatrix*> children;
  SparseSymMatrix* diag;
  SparseGenMatrix* border;
  int id;
  long long n;
  MPI_Comm mpiComm;
  int iAmDistrib;
  
  virtual void AddChild(StochSymMatrix* child);

  virtual StochSymMatrix* clone() const;

  virtual int isKindOf( int type ) const;
  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );

  virtual void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info );

  virtual void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info );

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;

  virtual long long size();

  virtual void symAtPutSubmatrix( int destRow, int destCol,
				  DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent );

  virtual void fromGetSpRow( int row, int col,
                             double A[], int lenA, int irowA[], int& nnz,
                             int rowExtent, int& info );

  virtual void atPutZeros( int row, int col,
			   int rowExtent, int colExtent );
  void mult ( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const override;
  void transMult ( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const override;
  
  double abmaxnorm() const override;
  
  void writeToStream(ostream& out) const override;

  void writeToStreamDense(std::ostream& out) const override;

  virtual void randomizePSD(double * seed);
  
  virtual void getDiagonal( OoqpVector& vec );
  virtual void setToDiagonal( OoqpVector& vec );
  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& x );

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info );

  virtual void SymmetricScale ( OoqpVector& vec );
  virtual void ColumnScale ( OoqpVector& vec );
  virtual void RowScale ( OoqpVector& vec );
  virtual void scalarMult( double num );

  // note: also used for dummy class!
  virtual void deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec)
  {
     deleteEmptyRowsCols(nnzVec, nullptr);
  }


 protected:
  StochSymMatrix* parent;
};

/** 
 * Dummy stochastic symmetric matrix
 */

class StochSymDummyMatrix : public StochSymMatrix
{


private:
   void writeToStreamDenseChild(stringstream& out, int offset) const override {};


protected:

public:

  StochSymDummyMatrix(int id_)
    : StochSymMatrix(id_, 0, 0, 0, MPI_COMM_NULL) {};

  virtual ~StochSymDummyMatrix(){};

  StochSymDummyMatrix* clone() const override { return new StochSymDummyMatrix(id); };

  void AddChild(StochSymMatrix* child) override {};

  int isKindOf( int type ) const override;

  void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override {};
  void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override {};

  void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info ) override {};

  void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info ) override {};

  void getSize( long long& m, long long& n ) const override { m = 0; n = 0; }
  void getSize( int& m, int& n ) const override { m = 0; n = 0; }

  long long size() override { return 0; }

  void symAtPutSubmatrix( int destRow, int destCol,
				  DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent ) override {};

  void fromGetSpRow( int row, int col,
                             double A[], int lenA, int irowA[], int& nnz,
                             int rowExtent, int& info ) override {};

  void atPutZeros( int row, int col,
			   int rowExtent, int colExtent ) override {};
  void mult ( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const override {};
  void transMult ( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const override {};
  
  double abmaxnorm() const override { return 0.0; }
  
  void writeToStream(ostream& out) const override {};
  void writeToStreamDense(std::ostream& out) const override {};

  void randomizePSD(double * seed) override {};
  
  void getDiagonal( OoqpVector& vec ) override {};
  void setToDiagonal( OoqpVector& vec ) override {};
  void atPutDiagonal( int idiag, OoqpVector& v ) override {};
  void fromGetDiagonal( int idiag, OoqpVector& x ) override {};

  void putSparseTriple( int irow[], int len, int jcol[], double A[],
				int& info ) override {};

  void SymmetricScale ( OoqpVector& vec ) override {};
  void ColumnScale ( OoqpVector& vec ) override {};
  void RowScale ( OoqpVector& vec ) override {};
  void scalarMult( double num ) override {};
};

typedef SmartPointer<StochSymMatrix> StochSymMatrixHandle;

#endif
