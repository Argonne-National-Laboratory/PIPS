#ifndef STOCHSYMMATRIX_H
#define STOCHSYMMATRIX_H

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "mpi.h"

class StochSymMatrix : public SymMatrix {
protected:

public:
  /** Constructs a matrix with local size 'local_n' having 'local_nnz' local nonzeros
      and set the global size and the id to to 'global_n' and 'id', respectively.
      The parameter 'id' is used for output/debug purposes only.
      The created matrix will have no children.*/
  StochSymMatrix(int id, int global_n, int local_n, int local_nnz, MPI_Comm mpiComm);
  StochSymMatrix(const vector<StochSymMatrix*> &blocks);
  virtual ~StochSymMatrix();

  std::vector<StochSymMatrix*> children;
  SparseSymMatrix* mat;
  int id;
  int n;
  MPI_Comm mpiComm;
  
  virtual void AddChild(StochSymMatrix* child);


  virtual int isKindOf( int type );
  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );

  virtual void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info );

  virtual void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info );

  virtual void getSize( int& m, int& n );

  virtual int size();

  virtual void symAtPutSubmatrix( int destRow, int destCol,
				  DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent );

  virtual void fromGetSpRow( int row, int col,
                             double A[], int lenA, int irowA[], int& nnz,
                             int rowExtent, int& info );

  virtual void atPutZeros( int row, int col,
			   int rowExtent, int colExtent );
  virtual void mult ( double beta,  OoqpVector& y,
		      double alpha, OoqpVector& x );
  virtual void transMult ( double beta,  OoqpVector& y,
			   double alpha, OoqpVector& x );
  
  virtual double abmaxnorm();
  
  virtual void writeToStream(ostream& out) const;

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
  
  
};

/** 
 * Dummy stochastics symmetric matrix
 */

class StochSymDummyMatrix : public StochSymMatrix {
protected:

public:

  StochSymDummyMatrix(int id_)
    : StochSymMatrix(id_, 0, 0, 0, MPI_COMM_NULL) {};

  virtual ~StochSymDummyMatrix(){};

  virtual void AddChild(StochSymMatrix* child){};

  virtual int isKindOf( int type );

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ){};
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ){};

  virtual void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info ){};

  virtual void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info ){};

  virtual void getSize( int& m, int& n ){m=0;n=0;}

  virtual int size(){return 0;}

  virtual void symAtPutSubmatrix( int destRow, int destCol,
				  DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent ){};

  virtual void fromGetSpRow( int row, int col,
                             double A[], int lenA, int irowA[], int& nnz,
                             int rowExtent, int& info ){};

  virtual void atPutZeros( int row, int col,
			   int rowExtent, int colExtent ){};
  virtual void mult ( double beta,  OoqpVector& y,
		      double alpha, OoqpVector& x ){};
  virtual void transMult ( double beta,  OoqpVector& y,
			   double alpha, OoqpVector& x ){};
  
  virtual double abmaxnorm(){return 0.0;}
  
  virtual void writeToStream(ostream& out) const{};

  virtual void randomizePSD(double * seed){};
  
  virtual void getDiagonal( OoqpVector& vec ){};
  virtual void setToDiagonal( OoqpVector& vec ){};
  virtual void atPutDiagonal( int idiag, OoqpVector& v ){};
  virtual void fromGetDiagonal( int idiag, OoqpVector& x ){};

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info ){};

  virtual void SymmetricScale ( OoqpVector& vec ){};
  virtual void ColumnScale ( OoqpVector& vec ){};
  virtual void RowScale ( OoqpVector& vec ){};
  virtual void scalarMult( double num ){};
  
  
};

typedef SmartPointer<StochSymMatrix> StochSymMatrixHandle;

#endif
