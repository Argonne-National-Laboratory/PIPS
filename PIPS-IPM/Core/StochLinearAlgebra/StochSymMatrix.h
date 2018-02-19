#ifndef STOCHSYMMATRIX_H
#define STOCHSYMMATRIX_H

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "mpi.h"

class StochSymMatrix : public SymMatrix {

private:

  virtual void deleteEmptyRowsCols(const OoqpVector& nnzVec, const OoqpVector* linkParent);

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

  virtual int isKindOf( int type );
  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );

  virtual void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info );

  virtual void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info );

  virtual void getSize( long long& m, long long& n );
  virtual void getSize( int& m, int& n );

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

  virtual void deleteEmptyRowsCols(const OoqpVector& nnzVec)
  {
     deleteEmptyRowsCols(nnzVec, NULL);
  }


 protected:
  StochSymMatrix* parent;
};

/** 
 * Dummy stochastic symmetric matrix
 */

class StochSymDummyMatrix : public StochSymMatrix {
protected:

public:

  StochSymDummyMatrix(int id_)
    : StochSymMatrix(id_, 0, 0, 0, MPI_COMM_NULL) {};

  virtual ~StochSymDummyMatrix(){};

  virtual StochSymDummyMatrix* clone() const { return new StochSymDummyMatrix(id); };

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

  virtual void getSize( long long& m, long long& n ){m=0;n=0;}
  virtual void getSize( int& m, int& n ){m=0;n=0;}

  virtual long long size(){return 0;}

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
  
  virtual void deleteEmptyRowsCols(const OoqpVector& nnzVec, const OoqpVector* linkParent) {};
  virtual void deleteEmptyRowsCols(const OoqpVector& nnzVec) {};
};

typedef SmartPointer<StochSymMatrix> StochSymMatrixHandle;

#endif
