/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#ifndef SPARSEGENMATRIX_SEQ_H
#define SPARSEGENMATRIX_SEQ_H

#include "OoqpVectorHandle.h"
#include "DoubleMatrix.h"
#include "SparseStorage.h"
#include "SparseGenMatrixHandle.h"

#include "SimpleVector.h"
#include "StochVector.h"

#include "mpi.h"

//#include "sTreeMultiStage.h"

#include "SparseGenMatrix.h"


#include "multiStageInputTree.hpp"



class sTreeMultiStage;


class Variables;

/** Represents sparse non-symmetric, possibly non-square matrices stored in
 *  row-major Harwell-Boeing format.
 *  @ingroup SparseLinearAlgebra
 */
class SeqSparseGenMatrix : public SparseGenMatrix {
protected:
 // SparseStorageHandle mStorage;
  int size;
  sTreeMultiStage* currTree;
  int iAmDistrib;

  bool isEmpty;
  
public:

  int id;
  MPI_Comm mpiComm;
	
//  SeqSparseGenMatrix( int id, int rows, int cols, int nnz );
  SeqSparseGenMatrix( int id_in, int rows, int cols, int nnz, MPI_Comm mpiComm_,
  							multiStageInputTree* in_);
  SeqSparseGenMatrix( int id_in, int rows, int cols, int nnz, sTreeMultiStage* tree_in, MPI_Comm mpiComm_,
  							multiStageInputTree* in_);
  ~SeqSparseGenMatrix();

  SeqSparseGenMatrix* parent;

  SeqSparseGenMatrix* child;

  std::vector<SeqSparseGenMatrix*> parents;

//  virtual void getSize( long long& m, long long& n );
//  virtual void getSize( int& m, int& n );
  virtual int numberOfNonZeros();

  virtual int isKindOf( int matType );

  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ){assert( "Not implemented" && 0 );}
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ){assert( "Not implemented" && 0 );}
  virtual void ColumnScale( OoqpVector& vec ) {assert( "Not implemented" && 0 );}
  virtual void RowScale( OoqpVector& vec ){assert( "Not implemented" && 0 );}
  virtual void SymmetricScale( OoqpVector &vec){assert( "Not implemented" && 0 );}
  virtual void scalarMult( double num);
  virtual void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info ){assert( "Not implemented" && 0 );}
  virtual void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent ){assert( "Not implemented" && 0 );}
  virtual void atPutSpRow( int col, double A[], int lenA, int jcolA[],
			   int& info ){assert( "Not implemented" && 0 );}

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info ){assert( "Not implemented" && 0 );}

  virtual void getDiagonal( OoqpVector& vec ){assert( "Not implemented" && 0 );}
  virtual void setToDiagonal( OoqpVector& vec ){assert( "Not implemented" && 0 );}

  virtual void mult ( double beta,  OoqpVector& y,
                      double alpha, OoqpVector& x );
  virtual void mult ( double beta,  double y[], int incy,
                      double alpha, double x[], int incx ){assert( "Not implemented" && 0 );}

  virtual void transMult( double beta,   OoqpVector& y,
			  double alpha,  OoqpVector& x );
  virtual void transMult( double beta,  OoqpVector& y_in, int incy,
			  double alpha, OoqpVector& x_in, int incx ){assert( "Not implemented" && 0 );}
  virtual void transMult( double beta,  double y_in[], int incy,
			  double alpha, double x_in[], int incx ){assert( "Not implemented" && 0 );}

  /** C = this^T * D * this where D=diag(d) is a diagonal matrix. */
  virtual void matTransDMultMat(OoqpVector& d, SymMatrix** res){assert( "Not implemented" && 0 );}
  /** C = this^T * inv(D) * this where D=diag(d) is a diagonal matrix. */
  virtual void matTransDinvMultMat(OoqpVector& d, SymMatrix** res){assert( "Not implemented" && 0 );}

  /** C = this * this^T */
  virtual void matMultTrans(SymMatrix** res){assert( "Not implemented" && 0 );}
  
  virtual double abmaxnorm();
  
  virtual void writeToStream(ostream& out) const {assert( "Not implemented" && 0 );}

  /** Make the elements in this matrix symmetric. The elements of interest
   *  must be in the lower triangle, and the upper triangle must be empty.
   *  @param info zero if the operation succeeded. Otherwise, insufficient
   *  space was allocated to symmetrize the matrix.
   */
  virtual void symmetrize( int& info ){assert( "Not implemented" && 0 );}

  virtual void randomize( double alpha, double beta, double * seed ){assert( "Not implemented" && 0 );}

  virtual void atPutDiagonal( int idiag, OoqpVector& v ){assert( "Not implemented" && 0 );}
  virtual void fromGetDiagonal( int idiag, OoqpVector& v ){assert( "Not implemented" && 0 );}

  SparseStorage * getStorage() { return mStorage.ptr(); }
  SparseStorage& getStorageRef() { return *mStorage; }
  int * krowM() { return mStorage->krowM; }
  int * jcolM() { return mStorage->jcolM; }
  double * M() { return mStorage->M; }

  virtual void copyMtxFromDouble(int copyLength,double *values){assert( "Not implemented" && 0 );}


  virtual void setAdditiveDiagonal(OoqpVector& diag){assert( "Not implemented" && 0 );}

  virtual void fromGetSpRow_WithRowStart( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int colExtent, int& info, int & rowStart){assert( "Not implemented" && 0 );}

  virtual double * getMatVal(){return mStorage->M; }

  virtual void fromGetDense_withMap( int row, int col, double * A, int lda,
					int rowExtent, int colExtent, int const FirstCall, std::map<int,int> &ValIdxMap ){assert( "Not implemented" && 0 );}

  virtual void printMatrixInMatlab(char *name){assert( "Not implemented" && 0 );}


  virtual void mult ( double beta,  StochVector& y,
                      double alpha, StochVector& x );	

  virtual void transMult( double beta,	StochVector& y,
			 double alpha,	StochVector& x );

  virtual void transMult( const int setA, double beta,	StochVector& y,
			 double alpha,	StochVector& x, StochVector* End_Par_Pos );

  virtual void addtransMultAtTarget ( const int setA,
				   double alpha,  StochVector& x_, SimpleVector* goal_Par );


  void AddParent(SeqSparseGenMatrix* parent_);


  bool parentIsEmpty;


 protected:
  // in the case of A'*A we internally form the transpose only once
  SparseGenMatrix* m_Mt;

  multiStageInputTree* in;

  virtual bool checkIfParentEmpty();
};

#endif

