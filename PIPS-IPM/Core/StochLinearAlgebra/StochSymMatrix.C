#include "StochSymMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixTypes.h"

#include <cassert>

using namespace std;


StochSymMatrix::StochSymMatrix(int id, int global_n, int local_n, int local_nnz, 
			       MPI_Comm mpiComm_)
  :id(id), n(global_n), mpiComm(mpiComm_)
{
  mat = new SparseSymMatrix(local_n, local_nnz);
}


StochSymMatrix::StochSymMatrix(const vector<StochSymMatrix*> &blocks) 
{
  mpiComm = blocks[0]->mpiComm;
  n = blocks[0]->n;
  id = blocks[0]->id;

  vector<SparseSymMatrix*> v(blocks.size());
  for(int i = 0; i < blocks.size(); i++) v[i] = blocks[i]->mat;
  mat = new SparseSymMatrix(v);

}


void StochSymMatrix::AddChild(StochSymMatrix* child)
{
  children.push_back(child);
}

StochSymMatrix::~StochSymMatrix()
{
  for(int it=0; it<children.size(); it++)
    delete children[it];

  if (mat) delete mat;
}

int StochSymMatrix::isKindOf( int type )
{
  return type == kStochSymMatrix || type == kSymMatrix;
}

void StochSymMatrix::atPutDense( int /* row */, int /* col */,
				 double * /* A */, int /* lda */,
				 int /* rowExtent */,
				 int /* colExtent */ )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::fromGetDense( int row, int col, double * A, int lda,
		   int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::symAtPutSpRow( int row, 
				    double A[], int lenA, int jcolA[],
				    int& info )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::fsymAtPutSpRow( int row, 
				     double A[], int lenA, int jcolA[],
				     int& info )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::getSize( int& m_, int& n_ )
{
  m_=n; n_=n;
}

int StochSymMatrix::size()
{
  return n;
}

void StochSymMatrix::symAtPutSubmatrix( int destRow, int destCol,
					DoubleMatrix& M,
					int srcRow, int srcCol,
					int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::fromGetSpRow( int row, int col,
				   double A[], int lenA, int irowA[], int& nnz,
				   int rowExtent, int& info )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::atPutZeros( int row, int col, int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

/** y = beta * y + alpha * this * x */
void StochSymMatrix::mult ( double beta,  OoqpVector& y_,
			    double alpha, OoqpVector& x_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  //check the tree compatibility
  int nChildren = children.size();
  assert(y.children.size() == nChildren);
  assert(x.children.size() == nChildren);

  //check the node size compatibility
  assert(this->mat->size() == y.vec->length());
  assert(this->mat->size() == x.vec->length());

  if (0.0 == alpha) {
    y.vec->scale( beta );
    return;
  } else {
    //if( alpha != 1.0 || beta != 1.0 ) {
    //  y.vec->scale( beta/alpha );
    //}

    mat->mult( beta, *y.vec, alpha, *x.vec );

    //if( 1.0 != alpha ) {
    //  y.vec->scale(alpha);
    //}
  }


  // reccursively multiply the children
  for (int it=0; it<nChildren; it++) {
    children[it]->mult(beta, *(y.children[it]), alpha, *(x.children[it]));
  }
}

/** y = beta * y + alpha * this^T * x */
void StochSymMatrix::transMult ( double beta,  OoqpVector& y_,
				 double alpha, OoqpVector& x_)
{
  // We are symmetric, this^T = this, therefore call 'mult' method
  this->mult(beta, y_, alpha, x_);
}
  
/** the magnitude of the element in this matrix with largest absolute value.
   */
double StochSymMatrix::abmaxnorm()
{
  double maxNorm=0.0, localMaxNorm, childMaxNorm;

  //!parallel stuff needed

  localMaxNorm = mat->abmaxnorm();
  maxNorm = max(localMaxNorm, maxNorm);
  
  for (int it=0; it<children.size(); it++) {
    childMaxNorm = children[it]->abmaxnorm();
    maxNorm = max(childMaxNorm, maxNorm);
  }

  return maxNorm;
}

void StochSymMatrix::writeToStream(ostream& out) const
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::randomizePSD(double * seed)
{
  assert( "Not implemented" && 0 );
}


void StochSymMatrix::getDiagonal( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  mat->getDiagonal( *vec.vec);

  for(int it=0; it<children.size(); it++)
    children[it]->getDiagonal(*vec.children[it]);
}

void StochSymMatrix::setToDiagonal( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  mat->setToDiagonal( *vec.vec);

  for(int it=0; it<children.size(); it++)
    children[it]->setToDiagonal(*vec.children[it]);
}

void StochSymMatrix::atPutDiagonal( int idiag, OoqpVector& v_ )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);

  //check the tree compatibility
  int nChildren = children.size();
  assert(v.children.size() == nChildren);

  //check the node size compatibility
  assert(this->mat->size() == v.vec->length());

  mat->atPutDiagonal ( idiag, *v.vec);

  for (int it=0; it<nChildren; it++) 
    children[it]->atPutDiagonal( idiag, *v.children[it]);
}

void StochSymMatrix::fromGetDiagonal( int idiag, OoqpVector& x_ )
{
  assert("The value of the parameter is not supported!" && idiag==0);

  StochVector& x = dynamic_cast<StochVector&>(x_);
  assert(x.children.size() == children.size());

  mat->getDiagonal(*x.vec);

  for (int it=0; it<children.size(); it++) 
    children[it]->getDiagonal(*x.children[it]);
}

void StochSymMatrix::putSparseTriple( int irow[], int len, int jcol[], 
				      double A[], int& info )
{
  assert("Not implemented!" && 0);
}

void StochSymMatrix::SymmetricScale( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  mat->SymmetricScale(*vec.vec);

  for (int it=0; it<children.size(); it++) 
    children[it]->SymmetricScale(*vec.children[it]);
}

void StochSymMatrix::ColumnScale( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  mat->ColumnScale(*vec.vec);

  for (int it=0; it<children.size(); it++) 
    children[it]->ColumnScale(*vec.children[it]);
}

void StochSymMatrix::RowScale ( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  mat->RowScale(*vec.vec);

  for (int it=0; it<children.size(); it++) 
    children[it]->RowScale(*vec.children[it]);
}


void StochSymMatrix::scalarMult( double num )
{
  mat->scalarMult(num);
  for (int it=0; it<children.size(); it++) 
    children[it]->scalarMult(num);
}



int StochSymDummyMatrix::isKindOf( int type ) 
{ 
  return type==kStochSymDummyMatrix;
}
