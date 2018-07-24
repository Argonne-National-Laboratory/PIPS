#include "StochGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "StochVector.h"
#include "SimpleVector.h"

using namespace std;

StochGenMatrix::StochGenMatrix(int id, 
			       long long global_m, long long global_n,
			       int A_m, int A_n, int A_nnz,
			       int B_m, int B_n, int B_nnz,
			       MPI_Comm mpiComm_)
  : id(id), m(global_m), n(global_n), 
    mpiComm(mpiComm_), iAmDistrib(0),
    workPrimalVec(NULL)
{
  //cout << "StochGenMatrix-> " << A_m << " " << A_n << " " << B_m << " " << B_n << " " << endl;

  Amat = new SparseGenMatrix(A_m, A_n, A_nnz);
  Bmat = new SparseGenMatrix(B_m, B_n, B_nnz);

  if(mpiComm!=MPI_COMM_NULL) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) iAmDistrib=1;
  }
}

/*StochGenMatrix::StochGenMatrix(const vector<StochGenMatrix*> &blocks) 
  : iAmDistrib(0), workPrimalVec(NULL)
{
  mpiComm = blocks[0]->mpiComm;
  n = blocks[0]->n;
  m = blocks[0]->m;
  id = blocks[0]->id;

  if(mpiComm!=MPI_COMM_NULL) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) iAmDistrib=1;
  }

  vector<SparseGenMatrix*> v(blocks.size());

  for (size_t i = 0; i < blocks.size(); i++) v[i] = blocks[i]->Amat;
  Amat = new SparseGenMatrix(v, false);
  for (size_t i = 0; i < blocks.size(); i++) v[i] = blocks[i]->Bmat;
  Bmat = new SparseGenMatrix(v, true);
}
*/

StochGenMatrix::~StochGenMatrix()
{
  //cout << "~~~~~~~~StochGenMatrix" << endl;
  for(size_t it=0; it<children.size(); it++)
    delete children[it];

  if (Amat)
    delete Amat;

  if (Bmat)
    delete Bmat;

  if(workPrimalVec)
    delete workPrimalVec;
}

void StochGenMatrix::AddChild(StochGenMatrix* child)
{
  children.push_back(child);
}

OoqpVector* StochGenMatrix::getWorkPrimalVec(const StochVector& origin)
{
  if(NULL == workPrimalVec)
    workPrimalVec = origin.dataClone();
  else
    assert(workPrimalVec->length() == origin.vec->length());
  return workPrimalVec;
}

int StochGenMatrix::isKindOf( int type )
{
  return type == kStochGenMatrix || type == kGenMatrix;
}

int StochGenDummyMatrix::isKindOf( int type )
{
  return type == kStochGenDummyMatrix;
}
void StochGenMatrix::getSize( long long& m_out, long long& n_out )
{
  m_out = m; n_out=n;
}
void StochGenMatrix::getSize( int& m_out, int& n_out )
{
  m_out = m; n_out=n;
}


void StochGenMatrix::atPutDense( int row, int col, double * A, int lda,
				 int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::fromGetDense( int row, int col, double * A, int lda,
				   int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::ColumnScale( OoqpVector& vec )
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::RowScale( OoqpVector& vec )
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::SymmetricScale( OoqpVector &vec)
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::scalarMult( double num)
{
  Amat->scalarMult(num);
  Bmat->scalarMult(num);
  for (size_t it=0; it<children.size(); it++) 
    children[it]->scalarMult(num);
}

void StochGenMatrix::fromGetSpRow( int row, int col,
				   double A[], int lenA, int jcolA[], int& nnz,
				   int colExtent, int& info )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
				     int srcRow, int srcCol,
				     int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::atPutSpRow( int col, 
				 double A[], int lenA, int jcolA[],
				 int& info )
{
  assert( "Not implemented" && 0 );
}


void StochGenMatrix::putSparseTriple( int irow[], int len, int jcol[], double A[], 
				      int& info )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::getDiagonal( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);

  //check compatibility
  assert( children.size() == vec.children.size() );

  // only local B gives the diagonal
  Bmat->getDiagonal(*vec.vec);

  //do it recursively
  for(size_t it=0; it<children.size(); it++)
    children[it]->getDiagonal(*vec.children[it]);
}
 
void StochGenMatrix::setToDiagonal( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  Bmat->setToDiagonal( *vec.vec);

  for(size_t it=0; it<children.size(); it++)
    children[it]->setToDiagonal(*vec.children[it]);
}

/* y = beta * y + alpha * this * x */
void StochGenMatrix::mult( double beta,  OoqpVector& y_,
			   double alpha, OoqpVector& x_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  //check the tree compatibility
  int nChildren = children.size();
  assert(y.children.size() - nChildren == 0);
  assert(x.children.size() - nChildren == 0);

  if (0.0 == alpha) {
    y.vec->scale( beta );
    return;
  } else {
    //if( alpha != 1.0 || beta != 1.0 ) {
    //  y.vec->scale( beta/alpha );
    //}

    Bmat->mult(beta, *y.vec, alpha, *x.vec);

    long long mA, nA; 
    Amat->getSize(mA,nA);
    if(nA>0) {
      //not the root
      Amat->mult(1.0, *y.vec, alpha, *x.parent->vec);
    }

    //if( 1.0 != alpha ) {
    //  y.vec->scale(alpha);
    //}
  }

  for(size_t it=0; it<children.size(); it++)
    children[it]->mult(beta, *y.children[it], alpha, *x.children[it]);
}

void StochGenMatrix::transMult ( double beta,   OoqpVector& y_,
				 double alpha,  OoqpVector& x_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

#ifdef DEBUG
  //check the tree compatibility
  int nChildren = children.size();
  assert(y.children.size() == nChildren);
  assert(x.children.size() == nChildren);
#endif

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);
  
  //preparations for the parallel case
  int iAmSpecial = 1;
  if(iAmDistrib) {
    int rank; MPI_Comm_rank(mpiComm, &rank);
    if(rank>0) iAmSpecial = 0;
  }

  if(iAmSpecial)
    //y_i = beta* y_i  +  alpha* B_i^T* x_i
    Bmat->transMult(beta, yvec, alpha, xvec);
  else
    yvec.setToZero();


  //!opt alloc buffer here and send it through the tree to be used by
  //!children when MPI_Allreduce

  //let the children compute their contribution
  for(size_t it=0; it<children.size(); it++) {
      children[it]->transMult2(beta, *y.children[it], alpha, *x.children[it], yvec);
  }

  if(iAmDistrib) {
    int locn=yvec.length();
    double* buffer = new double[locn];

    MPI_Allreduce(yvec.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
    yvec.copyFromArray(buffer);

    delete[] buffer;
  }
}

void StochGenMatrix::transMult2 ( double beta,   StochVector& y,
				  double alpha,  StochVector& x,
				  OoqpVector& yvecParent)
{
  //check the tree compatibility
  int nChildren = children.size();
  assert(y.children.size() - nChildren == 0);
  assert(x.children.size() - nChildren == 0);

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

#ifdef DEBUG
  int nA, mA;
  Amat->getSize(mA, nA);
  // this should NOT be the root
  assert(nA>0);
#endif

  //do A_i^T x_i and add it to yvecParent which already contains B_0^T x_0
  Amat->transMult(1.0, yvecParent, alpha, *x.vec);
  

  //preparations for the parallel case
  int iAmSpecial = 1;
  if(iAmDistrib) {
    int rank; MPI_Comm_rank(mpiComm, &rank);
    if(rank>0) iAmSpecial = 0;
  }

  if(iAmSpecial)
    //do A_i^T x_i and add it to yvecParent which already contains
    //B_0^T x_0
    Bmat->transMult(beta, yvec, alpha, xvec);
  else
    yvec.setToZero();
  
  for(size_t it=0; it<children.size(); it++)
    children[it]->transMult2(beta, *y.children[it], 
			     alpha, *x.children[it],
			     yvec);

  if(iAmDistrib) {
    int locn=yvec.length();
    double* buffer = new double[locn];

    MPI_Allreduce(yvec.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
    yvec.copyFromArray(buffer);

    delete[] buffer;
  }
}




double StochGenMatrix::abmaxnorm()
{
  double nrm = 0.0;
  
  for(size_t it=0; it<children.size(); it++)
    nrm = max(nrm, children[it]->abmaxnorm());

  if(iAmDistrib) {
    double nrmG=0;
    MPI_Allreduce(&nrm, &nrmG, 1, MPI_DOUBLE, MPI_MAX, mpiComm);
    nrm=nrmG;
  }

  nrm = max(nrm, max(Amat->abmaxnorm(), Bmat->abmaxnorm()));
  return nrm;
}

void StochGenMatrix::writeToStream(ostream& out) const
{
  //assert( "Has not been yet implemented" && 0 );
  printf("Amat:\n");
  Amat->writeToStream(out);
  printf("Bmat:\n");
  Bmat->writeToStream(out);
  for (int i = 0; i < children.size(); ++i) {
    printf("children[%d]->Amat:\n", i);
    children[i]->Amat->writeToStream(out);
    printf("children[%d]->Bmat:\n", i);
    children[i]->Bmat->writeToStream(out);
  }
}


/* Make the elements in this matrix symmetric. The elements of interest
 *  must be in the lower triangle, and the upper triangle must be empty.
 *  @param info zero if the operation succeeded. Otherwise, insufficient
 *  space was allocated to symmetrize the matrix.
 */
void StochGenMatrix::symmetrize( int& info )
{
  assert( "Has not been yet implemented" && 0 );
}


void StochGenMatrix::randomize( double alpha, double beta, double * seed )
{
  assert( "Has not been yet implemented" && 0 );
}


void StochGenMatrix::atPutDiagonal( int idiag, OoqpVector& v )
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::fromGetDiagonal( int idiag, OoqpVector& v )
{
  assert( "Has not been yet implemented" && 0 );
}

int StochGenMatrix::numberOfNonZeros()
{
  int nnz = 0;

  for(size_t it=0; it<children.size(); it++)
    nnz += children[it]->numberOfNonZeros();

  if(iAmDistrib) {
    int nnzG = 0;
    MPI_Allreduce(&nnz, &nnzG, 1, MPI_INT, MPI_SUM, mpiComm);
    nnz=nnzG;
  }

  nnz += Amat->numberOfNonZeros() + Bmat->numberOfNonZeros();
  return nnz;
}

void StochGenMatrix::matTransDMultMat(OoqpVector& d, SymMatrix** res)
{
  assert( "Has not been yet implemented" && 0 );
}
void StochGenMatrix::matTransDinvMultMat(OoqpVector& d, SymMatrix** res)
{
  assert( "Has not been yet implemented" && 0 );
}
