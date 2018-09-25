#include "StochGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "pipsdef.h"
#include <limits>
#include <algorithm>
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
  Blmat = new SparseGenMatrix(0, 0, 0);

  if(mpiComm!=MPI_COMM_NULL) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) iAmDistrib=1;
  }
}

StochGenMatrix::StochGenMatrix(int id,
			       long long global_m, long long global_n,
			       int A_m, int A_n, int A_nnz,
			       int B_m, int B_n, int B_nnz,
				   int Bl_m, int Bl_n, int Bl_nnz,
			       MPI_Comm mpiComm_)
  : id(id), m(global_m), n(global_n),
    mpiComm(mpiComm_), iAmDistrib(0),
    workPrimalVec(NULL)
{
  Amat = new SparseGenMatrix(A_m, A_n, A_nnz);
  Bmat = new SparseGenMatrix(B_m, B_n, B_nnz);
  Blmat = new SparseGenMatrix(Bl_m, Bl_n, Bl_nnz);

  if(mpiComm!=MPI_COMM_NULL) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) iAmDistrib=1;
  }
}

StochGenMatrix::StochGenMatrix(int id,
                long long global_m, long long global_n,
                MPI_Comm mpiComm_)
  : id(id), m(global_m), n(global_n),
    mpiComm(mpiComm_), iAmDistrib(0),
    workPrimalVec(NULL)
{
  Amat = NULL;
  Bmat = NULL;
  Blmat = NULL;

  if(mpiComm!=MPI_COMM_NULL) {
    int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) iAmDistrib=1;
  }
}

StochGenMatrix::~StochGenMatrix()
{
  //cout << "~~~~~~~~StochGenMatrix" << endl;
  for(size_t it=0; it<children.size(); it++)
    delete children[it];

  if (Amat)
    delete Amat;

  if (Bmat)
    delete Bmat;

  if (Blmat)
    delete Blmat;

  if(workPrimalVec)
    delete workPrimalVec;
}


StochGenMatrix* StochGenMatrix::cloneFull(bool switchToDynamicStorage) const
{
   StochGenMatrix* clone = new StochGenMatrix(id, m, n, mpiComm);

   // clone submatrices
   clone->Amat = Amat->cloneFull(switchToDynamicStorage);
   clone->Bmat = Bmat->cloneFull(switchToDynamicStorage);
   clone->Blmat = Blmat->cloneFull(switchToDynamicStorage);

   for( size_t it = 0; it < children.size(); it++ )
      clone->children.push_back(children[it]->cloneFull(switchToDynamicStorage));

   return clone;
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

void StochGenMatrix::ColumnScale2( OoqpVector& vec, OoqpVector& parentvec )
{
   StochVector& scalevec = dynamic_cast<StochVector&>(vec);
   SimpleVector& scalevecparent = dynamic_cast<SimpleVector&>(parentvec);

   assert(scalevec.children.size() == 0 && children.size() == 0);

   Amat->ColumnScale(scalevecparent);
   Bmat->ColumnScale(*scalevec.vec);
   Blmat->ColumnScale(*scalevec.vec);
}

void StochGenMatrix::ColumnScale( OoqpVector& vec )
{
   StochVector& scalevec = dynamic_cast<StochVector&>(vec);

   assert(children.size() == scalevec.children.size());

   Bmat->ColumnScale(*scalevec.vec);
   Blmat->ColumnScale(*scalevec.vec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->ColumnScale2(*(scalevec.children[it]), *scalevec.vec);
}

void StochGenMatrix::RowScale2( OoqpVector& vec, OoqpVector* linkingvec )
{
   StochVector& scalevec = dynamic_cast<StochVector&>(vec);

   assert(scalevec.children.size() == 0 && children.size() == 0);

   Amat->RowScale(*scalevec.vec);
   Bmat->RowScale(*scalevec.vec);

   if( linkingvec )
   {
      SimpleVector* vecl = dynamic_cast<SimpleVector*>(linkingvec);
      Blmat->RowScale(*vecl);
   }
}

void StochGenMatrix::RowScale( OoqpVector& vec )
{
   StochVector& scalevec = dynamic_cast<StochVector&>(vec);

   assert(children.size() == scalevec.children.size());

   Bmat->RowScale(*scalevec.vec);

   SimpleVector* vecl = dynamic_cast<SimpleVector*>(scalevec.vecl);

   if( vecl )
      Blmat->RowScale(*vecl);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->RowScale2(*(scalevec.children[it]), vecl);
}

void StochGenMatrix::SymmetricScale( OoqpVector &vec)
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::scalarMult( double num)
{
  Amat->scalarMult(num);
  Bmat->scalarMult(num);
  Blmat->scalarMult(num);

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

  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

  if (0.0 == alpha) {
    y.vec->scale( beta );

    for(size_t it=0; it<children.size(); it++)
      children[it]->mult(beta, *y.children[it], alpha, *x.children[it]);

    return;
  } else {
    //if( alpha != 1.0 || beta != 1.0 ) {
    //  y.vec->scale( beta/alpha );
    //}

    Bmat->mult(beta, yvec, alpha, xvec);

    long long mA, nA; 
    Amat->getSize(mA,nA);

    // not at root?
    if(nA>0) {
      Amat->mult(1.0, yvec, alpha, *x.parent->vec);

    }

    //if( 1.0 != alpha ) {
    //  y.vec->scale(alpha);
    //}
  }

  int blm, bln;
  Blmat->getSize(blm, bln);

  /* linking constraints present? */
  if( blm > 0 )
  {
	SimpleVector& yvecl = dynamic_cast<SimpleVector&>(*y.vecl);

	if( iAmSpecial(iAmDistrib, mpiComm) )
	  Blmat->mult(beta, yvecl, alpha, xvec);
	else
	  yvecl.setToZero();

    for(size_t it=0; it<children.size(); it++)
	   children[it]->mult2(beta, *y.children[it], alpha, *x.children[it], yvecl);

	if(iAmDistrib) {
	  // sum up linking constraints vectors
	  int locn=yvecl.length();
	  double* buffer = new double[locn];

	  MPI_Allreduce(yvecl.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
	  yvecl.copyFromArray(buffer);

	  delete[] buffer;
    }
  }
  else
  {
    for(size_t it=0; it<children.size(); it++)
      children[it]->mult(beta, *y.children[it], alpha, *x.children[it]);
  }
}


/* mult method for children; needed only for linking constraints */
void StochGenMatrix::mult2( double beta,  OoqpVector& y_,
			   double alpha, OoqpVector& x_, OoqpVector& yparentl_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  if( 0.0 == alpha ) {
    y.vec->scale( beta );
    return;
  }

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

  Bmat->mult(beta, yvec, alpha, xvec);
  Blmat->mult(1.0, yparentl_, alpha, xvec);
  Amat->mult(1.0, *y.vec, alpha, *x.parent->vec);

  // not implemented
  assert(children.size() == 0);
}


void StochGenMatrix::transMult ( double beta,   OoqpVector& y_,
				 double alpha,  OoqpVector& x_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  // assert tree compatibility
  assert(y.children.size() == children.size());
  assert(x.children.size() == children.size());

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

  int blm, bln;
  Blmat->getSize(blm, bln);

  // with linking constraints?
  if( blm > 0 )
  {
    assert(x.vecl);
    SimpleVector& xvecl = dynamic_cast<SimpleVector&>(*x.vecl);

    if( iAmSpecial(iAmDistrib, mpiComm) )
    {
      //y_i = beta* y_i  +  alpha* B_i^T* x_i
      Bmat->transMult(beta, yvec, alpha, xvec);

	   //y_i = y_i  +  alpha* Bl_0^T* xl_i
	   Blmat->transMult(1.0, yvec, alpha, xvecl);
    }
    else
    {
      yvec.setToZero();
    }

    //!opt alloc buffer here and send it through the tree to be used by
    //!children when MPI_Allreduce
    //let the children compute their contribution
    for(size_t it=0; it<children.size(); it++) {
      children[it]->transMult2(beta, *y.children[it], alpha, *x.children[it], yvec, xvecl);
    }
  }
  else // no linking constraints
  {
    if( iAmSpecial(iAmDistrib, mpiComm) )
      Bmat->transMult(beta, yvec, alpha, xvec);
    else
      yvec.setToZero();

    for(size_t it=0; it<children.size(); it++) {
      children[it]->transMult2(beta, *y.children[it], alpha, *x.children[it], yvec);
    }
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
				  OoqpVector& yvecParent, OoqpVector& xvecl)
{
  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

#ifdef STOCH_TESTING
  int nA, mA;
  Amat->getSize(mA, nA);
  // this should NOT be the root
  assert(nA>0);
#endif

  //do A_i^T x_i and add it to yvecParent which already contains B_0^T x_0
  Amat->transMult(1.0, yvecParent, alpha, *x.vec);

  if( iAmSpecial(iAmDistrib, mpiComm) )
  {
    Bmat->transMult(beta, yvec, alpha, xvec);

	 //y_i = y_i  +  alpha* Bl_i^T* xl_i
	 Blmat->transMult(1.0, yvec, alpha, xvecl);
  }
  else
    yvec.setToZero();
  
  assert(children.size() == 0);
#if 0
  for(size_t it=0; it<children.size(); it++)
    children[it]->transMult2(beta, *y.children[it], 
			     alpha, *x.children[it],
			     yvec);
#endif


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
  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

#ifdef STOCH_TESTING
  int nA, mA;
  Amat->getSize(mA, nA);
  // this should NOT be the root
  assert(nA>0);
#endif

  //do A_i^T x_i and add it to yvecParent which already contains B_0^T x_0
  Amat->transMult(1.0, yvecParent, alpha, *x.vec);

  if( iAmSpecial(iAmDistrib, mpiComm) )
  {
    //do A_i^T x_i and add it to yvecParent which already contains
    //B_0^T x_0
    Bmat->transMult(beta, yvec, alpha, xvec);
  }
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
  nrm = max(nrm, Blmat->abmaxnorm());

  return nrm;
}

void StochGenMatrix::getLinkVarsNnz(std::vector<int>& vec) const
{
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->getLinkVarsNnzChild(vec);

   if( iAmDistrib )
   {
      int* buffer = new int[vec.size()];
      MPI_Allreduce(&vec[0], buffer, vec.size(), MPI_INT, MPI_SUM, mpiComm);

      std::memcpy(&vec[0], buffer, vec.size() * sizeof(int));

      delete[] buffer;
   }
}

void StochGenMatrix::getLinkVarsNnzChild(std::vector<int>& vec) const
{
   assert(children.size() == 0);

   Amat->getLinkVarsNnz(vec);
}

void StochGenMatrix::writeToStream(ostream& out) const
{
  assert( "Has not been yet implemented" && 0 );
}

void
StochGenMatrix::writeToStreamDense(ostream& out) const
{
   int rank;
   MPI_Comm_rank(mpiComm, &rank);
   int world_size;
   MPI_Comm_size(mpiComm, &world_size);
   int m, n;
   int offset = 0;
   stringstream sout;
   MPI_Status status;
   int l;

   if( iAmDistrib )
            MPI_Barrier(mpiComm);

   if( iAmDistrib && rank > 0 )  // receive offset from previous process
      MPI_Recv(&offset, 1, MPI_INT, (rank - 1), 0, mpiComm, MPI_STATUS_IGNORE);

   else  //  !iAmDistrib || (iAmDistrib && rank == 0)
      this->Bmat->writeToStreamDense(out);

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->writeToStreamDenseChild(sout, offset);
      children[it]->Bmat->getSize(m, n);
      offset += n;
   }

   if( iAmDistrib && rank > 0 )
   {
      std::string str = sout.str();
      // send string to rank ZERO to print it there:
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, rank, mpiComm);
      // send offset to next process:
      if( rank < world_size - 1 )
         MPI_Ssend(&offset, 1, MPI_INT, rank + 1, 0, mpiComm);
   }

   else if( !iAmDistrib )
      out << sout.str();

   else if( iAmDistrib && rank == 0 )
   {
      out << sout.str();
      MPI_Ssend(&offset, 1, MPI_INT, rank + 1, 0, mpiComm);

      for( int p = 1; p < world_size; p++ )
      {
         MPI_Probe(p, p, mpiComm, &status);
         MPI_Get_count(&status, MPI_CHAR, &l);
         char *buf = new char[l];
         MPI_Recv(buf, l, MPI_CHAR, p, p, mpiComm, &status);
         string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }

   }

   // linking contraints ?
   int mlink, nlink;
   this->Blmat->getSize(mlink, nlink);
   if( mlink > 0 )
   {
      if( iAmDistrib )
         MPI_Barrier(mpiComm);

      // for each row r do:
      for( int r = 0; r < mlink; r++ )
      {
         if( iAmDistrib )
         {
            MPI_Barrier(mpiComm);

            // process Zero collects all the information and then prints it.
            if( rank == 0 )
            {
               out << this->Blmat->writeToStreamDenseRow(r);

               out << writeToStreamDenseRowLink(r);

               for( int p = 1; p < world_size; p++ )
               {
                  MPI_Probe(p, r + 1, mpiComm, &status);
                  MPI_Get_count(&status, MPI_CHAR, &l);
                  char *buf = new char[l];
                  MPI_Recv(buf, l, MPI_CHAR, p, r + 1, mpiComm, &status);
                  string rowPartFromP(buf, l);
                  out << rowPartFromP;
                  delete[] buf;
               }
               out << endl;

            }
            else // rank != 0
            {
               std::string str = writeToStreamDenseRowLink(r);
               MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, r + 1, mpiComm);
            }
         }
         else // not distributed
         {
            stringstream sout;
            this->Blmat->writeToStreamDenseRow(sout, r);
            for( size_t it = 0; it < children.size(); it++ )
               children[it]->Blmat->writeToStreamDenseRow(sout, r);

            out << sout.rdbuf() << endl;
         }
      }
   }
   if( iAmDistrib )
      MPI_Barrier(mpiComm);
}

/** writes child matrix blocks, offset indicates the offset between A and B block. */
void StochGenMatrix::writeToStreamDenseChild(stringstream& out, int offset) const
{
   int mA, mB, n;
   this->Amat->getSize(mA, n);
   this->Bmat->getSize(mB, n);
   assert(mA == mB );
   for(int r=0; r < mA; r++)
   {
      this->Amat->writeToStreamDenseRow(out, r);
      for(int i=0; i<offset; i++)
         out <<'\t';
      this->Bmat->writeToStreamDenseRow(out, r);
      out << endl;
   }
}

/** returns a string containing the linking-row rowidx of the children. */
std::string StochGenMatrix::writeToStreamDenseRowLink(int rowidx) const
{
   std::string str_all;
   for( size_t it = 0; it < children.size(); it++ )
   {
      std::string str = children[it]->Blmat->writeToStreamDenseRow(rowidx);
      str_all.append(str);
   }
   return str_all;
}

void StochGenMatrix::writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   string rt;
   if( rowType == 0 )
      rt = "E";
   else if( rowType == 1 )
      rt = "L";
   else if( rowType == 2 )
      rt = "G";
   else
      assert(0);

   StochVector* irhsStoch;
   if( irhs )
      irhsStoch = dynamic_cast<StochVector*>(irhs);

   int m, n;
   if( myRank == 0)
   {
      // A_0 block:
      this->Bmat->getSize(m, n);
      for(int i=0; i<m; i++)
      {
         if( !irhs || (irhs && dynamic_cast<SimpleVector*>(irhsStoch->vec)->elements()[i] != 0.0) )
            out<< " "<<rt<<" row_"<<rt<<"_"<<"R" <<"_"<<i <<endl;
      }
      // linking rows:
      if( Blmat )
      {
         this->Blmat->getSize(m, n);
         for(int i=0; i<m; i++)
         {
            if( !irhs || (irhs && dynamic_cast<SimpleVector*>(irhsStoch->vecl)->elements()[i] != 0.0) )
               out<<" "<< rt<<" row_"<<rt<<"_"<<"L" <<"_"<<i <<endl;
         }
      }
   }
   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->Amat->getSize(m, n);
      for(int i=0; i<m; i++)
      {
         if( !irhs || (irhs && dynamic_cast<SimpleVector*>(irhsStoch->children[it]->vec)->elements()[i] != 0.0) )
            out<<" "<< rt<<" row_"<<rt<<"_"<<it <<"_"<<i <<endl;
      }
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

void StochGenMatrix::initTransposed(bool dynamic)
{
   Bmat->initTransposed(dynamic);
   Blmat->initTransposed(dynamic);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->initTransposedChild(dynamic);
}

void StochGenMatrix::deleteTransposed()
{
   Amat->deleteTransposed();
   Bmat->deleteTransposed();
   Blmat->deleteTransposed();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->deleteTransposed();
}

void StochGenMatrix::initTransposedChild(bool dynamic)
{
   Amat->initTransposed(dynamic);
   Bmat->initTransposed(dynamic);

   if( Blmat != NULL )
      Blmat->initTransposed(dynamic);
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

  nnz += Amat->numberOfNonZeros() + Bmat->numberOfNonZeros() + Blmat->numberOfNonZeros();

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


void StochGenMatrix::getNnzPerRow(OoqpVector& nnzVec, OoqpVector* linkParent)
{
   StochVector& nnzVecStoch = dynamic_cast<StochVector&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   SimpleVector* nnzvecl = NULL;

   Bmat->addNnzPerRow(*(nnzVecStoch.vec));

   if( linkParent != NULL )
      Amat->addNnzPerRow(*(nnzVecStoch.vec));

   /* with linking constraints? */
   if( nnzVecStoch.vecl || linkParent )
   {
      assert(nnzVecStoch.vecl == NULL || linkParent == NULL);

      if( linkParent )
         nnzvecl = dynamic_cast<SimpleVector*>(linkParent);
      else
         nnzvecl = dynamic_cast<SimpleVector*>(nnzVecStoch.vecl);

      if( linkParent != NULL || iAmSpecial(iAmDistrib, mpiComm) )
         Blmat->addNnzPerRow(*nnzvecl);
   }


   for( size_t it = 0; it < children.size(); it++ )
     children[it]->getNnzPerRow(*(nnzVecStoch.children[it]), nnzvecl);

   // distributed, with linking constraints, and at root?
   if( iAmDistrib && nnzVecStoch.vecl != NULL && linkParent == NULL )
   {
      // sum up linking constraints vectors
      const int locn = nnzvecl->length();
      double* buffer = new double[locn];

      MPI_Allreduce(nnzvecl->elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      nnzvecl->copyFromArray(buffer);

      delete[] buffer;
   }
}

void StochGenMatrix::getNnzPerCol(OoqpVector& nnzVec, OoqpVector* linkParent)
{
   StochVector& nnzVecStoch = dynamic_cast<StochVector&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   SimpleVector* const vec = dynamic_cast<SimpleVector*>(nnzVecStoch.vec);

   if( iAmSpecial(iAmDistrib, mpiComm) || linkParent != NULL )
   {
      Bmat->addNnzPerCol(*(vec));

      int blm, bln;
      Blmat->getSize(blm, bln);

      /* with linking constraints? */
      if( blm > 0 )
         Blmat->addNnzPerCol(*vec);
   }

   // not at root?
   if( linkParent != NULL )
      Amat->addNnzPerCol(*linkParent);
   else
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->getNnzPerCol(*(nnzVecStoch.children[it]), vec);
   }

   // distributed and at root?
   if( iAmDistrib && linkParent == NULL )
   {
      const int locn = vec->length();
      double* const entries = vec->elements();
      double* buffer = new double[locn];

      MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      vec->copyFromArray(buffer);

      delete[] buffer;
   }
}

void StochGenMatrix::getRowMinMaxVec(bool getMin, bool initializeVec,
      const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent)
{
   StochVector& minmaxVecStoch = dynamic_cast<StochVector&>(minmaxVec);
   const StochVector* const colScaleVecStoch = dynamic_cast<const StochVector*>(colScaleVec);

   SimpleVector* mvecl = NULL;
   const SimpleVector* const covecparent = dynamic_cast<const SimpleVector*>(colScaleParent);
   const SimpleVector* const covec = colScaleVecStoch != NULL ?
         dynamic_cast<SimpleVector*>(colScaleVecStoch->vec) : NULL;

   // assert tree compatibility
   assert(minmaxVecStoch.children.size() == children.size());

   Bmat->getRowMinMaxVec(getMin, initializeVec, covec, *(minmaxVecStoch.vec));

   // not at root?
   if( linkParent != NULL )
      Amat->getRowMinMaxVec(getMin, false, covecparent, *(minmaxVecStoch.vec));

   /* with linking constraints? */
   if( minmaxVecStoch.vecl || linkParent )
   {
      assert(minmaxVecStoch.vecl == NULL || linkParent == NULL);

      // at root?
      if( linkParent == NULL )
      {
         mvecl = dynamic_cast<SimpleVector*>(minmaxVecStoch.vecl);

         if( initializeVec )
         {
            if( getMin )
               mvecl->setToConstant(std::numeric_limits<double>::max());
            else
               mvecl->setToZero();
         }
      }
      else
      {
         mvecl = dynamic_cast<SimpleVector*>(linkParent);
      }

      if( linkParent != NULL || iAmSpecial(iAmDistrib, mpiComm) )
         Blmat->getRowMinMaxVec(getMin, false, covec, *mvecl);
   }

   if( colScaleVec != NULL )
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->getRowMinMaxVec(getMin, initializeVec, colScaleVecStoch->children[it], covec,
               *(minmaxVecStoch.children[it]), mvecl);
   }
   else
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->getRowMinMaxVec(getMin, initializeVec, NULL, NULL,
               *(minmaxVecStoch.children[it]), mvecl);
   }

   // distributed, with linking constraints, and at root?
   if( iAmDistrib && minmaxVecStoch.vecl != NULL && linkParent == NULL )
   {
      assert(mvecl != NULL);

      // sum up linking constraints vectors
      const int locn = mvecl->length();
      double* buffer = new double[locn];

      if( getMin )
         MPI_Allreduce(mvecl->elements(), buffer, locn, MPI_DOUBLE, MPI_MIN, mpiComm);
      else
         MPI_Allreduce(mvecl->elements(), buffer, locn, MPI_DOUBLE, MPI_MAX, mpiComm);

      mvecl->copyFromArray(buffer);

      delete[] buffer;
   }
}

void StochGenMatrix::getColMinMaxVec(bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleLink, OoqpVector& minmaxVec, OoqpVector* minmaxParent)
{
   StochVector& minmaxVecStoch = dynamic_cast<StochVector&>(minmaxVec);
   const StochVector* rowScaleVecStoch = dynamic_cast<const StochVector*>(rowScaleVec);

   // assert tree compatibility
   assert(minmaxVecStoch.children.size() == children.size());

   SimpleVector* const mvec = dynamic_cast<SimpleVector*>(minmaxVecStoch.vec);
   const SimpleVector* const covec = rowScaleVecStoch != NULL ?
         dynamic_cast<SimpleVector*>(rowScaleVecStoch->vec) : NULL;

   Bmat->getColMinMaxVec(getMin, initializeVec, covec, *(mvec));

   int blm, bln;
   Blmat->getSize(blm, bln);

   const SimpleVector* covecl = NULL;

   /* with linking constraints? */
   if( blm > 0 )
   {
      // with rowScale vector?
      if( rowScaleVecStoch != NULL )
      {
         // at root?
         if( minmaxParent == NULL )
            covecl = dynamic_cast<SimpleVector*>(rowScaleVecStoch->vecl);
         else
            covecl = dynamic_cast<const SimpleVector*>(rowScaleLink);

         assert(covecl != NULL);
      }

      if( iAmSpecial(iAmDistrib, mpiComm) || minmaxParent != NULL )
         Blmat->getColMinMaxVec(getMin, false, covecl, *mvec);
   }

   // not at root?
   if( minmaxParent != NULL )
      Amat->getColMinMaxVec(getMin, false, covec, *(minmaxParent));
   else
   {
      if( rowScaleVecStoch )
      {
         for( size_t it = 0; it < children.size(); it++ )
            children[it]->getColMinMaxVec(getMin, initializeVec, rowScaleVecStoch->children[it], covecl, *(minmaxVecStoch.children[it]), mvec);
      }
      else
      {
         for( size_t it = 0; it < children.size(); it++ )
            children[it]->getColMinMaxVec(getMin, initializeVec, NULL, NULL, *(minmaxVecStoch.children[it]), mvec);
      }
   }

   // distributed and at root?
   if( iAmDistrib && minmaxParent == NULL )
   {
      const int locn = mvec->length();
      double* const entries = mvec->elements();
      double* buffer = new double[locn];

      if( getMin )
         MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_MIN, mpiComm);
      else
         MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_MAX, mpiComm);

      mvec->copyFromArray(buffer);

      delete[] buffer;
   }
}


void StochGenMatrix::addRowSums( OoqpVector& sumVec, OoqpVector* linkParent )
{
   StochVector& sumVecStoch = dynamic_cast<StochVector&>(sumVec);
   SimpleVector* mvecl = NULL;

   // assert tree compatibility
   assert(sumVecStoch.children.size() == children.size());

   Bmat->addRowSums(*sumVecStoch.vec);

   // not at root?
   if( linkParent != NULL )
      Amat->addRowSums(*sumVecStoch.vec);

   /* with linking constraints? */
   if( sumVecStoch.vecl || linkParent )
   {
      assert(sumVecStoch.vecl == NULL || linkParent == NULL);

      // at root?
      if( linkParent == NULL )
         mvecl = dynamic_cast<SimpleVector*>(sumVecStoch.vecl);
      else
         mvecl = dynamic_cast<SimpleVector*>(linkParent);

      if( linkParent != NULL || iAmSpecial(iAmDistrib, mpiComm) )
         Blmat->addRowSums(*mvecl);
   }

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->addRowSums(*(sumVecStoch.children[it]), mvecl);

   // distributed, with linking constraints, and at root?
   if( iAmDistrib && sumVecStoch.vecl != NULL && linkParent == NULL )
   {
      assert(mvecl != NULL);

      // sum up linking constraints vectors
      const int locn = mvecl->length();
      double* buffer = new double[locn];

      MPI_Allreduce(mvecl->elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      mvecl->copyFromArray(buffer);

      delete[] buffer;
   }
}

void StochGenMatrix::addColSums( OoqpVector& sumVec, OoqpVector* linkParent )
{
   StochVector& sumVecStoch = dynamic_cast<StochVector&>(sumVec);

   // assert tree compatibility
   assert(sumVecStoch.children.size() == children.size());

   SimpleVector* const mvec = dynamic_cast<SimpleVector*>(sumVecStoch.vec);

   if( iAmSpecial(iAmDistrib, mpiComm) || linkParent != NULL )
      Bmat->addColSums(*mvec);

   int blm, bln;
   Blmat->getSize(blm, bln);

   /* with linking constraints? */
   if( blm > 0 && (iAmSpecial(iAmDistrib, mpiComm) || linkParent != NULL) )
      Blmat->addColSums(*mvec);

   // not at root?
   if( linkParent != NULL )
      Amat->addColSums(*linkParent);
   else
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->addColSums(*(sumVecStoch.children[it]), mvec);
   }

   // distributed and at root?
   if( iAmDistrib && linkParent == NULL )
   {
      const int locn = mvec->length();
      double* const entries = mvec->elements();
      double* buffer = new double[locn];

      MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      mvec->copyFromArray(buffer);

      delete[] buffer;
   }
}

void StochGenMatrix::initStaticStorageFromDynamic(const OoqpVector& rowNnzVec, const OoqpVector& colNnzVec, const OoqpVector* rowLinkVec, const OoqpVector* colParentVec)
{
   const StochVector& rowNnzVecStoch = dynamic_cast<const StochVector&>(rowNnzVec);
   const StochVector& colNnzVecStoch = dynamic_cast<const StochVector&>(colNnzVec);

   assert(rowNnzVecStoch.children.size() == colNnzVecStoch.children.size());

   const SimpleVector* const rowvec = dynamic_cast<const SimpleVector*>(rowNnzVecStoch.vec);
   const SimpleVector* const colvec = dynamic_cast<const SimpleVector*>(colNnzVecStoch.vec);

   const SimpleVector* const rowlink = dynamic_cast<const SimpleVector*>(rowNnzVecStoch.vecl);

   Amat->initStaticStorageFromDynamic(*rowvec, colParentVec); // initialized with colVec == NULL for parent
   Bmat->initStaticStorageFromDynamic(*rowvec, colvec);

   // at root?
   if( colParentVec == NULL )
   {
      assert(rowLinkVec == NULL);

      if( rowlink != NULL )
         Blmat->initStaticStorageFromDynamic(*rowlink, colvec);

      for( size_t it = 0; it < children.size(); it++ )
         children[it]->initStaticStorageFromDynamic(*(rowNnzVecStoch.children[it]), *(colNnzVecStoch.children[it]), rowlink, colvec);
   }
   else
   {
      if( rowLinkVec != NULL)
         Blmat->initStaticStorageFromDynamic(*rowLinkVec, colvec);
   }

}

void StochGenMatrix::freeDynamicStorage()
{
   Amat->freeDynamicStorage();
   Bmat->freeDynamicStorage();
   Blmat->freeDynamicStorage();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->freeDynamicStorage();
}

void StochGenMatrix::updateKLinkConsCount(std::vector<int>& linkCount) const
{
   if( Blmat == NULL )
      return;

   int m, n;

   Blmat->getSize(m, n);
   assert(m > 0);
   assert(linkCount.size() == size_t(m));

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         assert(children[it]->Blmat);
         children[it]->Blmat->updateNonEmptyRowsCount(linkCount);
      }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &linkCount[0], m, MPI_INT, MPI_SUM, mpiComm);
}

void StochGenMatrix::updateKLinkVarsCount(std::vector<int>& linkCount) const
{
   int m, n;

   Bmat->getSize(m, n);

   if( n == 0 )
      return;

   assert(linkCount.size() == size_t(n));

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         children[it]->Amat->getTranspose().updateNonEmptyRowsCount(linkCount);
         children[it]->Amat->deleteTransposed();
      }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &linkCount[0], n, MPI_INT, MPI_SUM, mpiComm);
}

std::vector<int> StochGenMatrix::get2LinkStartBlocks() const
{
   if( Blmat == NULL )
      return std::vector<int>();

   int m, n;
   Blmat->getSize(m, n);

   if( m == 0 )
      return std::vector<int>();

   assert(m > 0);

   std::vector<int> linkCount(m, 0);
   std::vector<int> linkBlockStart(m, -1);
   std::vector<int> linkBlockEnd(m, -1);
   std::vector<bool> is2link(m, false);

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         assert(children[it]->Blmat);
         children[it]->Blmat->updateNonEmptyRowsCount(it, linkCount, linkBlockStart, linkBlockEnd);
      }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &linkCount[0], m, MPI_INT, MPI_SUM, mpiComm);

   // set block identifier
   for( int i = 0; i < m; i++ )
   {
      assert(linkBlockEnd[i] == -1 || linkBlockStart[i] <= linkBlockEnd[i]);

      if( linkCount[i] == 2 && (linkBlockEnd[i] - linkBlockStart[i]) == 1 )
      {
         assert(linkBlockStart[i] >= 0 && linkBlockEnd[i] >= 0);
         is2link[i] = true;
      }
   }

   if( iAmDistrib )
   {
      // find 2-links between different processes

      int myrank, size;
      MPI_Comm_rank(mpiComm, &myrank);
      MPI_Comm_size(mpiComm, &size);

      assert(size > 0);

      // 1. allgather number of local 2-link candidates
      std::vector<int> localCandsIdx;
      std::vector<int> localCandsBlock;
      std::vector<int> candsPerProc(size, -1);

      for( int i = 0; i < m; i++ )
         if( linkCount[i] == 2 && linkBlockStart[i] >= 0 && linkBlockEnd[i] == -1 )
         {
            assert(!is2link[i]);

            localCandsIdx.push_back(i);
            assert(unsigned(linkBlockStart[i]) < children.size());

            localCandsBlock.push_back(linkBlockStart[i]);
         }

      int localcount = localCandsIdx.size();

      MPI_Allgather(&localcount, 1, MPI_INT, &candsPerProc[0], 1, MPI_INT, mpiComm);

#ifndef NDEBUG
      for( size_t i = 0; i < candsPerProc.size(); i++ )
         assert(candsPerProc[i] >= 0);
#endif

      // 2. allgatherv 2-link candidates
      std::vector<int> displacements(size + 1, 0);
      for( int i = 1; i < size; i++ )
         displacements[i] = candsPerProc[i - 1] + displacements[i - 1];

      const int nAllCands = displacements[size - 1] + candsPerProc[size - 1];

      displacements[size] = nAllCands;
      std::vector<int> allCandsRow(nAllCands, -1);
      std::vector<int> allCandsBlock(nAllCands, -1);

      MPI_Allgatherv(&localCandsIdx[0], localcount, MPI_INT, &allCandsRow[0], &candsPerProc[0],
            &displacements[0], MPI_INT, mpiComm);

      MPI_Allgatherv(&localCandsBlock[0], localcount, MPI_INT, &allCandsBlock[0], &candsPerProc[0],
            &displacements[0], MPI_INT, mpiComm);

#ifndef NDEBUG
      for( size_t i = 0; i < allCandsRow.size(); i++ )
         assert(allCandsRow[i] >= 0 && allCandsRow[i] < m && allCandsBlock[i] >= 0 && allCandsBlock[i] < int(children.size()));
#endif


      // 3. check which candidates are indeed 2-links
      std::vector<int> blocksHash(m, -1);
      for( int i = 0; i < size - 1; i++ )
      {
         // hash
         for( int j = displacements[i]; j < displacements[i + 1]; j++ )
            blocksHash[allCandsRow[j]] = allCandsBlock[j];

         // compare with next
         for( int j = displacements[i + 1]; j < displacements[i + 2]; j++ )
         {
            assert(allCandsBlock[j] > 0);
            const int candRow = allCandsRow[j];
            if( blocksHash[candRow] >= 0 )
            {
               assert(blocksHash[candRow] != allCandsBlock[j]);

               const int startBlock = std::min(blocksHash[candRow], allCandsBlock[j]);
               const int endBlock = std::max(blocksHash[candRow], allCandsBlock[j]);

               assert(startBlock >= 0 && endBlock >= 0);

               if( endBlock - startBlock != 1 )
                  continue;

               is2link[candRow] = true;

               // start block owned by this MPI process?
               if( !children[startBlock]->isKindOf(kStochGenDummyMatrix) )
               {
                  assert(children[endBlock]->isKindOf(kStochGenDummyMatrix));
                  linkBlockStart[candRow] = startBlock;
               }
               else
               {
                  assert(children[startBlock]->isKindOf(kStochGenDummyMatrix));
                  linkBlockStart[candRow] = -1;
               }
            }
         }

         // un-hash
         for( int j = displacements[i]; j < displacements[i + 1]; j++ )
            blocksHash[allCandsRow[j]] = -1;
      }
   }

   // correct block identifier
   for( int i = 0; i < m; i++ )
      if( !is2link[i] )
         linkBlockStart[i] = -1;

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &linkBlockStart[0], m, MPI_INT, MPI_MAX, mpiComm);

   return linkBlockStart;
}


void StochGenMatrix::permuteLinkingVars(const std::vector<unsigned int>& permvec)
{
   if( Blmat )
      Blmat->permuteCols(permvec);

   Bmat->permuteCols(permvec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->permuteLinkingVarsChild(permvec);
}

void StochGenMatrix::permuteLinkingVarsChild(const std::vector<unsigned int>& permvec)
{
   Amat->permuteCols(permvec);

   assert(children.size() == 0);
}

void StochGenMatrix::permuteLinkingCons(const std::vector<unsigned int>& permvec)
{
   if( Blmat )
      Blmat->permuteRows(permvec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->permuteLinkingCons(permvec);
}


void StochGenMatrix::updateTransposed()
{
  Amat->updateTransposed();
  Bmat->updateTransposed();
  Blmat->updateTransposed();

  for( size_t it = 0; it < children.size(); it++ )
     children[it]->updateTransposed();
}

