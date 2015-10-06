/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "SparseGenMatrix.h"
#include "SparseStorage.h"
#include <cassert>
#include "SimpleVector.h"

#include "DoubleMatrixTypes.h"

#include "SparseSymMatrix.h"


void SparseGenMatrix::printMatrixInMatlab(char *name)
{
  int *krowM = this->krowM();
  int *jcolM = this->jcolM();
  double *M = this->M();
  int nnz = this->numberOfNonZeros();
 
  int nn,mm;
  this->getSize(mm,nn);

   FILE* outfile = fopen( name, "wr" );
   
		fprintf(outfile,"n_dim_%s = %d;\n\n", name,nn); 
		fprintf(outfile,"m_dim_%s = %d;\n\n", name,mm); 
		fprintf(outfile,"nnz_%s = %d;\n\n", name,nnz); 
		
		int findkk=0;
		fprintf(outfile,"rowId_%s = [ ",name);
		for(int ii=0;ii<mm;ii++)
		  for(int kk=krowM[ii];kk<krowM[ii+1];kk++)
			if(findkk%10==0){
			 fprintf(outfile,"%d, ... \n",ii+1);
			 findkk++;
			}
			else{
			 fprintf(outfile,"%d,", ii+1);
			 findkk++;
			}
		fprintf(outfile,"]; \n\n");
	
		findkk=0;
		fprintf(outfile,"colId_%s = [ ",name);
		for(int kk=0;kk<nnz;kk++)
			if(findkk%10==0){
			 fprintf(outfile,"%d, ... \n", jcolM[kk]+1);
					 findkk++;
			}
			else{
			 fprintf(outfile,"%d,", jcolM[kk]+1);
					 findkk++;
			}
		fprintf(outfile,"]; \n\n");
	
		findkk=0;
		fprintf(outfile,"elts_%s = [ ",name);
		for(int kk=0;kk<nnz;kk++)
			if(findkk%10==0){
			  fprintf(outfile,"%5.17g, ... \n", M[kk]);
					 findkk++;
			}
			else{
			  fprintf(outfile,"%5.17g,", M[kk]);
					 findkk++;
			}
		fprintf(outfile,"]; \n\n"); 


  fclose(outfile);

}






int SparseGenMatrix::isKindOf( int type )
{
  return type == kSparseGenMatrix || type == kGenMatrix;
}

SparseGenMatrix::SparseGenMatrix()
    : m_Mt(NULL)
{
}

SparseGenMatrix::SparseGenMatrix( int rows, int cols, int nnz )
  : m_Mt(NULL)
{
  mStorage = SparseStorageHandle( new SparseStorage( rows, cols, nnz ) );
}


SparseGenMatrix::SparseGenMatrix( int rows, int cols, int nnz,
				  int krowM[], int jcolM[],  double M[],
				  int deleteElts)
  : m_Mt(NULL)
{
  //cout << "SparseGenMatrix1  " << rows << " " << cols << " " << nnz << endl;
  mStorage = SparseStorageHandle( new SparseStorage( rows, cols,nnz, 
						     krowM, jcolM, M,
						     deleteElts) );
}

/*SparseGenMatrix::SparseGenMatrix(const std::vector<SparseGenMatrix*> &blocks, bool diagonal)
  : m_Mt(NULL)
{
  vector<SparseStorage*> v(blocks.size());
  for (size_t i = 0; i < blocks.size(); i++) v[i] = blocks[i]->mStorage;
  mStorage = SparseStorageHandle(new SparseStorage(v,diagonal));
}
*/

SparseGenMatrix::~SparseGenMatrix()
{
//  cout << "~~~~~~~~~SparseGenMatrix " << mStorage->refs()  << endl;
  if(m_Mt) delete m_Mt;
  
}


void SparseGenMatrix::atPutDense( int row, int col, double * A, int lda,
				      int rowExtent, int colExtent )
{
  mStorage->atPutDense( row, col, A, lda, rowExtent, colExtent );
}


void SparseGenMatrix::fromGetDense( int row, int col, double * A, int lda,
					int rowExtent, int colExtent )
{
  mStorage->fromGetDense( row, col, A, lda, rowExtent, colExtent );
}
  

void SparseGenMatrix::fromGetSpRow( int row, int col,
				    double A[], int lenA,
				    int jcolA[], int& nnz,
				    int colExtent, int& info )
{
  mStorage->fromGetSpRow( row, col, A, lenA, jcolA, nnz,
			  colExtent, info );
}


void SparseGenMatrix::putSparseTriple( int irow[], int len,
					   int jcol[], double A[], 
					   int& info )
{
  mStorage->putSparseTriple( irow, len, jcol, A, info );
}


void SparseGenMatrix::writeToStream(ostream& out) const
{
  mStorage->writeToStream( out );
}


void SparseGenMatrix::randomize( double alpha, double beta, double * seed )
{
  mStorage->randomize( alpha, beta, seed );
}


void SparseGenMatrix::getDiagonal( OoqpVector& vec )
{
  mStorage->getDiagonal( vec );
}


void SparseGenMatrix::setToDiagonal( OoqpVector& vec )
{
  mStorage->setToDiagonal( vec );
}


void SparseGenMatrix::atPutSpRow( int row, double A[],
				      int lenA, int jcolA[], int& info )
{
  mStorage->atPutSpRow( row, A, lenA, jcolA, info );
}


int SparseGenMatrix::numberOfNonZeros()
{
  return mStorage->numberOfNonZeros();
}


void SparseGenMatrix::symmetrize( int& info ) 
{
  mStorage->symmetrize( info );
}

int* SparseGenMatrix::symmetrize_set( int& info) 
{
  return mStorage->symmetrize_set( info);
}

void SparseGenMatrix::symmetrize_valonly( double *val_lower,int *goffIDX) 
{
  mStorage->symmetrize_valonly( val_lower,goffIDX);
}




void SparseGenMatrix::getSize( long long& m, long long& n )
{
  m = mStorage->m;
  n = mStorage->n;
}
void SparseGenMatrix::getSize( int& m, int& n )
{
  m = mStorage->m;
  n = mStorage->n;
}


void SparseGenMatrix::atPutSubmatrix( int destRow, int destCol,
					  DoubleMatrix& M,
					  int srcRow, int srcCol,
					  int rowExtent, int colExtent )
{
  int i, k;
  int info, nnz;

  int *    ja = new int[colExtent];
  double * a  = new double[colExtent];

  nnz = 0;
  for ( i = 0; i < rowExtent; i++ ) {
    M.fromGetSpRow( srcRow + i, srcCol, a, colExtent, ja,
		     nnz, colExtent, info );
    for( k = 0; k < nnz; k++ ) {
      ja[k] += (destCol - srcCol);
    }
    mStorage->atPutSpRow( destRow + i, a, nnz, ja, info );
  }

  delete [] ja;
  delete [] a;
}


void SparseGenMatrix::mult ( double beta,  OoqpVector& y_in,
				 double alpha, OoqpVector& x_in )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert( x.n == mStorage->n && y.n == mStorage->m );

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->mult( beta, yv, 1, alpha, xv, 1 );
}

void SparseGenMatrix::mult ( double beta,  double y[], int incy,
			     double alpha, double x[], int incx )
{
  mStorage->mult( beta, y, incy, alpha, x, incx);
}


void SparseGenMatrix::transMult ( double beta,   OoqpVector& y_in,
				  double alpha,  OoqpVector& x_in )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);

  assert( x.n == mStorage->m && y.n == mStorage->n );

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->transMult( beta, yv, 1, alpha, xv, 1 );
}

// wrapper added by cpetra
void SparseGenMatrix::transMult( double beta,  OoqpVector& y_in, int incy,
			         double alpha, OoqpVector& x_in, int incx )
{
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  assert(x.n>0 && y.n>0);
  assert(x.n>=incx*mStorage->m);
  assert(y.n>=incy*mStorage->n);

  double *xv = 0, *yv = 0;

  if( x.n > 0 ) xv = &x[0];
  if( y.n > 0 ) yv = &y[0];

  mStorage->transMult( beta, yv, incy, alpha, xv, incx );
}

// wrapper added by cpetra
void SparseGenMatrix::transMult( double beta,  double yv[], int incy,
				 double alpha, double xv[], int incx )
{
  mStorage->transMult( beta, yv, incy, alpha, xv, incx );
}


double SparseGenMatrix::abmaxnorm()
{
  return mStorage->abmaxnorm();
}


void SparseGenMatrix::atPutDiagonal( int idiag, OoqpVector& vvec )
{
  SimpleVector & v = dynamic_cast<SimpleVector &>(vvec);

  mStorage->atPutDiagonal( idiag, &v[0], 1, v.n );
}


void SparseGenMatrix::fromGetDiagonal( int idiag, OoqpVector& vvec )
{
  mStorage->fromGetDiagonal( idiag, vvec );
}

void SparseGenMatrix::ColumnScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );
}

void SparseGenMatrix::SymmetricScale( OoqpVector& vec )
{
  mStorage->SymmetricScale( vec );
}

void SparseGenMatrix::RowScale( OoqpVector& vec )
{
  mStorage->ColumnScale( vec );
}

void SparseGenMatrix::scalarMult( double num )
{
  mStorage->scalarMult( num );
}

void SparseGenMatrix::matTransDMultMat(OoqpVector& d_, SymMatrix** res)
{
  SimpleVector& d = dynamic_cast<SimpleVector &>(d_);

  int m=mStorage->m; int n=mStorage->n; int nnz=mStorage->numberOfNonZeros();

  if(*res==NULL) {
    assert(m_Mt==NULL);
    //we need to form the transpose
    m_Mt=new SparseGenMatrix(n,m,nnz);
    mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());

    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=NULL; int* jcolMtM=NULL; double* dMtM=NULL;
    mStorage->matTransDSymbMultMat(&d[0], 
				   m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
				   &krowMtM, &jcolMtM, &dMtM);

    *res = new SparseSymMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
  }

  assert(res); 
  assert(m_Mt);

  SparseSymMatrix* MtDM = dynamic_cast<SparseSymMatrix*>(*res);

  mStorage->matTransDMultMat(&d[0], 
			     m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
			     MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}

void SparseGenMatrix::matTransDinvMultMat(OoqpVector& d_, SymMatrix** res)
{
  SimpleVector& d = dynamic_cast<SimpleVector &>(d_);

  int m=mStorage->m; int n=mStorage->n; int nnz=mStorage->numberOfNonZeros();

  if(*res==NULL) {
    assert(m_Mt==NULL);
    //we need to form the transpose
    m_Mt=new SparseGenMatrix(n,m,nnz);
    mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());

    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=NULL; int* jcolMtM=NULL; double* dMtM=NULL;
    mStorage->matTransDSymbMultMat(&d[0], 
				   m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
				   &krowMtM, &jcolMtM, &dMtM);

    *res = new SparseSymMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
  }

  assert(res); 
  assert(m_Mt);

  SparseSymMatrix* MtDM = dynamic_cast<SparseSymMatrix*>(*res);

  mStorage->matTransDinvMultMat(&d[0], 
			     m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(),
			     MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}
void SparseGenMatrix::matMultTrans(SymMatrix** res)
{
  int m=mStorage->m; int n=mStorage->n; 
  int nnz=mStorage->numberOfNonZeros();

  SimpleVector d(n); d.setToConstant(1.0);

  if(*res==NULL) {
    //assert(m_Mt==NULL);

    //we need to form the transpose
    if(!m_Mt) {
      m_Mt = new SparseGenMatrix(n,m,nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
    }
    //find the sparsity pattern of the product -> the buffers for result will be allocated
    int* krowMtM=NULL; int* jcolMtM=NULL; double* dMtM=NULL;

    m_Mt->mStorage->matTransDSymbMultMat(&d[0], 
					 krowM(), jcolM(), M(),
					 &krowMtM, &jcolMtM, &dMtM);
    *res = new SparseSymMatrix(m, krowMtM[m], krowMtM, jcolMtM, dMtM, 1);
  }
  SparseSymMatrix* MMt = dynamic_cast<SparseSymMatrix*>(*res);
  m_Mt->mStorage->matTransDMultMat(&d[0],
				   krowM(), jcolM(), M(),
				   MMt->krowM(), MMt->jcolM(), MMt->M());
}

void SparseGenMatrix::copyMtxFromDouble(int copyLength, double *values)
{
  mStorage->copyMtxFromDouble(copyLength,values);
}

/* REMOVE_NY
void SparseGenMatrix::changeSlackPartQP(Variables * vars_, OoqpVector & vecS_,
			int nSlack_in, int nxL_in, int nxU_in, int nsL_in, int nsU_in,  int const caseFlag)
{
  SimpleVector& vSlack = dynamic_cast<SimpleVector &>(vecS_);
  double *valSlack = vSlack.elements();
  
  int nx = vSlack.n;
  int m = mStorage->m; 
  assert ( nx = mStorage->n); 
//  int nnz = mStorage->numberOfNonZeros();

  int i;
  if(caseFlag == 1)
    for(i= 0; i<nSlack_in; i++){
	  M()[krowM()[m-i]-1] = -1.0 * ((1<valSlack[nx-i-1])?1:valSlack[nx-i-1]);
    }
  else if (caseFlag == 0)
    for(i= 0; i<nSlack_in; i++){
	  M()[krowM()[m-i]-1] = -1.0;
    }
}
*/

void SparseGenMatrix::setAdditiveDiagonal(OoqpVector& v )
{
  mStorage->setAdditiveDiagonal( v );
}


  

void SparseGenMatrix::fromGetSpRow_WithRowStart( int row, int col,
				    double A[], int lenA,
				    int jcolA[], int& nnz,
				    int colExtent, int& info, int & rowStart)
{
  mStorage->fromGetSpRow_WithRowStart( row, col, A, lenA, jcolA, nnz,
			  colExtent, info, rowStart );
}

void SparseGenMatrix::fromGetDense_withMap( int row, int col, double * A, int lda,
					int rowExtent, int colExtent, int const FirstCall, std::map<int,int> &ValIdxMap )
{
  mStorage->fromGetDense_withMap( row, col, A, lda, rowExtent, colExtent, FirstCall, ValIdxMap);
}

