#include "SparseSymMatrixRowMajList.h"
#include <cassert>
#include <cmath>
#include "SimpleVector.h"

#include "DoubleMatrixTypes.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include <cstdlib>
using namespace std;



SparseSymMatrixRowMajList::SparseSymMatrixRowMajList( int size )
  : nnz(0)
{
  //mStorage = SparseStorageHandle( new SparseStorage(size, size, nnz) );
  vlmat.reserve(size);
  for(int i=0; i<size; i++) {
    std::list<ColVal> row; row.push_back(ColVal(i, 0.));
    vlmat.push_back(row);
  }
  nnz += size;
}


void SparseSymMatrixRowMajList::putSparseTriple( int irow[], int len,
					   int jcol[], double A[], 
					   int& info )
{
  assert(false && "not implemented"); //mStorage->putSparseTriple( irow, len, jcol, A, info );
}


void SparseSymMatrixRowMajList::fromGetDense( int row, int col, double * A, int lda,
					int rowExtent, int colExtent )
{
  assert(false && "not implemented"); //mStorage->fromGetDense( row, col, A, lda, rowExtent, colExtent );
}


void SparseSymMatrixRowMajList::getDiagonal( OoqpVector& vec )
{
  assert(false && "not implemented"); //mStorage->getDiagonal( vec );
}

void SparseSymMatrixRowMajList::setToDiagonal( OoqpVector& vec )
{
  assert(false && "not implemented");// mStorage->setToDiagonal( vec );
}

void SparseSymMatrixRowMajList::symAtPutSpRow( int row, double A[],
					 int lenA, int jcolA[],
					 int& info )
{
  assert(false && "not implemented"); //
  // // Lower triangular put
  // int lA = lenA;
  // while( lA > 0 && jcolA[lA - 1] > row ) lA--;
  // if( lA > 0 ) {
  //   mStorage->atPutSpRow( row, A, lA, jcolA, info );
  // } else {
  //   info = 0;
  // }
}

void SparseSymMatrixRowMajList::fromGetSpRow( int row, int col,
					double A[], int lenA,
					int jcolA[], int& nnz,
					int colExtent, int& info )
{
   assert(false && "not implemented"); // mStorage->fromGetSpRow( row, col, A, lenA, jcolA, nnz, colExtent, info );
}


void SparseSymMatrixRowMajList::randomizePSD(double * seed)
{
  assert(false && "not implemented"); //
}


void SparseSymMatrixRowMajList::symAtPutSubmatrix( int destRow, int destCol,
					     DoubleMatrix& M,
					     int srcRow, int srcCol,
					     int rowExtent, int colExtent )
{
  assert(false && "not implemented"); //

//   int i, k;
//   int info, nnz, nnzR;

// //  int *    ja = new int[colExtent];
// //  double * a = new double[colExtent];

//   nnz = 0;
//   for ( i = 0; i < rowExtent; i++ ) {

// 	nnzR	  = (M.krowM()[i+1])-(M.krowM()[i]);
// 	int *	 ja = new int[nnzR];
// 	double * a = new double[nnzR];
// 	M.fromGetSpRow( srcRow + i, srcCol, a, nnzR, ja, nnz, colExtent, info );
	
// //    M.fromGetSpRow( srcRow + i, srcCol, a, colExtent, ja, nnz, colExtent, info );
//     for( k = 0; k < nnz; k++ ) {
//       ja[k] += (destCol - srcCol);
//     }
//     this->symAtPutSpRow( destRow + i, a, nnz, ja, info );
//     assert( info == 0 );

// 	assert(nnzR==nnz);
//     delete [] a;
//     delete [] ja;	
//   }

// //  delete [] a;
// //  delete [] ja;
}

// Pass these to storage
void SparseSymMatrixRowMajList::getSize( long long& m, long long& n )
{
  m = n = vlmat.size();
  // int mint, nint;
  // mStorage->getSize( mint, nint );
  // m=mint; n=nint;
}

void SparseSymMatrixRowMajList::getSize( int& m, int& n )
{
  m = n = vlmat.size();
}


long long SparseSymMatrixRowMajList::size()
{
  return vlmat.size();
}
void SparseSymMatrixRowMajList::mult ( double beta,  OoqpVector& y_in,
				 double alpha, OoqpVector& x_in )
{
  assert(false && "not implemented"); //

  // SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  // SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  // assert( x.n == mStorage->n &&  y.n == mStorage->m );
  
  // double *xv = 0, *yv = 0;
  // if( x.n > 0 ) xv = &x[0];
  // if( y.n > 0 ) yv = &y[0];

  // this->mult( beta, yv, 1, alpha, xv, 1 );
}

void SparseSymMatrixRowMajList::transMult ( double beta,   OoqpVector& y_in,
				      double alpha,  OoqpVector& x_in )
{
  assert(false && "not implemented"); //
  // SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  // SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  
  // assert( x.n == mStorage->n &&  y.n == mStorage->m );
  
  // double *xv = 0, *yv = 0;
  // if( x.n > 0 ) xv = &x[0];
  // if( y.n > 0 ) yv = &y[0];

  // this->mult( beta, yv, 1, alpha, xv, 1 );
}
  
void SparseSymMatrixRowMajList::transMult ( double beta,  double y[], int incy,
				      double alpha, double x[], int incx )
{
  assert(false && "not implemented"); //this->mult( beta, y, incy, alpha, x, incx );
}

double SparseSymMatrixRowMajList::abmaxnorm()
{
  assert(false && "not implemented"); 
  return 1.;
  //return mStorage->abmaxnorm();
}

void SparseSymMatrixRowMajList::writeToStream(ostream& out) const
{
  assert(false && "not implemented"); //mStorage->writeToStream( out );
}

void SparseSymMatrixRowMajList::mult ( double beta,  double y[], int incy,
				 double alpha, double x[], int incx )
{
  assert(false && "not implemented"); //

  // int m, n, i, j, k;
  // this->getSize(m, n); 

  // int * jcolM = mStorage->jcolM;
  // int * krowM = mStorage->krowM;
  // double * M  = mStorage->M;

  // for ( i = 0; i < m; i++ ) {
  //   y[i * incy] *= beta;
  // }
  // for ( i = 0; i < n; i++ ) {
  //   for( k = krowM[i]; k < krowM[i+1]; k++ ) {
  //     j = jcolM[k];
	  
  //     y[i * incy] += alpha * M[k] * x[j * incx];
  //     if ( i != j ) {
  // 	y[j * incy] += alpha * M[k] * x[i * incx];
  //     }
  //   }
  // }
}


void SparseSymMatrixRowMajList::atPutDiagonal( int idiag, OoqpVector& v_ )
{
  assert(idiag<vlmat.size());
  SimpleVector& v = dynamic_cast<SimpleVector&>(v_);
  for(int i=idiag; i<idiag+v.n; i++) {
    if(vlmat[i].empty()) {
      vlmat[i].push_back(ColVal(i,v[i-idiag]));
      nnz++;
    } else {
      assert(vlmat[i].front().jcol>=i); //see atAddDiagonal for how to properly implement this case
      vlmat[i].front().M = v[i-idiag];
    }
  }
}

void SparseSymMatrixRowMajList::atAddDiagonal( int idiag, OoqpVector& v_ )
{
  assert(idiag<vlmat.size());
  SimpleVector& v = dynamic_cast<SimpleVector&>(v_);
  for(int i=idiag; i<idiag+v.n; i++) {
    if(vlmat[i].empty()) {
      vlmat[i].push_back(ColVal(i,v[i-idiag]));
      nnz++;
    } else {
      assert(vlmat[i].back().jcol<=i);
      if(vlmat[i].back().jcol<i) {
	vlmat[i].push_back(ColVal(i,v[i-idiag]));
	nnz++;
      } else {
	vlmat[i].back().M += v[i-idiag];
      }
    }
  }
}


void SparseSymMatrixRowMajList::fromGetDiagonal( int idiag, OoqpVector& v )
{
  assert(false && "not implemented"); // mStorage->fromGetDiagonal( idiag, v );
}

void SparseSymMatrixRowMajList::SymmetricScale( OoqpVector& vec )
{
  assert(false && "not implemented"); //mStorage->SymmetricScale( vec );
}

void SparseSymMatrixRowMajList::ColumnScale( OoqpVector& vec )
{
  assert(false && "not implemented"); //mStorage->ColumnScale( vec );
}

void SparseSymMatrixRowMajList::RowScale( OoqpVector& vec )
{
  assert(false && "not implemented"); //mStorage->RowScale( vec );
}

void SparseSymMatrixRowMajList::scalarMult( double num )
{
  assert(false && "not implemented"); //mStorage->scalarMult( num );
}

void SparseSymMatrixRowMajList::reduceToLower()
{
  assert(false && "not implemented"); //mStorage->reduceToLower();
}

void SparseSymMatrixRowMajList::copyMtxFromDouble(int copyLength,double *values)
{
  assert(false && "not implemented"); //mStorage->copyMtxFromDouble(copyLength,values);
}

void SparseSymMatrixRowMajList::setAdditiveDiagonal(OoqpVector& v )
{
  assert(false && "not implemented"); //mStorage->setAdditiveDiagonal( v );
}

void SparseSymMatrixRowMajList::copyDiagonalVal_From( int idiag, OoqpVector& v, bool firstCall,std::map<int,int> &ValIdxMap )
{
  assert(false && "not implemented"); //mStorage->copyDiagonalVal_From( idiag, v, firstCall,ValIdxMap);
}


void SparseSymMatrixRowMajList::symAtPutSpRow_CorrectMap( int row, double A[],
					 int lenA, int jcolA[],
					 int& info, std::map<int,int> &ValIdxMap, int const constIDX)
{
  assert(false && "not implemented"); //
  // // Lower triangular put
  // int lA = lenA;
  // while( lA > 0 && jcolA[lA - 1] > row ) lA--;
  // if( lA > 0 ) {
  //   mStorage->atPutSpRow_CorrectMap( row, A, lA, jcolA, info,ValIdxMap, constIDX);
  // } else {
  //   info = 0;
  // }
}

void SparseSymMatrixRowMajList::fromGetSpRow_WithRowStart( int row, int col,
					double A[], int lenA,
					int jcolA[], int& nnz,
					int colExtent, int& info, int & rowStart )
{
  assert(false && "not implemented"); //mStorage->fromGetSpRow_WithRowStart( row, col, A, lenA, jcolA, nnz, colExtent, info, rowStart);
}


void SparseSymMatrixRowMajList::atAddSpRow(const int& row, int* jcolSrc, double* M, const int& nelems)
{
  int itSrc=0;
  list<ColVal>::iterator itDest = vlmat[row].begin();

  while(itSrc<nelems) {
    while(itDest!=vlmat[row].end() && itDest->jcol < jcolSrc[itSrc]) {
      ++itDest;
    }
    if(itDest==vlmat[row].end())
      break;

    if(itDest->jcol == jcolSrc[itSrc]) {
      if(itDest->jcol > row) {
	assert(false && "should not find this jcol here"); //this->putElem(itDest->jcol, row, M[itSrc]); 
      } else {
	itDest->M += M[itSrc];
      }
      itSrc++;
    } else {
      if(jcolSrc[itSrc]>row) {
	assert(false && "remove me when testing with the small example");
	//the addElem below is highly inefficient for large-problems and it is here
	//to ensure that the small test examples work despite the fact that they 
	//the upper triangular of Q.
	this->addElem(jcolSrc[itSrc], row, M[itSrc]);
      } else {
	vlmat[row].insert(itDest, ColVal(jcolSrc[itSrc], M[itSrc]));
	nnz++;
      }
      itSrc++;
    }
  }

  while(itSrc<nelems) {
    assert(itDest==vlmat[row].end());
    if(jcolSrc[itSrc]>row) {
      //assert(false && "remove me when testing with the small example");
      //the addElem below is highly inefficient for large-problems and it is here
      //to ensure that the small test examples work despite the fact that they 
      //the upper triangular of Q.
      this->addElem(jcolSrc[itSrc], row, M[itSrc]);
    } else {
      vlmat[row].push_back(ColVal(jcolSrc[itSrc], M[itSrc]));
      nnz++;
    }
    itSrc++;
  }
}
void SparseSymMatrixRowMajList::atAddSpRow(const int& row, std::list<ColVal>& colvalSrc)
{
  int extrannzdest = 0;
  list<ColVal>::iterator itSrc  = colvalSrc.begin();
  list<ColVal>::iterator itDest = vlmat[row].begin();
  while(itSrc!=colvalSrc.end()) {
    while(itDest!=vlmat[row].end() && itDest->jcol < itSrc->jcol) {
      ++itDest;
    }

    if(itDest==vlmat[row].end())
      break;

    if(itDest->jcol == itSrc->jcol) {
      assert(itDest->jcol <= row && "should not find this jcol here"); //this->putElem(itDest->jcol, row, M[itSrc]); 
      itDest->M += itSrc->M;
      ++itSrc;
    } else {
      vlmat[row].insert(itDest, *itSrc);
      //assert(false);
      extrannzdest++;
      ++itSrc;
    }
  }

  while(itSrc!=colvalSrc.end()) {
    assert(itDest==vlmat[row].end());
    assert(itSrc->jcol<=row && "lower triangle elements only !?!");

    vlmat[row].push_back(*itSrc);
    extrannzdest++;
    ++itSrc;
  }
  nnz += extrannzdest;
}

void SparseSymMatrixRowMajList::setToConstant(const double& val)
{
  for(int i=0; i<vlmat.size(); i++) {
    for(list<ColVal>::iterator it=vlmat[i].begin(); it!=vlmat[i].end(); ++it) {
      it->M=val;
    }
  }
}


/* Precondition: both colDest and jcolSrc are supposed to be ordered */
void SparseSymMatrixRowMajList::
mergeSetColumn(const int& rowDestIdx, const int& colDestIdxOffset,
	       int* jcolSrc, double* Msrc, const int& nElemsSrc, const int& colSrcIdxOffset,
	       const int& colExtent)
{
  std::list<ColVal>& colDest = vlmat[rowDestIdx];
  std::list<ColVal>::iterator itDest = colDest.begin();
  int itSrc=0;
  
  while(itSrc<nElemsSrc) {
#ifdef DEBUG
    if(itSrc+1<nElemsSrc)
      assert(jcolSrc[itSrc]<jcolSrc[itSrc+1]);
#endif
    //skip elements not in colSrcIdxOffset+0, colSrcIdxOffset+1,...,colSrcIdxOffset+colExtent-1
    if(jcolSrc[itSrc]-colSrcIdxOffset>=colExtent) {
      itSrc++;
      continue;
    }

    while(itDest!=colDest.end() && (itDest->jcol-colDestIdxOffset < jcolSrc[itSrc]-colSrcIdxOffset)) {
      ++itDest;
    }
    if(itDest==colDest.end())
      break;
    if(itDest->jcol-colDestIdxOffset == jcolSrc[itSrc]-colSrcIdxOffset) {
      assert(jcolSrc[itSrc]-colSrcIdxOffset+colDestIdxOffset<=rowDestIdx);
      itDest->M = Msrc[itSrc];
      itSrc++;
    } else {
      int colDestIdx = jcolSrc[itSrc]-colSrcIdxOffset+colDestIdxOffset;
      if(colDestIdx>rowDestIdx) {
	//if(!bSymUpdate) {
	//  //this code should not run given the conventions for storing the matrices in PIPS.
	//  assert(!bSymUpdate &&  "careful here: sym update not enforced and you are trying "
	//	 "to put entries in the upper triangular; do you know what you're doing? (2)");
	//  //nothing needs to be done if the symmetric update is not enforced. 
	//} else {
	assert(false && "remove me after testing with the small examples");
	//the addElem below is highly inefficient for large-problems and it is here
	//to ensure that the small test examples work despite the fact that they 
	//the upper triangular of Q. 
	
	//add the element in the upper triangular part, that is, instead of inserting at
	//(rowDestIdx, colDestIdx) we insert at (colDestIdx, rowDestIdx)
	this->putElem(colDestIdx,rowDestIdx,Msrc[itSrc]);
	itSrc++;
	//}
      } else {
	//this case is for itDest->jcol-colDestIdxOffset > jcolSrc[itSrc]-colSrcIdxOffset
	colDest.insert(itDest, ColVal(jcolSrc[itSrc]-colSrcIdxOffset+colDestIdxOffset,Msrc[itSrc]));
	nnz++;
	itSrc++;
      }
    }
  }

  //add the elements that have col indexes greater than anything in colDest list.
  while(itSrc<nElemsSrc) {
    if(jcolSrc[itSrc]-colSrcIdxOffset>=colExtent) {
      itSrc++;
      continue;
    }

    assert(itDest==colDest.end());
    int colDestIdx = jcolSrc[itSrc]-colSrcIdxOffset+colDestIdxOffset;
    if(colDestIdx>rowDestIdx) {
      //if(!bSymUpdate) {
	//this code should not run given the conventions for storing the matrices in PIPS.
	//assert(!bSymUpdate &&  "careful here: sym update not enforced and you are trying "
	//       "to put entries in the upper triangular; do you know what you're doing? (1)");
	//this->putElem(colDestIdx,rowDestIdx,Msrc[itSrc]);
      //nothing needs to be done if the symmetric update is not enforced. 
      //} else {
	//assert(false && "remove me after testing with the small examples");
      //the addElem below is highly inefficient for large-problems and it is here
      //to ensure that the small test examples work despite the fact that they provide
      //the upper triangular of Q.
      
      //add the element in the upper triangular part, that is, instead of inserting at
      //(rowDestIdx, colDestIdx) we insert at (colDestIdx, rowDestIdx)
	//assert(false);
	this->putElem(colDestIdx,rowDestIdx,Msrc[itSrc]);
	//}
    } else  {
      colDest.push_back(ColVal(jcolSrc[itSrc]-colSrcIdxOffset+colDestIdxOffset,Msrc[itSrc]));
      nnz++;
    }
    itSrc++;
  }
}
void SparseSymMatrixRowMajList::putElem(const int& row, const int& col, const double& M)
{
  std::list<ColVal>& colDest = vlmat[row];
  std::list<ColVal>::iterator itDest = colDest.begin();
  while(itDest!=colDest.end()) {
    if(itDest->jcol > col) {
      colDest.insert(itDest, ColVal(col, M));
      nnz++;
      return;
    } else { 
      if(itDest->jcol == col) {
	itDest->M = M;
	return;
      }
    }
    ++itDest;
  }
  //here itDest==colDest.end(), which means the element needs to be inserted at the end of the list
  colDest.push_back(ColVal(col, M));
  nnz++;
}

void SparseSymMatrixRowMajList::addElem(const int& row, const int& col, const double& M)
{
  std::list<ColVal>& colDest = vlmat[row];
  std::list<ColVal>::iterator itDest = colDest.begin();
  while(itDest!=colDest.end()) {
    if(itDest->jcol > col) {
      colDest.insert(itDest, ColVal(col, M));
      nnz++;
      return;
    } else { 
      if(itDest->jcol == col) {
	itDest->M += M;
	return;
      }
    }
    ++itDest;
  }
  //here itDest==colDest.end(), which means the element needs to be inserted at the end of the list
  colDest.push_back(ColVal(col, M));
  nnz++;
}

void SparseSymMatrixRowMajList::symAtSetSubmatrix( int destRow, int destCol, DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent)
{
  int i, k, colIdx, colCount, tmp;
  //bool srcIsSym = false;
  int* krowSrc = NULL;
  int* jcolSrc = NULL;
  double*  MSrc = NULL;


  SparseSymMatrix* pM = dynamic_cast<SparseSymMatrix*>(&M);
  if(pM==NULL) {
    SparseGenMatrix* pM = dynamic_cast<SparseGenMatrix*>(&M);
    if(pM==NULL) {
      assert(false && "this code expects a sparse matrix type");
      exit(-1);
    } else {
      krowSrc = pM->krowM(); jcolSrc = pM->jcolM(); MSrc = pM->M();
    }
  } else {
    //srcIsSym = true;
    krowSrc = pM->krowM(); jcolSrc = pM->jcolM(); MSrc = pM->M(); 
  }
 
  for ( i = 0; i < rowExtent; i++ ) {
    tmp = i+srcRow;
    colIdx = krowSrc[tmp];
    colCount = krowSrc[tmp+1]-colIdx;
    
    //if(srcIsSym) {
      //if(colCount>0) {
      //assert(destCol+jcolSrc[krowSrc[i+srcRow]]>=i+destRow && 
      //       "symmetric matrices need to have only elements in the upper triangle");
      //}
    //}

    mergeSetColumn(i+destRow, destCol,
		   jcolSrc+colIdx, MSrc+colIdx, colCount, srcCol,
		   colExtent);
  }
}

void SparseSymMatrixRowMajList::
symAtAddSubmatrix( int destRow, int destCol,
		   DoubleMatrix& Mat,
		   int srcRow, int srcCol,
		   int rowExtent, int colExtent )
{
  int i, k, colIdx, colCount, tmp;

  int* krowSrc = NULL;
  int* jcolSrc = NULL;
  double*  MSrc = NULL;

  assert(destCol==0 && "code not supporting/not tested with nonzero destCol");

  SparseSymMatrix* pM = dynamic_cast<SparseSymMatrix*>(&Mat);
  if(pM==NULL) {
    SparseGenMatrix* pM = dynamic_cast<SparseGenMatrix*>(&Mat);
    if(pM==NULL) {
      assert(false && "this code expects a sparse matrix type");
      exit(-1);
    } else {
      krowSrc = pM->krowM(); jcolSrc = pM->jcolM(); MSrc = pM->M();
    }
  } else {
    krowSrc = pM->krowM(); jcolSrc = pM->jcolM(); MSrc = pM->M(); 
  }
 
  for ( i = 0; i < rowExtent; i++ ) {
    tmp = i+srcRow;
    colIdx = krowSrc[tmp];
    colCount = krowSrc[tmp+1]-colIdx;

    //assert(jcolSrc[colIdx]>=i);

    this->atAddSpRow(i+destRow, jcolSrc+colIdx, MSrc+colIdx, colCount);
  }
}

void SparseSymMatrixRowMajList::atGetSparseTriplet(int* ii, int* jj, double* MM, bool fortran/*=true*/)
{
  int o=fortran?1:0;
  int itnz=0;
  for(int i=0; i<vlmat.size(); i++) {
    for(list<ColVal>::const_iterator it=vlmat[i].begin(); it!=vlmat[i].end(); ++it) {
      //printf("i=%d j=%d M=%g\n", i+o, it->jcol+o,it->M);
      assert(it->jcol<=i);
      ii[itnz]=i+o; 
      jj[itnz]=it->jcol+o;
      MM[itnz]=it->M;
      itnz++;
    }
  }
}

void SparseSymMatrixRowMajList::printMatrixInMatlab( char *name)
{
  assert(false && "not implemented");
}

bool SparseSymMatrixRowMajList::
fromGetSparseTriplet_w_patternMatch(const int* irow, const int* jcol, const int& nnz_in, double* M_out)
{
  if(nnz_in!=nnz) return false;
  if(nnz_in==0) return true;

  int it_in=0; 
  int row=irow[it_in];
  list<ColVal>::iterator it=vlmat[row].begin(); 
  for(int it_in=0; it_in<nnz_in; ) {
    assert(row<=irow[it_in]);    

    if(row<irow[it_in]) {
      row = irow[it_in];
      it=vlmat[row].begin(); 
    }
    while(it!=vlmat[row].end() && it->jcol<jcol[it_in]) {
      ++it;
    }
    if(it==vlmat[row].end()) {
      return false; //could not find irow[it_in] and jcol[it_in] in 'this'
    }
    //assert(it->jcol==jcol[it_in]);
    if(it->jcol > jcol[it_in]) return false;

    M_out[it_in] = it->M;

    ++it; it_in++;
  }
  
  return true;
}

// Performs the following 
//
// M[ this.(i,j) \setintersection (irn,jcn)] = this.M
// M[(irn,jcn) \setdiff this.(i,j) ]         = 0.
// (irn_diff,jcn_diff) = this.(i,j) \setdiff (irn,jcn) and M_diff = this[irn_diff,jcn_diff]
//
// Input:
//  - irn, jcn: triplet indexes, all of size nnz_
// 
// Output:
//  - M: values corresponding to entries in (irn,jcn), of size nnz_. The method copies in M the values 
// stored in 'this' that have indexes present in (irn,jcn). When (i,j) is in (irn,jcn) but not in 'this',
// M[i,j] is set to 0. The indexes in 'this' that are not in (irn,jcn) are returned in (irn_diff, jcn_diff)
//  - irn_diff, jcn_diff: allocated by this function or returned as NULL; caller is responsible for freeing
// these; contain the indexes in 'this' that are not in (irn,jcn)
//  - nnz_diff: number of nz in (irn_diff,jcn_diff)
//  - M_diff: entries of the difference
// 
// Return: false if diff is non-empty, otherwise true
bool SparseSymMatrixRowMajList::
fromGetIntersectionSparseTriplet_w_diff(const int* irow, const int* jcol, const int& nnz_in, 
					double* M_out,
					int** irow_diff, int** jcol_diff, double**M_diff, int& nnz_diff)

{
  *irow_diff = *jcol_diff = NULL; nnz_diff=0;
  if(nnz_in==0) {
    if(this->nnz == 0) return true;
  }
  
  list<ColVal>::iterator it_this;
  std::vector<int> virow_diff, vjcol_diff;
  std::vector<double> vM_diff;
  virow_diff.reserve(1024); vjcol_diff.reserve(1024); vM_diff.reserve(1024);//should cover most power grid problems

  int row_in_start = nnz_in>0?irow[0]:0;

  //elements in this' rows up to row_in_start will be in diff
  for(int row=0; row<row_in_start; row++) {
    
    for(it_this = this->vlmat[row].begin(); it_this!=this->vlmat[row].end(); ++it_this) {
      virow_diff.push_back(row);
      vjcol_diff.push_back(it_this->jcol);
      vM_diff.push_back(it_this->M);
    }
  }

  int row=-1; 
  // go over rows that are present both in this and (irow,jcol)
  for(int it_in=0; it_in<nnz_in; ) {
    assert(row<=irow[it_in]);
    if(row<irow[it_in]) {
      row=irow[it_in];
      it_this=vlmat[row].begin(); 
    }


    //iterate in this' row till jcol[it_in] is reached or not found
    while(it_this!=vlmat[row].end() && it_this->jcol<jcol[it_in]) {
      // current element in this not in (irow,jcol)
      virow_diff.push_back(row);
      vjcol_diff.push_back(it_this->jcol);
      vM_diff.push_back(it_this->M);
      ++it_this;
    }

    if(it_this==vlmat[row].end()) {
      //could not find irow[it_in] and jcol[it_in] in 'this'
      M_out[it_in] = 0.;
      it_in++;
    } else {
      // we have that it_this->jcol >= jcol[it_in])
      if(it_this->jcol==jcol[it_in]) {
	//element found
	M_out[it_in] = it_this->M;
	++it_this; it_in++;
      } else {
	assert(it_this->jcol>jcol[it_in]);
	//element in 'this' not in (irow,jcol)
	virow_diff.push_back(row);
	vjcol_diff.push_back(it_this->jcol);
	vM_diff.push_back(it_this->M);
	++it_this; it_in++;
      }
    }
  } //end for over elements in (irow,jcol)

  // go over rows that are present in this but not in (irow,jcol) 
  for(; row<this->vlmat.size(); row++) {
    for(it_this = this->vlmat[row].begin(); it_this!=this->vlmat[row].end(); ++it_this) {
      //element in diff
      virow_diff.push_back(row);
      vjcol_diff.push_back(it_this->jcol);
      vM_diff.push_back(it_this->M);
    }
  }

  nnz_diff = virow_diff.size();
  assert(nnz_diff == vjcol_diff.size());
  assert(nnz_diff == vM_diff.size());
  if(nnz_diff==0) 
    return true;
  //else

  *irow_diff = new int[nnz_diff];
  assert(NULL != *irow_diff && "insufficient memory !?!");
  std::copy(virow_diff.begin(), virow_diff.end(), *irow_diff);
  std::vector<int>().swap(virow_diff);//deallocate before the next allocation

  *jcol_diff = new int[nnz_diff];
  assert(NULL != *jcol_diff && "insufficient memory !?!");
  std::copy(vjcol_diff.begin(), vjcol_diff.end(), *jcol_diff);
  std::vector<int>().swap(vjcol_diff);//deallocate

  *M_diff = new double[nnz_diff];
  assert(NULL != *M_diff && "insufficient memory !?!");
  std::copy(vM_diff.begin(), vM_diff.end(), *M_diff);

  return false;
}

//copies (i,j,M) to 'this'; returns false if an entry of (i,j,M) not found in 'this', otherwise true.
bool SparseSymMatrixRowMajList::
atPutSparseTriplet(const int* irow, const int* jcol, const double* M, const int& nnz_in)
{
  if(nnz_in!=nnz) return false;
  if(nnz_in==0) return true;

  int it_in=0;
  int row=irow[it_in];
  list<ColVal>::iterator it=vlmat[row].begin();
  while(it_in<nnz_in) {
    assert(row<=irow[it_in]);

    if(row<irow[it_in]) {
      row = irow[it_in];
      it=vlmat[row].begin(); 
    }
    while(it!=vlmat[row].end() && it->jcol<jcol[it_in]) {
      ++it;
    }
    if(it==vlmat[row].end()) {
      return false; //could not find irow[it_in] and jcol[it_in] in 'this'
    }
    assert(it->jcol==jcol[it_in]);
    it->M = M[it_in];

    ++it; it_in++;
  }
  
  return true;
}
