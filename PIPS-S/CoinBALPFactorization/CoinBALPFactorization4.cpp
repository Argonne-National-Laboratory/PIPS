/* $Id: CoinFactorization.hpp 1767 2015-01-05 12:36:13Z forrest $ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).
/* 
Authors

John Forrest
*/

/* This PIPS-S file is under the EPL license since it is derivative work of 
   COIN-OR
*/
#include "CoinBALPFactorization.hpp"



// For semi-sparse
#define BITS_PER_CHECK 8
#define CHECK_SHIFT 3
typedef unsigned char CoinCheckZero;

// routines for BTRAN

void CoinBALPFactorization::updateColumnTransposeUQ ( CoinIndexedVector * region,
				      CoinIndexedVector * vec) const {
	assert(!vec->packedMode());
	int * regionIndex = region->getIndices();
	int numberNonZero;
	const int *back = pivotColumnBack_.array();
	double * regionElts = region->denseVector();

	numberNonZero = vec->getNumElements();
	int * vecIdx = vec->getIndices();
	double * vecElts = vec->denseVector();
	int nnz = 0;
	for (int i = 0; i < numberNonZero; i++) {
		int row = vecIdx[i];
		double value = vecElts[row];
		//if (row >= numberColumns_) continue; // applying ETAs leaves some of these
		vecElts[row] = 0.0;
		if (value == COIN_INDEXED_REALLY_TINY_ELEMENT) continue;
		assert(row < numberColumns_);
		row = back[row];
		regionIndex[nnz++] = row;
		regionElts[row] = value;
	}
	region->setNumElements(nnz);
	vec->setNumElements(0);

	CoinFactorizationDouble * COIN_RESTRICT pivotRegion = pivotRegion_.array();
	int smallestIndex=numberRowsExtra_;
	for (int j = 0; j < nnz; j++ ) {
		int iRow = regionIndex[j];
		smallestIndex = CoinMin(smallestIndex,iRow);
		regionElts[iRow] *= pivotRegion[iRow];
	}
	updateColumnTransposeU ( region,smallestIndex );

}

void CoinBALPFactorization::updateColumnTransposeG ( CoinIndexedVector * region,
				      CoinIndexedVector * vec) const {
	assert(!vec->packedMode());
	int * regionIndex = region->getIndices();
	int numberNonZero;
	const int *back = permuteBack_.array();
	double * regionElts = region->denseVector();

	updateColumnTransposeL( region );

	numberNonZero = region->getNumElements();
	int * vecIdx = vec->getIndices();
	double * vecElts = vec->denseVector();
	for (int j = 0; j < numberNonZero; j++) {
		int row = regionIndex[j];
		double value = regionElts[row];
		regionElts[row] = 0.0;
		assert(row < numberRows_);
		row = back[row];
		vecElts[row] = value;
		vecIdx[j] = row;
	}
	region->setNumElements(0);
	vec->setNumElements(numberNonZero);

}

int
CoinBALPFactorization::updateColumnTranspose ( CoinIndexedVector * region,
                                          CoinIndexedVector * vec ) 
  const
{
	 updateColumnTransposeUQ(region,vec);
	 updateColumnTransposeG(region,vec);
	 return vec->getNumElements();
}



//  updateColumnTransposeU.  Updates part of column transpose (BTRANU)
//assumes index is sorted i.e. region is correct
//does not sort by sign
void
CoinBALPFactorization::updateColumnTransposeU ( CoinIndexedVector * regionSparse,
					    int smallestIndex) const
{
  int number = regionSparse->getNumElements (  );
  int goSparse;
  if (collectStatistics_) btranCountInput_ += number;
  // Guess at number at end
  if (sparseThreshold_>0) {
    if (btranAverageAfterU_) {
      int newNumber = static_cast<int> (number*btranAverageAfterU_);
      if (newNumber< sparseThreshold_)
	goSparse = 2;
      else if (newNumber< sparseThreshold2_)
	goSparse = 1;
      else
	goSparse = 0;
    } else {
      if (number<sparseThreshold_) 
	goSparse = 2;
      else
	goSparse = 0;
    }
  } else {
    goSparse=0;
  }
  switch (goSparse) {
  case 0: // densish
    updateColumnTransposeUDensish(regionSparse,smallestIndex);
    break;
  case 1: // middling
    updateColumnTransposeUSparsish(regionSparse,smallestIndex);
    break;
  case 2: // sparse
    updateColumnTransposeUSparse(regionSparse);
    break;
  }
  if (collectStatistics_) { 
    btranCountAfterU_ += regionSparse->getNumElements (  );
    numberBtranCounts_++;
  }
}


/* Updates part of column transpose (BTRANU) when densish,
   assumes index is sorted i.e. region is correct */
void 
CoinBALPFactorization::updateColumnTransposeUDensish 
                        ( CoinIndexedVector * regionSparse,
			  int smallestIndex) const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int numberNonZero = regionSparse->getNumElements (  );
  double tolerance = zeroTolerance_;
  
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  
  const CoinBigIndex *startRow = startRowU_.array();
  
  const CoinBigIndex *convertRowToColumn = convertRowToColumnU_.array();
  const int *indexColumn = indexColumnU_.array();
  
  const CoinFactorizationDouble * element = elementU_.array();
  int last = numberU_;
  
  const int *numberInRow = numberInRow_.array();
  numberNonZero = 0;
  for (int i=smallestIndex ; i < last; i++ ) {
    CoinFactorizationDouble pivotValue = region[i];
    if ( fabs ( pivotValue ) > tolerance ) {
      CoinBigIndex start = startRow[i];
      int numberIn = numberInRow[i];
      CoinBigIndex end = start + numberIn;
      for (CoinBigIndex j = start ; j < end; j ++ ) {
	int iRow = indexColumn[j];
	CoinBigIndex getElement = convertRowToColumn[j];
	CoinFactorizationDouble value = element[getElement];
	region[iRow] -=  value * pivotValue;
      }     
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }       
  }
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}
/* Updates part of column transpose (BTRANU) when sparsish,
      assumes index is sorted i.e. region is correct */
void 
CoinBALPFactorization::updateColumnTransposeUSparsish 
                        ( CoinIndexedVector * regionSparse,
			  int smallestIndex) const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int numberNonZero = regionSparse->getNumElements (  );
  double tolerance = zeroTolerance_;
  
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  
  const CoinBigIndex *startRow = startRowU_.array();
  
  const CoinBigIndex *convertRowToColumn = convertRowToColumnU_.array();
  const int *indexColumn = indexColumnU_.array();
  
  const CoinFactorizationDouble * element = elementU_.array();
  int last = numberU_;
  
  const int *numberInRow = numberInRow_.array();
  
  // mark known to be zero
  int nInBig = sizeof(CoinBigIndex)/sizeof(int);
  CoinCheckZero * COIN_RESTRICT mark = reinterpret_cast<CoinCheckZero *> (sparse_.array()+(2+nInBig)*maximumRowsExtra_);

  for (int i=0;i<numberNonZero;i++) {
    int iPivot=regionIndex[i];
    int iWord = iPivot>>CHECK_SHIFT;
    int iBit = iPivot-(iWord<<CHECK_SHIFT);
    if (mark[iWord]) {
      mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
    } else {
      mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
    }
  }

  numberNonZero = 0;
  // Find convenient power of 2
  smallestIndex = smallestIndex >> CHECK_SHIFT;
  int kLast = last>>CHECK_SHIFT;
  // do in chunks

  for (int k=smallestIndex;k<kLast;k++) {
    unsigned int iMark = mark[k];
    if (iMark) {
      // something in chunk - do all (as imark may change)
      int i = k<<CHECK_SHIFT;
      int iLast = i+BITS_PER_CHECK;
      for ( ; i < iLast; i++ ) {
	CoinFactorizationDouble pivotValue = region[i];
	if ( fabs ( pivotValue ) > tolerance ) {
	  CoinBigIndex start = startRow[i];
	  int numberIn = numberInRow[i];
	  CoinBigIndex end = start + numberIn;
	  for (CoinBigIndex j = start ; j < end; j ++ ) {
	    int iRow = indexColumn[j];
	    CoinBigIndex getElement = convertRowToColumn[j];
	    CoinFactorizationDouble value = element[getElement];
	    int iWord = iRow>>CHECK_SHIFT;
	    int iBit = iRow-(iWord<<CHECK_SHIFT);
	    if (mark[iWord]) {
	      mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
	    } else {
	      mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
	    }
	    region[iRow] -=  value * pivotValue;
	  }     
	  regionIndex[numberNonZero++] = i;
	} else {
	  region[i] = 0.0;
	}       
      }
      mark[k]=0;
    }
  }
  mark[kLast]=0;
  for (int i = kLast<<CHECK_SHIFT; i < last; i++ ) {
    CoinFactorizationDouble pivotValue = region[i];
    if ( fabs ( pivotValue ) > tolerance ) {
      CoinBigIndex start = startRow[i];
      int numberIn = numberInRow[i];
      CoinBigIndex end = start + numberIn;
      for (CoinBigIndex j = start ; j < end; j ++ ) {
	int iRow = indexColumn[j];
	CoinBigIndex getElement = convertRowToColumn[j];
	CoinFactorizationDouble value = element[getElement];

	region[iRow] -=  value * pivotValue;
      }     
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }       
  }
#ifdef COIN_DEBUG
  for (int i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}
/* Updates part of column transpose (BTRANU) when sparse,
   assumes index is sorted i.e. region is correct */
void 
CoinBALPFactorization::updateColumnTransposeUSparse ( 
		   CoinIndexedVector * regionSparse) const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int numberNonZero = regionSparse->getNumElements (  );
  double tolerance = zeroTolerance_;
  
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  const CoinBigIndex *startRow = startRowU_.array();
  
  const CoinBigIndex *convertRowToColumn = convertRowToColumnU_.array();
  const int *indexColumn = indexColumnU_.array();
  
  const CoinFactorizationDouble * element = elementU_.array();
  
  const int *numberInRow = numberInRow_.array();
  
  // use sparse_ as temporary area
  // mark known to be zero
  int * COIN_RESTRICT stack = sparse_.array();  /* pivot */
  int * COIN_RESTRICT list = stack + maximumRowsExtra_;  /* final list */
  CoinBigIndex * COIN_RESTRICT next = reinterpret_cast<CoinBigIndex *> (list + maximumRowsExtra_);  /* jnext */
  char * COIN_RESTRICT mark = reinterpret_cast<char *> (next + maximumRowsExtra_);
  int nList;
#ifdef COIN_DEBUG
  for (int i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif
  nList=0;
  for (int k=0;k<numberNonZero;k++) {
    int kPivot=regionIndex[k];
    stack[0]=kPivot;
    CoinBigIndex j=startRow[kPivot]+numberInRow[kPivot]-1;
    next[0]=j;
    int nStack=1;
    while (nStack) {
      /* take off stack */
      kPivot=stack[--nStack];
      if (mark[kPivot]!=1) {
	j=next[nStack];
	if (j>=startRow[kPivot]) {
	  kPivot=indexColumn[j--];
	  /* put back on stack */
	  next[nStack++] =j;
	  if (!mark[kPivot]) {
	    /* and new one */
	    j=startRow[kPivot]+numberInRow[kPivot]-1;
	    stack[nStack]=kPivot;
	    mark[kPivot]=2;
	    next[nStack++]=j;
	  }
	} else {
	  // finished
	  list[nList++]=kPivot;
	  mark[kPivot]=1;
	}
      }
    }
  }
  numberNonZero=0;
  for (int i=nList-1;i>=0;i--) {
    int iPivot = list[i];
    mark[iPivot]=0;
    CoinFactorizationDouble pivotValue = region[iPivot];
    if ( fabs ( pivotValue ) > tolerance ) {
      CoinBigIndex start = startRow[iPivot];
      int numberIn = numberInRow[iPivot];
      CoinBigIndex end = start + numberIn;
      for (CoinBigIndex j=start ; j < end; j ++ ) {
	int iRow = indexColumn[j];
	CoinBigIndex getElement = convertRowToColumn[j];
	CoinFactorizationDouble value = element[getElement];
	region[iRow] -= value * pivotValue;
      }     
      regionIndex[numberNonZero++] = iPivot;
    } else {
      region[iPivot] = 0.0;
    }       
  }       
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}



//  updateColumnTransposeL.  Updates part of column transpose (BTRANL)
void
CoinBALPFactorization::updateColumnTransposeL ( CoinIndexedVector * regionSparse ) const
{
  int number = regionSparse->getNumElements (  );
  if (!numberL_&&!numberDense_) {
    if (sparse_.array()||number<numberRows_)
      return;
  }
  int goSparse;
  /*if (sparseThreshold_>0) {
    if (btranAverageAfterL_) {
      int newNumber = static_cast<int> (number*btranAverageAfterL_);
      if (newNumber< sparseThreshold_)
	goSparse = 2;
      else if (newNumber< sparseThreshold2_)
	goSparse = 1;
      else
	goSparse = 0;
    } else {
      if (number<sparseThreshold_) 
	goSparse = 2;
      else
	goSparse = 0;
    }
  } else {
    goSparse=-1;
  }*/
  if (sparseThreshold_>0) {
  	goSparse = 0;
  } else {
  	goSparse = -1;
  }
  switch (goSparse) {
  case -1: // No row copy
    updateColumnTransposeLDensish(regionSparse);
    break;
  case 0: // densish but by row
    updateColumnTransposeLByRow(regionSparse);
    break;
  case 1: // middling(and by row)
    updateColumnTransposeLSparsish(regionSparse);
    break;
  case 2: // sparse
    updateColumnTransposeLSparse(regionSparse);
    break;
  }
}


/*  updateColumnTransposeLDensish.  
    Updates part of column transpose (BTRANL) dense by column */
// uses a plain triangular solve approach
void
CoinBALPFactorization::updateColumnTransposeLDensish 
     ( CoinIndexedVector * regionSparse ) const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  int numberNonZero;
  double tolerance = zeroTolerance_;
  int base;
  int first = -1;
  
  numberNonZero=0;
  //scan
  for (first=numberRows_-1;first>=0;first--) {
    if (region[first]) 
      break;
  }
  if ( first >= 0 ) {
    base = baseL_;
    const CoinBigIndex * COIN_RESTRICT startColumn = startColumnL_.array();
    const int * COIN_RESTRICT indexRow = indexRowL_.array();
    const CoinFactorizationDouble * COIN_RESTRICT element = elementL_.array();
    //int last = baseL_ + numberL_;
    int last = numberRows_;

    if ( first >= last ) {
      first = last - 1;
    }       
    for (int i = first ; i >= base; i-- ) {
      CoinBigIndex j;
      CoinFactorizationDouble pivotValue = region[i];
      for ( j= startColumn[i] ; j < startColumn[i+1]; j++ ) {
	int iRow = indexRow[j];
	CoinFactorizationDouble value = element[j];
	pivotValue -= value * region[iRow];
      }       
      if ( fabs ( pivotValue ) > tolerance ) {
	region[i] = pivotValue;
	regionIndex[numberNonZero++] = i;
      } else { 
	region[i] = 0.0;
      }       
    }       
    //may have stopped early
    if ( first < base ) {
      base = first + 1;
    }
    if (base > 5) { // why 5?
      int i=base-1;
      CoinFactorizationDouble pivotValue=region[i];
      bool store = fabs(pivotValue) > tolerance;
      for (; i > 0; i-- ) {
	bool oldStore = store;
	CoinFactorizationDouble oldValue = pivotValue;
	pivotValue = region[i-1];
	store = fabs ( pivotValue ) > tolerance ;
	if (!oldStore) {
	  region[i] = 0.0;
	} else {
	  region[i] = oldValue;
	  regionIndex[numberNonZero++] = i;
	}       
      }     
      if (store) {
	region[0] = pivotValue;
	regionIndex[numberNonZero++] = 0;
      } else {
	region[0] = 0.0;
      }       
    } else {
      for (int i = base -1 ; i >= 0; i-- ) {
	CoinFactorizationDouble pivotValue = region[i];
	if ( fabs ( pivotValue ) > tolerance ) {
	  region[i] = pivotValue;
	  regionIndex[numberNonZero++] = i;
	} else {
	  region[i] = 0.0;
	}       
      }     
    }
  } 
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}
/*  updateColumnTransposeLByRow. 
    Updates part of column transpose (BTRANL) densish but by row */
void
CoinBALPFactorization::updateColumnTransposeLByRow 
    ( CoinIndexedVector * regionSparse ) const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  int numberNonZero;
  double tolerance = zeroTolerance_;
  int first = -1;
  
  // use row copy of L
  const CoinFactorizationDouble * element = elementByRowL_.array();
  const CoinBigIndex * startRow = startRowL_.array();
  const int * column = indexColumnL_.array();
  for (first=numberRows_-1;first>=0;first--) {
    if (region[first]) 
      break;
  }
  numberNonZero=0;
  for (int i=first;i>=0;i--) {
    CoinFactorizationDouble pivotValue = region[i];
    if ( fabs ( pivotValue ) > tolerance ) {
      regionIndex[numberNonZero++] = i;
      CoinBigIndex j;
      for (j = startRow[i + 1]-1;j >= startRow[i]; j--) {
	int iRow = column[j];
	CoinFactorizationDouble value = element[j];
	region[iRow] -= pivotValue*value;
      }
    } else {
      region[i] = 0.0;
    }     
  }
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}
// Updates part of column transpose (BTRANL) when sparsish by row
void
CoinBALPFactorization::updateColumnTransposeLSparsish 
    ( CoinIndexedVector * regionSparse ) const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  int numberNonZero = regionSparse->getNumElements();
  double tolerance = zeroTolerance_;
  
  // use row copy of L
  const CoinFactorizationDouble * element = elementByRowL_.array();
  const CoinBigIndex * startRow = startRowL_.array();
  const int * column = indexColumnL_.array();
  // mark known to be zero
  int nInBig = sizeof(CoinBigIndex)/sizeof(int);
  CoinCheckZero * COIN_RESTRICT mark = reinterpret_cast<CoinCheckZero *> (sparse_.array()+(2+nInBig)*maximumRowsExtra_);
  for (int i=0;i<numberNonZero;i++) {
    int iPivot=regionIndex[i];
    int iWord = iPivot>>CHECK_SHIFT;
    int iBit = iPivot-(iWord<<CHECK_SHIFT);
    if (mark[iWord]) {
      mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
    } else {
      mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
    }
  }
  numberNonZero = 0;
  // First do down to convenient power of 2
  int jLast = (numberRows_-1)>>CHECK_SHIFT;
  jLast = (jLast<<CHECK_SHIFT);
  for (int i=numberRows_-1;i>=jLast;i--) {
    CoinFactorizationDouble pivotValue = region[i];
    if ( fabs ( pivotValue ) > tolerance ) {
      regionIndex[numberNonZero++] = i;
      CoinBigIndex j;
      for (j = startRow[i + 1]-1;j >= startRow[i]; j--) {
	int iRow = column[j];
	CoinFactorizationDouble value = element[j];
	int iWord = iRow>>CHECK_SHIFT;
	int iBit = iRow-(iWord<<CHECK_SHIFT);
	if (mark[iWord]) {
	  mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
	} else {
	  mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
	}
	region[iRow] -= pivotValue*value;
      }
    } else {
      region[i] = 0.0;
    }     
  }
  // and in chunks
  jLast = jLast>>CHECK_SHIFT;
  mark[jLast]=0;
  for (int k=jLast-1;k>=0;k--) {
    unsigned int iMark = mark[k];
    if (iMark) {
      // something in chunk - do all (as imark may change)
      int iLast = k<<CHECK_SHIFT;
      for (int i = iLast+BITS_PER_CHECK-1; i >= iLast; i-- ) {
	CoinFactorizationDouble pivotValue = region[i];
	if ( fabs ( pivotValue ) > tolerance ) {
	  regionIndex[numberNonZero++] = i;
	  CoinBigIndex j;
	  for (j = startRow[i + 1]-1;j >= startRow[i]; j--) {
	    int iRow = column[j];
	    CoinFactorizationDouble value = element[j];
	    int iWord = iRow>>CHECK_SHIFT;
	    int iBit = iRow-(iWord<<CHECK_SHIFT);
	    if (mark[iWord]) {
	      mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
	    } else {
	      mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
	    }
	    region[iRow] -= pivotValue*value;
	  }
	} else {
	  region[i] = 0.0;
	}     
      }
      mark[k]=0;
    }
  }
#ifdef COIN_DEBUG
  for (int i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}
/*  updateColumnTransposeLSparse. 
    Updates part of column transpose (BTRANL) sparse */
void
CoinBALPFactorization::updateColumnTransposeLSparse 
    ( CoinIndexedVector * regionSparse ) const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  int numberNonZero = regionSparse->getNumElements (  );
  double tolerance = zeroTolerance_;
  
  // use row copy of L
  const CoinFactorizationDouble * element = elementByRowL_.array();
  const CoinBigIndex * startRow = startRowL_.array();
  const int * column = indexColumnL_.array();
  // use sparse_ as temporary area
  // mark known to be zero
  int * COIN_RESTRICT stack = sparse_.array();  /* pivot */
  int * COIN_RESTRICT list = stack + maximumRowsExtra_;  /* final list */
  CoinBigIndex * COIN_RESTRICT next = reinterpret_cast<CoinBigIndex *> (list + maximumRowsExtra_);  /* jnext */
  char * COIN_RESTRICT mark = reinterpret_cast<char *> (next + maximumRowsExtra_);
  int nList;
  int number = numberNonZero;
#ifdef COIN_DEBUG
  for (i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif
  nList=0;
  for (int k=0;k<number;k++) {
    int kPivot=regionIndex[k];
    if(!mark[kPivot]&&region[kPivot]) {
      stack[0]=kPivot;
      CoinBigIndex j=startRow[kPivot+1]-1;
      int nStack=0;
      while (nStack>=0) {
	/* take off stack */
	if (j>=startRow[kPivot]) {
	  int jPivot=column[j--];
	  /* put back on stack */
	  next[nStack] =j;
	  if (!mark[jPivot]) {
	    /* and new one */
	    kPivot=jPivot;
	    j = startRow[kPivot+1]-1;
	    stack[++nStack]=kPivot;
	    mark[kPivot]=1;
	    next[nStack]=j;
	  }
	} else {
	  /* finished so mark */
	  list[nList++]=kPivot;
	  mark[kPivot]=1;
	  --nStack;
	  if (nStack>=0) {
	    kPivot=stack[nStack];
	    j=next[nStack];
	  }
	}
      }
    }
  }
  numberNonZero=0;
  for (int i=nList-1;i>=0;i--) {
    int iPivot = list[i];
    mark[iPivot]=0;
    CoinFactorizationDouble pivotValue = region[iPivot];
    if ( fabs ( pivotValue ) > tolerance ) {
      regionIndex[numberNonZero++] = iPivot;
      CoinBigIndex j;
      for ( j = startRow[iPivot]; j < startRow[iPivot+1]; j ++ ) {
	int iRow = column[j];
	CoinFactorizationDouble value = element[j];
	region[iRow] -= value * pivotValue;
      }
    } else {
      region[iPivot]=0.0;
    }
  }
  //set counts
  regionSparse->setNumElements ( numberNonZero );
}
