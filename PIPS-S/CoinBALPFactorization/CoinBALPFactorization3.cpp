#include "CoinBALPFactorization.hpp"

// routines for T and FTRAN



void CoinBALPFactorization::getRowOfT(CoinIndexedVector &v, int row) {
	assert(v.getNumElements() == 0);
	assert(!v.packedMode());
	assert(v.capacity() >= numberColumnsT_);
	if (row < numberColumns_) printf("asked for XT row, why?\n");
	int nElts = numberInRowT_.array()[row];
	double *elts = v.denseVector();
	int *idx = v.getIndices();
	
	CoinBigIndex *convertRowToColumnT = convertRowToColumnT_.array();
	int *indexColumn = indexColumnT_.array();
	CoinFactorizationDouble *elementT = elementT_.array();
	CoinBigIndex start = startRowT_.array()[row];
	CoinBigIndex end = start + nElts;
	assert(end == startRowT_.array()[row+1]);

	nElts = 0;
	for (CoinBigIndex j = start; j < end; j++) {
		int col = indexColumn[j];
		idx[nElts++] = col;
		elts[col] = elementT[convertRowToColumnT[j]];
	}
	assert(nElts == numberInRowT_.array()[row]);
	v.setNumElements(nElts);
	
}

// copy a row of T into packed form
// assume buffers have enough room
// user must call getNumberInRowOfT first
void CoinBALPFactorization::getRowOfT(int *colIdx,double *elts, int row) {
	assert(row >= 0 && row <= numberRows_);
	int nElts = numberInRowT_.array()[row];
	
	CoinBigIndex *convertRowToColumnT = convertRowToColumnT_.array();
	int *indexColumn = indexColumnT_.array();
	CoinFactorizationDouble *elementT = elementT_.array();
	CoinBigIndex start = startRowT_.array()[row];
	CoinBigIndex end = start + nElts;
	assert(end == startRowT_.array()[row+1]);

	nElts = 0;
	for (CoinBigIndex j = start; j < end; j++) {
		int col = indexColumn[j];
		colIdx[nElts] = col;
		elts[nElts++] = elementT[convertRowToColumnT[j]];
	}
}

int CoinBALPFactorization::getNumberInRowOfT(int row) const {
	if (row <= numberRows_ - numberColumns_) printf("asked for XT row, why?\n");
	assert(row >= 0 && row <= numberRows_);
	return numberInRowT_.array()[row];
}


void CoinBALPFactorization::multXT(const CoinIndexedVector &in, CoinIndexedVector &out) {

	assert(!in.packedMode() && !out.packedMode());
	assert(in.capacity() >= numberColumnsT_);
	assert(out.capacity() >= numberColumns_);

	const double *inElts = in.denseVector();
	double *outElts = out.denseVector();
	const int *inIdx = in.getIndices();
	int *outIdx = out.getIndices();
	const CoinFactorizationDouble *element = elementT_.array();

	int nEltsIn = in.getNumElements();
	
	if (nEltsIn < 0.1*numberColumnsT_) {
		//printf("did hypersparse multXT\n");
		// take linear combinations of the columns, ignoring bottom rows
		// it is easier to check here than in PRICE if to include the row

		CoinBigIndex *startColumn = startColumnT_.array();
		int *indexRow = indexRowT_.array();
		int *numberInColumn = numberInColumnT_.array();

		for (int j = 0; j < nEltsIn; j++) {
			int i = inIdx[j];
			double mult = -inElts[i]; // note negative sign
			CoinBigIndex start = startColumn[i];
			CoinBigIndex end = start + numberInColumn[i];
			for (CoinBigIndex q = start; q < end; q++) {
				int row = indexRow[q];
				if (row >= numberColumns_) continue;
				out.quickAdd(row,mult*element[q]);
			}
		}
	} else {
		//printf("did sparse multXT\n");
		// take dot products with the rows

		int nEltsOut = out.getNumElements();
		CoinBigIndex *startRow = startRowT_.array();
		int *indexColumn = indexColumnT_.array();
		int *numberInRow = numberInRowT_.array();
		CoinBigIndex *convertRowToColumn = convertRowToColumnT_.array();

		for (int row = 0; row < numberColumns_; row++) {
			CoinBigIndex start = startRow[row];
			CoinBigIndex end = start + numberInRow[row];
			double val = 0;
			for (CoinBigIndex q = start; q < end; q++) {
				int col = indexColumn[q];
				val += inElts[col]*element[convertRowToColumn[q]];
			}
			if (fabs(val) > zeroTolerance_) {
				if (!outElts[row]) outIdx[nEltsOut++] = row;
				outElts[row] -= val; // note negative sign
			}
		}
		out.setNumElements(nEltsOut);
	}

	//printf("multXT out density: %f\n", out.getNumElements()/static_cast<double>(numberColumns_));

}


void CoinBALPFactorization::multXTTranspose(const CoinIndexedVector &in, CoinIndexedVector &out) {

	assert(!in.packedMode() && !out.packedMode());
	assert(out.capacity() >= numberColumnsT_);
	assert(in.capacity() >= numberColumns_);

	const double *inElts = in.denseVector();
	const int *inIdx = in.getIndices();
	const CoinFactorizationDouble *element = elementT_.array();

	int nEltsIn = in.getNumElements();

	CoinBigIndex *startRow = startRowT_.array();
	int *indexColumn = indexColumnT_.array();
	int *numberInRow = numberInRowT_.array();
	CoinBigIndex *convertRowToColumn = convertRowToColumnT_.array();

	// take linear combination of the rows
	// given the position of this operation, input should be hyper-sparse
	for (int j = 0; j < nEltsIn; j++) {
		int row = inIdx[j];
		assert(row < numberColumns_);
		double mult = -inElts[row]; // note negative sign
		CoinBigIndex start = startRow[row];
		CoinBigIndex end = start+numberInRow[row];
		for (CoinBigIndex q = start; q < end; q++) {
			int col = indexColumn[q];
			out.quickAdd(col,mult*element[convertRowToColumn[q]]);
		}
	}
}

// add *negative* answer to out
void CoinBALPFactorization::multXTTranspose(const CoinIndexedVector &in, double *out) {

	assert(!in.packedMode());
	assert(in.capacity() >= numberColumns_);

	const double *inElts = in.denseVector();
	const int *inIdx = in.getIndices();
	const CoinFactorizationDouble *element = elementT_.array();

	int nEltsIn = in.getNumElements();

	CoinBigIndex *startRow = startRowT_.array();
	int *indexColumn = indexColumnT_.array();
	int *numberInRow = numberInRowT_.array();
	CoinBigIndex *convertRowToColumn = convertRowToColumnT_.array();

	// take linear combination of the rows
	// given the position of this operation, input should be hyper-sparse
	for (int j = 0; j < nEltsIn; j++) {
		int row = inIdx[j];
		assert(row < numberColumns_);
		double mult = inElts[row];
		CoinBigIndex start = startRow[row];
		CoinBigIndex end = start+numberInRow[row];
		for (CoinBigIndex q = start; q < end; q++) {
			int col = indexColumn[q];
			out[col] -= mult*element[convertRowToColumn[q]];
		}
	}
}

// For semi-sparse
#define BITS_PER_CHECK 8
#define CHECK_SHIFT 3
typedef unsigned char CoinCheckZero;


/* Updates one column (FTRAN) from region2 and permutes.
   region1 starts as zero
   Note - if regionSparse2 packed on input - will be packed on output
   - returns un-permuted result in region2 and region1 is zero */
int CoinBALPFactorization::updateColumn ( CoinIndexedVector * region,
				      CoinIndexedVector * vec) 
  const
{

	assert(numberRows_ == numberColumns_); // updateColumn doesn't make sense for BALP
	updateColumnG(region, vec);
	updateColumnUQ(region,vec);
	return vec->getNumElements();

}

void CoinBALPFactorization::updateColumnUQ ( CoinIndexedVector * region,
				      CoinIndexedVector * vec) const {
	assert(!vec->packedMode());
	int * regionIndex = region->getIndices();
	int numberNonZero;
	double * regionElts = region->denseVector();

	updateColumnU(region, regionIndex);

	const int *pivotColumn = pivotColumn_.array();
	numberNonZero = region->getNumElements();
	int * vecIdx = vec->getIndices();
	double * vecElts = vec->denseVector();
	int number = 0;
	for (int j = 0; j < numberNonZero; j++) {
		int row = regionIndex[j];
		double value = regionElts[row];
		regionElts[row] = 0.0;
		if (fabs(value) > zeroTolerance_) {
			row = pivotColumn[row];
			vecElts[row] = value;
			vecIdx[number++] = row;
		}
	}
	vec->setNumElements(number);
	region->setNumElements(0);

}

void CoinBALPFactorization::updateColumnG ( CoinIndexedVector * region,
				      CoinIndexedVector * vec) const {
	assert(!vec->packedMode());
	int * regionIndex = region->getIndices();
	int numberNonZero;
	const int *permute = permute_.array();
	double * regionElts = region->denseVector();

	numberNonZero = vec->getNumElements();
	int * vecIdx = vec->getIndices();
	double * vecElts = vec->denseVector();
	int nnz = 0;
	for (int j = 0; j < numberNonZero; j++) {
		int row = vecIdx[j];
		double value = vecElts[row];
		vecElts[row] = 0.0;
		if (value == COIN_INDEXED_REALLY_TINY_ELEMENT) continue; // applying ETAs leaves some of these
		assert(row < numberRows_);
		row = permute[row];
		regionElts[row] = value;
		regionIndex[nnz++] = row;
	}
	region->setNumElements(nnz);
	updateColumnL( region, regionIndex );

}

//  updateColumnL.  Updates part of column (FTRANL)
void
CoinBALPFactorization::updateColumnL ( CoinIndexedVector * regionSparse,
					   int * COIN_RESTRICT regionIndex) const
{
  if (collectStatistics_) ftranCountInput_ += regionSparse->getNumElements();
  if (numberL_) {
    int number = regionSparse->getNumElements (  );
    int goSparse;
    // Guess at number at end
    if (sparseThreshold_>0) {
      if (ftranAverageAfterL_) {
	int newNumber = static_cast<int> (number*ftranAverageAfterL_);
	if (newNumber< sparseThreshold_&&(numberL_<<2)>newNumber)
	  goSparse = 2;
	else if (newNumber< sparseThreshold2_&&(numberL_<<1)>newNumber)
	  goSparse = 1;
	else
	  goSparse = 0;
      } else {
	if (number<sparseThreshold_&&(numberL_<<2)>number) 
	  goSparse = 2;
	else
	  goSparse = 0;
      }
    } else {
      goSparse=0;
    }
    switch (goSparse) {
    case 0: // densish
      updateColumnLDensish(regionSparse,regionIndex);
      break;
    case 1: // middling
      updateColumnLSparsish(regionSparse,regionIndex);
      break;
    case 2: // sparse
      updateColumnLSparse(regionSparse,regionIndex);
      break;
    }
  }
  if (collectStatistics_)
  	ftranCountAfterL_ += regionSparse->getNumElements();
}
// Updates part of column (FTRANL) when densish
void 
CoinBALPFactorization::updateColumnLDensish ( CoinIndexedVector * regionSparse ,
					  int * COIN_RESTRICT regionIndex)
  const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int number = regionSparse->getNumElements (  );
  int numberNonZero;
  double tolerance = zeroTolerance_;
  
  numberNonZero = 0;
  
  const CoinBigIndex * COIN_RESTRICT startColumn = startColumnL_.array();
  const int * COIN_RESTRICT indexRow = indexRowL_.array();
  const CoinFactorizationDouble * COIN_RESTRICT element = elementL_.array();
  int last = numberColumns1_;
  assert ( last == baseL_ + numberL_);
  int smallestIndex = numberRowsExtra_;
  // do easy ones
  for (int k=0;k<number;k++) {
    int iPivot=regionIndex[k];
    if (iPivot>=baseL_) 
      smallestIndex = CoinMin(iPivot,smallestIndex);
    else
      regionIndex[numberNonZero++]=iPivot;
  }
  // now others
  for (int i = smallestIndex; i < last; i++ ) {
    CoinFactorizationDouble pivotValue = region[i];
    
    if ( fabs(pivotValue) > tolerance ) {
      CoinBigIndex start = startColumn[i];
      CoinBigIndex end = startColumn[i + 1];
      for (CoinBigIndex j = start; j < end; j ++ ) {
	int iRow = indexRow[j];
	CoinFactorizationDouble result = region[iRow];
	CoinFactorizationDouble value = element[j];

	region[iRow] = result - value * pivotValue;
      }     
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }       
  }
  // catch nonzeros in lower part
  for (int i = last; i < numberRows_; i++) {
	  CoinFactorizationDouble value = region[i];
	  if ( fabs(value) > tolerance) {
		  regionIndex[numberNonZero++] = i;
	  }
	  else {
		  region[i] = 0.0;
	  }
  }
   
  regionSparse->setNumElements ( numberNonZero );
} 
// Updates part of column (FTRANL) when sparsish
void 
CoinBALPFactorization::updateColumnLSparsish ( CoinIndexedVector * regionSparse,
					   int * COIN_RESTRICT regionIndex)
  const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int number = regionSparse->getNumElements (  );
  int numberNonZero;
  double tolerance = zeroTolerance_;
  
  numberNonZero = 0;
  
  const CoinBigIndex *startColumn = startColumnL_.array();
  const int *indexRow = indexRowL_.array();
  const CoinFactorizationDouble *element = elementL_.array();
  int last = numberColumns1_;
  assert ( last == baseL_ + numberL_);
#if DENSE_CODE==1
  //can take out last bit of sparse L as empty
  last -= numberDense_;
#endif
  // mark known to be zero
  int nInBig = sizeof(CoinBigIndex)/sizeof(int);
  CoinCheckZero * COIN_RESTRICT mark = reinterpret_cast<CoinCheckZero *> (sparse_.array()+(2+nInBig)*maximumRowsExtra_);
  int smallestIndex = numberRowsExtra_;
  // do easy ones
  for (int k=0;k<number;k++) {
    int iPivot=regionIndex[k];
    if (iPivot<baseL_) { 
      regionIndex[numberNonZero++]=iPivot;
    } else {
      smallestIndex = CoinMin(iPivot,smallestIndex);
      int iWord = iPivot>>CHECK_SHIFT;
      int iBit = iPivot-(iWord<<CHECK_SHIFT);
      if (mark[iWord]) {
	mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
      } else {
	mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
      }
    }
  }
  // now others
  // First do up to convenient power of 2
  int jLast = (smallestIndex+BITS_PER_CHECK-1)>>CHECK_SHIFT;
  jLast = CoinMin((jLast<<CHECK_SHIFT),last);
  int i;
  for ( i = smallestIndex; i < jLast; i++ ) {
    CoinFactorizationDouble pivotValue = region[i];
    CoinBigIndex start = startColumn[i];
    CoinBigIndex end = startColumn[i + 1];
    
    if ( fabs(pivotValue) > tolerance ) {
      for (CoinBigIndex j = start; j < end; j ++ ) {
	int iRow = indexRow[j];
	CoinFactorizationDouble result = region[iRow];
	CoinFactorizationDouble value = element[j];
	region[iRow] = result - value * pivotValue;
	int iWord = iRow>>CHECK_SHIFT;
	int iBit = iRow-(iWord<<CHECK_SHIFT);
	if (mark[iWord]) {
	  mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
	} else {
	  mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
	}
      }     
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }       
  }
  
  int kLast = last>>CHECK_SHIFT;
  if (jLast<last) {
    // now do in chunks
    for (int k=(jLast>>CHECK_SHIFT);k<kLast;k++) {
      unsigned int iMark = mark[k];
      if (iMark) {
	// something in chunk - do all (as imark may change)
	i = k<<CHECK_SHIFT;
	int iLast = i+BITS_PER_CHECK;
	for ( ; i < iLast; i++ ) {
	  CoinFactorizationDouble pivotValue = region[i];
	  CoinBigIndex start = startColumn[i];
	  CoinBigIndex end = startColumn[i + 1];
	  
	  if ( fabs(pivotValue) > tolerance ) {
	    CoinBigIndex j;
	    for ( j = start; j < end; j ++ ) {
	      int iRow = indexRow[j];
	      CoinFactorizationDouble result = region[iRow];
	      CoinFactorizationDouble value = element[j];
	      region[iRow] = result - value * pivotValue;
	      int iWord = iRow>>CHECK_SHIFT;
	      int iBit = iRow-(iWord<<CHECK_SHIFT);
	      if (mark[iWord]) {
		mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
	      } else {
		mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
	      }
	    }     
	    regionIndex[numberNonZero++] = i;
	  } else {
	    region[i] = 0.0;
	  }       
	}
	mark[k]=0; // zero out marked
      }
    }
    i = kLast<<CHECK_SHIFT;
  }
  for ( ; i < last; i++ ) {
    CoinFactorizationDouble pivotValue = region[i];
    CoinBigIndex start = startColumn[i];
    CoinBigIndex end = startColumn[i + 1];
    
    if ( fabs(pivotValue) > tolerance ) {
      for (CoinBigIndex j = start; j < end; j ++ ) {
	int iRow = indexRow[j];
	CoinFactorizationDouble result = region[iRow];
	CoinFactorizationDouble value = element[j];
	region[iRow] = result - value * pivotValue;
      }     
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }       
  }
  // Now trivial part
  for ( ; i < numberRows_; i++ ) {
    double pivotValue = region[i];
    if ( fabs(pivotValue) > tolerance ) {
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }       
  }
  // zero out ones that might have been skipped
  mark[smallestIndex>>CHECK_SHIFT]=0;
  int kkLast = (numberRows_+BITS_PER_CHECK-1)>>CHECK_SHIFT;
  CoinZeroN(mark+kLast,kkLast-kLast);
  regionSparse->setNumElements ( numberNonZero );
}
// Updates part of column (FTRANL) when sparse
void 
CoinBALPFactorization::updateColumnLSparse ( CoinIndexedVector * regionSparse ,
					   int * COIN_RESTRICT regionIndex)
  const
{
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  int number = regionSparse->getNumElements (  );
  int numberNonZero;
  double tolerance = zeroTolerance_;
  
  numberNonZero = 0;
  
  const CoinBigIndex *startColumn = startColumnL_.array();
  const int *indexRow = indexRowL_.array();
  const CoinFactorizationDouble *element = elementL_.array();
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
  for (int k=0;k<number;k++) {
    int kPivot=regionIndex[k];
    if (kPivot>=baseL_ && kPivot < numberColumns1_) {
      assert (kPivot<numberRowsExtra_);
      //if (kPivot>=numberRowsExtra_) abort();
      if(!mark[kPivot]) {
	stack[0]=kPivot;
	CoinBigIndex j=startColumn[kPivot+1]-1;
        int nStack=0;
	while (nStack>=0) {
	  /* take off stack */
	  if (j>=startColumn[kPivot]) {
	    int jPivot=indexRow[j--];
	    assert (jPivot>=baseL_&&jPivot<numberRowsExtra_);
	    //if (jPivot<baseL_||jPivot>=numberRowsExtra_) abort();
	    /* put back on stack */
	    next[nStack] =j;
	    if (!mark[jPivot]) {
	      /* and new one */
	      kPivot=jPivot;
	      j = startColumn[kPivot+1]-1;
	      stack[++nStack]=kPivot;
	      assert (kPivot<numberRowsExtra_);
	      //if (kPivot>=numberRowsExtra_) abort();
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
	      assert (kPivot<numberRowsExtra_);
	      j=next[nStack];
	    }
	  }
	}
      }
    } else {
      // just put on list
      regionIndex[numberNonZero++]=kPivot;
    }
  }
  for (int i=nList-1;i>=0;i--) {
    int iPivot = list[i];
    mark[iPivot]=0;
    CoinFactorizationDouble pivotValue = region[iPivot];
    if ( fabs ( pivotValue ) > tolerance ) {
      regionIndex[numberNonZero++]=iPivot;
      for (CoinBigIndex j = startColumn[iPivot]; 
	   j < startColumn[iPivot+1]; j ++ ) {
	int iRow = indexRow[j];
	CoinFactorizationDouble value = element[j];
	region[iRow] -= value * pivotValue;
      }
    } else {
      region[iPivot]=0.0;
    }
  }
  regionSparse->setNumElements ( numberNonZero );
}



//  updateColumnU.  Updates part of column (FTRANU)
void
CoinBALPFactorization::updateColumnU ( CoinIndexedVector * regionSparse,
				   int * indexIn) const
{
  int numberNonZero = regionSparse->getNumElements (  );
	if (collectStatistics_)
		ftranCountAfterR_ += numberNonZero;
	

  int goSparse;
  // Guess at number at end
  if (sparseThreshold_>0) {
      int newNumber = static_cast<int> (numberNonZero*ftranAverageAfterU_);
      if (newNumber< sparseThreshold_)
	goSparse = 2;
      else if (newNumber< sparseThreshold2_)
	goSparse = 1;
      else
	goSparse = 0;
    
  } else {
    goSparse=0;
  }
  // override because only densish solve works for now
  //goSparse = 0;
  switch (goSparse) {
  case 0: // densish
    {
      double *region = regionSparse->denseVector (  );
      int * regionIndex = regionSparse->getIndices();
      int numberNonZero=updateColumnUDensish(region,regionIndex);
      regionSparse->setNumElements ( numberNonZero );
    }
    break;
  case 1: // middling
    updateColumnUSparsish(regionSparse,indexIn);
    break;
  case 2: // sparse
    updateColumnUSparse(regionSparse,indexIn);
    break;
  }
  if (collectStatistics_) { 
    ftranCountAfterU_ += regionSparse->getNumElements (  );
    numberFtranCounts_++;
  }
}

// Updates part of column (FTRANU) real work
int 
CoinBALPFactorization::updateColumnUDensish ( double * COIN_RESTRICT region, 
					  int * COIN_RESTRICT regionIndex) const
{
  double tolerance = zeroTolerance_;
  const CoinBigIndex *startColumn = startColumnU_.array();
  const int *indexRow = indexRowU_.array();
  const CoinFactorizationDouble *element = elementU_.array();
  int numberNonZero = 0;
  const int *numberInColumn = numberInColumn_.array();
  const CoinFactorizationDouble *pivotRegion = pivotRegion_.array();

  
  for (int i = numberU_-1 ; i >= numberSlacks_; i-- ) {
    CoinFactorizationDouble pivotValue = region[i];
    if (pivotValue) {

      region[i] = 0.0;
      if ( fabs ( pivotValue ) > tolerance ) {
	CoinBigIndex start = startColumn[i];
	const CoinFactorizationDouble * thisElement = element+start;
	const int * thisIndex = indexRow+start;

	for (CoinBigIndex j=numberInColumn[i]-1 ; j >=0; j-- ) {
	  int iRow = thisIndex[j];
	  CoinFactorizationDouble regionValue = region[iRow];
	  CoinFactorizationDouble value = thisElement[j];
	  region[iRow] = regionValue - value * pivotValue;
	}
	pivotValue *= pivotRegion[i];
	region[i]=pivotValue;
	regionIndex[numberNonZero++]=i;
      }
    }
  }
    
  // now do slacks
  assert(numberSlacks_ == 0);
  if (slackValue_==-1.0) {
    for (int i = numberSlacks_-1; i>=0;i--) {
      double value = region[i];
      if ( value ) {
	region[i]=-value;
	regionIndex[numberNonZero]=i;
	if ( fabs(value) > tolerance ) 
	  numberNonZero++;
	else 
	  region[i]=0.0;
      }
    }
  } else {
    assert (slackValue_==1.0);
    for (int i = numberSlacks_-1; i>=0;i--) {
      double value = region[i];
      double absValue = fabs ( value );
      if ( value ) {
	region[i]=0.0;
	if ( absValue > tolerance ) {
	  region[i]=value;
	  regionIndex[numberNonZero++]=i;
	}
      }
    }
  }
  return numberNonZero;
}
//  updateColumnU.  Updates part of column (FTRANU)
/*
  Since everything is in order I should be able to do a better job of
  marking stuff - think.  Also as L is static maybe I can do something
  better there (I know I could if I marked the depth of every element
  but that would lead to other inefficiencies.
*/
void
CoinBALPFactorization::updateColumnUSparse ( CoinIndexedVector * regionSparse,
					 int * COIN_RESTRICT indexIn) const
{
  int numberNonZero = regionSparse->getNumElements (  );
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  double tolerance = zeroTolerance_;
  const CoinBigIndex *startColumn = startColumnU_.array();
  const int *indexRow = indexRowU_.array();
  const CoinFactorizationDouble *element = elementU_.array();
  const CoinFactorizationDouble *pivotRegion = pivotRegion_.array();
  // use sparse_ as temporary area
  // mark known to be zero
  int * COIN_RESTRICT stack = sparse_.array();  /* pivot */
  int * COIN_RESTRICT list = stack + maximumRowsExtra_;  /* final list */
  CoinBigIndex * COIN_RESTRICT next = reinterpret_cast<CoinBigIndex *> (list + maximumRowsExtra_);  /* jnext */
  char * COIN_RESTRICT mark = reinterpret_cast<char *> (next + maximumRowsExtra_);
#ifdef COIN_DEBUG
  for (int i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif

  // move slacks to end of stack list
  int * COIN_RESTRICT putLast = stack+maximumRowsExtra_;
  int * COIN_RESTRICT put = putLast;

  const int *numberInColumn = numberInColumn_.array();
  int nList = 0;
  for (int i=0;i<numberNonZero;i++) {
    int kPivot=indexIn[i];
    stack[0]=kPivot;
    CoinBigIndex j=startColumn[kPivot]+numberInColumn[kPivot]-1;
    int nStack=1;
    next[0]=j;
    while (nStack) {
      /* take off stack */
      int kPivot=stack[--nStack];
      if (mark[kPivot]!=1) {
	j=next[nStack];
	if (j>=startColumn[kPivot]) {
	  kPivot=indexRow[j--];
	  /* put back on stack */
	  next[nStack++] =j;
	  if (!mark[kPivot]) {
	    /* and new one */
	    int numberIn = numberInColumn[kPivot];
	    if (numberIn) {
	      j = startColumn[kPivot]+numberIn-1;
	      stack[nStack]=kPivot;
	      mark[kPivot]=2;
	      next[nStack++]=j;
	    } else {
	      // can do immediately
	      /* finished so mark */
	      mark[kPivot]=1;
	      if (kPivot>=numberSlacks_) {
		list[nList++]=kPivot;
	      } else {
		// slack - put at end
		--put;
		*put=kPivot;
	      }
	    }
	  }
	} else {
	  /* finished so mark */
	  mark[kPivot]=1;
	  if (kPivot>=numberSlacks_) {
	    list[nList++]=kPivot;
	  } else {
	    // slack - put at end
	    assert (!numberInColumn[kPivot]);
	    --put;
	    *put=kPivot;
	  }
	}
      }
    }
  }
#if 0
  {
    std::sort(list,list+nList);
    int i;
    int last;
    last =-1;
    for (i=0;i<nList;i++) {
      int k = list[i];
      assert (k>last);
      last=k;
    }
    std::sort(put,putLast);
    int n = putLast-put;
    last =-1;
    for (i=0;i<n;i++) {
      int k = put[i];
      assert (k>last);
      last=k;
    }
  }
#endif
  numberNonZero=0;
  for (int i=nList-1;i>=0;i--) {
    int iPivot = list[i];
    mark[iPivot]=0;
    CoinFactorizationDouble pivotValue = region[iPivot];
    region[iPivot]=0.0;
    if ( fabs ( pivotValue ) > tolerance ) {
      CoinBigIndex start = startColumn[iPivot];
      int number = numberInColumn[iPivot];
      
      CoinBigIndex j;
      for ( j = start; j < start+number; j++ ) {
	CoinFactorizationDouble value = element[j];
	int iRow = indexRow[j];
	region[iRow] -=  value * pivotValue;
      }
      pivotValue *= pivotRegion[iPivot];
      region[iPivot]=pivotValue;
      regionIndex[numberNonZero++]=iPivot;
    }
  }
  // slacks
#ifndef COIN_FAST_CODE
  if (slackValue_==1.0) {
    for (;put<putLast;put++) {
      int iPivot = *put;
      mark[iPivot]=0;
      CoinFactorizationDouble pivotValue = region[iPivot];
      region[iPivot]=0.0;
      if ( fabs ( pivotValue ) > tolerance ) {
	region[iPivot]=pivotValue;
	regionIndex[numberNonZero++]=iPivot;
      }
    }
  } else {
#endif
    for (;put<putLast;put++) {
      int iPivot = *put;
      mark[iPivot]=0;
      CoinFactorizationDouble pivotValue = region[iPivot];
      region[iPivot]=0.0;
      if ( fabs ( pivotValue ) > tolerance ) {
	region[iPivot]=-pivotValue;
	regionIndex[numberNonZero++]=iPivot;
      }
    }
#ifndef COIN_FAST_CODE
  }
#endif
  regionSparse->setNumElements ( numberNonZero );
}
//  updateColumnU.  Updates part of column (FTRANU)
/*
  Since everything is in order I should be able to do a better job of
  marking stuff - think.  Also as L is static maybe I can do something
  better there (I know I could if I marked the depth of every element
  but that would lead to other inefficiencies.
*/
#ifdef COIN_DEVELOP
double ncall_SZ=0.0;
double nrow_SZ=0.0;
double nslack_SZ=0.0;
double nU_SZ=0.0;
double nnz_SZ=0.0;
double nDone_SZ=0.0;
#endif
void
CoinBALPFactorization::updateColumnUSparsish ( CoinIndexedVector * regionSparse,
					   int * COIN_RESTRICT indexIn) const
{
  int * COIN_RESTRICT regionIndex = regionSparse->getIndices (  );
  // mark known to be zero
  int * COIN_RESTRICT stack = sparse_.array();  /* pivot */
  int * COIN_RESTRICT list = stack + maximumRowsExtra_;  /* final list */
  CoinBigIndex * COIN_RESTRICT next = reinterpret_cast<CoinBigIndex *> (list + maximumRowsExtra_);  /* jnext */
  CoinCheckZero * COIN_RESTRICT mark = reinterpret_cast<CoinCheckZero *> (next + maximumRowsExtra_);
  const int *numberInColumn = numberInColumn_.array();
#ifdef COIN_DEBUG
  for (int i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif

  int nMarked=0;
  int numberNonZero = regionSparse->getNumElements (  );
  double * COIN_RESTRICT region = regionSparse->denseVector (  );
  double tolerance = zeroTolerance_;
  const CoinBigIndex *startColumn = startColumnU_.array();
  const int *indexRow = indexRowU_.array();
  const CoinFactorizationDouble *element = elementU_.array();
  const CoinFactorizationDouble *pivotRegion = pivotRegion_.array();
#ifdef COIN_DEVELOP
  ncall_SZ++;
  nrow_SZ += numberRows_;
  nslack_SZ += numberSlacks_;
  nU_SZ += numberU_;
#endif

  for (int ii=0;ii<numberNonZero;ii++) {
    int iPivot=indexIn[ii];
    int iWord = iPivot>>CHECK_SHIFT;
    int iBit = iPivot-(iWord<<CHECK_SHIFT);
    if (mark[iWord]) {
      mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
    } else {
      mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
      stack[nMarked++]=iWord;
    }
  }
  numberNonZero = 0;
  // First do down to convenient power of 2
  CoinBigIndex jLast = (numberU_-1)>>CHECK_SHIFT;
  jLast = CoinMax((jLast<<CHECK_SHIFT),static_cast<CoinBigIndex> (numberSlacks_));
  int i;
  for ( i = numberU_-1 ; i >= jLast; i-- ) {
    CoinFactorizationDouble pivotValue = region[i];
    region[i] = 0.0;
    if ( fabs ( pivotValue ) > tolerance ) {
#ifdef COIN_DEVELOP
      nnz_SZ ++;
#endif
      CoinBigIndex start = startColumn[i];
      const CoinFactorizationDouble * thisElement = element+start;
      const int * thisIndex = indexRow+start;
      
#ifdef COIN_DEVELOP
      nDone_SZ += numberInColumn[i];
#endif
      for (int j=numberInColumn[i]-1 ; j >=0; j-- ) {
	int iRow0 = thisIndex[j];
	CoinFactorizationDouble regionValue0 = region[iRow0];
	CoinFactorizationDouble value0 = thisElement[j];
	int iWord = iRow0>>CHECK_SHIFT;
	int iBit = iRow0-(iWord<<CHECK_SHIFT);
	if (mark[iWord]) {
	  mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
	} else {
	  mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
	  stack[nMarked++]=iWord;
	}
	region[iRow0] = regionValue0 - value0 * pivotValue;
      }
      pivotValue *= pivotRegion[i];
      region[i]=pivotValue;
      regionIndex[numberNonZero++]=i;
    }
  }
  int kLast = (numberSlacks_+BITS_PER_CHECK-1)>>CHECK_SHIFT;
  if (jLast>numberSlacks_) {
    // now do in chunks
    for (int k=(jLast>>CHECK_SHIFT)-1;k>=kLast;k--) {
      unsigned int iMark = mark[k];
      if (iMark) {
	// something in chunk - do all (as imark may change)
	int iLast = k<<CHECK_SHIFT;
	for ( i = iLast+BITS_PER_CHECK-1 ; i >= iLast; i-- ) {
	  CoinFactorizationDouble pivotValue = region[i];
	  if (pivotValue) {
#ifdef COIN_DEVELOP
	    nnz_SZ ++;
#endif
	    region[i] = 0.0;
	    if ( fabs ( pivotValue ) > tolerance ) {
	      CoinBigIndex start = startColumn[i];
	      const CoinFactorizationDouble * thisElement = element+start;
	      const int * thisIndex = indexRow+start;
#ifdef COIN_DEVELOP
	      nDone_SZ += numberInColumn[i];
#endif
	      for (int j=numberInColumn[i]-1 ; j >=0; j-- ) {
		int iRow0 = thisIndex[j];
		CoinFactorizationDouble regionValue0 = region[iRow0];
		CoinFactorizationDouble value0 = thisElement[j];
		int iWord = iRow0>>CHECK_SHIFT;
		int iBit = iRow0-(iWord<<CHECK_SHIFT);
		if (mark[iWord]) {
		  mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
		} else {
		  mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
		  stack[nMarked++]=iWord;
		}
		region[iRow0] = regionValue0 - value0 * pivotValue;
	      }
	      pivotValue *= pivotRegion[i];
	      region[i]=pivotValue;
	      regionIndex[numberNonZero++]=i;
	    }
	  }
	}
	mark[k]=0;
      }
    }
    i = (kLast<<CHECK_SHIFT)-1;
  }
  for ( ; i >= numberSlacks_; i-- ) {
    CoinFactorizationDouble pivotValue = region[i];
    region[i] = 0.0;
    if ( fabs ( pivotValue ) > tolerance ) {
#ifdef COIN_DEVELOP
      nnz_SZ ++;
#endif
      CoinBigIndex start = startColumn[i];
      const CoinFactorizationDouble * thisElement = element+start;
      const int * thisIndex = indexRow+start;
#ifdef COIN_DEVELOP
      nDone_SZ += numberInColumn[i];
#endif
      for (int j=numberInColumn[i]-1 ; j >=0; j-- ) {
	int iRow0 = thisIndex[j];
	CoinFactorizationDouble regionValue0 = region[iRow0];
	CoinFactorizationDouble value0 = thisElement[j];
	int iWord = iRow0>>CHECK_SHIFT;
	int iBit = iRow0-(iWord<<CHECK_SHIFT);
	if (mark[iWord]) {
	  mark[iWord] = static_cast<CoinCheckZero>(mark[iWord] | (1<<iBit));
	} else {
	  mark[iWord] = static_cast<CoinCheckZero>(1<<iBit);
	  stack[nMarked++]=iWord;
	}
	region[iRow0] = regionValue0 - value0 * pivotValue;
      }
      pivotValue *= pivotRegion[i];
      region[i]=pivotValue;
      regionIndex[numberNonZero++]=i;
    }
  }
  
  if (numberSlacks_) {
    // now do slacks
#ifndef COIN_FAST_CODE
    double factor = slackValue_;
    if (factor==1.0) {
      // First do down to convenient power of 2
      CoinBigIndex jLast = (numberSlacks_-1)>>CHECK_SHIFT;
      jLast = jLast<<CHECK_SHIFT;
      for ( i = numberSlacks_-1; i>=jLast;i--) {
	double value = region[i];
	double absValue = fabs ( value );
	if ( value ) {
	  region[i]=0.0;
	  if ( absValue > tolerance ) {
	  region[i]=value;
	  regionIndex[numberNonZero++]=i;
	  }
	}
      }
      mark[jLast]=0;
      // now do in chunks
      for (int k=(jLast>>CHECK_SHIFT)-1;k>=0;k--) {
	unsigned int iMark = mark[k];
	if (iMark) {
	  // something in chunk - do all (as imark may change)
	  int iLast = k<<CHECK_SHIFT;
	  i = iLast+BITS_PER_CHECK-1;
	  for ( ; i >= iLast; i-- ) {
	    double value = region[i];
	    double absValue = fabs ( value );
	    if ( value ) {
	      region[i]=0.0;
	      if ( absValue > tolerance ) {
		region[i]=value;
		regionIndex[numberNonZero++]=i;
	      }
	    }
	  }
	  mark[k]=0;
	}
      }
    } else {
      assert (factor==-1.0);
#endif
      // First do down to convenient power of 2
      CoinBigIndex jLast = (numberSlacks_-1)>>CHECK_SHIFT;
      jLast = jLast<<CHECK_SHIFT;
      for ( i = numberSlacks_-1; i>=jLast;i--) {
	double value = region[i];
	double absValue = fabs ( value );
	if ( value ) {
	  region[i]=0.0;
	  if ( absValue > tolerance ) {
	    region[i]=-value;
	    regionIndex[numberNonZero++]=i;
	  }
	}
      }
      mark[jLast]=0;
      // now do in chunks
      for (int k=(jLast>>CHECK_SHIFT)-1;k>=0;k--) {
	unsigned int iMark = mark[k];
	if (iMark) {
	  // something in chunk - do all (as imark may change)
	  int iLast = k<<CHECK_SHIFT;
	  i = iLast+BITS_PER_CHECK-1;
	  for ( ; i >= iLast; i-- ) {
	    double value = region[i];
	    double absValue = fabs ( value );
	    if ( value ) {
	      region[i]=0.0;
	      if ( absValue > tolerance ) {
		region[i]=-value;
		regionIndex[numberNonZero++]=i;
	      }
	    }
	  }
	  mark[k]=0;
	}
      }
#ifndef COIN_FAST_CODE
    }
#endif
  }
  regionSparse->setNumElements ( numberNonZero );
  mark[(numberU_-1)>>CHECK_SHIFT]=0;
  mark[numberSlacks_>>CHECK_SHIFT]=0;
  if (numberSlacks_)
    mark[(numberSlacks_-1)>>CHECK_SHIFT]=0;
#ifdef COIN_DEBUG
  for (i=0;i<maximumRowsExtra_;i++) {
    assert (!mark[i]);
  }
#endif
}


