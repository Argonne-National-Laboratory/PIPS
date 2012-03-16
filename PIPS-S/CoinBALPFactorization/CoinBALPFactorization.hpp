#ifndef COINBALP_HPP
#define COINBALP_HPP

#include "CoinFactorization.hpp"

// modification of CoinFactorization to perform "tall" LU and handle first-stage block 
// intended as a stand-alone code, doesn't depend on any code from simplex project
// idea is to keep modifications of CoinUtils separate so CoinUtils can be updated
// however, because we access internals, changes in CoinUtils could break this code.
// designed with CoinUtils-2.8.0


class CoinBALPFactorization : protected CoinFactorization {
// protected inheritance so we don't expose any functions that
// might not have the desired effect

public:
	CoinBALPFactorization();
	~CoinBALPFactorization();

	void show_self() const;

	/** When part of LP - given by basic variables.
	Actually does factorization.
	Arrays passed in have non negative value to say basic.
	Basic row refers to corresponding slack variable.
	Takes rectangular matrix in two parts: [ M1 M2 ]
	Only eliminates lower triangle of M1, but applies operations
	to M2. Column permutations are restricted to M1.
	Stores result in T matrix
	returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
	int factorizeRect ( const CoinPackedMatrix & M1,
		  const CoinPackedMatrix & M2,
		  int *M1columnIsBasic,
		  int *M1rowIsBasic,
		  int *M2columnIsBasic,
		  int M2extracolumns,
		  double areaFactor = 0.0 );

	// like CoinFactorization, but keep the permutations internal
	// matrix given as unsorted triplets
	int factorizeSquare(int numberOfColumns, CoinBigIndex numberOfElements, const int indicesRow[],
				const int indicesColumn[], const double elements[], double areaFactor = 0.0);
	

	void getRowOfT(CoinIndexedVector&,int row);
	void getRowOfT(int *idx,double *elts, int row);
	int getNumberInRowOfT(int row) const;
	
	// multiply a vector times the first numberColumns_ rows of T
	// and *subtract* it from out vector
	void multXT(const CoinIndexedVector &in, CoinIndexedVector &out);
	void multXTTranspose(const CoinIndexedVector &in, CoinIndexedVector &out);
	void multXTTranspose(const CoinIndexedVector &in, double *out);

	// FTRAN
	int updateColumn(CoinIndexedVector * regionSparse, 
			  CoinIndexedVector * regionSparse2) const;

	/// Updates part of column (FTRANL)
	void updateColumnL ( CoinIndexedVector * region, int * indexIn ) const;

	/// Updates part of column (FTRANU)
	void updateColumnU ( CoinIndexedVector * region, int * indexIn) const;

	// applies L^{-1}P (FTRANG)
	// input vector is "vec"
	// output is in "region"
	void updateColumnG (CoinIndexedVector *region, CoinIndexedVector*vec) const;

	// applies QU^{-1}
	// input vector is "region"
	// output vector is "vec"
	void updateColumnUQ(CoinIndexedVector *region, CoinIndexedVector*vec) const;

	// applies P^{-1}L^{-T} (BTRANG)
	// input vector is "region"
	// output is in "vec"
	void updateColumnTransposeG (CoinIndexedVector *region, CoinIndexedVector*vec) const;

	// applies U^{-T}Q^{-1}
	// input vector is "vec"
	// output vector is "region"
	void updateColumnTransposeUQ(CoinIndexedVector *region, CoinIndexedVector*vec) const;
	
	/** Updates one column (BTRAN) from regionSparse2
	regionSparse starts as zero and is zero at end 
	Note - if regionSparse2 packed on input - will be packed on output
	*/
	int updateColumnTranspose ( CoinIndexedVector * regionSparse,
			      CoinIndexedVector * regionSparse2) const;

	/** Updates part of column transpose (BTRANU),
	assumes index is sorted i.e. region is correct */	
	void updateColumnTransposeU ( CoinIndexedVector * region,
				int smallestIndex) const;
	
	/// Updates part of column transpose (BTRANL)
	void updateColumnTransposeL ( CoinIndexedVector * region ) const;

	int numberColumns() const { return CoinFactorization::numberColumns(); }
	
	// make a row copy of L
	void goSparse() { sparseThreshold(0); CoinFactorization::goSparse(); }
	// these should be called also for hyper-sparse solves
	void checkSparse() { CoinFactorization::checkSparse(); }
	void setCollectStatistics(bool b) { CoinFactorization::setCollectStatistics(b); }

	// this is only a hint to avoid allocating too much space for factors
	// if we have many scenarios
	void setNumScenariosPerProc(int n) { scenariosPerProc_ = n; }

protected:
	void preProcessRect(int state);
	int factorRect();
	int factorSparseRect();
	int factorSparseSmallRect();
	int factorSparseLargeRect();
	void cleanupRect();

	void gutsOfDestructor();
	//void gutsOfInitialize(int);
	void getAreas(int nrows, int ncols1, int ncols2, CoinBigIndex, CoinBigIndex);

	bool pivotColumnSingleton ( int pivotRow, int pivotColumn );
	bool pivotOneOtherRow ( int pivotRow, int pivotColumn );

	
	/// Updates part of column (FTRANL) when densish
	void updateColumnLDensish ( CoinIndexedVector * region, int * indexIn ) const;
	/// Updates part of column (FTRANL) when sparse
	void updateColumnLSparse ( CoinIndexedVector * region, int * indexIn ) const;
	/// Updates part of column (FTRANL) when sparsish
	void updateColumnLSparsish ( CoinIndexedVector * region, int * indexIn ) const;

	
	/// Updates part of column (FTRANU) when sparse
	void updateColumnUSparse ( CoinIndexedVector * regionSparse, 
			     int * indexIn) const;
	/// Updates part of column (FTRANU) when sparsish
	void updateColumnUSparsish ( CoinIndexedVector * regionSparse, 
			       int * indexIn) const;
	/// Updates part of column (FTRANU)
	int updateColumnUDensish ( double * COIN_RESTRICT region, 
			     int * COIN_RESTRICT regionIndex) const;
	
	
	/** Updates part of column transpose (BTRANU) when sparsish,
	assumes index is sorted i.e. region is correct */
	void updateColumnTransposeUSparsish ( CoinIndexedVector * region,
					int smallestIndex) const;
	/** Updates part of column transpose (BTRANU) when densish,
	assumes index is sorted i.e. region is correct */
	void updateColumnTransposeUDensish ( CoinIndexedVector * region,
				       int smallestIndex) const;
	/** Updates part of column transpose (BTRANU) when sparse,
	assumes index is sorted i.e. region is correct */
	void updateColumnTransposeUSparse ( CoinIndexedVector * region) const;

	
	/// Updates part of column transpose (BTRANL) when densish by column
	void updateColumnTransposeLDensish ( CoinIndexedVector * region ) const;
	/// Updates part of column transpose (BTRANL) when densish by row
	void updateColumnTransposeLByRow ( CoinIndexedVector * region ) const;
	/// Updates part of column transpose (BTRANL) when sparsish by row
	void updateColumnTransposeLSparsish ( CoinIndexedVector * region ) const;
	/// Updates part of column transpose (BTRANL) when sparse (by Row)
	void updateColumnTransposeLSparse ( CoinIndexedVector * region ) const;



	
	// number of columns in first block
	// this is equal to numberColumns_ after factorization
	int numberColumns1_;
	
	// number of columns in T block
	int numberColumnsT_;

	/// For rectangular LU, flag which rows were pivoted on
	// 0 - no, 1 - yes
	CoinIntArrayWithLength rowPivotFlag_;

	// row copy
	CoinBigIndexArrayWithLength startRowT_;	
	
	/// Length of T
	CoinBigIndex lengthT_;

	/// Length of area reserved for T
	CoinBigIndex lengthAreaT_;

	/// Elements of T
	CoinFactorizationDoubleArrayWithLength elementT_;

	/// Row indices of T
	CoinIntArrayWithLength indexRowT_;

	/// Start of each column in T
	CoinBigIndexArrayWithLength startColumnT_;

	/// Converts rows to columns in T
	// i.e. map from indices in row copy to index in elementsU
	CoinBigIndexArrayWithLength convertRowToColumnT_;

	CoinIntArrayWithLength indexColumnT_;

	/// Number in each Row of T
	CoinIntArrayWithLength numberInRowT_;

	/// Number in each Column of T
	CoinIntArrayWithLength numberInColumnT_;

	/// Hint to help determine how much space to allocate for the factors
	int scenariosPerProc_;



  template <class T>  inline bool
  pivot ( int pivotRow,
	  int pivotColumn,
	  CoinBigIndex pivotRowPosition,
	  CoinBigIndex pivotColumnPosition,
	  CoinFactorizationDouble work[],
	  unsigned int workArea2[],
	  int increment2,
	  T markRow[] ,
	  int largeInteger)
{
  int *indexColumnU = indexColumnU_.array();
  CoinBigIndex *startColumnU = startColumnU_.array();
  int *numberInColumn = numberInColumn_.array();
  CoinFactorizationDouble *elementU = elementU_.array();
  int *indexRowU = indexRowU_.array();
  CoinBigIndex *startRowU = startRowU_.array();
  int *numberInRow = numberInRow_.array();
  CoinFactorizationDouble *elementL = elementL_.array();
  int *indexRowL = indexRowL_.array();
  int *saveColumn = saveColumn_.array();
  int *nextRow = nextRow_.array();
  int *lastRow = lastRow_.array() ;

  //store pivot columns (so can easily compress)
  int numberInPivotRow = numberInRow[pivotRow] - 1;
  CoinBigIndex startColumn = startColumnU[pivotColumn];
  int numberInPivotColumn = numberInColumn[pivotColumn] - 1;
  CoinBigIndex endColumn = startColumn + numberInPivotColumn + 1;
  int put = 0;
  CoinBigIndex startRow = startRowU[pivotRow];
  CoinBigIndex endRow = startRow + numberInPivotRow + 1;

  if ( pivotColumnPosition < 0 ) {
    for ( pivotColumnPosition = startRow; pivotColumnPosition < endRow; pivotColumnPosition++ ) {
      int iColumn = indexColumnU[pivotColumnPosition];
      if ( iColumn != pivotColumn ) {
	saveColumn[put++] = iColumn;
      } else {
        break;
      }
    }
  } else {
    for (CoinBigIndex i = startRow ; i < pivotColumnPosition ; i++ ) {
      saveColumn[put++] = indexColumnU[i];
    }
  }
  assert (pivotColumnPosition<endRow);
  assert (indexColumnU[pivotColumnPosition]==pivotColumn);
  pivotColumnPosition++;
  for ( ; pivotColumnPosition < endRow; pivotColumnPosition++ ) {
    saveColumn[put++] = indexColumnU[pivotColumnPosition];
  }
  //take out this bit of indexColumnU
  int next = nextRow[pivotRow];
  int last = lastRow[pivotRow];

  nextRow[last] = next;
  lastRow[next] = last;
  nextRow[pivotRow] = numberGoodU_;	//use for permute
  lastRow[pivotRow] = -2;
  numberInRow[pivotRow] = 0;
  //store column in L, compress in U and take column out
  CoinBigIndex l = lengthL_;

  if ( l + numberInPivotColumn > lengthAreaL_ ) {
    //need more memory
    if ((messageLevel_&4)!=0) 
      printf("more memory needed in middle of invert\n");
    return false;
  }
  //l+=currentAreaL_->elementByColumn-elementL;
  CoinBigIndex lSave = l;

  CoinBigIndex * startColumnL = startColumnL_.array();
  startColumnL[numberGoodL_] = l;	//for luck and first time
  numberGoodL_++;
  startColumnL[numberGoodL_] = l + numberInPivotColumn;
  lengthL_ += numberInPivotColumn;
  if ( pivotRowPosition < 0 ) {
    for ( pivotRowPosition = startColumn; pivotRowPosition < endColumn; pivotRowPosition++ ) {
      int iRow = indexRowU[pivotRowPosition];
      if ( iRow != pivotRow ) {
	indexRowL[l] = iRow;
	elementL[l] = elementU[pivotRowPosition];
	markRow[iRow] = static_cast<T>(l - lSave);
	l++;
	//take out of row list
	CoinBigIndex start = startRowU[iRow];
	CoinBigIndex end = start + numberInRow[iRow];
	CoinBigIndex where = start;

	while ( indexColumnU[where] != pivotColumn ) {
	  where++;
	}			/* endwhile */
#if DEBUG_COIN
	if ( where >= end ) {
	  abort (  );
	}
#endif
	indexColumnU[where] = indexColumnU[end - 1];
	numberInRow[iRow]--;
      } else {
	break;
      }
    }
  } else {
    CoinBigIndex i;

    for ( i = startColumn; i < pivotRowPosition; i++ ) {
      int iRow = indexRowU[i];

      markRow[iRow] = static_cast<T>(l - lSave);
      indexRowL[l] = iRow;
      elementL[l] = elementU[i];
      l++;
      //take out of row list
      CoinBigIndex start = startRowU[iRow];
      CoinBigIndex end = start + numberInRow[iRow];
      CoinBigIndex where = start;

      while ( indexColumnU[where] != pivotColumn ) {
	where++;
      }				/* endwhile */
#if DEBUG_COIN
      if ( where >= end ) {
	abort (  );
      }
#endif
      indexColumnU[where] = indexColumnU[end - 1];
      numberInRow[iRow]--;
      assert (numberInRow[iRow]>=0);
    }
  }
  assert (pivotRowPosition<endColumn);
  assert (indexRowU[pivotRowPosition]==pivotRow);
  CoinFactorizationDouble pivotElement = elementU[pivotRowPosition];
  CoinFactorizationDouble pivotMultiplier = 1.0 / pivotElement;

  pivotRegion_.array()[numberGoodU_] = pivotMultiplier;
  pivotRowPosition++;
  for ( ; pivotRowPosition < endColumn; pivotRowPosition++ ) {
    int iRow = indexRowU[pivotRowPosition];
    
    markRow[iRow] = static_cast<T>(l - lSave);
    indexRowL[l] = iRow;
    elementL[l] = elementU[pivotRowPosition];
    l++;
    //take out of row list
    CoinBigIndex start = startRowU[iRow];
    CoinBigIndex end = start + numberInRow[iRow];
    CoinBigIndex where = start;
    
    while ( indexColumnU[where] != pivotColumn ) {
      where++;
    }				/* endwhile */
#if DEBUG_COIN
    if ( where >= end ) {
      abort (  );
    }
#endif
    indexColumnU[where] = indexColumnU[end - 1];
    numberInRow[iRow]--;
    assert (numberInRow[iRow]>=0);
  }
  markRow[pivotRow] = static_cast<T>(largeInteger);
  //compress pivot column (move pivot to front including saved)
  numberInColumn[pivotColumn] = 0;
  //use end of L for temporary space
  int *indexL = &indexRowL[lSave];
  CoinFactorizationDouble *multipliersL = &elementL[lSave];

  //adjust
  int j;

  for ( j = 0; j < numberInPivotColumn; j++ ) {
    multipliersL[j] *= pivotMultiplier;
  }
  //zero out fill
  CoinBigIndex iErase;
  for ( iErase = 0; iErase < increment2 * numberInPivotRow;
	iErase++ ) {
    workArea2[iErase] = 0;
  }
  CoinBigIndex added = numberInPivotRow * numberInPivotColumn;
  unsigned int *temp2 = workArea2;
  int * nextColumn = nextColumn_.array();

  //pack down and move to work
  int jColumn;
  for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
    int iColumn = saveColumn[jColumn];
    CoinBigIndex startColumn = startColumnU[iColumn];
    CoinBigIndex endColumn = startColumn + numberInColumn[iColumn];
    int iRow = indexRowU[startColumn];
    CoinFactorizationDouble value = elementU[startColumn];
    double largest;
    CoinBigIndex put = startColumn;
    CoinBigIndex positionLargest = -1;
    CoinFactorizationDouble thisPivotValue = 0.0;

    //compress column and find largest not updated
    bool checkLargest;
    int mark = markRow[iRow];

    if ( mark == largeInteger+1 ) {
      largest = fabs ( value );
      positionLargest = put;
      put++;
      checkLargest = false;
    } else {
      //need to find largest
      largest = 0.0;
      checkLargest = true;
      if ( mark != largeInteger ) {
	//will be updated
	work[mark] = value;
	int word = mark >> COINFACTORIZATION_SHIFT_PER_INT;
	int bit = mark & COINFACTORIZATION_MASK_PER_INT;

	temp2[word] = temp2[word] | ( 1 << bit );	//say already in counts
	added--;
      } else {
	thisPivotValue = value;
      }
    }
    CoinBigIndex i;
    for ( i = startColumn + 1; i < endColumn; i++ ) {
      iRow = indexRowU[i];
      value = elementU[i];
      int mark = markRow[iRow];
      assert(iRow >= 0 && iRow < numberRows_);

      if ( mark == largeInteger+1 ) {
	//keep
	indexRowU[put] = iRow;
	elementU[put] = value;
	if ( checkLargest ) {
	  double absValue = fabs ( value );

	  if ( absValue > largest ) {
	    largest = absValue;
	    positionLargest = put;
	  }
	}
	put++;
      } else if ( mark != largeInteger ) {
	//will be updated
	work[mark] = value;
	int word = mark >> COINFACTORIZATION_SHIFT_PER_INT;
	int bit = mark & COINFACTORIZATION_MASK_PER_INT;

	temp2[word] = temp2[word] | ( 1 << bit );	//say already in counts
	added--;
      } else {
	thisPivotValue = value;
      }
    }
    //slot in pivot
    elementU[put] = elementU[startColumn];
    indexRowU[put] = indexRowU[startColumn];
    if ( positionLargest == startColumn ) {
      positionLargest = put;	//follow if was largest
    }
    put++;
    elementU[startColumn] = thisPivotValue;
    indexRowU[startColumn] = pivotRow;
    //clean up counts
    startColumn++;
    numberInColumn[iColumn] = put - startColumn;
    int * numberInColumnPlus = numberInColumnPlus_.array();
    numberInColumnPlus[iColumn]++;
    startColumnU[iColumn]++;
    //how much space have we got
    int next = nextColumn[iColumn];
    CoinBigIndex space;

    space = startColumnU[next] - put - numberInColumnPlus[next];
    //assume no zero elements
    if ( numberInPivotColumn > space ) {
      //getColumnSpace also moves fixed part
      if ( !getColumnSpace ( iColumn, numberInPivotColumn ) ) {
	return false;
      }
      //redo starts
      // added positionLargest >=0 condition
      // seems like -1 indicates no largest, so keep this indication
      // was causing memory corruption (very nasty to debug) - Miles
      if (positionLargest>=0) positionLargest = positionLargest + startColumnU[iColumn] - startColumn;
      startColumn = startColumnU[iColumn];
      put = startColumn + numberInColumn[iColumn];
    }
    double tolerance = zeroTolerance_;

    int *nextCount = nextCount_.array();
    for ( j = 0; j < numberInPivotColumn; j++ ) {
      value = work[j] - thisPivotValue * multipliersL[j];
      double absValue = fabs ( value );

      if ( absValue > tolerance ) {
	work[j] = 0.0;
	assert (put<lengthAreaU_); 
	elementU[put] = value;
	indexRowU[put] = indexL[j];
	assert(indexRowU[put] >= 0 && indexRowU[put] < numberRows_);
	if ( absValue > largest ) {
	  largest = absValue;
	  positionLargest = put;
	}
	put++;
      } else {
	work[j] = 0.0;
	added--;
	int word = j >> COINFACTORIZATION_SHIFT_PER_INT;
	int bit = j & COINFACTORIZATION_MASK_PER_INT;

	if ( temp2[word] & ( 1 << bit ) ) {
	  //take out of row list
	  iRow = indexL[j];
	  CoinBigIndex start = startRowU[iRow];
	  CoinBigIndex end = start + numberInRow[iRow];
	  CoinBigIndex where = start;

	  while ( indexColumnU[where] != iColumn ) {
	    where++;
	  }			/* endwhile */
#if DEBUG_COIN
	  if ( where >= end ) {
	    abort (  );
	  }
#endif
	  indexColumnU[where] = indexColumnU[end - 1];
	  numberInRow[iRow]--;
	} else {
	  //make sure won't be added
	  int word = j >> COINFACTORIZATION_SHIFT_PER_INT;
	  int bit = j & COINFACTORIZATION_MASK_PER_INT;

	  temp2[word] = temp2[word] | ( 1 << bit );	//say already in counts
	}
      }
    }
    numberInColumn[iColumn] = put - startColumn;
    //move largest
    // without earlier change, this can trigger with numberInColumn[iColumn] == 0, not good
    if ( positionLargest >= 0 ) {
      value = elementU[positionLargest];
      iRow = indexRowU[positionLargest];
      elementU[positionLargest] = elementU[startColumn];
      indexRowU[positionLargest] = indexRowU[startColumn];
      elementU[startColumn] = value;
      indexRowU[startColumn] = iRow;
    }
    //linked list for column
    if ( iColumn < numberColumns1_ && nextCount[iColumn + numberRows_] != -2 ) {
      //modify linked list
      deleteLink ( iColumn + numberRows_ );
      addLink ( iColumn + numberRows_, numberInColumn[iColumn] );
    }
    temp2 += increment2;
  }
  //get space for row list
  unsigned int *putBase = workArea2;
  int bigLoops = numberInPivotColumn >> COINFACTORIZATION_SHIFT_PER_INT;
  int i = 0;

  // do linked lists and update counts
  while ( bigLoops ) {
    bigLoops--;
    int bit;
    for ( bit = 0; bit < COINFACTORIZATION_BITS_PER_INT; i++, bit++ ) {
      unsigned int *putThis = putBase;
      int iRow = indexL[i];

      //get space
      int number = 0;
      int jColumn;

      for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
	unsigned int test = *putThis;

	putThis += increment2;
	test = 1 - ( ( test >> bit ) & 1 );
	number += test;
      }
      int next = nextRow[iRow];
      CoinBigIndex space;

      space = startRowU[next] - startRowU[iRow];
      number += numberInRow[iRow];
      if ( space < number ) {
	if ( !getRowSpace ( iRow, number ) ) {
	  return false;
	}
      }
      // now do
      putThis = putBase;
      next = nextRow[iRow];
      number = numberInRow[iRow];
      CoinBigIndex end = startRowU[iRow] + number;
      int saveIndex = indexColumnU[startRowU[next]];

      //add in
      for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
	unsigned int test = *putThis;

	putThis += increment2;
	test = 1 - ( ( test >> bit ) & 1 );
	indexColumnU[end] = saveColumn[jColumn];
	end += test;
      }
      //put back next one in case zapped
      indexColumnU[startRowU[next]] = saveIndex;
      markRow[iRow] = static_cast<T>(largeInteger+1);
      number = end - startRowU[iRow];
      numberInRow[iRow] = number;
      deleteLink ( iRow );
      addLink ( iRow, number );
    }
    putBase++;
  }				/* endwhile */
  int bit;

  for ( bit = 0; i < numberInPivotColumn; i++, bit++ ) {
    unsigned int *putThis = putBase;
    int iRow = indexL[i];

    //get space
    int number = 0;
    int jColumn;

    for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
      unsigned int test = *putThis;

      putThis += increment2;
      test = 1 - ( ( test >> bit ) & 1 );
      number += test;
    }
    int next = nextRow[iRow];
    CoinBigIndex space;

    space = startRowU[next] - startRowU[iRow];
    number += numberInRow[iRow];
    if ( space < number ) {
      if ( !getRowSpace ( iRow, number ) ) {
	return false;
      }
    }
    // now do
    putThis = putBase;
    next = nextRow[iRow];
    number = numberInRow[iRow];
    CoinBigIndex end = startRowU[iRow] + number;
    int saveIndex;

    saveIndex = indexColumnU[startRowU[next]];

    //add in
    for ( jColumn = 0; jColumn < numberInPivotRow; jColumn++ ) {
      unsigned int test = *putThis;

      putThis += increment2;
      test = 1 - ( ( test >> bit ) & 1 );

      indexColumnU[end] = saveColumn[jColumn];
      end += test;
    }
    indexColumnU[startRowU[next]] = saveIndex;
    markRow[iRow] = static_cast<T>(largeInteger+1);
    number = end - startRowU[iRow];
    numberInRow[iRow] = number;
    deleteLink ( iRow );
    addLink ( iRow, number );
  }
  markRow[pivotRow] = static_cast<T>(largeInteger+1);
  //modify linked list for pivots
  deleteLink ( pivotRow );
  deleteLink ( pivotColumn + numberRows_ );
  totalElements_ += added;
  return true;
}



};



#endif
