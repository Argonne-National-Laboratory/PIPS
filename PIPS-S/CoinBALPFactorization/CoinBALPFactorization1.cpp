#include "CoinBALPFactorization.hpp"

#include "PIPSLogging.hpp"

// changes to existing CoinFactorization auxiliary functions

CoinBALPFactorization::CoinBALPFactorization() {
	// parent constructor calls its version of gutsOfInitialize
	// just do what's missing
	// would be more convenient if CoinFactorization functions were virtual
	//rowPivotFlag_.conditionalNew(1);
	scenariosPerProc_ = 1;
}

CoinBALPFactorization::~CoinBALPFactorization() {
	rowPivotFlag_.conditionalDelete();
}

void CoinBALPFactorization::gutsOfDestructor() {
	CoinFactorization::gutsOfDestructor();
	rowPivotFlag_.conditionalDelete();

	numberColumns1_ = 0;
	numberColumnsT_ = 0;
	lengthT_ = 0;
	lengthAreaT_ = 0;
	elementT_.conditionalDelete();
	indexRowT_.conditionalDelete();
	indexColumnT_.conditionalDelete();
	//startRowU_.conditionalDelete();
	//startColumnU_.conditionalDelete();
	numberInColumnT_.conditionalDelete();
	numberInRowT_.conditionalDelete();
	convertRowToColumnT_.conditionalDelete();
	startColumnT_.conditionalDelete();
	startRowT_.conditionalDelete();

}

/*
void CoinBALPFactorization::gutsOfInitialize(int type) {
	CoinFactorization::gutsOfInitialize(type);
	if ( (type&4) != 0 ) {
		rowPivotFlag_.conditionalNew(1);
	}
}*/

void CoinBALPFactorization::getAreas(int numberOfRows, int numberOfColumns1, 
	int numberOfColumns2, CoinBigIndex maximumL, CoinBigIndex maximumU) {

	CoinFactorization::getAreas(numberOfRows,numberOfColumns1+numberOfColumns2,maximumL,maximumU);
	rowPivotFlag_.conditionalNew( numberRows_);
	memset(rowPivotFlag_.array(),0,numberRows_*sizeof(int));

	numberColumnsT_ = numberOfColumns2;
	// calculate this exactly when we're done factorizing
	//lengthAreaT_ = static_cast<double>(maximumU*numberOfColumns2)/(numberOfColumns1+numberOfColumns2); // no specific reason to choose this

	//elementT_.conditionalNew( lengthAreaT_ );
	//indexRowT_.conditionalNew( lengthAreaT_ );
	//indexColumnT_.conditionalNew( lengthAreaT_ );
	//convertRowToColumnT_.conditionalNew(lengthAreaT_);
	numberInColumnT_.conditionalNew( maximumColumnsExtra_ + 1 );
	numberInRowT_.conditionalNew( maximumRowsExtra_ + 1 );
	startColumnT_.conditionalNew( maximumColumnsExtra_ + 1 );
	startRowT_.conditionalNew( maximumRowsExtra_ + 1 );

}


bool
CoinBALPFactorization::pivotColumnSingleton ( int pivotRow,
				     int pivotColumn )
{
  int * numberInRow = numberInRow_.array();
  int * numberInColumn = numberInColumn_.array();
  int * numberInColumnPlus = numberInColumnPlus_.array();
  //store pivot columns (so can easily compress)
  int numberDoRow = numberInRow[pivotRow] - 1;
  CoinBigIndex * startColumnU = startColumnU_.array();
  CoinBigIndex startColumn = startColumnU[pivotColumn];
  int put = 0;
  CoinBigIndex * startRowU = startRowU_.array();
  CoinBigIndex startRow = startRowU[pivotRow];
  CoinBigIndex endRow = startRow + numberDoRow + 1;
  int * indexColumnU = indexColumnU_.array();
  int * indexRowU = indexRowU_.array();
  int * saveColumn = saveColumn_.array();
  CoinBigIndex i;

  for ( i = startRow; i < endRow; i++ ) {
    int iColumn = indexColumnU[i];

    if ( iColumn != pivotColumn ) {
      saveColumn[put++] = iColumn;
    }
  }
  int * nextRow = nextRow_.array();
  int * lastRow = lastRow_.array();
  //take out this bit of indexColumnU
  int next = nextRow[pivotRow];
  int last = lastRow[pivotRow];

  nextRow[last] = next;
  lastRow[next] = last;
  nextRow[pivotRow] = numberGoodU_;	//use for permute
  lastRow[pivotRow] =-2; //mark
  //clean up counts
  CoinFactorizationDouble *elementU = elementU_.array();
  CoinFactorizationDouble pivotElement = elementU[startColumn];

  pivotRegion_.array()[numberGoodU_] = 1.0 / pivotElement;
  numberInColumn[pivotColumn] = 0;
  //totalElements_ --;
  //numberInColumnPlus[pivotColumn]++; // don't add b/c pivot elements don't stay in U
  //move pivot row in other columns to safe zone
  for ( i = 0; i < numberDoRow; i++ ) {
    int iColumn = saveColumn[i];

    if ( numberInColumn[iColumn] ) {
      int number = numberInColumn[iColumn] - 1;

			if (iColumn < numberColumns1_) {
				//modify linked list
				deleteLink ( iColumn + numberRows_ );
				addLink ( iColumn + numberRows_, number );
      }
			//move pivot row element
      if ( number ) {
	CoinBigIndex start = startColumnU[iColumn];
	CoinBigIndex pivot = start;
	int iRow = indexRowU[pivot];
	while ( iRow != pivotRow ) {
	  pivot++;
	  iRow = indexRowU[pivot];
	}
        assert (pivot < startColumnU[iColumn] +
                numberInColumn[iColumn]);
	if ( pivot != start ) {
	  //move largest one up
	  CoinFactorizationDouble value = elementU[start];

	  iRow = indexRowU[start];
	  elementU[start] = elementU[pivot];
	  indexRowU[start] = indexRowU[pivot];
	  elementU[pivot] = elementU[start + 1];
	  indexRowU[pivot] = indexRowU[start + 1];
	  elementU[start + 1] = value;
	  indexRowU[start + 1] = iRow;
	} else {
	  //find new largest element
	  int iRowSave = indexRowU[start + 1];
	  CoinFactorizationDouble valueSave = elementU[start + 1];
	  double valueLargest = fabs ( valueSave );
	  CoinBigIndex end = start + numberInColumn[iColumn];
	  CoinBigIndex largest = start + 1;

	  CoinBigIndex k;
	  for ( k = start + 2; k < end; k++ ) {
	    CoinFactorizationDouble value = elementU[k];
	    double valueAbs = fabs ( value );

	    if ( valueAbs > valueLargest ) {
	      valueLargest = valueAbs;
	      largest = k;
	    }
	  }
	  indexRowU[start + 1] = indexRowU[largest];
	  elementU[start + 1] = elementU[largest];
	  indexRowU[largest] = iRowSave;
	  elementU[largest] = valueSave;
	}
      }
      //clean up counts
      numberInColumn[iColumn]--;
      numberInColumnPlus[iColumn]++;
      startColumnU[iColumn]++;
      //totalElements_--;
    }
  }
  //modify linked list for pivots
  deleteLink ( pivotRow );
  deleteLink ( pivotColumn + numberRows_ );
  numberInRow[pivotRow] = 0;
  //put in dummy pivot in L
  CoinBigIndex l = lengthL_;

  CoinBigIndex * startColumnL = startColumnL_.array();
  startColumnL[numberGoodL_] = l;	//for luck and first time
  numberGoodL_++;
  startColumnL[numberGoodL_] = l;
  return true;
}


bool 
CoinBALPFactorization::pivotOneOtherRow ( int pivotRow,
					   int pivotColumn )
{
  int * numberInRow = numberInRow_.array();
  int * numberInColumn = numberInColumn_.array();
  int * numberInColumnPlus = numberInColumnPlus_.array();
  int numberInPivotRow = numberInRow[pivotRow] - 1;
  CoinBigIndex * startRowU = startRowU_.array();
  CoinBigIndex * startColumnU = startColumnU_.array();
  CoinBigIndex startColumn = startColumnU[pivotColumn];
  CoinBigIndex startRow = startRowU[pivotRow];
  CoinBigIndex endRow = startRow + numberInPivotRow + 1;

  //take out this bit of indexColumnU
  int * nextRow = nextRow_.array();
  int * lastRow = lastRow_.array();
  int next = nextRow[pivotRow];
  int last = lastRow[pivotRow];

  nextRow[last] = next;
  lastRow[next] = last;
  nextRow[pivotRow] = numberGoodU_;	//use for permute
  lastRow[pivotRow] = -2;
  numberInRow[pivotRow] = 0;
  //store column in L, compress in U and take column out
  CoinBigIndex l = lengthL_;

  if ( l + 1 > lengthAreaL_ ) {
    //need more memory
    //if ((messageLevel_&4)!=0) 
    PIPS_APP_LOG_SEV(warning)  << "more memory needed in middle of invert";
    return false;
  }
  //l+=currentAreaL_->elementByColumn-elementL_;
  //CoinBigIndex lSave=l;
  CoinBigIndex * startColumnL = startColumnL_.array();
  CoinFactorizationDouble * elementL = elementL_.array();
  int * indexRowL = indexRowL_.array();
  startColumnL[numberGoodL_] = l;	//for luck and first time
  numberGoodL_++;
  startColumnL[numberGoodL_] = l + 1;
  lengthL_++;
  CoinFactorizationDouble pivotElement;
  CoinFactorizationDouble otherMultiplier;
  int otherRow;
  int * saveColumn = saveColumn_.array();
  CoinFactorizationDouble *elementU = elementU_.array();
  int * indexRowU = indexRowU_.array();

  if ( indexRowU[startColumn] == pivotRow ) {
    pivotElement = elementU[startColumn];
    otherMultiplier = elementU[startColumn + 1];
    otherRow = indexRowU[startColumn + 1];
  } else {
    pivotElement = elementU[startColumn + 1];
    otherMultiplier = elementU[startColumn];
    otherRow = indexRowU[startColumn];
  }
  int numberSave = numberInRow[otherRow];
  CoinFactorizationDouble pivotMultiplier = 1.0 / pivotElement;

  CoinFactorizationDouble * pivotRegion = pivotRegion_.array();
  pivotRegion[numberGoodU_] = pivotMultiplier;
  numberInColumn[pivotColumn] = 0;
  otherMultiplier = otherMultiplier * pivotMultiplier;
  indexRowL[l] = otherRow;
  elementL[l] = otherMultiplier;
  //take pivot column out of row list for otherRow
  CoinBigIndex start = startRowU[otherRow];
  CoinBigIndex end = start + numberSave;
  CoinBigIndex where = start;
  int * indexColumnU = indexColumnU_.array();

  while ( indexColumnU[where] != pivotColumn ) {
    where++;
  }				/* endwhile */
  assert ( where < end );
  end--;
  indexColumnU[where] = indexColumnU[end]; //this works b/c columns aren't ordered within rows
  int numberAdded = 0;
  int numberDeleted = 0;

  //pack down and move to work
  int j;
  const int * nextCount = nextCount_.array();
  int * nextColumn = nextColumn_.array();
	
  // forevery nonzero column in pivot row
  // e.g. every nonzero entry in pivot row, which can cause fill-in
  // but fill-in only occurs in otherRow
  for ( j = startRow; j < endRow; j++ ) {
    int iColumn = indexColumnU[j];
    if ( iColumn != pivotColumn ) {
      // need to update U[otherRow,iColumn]
      CoinBigIndex startColumn = startColumnU[iColumn];
      CoinBigIndex endColumn = startColumn + numberInColumn[iColumn];
      int iRow = indexRowU[startColumn];
      CoinFactorizationDouble value = elementU[startColumn];
      double largest;
      bool foundOther = false;

      //leave room for pivot
      //pivot element goes to beginning of iColumn
      //so push everything back by one element
      CoinBigIndex put = startColumn + 1;
      CoinBigIndex positionLargest = -1;
      CoinFactorizationDouble thisPivotValue = 0.0;
      CoinFactorizationDouble otherElement = 0.0;
      CoinFactorizationDouble nextValue = elementU[put];;
      int nextIRow = indexRowU[put];

      //compress column and find largest not updated
      if ( iRow != pivotRow ) {
	if ( iRow != otherRow ) {
	  largest = fabs ( value );
	  elementU[put] = value;
	  indexRowU[put] = iRow;
	  positionLargest = put;
	  put++;
	  CoinBigIndex i;
	  for ( i = startColumn + 1; i < endColumn; i++ ) {
	    iRow = nextIRow;
	    value = nextValue;
#ifdef ZEROFAULT
	    // doesn't matter reading uninitialized but annoys checking
	    if ( i + 1 < endColumn ) {
#endif
	      nextIRow = indexRowU[i + 1];
	      nextValue = elementU[i + 1];
#ifdef ZEROFAULT
	    }
#endif
	    if ( iRow != pivotRow ) {
	      if ( iRow != otherRow ) {
		//keep
		indexRowU[put] = iRow;
		elementU[put] = value;;
		put++;
	      } else {
		otherElement = value;
		foundOther = true;
	      }
	    } else {
	      thisPivotValue = value;
	    }
	  }
	} else {
	  otherElement = value;
	  foundOther = true;
	  //need to find largest
	  largest = 0.0;
	  CoinBigIndex i;
	  for ( i = startColumn + 1; i < endColumn; i++ ) {
	    iRow = nextIRow;
	    value = nextValue;
#ifdef ZEROFAULT
	    // doesn't matter reading uninitialized but annoys checking
	    if ( i + 1 < endColumn ) {
#endif
	      nextIRow = indexRowU[i + 1];
	      nextValue = elementU[i + 1];
#ifdef ZEROFAULT
	    }
#endif
	    if ( iRow != pivotRow ) {
	      //keep
	      indexRowU[put] = iRow;
	      elementU[put] = value;;
	      double absValue = fabs ( value );

	      if ( absValue > largest ) {
		largest = absValue;
		positionLargest = put;
	      }
	      put++;
	    } else {
	      thisPivotValue = value;
	    }
	  }
	}
      } else {
	//need to find largest
	largest = 0.0;
	thisPivotValue = value;
	CoinBigIndex i;
	for ( i = startColumn + 1; i < endColumn; i++ ) {
	  iRow = nextIRow;
	  value = nextValue;
#ifdef ZEROFAULT
	  // doesn't matter reading uninitialized but annoys checking
	  if ( i + 1 < endColumn ) {
#endif
	    nextIRow = indexRowU[i + 1];
	    nextValue = elementU[i + 1];
#ifdef ZEROFAULT
	  }
#endif
	  if ( iRow != otherRow ) {
	    //keep
	    indexRowU[put] = iRow;
	    elementU[put] = value;;
	    double absValue = fabs ( value );

	    if ( absValue > largest ) {
	      largest = absValue;
	      positionLargest = put;
	    }
	    put++;
	  } else {
	    otherElement = value;
	    foundOther = true;
	  }
	}
      }
      //slot in pivot
      elementU[startColumn] = thisPivotValue; // this is U[pivotRow,iColumn]
      indexRowU[startColumn] = pivotRow;
      //clean up counts
      startColumn++;
      numberInColumn[iColumn] = put - startColumn;
      numberInColumnPlus[iColumn]++;
      startColumnU[iColumn]++;
      otherElement = otherElement - thisPivotValue * otherMultiplier;
      double absValue = fabs ( otherElement );

      if ( absValue > zeroTolerance_ ) {
	if ( !foundOther ) {
	  //otherRow had zero here, may need more space for fill-in
	  saveColumn[numberAdded++] = iColumn; //remember which column we added element to
	  int next = nextColumn[iColumn];
	  CoinBigIndex space;

	  space = startColumnU[next] - put - numberInColumnPlus[next];
	  if ( space <= 0 ) {
	    //getColumnSpace also moves fixed part
	    int number = numberInColumn[iColumn];

	    if ( !getColumnSpace ( iColumn, number + 1 ) ) {
	      return false;
	    }
	    //redo starts
	    positionLargest =
	      positionLargest + startColumnU[iColumn] - startColumn;
	    startColumn = startColumnU[iColumn];
	    put = startColumn + number;
	  }
	}
	elementU[put] = otherElement;
	indexRowU[put] = otherRow;
	if ( absValue > largest ) {
	  largest = absValue;
	  positionLargest = put;
	}
	put++;
      } else {
	if ( foundOther ) {
	  // row operation turned nonzero into zero!
	  numberDeleted++;
	  //take out of row list
	  CoinBigIndex where = start;

	  while ( indexColumnU[where] != iColumn ) {
	    where++;
	  }			/* endwhile */
	  assert ( where < end );
	  end--;
	  indexColumnU[where] = indexColumnU[end];
	}
      }
      numberInColumn[iColumn] = put - startColumn;
      //move largest
      //largest always at beginning of column
      if ( positionLargest >= 0 ) {
	value = elementU[positionLargest];
	iRow = indexRowU[positionLargest];
	elementU[positionLargest] = elementU[startColumn];
	indexRowU[positionLargest] = indexRowU[startColumn];
	elementU[startColumn] = value;
	indexRowU[startColumn] = iRow;
      }
      //linked list for column
      if (iColumn < numberColumns1_ && nextCount[iColumn + numberRows_] != -2) {
	//modify linked list
	deleteLink ( iColumn + numberRows_ );
	addLink ( iColumn + numberRows_, numberInColumn[iColumn] );
      }
    }
  }
  //get space for row list
  next = nextRow[otherRow];
  CoinBigIndex space;

  space = startRowU[next] - end;
  totalElements_ += numberAdded - numberDeleted;
  int number = numberAdded + ( end - start );

  if ( space < numberAdded ) {
    numberInRow[otherRow] = end - start;
    if ( !getRowSpace ( otherRow, number ) ) {
      return false;
    }
    end = startRowU[otherRow] + end - start;
  }
  // do linked lists and update counts
  numberInRow[otherRow] = number;
  if ( number != numberSave ) {
    deleteLink ( otherRow );
    addLink ( otherRow, number );
  }
  for ( j = 0; j < numberAdded; j++ ) {
    indexColumnU[end++] = saveColumn[j];
  }
  //modify linked list for pivots
  deleteLink ( pivotRow );
  deleteLink ( pivotColumn + numberRows_ );
  return true;
}

