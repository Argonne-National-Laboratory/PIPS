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
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CoinUtilsConfig.h"

#include <cassert>
#include <cstdio>

#include "CoinBALPFactorization.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinFinite.hpp"
#include <iostream>

#include "PIPSLogging.hpp"

int CoinBALPFactorization::factorizeRect (
    const CoinPackedMatrix & M1,
    const CoinPackedMatrix & M2,
    int *M1columnIsBasic,
    int *M1rowIsBasic,
    int *M2columnIsBasic,
    int M2extracolumns,
    double areaFactor )
{
  // maybe for speed will be better to leave as many regions as possible
  gutsOfDestructor();
  gutsOfInitialize(2);
  // ? is this correct
  //if (biasLU_==2)
  //biasLU_=3;
  if (areaFactor)
    areaFactor_ = areaFactor;
  const int * row1 = M1.getIndices();
  const CoinBigIndex * columnStart1 = M1.getVectorStarts();
  const int * columnLength1 = M1.getVectorLengths(); 
  const double * element1 = M1.getElements();
  int numberRows1=M1.getNumRows();
  int numberColumns1=M1.getNumCols();

  const int * row2 = M2.getIndices();
  const CoinBigIndex * columnStart2 = M2.getVectorStarts();
  const int * columnLength2 = M2.getVectorLengths(); 
  const double * element2 = M2.getElements();
  int numberRows2=M2.getNumRows();
  int numberColumns2=M2.getNumCols() + M2extracolumns;

  assert(numberRows1 == numberRows2);

  int numberRows = numberRows1;

  int numberBasic1 = 0, numberBasic2 = 0, numberBasic;
  CoinBigIndex numberElements=0;

  // compute how much in basis

  int i;

  for (i=0;i<numberColumns1;i++) {
    if (M1columnIsBasic[i]>=0) {
      numberBasic1++;
      numberElements += columnLength1[i];
    }
  }
  for (i=0;i<numberRows;i++) {
    if (M1rowIsBasic[i]>=0) {
      numberBasic1++;
      numberElements += 1;
    }
  }
  //printf("numbercolumns2: %d\n", numberColumns2);
  for (i=0;i<numberColumns2;i++) {
    //printf("%d %d %d %d\n",i, numberColumns2, M2columnIsBasic[i],numberBasic2);
    if (M2columnIsBasic[i]>=0) {
      numberBasic2++;
      if (i < M2.getNumCols()) numberElements += columnLength2[i];
    }
  }
  if ( numberBasic1 > numberRows1 ) {
    return -2; // say too many in basis -- M1 should be tall
  }
  numberBasic = numberBasic1 + numberBasic2;
  //printf("nbasic1: %d nbasic2: %d\n",numberBasic1, numberBasic2);

  // this can be tuned per application, it depends highly on the amount of fill-in
  // you can decrease this until you start to see "Factorization did ... compressions"
  numberElements = 4 * numberBasic + 3 * numberElements + 60000/scenariosPerProc_;
  getAreas ( numberRows, numberBasic1, numberBasic2, numberElements,
      2 * numberElements );
  //fill
  //copy
  numberBasic=0;
  numberElements=0;
  int * indexColumnU = indexColumnU_.array();
  int * indexRowU = indexRowU_.array();
  CoinFactorizationDouble * elementU = elementU_.array();
  for (i=0;i<numberColumns1;i++) {
    if (M1columnIsBasic[i]>=0) {
      CoinBigIndex j;
      for (j=columnStart1[i];j<columnStart1[i]+columnLength1[i];j++) {
	indexRowU[numberElements]=row1[j];
	indexColumnU[numberElements]=numberBasic;
	elementU[numberElements++]=element1[j];
      }
      numberBasic++;
    }
  }
  for (i=0;i<numberRows;i++) {
     if (M1rowIsBasic[i]>=0) {
       indexRowU[numberElements] = i;
       indexColumnU[numberElements] = numberBasic++;
       elementU[numberElements++] = -1.0; // slack
     }
  }
  for (i=0;i<numberColumns2;i++) {
    if (M2columnIsBasic[i]>=0) {
      if (i < M2.getNumCols()) {
        CoinBigIndex j;
        for (j=columnStart2[i];j<columnStart2[i]+columnLength2[i];j++) {
	  indexRowU[numberElements]=row2[j];
	  indexColumnU[numberElements]=numberBasic;
	  elementU[numberElements++]=element2[j];
        }
      }
      numberBasic++;
    }
  }

  lengthU_ = numberElements;
  maximumU_ = numberElements;
  numberColumns1_ = numberBasic1;

	//printf("ncolsT: %d ncols1: %d ncols: %d nrows: %d\n",numberColumnsT_,numberColumns1_, numberColumns_, numberRows_);
  preProcessRect ( 0);
  factorRect ();
  numberBasic=0;
  //printf("done in factorize, status: %d\n", status_);
  if (status_ == 0) {
  	//printf("got it -- done\n");

  } else if (status_ == -1) {
    const int * pivotColumn = pivotColumn_.array();

    for (i=0;i<numberColumns1_;i++) {
      if (M1columnIsBasic[i]>=0) {
	if (pivotColumn[numberBasic]>=0) 
	  M1columnIsBasic[i]=pivotColumn[numberBasic];
	else 
	  M1columnIsBasic[i]=-1;
	numberBasic++;
      }
    }
  }
	//printf("returning\n");
  return status_;
}



//Does most of factorization
  int
CoinBALPFactorization::factorRect ()
{
  int * lastColumn = lastColumn_.array();
  int * lastRow = lastRow_.array();
  //sparse
  status_ = factorSparseRect ();
  //printf("status: %d\n",status_);
  switch ( status_ ) {
    case 0:			//finished
      totalElements_ = 0;
      {
	int * pivotColumn = pivotColumn_.array();
	//printf("goodU: %d ncols1: %d ncols: %d nrows: %d\n",numberGoodU_,numberColumns1_, numberColumns_, numberRows_);
	if ( numberGoodU_ < numberColumns1_ ) {
	  int i, k;
	  // Clean out unset nextRow
	  int * nextRow = nextRow_.array();
	  //int nSing =0;
	  k=nextRow[maximumRowsExtra_];
	  while (k!=maximumRowsExtra_) {
	    int iRow = k;
	    k=nextRow[k];
	    //nSing++;
	    nextRow[iRow]=-1;
	  }
	  //assert (nSing);
	  //printf("%d singularities - good %d rows %d\n",nSing,
	  //     numberGoodU_,numberRows_);
	  // Now nextRow has -1 or sequence into numberGoodU_;
	  int * permuteA = permute_.array();
	  for ( i = 0; i < numberRows_; i++ ) {
	    int iGood= nextRow[i];
	    if (iGood>=0)
	      permuteA[iGood]=i;
	  }

	  // swap arrays
	  permute_.swap(nextRow_);
	  int * permute = permute_.array();
	  for ( i = 0; i < numberRows_; i++ ) {
	    lastRow[i] = -1;
	  }
	  for ( i = 0; i < numberColumns_; i++ ) {
	    lastColumn[i] = -1;
	  }
	  for ( i = 0; i < numberGoodU_; i++ ) {
	    int goodRow = permuteA[i];	//valid pivot row
	    int goodColumn = pivotColumn[i];

	    lastRow[goodRow] = goodColumn;	//will now have -1 or column sequence
	    lastColumn[goodColumn] = goodRow;	//will now have -1 or row sequence
	  }
	  nextRow_.conditionalDelete();
	  k = 0;
	  //copy back and count
	  for ( i = 0; i < numberRows_; i++ ) {
	    permute[i] = lastRow[i];
	    if ( permute[i] < 0 ) {
	      //std::cout << i << " " <<permute[i] << std::endl;
	    } else {
	      k++;
	    }
	  }
	  for ( i = 0; i < numberColumns_; i++ ) {
	    pivotColumn[i] = lastColumn[i];
	  }
	  //if ((messageLevel_&4)!=0) 
	  PIPS_APP_LOG_SEV(warning) <<"Factorization has "<<numberRows_-k
				    <<" singularities";
	  status_ = -1;
	}
      }
      break;
      // dense
    case 2:
      status_=factorDense();
      if(!status_) 
	break;
    default:
      //singular ? or some error
      //if ((messageLevel_&4)!=0) 
      PIPS_APP_LOG_SEV(error) << "CoinBALPFactorization::factorRect() Error " << status_ ;
      break;
  }				/* endswitch */
  //clean up
  //printf("end switch, status: %d\n",status_);
  if ( !status_ ) {
    if ( numberCompressions_)
      PIPS_APP_LOG_SEV(info) << "        Factorization did "<<numberCompressions_
			     << " compressions";
    if ( numberCompressions_ > 10 ) {
      areaFactor_ *= 1.1;
    }
    numberCompressions_=0;
    //printf("calling cleanup\n");
    cleanupRect (  );
  }
  //printf("returning\n");
  return status_;
}

static int counter1=0;


  int
CoinBALPFactorization::factorSparseRect ()
{
  int larger;

  if ( numberRows_ < numberColumns_ ) {
    larger = numberColumns_;
  } else {
    larger = numberRows_;
  }
  int returnCode;
#define LARGELIMIT 65530
#define SMALL_SET 65531
#define SMALL_UNSET (SMALL_SET+1)
#define LARGE_SET COIN_INT_MAX-10
#define LARGE_UNSET (LARGE_SET+1)
  if ( larger < LARGELIMIT )
    returnCode = factorSparseSmallRect();
  else {
    returnCode = factorSparseLargeRect();
  }
  return returnCode;
}
//  factorSparse.  Does sparse phase of factorization
//return code is <0 error, 0= finished
  int
CoinBALPFactorization::factorSparseSmallRect ()
{
  int *indexRow = indexRowU_.array();
  int *indexColumn = indexColumnU_.array();
  CoinFactorizationDouble *element = elementU_.array();
  int count = 1;
  workArea_.conditionalNew(numberRows_);
  CoinFactorizationDouble * workArea = workArea_.array();
#ifndef NDEBUG
  counter1++;
#endif
  CoinZeroN ( workArea, numberRows_ );
  //get space for bit work area
  CoinBigIndex workSize = 1000;
  workArea2_.conditionalNew(workSize);
  unsigned int * workArea2 = workArea2_.array();

  //set markRow so no rows updated
  unsigned short * markRow = reinterpret_cast<unsigned short *> (markRow_.array());
  CoinFillN (  markRow, numberRows_, static_cast<unsigned short> (SMALL_UNSET));
  int status = 0;
  int * numberInRow = numberInRow_.array();
  int * numberInColumn = numberInColumn_.array();
  CoinBigIndex * startColumnU = startColumnU_.array();

  /*
  int * numberInColumnPlus = numberInColumnPlus_.array();
  //do slacks first
   CoinBigIndex * startColumnL = startColumnL_.array();
  // not sure why needs to be square
  // TODO: make this work for rectangular
	
	if (biasLU_<3&&numberColumns_==numberRows_) {
	//if (biasLU_<3) {
		printf("in slacks\n");
		int iPivotColumn;
		int * pivotColumn = pivotColumn_.array();
		int * nextRow = nextRow_.array();
		int * lastRow = lastRow_.array();
		for ( iPivotColumn = 0; iPivotColumn < numberColumns1_;
				iPivotColumn++ ) {
			if ( numberInColumn[iPivotColumn] == 1 ) {
				CoinBigIndex start = startColumnU[iPivotColumn];
				CoinFactorizationDouble value = element[start];
				if ( value == slackValue_ && numberInColumnPlus[iPivotColumn] == 0 ) {
					// treat as slack
					int iRow = indexRow[start];
					// but only if row not marked
					if (numberInRow[iRow]>0) {
						totalElements_ -= numberInRow[iRow];
						//take out this bit of indexColumnU
						int next = nextRow[iRow];
						int last = lastRow[iRow];

						nextRow[last] = next;
						lastRow[next] = last;
						nextRow[iRow] = numberGoodU_;	//use for permute
						lastRow[iRow] = -2; //mark
						//modify linked list for pivots
						deleteLink ( iRow );
						numberInRow[iRow]=-1;
						numberInColumn[iPivotColumn]=0;
						numberGoodL_++;
						startColumnL[numberGoodL_] = 0;
						pivotColumn[numberGoodU_] = iPivotColumn;
						numberGoodU_++;
					}
				}
			}
		}
		// redo
		preProcessRect(4);
		CoinFillN (  markRow, numberRows_, static_cast<unsigned short> (SMALL_UNSET));
	}
  */
  numberSlacks_ = numberGoodU_;
  int *nextCount = nextCount_.array();
  int *firstCount = firstCount_.array();
  CoinBigIndex *startRow = startRowU_.array();
  CoinBigIndex *startColumn = startColumnU;

  double pivotTolerance = pivotTolerance_;
  int numberTrials = numberTrials_;
  int numberRows = numberRows_;
  // Put column singletons first - (if false)
  separateLinks(1,(biasLU_>1));
#ifndef NDEBUG
  int counter2=0;
#endif
  while ( count <= biggerDimension_ ) {
    //printf("count = %d\n",count);
#ifndef NDEBUG
		counter2++;
		int badRow=-1;
		if (counter1==-1&&counter2>=0) {
			// check counts consistent
			for (int iCount=1;iCount<numberRows_;iCount++) {
				int look = firstCount[iCount];
				while ( look >= 0 ) {
					if ( look < numberRows_ ) {
						int iRow = look;
						if (iRow==badRow)
							printf("row count for row %d is %d\n",iCount,iRow);
						if ( numberInRow[iRow] != iCount ) {
							printf("failed debug on %d entry to factorSparse and %d try\n",
									counter1,counter2);
							printf("row %d - count %d number %d\n",iRow,iCount,numberInRow[iRow]);
							abort();
						}
						look = nextCount[look];
					} else {
						int iColumn = look - numberRows;
						if ( numberInColumn[iColumn] != iCount ) {
							printf("failed debug on %d entry to factorSparse and %d try\n",
									counter1,counter2);
							printf("column %d - count %d number %d\n",iColumn,iCount,numberInColumn[iColumn]);
							abort();
						}
						look = nextCount[look];
					}
				}
			}
		}
#endif
    CoinBigIndex minimumCount = COIN_INT_MAX;
		double minimumCost = COIN_DBL_MAX;

		int iPivotRow = -1;
		int iPivotColumn = -1;
		int pivotRowPosition = -1;
		int pivotColumnPosition = -1;
		int look = firstCount[count];
		int trials = 0;
		int * pivotColumn = pivotColumn_.array();

		if ( count == 1 && firstCount[1] >= 0 &&!biasLU_) {
			//do column singletons first to put more in U
			while ( look >= 0 ) {
				if ( look < numberRows_ ) {
					look = nextCount[look];
				} else {
					int iColumn = look - numberRows_;

					assert ( numberInColumn[iColumn] == count );
					CoinBigIndex start = startColumnU[iColumn];
					int iRow = indexRow[start];

					iPivotRow = iRow;
					pivotRowPosition = start;
					iPivotColumn = iColumn;
					assert (iPivotRow>=0&&iPivotColumn>=0);
					assert (iPivotColumn < numberColumns1_);
					pivotColumnPosition = -1;
					look = -1;
					break;
				}
			}			/* endwhile */
			if ( iPivotRow < 0 ) {
				//back to singletons
				look = firstCount[1];
			}
		}
		while ( look >= 0 ) {
			//printf("look = %d\n",look);
			if ( look < numberRows_ ) {
				int iRow = look;
#ifndef NDEBUG        
				if ( numberInRow[iRow] != count ) {
					printf("failed on %d entry to factorSparse and %d try\n",
							counter1,counter2);
					printf("row %d - count %d number %d\n",iRow,count,numberInRow[iRow]);
					abort();
				}
#endif
				look = nextCount[look];
				bool rejected = false;
				CoinBigIndex start = startRow[iRow];
				CoinBigIndex end = start + count;

				CoinBigIndex i;
				for ( i = start; i < end; i++ ) {
					int iColumn = indexColumn[i];
					assert (numberInColumn[iColumn]>0);
					if (iColumn >= numberColumns1_) continue;
					double cost = ( count - 1 ) * numberInColumn[iColumn];

					if ( cost < minimumCost ) {
						CoinBigIndex where = startColumn[iColumn];
						double minimumValue = element[where];

						minimumValue = fabs ( minimumValue ) * pivotTolerance;
						while ( indexRow[where] != iRow ) {
							where++;
						}			/* endwhile */
						assert ( where < startColumn[iColumn] +
								numberInColumn[iColumn]);
						CoinFactorizationDouble value = element[where];

						value = fabs ( value );
						if ( value >= minimumValue ) {
							minimumCost = cost;
							minimumCount = numberInColumn[iColumn];
							iPivotRow = iRow;
							pivotRowPosition = -1;
							iPivotColumn = iColumn;
							assert (iPivotRow>=0&&iPivotColumn>=0);
							pivotColumnPosition = i;
							rejected=false;
							if ( minimumCount < count ) {
								look = -1;
								break;
							}
						} else if ( iPivotRow == -1 ) {
							rejected = true;
						}
					}
				}
				trials++;
				if ( trials >= numberTrials && iPivotRow >= 0 ) {
					look = -1;
					break;
				}
				if ( rejected ) {
					//take out for moment
					//eligible when row changes
					deleteLink ( iRow );
					addLink ( iRow, biggerDimension_ + 1 );
				}
			} else {
				int iColumn = look - numberRows;

				assert ( numberInColumn[iColumn] == count );
				look = nextCount[look];
				CoinBigIndex start = startColumn[iColumn];
				CoinBigIndex end = start + numberInColumn[iColumn];
				CoinFactorizationDouble minimumValue = element[start];

				minimumValue = fabs ( minimumValue ) * pivotTolerance;
				CoinBigIndex i;
				for ( i = start; i < end; i++ ) {
					CoinFactorizationDouble value = element[i];

					value = fabs ( value );
					if ( value >= minimumValue ) {
						int iRow = indexRow[i];
						int nInRow = numberInRow[iRow];
						assert (nInRow>0);
						double cost = ( count - 1 ) * nInRow;

						if ( cost < minimumCost ) {
							minimumCost = cost;
							minimumCount = nInRow;
							iPivotRow = iRow;
							pivotRowPosition = i;
							iPivotColumn = iColumn;
							assert(iPivotColumn < numberColumns1_);
							assert (iPivotRow>=0&&iPivotColumn>=0);
							pivotColumnPosition = -1;
							if ( minimumCount <= count + 1 ) {
								look = -1;
								break;
							}
						}
					}
				}
				trials++;
				if ( trials >= numberTrials && iPivotRow >= 0 ) {
					look = -1;
					break;
				}
			}
		} /* endwhile */
    if (iPivotRow>=0) {
      assert (iPivotRow<numberRows_);
      assert (iPivotColumn < numberColumns1_);
      int numberDoRow = numberInRow[iPivotRow] - 1;
      int numberDoColumn = numberInColumn[iPivotColumn] - 1;
    //printf("prow: %d pcol: %d doRow: %d doCol %d\n",iPivotRow, iPivotColumn,numberDoRow,numberDoColumn);

      totalElements_ -= ( numberDoRow + numberDoColumn + 1 );
      if ( numberDoColumn > 0 ) {
	if ( numberDoRow > 0 ) {
	  if ( numberDoColumn > 1 ) {
	    //  if (1) {
	    //need to adjust more for cache and SMP
	    //allow at least 4 extra
	    int increment = numberDoColumn + 1 + 4;

	    if ( increment & 15 ) {
	      increment = increment & ( ~15 );
	      increment += 16;
	    }
	    int increment2 =

	      ( increment + COINFACTORIZATION_BITS_PER_INT - 1 ) >> COINFACTORIZATION_SHIFT_PER_INT;
	    CoinBigIndex size = increment2 * numberDoRow;

	    if ( size > workSize ) {
	      workSize = size;
	      workArea2_.conditionalNew(workSize);
	      workArea2 = workArea2_.array();
	    }
	    bool goodPivot;
	    //branch out to best pivot routine
	    //printf("call pivot\n");
	    goodPivot = pivot ( iPivotRow, iPivotColumn,
		pivotRowPosition, pivotColumnPosition,
		workArea, workArea2, 
		increment2,  markRow ,
		SMALL_SET);

	    if ( !goodPivot ) {
	      status = -99;
	      count=biggerDimension_+1;
	      break;
	    }
	  } else {
	  	//printf("call pivotoneotherrow\n");
	    if ( !pivotOneOtherRow ( iPivotRow, iPivotColumn ) ) {
	      status = -99;
	      count=biggerDimension_+1;
	      break;
	    }
	  }
	  } else {
	  	//printf("call pivotrowsingleton\n");
	    assert (!numberDoRow);
	    if ( !pivotRowSingleton ( iPivotRow, iPivotColumn ) ) {
	      status = -99;
	      count=biggerDimension_+1;
	      break;
	    }
	  }
	} else {
	  	//printf("call pivotcolumnsingleton\n");
	  assert (!numberDoColumn);
	  if ( !pivotColumnSingleton ( iPivotRow, iPivotColumn ) ) {
	    status = -99;
	    count=biggerDimension_+1;
	    break;
	  }
	}
	assert (nextRow_.array()[iPivotRow]==numberGoodU_);
	rowPivotFlag_.array()[iPivotRow] = 1;
	pivotColumn[numberGoodU_] = iPivotColumn;
	numberGoodU_++;
	// This should not need to be trapped here - but be safe
	if (numberGoodU_==numberColumns1_) 
	  break;
#if COIN_DEBUG==2
	checkConsistency (  );
#endif


	// start at 1 again
	count = 1;
      } else {
	//end of this - onto next
	count++;
      } 
    }				/* endwhile */
    workArea_.conditionalDelete() ;
    workArea2_.conditionalDelete() ;
    return status;
  }


// seems like the main difference here is with markRow as short vs. int
int
CoinBALPFactorization::factorSparseLargeRect (  )
{
  int *indexRow = indexRowU_.array();
  int *indexColumn = indexColumnU_.array();
  CoinFactorizationDouble *element = elementU_.array();
  int count = 1;
  workArea_.conditionalNew(numberRows_);
  CoinFactorizationDouble * workArea = workArea_.array();
#ifndef NDEBUG
  counter1++;
#endif
  CoinZeroN ( workArea, numberRows_ );
  //get space for bit work area
  CoinBigIndex workSize = 1000;
  workArea2_.conditionalNew(workSize);
  unsigned int * workArea2 = workArea2_.array();

  //set markRow so no rows updated
  int * markRow = markRow_.array();
  CoinFillN ( markRow, numberRows_, COIN_INT_MAX-10+1);
  int status = 0;
  int * numberInRow = numberInRow_.array();
  int * numberInColumn = numberInColumn_.array();
  CoinBigIndex * startColumnU = startColumnU_.array();
  /*
  //do slacks first
  int * numberInColumnPlus = numberInColumnPlus_.array();
  CoinBigIndex * startColumnL = startColumnL_.array();
	if (biasLU_<3&&numberColumns_==numberRows_) {
		int iPivotColumn;
		int * pivotColumn = pivotColumn_.array();
		int * nextRow = nextRow_.array();
		int * lastRow = lastRow_.array();
		for ( iPivotColumn = 0; iPivotColumn < numberColumns1_;
				iPivotColumn++ ) {
			if ( numberInColumn[iPivotColumn] == 1 ) {
				CoinBigIndex start = startColumnU[iPivotColumn];
				CoinFactorizationDouble value = element[start];
				if ( value == slackValue_ && numberInColumnPlus[iPivotColumn] == 0 ) {
					// treat as slack
					int iRow = indexRow[start];
					// but only if row not marked
					if (numberInRow[iRow]>0) {
						totalElements_ -= numberInRow[iRow];
						//take out this bit of indexColumnU
						int next = nextRow[iRow];
						int last = lastRow[iRow];

						nextRow[last] = next;
						lastRow[next] = last;
						nextRow[iRow] = numberGoodU_;	//use for permute
						lastRow[iRow] = -2; //mark
						//modify linked list for pivots
						deleteLink ( iRow );
						numberInRow[iRow]=-1;
						numberInColumn[iPivotColumn]=0;
						numberGoodL_++;
						startColumnL[numberGoodL_] = 0;
						pivotColumn[numberGoodU_] = iPivotColumn;
						numberGoodU_++;
					}
				}
      }
    }
    // redo
    preProcessRect(4);
    CoinFillN ( markRow, numberRows_, COIN_INT_MAX-10+1);
  }*/
  numberSlacks_ = numberGoodU_;
  int *nextCount = nextCount_.array();
  int *firstCount = firstCount_.array();
  CoinBigIndex *startRow = startRowU_.array();
  CoinBigIndex *startColumn = startColumnU;
  //double *elementL = elementL_.array();
  //int *indexRowL = indexRowL_.array();
  //int *saveColumn = saveColumn_.array();
  //int *nextRow = nextRow_.array();
  //int *lastRow = lastRow_.array() ;
  double pivotTolerance = pivotTolerance_;
  int numberTrials = numberTrials_;
  int numberRows = numberRows_;
  // Put column singletons first - (if false)
  separateLinks(1,(biasLU_>1));
#ifndef NDEBUG
  int counter2=0;
#endif
  while ( count <= biggerDimension_ ) {
#ifndef NDEBUG
    counter2++;
    int badRow=-1;
    if (counter1==-1&&counter2>=0) {
      // check counts consistent
      for (int iCount=1;iCount<numberRows_;iCount++) {
        int look = firstCount[iCount];
        while ( look >= 0 ) {
          if ( look < numberRows_ ) {
            int iRow = look;
            if (iRow==badRow)
              printf("row count for row %d is %d\n",iCount,iRow);
            if ( numberInRow[iRow] != iCount ) {
              printf("failed debug on %d entry to factorSparse and %d try\n",
                     counter1,counter2);
              printf("row %d - count %d number %d\n",iRow,iCount,numberInRow[iRow]);
              abort();
            }
            look = nextCount[look];
          } else {
            int iColumn = look - numberRows;
            if ( numberInColumn[iColumn] != iCount ) {
              printf("failed debug on %d entry to factorSparse and %d try\n",
                     counter1,counter2);
              printf("column %d - count %d number %d\n",iColumn,iCount,numberInColumn[iColumn]);
              abort();
            }
            look = nextCount[look];
          }
        }
      }
    }
#endif
    CoinBigIndex minimumCount = COIN_INT_MAX;
    double minimumCost = COIN_DBL_MAX;

    int iPivotRow = -1;
    int iPivotColumn = -1;
    int pivotRowPosition = -1;
    int pivotColumnPosition = -1;
    int look = firstCount[count];
    int trials = 0;
    int * pivotColumn = pivotColumn_.array();

    if ( count == 1 && firstCount[1] >= 0 &&!biasLU_) {
      //do column singletons first to put more in U
      while ( look >= 0 ) {
        if ( look < numberRows_ ) {
          look = nextCount[look];
        } else {
          int iColumn = look - numberRows_;
          
          assert ( numberInColumn[iColumn] == count );
          CoinBigIndex start = startColumnU[iColumn];
          int iRow = indexRow[start];
          
          iPivotRow = iRow;
          pivotRowPosition = start;
          iPivotColumn = iColumn;
          assert (iPivotRow>=0&&iPivotColumn>=0);
	  assert (iPivotColumn < numberColumns1_);
          pivotColumnPosition = -1;
          look = -1;
          break;
        }
      }			/* endwhile */
      if ( iPivotRow < 0 ) {
        //back to singletons
        look = firstCount[1];
      }
    }
    while ( look >= 0 ) {
      if ( look < numberRows_ ) {
        int iRow = look;
#ifndef NDEBUG        
        if ( numberInRow[iRow] != count ) {
          printf("failed on %d entry to factorSparse and %d try\n",
                 counter1,counter2);
          printf("row %d - count %d number %d\n",iRow,count,numberInRow[iRow]);
          abort();
        }
#endif
        look = nextCount[look];
        bool rejected = false;
        CoinBigIndex start = startRow[iRow];
        CoinBigIndex end = start + count;
        
        CoinBigIndex i;
        for ( i = start; i < end; i++ ) {
          int iColumn = indexColumn[i];
          assert (numberInColumn[iColumn]>0);
	  if (iColumn >= numberColumns1_) continue;
          double cost = ( count - 1 ) * numberInColumn[iColumn];
          
          if ( cost < minimumCost ) {
            CoinBigIndex where = startColumn[iColumn];
            CoinFactorizationDouble minimumValue = element[where];
            
            minimumValue = fabs ( minimumValue ) * pivotTolerance;
            while ( indexRow[where] != iRow ) {
              where++;
            }			/* endwhile */
            assert ( where < startColumn[iColumn] +
                     numberInColumn[iColumn]);
            CoinFactorizationDouble value = element[where];
            
            value = fabs ( value );
            if ( value >= minimumValue ) {
              minimumCost = cost;
              minimumCount = numberInColumn[iColumn];
              iPivotRow = iRow;
              pivotRowPosition = -1;
              iPivotColumn = iColumn;
              assert (iPivotRow>=0&&iPivotColumn>=0);
              pivotColumnPosition = i;
              rejected=false;
              if ( minimumCount < count ) {
                look = -1;
                break;
              }
            } else if ( iPivotRow == -1 ) {
              rejected = true;
            }
          }
        }
        trials++;
        if ( trials >= numberTrials && iPivotRow >= 0 ) {
          look = -1;
          break;
        }
        if ( rejected ) {
          //take out for moment
          //eligible when row changes
          deleteLink ( iRow );
          addLink ( iRow, biggerDimension_ + 1 );
        }
      } else {
        int iColumn = look - numberRows;
        
        assert ( numberInColumn[iColumn] == count );
        look = nextCount[look];
        CoinBigIndex start = startColumn[iColumn];
        CoinBigIndex end = start + numberInColumn[iColumn];
        CoinFactorizationDouble minimumValue = element[start];
        
        minimumValue = fabs ( minimumValue ) * pivotTolerance;
        CoinBigIndex i;
        for ( i = start; i < end; i++ ) {
          CoinFactorizationDouble value = element[i];
          
          value = fabs ( value );
          if ( value >= minimumValue ) {
            int iRow = indexRow[i];
            int nInRow = numberInRow[iRow];
            assert (nInRow>0);
            double cost = ( count - 1 ) * nInRow;
            
            if ( cost < minimumCost ) {
              minimumCost = cost;
              minimumCount = nInRow;
              iPivotRow = iRow;
              pivotRowPosition = i;
              iPivotColumn = iColumn;
	      assert(iPivotColumn < numberColumns1_);
              assert (iPivotRow>=0&&iPivotColumn>=0);
              pivotColumnPosition = -1;
              if ( minimumCount <= count + 1 ) {
                look = -1;
                break;
              }
            }
          }
        }
        trials++;
        if ( trials >= numberTrials && iPivotRow >= 0 ) {
          look = -1;
          break;
        }
      }
    }				/* endwhile */
    if (iPivotRow>=0) {
      if ( iPivotRow >= 0 ) {
	assert (iPivotRow<numberRows_);
        assert (iPivotColumn < numberColumns1_);
        int numberDoRow = numberInRow[iPivotRow] - 1;
        int numberDoColumn = numberInColumn[iPivotColumn] - 1;
        
        totalElements_ -= ( numberDoRow + numberDoColumn + 1 );
        if ( numberDoColumn > 0 ) {
          if ( numberDoRow > 0 ) {
            if ( numberDoColumn > 1 ) {
              //  if (1) {
              //need to adjust more for cache and SMP
              //allow at least 4 extra
              int increment = numberDoColumn + 1 + 4;
              
              if ( increment & 15 ) {
                increment = increment & ( ~15 );
                increment += 16;
              }
              int increment2 =
                
                ( increment + COINFACTORIZATION_BITS_PER_INT - 1 ) >> COINFACTORIZATION_SHIFT_PER_INT;
              CoinBigIndex size = increment2 * numberDoRow;
              
              if ( size > workSize ) {
                workSize = size;
		workArea2_.conditionalNew(workSize);
		workArea2 = workArea2_.array();
              }
              bool goodPivot;
              
	      //might be able to do better by permuting
	      //branch out to best pivot routine 
	      goodPivot = pivot ( iPivotRow, iPivotColumn,
				  pivotRowPosition, pivotColumnPosition,
				  workArea, workArea2, 
				  increment2, markRow ,
				  LARGE_SET);

              if ( !goodPivot ) {
                status = -99;
                count=biggerDimension_+1;
                break;
              }
            } else {
              if ( !pivotOneOtherRow ( iPivotRow, iPivotColumn ) ) {
                status = -99;
                count=biggerDimension_+1;
                break;
              }
            }
          } else {
            assert (!numberDoRow);
            if ( !pivotRowSingleton ( iPivotRow, iPivotColumn ) ) {
              status = -99;
              count=biggerDimension_+1;
              break;
            }
          }
        } else {
          assert (!numberDoColumn);
          if ( !pivotColumnSingleton ( iPivotRow, iPivotColumn ) ) {
            status = -99;
            count=biggerDimension_+1;
            break;
          }
        }
	assert (nextRow_.array()[iPivotRow]==numberGoodU_);
	rowPivotFlag_.array()[iPivotRow] = 1;
        pivotColumn[numberGoodU_] = iPivotColumn;
        numberGoodU_++;
        // This should not need to be trapped here - but be safe
        if (numberGoodU_==numberColumns1_) 
          break;
      }
#if COIN_DEBUG==2
      checkConsistency (  );
#endif


      // start at 1 again
      count = 1;
    } else {
      //end of this - onto next
      count++;
    } 
  }				/* endwhile */
  workArea_.conditionalDelete() ;
  workArea2_.conditionalDelete() ;
  return status;
}

  //  show_self.  Debug show object
  void
    CoinBALPFactorization::show_self () const
    {
      int i;

      const int * pivotColumn = pivotColumn_.array();
      std::cout.precision(14);
      for ( i = 0; i < numberColumns1_; i++ ) {
	std::cout << "c " << i << " " << pivotColumn[i];
	if (pivotColumnBack()) std::cout<< " " << pivotColumnBack()[i];
	std::cout<< " " << pivotRegion_.array()[i];
	std::cout << std::endl;
      }
      for ( i = 0; i < numberRows_; i++) {
	std::cout << "r " << i << " " << permute_.array()[i] << " " << permuteBack_.array()[i] << std::endl;
      }
      for ( i = 0; i < numberColumns_; i++ ) {
	std::cout << "u " << i << " " << numberInColumn_.array()[i] << std::endl;
	int j;
	CoinSort_2(indexRowU_.array()+startColumnU_.array()[i],
	    indexRowU_.array()+startColumnU_.array()[i]+numberInColumn_.array()[i],
	    elementU_.array()+startColumnU_.array()[i]);
	for ( j = startColumnU_.array()[i]; j < startColumnU_.array()[i] + numberInColumn_.array()[i];
	    j++ ) {
	  assert (indexRowU_.array()[j]>=0&&indexRowU_.array()[j]<numberRows_);
	  assert (elementU_.array()[j]>-1.0e50&&elementU_.array()[j]<1.0e50);
	  std::cout << indexRowU_.array()[j] << " " << elementU_.array()[j] << std::endl;
	}
      }
      for ( i = 0; i < numberColumns1_; i++ ) {
	std::cout << "l " << i << " " << startColumnL_.array()[i + 1] -
	  startColumnL_.array()[i] << std::endl;
	CoinSort_2(indexRowL_.array()+startColumnL_.array()[i],
	    indexRowL_.array()+startColumnL_.array()[i+1],
	    elementL_.array()+startColumnL_.array()[i]);
	int j;
	for ( j = startColumnL_.array()[i]; j < startColumnL_.array()[i + 1]; j++ ) {
	  std::cout << indexRowL_.array()[j] << " " << elementL_.array()[j] << std::endl;
	}
      }

      for ( i = 0; i < numberColumnsT_; i++ ) {
	std::cout << "t " << i << " " << startColumnT_.array()[i + 1] -
	  startColumnT_.array()[i] << std::endl;
	CoinSort_2(indexRowT_.array()+startColumnT_.array()[i],
	    indexRowT_.array()+startColumnT_.array()[i+1],
	    elementT_.array()+startColumnT_.array()[i]);
	int j;
	for ( j = startColumnT_.array()[i]; j < startColumnT_.array()[i + 1]; j++ ) {
	  std::cout << indexRowT_.array()[j] << " " << elementT_.array()[j] << std::endl;
	}
      }
    }

//  cleanup.  End of factorization
void
CoinBALPFactorization::cleanupRect (  )
{

  getColumnSpace ( 0, COIN_INT_MAX >> 1 );	//compress
  // swap arrays
  //  numberInColumn_.swap(numberInColumnPlus_); // why??
  CoinBigIndex * startColumnU = startColumnU_.array();
  CoinBigIndex lastU = startColumnU[maximumColumnsExtra_];

  //free some memory here
  saveColumn_.conditionalDelete();
  markRow_.conditionalDelete() ;
  //firstCount_.conditionalDelete() ;
  nextCount_.conditionalDelete() ;
  lastCount_.conditionalDelete() ;
  int * numberInRow = numberInRow_.array();
  int * numberInColumn = numberInColumn_.array();
  int * numberInColumnPlus = numberInColumnPlus_.array();
  //make column starts OK
  //for best cache behavior get in order (last pivot at bottom of space)
  //that will need thinking about
  //use nextRow for permutation  (as that is what it is)
  int i;


  // swap arrays
  permute_.swap(nextRow_);
  //safety feature
  int * permute = permute_.array();
  int * rowPivotFlag = rowPivotFlag_.array();
  permute[numberRows_] = 0;
  permuteBack_.conditionalNew( maximumRowsExtra_ + 1);
  // handle rows that weren't pivoted
  // -- put them at the end
  int nonPivotCount = numberGoodU_;
  for ( i = 0; i < numberRows_; i++) {
	if (rowPivotFlag[i] == 0) { // didn't pivot
		permute[i] = nonPivotCount++;
	}
  }
  assert(nonPivotCount == numberRows_);

  int * permuteBack = permuteBack_.array();
#ifdef ZEROFAULT
  memset(permuteBack_.array(),'w',(maximumRowsExtra_+1)*sizeof(int));
#endif
  for ( i = 0; i < numberRows_; i++ ) {
    int iRow = permute[i];

    permuteBack[iRow] = i;
  }
  //redo nextRow_

#ifndef NDEBUG
  for ( i = 0; i < numberRows_; i++ ) {
    assert (permute[i]>=0&&permute[i]<numberRows_);
    assert (permuteBack[i]>=0&&permuteBack[i]<numberRows_);
  }
#endif

    // Redo total elements
    totalElements_=0;
    for ( i = 0; i < numberColumns_; i++ ) {
      int number1 = numberInColumnPlus[i];
      int number2 = numberInColumn[i];
      totalElements_ += number1;
      totalElements_ += number2;
      startColumnU[i] -= number1;
      numberInColumnPlus[i] = 0;
      numberInColumn[i] += number1;
    }

  int numberU = 0;
  pivotColumnBack_.conditionalNew( maximumRowsExtra_ + 1);
#ifdef ZEROFAULT
  memset(pivotColumnBack(),'q',(maximumRowsExtra_+1)*sizeof(int));
#endif
  const int * pivotColumn = pivotColumn_.array();
  int * pivotColumnB = pivotColumnBack();
  int *indexColumnU = indexColumnU_.array();
  CoinBigIndex *startColumn = startColumnU;
  int *indexRowU = indexRowU_.array();
  CoinFactorizationDouble *elementU = elementU_.array();

	assert(numberGoodU_ == numberColumns1_);
    for ( i = 0; i < numberColumns1_; i++ ) {
      int iColumn = pivotColumn[i]; // ith column pivoted
      
      pivotColumnB[iColumn] = i; // iColumn is in ith position
      if ( iColumn >= 0 ) {
	//wanted
	if ( numberU != iColumn ) {
	  numberInColumnPlus[iColumn] = numberU;
	} else {
	  numberInColumnPlus[iColumn] = -1;	//already in correct place
	}
	numberU++;
      }
    }
    // looks like this loop just relabels the columns
    // according to their pivot order
    for ( i = 0; i < numberColumns1_; i++ ) {
      int number = numberInColumn[i];	//always 0? -- no
      int where = numberInColumnPlus[i];
      
      numberInColumnPlus[i] = -1;
      CoinBigIndex start = startColumnU[i];
      
      while ( where >= 0 ) {
	//put where it should be
	int numberNext = numberInColumn[where];	//always 0? -- no
	int whereNext = numberInColumnPlus[where];
	CoinBigIndex startNext = startColumnU[where];
	
	numberInColumn[where] = number;
	numberInColumnPlus[where] = -1;
	startColumnU[where] = start;
	number = numberNext;
	where = whereNext;
	start = startNext;
      }				/* endwhile */
    }
    //sort - using indexColumn
    //reorder factor so columns are in order in memory
    //indexColumnU used as map from old indices to new indices
    // also serves to pack down
    CoinFillN ( indexColumnU_.array(), lastU, -1 );
    CoinBigIndex k = 0;
    
    for ( i = numberSlacks_; i < numberColumns_; i++ ) {
      CoinBigIndex start = startColumn[i];
      CoinBigIndex end = start + numberInColumn[i];
      
      CoinBigIndex j;
      for ( j = start; j < end; j++ ) {
	indexColumnU[j] = k++;
      }
    }
    for ( i = numberSlacks_; i < numberColumns_; i++ ) {
      CoinBigIndex start = startColumn[i];
      CoinBigIndex end = start + numberInColumn[i];
      
      CoinBigIndex j;
      for ( j = start; j < end; j++ ) {
	CoinBigIndex k = indexColumnU[j];
	int iRow = indexRowU[j];
	assert(iRow >= 0 && iRow < numberRows_);
	CoinFactorizationDouble element = elementU[j];
	
	while ( k != -1 ) {
	  CoinBigIndex kNext = indexColumnU[k];
	  int iRowNext = indexRowU[k];
	  assert(iRowNext >= 0 && iRowNext < numberRows_);
	  CoinFactorizationDouble elementNext = elementU[k];
	  
	  indexColumnU[k] = -1;
	  indexRowU[k] = iRow;
	  elementU[k] = element;
	  k = kNext;
	  iRow = iRowNext;
	  element = elementNext;
	}				/* endwhile */
      }
    }
    CoinZeroN ( startColumnU, numberSlacks_ );
    k = 0;
    for ( i = numberSlacks_; i < numberColumns1_; i++ ) {
      startColumnU[i] = k;
      k += numberInColumn[i];
    }
    maximumU_=k; // where is this used?

	// extract T
	CoinBigIndex *startColumnT = startColumnT_.array();
	CoinBigIndex *startRowT = startRowT_.array();
	int *numberInColumnT = numberInColumnT_.array();
	int *numberInRowT = numberInRowT_.array();

	int offsetT = k;
	k = 0;
	for (i = numberColumns1_; i < numberColumns_; i++) {
		int iT = i - numberColumns1_;
		startColumnT[iT] = k;
		k += (numberInColumnT[iT] = numberInColumn[i]);
	}
	// k contains total number of elements in T.
	
	startColumnT[numberColumnsT_] = k;
	// no space for updating factor
	lengthAreaT_ = k;
	//if (numberColumnsT_) printf("Density of T: %f\n",static_cast<double>(k)/(static_cast<double>(numberColumnsT_)*numberRows_));

	elementT_.conditionalNew( lengthAreaT_ );
	indexRowT_.conditionalNew( lengthAreaT_ );
	indexColumnT_.conditionalNew( lengthAreaT_ );
	convertRowToColumnT_.conditionalNew(lengthAreaT_);
	
	CoinFactorizationDouble *elementT = elementT_.array();
	int *indexRowT = indexRowT_.array();
	int *indexColumnT = indexColumnT_.array();
	
	memcpy(elementT,elementU+offsetT,k*sizeof(CoinFactorizationDouble));
	memcpy(indexRowT,indexRowU+offsetT,k*sizeof(int));

	assert(numberColumns_ == numberColumns1_ + numberColumnsT_);
	numberColumns_ = numberColumns1_;




	//if ( (messageLevel_ & 8)) {
	//PIPS_APP_LOG_SEV(info)
	std::stringstream ss; ss <<"        length of U "<<totalElements_<<", length of L "<<lengthL_;
	if (numberDense_)  ss<<" plus "<<numberDense_*numberDense_<<" from "<<numberDense_<<" dense rows";
	PIPS_ALG_LOG_SEV(debug) << ss.str();
	//std::cout<<std::endl;
	//}
  // and add L and dense
  totalElements_ += numberDense_*numberDense_+lengthL_;
  int * nextColumn = nextColumn_.array();
  int * lastColumn = lastColumn_.array();
  // See whether to have extra copy of R
  if (maximumU_>10*numberRows_||numberRows_<200) {
    // NO
    //printf("no extra copy of R\n");
    numberInColumnPlus_.conditionalDelete() ;
  } else {
  	//printf("extra copy of R\n");
    for ( i = 0; i < numberColumns_; i++ ) {
      lastColumn[i] = i - 1;
      nextColumn[i] = i + 1;
      numberInColumnPlus[i]=0;
    }
    nextColumn[numberColumns_ - 1] = maximumColumnsExtra_;
    lastColumn[maximumColumnsExtra_] = numberColumns_ - 1;
    nextColumn[maximumColumnsExtra_] = 0;
    lastColumn[0] = maximumColumnsExtra_;
  }
  numberU_ = numberU;
  numberGoodU_ = numberU;
  numberL_ = numberGoodL_;
#if COIN_DEBUG
  for ( i = 0; i < numberRows_; i++ ) {
    if ( permute[i] < 0 ) {
      std::cout << i << std::endl;
      abort (  );
    }
  }
#endif
  // relabel rows
  CoinZeroN(numberInRow,numberRows_);
  for ( i = numberSlacks_; i < numberColumns_; i++ ) {
    CoinBigIndex start = startColumnU[i];
    CoinBigIndex end = start + numberInColumn[i];

    totalElements_ += numberInColumn[i]; // didn't we set this already? may be a bug
    CoinBigIndex j;

    for ( j = start; j < end; j++ ) {
      int iRow = indexRowU[j];
      iRow = permute[iRow];
      indexRowU[j] = iRow;
      numberInRow[iRow]++; 
    }
  }

  CoinZeroN(numberInRowT,numberRows_);
  for ( i = 0; i < numberColumnsT_; i++ ) {
    CoinBigIndex start = startColumnT[i];
    CoinBigIndex end = start + numberInColumnT[i];

    totalElements_ += numberInColumnT[i];
    CoinBigIndex j;

    for ( j = start; j < end; j++ ) {
      int iRow = indexRowT[j];
      iRow = permute[iRow];
      indexRowT[j] = iRow;
      numberInRowT[iRow]++; 
    }
  }

  
  CoinFactorizationDouble * pivotRegion = pivotRegion_.array();
    //space for cross reference
    convertRowToColumnU_.conditionalNew( lengthAreaU_ );
    CoinBigIndex *convertRowToColumn = convertRowToColumnU_.array();
    CoinBigIndex j = 0;
    CoinBigIndex *startRow = startRowU_.array();
    
    int iRow;
    for ( iRow = 0; iRow < numberRows_; iRow++ ) {
      startRow[iRow] = j;
      j += numberInRow[iRow];
    }
    CoinBigIndex numberInU = j;
    
    CoinZeroN ( numberInRow_.array(), numberRows_ );
    
    for ( i = numberSlacks_; i < numberColumns_; i++ ) {
      CoinBigIndex start = startColumnU[i];
      CoinBigIndex end = start + numberInColumn[i];
      
      //CoinFactorizationDouble pivotValue = pivotRegion[i];

      CoinBigIndex j;
      for ( j = start; j < end; j++ ) {
	int iRow = indexRowU[j];
	int iLook = numberInRow[iRow];
	
	numberInRow[iRow] = iLook + 1;
	CoinBigIndex k = startRow[iRow] + iLook;
	
	indexColumnU[k] = i;
	convertRowToColumn[k] = j;
	//multiply by pivot
        elementU[j] *= pivotRegion[i];
      }
    }

    // row copy of T
    CoinBigIndex *convertRowToColumnT = convertRowToColumnT_.array();
    j = 0;
    
    for ( iRow = 0; iRow < numberRows_; iRow++ ) {
      startRowT[iRow] = j;
      j += numberInRowT[iRow];
    }
    startRowT[numberRows_] = j;
    
    //CoinBigIndex numberInT = j;
    
    CoinZeroN ( numberInRowT, numberRows_ );
    
    for ( i = 0; i < numberColumnsT_; i++ ) {
      CoinBigIndex start = startColumnT[i];
      CoinBigIndex end = start + numberInColumnT[i];
      
      CoinBigIndex j;
      for ( j = start; j < end; j++ ) {
	int iRow = indexRowT[j];
	int iLook = numberInRowT[iRow];
	
	numberInRowT[iRow] = iLook + 1;
	CoinBigIndex k = startRowT[iRow] + iLook;
	
	indexColumnT[k] = i;
	convertRowToColumnT[k] = j;
      }
    }
    int * nextRow = nextRow_.array();
    int * lastRow = lastRow_.array();
    for ( j = 0; j < numberRows_; j++ ) {
      lastRow[j] = j - 1;
      nextRow[j] = j + 1;
    }
    nextRow[numberRows_ - 1] = maximumRowsExtra_;
    lastRow[maximumRowsExtra_] = numberRows_ - 1;
    nextRow[maximumRowsExtra_] = 0;
    lastRow[0] = maximumRowsExtra_;
    startRow[maximumRowsExtra_] = numberInU;


  int firstReal = numberRows_;

  CoinBigIndex * startColumnL = startColumnL_.array();
  int * indexRowL = indexRowL_.array();
  // apply row permutation to L
  for ( i = numberGoodU_ - 1; i >= 0; i-- ) {
    CoinBigIndex start = startColumnL[i];
    CoinBigIndex end = startColumnL[i + 1];

    totalElements_ += end - start;
    if ( end > start ) {
      firstReal = i;
      CoinBigIndex j;
      for ( j = start; j < end; j++ ) {
	int iRow = indexRowL[j];
	iRow = permute[iRow];
	assert (iRow>firstReal);
	indexRowL[j] = iRow;
      }
    }
  }
  // add in columns of L at the end that are trivial
  // because nothing was eliminated.
  // would be more efficient to put these at the beginning,
  // because already support for these (see baseL)
  for ( i = numberGoodU_+1; i <= numberRows_; i++) {
	  startColumnL[i] = startColumnL[numberGoodU_];
  }

  baseL_ = firstReal;
  numberL_ -= firstReal;
  factorElements_ = totalElements_;
  //can delete pivotRowL_ as not used
  pivotRowL_.conditionalDelete() ;
  //use L for R if room
  CoinBigIndex space = lengthAreaL_ - lengthL_;
  CoinBigIndex spaceUsed = lengthL_ + lengthU_;

  int needed = ( spaceUsed + numberRows_ - 1 ) / numberRows_;

  needed = needed * 2 * maximumPivots_;
  if ( needed < 2 * numberRows_ ) {
    needed = 2 * numberRows_;
  }
  if (numberInColumnPlus_.array()) {
    // Need double the space for R
    space = space/2;
    startColumnR_.conditionalNew( maximumPivots_ + 1 + maximumColumnsExtra_ + 1 );
    CoinBigIndex * startR = startColumnR_.array() + maximumPivots_+1;
    CoinZeroN (startR,(maximumColumnsExtra_+1));
  } else {
    startColumnR_.conditionalNew(maximumPivots_ + 1 );
  }
#ifdef ZEROFAULT
  memset(startColumnR_.array(),'z',(maximumPivots_ + 1)*sizeof(CoinBigIndex));
#endif
  if ( space >= needed ) {
    lengthR_ = 0;
    lengthAreaR_ = space;
    elementR_ = elementL_.array() + lengthL_;
    indexRowR_ = indexRowL_.array() + lengthL_;
  } else {
    lengthR_ = 0;
    lengthAreaR_ = space;
    elementR_ = elementL_.array() + lengthL_;
    indexRowR_ = indexRowL_.array() + lengthL_;
    if ((messageLevel_&4))
      PIPS_APP_LOG_SEV(warning)<<"Factorization may need some increasing area space";
    if ( areaFactor_ ) {
      areaFactor_ *= 1.1;
    } else {
      areaFactor_ = 1.1;
    }
  }
  numberR_ = 0;
}


  //  preProcess.  PreProcesses raw triplet data
  //state is 0 - triplets, 1 - some counts etc , 2 - ..
  // modifications to only pivot first numberColumns1 columns
  // not used
	void
CoinBALPFactorization::preProcessRect ( int state )
{
	int *indexRow = indexRowU_.array();
	int *indexColumn = indexColumnU_.array();
	CoinFactorizationDouble *element = elementU_.array();
	CoinBigIndex numberElements = lengthU_;
	int *numberInRow = numberInRow_.array();
	int *numberInColumn = numberInColumn_.array();
	int *numberInColumnPlus = numberInColumnPlus_.array();
	CoinBigIndex *startRow = startRowU_.array();
	CoinBigIndex *startColumn = startColumnU_.array();
	int numberRows = numberRows_;
	int numberColumns = numberColumns_;

	if (state<4)
		totalElements_ = numberElements;
	//state falls through to next state
	switch ( state ) {
		case 0:			//counts
			{
				CoinZeroN ( numberInRow, numberRows + 1 );
				CoinZeroN ( numberInColumn, maximumColumnsExtra_ + 1 );
				CoinBigIndex i;
				for ( i = 0; i < numberElements; i++ ) {
					int iRow = indexRow[i];
					int iColumn = indexColumn[i];

					numberInRow[iRow]++;
					numberInColumn[iColumn]++;
				}
			}
		case 1:			//sort
			{
				CoinBigIndex i, k;

				i = 0;
				int iColumn;
				// startColumn is initialized
				for ( iColumn = 0; iColumn < numberColumns; iColumn++ ) {
					//position after end of Column
					i += numberInColumn[iColumn];
					startColumn[iColumn] = i;
				}
				// startColumn[0] = numberInColumn[0]
				// but should have startColumn[0] == 0
				// following code reorders triplet elements into column-major form
				for ( k = numberElements - 1; k >= 0; k-- ) {
					int iColumn = indexColumn[k];

					if ( iColumn >= 0 ) {
						CoinFactorizationDouble value = element[k];
						int iRow = indexRow[k];

						indexColumn[k] = -1;
						while ( true ) {
							CoinBigIndex iLook = startColumn[iColumn] - 1;

							startColumn[iColumn] = iLook;
							CoinFactorizationDouble valueSave = element[iLook];
							int iColumnSave = indexColumn[iLook];
							int iRowSave = indexRow[iLook];

							element[iLook] = value;
							indexRow[iLook] = iRow;
							indexColumn[iLook] = -1;
							if ( iColumnSave >= 0 ) {
								iColumn = iColumnSave;
								value = valueSave;
								iRow = iRowSave;
							} else {
								break;
							}
						}			/* endwhile */
					}
				}
			}
		case 2:			//move largest in column to beginning
			//and do row part
			{
				CoinBigIndex i, k;

				i = 0;
				int iRow;
				for ( iRow = 0; iRow < numberRows; iRow++ ) {
					startRow[iRow] = i;
					i += numberInRow[iRow];
				}
				CoinZeroN ( numberInRow, numberRows );
				int iColumn;
				// TODO: we don't need to move largest in column to beginning
				// for non-pivot columns
				for ( iColumn = 0; iColumn < numberColumns; iColumn++ ) {
					int number = numberInColumn[iColumn];

					if ( number ) {
						CoinBigIndex first = startColumn[iColumn];
						CoinBigIndex largest = first;
						int iRowSave = indexRow[first];
						CoinFactorizationDouble valueSave = element[first];
						double valueLargest = fabs ( valueSave );
						int iLook = numberInRow[iRowSave];

						numberInRow[iRowSave] = iLook + 1;
						indexColumn[startRow[iRowSave] + iLook] = iColumn;
						for ( k = first + 1; k < first + number; k++ ) {
							int iRow = indexRow[k];
							int iLook = numberInRow[iRow];

							numberInRow[iRow] = iLook + 1;
							indexColumn[startRow[iRow] + iLook] = iColumn;
							CoinFactorizationDouble value = element[k];
							double valueAbs = fabs ( value );

							if ( valueAbs > valueLargest ) {
								valueLargest = valueAbs;
								largest = k;
							}
						}
						indexRow[first] = indexRow[largest];
						element[first] = element[largest];
						indexRow[largest] = iRowSave;
						element[largest] = valueSave;
					}
				}
			}
		case 3:			//links and initialize pivots
			{
				int *lastRow = lastRow_.array();
				int *nextRow = nextRow_.array();
				int *lastColumn = lastColumn_.array();
				int *nextColumn = nextColumn_.array();

				CoinFillN ( firstCount_.array(), biggerDimension_ + 2, -1 );
				CoinFillN ( pivotColumn_.array(), numberColumns_, -1 );
				CoinZeroN ( numberInColumnPlus, maximumColumnsExtra_ + 1 );
				int iRow;
				for ( iRow = 0; iRow < numberRows; iRow++ ) {
					lastRow[iRow] = iRow - 1;
					nextRow[iRow] = iRow + 1;
					int number = numberInRow[iRow];

					addLink ( iRow, number );
				}
				lastRow[maximumRowsExtra_] = numberRows - 1;
				nextRow[maximumRowsExtra_] = 0;
				lastRow[0] = maximumRowsExtra_;
				nextRow[numberRows - 1] = maximumRowsExtra_;
				startRow[maximumRowsExtra_] = numberElements;
				int iColumn;
				for ( iColumn = 0; iColumn < numberColumns_; iColumn++ ) {
					lastColumn[iColumn] = iColumn - 1;
					nextColumn[iColumn] = iColumn + 1;
					int number = numberInColumn[iColumn];
					if (iColumn < numberColumns1_) 
						addLink ( iColumn + numberRows, number );
				}
				lastColumn[maximumColumnsExtra_] = numberColumns - 1;
				nextColumn[maximumColumnsExtra_] = 0;
				lastColumn[0] = maximumColumnsExtra_;
				if (numberColumns)
					nextColumn[numberColumns - 1] = maximumColumnsExtra_;
				startColumn[maximumColumnsExtra_] = numberElements;
			}
			break;
		case 4:			//move largest in column to beginning
			{
				printf("CASE 4!!\n");
				CoinBigIndex i, k;
				CoinFactorizationDouble * pivotRegion = pivotRegion_.array();
				int iColumn;
				int iRow;
				for ( iRow = 0; iRow < numberRows; iRow++ ) {
					if( numberInRow[iRow]>=0) {
						// zero count
						numberInRow[iRow]=0;
					} else {
						// empty
						//numberInRow[iRow]=-1; already that
					}
				}
				//CoinZeroN ( numberInColumnPlus, maximumColumnsExtra_ + 1 );
				for ( iColumn = 0; iColumn < numberColumns; iColumn++ ) {
					int number = numberInColumn[iColumn];

					if ( number ) {
						// use pivotRegion and startRow for remaining elements
						CoinBigIndex first = startColumn[iColumn];
						CoinBigIndex largest = -1;

						double valueLargest = -1.0;
						int nOther=0;
						k = first;
						CoinBigIndex end = first+number;
						for (  ; k < end; k++ ) {
							int iRow = indexRow[k];
							assert (iRow<numberRows_);
							CoinFactorizationDouble value = element[k];
							if (numberInRow[iRow]>=0) {
								numberInRow[iRow]++;
								double valueAbs = fabs ( value );
								if ( valueAbs > valueLargest ) {
									valueLargest = valueAbs;
									largest = nOther;
								}
								startRow[nOther]=iRow;
								pivotRegion[nOther++]=value;
							} else {
								indexRow[first] = iRow;
								element[first++] = value;
							}
						}
						numberInColumnPlus[iColumn]=first-startColumn[iColumn];
						startColumn[iColumn]=first;
						//largest
						if (largest>=0) {
							indexRow[first] = startRow[largest];
							element[first++] = pivotRegion[largest];
						}
						for (k=0;k<nOther;k++) {
							if (k!=largest) {
								indexRow[first] = startRow[k];
								element[first++] = pivotRegion[k];
							}
						}
						numberInColumn[iColumn]=first-startColumn[iColumn];
					}
				}
				//and do row part
				i = 0;
				for ( iRow = 0; iRow < numberRows; iRow++ ) {
					startRow[iRow] = i;
					int n=numberInRow[iRow];
					if (n>0) {
						numberInRow[iRow]=0;
						i += n;
					}
				}
				for ( iColumn = 0; iColumn < numberColumns; iColumn++ ) {
					int number = numberInColumn[iColumn];

					if ( number ) {
						CoinBigIndex first = startColumn[iColumn];
						for ( k = first; k < first + number; k++ ) {
							int iRow = indexRow[k];
							int iLook = numberInRow[iRow];

							numberInRow[iRow] = iLook + 1;
							indexColumn[startRow[iRow] + iLook] = iColumn;
						}
					}
				}
			}
			// modified 3
			{
				//set markRow so no rows updated
				//CoinFillN ( markRow_.array(), numberRows_, -1 );
				int *lastColumn = lastColumn_.array();
				int *nextColumn = nextColumn_.array();
				CoinFactorizationDouble * pivotRegion = pivotRegion_.array();

				int iRow;
				int numberGood=0;
				startColumnL_.array()[0] = 0;	//for luck
				for ( iRow = 0; iRow < numberRows; iRow++ ) {
					int number = numberInRow[iRow];
					if (number<0) {
						numberInRow[iRow]=0;
						pivotRegion[numberGood++]=slackValue_;
					}
				}
				int iColumn;
				for ( iColumn = 0 ; iColumn < numberColumns_; iColumn++ ) {
					lastColumn[iColumn] = iColumn - 1;
					nextColumn[iColumn] = iColumn + 1;
					int number = numberInColumn[iColumn];
					if ( iColumn < numberColumns1_ ) {
						deleteLink(iColumn+numberRows);
						addLink ( iColumn + numberRows, number );
					}
				}
				lastColumn[maximumColumnsExtra_] = numberColumns - 1;
				nextColumn[maximumColumnsExtra_] = 0;
				lastColumn[0] = maximumColumnsExtra_;
				if (numberColumns)
					nextColumn[numberColumns - 1] = maximumColumnsExtra_;
				startColumn[maximumColumnsExtra_] = numberElements;
			}
	}				/* endswitch */
}



int CoinBALPFactorization::factorizeSquare(int numberOfColumns, 
	CoinBigIndex numberOfElements, const int indicesRow[],
	const int indicesColumn[], const double elements[], double areaFactor)
{
 
 gutsOfDestructor();
  gutsOfInitialize(2);
  if (areaFactor)
    areaFactor_ = areaFactor;
  CoinBigIndex numberElements = 3*numberOfColumns+3*numberOfElements+20000;
  getAreas ( numberOfColumns, numberOfColumns, 0, numberElements, 2*numberElements );
  //copy
  CoinMemcpyN ( indicesRow, numberOfElements, indexRowU_.array() );
  CoinMemcpyN ( indicesColumn, numberOfElements, indexColumnU_.array() );
  int i;
  CoinFactorizationDouble * elementU = elementU_.array();
  for (i=0;i<numberOfElements;i++)
    elementU[i] = elements[i];
  lengthU_ = numberOfElements;
  maximumU_ = numberOfElements;
  numberColumns1_ = numberOfColumns;
  preProcessRect ( 0 );
  factorRect (  );
 
  return status_;
}

