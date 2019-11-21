/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonColumns.h"


StochPresolverSingletonColumns::StochPresolverSingletonColumns(PresolveData& presData, const sData& origProb, StochPostsolver* postsolver)
   : StochPresolverBase(presData, origProb, postsolver)
{
 // todo
}

StochPresolverSingletonColumns::~StochPresolverSingletonColumns()
{
 // todo
}


void StochPresolverSingletonColumns::applyPresolving()
{

}

void StochPresolverSingletonColumns::countSingletonColumns()
{
   
}

void StochPresolverSingletonColumns::initSingletonColumns(int& nSColEq, int& nSColIneq, int& nZeroCostSc,
      int& nLowerBound, int& nUpperBound, int& nBothBounds, int& nNoBounds, int& nSColEqLinkRow, int& nSColIneqLinkrow, int& nSC)
{
  
}

void StochPresolverSingletonColumns::initSingletonColsBlock(int it, SimpleVector const * nnzColSimple,
      int& nSColEq, int& nSColIneq, int& nZeroCostSC, int& nLowerBound, int& nUpperBound, int& nBothBounds, int& nNoBounds,
      int& nSColEqLinkRow, int& nSColIneqLinkrow, int& nSC)
{
  
}

void StochPresolverSingletonColumns::synchronizeSumSeveral(int& val0, int& val1, int& val2, int& val3, int& val4, int& val5, int& val6,
      int&val7, int& val8, int& val9 )
{
   
}

