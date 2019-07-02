/*
 * StochPresolverSingletonColumns.h
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_

#include "StochPresolverBase.h"
#include <limits>


class StochPresolverSingletonColumns : public StochPresolverBase
{
public:
   StochPresolverSingletonColumns(PresolveData& presData, const sData& origProb, StochPostsolver* postsolver);

   ~StochPresolverSingletonColumns();

   // remove singleton columns
   virtual void applyPresolving();

private:
   void countSingletonColumns();
   void initSingletonColumns(int& nSColEq, int& nSCIneq, int& nZeroCostSc, int& nLowerBound, int& nUpperBound,
         int& nBothBounds, int& nNoBounds, int& nSColEqLinkRow, int& nSColIneqLinkrow, int& nSC);
   void initSingletonColsBlock(int it, SimpleVector const * nnzColSimple, int& nSColEq, int& nSCIneq, int& nZeroCostSc,
         int& nLowerBound, int& nUpperBound, int& nBothBounds, int& nNoBounds, int& nSColEqLinkRow, int& nSColIneqLinkrow, int& nSC);
   void synchronizeSumSeveral(int& val0, int& val1, int& val2, int& val3, int& val4, int& val5, int& val6, int&val7, int& val8, int& val9 );

};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_ */
