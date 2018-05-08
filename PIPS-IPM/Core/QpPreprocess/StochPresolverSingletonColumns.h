/*
 * StochPresolverSingletonColumns.h
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_

#include "StochPresolverBase.h"

class StochPresolverSingletonColumns : public StochPresolverBase
{
public:
   StochPresolverSingletonColumns(PresolveData& presData);

   ~StochPresolverSingletonColumns();

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

private:
   void initSingletonColumns(int& nSColEq, int& nSCIneq, int& nZeroCostSc);
   void initSingletonColsBlock(int it, SimpleVector const * nnzColSimple, int& nSColEq, int& nSCIneq, int& nZeroCostSc);

};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_ */
