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
   StochPresolverSingletonColumns(PresolveData& presData, const sData& origProb);

   ~StochPresolverSingletonColumns();

   // remove singleton columns
   virtual void applyPresolving();

private:
   long long removed_cols;

   bool removeSingletonColumn(int node, int col);

};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_ */
