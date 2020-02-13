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

   virtual void applyPresolving();

private:
   long long removed_cols;

   bool removeSingletonColumn( const INDEX& col );
   INDEX findRowForColumnSingleton( const INDEX& col, bool& found );
   INDEX findRowForLinkingSingleton( int col, bool& found );
   INDEX findRowForLinkingSingletonInSystem( int col, SystemType system_type, bool& found );
   INDEX findRowForNonlinkingSingelton( const INDEX& col, bool& found );
   bool findRowForSingletonColumnInMatrix( const SparseStorageDynamic& mat, int& row, const int& col );

   void checkColImpliedFree( const INDEX& col, const INDEX& row, bool& lb_implied_free, bool& ub_implied_free );
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_ */
