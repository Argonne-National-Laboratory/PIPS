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
   bool findRowForColumnSingleton( SystemType& system_type, int& node_row, int& row, bool& linking, const int& node_col, const int& col );
   bool findRowForLinkingSingleton( SystemType& system_type, int& node_row, int& row, bool& linking, const int& col);
   bool findRowForLinkingSingletonInSystem( SystemType system_type, int& node_row, int& row, bool& linking, const int& col);
   bool findRowForNonlinkingSingelton( SystemType& system_type, int& row, bool& linking, const int& node_col, const int& col);
   bool findRowForSingletonColumnInMatrix( const SparseStorageDynamic& mat, int& row, const int& col);


   void checkColImpliedFree( SystemType system_type, int node_row, int row, bool linking_row, int node_col, int col, bool& lb_implied_free, bool& ub_implied_free);
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONCOLUMNS_H_ */
