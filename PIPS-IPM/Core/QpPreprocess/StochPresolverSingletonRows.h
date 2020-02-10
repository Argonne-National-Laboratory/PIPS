/*
 * StochPresolverSingletonRows.h
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_

#include "StochPresolverBase.h"
#include "SparseStorageDynamic.h"
#include "PresolveData.h"

class StochPresolverSingletonRows : public StochPresolverBase
{
public:
   StochPresolverSingletonRows(PresolveData& presData, const sData& origProb);

   ~StochPresolverSingletonRows();

   // remove singleton rows
   virtual void applyPresolving();

private:
   long long removed_rows; 

   bool removeSingletonRow( const INDEX& row );
   void getBoundsAndColFromSingletonRow( const INDEX& row, int& node_col, int& col_idx, double& xlow_new, double& xupp_new);

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
