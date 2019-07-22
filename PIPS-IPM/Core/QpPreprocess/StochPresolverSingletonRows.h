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
#include <fstream>

class StochPresolverSingletonRows : public StochPresolverBase
{
public:
   StochPresolverSingletonRows(PresolveData& presData, const sData& origProb);

   ~StochPresolverSingletonRows();

   // remove singleton rows
   virtual void applyPresolving();

private:
   long long removed_rows; 

   bool removeSingletonRow(SystemType system_type, int node, int row_idx);
   void getBoundsAndColFromSingletonRow( SystemType system_type, int& node, int row_idx, BlockType block_type, int& col_idx, double& ubx, double& lbx );

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
