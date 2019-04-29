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
   StochPresolverSingletonRows(PresolveData& presData);

   ~StochPresolverSingletonRows();

   // remove singleton rows
   virtual void applyPresolving();

private:
   void doSingletonRows(int& n_sing_sys, int& n_sing_other_sys, SystemType system_type);

   void procSingletonRowRoot(SystemType system_type);
   void procSingletonRowChild(int it, int& n_singleton_sys, int& n_singleton_other_sys, SystemType system_type);

   void processSingletonBlock(SystemType system_type, BlockType block_type, int node);

   void calculateNewBoundsOnVariable(double& new_xlow, double& new_xupp, const double& iclow, const double& clow,
         const double& icupp, const double& cupp, double aik) const;

   void getValuesForSR(SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const;

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
