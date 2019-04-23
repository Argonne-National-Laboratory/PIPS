/*
 * StochPresolverSingletonRows.h
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_

#include "StochPresolverBase.h"
#include <vector>

class StochPresolverSingletonRows : public StochPresolverBase
{
public:
   StochPresolverSingletonRows(PresolveData& presData);

   ~StochPresolverSingletonRows();

   // remove singleton rows
   virtual void applyPresolving();

private:
   void doSingletonRows(int& n_sing_sys, int& n_sing_other_sys, SystemType system_type);
   void doSingletonLinkRows(int& newSREq, int& newSRIneq);

   void procSingletonRowRoot(SystemType system_type);
   void procSingletonRowChild(int it, int& n_singleton_sys, int& n_singleton_other_sys, SystemType system_type);

   void procSingletonRowChildBmat(int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR, SystemType system_type);
   void removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR);

   void processSingletonBlock(SystemType system_type, BlockType block_type, int node);

   void calculateNewBoundsOnVariable(double& new_xlow, double& new_xupp, const double& iclow, const double& clow,
         const double& icupp, const double& cupp, double aik) const;

   void updateLinkingVarsBounds();
   void getValuesForSR(SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const;

   // todo move
   bool tightenBounds(double new_xlow, double new_xupp, double& ixlow, double& old_xlow, double& ixupp, double& old_xupp) const;
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
