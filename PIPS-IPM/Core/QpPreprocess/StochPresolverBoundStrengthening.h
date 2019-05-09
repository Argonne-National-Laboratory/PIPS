/*
 * StochPresolverBoundStrengthening.h
 *
 *  Created on: 28.05.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_

#include "StochPresolverBase.h"

class StochPresolverBoundStrengthening : public StochPresolverBase
{
public:
   StochPresolverBoundStrengthening(PresolveData& presData, const sData& origProb);

   virtual ~StochPresolverBoundStrengthening();

   virtual void applyPresolving();

private:
   void doBoundStrengthParent(SystemType system_type);
   void doBoundStrengthChild(int it, SystemType system_type);
   void strenghtenBoundsInBlock( SparseStorageDynamic& matrix, BlockType block_type, int rowIdx, double partMinActivity, double partMaxActivity,
         SystemType system_type, int node);
   double computeNewBound(const SimpleVector& bounds, double activity, double matrixEntry, int rowIdx) const;
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_ */
