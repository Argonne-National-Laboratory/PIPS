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
   StochPresolverBoundStrengthening(PresolveData& presData);

   virtual ~StochPresolverBoundStrengthening();

   virtual void applyPresolving();

private:
   bool setCPforBounds(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   void doBoundStrengthParent(SystemType system_type);
   void doBoundStrengthChild(int it, SystemType system_type);
   double computeNewBound(bool rhs, double activity, double matrixEntry, int rowIdx, SystemType system_type) const;
   void strenghtenBoundsInBlock( SparseStorageDynamic& matrix, bool childBlock, int rowIdx, double partMinActivity, double partMaxActivity,
         SystemType system_type, bool atRoot, std::vector<COLUMNFORDELETION>* colAdaptLinkBlock);
   void setNewBoundsIfTighter(int index, double new_low, double new_upp, SimpleVector& ilow, SimpleVector& low, SimpleVector& iupp, SimpleVector& upp) const;
   bool checkNewBoundTightens(bool uppLow, int colIdx, double newBound, const SimpleVector& ixbound, const SimpleVector& xbound ) const;
   void strengthenLinkingVarsBounds();
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_ */
