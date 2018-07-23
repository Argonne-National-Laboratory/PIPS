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

   ~StochPresolverBoundStrengthening();

   virtual bool applyPresolving(int& nelims);

private:
   bool setCPforBounds(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   void computeActivityBlockwise( SparseStorageDynamic& matrix, int rowIdx, int colIdx,
         double& infRow, double& supRow,
         SimpleVector& xlow, SimpleVector& ixlow, SimpleVector& xupp, SimpleVector& ixupp);
   bool doBoundStrengthParent(SystemType system_type);
   bool doBoundStrengthChild(SystemType system_type);
   double computeNewBound(bool rhs, double activity, double matrixEntry, int rowIdx, SystemType system_type);
   bool strenghtenBoundsInBlock( SparseStorageDynamic& matrix, bool childBlock,
         int rowIdx, double partMinActivity, double partMaxActivity, SystemType system_type, bool atRoot);
   void setNewBoundsIfTighter(int index, double new_low, double new_upp,
         SimpleVector& ilow, SimpleVector& low, SimpleVector& iupp, SimpleVector& upp);
   bool checkNewBoundTightens(bool uppLow, int colIdx, double newBound,
         SimpleVector& ixbound, SimpleVector& xbound ) const;
   void strengthenLinkingVarsBounds();
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_ */
