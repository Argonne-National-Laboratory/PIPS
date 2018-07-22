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
   void doBoundStrengthParent(SystemType system_type);
   void doBoundStrengthChild(SystemType system_type);
   double computeNewBound(bool rhs, double activity, double matrixEntry, int rowIdx, SystemType system_type);
   void strenghtenBoundsInBlock( SparseStorageDynamic& matrix, bool childBlock,
         int rowIdx, double partMinActivity, double partMaxActivity, SystemType system_type);
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_ */
