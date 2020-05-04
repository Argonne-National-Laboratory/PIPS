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

   bool applyPresolving() override;

private:
   /** limit for rounds of bound strengthening per call of presolver */
   const int limit_iter;
   /** min entry to devide by in order to derive a bound */
   const double limit_entry;
   /** max activity to be devided */
   const double limit_partial_activity;
   /** max bounds proposed from bounds strengthening presolver */
   const double limit_bounds;

   long long tightenings;

   bool strenghtenBoundsInNode(SystemType system_type, int node);
   bool strenghtenBoundsInBlock( SystemType system_type, int node, BlockType block_type);
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_ */
