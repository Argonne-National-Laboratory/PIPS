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
   long long tightenings;

   bool local_bound_tightenings;
   const unsigned int n_linking_vars;

   std::vector<double> ub_linking_var;
   std::vector<double> lb_linking_var;
   std::vector<INDEX> rows_ub;
   std::vector<INDEX> rows_lb;

   void resetArrays();
   void communicateLinkingVarBounds();

   bool strenghtenBoundsInNode(SystemType system_type, int node);
   bool strenghtenBoundsInBlock( SystemType system_type, int node, BlockType block_type);
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_ */
