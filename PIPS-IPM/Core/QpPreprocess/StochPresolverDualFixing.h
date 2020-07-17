/*
 * StochPresolverDualFixing.h
 *
 *  Created on: 02.05.2018
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUALFIXING_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUALFIXING_H_

#include "StochPresolverBase.h"

class StochPresolverDualFixing : public StochPresolverBase
{
   public:
      StochPresolverDualFixing(PresolveData& pres_data, const sData& orig_prob);

      virtual ~StochPresolverDualFixing();

      // check for dual fixing in all columns
      bool applyPresolving() override;

   private:
      int fixed_columns;


      int applyDualFixingNode(int node);
};


#endif /* PIPS_IPM_CORE_STOCHPRESOLVERDUALFIXING_H */
