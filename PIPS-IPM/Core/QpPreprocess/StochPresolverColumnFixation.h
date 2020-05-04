/*
 * StochPresolverColumnFixation.h
 *
 *  Created on: 15.05.2019
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERCOLUMNFIXATION_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERCOLUMNFIXATION_H_

#include "StochPresolverBase.h"

class StochPresolverColumnFixation: public StochPresolverBase
{
   public:

      StochPresolverColumnFixation(PresolveData& presData, const sData& origProb);

      virtual ~StochPresolverColumnFixation();

      void applyPresolving() override;

   private:
      /** limit on the possible impact a column can have on the problem */
      const double limit_fixation_max_fixing_impact;

      int fixed_columns;

};








#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERCOLUMNFIXATION_H_ */
