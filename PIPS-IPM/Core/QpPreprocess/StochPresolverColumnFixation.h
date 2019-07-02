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

      virtual void applyPresolving();

   private:
      int fixed_columns;

};








#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERCOLUMNFIXATION_H_ */
