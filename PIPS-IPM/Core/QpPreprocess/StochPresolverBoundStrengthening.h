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

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

private:

};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBOUNDSTRENGTHENING_H_ */
