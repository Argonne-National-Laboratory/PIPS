/*
 * StochPresolverDuplicateRows.h
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUPLICATEROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUPLICATEROWS_H_

#include "StochPresolverBase.h"

#include <boost/unordered_set.hpp>


class StochPresolverDuplicateRows : public StochPresolverBase
{
public:
   StochPresolverDuplicateRows(PresolveData& presData);

   ~StochPresolverDuplicateRows();

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

private:
   void countDuplicateRows(StochGenMatrix& matrix, SystemType system_type);
   bool compareCoefficients(SparseStorageDynamic& matrix, int i, int j) const;

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUPLICATEROWS_H_ */
