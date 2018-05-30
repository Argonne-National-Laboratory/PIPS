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
   // pointers to the normalized and copied matrix blocks
   SparseStorageDynamic* norm_Amat;
   SparseStorageDynamic* norm_Bmat;
   SparseStorageDynamic* norm_Cmat;
   SparseStorageDynamic* norm_Dmat;
   SimpleVector* norm_b;
   SimpleVector* norm_c;
   SimpleVector* norm_d;

   // number of rows of the A or C block
   int mA;

   // number of columns of the A or C block
   int nA;

   bool setNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC);
   void countDuplicateRows(StochGenMatrix& matrix, SystemType system_type);
   bool compareCoefficients(SparseStorageDynamic& matrix, int i, int j) const;

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUPLICATEROWS_H_ */
