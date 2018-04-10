/*
 * StochPresolverSingletonRows.h
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_

#include "StochPresolverBase.h"
#include <vector>


class StochPresolverSingletonRows : public StochPresolverBase
{
public:
   StochPresolverSingletonRows(PresolveData& presData);

   ~StochPresolverSingletonRows();

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

   // data

private:

   // private methods:
   int initSingletonRows(SystemType system_type);
   int initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple);
   bool doSingletonRowsA(int& newSREq, int& newSRIneq);

   bool procSingletonRowRoot(StochGenMatrix& stochMatrix);
   bool procSingletonRowChild(StochGenMatrix& stochMatrix, int it, int& newSR, int& newSRIneq);
   bool procSingletonRowChildAmat(SparseStorageDynamic& A_mat, int it);
   bool procSingletonRowChildBmat(SparseStorageDynamic& B_mat, int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR);
   bool removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR);

   bool removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx);

   int doSingletonRowsC();
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
