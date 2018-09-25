/*
 * StochPresolverSingletonRows.h
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_

#include "StochPresolverBase.h"

class StochPresolverSingletonRows : public StochPresolverBase
{
public:
   StochPresolverSingletonRows(PresolveData& presData);

   ~StochPresolverSingletonRows();

   // remove singleton rows
   virtual void applyPresolving();

private:

   // private methods:
   int initSingletonRows(SystemType system_type);
   int initSingletonRowsBlock(int it, SimpleVector const * nnzRowSimple);
   void doSingletonRowsA(int& newSREq, int& newSRIneq);
   void doSingletonRowsC(int& newSREq, int& newSRIneq);
   void doSingletonLinkRows(int& newSREq, int& newSRIneq);

   void procSingletonRowRoot(StochGenMatrix& stochMatrix, SystemType system_type);
   void procSingletonRowChildEquality(int it, int& newSR, int& newSRIneq);
   void procSingletonRowChildAmat(int it, SystemType system_type);
   void procSingletonRowChildBmat(int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR,
         SystemType system_type);
   void removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR);

   void removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx);
   void removeSingleRowEntryB0Inequality(SparseStorageDynamic& storage, SparseStorageDynamic& storageTransposed, int rowIdx);
   void procSingletonRowChildInequality(int it, int& newSREq, int& newSRIneq);

   void calculateNewBoundsOnVariable(double& newxlow, double& newxupp, int rowIdx, double aik) const;

   void updateLinkingVarsBounds();
   void getValuesForSR(SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const;

   /** initialize current pointer for matrices and vectors.
    * If it==-1, we are at parent and want block B_0 (Bmat).
    * Returns false if it is a dummy child. */
   bool updateCPForSingletonRow(int it, SystemType system_type);

};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
