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

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

private:

   // private methods:
   int initSingletonRows(SystemType system_type);
   int initSingletonRowsBlock(int it, SimpleVector const * nnzRowSimple);
   bool doSingletonRowsA(int& newSREq, int& newSRIneq);
   bool doSingletonRowsC(int& newSREq, int& newSRIneq);

   bool procSingletonRowRoot(StochGenMatrix& stochMatrix, SystemType system_type);
   bool procSingletonRowChild(int it, int& newSR, int& newSRIneq);
   bool procSingletonRowChildAmat(int it, SystemType system_type);
   bool procSingletonRowChildBmat(int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR,
         SystemType system_type);
   bool removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR);

   bool removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx);
   bool removeSingleRowEntryB0Inequality(SparseStorageDynamic& storage, SparseStorageDynamic& storageTransposed, int rowIdx);
   bool procSingletonRowChildInequality(int it, int& newSREq, int& newSRIneq);

   void calculateNewBoundsOnVariable(double& newxlow, double& newxupp, int rowIdx, double aik) const;
   bool storeColValInColAdaptParentAndAdaptOffset(int colIdx, double value, double* g);

   void updateLinkingVarsBounds();
   void getValuesForSR(SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const;

   /** initialize current pointer for matrices and vectors.
    * If it==-1, we are at parent and want block B_0 (Bmat).
    * Returns false if it is a dummy child. */
   bool updateCPForSingletonRow(int it, SystemType system_type);
   bool updateCPForSingletonRowInequalityBChild( int it );
   bool updateCPForSingletonRowEqualityBChild( int it );

};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
