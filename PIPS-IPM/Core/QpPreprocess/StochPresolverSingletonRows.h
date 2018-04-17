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
#include <limits>

typedef struct
{
   int colIdx;
   double newxlow;
   double newxupp;
} XBOUNDS;

struct xbounds_col_is_smaller
{
    bool operator()(const XBOUNDS& x, const XBOUNDS& y) const
    {
        return x.colIdx < y.colIdx;
    }
};

class StochPresolverSingletonRows : public StochPresolverBase
{
public:
   StochPresolverSingletonRows(PresolveData& presData);

   ~StochPresolverSingletonRows();

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

   // data
   std::vector<XBOUNDS> newBoundsParent;

private:

   // private methods:
   int initSingletonRows(SystemType system_type);
   int initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple);
   bool doSingletonRowsA(int& newSREq, int& newSRIneq);
   bool doSingletonRowsC(int& newSREq, int& newSRIneq);

   bool procSingletonRowRoot(StochGenMatrix& stochMatrix, SystemType system_type);
   bool procSingletonRowChild(StochGenMatrix& stochMatrix, int it, int& newSR, int& newSRIneq);
   bool procSingletonRowChildAmat(int it, SystemType system_type);
   bool procSingletonRowChildBmat(SparseStorageDynamic& B_mat, int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR,
         SystemType system_type);
   bool removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR);

   bool removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx);
   bool removeSingleRowEntryB0Inequality(SparseStorageDynamic& storage, SparseStorageDynamic& storageTransposed, int rowIdx);
   bool procSingletonRowChildInequality(StochGenMatrix& stochMatrix, int it);

   void calculateNewBoundsOnVariable(double& newxlow, double& newxupp, int rowIdx, double aik);
   bool newBoundsImplyInfeasible(double newxlow, double newxupp, int colIdx,
         double* ixlow, double* ixupp, double* xlow, double* xupp);
   bool newBoundsFixVariable(double& value, double newxlow, double newxupp, int colIdx,
         double* ixlow, double* ixupp, double* xlow, double* xupp);
   bool newBoundsTightenOldBounds(double newxlow, double newxupp, int colIdx,
         double* ixlow, double* ixupp, double* xlow, double* xupp);
   bool storeColValInColAdaptParentAndAdaptOffset(int colIdx, double value, double* g);
   bool storeNewBoundsParent(int colIdx, double newxlow, double newxupp);

   bool combineNewBoundsParent();
   void updateLinkingVarsBounds();
   void getValuesForSR(SparseStorageDynamic& storage, int rowIdx, int& colIdx, double& aik);
   XBOUNDS getNewBoundsParent(int i);
   void setNewBoundsParent(int i, int colIdx, double newxlow, double newxupp);
   int getNumberNewBoundsParent();
   void addNewBoundsParent(XBOUNDS newXBounds);
   void clearNewBoundsParent();

   void setNewXBounds(int colIdx, double newxlow, double newxupp, double* ixlow, double* xlow, double* ixupp, double* xupp);
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
