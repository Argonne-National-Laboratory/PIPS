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
   bool newBoundsImplyInfeasible(double newxlow, double newxupp, int colIdx,
         double* ixlow, double* ixupp, double* xlow, double* xupp) const;
   bool newBoundsFixVariable(double& value, double newxlow, double newxupp, int colIdx,
         double* ixlow, double* ixupp, double* xlow, double* xupp) const;
   bool newBoundsTightenOldBounds(double newxlow, double newxupp, int colIdx,
         double* ixlow, double* ixupp, double* xlow, double* xupp) const;
   bool storeColValInColAdaptParentAndAdaptOffset(int colIdx, double value, double* g);
   bool storeNewBoundsParent(int colIdx, double newxlow, double newxupp);

   bool combineNewBoundsParent();
   void updateLinkingVarsBounds();
   void getValuesForSR(SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const;
   XBOUNDS getNewBoundsParent(int i) const;
   void setNewBoundsParent(int i, int colIdx, double newxlow, double newxupp);
   int getNumberNewBoundsParent() const;
   void addNewBoundsParent(XBOUNDS newXBounds);
   void clearNewBoundsParent();

   void setNewXBounds(int colIdx, double newxlow, double newxupp, double* ixlow, double* xlow, double* ixupp, double* xupp) const;
   void synchronizeNumberSR(int& newSREq, int& newSRIneq) const;
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERSINGLETONROWS_H_ */
