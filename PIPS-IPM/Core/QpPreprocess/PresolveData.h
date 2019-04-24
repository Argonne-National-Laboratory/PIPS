/*
 * PresolveData.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_

#include "sData.h"
#include <algorithm>


typedef struct
{
   int colIdx;
   double val;
} COLUMNFORDELETION;

struct col_is_smaller
{
    bool operator()(const COLUMNFORDELETION& x, const COLUMNFORDELETION& y) const
    {
        return x.colIdx < y.colIdx;
    }
};

class PresolveData
{

   public:
      sData* presProb;


      // number of non-zero elements of each row / column
      StochVectorHandle nRowElemsA;
      StochVectorHandle nRowElemsC;
      StochVectorHandle nColElems;

      // number of removed elements of each row / column
      StochVectorHandle redRowA; // todo maybe rename ?
      StochVectorHandle redRowC;
      StochVectorHandle redCol;


   public:
      PresolveData(const sData* sorigprob);
      ~PresolveData();

      sData* finalize();

      bool combineColAdaptParent();

      bool reductionsEmpty();
      void resetRedCounters();

      // todo getter, setter for element access of nnz counter???
      int getNChildren() const;
      double getObjOffset() const;
      double addObjOffset(double addOffset);
      void setObjOffset(double offset);

      COLUMNFORDELETION getColAdaptParent(int i) const;
      int getNumberColAdParent() const;
      void addColToAdaptParent(COLUMNFORDELETION colToAdapt);
      void clearColAdaptParent();

   private:
      // number of children
      int nChildren;
      // objective offset created by presolving
      double objOffset;

      std::vector<COLUMNFORDELETION> linkingVariablesMarkedForDeletion;

      // initialize row and column nnz counter
      void initNnzCounter();
      void initialize();

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */

