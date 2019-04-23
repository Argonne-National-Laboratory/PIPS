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
} COLUMNTOADAPT;

struct col_is_smaller
{
    bool operator()(const COLUMNTOADAPT& x, const COLUMNTOADAPT& y) const
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

      //StochGenMatrix& Apres;
      //StochGenMatrix& Cpres;

      PresolveData(const sData* sorigprob);
      ~PresolveData();

      sData* finalize();

      bool combineColAdaptParent();

      bool reductionsEmpty();
      void resetRedCounters();

      void resetBlocks();
      // todo getter, setter for element access of nnz counter???
      int getNChildren() const;
      double getObjOffset() const;
      double addObjOffset(double addOffset);
      void setObjOffset(double offset);

      int getSingletonRow(int i) const;
      int getNumberSR() const;
      void addSingletonRow(int value);
      void setSingletonRow(int i, int value);
      void clearSingletonRows();
      int getSingletonRowIneq(int i) const;
      int getNumberSRIneq() const;
      void addSingletonRowIneq(int value);
      void setSingletonRowIneq(int i, int value);
      void clearSingletonRowsIneq();

      void setBlocks(int i, double value);
      double getBlocks(int i) const;
      void setBlocksIneq(int i, double value);
      double getBlocksIneq(int i) const;

      COLUMNTOADAPT getColAdaptParent(int i) const;
      int getNumberColAdParent() const;
      void addColToAdaptParent(COLUMNTOADAPT colToAdapt);
      void clearColAdaptParent();

   private:
      // number of children
      int nChildren;
      // objective offset created by presolving
      double objOffset;

      // variables used for singleton row elimination:
      /** vector containing the row indices of singleton rows */
      std::vector<int> singletonRows;
      std::vector<int> singletonRowsIneq;

      /** array of length nChildren+3 to store start indices for singletonRows
       * that correspond to the correct block. As blocks[0] represents the parent block,
       * the child block 'it' is accessed using the index 'it+1'.
       * The linking-row block is accessed using the index nChildren+2. */
      int* blocks;
      int* blocksIneq;

      /** vector containing the column indices of entries that were found during the
       * singleton row routine. Along with the column index, the value needed for
       * adaptation is stored. */
      std::vector<COLUMNTOADAPT> colAdaptParent;

      // initialize row and column nnz counter
      void initNnzCounter();
      void initialize();

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */

