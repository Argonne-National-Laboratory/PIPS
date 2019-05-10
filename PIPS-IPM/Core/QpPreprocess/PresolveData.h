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


      /* stuff for handling the update and changes of activities of certain rows */

      StochVectorHandle actmax_eq;
      StochVectorHandle actmin_eq;

      StochVectorHandle actmax_ineq;
      StochVectorHandle actmin_ineq;

      SimpleVector actmax_eq_chgs;
      SimpleVector actmin_eq_chgs;
      SimpleVector actmax_ineq_chgs;
      SimpleVector actmin_ineq_chgs;


   public:
      PresolveData(const sData* sorigprob);
      ~PresolveData();

      sData* finalize();


      /* compute and update activities */
      void recomputeActivities() { recomputeActivities(false); }
      void recomputeActivities(bool linking_only);
      void updateLinkingRowsActivities();
   private:
      void addActivityOfBlock( const SparseStorageDynamic& matrix, SimpleVector& min_activities, SimpleVector& max_activities,
            const SimpleVector& xlow, const SimpleVector& ixlow, const SimpleVector& xupp, const SimpleVector& ixupp) const;
   public:

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

