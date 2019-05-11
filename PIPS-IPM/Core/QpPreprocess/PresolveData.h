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
      bool outdated_activities;

      StochVectorHandle actmax_eq;
      StochVectorHandle actmin_eq;

      StochVectorHandle actmax_ineq;
      StochVectorHandle actmin_ineq;

      /* for better MPI communication we allocate one contiguous array in storage and
       * make four SimppleVectors each pointing to a part of it */
      int lenght_array_act_chgs;
      double* array_act_chgs;
      SimpleVector actmax_eq_chgs;
      SimpleVector actmin_eq_chgs;
      SimpleVector actmax_ineq_chgs;
      SimpleVector actmin_ineq_chgs;

   private:
      int my_rank;
      bool distributed;

      // number of children
      int nChildren;
      // objective offset created by presolving
      double objOffset;

      std::vector<COLUMNFORDELETION> linkingVariablesMarkedForDeletion;

   public:
      PresolveData(const sData* sorigprob);
      ~PresolveData();

      sData* finalize();


      /* compute and update activities */
      void recomputeActivities() { recomputeActivities(false); }
      void recomputeActivities(bool linking_only);
      void updateLinkingRowActivities();
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
      // initialize row and column nnz counter
      void initNnzCounter();
      void initialize();
};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */

