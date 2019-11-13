/*
 * StochPresolverBase.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_

#include "StochVector.h"
#include "StochGenMatrix.h"
#include "PresolveData.h"
#include "sData.h"
#include "SystemType.h"
#include "StochPostsolver.h"
#include "pipsport.h"

#include <vector>

class StochPresolverBase
{
public:
   StochPresolverBase(PresolveData& presData, const sData& origProb);
   virtual ~StochPresolverBase();

   // todo return bool whether enough eliminations
   virtual void applyPresolving() = 0;
   void countRowsCols(); // theoretically const but sets pointers

protected:
   void updatePointersForCurrentNode(int node, SystemType system_type);
private:
   void countRowsBlock(int& n_rows_total, int& n_rows_empty, int& n_rows_onesided, int& n_rows_boxed, int& n_rows_fixed, int& n_rows_singleton, SystemType system_type, 
      BlockType block_type) const; 
   void countBoxedColumns( int& n_cols_total, int& n_cols_empty, int& n_cols_free, int& n_cols_onesided, int& n_cols_boxed, int& n_cols_singleton, 
      int& n_cols_orig_free, int& n_cols_orig_free_removed, const SimpleVector& ixlow_orig, const SimpleVector& ixupp_orig,
      BlockType block_type) const;

   void setPointersMatrices(GenMatrixHandle mat, int node);
   void setPointersMatrixBoundsActivities(SystemType system_type, int node);
   void setPointersVarBounds(int node);
   void setPointersObjective(int node);
   void setReductionPointers(SystemType system_type, int node);

protected:
   int my_rank;
   bool distributed;

   /* not owned by the class itself - given from the outside */
   PresolveData& presData;

   // pointers to the currently needed matrices and vectors for presolving
   const sData& origProb;

   const SparseStorageDynamic* currAmat;
   const SparseStorageDynamic* currAmatTrans;
   const SparseStorageDynamic* currBmat;
   const SparseStorageDynamic* currBmatTrans;
   const SparseStorageDynamic* currBlmat;
   const SparseStorageDynamic* currBlmatTrans;

   const SimpleVector* currxlowParent;
   const SimpleVector* currIxlowParent;
   const SimpleVector* currxuppParent;
   const SimpleVector* currIxuppParent;
   const SimpleVector* currxlowChild;
   const SimpleVector* currIxlowChild;
   const SimpleVector* currxuppChild ;
   const SimpleVector* currIxuppChild;

   const SimpleVector* currEqRhs;
   const SimpleVector* currIneqLhs;
   const SimpleVector* currIclow;
   const SimpleVector* currIneqRhs;
   const SimpleVector* currIcupp;
   const SimpleVector* currEqRhsLink;
   const SimpleVector* currIneqLhsLink;
   const SimpleVector* currIclowLink;
   const SimpleVector* currIneqRhsLink;
   const SimpleVector* currIcuppLink;

   const SimpleVector* currgParent;
   const SimpleVector* currgChild;

   const SimpleVectorBase<int>* currNnzRow;
   const SimpleVectorBase<int>* currNnzRowLink;

   const SimpleVectorBase<int>* currNnzColParent;
   const SimpleVectorBase<int>* currNnzColChild;

   /** the number of children */
   int nChildren;
   /** number of entry eliminations on this process in the current elimination routine */
   int localNelims;
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
