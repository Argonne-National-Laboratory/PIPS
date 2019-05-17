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
#include <vector>

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
   void countRowsBlock(int& n_rows, int& n_ranged_rows, int& n_fixed_rows, int& n_singleton_rows, SystemType system_type, BlockType block_type) const;
   void countBoxedColumns(int& nBoxCols, int& nColsTotal, int& nFreeVars, int& nOnesidedVars, int& nSingletonVars, int& nSingletonVarsImpliedFree,
         const SimpleVector& ixlow_orig, const SimpleVector& xlow_orig, const SimpleVector& ixupp_orig, const SimpleVector& xupp_orig, BlockType block_type) const;
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

   const SimpleVector* currActMax;
   const SimpleVector* currActMin;
   const SimpleVector* currActMaxLink;
   const SimpleVector* currActMinLink;

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

   const SimpleVector* currNnzRow;
//   const SimpleVector* currRedRowLink;
   const SimpleVector* currNnzRowLink;

//   const SimpleVector* currRedColParent;
   const SimpleVector* currNnzColParent;
   const SimpleVector* currNnzColChild;

   /** the number of children */
   int nChildren;
   /** number of entry eliminations on this process in the current elimination routine */
   int localNelims;

   /* objective offset resulting from local presolving */
   double indivObjOffset;
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
