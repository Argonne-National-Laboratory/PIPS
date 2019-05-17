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
private:
   void countRowsBlock(int& n_rows, int& n_ranged_rows, int& n_fixed_rows, int& n_singleton_rows, SystemType system_type, BlockType block_type) const;
   void countBoxedColumns(int& nBoxCols, int& nColsTotal, int& nFreeVars, int& nOnesidedVars, int& nSingletonVars, int& nSingletonVarsImpliedFree,
         const SimpleVector& ixlow_orig, const SimpleVector& xlow_orig, const SimpleVector& ixupp_orig, const SimpleVector& xupp_orig, BlockType block_type) const;
protected:
   void updatePointersForCurrentNode(int node, SystemType system_type);
private:
   void setPointersMatrices(GenMatrixHandle mat, int node);
   void setPointersMatrixBoundsActivities(SystemType system_type, int node);
   void setPointersVarBounds(int node);
   void setPointersObjective(int node);
   void setReductionPointers(SystemType system_type, int node);
protected:
public:
   /** checks if all processes have the same root node data.
    *
    *  Applies a MIP_Reduce on all root node data (e.g. A_0, C_0, b_0, f_0...) and compares that
    *  to the root node data in process with id = 0.
    */
//   bool checkRootNodeDataInSync(const sData& sData) const;

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

//   std::vector<XBOUNDS> newBoundsParent;

   /* objective offset resulting from local presolving */
   double indivObjOffset;


   /* operations that actually change the presolved problem - so all reductions etc */

   /* swap two entries in the SparseStorageDynamic format */
   // todo move to presData
//   void updateAndSwap( SparseStorageDynamic* storage, int rowidx, int& indexK, int& rowEnd, double* redCol, int& nelims, bool linking = false);
//   // methods to update the transposed matrix:
//   void updateTransposedSubmatrix(SparseStorageDynamic* transStorage, std::vector<std::pair<int,int> >& elements) const;
//   void updateLinkingVarsBlocks(int& newSREq, int& newSRIneq); // todo check
//   void setNewBounds(int index, double new_low, double new_upp,
//         double* ilow, double* low, double* iupp, double* upp) const;
//   void setNewBound(int index, double new_bound,
//         SimpleVector* bound_vector, SimpleVector* i_bound_vector) const;
//   void clearRow(SparseStorageDynamic& storage, const int rowIdx) const;
//   void removeRow(int rowIdx, SparseStorageDynamic& Ablock, SparseStorageDynamic& AblockTrans,
//         SparseStorageDynamic* Bblock, SparseStorageDynamic* BblockTrans, SimpleVector& nnzRow,
//         SimpleVector& redColParent, SimpleVector* nnzColChild);
//   void removeRowInBblock(int rowIdx, SparseStorageDynamic* Bblock,
//         SparseStorageDynamic* BblockTrans, SimpleVector* nnzColChild);
//private:
//   void adaptChildBmat( std::vector<COLUMNFORDELETION> const & colAdaptBlock, SystemType system_type, int& newSR);
//
//
//
//   bool newBoundsTightenOldBounds(double new_low, double new_upp, int index,
//         double* ilow, double* iupp, double* low, double* upp) const;
//


//
//private:
//   void countBoxedColumns(int& nBoxCols, int& nColsTotal, int& nFreeVars, int& nOnesidedVars, int& nSingletonVars, int& nSingletonVarsImpliedFree,
//         const SimpleVector& ixlow_orig, const SimpleVector& xlow_orig, const SimpleVector& ixupp_orig, const SimpleVector& xupp_orig, BlockType block_type) const;
//
//protected:
//   bool variableFixationValid(double fixation_value, const double& ixlow, const double& xlow, const double& ixupp, const double& xupp, bool print_message = false) const;
//   bool tightenBounds(double new_xlow, double new_xupp, double& ixlow, double& old_xlow, double& ixupp, double& old_xupp) const;
//
//
//
//   void synchronize(int& value) const; // todo move
//
//public:
//
//   int colAdaptLinkVars(int it, SystemType system_type);
//   int colAdaptBl0(SystemType system_type);
//
//   void storeColValInColAdaptParent(int colIdx, double value);
//   bool newBoundsImplyInfeasible(double newxlow, double newxupp, int colIdx,
//      const double* ixlow, const double* ixupp, const double* xlow, const double* xupp) const;
//
//
//   // todo use these instead of all reduce min and max? what is more expensive on average?
//   void storeNewBoundsParent(int colIdx, double newxlow, double newxupp);
//   void combineNewBoundsParent();
//   XBOUNDS getNewBoundsParent(int i) const;
//   void setNewBoundsParent(int i, int colIdx, double newxlow, double newxupp);
//   int getNumberNewBoundsParent() const;
//   void addNewBoundsParent(XBOUNDS newXBounds);
//   void clearNewBoundsParent();
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
