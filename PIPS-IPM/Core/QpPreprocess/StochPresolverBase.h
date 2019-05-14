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

   void abortInfeasible(MPI_Comm comm) const;
   void countRowsCols(); // theoretically const but sets pointers
private:
   void countRowsBlock(int& n_rows, int& n_ranged_rows, int& n_fixed_rows, int& n_singleton_rows, SystemType system_type, BlockType block_type) const;
   void countBoxedColumns(int& nBoxCols, int& nColsTotal, int& nFreeVars, int& nOnesidedVars, int& nSingletonVars, int& nSingletonVarsImpliedFree,
         const SimpleVector& ixlow_orig, const SimpleVector& xlow_orig, const SimpleVector& ixupp_orig, const SimpleVector& xupp_orig, BlockType block_type) const;
protected:
   void updatePointersForCurrentNode(int node, SystemType system_type);
private:
   void setPointersMatrices(GenMatrixHandle mat, int node);
   void setPointersMatrixBounds(SystemType system_type, int node);
   void setPointersVarBounds(int node);
   void setPointersObjective(int node);
   void setReductionPointers(SystemType system_type, int node);
protected:
   bool hasLinking(SystemType system_type) const;

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

   // todo do we want to make these adjustable?
   static const double feastol = 1.0e-6; // was 1.0e-6
   static const double infinity = 1.0e30;
   // todo rename for more clarity
   static const double tolerance1 = 1.0e-3;  // for model cleanup // was 1.0e-3
   static const double tolerance2 = 1.0e-2;  // for model cleanup // was 1.0e-2
   static const double tol_matrix_entry = 1.0e-10;//1.0e-10; // for model cleanup // was 1.0e-10
   static const double tolerance4 = 1.0e-12; // for variable fixing
   static const double limit1 = 1.0e3;   // for bound strengthening
   static const double limit2 = 1.0e8;   // for bound strengthening
   static const int maxIterSR = 10;
   static const double tol_compare_double = 1.0e-8;

   /* not owned by the class itself - given from the outside */
   PresolveData& presData;

   // pointers to the currently needed matrices and vectors for presolving
   const sData& origProb;

   SparseStorageDynamic* currAmat;
   SparseStorageDynamic* currAmatTrans;
   SparseStorageDynamic* currBmat;
   SparseStorageDynamic* currBmatTrans;
   SparseStorageDynamic* currBlmat;
   SparseStorageDynamic* currBlmatTrans;

   SimpleVector* currxlowParent;
   SimpleVector* currIxlowParent;
   SimpleVector* currxuppParent;
   SimpleVector* currIxuppParent;
   SimpleVector* currxlowChild;
   SimpleVector* currIxlowChild;
   SimpleVector* currxuppChild ;
   SimpleVector* currIxuppChild;

   SimpleVector* currEqRhs;
   SimpleVector* currIneqLhs;
   SimpleVector* currIclow;
   SimpleVector* currIneqRhs;
   SimpleVector* currIcupp;
   SimpleVector* currEqRhsLink;
   SimpleVector* currIneqLhsLink;
   SimpleVector* currIclowLink;
   SimpleVector* currIneqRhsLink;
   SimpleVector* currIcuppLink;

   SimpleVector* currgParent;
   SimpleVector* currgChild;

   SimpleVector* currNnzRow;
   SimpleVector* currRedRowLink;
   SimpleVector* currNnzRowLink;

   SimpleVector* currRedColParent;
   SimpleVector* currNnzColParent;
   SimpleVector* currNnzColChild;

   /** the number of children */
   int nChildren;
   /** number of entry eliminations on this process in the current elimination routine */
   int localNelims;

   std::vector<XBOUNDS> newBoundsParent;

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
//   void allreduceAndUpdate(MPI_Comm comm, SimpleVector& adaptionsVector, SimpleVector& baseVector);
//   bool removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx, double& m) const;
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
//   bool newBoundsFixVariable(double& value, double newxlow, double newxupp, int colIdx,
//         const double* ixlow, const double* ixupp, const double* xlow, const double* xupp) const;
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
//
//protected:
//   void countSingletonRows(int& n_singletons_equality, int& n_singletons_inequality) const;
//private:
//   void countSingletonRowsSystem(int& n_singletons, SystemType system_type) const;
//
//protected:
//   void allreduceAndApplyRhsLhsReductions(SystemType system_type);
//   void allreduceAndUpdateVarBounds();
//
//private:
//   void setVarboundsToInftyForAllreduce() const;
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
