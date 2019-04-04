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
#include <vector>

// todo maybe make these std::pair?
// eg typedef MTRXENTRY std::pair<int,int>
typedef struct
{
   int rowIdx;
   int colIdx;
} MTRXENTRY;

typedef struct
{
   int start;
   int end;
} BLOCKS;

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

enum SystemType {EQUALITY_SYSTEM, INEQUALITY_SYSTEM};
enum BlockType {LINKING_VARS_BLOCK, CHILD_BLOCK};


class StochPresolverBase
{
public:
   StochPresolverBase(PresolveData& presData);
   virtual ~StochPresolverBase();

   // todo return bool whether enough eliminations
   virtual void applyPresolving() = 0;

   bool verifyNnzcounters();
   void countRowsCols(); // theoretically const but sets pointers

   /** checks if all processes have the same root node data.
    *
    *  Applies a MIP_Reduce on all root node data (e.g. A_0, C_0, b_0, f_0...) and compares that
    *  to the root node data in process with id = 0.
    */
   bool checkRootNodeDataInSync(const sData& sData) const;

protected:
   // todo do we want to make these adjustable?
   static const double feastol = 1.0e-6;
   static const double infinity = 1.0e30;
   // todo rename for more clarity
   static const double tolerance1 = 1.0e-3;  // for model cleanup
   static const double tolerance2 = 1.0e-2;  // for model cleanup
   static const double tol_matrix_entry = 1.0e-10; // for model cleanup
   static const double tolerance4 = 1.0e-12; // for variable fixing
   static const double limit1 = 1.0e3;   // for bound strengthening
   static const double limit2 = 1.0e8;   // for bound strengthening
   static const int maxIterSR = 20;
   static const double tol_compare_double = 1.0e-8;

   /* not owned by the class itself - given from the outside */
   PresolveData& presData;

   // pointers to the currently needed matrices and vectors for presolving
   sData* presProb;


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
   SimpleVector* currIneqLhs; // todo rename ?
   SimpleVector* currIclow;
   SimpleVector* currIneqRhs; // TODO rename ?
   SimpleVector* currIcupp;
   SimpleVector* currEqRhsLink; // TODO rename ?
   SimpleVector* currIneqLhsLink; // TODO rename ?
   SimpleVector* currIclowLink;
   SimpleVector* currIneqRhsLink; // TODO rename ?
   SimpleVector* currIcuppLink;

   // todo ? not nulled in setpointerstonull - probably intended not sure why though
   double* currEqRhsAdaptionsLink;
   double* currInEqRhsAdaptionsLink;
   double* currInEqLhsAdaptionsLink;

   SimpleVector* currgParent;
   SimpleVector* currgChild;

   SimpleVector* currRedRow;
   SimpleVector* currNnzRow;
   SimpleVector* currRedRowLink;
   SimpleVector* currRedColParent;
   SimpleVector* currRedColChild;
   SimpleVector* currNnzColParent;
   SimpleVector* currNnzColChild;

   /** the number of children */
   int nChildren;
   /** number of entry eliminations on this process in the current elimination routine */
   int localNelims;

   /** vector containing the removed entries */
   std::vector<MTRXENTRY> removedEntries;

   /** array of length nChildren+1 to store start and end indices for removedEntries
    * that correspond to the linking-variable block (usually Amat).
    * As linkVarsBlocks[0] represents the parent block, the child block 'it' is accessed
    * using the index 'it+1'. */
   BLOCKS* linkVarsBlocks;
   /** array of length nChildren+1 to store start and end indices for removedEntries
    * that correspond to the child block (usually Bmat).
    * As childBlocks[0] represents the parent block which has no 'free' block,
    * the child block 'it' is accessed using the index 'it+1'. */
   BLOCKS* childBlocks;

   std::vector<XBOUNDS> newBoundsParent;

   /* objective offset resulting from local presolving*/
   double indivObjOffset;

   /* swap two entries in the SparseStorageDynamic format */
   void updateAndSwap( SparseStorageDynamic* storage, int rowidx, int& indexK, int& rowEnd, double* redCol, int& nelims);
   void updateRhsNRowLink();
   void updateNnzUsingReductions( OoqpVector* nnzVector, OoqpVector* redVector) const;
   void updateNnzColParent(MPI_Comm comm);
   void allreduceAndUpdate(MPI_Comm comm, SimpleVector& adaptionsVector, SimpleVector& baseVector);

   void storeRemovedEntryIndex(int rowidx, int colidx, int it, BlockType block_type);
   // methods to update the transposed matrix:
   void updateTransposed(StochGenMatrix& matrix);
   void updateTransposedSubmatrix(SparseStorageDynamic& transStorage, const int blockStart, const int blockEnd) const;

   void updateLinkingVarsBlocks(int& newSREq, int& newSRIneq);
   bool newBoundsTightenOldBounds(double new_low, double new_upp, int index,
         double* ilow, double* iupp, double* low, double* upp) const;
   void setNewBounds(int index, double new_low, double new_upp,
         double* ilow, double* low, double* iupp, double* upp) const;
   void setNewBound(int index, double new_bound,
         SimpleVector* bound_vector, SimpleVector* i_bound_vector) const;

   // methods to update the current pointers:
   void setCurrentPointersToNull();
   /** initialize current pointer for matrices and vectors.
    * If it==-1, we are at parent and want block B_0 (Bmat).
    * Returns false if it is a dummy child. */
   bool updateCurrentPointersForColAdapt(int it, SystemType system_type);
   bool updateCPforColAdaptF0( SystemType system_type );
   bool updateCPforAdaptFixationsBChild( int it, SystemType system_type);

   void setCPAmatsRoot(GenMatrixHandle matrixHandle);
   bool setCPAmatsChild(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   bool setCPBmatsChild(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   bool setPointersForAmatBmat(int node, SystemType system_type);
   bool setCPLinkConstraint(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   void setCPBlmatsRoot(GenMatrixHandle matrixHandle);
   void setCPBlmatsChild(GenMatrixHandle matrixHandle, int it);
   void setCPColumnRoot();
   void setCPColumnChild(int it);

   void setPointersRowRoot(SystemType system_type);
   void setCPRowRootEquality();
   void setCPRowRootInequality();

   void setCPRowRootIneqOnlyLhsRhs();

   void setPointersRowChildNode(SystemType system_type, int node);
   void setCPRowChildEquality(int it);
   void setCPRowChildInequality(int it);

   void setCPRowChildIneqOnlyLhsRhs(int it);
   void setCPRowLinkEquality();
   void setCPRhsLinkInequality();

   void resetLinkvarsAndChildBlocks();
   void resetEqRhsAdaptionsLink();
   void resetIneqRhsAdaptionsLink();

   bool removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx, double& m) const;
   void clearRow(SparseStorageDynamic& storage, const int rowIdx) const;
   void removeRow(int rowIdx, SparseStorageDynamic& Ablock, SparseStorageDynamic& AblockTrans,
         SparseStorageDynamic* Bblock, SparseStorageDynamic* BblockTrans, SimpleVector& nnzRow,
         SimpleVector& redColParent, SimpleVector* nnzColChild);
   void removeRowInBblock(int rowIdx, SparseStorageDynamic* Bblock,
         SparseStorageDynamic* BblockTrans, SimpleVector* nnzColChild);

   bool nodeIsDummy(int it, SystemType system_type);
   bool hasLinking(SystemType system_type) const;
   void getRankDistributed(MPI_Comm comm, int& myRank, bool& iAmDistrib) const;
   void abortInfeasible(MPI_Comm comm) const;
   void synchronize(int& value) const;
   void synchronizeSum(int& first, int& second) const;


   void adaptChildBmat( std::vector<COLUMNTOADAPT> const & colAdaptBlock, SystemType system_type, int& newSR);
   void adaptChildBlmat( std::vector<COLUMNTOADAPT> const & colAdaptBlock, SystemType system_type);
   int adaptChildBmatCol(int colIdx, double val, SystemType system_type);
   void adaptOtherSystemChildB(SystemType system_type, std::vector<COLUMNTOADAPT> const & colAdaptBblock, int& newSR);

   int colAdaptLinkVars(int it, SystemType system_type);
   int colAdaptF0(SystemType system_type);

   bool newBoundsFixVariable(double& value, double newxlow, double newxupp, int colIdx,
      double* ixlow, double* ixupp, double* xlow, double* xupp) const;
   int fixVarInChildBlockAndStore(int colIdx, double val, SystemType system_type,
         std::vector<COLUMNTOADAPT> & colAdaptLinkBlock);
   void storeColValInColAdaptParent(int colIdx, double value);
   bool newBoundsImplyInfeasible(double newxlow, double newxupp, int colIdx,
      double* ixlow, double* ixupp, double* xlow, double* xupp) const;
   void storeNewBoundsParent(int colIdx, double newxlow, double newxupp);
   void combineNewBoundsParent();
   XBOUNDS getNewBoundsParent(int i) const;
   void setNewBoundsParent(int i, int colIdx, double newxlow, double newxupp);
   int getNumberNewBoundsParent() const;
   void addNewBoundsParent(XBOUNDS newXBounds);
   void clearNewBoundsParent();

   void sumIndivObjOffset();
   void computeActivityBlockwise( SparseStorageDynamic& matrix, int rowIdx, int colIdx,
         double& infRow, double& supRow,
         SimpleVector& xlow, SimpleVector& ixlow, SimpleVector& xupp, SimpleVector& ixupp);

   void countRangedRowsBlock(int& nRangedRows, int& nRowsIneq) const;
   void countEqualityRowsBlock(int& nRowsEq) const;
   void countSingletonRowsBlock(int& nSingletonRows) const;
   void countBoxedColumns(int& nBoxCols, int& nColsTotal, int& nFreeVars) const;
};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
