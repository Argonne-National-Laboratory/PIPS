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
#include "DoubleMatrixTypes.h"
#include "SmartPointer.h"
#include "PresolveData.h"
#include "sData.h"
#include <vector>
#include <cassert>

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

enum SystemType {EQUALITY_SYSTEM, INEQUALITY_SYSTEM};
enum BlockType {LINKING_VARS_BLOCK, CHILD_BLOCK};


class StochPresolverBase
{
public:
   StochPresolverBase(PresolveData& presData);

   virtual ~StochPresolverBase();

   virtual bool applyPresolving(int& nelims) = 0;


protected:
   static const double feastol = 1.0e-6;
   static const double infinity = 10e30;
   static const double tolerance1 = 1.0e-3;
   static const double tolerance2 = 1.0e-2;
   static const double tolerance3 = 1.0e-10;
   static const int maxIterSR = 20;


   // pointers to the currently needed matrices and vectors for presolving
   sData* presProb;
   PresolveData& presData;

   SparseStorageDynamic* currAmat;
   SparseStorageDynamic* currAmatTrans;
   SparseStorageDynamic* currBmat;
   SparseStorageDynamic* currBmatTrans;
   SparseStorageDynamic* currBlmat;
   SparseStorageDynamic* currBlmatTrans;
   SimpleVector* currxlowParent;
   SimpleVector* currxlowChild;
   SimpleVector* currxuppParent;
   SimpleVector* currxuppChild;
   SimpleVector* currIxlowParent;
   SimpleVector* currIxlowChild;
   SimpleVector* currIxuppParent;
   SimpleVector* currIxuppChild;
   SimpleVector* currEqRhs;
   SimpleVector* currIneqRhs;
   SimpleVector* currIneqLhs;
   SimpleVector* currIcupp;
   SimpleVector* currIclow;
   SimpleVector* currEqRhsLink;
   SimpleVector* currIneqRhsLink;
   SimpleVector* currIneqLhsLink;
   SimpleVector* currIcuppLink;
   SimpleVector* currIclowLink;
   double* currEqRhsAdaptionsLink;
   double* currInEqRhsAdaptionsLink;
   double* currInEqLhsAdaptionsLink;

   SimpleVector* currgParent;
   SimpleVector* currgChild;

   SimpleVector* currNnzRow;
   SimpleVector* currRedRow;
   SimpleVector* currRedRowLink;
   SimpleVector* currRedColParent;
   SimpleVector* currNnzColParent;
   SimpleVector* currRedColChild;
   SimpleVector* currNnzColChild;

   /** the number of children */
   int nChildren;
   /** number of eliminations on this process in the current elimination routine */
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

   /* swap two entries in the SparseStorageDynamic format */
   void updateAndSwap( SparseStorageDynamic* storage, int rowidx, int& indexK, int& rowEnd, double* redCol, int& nelims);
   void updateRhsNRowLink();
   void updateNnzUsingReductions( OoqpVector* nnzVector, OoqpVector* redVector);

   void storeRemovedEntryIndex(int rowidx, int colidx, int it, BlockType block_type);
   // methods to update the transposed matrix:
   void updateTransposed(StochGenMatrix& matrix);
   void updateTransposedSubmatrix(SparseStorageDynamic& transStorage, const int blockStart, const int blockEnd);

   void updateLinkingVarsBlocks(int& newSREq, int& newSRIneq);

   // methods to update the current pointers:
   void setCurrentPointersToNull();
   /** initialize current pointer for matrices and vectors.
    * If it==-1, we are at parent and want block B_0 (Bmat).
    * Returns false if it is a dummy child. */
   bool updateCPforTinyEntry(int it, SystemType system_type);
   bool updateCPForSingletonRow(int it, SystemType system_type);
   bool updateCPForSingletonRowInequalityBChild( int it );
   bool updateCurrentPointersForColAdapt(int it, SystemType system_type);
   bool updateCPforColAdaptF0( SystemType system_type );

   void setCPAmatsRoot(GenMatrixHandle matrixHandle);
   bool setCPAmatsChild(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   bool setCPBmatsChild(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   bool setCPAmatBmat(GenMatrixHandle matrixHandle, int it, SystemType system_type);
   void setCPBlmatsRoot(GenMatrixHandle matrixHandle);
   void setCPBlmatsChild(GenMatrixHandle matrixHandle, int it);
   void setCPColumnRoot();
   void setCPColumnChild(int it);
   void setCPRowRootEquality();
   void setCPRowRootInequality();
   void setCPRowChildEquality(int it);
   void setCPRowChildInequality(int it);
   void setCPRowLinkEquality();
   void setCPRhsLinkInequality();

   void resetLinkvarsAndChildBlocks();
   void resetEqRhsAdaptionsLink();
   void resetIneqRhsAdaptionsLink();

   double removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx);
   void clearRow(SparseStorageDynamic& storage, const int rowIdx);

   bool childIsDummy(StochGenMatrix& matrix, int it, SystemType system_type);
   bool hasLinking(SystemType system_type);
   void getRankDistributed(MPI_Comm comm, int& myRank, bool& iAmDistrib);

   bool adaptChildBmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type, int& newSR);
   bool adaptChildBlmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type);
   int adaptChildBmatCol(int colIdx, double val, SystemType system_type);
   bool adaptInequalityChildB(std::vector<COLUMNTOADAPT> & colAdaptBblock, int& newSRIneq);

   int colAdaptLinkVars(int it, SystemType system_type);
   int colAdaptF0(SystemType system_type);

};



#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
