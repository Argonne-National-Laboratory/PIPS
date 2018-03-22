/*
 * StochPresolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_

#include "QpPresolver.h"
#include "StochVector.h"
#include "StochGenMatrix.h"
#include "SmartPointer.h"
#include "sData.h"
#include <vector>

class Data;

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
   double val;
} COLUMNTOADAPT;

struct col_is_smaller
{
    bool operator()(const COLUMNTOADAPT& x, const COLUMNTOADAPT& y) const
    {
        return x.colIdx < y.colIdx;
    }
};

enum SystemType {EQUALITY_SYSTEM, INEQUALITY_SYSTEM};
enum BlockType {LINKING_VARS_BLOCK, CHILD_BLOCK};

/**  * @defgroup QpPreprocess
 *
 * QP presolver
 * @{
 */

/**
 * Derived class for QP presolvers.
 */
class StochPresolver : public QpPresolver
{
private:
  static const double feastol = 1.0e-6;
  static const double infinity = 10e30;
  static const double tolerance1 = 1.0e-3;
  static const double tolerance2 = 1.0e-2;
  static const double tolerance3 = 1.0e-10;
  static const int maxIterSR = 20;

  // number of non-zero elements of each row
  StochVectorHandle nRowElemsA;
  StochVectorHandle nRowElemsC;

  // number of non-zero elements of each column
  StochVectorHandle nColElems;

  // number of removed elements of each row / column
  StochVectorHandle redRowA;
  StochVectorHandle redRowC;
  StochVectorHandle redCol;

  // pointers to the currently needed matrices and vectors for presolving
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

  // variables used for singleton row elimination:
  /** vector containing the row indices of singleton rows */
  std::vector<int> singletonRows;
  /** array of length nChildren+3 to store start indices for singletonRows
   * that correspond to the correct block. As blocks[0] represents the parent block,
   * the child block 'it' is accessed using the index 'it+1'.
   * The linking-row block is accessed using the index nChildren+2. */
  int* blocks;
  std::vector<int> singletonRowsIneq;
  int* blocksIneq;

  // variables used for singleton row elimination:
  /** vector containing the row indices of singleton rows */
  std::vector<int> singletonRows2;
  /** array of length nChildren+3 to store start indices for singletonRows
   * that correspond to the correct block. As blocks[0] represents the parent block,
   * the child block 'it' is accessed using the index 'it+1'.
   * The linking-row block is accessed using the index nChildren+2. */
  int* blocks2;
  std::vector<int> singletonRowsIneq2;
  int* blocksIneq2;

  /** vector containing the column indices of entries that were found during the
   * singleton row routine. Along with the column index, the value needed for
   * adaptation is stored. */
  std::vector<COLUMNTOADAPT> colAdaptParent;

  /** objective offset created by presolving*/
  double objOffset;

  // initialize row and column nnz counter
  void initNnzCounter();

  void setCurrentPointersToNull();

  /** initialize current pointer for matrices and vectors.
   * If it==-1, we are at parent and want block B_0 (Bmat).
   * Returns false if it is a dummy child. */
  bool updateCurrentPointers(int it, SystemType system_type);

  // remove small matrix entries and return number of eliminations
  int removeTinyEntries();

  int removeTinyEntriesSystemA();
  int removeTinyEntriesSystemC();

  int removeTinyChild( int it, SystemType system_type );

  int removeTinyInnerLoop( int it, SystemType system_type, BlockType block_type );
  void updateAndSwap( SparseStorageDynamic* storage, int rowidx, int& indexK, int& rowEnd, double* redCol, int& nelims);
  void storeRemovedEntryIndex(int rowidx, int colidx, int it, BlockType block_type);

  void updateNnzUsingReductions( OoqpVector* nnzVector, OoqpVector* redVector);

  void updateTransposed(StochGenMatrix& matrix);
  void updateTransposedSubmatrix(SparseStorageDynamic& transStorage, const int blockStart, const int blockEnd);

  int doSingletonRows();
  int initSingletonRows(SystemType system_type);
  int initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple);
  bool doSingletonRowsA(int& newSREq, int& newSRIneq);
  void updateLinkingVarsBlocks(int& newSREq, int& newSRIneq);
  bool updateCurrentPointersForSingletonRow(int it, SystemType system_type);
  bool updateCPForSingletonRowInequalityBChild( int it );
  bool procSingletonRowRoot(StochGenMatrix& stochMatrix);
  bool procSingletonRowChild(StochGenMatrix& stochMatrix, int it, int& newSR, int& newSRIneq);
  bool procSingletonRowChildAmat(SparseStorageDynamic& A_mat, int it);
  bool procSingletonRowChildBmat(SparseStorageDynamic& B_mat, int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR);
  bool removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR);
  int adaptChildBmatCol(int colIdx, double val, SystemType system_type);
  bool adaptInequalityChildB(std::vector<COLUMNTOADAPT> & colAdaptBblock, int& newSRIneq);
  bool removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx);
  bool adaptChildBmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type, int& newSR);
  bool adaptChildBlmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type);

  int colAdaptLinkVars(int it, SystemType system_type);
  bool updateCurrentPointersForColAdapt(int it, SystemType system_type);
  bool updateCPforColAdaptF0( SystemType system_type );
  int colAdaptF0(SystemType system_type);
  bool combineColAdaptParent();
  void updateRhsNRowLink();
  int doSingletonRowsC();

  void resetLinkvarsAndChildBlocks();
  void resetBlocks();
  void resetRedCounters();
  void resetEqRhsAdaptionsLink();
  void resetIneqRhsAdaptionsLink();

  double removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx);
  void clearRow(SparseStorageDynamic& storage, const int rowIdx);

  bool childIsDummy(StochGenMatrix& matrix, int it, SystemType system_type);
  bool hasLinking(SystemType system_type);
  void getRankDistributed( MPI_Comm comm, int& myRank, bool& iAmDistrib );
  void updateObjOffset();

  sData* presProb;

public:

  StochPresolver(const Data* prob);
  virtual ~StochPresolver();

  virtual Data* presolve();
};

//@}




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_ */
