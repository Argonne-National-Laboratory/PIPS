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
  static const int maxIterSR = 20;

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

  int doSingletonRows();
  int initSingletonRows(SystemType system_type);
  int initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple);
  bool doSingletonRowsA(int& newSREq, int& newSRIneq);

  bool procSingletonRowRoot(StochGenMatrix& stochMatrix);
  bool procSingletonRowChild(StochGenMatrix& stochMatrix, int it, int& newSR, int& newSRIneq);
  bool procSingletonRowChildAmat(SparseStorageDynamic& A_mat, int it);
  bool procSingletonRowChildBmat(SparseStorageDynamic& B_mat, int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR);
  bool removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR);
  int adaptChildBmatCol(int colIdx, double val, SystemType system_type);
  bool adaptInequalityChildB(std::vector<COLUMNTOADAPT> & colAdaptBblock, int& newSRIneq);
  bool removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx);

  int colAdaptLinkVars(int it, SystemType system_type);

  int colAdaptF0(SystemType system_type);
  bool combineColAdaptParent();
  int doSingletonRowsC();

  void updateObjOffset();

public:

  StochPresolver(const Data* prob);
  virtual ~StochPresolver();

  virtual Data* presolve();
};

//@}




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_ */
