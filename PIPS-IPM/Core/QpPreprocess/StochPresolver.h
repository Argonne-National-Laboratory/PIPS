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
  SparseStorageDynamic* currBmat;
  SparseStorageDynamic* currBlmat;
  SimpleVector* currxlowParent;
  SimpleVector* currxlowChild;
  SimpleVector* currxuppParent;
  SimpleVector* currxuppChild;
  SimpleVector* currEqRhs;
  SimpleVector* currIneqRhs;
  SimpleVector* currIneqLhs;
  SimpleVector* currIcupp;
  SimpleVector* currIclow;

  SimpleVector* currNnzRow;
  SimpleVector* currRedRow;
  SimpleVector* currRedColParent;
  SimpleVector* currRedColChild;

  int nChildren;
  int localNelims;
  std::vector<MTRXENTRY> removedEntries;
  BLOCKS* linkVarsBlocks;
  BLOCKS* childBlocks;

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
  void resetBlocks();

  sData* presProb;

public:

  StochPresolver(const Data* prob);
  virtual ~StochPresolver();

  virtual Data* presolve();
};

//@}




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_ */
