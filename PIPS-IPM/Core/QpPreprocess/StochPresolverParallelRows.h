/*
 * StochPresolverParallelRows.h
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERPARALLELROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERPARALLELROWS_H_

#include "StochPresolverBase.h"

#include <boost/unordered_set.hpp>

namespace rowlib
{
   static const double offset_hash_double = 0.127;

    struct rowWithColInd
    {
        int id;
        int offset_nA;
        int lengthA;
        int* colIndicesA;
        double* norm_entriesA;
        int lengthB;
        int* colIndicesB;
        double* norm_entriesB;

        rowWithColInd(int id, int offset, int lenA, int* colA, double* entA, int lenB, int* colB, double* entB)
            : id(id), offset_nA(offset), lengthA(lenA), colIndicesA(colA), norm_entriesA(entA),
              lengthB(lenB), colIndicesB(colB), norm_entriesB(entB)  {}
    };

    bool operator==(rowWithColInd const& a, rowWithColInd const& b);
    std::size_t hash_value(rowWithColInd const& b);

    struct rowWithEntries
    {
       int id;
       int offset_nA;
       int lengthA;
       int* colIndicesA;
       double* norm_entriesA;
       int lengthB;
       int* colIndicesB;
       double* norm_entriesB;

        rowWithEntries(int id, int offset, int lenA, int* colA, double* entA, int lenB, int* colB, double* entB)
            : id(id), offset_nA(offset), lengthA(lenA), colIndicesA(colA), norm_entriesA(entA),
              lengthB(lenB), colIndicesB(colB), norm_entriesB(entB) {}
    };

    bool operator==(rowWithEntries const& a, rowWithEntries const& b);
    std::size_t hash_value(rowWithEntries const& b);
}

class StochPresolverParallelRows : public StochPresolverBase
{
public:
   StochPresolverParallelRows(PresolveData& presData);

   ~StochPresolverParallelRows();

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

private:
   SparseStorageDynamic* currCmat;
   SparseStorageDynamic* currCmatTrans;
   SparseStorageDynamic* currDmat;
   SparseStorageDynamic* currDmatTrans;
   SimpleVector* currNnzRowC;

   // pointers to the normalized and copied matrix blocks
   SparseStorageDynamic* norm_Amat;
   SparseStorageDynamic* norm_Bmat;
   SparseStorageDynamic* norm_Cmat;
   SparseStorageDynamic* norm_Dmat;
   SimpleVector* norm_b;
   SimpleVector* norm_clow;
   SimpleVector* norm_cupp;
   SimpleVector* norm_iclow;
   SimpleVector* norm_icupp;
   SimpleVector* norm_factorC;
   SimpleVector* norm_factorA;

   // data for the nearly parallel row case
   SimpleVector* rowContainsSingletonVariableA;
   SimpleVector* rowContainsSingletonVariableC;
   SimpleVector* singletonCoeffsColParent;
   SimpleVector* singletonCoeffsColChild;
   SimpleVector* normNnzRowA;
   SimpleVector* normNnzRowC;
   SimpleVector* normNnzColParent;
   SimpleVector* normNnzColChild;
   SparseStorageDynamic* norm_AmatTrans;
   SparseStorageDynamic* norm_BmatTrans;
   SparseStorageDynamic* norm_CmatTrans;
   SparseStorageDynamic* norm_DmatTrans;

   // number of rows of the A or B block
   int mA;
   // number of columns of the A or C block
   int nA;

   boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> > rowsFirstHashTable;
   boost::unordered_set<rowlib::rowWithEntries, boost::hash<rowlib::rowWithEntries> > rowsSecondHashTable;

   bool setNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC);
   void deleteNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC);
   void removeSingletonVars();
   void removeEntry(int colIdx, SimpleVector& rowContainsSingletonVar,
         SparseStorageDynamic& matrix, SparseStorageDynamic& matrixTrans, SimpleVector& nnzRow, SimpleVector& nnzCol,
         BlockType block_type);
   void normalizeBlocksRowwise( SystemType system_type, SparseStorageDynamic& Ablock, SparseStorageDynamic* Bblock,
         SimpleVector& Rhs, SimpleVector* Lhs, SimpleVector* iRhs, SimpleVector* iLhs);
   void insertRowsIntoHashtable( boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> > &rows,
         SparseStorageDynamic& Ablock, SparseStorageDynamic* Bblock, SystemType system_type, SimpleVector& nnzRow );
   bool compareRowsInSecondHashTable(int& nRowElims, int it);
   bool checkRowsAreParallel( rowlib::rowWithEntries row1, rowlib::rowWithEntries row2);
   void eliminateOriginalRow(int rowId, int& nRowElims);
   void removeRow(int rowIdx, SparseStorageDynamic& Ablock, SparseStorageDynamic& AblockTrans,
         SparseStorageDynamic* Bblock, SparseStorageDynamic* BblockTrans, SimpleVector& nnzRow,
         SimpleVector& redColParent, SimpleVector* nnzColChild);
   bool tightenOriginalBoundsOfRow1(int rowId1, int rowId2);
   double getSingletonCoefficient(int singleColIdx);
   void tightenBoundsForSingleVar(int singleColIdx, double newxlow, double newxupp);
   bool doNearlyParallelRowCase3(int rowId1, int rowId2);
   double addCostToVariable1(int singleColIdx1, int singleColIdx2, double t);

   void countDuplicateRows(StochGenMatrix& matrix, SystemType system_type);
   bool compareCoefficients(SparseStorageDynamic& matrix, int i, int j) const;

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERPARALLELROWS_H_ */
