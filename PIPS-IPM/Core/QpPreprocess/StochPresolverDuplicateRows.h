/*
 * StochPresolverDuplicateRows.h
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUPLICATEROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUPLICATEROWS_H_

#include "StochPresolverBase.h"

#include <boost/unordered_set.hpp>

namespace rowlib
{
    struct row
    {
        int id;
        int length;
        int* colIndices;

        row(int i, int n, int* t)
            : id(i), length(n), colIndices(t) {}
    };

    bool operator==(row const& a, row const& b)
    {
        return a.id == b.id;
    }

    std::size_t hash_value(row const& b)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, b.length);
        for(int i=0; i<b.length; i++)
            boost::hash_combine(seed, b.colIndices[i]);
        return seed;
    }
}

class StochPresolverDuplicateRows : public StochPresolverBase
{
public:
   StochPresolverDuplicateRows(PresolveData& presData);

   ~StochPresolverDuplicateRows();

   // remove small matrix entries and return number of eliminations
   virtual bool applyPresolving(int& nelims);

private:
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

   // number of rows of the A or C block
   int mA;

   // number of columns of the A or C block
   int nA;

   bool setNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC);
   void deleteNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC);
   void normalizeBLocksRowwise( SystemType system_type, SparseStorageDynamic* Ablock, SparseStorageDynamic* Bblock,
         SimpleVector* Rhs, SimpleVector* Lhs, SimpleVector* iRhs, SimpleVector* iLhs);
   void insertRowsIntoHashtable( boost::unordered_set<rowlib::row, boost::hash<rowlib::row> > rows, int it,
         StochGenMatrix& matrixA, StochGenMatrix& matrixC );

   void countDuplicateRows(StochGenMatrix& matrix, SystemType system_type);
   bool compareCoefficients(SparseStorageDynamic& matrix, int i, int j) const;

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERDUPLICATEROWS_H_ */
