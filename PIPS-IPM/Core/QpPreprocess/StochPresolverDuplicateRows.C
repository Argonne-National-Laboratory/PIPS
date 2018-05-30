/*
 * StochPresolverDuplicateRows.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverDuplicateRows.h"

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

StochPresolverDuplicateRows::StochPresolverDuplicateRows(PresolveData& presData)
: StochPresolverBase(presData)
{
   norm_Amat = NULL;
   norm_Bmat= NULL;
   norm_Cmat = NULL;
   norm_Dmat = NULL;
   norm_b = NULL;
   norm_c = NULL;
   norm_d = NULL;

   mA = 0;
   nA = 0;
}

StochPresolverDuplicateRows::~StochPresolverDuplicateRows()
{
 // todo
}

bool StochPresolverDuplicateRows::applyPresolving(int& nelims)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

  /* boost::unordered_set<rowlib::row, boost::hash<rowlib::row> > rows;
   int col0[] = {1,2,3};
   rowlib::row row0(0, 3, col0);
   int col1[] = {1,2,3};
   rowlib::row row1(1, 3, col1);
   int col2[] = {1,2,4};
   rowlib::row row2(2, 3, col2);

   rows.insert(row0);
   rows.insert(row1);
   rows.insert(row2);
   int col3[] = {0,1,2,3};
   rows.emplace(3,4,col3);

   unsigned n = rows.bucket_count();
   std::cout << "unordered_set rows has size " << rows.size() << '\n';
   std::cout << "unordered_set rows has " << n << " buckets.\n";
   for (unsigned i=0; i<n; ++i)
   {
       std::cout << "bucket #" << i << " contains: ";
       for (boost::unordered_set<rowlib::row>::local_iterator it = rows.begin(i); it!=rows.end(i); ++it)
           std::cout << " row id#"<<it->id;
       std::cout << "\n";
   }
   assert( rows.find(row0) != rows.end());*/

  /*if( myRank == 0 ) cout<<"Before duplicate Row Presolving:"<<endl;
   countRowsCols();*/

   StochGenMatrix& matrixA = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& matrixC = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   // for children:
   for( size_t it = 0; it< matrixA.children.size(); it++)
   {
      // copy and normalize A,B,C,D and b,c,d
      setNormalizedPointers((int)it, matrixA, matrixC );

      // Initialize unordered set 'rows'
      // Per row, add row to 'rows'
      // Second Hashing: Per bucket, do ...
      // When two parallel rows are found, check if they are both =, both <=, or = and and <=
      // The action has to be applied to the original matrices (not the normalized copies)
   }

   SparseStorageDynamic norm_storage = matrixA.Bmat->getStorageDynamicRef();

   assert(0);


   //countDuplicateRows(matrixC, INEQUALITY_SYSTEM);
   //countDuplicateRows(matrixA, EQUALITY_SYSTEM);

   return true;
}

bool StochPresolverDuplicateRows::setNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC)
{
   // check if it is no dummy child
   // copy the matrices
   // normalize the rows

   if( !childIsDummy(matrixA, it, EQUALITY_SYSTEM) )
   {
      // todo: copy instead of setting pointers to original data
      norm_Amat = dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Amat)->getStorageDynamic();
      norm_Bmat = dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Bmat)->getStorageDynamic();
      norm_b = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
   }

   return true;
}

void StochPresolverDuplicateRows::countDuplicateRows(StochGenMatrix& matrix, SystemType system_type)
{
   //TODO: several duplicate rows are counted too often because each duplicate pair is counted as one
   // ie, if rows i,j,k,l are duplicate, they are counted (i,j) (i,k) (i,l) (j,k) (j,l)(k,l) as 6 instead of 4
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   int duplicRow = 0;

   // the B_0 block:
   if( myRank == 0)
   {
      currAmat = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamic();
      if(system_type==EQUALITY_SYSTEM)
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
      else
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec);
      for( int i=0; i<currNnzRow->n; i++)
      {
         for( int j=i+1; j<currNnzRow->n; j++)
         {
            double* nRow = currNnzRow->elements();
            if( nRow[i] != 0.0 && nRow[i] == nRow[j] &&
                  compareCoefficients(*currAmat, i, j))
            {
               duplicRow++;
               //cout<<"Row "<<i<<" and "<<j<<" of B_0 are duplicate rows."<<endl;
            }
         }
      }
      if(system_type==EQUALITY_SYSTEM)
         cout<<"There are "<<duplicRow<<" duplicate rows in the A_0 block of A."<<endl;
      else
         cout<<"There are "<<duplicRow<<" duplicate rows in the A_0 block of C."<<endl;
   }

   // the children:
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( (system_type==EQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, EQUALITY_SYSTEM)) ||
            (system_type == INEQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, INEQUALITY_SYSTEM)) )
      {
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
         if(system_type==EQUALITY_SYSTEM)
            currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
         else
            currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
         for( int i=0; i<currNnzRow->n; i++)
         {
            for( int j=i+1; j<currNnzRow->n; j++)
            {
               double* nRow = currNnzRow->elements();
               if( nRow[i] != 0.0 && nRow[i] == nRow[j] &&
                     compareCoefficients(*currAmat, i, j) && compareCoefficients(*currBmat, i, j) )
               {
                  duplicRow++;
                  //cout<<"Row "<<i<<" and "<<j<<" of child "<<it<<" are duplicate rows."<<endl;
               }
            }
         }
      }
   }

   // the linking constraints:
   //TODO: more efficient communication, maybe vector containing info about all rows?
   if( (system_type==EQUALITY_SYSTEM && hasLinking(EQUALITY_SYSTEM)) ||
         (system_type==INEQUALITY_SYSTEM && hasLinking(INEQUALITY_SYSTEM)))
   {
      int nDuplicLinkRow = 0;
      if(system_type==EQUALITY_SYSTEM)
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      else
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      for(int i=0; i<currNnzRow->n; i++)
      {
         for( int j=i+1; j<currNnzRow->n; j++)
         {
            int notParallel = 0;
            double* nRow = currNnzRow->elements();
            currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamic();
            if( nRow[i] != 0.0 && nRow[i] == nRow[j] &&
                  compareCoefficients(*currBlmat, i, j) )
            {
               // check other F_i blocks as well:
               for( size_t it = 0; it< matrix.children.size(); it++)
               {
                  if( (system_type==EQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, EQUALITY_SYSTEM)) ||
                        (system_type==INEQUALITY_SYSTEM && !childIsDummy( matrix, (int)it, INEQUALITY_SYSTEM)) )
                  {
                     currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamic();
                     if( !compareCoefficients(*currBlmat, i, j))
                        notParallel++;
                  }
                  if(system_type==EQUALITY_SYSTEM)
                     currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
                  else
                     currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
               }
               if( iAmDistrib )
                  MPI_Allreduce(MPI_IN_PLACE, &notParallel, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               if( myRank == 0 && notParallel == 0 )
               {
                  duplicRow++;
                  nDuplicLinkRow++;
                  //cout<<"Row "<<i<<" and "<<j<<" of Linking Constraints are duplicate rows."<<endl;
               }
            }
         }
      }
      if(myRank == 0)
      {
         assert(nDuplicLinkRow <= currNnzRow->n);
         //duplicRow += nDuplicLinkRow;
         if(system_type==EQUALITY_SYSTEM)
            cout<<"There are "<<nDuplicLinkRow<<" duplicate rows in the linking rows of A."<<endl;
         else
            cout<<"There are "<<nDuplicLinkRow<<" duplicate rows in the linking rows of C."<<endl;
      }
   }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &duplicRow, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if(myRank == 0)
   {
      if(system_type==EQUALITY_SYSTEM)
         cout<<"There are "<<duplicRow<<" duplicate rows in A."<<endl;
      else
         cout<<"There are "<<duplicRow<<" duplicate rows in C."<<endl;
   }
}

/* Compares the matrix coeffients in row i and j.
 * Returns true if they are the same or false if not.
 * The method assumes ordered rows, so it should only be called when it is
 * ensured that the order in each row is intact.
 * Empty rows can be considered as duplicate rows.
 */
bool StochPresolverDuplicateRows::compareCoefficients(SparseStorageDynamic& matrix, int i, int j) const
{
   assert( currNnzRow->elements()[i] == currNnzRow->elements()[j] );
   assert( i>=0 && i<matrix.m );
   assert( j>=0 && j<matrix.m );

   int rowLen = matrix.rowptr[i].end - matrix.rowptr[i].start;

   if( rowLen != matrix.rowptr[j].end - matrix.rowptr[j].start )
      return false;  // row entries can be in Amat or Bmat blocks
   if( rowLen == 0.0 )
      return true;  // empty rows

   bool allNegatedCoeffs = false;
   // compare first entry to determine the sign:
   int k_i = matrix.rowptr[i].start;
   int k_j = matrix.rowptr[j].start;
   if( matrix.jcolM[k_i] != matrix.jcolM[k_j] )
      return false;
   if( matrix.M[k_i] == - matrix.M[k_j] )
      allNegatedCoeffs = true;
   else if ( matrix.M[k_i] != matrix.M[k_j] )
      return false;

   for( int k=1; k<rowLen; k++ )
   {
      k_i = matrix.rowptr[i].start + k;
      k_j = matrix.rowptr[j].start + k;

      if( matrix.jcolM[k_i] != matrix.jcolM[k_j] )
         return false;
      if( (allNegatedCoeffs && matrix.M[k_i] != -matrix.M[k_j]) ||
            (!allNegatedCoeffs && matrix.M[k_i] != matrix.M[k_j]) )
         return false;
   }
   return true;
}
