/*
 * StochPresolverDuplicateRows.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverDuplicateRows.h"

StochPresolverDuplicateRows::StochPresolverDuplicateRows(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
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

   if( myRank == 0 )
      cout<<"Before duplicate Row Presolving:"<<endl;
   countRowsCols();

   int duplicRow = 0;
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));
   // for now, only look in C
   // the B_0 block:
   if( myRank == 0)
   {
      currAmat = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamic();
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
               cout<<"Row "<<i<<" and "<<j<<" of B_0 are duplicate rows."<<endl;
            }
         }
      }
   }

   // the children:
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( !childIsDummy( matrix, (int)it, INEQUALITY_SYSTEM) )
      {
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
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
                  cout<<"Row "<<i<<" and "<<j<<" of child "<<it<<" are duplicate rows."<<endl;
               }
            }
         }
      }
   }

   // the linking constraints:
   if( hasLinking( INEQUALITY_SYSTEM ))
   {
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      for(int i=0; i<currNnzRow->n; i++)
      {
         for( int j=i+1; j<currNnzRow->n; j++)
         {
            int counter = 0;
            double* nRow = currNnzRow->elements();
            currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamic();
            if( nRow[i] != 0.0 && nRow[i] == nRow[j] &&
                  compareCoefficients(*currBlmat, i, j) )
            {
               // check other F_i blocks as well:
               for( size_t it = 0; it< matrix.children.size(); it++)
               {
                  if( !childIsDummy( matrix, (int)it, INEQUALITY_SYSTEM) )
                  {
                     currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamic();
                     if( !compareCoefficients(*currBlmat, i, j))
                        counter++;
                  }
                  currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
               }
               if( iAmDistrib )
                  MPI_Allreduce(MPI_IN_PLACE, &counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
               if( myRank == 0 && counter == 0 )
               {
                  duplicRow++;
                  cout<<"Row "<<i<<" and "<<j<<" of Linking Constraints are duplicate rows."<<endl;
               }
            }
         }
      }
   }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &duplicRow, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if(myRank == 0) cout<<"There are "<<duplicRow<<" duplicate rows in C."<<endl;

   return true;
}

/* Compares the matrix coeffients in row i and j.
 * Returns true if they are the same or false if not.
 * The method assumes ordered rows, so it should only be called when it is
 * ensured that the order in each row is intact.
 * Empty rows can be considered as duplicate rows.
 */
bool StochPresolverDuplicateRows::compareCoefficients(SparseStorageDynamic& matrix, int i, int j)
{
   assert( currNnzRow->elements()[i] == currNnzRow->elements()[j] );

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
