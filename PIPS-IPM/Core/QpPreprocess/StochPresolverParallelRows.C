/*
 * StochPresolverParallelRows.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

//#define PIPS_DEBUG
#include "StochPresolverParallelRows.h"

namespace rowlib
{
   bool operator==(rowWithColInd const& row1, rowWithColInd const& row2)
   {
      return row1.id == row2.id;
   }

   std::size_t hash_value(rowWithColInd const& row)
   {
      std::size_t seed = 0;
      boost::hash_combine(seed, (row.lengthA+row.lengthB) );
      for( int i = 0; i < row.lengthA; i++ )
         boost::hash_combine(seed, row.colIndicesA[i]);
      for( int i = 0; i < row.lengthB; i++ )
         boost::hash_combine(seed, (row.colIndicesB[i] + row.offset_nA) );
      return seed;
   }

   bool operator==(rowWithEntries const& row1, rowWithEntries const& row2)
   {
      return row1.id == row2.id;
   }

   std::size_t hash_value(rowWithEntries const& row)
   {
      std::size_t seed = 0;
      for( int i = 0; i < row.lengthA; i++ )
      {
         // Instead of hashing the normalized double coefficient, use an integer representation:
         int value_to_hash;
         double mantisse = std::frexp( (row.norm_entriesA[i] + offset_hash_double), &value_to_hash);
         value_to_hash += 10*( (int)trunc(mantisse*10000) );
         boost::hash_combine(seed, value_to_hash);
      }
      for( int i = 0; i < row.lengthB; i++ )
      {
         int value_to_hash;
         double mantisse = std::frexp( (row.norm_entriesB[i] + offset_hash_double), &value_to_hash);
         value_to_hash += 10*( (int)trunc(mantisse*10000) );
         boost::hash_combine(seed, value_to_hash);
      }
      return seed;
   }
}

StochPresolverParallelRows::StochPresolverParallelRows(PresolveData& presData)
: StochPresolverBase(presData)
{
   currCmat = NULL;
   currCmatTrans = NULL;
   currDmat = NULL;
   currDmatTrans = NULL;
   currNnzRowC = NULL;
   norm_Amat = NULL;
   norm_Bmat= NULL;
   norm_Cmat = NULL;
   norm_Dmat = NULL;
   norm_b = NULL;
   norm_cupp = NULL;
   norm_clow = NULL;
   norm_icupp = NULL;
   norm_iclow = NULL;
   norm_factorA = NULL;
   rowContainsSingletonVariableA = NULL;
   norm_factorC = NULL;
   rowContainsSingletonVariableC = NULL;
   singletonCoeffsColParent = NULL;
   singletonCoeffsColChild = NULL;
   normNnzRowA = NULL;
   normNnzRowC = NULL;
   normNnzColParent = NULL;
   normNnzColChild = NULL;
   norm_AmatTrans = NULL;
   norm_BmatTrans = NULL;
   norm_CmatTrans = NULL;
   norm_DmatTrans = NULL;
   currgChild = NULL;
   currgParent = NULL;
   gParentAdaptions = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec)->cloneFull();
   gParentAdaptions->setToZero();

   mA = 0;
   nA = 0;
}

StochPresolverParallelRows::~StochPresolverParallelRows()
{
 // todo
   rowsFirstHashTable = boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> >();
   rowsSecondHashTable  = boost::unordered_set<rowlib::rowWithEntries, boost::hash<rowlib::rowWithEntries> >();
   delete gParentAdaptions;

}

void StochPresolverParallelRows::applyPresolving()
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- Before Parallel Row Presolving:"<<endl;
   countRowsCols();
#endif

   if( myRank == 0 ) cout<<"Start Parallel Row Presolving..."<<endl;

   StochGenMatrix& matrixA = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& matrixC = dynamic_cast<StochGenMatrix&>(*(presProb->C));
   presData.resetRedCounters();
   gParentAdaptions->setToZero();
   clearNewBoundsParent();
   int nRowElims = 0;
   indivObjOffset = 0.0;

   // for children:
   for( size_t child_it = 0; child_it< matrixA.children.size(); child_it++)
   {
      // copy and normalize A,B,C,D and b,c,d:
      if( setNormalizedPointers((int)child_it, matrixA, matrixC ) )
      {
         // Prepare unordered set 'rowsFirstHashTable':
         rowsFirstHashTable.clear();

         // Per row, add row to the set 'rowsFirstHashTable':
         if( norm_Amat )
            insertRowsIntoHashtable( rowsFirstHashTable, *norm_Amat, norm_Bmat, EQUALITY_SYSTEM, *normNnzRowA );
         assert( (int)rowsFirstHashTable.size() <= mA );
         if( norm_Cmat )
            insertRowsIntoHashtable( rowsFirstHashTable, *norm_Cmat, norm_Dmat, INEQUALITY_SYSTEM, *normNnzRowC );
         assert( (int)rowsFirstHashTable.size() <= mA + norm_Cmat->m);

         // Second Hashing: Per bucket, do Second Hashing:
         rowsSecondHashTable.clear();
         for( size_t i = 0; i < rowsFirstHashTable.bucket_count(); ++i )
         {
            // skip bins with less than 2 elements:
            if( rowsFirstHashTable.bucket_size(i)<2 )
               continue;
            // insert elements from first Hash-bin into the second Hash-table:
            for( boost::unordered_set<rowlib::rowWithColInd>::local_iterator it =
                  rowsFirstHashTable.begin(i); it != rowsFirstHashTable.end(i); ++it )
            {
               rowsSecondHashTable.emplace(it->id, it->offset_nA, it->lengthA, it->colIndicesA, it->norm_entriesA,
                     it->lengthB, it->colIndicesB, it->norm_entriesB);
            }

            // Compare the rows in the final (from second hash) bin:
            compareRowsInSecondHashTable(nRowElims, (int)child_it);

            rowsSecondHashTable.clear();
         }

         // Objects created with new have to be deleted at the end of each child (norm_Amat etc)
         deleteNormalizedPointers((int)child_it, matrixA, matrixC);
      }
   }
   // combine the bounds of linking-variables:
   combineNewBoundsParent();
   // update the bounds of linking-variables:
   tightenLinkingVarsBounds();

   // update NnzColParent:
   updateNnzColParent(MPI_COMM_WORLD);
   presData.resetRedCounters();

   rowsFirstHashTable.clear();
   rowsSecondHashTable.clear();

   // for the A_0 and C_0 blocks:
   setNormalizedPointers(-1, matrixA, matrixC );
   // First Hashing: Fill 'rowsFirstHashTable':
   if( norm_Amat )
      insertRowsIntoHashtable( rowsFirstHashTable, *norm_Amat, NULL, EQUALITY_SYSTEM, *normNnzRowA );
   assert( (int)rowsFirstHashTable.size() <= mA );
   if( norm_Cmat )
      insertRowsIntoHashtable( rowsFirstHashTable, *norm_Cmat, NULL, INEQUALITY_SYSTEM, *normNnzRowC );
   assert( (int)rowsFirstHashTable.size() <= mA + norm_Cmat->m);
   // Second Hashing: Per bucket, do Second Hashing:
   for( size_t i = 0; i < rowsFirstHashTable.bucket_count(); ++i )
   {
      if( rowsFirstHashTable.bucket_size(i)<2 )
         continue;
      // insert elements from first Hash-bin into the second Hash-table:
      for( boost::unordered_set<rowlib::rowWithColInd>::local_iterator it =
            rowsFirstHashTable.begin(i); it != rowsFirstHashTable.end(i); ++it )
      {
         rowsSecondHashTable.emplace(it->id, it->offset_nA, it->lengthA, it->colIndicesA, it->norm_entriesA,
               it->lengthB, it->colIndicesB, it->norm_entriesB);
      }

      // Compare the rows in the final (from second hash) bin:
      compareRowsInSecondHashTable(nRowElims, -1);
      rowsSecondHashTable.clear();
   }
   // update objective coefficients of linking variables (gParent):
   allreduceAndUpdate(MPI_COMM_WORLD, *gParentAdaptions, *currgParent);

   deleteNormalizedPointers(-1, matrixA, matrixC);

   // update NnzColParent: As all processses have the same information, no communication necessary
   updateNnzUsingReductions(presData.nColElems->vec, presData.redCol->vec);
   presData.resetRedCounters();
   // synchronize nRowElims:
   synchronize(nRowElims);


   if( myRank == 0 )
      cout<<"Removed "<<nRowElims<<" Rows in Parallel Row Presolving."<<endl;

   // Sum up individual objOffset and then add it to the global objOffset:
   sumIndivObjOffset();
   presData.addObjOffset(indivObjOffset);

#ifdef TIMING
   if( myRank == 0 )
      cout<<"Global objOffset is now: "<<presData.getObjOffset()<<endl;
#endif

#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- After Parallel Row Presolving:"<<endl;
   countRowsCols();
#endif
}

/** If it is no dummy child, sets normalized pointers:
 * For the matrix blocks A, B, C and D and their transposed matrices and for the lhs/rhs.
 * Sets the pointers to currNnzRow, currNnzRowC, currNnzColChild, currRedColParent.
 * Sets mA and nA correctly.
 */
bool StochPresolverParallelRows::setNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC)
{
   // At the root
   if( it == -1 )
   {
      norm_Amat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.Bmat)->getStorageDynamicRef());
      norm_AmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.Bmat)->getStorageDynamicTransposedRef());
      norm_b = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec)->cloneFull();
      setCPAmatsRoot( presProb->A );
      norm_factorA = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec)->cloneFull();
      norm_factorA->setToZero();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
      normNnzRowA = currNnzRow->cloneFull();
      rowContainsSingletonVariableA = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec)->cloneFull();
      rowContainsSingletonVariableA->setToConstant( -1.0 );

      norm_Bmat = NULL;
      norm_BmatTrans = NULL;
      currBmat = NULL;
      currBmatTrans = NULL;

      norm_Cmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamicRef());
      norm_CmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamicTransposedRef());
      norm_cupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec)->cloneFull();
      norm_clow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vec)->cloneFull();
      norm_icupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vec)->cloneFull();
      norm_iclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vec)->cloneFull();
      norm_factorC = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec)->cloneFull();
      norm_factorC->setToZero();
      currCmat = dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamic();
      currCmatTrans = dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamicTransposed();
      currNnzRowC = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec);
      normNnzRowC = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec)->cloneFull();
      setCPRowRootIneqOnlyLhsRhs();
      rowContainsSingletonVariableC = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec)->cloneFull();
      rowContainsSingletonVariableC->setToConstant( -1.0 );

      norm_Dmat = NULL;
      norm_DmatTrans = NULL;
      currDmat = NULL;
      currDmatTrans = NULL;

      currNnzColParent = dynamic_cast<SimpleVector*>(presData.nColElems->vec);
      setCPColumnRoot();   // set pointers to currxlowParent etc. and to redColParent
      singletonCoeffsColParent = dynamic_cast<SimpleVector*>(presData.nColElems->vec)->cloneFull();
      singletonCoeffsColParent->setToZero();
      currNnzColChild = NULL;
      normNnzColParent = dynamic_cast<SimpleVector*>(presData.nColElems->vec)->cloneFull();
      normNnzColChild = NULL;

      currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);

      mA = norm_Amat->m;
      nA = norm_Amat->n;

      removeSingletonVars();

      if( norm_Amat )
         normalizeBlocksRowwise( EQUALITY_SYSTEM, *norm_Amat, NULL, *norm_b, NULL, NULL, NULL);
      if( norm_Cmat )
         normalizeBlocksRowwise( INEQUALITY_SYSTEM, *norm_Cmat, NULL, *norm_cupp, norm_clow, norm_icupp, norm_iclow);

      return true;
   }
   // else, check if it is no dummy child and copy the matrices:
   if( !childIsDummy(matrixA, it, EQUALITY_SYSTEM) )
   {
      norm_Amat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Amat)->getStorageDynamicRef());
      norm_AmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Amat)->getStorageDynamicTransposedRef());
      norm_Bmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Bmat)->getStorageDynamicRef());
      norm_BmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[it]->Bmat)->getStorageDynamicTransposedRef());
      norm_b = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec)->cloneFull();
      norm_factorA = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec)->cloneFull();
      norm_factorA->setToZero();

      setCPAmatsChild( presProb->A, it, EQUALITY_SYSTEM);
      setCPBmatsChild( presProb->A, it, EQUALITY_SYSTEM);
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
      normNnzRowA = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec)->cloneFull();

      // set all entries in rowContainsSingletonVariable to -1 to distinguish if there is a SV in a row and at which column index
      rowContainsSingletonVariableA = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec)->cloneFull();
      rowContainsSingletonVariableA->setToConstant( -1.0 );
   }
   else
   {
      norm_Amat = NULL;
      norm_AmatTrans = NULL;
      norm_Bmat = NULL;
      norm_BmatTrans = NULL;
      norm_b = NULL;
      norm_factorA = NULL;

      currAmat = NULL;
      currAmatTrans = NULL;
      currBmat = NULL;
      currBmatTrans = NULL;
      currNnzRow = NULL;
      normNnzRowA = NULL;

      rowContainsSingletonVariableA = NULL;
   }
   if( !childIsDummy(matrixC, it, INEQUALITY_SYSTEM) )
   {
      norm_Cmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Amat)->getStorageDynamicRef());
      norm_CmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Amat)->getStorageDynamicTransposedRef());
      norm_Dmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Bmat)->getStorageDynamicRef());
      norm_DmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Bmat)->getStorageDynamicTransposedRef());
      norm_cupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec)->cloneFull();
      norm_clow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec)->cloneFull();
      norm_icupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec)->cloneFull();
      norm_iclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec)->cloneFull();
      norm_factorC = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec)->cloneFull();
      norm_factorC->setToZero();

      currCmat = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Amat)->getStorageDynamic();
      currCmatTrans = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Amat)->getStorageDynamicTransposed();
      currDmat = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Bmat)->getStorageDynamic();
      currDmatTrans = dynamic_cast<SparseGenMatrix*>(matrixC.children[it]->Bmat)->getStorageDynamicTransposed();
      currNnzRowC = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
      normNnzRowC = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec)->cloneFull();

      setCPRowChildIneqOnlyLhsRhs(it); // set lhs and rhs

      rowContainsSingletonVariableC = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec)->cloneFull();
      rowContainsSingletonVariableC->setToConstant( -1.0 );
   }
   else
   {
      norm_Cmat = NULL;
      norm_CmatTrans = NULL;
      norm_Dmat = NULL;
      norm_DmatTrans = NULL;
      norm_cupp = NULL;
      norm_clow = NULL;
      norm_icupp = NULL;
      norm_iclow = NULL;
      norm_factorC = NULL;

      currCmat = NULL;
      currCmatTrans = NULL;
      currDmat = NULL;
      currDmatTrans = NULL;
      currNnzRowC = NULL;
      normNnzRowC = NULL;

      currIneqRhs = NULL;
      currIneqLhs = NULL;
      currIcupp = NULL;
      currIclow = NULL;

      rowContainsSingletonVariableC = NULL;
   }
   if( !norm_Amat && !norm_Cmat )   // case no child exists
   {
      currRedColParent = NULL;
      currNnzColChild = NULL;
      normNnzColChild = NULL;
      normNnzColParent = NULL;
      return false;
   }
   else  // set mA, nA correctly
   {
      mA = (norm_Amat) ? norm_Amat->m : 0;
      nA = (norm_Amat) ? norm_Amat->n : 0;

      setCPColumnChild(it);   // set pointers to currxlowChild etc.
      setCPColumnRoot();      // set pointers to currxlowParent etc. and to currRedColParent
      currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);

      normNnzColParent = dynamic_cast<SimpleVector*>(presData.nColElems->vec)->cloneFull();
      normNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec)->cloneFull();

      currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec);
      currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);

      singletonCoeffsColParent = dynamic_cast<SimpleVector*>(presData.nColElems->vec)->cloneFull();
      singletonCoeffsColParent->setToZero();
      singletonCoeffsColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec)->cloneFull();
      singletonCoeffsColChild->setToZero();
   }
   if(norm_Amat && norm_Cmat) // only for assertions
   {
      assert( norm_Amat->n == norm_Cmat->n );
      assert( norm_Bmat->n == norm_Dmat->n );
      assert( norm_Amat->m == norm_Bmat->m );
      assert( norm_Cmat->m == norm_Dmat->m );
   }

   // before normalizing, remove the singleton variables:
   removeSingletonVars();

   // normalize the rows:
   if( norm_Amat )
      normalizeBlocksRowwise( EQUALITY_SYSTEM, *norm_Amat, norm_Bmat, *norm_b, NULL, NULL, NULL);
   if( norm_Cmat )
      normalizeBlocksRowwise( INEQUALITY_SYSTEM, *norm_Cmat, norm_Dmat, *norm_cupp, norm_clow, norm_icupp, norm_iclow);

   return true;
}

void StochPresolverParallelRows::deleteNormalizedPointers(int it, StochGenMatrix& matrixA, StochGenMatrix& matrixC)
{
   delete normNnzColParent;
   delete singletonCoeffsColParent;

   if( it == -1 ) // Case at root
   {
      assert( norm_Amat && norm_b );
      delete norm_Amat;
      delete norm_AmatTrans;
      delete norm_b;
      assert( norm_Cmat && norm_cupp && norm_clow && norm_icupp && norm_iclow );
      delete norm_Cmat;
      delete norm_CmatTrans;
      delete norm_cupp;
      delete norm_clow;
      delete norm_icupp;
      delete norm_iclow;
      delete norm_factorA;
      delete norm_factorC;
      delete normNnzRowA;
      delete normNnzRowC;
      delete rowContainsSingletonVariableA;
      delete rowContainsSingletonVariableC;
      return;
   }
   bool childExists = false;
   if( !childIsDummy(matrixA, it, EQUALITY_SYSTEM) )
   {
      childExists = true;
      assert( norm_Amat && norm_Bmat && norm_b );
      delete norm_Amat;
      delete norm_AmatTrans;
      delete norm_Bmat;
      delete norm_BmatTrans;
      delete norm_b;
      delete norm_factorA;
      delete normNnzRowA;
      delete rowContainsSingletonVariableA;
   }
   if( !childIsDummy(matrixC, it, INEQUALITY_SYSTEM) )
   {
      childExists = true;
      assert( norm_Cmat && norm_Dmat && norm_cupp && norm_clow && norm_icupp && norm_iclow );
      delete norm_Cmat;
      delete norm_CmatTrans;
      delete norm_Dmat;
      delete norm_DmatTrans;
      delete norm_cupp;
      delete norm_clow;
      delete norm_icupp;
      delete norm_iclow;
      delete norm_factorC;
      delete normNnzRowC;
      delete rowContainsSingletonVariableC;
   }
   if( childExists )
   {
      delete normNnzColChild;
      delete singletonCoeffsColChild;
   }
}

void StochPresolverParallelRows::removeSingletonVars()
{
   // for the A_i blocks:
   for( int i = 0; i < normNnzColParent->n; i++ )
      {
         if( normNnzColParent->elements()[i] == 1.0 )
         {
            if( norm_AmatTrans &&
                  (norm_AmatTrans->rowptr[i].start +1 == norm_AmatTrans->rowptr[i].end) )
            {
               removeEntry(i, *rowContainsSingletonVariableA, *norm_Amat, *norm_AmatTrans,
                     *normNnzRowA, *normNnzColParent, LINKING_VARS_BLOCK);
            }
            else if(norm_CmatTrans &&
                  (norm_CmatTrans->rowptr[i].start +1 == norm_CmatTrans->rowptr[i].end) )
            {
               removeEntry(i, *rowContainsSingletonVariableC, *norm_Cmat, *norm_CmatTrans,
                     *normNnzRowC, *normNnzColParent, LINKING_VARS_BLOCK);
            }
            // else, the singleton entry is in one of the other A_i or C_i blocks
         }
      }

   // for the child block Bmat and Dmat:
   if( normNnzColChild )
   {
      for( int i = 0; i < normNnzColChild->n; i++ )
      {
         if( normNnzColChild->elements()[i] == 1.0 )
         {
            if( norm_BmatTrans &&
                  (norm_BmatTrans->rowptr[i].start +1 == norm_BmatTrans->rowptr[i].end) )
            {
               removeEntry(i, *rowContainsSingletonVariableA, *norm_Bmat, *norm_BmatTrans,
                     *normNnzRowA, *normNnzColChild, CHILD_BLOCK);
            }
            else if(norm_DmatTrans &&
                  (norm_DmatTrans->rowptr[i].start +1 == norm_DmatTrans->rowptr[i].end) )
            {
               removeEntry(i, *rowContainsSingletonVariableC, *norm_Dmat, *norm_DmatTrans,
                     *normNnzRowC, *normNnzColChild, CHILD_BLOCK);
            }
            // else, the singleton entry is in the linking block Blmat
         }
      }
   }
}

/** Removes a singleton entry in column colIdx in Bblock and BblockTrans and adapts
 * the nnz vectors nnzRow and nnzColChild accordingly. Sets the entries in
 * rowContainsSingletonVar to the corresponding column index in which the singleton entry occurs.
 */
void StochPresolverParallelRows::removeEntry(int colIdx, SimpleVector& rowContainsSingletonVar,
      SparseStorageDynamic& matrix, SparseStorageDynamic& matrixTrans, SimpleVector& nnzRow, SimpleVector& nnzCol,
      BlockType block_type)
{
   assert( colIdx >= 0 && colIdx < matrixTrans.m );
   assert( matrixTrans.rowptr[colIdx].start +1 == matrixTrans.rowptr[colIdx].end);
   assert( nnzRow.n == matrix.m );
   assert( matrix.n == nnzCol.n );
   assert( nnzCol.n == matrixTrans.m );
   assert( nnzCol.elements()[colIdx] == 1.0 );

   // First, find indices of the singleton entry:
   const int k = matrixTrans.rowptr[colIdx].start;
   const int rowIdx = matrixTrans.jcolM[k];
   assert( rowIdx < nnzRow.n );

   // check if there are no more than one singleton Entries in this row:
   if( rowContainsSingletonVar.elements()[rowIdx] >= 0 )
   {
      rowContainsSingletonVar.elements()[rowIdx] = -2.0;
      return;
   }
   // store the colIdx in rowContainsSingletonVar:
   // (possibly add offset to colIdx so that Amat and Bmat are distinct)
   if( block_type == LINKING_VARS_BLOCK )
      rowContainsSingletonVar.elements()[rowIdx] = colIdx;
   else // if( block_type == CHILD_BLOCK )
      rowContainsSingletonVar.elements()[rowIdx] = colIdx + nA;

   // Second, remove the entry from norm_Bmat and norm_BmatTrans:
   double coeff = 0;
   removeEntryInDynamicStorage(matrix, rowIdx, colIdx, coeff);
   matrixTrans.rowptr[colIdx].end --;
   nnzRow.elements()[rowIdx]--;
   nnzCol.elements()[colIdx] = 0.0;
   if( block_type == LINKING_VARS_BLOCK )
      singletonCoeffsColParent->elements()[colIdx] = coeff;
   else
      singletonCoeffsColChild->elements()[colIdx] = coeff;
}

void StochPresolverParallelRows::normalizeBlocksRowwise( SystemType system_type,
      SparseStorageDynamic& Ablock, SparseStorageDynamic* Bblock,
      SimpleVector& Rhs, SimpleVector* Lhs, SimpleVector* iRhs, SimpleVector* iLhs)
{
   int nRows = Ablock.m;
   if(Bblock) assert( nRows == Bblock->m);
   if( nRows != Rhs.n )
      cout<<"nRows="<<nRows<<", Rhs->n="<<Rhs.n<<endl;
   assert( nRows == Rhs.n );

   if( system_type == INEQUALITY_SYSTEM )
   {
      assert( Lhs && iRhs && iLhs );
      assert( Lhs->n == nRows && iRhs->n == nRows && iLhs->n == nRows );
   }

   for( int i=0; i<nRows; i++)   // row i
   {
      double maxValue = 0.0;
      bool negateRow = false;
      const int rowStartA = Ablock.rowptr[i].start;
      const int rowEndA = Ablock.rowptr[i].end;
      if( rowStartA < rowEndA )
      {
         if( Ablock.M[rowStartA] < 0)
            negateRow = true;
         for(int k=rowStartA; k<rowEndA; k++)
         {
            if( fabs(Ablock.M[k]) > maxValue )
               maxValue = fabs(Ablock.M[k]);
         }
      }
      int rowStartB = 0;
      int rowEndB = 0;
      if(Bblock)
      {
         rowStartB = Bblock->rowptr[i].start;
         rowEndB = Bblock->rowptr[i].end;
         if( rowStartB < rowEndB )
         {
            if( Bblock->M[rowStartB] < 0)
               negateRow = true;
            for(int k=rowStartB; k<rowEndB; k++)
            {
               if( fabs(Bblock->M[k]) > maxValue )
                  maxValue = fabs(Bblock->M[k]);
            }
         }
      }
      // normalize the row by dividing all entries by maxValue and possibly by -1, if negateRow.
      if(negateRow)
         maxValue *= -1.0;
      if( rowStartA < rowEndA )
      {
         for(int k=rowStartA; k<rowEndA; k++)
            Ablock.M[k] /= maxValue;
      }
      if( Bblock && (rowStartB < rowEndB) )
      {
         for(int k=rowStartB; k<rowEndB; k++)
            Bblock->M[k] /= maxValue;
      }
      if( system_type == EQUALITY_SYSTEM )
      {
         Rhs.elements()[i] /= maxValue;
         norm_factorA->elements()[i] = maxValue;
      }
      else
      {
         if( iLhs->elements()[i] != 0.0 )
            Lhs->elements()[i] /= maxValue;
         if( iRhs->elements()[i] != 0.0 )
            Rhs.elements()[i] /= maxValue;
         // if multiplied by a negative value, lhs and rhs have to be swapped:
         if(negateRow)
         {
            std::swap(Lhs->elements()[i], Rhs.elements()[i]);
            std::swap(iLhs->elements()[i], iRhs->elements()[i]);
         }
         norm_factorC->elements()[i] = maxValue;
      }
   }
}

/** Inserts all non-empty rows into the given unordered set 'rows'.
 * Ablock and Bblock are supposed to be the two matrix blocks of a child.
 * If at root, Bblock should be NULL.
 */
void StochPresolverParallelRows::insertRowsIntoHashtable( boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> > &rows,
      SparseStorageDynamic& Ablock, SparseStorageDynamic* Bblock, SystemType system_type, SimpleVector& nnzRow )
{
   if( Bblock) assert( Ablock.m == Bblock->m );
   if( system_type == EQUALITY_SYSTEM )
      assert( mA == Ablock.m );

   for(int i=0; i<Ablock.m; i++)
   {
      // ignore rows containing more than one singleton entry:
      if( system_type == EQUALITY_SYSTEM && rowContainsSingletonVariableA->elements()[i] == -2.0 )
         continue;
      if( system_type == INEQUALITY_SYSTEM && rowContainsSingletonVariableC->elements()[i] == -2.0 )
         continue;

      // calculate rowId including possible offset (for Inequality rows):
      int rowId = i;
      if( nnzRow.elements()[rowId] == 0.0 )
         continue;
      if( system_type == INEQUALITY_SYSTEM )
         rowId += mA;

      // calculate rowlength of A- and B-block
      const int rowStartA = Ablock.rowptr[i].start;
      const int rowlengthA = Ablock.rowptr[i].end - rowStartA;
      if( Bblock ){
         const int rowStartB = Bblock->rowptr[i].start;
         const int rowlengthB = Bblock->rowptr[i].end - rowStartB;

         // colIndices and normalized entries are set as pointers to the original data.
         // create and insert the new element:
         assert( rowlengthA != 0 || rowlengthB != 0 );
         rows.emplace(rowId, nA, rowlengthA, &(Ablock.jcolM[rowStartA]), &(Ablock.M[rowStartA]),
               rowlengthB, &(Bblock->jcolM[rowStartB]), &(Bblock->M[rowStartB]));
      }
      else
      {
         assert( rowlengthA != 0 );
         rows.emplace(rowId, nA, rowlengthA, &(Ablock.jcolM[rowStartA]),
               &(Ablock.M[rowStartA]), 0, NULL, NULL);
      }
   }
}

/*
 * Per bucket in rowsSecondHashTable, compare the containing rows and check if they are parallel.
 * If so, consider the different possible cases. If a row can be removed, the action is applied
 * to the original matrices (not the normalized copies).
 * @param it represents the child number (or -1 for parent blocks)
 */
void StochPresolverParallelRows::compareRowsInSecondHashTable(int& nRowElims, int it)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( rowsSecondHashTable.empty() )
      return;

   for (size_t i=0; i<rowsSecondHashTable.bucket_count(); ++i)
   {
      for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator it1 = rowsSecondHashTable.begin(i);
            it1!=rowsSecondHashTable.end(i); ++it1)
      {
         // either pairwise comparison OR lexicographical sorting and then compare only neighbors.
         // Here: parwise comparison:
         for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator it2 = it1;
                                    it2!=rowsSecondHashTable.end(i); ++it2)
         {
            // When two parallel rows are found, check if they are both =, both <=, or = and <=
            if( checkRowsAreParallel( *it1, *it2) )
            {
               //cout<<"Found two rows with parallel coefficients: Row "<<it1->id<<"("<<it1->lengthA+it1->lengthB<<")"<<" and "<<it2->id<<"("<<it2->lengthA+it2->lengthB<<")"<<endl;
               if( it1->id < mA && it2->id < mA )
               {
                  // Case both constraints are equalities:
                  // check if one of them contains a singleton variable:
                  if( rowContainsSingletonVariableA->elements()[it1->id] != -1.0
                        || rowContainsSingletonVariableA->elements()[it2->id] != -1.0 )
                  {
                     // nearly parallel case 1:
                     //PIPSdebugMessage("Nearly Parallel Rows, case 1. \n");
                     // make sure that a_r2 != 0, otherwise switch ids.
                     int id1 = it1->id;
                     int id2 = it2->id;
                     bool swappedId1Id2 = false;
                     if( rowContainsSingletonVariableA->elements()[id2] == -1.0 )
                     {   // swap ids so that id2 has a singleton entry:
                        std::swap(id1, id2);
                        swappedId1Id2 = true;
                     }
                     assert( rowContainsSingletonVariableA->elements()[id2] != -1.0 );
                     // if row id1 was already deleted, do not continue the procedure:
                     if( (id1 < mA && currNnzRow->elements()[id1] == 0.0)
                           || (id1 > mA && currNnzRowC->elements()[id1 - mA] == 0.0) )
                        continue;

                     doNearlyParallelRowCase1(id1, id2, it);
                     // delete the inequality constraint id2 (which could be either it1 or it2 now):
                     (swappedId1Id2) ?
                           eliminateOriginalRow((int)it1->id, nRowElims) : eliminateOriginalRow((int)it2->id, nRowElims);
                  }
                  else
                  {
                     //PIPSdebugMessage("Really Parallel Rows, case 1. \n");
                     if(  !PIPSisEQ( norm_b->elements()[it1->id], norm_b->elements()[it2->id]) )
                        abortInfeasible(MPI_COMM_WORLD);

                     // delete row2 in the original system:
                     eliminateOriginalRow((int) it2->id, nRowElims);
                  }
               }
               else if( it1->id >= mA && it2->id >= mA )
               {
                  // Case both constraints are inequalities
                  const int rowId1 = it1->id - mA;
                  const int rowId2 = it2->id - mA;

                  if( rowContainsSingletonVariableC->elements()[rowId1] != -1.0
                        || rowContainsSingletonVariableC->elements()[rowId2] != -1.0 )
                  {
                     // nearly parallel case 3:
                     //PIPSdebugMessage("Nearly Parallel Rows, case 3. \n");
                     if( rowContainsSingletonVariableC->elements()[rowId1] != -1.0
                        && rowContainsSingletonVariableC->elements()[rowId2] != -1.0 )
                     {
                        if( doNearlyParallelRowCase3(rowId1, rowId2, it) )
                           eliminateOriginalRow((int)it2->id, nRowElims);
                     }
                  }
                  else
                  {
                     //PIPSdebugMessage("Really Parallel Rows, case 3. \n");
                     // tighten bounds in original and normalized system:
                     tightenOriginalBoundsOfRow1(rowId1, rowId2);

                     // delete row2 in the original system:
                     eliminateOriginalRow((int)it2->id, nRowElims);
                  }
               }
               else
               {  // Case one constraint is an equality, one an inequality
                  int id1 = it1->id;
                  int id2 = it2->id;
                  bool swappedId1Id2 = false;
                  if( id1 >= mA )   // swap ids so that id2 is the inequality constraint.
                  {
                     std::swap(id1, id2);
                     swappedId1Id2 = true;
                  }
                  const int ineqRowId = id2 - mA;

                  if( rowContainsSingletonVariableA->elements()[id1] == -1.0
                        && rowContainsSingletonVariableC->elements()[ineqRowId] == -1.0 )
                  {
                     //PIPSdebugMessage("Really Parallel Rows, case 2. \n");
                     // check for infeasibility:
                     if( norm_iclow->elements()[ineqRowId] != 0.0
                           && norm_clow->elements()[ineqRowId] > norm_b->elements()[id1] )
                        abortInfeasible(MPI_COMM_WORLD);
                     if( norm_icupp->elements()[ineqRowId] != 0.0
                           && norm_cupp->elements()[ineqRowId] < norm_b->elements()[id1] )
                        abortInfeasible(MPI_COMM_WORLD);
                     // remove inequality row.
                  }
                  else if( rowContainsSingletonVariableA->elements()[id1] != -1.0
                        && rowContainsSingletonVariableC->elements()[ineqRowId] == -1.0 )
                  {
                     // nearly parallel case 2:
                     //PIPSdebugMessage("Nearly Parallel Rows, case 2. \n");
                     // compute the new variable bound for x_id1:
                     const int singleColIdx = rowContainsSingletonVariableA->elements()[id1];
                     double coeff_singleton = getSingletonCoefficient(singleColIdx);

                     double newxlow = - std::numeric_limits<double>::max();
                     double newxupp = std::numeric_limits<double>::max();

                     // check if f_q * a_q is positive or negative:
                     double faq = norm_factorA->elements()[id1] * coeff_singleton;
                     if( faq > 0 )
                     {
                        if( norm_iclow->elements()[ineqRowId] != 0.0 )
                           newxupp = ( norm_b->elements()[id1] - norm_clow->elements()[ineqRowId] ) * norm_factorA->elements()[id1] / coeff_singleton;
                        if( norm_icupp->elements()[ineqRowId] != 0.0 )
                           newxlow = ( norm_b->elements()[id1] - norm_cupp->elements()[ineqRowId] ) * norm_factorA->elements()[id1] / coeff_singleton;
                     }
                     else  //if( fq*aq < 0 )
                     {
                        if( norm_iclow->elements()[ineqRowId] != 0.0 )
                           newxlow = ( norm_b->elements()[id1] - norm_clow->elements()[ineqRowId] ) * norm_factorA->elements()[id1] / coeff_singleton;
                        if( norm_icupp->elements()[ineqRowId] != 0.0 )
                           newxupp = ( norm_b->elements()[id1] - norm_cupp->elements()[ineqRowId] ) * norm_factorA->elements()[id1] / coeff_singleton;
                     }

                     // effectively tighten bounds of variable x_id1:
                     tightenBoundsForSingleVar(singleColIdx, newxlow, newxupp);

                     // remove inequality row.

                  }
                  else
                     continue;   // no action if the inequality contains a singleton entry.

                  // in both the parallel row case and in the nearly parallel row case,
                  // delete the inequality constraint id2 (which could be either it1 or it2 now):
                  (swappedId1Id2) ?
                        eliminateOriginalRow((int)it1->id, nRowElims) : eliminateOriginalRow((int)it2->id, nRowElims);

               }
            }
         }
      }
   }
}

/**
 * Compare two rowWithEntries rows if the normalized coefficients are the same, at the same
 * columns indices. If yes, return true. Else, return false.
 */
bool StochPresolverParallelRows::checkRowsAreParallel( rowlib::rowWithEntries row1, rowlib::rowWithEntries row2)
{
   assert( row1.id >= 0 && row2.id >= 0 );
   if( row1.id == row2.id )
      return false;

   if( row1.lengthA != row2.lengthA || row1.lengthB != row2.lengthB )
      return false;
   for( int i=0; i<row1.lengthA; i++)
   {
      if( row1.colIndicesA[i] != row2.colIndicesA[i] )
         return false;
      if( !PIPSisEQ(row1.norm_entriesA[i], row2.norm_entriesA[i], tol_compare_double) )
         return false;
   }
   for( int i=0; i<row1.lengthB; i++)
   {
      if( row1.colIndicesB[i] != row2.colIndicesB[i] )
         return false;
      if( !PIPSisEQ(row1.norm_entriesB[i], row2.norm_entriesB[i], tol_compare_double) )
         return false;
   }
   return true;
}

/** Eliminate the row in the original system. */
void StochPresolverParallelRows::eliminateOriginalRow(int rowId, int& nRowElims)
{
   assert(rowId>=0);
   if(rowId < mA) // equality row
   {
      assert( norm_Amat );
      if(currNnzRow->elements()[rowId] != 0.0) // check if row was already removed
      {
//         cout<<"Delete row in A with id# "<<rowId<<endl;
         nRowElims++;
         removeRow(rowId, *currAmat, *currAmatTrans, currBmat, currBmatTrans,
               *currNnzRow, *currRedColParent, currNnzColChild);
      }
   }
   else  // inequality row
   {
      assert( norm_Cmat );
      assert(rowId < mA + norm_Cmat->m );
      const int rowIdC = rowId - mA;
      if(currNnzRowC->elements()[rowIdC] != 0.0) // check if row was already removed
      {
//         cout<<"Delete row in C with id# "<<rowIdC<<endl;
         nRowElims++;
         removeRow(rowIdC, *currCmat, *currCmatTrans, currDmat, currDmatTrans,
               *currNnzRowC, *currRedColParent, currNnzColChild);
      }
   }
}

/**
 * Tightens the original lower and upper bounds of the first row, given the lower
 * and upper bounds of the second row. The normalized bounds are compared and the
 * normalizing factor of row1 is used to determine which bound can be tightened to
 * which value. The normalized bounds of row1 are also updated.
 * Assumes that both rows are inequality constraints.
 */
void StochPresolverParallelRows::tightenOriginalBoundsOfRow1(int rowId1, int rowId2)
{
   assert( currCmat );
   assert( rowId1 >= 0 && rowId2 >= 0 );
   assert( rowId1 < currCmat->m && rowId2 < currCmat->m );
   assert( norm_factorC && norm_factorC->n == currCmat->m );

   // the normalizing factor of the first inequality:
   double factor = norm_factorC->elements()[rowId1];

   // Set norm_low/norm_upp to the bounds of the second inequality:
   double norm_low = -std::numeric_limits<double>::max();
   double norm_upp = std::numeric_limits<double>::max();
   if( norm_iclow->elements()[rowId2] != 0.0 )
      norm_low = norm_clow->elements()[rowId2];
   if( norm_icupp->elements()[rowId2] != 0.0 )
      norm_upp = norm_cupp->elements()[rowId2];

   // test for infeasibility:
   if( (norm_iclow->elements()[rowId1] != 0.0 && PIPSisLT(norm_upp, norm_clow->elements()[rowId1]) )
         || (norm_icupp->elements()[rowId1] != 0.0 && PIPSisLT(norm_cupp->elements()[rowId1], norm_low) ) )
   {
      cout<<"Detected infeasibility during parallel row presolving."<<endl;
      abortInfeasible(MPI_COMM_WORLD);
   }
   // test if the new bounds are tightening:
   if( ( norm_iclow->elements()[rowId1] != 0.0 && PIPSisLT(norm_clow->elements()[rowId1], norm_low) )
         || ( norm_iclow->elements()[rowId1] == 0.0 && norm_low > -std::numeric_limits<double>::max() ))
   {
      // tighten normalized lower bound:
      setNewBound( rowId1, norm_low, norm_clow, norm_iclow);

      // tighten original lower (factor>0) or upper ((factor<0) bound:
      (factor>0) ? setNewBound( rowId1, factor * norm_low, currIneqLhs, currIclow) :
            setNewBound( rowId1, factor * norm_low, currIneqRhs, currIcupp);
   }
   if( ( norm_icupp->elements()[rowId1] != 0.0 && PIPSisLT(norm_upp, norm_cupp->elements()[rowId1]) )
         || ( norm_icupp->elements()[rowId1] == 0.0 && norm_upp < std::numeric_limits<double>::max() ))
   {
      // tighten normalized upper bound:
      setNewBound( rowId1, norm_upp, norm_cupp, norm_icupp);

      // tighten original upper (factor>0) or lower (factor<0) bound:
      (factor>0) ? setNewBound( rowId1, factor * norm_upp, currIneqRhs, currIcupp) :
         setNewBound( rowId1, factor * norm_upp, currIneqLhs, currIclow);
   }
}

/** Returns the matrix coefficient of the singleton variable with index singleColIdx.
 * This index is possibly offset by nA, if it is in a B- or D-block.
 * If singleColIdx == -1, then 0.0 is returned.
 */
double StochPresolverParallelRows::getSingletonCoefficient(int singleColIdx)
{
   assert( singleColIdx >= -1.0);
   if( singleColIdx >= nA )
   {
      assert( singletonCoeffsColChild );
      assert( singleColIdx < nA + singletonCoeffsColChild->n );
   }
   else
      assert( singletonCoeffsColParent );

   if( singleColIdx == -1.0 )
      return 0.0;

   return ( singleColIdx >= nA ) ? singletonCoeffsColChild->elements()[singleColIdx - nA] :
         singletonCoeffsColParent->elements()[singleColIdx];
}

/** Tightens the bounds for a singleton variable with index singleColIdx.
 * @param newxlow and @param newxupp are the new candidate bounds.
 * This index is possibly offset by nA, if it is in a B- or D-block.
 */
void StochPresolverParallelRows::tightenBoundsForSingleVar(int singleColIdx, double newxlow, double newxupp)
{
   assert( singleColIdx >= 0 );

   if( singleColIdx < nA )   // variable x_1 is in the parent block Amat
   {
      setNewBounds(singleColIdx, newxlow, newxupp, currIxlowParent->elements(),
            currxlowParent->elements(), currIxuppParent->elements(), currxuppParent->elements());
   }
   else   // variable x_1 is in the child block Bmat
   {
      assert( singleColIdx < nA + currxuppChild->n );
      setNewBounds(singleColIdx-nA, newxlow, newxupp, currIxlowChild->elements(),
            currxlowChild->elements(), currIxuppChild->elements(), currxuppChild->elements());
   }
}

bool StochPresolverParallelRows::doNearlyParallelRowCase1(int rowId1, int rowId2, int it)
{
   assert( rowId1 >= 0 && rowId1 < mA );
   assert( rowId2 >= 0 && rowId2 < mA );
   assert( it >= -1 && it < nChildren );
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   // calculate t and d:
   const int singleColIdx1 = rowContainsSingletonVariableA->elements()[rowId1];
   const int singleColIdx2 = rowContainsSingletonVariableA->elements()[rowId2];

   const double coeff_singleton1 = getSingletonCoefficient(singleColIdx1);
   const double coeff_singleton2 = getSingletonCoefficient(singleColIdx2);
   const double s = norm_factorA->elements()[rowId1] / norm_factorA->elements()[rowId2];
   const double t = coeff_singleton1 / coeff_singleton2 / s;
   const double d = (norm_b->elements()[rowId2] - norm_b->elements()[rowId1])
            * norm_factorA->elements()[rowId2] / coeff_singleton2;

   if( singleColIdx1 != -1.0 )
   {
      // tighten the bounds of variable x_1:
      double newxlow = - std::numeric_limits<double>::max();
      double newxupp = std::numeric_limits<double>::max();

      // calculate new bounds depending on the sign of t:
      if( t > 0 )
      {
         if( singleColIdx2 >= nA )
         {
            computeXminusYdivZ( newxlow, *currIxlowChild, *currxlowChild, singleColIdx2 - nA, d, t);
            computeXminusYdivZ( newxupp, *currIxuppChild, *currxuppChild, singleColIdx2 - nA, d, t);
         }
         else
         {
            computeXminusYdivZ( newxlow, *currIxlowParent, *currxlowParent, singleColIdx2, d, t);
            computeXminusYdivZ( newxupp, *currIxuppParent, *currxuppParent, singleColIdx2, d, t);
         }
      }
      else if( t < 0 )
      {
         if( singleColIdx2 >= nA )
         {
            computeXminusYdivZ( newxlow, *currIxuppChild, *currxuppChild, singleColIdx2 - nA, d, t);
            computeXminusYdivZ( newxupp, *currIxlowChild, *currxlowChild, singleColIdx2 - nA, d, t);
         }
         else
         {
            computeXminusYdivZ( newxlow, *currIxuppParent, *currxuppParent, singleColIdx2, d, t);
            computeXminusYdivZ( newxupp, *currIxlowParent, *currxlowParent, singleColIdx2, d, t);
         }
      }

      // effectively tighten bounds of variable x_id1:
      if( singleColIdx1 >= nA || it == -1 )
         tightenBoundsForSingleVar(singleColIdx1, newxlow, newxupp);
      else  // if singleColIdx1 < nA && it > -1
         storeNewBoundsParent(singleColIdx1, newxlow, newxupp);

      // adapt objective function:
      adaptObjective( singleColIdx1, singleColIdx2, t, d, it);
   }
   else
   {
      adaptObjective( -1, singleColIdx2, t, d, it);
   }
   return true;
}
/**
 * Executes the Nearly Parallel Row Case 3: Both constraints are inequalities and both contain
 * a singleton variable entry.
 * The row indices rowId1, rowId2 are already de-offset, so they should be in the range [0,Cmat->m).
 */
bool StochPresolverParallelRows::doNearlyParallelRowCase3(int rowId1, int rowId2, int it)
{
   assert( rowId1 >= 0 && rowId2 >= 0 );
   assert( rowContainsSingletonVariableC );
   assert( rowId1 < rowContainsSingletonVariableC->n );
   assert( rowId2 < rowContainsSingletonVariableC->n );
   assert( norm_factorC );

   // First, do all the checks to verify if the third case applies and all conditions are met.
   // check s>0:
   const double s = norm_factorC->elements()[rowId1] / norm_factorC->elements()[rowId2];
   if( s <= 0)
      return false;

   // check a_q!=0.0, a_r!=0.0:
   const int singleColIdx1 = rowContainsSingletonVariableC->elements()[rowId1];
   const int singleColIdx2 = rowContainsSingletonVariableC->elements()[rowId2];
   const double coeff_singleton1 = getSingletonCoefficient(singleColIdx1);
   const double coeff_singleton2 = getSingletonCoefficient(singleColIdx2);
   if( coeff_singleton1 == 0.0 || coeff_singleton2 == 0.0 )
      return false;

   // check d_q == s * d_r (here s==1, so just d_q == d_r:
   if( (norm_iclow->elements()[rowId1] != 0.0 && norm_iclow->elements()[rowId1] == 0.0 )
      || (norm_iclow->elements()[rowId1] == 0.0 && norm_iclow->elements()[rowId1] != 0.0) )
      return false;
   if( (norm_iclow->elements()[rowId1] != 0.0 && norm_iclow->elements()[rowId1] != 0.0)
         && ( !PIPSisEQ(norm_clow->elements()[rowId1], norm_clow->elements()[rowId2], tol_compare_double) ) )
      return false;

   // check f_q == f_r
   if( (norm_icupp->elements()[rowId1] != 0.0 && norm_icupp->elements()[rowId1] == 0.0 )
      || (norm_icupp->elements()[rowId1] == 0.0 && norm_icupp->elements()[rowId1] != 0.0) )
      return false;
   if( (norm_icupp->elements()[rowId1] != 0.0 && norm_icupp->elements()[rowId1] != 0.0)
         && ( !PIPSisEQ(norm_cupp->elements()[rowId1], norm_cupp->elements()[rowId2], tol_compare_double) ) )
      return false;

   // check c_1*c_2 >= 0:
   const double c1 = (singleColIdx1 < nA) ? currgParent->elements()[singleColIdx1] : currgChild->elements()[singleColIdx1 - nA];
   const double c2 = (singleColIdx2 < nA) ? currgParent->elements()[singleColIdx2] : currgChild->elements()[singleColIdx2 - nA];
   if( c1 * c2 < 0 )
      return false;

   // check lower and upper bounds:
   double ixlow1 = (singleColIdx1 < nA) ? currIxlowParent->elements()[singleColIdx1] : currIxlowChild->elements()[singleColIdx1 - nA];
   double ixupp1 = (singleColIdx1 < nA) ? currIxuppParent->elements()[singleColIdx1] : currIxuppChild->elements()[singleColIdx1 - nA];
   double ixlow2 = (singleColIdx2 < nA) ? currIxlowParent->elements()[singleColIdx2] : currIxlowChild->elements()[singleColIdx2 - nA];
   double ixupp2 = (singleColIdx2 < nA) ? currIxuppParent->elements()[singleColIdx2] : currIxuppChild->elements()[singleColIdx2 - nA];

   double xlow1 = (singleColIdx1 < nA) ? currxlowParent->elements()[singleColIdx1] : currxlowChild->elements()[singleColIdx1 - nA];
   double xupp1 = (singleColIdx1 < nA) ? currxuppParent->elements()[singleColIdx1] : currxuppChild->elements()[singleColIdx1 - nA];
   double xlow2 = (singleColIdx2 < nA) ? currxlowParent->elements()[singleColIdx2] : currxlowChild->elements()[singleColIdx2 - nA];
   double xupp2 = (singleColIdx2 < nA) ? currxuppParent->elements()[singleColIdx2] : currxuppChild->elements()[singleColIdx2 - nA];

   // check aq * l1 = s * ar * l2:
   if( ixlow1 != 0.0 && ixlow2 != 0.0 &&
         !PIPSisEQ(coeff_singleton1 * xlow1, s* coeff_singleton2 * xlow2, tol_compare_double))
      return false;
   if( ixlow1 != 0.0 && ixlow2 == 0.0 && xlow1 > -std::numeric_limits<double>::max() )
      return false;
   if( ixlow1 == 0.0 && ixlow2 != 0.0 && xlow2 > -std::numeric_limits<double>::max() )
      return false;

   // check aq * u1 = s * ar * u2:
   if( ixupp1 != 0.0 && ixupp2 != 0.0 &&
         !PIPSisEQ(coeff_singleton1 * xupp1, s* coeff_singleton2 * xupp2, tol_compare_double))
      return false;
   if( ixupp1 != 0.0 && ixupp2 == 0.0 && xupp1 < std::numeric_limits<double>::max() )
      return false;
   if( ixupp1 == 0.0 && ixupp2 != 0.0 && xupp2 < std::numeric_limits<double>::max() )
      return false;

   // Second, aggregate x_2: adapt objectiveCost(x_1)
   const double t = coeff_singleton1 / coeff_singleton2 / s;

   adaptObjective( singleColIdx1, singleColIdx2, t, 0.0, it);

   return true;
}

/**
 * Adapt the objective function, at coefficient of colIdx1 and if @param offset, then also the objOffset is adapted.
 */
void StochPresolverParallelRows::adaptObjective( int colIdx1, int colIdx2, double t, double d, int it)
{
   assert( colIdx1 >= -1 && colIdx2 >= 0 );
   assert( it >= -1 && it < nChildren );
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   const double costOfVar2 = (colIdx2 < nA)
         ? currgParent->elements()[colIdx2] : currgChild->elements()[colIdx2 - nA];

   if( colIdx1 == -1 )
   {  // in case a_q1 == 0.0 ( singleColIdx1 == -1.0), only adapt objective offset:
      if( myRank == 0 || it > -1 ) // only add the objective offset for root as process ZERO:
         indivObjOffset += d * costOfVar2;
   }
   else if( colIdx1 >= nA )
   {
      currgChild->elements()[colIdx1 - nA] += t * costOfVar2;
      indivObjOffset += d * costOfVar2;
   }
   else if( it > -1 )   // case singleColIdx1 a linking variable, and currently working on a child
   {  // do not adapt currgParent directly, but save the adaption to communicate it at the end:
      gParentAdaptions->elements()[colIdx1] += t * costOfVar2;
      indivObjOffset += d * costOfVar2;
   }
   else  // case Parent block that all processes have
   {
      currgParent->elements()[colIdx1] += t * costOfVar2;
      // only add the objective offset for root as process ZERO:
      if( myRank == 0 ) indivObjOffset += d * costOfVar2;
   }
}

void StochPresolverParallelRows::computeXminusYdivZ( double& result,
      SimpleVector& ixvec, SimpleVector& xvec, int index, double y, double z)
{
   if( ixvec.elements()[index] != 0.0 )
      result = ( xvec.elements()[index] - y ) / z;
}

/** Should be called right after combineNewBoundsParent().
 * Updates the bounds on the linking variables.
 */
void StochPresolverParallelRows::tightenLinkingVarsBounds()
{
   setCPColumnRoot();

   // apply updated newBoundsParent to the variable bounds.
   for(int i=0; i<getNumberNewBoundsParent(); i++)
   {
      XBOUNDS newbounds = getNewBoundsParent(i);
      assert( newbounds.colIdx >= 0 && newbounds.colIdx < nA );
      tightenBoundsForSingleVar(newbounds.colIdx, newbounds.newxlow, newbounds.newxupp);
   }
   clearNewBoundsParent();
}
