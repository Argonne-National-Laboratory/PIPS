/*
 * StochPresolverParallelRows.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

// todo: should we ever decide to switch to a newer c++ standard - this can be optimized

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

StochPresolverParallelRows::StochPresolverParallelRows(PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb)
{
   mA = 0;
   nA = 0;
}

StochPresolverParallelRows::~StochPresolverParallelRows()
{
 // todo : check what this line does
   row_support_hashtable = boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> >();
   row_coefficients_hashtable  = boost::unordered_set<rowlib::rowWithEntries, boost::hash<rowlib::rowWithEntries> >();
}

/// presolve assumes that all rows are in their correct blocks -> linking rows are not pure local/linking rows with one singleton column cannot be local up to that singleton column
/// linking variables not in A0/C0 cannot be completely in the linking vars block etc.
void StochPresolverParallelRows::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before parallel row presolving:" << std::endl;
   }
   countRowsCols();
#endif

   /// first hash support of all rows then per hashbucket hash coeffs to find (nearly) parallel rows
   int nRowElims = 0;

   /// non-linking non-root part of matrices
   /// since we assume that all linking rows actually are pure linking rows we need not consider them here
   for(int node = 0; node < nChildren; ++node)
   {
      if( !presData.nodeIsDummy(node, EQUALITY_SYSTEM) || !presData.nodeIsDummy(node, INEQUALITY_SYSTEM) )
      {
         /// copy and normalize A_i, B_i, C_i, D_i and b_i, clow_i, cupp_i
         setNormalizedPointers(node);

         row_support_hashtable.clear();

         // Per row, add row to the set 'row_support_hashtable':
         if( !presData.nodeIsDummy(node, EQUALITY_SYSTEM) )
         {
            assert(norm_Amat);
            assert(norm_Bmat);
            assert(normNnzRowA);
            insertRowsIntoHashtable( row_support_hashtable, norm_Amat, norm_Bmat, EQUALITY_SYSTEM, normNnzRowA );
         }
         assert( static_cast<int>(row_support_hashtable.size()) <= mA );

         if( !presData.nodeIsDummy(node, INEQUALITY_SYSTEM) )
         {
            assert(norm_Cmat);
            assert(norm_Dmat);
            assert(normNnzRowC);
            insertRowsIntoHashtable( row_support_hashtable, norm_Cmat, norm_Dmat, INEQUALITY_SYSTEM, normNnzRowC );
         }
         assert( static_cast<int>(row_support_hashtable.size()) <= mA + norm_Cmat->m);

         // Second Hashing: Per bucket, do Second Hashing:
         row_coefficients_hashtable.clear();
         for( size_t i = 0; i < row_support_hashtable.bucket_count(); ++i )
         {
            // skip bins with less than 2 elements:
            if( row_support_hashtable.bucket_size(i) < 2 )
               continue;
            // insert elements from first Hash-bin into the second Hash-table:
            for( boost::unordered_set<rowlib::rowWithColInd>::local_iterator it =
                  row_support_hashtable.begin(i); it != row_support_hashtable.end(i); ++it )
            {
               row_coefficients_hashtable.emplace(it->id, it->offset_nA, it->lengthA, it->colIndicesA, it->norm_entriesA,
                     it->lengthB, it->colIndicesB, it->norm_entriesB);
            }

            // Compare the rows in the final (from second hash) bin:
            compareRowsInCoeffHashTable(nRowElims, node);

            row_coefficients_hashtable.clear();
         }

         deleteNormalizedPointers( node );
      }
   }

   presData.allreduceLinkingVarBounds();
   presData.allreduceAndApplyLinkingRowActivities();
   presData.allreduceAndApplyNnzChanges();
   presData.allreduceAndApplyBoundChanges();

   row_support_hashtable.clear();
   row_coefficients_hashtable.clear();

   // for the A_0 and C_0 blocks:
   setNormalizedPointers(-1);
   // First Hashing: Fill 'row_support_hashtable':
   if( norm_Amat )
      insertRowsIntoHashtable( row_support_hashtable, norm_Amat, NULL, EQUALITY_SYSTEM, normNnzRowA );
   assert( (int)row_support_hashtable.size() <= mA );
   if( norm_Cmat )
      insertRowsIntoHashtable( row_support_hashtable, norm_Cmat, NULL, INEQUALITY_SYSTEM, normNnzRowC );
   assert( (int)row_support_hashtable.size() <= mA + norm_Cmat->m);
   // Second Hashing: Per bucket, do Second Hashing:
   for( size_t i = 0; i < row_support_hashtable.bucket_count(); ++i )
   {
      if( row_support_hashtable.bucket_size(i) < 2 )
         continue;
      // insert elements from first Hash-bin into the second Hash-table:
      for( boost::unordered_set<rowlib::rowWithColInd>::local_iterator it =
            row_support_hashtable.begin(i); it != row_support_hashtable.end(i); ++it )
      {
         row_coefficients_hashtable.emplace(it->id, it->offset_nA, it->lengthA, it->colIndicesA, it->norm_entriesA,
               it->lengthB, it->colIndicesB, it->norm_entriesB);
      }

      // Compare the rows in the final (from second hash) bin:
      compareRowsInCoeffHashTable(nRowElims, -1);
      row_coefficients_hashtable.clear();
   }

   deleteNormalizedPointers(-1);

   // todo : not necessary?
   presData.allreduceLinkingVarBounds();
   presData.allreduceAndApplyLinkingRowActivities();
   presData.allreduceAndApplyNnzChanges();
   presData.allreduceAndApplyBoundChanges();
   presData.allreduceAndApplyObjVecChanges();

   // TODO:add detection for linking constraints


   synchronize(nRowElims);
   if( my_rank == 0 )
      std::cout << "Removed " << nRowElims << " Rows in Parallel Row Presolving." << std::endl;

#ifndef NDEBUG
   if( my_rank == 0 )
      std::cout << "--- After parallel row presolving:" << std::endl;
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());
}

/** If it is no dummy child, sets normalized pointers:
 * For the matrix blocks A, B, C and D and their transposed matrices and for the lhs/rhs.
 * Sets the pointers to currNnzRow, currNnzRowC, currNnzColChild, currRedColParent.
 * Sets mA and nA correctly.
 */
// todo does not set Bl mat - maybe not needed
void StochPresolverParallelRows::setNormalizedPointersMatrices(int node)
{
   assert(-1 <= node && node <= nChildren);

   const StochGenMatrix& matrixA = dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().A));
   const StochGenMatrix& matrixC = dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C));

   if(node == -1)
   {
      /* EQUALITY_SYSTEM */
      norm_Amat = NULL;
      norm_AmatTrans = NULL;

      norm_Bmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.Bmat)->getStorageDynamicRef());
      norm_BmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.Bmat)->getStorageDynamicTransposedRef());

      /* INEQUALITY_SYSTEM */
      norm_Cmat = NULL;
      norm_CmatTrans = NULL;

      norm_Dmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamicRef());
      norm_DmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamicTransposedRef());
   }
   else
   {
      /* EQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, EQUALITY_SYSTEM ) )
      {
         norm_Amat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Amat)->getStorageDynamicRef());
         norm_AmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Amat)->getStorageDynamicTransposedRef());

         norm_Bmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Bmat)->getStorageDynamicRef());
         norm_BmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Bmat)->getStorageDynamicTransposedRef());
      }
      else
         norm_Amat = norm_AmatTrans = norm_Bmat = norm_BmatTrans = NULL;

      /* INEQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM ) )
      {
         norm_Cmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Amat)->getStorageDynamicRef());
         norm_CmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Amat)->getStorageDynamicTransposedRef());

         norm_Dmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Bmat)->getStorageDynamicRef());
         norm_DmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Bmat)->getStorageDynamicTransposedRef());
      }
      else
         norm_Cmat = norm_CmatTrans = norm_Dmat = norm_DmatTrans = NULL;
   }
}

void StochPresolverParallelRows::setNormalizedPointersMatrixBounds(int node)
{
   assert(-1 <= node && node <= nChildren);

   if(node == -1)
   {
      norm_b = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bA)).vec->clone());

      norm_cupp = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bu)).vec->clone());
      norm_icupp = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().icupp)).vec->clone());
      norm_clow = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bl)).vec->clone());
      norm_iclow = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().iclow)).vec->clone());
   }
   else
   {
      /* EQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, EQUALITY_SYSTEM ) )
      {
         norm_b = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bA)).children[node]->vec->clone());
      }
      else
         norm_b = NULL;

      /* INEQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM ) )
      {
         norm_cupp = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bu)).children[node]->vec->clone());
         norm_clow = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bl)).children[node]->vec->clone());
         norm_icupp = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().icupp)).children[node]->vec->clone());
         norm_iclow = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().iclow)).children[node]->vec->clone());
      }
      else
         norm_cupp = norm_icupp = norm_clow = norm_iclow = NULL;
   }
}

// TODO : does not yet set any pointers for linking constraints of the other system - necessary?
/* sets an extended set of pointers for the current node*/
void StochPresolverParallelRows::updateExtendedPointersForCurrentNode(int node)
{
   assert(-1 <= node && node < nChildren);
   assert(!presData.nodeIsDummy(node, EQUALITY_SYSTEM) || !presData.nodeIsDummy(node, INEQUALITY_SYSTEM));

   if(node == -1)
   {
      updatePointersForCurrentNode(-1, EQUALITY_SYSTEM);

      /* INEQUALITY_SYSTEM */
      currCmat = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).Bmat)->getStorageDynamic();
      currCmatTrans = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).Bmat)->getStorageDynamicTransposed();

      currDmat = NULL;
      currDmatTrans = NULL;

      currIneqRhs = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bu)).vec);
      currIneqLhs = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bl)).vec);
      currIcupp = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().icupp)).vec);
      currIclow = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().iclow)).vec);

      currNnzRowC = dynamic_cast<const SimpleVectorBase<int>*>(presData.getNnzsRowC().vec);
   }
   else
   {
      /* EQUALITY_SYSTEM */
      if(!presData.nodeIsDummy(node, EQUALITY_SYSTEM))
      {
         updatePointersForCurrentNode(node, EQUALITY_SYSTEM);
      }
      else
      {
         updatePointersForCurrentNode(node, INEQUALITY_SYSTEM);
         currAmat = currAmatTrans = currBmat = currBmatTrans = currBlmat = currBlmatTrans = NULL;
         currNnzRow = currNnzRowLink = NULL;
         currEqRhs = currEqRhsLink = NULL;
      }

      /* INEQUALITY_SYSTEM */
      if(!presData.nodeIsDummy(-1, INEQUALITY_SYSTEM))
      {
         currCmat = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).children[node]->Amat)->getStorageDynamic();
         currCmatTrans = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).children[node]->Amat)->getStorageDynamicTransposed();

         currDmat = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).children[node]->Bmat)->getStorageDynamic();
         currDmatTrans = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).children[node]->Bmat)->getStorageDynamicTransposed();

         currIneqRhs = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bu)).children[node]->vec);
         currIneqLhs = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bl)).children[node]->vec);
         currIcupp = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().icupp)).children[node]->vec);
         currIclow = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().iclow)).children[node]->vec);

         currNnzRowC = dynamic_cast<const SimpleVectorBase<int>*>(presData.getNnzsRowC().children[node]->vec);
      }
      else
      {
         currCmat = currCmatTrans = currDmat = currDmatTrans = NULL;
         currNnzRowC = NULL;
      }
   }

}

void StochPresolverParallelRows::setNormalizedNormFactors(int node)
{
   assert(-1 <= node && node <= nChildren);

   if(node == -1)
   {
      norm_factorA = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bA)).vec->clone());
      norm_factorA->setToZero();
      norm_factorC = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bu)).vec->clone());
      norm_factorC->setToZero();
   }
   else
   {
      /* EQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, EQUALITY_SYSTEM ) )
      {
         norm_factorA = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bA)).children[node]->vec->clone());
         norm_factorA->setToZero();
      }
      else
         norm_factorA = NULL;

      /* INEQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM ) )
      {
         norm_factorC = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bu)).children[node]->vec->clone());
         norm_factorC->setToZero();
      }
      else
         norm_factorC = NULL;
   }
}

void StochPresolverParallelRows::setNormalizedSingletonFlags(int node)
{
   assert(-1 <= node && node <= nChildren);

   singletonCoeffsColParent = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*presData.getPresProb().g).vec->clone());
   singletonCoeffsColParent->setToZero();

   if(node == -1)
   {
      rowContainsSingletonVariableA = new SimpleVectorBase<int>(presData.getNnzsRowA().vec->length());
      rowContainsSingletonVariableA->setToConstant( -1.0 );
      rowContainsSingletonVariableC = new SimpleVectorBase<int>(presData.getNnzsRowC().vec->length());
      rowContainsSingletonVariableC->setToConstant( -1.0 );

      singletonCoeffsColChild = NULL;
   }
   else
   {
      singletonCoeffsColChild = dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*presData.getPresProb().g).children[node]->vec->clone());
      singletonCoeffsColChild->setToZero();

      /* EQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, EQUALITY_SYSTEM ) )
      {
         rowContainsSingletonVariableA = new SimpleVectorBase<int>(presData.getNnzsRowA().children[node]->vec->length());
         rowContainsSingletonVariableA->setToConstant( -1.0 );
      }
      else
         rowContainsSingletonVariableA = NULL;

      /* INEQUALITY_SYSTEM */
      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM ) )
      {
         rowContainsSingletonVariableC = new SimpleVectorBase<int>(presData.getNnzsRowC().children[node]->vec->length());
         rowContainsSingletonVariableC->setToConstant( -1.0 );
      }
      else
         rowContainsSingletonVariableC = NULL;
   }
}

void StochPresolverParallelRows::setNormalizedReductionPointers(int node)
{
   assert(-1 <= node && node <= nChildren);

   normNnzColParent = dynamic_cast<SimpleVectorBase<int>*>(presData.getNnzsCol().vec->cloneFull());
   normNnzColChild = (node == -1) ? NULL : 
      dynamic_cast<SimpleVectorBase<int>*>(presData.getNnzsCol().children[node]->vec->cloneFull());

   if(node == -1)
   {
      normNnzRowA = dynamic_cast<SimpleVectorBase<int>*>(currNnzRow->cloneFull());
      normNnzRowC = dynamic_cast<SimpleVectorBase<int>*>(presData.getNnzsRowC().vec->cloneFull());
   }
   else
   {
      /* EQUALITY_SYSTEM */
      normNnzRowA = (!presData.nodeIsDummy( node, EQUALITY_SYSTEM )) ? 
         dynamic_cast<SimpleVectorBase<int>*>(presData.getNnzsRowA().children[node]->vec->cloneFull()) : NULL;
      /* INEQUALITY_SYSTEM */
      normNnzRowC =(!presData.nodeIsDummy( node, INEQUALITY_SYSTEM )) ? 
         dynamic_cast<SimpleVectorBase<int>*>(presData.getNnzsRowC().children[node]->vec->cloneFull()) : NULL;
   }
}

void StochPresolverParallelRows::setNormalizedPointers(int node)
{
   assert( !presData.nodeIsDummy(node, EQUALITY_SYSTEM) || !presData.nodeIsDummy(node, INEQUALITY_SYSTEM) );
   assert(-1 <= node && node <= nChildren );

   updateExtendedPointersForCurrentNode(node);

   /* set normalized matrix pointers for A B (not Bl) */
   setNormalizedPointersMatrices(node);

   /* set normalized lhs rhs for equations */
   setNormalizedPointersMatrixBounds(node);

   /* set up pointers for normalization factors */
   setNormalizedNormFactors(node);

   /* set up singleton flag pointers */
   setNormalizedSingletonFlags(node);

   /* set reduction pointers columns and rows */
   setNormalizedReductionPointers(node);

   /* set mA, nA */ // todo what is this for?
   mA = (norm_Amat) ? norm_Amat->m : 0;
   nA = (norm_Amat) ? norm_Amat->n : 0;

   /* remove singleton columns before normalization */
   removeSingletonVars();

   /* normalization of all rows */
   if( !presData.nodeIsDummy(node, EQUALITY_SYSTEM) )
   {
      assert(norm_Bmat);
      if( node != -1)
         assert(norm_Amat);
      assert(norm_b);
      normalizeBlocksRowwise( EQUALITY_SYSTEM, norm_Amat, norm_Bmat, norm_b, NULL, NULL, NULL);
   }
   if( !presData.nodeIsDummy(node, INEQUALITY_SYSTEM) )
   {
      assert(norm_Dmat);
      if(node != -1)
         assert(norm_Cmat)
      assert(norm_cupp);
      assert(norm_clow);
      assert(norm_icupp);
      assert(norm_iclow);
      normalizeBlocksRowwise( INEQUALITY_SYSTEM, norm_Cmat, norm_Dmat, norm_cupp, norm_clow, norm_icupp, norm_iclow);
   }

   /* asserts */
   assert( norm_Bmat || norm_Dmat);

   if(node != -1 && norm_Amat && norm_Cmat)
   {
      assert( norm_Amat->n == norm_Cmat->n );
      assert( norm_Bmat->n == norm_Dmat->n );
      assert( norm_Amat->m == norm_Bmat->m );
      assert( norm_Cmat->m == norm_Dmat->m );
   }
}

// todo check once
void StochPresolverParallelRows::deleteNormalizedPointers(int node)
{
   delete normNnzColParent;
   delete singletonCoeffsColParent;

   if( node == -1 ) // Case at root
   {
      assert( norm_Bmat && norm_b );
      delete norm_Bmat;
      delete norm_BmatTrans;
      delete norm_b;
      assert( norm_Dmat && norm_cupp && norm_clow && norm_icupp && norm_iclow );
      delete norm_Dmat;
      delete norm_DmatTrans;
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

   /* node != -1 */
   bool childExists = false;
   if( !presData.nodeIsDummy(node, EQUALITY_SYSTEM) )
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

   if( !presData.nodeIsDummy(node, INEQUALITY_SYSTEM) )
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
   assert(normNnzColParent);

   for( int col = 0; col < normNnzColChild->n; col++ )
   {
      if( (*normNnzColChild)[col] == 1 )
      {
         // check if the singleton column is part of the current b_mat/d_mat
         // else, the singleton entry is in one of the other B_i or D_i blocks
         if( norm_BmatTrans && (norm_BmatTrans->rowptr[col].start + 1 == norm_BmatTrans->rowptr[col].end) )
         {
            removeEntry(col, *rowContainsSingletonVariableA, *norm_Bmat, *norm_BmatTrans,
                  *normNnzRowA, *normNnzColChild, (normNnzColChild == NULL) ? true : false);
         }
         else if(norm_DmatTrans && (norm_DmatTrans->rowptr[col].start + 1 == norm_DmatTrans->rowptr[col].end) )
         {
            removeEntry(col, *rowContainsSingletonVariableC, *norm_Dmat, *norm_DmatTrans,
                  *normNnzRowC, *normNnzColChild, (normNnzColChild == NULL) ? true : false);
         }
      }
   }

   // for the child block Bmat and Dmat:
   // if there is an a_mat == we are not in the root node
   if( normNnzColChild )
   {
   for( int col = 0; col < normNnzColParent->n; col++ )
      {
         if( (*normNnzColParent)[col] == 1.0 )
         {
            // check if the singleton column is part of the current a_mat/c_mat
            // else, the singleton entry is in one of the other A_i or C_i blocks
            if( norm_AmatTrans && (norm_AmatTrans->rowptr[col].start + 1 == norm_AmatTrans->rowptr[col].end) )
            {
               removeEntry(col, *rowContainsSingletonVariableA, *norm_Amat, *norm_AmatTrans,
                     *normNnzRowA, *normNnzColParent, true);
            }
            else if(norm_CmatTrans && (norm_CmatTrans->rowptr[col].start +1 == norm_CmatTrans->rowptr[col].end) )
            {
               removeEntry(col, *rowContainsSingletonVariableC, *norm_Cmat, *norm_CmatTrans,
                     *normNnzRowC, *normNnzColParent, true);
            }
         }
      }
   }
}

/** Removes a singleton entry in column colIdx in matrix and matrix_trans and adapts
 * the nnz vectors nnzRow and nnzColChild accordingly. Sets the entries in
 * rowContainsSingletonVar to the corresponding column index in which the singleton entry occurs.
 */
void StochPresolverParallelRows::removeEntry(int colIdx, SimpleVectorBase<int>& rowContainsSingletonVar,
      SparseStorageDynamic& matrix, SparseStorageDynamic& matrixTrans, SimpleVectorBase<int>& nnzRow, SimpleVectorBase<int>& nnzCol,
      bool parent)
{
   assert( 0 <= colIdx && colIdx < matrixTrans.m );
   assert( matrixTrans.rowptr[colIdx].start + 1 == matrixTrans.rowptr[colIdx].end);
   assert( nnzRow.n == matrix.m );
   assert( matrix.n == nnzCol.n );
   assert( nnzCol.n == matrixTrans.m );
   assert( nnzCol[colIdx] == 1 );

   // First, find indices of the singleton entry:
   const int k = matrixTrans.rowptr[colIdx].start;
   const int rowIdx = matrixTrans.jcolM[k];
   assert( rowIdx < nnzRow.n );

   // check if there are no more than one singleton entry in this row
   // if so the row is neither parallel nor nearly parallel to any other row
   if( 0 <= rowContainsSingletonVar[rowIdx] )
   {
      rowContainsSingletonVar[rowIdx] = -2;
      return;
   }

   // store the colIdx in rowContainsSingletonVar
   // (possibly add offset to colIdx so that Amat and Bmat are distinct)
   if( parent )
      rowContainsSingletonVar[rowIdx] = colIdx;
   else // if( !parent )
      rowContainsSingletonVar[rowIdx] = colIdx + nA;

   // Second, remove the entry from norm_Bmat and norm_BmatTrans:
   const double coeff = removeEntryInDynamicStorage(matrix, rowIdx, colIdx);

   matrixTrans.rowptr[colIdx].end--;
   nnzRow[rowIdx]--;
   nnzCol[colIdx] = 0;

   if( parent )
      (*singletonCoeffsColParent)[colIdx] = coeff;
   else
      (*singletonCoeffsColChild)[colIdx] = coeff;
}

/** removes row col from dynamic storage */
double StochPresolverParallelRows::removeEntryInDynamicStorage(SparseStorageDynamic& storage, int row, int col) const
{
   int i = -1;
   int end = storage.rowptr[row].end;
   int start = storage.rowptr[row].start;

   for( i = start; i < end; i++)
   {
      if( storage.jcolM[i] == col )
         break;
   }
   assert( storage.jcolM[i] == col);
   const double coeff = storage.M[i];
   std::swap(storage.M[i], storage.M[end-1]);

   std::swap(storage.jcolM[i], storage.jcolM[end-1]);
   storage.rowptr[row].end--;
   return coeff;
}


// todo there seems to be no numerical threshold for the normalization below .. this should be fixed - rows with fairly different coefficients could be regarded equal
/// cupp can be either the rhs for the equality system or upper bounds for inequalities
void StochPresolverParallelRows::normalizeBlocksRowwise( SystemType system_type,
      SparseStorageDynamic* a_mat, SparseStorageDynamic* b_mat,
      SimpleVector* cupp, SimpleVector* clow, SimpleVector* icupp, SimpleVector* iclow) const
{
   assert(b_mat);
   assert(cupp);
   assert( b_mat->m == cupp->n );
   if( a_mat )
      assert( a_mat->m == b_mat->m);

   if( system_type == INEQUALITY_SYSTEM )
   {
      assert( clow && iclow && icupp);
      assert( clow->n == cupp->n && iclow->n == clow->n && iclow->n == clow->n );
   }


   int n_rows = b_mat->m;

   /// for every row find the max value and normalize by that
   for( int row = 0; row < n_rows; row++)
   {
      double absmax = 0.0;
      bool negate_row = false;

      const int row_B_start = b_mat->rowptr[row].start;
      const int row_B_end = b_mat->rowptr[row].end;
      if( row_B_start < row_B_end )
      {
         if( b_mat->M[row_B_start] < 0)
            negate_row = true;

         for(int k = row_B_start; k < row_B_end; k++)
         {
            if( absmax < std::fabs(b_mat->M[k]) )
               absmax = std::fabs(b_mat->M[k]);
         }
      }

      int row_A_start = 0;
      int row_A_end = 0;
      if( A_mat )
      {
         row_A_start = a_mat->rowptr[row].start;
         row_A_end = a_mat->rowptr[row].end;
         if( row_A_start < row_A_end )
         {
            // todo: this seems wrong to me - I think here we can only decide about negation if the row in a_mat is empty EDIT: fixed it?
            //if( row_B_start == row_B_end && a_mat->M[row_A_start] < 0)
            if( a_mat->M[row_A_start] < 0)
               negate_row = true;
            for(int k = row_A_start; k < row_A_end; k++)
            {
               if( absmax < std::fabs(a_mat->M[k]) )
                  absmax = std::fabs(a_mat->M[k]);
            }
         }
      }

      // normalize the row by dividing all entries by abs_max and possibly by -1 if negate_row.
      if(negate_row)
         absmax = absmax * -1.0;
      
      for(int k = row_B_start; k < row_B_end; k++)
         b_mat->M[k] /= absmax;

      for(int k = row_A_start; k < row_A_end; k++)
         a_mat->M[k] /= absmax;

      if( system_type == EQUALITY_SYSTEM )
      {
         (*cupp)[row] /= absmax;
         (*norm_factorA)[row] = absmax;
      }
      else
      {
         if( !PIPSisZero((*iclow)[row]) )
            (*clow)[row] /= absmax;
         if( !PIPSisZero((*icupp)[row]) )
            (*cupp)[row] /= absmax;
         // if multiplied by a negative value, lhs and rhs have to be swapped:
         if(negate_row)
         {
            std::swap( (*clow)[row], (*cupp)[row] );
            std::swap( (*iclow)[row], (*icupp)[row] );
         }
         (*norm_factorC)[row] = absmax;
      }
   }
}

/** Inserts all non-empty rows into the given unordered set 'rows'.
 * Ablock and Bblock are supposed to be the two matrix blocks of a child.
 * If at root, Bblock should be NULL.
 */
// todo : I think this is wrong or at least not complete - in theory we should sort the rows first (according to the colindices)
void StochPresolverParallelRows::insertRowsIntoHashtable( boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> > &rows,
      SparseStorageDynamic* a_mat, SparseStorageDynamic* b_mat, SystemType system_type, SimpleVectorBase<int>* nnzRow )
{
   assert(b_mat);
   if(a_mat)
     assert( a_mat->m == b_mat->m );
   if( system_type == EQUALITY_SYSTEM )
      assert( mA == b_mat->m );

   for(int row = 0; row < b_mat->m; row++)
   {
      // ignore rows containing more than one singleton entry:
      if( system_type == EQUALITY_SYSTEM && (*rowContainsSingletonVariableA)[row] == -2 )
         continue;
      if( system_type == INEQUALITY_SYSTEM && (*rowContainsSingletonVariableC)[row] == -2 )
         continue;

      // calculate rowId including possible offset (for Inequality rows):
      int rowId = row;
      if( (*nnzRow)[rowId] == 0 )
         continue;
      if( system_type == INEQUALITY_SYSTEM )
         rowId += mA;

      // calculate rowlength of a_mat and b_mat
      const int row_B_start =  b_mat->rowptr[row].start;
      const int row_B_length = b_mat->rowptr[row].end - row_B_start;
      
      if( a_mat )
      {
         const int row_A_start = a_mat->rowptr[row].start;
         const int row_A_length = a_mat->rowptr[row].end - row_A_start;

         if( row_B_length == 0 && row_A_length == 0 )
            continue;
         // colIndices and normalized entries are set as pointers to the original data.
         // create and insert the new element:
         rows.emplace(rowId, nA, row_B_length, &(b_mat->jcolM[row_B_start]), &(b_mat->M[row_B_start]),
               row_A_length, &(a_mat->jcolM[row_A_start]), &(a_mat->M[row_A_start]));
      }
      else
      {
         if(row_B_length == 0)
            continue;

         rows.emplace(rowId, nA, row_B_length, &(b_mat->jcolM[row_B_start]),
               &(b_mat->M[row_B_start]), 0, NULL, NULL);
      }
   }
}

/*
 * Per bucket in row_coefficients_hashtable, compare the containing rows and check if they are parallel.
 * If so, consider the different possible cases. If a row can be removed, the action is applied
 * to the original matrices (not the normalized copies).
 * @param it represents the child number (or -1 for parent blocks)
 */
void StochPresolverParallelRows::compareRowsInCoeffHashTable(int& nRowElims, int node)
{
   if( row_coefficients_hashtable.empty() )
      return;

   for (size_t i = 0; i < row_coefficients_hashtable.bucket_count(); ++i)
   {
      for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator it1 = row_coefficients_hashtable.begin(i);
            it1 != row_coefficients_hashtable.end(i); ++it1)
      {
         // either pairwise comparison OR lexicographical sorting and then compare only neighbors.
         // Here: parwise comparison: // todo make lexicographical
         for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator it2 = it1;
                                    it2 != row_coefficients_hashtable.end(i); ++it2)
         {
            // When two parallel rows are found, check if they are both =, both <=, or = and <=
            if( checkRowsAreParallel( *it1, *it2) )
            {
               if( it1->id < mA && it2->id < mA )
               {
                  int id1 = it1->id;
                  int id2 = it2->id;
                  // Case both constraints are equalities:
                  // check if one of them contains a singleton variable:
                  if( (*rowContainsSingletonVariableA)[id1] != -1
                        || (*rowContainsSingletonVariableA)[id2] != -1 )
                  {
                     // nearly parallel case 1:
                     // make sure that a_r2 != 0, otherwise switch ids.
                     if( (*rowContainsSingletonVariableA)[id2] != -1 )
                        std::swap(id1, id2);

                     assert( (*rowContainsSingletonVariableA)[id2] != -1 );
                     // if row id1 was already deleted, do not continue the procedure:
                     if( (id1 < mA && (*currNnzRow)[id1]) == 0 )
                        continue;

                     // case two is basically case one
                     doNearlyParallelRowCase1(id1, id2, node);
                  }
                  else
                  {
                     if( !PIPSisEQ( (*norm_b)[id1], (*norm_b)[id2]) )
                        abortInfeasible(MPI_COMM_WORLD, "Found parallel equality rows with non-compatible right hand sides",
                           "StochPresolverParallelRows.C", "compareRowsInCoeffHashTable"); 
                     // delete row2 in the original system:
                     presData.removeRedundantRow( EQUALITY_SYSTEM, node, id2, false);
                  }
               }
               else if( it1->id >= mA && it2->id >= mA )
               {
                  // Case both constraints are inequalities
                  const int id1 = it1->id - mA;
                  const int id2 = it2->id - mA;

                  if( (*rowContainsSingletonVariableC)[id1] != -1
                        || (*rowContainsSingletonVariableC)[id2] != -1 )
                  {
                     // nearly parallel case 3:
                     if( (*rowContainsSingletonVariableC)[id1] != -1
                        && (*rowContainsSingletonVariableC)[id2] != -1 )
                     {
                        doNearlyParallelRowCase3(id1, id2, node);
                     }
                  }
                  else
                  {
                     // tighten bounds in original and normalized system:
                     tightenOriginalBoundsOfRow1( INEQUALITY_SYSTEM, node, id1, id2 );

                     // delete row2 in the original system:
                     presData.removeRedundantRow( INEQUALITY_SYSTEM, node, id2, false);
                  }
               }
               else
               {  // Case one constraint is an equality, one an inequality
                  int id1 = it1->id;
                  int id2 = it2->id;
                  if( id1 >= mA )   // swap ids so that id2 is the inequality constraint.
                     std::swap(id1, id2);

                  const int ineqRowId = id2 - mA;

                  if( (*rowContainsSingletonVariableA)[id1] == -1
                        && (*rowContainsSingletonVariableC)[ineqRowId] == -1 )
                  {
                     //PIPSdebugMessage("Really Parallel Rows, case 2. \n");
                     // check for infeasibility:
                     if( !PIPSisZero((*norm_iclow)[ineqRowId])
                           && PIPSisLT( (*norm_b)[id1], (*norm_clow)[ineqRowId]) )
                        abortInfeasible(MPI_COMM_WORLD, "Found parallel inequality and equality rows where rhs/lhs do not match",
                           "StochPresolverParallelRows.C", "compareRowsInCoeffHashTable"); 
                     if( !PIPSisZero((*norm_icupp)[ineqRowId])
                           && PIPSisLT((*norm_cupp)[ineqRowId], (*norm_b)[id1]) )
                        abortInfeasible(MPI_COMM_WORLD, "Found parallel inequality and equality rows where rhs/lhs do not match",
                           "StochPresolverParallelRows.C", "compareRowsInCoeffHashTable"); 
                     // remove inequality row.
                  }
                  else if( (*rowContainsSingletonVariableA)[id1] != -1
                        && (*rowContainsSingletonVariableC)[ineqRowId] == -1 )
                  {
                     // nearly parallel case 2:
                     PIPSdebugMessage("Nearly Parallel Rows, case 2. \n");
                     // compute the new variable bound for x_id1:
                     int singleColIdx = (*rowContainsSingletonVariableA)[id1];
                     double coeff_singleton = getSingletonCoefficient(singleColIdx);

                     double newxlow = - std::numeric_limits<double>::max();
                     double newxupp = std::numeric_limits<double>::max();

                     // check if f_q * a_q is positive or negative:
                     double faq = (*norm_factorA)[id1] * coeff_singleton;
                     if( faq > 0 )
                     {
                        if( !PIPSisZero((*norm_iclow)[ineqRowId]) )
                           newxupp = ( (*norm_b)[id1] - (*norm_clow)[ineqRowId] ) * (*norm_factorA)[id1] / coeff_singleton;
                        if( !PIPSisZero((*norm_icupp)[ineqRowId]) )
                           newxlow = ( (*norm_b)[id1] - (*norm_cupp)[ineqRowId] ) * (*norm_factorA)[id1] / coeff_singleton;
                     }
                     else  //if( fq*aq < 0 )
                     {
                        if( !PIPSisZero((*norm_iclow)[ineqRowId]) )
                           newxlow = ( (*norm_b)[id1] - (*norm_clow)[ineqRowId] ) * (*norm_factorA)[id1] / coeff_singleton;
                        if( !PIPSisZero((*norm_icupp)[ineqRowId]) )
                           newxupp = ( (*norm_b)[id1] - (*norm_cupp)[ineqRowId] ) * (*norm_factorA)[id1] / coeff_singleton;
                     }

                     BlockType block_type = A_MAT;
                     if(singleColIdx > nA || node == -1)
                        block_type = B_MAT;

                     presData.rowPropagatedBounds(EQUALITY_SYSTEM, node, block_type, id1, singleColIdx, newxupp, newxlow);
                     // remove inequality row -> after block
                  }
                  else
                     continue;   // no action if the inequality contains a singleton entry.

                  // in both the parallel row case and in the nearly parallel row case,
                  // delete the inequality constraint id2 (which could be either it1 or it2 now):
                  presData.removeRedundantRow(INEQUALITY_SYSTEM, node, ineqRowId, false);

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
   for( int i = 0; i < row1.lengthA; i++)
   {
      if( row1.colIndicesA[i] != row2.colIndicesA[i] )
         return false;
      if( !PIPSisEQ(row1.norm_entriesA[i], row2.norm_entriesA[i], tol_compare_double) )
         return false;
   }
   for( int i = 0; i < row1.lengthB; i++)
   {
      if( row1.colIndicesB[i] != row2.colIndicesB[i] )
         return false;
      if( !PIPSisEQ(row1.norm_entriesB[i], row2.norm_entriesB[i], tol_compare_double) )
         return false;
   }
   return true;
}

/**
 * Tightens the original lower and upper bounds of the first row, given the lower
 * and upper bounds of the second row. The normalized bounds are compared and the
 * normalizing factor of row1 is used to determine which bound can be tightened to
 * which value. The normalized bounds of row1 are also updated.
 * Assumes that both rows are inequality constraints.
 */
void StochPresolverParallelRows::tightenOriginalBoundsOfRow1(SystemType system_type, int node, int rowId1, int rowId2)
{
   assert( system_type == INEQUALITY_SYSTEM);
   assert( currCmat );
   assert( rowId1 >= 0 && rowId2 >= 0 );
   assert( rowId1 < currCmat->m && rowId2 < currCmat->m );
   assert( norm_factorC && norm_factorC->n == currCmat->m );

   // the normalizing factor of the first inequality:
   double factor = (*norm_factorC)[rowId1];

   // Set norm_low/norm_upp to the bounds of the second inequality:
   double norm_low_row2 = -std::numeric_limits<double>::infinity();
   double norm_upp_row2 = std::numeric_limits<double>::infinity();

   if( !PIPSisZero((*norm_iclow)[rowId2]) )
      norm_low_row2 = (*norm_clow)[rowId2];
   if( !PIPSisZero((*norm_icupp)[rowId2]) )
      norm_upp_row2 = (*norm_cupp)[rowId2];

   // test for infeasibility: // todo
   if( ( !PIPSisZero((*norm_iclow)[rowId1]) && PIPSisLT( norm_upp_row2, (*norm_clow)[rowId1]) )
         || ( !PIPSisZero((*norm_icupp)[rowId1]) && PIPSisLT( (*norm_cupp)[rowId1], norm_low_row2) ) )
   {
      abortInfeasible(MPI_COMM_WORLD, "Found incompatible row rhs/lhs", "StochPresolverParallelRows.C", "tightenOriginalBoundsOfRow1");
   }

   double new_lhs = -std::numeric_limits<double>::infinity();
   double new_rhs = std::numeric_limits<double>::infinity();

   // todo numeric safeguards // todo check if inifty * factor == infty?
   if( ( !PIPSisZero((*norm_iclow)[rowId1]) && PIPSisLT( (*norm_clow)[rowId1], norm_low_row2) ) ||
      ( PIPSisZero((*norm_iclow)[rowId1]) && norm_low_row2 > -std::numeric_limits<double>::infinity() ) )
   {
      (*norm_clow)[rowId1] = /*std::max(*/norm_low_row2;//, (*norm_clow)[rowId1]);
      (*norm_iclow)[rowId1] = 1.0;


      ( PIPSisLT( 0.0, factor) ) ? new_lhs = factor * norm_low_row2 : new_rhs = factor * norm_low_row2;
   }

   if( ( !PIPSisZero((*norm_icupp)[rowId1]) && PIPSisLT(norm_upp_row2, (*norm_cupp)[rowId1]) )
         || ( PIPSisZero((*norm_icupp)[rowId1]) && norm_upp_row2 < std::numeric_limits<double>::infinity() ))
   {
      (*norm_cupp)[rowId1] = /*std::min(*/norm_upp_row2;/*, (*norm_cupp)[rowId1]);*/ // todo!!!
      (*norm_icupp)[rowId1] = 1.0;

      ( PIPSisLT( 0.0, factor) ) ? new_rhs = factor * norm_upp_row2 : new_lhs = factor * norm_upp_row2;
   }

   assert( PIPSisLE( new_lhs, new_rhs) );
   presData.tightenRowBoundsParallelRow(INEQUALITY_SYSTEM, node, rowId1, new_lhs, new_rhs, false);
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

   return ( singleColIdx >= nA ) ? (*singletonCoeffsColChild)[singleColIdx - nA] :
         (*singletonCoeffsColParent)[singleColIdx];
}

void StochPresolverParallelRows::doNearlyParallelRowCase1(int rowId1, int rowId2, int it)
{
   /* assert both rows are equality constraints */
   assert( rowId1 >= 0 && rowId1 < mA );
   assert( rowId2 >= 0 && rowId2 < mA );
   assert( it >= -1 && it < nChildren );

   // calculate t and d:
   const int singleColIdx1 = (*rowContainsSingletonVariableA)[rowId1];
   const int singleColIdx2 = (*rowContainsSingletonVariableA)[rowId2];
   assert( !PIPSisEQ(singleColIdx2, -1) );

   const double coeff_singleton1 = getSingletonCoefficient(singleColIdx1);
   const double coeff_singleton2 = getSingletonCoefficient(singleColIdx2);
   assert( !PIPSisZero(coeff_singleton2) );

   const double s = (*norm_factorA)[rowId1] / (*norm_factorA)[rowId2];
   const double t = coeff_singleton1 / coeff_singleton2 / s;
   const double d = ( (*norm_b)[rowId2] - (*norm_b)[rowId1] )
            * (*norm_factorA)[rowId2] / coeff_singleton2;

   const int col1_idx = (singleColIdx1 < nA) ? singleColIdx1 : (singleColIdx1 - nA);
   const int col2_idx = (singleColIdx2 < nA) ? singleColIdx2 : (singleColIdx2 - nA);
   const int node_var1 = (singleColIdx1 < nA) ? -1 : it;
   const int node_var2 = (singleColIdx2 < nA) ? -1 : it;     
   const SimpleVector* ixlow = (singleColIdx2 < nA) ? currIxlowParent : currIxlowChild;
   const SimpleVector* ixupp = (singleColIdx2 < nA) ? currIxuppParent : currIxuppChild;
   const SimpleVector* xlow = (singleColIdx2 < nA) ? currxlowParent : currxlowChild;
   const SimpleVector* xupp = (singleColIdx2 < nA) ? currxuppParent : currxuppChild;

   if( !PIPSisEQ(singleColIdx1, -1.0) )
   {
      // tighten the bounds of variable x_1:
      double newxlow = -std::numeric_limits<double>::max();
      double newxupp = std::numeric_limits<double>::max();

      // calculate new bounds depending on the sign of t:
      if( PIPSisLT(t, 0) )
      {
         std::swap(ixlow, ixupp);
         std::swap(xlow, xupp);
      }

      if( PIPSisEQ( (*ixlow)[col2_idx], 1.0 ) )
         newxlow = ((*xlow)[col2_idx] - d) / t;
      if( PIPSisEQ( (*ixupp)[col2_idx], 1.0 ) )
         newxupp = ((*xupp)[col2_idx] - d) / t;

      presData.substituteVariableParallelRows(EQUALITY_SYSTEM, it, col1_idx, rowId1, node_var1, col2_idx, rowId2, node_var2, t, d);
      // effectively tighten bounds of variable x_id1:
      presData.rowPropagatedBounds(EQUALITY_SYSTEM, it, B_MAT, rowId1, col1_idx, newxupp, newxlow);
   }
   else
   {
      /* row1 does not have a singleton variable - row2 does */
      /* var2 = d */
      assert( PIPSisZero(t) );
      presData.fixColumn( node_var2, col2_idx, d );
   }
}
/**
 * Executes the Nearly Parallel Row Case 3: Both constraints are inequalities and both contain
 * a singleton variable entry.
 * The row indices rowId1, rowId2 are already de-offset, so they should be in the range [0,Cmat->m).
 */
void StochPresolverParallelRows::doNearlyParallelRowCase3(int rowId1, int rowId2, int it)
{
   assert( rowId1 >= 0 && rowId2 >= 0 );
   assert( rowContainsSingletonVariableC );
   assert( rowId1 < rowContainsSingletonVariableC->n );
   assert( rowId2 < rowContainsSingletonVariableC->n );
   assert( norm_factorC );
   assert( !PIPSisZero((*norm_factorC)[rowId2]) );

   // First, do all the checks to verify if the third case applies and all conditions are met.

   const int singleColIdx1 = (*rowContainsSingletonVariableC)[rowId1];
   const int singleColIdx2 = (*rowContainsSingletonVariableC)[rowId2];
   const int col1_idx = (singleColIdx1 < nA) ? singleColIdx1 : singleColIdx1 - nA;
   const int col2_idx = (singleColIdx2 < nA) ? singleColIdx2 : singleColIdx2 - nA;
   const double coeff_singleton1 = getSingletonCoefficient(singleColIdx1);
   const double coeff_singleton2 = getSingletonCoefficient(singleColIdx2);
   const double s = (*norm_factorC)[rowId1] / (*norm_factorC)[rowId2];

   assert(col1_idx < nA);
   assert(col2_idx < nA);
   
   // check s > 0:
   if( PIPSisLT(s, 0.0) )
      return;
   // check a_q != 0.0, a_r!= 0.0:
   if( PIPSisZero(coeff_singleton1) || PIPSisZero(coeff_singleton2) )
      return;

   // check d_q == s * d_r (here s == 1, so just d_q == d_r) and f_q == f_r
   if(   !PIPSisEQ( (*norm_iclow)[rowId1], (*norm_iclow)[rowId2] )
      || !PIPSisEQ( (*norm_icupp)[rowId1], (*norm_icupp)[rowId2] ) )
      return;

   if( !PIPSisZero( (*norm_iclow)[rowId1] ) && !PIPSisEQ( (*norm_clow)[rowId1], (*norm_clow)[rowId2] ) )
      return;
   if( !PIPSisZero( (*norm_icupp)[rowId1] ) && !PIPSisEQ( (*norm_cupp)[rowId1], (*norm_cupp)[rowId2] ) )
      return;

   // check c_1 * c_2 >= 0:
   const double c1 = (singleColIdx1 < nA) ? (*currgParent)[col1_idx] : (*currgChild)[col1_idx];
   const double c2 = (singleColIdx2 < nA) ? (*currgParent)[col2_idx] : (*currgChild)[col2_idx];
   if( c1 * c2 < 0 )
      return;

   // check lower and upper bounds:
   const double ixlow1 = (singleColIdx1 < nA) ? (*currIxlowParent)[col1_idx] : (*currIxlowChild)[col1_idx];
   const double ixupp1 = (singleColIdx1 < nA) ? (*currIxuppParent)[col1_idx] : (*currIxuppChild)[col1_idx];
   const double ixlow2 = (singleColIdx2 < nA) ? (*currIxlowParent)[col2_idx] : (*currIxlowChild)[col2_idx];
   const double ixupp2 = (singleColIdx2 < nA) ? (*currIxuppParent)[col2_idx] : (*currIxuppChild)[col2_idx];

   const double xlow1 = (singleColIdx1 < nA) ? (*currxlowParent)[col1_idx] : (*currxlowChild)[col1_idx];
   const double xupp1 = (singleColIdx1 < nA) ? (*currxuppParent)[col1_idx] : (*currxuppChild)[col1_idx];
   const double xlow2 = (singleColIdx2 < nA) ? (*currxlowParent)[col2_idx] : (*currxlowChild)[col2_idx];
   const double xupp2 = (singleColIdx2 < nA) ? (*currxuppParent)[col2_idx] : (*currxuppChild)[col2_idx];

   // check aq * l1 = s * ar * l2:
   if( !PIPSisEQ(ixlow1, ixlow2) )
      return;
   if( !PIPSisZero(ixlow1) && !PIPSisZero(ixlow2) &&
         !PIPSisEQ(coeff_singleton1 * xlow1, s * coeff_singleton2 * xlow2) )
      return;

   // check aq * u1 = s * ar * u2:
   if( !PIPSisEQ(ixupp1, ixupp2) )
      return;
   if( !PIPSisZero(ixupp1) && !PIPSisZero(ixupp2) &&
         !PIPSisEQ(coeff_singleton1 * xupp1, s * coeff_singleton2 * xupp2) )
      return;

   // Second, aggregate x_2: adapt objectiveCost(x_1)
   const int node_var1 = (singleColIdx1 < nA) ? -1 : it;
   const int node_var2 = (singleColIdx2 < nA) ? -1 : it;
   const double t = coeff_singleton1 / coeff_singleton2 / s;

   presData.substituteVariableParallelRows(INEQUALITY_SYSTEM, it, col1_idx, rowId1, node_var1, col2_idx, rowId2, node_var2, t, 0.0);
}
