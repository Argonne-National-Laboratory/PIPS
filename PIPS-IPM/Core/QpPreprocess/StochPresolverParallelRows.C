/*
 * StochPresolverParallelRows.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

// todo: should we ever decide to switch to a newer c++ standard - this can be optimized

//#define PIPS_DEBUG
#include "pipsport.h"
#include "StochPresolverParallelRows.h"
#include "pipsport.h"
#include "StochVectorUtilities.h"

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
   setExtendedPointersToNull();
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
      if( !presData.nodeIsDummy(node) )
      {
         /// copy and normalize A_i, B_i, C_i, D_i and b_i, clow_i, cupp_i
         setNormalizedPointers(node);

         row_support_hashtable.clear();

         // Per row, add row to the set 'row_support_hashtable':
         assert(norm_Amat);
         assert(norm_Bmat);
         assert(normNnzRowA);

         insertRowsIntoHashtable( row_support_hashtable, norm_Amat, norm_Bmat, EQUALITY_SYSTEM, normNnzRowA );
         assert( static_cast<int>(row_support_hashtable.size()) <= mA );

         assert(norm_Cmat);
         assert(norm_Dmat);
         assert(normNnzRowC);

         insertRowsIntoHashtable( row_support_hashtable, norm_Cmat, norm_Dmat, INEQUALITY_SYSTEM, normNnzRowC );

         assert( static_cast<int>(row_support_hashtable.size()) <= mA + norm_Cmat->getM() );

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
   assert(norm_Bmat); assert(norm_Dmat);
   // First Hashing: Fill 'row_support_hashtable':
   insertRowsIntoHashtable( row_support_hashtable, nullptr, norm_Bmat, EQUALITY_SYSTEM, normNnzRowA );
   assert( static_cast<int>(row_support_hashtable.size()) <= mA );

   insertRowsIntoHashtable( row_support_hashtable, nullptr, norm_Dmat, INEQUALITY_SYSTEM, normNnzRowC );
   assert( static_cast<int>(row_support_hashtable.size()) <= mA + norm_Dmat->getM());
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


   nRowElims = PIPS_MPIgetMax(nRowElims, MPI_COMM_WORLD);
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
   assert(-1 <= node && node < nChildren);

   const StochGenMatrix& matrixA = dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().A));
   const StochGenMatrix& matrixC = dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C));

   if(node == -1)
   {
      /* EQUALITY_SYSTEM */
      norm_Amat = nullptr;
      norm_AmatTrans = nullptr;

      norm_Bmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.Bmat)->getStorageDynamicRef());
      norm_BmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.Bmat)->getStorageDynamicTransposedRef());

      /* INEQUALITY_SYSTEM */
      norm_Cmat = nullptr;
      norm_CmatTrans = nullptr;

      norm_Dmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamicRef());
      norm_DmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.Bmat)->getStorageDynamicTransposedRef());
   }
   else
   {
      if( !presData.nodeIsDummy( node ) )
      {
         /* EQUALITY_SYSTEM */
         norm_Amat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Amat)->getStorageDynamicRef());
         norm_AmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Amat)->getStorageDynamicTransposedRef());

         norm_Bmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Bmat)->getStorageDynamicRef());
         norm_BmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixA.children[node]->Bmat)->getStorageDynamicTransposedRef());

         /* INEQUALITY_SYSTEM */
         norm_Cmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Amat)->getStorageDynamicRef());
         norm_CmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Amat)->getStorageDynamicTransposedRef());

         norm_Dmat = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Bmat)->getStorageDynamicRef());
         norm_DmatTrans = new SparseStorageDynamic(dynamic_cast<SparseGenMatrix*>(matrixC.children[node]->Bmat)->getStorageDynamicTransposedRef());
      }
      else
      {
         norm_Amat = norm_AmatTrans = norm_Bmat = norm_BmatTrans = nullptr;
         norm_Cmat = norm_CmatTrans = norm_Dmat = norm_DmatTrans = nullptr;
      }

   }
}

void StochPresolverParallelRows::setNormalizedPointersMatrixBounds(int node)
{
   assert(-1 <= node && node < nChildren);
   assert(!presData.nodeIsDummy(node));

   norm_b = dynamic_cast<SimpleVector*>(getSimpleVecFromRowStochVec(*presData.getPresProb().bA, node, false).cloneFull());

   norm_cupp = dynamic_cast<SimpleVector*>(getSimpleVecFromRowStochVec(*presData.getPresProb().bu, node, false).cloneFull());
   norm_icupp = dynamic_cast<SimpleVector*>(getSimpleVecFromRowStochVec(*presData.getPresProb().icupp, node, false).cloneFull());
   norm_clow = dynamic_cast<SimpleVector*>(getSimpleVecFromRowStochVec(*presData.getPresProb().bl, node, false).cloneFull());
   norm_iclow = dynamic_cast<SimpleVector*>(getSimpleVecFromRowStochVec(*presData.getPresProb().iclow, node, false).cloneFull());
}

// TODO : does not yet set any pointers for linking constraints of the other system - necessary?
/* sets an extended set of pointers for the current node*/
void StochPresolverParallelRows::updateExtendedPointersForCurrentNode(int node)
{
   assert(-1 <= node && node < nChildren);
   assert(!presData.nodeIsDummy(node));

   updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

   if(node == -1)
   {

      /* INEQUALITY_SYSTEM */
      currCmat = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).Bmat)->getStorageDynamic();
      currCmatTrans = dynamic_cast<SparseGenMatrix*>(dynamic_cast<const StochGenMatrix&>(*(presData.getPresProb().C)).Bmat)->getStorageDynamicTransposed();

      currDmat = nullptr;
      currDmatTrans = nullptr;

      currIneqRhs = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bu)).vec);
      currIneqLhs = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bl)).vec);
      currIcupp = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().icupp)).vec);
      currIclow = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().iclow)).vec);

      currNnzRowC = dynamic_cast<const SimpleVectorBase<int>*>(presData.getNnzsRowC().vec);
   }
   else
   {

      /* INEQUALITY_SYSTEM */
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

}

void StochPresolverParallelRows::setNormalizedNormFactors(int node)
{
   assert(-1 <= node && node < nChildren);
   assert(!presData.nodeIsDummy(node));

   norm_factorA = dynamic_cast<SimpleVector*>(getSimpleVecFromRowStochVec(*presData.getPresProb().bA, node, false).clone());
   norm_factorA->setToZero();
   norm_factorC = dynamic_cast<SimpleVector*>(getSimpleVecFromRowStochVec(*presData.getPresProb().bu, node, false).clone());
   norm_factorC->setToZero();
}

void StochPresolverParallelRows::setNormalizedSingletonFlags(int node)
{
   assert(-1 <= node && node < nChildren);
   assert(!presData.nodeIsDummy(node));

   singletonCoeffsColParent = dynamic_cast<SimpleVector*>(getSimpleVecFromColStochVec(*presData.getPresProb().g, -1).clone());
   singletonCoeffsColParent->setToZero();

   rowContainsSingletonVariableA = new SimpleVectorBase<int>(getSimpleVecFromRowStochVec(presData.getNnzsRowA(), node, false).length());
   rowContainsSingletonVariableA->setToConstant( -1 );
   rowContainsSingletonVariableC = new SimpleVectorBase<int>(getSimpleVecFromRowStochVec(presData.getNnzsRowC(), node, false).length());
   rowContainsSingletonVariableC->setToConstant( -1 );
   
   if(node == -1)
      singletonCoeffsColChild = nullptr;
   else
   {
      singletonCoeffsColChild = dynamic_cast<SimpleVector*>(getSimpleVecFromColStochVec(*presData.getPresProb().g, node).clone());
      singletonCoeffsColChild->setToZero();
   }
}

void StochPresolverParallelRows::setNormalizedReductionPointers(int node)
{
   assert(-1 <= node && node < nChildren);

   normNnzColParent = dynamic_cast<SimpleVectorBase<int>*>(getSimpleVecFromColStochVec(presData.getNnzsCol(), -1).cloneFull());
   normNnzColChild = (node == -1) ? nullptr : 
      dynamic_cast<SimpleVectorBase<int>*>(getSimpleVecFromColStochVec(presData.getNnzsCol(), node).cloneFull());

   normNnzRowA = dynamic_cast<SimpleVectorBase<int>*>(getSimpleVecFromRowStochVec(presData.getNnzsRowA(), node, false).cloneFull());
   normNnzRowC = dynamic_cast<SimpleVectorBase<int>*>(getSimpleVecFromRowStochVec(presData.getNnzsRowC(), node, false).cloneFull());
}

void StochPresolverParallelRows::setNormalizedPointers(int node)
{
   assert( !presData.nodeIsDummy(node) );
   assert(-1 <= node && node < nChildren );

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

   assert(norm_Bmat);
   /* set mA, nA */ // todo what is this for?
   mA = (norm_Amat) ? norm_Amat->getM() : norm_Bmat->getM();
   nA = (norm_Amat) ? norm_Amat->getN() : norm_Bmat->getN();

   /* remove singleton columns before normalization */
   removeSingletonVars();

   /* normalization of all rows */
   if( !presData.nodeIsDummy(node) )
   {
      assert(norm_Bmat);
      if( node != -1)
         assert(norm_Amat);
      assert(norm_b);
      normalizeBlocksRowwise( EQUALITY_SYSTEM, norm_Amat, norm_Bmat, norm_b, nullptr, nullptr, nullptr);

      assert(norm_Dmat);
      if(node != -1)
         assert(norm_Cmat);
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
      assert( norm_Amat->getN() == norm_Cmat->getN() );
      assert( norm_Bmat->getN() == norm_Dmat->getN() );
      assert( norm_Amat->getM() == norm_Bmat->getM() );
      assert( norm_Cmat->getM() == norm_Dmat->getM() );
   }
}

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
   if( !presData.nodeIsDummy(node) )
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

void StochPresolverParallelRows::setExtendedPointersToNull()
{
   currCmat = nullptr;
   currCmatTrans = nullptr;
   currDmat = nullptr;
   currDmatTrans = nullptr;
   currNnzRowC = nullptr;
   norm_Amat = nullptr;
   norm_Bmat= nullptr;
   norm_Cmat = nullptr;
   norm_Dmat = nullptr;
   norm_b = nullptr;
   norm_cupp = nullptr;
   norm_clow = nullptr;
   norm_icupp = nullptr;
   norm_iclow = nullptr;
   norm_factorA = nullptr;
   rowContainsSingletonVariableA = nullptr;
   norm_factorC = nullptr;
   rowContainsSingletonVariableC = nullptr;
   singletonCoeffsColParent = nullptr;
   singletonCoeffsColChild = nullptr;
   normNnzRowA = nullptr;
   normNnzRowC = nullptr;
   normNnzColParent = nullptr;
   normNnzColChild = nullptr;
   norm_AmatTrans = nullptr;
   norm_BmatTrans = nullptr;
   norm_CmatTrans = nullptr;
   norm_DmatTrans = nullptr;
   currgChild = nullptr;
   currgParent = nullptr;
}

void StochPresolverParallelRows::removeSingletonVars()
{
   assert(normNnzColChild || normNnzColParent);
   const bool at_root_node = (normNnzColChild == nullptr);
   
   SimpleVectorBase<int>* nnzs_bmat = (at_root_node) ? normNnzColParent :
      normNnzColChild;
   
   /* Bmat */
   for( int col = 0; col < nnzs_bmat->n; col++ )
   {
      if( (*nnzs_bmat)[col] == 1 )
      {
         // check if the singleton column is part of the current b_mat/d_mat
         // else, the singleton entry is in one of the other B_i or D_i blocks
         if( norm_BmatTrans && (norm_BmatTrans->getRowPtr(col).start + 1 == norm_BmatTrans->getRowPtr(col).end) )
         {
            removeEntry(col, *rowContainsSingletonVariableA, *norm_Bmat, *norm_BmatTrans,
                  *normNnzRowA, *nnzs_bmat, at_root_node);
         }
         else if(norm_DmatTrans && (norm_DmatTrans->getRowPtr(col).start + 1 == norm_DmatTrans->getRowPtr(col).end) )
         {
            removeEntry(col, *rowContainsSingletonVariableC, *norm_Dmat, *norm_DmatTrans,
                  *normNnzRowC, *nnzs_bmat, at_root_node);
         }
      }
   }

   // for the child block Bmat and Dmat:
   // if there is an a_mat == we are not in the root node
   if( !at_root_node )
   {
   for( int col = 0; col < normNnzColParent->n; col++ )
      {
         if( (*normNnzColParent)[col] == 1 )
         {
            // check if the singleton column is part of the current a_mat/c_mat
            // else, the singleton entry is in one of the other A_i or C_i blocks
            if( norm_AmatTrans && (norm_AmatTrans->getRowPtr(col).start + 1 == norm_AmatTrans->getRowPtr(col).end) )
            {
               removeEntry(col, *rowContainsSingletonVariableA, *norm_Amat, *norm_AmatTrans,
                     *normNnzRowA, *normNnzColParent, true);
            }
            else if(norm_CmatTrans && (norm_CmatTrans->getRowPtr(col).start + 1 == norm_CmatTrans->getRowPtr(col).end) )
            {
               removeEntry(col, *rowContainsSingletonVariableC, *norm_Cmat, *norm_CmatTrans,
                     *normNnzRowC, *normNnzColParent, true);
            }
         }
      }
   }
}

/** Removes a singleton entry in column col in matrix and matrix_trans and adapts
 * the nnz vectors nnzRow and nnzColChild accordingly. Sets the entries in
 * rowContainsSingletonVar to the corresponding column index in which the singleton entry occurs.
 */
void StochPresolverParallelRows::removeEntry(int col, SimpleVectorBase<int>& rowContainsSingletonVar,
      SparseStorageDynamic& matrix, SparseStorageDynamic& matrixTrans, SimpleVectorBase<int>& nnzRow, SimpleVectorBase<int>& nnzCol,
      bool parent)
{
   assert( 0 <= col && col < matrixTrans.getM() );
   assert( matrixTrans.getRowPtr(col).start + 1 == matrixTrans.getRowPtr(col).end);
   assert( nnzRow.n == matrix.getM() );
   assert( matrix.getN() == nnzCol.n );
   assert( nnzCol.n == matrixTrans.getM() );
   assert( PIPSisEQ( nnzCol[col], 1.0 ) );

   // get indices of the singleton entry int mat_trans
   const int row = matrixTrans.getJcolM(matrixTrans.getRowPtr(col).start);

   assert( row < nnzRow.n );

   // check if there are no more than one singleton entry in this row
   // if so the row is neither parallel nor nearly parallel to any other row
   if( 0 <= rowContainsSingletonVar[row] )
   {
      rowContainsSingletonVar[row] = -2.0;
      return;
   }

   // store the col in rowContainsSingletonVar
   // (possibly add offset to col so that Amat and Bmat are distinct)
   if( parent )
      rowContainsSingletonVar[row] = col;
   else 
      rowContainsSingletonVar[row] = col + nA;

   // find row in matrix
   assert(0 <= row && row < matrix.getM());
   int i = -1;
   int end = matrix.getRowPtr(row).end;
   int start = matrix.getRowPtr(row).start;

   for( i = start; i < end; i++)
   {
      if( matrix.getJcolM(i) == col )
         break;
   }
   assert( matrix.getJcolM(i) == col);
   const double val = matrix.getMat(i);

   matrix.removeEntryAtIndex( row, i);
   matrixTrans.removeEntryAtIndex( col, matrixTrans.getRowPtr(col).start );

   nnzRow[row]--;
   nnzCol[col] = 0.0;

   if( parent )
      (*singletonCoeffsColParent)[col] = val;
   else
      (*singletonCoeffsColChild)[col] = val;
}


// todo there seems to be no numerical threshold for the normalization below .. this should be fixed - rows with fairly different coefficients could be regarded equal
/// cupp can be either the rhs for the equality system or upper bounds for inequalities
void StochPresolverParallelRows::normalizeBlocksRowwise( SystemType system_type,
      SparseStorageDynamic* a_mat, SparseStorageDynamic* b_mat,
      SimpleVector* cupp, SimpleVector* clow, SimpleVector* icupp, SimpleVector* iclow) const
{
   assert(b_mat);
   assert(cupp);
   assert( b_mat->getM() == cupp->n );
   if( a_mat )
      assert( a_mat->getM() == b_mat->getM());

   if( system_type == INEQUALITY_SYSTEM )
   {
      assert( clow && iclow && icupp);
      assert( clow->n == cupp->n && iclow->n == clow->n && iclow->n == clow->n );
   }

   const int n_rows = b_mat->getM();

   /// for every row find the max value and normalize by that
   for( int row = 0; row < n_rows; row++)
   {
      double absmax = 0.0;
      bool negate_row = false;

      const int row_B_start = b_mat->getRowPtr(row).start;
      const int row_B_end = b_mat->getRowPtr(row).end;

      if( row_B_start < row_B_end )
      {
         if( b_mat->getMat(row_B_start) < 0)
            negate_row = true;

         for(int k = row_B_start; k < row_B_end; k++)
         {
            if( absmax < std::fabs(b_mat->getMat(k)) )
               absmax = std::fabs(b_mat->getMat(k));
         }
      }

      int row_A_start = 0;
      int row_A_end = 0;
      if( a_mat )
      {
    	 negate_row = false;

         row_A_start = a_mat->getRowPtr(row).start;
         row_A_end = a_mat->getRowPtr(row).end;
         if( row_A_start < row_A_end )
         {
            if( a_mat->getMat(row_A_start) < 0)
               negate_row = true;

            for(int k = row_A_start; k < row_A_end; k++)
            {
               if( absmax < std::fabs(a_mat->getMat(k)) )
                  absmax = std::fabs(a_mat->getMat(k));
            }
         }
      }

      if( PIPSisZero(absmax) )
      {
         if(system_type == EQUALITY_SYSTEM)
            (*norm_factorA)[row] = absmax;
         else
            (*norm_factorC)[row] = absmax;

         continue;
      }

      // normalize the row by dividing all entries by abs_max and possibly by -1 if negate_row.
      if(negate_row)
         absmax = absmax * -1.0;
      
      for(int k = row_B_start; k < row_B_end; k++)
         b_mat->setMat(k, b_mat->getMat(k) / absmax);

      for(int k = row_A_start; k < row_A_end; k++)
         a_mat->setMat(k, b_mat->getMat(k) / absmax);

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
         // if multiplied by a negative value, lhs and rhs have to be swapped:row_B_start
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
 * If at root, Bblock should be nullptr.
 */
// todo : I think this is wrong or at least not complete - in theory we should sort the rows first (according to the colindices)
void StochPresolverParallelRows::insertRowsIntoHashtable( boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> > &rows,
      SparseStorageDynamic* a_mat, SparseStorageDynamic* b_mat, SystemType system_type, SimpleVectorBase<int>* nnzRow )
{
   assert(b_mat);
   if(a_mat)
     assert( a_mat->getM() == b_mat->getM() );
   if( system_type == EQUALITY_SYSTEM && (b_mat != nullptr && a_mat == nullptr) )
      assert( mA == b_mat->getM() );
   if( system_type == EQUALITY_SYSTEM && a_mat != nullptr)
      assert( mA == a_mat->getM() );
      
   for(int row = 0; row < b_mat->getM(); row++)
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
      const int row_B_start =  b_mat->getRowPtr(row).start;
      const int row_B_length = b_mat->getRowPtr(row).end - row_B_start;
      
      if( a_mat )
      {
         const int row_A_start = a_mat->getRowPtr(row).start;
         const int row_A_length = a_mat->getRowPtr(row).end - row_A_start;

         if( row_B_length == 0 && row_A_length == 0 )
            continue;
         // colIndices and normalized entries are set as pointers to the original data.
         // create and insert the new element:
         rows.emplace(rowId, nA, row_B_length, &(b_mat->getJcolM()[row_B_start]), &(b_mat->getMat()[row_B_start]),
               row_A_length, &(a_mat->getJcolM()[row_A_start]), &(a_mat->getMat()[row_A_start]));
      }
      else
      {
         if(row_B_length == 0)
            continue;

         rows.emplace(rowId, nA, row_B_length, &(b_mat->getJcolM()[row_B_start]),
               &(b_mat->getMat()[row_B_start]), 0, nullptr, nullptr);
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
      for (boost::unordered_set<rowlib::rowWithEntries>::local_iterator row_one_iter = row_coefficients_hashtable.begin(i);
         row_one_iter != row_coefficients_hashtable.end(i); ++row_one_iter)
      {
         const SystemType row1_system = (row_one_iter->id < mA) ? EQUALITY_SYSTEM : INEQUALITY_SYSTEM;
         const int row1 = (row1_system == EQUALITY_SYSTEM) ? row_one_iter->id : row_one_iter->id - mA;
         const int row1_id = row_one_iter->id;
         assert(0 <= row1);

         /* if row1 has been removed in the meanwhile do not continue with it */
         if( presData.wasRowRemoved( INDEX(ROW, node, row1, false, row1_system) ) )
            continue;

         // either pairwise comparison OR lexicographical sorting and then compare only neighbors.
         // Here: pairwise comparison: // todo make lexicographical
         boost::unordered_set<rowlib::rowWithEntries>::local_iterator row_two_iter = row_one_iter;
         while ( ++row_two_iter != row_coefficients_hashtable.end(i) )
         {
            /* if at some point row1 was removed we have to get a new row1 */
            if( presData.wasRowRemoved( INDEX(ROW, node, row1, false, row1_system) ) )
               break;

            assert( row_two_iter->id != row_one_iter->id );
            const SystemType row2_system = (row_two_iter->id < mA) ? EQUALITY_SYSTEM : INEQUALITY_SYSTEM;
            const int row2 = (row2_system == EQUALITY_SYSTEM) ? row_two_iter->id : row_two_iter->id - mA;
            const int row2_id = row_two_iter->id;
            assert(0 <= row2);

            /* if row2 has been removed in the meanwhile do not continue with it */
            if( presData.wasRowRemoved( INDEX(ROW, node, row2, false, row2_system) ) )
               continue;

            bool removed = false;
            /* When two parallel rows are found, check if they are both =, both <=, or = and <= */
            if( checkRowsAreParallel( *row_one_iter, *row_two_iter) )
            {
               if( row1_system == EQUALITY_SYSTEM && row2_system == EQUALITY_SYSTEM )
               {
                  assert(row1_id == row1);
                  assert(row2_id == row2);

                  /* check if one constraint contains a singleton variable */
                  if( rowContainsSingletonVariable(row1_id) || rowContainsSingletonVariable(row2_id) )
                     removed = twoNearlyParallelEqualityRows(row1_id, row2_id, node);
                  else
                     removed = twoParallelEqualityRows(row1_id, row2_id, node);
               }
               else if( row1_system == INEQUALITY_SYSTEM && row2_system == INEQUALITY_SYSTEM )
               {
                  assert( row1 == row1_id - mA);
                  assert( row2 == row2_id - mA);

                  if( rowContainsSingletonVariable(row1_id) && rowContainsSingletonVariable(row2_id) )
                     removed = twoNearlyParallelInequalityRows(row1_id, row2_id, node);
                  else if( !rowContainsSingletonVariable(row1_id) && !rowContainsSingletonVariable(row2_id) )
                     removed = twoParallelInequalityRows(row1_id, row2_id, node);
               }
               else
               {
                  assert( (row1_system == EQUALITY_SYSTEM && row2_system == INEQUALITY_SYSTEM) ||
                     (row1_system == INEQUALITY_SYSTEM && row2_system == EQUALITY_SYSTEM) );
                  assert( (row1 == row1_id && row2 == row2_id - mA) ||
                     (row1 == row1_id - mA && row2 == row2_id) );

                  const int row_ineq_id = (row1_system == INEQUALITY_SYSTEM) ? row1_id : row2_id;
                  const int row_ineq = (row1_system == INEQUALITY_SYSTEM) ? row1 : row2;

                  const int row_eq_id = (row1_system == EQUALITY_SYSTEM) ? row1_id : row2_id;
                  const int row_eq = (row1_system == EQUALITY_SYSTEM) ? row1 : row2;


                  if( !rowContainsSingletonVariable(row_eq_id) && !rowContainsSingletonVariable(row_ineq_id) )
                     removed = parallelEqualityAndInequalityRow(row_eq, row_ineq, node);
                  else if( rowContainsSingletonVariable(row_eq_id) && !rowContainsSingletonVariable(row_ineq_id) )
                     removed = nearlyParallelEqualityAndInequalityRow(row_eq, row_ineq, node);
               }
            }

            if( removed )
               ++nRowElims;
         }
      }
   }
}

/**
 * Compare two rowWithEntries rows if the normalized coefficients are the same, at the same
 * columns indices. If yes, return true. Else, return false.
 */
bool StochPresolverParallelRows::checkRowsAreParallel( const rowlib::rowWithEntries& row1, const rowlib::rowWithEntries& row2)
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

bool StochPresolverParallelRows::twoParallelEqualityRows(int row1_id, int row2_id, int node) const
{
   const int row1 = row1_id;
   const int row2 = row2_id;
   assert( -1 <= node && node < nChildren );
   assert( row1 < mA && row2 < mA);

   if( !PIPSisEQ( (*norm_b)[row1], (*norm_b)[row2]) )
      PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found parallel equality rows with non-compatible right hand sides",
         "StochPresolverParallelRows.C", "compareRowsInCoeffHashTable");

   /* one of the rows can be discarded */
   presData.removeRedundantParallelRow( INDEX(ROW, node, row2, false, EQUALITY_SYSTEM), INDEX(ROW, node, row1, false, EQUALITY_SYSTEM) );
   return true;
}

bool StochPresolverParallelRows::twoNearlyParallelEqualityRows(int row1_id, int row2_id, int node ) const
{
   int row1 = row1_id;
   int row2 = row2_id;
   assert( -1 <= node && node < nChildren );
   assert( row1 < mA && row2 < mA );
   assert( 0 <= row1 && 0 <= row2 );

   /* make sure that row2 contains the singleton variable, otherwise switch row1 and row2 */
   if( (*rowContainsSingletonVariableA)[row2] != -1 )
      std::swap(row1, row2);

   assert((*rowContainsSingletonVariableA)[row2] != -1);

   /* get the singleton columns (if existent) */
   const int col1_index = (*rowContainsSingletonVariableA)[row1];
   const int col2_index = (*rowContainsSingletonVariableA)[row2];
   const int node_var1 = (col1_index < nA) ? -1 : node;
   const int node_var2 = (col2_index < nA) ? -1 : node;

   const double a_col1 = getSingletonCoefficient(col1_index);
   const double a_col2 = getSingletonCoefficient(col2_index);
   /* we swapped the rows so this would hold */
   assert( !PIPSisZero(a_col2) );

   const int col1 = ( node_var1 == -1 ) ? col1_index : col1_index - nA;
   const int col2 = ( node_var2 == -1 ) ? col2_index : col2_index - nA;
   assert( -1 <= col1 );
   assert( 0 <= col2 );
   if( node_var1 == -1 )
      assert( col1 < nA );
   if( node_var2 == -1 )
      assert( col2 < nA );

   /* we now can substitute the singleton column col2 out of the problem with x2 = t x1 + d where
    *    t = a_col1 / (s * a_col2)
    *    d = (b_row2 - b_row1/s)/a_col2
    */

   /* calculate t and d */
   const double s = (*norm_factorA)[row1] / (*norm_factorA)[row2];
   const double t = a_col1 / (a_col2 * s);
   const double d = ( (*norm_b)[row2] - (*norm_b)[row1] )
            * (*norm_factorA)[row2] / a_col2;

   double ixlow_col2 = (node_var2 == -1) ? (*currIxlowParent)[col2] : (*currIxlowChild)[col2];
   double ixupp_col2 = (node_var2 == -1) ? (*currIxuppParent)[col2] : (*currIxuppChild)[col2];
   double xlow_col2 = (node_var2 == -1) ? (*currxlowParent)[col2] : (*currxlowChild)[col2];
   double xupp_col2 = (node_var2 == -1) ? (*currxuppParent)[col2] : (*currxuppChild)[col2];

   /* effectively tighten bounds of variable col2 */
   double xlow_new = INF_NEG_PRES;
   double xupp_new = INF_POS_PRES;

   /* no col1 singleton */
   if( col1 == -1 )
      assert( PIPSisZero(t) );
   else
   {
      /* tighten the bounds of singleton variable x_1 (if existent) */

      /* calculate new bounds depending on the sign of t */
      if( PIPSisLT(t, 0) )
      {
         std::swap(ixlow_col2, ixupp_col2);
         std::swap(xlow_col2, xupp_col2);
      }

      if( !PIPSisZero( ixlow_col2 ) )
         xlow_new = (xlow_col2 - d) / t;
      if( !PIPSisZero( ixupp_col2 ) )
         xupp_new = (xupp_col2 - d) / t;

      const double ixlow_col1 = (node_var1 == -1) ? (*currIxlowParent)[col1] : (*currIxlowChild)[col1];
      const double ixupp_col1 = (node_var1 == -1) ? (*currIxuppParent)[col1] : (*currIxuppChild)[col1];
      const double xlow_col1 = (node_var1 == -1) ? (*currxlowParent)[col1] : (*currxlowChild)[col1];
      const double xupp_col1 = (node_var1 == -1) ? (*currxuppParent)[col1] : (*currxuppChild)[col1];

      if( !PIPSisZero( ixlow_col1 ) )
         xlow_new = std::max(xlow_new, xlow_col1);
      if( !PIPSisZero( ixupp_col1 ) )
         xupp_new = std::min(xupp_new, xupp_col1);
   }

   const INDEX row1_INDEX(ROW, node, row1, false, EQUALITY_SYSTEM);
   const INDEX row2_INDEX(ROW, node, row2, false, EQUALITY_SYSTEM);

   const INDEX col1_INDEX = (col1 == -1) ? INDEX() : INDEX(COL, node_var1, col1);
   const INDEX col2_INDEX(COL, node_var2, col2);
   /* tighten bounds of x1 */
   presData.tightenBoundsNearlyParallelRows( row1_INDEX, row2_INDEX, col1_INDEX, col2_INDEX, xlow_new, xupp_new, t, d, s);

   /* now variable substitution is possible */
   presData.substituteVariableNearlyParallelRows( row1_INDEX, row2_INDEX, col1_INDEX, col2_INDEX, t, d, s );

   /* now row2 is redundant and can be discarded */
   presData.removeRedundantParallelRow( row2_INDEX, row1_INDEX );

   return true;
}

bool StochPresolverParallelRows::twoParallelInequalityRows(int row1, int row2, int node) const
{
    /* tighten bounds in original and normalized system */
   tightenOriginalBoundsOfRow1( INEQUALITY_SYSTEM, node, row1, row2 );

    /* delete row2 in the original system */
   presData.removeRedundantParallelRow( INDEX(ROW, node, row2, false, INEQUALITY_SYSTEM), INDEX(ROW, node, row1, false, INEQUALITY_SYSTEM) );

   return true;
}

/**
 * Tightens the original lower and upper bounds of the first row, given the lower
 * and upper bounds of the second row. The normalized bounds are compared and the
 * normalizing factor of row1 is used to determine which bound can be tightened to
 * which value. The normalized bounds of row1 are also updated.
 * Assumes that both rows are inequality constraints.
 */
void StochPresolverParallelRows::tightenOriginalBoundsOfRow1(SystemType system_type, int node, int row1, int row2) const
{
   assert( system_type == INEQUALITY_SYSTEM);
   assert( currCmat );
   assert( 0 <= row1 && 0 <= row2 );
   assert( row1 < currCmat->getM() && row2 < currCmat->getM() );
   assert( norm_factorC && norm_factorC->n == currCmat->getM() );

   const double norm_factor_row1 = (*norm_factorC)[row1];
   const double norm_factor_row2 = (*norm_factorC)[row2];

   const double norm_clow_row2 = PIPSisZero( (*norm_iclow)[row2] ) ? INF_NEG_PRES : (*norm_clow)[row2];
   const double norm_cupp_row2 = PIPSisZero( (*norm_icupp)[row2] ) ? INF_POS_PRES : (*norm_cupp)[row2];

   double& norm_clow_row1 = (*norm_clow)[row1];
   double& norm_cupp_row1 = (*norm_cupp)[row1];
   double& iclow_row1 = (*norm_iclow)[row1];
   double& icupp_row1 = (*norm_icupp)[row1];

   /* test for infeasibility */
   if( ( !PIPSisZero(iclow_row1) && PIPSisLT( norm_cupp_row2, norm_clow_row1 ) )
         || ( !PIPSisZero(iclow_row1) && PIPSisLT( norm_cupp_row1, norm_clow_row2) ) )
   {
      PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found incompatible row rhs/lhs", "StochPresolverParallelRows.C", "tightenOriginalBoundsOfRow1");
   }

   double new_lhs = INF_NEG_PRES;
   double new_rhs = INF_POS_PRES;

   if( PIPSisLT( norm_clow_row1, norm_clow_row2) )
   {
      assert(norm_clow_row2 != INF_NEG_PRES);
      norm_clow_row1 = norm_clow_row2;
      iclow_row1 = 1.0;

      ( PIPSisLT( 0.0, norm_factor_row1) ) ? new_lhs = norm_factor_row1 * norm_clow_row2 : new_rhs = norm_factor_row1 * norm_clow_row2;
   }

   if( PIPSisLT( norm_cupp_row2, norm_cupp_row1) )
   {
      assert(norm_cupp_row2 != INF_POS_PRES);
      norm_cupp_row1 = norm_cupp_row2;
      icupp_row1 = 1.0;

      ( PIPSisLT( 0.0, norm_factor_row1) ) ? new_rhs = norm_factor_row1 * norm_cupp_row2 : new_lhs = norm_factor_row1 * norm_cupp_row2;
   }

   assert( PIPSisLE( new_lhs, new_rhs) );
   assert( !PIPSisZero( norm_factor_row1 / norm_factor_row2 ) );

   presData.tightenRowBoundsParallelRow( INDEX(ROW, node, row1, false, INEQUALITY_SYSTEM), INDEX(ROW, node, row2, false, INEQUALITY_SYSTEM), new_lhs, new_rhs, norm_factor_row1 / norm_factor_row2);
}

/** Returns the matrix coefficient of the singleton variable with index singleColIdx.
 * This index is possibly offset by nA, if it is in a B- or D-block.
 * If singleColIdx == -1, then 0.0 is returned.
 */
double StochPresolverParallelRows::getSingletonCoefficient(int singleton_index) const
{
   assert( singleton_index >= -1);
   if( singleton_index >= nA )
   {
      assert( singletonCoeffsColChild );
      assert( singleton_index < nA + singletonCoeffsColChild->n );
   }
   else
      assert( singletonCoeffsColParent );

   if( singleton_index == -1 )
      return 0.0;

   return ( singleton_index >= nA ) ? (*singletonCoeffsColChild)[singleton_index - nA] :
         (*singletonCoeffsColParent)[singleton_index];
}

bool StochPresolverParallelRows::rowContainsSingletonVariable( int row_index ) const
{
   assert(0 <= row_index);

   if( row_index < mA )
   {
      assert( rowContainsSingletonVariableA );
      return (*rowContainsSingletonVariableA)[row_index] != -1;
   }
   else
   {
      assert( rowContainsSingletonVariableC );
      assert( row_index - mA < rowContainsSingletonVariableC->n );
      return (*rowContainsSingletonVariableC)[row_index - mA] != -1;
   }
}

bool StochPresolverParallelRows::parallelEqualityAndInequalityRow(int row_eq, int row_ineq, int node) const
{
   /* check for infeasibility */
   if( !PIPSisZero( (*norm_iclow)[row_ineq] ) && PIPSisLT( (*norm_b)[row_eq], (*norm_clow)[row_ineq] ) )
      PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found parallel inequality and equality rows where rhs/lhs do not match",
         "StochPresolverParallelRows.C", "compareRowsInCoeffHashTable");
   if( !PIPSisZero( (*norm_icupp)[row_ineq] ) && PIPSisLT( (*norm_cupp)[row_ineq], (*norm_b)[row_eq] ) )
      PIPS_MPIabortInfeasible(MPI_COMM_WORLD, "Found parallel inequality and equality rows where rhs/lhs do not match",
         "StochPresolverParallelRows.C", "compareRowsInCoeffHashTable");

   /* remove the inequality row from the system */
   presData.removeRedundantParallelRow( INDEX(ROW, node, row_ineq, false, INEQUALITY_SYSTEM), INDEX(ROW, node, row_eq, false, EQUALITY_SYSTEM) );

   return true;
}


/**
 * Executes the Nearly Parallel Row Case 3: Both constraints are inequalities and both contain
 * a singleton variable entry.
 * The row indices rowId1, rowId2 are already de-offset, so they should be in the range [0,Cmat->getM()).
 */
bool StochPresolverParallelRows::twoNearlyParallelInequalityRows(int row1, int row2, int node) const
{
   assert( 0 <= row1 && 0 <= row2 );
   assert( rowContainsSingletonVariableC );
   assert( row1 < rowContainsSingletonVariableC->n );
   assert( row2 < rowContainsSingletonVariableC->n );
   assert( norm_factorC );
   assert( !PIPSisZero((*norm_factorC)[row2]) );

   assert( rowContainsSingletonVariable( row1 + mA ) );
   assert( rowContainsSingletonVariable( row2 + mA ) );

   /* s > 0 */
   const double s = (*norm_factorC)[row1] / (*norm_factorC)[row2];

   if( PIPSisLT(s, 0.0) )
      return false;

   const int col1_idx = (*rowContainsSingletonVariableC)[row1];
   const int col2_idx = (*rowContainsSingletonVariableC)[row2];
   const int node_col1 = (col1_idx < nA) ? -1 : node;
   const int node_col2 = (col2_idx < nA) ? -1 : node;

   /* a_col1 != 0.0, a_col2 != 0.0 */
   const double a_col1 = getSingletonCoefficient(col1_idx);
   const double a_col2 = getSingletonCoefficient(col2_idx);

   if( PIPSisZero(a_col1) || PIPSisZero(a_col2) )
      return false;

   const int col1 = (col1_idx < nA) ? col1_idx : col1_idx - nA;
   const int col2 = (col2_idx < nA) ? col2_idx : col2_idx - nA;

   /* clow_row1 = s * clow_row2 && clow_row1 = s * clow_row2 */
   if( !PIPSisEQ( (*norm_iclow)[row1], (*norm_iclow)[row2] ) || !PIPSisEQ( (*norm_icupp)[row1], (*norm_icupp)[row2] ) )
      return false;
   if( !PIPSisZero( (*norm_iclow)[row1] ) && !PIPSisEQ( (*norm_clow)[row1], s * (*norm_clow)[row2] ) )
      return false;
   if( !PIPSisZero( (*norm_icupp)[row1] ) && !PIPSisEQ( (*norm_cupp)[row1], s * (*norm_cupp)[row2] ) )
      return false;

   /* c_1 * c_2 >= 0 */
   const double c1 = (node_col1 == -1) ? (*currgParent)[col1] : (*currgChild)[col1];
   const double c2 = (node_col2 == -1) ? (*currgParent)[col2] : (*currgChild)[col2];

   if( PIPSisLT(c1 * c2, 0.0) )
      return false;

   /* lower and upper bounds matching */
   const double ixlow_col1 = (node_col1 == -1) ? (*currIxlowParent)[col1] : (*currIxlowChild)[col1];
   const double ixupp_col1 = (node_col1 == -1) ? (*currIxuppParent)[col1] : (*currIxuppChild)[col1];
   const double ixlow_col2 = (node_col2 == -1) ? (*currIxlowParent)[col2] : (*currIxlowChild)[col2];
   const double ixupp_col2 = (node_col2 == -1) ? (*currIxuppParent)[col2] : (*currIxuppChild)[col2];

   const double xlow_col1 = (node_col1 == -1) ? (*currxlowParent)[col1] : (*currxlowChild)[col1];
   const double xupp_col1 = (node_col1 == -1) ? (*currxuppParent)[col1] : (*currxuppChild)[col1];
   const double xlow_col2 = (node_col2 == -1) ? (*currxlowParent)[col2] : (*currxlowChild)[col2];
   const double xupp_col2 = (node_col2 == -1) ? (*currxuppParent)[col2] : (*currxuppChild)[col2];

   /* a_col1 * xlow_col1 = s * a_col2 * xlow_col2 */
   if( !PIPSisEQ(ixlow_col1, ixlow_col2) )
      return false;
   if( !PIPSisZero(ixlow_col1) && !PIPSisZero(ixlow_col2) &&
         !PIPSisEQ(a_col1 * xlow_col1, s * a_col2 * xlow_col2) )
      return false;

   /* a_col1 * xupp_col1 = s * a_col2 * xupp_col2 */
   if( !PIPSisEQ(ixupp_col1, ixupp_col2) )
      return false;
   if( !PIPSisZero(ixupp_col1) && !PIPSisZero(ixupp_col2) &&
         !PIPSisEQ(a_col1 * xupp_col1, s * a_col2 * xupp_col2) )
      return false;

   /* aggregate x_2: adapt objectiveCost(x_1) */
   assert(!PIPSisZero(s * a_col2));
   const double t = a_col1 / (s * a_col2);

   const INDEX row1_INDEX(ROW, node, row1, false, INEQUALITY_SYSTEM);
   const INDEX row2_INDEX(ROW, node, row2, false, INEQUALITY_SYSTEM);

   const INDEX col1_INDEX(COL, node_col1, col1);
   const INDEX col2_INDEX(COL, node_col2, col2);

   /* variable can be substituted */
   presData.substituteVariableNearlyParallelRows( row1_INDEX, row2_INDEX, col1_INDEX, col2_INDEX, t, 0.0, s);

   /* row2 is now redundant */
   presData.removeRedundantParallelRow( row2_INDEX, row1_INDEX );

   return true;
}

bool StochPresolverParallelRows::nearlyParallelEqualityAndInequalityRow(int row_eq, int row_ineq, int node) const
{
   assert( 0 <= row_eq && row_eq < mA );
   assert( rowContainsSingletonVariable(row_eq) );
   assert( !rowContainsSingletonVariable( row_ineq + mA) );

   /* compute the new variable bounds for row_eq */
   const int col_idx = (*rowContainsSingletonVariableA)[row_eq];
   const double a_col = getSingletonCoefficient(col_idx);
   assert( !PIPSisZero(a_col) );

   double xlow_new = INF_NEG_PRES;
   double xupp_new = INF_POS_PRES;

   const double s = (*norm_factorA)[row_eq] / (*norm_factorC)[row_ineq];
   const double faq =  s * a_col;

   if( PIPSisLT(0, faq) )
   {
      if( !PIPSisZero((*norm_iclow)[row_ineq]) )
         xupp_new = ( (*norm_b)[row_eq] - (*norm_clow)[row_ineq] ) * (*norm_factorA)[row_eq] / a_col;
      if( !PIPSisZero((*norm_icupp)[row_ineq]) )
         xlow_new = ( (*norm_b)[row_eq] - (*norm_cupp)[row_ineq] ) * (*norm_factorA)[row_eq] / a_col ;
   }
   else if( PIPSisLT(faq, 0.0) )
   {
      if( !PIPSisZero((*norm_iclow)[row_ineq]) )
         xlow_new = ( (*norm_b)[row_eq] - (*norm_clow)[row_ineq] ) * (*norm_factorA)[row_eq] / a_col;
      if( !PIPSisZero((*norm_icupp)[row_ineq]) )
         xupp_new = ( (*norm_b)[row_eq] - (*norm_cupp)[row_ineq] ) * (*norm_factorA)[row_eq] / a_col;
   }

   const int col = (col_idx > nA ) ? (col_idx - nA) : col_idx;
   const int node_col = (col_idx > nA ) ? node : -1;
   const bool linking_row = false;

   const INDEX row1_INDEX(ROW, node, linking_row, row_eq, EQUALITY_SYSTEM);
   const INDEX row2_INDEX(ROW, node, linking_row, row_ineq, INEQUALITY_SYSTEM);

   presData.tightenBoundsNearlyParallelRows( row1_INDEX, row2_INDEX, INDEX(COL, node_col, col), INDEX(), xlow_new, xupp_new, INF_POS_PRES, INF_POS_PRES, s ) ;

   presData.removeRedundantParallelRow( row2_INDEX, row1_INDEX );

   return true;
}
