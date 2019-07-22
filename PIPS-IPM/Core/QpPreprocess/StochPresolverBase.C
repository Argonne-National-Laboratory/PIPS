/*
 * StochPresolverBase.cpp
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

//#define PIPS_DEBUG
#include "StochPresolverBase.h"
#include "DoubleMatrixTypes.h"
#include "SmartPointer.h"
#include "pipsdef.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <cmath> // std::isfinite


StochPresolverBase::StochPresolverBase(PresolveData& presData, const sData& origProb) :
      presData(presData), origProb(origProb)
{
   getRankDistributed(MPI_COMM_WORLD, my_rank, distributed);

   localNelims = 0;
   nChildren = presData.getNChildren();
}

StochPresolverBase::~StochPresolverBase()
{
}

void StochPresolverBase::countRowsCols()// method is const but changes pointers
{
   int zero_dummy = 0; 

   int n_rows_eq_empty = 0;
   int n_rows_ineq_empty = 0;
   int n_rows_eq = 0;
   int n_rows_ineq = 0;

   int n_rows_eq_fixed = 0;
   int n_rows_fixed_ineq = 0;
   int n_rows_boxed_ineq = 0;
   int n_rows_onsided_ineq = 0;

   int n_rows_singleton_eq = 0;
   int n_rows_singleton_ineq = 0;
   
   int n_cols = 0;
   int n_cols_empty = 0;
   int n_cols_free = 0;
   int n_cols_onesided = 0;
   int n_cols_boxed = 0;
   int n_cols_singleton = 0;

   /* root nodes of equality and inequality system - linking and non linking */
   if( my_rank == 0 )
   {
      updatePointersForCurrentNode(-1, EQUALITY_SYSTEM);

      int n_rows_eq_linking = 0;
      int n_rows_ineq_linking = 0;
      int n_rows_eq_linking_empty = 0;
      int n_rows_ineq_linking_empty = 0;
      int n_rows_singleton_eq_linking = 0;
      int n_rows_singleton_ineq_linking = 0;


      countRowsBlock(n_rows_eq, n_rows_eq_empty, zero_dummy, zero_dummy, n_rows_eq_fixed, n_rows_singleton_eq, EQUALITY_SYSTEM, LINKING_VARS_BLOCK);
      assert( n_rows_eq - n_rows_eq_empty == n_rows_eq_fixed );
      countRowsBlock(n_rows_eq_linking, n_rows_eq_linking_empty, zero_dummy, zero_dummy, n_rows_eq_fixed, n_rows_singleton_eq_linking, EQUALITY_SYSTEM, 
         LINKING_CONS_BLOCK); 
      assert( zero_dummy == 0 );
      updatePointersForCurrentNode(-1, INEQUALITY_SYSTEM);

      countRowsBlock(n_rows_ineq, n_rows_ineq_empty, n_rows_onsided_ineq, n_rows_boxed_ineq, n_rows_fixed_ineq, n_rows_singleton_ineq, INEQUALITY_SYSTEM,
         LINKING_VARS_BLOCK);
      countRowsBlock(n_rows_ineq_linking, n_rows_ineq_linking_empty, n_rows_onsided_ineq, n_rows_boxed_ineq, n_rows_fixed_ineq, n_rows_singleton_ineq_linking,
         INEQUALITY_SYSTEM, LINKING_CONS_BLOCK);

      const SimpleVector& xlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.blx).vec);
      const SimpleVector& xupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.bux).vec);
      const SimpleVector& ixlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixlow).vec);
      const SimpleVector& ixupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixupp).vec);

      countBoxedColumns(n_cols, n_cols_empty, n_cols_free, n_cols_onesided, n_cols_boxed, n_cols_singleton, ixlow_orig, xlow_orig, ixupp_orig, xupp_orig,
         LINKING_VARS_BLOCK);
      assert( n_cols - n_cols_empty == n_cols_free + n_cols_onesided + n_cols_boxed);

      std::cout << "#linking_vars:\t" << n_cols << "\t(#empty: n_cols_empty " << n_cols_empty << ", #singleton: " << n_cols_singleton << 
         ", #free: " << n_cols_free << ", #onesided: " << n_cols_onesided << ", #boxed: " << n_cols_boxed << ")" << std::endl;

      std::cout << "#rows B0:\t" << n_rows_eq << "\t(#empty: n_rows_empty " << n_rows_eq_empty << ", #singleton: " << n_rows_singleton_eq << ")" << std::endl;
      std::cout << "#rows Bl_0:\t" << n_rows_eq_linking << "\t(#empty: n_rows_empty " << n_rows_eq_linking_empty << ", #singleton: " <<
         n_rows_singleton_eq_linking << ")" << std::endl;
      std::cout << "#rows D0:\t" << n_rows_ineq << "\t(#empty: n_rows_empty " << n_rows_ineq_empty << ", #singleton: " << n_rows_singleton_ineq << ")" << std::endl;
      std::cout << "#rows Dl_0:\t" << n_rows_ineq_linking << "\t(#empty: n_rows_empty " << n_rows_ineq_linking_empty << ", #singleton: " << 
         n_rows_singleton_ineq_linking << ")" << std::endl;

      n_rows_eq += n_rows_eq_linking;
      n_rows_ineq += n_rows_ineq_linking;
      n_rows_eq_empty += n_rows_eq_linking_empty;
      n_rows_ineq_empty += n_rows_ineq_linking_empty;
      n_rows_singleton_eq += n_rows_singleton_eq_linking;
      n_rows_singleton_ineq += n_rows_singleton_ineq_linking;

      assert( n_rows_eq - n_rows_eq_empty == n_rows_eq_fixed );
      assert( n_rows_ineq - n_rows_ineq_empty == n_rows_onsided_ineq + n_rows_boxed_ineq + n_rows_fixed_ineq );
   }

   /* child nodes in both systems */
   for( int node = 0; node < nChildren; node++)
   {
      assert( (presData.nodeIsDummy( node, EQUALITY_SYSTEM) && presData.nodeIsDummy( node, INEQUALITY_SYSTEM)) ||
            (!presData.nodeIsDummy( node, EQUALITY_SYSTEM) && !presData.nodeIsDummy( node, INEQUALITY_SYSTEM) ));

      /* equality system */
      if(!presData.nodeIsDummy( node, EQUALITY_SYSTEM))
      {
         updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

         if( n_rows_eq - n_rows_eq_empty != n_rows_eq_fixed  )
         {
            std::cout << my_rank << " " << node << " " << n_rows_eq << " " << n_rows_eq_empty << " " << n_rows_eq_fixed << std::endl;
         }

         assert( n_rows_eq - n_rows_eq_empty == n_rows_eq_fixed );

         countRowsBlock(n_rows_eq, n_rows_eq_empty, zero_dummy, zero_dummy, n_rows_eq_fixed, n_rows_singleton_eq, EQUALITY_SYSTEM, CHILD_BLOCK);

         assert( zero_dummy == 0 );
         assert( n_rows_eq - n_rows_eq_empty == n_rows_eq_fixed );
      }

      /* inequality system */
      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM) )
      {
         updatePointersForCurrentNode(node, INEQUALITY_SYSTEM);
         countRowsBlock(n_rows_ineq, n_rows_ineq_empty, n_rows_onsided_ineq, n_rows_boxed_ineq, n_rows_fixed_ineq, n_rows_singleton_ineq, INEQUALITY_SYSTEM, 
            CHILD_BLOCK);
         assert( n_rows_ineq - n_rows_ineq_empty == n_rows_onsided_ineq + n_rows_boxed_ineq + n_rows_fixed_ineq );
      }

      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM) )
      {
         const SimpleVector& xlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.blx).children[node]->vec);
         const SimpleVector& xupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.bux).children[node]->vec);
         const SimpleVector& ixlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixlow).children[node]->vec);
         const SimpleVector& ixupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixupp).children[node]->vec);

         countBoxedColumns(n_cols, n_cols_empty, n_cols_free, n_cols_onesided, n_cols_boxed, n_cols_singleton, ixlow_orig, xlow_orig, ixupp_orig, xupp_orig,
            CHILD_BLOCK);
         assert( n_cols - n_cols_empty == n_cols_free + n_cols_onesided + n_cols_boxed);
      }
      else if( !presData.nodeIsDummy( node, EQUALITY_SYSTEM) )
      {
         updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

         const SimpleVector& xlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.blx).children[node]->vec);
         const SimpleVector& xupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.bux).children[node]->vec);
         const SimpleVector& ixlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixlow).children[node]->vec);
         const SimpleVector& ixupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixupp).children[node]->vec);

         countBoxedColumns(n_cols, n_cols_empty, n_cols_free, n_cols_onesided, n_cols_boxed, n_cols_singleton, ixlow_orig, xlow_orig, ixupp_orig, xupp_orig,
            CHILD_BLOCK);
         assert( n_cols - n_cols_empty == n_cols_free + n_cols_onesided + n_cols_boxed);
      }
   }

#if 0//TIMING // TODO
   // count how many linking rows do not really link two blocks:
#endif

   /* sync data */
   if( distributed )
   {
      int* count = new int[15];

      count[0] = n_rows_eq;
      count[1] = n_rows_ineq;
      count[2] = n_rows_eq_empty;
      count[3] = n_rows_ineq_empty;
      count[4] = n_rows_fixed_ineq;
      count[5] = n_rows_boxed_ineq;
      count[6] = n_rows_onsided_ineq;
      count[7] = n_rows_singleton_eq;
      count[8] = n_rows_singleton_ineq;
      count[9] = n_cols;
      count[10] = n_cols_empty;
      count[11] = n_cols_free;
      count[12] = n_cols_onesided;
      count[13] = n_cols_boxed;
      count[14] = n_cols_singleton;

      MPI_Allreduce(MPI_IN_PLACE, count, 15, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      n_rows_eq = count[0];
      n_rows_ineq = count[1];
      n_rows_eq_empty = count[2];
      n_rows_ineq_empty = count[3];
      n_rows_fixed_ineq = count[4];
      n_rows_boxed_ineq = count[5];
      n_rows_onsided_ineq = count[6];
      n_rows_singleton_eq = count[7];
      n_rows_singleton_ineq = count[8];

      n_cols = count[9];
      n_cols_empty = count[10];
      n_cols_free = count[11];
      n_cols_onesided = count[12];
      n_cols_boxed = count[13];
      n_cols_singleton = count[14];

      delete[] count;
   }
   n_rows_eq_fixed = n_rows_eq - n_rows_eq_empty; 

   if( my_rank == 0 )
   {
      std::cout << "#rows_total:\t" << n_rows_eq + n_rows_ineq << "\t(#empty: " << n_rows_eq_empty + n_rows_ineq_empty << ", #fixed: " << 
         n_rows_eq_fixed + n_rows_fixed_ineq << ", #boxed: " << n_rows_boxed_ineq << ", #onesided: " << n_rows_onsided_ineq << ", #singleton: " << 
         n_rows_singleton_eq + n_rows_singleton_ineq << ")" << std::endl;
      std::cout << "#rows non-empty A:\t" << n_rows_eq - n_rows_eq_empty << "\t(#singleton: " << n_rows_singleton_eq << ")" << std::endl;
      std::cout << "#rows non-empty C:\t" << n_rows_ineq - n_rows_ineq_empty << "\t(#singleton: " << n_rows_singleton_ineq << ")" << std::endl;
      
      std::cout << "#vars_total:\t" << n_cols << "\t(#empty: " << n_cols_empty << ", #free: " << n_cols_free << ", #onesided: " << n_cols_onesided << 
         ", #boxed: " << n_cols_boxed << ", #singleton: " << n_cols_singleton << ")" << std::endl;
   }
}

void StochPresolverBase::countRowsBlock(int& n_rows, int& n_rows_empty, int& n_rows_onesided, int& n_rows_boxed, int& n_rows_fixed, int& n_rows_singleton, 
   SystemType system_type, BlockType block_type) const
{
   if(block_type == LINKING_CONS_BLOCK)
      if( !presData.hasLinking(system_type) )
         return;

   const SimpleVector* nnz_row = (block_type != LINKING_CONS_BLOCK) ? currNnzRow : currNnzRowLink;
   const SimpleVector* iclow = (block_type != LINKING_CONS_BLOCK) ? currIclow : currIclowLink;
   const SimpleVector* lhs = (block_type != LINKING_CONS_BLOCK) ? currIneqLhs : currIneqLhsLink;
   const SimpleVector* icupp = (block_type != LINKING_CONS_BLOCK) ? currIcupp : currIcuppLink;
   const SimpleVector* rhs = (block_type != LINKING_CONS_BLOCK) ? currIneqRhs : currIneqRhsLink;

   if(system_type == EQUALITY_SYSTEM)
      rhs = (block_type != LINKING_CONS_BLOCK) ? currEqRhs : currEqRhsLink;

   assert(nnz_row);
   if(system_type == EQUALITY_SYSTEM)
   {
      assert(rhs); assert(lhs == NULL); assert(iclow == NULL); assert(icupp == NULL);
      assert( nnz_row->n == rhs->n );
   }
   else
   {
      assert(lhs); assert(rhs); assert(iclow); assert(icupp);
      assert( lhs->n == rhs->n );
      assert( iclow->n == icupp->n );
      assert( lhs->n == nnz_row->n );
      assert( iclow->n == nnz_row->n );
   }

   n_rows += rhs->n;

   for(int i = 0; i < rhs->n; ++i)
   {
      if( (*nnz_row)[i] != 0.0)
      {
         if( (*nnz_row)[i] == 1.0)
            ++n_rows_singleton;

         if(system_type == EQUALITY_SYSTEM)
         {
            ++n_rows_fixed;
         }
         else
         {
            if( (*iclow)[i] != 0.0 && (*icupp)[i] != 0.0 )
            {
               if( PIPSisEQ( (*lhs)[i], (*rhs)[i]) )
                  ++n_rows_fixed;
               else
                  ++n_rows_boxed;
            }
            else
            {
               assert( (*iclow)[i] != 0.0 || (*icupp)[i] != 0.0);
               ++n_rows_onesided;
            }
         }
      }
      else
      {
         ++n_rows_empty;
      }
   }
}

void StochPresolverBase::countBoxedColumns( int& n_cols, int& n_cols_empty, int& n_cols_free, int& n_cols_onesided, int& n_cols_boxed, int& n_cols_singleton,
   const SimpleVector& ixlow_orig, const SimpleVector& xlow_orig, const SimpleVector& ixupp_orig, const SimpleVector& xupp_orig, BlockType block_type) const
{
   const SimpleVector& ixlow = (block_type == LINKING_VARS_BLOCK) ? *currIxlowParent : *currIxlowChild;
   const SimpleVector& ixupp = (block_type == LINKING_VARS_BLOCK) ? *currIxuppParent : *currIxuppChild;
   const SimpleVector& curr_nnz = (block_type == LINKING_VARS_BLOCK) ? *currNnzColParent : *currNnzColChild;

#ifndef NDEBUG
   const SimpleVector& xupp = (block_type == LINKING_VARS_BLOCK) ? *currxuppParent : *currxuppChild;
   const SimpleVector& xlow = (block_type == LINKING_VARS_BLOCK) ? *currxlowParent : *currxlowChild;
#endif
   assert( ixlow.n == ixupp.n );
   assert( ixlow.n == xlow.n );
   assert( xlow.n == xupp.n );
   assert( ixlow.n == curr_nnz.n );

   n_cols += ixlow.n;
   for( int i = 0; i < ixlow.n; i++ )
   {
      if( curr_nnz[i] != 0.0 )
      {
         if(curr_nnz[i] == 1.0)
         {
            ++n_cols_singleton;
         }

         if( ixlow[i] != 0.0 && ixupp[i] != 0.0 )
            ++n_cols_boxed;
         else if( ixlow[i] == 0.0 && ixupp[i] == 0.0)
            ++n_cols_free;
         else
            ++n_cols_onesided;
      }
      else
      {
         ++n_cols_empty;
      }
   }
}

/**
 * set all pointers to the currently necessary data
 * If node == -1 we are in the root node
 */
void StochPresolverBase::updatePointersForCurrentNode(int node, SystemType system_type)
{
   assert( !presData.nodeIsDummy(node, system_type) );
   assert(-1 <= node && node <= nChildren );
   assert(system_type == EQUALITY_SYSTEM || system_type == INEQUALITY_SYSTEM);

   const GenMatrixHandle matrix = (system_type == EQUALITY_SYSTEM) ? presData.getPresProb().A : presData.getPresProb().C;

   /* set matrix pointers for A B and Bl */
   setPointersMatrices(matrix, node);

   /* set lhs rhs for equations */
   setPointersMatrixBoundsActivities(system_type, node);

   /* set x lower upper bounds */
   setPointersVarBounds(node);

   /* set adaptions ? todo */

   /* set objective function pointers */
   setPointersObjective(node);

   /* set reduction pointers columns and rows */
   setReductionPointers(system_type, node);
}

// todo : set pointers NULL if no linking constraints?
void StochPresolverBase::setPointersMatrices(const GenMatrixHandle mat, int node)
{
   assert(-1 <= node && node < nChildren);
   const StochGenMatrix& smat = dynamic_cast<const StochGenMatrix&>(*mat);

   /* in root node only A0 and Bl0 are present */
   if( node == -1 )
   {
      currAmat = dynamic_cast<const SparseGenMatrix*>(smat.Bmat)->getStorageDynamic();
      currAmatTrans = dynamic_cast<const SparseGenMatrix*>(smat.Bmat)->getStorageDynamicTransposed();

      currBmat = NULL;
      currBmatTrans = NULL;

      currBlmat =
            dynamic_cast<const SparseGenMatrix*>(smat.Blmat)->getStorageDynamic();
      currBlmatTrans =
            dynamic_cast<const SparseGenMatrix*>(smat.Blmat)->getStorageDynamicTransposed();
   }
   else
   {
      currAmat =
            dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Amat)->getStorageDynamic();
      currAmatTrans =
            dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Amat)->getStorageDynamicTransposed();
      currBmat =
            dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Bmat)->getStorageDynamic();
      currBmatTrans =
            dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Bmat)->getStorageDynamicTransposed();
      currBlmat =
            dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Blmat)->getStorageDynamic();
      currBlmatTrans =
            dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Blmat)->getStorageDynamicTransposed();
   }
}

// todo make one for bounds one for activities
void StochPresolverBase::setPointersMatrixBoundsActivities(SystemType system_type, int node)
{
   assert(-1 <= node && node < nChildren);

   /* non-linking constraints */
   const StochVector& lhs = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<const StochVector&>(*(presData.getPresProb().bA))
         : dynamic_cast<const StochVector&>(*(presData.getPresProb().bl));
   const StochVector& lhs_idx = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<const StochVector&>(*(presData.getPresProb().bA))
         : dynamic_cast<const StochVector&>(*(presData.getPresProb().iclow));
   const StochVector& rhs = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<const StochVector&>(*(presData.getPresProb().bA))
         : dynamic_cast<const StochVector&>(*(presData.getPresProb().bu));
   const StochVector& rhs_idx = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<const StochVector&>(*(presData.getPresProb().bA))
         : dynamic_cast<const StochVector&>(*(presData.getPresProb().icupp));

   if( system_type == EQUALITY_SYSTEM )
   {
      if( node == -1 )
      {
         currEqRhs = dynamic_cast<const SimpleVector*>(rhs.vec);
      }
      else
      {
         currEqRhs = dynamic_cast<const SimpleVector*>(rhs.children[node]->vec);
         assert(rhs.children[node]->vecl == NULL);
      }
      currIneqLhs = currIclow = currIneqRhs = currIcupp = currIneqLhsLink =
            currIclowLink = currIneqRhsLink = currIcuppLink = NULL;

      if( presData.hasLinking(system_type) )
         currEqRhsLink = dynamic_cast<const SimpleVector*>(rhs.vecl);
      else
         currEqRhsLink = NULL;
   }
   else
   {

      if( node == -1 )
      {
         currIneqLhs = dynamic_cast<const SimpleVector*>(lhs.vec);
         currIclow = dynamic_cast<const SimpleVector*>(lhs_idx.vec);
         currIneqRhs = dynamic_cast<const SimpleVector*>(rhs.vec);
         currIcupp = dynamic_cast<const SimpleVector*>(rhs_idx.vec);
      }
      else
      {
         currIneqLhs = dynamic_cast<const SimpleVector*>(lhs.children[node]->vec);
         currIclow = dynamic_cast<const SimpleVector*>(lhs_idx.children[node]->vec);
         currIneqRhs = dynamic_cast<const SimpleVector*>(rhs.children[node]->vec);
         currIcupp = dynamic_cast<const SimpleVector*>(rhs_idx.children[node]->vec);

         assert(lhs.children[node]->vecl == NULL);
         assert(lhs_idx.children[node]->vecl == NULL);
         assert(rhs.children[node]->vecl == NULL);
         assert(rhs_idx.children[node]->vecl == NULL);
      }

      currEqRhs = currEqRhsLink = NULL;

      if( presData.hasLinking(system_type) )
      {
         currIneqLhsLink = dynamic_cast<const SimpleVector*>(lhs.vecl);
         currIclowLink = dynamic_cast<const SimpleVector*>(lhs_idx.vecl);
         currIneqRhsLink = dynamic_cast<const SimpleVector*>(rhs.vecl);
         currIcuppLink = dynamic_cast<const SimpleVector*>(rhs_idx.vecl);
      }
      else
      {
         currIneqLhsLink = currIclowLink = currIneqRhsLink = currIcuppLink = NULL;
      }
   }
}

void StochPresolverBase::setPointersVarBounds(int node)
{
   assert(-1 <= node && node <= nChildren);

   currxlowParent = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().blx)).vec);
   currIxlowParent = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixlow)).vec);
   currxuppParent = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bux)).vec);
   currIxuppParent = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixupp)).vec);

   assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().blx)).vecl == NULL);
   assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixlow)).vecl == NULL);
   assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().bux)).vecl == NULL);
   assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixupp)).vecl == NULL);

   if(node != -1)
   {
      currxlowChild = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().blx)).children[node]->vec);
      currxuppChild = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().bux)).children[node]->vec);
      currIxlowChild = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixlow)).children[node]->vec);
      currIxuppChild = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixupp)).children[node]->vec);

      assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().blx)).children[node]->vecl == NULL);
      assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().bux)).children[node]->vecl == NULL);
      assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixlow)).children[node]->vecl == NULL);
      assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().ixupp)).children[node]->vecl == NULL);
   }
   else
      currxlowChild = currxuppChild = currIxlowChild = currIxuppChild = NULL;
}

void StochPresolverBase::setPointersObjective(int node)
{
   currgParent = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().g)).vec);
   if(node != -1)
   {
      currgChild = dynamic_cast<const SimpleVector*>(dynamic_cast<const StochVector&>(*(presData.getPresProb().g)).children[node]->vec);
      assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().g)).children[node]->vecl == NULL);
   }
   else
      currgChild = NULL;

   assert(dynamic_cast<const StochVector&>(*(presData.getPresProb().g)).vecl == NULL);

}

void StochPresolverBase::setReductionPointers(SystemType system_type, int node){
   assert(-1 <= node && node <= nChildren);

//   const SimpleVectorHandle row_red = (system_type == EQUALITY_SYSTEM) ? presData.nnzs_row_A_chgs : presData.nnzs_row_C_chgs;
   const StochVector& row_nnz = (system_type == EQUALITY_SYSTEM) ? presData.getNnzsRowA() : presData.getNnzsRowC();

   /* rows */
   if( node == -1)
   {
      currNnzRow = dynamic_cast<const SimpleVector*>(row_nnz.vec);
   }
   else
   {
      assert(row_nnz.children[node]->vec != NULL);

      currNnzRow = dynamic_cast<const SimpleVector*>(row_nnz.children[node]->vec);
   }

   if( presData.hasLinking(system_type) )
   {
//      currRedRowLink = &(*row_red);
      currNnzRowLink = dynamic_cast<const SimpleVector*>(row_nnz.vecl);;
   }
   else
//      currRedRowLink = currNnzRowLink = NULL;
      ;
   /* colums */
//   currRedColParent = &(*presData.nnzs_col_chgs);
   currNnzColParent = dynamic_cast<const SimpleVector*>(presData.getNnzsCol().vec);;

   assert(presData.getNnzsCol().vecl == NULL);

   if(node != -1)
   {
      currNnzColChild = dynamic_cast<const SimpleVector*>(presData.getNnzsCol().children[node]->vec);

      assert(presData.getNnzsCol().children[node]->vecl == NULL);
   }
   else
      currNnzColChild = NULL;
}

//void StochPresolverBase::setCurrentPointersToNull()
//{
//   currAmat = NULL;
//   currAmatTrans = NULL;
//   currBmat = NULL;
//   currBmatTrans = NULL;
//   currBlmat = NULL;
//   currBlmatTrans = NULL;
//
//   currxlowParent = NULL;
//   currIxlowParent = NULL;
//   currxuppParent = NULL;
//   currIxuppParent = NULL;
//   currxlowChild = NULL;
//   currIxlowChild = NULL;
//   currxuppChild = NULL;
//   currIxuppChild = NULL;
//   currEqRhs = NULL;
//   currIneqLhs = NULL;
//   currIclow = NULL;
//   currIneqRhs = NULL;
//   currIcupp = NULL;
//   currEqRhsLink = NULL;
//   currIneqLhsLink = NULL;
//   currIclowLink = NULL;
//   currIneqRhsLink = NULL;
//   currIcuppLink = NULL;
//
//   currgParent = NULL;
//   currgChild = NULL;
//
//   currRedRow = NULL;
//   currNnzRow = NULL;
//   currRedRowLink = NULL;
//   currNnzColParent = NULL;
//   currNnzColChild = NULL;
//}
//
//bool StochPresolverBase::variableFixationValid(double fixation_value, const double& ixlow, const double& xlow, const double& ixupp, const double& xupp, bool print_message) const
//{
//   if( (ixlow != 0.0 && PIPSisLT(fixation_value, xlow)) || (ixupp != 0.0 && PIPSisLT(xupp, fixation_value)) )
//   {
//      std::cout << "Presolve detected infeasibility! Fixation of variable to invalid value - value: " << fixation_value << "\t bounds: x € ["
//            << ((ixlow == 0.0) ? -std::numeric_limits<double>::infinity() : xlow) << ", " << ((ixupp == 0.0) ? std::numeric_limits<double>::infinity() : xupp) << "]" << std::endl;
//
//      return false;
//   }
//   return true;
//}