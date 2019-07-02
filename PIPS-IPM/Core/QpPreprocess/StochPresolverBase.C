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
   int n_rows_eq = 0;
   int n_rows_eq_linking = 0;
   int n_rows_ineq = 0;
   int n_rows_ineq_linking = 0;
   int n_fixed_rows = 0;
   int n_ranged_rows = 0;
   int n_singleton_rows_eq = 0;
   int n_singleton_rows_ineq = 0;
   int n_cols = 0;
   int n_boxed_cols = 0;
   int n_onesided_cols = 0;
   int n_free_cols = 0;
   int n_singleton_cols = 0;
   int n_singleton_implied_free = 0;

   /* root nodes of equality and inequality system - linking and non linking */
   if( my_rank == 0 )
   {
      updatePointersForCurrentNode(-1, EQUALITY_SYSTEM);

      int n_rows_linking_ranged = 0;
      int n_rows_eq_linking_singleton = 0;
      int n_rows_ineq_linking_singleton = 0;

      countRowsBlock(n_rows_eq, n_ranged_rows, n_fixed_rows, n_singleton_rows_eq, EQUALITY_SYSTEM, LINKING_VARS_BLOCK);
      countRowsBlock(n_rows_eq_linking, n_rows_linking_ranged, n_fixed_rows, n_rows_eq_linking_singleton, EQUALITY_SYSTEM, LINKING_CONS_BLOCK);
      n_singleton_rows_eq += n_rows_eq_linking_singleton;
      assert(n_rows_linking_ranged == 0);

      updatePointersForCurrentNode(-1, INEQUALITY_SYSTEM);

      countRowsBlock(n_rows_ineq, n_ranged_rows, n_fixed_rows, n_singleton_rows_ineq, INEQUALITY_SYSTEM, LINKING_VARS_BLOCK);
      countRowsBlock(n_rows_ineq_linking, n_rows_linking_ranged, n_fixed_rows, n_rows_ineq_linking_singleton, INEQUALITY_SYSTEM, LINKING_CONS_BLOCK);
      n_singleton_rows_ineq += n_rows_ineq_linking_singleton;
      n_ranged_rows += n_rows_linking_ranged;

      const SimpleVector& xlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.blx).vec);
      const SimpleVector& xupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.bux).vec);
      const SimpleVector& ixlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixlow).vec);
      const SimpleVector& ixupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixupp).vec);

      countBoxedColumns( n_boxed_cols, n_cols, n_free_cols, n_onesided_cols, n_singleton_cols, n_singleton_implied_free, ixlow_orig, xlow_orig, ixupp_orig, xupp_orig, LINKING_VARS_BLOCK);

      std::cout << "#linking_vars:\t" << n_cols << "\t(#singleton: " << n_singleton_cols << ", #free: " << n_free_cols << ", #impliedFree: "
            << n_singleton_implied_free << ", #onesided: " << n_onesided_cols << ", #boxed: " << n_boxed_cols << ")" << std::endl;

      std::cout << "#rows B0:\t" << n_rows_eq << "\t(#singleton: " << n_singleton_rows_eq << ")" << std::endl;
      std::cout << "#rows Bl_0:\t" << n_rows_eq_linking << "\t(#singleton: " << n_rows_eq_linking_singleton << ")" << std::endl;
      std::cout << "#rows D0:\t" << n_rows_ineq << "\t(#singleton: " << n_singleton_rows_ineq << ")" << std::endl;
      std::cout << "#rows Dl_0:\t" << n_rows_ineq_linking << "\t(#singleton: " << n_rows_ineq_linking_singleton << ", #ranged: "
            << n_rows_linking_ranged << ", #fixed: " << n_fixed_rows << ")" << std::endl;
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
         countRowsBlock(n_rows_eq, n_ranged_rows, n_fixed_rows, n_singleton_rows_eq, EQUALITY_SYSTEM, CHILD_BLOCK);
      }

      /* inequality system */
      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM) )
      {
         updatePointersForCurrentNode(node, INEQUALITY_SYSTEM);
         countRowsBlock(n_rows_ineq, n_ranged_rows, n_fixed_rows, n_singleton_rows_ineq, INEQUALITY_SYSTEM, CHILD_BLOCK);
      }

      if( !presData.nodeIsDummy( node, INEQUALITY_SYSTEM) )
      {
         const SimpleVector& xlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.blx).children[node]->vec);
         const SimpleVector& xupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.bux).children[node]->vec);
         const SimpleVector& ixlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixlow).children[node]->vec);
         const SimpleVector& ixupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixupp).children[node]->vec);

         countBoxedColumns( n_boxed_cols, n_cols, n_free_cols, n_onesided_cols, n_singleton_cols, n_singleton_implied_free, ixlow_orig, xlow_orig, ixupp_orig, xupp_orig, CHILD_BLOCK);
      }
      else if( !presData.nodeIsDummy( node, EQUALITY_SYSTEM) )
      {
         updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

         const SimpleVector& xlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.blx).children[node]->vec);
         const SimpleVector& xupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.bux).children[node]->vec);
         const SimpleVector& ixlow_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixlow).children[node]->vec);
         const SimpleVector& ixupp_orig = dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector& >(*origProb.ixupp).children[node]->vec);

         countBoxedColumns( n_boxed_cols, n_cols, n_free_cols, n_onesided_cols, n_singleton_cols, n_singleton_implied_free, ixlow_orig, xlow_orig, ixupp_orig, xupp_orig, CHILD_BLOCK);
      }
   }

#if 0//TIMING // TODO
   // count how many linking rows do not really link two blocks:
   if( hasLinking(EQUALITY_SYSTEM) )
   {
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      int* rowHasEntryInBlocks = new int[currNnzRow->n];
      for( int i = 0; i < currNnzRow->n; i++ )
         rowHasEntryInBlocks[i] = 0;
      for( size_t it = 0; it < presData.nRowElemsA->children.size(); it++)
      {
         if( !presData.nodeIsDummy( it, EQUALITY_SYSTEM))
         {
//            setCPBlmatsChild(presProb->A, (int)it); // todo
            for( int i = 0; i < currNnzRow->n; i++ )
            {
               if( currNnzRow->elements()[i] != 0.0 )
                  if( currBlmat->rowptr[i].start != currBlmat->rowptr[i].end)
                     rowHasEntryInBlocks[i]++;
            }
         }
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);

      }
      MPI_Allreduce(MPI_IN_PLACE, rowHasEntryInBlocks, currNnzRow->n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      int linkRows1Blocks = 0;
      int linkRows2Blocks = 0;
      for( int i = 0; i < currNnzRow->n; i++ )
      {
         if(rowHasEntryInBlocks[i] == 1)
            linkRows1Blocks++;
         if(rowHasEntryInBlocks[i] == 2)
            linkRows2Blocks++;
      }
      if( myRank == 0 )
      {
         cout<<"1-link rows in A: "<<linkRows1Blocks<<endl;
         cout<<"2-link rows in A: "<<linkRows2Blocks<<endl;
      }
   }
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      int* rowHasEntryInBlocks = new int[currNnzRow->n];
      for( int i = 0; i < currNnzRow->n; i++ )
         rowHasEntryInBlocks[i] = 0;
      for( size_t it = 0; it < presData.nRowElemsC->children.size(); it++)
      {
         if( !presData.nodeIsDummy( it, INEQUALITY_SYSTEM))
         {
//            setCPBlmatsChild(presProb->C, (int)it); // todo
            for( int i = 0; i < currNnzRow->n; i++ )
            {
               if( currNnzRow->elements()[i] != 0.0 )
                  if( currBlmat->rowptr[i].start != currBlmat->rowptr[i].end)
                     rowHasEntryInBlocks[i]++;
            }
         }
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      }

      MPI_Allreduce(MPI_IN_PLACE, rowHasEntryInBlocks, currNnzRow->n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      int linkRows1Blocks = 0;
      int linkRows2Blocks = 0;

      for( int i = 0; i < currNnzRow->n; i++ )
      {
         if(rowHasEntryInBlocks[i] == 1)
            linkRows1Blocks++;
         if(rowHasEntryInBlocks[i] == 2)
            linkRows2Blocks++;
      }
      if( myRank == 0 )
      {
         cout<<"1-link rows in C: "<<linkRows1Blocks<<endl;
         cout<<"2-link rows in C: "<<linkRows2Blocks<<endl;
      }
   }
#endif

   /* sync data */
   if( distributed )
   {

      int* count = new int[12];
      count[0] = n_rows_eq;
      count[1] = n_rows_ineq;
      count[2] = n_fixed_rows;
      count[3] = n_ranged_rows;
      count[4] = n_singleton_rows_eq;
      count[5] = n_singleton_rows_ineq;
      count[6] = n_cols;
      count[7] = n_boxed_cols;
      count[8] = n_onesided_cols;
      count[9] = n_free_cols;
      count[10] = n_singleton_cols;
      count[11] = n_singleton_implied_free;

      MPI_Allreduce(MPI_IN_PLACE, count, 9, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      n_rows_eq = count[0];
      n_rows_ineq = count[1];
      n_fixed_rows = count[2];
      n_ranged_rows = count[3];
      n_singleton_rows_eq = count[4];
      n_singleton_rows_ineq = count[5];
      n_cols = count[6];
      n_boxed_cols = count[7];
      n_onesided_cols = count[8];
      n_free_cols = count[9];
      n_singleton_cols = count[10];
      n_singleton_implied_free = count[11];

      delete[] count;
   }

   if( my_rank == 0 )
   {
      std::cout << "#rows_total:\t" << n_rows_eq + n_rows_ineq << "\t(#fixed: " << n_fixed_rows << ", #ranged: " << n_ranged_rows << ", #singleton: " << n_singleton_rows_eq + n_singleton_rows_ineq<< ")" << std::endl;
      std::cout << "#rows A:\t" << n_rows_eq << "\t(#singleton: " << n_singleton_rows_eq << ")" << std::endl;
      std::cout << "#rows C:\t" << n_rows_ineq << "\t(#singleton: " << n_singleton_rows_ineq << ")" << std::endl;
      std::cout << "#vars_total:\t" << n_cols << "\t(#singleton: " << n_singleton_cols << ", #free: " << n_free_cols << ", #impliedFree: "
            << n_singleton_implied_free << ", #onesided: " << n_onesided_cols << ", #boxed: " << n_boxed_cols << ")" << std::endl;
   }
}

void StochPresolverBase::countRowsBlock(int& n_rows, int& n_ranged_rows, int& n_fixed_rows, int& n_singleton_rows, SystemType system_type, BlockType block_type) const
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
   }
   else
   {
      assert(lhs); assert(rhs); assert(iclow); assert(icupp); assert( iclow->n == icupp->n );
   }

   for(int i = 0; i < rhs->n; ++i)
   {
      if(nnz_row->elements()[i] != 0.0)
      {
         n_rows++;
         if(nnz_row->elements()[i] == 1.0)
            n_singleton_rows++;

         if(system_type == EQUALITY_SYSTEM)
         {
            n_fixed_rows++;
         }
         else
         {
            if( iclow->elements()[i] != 0.0 && icupp->elements()[i] != 0.0 )
            {
               if( PIPSisEQ(lhs->elements()[i], rhs->elements()[i]))
                  n_fixed_rows++;
               else
                  n_ranged_rows++;
            }
            else
               assert(iclow->elements()[i] != 0.0 || icupp->elements()[i] != 0.0);
         }
      }
   }
}

void StochPresolverBase::countBoxedColumns(int& nBoxCols, int& nColsTotal, int& nFreeVars, int& nOnesidedVars, int& nSingletonVars, int& nSingletonVarsImpliedFree,
      const SimpleVector& ixlow_orig, const SimpleVector& xlow_orig, const SimpleVector& ixupp_orig, const SimpleVector& xupp_orig, BlockType block_type) const
{
   const SimpleVector& ixlow = (block_type == LINKING_VARS_BLOCK) ? *currIxlowParent : *currIxlowChild;
   const SimpleVector& ixupp = (block_type == LINKING_VARS_BLOCK) ? *currIxuppParent : *currIxuppChild;
   const SimpleVector& xlow = (block_type == LINKING_VARS_BLOCK) ? *currxlowParent : *currxlowChild;
   const SimpleVector& xupp = (block_type == LINKING_VARS_BLOCK) ? *currxuppParent : *currxuppChild;
   const SimpleVector& curr_nnz = (block_type == LINKING_VARS_BLOCK) ? *currNnzColParent : *currNnzColChild;
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );
   assert( ixlow.n == ixupp.n );

   for( int i = 0; i < ixlow.n; i++ )
   {
      if( curr_nnz[i] != 0.0 )
      {
         if(curr_nnz[i] == 1.0)
         {
            ++nSingletonVars;

            /* check whether or not implied free */
            if( ( ixlow_orig[i] == 0.0 || PIPSisLT(xlow_orig[i], xlow[i]) ) && ( ixupp[i] == 0.0 && PIPSisLT(xupp[i], xupp_orig[i]) ) )
            {
               ++nSingletonVarsImpliedFree;
            }
         }

         nColsTotal++;
         if( ixlow[i] != 0.0 && ixupp[i] != 0.0 )
            nBoxCols++;
         else if( ixlow[i] == 0.0 && ixupp[i] == 0.0)
            nFreeVars++;
         else
         {
            assert(ixlow[i] != 0.0 || ixupp[i] != 0.0);
            ++nOnesidedVars;
         }
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