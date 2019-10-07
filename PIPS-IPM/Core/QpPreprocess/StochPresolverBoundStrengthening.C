/*
 * StochPresolverBoundStrengthening.C
 *
 *  Created on: 28.05.2018
 *      Author: Svenja Uslu
 */

#include "StochPresolverBoundStrengthening.h"
#include <limits>
#include <cmath>
#include "pipsdef.h"

// todo exhaustive

StochPresolverBoundStrengthening::StochPresolverBoundStrengthening(
      PresolveData& presData, StochPostsolver* postsolver) :
      StochPresolverBase(presData, postsolver)
{
   // todo
}

StochPresolverBoundStrengthening::~StochPresolverBoundStrengthening()
{
   // todo
}

// todo print variables bounds strengthened
void StochPresolverBoundStrengthening::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

#ifndef NDEBUG
   if( myRank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before Bound Strengthening Presolving:" << std::endl;
   }
   countRowsCols();
#endif

   if( myRank == 0 )
      std::cout << "Start Bound Strengthening Presolving..." << std::endl;


#ifndef NDEBUG
   StochVector* blx_clone_before = dynamic_cast<StochVector*>(presData.presProb->blx->cloneFull());
   StochVector* bux_clone_before = dynamic_cast<StochVector*>(presData.presProb->bux->cloneFull());

   blx_clone_before->componentMult(*presData.presProb->ixlow);
   bux_clone_before->componentMult(*presData.presProb->ixupp);

   double lower_inf_before = blx_clone_before->infnorm();
   double lower_two_before = blx_clone_before->twonorm();
   double upper_inf_before = bux_clone_before->infnorm();
   double upper_two_before = bux_clone_before->twonorm();

   if( myRank == 0 )
   {
      std::cout << "lower bound infnorm: " << lower_inf_before << "\t\tupper bound infnorm: " << upper_inf_before << std::endl;
      std::cout << "lower bound twonorm: " << lower_two_before << "\t\tupper bound twonorm: " << upper_two_before << std::endl;
   }
#endif

   // root:
   doBoundStrengthParent( EQUALITY_SYSTEM );
   doBoundStrengthParent( INEQUALITY_SYSTEM );

   // children:
   for( int child_it = 0; child_it < nChildren; child_it++)
   {
      // dummy child?
      if( !nodeIsDummy(child_it, EQUALITY_SYSTEM) )
         doBoundStrengthChild(child_it, EQUALITY_SYSTEM);

      if( !nodeIsDummy(child_it, INEQUALITY_SYSTEM) )
         doBoundStrengthChild(child_it, INEQUALITY_SYSTEM);
   }

   /* allreduce found deletions for linking variables */
   if( !presData.combineColAdaptParent() )
   {
      abortInfeasible(MPI_COMM_WORLD );
   }

   int a = 0;
   int b = 0;
   updateLinkingVarsBlocks(a, b);

   /* allreduce rhs lhs changes and rdeuction counters and apply them to non-zero counters */
   allreduceAndApplyNnzReductions(EQUALITY_SYSTEM);
   allreduceAndApplyNnzReductions(INEQUALITY_SYSTEM);

   allreduceAndApplyRhsLhsReductions(EQUALITY_SYSTEM);
   allreduceAndApplyRhsLhsReductions(INEQUALITY_SYSTEM);

   /* allreduce new var bounds */
   allreduceAndUpdateVarBounds();

   // Sum up individual objOffset and then add it to the global objOffset:
   sumIndivObjOffset();
   presData.addObjOffset(indivObjOffset);
   indivObjOffset = 0;

   if( myRank == 0 )
   {
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;
   }

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- After bound strengthening presolving:" << std::endl;
   countRowsCols();
   if( myRank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif


#ifndef NDEBUG
   StochVector* blx_clone_after = dynamic_cast<StochVector*>(presData.presProb->blx->cloneFull());
   StochVector* bux_clone_after = dynamic_cast<StochVector*>(presData.presProb->bux->cloneFull());

   blx_clone_after->componentMult(*presData.presProb->ixlow);
   bux_clone_after->componentMult(*presData.presProb->ixupp);

   double lower_inf_after = blx_clone_after->infnorm();
   double lower_two_after = blx_clone_after->twonorm();
   double upper_inf_after = bux_clone_after->infnorm();
   double upper_two_after = bux_clone_after->twonorm();

   if( myRank == 0 )
   {
      std::cout << "lower bound infnorm: " << lower_inf_after << "\t\tupper bound infnorm: " << upper_inf_after << std::endl;
      std::cout << "lower bound infnorm: " << lower_two_after << "\t\tupper bound infnorm: " << upper_two_after << std::endl;
   }
#endif

   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);
}

// todo no Bl0block?
void StochPresolverBoundStrengthening::doBoundStrengthParent(SystemType system_type)
{
   updatePointersForCurrentNode(-1, system_type);

   for(int rowIdx = 0; rowIdx < currAmat->m; rowIdx++)
      strenghtenBoundsInBlock( *currAmat, LINKING_VARS_BLOCK, rowIdx, 0.0, 0.0, system_type, -1);
}

// todo no Bl block
void StochPresolverBoundStrengthening::doBoundStrengthChild(int it, SystemType system_type)
{
   assert(it > -1);

   updatePointersForCurrentNode(it, system_type);

   for(int i=0; i < currAmat->m; i++) // i ~ rowIndex
   {
      // First, compute minActivity and maxActivity of the whole row (per block):
      double partMinActivityA = 0.0, partMaxActivityA = 0.0, partMinActivityB = 0.0, partMaxActivityB = 0.0;
      computeActivityBlockwise(*currAmat, i, -1, partMinActivityA, partMaxActivityA, *currxlowParent, *currIxlowParent, *currxuppParent, *currIxuppParent);
      computeActivityBlockwise(*currBmat, i, -1, partMinActivityB, partMaxActivityB, *currxlowChild, *currIxlowChild, *currxuppChild, *currIxuppChild);

      strenghtenBoundsInBlock( *currAmat, LINKING_VARS_BLOCK, i, partMinActivityB, partMaxActivityB, system_type, it);
      strenghtenBoundsInBlock( *currBmat, CHILD_BLOCK, i, partMinActivityA, partMaxActivityA, system_type, it);
   }
}

double StochPresolverBoundStrengthening::computeNewBound(const SimpleVector& bounds, double activity, double matrixEntry, int rowIdx) const
{
   assert( matrixEntry != 0.0 );
   assert( 0 <= rowIdx && rowIdx <= bounds.n);

   return (bounds[rowIdx] - activity) / matrixEntry;
}

/**
 * Strengthen the variable bounds in the block matrix. If childBlock==true, then a block B_i or D_i is considered.
 * partMinActivity and partMaxActivity represent the partial row activity of the respective other block.
 */
void StochPresolverBoundStrengthening::strenghtenBoundsInBlock( SparseStorageDynamic& matrix, BlockType block_type, int rowIdx, double partMinActivity,
      double partMaxActivity, SystemType system_type, int node)
{
   updatePointersForCurrentNode(node, system_type); // todo
   assert( 0 <= rowIdx && rowIdx < matrix.m );

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( (partMinActivity == -std::numeric_limits<double>::max() || partMinActivity == -std::numeric_limits<double>::infinity())
         && (partMaxActivity == std::numeric_limits<double>::max() || partMaxActivity == std::numeric_limits<double>::infinity()) )
      return;

   double row_activity_min = partMinActivity;
   double row_activity_max = partMaxActivity;

   SimpleVector& xlow = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currxlowParent : *currxlowChild;
   SimpleVector& ixlow = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currIxlowParent : *currIxlowChild;
   SimpleVector& xupp = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currxuppParent : *currxuppChild;
   SimpleVector& ixupp = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currIxuppParent : *currIxuppChild;

   SimpleVector& iclow = (block_type == LINKING_CONS_BLOCK) ? *currIclowLink : *currIclow;
   SimpleVector& clow = (block_type == LINKING_CONS_BLOCK) ? *currIneqLhsLink : *currIneqLhs;
   SimpleVector& icupp = (block_type == LINKING_CONS_BLOCK) ? *currIcuppLink : *currIcupp;
   SimpleVector& cupp = (block_type == LINKING_CONS_BLOCK) ? *currIneqRhsLink : *currIneqRhs;

   SimpleVector& rhs = (block_type == LINKING_CONS_BLOCK) ? *currEqRhsLink : *currEqRhs;

   // Compute remaining activity of the (complete) row
   computeActivityBlockwise(matrix, rowIdx, -1, row_activity_min, row_activity_max, xlow, ixlow, xupp, ixupp);

   if( (row_activity_min == -std::numeric_limits<double>::max() || row_activity_min == -std::numeric_limits<double>::infinity())
         && (row_activity_max == std::numeric_limits<double>::max() || row_activity_max == std::numeric_limits<double>::infinity()) )
      return;

   for( int j = matrix.rowptr[rowIdx].start; j < matrix.rowptr[rowIdx].end; j++ )
   {
      // compute the possible new bounds on variable x_colIdx:
      const int colIdx = matrix.jcolM[j];
      const double a_ik = matrix.M[j];

      assert( !PIPSisZero(a_ik) );
      assert( ixlow[colIdx] != 0.0 || ixupp[colIdx] != 0.0 );

      /* subtract current entry from row activity */
      double row_activity_min_without_curr = ( PIPSisLE(a_ik, 0) ) ? row_activity_min - a_ik * xupp[colIdx] : row_activity_min - a_ik * xlow[colIdx];
      double row_activity_max_without_curr = ( PIPSisLE(a_ik, 0) ) ? row_activity_max - a_ik * xlow[colIdx] : row_activity_max - a_ik * xupp[colIdx];

      /** Computes the new bound (bA - activity)/matrixEntry.
       * If systemType==EQUALITY_SYSTEM, then bA is the rhs currEqRhs.
       * Else, the boolean rhs defines if the lhs or the rhs of the inequality should be used as bA in the computation:
       * If rhs is true, then the right hand side cupp is used. If false, then the lhs clow is used.
       */

      double newBoundLow = -std::numeric_limits<double>::max();
      double newBoundUpp = std::numeric_limits<double>::max();

      if(system_type == EQUALITY_SYSTEM)
      {
         if( PIPSisLT(0.0, a_ik) )
         {
            if(row_activity_min_without_curr > -std::numeric_limits<double>::max())
               newBoundUpp = computeNewBound(rhs, row_activity_min_without_curr, a_ik, rowIdx);
            if(row_activity_max_without_curr < std::numeric_limits<double>::max())
               newBoundLow = computeNewBound(rhs, row_activity_max_without_curr, a_ik, rowIdx);
         }
         else
         {
            if(row_activity_min_without_curr > -std::numeric_limits<double>::max())
               newBoundLow = computeNewBound(rhs, row_activity_min_without_curr, a_ik, rowIdx);
            if(row_activity_max_without_curr < std::numeric_limits<double>::max())
               newBoundUpp = computeNewBound(rhs, row_activity_max_without_curr, a_ik, rowIdx);
         }
      }
      else
      {
         if( PIPSisLT(0.0, a_ik) )
         {
            if(row_activity_min_without_curr > -std::numeric_limits<double>::max() && icupp[rowIdx] != 0.0)
               newBoundUpp = computeNewBound(cupp, row_activity_min_without_curr, a_ik, rowIdx);
            if(row_activity_max_without_curr < std::numeric_limits<double>::max() && iclow[rowIdx] != 0.0)
               newBoundLow = computeNewBound(clow, row_activity_max_without_curr, a_ik, rowIdx);
         }
         else
         {
            if(row_activity_min_without_curr > -std::numeric_limits<double>::max() && icupp[rowIdx] != 0.0)
               newBoundLow = computeNewBound(cupp, row_activity_min_without_curr, a_ik, rowIdx);
            if(row_activity_max_without_curr < std::numeric_limits<double>::max() && iclow[rowIdx] != 0.0)
               newBoundUpp = computeNewBound(clow, row_activity_max_without_curr, a_ik, rowIdx);
         }
      }

      /* check if new bounds valid */ // todo
      if( newBoundsImplyInfeasible(newBoundLow, newBoundUpp, colIdx, ixlow.elements(), ixupp.elements(), xlow.elements(), xupp.elements()) )
         abortInfeasible(MPI_COMM_WORLD);

      /* check if new bounds imply fixation of the variable */
      double fixation_value = 0.0;
      if( newBoundsFixVariable(fixation_value, newBoundLow, newBoundUpp, colIdx, ixlow.elements(),
            ixupp.elements(), xlow.elements(), xupp.elements()) )
      {
         /* for Amat we store deletions - collect them and apply them later */
          if( node == -1 || block_type == LINKING_VARS_BLOCK )
          {
             /* only 0 process stores fixations in root node B0/Bl0 - this is to reduce MPI communications - it is not necessary */
             if( myRank == 0 && node == -1 )
                storeColValInColAdaptParent(colIdx, fixation_value);
             else if( block_type == LINKING_VARS_BLOCK )
                storeColValInColAdaptParent(colIdx, fixation_value);
          }
          else
          {
             /* delete variable */
             deleteNonlinkColumnFromSystem(node, colIdx, fixation_value);
          }
      }
      else
      {
         tightenBounds(newBoundLow, newBoundUpp, ixlow[colIdx], xlow[colIdx], ixupp[colIdx], xupp[colIdx]);
      }
   }
}
