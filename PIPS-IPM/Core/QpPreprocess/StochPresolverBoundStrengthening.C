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

// todo : make a list of fixed vars ? then variable fixing does not need to iterate the full matrix..
StochPresolverBoundStrengthening::StochPresolverBoundStrengthening(
      PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb), tightenings(0)
{
}

StochPresolverBoundStrengthening::~StochPresolverBoundStrengthening()
{
}

// todo print variables bounds strengthened
// todo print removed fixed cols
void StochPresolverBoundStrengthening::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before Bound Strengthening Presolving:" << std::endl;
   }
   countRowsCols();
#endif

   if( my_rank == 0 )
      std::cout << "Start Bound Strengthening Presolving..." << std::endl;

   int max_iter = 1; // todo
   int iter = 0;
   bool tightened;

   do
   {
      ++iter;
      tightened = false;
      /* root nodes */
      if( strenghtenBoundsInNode( EQUALITY_SYSTEM, -1) )
        tightened = true;
      if( strenghtenBoundsInNode( INEQUALITY_SYSTEM, -1) )
        tightened = true;

      // children:
      for( int node = 0; node < nChildren; node++)
      {
         // dummy child?
         if( !presData.nodeIsDummy(node, EQUALITY_SYSTEM) )
         {
            if( strenghtenBoundsInNode(EQUALITY_SYSTEM, node) )
              tightened = true;
         }

         if( !presData.nodeIsDummy(node, INEQUALITY_SYSTEM) )
         {
            if( strenghtenBoundsInNode(INEQUALITY_SYSTEM, node) )
              tightened = true;
         }
      }
   /* update bounds on all processors */
   }
   while( tightened && iter < max_iter );


   presData.allreduceLinkingVarBounds();
   presData.allreduceAndApplyLinkingRowActivities();

#ifndef NDEBUG
   MPI_Allreduce(MPI_IN_PLACE, &tightenings, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   if( my_rank == 0 )
      std::cout << "--- After " << iter << " rounds of bound strengthening and " << tightenings << " times of tightening bounds:" << std::endl;
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyActivities());
   assert(presData.verifyNnzcounters());
}

bool StochPresolverBoundStrengthening::strenghtenBoundsInNode(SystemType system_type, int node)
{
   assert( -1 <= node && node < nChildren );

   bool tightened = false;

   if( strenghtenBoundsInBlock(system_type, node, LINKING_VARS_BLOCK) )
      tightened = true;

   if( presData.hasLinking(system_type) )
   {
      //if( strenghtenBoundsInBlock(system_type, node, LINKING_CONS_BLOCK) )
      //   tightened = true;
   }

   if(node != -1)
   {
      if( strenghtenBoundsInBlock(system_type, node, CHILD_BLOCK) )
         tightened = true;
   }

   return tightened;
}

/**
 * Strengthen the variable bounds in the block matrix. If childBlock==true, then a block B_i or D_i is considered.
 * partMinActivity and partMaxActivity represent the partial row activity of the respective other block.
 */
bool StochPresolverBoundStrengthening::strenghtenBoundsInBlock( SystemType system_type, int node, BlockType block_type)
{
   assert( -1 <= node && node < nChildren );
   updatePointersForCurrentNode(node, system_type);

   // todo : if bound is too big we do not accept it
   // todo : if current entry is too small we assume that dividing by it is not numerically stable and skip it
   const double numeric_limit_entry = 1e-7;
   const double numeric_limit_bounds = 1e12; // todo
   bool tightened = false;

   const SimpleVector& xlow = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currxlowParent : *currxlowChild;
   const SimpleVector& ixlow = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currIxlowParent : *currIxlowChild;
   const SimpleVector& xupp = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currxuppParent : *currxuppChild;
   const SimpleVector& ixupp = (node == -1 || block_type == LINKING_VARS_BLOCK) ? *currIxuppParent : *currIxuppChild;

   const SimpleVector& iclow = (block_type == LINKING_CONS_BLOCK) ? *currIclowLink : *currIclow;
   const SimpleVector& clow = (block_type == LINKING_CONS_BLOCK) ? *currIneqLhsLink : *currIneqLhs;
   const SimpleVector& icupp = (block_type == LINKING_CONS_BLOCK) ? *currIcuppLink : *currIcupp;
   const SimpleVector& cupp = (block_type == LINKING_CONS_BLOCK) ? *currIneqRhsLink : *currIneqRhs;
   const SimpleVector& rhs = (block_type == LINKING_CONS_BLOCK) ? *currEqRhsLink : *currEqRhs;

   const SimpleVector& nnzs_row = (block_type == LINKING_CONS_BLOCK) ? *currNnzRowLink : *currNnzRow;

   const SparseStorageDynamic* mat;

   if(block_type == LINKING_CONS_BLOCK)
      mat = currBlmat;
   else if(block_type == LINKING_VARS_BLOCK)
      mat = currAmat;
   else
   {
      assert(node != -1);
      mat = currBmat;
   }
   assert(mat);

   /* for every row in the current block and every entry in said row check if we can improve on the currently known bounds */
   for(int row = 0; row < mat->getM(); ++row)
   {
      double actmin_part, actmax_part;
      int actmin_ubndd, actmax_ubndd;

      presData.getRowActivities(system_type, node, block_type, row, actmax_part, actmin_part, actmax_ubndd, actmin_ubndd);

      /* two or more unbounded variables make it impossible to derive new bounds so skip the row completely */
      if( actmin_ubndd >= 2 && actmax_ubndd >= 2)
         continue;

      /* if the partial row activities (so the activities of all bounded variables) exceed some limit we skip the row since no useful
       * and numerically stable bounds will be obtained here
       */
      if( (actmin_part <= -std::numeric_limits<double>::max() ) && (actmax_part >= std::numeric_limits<double>::max()) )
         continue;

      for( int j = mat->getRowPtr(row).start; j < mat->getRowPtr(row).end; j++ )
      {
         assert( mat->getRowPtr(row).end - mat->getRowPtr(row).end < nnzs_row[row] );

         // compute the possible new bounds on variable x_colIdx:
         const int col = mat->getJcolM(j);
         const double a_ik = mat->getMat(j);

         assert( !PIPSisZero(a_ik) );
         if( PIPSisLT(std::fabs(a_ik), numeric_limit_entry) )
            continue;

         /* row activities without the entry currently in focus */
         double actmin_row_without_curr = -std::numeric_limits<double>::infinity();
         double actmax_row_without_curr = std::numeric_limits<double>::infinity();

         /* subtract current entry from row activity */
         if(actmin_ubndd == 0.0)
            actmin_row_without_curr = ( PIPSisLE(a_ik, 0) ) ? actmin_part - a_ik * xupp[col] : actmin_part - a_ik * xlow[col];
         else if(actmin_ubndd == 1.0)
         {
            /* if the current entry is the unbounded one we can deduce bounds and the partial activity is the row activity excluding the current col */
            if( (PIPSisLE(a_ik, 0.0) && ixupp[col] == 0.0) || (PIPSisLE(0.0, a_ik) && ixlow[col] == 0.0) )
               actmin_row_without_curr = actmin_part;
         }

         if(actmax_ubndd == 0.0)
            actmax_row_without_curr = ( PIPSisLE(a_ik, 0) ) ? actmax_part - a_ik * xlow[col] : actmax_part - a_ik * xupp[col];
         else if(actmax_ubndd == 1.0)
         {
            /* if the current entry is the unbounded one we can deduce bounds and the partial activity is the row activity excluding the current col */
            if( (PIPSisLE(a_ik, 0.0) && ixlow[col] == 0.0) || (PIPSisLE(0.0, a_ik) && ixupp[col] == 0.0) )
               actmax_row_without_curr = actmax_part;
         }

         /* a singleton row has zero activity without the current column */
         if(nnzs_row[row] == 1)
         {
            assert(actmin_ubndd <= 1);
            assert(actmax_ubndd <= 1);
            assert( PIPSisZero(actmin_row_without_curr, feastol) );
            assert( PIPSisZero(actmax_row_without_curr, feastol) );

            actmin_row_without_curr = actmax_row_without_curr = 0;
         }

         double lbx_new = -std::numeric_limits<double>::infinity();
         double ubx_new = std::numeric_limits<double>::infinity();

         if(system_type == EQUALITY_SYSTEM)
         {
            if( PIPSisLT(0.0, a_ik) )
            {
               ubx_new = (rhs[row] - actmin_row_without_curr) / a_ik;
               lbx_new = (rhs[row] - actmax_row_without_curr) / a_ik;
            }
            else
            {
               lbx_new = (rhs[row] - actmin_row_without_curr) / a_ik;
               ubx_new = (rhs[row] - actmax_row_without_curr) / a_ik;
            }
         }
         else
         {
            if( PIPSisLT(0.0, a_ik) )
            {
               if( icupp[row] != 0.0 )
                  ubx_new = (cupp[row] - actmin_row_without_curr) / a_ik;
               if( iclow[row] != 0.0 )
                  lbx_new = (clow[row] - actmax_row_without_curr) / a_ik;
            }
            else
            {
               if( icupp[row] != 0.0 )
                  lbx_new = (cupp[row] - actmin_row_without_curr) / a_ik;
               if( iclow[row] != 0.0 )
                  ubx_new = (clow[row] - actmax_row_without_curr) / a_ik;
            }
         }

         if( std::fabs(ubx_new) > numeric_limit_bounds )
            ubx_new = std::numeric_limits<double>::infinity();
         if( std::fabs(lbx_new) > numeric_limit_bounds )
            lbx_new = -std::numeric_limits<double>::infinity();

         bool row_propagated = presData.rowPropagatedBounds(system_type, node, block_type, row, col, ubx_new, lbx_new);

         if(row_propagated && (node != -1 || my_rank == 0))
            ++tightenings;
         tightened = tightened || row_propagated;
      }
   }

   return tightened;
}
