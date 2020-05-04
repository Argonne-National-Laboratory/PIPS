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
#include "StochOptions.h"


StochPresolverBoundStrengthening::StochPresolverBoundStrengthening(
      PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb),
      limit_iter( pips_options::getIntParameter("PRESOLVE_BOUND_STR_MAX_ITER") ),
      limit_entry( pips_options::getDoubleParameter("PRESOLVE_BOUND_STR_NUMERIC_LIMIT_ENTRY") ),
      limit_partial_activity( pips_options::getDoubleParameter("PRESOLVE_BOUND_STR_MAX_PARTIAL_ACTIVITY") ),
      limit_bounds( pips_options::getDoubleParameter("PRESOLVE_BOUND_STR_NUMERIC_LIMIT_BOUNDS") ),
      tightenings(0)
{
}

StochPresolverBoundStrengthening::~StochPresolverBoundStrengthening()
{
}

bool StochPresolverBoundStrengthening::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());

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
         if( !presData.nodeIsDummy(node) )
         {
            if( strenghtenBoundsInNode(EQUALITY_SYSTEM, node) )
              tightened = true;

            if( strenghtenBoundsInNode(INEQUALITY_SYSTEM, node) )
              tightened = true;
         }
      }
   /* update bounds on all processors */
   }
   while( tightened && iter < limit_iter );

   presData.allreduceLinkingVarBounds();
   presData.allreduceAndApplyLinkingRowActivities();

#ifndef NDEBUG
   tightenings = PIPS_MPIgetSum(tightenings, MPI_COMM_WORLD);
   if( my_rank == 0 )
      std::cout << "--- After " << iter << " rounds of bound strengthening and " << tightenings << " times of tightening bounds:" << std::endl;
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   return tightened;

   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());
}

bool StochPresolverBoundStrengthening::strenghtenBoundsInNode(SystemType system_type, int node)
{
   assert( -1 <= node && node < nChildren );

   bool tightened = false;

   if( strenghtenBoundsInBlock(system_type, node, B_MAT) )
      tightened = true;

   if( presData.hasLinking(system_type) )
   {
      if( strenghtenBoundsInBlock(system_type, node, BL_MAT) )
        tightened = true;
   }

   if(node != -1)
   {
      if( strenghtenBoundsInBlock(system_type, node, A_MAT) )
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

   bool tightened = false;

   const SimpleVector& xlow = (node == -1 || block_type == A_MAT) ? *currxlowParent : *currxlowChild;
   const SimpleVector& ixlow = (node == -1 || block_type == A_MAT) ? *currIxlowParent : *currIxlowChild;
   const SimpleVector& xupp = (node == -1 || block_type == A_MAT) ? *currxuppParent : *currxuppChild;
   const SimpleVector& ixupp = (node == -1 || block_type == A_MAT) ? *currIxuppParent : *currIxuppChild;

   const SimpleVector& iclow = (block_type == BL_MAT) ? *currIclowLink : *currIclow;
   const SimpleVector& clow = (block_type == BL_MAT) ? *currIneqLhsLink : *currIneqLhs;
   const SimpleVector& icupp = (block_type == BL_MAT) ? *currIcuppLink : *currIcupp;
   const SimpleVector& cupp = (block_type == BL_MAT) ? *currIneqRhsLink : *currIneqRhs;
   const SimpleVector& rhs = (block_type == BL_MAT) ? *currEqRhsLink : *currEqRhs;

   const SimpleVectorBase<int>& nnzs_row = (block_type == BL_MAT) ? *currNnzRowLink : *currNnzRow;

   const SparseStorageDynamic* mat;

   if(block_type == BL_MAT)
      mat = currBlmat;
   else if(block_type == A_MAT)
   {
      assert(node != -1);
      mat = currAmat;
   }
   else
      mat = currBmat;

   const bool linking = (block_type == BL_MAT);
   assert(mat);

   /* for every row in the current block and every entry in said row check if we can improve on the currently known bounds */
   for(int row = 0; row < mat->getM(); ++row)
   {
      const INDEX row_INDEX(ROW, linking ? -1 : node, row, linking, system_type);

      double actmin_part, actmax_part;
      int actmin_ubndd, actmax_ubndd;

      presData.getRowActivities( row_INDEX, actmax_part, actmin_part, actmax_ubndd, actmin_ubndd);

      /* two or more unbounded variables make it impossible to derive new bounds so skip the row completely */
      if( actmin_ubndd >= 2 && actmax_ubndd >= 2)
         continue;

      /* if the partial row activities (so the activities of all bounded variables) exceed some limit we skip the row since no useful
       * and numerically stable bounds will be obtained here
       */
      if( std::fabs(actmin_part) >= limit_partial_activity && std::fabs(actmax_part) >= limit_partial_activity )
         continue;

      for( int j = mat->getRowPtr(row).start; j < mat->getRowPtr(row).end; j++ )
      {
         assert( mat->getRowPtr(row).end - mat->getRowPtr(row).end < nnzs_row[row] );

         // compute the possible new bounds on variable x_colIdx:
         const int col = mat->getJcolM(j);
         const double a_ik = mat->getMat(j);

         assert( !PIPSisZero(a_ik) );
         if( PIPSisLT(std::fabs(a_ik), limit_entry ) )
            continue;

         /* row activities without the entry currently in focus */
         double actmin_row_without_curr = -std::numeric_limits<double>::infinity();
         double actmax_row_without_curr = std::numeric_limits<double>::infinity();

         /* subtract current entry from row activity */
         if(actmin_ubndd == 0)
            actmin_row_without_curr = ( PIPSisLE(a_ik, 0) ) ? actmin_part - a_ik * xupp[col] : actmin_part - a_ik * xlow[col];
         else if(actmin_ubndd == 1)
         {
            /* if the current entry is the unbounded one we can deduce bounds and the partial activity is the row activity excluding the current col */
            if( (PIPSisLE(a_ik, 0.0) && PIPSisZero(ixupp[col])) || (PIPSisLE(0.0, a_ik) && PIPSisZero(ixlow[col])) )
               actmin_row_without_curr = actmin_part;
         }

         if(actmax_ubndd == 0)
            actmax_row_without_curr = ( PIPSisLE(a_ik, 0) ) ? actmax_part - a_ik * xlow[col] : actmax_part - a_ik * xupp[col];
         else if(actmax_ubndd == 1)
         {
            /* if the current entry is the unbounded one we can deduce bounds and the partial activity is the row activity excluding the current col */
            if( (PIPSisLE(a_ik, 0.0) && PIPSisZero(ixlow[col])) || (PIPSisLE(0.0, a_ik) && PIPSisZero(ixupp[col])) )
               actmax_row_without_curr = actmax_part;
         }

         /* a singleton row has zero activity without the current column - we skip it here though - it will be left for the singleton row presolver */
         if(nnzs_row[row] == 1)
         {
            assert(actmin_ubndd <= 1);
            assert(actmax_ubndd <= 1);
            assert( PIPSisZero(actmin_row_without_curr, feastol) );
            assert( PIPSisZero(actmax_row_without_curr, feastol) );

            actmin_row_without_curr = actmax_row_without_curr = 0;
            continue;
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
               if( !PIPSisZero(icupp[row]) )
                  ubx_new = (cupp[row] - actmin_row_without_curr) / a_ik;
               if( !PIPSisZero(iclow[row]) )
                  lbx_new = (clow[row] - actmax_row_without_curr) / a_ik;
            }
            else
            {
               if( !PIPSisZero(icupp[row]) )
                  lbx_new = (cupp[row] - actmin_row_without_curr) / a_ik;
               if( !PIPSisZero(iclow[row]) )
                  ubx_new = (clow[row] - actmax_row_without_curr) / a_ik;
            }
         }

         if( std::fabs(ubx_new) > limit_bounds )
            ubx_new = INF_POS_PRES;
         if( std::fabs(lbx_new) > limit_bounds )
            lbx_new = INF_NEG_PRES;

         const int node_col = (block_type == A_MAT || node == -1) ? -1 : node;

//         bool row_propagated = presData.rowPropagatedBounds( INDEX(ROW, node, row, linking, system_type), INDEX(COL, node_col, col), lbx_new, ubx_new);
         bool row_propagated = presData.rowPropagatedBoundsNonTight( row_INDEX, INDEX(COL, node_col, col), lbx_new, ubx_new);

         if(row_propagated && (node != -1 || my_rank == 0))
            ++tightenings;
         tightened = tightened || row_propagated;
      }
   }

   return tightened;
}
