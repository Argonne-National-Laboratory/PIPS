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

StochPresolverBoundStrengthening::StochPresolverBoundStrengthening(
      PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb),
      tightenings(0),
      local_bound_tightenings(false),
      n_linking_vars( dynamic_cast<const StochVector&>(*origProb.g).vec->length() ),
      n_eq_linking_rows( dynamic_cast<const StochVector&>(*origProb.bA).vecl->length() ),
      n_ineq_linking_rows( dynamic_cast<const StochVector&>(*origProb.bl).vecl->length() ),
      ub_linking_var(n_linking_vars),
      lb_linking_var(n_linking_vars),
      rows_ub(n_linking_vars),
      rows_lb(n_linking_vars),
      used_linking_eq_row(n_eq_linking_rows),
      used_linking_ineq_row(n_ineq_linking_rows)
{
}

StochPresolverBoundStrengthening::~StochPresolverBoundStrengthening()
{
}

void StochPresolverBoundStrengthening::resetArrays()
{
   /* upper bounds are stored inverted */
   std::fill( ub_linking_var.begin(), ub_linking_var.end(), INF_NEG_PRES );
   std::fill( lb_linking_var.begin(), lb_linking_var.end(), INF_NEG_PRES );
   std::fill( rows_ub.begin(), rows_ub.end(), INDEX() );
   std::fill( rows_lb.begin(), rows_lb.end(), INDEX() );
   std::fill( used_linking_eq_row.begin(), used_linking_eq_row.end(), false );
   std::fill( used_linking_eq_row.begin(), used_linking_eq_row.end(), false );
}

void StochPresolverBoundStrengthening::applyPresolving()
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

   resetArrays();
   presData.startBoundTightening();

   int iter = 0;
   bool tightened = false;

   do
   {
      ++iter;
      tightened = false;
      /* root nodes */
      /* it is important to do the root nodes first - we later communicate around bounds we found in A_MAT and want to get all the tightenings we can get from the root node before */
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
      communicateLinkingVarBounds();
      resetArrays();
   }
   while( tightened && iter < PRESOLVE_BOUND_STR_MAX_ITER );

   // TODO : bounds found with linking_rows: mark all rows that need to be stored by all procs and store them after tightening in the postsolver
   // later when postsolving and syncing the linking row multipliers use these rows to adjust all var bounds duals

   // TODO : can one only always undo the last bound-tightening of a variable? seems cheapest..
   // if I found some bound on a var, then found a better one which is tight -> i hope that all bounds found with the earlier are not tight - but that will not help for the comp slackness conditions..

//   presData.allreduceLinkingVarBounds();
   presData.allreduceAndApplyLinkingRowActivities();

#ifndef NDEBUG
   tightenings = PIPS_MPIgetSum(tightenings, MPI_COMM_WORLD);
   if( my_rank == 0 )
      std::cout << "--- After " << iter << " rounds of bound strengthening and " << tightenings << " times of tightening bounds:" << std::endl;
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif
   presData.endBoundTightening();

   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());
}

void StochPresolverBoundStrengthening::communicateLinkingVarBounds()
{
   PIPS_MPIgetLogicOrInPlace(local_bound_tightenings, MPI_COMM_WORLD);

   if( !local_bound_tightenings )
      return;

   // TODO : make vector bool
   std::vector<double> lbx_ubx( 2 * n_linking_vars );
   std::copy( lb_linking_var.begin(), lb_linking_var.end(), lbx_ubx.begin() );
   std::copy( ub_linking_var.begin(), ub_linking_var.end(), lbx_ubx.begin() + n_linking_vars );

   PIPS_MPIminArrayInPlace(lbx_ubx, MPI_COMM_WORLD);

   std::vector<int> best_proc_lbx_ubx( 2 * n_linking_vars, -1 );
   /* check lower bounds */
   for( unsigned int i = 0; i < n_linking_vars; ++i )
   {
      if( PIPSisEQ( lbx_ubx[i], lb_linking_var[i] ) )
      {
         best_proc_lbx_ubx[i] = my_rank;
      }
   }

   for( unsigned int i = 0; i < n_linking_vars; ++i )
   {
      if( PIPSisEQ( lbx_ubx[ n_linking_vars + i], ub_linking_var[i] ) )
      {
         best_proc_lbx_ubx[ n_linking_vars + i] = my_rank;
      }
   }

   PIPS_MPIminArrayInPlace(best_proc_lbx_ubx, MPI_COMM_WORLD);

   /* now that the best bounds and the respective processes are determined actually tighten the bounds */
   for( unsigned int i = 0; i < n_linking_vars; ++i )
   {
      if( best_proc_lbx_ubx[i] == my_rank )
      {
         const bool propagated = presData.rowPropagatedBounds(rows_lb[i], INDEX(COL, -1, i), lb_linking_var[i], INF_POS_PRES);
         if( propagated )
            ++tightenings;
      }
      else if( best_proc_lbx_ubx[i] != -1 )
      {
         const bool propagated = presData.rowPropagatedBounds(INDEX(), INDEX(COL, -1, i), lb_linking_var[i], INF_POS_PRES);
         if( propagated )
            ++tightenings;
      }
   }

   for( unsigned int i = 0; i < n_linking_vars; ++i )
   {
      if( best_proc_lbx_ubx[ n_linking_vars + i ] == my_rank )
      {
         const bool propagated = presData.rowPropagatedBounds(rows_ub[i], INDEX(COL, -1, i), INF_NEG_PRES, ub_linking_var[i]);
         if( propagated )
            ++tightenings;
      }
      else if( best_proc_lbx_ubx[ n_linking_vars + i ] != -1 )
      {
         const bool propagated = presData.rowPropagatedBounds(INDEX(), INDEX(COL, -1, i), INF_NEG_PRES, ub_linking_var[i]);
         if( propagated )
            ++tightenings;
      }
   }

   local_bound_tightenings = false;
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
   if( node == -1 || block_type == A_MAT )
      return false;

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
      if( std::fabs(actmin_part) >= PRESOLVE_BOUND_STR_MAX_PARTIAL_ACTIVITY && std::fabs(actmax_part) >= PRESOLVE_BOUND_STR_MAX_PARTIAL_ACTIVITY )
         continue;

      for( int j = mat->getRowPtr(row).start; j < mat->getRowPtr(row).end; j++ )
      {
         assert( mat->getRowPtr(row).end - mat->getRowPtr(row).end < nnzs_row[row] );

         // compute the possible new bounds on variable x_colIdx:
         const int col = mat->getJcolM(j);
         const double a_ik = mat->getMat(j);

         assert( !PIPSisZero(a_ik) );
         if( PIPSisLT(std::fabs(a_ik), PRESOLVE_BOUND_STR_NUMERIC_LIMIT_ENTRY) )
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

         if( std::fabs(ubx_new) > PRESOLVE_BOUND_STR_NUMERIC_LIMIT_BOUNDS )
            ubx_new = INF_POS_PRES;
         if( std::fabs(lbx_new) > PRESOLVE_BOUND_STR_NUMERIC_LIMIT_BOUNDS )
            lbx_new = INF_NEG_PRES;

         const int node_col = (block_type == A_MAT || node == -1) ? -1 : node;

         bool row_propagated = false;
         /* here we have to store and communicate found changes - only one process should change an upper / lower bound of a variable */
         if( block_type == A_MAT )
         {
            assert( node_col == -1 );

            /* store found upper bound if better */
            if( ubx_new != INF_POS_PRES )
            {
               /* we store upper bound inverted so that MPI communication can simple do a min over all entries */
               if( (PIPSisLT(ubx_new, cupp[row]) || PIPSisZero(icupp[row])) && PIPSisLT(ub_linking_var[col], -ubx_new) )
               {
                  local_bound_tightenings = true;
                  ub_linking_var[col] = -ubx_new;
                  rows_ub[col] = INDEX(COL, node_col, col);
               }
            }

            /* store found lower bound if better */
            if( lbx_new != INF_NEG_PRES )
            {
               if( (PIPSisLT(clow[row], lbx_new) || PIPSisZero(iclow[row])) && PIPSisLT(lbx_new, lb_linking_var[col]) )
               {
                  local_bound_tightenings = true;
                  lb_linking_var[col] = lbx_new;
                  rows_lb[col] = INDEX(COL, node_col, col);
               }
            }
         }
         else
            row_propagated = presData.rowPropagatedBounds(row_INDEX, INDEX(COL, node_col, col), lbx_new, ubx_new);
         //rowPropagatedBoundsNonTight( row_INDEX, INDEX(COL, node_col, col), lbx_new, ubx_new); // TODO: remove

         if(row_propagated && (node != -1 || my_rank == 0))
            ++tightenings;
         tightened = tightened || row_propagated;
      }
   }

   return tightened;
}
