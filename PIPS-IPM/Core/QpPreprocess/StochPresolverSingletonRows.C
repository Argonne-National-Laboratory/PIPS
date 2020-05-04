/*
 * StochPresolverSingletonRows.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

//#define PIPS_DEBUG
#include "StochPresolverSingletonRows.h"
#include "StochVectorUtilities.h"

#include <limits>
#include <iostream>
#include <algorithm>

StochPresolverSingletonRows::StochPresolverSingletonRows(PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb), removed_rows(0)
{
   const int n_linking_vars = dynamic_cast<const StochVector&>(*origProb.g).vec->length();

   buffer_found_singleton_equality.reserve( n_linking_vars );

   buffer_rows_lower.reserve( n_linking_vars );
   buffer_rows_upper.reserve( n_linking_vars );
   buffer_coeffs_lower.reserve( n_linking_vars );
   buffer_coeffs_upper.reserve( n_linking_vars );

   buffer_xlows.reserve( n_linking_vars );
   buffer_xupps.reserve( n_linking_vars );
}

StochPresolverSingletonRows::~StochPresolverSingletonRows()
{
}
 
bool StochPresolverSingletonRows::applyPresolving()
{
   assert(presData.presDataInSync());
   assert(presData.reductionsEmpty());

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before singleton Row Presolving:" << std::endl;
   }
   countRowsCols();
#endif

   int removed_rows_local = 0;
   resetBuffers();
   // main loop:
   while( !presData.getSingletonRows().empty() )
   {
      bool removed = false;
      const INDEX& row = presData.getSingletonRows().front();
      assert(row.isRow());

      removed = removeSingletonRow( row );
      
      if(removed && (row.getNode() > 0 || my_rank == 0))
         ++removed_rows_local;
      presData.getSingletonRows().pop();
   }

   PIPS_MPIgetSumInPlace(removed_rows_local, MPI_COMM_WORLD);

   removed_rows += removed_rows_local;

   assert( presData.getSingletonRows().empty() );

   /* sync the removal of singleton linking rows in the Ai/Ci blocks */
   removeSingletonLinkingColsSynced();

   presData.allreduceAndApplyNnzChanges();
   presData.allreduceAndApplyBoundChanges();
   presData.allreduceAndApplyLinkingRowActivities();
   presData.allreduceLinkingVarBounds();


#ifndef NDEBUG
   if(my_rank == 0)
      std::cout << "\tRemoved singleton rows during singleton row elimination: " << removed_rows << std::endl;

   if( my_rank == 0 )
      std::cout << "--- After singleton row presolving:" << std::endl;
   countRowsCols();
   if(my_rank == 0)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());

   if( removed_rows_local != 0 )
      return true;
   else
      return false;
}

/** Does one round of singleton row presolving for system A or C
 *
 * the blocks B,D,Fi,Gi (the blocks Bmat and Blmat of both A and C), the fixation and updating
 * of the columns is done. The fixed variables in one of the Amat blocks are stored in the
 * member variable colAdaptParent. Updating the blocks A,C,F0,G0 using colAdaptParent happens
 * in updateLinkingVarsBlocks() which should be called after this method.
 * Returns the number of newly found singleton rows (equality/inequality system) during adaption of B,D,Fi,Gi.
 */
bool StochPresolverSingletonRows::removeSingletonRow( const INDEX& row )
{
   assert( !presData.nodeIsDummy(row.getNode()) );

   if( presData.getNnzsRow(row) != 1 )
      return false;

   double xlow_new = INF_NEG;
   double xupp_new = INF_POS;

   double coeff = 0.0;

   int col_idx = -1;
   int node_col = -2;

   getBoundsAndColFromSingletonRow( row, node_col, col_idx, xlow_new, xupp_new, coeff );

   /* if the singleton entry was not found it was probably linking and someone else will remove it */
   if(node_col == -2 || col_idx == -1)
      assert( row.getLinking() );

   /* assert that singleton linking rows not in Bl0 are processed at the same time by exactly one process */
   if( row.getLinking() && node_col != -1 )
   {
      if(node_col == -2 || col_idx == -1)
      {
         assert( PIPS_MPIgetSum(0, MPI_COMM_WORLD) == 1 );
      }
      else
      {
         assert( PIPS_MPIgetSum(1, MPI_COMM_WORLD) == 1 );
      }
   }
   else if( row.getLinking() )
   {
      assert( PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD) );
   }

   /* because of postsolve here we only remove correctly placed singleton rows
    * singletons in Ai blocks need to be stored here and communicated later on
    */
   if( row.getNode() != -1 && node_col == -1 )
   {
      assert( !row.isLinkingRow() );

      if( PIPSisLT(xupp_new, buffer_xlows[col_idx]) || PIPSisLT( xlow_new, buffer_xupps[col_idx] ) )
         PIPS_MPIabortInfeasible("Found non-matching bounds on linking variables", "StochPresolverSingletonRows.C", "removeSingletonVar");

      /* if we already found a row - keep the better one */
      if( row.inEqSys() )
      {
         /* equality rows must be better or we are infeasible */
         if( !buffer_rows_lower[col_idx].isEmpty() )
            presData.removeRedundantRow(buffer_rows_lower[col_idx]);
         if( !buffer_rows_upper[col_idx].isEmpty() && buffer_rows_upper[col_idx].inInEqSys() )
            presData.removeRedundantRow(buffer_rows_upper[col_idx]);

         buffer_found_singleton_equality[col_idx] = 1;
         buffer_rows_lower[col_idx] = row;
         buffer_rows_upper[col_idx] = row;

         buffer_coeffs_lower[col_idx] = coeff;
         buffer_coeffs_upper[col_idx] = coeff;

         buffer_xlows[col_idx] = xlow_new;
         buffer_xupps[col_idx] = xupp_new;
      }
      else
      {
         if( buffer_found_singleton_equality[col_idx] )
            presData.removeRedundantRow(row);
         else
         {
            if( PIPSisLT(buffer_xlows[col_idx], xlow_new) )
            {
               if( !buffer_rows_lower[col_idx].isEmpty() )
                  presData.removeRedundantRow(buffer_rows_lower[col_idx]);

               buffer_rows_lower[col_idx] = row;
               buffer_coeffs_lower[col_idx] = coeff;
               buffer_xlows[col_idx] = xlow_new;
            }
            else
               presData.removeRedundantRow(row);

            if( PIPSisLT(xupp_new, buffer_xupps[col_idx]) )
            {
               if( !buffer_rows_upper[col_idx].isEmpty() )
                  presData.removeRedundantRow(buffer_rows_upper[col_idx]);

               buffer_rows_upper[col_idx] = row;
               buffer_coeffs_upper[col_idx] = coeff;
               buffer_xupps[col_idx] = xupp_new;
            }
            else
               presData.removeRedundantRow(row);
         }
      }
      return true;
   }

   assert( !PIPSisZero(coeff) );
   presData.removeSingletonRow(row, INDEX(COL, node_col, col_idx), xlow_new, xupp_new, coeff);

   if( my_rank == 0 || row.getNode() != -1 )
      return true;
   else
      return false;
}

void StochPresolverSingletonRows::getBoundsAndColFromSingletonRow(const INDEX& row, int& node_col, int& col_idx, double& xlow_new, double& xupp_new, double& coeff)
{
   double coeff_singleton = 0.0;

   const int node_row = row.getNode();
   const int row_index = row.getIndex();
   const SystemType system_type = row.getSystemType();
   const bool linking = row.getLinking();

   node_col = -2;
   col_idx = -1;

   if( !linking )
   {
      updatePointersForCurrentNode(node_row, system_type);

      if(node_row == -1)
      {
         assert(currBmat);
         assert(currBmat->getRowPtr(row_index).end - currBmat->getRowPtr(row_index).start == 1);
         coeff_singleton = currBmat->getMat(currBmat->getRowPtr(row_index).start);
         col_idx = currBmat->getJcolM(currBmat->getRowPtr(row_index).start);
         node_col = -1;
      }
      else
      {
         assert(currAmat); assert(currBmat); assert( row_index < currAmat->getM() ); assert( row_index < currBmat->getM() );
         assert( (currAmat->getRowPtr(row_index).end - currAmat->getRowPtr(row_index).start == 1) ||
            (currBmat->getRowPtr(row_index).end - currBmat->getRowPtr(row_index).start == 1) );
         
         if(currAmat->getRowPtr(row_index).end - currAmat->getRowPtr(row_index).start == 1)
         {
            coeff_singleton = currAmat->getMat(currAmat->getRowPtr(row_index).start);
            col_idx = currAmat->getJcolM(currAmat->getRowPtr(row_index).start);
            node_col = -1;
         }
         else
         {
            coeff_singleton = currBmat->getMat(currBmat->getRowPtr(row_index).start);
            col_idx = currBmat->getJcolM(currBmat->getRowPtr(row_index).start);
            node_col = node_row;
         }
      }
   }
   else
   {
      // todo : implement this more efficiently - we don't want to go through all our children to check whether a singleton entry is on our process or not
      // ideally we already know and also know the child
      // same should hold for singelton columns..

      for( int i = -1; i < presData.getNChildren(); ++i)
      {
         updatePointersForCurrentNode(i, system_type);

         assert( currBlmat->getRowPtr(row_index).end - currBlmat->getRowPtr(row_index).start == 1 ||
            currBlmat->getRowPtr(row_index).end - currBlmat->getRowPtr(row_index).start == 0 );
      
         if( currBlmat->getRowPtr(row_index).end - currBlmat->getRowPtr(row_index).start == 1 )
         {
            coeff_singleton = currBlmat->getMat(currBlmat->getRowPtr(row_index).start);
            col_idx = currBlmat->getJcolM(currBlmat->getRowPtr(row_index).start);
            node_col = i;
            break;
         }
      }

      /* singleton is not stored on our process */
      if( node_col == -2 )
         return;
   }

   assert( !PIPSisEQ(coeff_singleton, 0) );

   if(system_type == EQUALITY_SYSTEM)
   {
      const double& rhs = (linking) ? (*currEqRhsLink)[row_index] : (*currEqRhs)[row_index];

      xlow_new = rhs / coeff_singleton;
      xupp_new = rhs / coeff_singleton;
   }
   else
   {
      const double& iclow = (linking) ? (*currIclowLink)[row_index] : (*currIclow)[row_index];
      const double& icupp = (linking) ? (*currIcuppLink)[row_index] : (*currIcupp)[row_index];
      const double& clow = (linking) ? (*currIneqLhsLink)[row_index] : (*currIneqLhs)[row_index];
      const double& cupp = (linking) ? (*currIneqRhsLink)[row_index] : (*currIneqRhs)[row_index];


      assert( PIPSisEQ((*currIclow)[row_index], 1.0) || PIPSisEQ((*currIcupp)[row_index], 1.0));

      if( PIPSisLT(coeff_singleton, 0.0) )
      {
         if( PIPSisEQ( icupp, 1.0) )
            xlow_new = cupp / coeff_singleton;
         if( PIPSisEQ( iclow, 1.0) )
            xupp_new = clow / coeff_singleton;
      }
      else
      {
         if( PIPSisEQ( iclow , 1.0) )
            xlow_new = clow / coeff_singleton;
         if( PIPSisEQ( icupp, 1.0) )
            xupp_new = cupp / coeff_singleton;
      }
   }
   coeff = coeff_singleton;
   assert( PIPSisLE(xlow_new, xupp_new) );
}


void StochPresolverSingletonRows::removeSingletonLinkingColsSynced()
{
   const int n_linking_vars = buffer_xlows.size();
   /* if two procs found bounds on a variable then the better bound and after that the lower process index gets to tighten the variable */

   /* if some processes found fixations for a variable only these will compete (with their rank) for applying the tightening */
   std::vector<int> was_singleton_equality_found(n_linking_vars, 0);
   PIPS_MPImaxArray(buffer_found_singleton_equality, was_singleton_equality_found, MPI_COMM_WORLD);

   /* get best found bounds */
   std::vector<double> best_xlows(n_linking_vars, INF_NEG);
   std::vector<double> best_xupps(n_linking_vars, INF_POS);

   PIPS_MPImaxArray(buffer_xlows, best_xlows, MPI_COMM_WORLD);
   PIPS_MPIminArray(buffer_xupps, best_xupps, MPI_COMM_WORLD);

   std::vector<int> proc_lower(n_linking_vars, -1);
   std::vector<int> proc_upper(n_linking_vars, -1);

   for( int i = 0; i < n_linking_vars; ++i )
   {
      const double best_xlow = best_xlows[i];
      const double best_xupp = best_xupps[i];

      const double my_xlow = buffer_xlows[i];
      const double my_xupp = buffer_xupps[i];

      if( PIPSisLT(best_xupp, best_xlow) )
         PIPS_MPIabortInfeasible("Found non-matching bounds on linking variables", "StochPresolverSingletonRows.C", "removeSingletonLinkingColssSynced");

      const int i_found_eq = buffer_found_singleton_equality.at(i);
      const int was_eq_found = was_singleton_equality_found.at(i);

      if( was_eq_found )
         assert( best_xlow == best_xupp );
      if( i_found_eq )
      {
         assert( best_xlow == best_xupp );
         assert( my_xlow == my_xupp );
         assert( best_xlow == my_xlow );
      }

      if( i_found_eq == was_eq_found )
      {
         if( my_xlow == best_xlow )
            proc_lower[i] = my_rank;

         if( my_xupp == best_xupp )
            proc_upper[i] = my_rank;
      }
   }

   /* find the best ranking process to use its row */
   PIPS_MPIminArrayInPlace(proc_lower, MPI_COMM_WORLD);
   PIPS_MPIminArrayInPlace(proc_upper, MPI_COMM_WORLD);

   for( int i = 0; i < n_linking_vars; ++i )
   {
      const int proc_l = proc_lower[i];
      const int proc_u = proc_upper[i];

      if( proc_l == -1 && proc_u == -1 )
         continue;

      const INDEX col(COL, -1, i);
      /* equality row */
      if( proc_l == proc_u && was_singleton_equality_found[i] == 1 )
      {
         if( proc_l == my_rank )
         {
            assert( buffer_rows_lower[i] == buffer_rows_upper[i] );

            const INDEX& row = buffer_rows_lower[i];

            assert( buffer_coeffs_lower[i] == buffer_coeffs_upper[i] );
            assert( buffer_xlows[i] == best_xlows[i]);
            assert( buffer_xupps[i] == best_xupps[i]);
            assert( best_xlows[i] == best_xupps[i]);

            presData.removeSingletonRowSynced(row, col, best_xlows[i], best_xupps[i], buffer_coeffs_upper[i]);
         }
         else
         {
            presData.removeSingletonRowSynced(INDEX(), col, best_xlows[i], best_xupps[i], NAN);

            // if i found a row that is now redundant - remove it as redundant
            if( !buffer_rows_lower[i].isEmpty() )
               presData.removeRedundantRow( buffer_rows_lower[i] );
            if( !buffer_rows_upper[i].isEmpty() && buffer_rows_upper[i].inInEqSys() )
               presData.removeRedundantRow( buffer_rows_upper[i] );
         }
      }
      else
      {
         /* inequality for lower bound */
         if( proc_l != -1 )
         {
            if( proc_l == my_rank )
            {
               assert( buffer_xlows[i] == best_xlows[i] );

               presData.removeSingletonRowSynced(buffer_rows_lower[i], col, best_xlows[i], INF_POS, buffer_coeffs_lower[i]);
            }
            else
            {
               presData.removeSingletonRowSynced(INDEX(), col, best_xlows[i], INF_POS, NAN);

               // if i found a row that is now redundant - remove it as redundant
               if( !buffer_rows_lower[i].isEmpty() )
                  presData.removeRedundantRow( buffer_rows_lower[i] );
            }
         }
         /* inequality for upper bound */
         if( proc_u != -1 )
         {
            if( proc_u == my_rank )
            {
               assert( buffer_xupps[i] == best_xupps[i] );

               presData.removeSingletonRowSynced(buffer_rows_upper[i], col, INF_NEG, best_xupps[i], buffer_coeffs_upper[i]);
            }
            else
            {
               presData.removeSingletonRowSynced(INDEX(), col, INF_NEG, best_xupps[i], NAN);

               // if i found a row that is now redundant - remove it as redundant
               if( !buffer_rows_upper[i].isEmpty() )
                  presData.removeRedundantRow( buffer_rows_upper[i] );
            }
         }
      }
   }
}

void StochPresolverSingletonRows::resetBuffers()
{
   std::fill(buffer_found_singleton_equality.begin(), buffer_found_singleton_equality.end(), 0);

   std::fill(buffer_rows_lower.begin(), buffer_rows_lower.end(), INDEX());
   std::fill(buffer_rows_upper.begin(), buffer_rows_upper.end(), INDEX());

   std::fill(buffer_xlows.begin(), buffer_xlows.end(), INF_NEG);
   std::fill(buffer_xupps.begin(), buffer_xupps.end(), INF_NEG);

   std::fill(buffer_coeffs_lower.begin(), buffer_coeffs_lower.end(), NAN);
   std::fill(buffer_coeffs_upper.begin(), buffer_coeffs_upper.end(), NAN);
}
