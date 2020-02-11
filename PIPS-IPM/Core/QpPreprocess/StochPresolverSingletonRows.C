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
}

StochPresolverSingletonRows::~StochPresolverSingletonRows()
{
}
 
void StochPresolverSingletonRows::applyPresolving()
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

   if( presData.getSingletonRows().size() == 0 )
   {
      if( my_rank == 0)
         std::cout << "No more singletons left - exiting" << std::endl;
   }
#endif

   int removed_rows_local = 0;
   // main loop:
   while( !presData.getSingletonRows().empty() )
   {
      bool removed = false;
      const INDEX& row = presData.getSingletonRows().front();
      assert(row.isRow());

      removed = removeSingletonRow( row );
      
      if(removed && (row.node > 0 || my_rank == 0))
         ++removed_rows_local;
      presData.getSingletonRows().pop();
   }

   PIPS_MPIgetSumInPlace(removed_rows_local, MPI_COMM_WORLD);

   removed_rows += removed_rows_local;

   assert( presData.getSingletonRows().empty() );

   presData.syncPostsolveOfBoundsPropagatedByLinkingRows();

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
   assert( !presData.nodeIsDummy(row.node) );

   if(presData.getNnzsRow(row) != 1)
      return false;

   double xlow_new = INF_NEG_PRES;
   double xupp_new = INF_POS_PRES;

   double coeff = 0.0;

   int col_idx = -1;
   int node_col = -2;

   getBoundsAndColFromSingletonRow( row, node_col, col_idx, xlow_new, xupp_new, coeff );

   /* if the singleton entry was not found it was probably linking and someone else will remove it */
   if(node_col == -2 || col_idx == -1)
   {
      assert(row.linking);
      return false;
   }

   /* because of postsolve here we only remove correctly placed singleton rows */
   if( row.node != -1 && node_col == -1)
      return false;

   assert(!PIPSisZero(coeff));
   presData.removeSingletonRow(row, INDEX(COL, node_col, col_idx), xlow_new, xupp_new, coeff);

   if( my_rank == 0 || !row.linking )
      return true;
   else
      return false;
}

void StochPresolverSingletonRows::getBoundsAndColFromSingletonRow(const INDEX& row, int& node_col, int& col_idx, double& xlow_new, double& xupp_new, double& coeff)
{
   double coeff_singleton = 0.0;

   const int node_row = row.node;
   const int row_index = row.index;
   const SystemType system_type = row.system_type;
   const bool linking = row.linking;

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
