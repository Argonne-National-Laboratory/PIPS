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
      const sROWINDEX& row_idx = presData.getSingletonRows().front();

      removed = removeSingletonRow( row_idx.system_type, row_idx.node, row_idx.index );
      
      if(removed && (row_idx.node > 0 || my_rank == 0))
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
bool StochPresolverSingletonRows::removeSingletonRow(SystemType system_type, int node, int row_idx)
{
   assert( node == -2 || !presData.nodeIsDummy(node) );

   // todo: redesign - actually done twice
   updatePointersForCurrentNode(node, system_type);
   if( (*currNnzRow)[row_idx] != 1 )
      return false;

   double ubx = std::numeric_limits<double>::infinity();
   double lbx = -std::numeric_limits<double>::infinity();
   int col_idx = -1;
   BlockType block_type = BL_MAT;

   getBoundsAndColFromSingletonRow( system_type, node, row_idx, block_type, col_idx, ubx, lbx );

   const bool linking_row = (block_type == BL_MAT);

   /* because of postsolve here we only remove correctly placed singelton rows */
   if( node != -1 && block_type == A_MAT)
      return false;

   presData.rowPropagatedBounds( system_type, node, block_type, row_idx, col_idx, ubx, lbx );

   /* singleton linking rows will not get deleted here but later by model cleanup for synchronization reasons */
   if(!linking_row)
      presData.removeRedundantRow( system_type, node, row_idx, linking_row );

   updatePointersForCurrentNode(node, system_type);
   double ubx_new = (node == -1) ? (*currxuppParent)[col_idx] : (*currxuppChild)[col_idx];
   double lbx_new = (node == -1) ? (*currxlowParent)[col_idx] : (*currxlowChild)[col_idx];
   if( PIPSisEQ(ubx_new, lbx_new) )
      presData.fixColumn( node, col_idx, ubx_new);


   if( my_rank == 0 || !linking_row )
      return true;
   else
      return false;
}

void StochPresolverSingletonRows::getBoundsAndColFromSingletonRow( SystemType system_type, int& node, int row_idx, BlockType& block_type, 
   int& col_idx, double& ubx, double& lbx )
{
   double value = 0.0; 

   if(node != -2)
   {
      updatePointersForCurrentNode(node, system_type);

      if(node == -1)
      {
         assert(currBmat);
         assert( currBmat->getRowPtr(row_idx).end - currBmat->getRowPtr(row_idx).start == 1 );
         col_idx = currBmat->getJcolM(currBmat->getRowPtr(row_idx).start);
         value = currBmat->getMat(currBmat->getRowPtr(row_idx).start);
         block_type = B_MAT;
      }
      else
      {
         assert(currAmat); assert(currBmat); assert( row_idx < currAmat->getM() ); assert( row_idx < currBmat->getM() );
         assert( (currAmat->getRowPtr(row_idx).end - currAmat->getRowPtr(row_idx).start == 1) ||
            (currBmat->getRowPtr(row_idx).end - currBmat->getRowPtr(row_idx).start == 1) );
         
         if(currAmat->getRowPtr(row_idx).end - currAmat->getRowPtr(row_idx).start == 1)
         {
            col_idx = currAmat->getJcolM(currAmat->getRowPtr(row_idx).start);
            value = currAmat->getMat(currAmat->getRowPtr(row_idx).start);
            block_type = A_MAT;
         }
         else
         {
            col_idx = currBmat->getJcolM(currBmat->getRowPtr(row_idx).start);
            value = currBmat->getMat(currBmat->getRowPtr(row_idx).start);
            block_type = B_MAT;
         }
      }
   }
   else
   {
      block_type = BL_MAT;
      // todo : implement this more efficiently - we don't want to go through all our children to check whether a singleton entry is on our process or not
      // ideally we already know and also know the child
      // same should hold for singelton columns..
      assert( node == -2 );
      for( int i = -1; i < presData.getNChildren(); ++i)
      {
         updatePointersForCurrentNode(node, system_type);

         assert( currBlmat->getRowPtr(row_idx).end - currBlmat->getRowPtr(row_idx).start == 1 ||
            currBlmat->getRowPtr(row_idx).end - currBlmat->getRowPtr(row_idx).start == 0 );
      
         if( currBlmat->getRowPtr(row_idx).end - currBlmat->getRowPtr(row_idx).start == 1 )
         {
            assert(node == -2);
            col_idx = currBlmat->getJcolM(currBlmat->getRowPtr(row_idx).start);
            value = currBlmat->getMat(currBlmat->getRowPtr(row_idx).start);
            node = i;
            break;
         }
      }
      if( node == -2 )
         return;
   }

   assert( !PIPSisEQ(value, 0) );

   if(system_type == EQUALITY_SYSTEM)
   {
      if(block_type != BL_MAT)
      {
         ubx = (*currEqRhs)[row_idx] / value;
         lbx = (*currEqRhs)[row_idx] / value;
      }
      if(block_type == BL_MAT)
      {
         ubx = (*currEqRhsLink)[row_idx] / value;
         lbx = (*currEqRhsLink)[row_idx] / value;
      }  
   }
   else
   {
      if(block_type != BL_MAT)
      {
         assert( PIPSisEQ((*currIclow)[row_idx], 1.0) || PIPSisEQ((*currIcupp)[row_idx], 1.0));

         if( PIPSisLT(value, 0.0) )
         {
            if( PIPSisEQ((*currIcupp)[row_idx], 1.0) )
               lbx = (*currIneqRhs)[row_idx] / value;
            if( PIPSisEQ((*currIclow)[row_idx], 1.0) )
               ubx = (*currIneqLhs)[row_idx] / value;
         }
         else
         {
            if( PIPSisEQ((*currIclow)[row_idx], 1.0) )
               lbx = (*currIneqLhs)[row_idx] / value;
            if( PIPSisEQ((*currIcupp)[row_idx], 1.0) )
               ubx = (*currIneqRhs)[row_idx] / value;
         }
      }
      else
         assert(block_type != BL_MAT); // todo : can is possible

   }
   assert( PIPSisLE(lbx ,ubx) );
}
