/*
 * StochPresolverSingletonRows.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

//#define PIPS_DEBUG
#include "StochPresolverSingletonRows.h"
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
   // presData.allreduceAndApplyNnzChanges();
   // presData.allreduceAndApplyBoundChanges();
   // presData.allreduceAndApplyLinkingRowActivities();
   // presData.allreduceLinkingVarBounds();

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());  

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before singleton Row Presolving:" << std::endl;
   }
   countRowsCols();
#endif

   if( presData.getSingletonRows().size() == 0 )
   {
#ifndef NDEBUG
      if( my_rank == 0)
         std::cout << "No more singletons left - exiting" << std::endl;
#endif
   }


   // main loop:
   while( !presData.getSingletonRows().empty() )
   {
      bool removed = false;
      removed = removeSingletonRow( presData.getSingletonRows().back().system_type, presData.getSingletonRows().back().node, 
         presData.getSingletonRows().back().index );      
      
      if(removed)
         ++removed_rows;
      presData.getSingletonRows().pop_back();
   }

   // todo : should be enough to allreduce the offset once after all of presolving
   presData.allreduceObjOffset();
   if( my_rank == 0 )
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;

#ifndef NDEBUG
   if( my_rank == 0 )
      std::cout << "--- After singleton row presolving:" << std::endl;
   countRowsCols();
   if(my_rank == 0)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   // todo
   presData.allreduceAndApplyBoundChanges();
   presData.allreduceAndApplyNnzChanges();
   presData.allreduceAndApplyLinkingRowActivities();
   presData.allreduceLinkingVarBounds();

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());
}

/** Does one round of singleton rows presolving for system A or C
 *
 * the blocks B,D,Fi,Gi (the blocks Bmat and Blmat of both A and C), the fixation and updating
 * of the columns is done. The fixed variables in one of the Amat blocks are stored in the
 * member variable colAdaptParent. Updating the blocks A,C,F0,G0 using colAdaptParent happens
 * in updateLinkingVarsBlocks() which should be called after this method.
 * Returns the number of newly found singleton rows (equality/inequality system) during adaption of B,D,Fi,Gi.
 */
bool StochPresolverSingletonRows::removeSingletonRow(SystemType system_type, int node, int row_idx)
{
   assert( !presData.nodeIsDummy(node, system_type) || node == -2 );
   double ubx = std::numeric_limits<double>::infinity();
   double lbx = -std::numeric_limits<double>::infinity();
   int col_idx = -1;
   BlockType block_type = LINKING_CONS_BLOCK;

   getBoundsAndColFromSingletonRow( system_type, node, row_idx, block_type, col_idx, ubx, lbx );

   // todo : not sure whether linking conss get deleted even when block is dummy - they should..

   presData.rowPropagatedBounds( system_type, node, block_type, row_idx, col_idx, ubx, lbx );
   presData.removeRedundantRow( system_type, node, row_idx, (block_type == LINKING_CONS_BLOCK) );

   if( my_rank == 0 || block_type != LINKING_CONS_BLOCK )
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
         assert( currAmat->rowptr[row_idx].end - currAmat->rowptr[row_idx].start == 1.0 );
         col_idx = currAmat->jcolM[currAmat->rowptr[row_idx].start];
         value = currAmat->M[currAmat->rowptr[row_idx].start];
         block_type = LINKING_VARS_BLOCK;
      }
      else
      {
         assert(currAmat);
         assert(currBmat);
         assert( row_idx < currAmat->m );
         assert( row_idx < currBmat->m );
         assert( (currAmat->rowptr[row_idx].end - currAmat->rowptr[row_idx].start == 1.0) ||
            (currBmat->rowptr[row_idx].end - currBmat->rowptr[row_idx].start == 1.0) );
         
         if(currAmat->rowptr[row_idx].end - currAmat->rowptr[row_idx].start == 1.0)
         {
            col_idx = currAmat->jcolM[currAmat->rowptr[row_idx].start];
            value = currAmat->M[currAmat->rowptr[row_idx].start];
            block_type = LINKING_VARS_BLOCK;
         }
         else
         {
            col_idx = currBmat->jcolM[currBmat->rowptr[row_idx].start];
            value = currBmat->M[currBmat->rowptr[row_idx].start];
            block_type = CHILD_BLOCK;
         }
      }
   }
   else
   {
      block_type = LINKING_CONS_BLOCK;
      // todo : implement this more efficiently - we don't want to go through all our children to check wether a singlton entry is on our process or not - ideally we
      // already know and also know the child
      assert( node == -2 );
      for( int i = -1; i < presData.getNChildren(); ++i)
      {
         updatePointersForCurrentNode(node, system_type);

         assert( currBlmat->rowptr[row_idx].end - currBlmat->rowptr[row_idx].start == 1.0 ||
            currBlmat->rowptr[row_idx].end - currBlmat->rowptr[row_idx].start == 0.0 );
      
         if( currBlmat->rowptr[row_idx].end - currBlmat->rowptr[row_idx].start == 1.0 )
         {
            assert(node == -2);
            col_idx = currBlmat->jcolM[currBlmat->rowptr[row_idx].start];
            value = currBlmat->M[currBlmat->rowptr[row_idx].start];
            node = i;
            break;
         }
      }
      if( node == -2 )
         return;
   }

   assert(value != 0.0);

   if(system_type == EQUALITY_SYSTEM)
   {
      if(block_type != LINKING_CONS_BLOCK)
         ubx = lbx = (*currEqRhs)[row_idx] / value;
      if(block_type == LINKING_CONS_BLOCK)
         ubx = lbx = (*currEqRhsLink)[row_idx] / value;
   }
   else
   {

      if(block_type != LINKING_CONS_BLOCK)
      {
         assert( (*currIclow)[row_idx] == 1.0 || (*currIcupp)[row_idx] == 1.0 );

         if( (*currIclow)[row_idx] == 1.0 )
            lbx = (*currIneqLhs)[row_idx] / value;
         if( (*currIcupp)[row_idx] == 1.0 )
            ubx = (*currIneqRhs)[row_idx] / value;

         if( PIPSisLT( value, 0.0) )
            std::swap( lbx, ubx );
      }
   }
}