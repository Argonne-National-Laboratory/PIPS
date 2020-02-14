/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonColumns.h"
#include "StochVectorUtilities.h"

StochPresolverSingletonColumns::StochPresolverSingletonColumns(PresolveData &presData, const sData &origProb) :
   StochPresolverBase(presData, origProb), removed_cols(0)
{
}

StochPresolverSingletonColumns::~StochPresolverSingletonColumns()
{
}

void StochPresolverSingletonColumns::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());

   presData.putLinkingVarsSyncEvent();

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before singleton columns presolving:" << std::endl;
   }
   countRowsCols();
#endif

   if( presData.getSingletonCols().size() == 0 )
   {
#ifndef NDEBUG
      if( my_rank == 0 )
         std::cout << "No more singletons left - exiting" << std::endl;
#endif
   }

   int removed_cols_run = 0;

   // main loop:
   while( !presData.getSingletonCols().empty() )
   {
      bool removed = false;
      const INDEX& col = presData.getSingletonCols().front();
      removed = removeSingletonColumn(col);

      if( removed && (col.node != -1 || my_rank == 0) )
         ++removed_cols_run;

      presData.getSingletonCols().pop();
   }

#ifndef NDEBUG
   PIPS_MPIgetSumInPlace(removed_cols_run, MPI_COMM_WORLD);
   removed_cols += removed_cols_run;
   if( my_rank == 0 )
   {
      std::cout << "--- After singleton columns presolving:" << std::endl;
      std::cout << "\tRemoved columns during singleton column elimination: " << removed_cols << std::endl;
   }
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   //todo
   presData.allreduceAndApplyBoundChanges();
   presData.allreduceAndApplyNnzChanges();
   presData.allreduceAndApplyLinkingRowActivities();
   presData.allreduceLinkingVarBounds();
   presData.allreduceAndApplyObjVecChanges();
   presData.allreduceObjOffset();

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());
}

bool StochPresolverSingletonColumns::removeSingletonColumn(const INDEX& col)
{
   assert(col.isCol());

   const int node = col.node;
   const int col_index = col.index;

   assert( -1 <= node && node < nChildren );
   assert( !presData.nodeIsDummy(node) );

   updatePointersForCurrentNode(node, EQUALITY_SYSTEM);
   const SimpleVectorBase<int> &nnzs_col = (node == -1) ? *currNnzColParent : *currNnzColChild;

   if( nnzs_col[col_index] == 0 )
      return false;

   assert(nnzs_col[col_index] == 1);

   bool found = false;
   /* find the associated row via checking transposed matrices */
   INDEX row = findRowForColumnSingleton(col, found);

   if( !found )
   {
      // TODO : singleton columns in linking parts need communication
      assert(row.index_type == EMPTY_INDEX);
      assert( node == -1 );
      return false;
   }
   assert(row.isRow());

   // TODO : why - would this not be great?
   if( row.linking )
      return false;

   assert(-1 <= row.node && row.node < nChildren);
   updatePointersForCurrentNode(row.node, row.system_type);

   /* check whether col is free / implied free */
   bool lb_implied_free = false;
   bool ub_implied_free = false;

   // todo : pivot numerical limits

   checkColImpliedFree(col, row, lb_implied_free, ub_implied_free);
   bool implied_free = lb_implied_free && ub_implied_free;

   /* if objective of variable is zero we can just remove it from the problem together with the containing row */
   double obj = (node == -1) ? (*currgParent)[col_index] : (*currgChild)[col_index];

   if( implied_free && PIPSisEQ(obj, 0.0) )
   {
      presData.removeImpliedFreeColumnSingletonEqualityRow( row, col );
      return true;
   }

   /* equalitiy singleton variables */
   if( row.system_type == EQUALITY_SYSTEM )
   {
      /* (originally) free singleton columns just get deleted together with their row */
      if( implied_free )
      {
         presData.removeImpliedFreeColumnSingletonEqualityRow( row, col );
         return true;
      }
   }

   /* inequality singleton variables */
   // TODO
   return false;
}

INDEX StochPresolverSingletonColumns::findRowForColumnSingleton( const INDEX& col, bool& found )
{
   assert( col.isCol() );
   assert(-1 <= col.node && col.node < nChildren);

   if( col.node == -1 )
      return findRowForLinkingSingleton(col.index, found);
   else
      return findRowForNonlinkingSingelton(col, found);
}

INDEX StochPresolverSingletonColumns::findRowForLinkingSingleton(int col, bool& found)
{
   /* go through all children and check linking variables for the singleton row */

   /* equality part */
   INDEX row = findRowForLinkingSingletonInSystem(col, EQUALITY_SYSTEM, found);
   if( found )
      return row;

   /* inequality part */
   return findRowForLinkingSingletonInSystem(col, INEQUALITY_SYSTEM, found);
}

INDEX StochPresolverSingletonColumns::findRowForNonlinkingSingelton(const INDEX& col, bool& found)
{
   assert(col.isCol());
   assert(col.node != -1);

   if( presData.nodeIsDummy(col.node) )
   {
      found = false;
      return INDEX();
   }

   /* check Bmat and Blmat for the singleton row */
   SystemType row_system_type = EQUALITY_SYSTEM;
   bool row_linking = false;
   const int row_node = col.node;
   int row_index = -1;

   /* equality part */
   updatePointersForCurrentNode(row_node, row_system_type);

   /* Bmat */
   assert(currBmatTrans);

   found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col.index);
   if( found )
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);

   /* Blmat */
   assert(currBlmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col.index);
   if( found )
   {
      row_linking = true;
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);
   }

   /* inequality part */
   row_system_type = INEQUALITY_SYSTEM;
   updatePointersForCurrentNode(row_node, row_system_type);

   /* Bmat */
   assert(currBmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col.index);
   if( found )
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);

   /* Blmat */
   assert(currBlmatTrans);
   found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col.index);
   if( found )
   {
      row_linking = true;
      return INDEX(ROW, row_node, row_index, row_linking, row_system_type);
   }

   return INDEX();
}

INDEX StochPresolverSingletonColumns::findRowForLinkingSingletonInSystem(int col, SystemType system_type, bool& found)
{

   for( int row_node = -1; row_node < nChildren; ++row_node )
   {
      if( !presData.nodeIsDummy(row_node) )
      {
         int row_index = -1;
         updatePointersForCurrentNode(row_node, system_type);

         /* check transposed for entry */
         if( row_node == -1 )
         {
            /* Bmat */
            assert(currBmatTrans);
            found = findRowForSingletonColumnInMatrix(*currBmatTrans, row_index, col);
            if( found )
               return INDEX(ROW, row_node, row_index, false, system_type);

            /* Blmat */
            assert(currBlmatTrans);
            found = findRowForSingletonColumnInMatrix(*currBlmatTrans, row_index, col);
            if( found )
               return INDEX(ROW, row_node, row_index, true, system_type);
         }
         else
         {
            /* Amat */
            assert(currAmatTrans);
            found = findRowForSingletonColumnInMatrix(*currAmatTrans, row_index, col);
            if( found )
               return INDEX(ROW, row_node, row_index, false, system_type);
         }
      }
   }
   assert(found == false);
   return INDEX();
}

bool StochPresolverSingletonColumns::findRowForSingletonColumnInMatrix(
   const SparseStorageDynamic &mat,
   int& row,
   const int& col)
{
   assert(col < mat.getM());
   if( mat.getRowPtr(col).start != mat.getRowPtr(col).end )
   {
      assert((mat.getRowPtr(col).end - mat.getRowPtr(col).start) == 1);
      row = mat.getJcolM(mat.getRowPtr(col).start);
      return true;
   }
   return false;
}

void StochPresolverSingletonColumns::checkColImpliedFree(const INDEX& col, const INDEX& row,
   bool &lb_implied_free,
   bool &ub_implied_free)
{
   assert(col.isCol());
   assert(row.isRow());

   ub_implied_free = lb_implied_free = false;

   updatePointersForCurrentNode(row.node, row.system_type);

//   /* check whether free */
   const double ixupp = getSimpleVecFromColStochVec(*presData.getPresProb().ixupp, col.node)[col.index];
   const double ixlow = getSimpleVecFromColStochVec(*presData.getPresProb().ixlow, col.node)[col.index];

   if( PIPSisZero(ixupp) )
      ub_implied_free = true;
   if( PIPSisZero(ixlow) )
      lb_implied_free = true;

   /* check whether bound tightening found bounds from the variables row that make it implied free */
   ub_implied_free = ub_implied_free || presData.varBoundImpliedFreeBy(true, col, row);
   lb_implied_free = lb_implied_free || presData.varBoundImpliedFreeBy(false, col, row);

}
