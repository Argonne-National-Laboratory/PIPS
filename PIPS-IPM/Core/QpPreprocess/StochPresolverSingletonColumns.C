/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonColumns.h"


StochPresolverSingletonColumns::StochPresolverSingletonColumns(PresolveData& presData, const sData& origProb)
   : StochPresolverBase(presData, origProb), removed_cols(0)
{
 // todo
}

StochPresolverSingletonColumns::~StochPresolverSingletonColumns()
{
 // todo
}


void StochPresolverSingletonColumns::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

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
      if( my_rank == 0)
         std::cout << "No more singletons left - exiting" << std::endl;
#endif
   }

   while( !presData.getSingletonCols().empty() )
   {
      bool removed = false;
      removed = removeSingletonColumn( presData.getSingletonCols().front().node, 
         presData.getSingletonCols().front().index );      
      
      if(removed && ( presData.getSingletonCols().front().node != -1 || my_rank == 0) )
         ++removed_cols;
      presData.getSingletonCols().pop();
   }

   assert( presData.getSingletonCols().empty() );

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "--- After singleton columns presolving:" << std::endl;
      std::cout << "--- Removed " << removed_cols << " singleton columns" << std::endl;
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

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());  
}

bool StochPresolverSingletonColumns::removeSingletonColumn(int node_col, int col)
{
   assert( -1 <= node_col && node_col < nChildren );
   updatePointersForCurrentNode(node_col, EQUALITY_SYSTEM);

   /* this should not happen - singeltons are collected locally */
   if( presData.nodeIsDummy(node_col, EQUALITY_SYSTEM) )
      assert( false );

   const SimpleVector& nnzs_col = (node_col == -1) ? *currNnzColParent : *currNnzColChild;
   if(!(nnzs_col[col] <= 1))
      std::cout << "sing col with " << nnzs_col[col] << " entries" << std::endl;
   assert( nnzs_col[col] <= 1 );
   if( nnzs_col[col] == 0 )
      return false;

   int node_row = -3;
   int row = -3;
   bool linking_row = false; 
   SystemType system_type = EQUALITY_SYSTEM;
   /* find the associated row via checking the transposed matices */
   bool found = findRowForColumnSingleton( system_type, node_row, row, linking_row, node_col, col );

   if( !found )
   {  // TODO
      std::cout << "node_col, col: " << node_col << ", " << col << std::endl; 
      assert( node_col == -1 );
      return false;
   }

   assert( -1 <= node_row && node_row < nChildren );
   updatePointersForCurrentNode( node_row, system_type );
   
   /* check whether col is free / implied free */
   bool lb_implied_free = false;
   bool ub_implied_free = false;

   checkColImpliedFree( system_type, node_row, row, linking_row, node_col, col, lb_implied_free, ub_implied_free);
   bool implied_free = lb_implied_free && ub_implied_free;

   /* if objective of variable is zero we can just remove it from the problem together with the containing row */
   double obj = (node_col == -1) ? (*currgParent)[col] : (*currgChild)[col];
   if( implied_free && PIPSisEQ(obj, 0.0) )
   {
      presData.removeImpliedFreeColumnSingleton( system_type, node_row, row, linking_row, node_col, col );
      return true;
   }  

   /* equalitiy singleton variables */
   if( system_type == EQUALITY_SYSTEM )
   {
      /* (originally) free singleton columns just get deleted together with their row */
      if( implied_free )
         presData.removeImpliedFreeColumnSingleton( system_type, node_row, row, linking_row, node_col, col );
      return true;
   }

   /* inequality singleton variables */
   // TODO

   return false;
}

bool StochPresolverSingletonColumns::findRowForColumnSingleton( SystemType& system_type, int& node_row, int& row, bool& linking,
   const int& node_col, const int& col )
{
   if( node_col == -1 )
   {
      /* go through all children and check linking variables for the singleton row */
      
      /* equality part */
      system_type = EQUALITY_SYSTEM;
      linking = false;

      for( node_row = -1; node_row < nChildren; ++node_row)
      {
         if( !presData.nodeIsDummy( node_row, system_type) )
         {
            updatePointersForCurrentNode( node_row, system_type );

            /* check transposed for entry */
            if( node_row == -1 )
            {
               /* Amat */
               assert(currAmatTrans);
               if(currAmatTrans->getRowPtr()[col].start != currAmatTrans->getRowPtr()[col].end)
               {
                  assert( (currAmatTrans->getRowPtr()[col].end - currAmatTrans->getRowPtr()[col].start) == 1);
                  row = currAmatTrans->getJcolM()[currAmatTrans->getRowPtr()[col].start];
                  return true;
               }
               /* Blmat */
               assert(currBlmatTrans);
               if(currBlmatTrans->getRowPtr()[col].start != currBlmatTrans->getRowPtr()[col].end)
               {
                  assert( (currBlmatTrans->getRowPtr()[col].end - currBlmatTrans->getRowPtr()[col].start) == 1);
                  row = currBlmatTrans->getJcolM()[currBlmatTrans->getRowPtr()[col].start];
                  linking = true;
                  return true;
               }
            }
            else
            {
               /* Amat */
               assert(currAmatTrans);
               if(currAmatTrans->getRowPtr()[col].start != currAmatTrans->getRowPtr()[col].end)
               {
                  assert( (currAmatTrans->getRowPtr()[col].end - currAmatTrans->getRowPtr()[col].start) == 1);
                  row = currAmatTrans->getJcolM()[currAmatTrans->getRowPtr()[col].start];
                  return true;
               }
            }
         }
      }

      /* inequality part */
      system_type = INEQUALITY_SYSTEM;

      for( node_row = -1; node_row < nChildren; ++ node_row)
      {
         if( !presData.nodeIsDummy( node_row, system_type) )
         {
            updatePointersForCurrentNode( node_row, system_type );

            /* check transposed for entry */
            if( node_row == -1 )
            {
               /* Amat */
               assert(currAmatTrans);
               if(currAmatTrans->getRowPtr()[col].start != currAmatTrans->getRowPtr()[col].end)
               {
                  assert( (currAmatTrans->getRowPtr()[col].end - currAmatTrans->getRowPtr()[col].start) == 1);
                  row = currAmatTrans->getJcolM()[currAmatTrans->getRowPtr()[col].start];
                  return true;
               }
               /* Blmat */
               assert(currBlmatTrans);
               if(currBlmatTrans->getRowPtr()[col].start != currBlmatTrans->getRowPtr()[col].end)
               {
                  assert( (currBlmatTrans->getRowPtr()[col].end - currBlmatTrans->getRowPtr()[col].start) == 1);
                  row = currBlmatTrans->getJcolM()[currBlmatTrans->getRowPtr()[col].start];
                  linking = true;
                  return true;
               }
            }
            else
            {
               /* Amat */
               assert(currAmatTrans);
               if(currAmatTrans->getRowPtr()[col].start != currAmatTrans->getRowPtr()[col].end)
               {
                  assert( (currAmatTrans->getRowPtr()[col].end - currAmatTrans->getRowPtr()[col].start) == 1);
                  row = currAmatTrans->getJcolM()[currAmatTrans->getRowPtr()[col].start];
                  return true;
               }
            }
         }
      }
   }
   else
   {
      assert( 0 <= node_col && node_col < nChildren );
      /* check Bmat and Blmat for the singleton row */

      system_type = EQUALITY_SYSTEM;
      node_row = node_col;
      /* equality part */
      if( !presData.nodeIsDummy( node_row, system_type) )
      {
         updatePointersForCurrentNode( node_row, system_type );

         /* Bmat */
         assert(currBmatTrans);
         if(currBmatTrans->getRowPtr()[col].start != currBmatTrans->getRowPtr()[col].end)
         {
            assert( (currBmatTrans->getRowPtr()[col].end - currBmatTrans->getRowPtr()[col].start) == 1);
            row = currBmatTrans->getJcolM()[currBmatTrans->getRowPtr()[col].start];
            return true;
         }
         /* Blmat */
         assert(currBlmatTrans);
         if(currBlmatTrans->getRowPtr()[col].start != currBlmatTrans->getRowPtr()[col].end)
         {
            assert( (currBlmatTrans->getRowPtr()[col].end - currBlmatTrans->getRowPtr()[col].start) == 1);
            row = currBlmatTrans->getJcolM()[currBlmatTrans->getRowPtr()[col].start];
            linking = true;
            return true;
         }
      }

      /* inequality part */
      system_type = INEQUALITY_SYSTEM;
      if( !presData.nodeIsDummy( node_row, system_type) )
      {
         updatePointersForCurrentNode( node_row, system_type );

         /* Bmat */
         assert(currBmatTrans);
         if(currBmatTrans->getRowPtr()[col].start != currBmatTrans->getRowPtr()[col].end)
         {
            assert( (currBmatTrans->getRowPtr()[col].end - currBmatTrans->getRowPtr()[col].start) == 1);
            row = currBmatTrans->getJcolM()[currBmatTrans->getRowPtr()[col].start];
            return true;
         }
         /* Blmat */
         assert(currBlmatTrans);
         if(currBlmatTrans->getRowPtr()[col].start != currBlmatTrans->getRowPtr()[col].end)
         {
            assert( (currBlmatTrans->getRowPtr()[col].end - currBlmatTrans->getRowPtr()[col].start) == 1);
            row = currBlmatTrans->getJcolM()[currBlmatTrans->getRowPtr()[col].start];
            linking = true;
            return true;
         }
      }
   }

   return false;
}

void StochPresolverSingletonColumns::checkColImpliedFree( SystemType system_type, int node_row, int row, bool linking_row, int node_col, int col, 
   bool& lb_implied_free, bool& ub_implied_free)
{
   ub_implied_free = lb_implied_free = false;

   updatePointersForCurrentNode( node_row, system_type );

   /* check whether originally free */
   const SimpleVector& ixupp_orig = (node_col == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*origProb.ixupp).vec) 
      : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*origProb.ixupp).children[node_col]->vec);
   const SimpleVector& ixlow_orig = (node_col == -1) ? dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*origProb.ixlow).vec) 
      : dynamic_cast<const SimpleVector&>(*dynamic_cast<const StochVector&>(*origProb.ixlow).children[node_col]->vec);

   if( ixupp_orig[col] == 0.0 )
      ub_implied_free = true;
   if( ixlow_orig[col] == 0.0 )
      lb_implied_free = true;

   /* check whether bound tightening found bounds from the variables row that make it implied free */
   ub_implied_free = ub_implied_free || presData.varBoundImpliedFreeBy( true, node_col, col, system_type, node_row, row, linking_row );
   lb_implied_free = lb_implied_free || presData.varBoundImpliedFreeBy( false, node_col, col, system_type, node_row, row, linking_row );
}