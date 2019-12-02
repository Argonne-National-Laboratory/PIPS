/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonColumns.h"
#include "StochVectorUtilities.h"

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
      if( my_rank == 0)
         std::cout << "No more singletons left - exiting" << std::endl;
#endif
   }

   // main loop:
   while( !presData.getSingletonCols().empty() )
   {
      bool removed = false;
      removed = removeSingletonColumn( presData.getSingletonCols().front().node, 
         presData.getSingletonCols().front().index );      
      
      if(removed && (presData.getSingletonCols().front().node != -1 || my_rank == 0) )
         ++removed_cols;
      presData.getSingletonCols().pop();
   }

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
   presData.allreduceAndApplyObjVecChanges();
   presData.allreduceObjOffset();

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyNnzcounters());
   assert(presData.verifyActivities());  
}

bool StochPresolverSingletonColumns::removeSingletonColumn(const int& node_col, int col)
{
   assert( -1 <= node_col && node_col < nChildren );
   assert( !presData.nodeIsDummy(node_col) );

   updatePointersForCurrentNode(node_col, EQUALITY_SYSTEM);
   const SimpleVectorBase<int>& nnzs_col = (node_col == -1) ? *currNnzColParent : *currNnzColChild;

   if( nnzs_col[col] == 0 )
      return false;
   assert( nnzs_col[col] == 1 );

   int node_row = -3;
   int row = -3;
   bool linking_row = false; 
   SystemType system_type = EQUALITY_SYSTEM;
   
   /* find the associated row via checking the transposed matices */
   bool found = findRowForColumnSingleton( system_type, node_row, row, linking_row, node_col, col );

   if( linking_row )
      return false;

   if( !found )
   {  
      // TODO : singleton columns in linking parts need communication
      assert( node_col == -1 );
      return false;
   }

   assert( -1 <= node_row && node_row < nChildren );
   updatePointersForCurrentNode( node_row, system_type );
   
   /* check whether col is free / implied free */
   bool lb_implied_free = false;
   bool ub_implied_free = false;

   // todo : pivot numerical limits

   checkColImpliedFree( system_type, node_row, row, linking_row, node_col, col, lb_implied_free, ub_implied_free);
   bool implied_free = lb_implied_free && ub_implied_free;

   /* if objective of variable is zero we can just remove it from the problem together with the containing row */
   double obj = (node_col == -1) ? (*currgParent)[col] : (*currgChild)[col];
   if( implied_free && PIPSisEQ(obj, 0.0) )
   {
      presData.writeRowLocalToStreamDense(std::cout, system_type, node_row, linking_row, row);
      std::cout << node_col << "\taa\t" << col << std::endl;
      presData.removeImpliedFreeColumnSingleton( system_type, node_row, row, linking_row, node_col, col );
      return true;
   }  

   /* equalitiy singleton variables */
   if( system_type == EQUALITY_SYSTEM )
   {
      /* (originally) free singleton columns just get deleted together with their row */
      if( implied_free )
      {
         presData.writeRowLocalToStreamDense(std::cout, system_type, node_row, linking_row, row);
         std::cout << node_col << "\tbb\t" << col << std::endl;
         presData.removeImpliedFreeColumnSingleton( system_type, node_row, row, linking_row, node_col, col );
         return true;
      }
   }

   /* inequality singleton variables */
   // TODO

   return false;
}

bool StochPresolverSingletonColumns::findRowForColumnSingleton( SystemType& system_type, int& node_row, int& row, bool& linking,
   const int& node_col, const int& col )
{
   assert( -1 <= node_col && node_col < nChildren );
   if( node_col == -1 )
   {
      if( findRowForLinkingSingleton( system_type, node_row, row, linking, col) )
         return true;
   }
   else
   {
      node_row = node_col;
      if( findRowForNonlinkingSingelton( system_type, row, linking, node_col, col ) )
         return true;      
   }
   
   return false;
}

bool StochPresolverSingletonColumns::findRowForLinkingSingleton( SystemType& system_type, int& node_row, int& row, bool& linking, const int& col)
{
      /* go through all children and check linking variables for the singleton row */
      
      /* equality part */
      system_type = EQUALITY_SYSTEM;
      if( findRowForLinkingSingletonInSystem( system_type, node_row, row, linking, col) )
         return true;
      
      /* inequality part */
      system_type = INEQUALITY_SYSTEM;
      if( findRowForLinkingSingletonInSystem( system_type, node_row, row, linking, col) )
         return true;

      return false;
}

bool StochPresolverSingletonColumns::findRowForNonlinkingSingelton( SystemType& system_type, int& row, bool& linking, const int& node_col, const int& col)
{
   if( presData.nodeIsDummy( node_col ) )
      return false;

   /* check Bmat and Blmat for the singleton row */
   system_type = EQUALITY_SYSTEM;
   linking = false;

   /* equality part */
   updatePointersForCurrentNode( node_col, system_type );
   
   /* Bmat */
   assert(currBmatTrans);
   if( findRowForSingletonColumnInMatrix( *currBmatTrans, row, col) )
      return true;

   /* Blmat */
   assert(currBlmatTrans);
   if( findRowForSingletonColumnInMatrix( *currBlmatTrans, row, col) )
   {
      linking = true;
      return true;
   }
   

   /* inequality part */
   system_type = INEQUALITY_SYSTEM;
   updatePointersForCurrentNode( node_col, system_type );

   /* Bmat */
   assert(currBmatTrans);
   if( findRowForSingletonColumnInMatrix( *currBmatTrans, row, col) )
      return true;

   /* Blmat */
   assert(currBlmatTrans);
   if( findRowForSingletonColumnInMatrix( *currBlmatTrans, row, col) )
   {
      linking = true;
      return true;
   }

   return false;
}

bool StochPresolverSingletonColumns::findRowForLinkingSingletonInSystem( SystemType system_type, int& node_row, int& row, bool& linking, const int& col)
{
   linking = false;

   for( node_row = -1; node_row < nChildren; ++node_row)
   {
      if( !presData.nodeIsDummy( node_row ) )
      {
         updatePointersForCurrentNode( node_row, system_type );

         /* check transposed for entry */
         if( node_row == -1 )
         {
            /* Bmat */
            assert(currBmatTrans);
            if( findRowForSingletonColumnInMatrix( *currBmatTrans, row, col) )
               return true;

            /* Blmat */
            assert(currBlmatTrans);
            if( findRowForSingletonColumnInMatrix( *currBlmatTrans, row, col) )
            {
               linking = true;
               return true;
            }
         }
         else
         {
            /* Amat */
            assert(currAmatTrans);
            if( findRowForSingletonColumnInMatrix( *currAmatTrans, row, col) )
               return true;
         }
      }
   }
   return false;
}


bool StochPresolverSingletonColumns::findRowForSingletonColumnInMatrix( const SparseStorageDynamic& mat, int& row, const int& col)
{
   assert( col < mat.getM() );
   if( mat.getRowPtr(col).start != mat.getRowPtr(col).end )
   {
      assert( (mat.getRowPtr(col).end - mat.getRowPtr(col).start) == 1);
      row = mat.getJcolM(mat.getRowPtr(col).start);
      return true;
   }
   return false;
}


void StochPresolverSingletonColumns::checkColImpliedFree( SystemType system_type, int node_row, int row, bool linking_row, int node_col, int col, 
   bool& lb_implied_free, bool& ub_implied_free)
{
   ub_implied_free = lb_implied_free = false;

   updatePointersForCurrentNode( node_row, system_type );

   /* check whether originally free */
   const double ixupp_orig = getSimpleVecFromColStochVec(*origProb.ixupp, node_col)[col];
   const double ixlow_orig = getSimpleVecFromColStochVec(*origProb.ixlow, node_col)[col];

   if( PIPSisZero(ixupp_orig) )
      ub_implied_free = true;
   if( PIPSisZero(ixlow_orig) )
      lb_implied_free = true;

   /* check whether bound tightening found bounds from the variables row that make it implied free */
   // ub_implied_free = ub_implied_free || presData.varBoundImpliedFreeBy( true, node_col, col, system_type, node_row, row, linking_row );
   // lb_implied_free = lb_implied_free || presData.varBoundImpliedFreeBy( false, node_col, col, system_type, node_row, row, linking_row );
}