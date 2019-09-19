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

   if( presData.getSingletonRows().size() == 0 )
   {
#ifndef NDEBUG
      if( my_rank == 0)
         std::cout << "No more singletons left - exiting" << std::endl;
#endif
   }


   while( !presData.getSingletonCols().empty() )
   {
      //todo
      bool removed = false;
      //removed = removeSingletonCol( presData.getSingletonRows().back().system_type, presData.getSingletonRows().back().node, 
      //   presData.getSingletonRows().back().index );      
      
      if(removed)
         ++removed_cols;
      presData.getSingletonCols().pop_back();
   }

   presData.allreduceObjOffset();
   if( my_rank == 0 )
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;

#ifndef NDEBUG
   if( my_rank == 0 )
      std::cout << "--- After singleton columns presolving:" << std::endl;
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

bool StochPresolverSingletonColumns::removeSingletonColumn(int node, int col)
{


   return false;

}


