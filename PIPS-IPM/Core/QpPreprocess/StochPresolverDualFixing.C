/*
 * StochPresolverDualFixing.C
 *
 *  Created on: 02.05.2018
 *      Author: bzfkempk
 */

#include "StochPresolverDualFixing.h"


StochPresolverDualFixing::StochPresolverDualFixing(PresolveData& pres_data, const sData& orig_prob) :
   StochPresolverBase(pres_data, orig_prob), fixed_columns(0)
{
}

StochPresolverDualFixing::~StochPresolverDualFixing()
{
}

bool StochPresolverDualFixing::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());

   int n_fixed_cols_run = 0;

#ifndef NDEBUG
   if( my_rank == 0 && verbosity > 1 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before dual fixing:" << std::endl;
   }
   countRowsCols();
#endif

   /* check root node */
   /* allreduce findings and remove synced */
   n_fixed_cols_run = applyDualFixingNode(-1);


   /* check all other nodes if not dummy */

   for( int i = 0; i < nChildren; ++i )
   {
      if( !presData.nodeIsDummy(i) )
      {
         n_fixed_cols_run += applyDualFixingNode(i);
      }
   }


#ifndef NDEBUG
   if( my_rank == 0 && verbosity > 1 )
   {
      std::cout << "\tFixed columns in dual fixing: " << fixed_columns << std::endl;
      std::cout << "--- After dual fixing:" << std::endl;
   }
   else if( my_rank == 0 && verbosity == 1)
      std::cout << "DuFix:\t removed " << fixed_columns << " cols" << std::endl;

   countRowsCols();
   if( my_rank == 0 && verbosity > 1 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.presDataInSync());

   if( n_fixed_cols_run != 0 )
      return true;
   else
      return false;

}



int StochPresolverDualFixing::applyDualFixingNode(int node)
{


   return 0;
}
