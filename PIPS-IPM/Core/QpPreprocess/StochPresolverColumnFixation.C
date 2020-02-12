/*
 * StochPresolverColumnFixation.h
 *
 *  Created on: 15.05.2019
 *      Author: bzfkempk
 */

#include "pipsdef.h"
#include "StochPresolverColumnFixation.h"

#include <cmath>
#include <limits>

StochPresolverColumnFixation::StochPresolverColumnFixation(
      PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb), fixed_columns(0)
{
}

StochPresolverColumnFixation::~StochPresolverColumnFixation()
{
}

// todo atm can only be called once : does not recognize already fixed columns

/* scan through columns and fix those that have tight bounds */
void StochPresolverColumnFixation::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyActivities());
   assert(presData.verifyNnzcounters());

#ifndef NDEBUG
   if( my_rank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before column fixation presolving:" << std::endl;
   }
   countRowsCols();
#endif

   int fixed_columns_run = 0;

   /* remove fixed columns from system */
   updatePointersForCurrentNode(-1, EQUALITY_SYSTEM);

   /* linking variables */
   for( int col = 0; col < currgParent->n; ++col )
   {
      const INDEX col_index = INDEX(COL, -1, col);

      if(presData.wasColumnRemoved(col_index))
         continue;

      // todo : make criterion numerically more stable
      // if the column is fixed to a value I want abs(xlow - xupp) to be low enough such that
      // in the whole column for each aij != 0 abs(xlow - xupp) * abs(aij) < epsilon <= feastol
      // this way we ensure that even if the fixing is not quite accurate the column stays valid
      // todo : include also objective coefficient into this criterion ? - so that fixing the column somehow inaccurately will not majorly impact the objective value
      double absmax_row = 1.0;

      if( !PIPSisZero((*currIxlowParent)[col]) && !PIPSisZero((*currIxuppParent)[col]))
      {

         assert(PIPSisLE(0.0, (*currxuppParent)[col] - (*currxlowParent)[col], tolerance4));

         if( PIPSisLT( fabs((*currxuppParent)[col] - (*currxlowParent)[col]) * absmax_row, tolerance4) )
         {
            // verify if one of the bounds is integer:
            double intpart = 0.0;
            double value = 0.0;
            if( std::modf( (*currxlowParent)[col], &intpart) == 0.0 )
               value = (*currxlowParent)[col];
            else if( std::modf( (*currxuppParent)[col], &intpart) == 0.0 )
               value = (*currxuppParent)[col];
            else  // set the variable to the arithmetic mean:
               value = ( (*currxlowParent)[col] + (*currxuppParent)[col] ) / 2.0;

            presData.fixColumn( col_index, value);
            if(my_rank == 0)
               ++fixed_columns_run;
         }
      }
   }

   /* child nodes */
   for(int node = 0; node < nChildren; ++node)
   {
      if(presData.nodeIsDummy(node) )
         continue;

      updatePointersForCurrentNode(node, EQUALITY_SYSTEM);

      /* linking variables */
      for( int col = 0; col < currgChild->n; ++col )
      {
         // todo : make criterion numerically more stable
         // if the column is fixed to a value I want abs(xlow - xupp) to be low enough such that
         // in the whole column for each aij != 0 abs(xlow - xupp) * abs(aij) < epsilon <= feastol
         // this way we ensure that even if the fixing is not quite accurate the column stays valid
         // todo : include also objective coefficient into this criterion ?
         double absmax_row = 1.0;

         if( !PIPSisZero((*currIxlowChild)[col]) && !PIPSisZero((*currIxuppChild)[col]) )
         {
            assert(PIPSisLE(0.0, (*currxuppChild)[col] - (*currxlowChild)[col]));

            if( PIPSisLT( ((*currxuppChild)[col] - (*currxlowChild)[col]) * absmax_row, tolerance4) )
            {
               // verify if one of the bounds is integer:
               double intpart = 0.0;
               double value = 0.0;
               if( std::modf( (*currxlowChild)[col], &intpart) == 0.0 )
                  value = (*currxlowChild)[col];
               else if( std::modf( (*currxuppChild)[col], &intpart) == 0.0 )
                  value = (*currxuppChild)[col];
               else  // set the variable to the arithmetic mean:
                  value = ( (*currxlowChild)[col] + (*currxuppChild)[col] ) / 2.0;
               presData.fixColumn( INDEX(COL, node, col), value);
               ++fixed_columns_run;
            }
         }
      }
   }

   /* communicate the local changes */
   presData.allreduceAndApplyLinkingRowActivities();
   presData.allreduceAndApplyNnzChanges();
   presData.allreduceObjOffset();


#ifndef NDEBUG
   PIPS_MPIgetSumInPlace(fixed_columns_run, MPI_COMM_WORLD);
   fixed_columns += fixed_columns_run;
   if( my_rank == 0 )
      std::cout << "\tFixed columns during column fixation: " << fixed_columns << std::endl;
   if( my_rank == 0 )
      std::cout << "--- After column fixation presolving:" << std::endl;
   countRowsCols();
   if( my_rank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.getPresProb().isRootNodeInSync());
   assert(presData.verifyActivities());
   assert(presData.verifyNnzcounters());
}
