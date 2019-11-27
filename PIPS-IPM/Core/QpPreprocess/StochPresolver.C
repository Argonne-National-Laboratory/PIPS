/*
 * StochPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#include "StochPresolver.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctype.h>

#include "StochVector.h"
#include "StochGenMatrix.h"
#include "SmartPointer.h"
#include "sData.h"
#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include "StochVectorHandle.h"
#include "OoqpVector.h"
#include "StochGenMatrix.h"
#include "sTreeCallbacks.h"
#include "DoubleMatrixTypes.h"
#include "PresolveData.h"
#include "StochPostsolver.h"
#include "StochPresolverBoundStrengthening.h"
#include "StochPresolverModelCleanup.h"
#include "StochPresolverColumnFixation.h"
#include "StochPresolverSingletonRows.h"
//#include "StochPresolverSingletonColumns.h"
#include "StochPresolverParallelRows.h"
#include "pipschecks.h"


StochPresolver::StochPresolver(const Data* prob, Postsolver* postsolver = NULL)
 : QpPresolver(prob, postsolver)
{
   // todo
}

StochPresolver::~StochPresolver()
{
   // todo
}

Data* StochPresolver::presolve()
{
   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   if( my_rank == 0 )
      std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

   /* initialize presolve data */
   PresolveData presData(sorigprob, dynamic_cast<StochPostsolver*>(postsolver));

   assert( sorigprob->isRootNodeInSync() );
   assert( presData.getPresProb().isRootNodeInSync() );

   /* initialize all presolvers */
   StochPresolverBoundStrengthening presolverBS(presData, *sorigprob);
   StochPresolverParallelRows presolverParallelRow(presData, *sorigprob);
   StochPresolverModelCleanup presolverCleanup(presData, *sorigprob);
   StochPresolverColumnFixation presolverColFix(presData, *sorigprob);
   StochPresolverSingletonRows presolverSR(presData, *sorigprob);

   if( my_rank == 0 )
      std::cout <<"--- Before Presolving: " << std::endl;
   presolverCleanup.countRowsCols();

   // todo loop, and exhaustive
   // some list holding all presolvers - eg one presolving run
   // some while iterating over the list over and over until either every presolver says im done or some iterlimit is reached?
   presolverCleanup.applyPresolving();
   
   for( int i = 0; i < 1; ++i )
   {
      /* singleton rows */
      presolverSR.applyPresolving();
      presolverBS.applyPresolving();
      presolverParallelRow.applyPresolving();
      presolverColFix.applyPresolving();
   }

   // before the finalize call fix all empty rows and columns not yet fixed
   presolverCleanup.applyPresolving();
   
   if( my_rank == 0 )
      std::cout << "--- After Presolving:" << std::endl;
   presolverCleanup.countRowsCols();
   assert( presData.getPresProb().isRootNodeInSync() );

   // exit(1);


// todo : tell postsolver aboud released variables
   char* env = getenv("PIPS_RESET_FREE_VARIABLES");
   if( env != NULL )
   {
      std::string reset_vars(env);
      for(unsigned int i = 0; i < reset_vars.length(); ++i)
         reset_vars[i] = std::tolower(reset_vars[i]);
      //std::transform(reset_vars.begin(), reset_vars.end(), reset_vars.begin(), [](unsigned char c){ return std::tolower(c); }); 
      if(reset_vars == "true")
      {
          if( my_rank == 0 )
             std::cout << "Resetting bounds found in bound strengthening" << std::endl;

          presData.resetOriginallyFreeVarsBounds(*sorigprob);
          presolverCleanup.countRowsCols();
      }
   }   

   sData* finalPresData = presData.finalize();

   assert( finalPresData );
   assert( finalPresData->isRootNodeInSync() );

   // todo : verify presolver and postsolver have same amount of deleted rows and cols


   return finalPresData;
}
