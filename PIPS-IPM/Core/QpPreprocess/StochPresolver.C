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
#include "StochPresolverSingletonRows.h"
#include "StochPresolverSingletonColumns.h"
#include "PresolveData.h"
#include "StochPresolverParallelRows.h"
#include "StochPresolverBoundStrengthening.h"
#include "StochPresolverModelCleanup.h"
#include "pipschecks.h"

StochPresolver::StochPresolver(const Data* prob)
 : QpPresolver(prob)
{
   // todo
}

StochPresolver::~StochPresolver()
{
   // todo
}

Data* StochPresolver::presolve()
{
   int myRank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( myRank == 0 )
      std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

#if 0 // todo add flag
   ofstream myfile;
   myfile.open ("before_presolving.txt");
   sorigprob->writeToStreamDense(myfile);
   myfile.close();
#endif

   /* initialize presolve data */
   PresolveData presData(sorigprob);

   assert( sorigprob->isRootNodeInSync());
   assert( presData.presProb->isRootNodeInSync() );

   /* initialize all presolvers */
   StochPresolverBoundStrengthening presolverBS(presData);
   StochPresolverParallelRows presolverParallelRow(presData);
   StochPresolverModelCleanup presolverCleanup(presData);
   StochPresolverSingletonRows presolverSR(presData);

   if( myRank == 0 )
      std::cout <<"--- Before Presolving: " << std::endl;
   presolverSR.countRowsCols();

   // todo loop, and not exhaustive
   // some list holding all presolvers - eg one presolving run
   // some while iterating over the list over and over until either every presolver says im done or some iterlimit is reached?
   for( int i = 0; i < 1; ++i )
   {
      /* singleton rows */
//      presolverCleanup.applyPresolving();
      presolverSR.applyPresolving();
//      presolverBS.applyPresolving();
//      presolverParallelRow.applyPresolving();
//      presolverCleanup.applyPresolving();
   }

   if( myRank == 0 )
      std::cout << "--- After Presolving:" << std::endl;
   presolverCleanup.countRowsCols();

   {
   StochVector& g  = dynamic_cast<StochVector&>(*dynamic_cast<sData*>(presData.presProb)->g);
   double obj_inf = g.infnorm();
   double obj_min = 0;
   int obj_min_idx = -1;
   double obj_max = 0;
   int obj_max_idx = -1;
   g.min(obj_min, obj_min_idx);
   g.max(obj_max, obj_max_idx);

   if( myRank == 0 )
      std::cout << "obj_inf: " << obj_inf << "\tobj_min: " << obj_min << " at " << obj_min_idx << "\tobj_max: " << obj_max << " at " << obj_max_idx << std::endl;
   }

   // todo not yet available for dynamic storage
   assert( presData.presProb->isRootNodeInSync() );

//      presData.presProb->writeToStreamDense(std::cout);
   // i assume we actually apply our changes here and then return a valid sData object to the caller
   sData* finalPresData = presData.finalize();

//   finalPresData->writeToStreamDense(std::cout);

   assert( finalPresData->isRootNodeInSync() );

//   exit(1);

   {
      StochVector& g  = dynamic_cast<StochVector&>(*dynamic_cast<sData*>(finalPresData)->g);
      double obj_inf = g.infnorm();
      double obj_min = 0;
      int obj_min_idx = -1;
      double obj_max = 0;
      int obj_max_idx = -1;
      g.min(obj_min, obj_min_idx);
      g.max(obj_max, obj_max_idx);

      if( myRank == 0 )
         std::cout << "obj_inf: " << obj_inf << "\tobj_min: " << obj_min << " at " << obj_min_idx << "\tobj_max: " << obj_max << " at " << obj_max_idx << std::endl;
   }

   return finalPresData;
}


