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
#include "PresolveData.h"
#include "StochPostsolver.h"
#include "StochPresolverBoundStrengthening.h"
#include "StochPresolverModelCleanup.h"
#include "StochPresolverColumnFixation.h"
// #include "StochPresolverSingletonRows.h"
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
   PresolveData presData(sorigprob, dynamic_cast<StochPostsolver*>(postsolver));

   assert( sorigprob->isRootNodeInSync());
   assert( presData.getPresProb().isRootNodeInSync() );

   /* initialize all presolvers */
   StochPresolverBoundStrengthening presolverBS(presData, *sorigprob);
   StochPresolverParallelRows presolverParallelRow(presData, *sorigprob);
   StochPresolverModelCleanup presolverCleanup(presData, *sorigprob);
   StochPresolverColumnFixation presolverColFix(presData, *sorigprob);
//   StochPresolverSingletonRows presolverSR(presData, *sorigprob);

   if( myRank == 0 )
      std::cout <<"--- Before Presolving: " << std::endl;
   presolverCleanup.countRowsCols();

   // todo loop, and exhaustive
   // some list holding all presolvers - eg one presolving run
   // some while iterating over the list over and over until either every presolver says im done or some iterlimit is reached?
   for( int i = 0; i < 1; ++i )
   {
      /* singleton rows */
     // presolverSR.applyPresolving();
      // presolverBS.applyPresolving();
      // presolverColFix.applyPresolving();
      // presolverCleanup.applyPresolving();
//      presolverColFix.applyPresolving();
//      presolverBS.applyPresolving();

//      presolverBS.applyPresolving();

     presolverParallelRow.applyPresolving();
   }

   if( myRank == 0 )
      std::cout << "--- After Presolving:" << std::endl;
   presolverCleanup.countRowsCols();
   assert( presData.getPresProb().isRootNodeInSync() );
//      presData.presProb->writeToStreamDense(std::cout);

   sData* finalPresData = presData.finalize();
//   finalPresData->writeToStreamDense(std::cout);
   assert( finalPresData->isRootNodeInSync() );
   // exit(1);

   return finalPresData;
}
