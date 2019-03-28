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
//#include <math.h>
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
}

StochPresolver::~StochPresolver()
{
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

   assert( rootNodeInSyncSData(*sorigprob) );


   /* initialize all presolvers */
   StochPresolverBoundStrengthening presolverBS(presData);
   StochPresolverParallelRows presolverParallelRow(presData);
   StochPresolverModelCleanup presolverCleanup(presData);
   StochPresolverSingletonRows presolverSR(presData);

#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- Before Presolving:"<<endl;
   presolverSR.countRowsCols();
#endif

   // todo loop, and not exhaustive
   // some list holding all presolvers - eg one presolving run
   // some while iterating over the list over and over until either every presolver says im done or some iterlimit is reached?
   for( int i = 0; i < 1; ++i )
   {
      presolverSR.applyPresolving();
      //presolverBS.applyPresolving();
      // TODO bugged
      presolverCleanup.applyPresolving();
      //presolverParallelRow.applyPresolving();
      //presolverSR.applyPresolving();
   }

#ifndef NDEBUG
   assert( presolverSR.verifyNnzcounters() );
   if( myRank == 0 )
      std::cout << "--- After Presolving:" << std::endl;
   presolverSR.countRowsCols();
#endif


   // i assume we actually apply ur changes here and then return a valid sData object to the caller
   sData* finalPresData = presData.finalize();


   assert( rootNodeInSyncSData(*finalPresData) );

#if 0
   myfile.open("after_presolving.txt");
   finalPresData->writeToStreamDense(myfile);
   myfile.close();
#endif

   if( myRank==0 )
   {
      std::cout << "original problem:\t" << sorigprob->nx << " variables\t" << sorigprob->my << " equ. conss\t" << sorigprob->mz << " ineq. conss" << std::endl;
      std::cout << "presolved problem:\t" << finalPresData->nx << " variables\t" << finalPresData->my << " equ. conss\t" << finalPresData->mz << " ineq. conss" << std::endl;
   }
#ifdef TIMING
   std::cout << "sorigprob nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
   std::cout << "finalPresData nx, my, mz" << finalPresData->nx << " " << finalPresData->my << " " << finalPresData->mz << std::endl;
#endif

   return finalPresData;
}


