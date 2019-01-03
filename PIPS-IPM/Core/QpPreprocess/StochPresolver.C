/*
 * StochPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */


#include "StochPresolver.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <utility>
#include <math.h>
#include <algorithm>

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

   if( myRank == 0) std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

   // clone and initialize dynamic storage

   ofstream myfile;
   /*myfile.open ("before.txt");
   sorigprob->writeToStreamDense(myfile);
   myfile.close();*/

   // init presolve data presData:
   PresolveData presData(sorigprob);
   presData.initialize();

   // init all presolvers:

   StochPresolverBoundStrengthening presolverBS(presData);
   StochPresolverParallelRows presolverParallelRow(presData);
   StochPresolverModelCleanup presolverCleanup(presData);
   StochPresolverSingletonRows presolverSR(presData);

#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- Before Presolving:"<<endl;
   presolverSR.countRowsCols();
#endif

   presolverSR.applyPresolving();
   presolverBS.applyPresolving();
   presolverCleanup.applyPresolving();
   presolverParallelRow.applyPresolving();
   presolverSR.applyPresolving();


/*   cout<<"nRowElemsA "<<endl;
   nRowElemsA->writeToStreamAll(cout);
   cout<<"nRowElemsC "<<endl;
   nRowElemsC->writeToStreamAll(cout);
   cout<<"nColElems "<<endl;
   nColElems->writeToStreamAll(cout);
*/
#ifndef NDEBUG
   assert( presolverSR.verifyNnzcounters() );
   if( myRank == 0 )
      cout<<"--- After Presolving:"<<endl;
   presolverSR.countRowsCols();
#endif
   //if( myRank == 0) cout<<"Finalizing presolved Data."<<endl;
   sData* finalPresData = presData.finalize();

   /*myfile.open("after.txt");
   finalPresData->writeToStreamDense(myfile);
   myfile.close();*/

   if( myRank==0 )
   {
      std::cout << "original problem: variables, equ. constraints., inequ. constraints" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
      std::cout << "presolved problem: variables, equ. constraints., inequ. constraints " << finalPresData->nx << " " << finalPresData->my << " " << finalPresData->mz << std::endl;
   }
#ifdef TIMING
   std::cout << "sorigprob nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
   std::cout << "finalPresData nx, my, mz" << finalPresData->nx << " " << finalPresData->my << " " << finalPresData->mz << std::endl;
#endif

   return finalPresData;
}


