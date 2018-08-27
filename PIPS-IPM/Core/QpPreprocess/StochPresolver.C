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
   //StochPresolverSingletonColumns presolverSC(presData);

#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- Before Presolving:"<<endl;
   presolverSR.countRowsCols();
#endif

   // main presolving loop
   // while ( rerun )
      // presolverTiny.applyPresolving
      // presolverSingletonRow.applyPresolving

   presolverBS.applyPresolving();
   presolverParallelRow.applyPresolving();
   presolverCleanup.applyPresolving();
   presolverSR.applyPresolving();
   //presolverSC.applyPresolving();

/*   cout<<"nRowElemsA "<<endl;
   nRowElemsA->writeToStreamAll(cout);
   cout<<"nRowElemsC "<<endl;
   nRowElemsC->writeToStreamAll(cout);
   cout<<"nColElems "<<endl;
   nColElems->writeToStreamAll(cout);
*/
#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- After Presolving:"<<endl;
   presolverSR.countRowsCols();
#endif
   if( myRank == 0) cout<<"Finalizing presolved Data."<<endl;
   sData* finalPresData = presData.finalize();

   /*myfile.open("after.txt");
   finalPresData->writeToStreamDense(myfile);
   myfile.close();*/

   //presProb->writeToStreamDense(std::cout);

   std::cout << "sorigprob nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
   std::cout << "finalPresData nx, my, mz" << finalPresData->nx << " " << finalPresData->my << " " << finalPresData->mz << std::endl;

   //MPI_Barrier( MPI_COMM_WORLD);
   //assert(0);

   return finalPresData;
}


