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
#include "StochPresolverTinyEntries.h"
#include "StochPresolverSingletonRows.h"
#include "StochPresolverDuplicateRows.h"
#include "StochPresolverSingletonColumns.h"
#include "PresolveData.h"

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
   myfile.open ("before.txt");
   sorigprob->writeToStreamDense(myfile);
   myfile.close();

   int nelims = 0;

   // init presolve data presData:
   PresolveData presData(sorigprob);
   presData.initialize();

   // init all presolvers:

   //StochPresolverDuplicateRows presolverDuplicateRow(presData);
   //StochPresolverTinyEntries presolverTiny(presData);
   StochPresolverSingletonRows presolverSR(presData);
   //StochPresolverSingletonColumns presolverSC(presData);


   // main presolving loop
   // while ( rerun )
      // presolverTiny.applyPresolving
      // presolverSingletonRow.applyPresolving

   bool possfeas;
   //possfeas = presolverSC.applyPresolving(nelims);
   //possfeas = presolverDuplicateRow.applyPresolving(nelims);
   //possfeas = presolverTiny.applyPresolving(nelims);
   possfeas = presolverSR.applyPresolving(nelims);
   //possfeas = presolverSC.applyPresolving(nelims);

  // if( !possfeas )
  //    break;

/*   cout<<"nRowElemsA "<<endl;
   nRowElemsA->writeToStreamAll(cout);
   cout<<"nRowElemsC "<<endl;
   nRowElemsC->writeToStreamAll(cout);
   cout<<"nColElems "<<endl;
   nColElems->writeToStreamAll(cout);
*/
   if( myRank == 0) cout<<"Finalizing presolved Data."<<endl;
   sData* finalPresData = presData.finalize();

   myfile.open("after.txt");
   finalPresData->writeToStreamDense(myfile);
   myfile.close();

   //presProb->writeToStreamDense(std::cout);

   std::cout << "sorigprob nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
   std::cout << "finalPresData nx, my, mz" << finalPresData->nx << " " << finalPresData->my << " " << finalPresData->mz << std::endl;

   //MPI_Barrier( MPI_COMM_WORLD);
   //assert(0);

   return finalPresData;
}


