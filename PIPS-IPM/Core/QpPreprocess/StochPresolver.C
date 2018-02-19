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

#include "../Abstract/DoubleMatrix.h"
#include "../SparseLinearAlgebra/SparseGenMatrix.h"
#include "../StochLinearAlgebra/StochVectorHandle.h"
#include "../Vector/OoqpVector.h"
#include "StochGenMatrix.h"
#include "sTreeCallbacks.h"

StochPresolver::StochPresolver(const Data* prob)
 : QpPresolver(prob)
{
   const sData* sorigprob = dynamic_cast<const sData*>(prob);

   StochVectorHandle gclone(dynamic_cast<StochVector*>(sorigprob->g->clone()));
   nColElems = gclone;

   StochVectorHandle bAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   nRowElemsA = bAclone;

   StochVectorHandle icuppclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   nRowElemsC = icuppclone;

   presProb = NULL;
}


StochPresolver::~StochPresolver()
{
}

void
StochPresolver::initNnzCounter()
{
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   StochVectorHandle colClone(dynamic_cast<StochVector*>(nColElems->clone()));

   A.getNnzPerRow(*nRowElemsA);
   C.getNnzPerRow(*nRowElemsC);

   A.getNnzPerCol(*nColElems);
   C.getNnzPerCol(*colClone);

   nColElems->axpy(1.0, *colClone);

#if 0
   int rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if( rank == 0 ) std::cout << "write Cols with all cols" << std::endl;
   nColElems->writeToStreamAll(std::cout);

   if( rank == 0 ) std::cout << "write Rows A" << std::endl;
   nRowElemsA->writeToStreamAll(std::cout);

   if( rank == 0 ) std::cout << "write Rows C " << std::endl;
   nRowElemsC->writeToStreamAll(std::cout);
#endif


}

Data* StochPresolver::presolve()
{

   std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

   dynamic_cast<sTreeCallbacks*>(sorigprob->stochNode)->writeSizes(std::cout);


   // clone and initialize dynamic storage
   presProb = sorigprob->cloneFull(true);

   int rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   ofstream myfile;
   myfile.open ("before.txt");
   sorigprob->writeToStreamDense(myfile);
   myfile.close();

   StochGenMatrix& Apres = dynamic_cast<StochGenMatrix&>(*presProb->A);
   StochGenMatrix& Cpres = dynamic_cast<StochGenMatrix&>(*presProb->C);

   // initialized all dynamic transposed sub matrices
   Apres.initTransposed(true);
   Cpres.initTransposed(true);


   initNnzCounter();


   // todo: presolve
   //removeTinyEntries();


   presProb->cleanUpPresolvedData(*nRowElemsA, *nRowElemsC, *nColElems);


   std::cout << "AFTER \n" << std::endl;

   dynamic_cast<sTreeCallbacks*>(presProb->stochNode)->writeSizes(std::cout);

   Apres.deleteTransposed();
   Cpres.deleteTransposed();

   myfile.open("after.txt");
   presProb->writeToStreamDense(myfile);
   myfile.close();

   presProb->writeToStreamDense(std::cout);

   std::cout << "nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;

   //assert(0);

   return presProb;
}

void StochPresolver::removeTinyEntries()
{
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochVector& bA = dynamic_cast<StochVector&>(*(presProb->bA));

   StochVector& reductionsCol = *nColElems->clone();

   // remove tiny entries in A row-wise
   //removeTinyEntriesA(reductionsCol);

   // remove tiny entries in C row-wise
   removeTinyEntriesC(reductionsCol);


   // todo StochGenDummyMatrix
}


void StochPresolver::removeTinyEntriesC( StochVector& reductionsCol)
{
   double* cuppElems = (dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec))->elements();
   double* clowElems = (dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vec))->elements();
   double* icuppElems = (dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vec))->elements();
   double* iclowElems = (dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vec))->elements();

   double* const xuppElems = (dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presProb->bux)).vec))->elements();
   double* const xlowElems = (dynamic_cast<SimpleVector*>(dynamic_cast<const StochVector&>(*(presProb->blx)).vec))->elements();

   SimpleVector* nnzPerRowC = dynamic_cast<SimpleVector*>(nRowElemsC->vec);
   SimpleVector* nnzPerColC = dynamic_cast<SimpleVector*>(nColElems->vec);
   SimpleVector* reductionsRowSimple(nnzPerRowC);
   reductionsRowSimple->setToZero();
   SimpleVector* reductionsColSimple(nnzPerColC);
   reductionsColSimple->setToZero();
   double* reductRow = reductionsRowSimple->elements();
   double* reductCol = reductionsColSimple->elements();

   SparseStorageDynamic& storage = dynamic_cast<StochGenMatrix&>(*(presProb->C)).Bmat->getStorageDynamicRef();

   // todo: check correct sizes

   for( int r = 0; r < storage.m; r++ )
   {
      const int start = storage.rowptr[r].start;
      const int end = storage.rowptr[r].end;
      for( int k = start; k < end; k++ )
      {
         cout<<"consider entry ( "<<r<<", "<<storage.jcolM[k]<<" ) = "<<storage.M[k]<<endl;
         if( fabs(storage.M[k])<1e-3 && fabs(storage.M[k]) * (xuppElems[k]-xlowElems[k]) * nnzPerRowC->elements()[r] < 1e-2 * feastol )
         {
            if( icuppElems[r] != 0.0)
               cuppElems[r] -= storage.M[k] * xlowElems[k];
            if( iclowElems[r] != 0.0 )
               clowElems[r] -= storage.M[k] * xlowElems[k];
            storage.M[k] = 0.0;
            reductRow[r]++;
            reductCol[storage.jcolM[k]]++;
            cout<<"entry ( "<<r<<", "<<storage.jcolM[k]<<" ) removed"<<endl;
         }
         if(fabs(storage.M[k]) < 1e-10 ){
            storage.M[k] = 0.0;
            reductRow[r]++;
            reductCol[storage.jcolM[k]]++;
            cout<<"entry ( "<<r<<", "<<storage.jcolM[k]<<" ) removed"<<endl;
         }
      }
   }
}

