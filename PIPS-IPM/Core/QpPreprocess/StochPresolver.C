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
   // remove tiny entries in A

   // remove tiny entries in C
   removeTinyEntriesC();
}


void StochPresolver::removeTinyEntriesC()
{
   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1 )
      iAmDistrib = true;

   if( iAmDistrib && rank > 0 )
   {
      // removeTinyEntriesInnerLoop für Amat und adapt rhs. Row reductions selber. col reductions allreduce
      // Bmat und adapt rhs. row reductions selber. col reductions selber
      // Blmat und send adaptions to Zero. row reductions allreduce. col reductions selber
   }
   if(iAmDistrib && rank == 0 )
   {
      // Bmat und adapt rhs. row reductions selber. col reductions allreduce
      // Blmat und adapt rhs. row reductions allreduce. col reductions allreduce
      // receive all adaptations for rhs linking part and adapt them
   }
   if( !iAmDistrib )
   {
      // todo:
      // Bmat und adapt rhs. row and col reductions selber
      // Blmat und adapt rhs. row and col reductions selber
      // for children:
      //    Amat, Bmat, Blmat, adapt rhs, row, col reductions

      SimpleVector* nnzPerRowC = dynamic_cast<SimpleVector*>(nRowElemsC->vec);

      StochGenMatrix& matrixC = dynamic_cast<StochGenMatrix&>(*(presProb->C));
      const StochVector& xlow = dynamic_cast<const StochVector&>(*(presProb->blx));
      const StochVector& xupp = dynamic_cast<const StochVector&>(*(presProb->bux));

      // create empty vectors redRow, redCol, adaptionsRhs with the correct sizes:
      StochVector& redRow(*nRowElemsC);
      redRow.setToZero();
      SimpleVector* reductionsRowSimple = redRow.vec;
      StochVector& redCol(*nColElems);
      redCol.setToZero();
      SimpleVector* reductionsColSimple = redCol.vec;
      StochVector& adaptionsRhs(dynamic_cast<StochVector&>(*(presProb->bu)));
      adaptionsRhs.setToZero();

      StochVector& cupp = dynamic_cast<StochVector&>(*(presProb->bu));
      StochVector& clow = dynamic_cast<StochVector&>(*(presProb->bl));
      StochVector& icupp = dynamic_cast<StochVector&>(*(presProb->icupp));
      StochVector& iclow = dynamic_cast<StochVector&>(*(presProb->iclow));

      SparseStorageDynamic& storageB = matrixC.Bmat->getStorageDynamicRef();
      double* const xlowElemsB = (dynamic_cast<SimpleVector*>(xlow.vec))->elements();
      double* const xuppElemsB = (dynamic_cast<SimpleVector*>(xupp.vec))->elements();
      SimpleVector* adaptionsRhsVec(dynamic_cast<SimpleVector*>(adaptionsRhs.vec));

      // for Bmat:
      removeTinyEntriesInnerLoop( storageB, xlowElemsB, xuppElemsB, nnzPerRowC, reductionsRowSimple, reductionsColSimple, adaptionsRhsVec);
      adaptRhsC(adaptionsRhs, cupp, clow, icupp, iclow);

      SimpleVector* nnzPerRowClink = dynamic_cast<SimpleVector*>(nRowElemsC->vecl);
      SimpleVector* nnzPerColClink = dynamic_cast<SimpleVector*>(nColElems->vecl);
      SparseStorageDynamic& storageBlink = matrixC.Blmat->getStorageDynamicRef();
      double* const xlowElemsBlink = (dynamic_cast<SimpleVector*>(xlow.vecl))->elements();
      double* const xuppElemsBlink = (dynamic_cast<SimpleVector*>(xupp.vecl))->elements();
      SimpleVector* adaptionsRhsVecLink(dynamic_cast<SimpleVector*>(adaptionsRhs.vecl));

      // for Blmat:
      removeTinyEntriesInnerLoop( storageBlink, xlowElemsBlink, xuppElemsBlink, nnzPerRowClink, reductionsRowSimple, reductionsColSimple, adaptionsRhsVecLink);
      adaptRhsC(adaptionsRhs, cupp, clow, icupp, iclow);

      // todo: assert children size
      for( int it = 0; it< matrixC.children.size(); it++){
         //removeTinyEntriesCChild(matrixC.children[it], xlow.children[it], xupp.children[it], redRow.children[it], redCol.children[it], adaptionsRhs.children[it]);
         //adaptRhsC(adaptionsRhs.children[it], cupp.children[it], clow.children[it], icupp.children[it], iclow.children[it]);
         //(*nRowElemsC).children[it]->axpy(1.0, redRow.children[it]);
         //(*nColElems).children[it]->axpy(1.0, redCol.children[it]);
      }

      //(*nRowElemsC).axpy(1.0, reductionsRowSimple);
      //(*nColElems).axpy(1.0, reductionsColSimple);

   }

}

/** Calls removeTinyEntriesInnerLoop for a Child on the matrices Amat, Bmat, Blmat and collects the returned information in redRow, redCol
 * and adaptionsRhs. */
void StochPresolver::removeTinyEntriesCChild(StochGenMatrix* matrix, StochVector* xlow, StochVector* xupp, StochVector* redRow,
                                             StochVector* redCol, StochVector* adaptionsRhs)
{
   SparseStorageDynamic& storageA = dynamic_cast<SparseGenMatrix*>(matrix->Amat)->getStorageDynamicRef();
   SparseStorageDynamic& storageB = dynamic_cast<SparseGenMatrix*>(matrix->Bmat)->getStorageDynamicRef();
   SparseStorageDynamic& storageBl = dynamic_cast<SparseGenMatrix*>(matrix->Blmat)->getStorageDynamicRef();

   // todo: remove entries in Amat, Bmat, Blmat using innerLoop and return combined redRow, redCol and adaptionsRhs


}

/** Removes tiny entries in storage and stores the number of reductions in reductionsRow and reductionsCol respectively. The adaption
 * that have to be made to the rhs (and lhs) are stored in adaptionsRhs. */
void StochPresolver::removeTinyEntriesInnerLoop(SparseStorageDynamic& storage, double* const xlowElems, double* const xuppElems, SimpleVector* nnzPerRow,
                                SimpleVector* reductionsRow, SimpleVector* reductionsCol, SimpleVector* adaptionsRhs)
{
   for( int r = 0; r < storage.m; r++ )
   {
      const int start = storage.rowptr[r].start;
      const int end = storage.rowptr[r].end;
      for( int k = start; k < end; k++ )
      {
         cout << "consider entry ( " << r << ", " << storage.jcolM[k] << " ) = "<< storage.M[k] << endl;
         if( fabs(storage.M[k]) < 1e-3
               && fabs(storage.M[k]) * (xuppElems[k] - xlowElems[k]) * nnzPerRow->elements()[r] < 1e-2 * feastol )
         {
            adaptionsRhs->elements()[r] = -storage.M[k] * xlowElems[k];

            storage.M[k] = 0.0;
            reductionsRow->elements()[r]++;
            reductionsCol->elements()[storage.jcolM[k]]++;
            cout << "entry ( " << r << ", " << storage.jcolM[k] << " ) removed"<< endl;
         }
         if( fabs(storage.M[k]) < 1e-10 )
         {
            storage.M[k] = 0.0;
            reductionsRow->elements()[r]++;
            reductionsCol->elements()[storage.jcolM[k]]++;
            cout << "entry ( " << r << ", " << storage.jcolM[k] << " ) removed"<< endl;
         }
      }
   }
}

/** Adapt the rhs cupp and lhs clow by adding the values in adaptionsRhsStoch to them. */
void StochPresolver::adaptRhsC(StochVector& adaptionsRhsStoch, StochVector& cupp, StochVector& clow, StochVector& icupp, StochVector& iclow )
{
   SimpleVector* cuppS = dynamic_cast<SimpleVector*>(cupp.vec);
   SimpleVector* clowS = dynamic_cast<SimpleVector*>(clow.vec);
   SimpleVector* icuppS = dynamic_cast<SimpleVector*>(icupp.vec);
   SimpleVector* iclowS = dynamic_cast<SimpleVector*>(iclow.vec);
   SimpleVector* adaptionsRhs = dynamic_cast<SimpleVector*>(adaptionsRhsStoch.vec);

   assert(adaptionsRhs->n == cuppS->n);
   assert(adaptionsRhs->n == clowS->n);

   for( int i = 0; i < adaptionsRhs->n; i++ )
   {
      if( adaptionsRhs->elements()[i] != 0.0 )
      {
         if( icuppS->elements()[i] != 0.0 )
            cuppS->elements()[i] += adaptionsRhs->elements()[i];
         if( iclowS->elements()[i] != 0.0 )
            clowS->elements()[i] += adaptionsRhs->elements()[i];
      }
   }
}

/** Adapt the rhs bA by adding the values in adaptionsRhsStoch to it. */
void StochPresolver::adaptRhsA(StochVector* adaptionsRhsStoch, StochVector* b)
{
      SimpleVector* bA = dynamic_cast<SimpleVector*>(b->vec);
      SimpleVector* adaptionsRhs = dynamic_cast<SimpleVector*>(adaptionsRhsStoch->vec);

      assert(adaptionsRhs->n == bA->n);

      for( int i = 0; i < adaptionsRhs->n; i++ )
         if( adaptionsRhs->elements()[i] != 0.0 )
            bA->elements()[i] += adaptionsRhs->elements()[i];
}

