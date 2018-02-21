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
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmSpecial = true;
   if( world_size > 1 && myRank > 0)
      iAmSpecial = false;

   SimpleVector* nnzPerRowC = dynamic_cast<SimpleVector*>(nRowElemsC->vec);

   StochGenMatrix& matrixC = dynamic_cast<StochGenMatrix&>(*(presProb->C));
   const StochVector& xlow = dynamic_cast<const StochVector&>(*(presProb->blx));
   const StochVector& xupp = dynamic_cast<const StochVector&>(*(presProb->bux));

   // create empty vectors redRow, redCol, adaptionsRhs with the correct sizes:
   StochVector* redRow = (*nRowElemsC).clone();
   redRow->setToZero();
   SimpleVector* reductionsRowSimple = dynamic_cast<SimpleVector*>(redRow->vec);
   StochVector* redCol = (*nColElems).clone();
   redCol->setToZero();
   SimpleVector* reductionsColVec = dynamic_cast<SimpleVector*>(redCol->vec);
   StochVector* adaptionsRhs = (dynamic_cast<StochVector&>(*(presProb->bu))).clone();
   adaptionsRhs->setToZero();

   StochVector& cupp = dynamic_cast<StochVector&>(*(presProb->bu));
   StochVector& clow = dynamic_cast<StochVector&>(*(presProb->bl));
   StochVector& icupp = dynamic_cast<StochVector&>(*(presProb->icupp));
   StochVector& iclow = dynamic_cast<StochVector&>(*(presProb->iclow));

   SparseStorageDynamic& storageB = matrixC.Bmat->getStorageDynamicRef();
   double* const xlowElemsB = (dynamic_cast<SimpleVector*>(xlow.vec))->elements();
   double* const xuppElemsB = (dynamic_cast<SimpleVector*>(xupp.vec))->elements();
   SimpleVector* adaptionsRhsVec(dynamic_cast<SimpleVector*>(adaptionsRhs->vec));

   SimpleVector* nnzPerRowClink = dynamic_cast<SimpleVector*>(nRowElemsC->vecl);
   SimpleVector* reductionsRowSimpleLink = dynamic_cast<SimpleVector*>(redRow->vecl);
   SparseStorageDynamic& storageBlink = matrixC.Blmat->getStorageDynamicRef();
   double* const xlowElemsBlink = (dynamic_cast<SimpleVector*>(xlow.vecl))->elements();
   double* const xuppElemsBlink = (dynamic_cast<SimpleVector*>(xupp.vecl))->elements();
   SimpleVector* adaptionsRhsVecLink(dynamic_cast<SimpleVector*>(adaptionsRhs->vecl));

   if( iAmSpecial )
   {
      // for Bmat:
      removeTinyEntriesInnerLoop(storageB, xlowElemsB, xuppElemsB, nnzPerRowC, reductionsRowSimple, reductionsColVec, adaptionsRhsVec);
      adaptRhsC(dynamic_cast<SimpleVector*>(adaptionsRhs->vec), dynamic_cast<SimpleVector*>(cupp.vec), dynamic_cast<SimpleVector*>(clow.vec),
            dynamic_cast<SimpleVector*>(icupp.vec), dynamic_cast<SimpleVector*>(iclow.vec));
      (*nRowElemsC).vec->axpy(1.0, *redRow->vec);

      // for Blmat:
      removeTinyEntriesInnerLoop(storageBlink, xlowElemsBlink, xuppElemsBlink,
            nnzPerRowClink, reductionsRowSimpleLink, reductionsColVec, adaptionsRhsVecLink);
   }

   // todo: assert children size

   for( size_t it = 0; it< matrixC.children.size(); it++){

      removeTinyEntriesCChild(it, matrixC, xlow, xupp, nRowElemsC, redRow, redCol, adaptionsRhs);
      adaptRhsC(dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>((*adaptionsRhs).children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(cupp.children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(clow.children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(icupp.children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(iclow.children[it])->vec));
      (*nRowElemsC).children[it]->axpy(1.0, *dynamic_cast<OoqpVector*>(redRow->children[it]));
      (*nColElems).children[it]->axpy(1.0, *dynamic_cast<OoqpVector*>(redCol->children[it]));
   }

   // update nColElems.vec via AllReduce
   if( world_size > 1 )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      double* redColParentG = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(redColParent, redColParentG, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      redColParent = redColParentG;
   }
   (*nColElems).vec->axpy(1.0, *redCol->vec);

   // update nRowElemsC.vecl via AllReduce
   if( world_size > 1 )
   {
      double* redRowLink = dynamic_cast<SimpleVector*>(redRow->vecl)->elements();
      double* redRowLinkG = dynamic_cast<SimpleVector*>(redRow->vecl)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redRow->vecl)->length();
      MPI_Allreduce(redRowLink, redRowLinkG, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      redRowLink = redRowLinkG; //todo: check if only redRowLink is changed or original redRow ?
   }
   (*nRowElemsC).vecl->axpy(1.0, *redRow->vecl);

   // todo: update adaptionsRhs.vecl via AllReduce
   if( world_size > 1 )
   {
      double* adaptionsRhsLink = dynamic_cast<SimpleVector*>(adaptionsRhs->vecl)->elements();
      double* adaptionsRhsLinkG = dynamic_cast<SimpleVector*>(adaptionsRhs->vecl)->elements();
      int message_size = dynamic_cast<SimpleVector*>(adaptionsRhs->vecl)->length();
      MPI_Allreduce(adaptionsRhsLink, adaptionsRhsLinkG, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      adaptionsRhsLink = adaptionsRhsLinkG;
   }
   (*nRowElemsC).vecl->axpy(1.0, *redRow->vecl);
   adaptRhsC(dynamic_cast<SimpleVector*>(adaptionsRhs->vecl), dynamic_cast<SimpleVector*>(cupp.vecl), dynamic_cast<SimpleVector*>(clow.vecl),
         dynamic_cast<SimpleVector*>(icupp.vecl), dynamic_cast<SimpleVector*>(iclow.vecl));


}

/** Calls removeTinyEntriesInnerLoop for a Child on the matrices Amat, Bmat, Blmat and collects the returned information in redRow, redCol
 * and adaptionsRhs. */
void StochPresolver::removeTinyEntriesCChild(size_t it, StochGenMatrix& matrix, const StochVector& xlow, const StochVector& xupp, StochVector* nnzPerRow,
      StochVector* redRow, StochVector* redCol, StochVector* adaptionsRhs)
{
   SparseStorageDynamic& storageA = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamicRef();
   SparseStorageDynamic& storageB = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamicRef();
   SparseStorageDynamic& storageBl = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamicRef();

   // todo: remove entries in Amat, Bmat, Blmat using innerLoop and return combined redRow, redCol and adaptionsRhs

   double* const xlowElemsParent = (dynamic_cast<const SimpleVector*>(xlow.vec))->elements();
   double* const xuppElemsParent = (dynamic_cast<const SimpleVector*>(xupp.vec))->elements();
   double* const xlowElemsChild = (dynamic_cast<const SimpleVector*>(xlow.children[it]->vec))->elements();
   double* const xuppElemsChild = (dynamic_cast<const SimpleVector*>(xupp.children[it]->vec))->elements();
   SimpleVector* nnzRowChild = (dynamic_cast<SimpleVector*>(nnzPerRow->children[it]->vec));
   SimpleVector* nnzRowLink = (dynamic_cast<SimpleVector*>(nnzPerRow->vecl));
   SimpleVector* redRowChild = (dynamic_cast<SimpleVector*>(redRow->children[it]->vec));
   SimpleVector* redRowLink = (dynamic_cast<SimpleVector*>(redRow->vecl));
   SimpleVector* redColParent = (dynamic_cast<SimpleVector*>(redCol->vec));
   SimpleVector* redColChild = (dynamic_cast<SimpleVector*>(redCol->children[it]->vec));
   SimpleVector* rhsLink = (dynamic_cast<SimpleVector*>(adaptionsRhs->vecl));
   SimpleVector* rhsChild = (dynamic_cast<SimpleVector*>(adaptionsRhs->children[it]->vec));

   // for Amat:
   removeTinyEntriesInnerLoop( storageA, xlowElemsParent, xuppElemsParent, nnzRowChild, redRowChild, redColParent, rhsChild);

   // for Bmat:
   removeTinyEntriesInnerLoop( storageB, xlowElemsChild, xuppElemsChild, nnzRowChild, redRowChild, redColChild, rhsChild);

   // for Blmat:
   removeTinyEntriesInnerLoop( storageBl, xlowElemsParent, xuppElemsParent, nnzRowLink, redRowLink, redColChild, rhsLink);


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
void StochPresolver::adaptRhsC(SimpleVector* adaptionsRhs, SimpleVector* cupp, SimpleVector* clow, SimpleVector* icupp, SimpleVector* iclow )
{
  /* SimpleVector* cuppS = dynamic_cast<SimpleVector*>(cupp.vec);
   SimpleVector* clowS = dynamic_cast<SimpleVector*>(clow.vec);
   SimpleVector* icuppS = dynamic_cast<SimpleVector*>(icupp.vec);
   SimpleVector* iclowS = dynamic_cast<SimpleVector*>(iclow.vec);
   SimpleVector* adaptionsRhs = dynamic_cast<SimpleVector*>(adaptionsRhsStoch.vec);*/

   assert(adaptionsRhs->n == cupp->n);
   assert(adaptionsRhs->n == clow->n);

   for( int i = 0; i < adaptionsRhs->n; i++ )
   {
      if( adaptionsRhs->elements()[i] != 0.0 )
      {
         if( icupp->elements()[i] != 0.0 )
            cupp->elements()[i] += adaptionsRhs->elements()[i];
         if( iclow->elements()[i] != 0.0 )
            clow->elements()[i] += adaptionsRhs->elements()[i];
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

      // todo: children
}

