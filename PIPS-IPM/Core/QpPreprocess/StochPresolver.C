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
#include "DoubleMatrixTypes.h"

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


   // clone and initialize dynamic storage
   presProb = sorigprob->cloneFull(true);

   int myRank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

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
   int nElimsTinyEntries = removeTinyEntries();
   if( myRank == 0)
      std::cout << "In total, "<<nElimsTinyEntries<<" tiny entries were removed." << std::endl;


   presProb->cleanUpPresolvedData(*nRowElemsA, *nRowElemsC, *nColElems);


   std::cout << "AFTER \n" << std::endl;

   Apres.deleteTransposed();
   Cpres.deleteTransposed();

   myfile.open("after.txt");
   presProb->writeToStreamDense(myfile);
   myfile.close();

   //presProb->writeToStreamDense(std::cout);

   std::cout << "nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;

   //assert(0);

   return presProb;
}

int StochPresolver::removeTinyEntries()
{
   int nelims = 0;
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( myRank == 0)
      std::cout << "removing tiny entries..." << std::endl;
   // remove tiny entries in A

   // remove tiny entries in C
   nelims += removeTinyEntriesC();

   if( myRank == 0)
      std::cout << "removing tiny entries finished. Removed "<< nelims <<" entries in total." << std::endl;

   return nelims;
}


int StochPresolver::removeTinyEntriesC()
{
   int nelims = 0;

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1)
      iAmDistrib = true;

   double minval = -1.0;
   int index = -1;

   SimpleVector* nnzPerRowC = dynamic_cast<SimpleVector*>(nRowElemsC->vec);

   StochGenMatrix& matrixC = dynamic_cast<StochGenMatrix&>(*(presProb->C));
   const StochVector& xlow = dynamic_cast<const StochVector&>(*(presProb->blx));
   const StochVector& xupp = dynamic_cast<const StochVector&>(*(presProb->bux));

   // create empty vectors redRow, redCol, adaptionsRhs with the correct sizes:
   StochVector* redRow = (*nRowElemsC).clone();
   SimpleVector* reductionsRowSimple = dynamic_cast<SimpleVector*>(redRow->vec);
   StochVector* redCol = (*nColElems).clone();

   StochVector* adaptionsRhs = (dynamic_cast<StochVector&>(*(presProb->bu))).clone();

   StochVector& cupp = dynamic_cast<StochVector&>(*(presProb->bu));
   StochVector& clow = dynamic_cast<StochVector&>(*(presProb->bl));
   StochVector& icupp = dynamic_cast<StochVector&>(*(presProb->icupp));
   StochVector& iclow = dynamic_cast<StochVector&>(*(presProb->iclow));

   SparseStorageDynamic& storageB = matrixC.Bmat->getStorageDynamicRef();
   double* const xlowElemsB = (dynamic_cast<SimpleVector*>(xlow.vec))->elements();
   double* const xuppElemsB = (dynamic_cast<SimpleVector*>(xupp.vec))->elements();
   SimpleVector* adaptionsRhsVec(dynamic_cast<SimpleVector*>(adaptionsRhs->vec));


   SimpleVector* reductionsColVec = dynamic_cast<SimpleVector*>(redCol->vec);

   // for Bmat:
   int nelimsB0 = removeTinyEntriesInnerLoop(storageB, xlowElemsB, xuppElemsB, nnzPerRowC, reductionsRowSimple, reductionsColVec, adaptionsRhsVec);
   adaptRhsC(dynamic_cast<SimpleVector*>(adaptionsRhs->vec), dynamic_cast<SimpleVector*>(cupp.vec), dynamic_cast<SimpleVector*>(clow.vec),
         dynamic_cast<SimpleVector*>(icupp.vec), dynamic_cast<SimpleVector*>(iclow.vec));

   if( myRank == 0 )
      nelims += nelimsB0;

   // update nRowElemsC.vec
   (*nRowElemsC).vec->axpy(-1.0, *redRow->vec);

   (*nRowElemsC).vec->min(minval, index);
   assert( minval >= 0.0 );

   // assert children sizes
   assert( matrixC.children.size() == xlow.children.size() );
   assert( matrixC.children.size() == xupp.children.size() );
   assert( matrixC.children.size() == nRowElemsC->children.size() );
   assert( matrixC.children.size() == redRow->children.size() );
   assert( matrixC.children.size() == redCol->children.size() );
   assert( matrixC.children.size() == adaptionsRhs->children.size() );

   for( size_t it = 0; it< matrixC.children.size(); it++){

      if( matrixC.children[it]->isKindOf(kStochGenDummyMatrix) )
         continue;

      nelims += removeTinyEntriesCChild(it, matrixC, xlow, xupp, nRowElemsC, redRow, redCol, adaptionsRhs);
      adaptRhsC(dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>((*adaptionsRhs).children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(cupp.children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(clow.children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(icupp.children[it])->vec),
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>(iclow.children[it])->vec));

      dynamic_cast<StochVector*>((*nRowElemsC).children[it])->axpy(-1.0, *dynamic_cast<OoqpVector*>(redRow->children[it]));
      dynamic_cast<StochVector*>((*nColElems).children[it])->axpy(-1.0, *dynamic_cast<OoqpVector*>(redCol->children[it]));

      minval = -1.0;
      dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vec->min(minval, index);
      assert( minval >= 0.0 );
      if(dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vecl )
      {
         minval = -1.0;
         dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vecl->min(minval, index);
         assert( minval >= 0.0 );
      }
      minval = -1.0;
      dynamic_cast<StochVector*>((*nColElems).children[it])->vec->min(minval, index);
      assert( minval >= 0.0 );
   }

   // update nColElems.vec via AllReduce
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   (*nColElems).vec->axpy(-1.0, *redCol->vec);

   minval = -1.0;
   index = -1;
   (*nColElems).vec->min(minval, index);
   assert( minval >= 0.0 );

   if( nRowElemsC->vecl )
   {
      // for Blmat:
      SimpleVector* nnzPerRowClink = dynamic_cast<SimpleVector*>(nRowElemsC->vecl);
      SimpleVector* reductionsRowSimpleLink = dynamic_cast<SimpleVector*>(redRow->vecl);
      SparseStorageDynamic& storageBlink = matrixC.Blmat->getStorageDynamicRef();

      SimpleVector* adaptionsRhsVecLink(dynamic_cast<SimpleVector*>(adaptionsRhs->vecl));
      int nelimsLink = removeTinyEntriesInnerLoop(storageBlink, xlowElemsB, xuppElemsB,
            nnzPerRowClink, reductionsRowSimpleLink, reductionsColVec, adaptionsRhsVecLink);

      if( myRank == 0)
         nelims += nelimsLink;

      // update nRowElemsC.vecl and adaptionsRhs->vecl via AllReduce:
      if( iAmDistrib )
      {
         double* redRowLink = dynamic_cast<SimpleVector*>(redRow->vecl)->elements();
         int message_size = dynamic_cast<SimpleVector*>(redRow->vecl)->length();
         MPI_Allreduce(MPI_IN_PLACE, redRowLink, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         double* adaptionsRhsLink = dynamic_cast<SimpleVector*>(adaptionsRhs->vecl)->elements();
         message_size = dynamic_cast<SimpleVector*>(adaptionsRhs->vecl)->length();
         MPI_Allreduce(MPI_IN_PLACE, adaptionsRhsLink, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }

      (*nRowElemsC).vecl->axpy(-1.0, *redRow->vecl);

      minval = -1.0;
      index = -1;
      (*nRowElemsC).vecl->min(minval, index);
      assert( minval >= 0.0 );

      adaptRhsC(dynamic_cast<SimpleVector*>(adaptionsRhs->vecl), dynamic_cast<SimpleVector*>(cupp.vecl), dynamic_cast<SimpleVector*>(clow.vecl),
            dynamic_cast<SimpleVector*>(icupp.vecl), dynamic_cast<SimpleVector*>(iclow.vecl));
   }

   if( iAmDistrib )
   {
      MPI_Allreduce(MPI_IN_PLACE, &nelims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   }

   return nelims;
}

/** Calls removeTinyEntriesInnerLoop for a Child on the matrices Amat, Bmat, Blmat and collects the returned information in redRow, redCol
 * and adaptionsRhs. */
int StochPresolver::removeTinyEntriesCChild(size_t it, StochGenMatrix& matrix, const StochVector& xlow, const StochVector& xupp, StochVector* nnzPerRow,
      StochVector* redRow, StochVector* redCol, StochVector* adaptionsRhs)
{
   int nelims = 0;

   SparseStorageDynamic& storageA = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamicRef();
   SparseStorageDynamic& storageB = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamicRef();

   double* const xlowElemsParent = (dynamic_cast<const SimpleVector*>(xlow.vec))->elements();
   double* const xuppElemsParent = (dynamic_cast<const SimpleVector*>(xupp.vec))->elements();
   double* const xlowElemsChild = (dynamic_cast<const SimpleVector*>(xlow.children[it]->vec))->elements();
   double* const xuppElemsChild = (dynamic_cast<const SimpleVector*>(xupp.children[it]->vec))->elements();
   SimpleVector* nnzRowChild = dynamic_cast<SimpleVector*>(nnzPerRow->children[it]->vec);
   SimpleVector* redRowChild = dynamic_cast<SimpleVector*>(redRow->children[it]->vec);
   SimpleVector* redColParent = dynamic_cast<SimpleVector*>(redCol->vec);
   SimpleVector* redColChild = dynamic_cast<SimpleVector*>(redCol->children[it]->vec);
   SimpleVector* rhsChild = dynamic_cast<SimpleVector*>(adaptionsRhs->children[it]->vec);

   // for Amat:
   nelims += removeTinyEntriesInnerLoop( storageA, xlowElemsParent, xuppElemsParent, nnzRowChild, redRowChild, redColParent, rhsChild);

   // for Bmat:
   nelims += removeTinyEntriesInnerLoop( storageB, xlowElemsChild, xuppElemsChild, nnzRowChild, redRowChild, redColChild, rhsChild);

   if( nnzPerRow->vecl )
   {
      // for Blmat:
      SparseStorageDynamic& storageBl = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamicRef();
      SimpleVector* nnzRowLink = dynamic_cast<SimpleVector*>(nnzPerRow->vecl);
      SimpleVector* redRowLink = dynamic_cast<SimpleVector*>(redRow->vecl);
      SimpleVector* rhsLink = dynamic_cast<SimpleVector*>(adaptionsRhs->vecl);

      nelims += removeTinyEntriesInnerLoop( storageBl, xlowElemsParent, xuppElemsParent, nnzRowLink, redRowLink, redColChild, rhsLink);
   }

   return nelims;
}

/** Removes tiny entries in storage and stores the number of reductions in reductionsRow and reductionsCol respectively. The adaption
 * that have to be made to the rhs (and lhs) are stored in adaptionsRhs. */
int StochPresolver::removeTinyEntriesInnerLoop(SparseStorageDynamic& storage, double* const xlowElems, double* const xuppElems, SimpleVector* nnzPerRow,
                                SimpleVector* reductionsRow, SimpleVector* reductionsCol, SimpleVector* adaptionsRhs)
{
   int nelims = 0;
   double* adaptions = adaptionsRhs->elements();
   double* redRow = reductionsRow->elements();
   double* redCol = reductionsCol->elements();
   double* nnzRow = nnzPerRow->elements();

   for( int r = 0; r < storage.m; r++ )
   {
      const int start = storage.rowptr[r].start;
      const int end = storage.rowptr[r].end;
      for( int k = start; k < end; k++ )
      {
         if( fabs(storage.M[k]) < tolerance1
               && fabs(storage.M[k]) * (xuppElems[k] - xlowElems[k]) * nnzRow[r] < tolerance2 * feastol )
         {
            adaptions[r] -= storage.M[k] * xlowElems[k];

            cout << "Remove entry M ( "<< r << ", " << storage.jcolM[k] << " ) = "<<storage.M[k]<<" (by first test)"<<endl;
            storage.M[k] = 0.0;
            redRow[r]++;
            redCol[storage.jcolM[k]]++;
            nelims ++;
         }
         else if( fabs(storage.M[k]) < tolerance3 )
         {
            cout << "Remove entry M ( "<< r << ", " << storage.jcolM[k] << " ) = "<<storage.M[k]<<" (by second test)"<<endl;
            storage.M[k] = 0.0;
            redRow[r]++;
            redCol[storage.jcolM[k]]++;
            nelims ++;
         }
      }
   }
   return nelims;
}

/** Adapt the rhs cupp and lhs clow by adding the values in adaptionsRhsStoch to them. */
void StochPresolver::adaptRhsC(SimpleVector* adaptionsRhs, SimpleVector* cupp, SimpleVector* clow, SimpleVector* icupp, SimpleVector* iclow )
{
   assert(adaptionsRhs->n == cupp->n);
   assert(adaptionsRhs->n == clow->n);

   double* adaptions = adaptionsRhs->elements();
   double* icuppElems = icupp->elements();
   double* iclowElems = iclow->elements();
   double* cuppElems = cupp->elements();
   double* clowElems = clow->elements();

   for( int i = 0; i < adaptionsRhs->n; i++ )
   {
      if( adaptions[i] != 0.0 )
      {
         if( icuppElems[i] != 0.0 )
            cuppElems[i] += adaptions[i];
         if( iclowElems[i] != 0.0 )
            clowElems[i] += adaptions[i];
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

