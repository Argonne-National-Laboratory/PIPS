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
   redCol = (*nColElems).clone();

   StochVectorHandle bAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   nRowElemsA = bAclone;
   redRowA = (*nRowElemsA).clone();

   StochVectorHandle icuppclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   nRowElemsC = icuppclone;
   redRowC = (*nRowElemsC).clone();

   presProb = NULL;

   setCurrentPointersToNull();

   removedEntries = NULL;
   objOffset = 0.0;
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
   removedEntries = new MTRXENTRY[(int)ceil(0.2 * ((*presProb).my + (*presProb).mz) * (*presProb).nx)];
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
   nelims += removeTinyEntriesSystemA();

   // remove tiny entries in C
   nelims += removeTinyEntriesSystemC();

   if( myRank == 0)
      std::cout << "removing tiny entries finished. Removed "<< nelims <<" entries in total." << std::endl;

   return nelims;
}

void StochPresolver::setCurrentPointersToNull()
{
   currAmat = NULL;
   currBmat = NULL;
   currBlmat = NULL;
   currxlowParent = NULL;
   currxuppParent = NULL;
   currxlowChild = NULL;
   currxuppChild = NULL;
   currEqRhs = NULL;
   currIneqRhs = NULL;
   currIneqLhs = NULL;
   currIcupp = NULL;
   currIclow = NULL;
   currNnzRow = NULL;
   currRedRow = NULL;
   currRedColParent = NULL;
   currRedColChild = NULL;
}

bool StochPresolver::updateCurrentPointers(int it, bool equalitySystem)
{
   StochGenMatrix* matrix;
   currRedColParent = dynamic_cast<SimpleVector*>(redCol->vec);
   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);

   if( equalitySystem )
   {
      matrix = dynamic_cast<StochGenMatrix*>(SpAsPointer(presProb->A));

      if( it == -1 ) // case at root
      {
         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsA->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowA->vec);
      }
      else  // at child it
      {
         if( matrix->children[it]->isKindOf(kStochGenDummyMatrix))
         {
            assert( dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->isKindOf(kStochDummy) );
            assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->isKindOf(kStochDummy) );
            assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->isKindOf(kStochDummy) );
            assert( nRowElemsA->children[it]->isKindOf(kStochDummy) );
            assert( redRowA->children[it]->isKindOf(kStochDummy) );
            assert( redCol->children[it]->isKindOf(kStochDummy) );
            setCurrentPointersToNull();
            return false;
         }
         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->children[it]->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowA->children[it]->vec);
      }
      currIneqRhs = NULL;
      currIneqLhs = NULL;
      currIcupp = NULL;
      currIclow = NULL;
   }
   else  // inequality
   {
      matrix = dynamic_cast<StochGenMatrix*>(SpAsPointer(presProb->C));
      currEqRhs = NULL;

      if( it == -1 ) // case at root
      {
         currIneqRhs =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec);
         currIneqLhs =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vec);
         currIcupp =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vec);
         currIclow =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowC->vec);
      }
      else  // at child it
      {
         if( matrix->children[it]->isKindOf(kStochGenDummyMatrix))
         {
            assert( dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->isKindOf(kStochDummy) );
            assert( dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->isKindOf(kStochDummy) );
            assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->isKindOf(kStochDummy) );
            assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->isKindOf(kStochDummy) );
            assert( dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->isKindOf(kStochDummy) );
            assert( dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->isKindOf(kStochDummy) );
            assert( nRowElemsC->children[it]->isKindOf(kStochDummy) );
            assert( redRowC->children[it]->isKindOf(kStochDummy) );
            assert( redCol->children[it]->isKindOf(kStochDummy) );
            setCurrentPointersToNull();
            return false;
         }
         currIneqRhs =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec);
         currIneqLhs =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec);
         currIcupp =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec);
         currIclow =
               dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->children[it]->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowC->children[it]->vec);

      }
   }
   if( it == -1 ) // case at root
   {
      currAmat = NULL;
      currBmat = matrix->Bmat->getStorageDynamic();
      //currBlmat = matrix->Blmat->getStorageDynamic();
      currRedColChild = NULL;

      currxlowChild = NULL;
      currxuppChild = NULL;
   }

   else  // at child it
   {
      currAmat = dynamic_cast<SparseGenMatrix*>(matrix->children[it]->Amat)->getStorageDynamic();
      currBmat = dynamic_cast<SparseGenMatrix*>(matrix->children[it]->Bmat)->getStorageDynamic();
      //currBlmat = dynamic_cast<SparseGenMatrix*>(matrix->children[it]->Blmat)->getStorageDynamic();
      currRedColChild = dynamic_cast<SimpleVector*>(redCol->children[it]->vec);

      currxlowChild =
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->vec);
      currxuppChild =
            dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->vec);
   }
   return true;
}

int StochPresolver::removeTinyEntriesSystemA()
{
   int nelims = 0;
   double minval = -1.0;
   int index = -1;

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1)
      iAmDistrib = true;

   redRowC->setToZero();
   redCol->setToZero();
   setCurrentPointersToNull();

   StochGenMatrix* matrix = dynamic_cast<StochGenMatrix*>(SpAsPointer(presProb->A));


   int nelimsB0 = 0;
   if( updateCurrentPointers( -1, true) )
   {
      nelimsB0 = removeTinyInnerLoop(*currBmat, currxlowParent, currxuppParent, currRedColParent);

      // update nRowElems.vec
      (*nRowElemsA).vec->axpy(-1.0, *currRedRow);

#ifndef NDEBUG
      minval = -1.0;
      (*nRowElemsA).vec->min(minval, index);
      assert( minval >= 0.0 );
#endif
   }

   if( myRank == 0 )
      nelims += nelimsB0;

   // assert children sizes
   assert( matrix->children.size() == nRowElemsC->children.size() );
   assert( matrix->children.size() == redCol->children.size() );

   // go through the children
   for( size_t it = 0; it< matrix->children.size(); it++)
   {
      if( updateCurrentPointers(int(it), true) )
      {
         nelims += removeTinyChild();

         // update nRowElems.children[it] and nColElems.children[it]
         dynamic_cast<SimpleVector*>(dynamic_cast<StochVector*>((*nRowElemsA).children[it])->vec)->axpy(-1.0, *currRedRow);
         assert( dynamic_cast<StochVector*>((*nRowElemsA).children[it])->vecl == NULL );

         dynamic_cast<StochVector*>((*nColElems).children[it])->vec->axpy(-1.0, *currRedColChild);

#ifndef NDEBUG
         minval = -1.0;
         dynamic_cast<StochVector*>((*nRowElemsA).children[it])->vec->min(minval, index);
         assert( minval >= 0.0 );

         minval = -1.0;
         dynamic_cast<StochVector*>((*nColElems).children[it])->vec->min(minval, index);
         assert( minval >= 0.0 );
#endif
      }
   }

   // update nColElems.vec via AllReduce
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   (*nColElems).vec->axpy(-1.0, *redCol->vec);

#ifndef NDEBUG
   minval = -1.0;
   (*nColElems).vec->min(minval, index);
   assert( minval >= 0.0 );
#endif

   // todo: special treatment for the linking rows

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nelims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   return nelims;
}

int StochPresolver::removeTinyEntriesSystemC()
{
   int nelims = 0;
   double minval = -1.0;
   int index = -1;

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1)
      iAmDistrib = true;

   redRowC->setToZero();
   redCol->setToZero();
   setCurrentPointersToNull();

   StochGenMatrix* matrix = dynamic_cast<StochGenMatrix*>(SpAsPointer(presProb->C));

   int nelimsB0 = 0;
   if( updateCurrentPointers( -1, false) )
   {
      nelimsB0 = removeTinyInnerLoop(*currBmat, currxlowParent, currxuppParent, currRedColParent);

      // update nRowElems.vec
      (*nRowElemsC).vec->axpy(-1.0, *currRedRow);

#ifndef NDEBUG
      minval = -1.0;
      (*nRowElemsC).vec->min(minval, index);
      assert( minval >= 0.0 );
#endif
   }

   if( myRank == 0 )
      nelims += nelimsB0;

   // assert children sizes
   assert( matrix->children.size() == nRowElemsC->children.size() );
   assert( matrix->children.size() == redCol->children.size() );

   // go through the children
   for( size_t it = 0; it< matrix->children.size(); it++)
   {
      if( updateCurrentPointers((int)it, false) )
      {
         nelims += removeTinyChild();

         // update nRowElems.children[it] and nColElems.children[it]
         dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vec->axpy(-1.0, *currRedRow);
         assert( dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vecl == NULL );

         dynamic_cast<StochVector*>((*nColElems).children[it])->vec->axpy(-1.0, *currRedColChild);

#ifndef NDEBUG
         minval = -1.0;
         dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vec->min(minval, index);
         assert( minval >= 0.0 );

         minval = -1.0;
         dynamic_cast<StochVector*>((*nColElems).children[it])->vec->min(minval, index);
         assert( minval >= 0.0 );
#endif
      }
   }

   // update nColElems.vec via AllReduce
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   (*nColElems).vec->axpy(-1.0, *redCol->vec);

#ifndef NDEBUG
   minval = -1.0;
   (*nColElems).vec->min(minval, index);
   assert( minval >= 0.0 );
#endif

   // todo: special treatment for the linking rows

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nelims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   return nelims;
}

/** Calls removeTinyInnerLoop for a Child on the matrices Amat, Bmat which also adapts the rhs. */
int StochPresolver::removeTinyChild()
{
   int nelims = 0;

   // for Amat:
   nelims += removeTinyInnerLoop( *currAmat, currxlowParent, currxuppParent, currRedColParent);

   // for Bmat:
   nelims += removeTinyInnerLoop( *currBmat, currxlowChild, currxuppChild, currRedColChild);

   // todo: special treatment for the linking rows

   return nelims;
}

/** Removes tiny entries in storage and adapts the rhs accordingly. */
int StochPresolver::removeTinyInnerLoop(SparseStorageDynamic& storage, SimpleVector* const xlow, SimpleVector* const xupp, SimpleVector* reductionsCol)
{
   int nelims = 0;
   double* redRow = currRedRow->elements();
   double* redCol = reductionsCol->elements();
   double* nnzRow = currNnzRow->elements();
   double* xuppElems = xupp->elements();
   double* xlowElems = xlow->elements();

   bool equality = false;
   if( currEqRhs != NULL)
      equality = true;

   double* rhsElems;
   double* icuppElems;
   double* iclowElems;
   double* cuppElems;
   double* clowElems;

   if( equality )
      rhsElems = currEqRhs->elements();
   else
   {
      icuppElems = currIcupp->elements();
      iclowElems = currIclow->elements();
      cuppElems = currIneqRhs->elements();
      clowElems = currIneqLhs->elements();
   }

   for( int r = 0; r < storage.m; r++ )
   {
      const int start = storage.rowptr[r].start;
      int end = storage.rowptr[r].end;
      int k = start;
      while( k < end )
      {
         const int col = storage.jcolM[k];
         if( fabs(storage.M[k]) < tolerance3 )
         {
            cout << "Remove entry M ( "<< r << ", " << storage.jcolM[k] << " ) = "<<storage.M[k]<<" (by first test)"<<endl;

            redRow[r]++;
            redCol[storage.jcolM[k]]++;
            nelims ++;

            //removedEntries[nelimsGlobal].rowIdx = r;
            //removedEntries[nelimsGlobal].colIdx = storage.jcolM[k];
            //storage.M[k] = 0.0;
            std::swap(storage.M[k],storage.M[end-1]);
            std::swap(storage.jcolM[k],storage.jcolM[end-1]);
            storage.rowptr[r].end --;
            end = storage.rowptr[r].end;
            k--;
         }
         else if( fabs(storage.M[k]) < tolerance1
               && fabs(storage.M[k]) * (xuppElems[col] - xlowElems[col]) * nnzRow[r] < tolerance2 * feastol )
         {
            if( equality )
               rhsElems[r] -= storage.M[k] * xlowElems[col];
            else
            {
               if( icuppElems[r] != 0.0 )
                  cuppElems[r] -= storage.M[k] * xlowElems[col];
               if( iclowElems[r] != 0.0 )
                  clowElems[r] -= storage.M[k] * xlowElems[col];
            }

            cout << "Remove entry M ( "<< r << ", " << col << " ) = "<<storage.M[k]<<" (by second test)"<<endl;
            redRow[r]++;
            redCol[storage.jcolM[k]]++;
            nelims ++;

            //removedEntries[nelimsGlobal]rowIdx = r;
            //removedEntries[nelimsGlobal].colIdx = storage.jcolM[k];
            // storage.M[k] = 0.0;
            std::swap(storage.M[k],storage.M[end-1]);
            std::swap(storage.jcolM[k],storage.jcolM[end-1]);
            storage.rowptr[r].end --;
            end = storage.rowptr[r].end;
            k--;
         }
         k++;
      }
   }
   return nelims;
}


