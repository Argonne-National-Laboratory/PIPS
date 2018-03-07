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
   StochVectorHandle colClone(dynamic_cast<StochVector*>(sorigprob->g->clone()));
   redCol = colClone;

   StochVectorHandle bAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   nRowElemsA = bAclone;
   StochVectorHandle rowAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   redRowA = rowAclone;

   StochVectorHandle icuppclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   nRowElemsC = icuppclone;
   StochVectorHandle rowCclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   redRowC = rowCclone;

   presProb = NULL;

   setCurrentPointersToNull();

   localNelims = 0;
   nChildren = nColElems->children.size();
   linkVarsBlocks = new BLOCKS[nChildren + 1];
   childBlocks = new BLOCKS[nChildren + 1];
   resetLinkvarsAndChildBlocks();

   blocks = new int[nChildren + 3];
   resetBlocks();
   colBlocks = new int[nChildren + 2];
//todo: reset colBlocks

   objOffset = 0.0;
}


StochPresolver::~StochPresolver()
{
   delete[] linkVarsBlocks;
   delete[] childBlocks;
   delete[] blocks;
   delete[] colBlocks;
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

   removedEntries.reserve((int)ceil(0.2 * ((*presProb).my + (*presProb).mz) * (*presProb).nx));
   int nElimsTinyEntries = removeTinyEntries();
   if( myRank == 0)
      std::cout << "In total, "<<nElimsTinyEntries<<" tiny entries were removed." << std::endl;

   int nElimsSingletonRow = doSingletonRows();
   if( myRank == 0)
         std::cout << "In total, "<<nElimsSingletonRow<<" singleton rows were considered." << std::endl;

   presProb->cleanUpPresolvedData(*nRowElemsA, *nRowElemsC, *nColElems);

   if( myRank == 0)
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

   nelims += removeTinyEntriesSystemA();
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*presProb->A);
   updateTransposed(A);

   if( myRank == 0)
         std::cout << "removing tiny entries in inequality system C..." << std::endl;

   nelims += removeTinyEntriesSystemC();
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*presProb->C);
   updateTransposed(C);

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

bool StochPresolver::updateCurrentPointers(int it, SystemType system_type)
{
   currRedColParent = dynamic_cast<SimpleVector*>(redCol->vec);
   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);

   if( system_type == EQUALITY_SYSTEM )
   {
      StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

      if( it == -1 ) // case at root
      {
         currAmat = matrix.Bmat->getStorageDynamic();   // save Bmat as currAmat for easy computation in TinyInnerLoop

         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsA->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowA->vec);
      }
      else  // at child it
      {
         if( matrix.children[it]->isKindOf(kStochGenDummyMatrix))
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
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
         //currBlmat = dynamic_cast<SparseGenMatrix*>(matrix->children[it]->Blmat)->getStorageDynamic();

         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->children[it]->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowA->children[it]->vec);
      }
      currIneqRhs = NULL;
      currIneqLhs = NULL;
      currIcupp = NULL;
      currIclow = NULL;
   }
   else  // system_type == INEQUALITY_SYSTEM
   {
      assert( system_type == INEQUALITY_SYSTEM );

      StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));
      currEqRhs = NULL;

      if( it == -1 ) // case at root
      {
         currAmat = matrix.Bmat->getStorageDynamic();   // save Bmat as currAmat for easy computation in TinyInnerLoop

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
         if( matrix.children[it]->isKindOf(kStochGenDummyMatrix))
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
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
         //currBlmat = dynamic_cast<SparseGenMatrix*>(matrix->children[it]->Blmat)->getStorageDynamic();

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
      currBmat = NULL;  // Bmat is saved as currAmat for easy computation in TinyInnerLoop
      //currBlmat = matrix->Blmat->getStorageDynamic();
      currRedColChild = NULL;

      currxlowChild = NULL;
      currxuppChild = NULL;
   }

   else  // at child it
   {
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
   localNelims = 0;

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

   int nelimsB0 = 0;
   if( updateCurrentPointers( -1, EQUALITY_SYSTEM) )
   {
      nelimsB0 = removeTinyInnerLoop( -1, EQUALITY_SYSTEM, LINKING_VARS_BLOCK );
      assert( nelimsB0 == localNelims );
      updateNnzUsingReductions( (*nRowElemsA).vec, currRedRow);
   }

   if( myRank == 0 )
      nelims += nelimsB0;

   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

   assert( matrix.children.size() == nRowElemsC->children.size() );
   assert( matrix.children.size() == redCol->children.size() );

   // go through the children
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( updateCurrentPointers(int(it), EQUALITY_SYSTEM) )
      {
         nelims += removeTinyChild((int)it, EQUALITY_SYSTEM);

         updateNnzUsingReductions( dynamic_cast<StochVector*>((*nRowElemsA).children[it])->vec, currRedRow);
         updateNnzUsingReductions( dynamic_cast<StochVector*>((*nColElems).children[it])->vec, currRedColChild);

         assert( dynamic_cast<StochVector*>((*nRowElemsA).children[it])->vecl == NULL );
      }
   }

   // update nColElems.vec via AllReduce
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions((*nColElems).vec, redCol->vec);

   // todo: special treatment for the linking rows

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nelims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   return nelims;
}

int StochPresolver::removeTinyEntriesSystemC()
{
   int nelims = 0;
   localNelims = 0;

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

   int nelimsB0 = 0;
   if( updateCurrentPointers( -1, INEQUALITY_SYSTEM) )
   {
      nelimsB0 = removeTinyInnerLoop( -1, INEQUALITY_SYSTEM, LINKING_VARS_BLOCK );
      assert( nelimsB0 == localNelims );
      updateNnzUsingReductions( (*nRowElemsC).vec, currRedRow);
   }

   if( myRank == 0 )
      nelims += nelimsB0;

   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   assert( matrix.children.size() == nRowElemsC->children.size() );
   assert( matrix.children.size() == redCol->children.size() );

   // go through the children
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( updateCurrentPointers((int)it, INEQUALITY_SYSTEM) )
      {
         nelims += removeTinyChild((int)it, INEQUALITY_SYSTEM);

         updateNnzUsingReductions( dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vec, currRedRow);
         updateNnzUsingReductions( dynamic_cast<StochVector*>((*nColElems).children[it])->vec, currRedColChild);

         assert( dynamic_cast<StochVector*>((*nRowElemsC).children[it])->vecl == NULL );
      }
   }

   // update nColElems.vec via AllReduce
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions((*nColElems).vec, redCol->vec);

   // todo: special treatment for the linking rows

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nelims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   return nelims;
}

/** Calls removeTinyInnerLoop for a Child on the matrices Amat, Bmat which also adapts the rhs. */
int StochPresolver::removeTinyChild( int it, SystemType system_type )
{
   int nelims = 0;

   // for Amat:
   nelims += removeTinyInnerLoop( it, system_type, LINKING_VARS_BLOCK );

   // for Bmat:
   nelims += removeTinyInnerLoop( it, system_type, CHILD_BLOCK );

   // todo: special treatment for the linking rows

   return nelims;
}

/** Removes tiny entries in storage and adapts the rhs accordingly.
 *  If block_type == LINKING_VARS_BLOCK, then block Amat is considered.
 *  If block_type == CHILD_BLOCK, then block Bmat is considered. */
int StochPresolver::removeTinyInnerLoop( int it, SystemType system_type, BlockType block_type )
{
   // Setting all the pointers correctly
   int nelims = 0;
   double* nnzRow = currNnzRow->elements();
   double* xlowElems;
   double* xuppElems;
   double* redCol;
   SparseStorageDynamic* storage;

   if( block_type == LINKING_VARS_BLOCK )
   {
      storage = currAmat;
      xlowElems = currxlowParent->elements();
      xuppElems = currxuppParent->elements();
      redCol = currRedColParent->elements();

      linkVarsBlocks[it+1].start = localNelims;
   }
   else if( block_type == CHILD_BLOCK )
   {
      storage = currBmat;
      xlowElems = currxlowChild->elements();
      xuppElems = currxuppChild->elements();
      redCol = currRedColChild->elements();

      childBlocks[it+1].start = localNelims;
   }

   double* rhsElems;
   double* icuppElems;
   double* iclowElems;
   double* cuppElems;
   double* clowElems;

   if( system_type == EQUALITY_SYSTEM )
      rhsElems = currEqRhs->elements();
   else
   {
      assert( system_type == INEQUALITY_SYSTEM );
      icuppElems = currIcupp->elements();
      iclowElems = currIclow->elements();
      cuppElems = currIneqRhs->elements();
      clowElems = currIneqLhs->elements();
   }

   // the actual work starts here:
   for( int r = 0; r < storage->m; r++ )
   {
      const int start = storage->rowptr[r].start;
      int end = storage->rowptr[r].end;
      int k = start;
      while( k < end )
      {
         const int col = storage->jcolM[k];
         if( fabs(storage->M[k]) < tolerance3 )
         {
            cout << "Remove entry M ( "<< r << ", " << storage->jcolM[k] << " ) = "<<storage->M[k]<<" (by first test)"<<endl;
            storeRemovedEntryIndex(r, storage->jcolM[k], it, block_type);
            updateAndSwap(storage, r, k, end, redCol, nelims);
         }
         else if( fabs(storage->M[k]) < tolerance1
               && fabs(storage->M[k]) * (xuppElems[col] - xlowElems[col]) * nnzRow[r] < tolerance2 * feastol )
         {
            if( system_type == EQUALITY_SYSTEM )
               rhsElems[r] -= storage->M[k] * xlowElems[col];
            else
            {
               if( icuppElems[r] != 0.0 )
                  cuppElems[r] -= storage->M[k] * xlowElems[col];
               if( iclowElems[r] != 0.0 )
                  clowElems[r] -= storage->M[k] * xlowElems[col];
            }

            cout << "Remove entry M ( "<< r << ", " << col << " ) = "<<storage->M[k]<<" (by second test)"<<endl;
            storeRemovedEntryIndex(r, storage->jcolM[k], it, block_type);
            updateAndSwap(storage, r, k, end, redCol, nelims);
         }
         k++;
      }
   }
   if( block_type == LINKING_VARS_BLOCK )
      linkVarsBlocks[it+1].end = localNelims;
   else if( block_type == CHILD_BLOCK )
      childBlocks[it+1].end = localNelims;

   return nelims;
}

void StochPresolver::updateAndSwap( SparseStorageDynamic* storage, int rowidx, int& indexK, int& rowEnd, double* redCol, int& nelims)
{
   double* redRow = currRedRow->elements();

   redRow[rowidx]++;
   redCol[storage->jcolM[indexK]]++;

   std::swap(storage->M[indexK],storage->M[rowEnd-1]);
   std::swap(storage->jcolM[indexK],storage->jcolM[rowEnd-1]);
   storage->rowptr[rowidx].end --;
   rowEnd = storage->rowptr[rowidx].end;
   indexK--;

   nelims++;
}

void StochPresolver::storeRemovedEntryIndex(int rowidx, int colidx, int it, BlockType block_type)
{
   assert( (int)removedEntries.size() == localNelims );

   if( block_type == LINKING_VARS_BLOCK )
      assert( linkVarsBlocks[it+1].start <= localNelims );

   else if( block_type == CHILD_BLOCK )
      assert( childBlocks[it+1].start <= localNelims );

   MTRXENTRY entry = {rowidx, colidx};
   removedEntries.push_back(entry);

   localNelims++;
}

/** Update the nnzVector by subtracting the reductions vector. */
void StochPresolver::updateNnzUsingReductions( OoqpVector* nnzVector, OoqpVector* redVector)
{
   SimpleVector* redSimple = dynamic_cast<SimpleVector*>(redVector);
   nnzVector->axpy(-1.0, *redSimple);

#ifndef NDEBUG
   double minval = -1.0;
   int index = -1;
   nnzVector->min(minval, index);
   assert( minval >= 0.0 );
#endif
}

void StochPresolver::updateTransposed(StochGenMatrix& matrix)
{
   // update matrix' using removedEntries, linkVarsBlocks, childBlocks

   if( linkVarsBlocks[0].start != linkVarsBlocks[0].end )
   {
      assert(linkVarsBlocks[0].start < linkVarsBlocks[0].end);
      SparseStorageDynamic& B0_trans = matrix.Bmat->getStorageDynamicTransposedRef();
      updateTransposedSubmatrix(B0_trans, linkVarsBlocks[0].start, linkVarsBlocks[0].end);
   }

   for( int i = 1; i < nChildren+1; i++)
   {
      if( linkVarsBlocks[i].start != linkVarsBlocks[i].end )
      {
         assert(linkVarsBlocks[i].start < linkVarsBlocks[i].end);
         SparseStorageDynamic& AChild_trans = matrix.children[i-1]->Amat->getStorageDynamicTransposedRef();
         updateTransposedSubmatrix(AChild_trans, linkVarsBlocks[i].start, linkVarsBlocks[i].end);
      }

      if( childBlocks[i].start != childBlocks[i].end )
      {
         assert(childBlocks[i].start < childBlocks[i].end);
         SparseStorageDynamic& BChild_trans = matrix.children[i-1]->Bmat->getStorageDynamicTransposedRef();
         updateTransposedSubmatrix(BChild_trans, childBlocks[i].start, childBlocks[i].end);
      }
   }

   // set localNelims, removedEntries, linkVarsBlocks, childBlocks to zero again.
   localNelims = 0;
   removedEntries.clear();
   resetLinkvarsAndChildBlocks();
}

int StochPresolver::doSingletonRows()
{
   int nelims = 0;
   int nSingleRows = initSingletonRows(EQUALITY_SYSTEM);
   cout<<"Found "<<nSingleRows<<" singleton rows in equality system A."<<endl;
   nelims += doSingletonRowsA();
   nelims += doSingletonRowsC();
   return nelims;
}

int StochPresolver::initSingletonRows(SystemType system_type)
{
   int nSingletonRows = 0;
   if( system_type == EQUALITY_SYSTEM )
   {
      assert(singletonRows.size() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(nRowElemsA->vec);
      nSingletonRows += initSingletonRowsBlock(-1, nRowASimple);

      assert((int)nRowElemsA->children.size() == nChildren);
      for( size_t it = 0; it < nRowElemsA->children.size(); it++)
      {
         SimpleVector* nRowASimpleChild = dynamic_cast<SimpleVector*>(nRowElemsA->children[it]->vec);
         nSingletonRows += initSingletonRowsBlock(int(it), nRowASimpleChild);
      }
      blocks[nChildren+1] = singletonRows.size();

      // todo: linking block nRowElemsA->vecl
      //blocks[nChildren+2] = singletonRows.size();
   }
   else
   {
      assert( system_type == INEQUALITY_SYSTEM );
      //todo: if rows from A and C should be stored, then another variable to store the indices is needed,
      // blocks is not enough.
      assert(singletonRows.size() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(nRowElemsC->vec);
      nSingletonRows += initSingletonRowsBlock(-1, nRowASimple);

      assert((int)nRowElemsC->children.size() == nChildren);
      for( size_t it = 0; it < nRowElemsC->children.size(); it++)
      {
         SimpleVector* nRowASimpleChild = dynamic_cast<SimpleVector*>(nRowElemsC->children[it]->vec);
         nSingletonRows += initSingletonRowsBlock(int(it), nRowASimpleChild);
      }
      blocks[nChildren+1] = singletonRows.size();

      // todo: linking block nRowElemsC->vecl
      //blocks[nChildren+2] = singletonRows.size();
   }
   return nSingletonRows;
}

int StochPresolver::initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple)
{
   int nSingletonRows = 0;

   blocks[it+1] = singletonRows.size();
   double* nnzRow = nnzRowSimple->elements();

   for( int i = 0; i < nnzRowSimple->n; i++)
      if( nnzRow[i] == 1.0 )
      {
         singletonRows.push_back(i);
         nSingletonRows++;
      }
   return nSingletonRows;
}

int StochPresolver::doSingletonRowsA()
{
   int nelims = 0;
   return nelims;
}

int StochPresolver::doSingletonRowsC()
{
   int nelims = 0;
   return nelims;
}

void StochPresolver::resetLinkvarsAndChildBlocks()
{
   for( int i = 0; i < nChildren+1; i++)
   {
      linkVarsBlocks[i].start = 0;
      linkVarsBlocks[i].end = 0;
      childBlocks[i].start = 0;
      childBlocks[i].end = 0;
   }
}

void StochPresolver::resetBlocks()
{
   for( int i = 0; i < nChildren+3; i++)
      blocks[i] = 0;
}

void StochPresolver::updateTransposedSubmatrix(SparseStorageDynamic& transStorage, const int blockStart, const int blockEnd)
{
   assert( blockEnd <= (int)removedEntries.size());

   for( int j = blockStart; j < blockEnd; j++)
   {
      int row_A = removedEntries[j].rowIdx;
      int row_At = removedEntries[j].colIdx;

      int start = transStorage.rowptr[row_At].start;
      int end = transStorage.rowptr[row_At].end;
      int col_At;

      for( col_At = start; col_At < end; col_At++)
      {
         if( transStorage.jcolM[col_At] == row_A )
            break;
      }

      std::swap(transStorage.M[col_At],transStorage.M[end-1]);
      std::swap(transStorage.jcolM[col_At],transStorage.jcolM[end-1]);
      transStorage.rowptr[row_At].end --;
   }
}


