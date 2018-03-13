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
   blocksIneq = new int[nChildren + 3];
   colBlocksChildren = new int[nChildren + 1];
   resetBlocks();
   resetColBlocks();

   objOffset = 0.0;
}


StochPresolver::~StochPresolver()
{
   delete[] linkVarsBlocks;
   delete[] childBlocks;
   delete[] blocks;
   delete[] blocksIneq;
   delete[] colBlocksChildren;
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
   //int nElimsTinyEntries = removeTinyEntries();
   //if( myRank == 0)
   //   std::cout << "In total, "<<nElimsTinyEntries<<" tiny entries were removed." << std::endl;

   doSingletonRows();

   myfile.open("middle.txt");
   presProb->writeToStreamDense(myfile);
   myfile.close();

/*   cout<<"nRowElemsA "<<endl;
   nRowElemsA->writeToStreamAll(cout);
   cout<<"nRowElemsC "<<endl;
   nRowElemsC->writeToStreamAll(cout);
   cout<<"nColElems "<<endl;
   nColElems->writeToStreamAll(cout);
*/
   presProb->cleanUpPresolvedData(*nRowElemsA, *nRowElemsC, *nColElems);

   if( myRank == 0)
      std::cout << "AFTER \n" << std::endl;

   Apres.deleteTransposed();
   Cpres.deleteTransposed();

   myfile.open("after.txt");
   presProb->writeToStreamDense(myfile);
   myfile.close();

   //presProb->writeToStreamDense(std::cout);

   std::cout << "sorigprob nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
   std::cout << "presProb nx, my, mz" << presProb->nx << " " << presProb->my << " " << presProb->mz << std::endl;

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
   currAmatTrans = NULL;
   currBmat = NULL;
   currBmatTrans = NULL;
   currBlmat = NULL;
   currxlowParent = NULL;
   currxuppParent = NULL;
   currxlowChild = NULL;
   currxuppChild = NULL;
   currIxlowParent = NULL;
   currIxlowChild = NULL;
   currIxuppParent = NULL;
   currIxuppChild = NULL;
   currEqRhs = NULL;
   currIneqRhs = NULL;
   currIneqLhs = NULL;
   currIcupp = NULL;
   currIclow = NULL;
   currNnzRow = NULL;
   currRedRow = NULL;
   currRedColParent = NULL;
   currRedColChild = NULL;
   currgParent = NULL;
   currgChild = NULL;
   currNnzColChild = NULL;
}

bool StochPresolver::updateCurrentPointers(int it, SystemType system_type)
{
   currRedColParent = dynamic_cast<SimpleVector*>(redCol->vec);
   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);
   currIxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).vec);
   currIxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).vec);

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
         if( childIsDummy(matrix, it, EQUALITY_SYSTEM) ) return false;

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

         currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec);
         currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vec);
         currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vec);
         currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowC->vec);
      }
      else  // at child it
      {
         if( childIsDummy(matrix, it, INEQUALITY_SYSTEM) ) return false;

         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
         //currBlmat = dynamic_cast<SparseGenMatrix*>(matrix->children[it]->Blmat)->getStorageDynamic();

         currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec);
         currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec);
         currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec);
         currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec);
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
      currIxlowChild = NULL;
      currIxuppChild = NULL;
   }
   else  // at child it
   {
      currRedColChild = dynamic_cast<SimpleVector*>(redCol->children[it]->vec);
      currxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->vec);
      currxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->vec);
      currIxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[it]->vec);
      currIxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[it]->vec);
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
   double* ixlowElems;
   double* ixuppElems;
   double* redCol;
   SparseStorageDynamic* storage;

   if( block_type == LINKING_VARS_BLOCK )
   {
      storage = currAmat;
      xlowElems = currxlowParent->elements();
      xuppElems = currxuppParent->elements();
      ixlowElems = currIxlowParent->elements();
      ixuppElems = currIxuppParent->elements();

      redCol = currRedColParent->elements();

      linkVarsBlocks[it+1].start = localNelims;
   }
   else if( block_type == CHILD_BLOCK )
   {
      storage = currBmat;
      xlowElems = currxlowChild->elements();
      xuppElems = currxuppChild->elements();
      ixlowElems = currIxlowChild->elements();
      ixuppElems = currIxuppChild->elements();

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
         else if( fabs(storage->M[k]) < tolerance1 && ixuppElems[col] != 0.0 && ixlowElems[col] != 0.0
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

   resetRedCounters();
   nelims += doSingletonRowsA();
   nelims += doSingletonRowsC();

   cout<<"Objective offset is "<<objOffset<<endl;
   return nelims;
}

int StochPresolver::initSingletonRows(SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1)
      iAmDistrib = true;

   int nSingletonRows = 0;

   if( system_type == EQUALITY_SYSTEM )
   {
      assert(singletonRows.size() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(nRowElemsA->vec);
      int nSingleRowsA0block =  initSingletonRowsBlock(-1, nRowASimple);
      if( myRank == 0 )
         nSingletonRows += nSingleRowsA0block;

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
   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nSingletonRows, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   cout<<"singletonRows: "<<endl;
   for(size_t i= 0; i<singletonRows.size(); i++)
      cout<<singletonRows[i]<<endl;

   cout<<"blocks: "<<endl;
   for(int i= 0; i<nChildren+3; i++)
      cout<<blocks[i]<<endl;


   return nSingletonRows;
}

int StochPresolver::initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple)
{
   int nSingletonRows = 0;

   blocks[it + 1] = singletonRows.size();
   double* nnzRow = nnzRowSimple->elements();

   for( int i = 0; i < nnzRowSimple->n; i++ )
      if( nnzRow[i] == 1.0 )
      {
         singletonRows.push_back(i);
         nSingletonRows++;
      }
   return nSingletonRows;
}

bool StochPresolver::doSingletonRowsA()
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1) iAmDistrib = true;

   // First, remove enries in the singleton rows
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

   updateCurrentPointersForSingletonRow(-1, EQUALITY_SYSTEM);

   bool possFeas = procSingletonRow(matrix, -1);
   if( !possFeas ) return false;

   assert(nChildren == (int)matrix.children.size());
   for( int it = 0; it < nChildren; it++ )
   {
      colBlocksChildren[it] = colAdaptChildren.size();
      if( updateCurrentPointersForSingletonRow(it, EQUALITY_SYSTEM) )
      {
         possFeas = procSingletonRow(matrix, it);
         if( !possFeas ) return false;
      }
   }
   colBlocksChildren[matrix.children.size()] = colAdaptChildren.size();

   // empty the singletonRow list
   for(int i = 0; i<(int)singletonRows.size(); i++)
      assert(singletonRows[i] == -1);
   singletonRows.clear();
   resetBlocks();

   // Update the nnz Vectors using the reduction Vectors
   updateNnzUsingReductions( nRowElemsA->vec, redRowA->vec);
   for(int i = 0; i<nChildren; i++)
   {
      updateNnzUsingReductions( dynamic_cast<StochVector*>((*nRowElemsA).children[i])->vec, dynamic_cast<SimpleVector*>(redRowA->children[i]->vec));
      updateNnzUsingReductions( dynamic_cast<StochVector*>((*nColElems).children[i])->vec, dynamic_cast<SimpleVector*>(redCol->children[i]->vec));
   }

   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions(nColElems->vec, redCol->vec);

   // Set the variable bounds to the corresponding value
   setRemovedVarsBoundsToValues();

   //cout<<"After setRemovedVarsBoundsToZero, before applyColAdapt:"<<endl;
   //presProb->writeToStreamDense(cout);

   // Then remove the corresponding entries to these variables (go through the columns)

   int newSREq = 0;
   int newSRIneq = 0;
   // todo AAAA Fehlerquelle hier drin:
   applyColAdapt(newSREq, newSRIneq);

   if( myRank == 0 )
      cout<<"Found "<<newSREq<<" new singleton rows in A and "<<newSRIneq<<" new singleton rows in C."<<endl;

   return true;
}

/** Update the current pointers for the singleton row routine.
 * If it==-1, we are at parent block. Else, et child[it].
 * Return false if child[it] is a dummy child. */
bool StochPresolver::updateCurrentPointersForSingletonRow(int it, SystemType system_type)
{
   //todo: for INEQUALITY_SYSTEM

   // current pointers the same for all children and parent (parent part of col vector)
   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);
   currIxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).vec);
   currIxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).vec);
   currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);

   if( it == -1 )
   {
      currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
      currgChild = NULL;
      currRedRow = dynamic_cast<SimpleVector*>(redRowA->vec);
      currRedColChild = NULL;
   }
   else  // at child it
   {
      if( nColElems->children[it]->isKindOf(kStochDummy))
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->isKindOf(kStochDummy) );
         assert( redRowA->children[it]->isKindOf(kStochDummy) );
         setCurrentPointersToNull();
         return false;
      }
      currxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->vec);
      currxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->vec);
      currIxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[it]->vec);
      currIxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[it]->vec);
      currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
      currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec);
      currRedRow = dynamic_cast<SimpleVector*>(redRowA->children[it]->vec);
   }
   return true;
}

bool StochPresolver::procSingletonRow(StochGenMatrix& stochMatrix, int it)
{
   bool possFeas = true;

   if( it == -1 ) // we are at parent. So only consider Bmat
   {
      SparseStorageDynamic& B0_mat = stochMatrix.Bmat->getStorageDynamicRef();

      assert( colAdaptParent.size() == 0 );

      for(int i = blocks[it+1]; i<blocks[it+2]; i++)
      {
         int rowIdx = singletonRows[i];
         singletonRows[i] = -1;  // for debugging purposes

         possFeas = removeSingleRowEntry(B0_mat, rowIdx, LINKING_VARS_BLOCK, true);
         if( !possFeas ) return false;
      }
   }
   else  // else, at child it. Consider Amat and Bmat
   {
      //colBlocksChildren[it] = colAdaptChildren.size();
      SparseStorageDynamic& A_mat = stochMatrix.children[it]->Amat->getStorageDynamicRef();
      SparseStorageDynamic& B_mat = stochMatrix.children[it]->Bmat->getStorageDynamicRef();

      for(int i = blocks[it+1]; i<blocks[it+2]; i++)
      {
         int rowIdx = singletonRows[i];
         singletonRows[i] = -1;  // for debugging purposes
         if( A_mat.rowptr[rowIdx].start +1 == A_mat.rowptr[rowIdx].end) // single entry in in the linking vars block Amat
            possFeas = removeSingleRowEntry(A_mat, rowIdx, LINKING_VARS_BLOCK, false);
         else   // single entry in in the child block Bmat
         {
            assert(B_mat.rowptr[rowIdx].start +1 == B_mat.rowptr[rowIdx].end);
            possFeas = removeSingleRowEntry(B_mat, rowIdx, CHILD_BLOCK, false);
         }
         if( !possFeas ) return false;
      }
   }

   // todo: consider Blmat
   return true;
}

/** Removes the single entry in row rowIdx, in either Amat or Bmat (depending on block_type).
 * Sets the rhs to 0, adds objective offset, increments redRow and redCol count.
 * Fixes the variable to 0 by setting both lower and uppder bound to 0.
 * Stores the corresponding column index in colAdaptParent or colAdaptChild respectively.
 * Return false if infeasibility is detected.*/
bool StochPresolver::removeSingleRowEntry(SparseStorageDynamic& storage, int rowIdx, BlockType block_type, bool parentZero)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   double* ixlow;
   double* ixupp;
   double* xlow;
   double* xupp;
   double* g;

   if( block_type == LINKING_VARS_BLOCK )
   {
      ixlow = currIxlowParent->elements();
      ixupp = currIxuppParent->elements();
      xlow = currxlowParent->elements();
      xupp = currxuppParent->elements();
      g = currgParent->elements();
   }
   else
   {
      assert(block_type == CHILD_BLOCK);
      ixlow = currIxlowChild->elements();
      ixupp = currIxuppChild->elements();
      xlow = currxlowChild->elements();
      xupp = currxuppChild->elements();
      g = currgChild->elements();
   }

   int indexK = storage.rowptr[rowIdx].start;
   int colIdx = storage.jcolM[indexK];

   assert(storage.rowptr[rowIdx].start +1 == storage.rowptr[rowIdx].end);

   double aik = storage.M[indexK];
   assert(aik != 0.0);
   cout<<"a_ik = "<<aik<<" at ("<<rowIdx<<" ,"<<colIdx<<" )" <<endl;

   double val = currEqRhs->elements()[rowIdx] / aik;

   if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
         || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
   {
      cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
      return false;
   }
   else{
      if( !parentZero || myRank == 0 )
         objOffset += g[colIdx] * val;

      storage.rowptr[rowIdx].end --;

      COLUMNTOADAPT colWithVal = {colIdx, val};
      if( block_type == LINKING_VARS_BLOCK )
         colAdaptParent.push_back(colWithVal);
      else
      {
         assert(block_type == CHILD_BLOCK);
         colAdaptChildren.push_back(colWithVal);
      }
   }

   return true;
}

/** For the removed entries in the singleton-row procedure, adapt the bounds for the concerned variables.
 * Set their lower and upper bounds (xlow/xupp) to 0.0 and their ixlow/ixupp to 1.0 . */
void StochPresolver::setRemovedVarsBoundsToValues()
{
   // todo: if there are same columns in the parent col block, might show infeasibility?

   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);
   currIxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).vec);
   currIxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).vec);

   for(int i=0; i<(int)colAdaptParent.size(); i++)
   {
      int colIdx = colAdaptParent[i].colIdx;
      double val = colAdaptParent[i].val;
      currxlowParent->elements()[colIdx] = val;
      currxuppParent->elements()[colIdx] = val;
      currIxlowParent->elements()[colIdx] = 1.0;
      currIxuppParent->elements()[colIdx] = 1.0;
   }

   for(int it=0; it<nChildren; it++)
   {
      if( nColElems->children[it]->isKindOf(kStochDummy))
         continue;
      currxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->vec);
      currxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->vec);
      currIxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[it]->vec);
      currIxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[it]->vec);

      for(int i = colBlocksChildren[it]; i < colBlocksChildren[it + 1]; i++)
      {
         int colIdx = colAdaptChildren[i].colIdx;
         int val = colAdaptChildren[i].val;
         currxlowParent->elements()[colIdx] = val;
         currxuppParent->elements()[colIdx] = val;
         currIxlowParent->elements()[colIdx] = 1.0;
         currIxuppParent->elements()[colIdx] = 1.0;
      }
   }
}

void StochPresolver::applyColAdapt(int& newSREq, int& newSRIneq)
{
   // go through colAdaptParent, then colAdaptChildren and for each column index j, go through A' to find all rows containing an entry in j.
   // for these rows, eliminate the entry (in A and in A'), adapt the rhs with the value, increase redRow and redCol, decrease nnzRow. If it is now 1, add
   // this row to the singletonRow list.

   // todo: verify that no two column indices inside are the same but with different values (->infeasible)
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1) iAmDistrib = true;

   if( colAdaptParent.size() > 0 )
   {
      updateCurrentPointersForColAdapt(-1, EQUALITY_SYSTEM);
      int newSREqRoot = colAdaptLinkVars(-1, EQUALITY_SYSTEM);

      updateCurrentPointersForColAdapt(-1, INEQUALITY_SYSTEM);
      int newSRIneqRoot = colAdaptLinkVars(-1, INEQUALITY_SYSTEM);
      if( myRank == 0)
      {
         newSREq += newSREqRoot;
         newSRIneq += newSRIneqRoot;
      }
   }

   for(int i= 0; i<nChildren; i++)
   {
      if( updateCurrentPointersForColAdapt( i, EQUALITY_SYSTEM) )
      {
         newSREq += colAdaptLinkVars(i, EQUALITY_SYSTEM);
         if( colBlocksChildren[i] < colBlocksChildren[i+1] )
            newSREq += colAdaptChild( i, EQUALITY_SYSTEM);
      }

      if( updateCurrentPointersForColAdapt( i, INEQUALITY_SYSTEM) )
      {
         newSRIneq += colAdaptLinkVars(i, INEQUALITY_SYSTEM);
         if( colBlocksChildren[i] < colBlocksChildren[i+1] )
            newSRIneq += colAdaptChild( i, INEQUALITY_SYSTEM);
      }
   }

   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      int* newSR = new int[2];
      newSR[0] = newSREq;
      newSR[1] = newSRIneq;
      MPI_Allreduce(MPI_IN_PLACE, newSR, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      newSREq = newSR[0];
      newSRIneq = newSR[1];
   }

   updateNnzUsingReductions(nColElems->vec, redCol->vec);

   // reset colAdaptParent, colAdaptChildren, colBlocksChildren
   colAdaptParent.clear();
   colAdaptChildren.clear();
   resetColBlocks();
}

/** Adapt the columns for the linking-variable-blocks (the A_i) blocks */
int StochPresolver::colAdaptLinkVars(int it, SystemType system_type)
{
   // todo: adapt linking block
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int newSingletonRows = 0;

   if( system_type == EQUALITY_SYSTEM )
   {
      if( it==-1) assert(singletonRows.size()==0);
      blocks[it+1] = singletonRows.size();
   }
   else
   {
      if( it==-1) assert(singletonRowsIneq.size()==0);
      blocksIneq[it+1] = singletonRowsIneq.size();
   }

   for( int i=0; i < (int)colAdaptParent.size(); i++)
   {
      int colIdxA = colAdaptParent[i].colIdx;
      double val = colAdaptParent[i].val;

      for( int j = currAmatTrans->rowptr[colIdxA].start; j < currAmatTrans->rowptr[colIdxA].end; j++ )
      {
         int rowIdxA = currAmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currAmat, rowIdxA, colIdxA);

         if( system_type == EQUALITY_SYSTEM )
            currEqRhs->elements()[rowIdxA] -= m * val;
         else
         {  // for inequality system, adapt both sides if they exist
            assert(system_type ==INEQUALITY_SYSTEM);
            if( currIcupp->elements()[rowIdxA] != 0.0 )
               currIneqRhs->elements()[rowIdxA] -= m * val;
            if( currIclow->elements()[rowIdxA] != 0.0 )
               currIneqLhs->elements()[rowIdxA] -=  m * val;
         }

         if( it > -1 || myRank == 0 )
            currRedColParent->elements()[colIdxA] ++;
         currNnzRow->elements()[rowIdxA] --;
         currRedRow->elements()[rowIdxA] ++;

         if(currNnzRow->elements()[rowIdxA] == 1)
         {
            newSingletonRows++;
            if( system_type == EQUALITY_SYSTEM )
               singletonRows.push_back(rowIdxA);
            else
               singletonRowsIneq.push_back(rowIdxA);
         }
      }
      clearRow(*currAmatTrans, colIdxA);
   }

   return newSingletonRows;
}

int StochPresolver::colAdaptChild( int it, SystemType system_type )
{
   // Bmat, BmatTrans, ncolchild, redcolchild, nrowchild, redrowchild, currEqRhs, icupp, cupp, iclow, clow are updated
   // todo: Linking Block Blmat ?
   int newSingletonRows = 0;

   for( int i = colBlocksChildren[it]; i<colBlocksChildren[it+1]; i++)
   {
      int colIdxB = colAdaptChildren[i].colIdx;
      double val = colAdaptChildren[i].val;

      // Block Bmat:
      for( int j = currBmatTrans->rowptr[colIdxB].start; j<currBmatTrans->rowptr[colIdxB].end; j++ )
      {
         int rowIdxB = currBmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currBmat, rowIdxB, colIdxB);

         if( system_type == EQUALITY_SYSTEM )
            currEqRhs->elements()[rowIdxB] -= m * val;
         else
         {  // for inequality system, adapt both sides if they exist
            assert(system_type ==INEQUALITY_SYSTEM);
            if( currIcupp->elements()[rowIdxB] != 0.0 )
               currIneqRhs->elements()[rowIdxB] -= m * val;
            if( currIclow->elements()[rowIdxB] != 0.0 )
               currIneqLhs->elements()[rowIdxB] -=  m * val;
         }

         currNnzColChild->elements()[colIdxB] --;
         currRedColChild->elements()[colIdxB] ++;
         currNnzRow->elements()[rowIdxB] --;
         currRedRow->elements()[rowIdxB] ++;
         assert( currNnzColChild->elements()[colIdxB] >= 0.0 && currNnzRow->elements()[rowIdxB] >= 0.0 );

         if(currNnzRow->elements()[rowIdxB] == 1)
         {
            newSingletonRows++;
            if( system_type == EQUALITY_SYSTEM )
               singletonRows.push_back(rowIdxB);
            else
               singletonRowsIneq.push_back(rowIdxB);
         }
      }
      clearRow(*currBmatTrans, colIdxB);
   }
   return newSingletonRows;
}

/** Return false if it is a dummy child. Else, set the current pointers to Amat, AmatTrans, Bmat, BmatTrans,
 * currNnzRow, currRedRow, currNnzColChild, currRedColChild.
 * currxlowChild, currxuppChild, currIxlowChild, currIxuppChild.
 * And depending on the systemType: either currEqRhs or currIneqRhs, currIneqLhs, currIcupp and currIclow. */
bool StochPresolver::updateCurrentPointersForColAdapt(int it, SystemType system_type)
{
   setCurrentPointersToNull();
   if( it == -1 )
   {
      if( system_type == EQUALITY_SYSTEM )
      {
         StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamic();
         currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamicTransposed();
         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsA->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowA->vec);
      }
      else  // system_type == INEQUALITY_SYSTEM
      {
         StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));

         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamic();
         currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamicTransposed();
         currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec);
         currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vec);
         currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vec);
         currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowC->vec);
      }
      currRedColParent = dynamic_cast<SimpleVector*>(redCol->vec);

   }
   else{ // at child it
      if( system_type == EQUALITY_SYSTEM )
      {
         StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

         if( childIsDummy(matrix, it, EQUALITY_SYSTEM) ) return false;

         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
         currBmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamicTransposed();
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamicTransposed();

         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsA->children[it]->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowA->children[it]->vec);
      }
      else  // system_type == INEQUALITY_SYSTEM
      {
         assert( system_type == INEQUALITY_SYSTEM );

         StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));

         if( childIsDummy(matrix, it, INEQUALITY_SYSTEM) ) return false;

         currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
         currBmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamicTransposed();
         currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
         currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamicTransposed();

         currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec);
         currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec);
         currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec);
         currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec);
         currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->children[it]->vec);
         currRedRow = dynamic_cast<SimpleVector*>(redRowC->children[it]->vec);
      }
      currNnzColChild = dynamic_cast<SimpleVector*>(nColElems->children[it]->vec);
      currRedColChild = dynamic_cast<SimpleVector*>(redCol->children[it]->vec);
      currRedColParent = dynamic_cast<SimpleVector*>(redCol->vec);
   }

   return true;
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
   {
      blocks[i] = 0;
      blocksIneq[i] = 0;
   }
}

void StochPresolver::resetColBlocks()
{
   for( int i = 0; i < nChildren+1; i++)
      colBlocksChildren[i] = 0;
}

void StochPresolver::resetRedCounters()
{
   redRowA->setToZero();
   redRowC->setToZero();
   redCol->setToZero();
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

double StochPresolver::removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx)
{
   int i;
   int end = storage.rowptr[rowIdx].end;
   int start = storage.rowptr[rowIdx].start;

   if( start == end )
      return 0.0;
   for( i=storage.rowptr[rowIdx].start; i<end; i++)
   {
      if( storage.jcolM[i] == colIdx )
                  break;
   }
   double m = storage.M[i];
   std::swap(storage.M[i],storage.M[end-1]);
   std::swap(storage.jcolM[i],storage.jcolM[end-1]);
   storage.rowptr[rowIdx].end --;

   return m;
}

void StochPresolver::clearRow(SparseStorageDynamic& storage, const int rowIdx)
{
   storage.rowptr[rowIdx].end = storage.rowptr[rowIdx].start;
}

bool StochPresolver::childIsDummy(StochGenMatrix& matrix, int it, SystemType system_type)
{
   if( matrix.children[it]->isKindOf(kStochGenDummyMatrix))
   {
      assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->isKindOf(kStochDummy) );
      assert( redCol->children[it]->isKindOf(kStochDummy) );

      if( system_type == EQUALITY_SYSTEM)
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->isKindOf(kStochDummy) );
         assert( nRowElemsA->children[it]->isKindOf(kStochDummy) );
         assert( redRowA->children[it]->isKindOf(kStochDummy) );
      }
      else
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->isKindOf(kStochDummy) );
         assert( nRowElemsC->children[it]->isKindOf(kStochDummy) );
         assert( redRowC->children[it]->isKindOf(kStochDummy) );
      }
      setCurrentPointersToNull();
      return true;
   }
   return false;
}


