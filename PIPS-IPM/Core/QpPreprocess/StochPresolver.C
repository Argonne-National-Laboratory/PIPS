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

   if(redRowA->vecl != NULL)  //todo correct test for linking?
      currEqRhsAdaptionsLink = new double[redRowA->vecl->n];
   if(redRowC->vecl != NULL)  //todo correct test for linking?
   {
      currInEqRhsAdaptionsLink = new double[redRowC->vecl->n];
      currInEqLhsAdaptionsLink = new double[redRowC->vecl->n];
   }

   localNelims = 0;
   nChildren = nColElems->children.size();
   linkVarsBlocks = new BLOCKS[nChildren + 1];
   childBlocks = new BLOCKS[nChildren + 1];
   resetLinkvarsAndChildBlocks();

   blocks = new int[nChildren + 3];
   blocksIneq = new int[nChildren + 3];
   resetBlocks();

   objOffset = 0.0;
}


StochPresolver::~StochPresolver()
{
   delete[] linkVarsBlocks;
   delete[] childBlocks;
   delete[] blocks;
   delete[] blocksIneq;
   delete[] currEqRhsAdaptionsLink;
   delete[] currInEqRhsAdaptionsLink;
   delete[] currInEqLhsAdaptionsLink;
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

   //doSingletonRows();

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
   currBlmatTrans = NULL;
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
   currEqRhsLink = NULL;
   currIneqRhsLink = NULL;
   currIneqLhsLink = NULL;
   currIcuppLink = NULL;
   currIclowLink = NULL;
   currNnzRow = NULL;
   currRedRow = NULL;
   currRedRowLink = NULL;
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
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1)
      iAmDistrib = true;

   int nelims = 0;
   int newSREq = initSingletonRows(EQUALITY_SYSTEM);
   int newSRIneq = 0;
   cout<<"Found "<<newSREq<<" singleton rows in equality system A."<<endl;

   int iter = 0;
   while( newSREq > 0 && iter < maxIterSR)
   {
      if( iter > 0 )
         newSREq = initSingletonRows(EQUALITY_SYSTEM);
      resetRedCounters();
      nelims += doSingletonRowsA(newSREq, newSRIneq);
      //nelims += doSingletonRowsC();
      if( myRank == 0 )
         cout<<"Found new singleton rows that were just created: "<<newSREq<<" in A, "<<newSRIneq<<" in C."<<endl;
      iter++;
   }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &objOffset, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if( myRank == 0 )
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

bool StochPresolver::doSingletonRowsA(int& newSREq, int& newSRIneq)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1) iAmDistrib = true;

   newSREq = 0;
   newSRIneq = 0;
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

   updateCurrentPointersForSingletonRow(-1, EQUALITY_SYSTEM);
   bool possFeas = procSingletonRowRoot(matrix);
   if( !possFeas ) return false;

   assert(nChildren == (int)matrix.children.size());
   for( int it = 0; it < nChildren; it++ )
   {
      if( updateCurrentPointersForSingletonRow(it, EQUALITY_SYSTEM) )
      {
         possFeas = procSingletonRowChild(matrix, it, newSREq, newSRIneq);
         if( !possFeas ) return false;
      }
   }

   // Update: redRowLink and lhs/rhs (Linking part) of both systems:
   updateRhsNRowLink();

   possFeas = combineColAdaptParent();
   if( !possFeas ) return false;

   // apply updated colAdaptParent to the Amat blocks
   for(int i= -1; i<nChildren; i++)
   {
      if( updateCurrentPointersForColAdapt( i, EQUALITY_SYSTEM) )
         newSREq += colAdaptLinkVars(i, EQUALITY_SYSTEM);
      if( updateCurrentPointersForColAdapt( i, INEQUALITY_SYSTEM) )
         newSRIneq += colAdaptLinkVars(i, INEQUALITY_SYSTEM);
   }

   // apply updated colAdaptParent to the F0 block (Blmat of parent):
   if( updateCPforColAdaptF0( EQUALITY_SYSTEM ) )
      colAdaptF0( EQUALITY_SYSTEM);
   if( updateCPforColAdaptF0( INEQUALITY_SYSTEM) )
      colAdaptF0( INEQUALITY_SYSTEM);
   colAdaptParent.clear();

   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions(nColElems->vec, redCol->vec);

   resetRedCounters();

   // empty the singletonRow list or todo: store the new singleton rows during the process
   for(int i = 0; i<(int)singletonRows.size(); i++)
      assert(singletonRows[i] == -1);
   singletonRows.clear();
   resetBlocks();

   if( iAmDistrib )
   {
      int* newSR = new int[2];
      newSR[0] = newSREq;
      newSR[1] = newSRIneq;
      MPI_Allreduce(MPI_IN_PLACE, newSR, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      newSREq = newSR[0];
      newSRIneq = newSR[1];
      delete[] newSR;
   }

   /*
   // print F_1 (blmat vom letzten child) for debugging
   SparseStorageDynamic& blmat_test = matrix.children[1]->Blmat->getStorageDynamicRef();
   cout<<"blmat_test: m,n,len: "<<blmat_test.m<<" "<<blmat_test.n<<" "<<blmat_test.len<<endl;
   blmat_test.writeToStreamDense(cout);
   cout<<"blmat_test jcolM:"<<endl;
   for(int i=0; i<blmat_test.len; i++)
      cout<<blmat_test.jcolM[i]<<endl;
   cout<<"blmat_test M:"<<endl;
   for(int i=0; i<blmat_test.len; i++)
      cout<<blmat_test.M[i]<<endl;
   cout<<"blmat_test krowptrs:"<<endl;
   for(int i=0; i<blmat_test.m; i++)
      cout<<"start: "<<blmat_test.rowptr[i].start<<" end: "<<blmat_test.rowptr[i].end<<endl;
      */

   return true;
}

/** Update the current pointers for the singleton row routine.
 * If it==-1, we are at parent block. Else, et child[it].
 * Return false if child[it] is a dummy child. */
bool StochPresolver::updateCurrentPointersForSingletonRow(int it, SystemType system_type)
{
   //todo: for INEQUALITY_SYSTEM
   setCurrentPointersToNull();

   // current pointers the same for all children and parent (parent part of col vector)
   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);
   currIxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).vec);
   currIxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).vec);
   currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);
   currRedColParent = dynamic_cast<SimpleVector*>(redCol->vec);

   if( it == -1 )
   {
      // todo: curr matrices
      currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
      currRedRow = dynamic_cast<SimpleVector*>(redRowA->vec);
      currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsA->vec);
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
      StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));
      currAmat = matrix.children[it]->Amat->getStorageDynamic();
      currAmatTrans = matrix.children[it]->Amat->getStorageDynamicTransposed();
      currBmat = matrix.children[it]->Bmat->getStorageDynamic();
      currBmatTrans = matrix.children[it]->Bmat->getStorageDynamicTransposed();

      if( hasLinking(system_type) )
      {
         currBlmat = matrix.children[it]->Blmat->getStorageDynamic();
         currBlmatTrans = matrix.children[it]->Blmat->getStorageDynamicTransposed();
         currEqRhsLink =  dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vecl);
         currRedRowLink = dynamic_cast<SimpleVector*>(redRowA->vecl);
      }

      currxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->vec);
      currxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->vec);
      currIxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[it]->vec);
      currIxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[it]->vec);
      currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
      currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec);
      currRedRow = dynamic_cast<SimpleVector*>(redRowA->children[it]->vec);
      currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsA->children[it]->vec);
      currNnzColChild = dynamic_cast<SimpleVector*>(nColElems->children[it]->vec);
      currRedColChild = dynamic_cast<SimpleVector*>(redCol->children[it]->vec);
   }
   return true;
}

bool StochPresolver::updateCPForSingletonRowInequalityBChild( int it )
{
   assert( it >= 0);
   setCurrentPointersToNull();

   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));
   if( childIsDummy(matrix, it, INEQUALITY_SYSTEM) ) return false;

   currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
   currBmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamicTransposed();

   currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec);
   currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec);
   currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec);
   currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec);
   currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->children[it]->vec);
   currRedRow = dynamic_cast<SimpleVector*>(redRowC->children[it]->vec);
   currRedColChild = dynamic_cast<SimpleVector*>(redCol->children[it]->vec);
   currNnzColChild = dynamic_cast<SimpleVector*>(nColElems->children[it]->vec);

   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      // todo: assert that all vectors and matrices have linking part
      assert(redRowC->vecl);
      currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamic();
      currBlmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamicTransposed();
      currIcuppLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vecl);
      currIclowLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vecl);
      currRedRowLink = dynamic_cast<SimpleVector*>(redRowC->vecl);
   }
   return true;
}

bool StochPresolver::procSingletonRowRoot(StochGenMatrix& stochMatrix)
{
   bool possFeas = true;

   SparseStorageDynamic& B0_mat = stochMatrix.Bmat->getStorageDynamicRef();
   assert( colAdaptParent.size() == 0 );

   for(int i = blocks[0]; i<blocks[1]; i++)
   {
      int rowIdx = singletonRows[i];
      singletonRows[i] = -1;  // for debugging purposes

      possFeas = removeSingleRowEntryB0(B0_mat, rowIdx);
      if( !possFeas ) return false;
   }

   return true;
}

/* Processing the singleton rows in child it, more precisely, goes through all singleton rows in Amat and Bmat.
 * Those in Amat are stored in colAdaptParent for later processing.
 * Those in Bmat are removed and stored in colAdaptLinkBlock. Furthermore, the corresponding fixed variables (columns)
 * in Bmat and in Blmat are removed.
 * Using this colAdaptLinkBlock, the variables (columns) are removed from the inequalities Bmat, Blmat as well.
 */
bool StochPresolver::procSingletonRowChild(StochGenMatrix& stochMatrix, int it, int& newSREq, int& newSRIneq)
{
   bool possFeas = true;

   SparseStorageDynamic& A_mat = stochMatrix.children[it]->Amat->getStorageDynamicRef();
   SparseStorageDynamic& B_mat = stochMatrix.children[it]->Bmat->getStorageDynamicRef();

   possFeas = procSingletonRowChildAmat( A_mat, it);
   if( !possFeas ) return false;

   std::vector<COLUMNTOADAPT> colAdaptLinkBlock;
   possFeas = procSingletonRowChildBmat( B_mat, it, colAdaptLinkBlock, newSREq);
   if( !possFeas ) return false;

   // using colAdaptLinkBlock, go through the columns in Blmat
   if( hasLinking(EQUALITY_SYSTEM) )
   {
      possFeas = adaptChildBlmat( colAdaptLinkBlock, EQUALITY_SYSTEM);
      if( !possFeas ) return false;
   }
   // and go through the columns in Bmat, Blmat of the inequality
   updateCPForSingletonRowInequalityBChild( it );
   possFeas = adaptInequalityChildB( colAdaptLinkBlock, newSRIneq );
   if( !possFeas ) return false;

   return true;
}

bool StochPresolver::procSingletonRowChildAmat(SparseStorageDynamic& A_mat, int it)
{
   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();
   double* g = currgParent->elements();

   for(int i = blocks[it+1]; i<blocks[it+2]; i++)
   {
      int rowIdx = singletonRows[i];
      if( A_mat.rowptr[rowIdx].start +1 == A_mat.rowptr[rowIdx].end )
      {
         singletonRows[i] = -1;  // for debugging purposes

         // store the column index with fixed value in colAdaptParent and adapt objOffset:
         int indexK = A_mat.rowptr[rowIdx].start;
         int colIdx = A_mat.jcolM[indexK];

         assert(A_mat.rowptr[rowIdx].start +1 == A_mat.rowptr[rowIdx].end);

         double aik = A_mat.M[indexK];
         assert(aik != 0.0);
         cout<<"a_ik = "<<aik<<" at ("<<rowIdx<<" ,"<<colIdx<<" )" <<endl;

         double val = currEqRhs->elements()[rowIdx] / aik;

         if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
               || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
         {
            cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
            return false;
         }
         objOffset += g[colIdx] * val;

         COLUMNTOADAPT colWithVal = {colIdx, val};
         colAdaptParent.push_back(colWithVal);
      }
   }
   return true;
}

bool StochPresolver::procSingletonRowChildBmat(SparseStorageDynamic& B_mat, int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR)
{
   for(int i = blocks[it+1]; i<blocks[it+2]; i++)
   {
      int rowIdx = singletonRows[i];
      if( rowIdx == -1 || B_mat.rowptr[rowIdx].start == B_mat.rowptr[rowIdx].end)
      {
         // entry was already in Amat or in a previous singleton row in Bmat
         continue;
      }
      else
      {
         assert( B_mat.rowptr[rowIdx].start +1 == B_mat.rowptr[rowIdx].end );
         singletonRows[i] = -1;  // for debugging purposes
         removeSingleRowEntryChildBmat(rowIdx, colAdaptLinkBlock, EQUALITY_SYSTEM, newSR);
      }
   }

   return true;
}

bool StochPresolver::removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR)
{
   double* ixlow = currIxlowChild->elements();
   double* ixupp = currIxuppChild->elements();
   double* xlow = currxlowChild->elements();
   double* xupp = currxuppChild->elements();
   double* g = currgChild->elements();

   assert(currBmat->rowptr[rowIdx].start +1 == currBmat->rowptr[rowIdx].end);

   int indexK = currBmat->rowptr[rowIdx].start;
   int colIdx = currBmat->jcolM[indexK];

   double aik = currBmat->M[indexK];
   assert(aik != 0.0);
   cout<<"a_ik = "<<aik<<" at ("<<rowIdx<<" ,"<<colIdx<<" )" <<endl;

   double val = currEqRhs->elements()[rowIdx] / aik;

   if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
         || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
   {
      cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
      return false;
   }
   objOffset += g[colIdx] * val;

   // adapt the col, val immediately in this block B_i and store them for the Blmat
   newSR += adaptChildBmatCol(colIdx, val, system_type);

   COLUMNTOADAPT colWithVal = {colIdx, val};
   colAdaptLinkBlock.push_back(colWithVal);

   return true;
}

/** For the given column index and the value to which this variable is to be fixed,
 * go through the current Bmat by going through the corresponding row of the transposed Bmat.
 * Adapt the rhs, the nnzRow, the nnzCol and all concerned entries in Bmat and BmatTransposed.
 * Return the number of newly found singleton rows.
 */
int StochPresolver::adaptChildBmatCol(int colIdx, double val, SystemType system_type)
{
   int newSingletonRows = 0;

   for( int j = currBmatTrans->rowptr[colIdx].start; j<currBmatTrans->rowptr[colIdx].end; j++ )
   {
      int rowIdxB = currBmatTrans->jcolM[j];
      double m = removeEntryInDynamicStorage(*currBmat, rowIdxB, colIdx);

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

      currNnzColChild->elements()[colIdx] --;
      currRedColChild->elements()[colIdx] ++;
      currNnzRow->elements()[rowIdxB] --;
      currRedRow->elements()[rowIdxB] ++;
      assert( currNnzColChild->elements()[colIdx] >= 0.0 && currNnzRow->elements()[rowIdxB] >= 0.0 );

      if(currNnzRow->elements()[rowIdxB] == 1)
         newSingletonRows++;
   }
   clearRow(*currBmatTrans, colIdx);

   return newSingletonRows;
}

/** Given the vector<COLUMNTOADAPT>, both the block Bmat and Blat (if existent) are updated accordingly.
 */
bool StochPresolver::adaptInequalityChildB( std::vector<COLUMNTOADAPT> & colAdaptBblock, int& newSRIneq )
{
   // Bmat blocks
   bool possFeas = adaptChildBmat( colAdaptBblock, INEQUALITY_SYSTEM, newSRIneq);
   if( !possFeas ) return false;

   // Blmat blocks
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      bool possFeas = adaptChildBlmat( colAdaptBblock, INEQUALITY_SYSTEM);
      if( !possFeas ) return false;
   }

   return true;
}
/** Removes the single entry in row rowIdx in Bmat_0.
 * Stores the corresponding column index in colAdaptParent.
 * Return false if infeasibility is detected.*/
bool StochPresolver::removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();
   double* g = currgParent->elements();

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
   else
   {
      if( myRank == 0 )
      {
         COLUMNTOADAPT colWithVal = {colIdx, val};
         bool uniqueAdditionToOffset = true;
         for(int i=0; i<(int)colAdaptParent.size(); i++)
         {
            if( colAdaptParent[i].colIdx == colIdx )
            {
               if( colAdaptParent[i].val != val )
               {
                  cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
                  return false;
               }
               uniqueAdditionToOffset = false;
            }

         }
         if( uniqueAdditionToOffset )
         {
            colAdaptParent.push_back(colWithVal);
            objOffset += g[colIdx] * val;
         }
      }
   }

   return true;
}

/*
 * Given a vector<COLUMNTOADAPT>, this routine goes through all columns inside and removes them
 * from the current Bmat block. Depending on the system_type, rhs (and lhs) are updated.
 */
bool StochPresolver::adaptChildBmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type, int& newSR )
{
   for(int i=0; i<(int)colAdaptBlock.size(); i++)
   {
      int colIdx = colAdaptBlock[i].colIdx;
      double val = colAdaptBlock[i].val;

      for( int j = currBmatTrans->rowptr[colIdx].start; j<currBmatTrans->rowptr[colIdx].end; j++ )
      {
         int rowIdx = currBmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currBmat, rowIdx, colIdx);

         if( system_type == EQUALITY_SYSTEM )
         {
            currEqRhs->elements()[rowIdx] -= m * val;
         }
         else
         {
            assert(system_type ==INEQUALITY_SYSTEM);
            if( currIcupp->elements()[rowIdx] != 0.0 )
               currIneqRhs->elements()[rowIdx] -= m * val;
            if( currIclow->elements()[rowIdx] != 0.0 )
               currIneqLhs->elements()[rowIdx] -=  m * val;

         }
         currRedRow->elements()[rowIdx] ++;
         currNnzRow->elements()[rowIdx] --;
         currNnzColChild->elements()[colIdx] --;
         currRedColChild->elements()[colIdx] ++;

         if( currNnzRow->elements()[rowIdx] == 1.0 )
            newSR++;

         assert( currNnzColChild->elements()[colIdx] >= 0.0 );
      }
      clearRow(*currBmatTrans, colIdx);
   }
   return true;
}

bool StochPresolver::adaptChildBlmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type)
{
   assert(currBlmat != NULL);
   if( system_type == EQUALITY_SYSTEM )
   {
      assert( redRowA->vecl->n == currRedRowLink->n );
      assert( redRowA->vecl->n == currBlmat->m );
   }
   else
   {
      assert(system_type == INEQUALITY_SYSTEM);
      assert( redRowC->vecl->n == currRedRowLink->n );
      assert( redRowC->vecl->n == currBlmat->m );
   }

   for(int i=0; i<(int)colAdaptBlock.size(); i++)
   {
      int colIdx = colAdaptBlock[i].colIdx;
      double val = colAdaptBlock[i].val;

      for( int j = currBlmatTrans->rowptr[colIdx].start; j<currBlmatTrans->rowptr[colIdx].end; j++ )
      {
         int rowIdx = currBlmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currBlmat, rowIdx, colIdx);

         if( system_type == EQUALITY_SYSTEM )
         {
            currEqRhsAdaptionsLink[rowIdx] -= m * val;
         }
         else
         {
            if( currIcuppLink->elements()[rowIdx] != 0.0 )
               currInEqRhsAdaptionsLink[rowIdx] -= m * val;
            if( currIclowLink->elements()[rowIdx] != 0.0 )
               currInEqLhsAdaptionsLink[rowIdx] -=  m * val;

         }
         currRedRowLink->elements()[rowIdx] ++;
         currNnzColChild->elements()[colIdx] --;
         currRedColChild->elements()[colIdx] ++;

         assert( currNnzColChild->elements()[colIdx] >= 0.0 );
      }
      clearRow(*currBlmatTrans, colIdx);
   }
   return true;
}

/** Adapt the columns for the linking-variable-blocks (the A_i) blocks */
int StochPresolver::colAdaptLinkVars(int it, SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int newSingletonRows = 0;

   for( int i=0; i < (int)colAdaptParent.size(); i++)
   {
      int colIdxA = colAdaptParent[i].colIdx;
      double val = colAdaptParent[i].val;

      for( int j = currAmatTrans->rowptr[colIdxA].start; j < currAmatTrans->rowptr[colIdxA].end; j++ )
      {
         int rowIdxA = currAmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currAmat, rowIdxA, colIdxA);
         cout<<"Removed entry "<<rowIdxA<<", "<<colIdxA<<" with value "<<m<<" in Amat of child "<<it<<endl;

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

         if( currNnzRow->elements()[rowIdxA] == 1.0 )
            if( it > -1 || myRank == 0 )  // damit die newSR von B_0 nicht doppelt gezählt werden
               newSingletonRows++;
      }
      clearRow(*currAmatTrans, colIdxA);
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

bool StochPresolver::updateCPforColAdaptF0( SystemType system_type )
{
   setCurrentPointersToNull();
   if( !hasLinking(system_type) )
      return false;

   if( system_type == EQUALITY_SYSTEM )
   {
      StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));
      currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamic();
      currBlmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamicTransposed();
      currEqRhsLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vecl);
      currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsA->vecl);
   }
   else  // system_type == INEQUALITY_SYSTEM
   {
      StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));

      currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamic();
      currBlmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamicTransposed();
      currIneqRhsLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vecl);
      currIneqLhsLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vecl);
      currIcuppLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vecl);
      currIclowLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vecl);
      currNnzRow = dynamic_cast<SimpleVector*>(nRowElemsC->vecl);
   }
   return true;
}

int StochPresolver::colAdaptF0(SystemType system_type)
{
   assert(currBlmat != NULL);
   assert( currNnzRow->n == currBlmat->m );

   int newSingletonRows = 0;

   for(int i=0; i<(int)colAdaptParent.size(); i++)
   {
      int colIdx = colAdaptParent[i].colIdx;
      double val = colAdaptParent[i].val;

      for( int j = currBlmatTrans->rowptr[colIdx].start; j<currBlmatTrans->rowptr[colIdx].end; j++ )
      {
         int rowIdx = currBlmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currBlmat, rowIdx, colIdx);

         if( system_type == EQUALITY_SYSTEM )
            currEqRhsLink->elements()[rowIdx] -= m * val;
         else
         {
            assert(system_type == INEQUALITY_SYSTEM);
            if( currIcuppLink->elements()[rowIdx] != 0.0 )
               currIneqRhsLink->elements()[rowIdx] -= m * val;
            if( currIclowLink->elements()[rowIdx] != 0.0 )
               currIneqLhsLink->elements()[rowIdx] -=  m * val;
         }
         dynamic_cast<SimpleVector*>(nColElems->vec)->elements()[colIdx] --;
         currNnzRow->elements()[rowIdx] --;

         assert( currNnzRow->elements()[rowIdx] >= 0.0 );
         assert( dynamic_cast<SimpleVector*>(nColElems->vec)->elements()[colIdx] >= 0.0 );

         if(currNnzRow->elements()[rowIdx] == 1.0)
            newSingletonRows++;
      }
      clearRow(*currBlmatTrans, colIdx);
   }
   return newSingletonRows;
}

bool StochPresolver::combineColAdaptParent()
{
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1) iAmDistrib = true;

   if( iAmDistrib )
   {
      // allgather the length of each colAdaptParent
      int mylen = colAdaptParent.size();
      int* recvcounts = new int[world_size];

      MPI_Allgather(&mylen, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

      // allgatherv the actual colAdaptParents
      // First, extract the colIdx and val into int* and double* arrays:
      int* colIndicesLocal = new int[mylen];
      double* valuesLocal = new double[mylen];
      for(int i=0; i<mylen; i++)
      {
         colIndicesLocal[i] = colAdaptParent[i].colIdx;
         valuesLocal[i] = colAdaptParent[i].val;
      }
      // Second, prepare the receive buffers:
      int lenghtGlobal = recvcounts[0];
      int* displs = new int[world_size];
      displs[0] = 0;
      for(int i=1; i<world_size; i++)
      {
         lenghtGlobal += recvcounts[i];
         displs[i] = displs[i-1] + recvcounts[i-1];
      }
      int* colIndicesGlobal = new int[lenghtGlobal];
      double* valuesGlobal = new double[lenghtGlobal];
      // Then, do the actual MPI communication:
      MPI_Allgatherv(colIndicesLocal, mylen, MPI_INT, colIndicesGlobal, recvcounts, displs , MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(valuesLocal, mylen, MPI_DOUBLE, valuesGlobal, recvcounts, displs , MPI_DOUBLE, MPI_COMM_WORLD);

      // Reconstruct a colAdaptParent:
      colAdaptParent.clear();
      for(int i=0; i<lenghtGlobal; i++)
      {
         COLUMNTOADAPT colWithVal = {colIndicesGlobal[i], valuesGlobal[i]};
         colAdaptParent.push_back(colWithVal);
      }

      delete[] recvcounts;
      delete[] colIndicesLocal;
      delete[] valuesLocal;
      delete[] colIndicesGlobal;
      delete[] valuesGlobal;
   }

   // Sort colIndicesGlobal (and valuesGlobal accordingly), remove duplicates and find infeasibilities
   std::sort(colAdaptParent.begin(), colAdaptParent.end(), col_is_smaller());

   if(colAdaptParent.size() > 0)
   {
      int colIdxCurrent = colAdaptParent[0].colIdx;
      double valCurrent = colAdaptParent[0].val;
      for(int i=1; i<(int)colAdaptParent.size(); i++)
      {
         if( colAdaptParent[i].colIdx == colIdxCurrent )
         {
            if( colAdaptParent[i].val != valCurrent )
               return false;
            else
               colAdaptParent.erase(colAdaptParent.begin()+i);   //todo: implement more efficiently
         }
         else{
            colIdxCurrent = colAdaptParent[i].colIdx;
            valCurrent = colAdaptParent[i].val;
         }
      }
   }
   assert( (int)colAdaptParent.size() <= nColElems->vec->n );

   return true;
}

void StochPresolver::updateRhsNRowLink()
{
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   bool iAmDistrib = false;
   if( world_size > 1) iAmDistrib = true;

   if( hasLinking(EQUALITY_SYSTEM) )
   {
      currEqRhsLink =  dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vecl);
      if( iAmDistrib )
      {
         // todo improve mpi communication, only one:
         double* redRowLink = dynamic_cast<SimpleVector*>(redRowA->vecl)->elements();
         int message_sizeA = dynamic_cast<SimpleVector*>(redRowA->vecl)->length();
         MPI_Allreduce(MPI_IN_PLACE, redRowLink, message_sizeA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, currEqRhsAdaptionsLink, message_sizeA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      // update nRowElemsA.vecl
      updateNnzUsingReductions(nRowElemsA->vecl, redRowA->vecl);
      // update rhs with += adaptionsRhsLink
      for(int i=0; i<currEqRhsLink->n; i++)
         currEqRhsLink->elements()[i] += currEqRhsAdaptionsLink[i];

   }
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      currIneqRhsLink =  dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vecl);
      currIneqLhsLink =  dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vecl);
      currIcuppLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vecl);
      currIclowLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vecl);

      if( iAmDistrib )
      {
         // todo improve mpi communication, only one:
         int message_sizeC = dynamic_cast<SimpleVector*>(redRowC->vecl)->length();
         double* redRowLinkC = dynamic_cast<SimpleVector*>(redRowC->vecl)->elements();
         MPI_Allreduce(MPI_IN_PLACE, redRowLinkC, message_sizeC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, currInEqRhsAdaptionsLink, message_sizeC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, currInEqLhsAdaptionsLink, message_sizeC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      // update nRowElemsC.vecl
      updateNnzUsingReductions(nRowElemsC->vecl, redRowC->vecl);
      // update rhs/lhs with += adaptionsRhsLink
      for(int i=0; i<currIneqRhsLink->n; i++)
      {
         if( currIcuppLink->elements()[i] != 0.0 )
            currIneqRhsLink->elements()[i] += currInEqRhsAdaptionsLink[i];
         if( currIclowLink->elements()[i] != 0.0 )
            currIneqLhsLink->elements()[i] += currInEqLhsAdaptionsLink[i];
      }
   }
   resetRhsAdaptionsLink();
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

void StochPresolver::resetRedCounters()
{
   redRowA->setToZero();
   redRowC->setToZero();
   redCol->setToZero();
}

void StochPresolver::resetRhsAdaptionsLink()
{
   for( int i = 0; i < redRowA->vecl->n; i++)
      currEqRhsAdaptionsLink[i] = 0.0;

   for( int i = 0; i < redRowC->vecl->n; i++ )
   {
      currInEqRhsAdaptionsLink[i] = 0.0;
      currInEqLhsAdaptionsLink[i] = 0.0;
   }
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

bool StochPresolver::hasLinking(SystemType system_type)
{
   int mlink, nlink;
   if( system_type == EQUALITY_SYSTEM )
   {
      dynamic_cast<StochGenMatrix&>(*(presProb->A)).Blmat->getSize(mlink, nlink);
      if( mlink > 0 )
         return true;
   }
   else
   {
      dynamic_cast<StochGenMatrix&>(*(presProb->C)).Blmat->getSize(mlink, nlink);
      if( mlink > 0 )
         return true;
   }
   return false;
}

