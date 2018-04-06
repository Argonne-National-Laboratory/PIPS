/*
 * StochPresolverBase.cpp
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#include "StochPresolverBase.h"



StochPresolverBase::StochPresolverBase(PresolveData& presData)
: presData(presData)
{


   presProb = presData.presProb;

   setCurrentPointersToNull();

   if(presData.redRowA->vecl != NULL)
      currEqRhsAdaptionsLink = new double[redRowA->vecl->n];
   else
      currEqRhsAdaptionsLink = NULL;
   if(redRowC->vecl != NULL)
   {
      currInEqRhsAdaptionsLink = new double[redRowC->vecl->n];
      currInEqLhsAdaptionsLink = new double[redRowC->vecl->n];
   }
   else
   {
      currInEqRhsAdaptionsLink = NULL;
      currInEqLhsAdaptionsLink = NULL;
   }

   localNelims = 0;
   nChildren = nColElems->children.size();
   linkVarsBlocks = new BLOCKS[nChildren + 1];
   childBlocks = new BLOCKS[nChildren + 1];
   resetLinkvarsAndChildBlocks();

   blocks = new int[nChildren + 3];
   blocksIneq = new int[nChildren + 3];
   resetBlocks();
}

StochPresolverBase::~StochPresolverBase()
{
   delete[] linkVarsBlocks;
   delete[] childBlocks;
   delete[] currEqRhsAdaptionsLink;
   delete[] currInEqRhsAdaptionsLink;
   delete[] currInEqLhsAdaptionsLink;
}

void StochPresolverBase::updateAndSwap( SparseStorageDynamic* storage, int rowidx, int& indexK, int& rowEnd, double* redCol, int& nelims)
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

void StochPresolverBase::updateRhsNRowLink()
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

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

      resetEqRhsAdaptionsLink();
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
      resetIneqRhsAdaptionsLink();
   }
}

void StochPresolverBase::setCurrentPointersToNull()
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


bool StochPresolverBase::updateCurrentPointers(int it, SystemType system_type)
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

/** Update the nnzVector by subtracting the reductions vector. */
void StochPresolverBase::updateNnzUsingReductions( OoqpVector* nnzVector, OoqpVector* redVector)
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

void StochPresolverBase::updateTransposed(StochGenMatrix& matrix)
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

void StochPresolverBase::updateTransposedSubmatrix(SparseStorageDynamic& transStorage, const int blockStart, const int blockEnd)
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


void StochPresolverBase::updateLinkingVarsBlocks(int& newSREq, int& newSRIneq)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   // apply updated colAdaptParent to the Amat blocks
   for(int i= -1; i<nChildren; i++)
   {
      // dummy child?
      if( updateCurrentPointersForColAdapt( i, EQUALITY_SYSTEM) )
         newSREq += colAdaptLinkVars(i, EQUALITY_SYSTEM);
      if( updateCurrentPointersForColAdapt( i, INEQUALITY_SYSTEM) )
         newSRIneq += colAdaptLinkVars(i, INEQUALITY_SYSTEM);
   }

   // apply updated colAdaptParent to the F0 block (Blmat of parent):
   // dummy child?
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

   // empty the singletonRow list
   for(int i = 0; i<(int)singletonRows.size(); i++)
      assert(singletonRows[i] == -1);
   singletonRows.clear();
   resetBlocks();

   if( iAmDistrib )
   {  // communicate newly found number of singleton rows so that all processes share this knowledge
      int newSR[2] = {newSREq, newSRIneq};
      MPI_Allreduce(MPI_IN_PLACE, newSR, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      newSREq = newSR[0];
      newSRIneq = newSR[1];
   }
}

/** Update the current pointers for the singleton row routine.
 * If it==-1, we are at parent block. Else, et child[it].
 * Return false if child[it] is a dummy child. */
bool StochPresolverBase::updateCurrentPointersForSingletonRow(int it, SystemType system_type)
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



bool StochPresolverBase::updateCPForSingletonRowInequalityBChild( int it )
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

bool StochPresolverBase::procSingletonRowRoot(StochGenMatrix& stochMatrix)
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


/** Return false if it is a dummy child. Else, set the current pointers to Amat, AmatTrans, Bmat, BmatTrans,
 * currNnzRow, currRedRow, currNnzColChild, currRedColChild.
 * currxlowChild, currxuppChild, currIxlowChild, currIxuppChild.
 * And depending on the systemType: either currEqRhs or currIneqRhs, currIneqLhs, currIcupp and currIclow. */
bool StochPresolverBase::updateCurrentPointersForColAdapt(int it, SystemType system_type)
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

bool StochPresolverBase::updateCPforColAdaptF0( SystemType system_type )
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


void StochPresolverBase::resetLinkvarsAndChildBlocks()
{
   for( int i = 0; i < nChildren+1; i++)
   {
      linkVarsBlocks[i].start = 0;
      linkVarsBlocks[i].end = 0;
      childBlocks[i].start = 0;
      childBlocks[i].end = 0;
   }
}


void StochPresolverBase::resetBlocks()
{
   for( int i = 0; i < nChildren+3; i++)
   {
      blocks[i] = 0;
      blocksIneq[i] = 0;
   }
}

void StochPresolverBase::resetRedCounters()
{
   redRowA->setToZero();
   redRowC->setToZero();
   redCol->setToZero();
}

void StochPresolverBase::resetEqRhsAdaptionsLink()
{
   assert(hasLinking(EQUALITY_SYSTEM));
   for( int i = 0; i < redRowA->vecl->n; i++)
      currEqRhsAdaptionsLink[i] = 0.0;
}

void StochPresolverBase::resetIneqRhsAdaptionsLink()
{
   assert(hasLinking(INEQUALITY_SYSTEM));
   for( int i = 0; i < redRowC->vecl->n; i++ )
   {
      currInEqRhsAdaptionsLink[i] = 0.0;
      currInEqLhsAdaptionsLink[i] = 0.0;
   }
}


double StochPresolverBase::removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx)
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

void StochPresolverBase::clearRow(SparseStorageDynamic& storage, const int rowIdx)
{
   storage.rowptr[rowIdx].end = storage.rowptr[rowIdx].start;
}

bool StochPresolverBase::childIsDummy(StochGenMatrix& matrix, int it, SystemType system_type)
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

bool StochPresolverBase::hasLinking(SystemType system_type)
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

void StochPresolverBase::getRankDistributed( MPI_Comm comm, int& myRank, bool& iAmDistrib )
{
   MPI_Comm_rank(comm, &myRank);
   int world_size;
   MPI_Comm_size(comm, &world_size);
   if( world_size > 1) iAmDistrib = true;
   else iAmDistrib = false;
}


/*
 * Given a vector<COLUMNTOADAPT>, this routine goes through all columns inside and removes them
 * from the current Bmat block. Depending on the system_type, rhs (and lhs) are updated.
 */
bool StochPresolverBase::adaptChildBmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type, int& newSR )
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

bool StochPresolverBase::adaptChildBlmat( std::vector<COLUMNTOADAPT> & colAdaptBlock, SystemType system_type)
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
