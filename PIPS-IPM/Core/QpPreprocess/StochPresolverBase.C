/*
 * StochPresolverBase.cpp
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#include "StochPresolverBase.h"
#include <algorithm>


StochPresolverBase::StochPresolverBase(PresolveData& presData)
: presData(presData)
{
   presProb = presData.presProb;

   setCurrentPointersToNull();

   if(presData.redRowA->vecl != NULL)
      currEqRhsAdaptionsLink = new double[presData.redRowA->vecl->n];
   else
      currEqRhsAdaptionsLink = NULL;
   if(presData.redRowC->vecl != NULL)
   {
      currInEqRhsAdaptionsLink = new double[presData.redRowC->vecl->n];
      currInEqLhsAdaptionsLink = new double[presData.redRowC->vecl->n];
   }
   else
   {
      currInEqRhsAdaptionsLink = NULL;
      currInEqLhsAdaptionsLink = NULL;
   }

   localNelims = 0;
   nChildren = presData.getNChildren();
   linkVarsBlocks = new BLOCKS[nChildren + 1];
   childBlocks = new BLOCKS[nChildren + 1];
   resetLinkvarsAndChildBlocks();

   indivObjOffset = 0.0;
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
         double* redRowLink = dynamic_cast<SimpleVector*>(presData.redRowA->vecl)->elements();
         int message_sizeA = dynamic_cast<SimpleVector*>(presData.redRowA->vecl)->length();
         MPI_Allreduce(MPI_IN_PLACE, redRowLink, message_sizeA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, currEqRhsAdaptionsLink, message_sizeA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      // update nRowElemsA.vecl
      updateNnzUsingReductions(presData.nRowElemsA->vecl, presData.redRowA->vecl);
      // update rhs with += adaptionsRhsLink
      for(int i=0; i<currEqRhsLink->n; i++)
         currEqRhsLink->elements()[i] += currEqRhsAdaptionsLink[i];

      resetEqRhsAdaptionsLink();
   }
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      setCPRhsLinkInequality();

      if( iAmDistrib )
      {
         // todo improve mpi communication, only one:
         int message_sizeC = dynamic_cast<SimpleVector*>(presData.redRowC->vecl)->length();
         double* redRowLinkC = dynamic_cast<SimpleVector*>(presData.redRowC->vecl)->elements();
         MPI_Allreduce(MPI_IN_PLACE, redRowLinkC, message_sizeC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, currInEqRhsAdaptionsLink, message_sizeC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, currInEqLhsAdaptionsLink, message_sizeC, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
      // update nRowElemsC.vecl
      updateNnzUsingReductions(presData.nRowElemsC->vecl, presData.redRowC->vecl);
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

void StochPresolverBase::storeRemovedEntryIndex(int rowidx, int colidx, int it, BlockType block_type)
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

/** Should be called right after doSingletonRowsA() or another method that stores
 * information to update in the member variable colAdaptParent.
 * Updates the blocks A,C,F0,G0 using colAdaptParent.
 * Returns the number of newly found singleton rows (equality/inequality system)
 * during adaption of A,C,F0,G0.
 */
void StochPresolverBase::updateLinkingVarsBlocks(int& newSREq, int& newSRIneq)
{
   cout<<"colAdaptParent has "<<presData.getNumberColAdParent()<<" entries."<<endl;
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
   presData.clearColAdaptParent();

   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(presData.redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions(presData.nColElems->vec, presData.redCol->vec);

   presData.resetRedCounters();

   // empty the singletonRow list
   presData.clearSingletonRows();

   if( iAmDistrib )
   {  // communicate newly found number of singleton rows so that all processes share this information
      int newSR[2] = {newSREq, newSRIneq};
      MPI_Allreduce(MPI_IN_PLACE, newSR, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      newSREq = newSR[0];
      newSRIneq = newSR[1];
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
   currNnzColParent = NULL;
}

bool StochPresolverBase::updateCPforTinyEntry(int it, SystemType system_type)
{
   setCurrentPointersToNull();
   setCPColumnRoot();

   if( system_type == EQUALITY_SYSTEM )
   {
      // dummy child? set currAmat and currBmat
      if( !setCPAmatBmat(presProb->A, it, EQUALITY_SYSTEM) ) return false;
      if( it == -1 ) // case at root
         setCPRowRootEquality();
      else  // at child it
         setCPRowChildEquality(it);
   }
   else  // system_type == INEQUALITY_SYSTEM
   {
      assert( system_type == INEQUALITY_SYSTEM );

      // dummy child? set currAmat and currBmat
      if( !setCPAmatBmat(presProb->C, it, INEQUALITY_SYSTEM) ) return false;

      if( it == -1 ) // case at root
         setCPRowRootInequality();
      else  // at child it
         setCPRowChildInequality(it);
   }
   if( it > -1 ) // at child it
      setCPColumnChild(it);

   return true;
}

/** Update the current pointers for the singleton row routine.
 * If it==-1, we are at parent block. Else, et child[it].
 * Return false if child[it] is a dummy child. */
bool StochPresolverBase::updateCPForSingletonRow(int it, SystemType system_type)
{
   setCurrentPointersToNull();

   setCPColumnRoot();
   currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);

   if( it == -1 )
   {
      if( system_type == EQUALITY_SYSTEM )
      {
         setCPAmatsRoot(presProb->C);
         setCPRowRootEquality();
      }
      else  // INEQUALITY_SYSTEM
      {
         assert( system_type == INEQUALITY_SYSTEM );
         setCPAmatsRoot(presProb->C);
         setCPRowRootInequality();
         currNnzColParent = dynamic_cast<SimpleVector*>(presData.nColElems->vec);
      }
   }
   else  // at child it
   {
      if( system_type == EQUALITY_SYSTEM )
      {
         // child is dummy? set currAmat, AmatTrans, Bmat, BmatTrans
         if( !setCPAmatsChild( presProb->A,  it, system_type)) return false;
         if( !setCPBmatsChild( presProb->A,  it, system_type)) return false;
         setCPRowChildEquality(it);

         if( hasLinking(system_type) )
         {
            setCPBlmatsChild( presProb->A, it);
            setCPRowLinkEquality();
         }
      }
      else  // INEQUALITY_SYSTEM
      {
         // child is dummy? set currAmat, AmatTrans, Bmat, BmatTrans
         if( !setCPAmatsChild( presProb->C,  it, system_type)) return false;
         if( !setCPBmatsChild( presProb->C,  it, system_type)) return false;
         setCPRowChildInequality(it);

         if( hasLinking(system_type) )
         {
            setCPBlmatsChild( presProb->C, it);
            setCPRhsLinkInequality();
            currRedRowLink = dynamic_cast<SimpleVector*>(presData.redRowC->vecl);
         }
      }
      setCPColumnChild(it);
      currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec);
      currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
   }
   return true;
}

/** Return false if it is a dummy child. */
bool StochPresolverBase::updateCPForSingletonRowInequalityBChild( int it )
{
   assert( it >= 0);
   setCurrentPointersToNull();

   // dummy child? set currBmat, currBmatTrans
   if( !setCPBmatsChild(presProb->C, it, INEQUALITY_SYSTEM)) return false;
   setCPRowChildInequality(it);

   currRedColChild = dynamic_cast<SimpleVector*>(presData.redCol->children[it]->vec);
   currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);

   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      setCPBlmatsChild(presProb->C, it);
      currIcuppLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vecl);
      currIclowLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vecl);
      currRedRowLink = dynamic_cast<SimpleVector*>(presData.redRowC->vecl);
   }
   return true;
}

/** Return false if it is a dummy child. */
bool StochPresolverBase::updateCPForSingletonRowEqualityBChild( int it )
{
   assert( it >= 0);
   setCurrentPointersToNull();

   // dummy child? set currBmat, currBmatTrans
   if( !setCPBmatsChild(presProb->A, it, EQUALITY_SYSTEM)) return false;
   setCPRowChildEquality(it);

   currRedColChild = dynamic_cast<SimpleVector*>(presData.redCol->children[it]->vec);
   currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);

   if( hasLinking(EQUALITY_SYSTEM) )
   {
      setCPBlmatsChild(presProb->A, it);
      currIcuppLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vecl);
      currIclowLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vecl);
      currRedRowLink = dynamic_cast<SimpleVector*>(presData.redRowA->vecl);
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
         setCPAmatsRoot(presProb->A);
         setCPRowRootEquality();
      }
      else  // system_type == INEQUALITY_SYSTEM
      {
         setCPAmatsRoot(presProb->C);
         setCPRowRootInequality();
      }
      currRedColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec);

   }
   else  // at child it
   {
      if( system_type == EQUALITY_SYSTEM )
      {
         if( !setCPAmatsChild( presProb->A, it, EQUALITY_SYSTEM)) return false;
         setCPRowChildEquality(it);
      }
      else  // system_type == INEQUALITY_SYSTEM
      {
         assert( system_type == INEQUALITY_SYSTEM );

         if( !setCPAmatsChild( presProb->C, it, INEQUALITY_SYSTEM)) return false;
         setCPRowChildInequality(it);
      }
      currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
      currRedColChild = dynamic_cast<SimpleVector*>(presData.redCol->children[it]->vec);
      currRedColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec);
   }
   return true;
}

/** Return false if it is a dummy child. */
bool StochPresolverBase::updateCPforColAdaptF0( SystemType system_type )
{
   setCurrentPointersToNull();
   if( !hasLinking(system_type) )
      return false;

   if( system_type == EQUALITY_SYSTEM )
   {
      setCPBlmatsRoot(presProb->A);
      setCPRowLinkEquality();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
   }
   else
   {
      setCPBlmatsRoot(presProb->C);
      setCPRhsLinkInequality();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      currRedRowLink = dynamic_cast<SimpleVector*>(presData.redRowC->vecl);
   }
   return true;
}

/** Set currAmat = root.Bmat */
void StochPresolverBase::setCPAmatsRoot(GenMatrixHandle matrixHandle)
{
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   currAmat = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamic();
   currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.Bmat)->getStorageDynamicTransposed();
}

bool StochPresolverBase::setCPAmatsChild(GenMatrixHandle matrixHandle, int it, SystemType system_type)
{
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   if( childIsDummy(matrix, it, system_type) )
      return false;
   currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
   currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamicTransposed();
   return true;
}

bool StochPresolverBase::setCPBmatsChild(GenMatrixHandle matrixHandle, int it, SystemType system_type)
{
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   if( childIsDummy(matrix, it, system_type) )
      return false;
   currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
   currBmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamicTransposed();
   return true;
}

/** Set the current pointers for currAmat and currBmat at the root.
 * Return false if it is a dummy child. */
bool StochPresolverBase::setCPAmatBmat(GenMatrixHandle matrixHandle, int it, SystemType system_type)
{
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   if( it == -1 ) // at root
   {
      // save Bmat as currAmat for easy computation in TinyInnerLoop
      currAmat = matrix.Bmat->getStorageDynamic();
   }
   else  // at child it
   {
      if( childIsDummy(matrix, it, system_type) ) return false;
      currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
      currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
   }
   return true;
}

void StochPresolverBase::setCPBlmatsRoot(GenMatrixHandle matrixHandle)
{
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamic();
   currBlmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.Blmat)->getStorageDynamicTransposed();
}

void StochPresolverBase::setCPBlmatsChild(GenMatrixHandle matrixHandle, int it)
{
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamic();
   currBlmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamicTransposed();
}

void StochPresolverBase::setCPColumnRoot()
{
   currRedColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec);
   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);
   currIxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).vec);
   currIxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).vec);
}

void StochPresolverBase::setCPColumnChild(int it)
{
   currRedColChild = dynamic_cast<SimpleVector*>(presData.redCol->children[it]->vec);
   currxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->vec);
   currxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->vec);
   currIxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[it]->vec);
   currIxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[it]->vec);
}

void StochPresolverBase::setCPRowRootEquality()
{
   currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowA->vec);
}

void StochPresolverBase::setCPRowRootInequality()
{
   currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec);
   currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vec);
   currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vec);
   currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vec);
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowC->vec);
}

void StochPresolverBase::setCPRowChildEquality(int it)
{
   currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowA->children[it]->vec);
}

void StochPresolverBase::setCPRowChildInequality(int it)
{
   currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec);
   currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec);
   currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec);
   currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec);
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowC->children[it]->vec);
}

void StochPresolverBase::setCPRowLinkEquality()
{
   currEqRhsLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vecl);
   currRedRowLink = dynamic_cast<SimpleVector*>(presData.redRowA->vecl);
}

void StochPresolverBase::setCPRhsLinkInequality()
{
   currIneqRhsLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vecl);
   currIneqLhsLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vecl);
   currIcuppLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vecl);
   currIclowLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vecl);
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

void StochPresolverBase::resetEqRhsAdaptionsLink()
{
   assert(hasLinking(EQUALITY_SYSTEM));
   for( int i = 0; i < presData.redRowA->vecl->n; i++)
      currEqRhsAdaptionsLink[i] = 0.0;
}

void StochPresolverBase::resetIneqRhsAdaptionsLink()
{
   assert(hasLinking(INEQUALITY_SYSTEM));
   for( int i = 0; i < presData.redRowC->vecl->n; i++ )
   {
      currInEqRhsAdaptionsLink[i] = 0.0;
      currInEqLhsAdaptionsLink[i] = 0.0;
   }
}

/** Removes the specified entry from storage and stores its value in m.
 * Returns false if the specified entry does not exist anymore in storage.
 * For example, if the entry was removed before because of redundancy.
 */
bool StochPresolverBase::removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx, double& m)
{
   int i = -1;
   int end = storage.rowptr[rowIdx].end;
   int start = storage.rowptr[rowIdx].start;

   for( i=start; i<end; i++)
   {
      if( storage.jcolM[i] == colIdx )
         break;
   }
   if( i < 0 || i == end )
      return false;
   m = storage.M[i];
   std::swap(storage.M[i],storage.M[end-1]);
   std::swap(storage.jcolM[i],storage.jcolM[end-1]);
   storage.rowptr[rowIdx].end --;

   return true;
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
      assert( presData.redCol->children[it]->isKindOf(kStochDummy) );

      if( system_type == EQUALITY_SYSTEM)
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->isKindOf(kStochDummy) );
         assert( presData.nRowElemsA->children[it]->isKindOf(kStochDummy) );
         assert( presData.redRowA->children[it]->isKindOf(kStochDummy) );
      }
      else
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->isKindOf(kStochDummy) );
         assert( presData.nRowElemsC->children[it]->isKindOf(kStochDummy) );
         assert( presData.redRowC->children[it]->isKindOf(kStochDummy) );
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
      {
         // todo: assert that all vectors and matrices have linking part
         assert(presData.redRowA->vecl);
         return true;
      }
   }
   else
   {
      dynamic_cast<StochGenMatrix&>(*(presProb->C)).Blmat->getSize(mlink, nlink);
      if( mlink > 0 )
      {
         // todo: assert that all vectors and matrices have linking part
         assert(presData.redRowC->vecl);
         return true;
      }
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
         double m = 0.0;
         bool entryExists = removeEntryInDynamicStorage(*currBmat, rowIdx, colIdx, m);
         if( !entryExists )
            continue;
         if( system_type == EQUALITY_SYSTEM )
            currEqRhs->elements()[rowIdx] -= m * val;
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
   cout<<"colAdaptBlock.size(): "<<colAdaptBlock.size()<<endl;
   assert(currBlmat != NULL);
   if( system_type == EQUALITY_SYSTEM )
   {
      assert( presData.redRowA->vecl->n == currRedRowLink->n );
      assert( presData.redRowA->vecl->n == currBlmat->m );
   }
   else
   {
      assert(system_type == INEQUALITY_SYSTEM);
      assert( presData.redRowC->vecl->n == currRedRowLink->n );
      assert( presData.redRowC->vecl->n == currBlmat->m );
   }

   for(int i=0; i<(int)colAdaptBlock.size(); i++)
   {
      int colIdx = colAdaptBlock[i].colIdx;
      double val = colAdaptBlock[i].val;

      for( int j = currBlmatTrans->rowptr[colIdx].start; j<currBlmatTrans->rowptr[colIdx].end; j++ )
      {
         int rowIdx = currBlmatTrans->jcolM[j];
         double m = 0.0;
         bool entryExists = removeEntryInDynamicStorage(*currBlmat, rowIdx, colIdx, m);
         if( !entryExists )
            continue;
         cout<<"Removed entry "<<rowIdx<<", "<<colIdx<<" with value "<<m<<" in F_i, system_type:"<<system_type<<endl;

         if( system_type == EQUALITY_SYSTEM )
            currEqRhsAdaptionsLink[rowIdx] -= m * val;
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

/** For the given column index and the value to which this variable is to be fixed,
 * go through the current Bmat by going through the corresponding row of the transposed Bmat.
 * Adapt the rhs, the nnzRow, the nnzCol and all concerned entries in Bmat and BmatTransposed.
 * Return the number of newly found singleton rows.
 */
int StochPresolverBase::adaptChildBmatCol(int colIdx, double val, SystemType system_type)
{
   int newSingletonRows = 0;

   for( int j = currBmatTrans->rowptr[colIdx].start; j<currBmatTrans->rowptr[colIdx].end; j++ )
   {
      int rowIdxB = currBmatTrans->jcolM[j];
      double m = 0.0;
      bool entryExists = removeEntryInDynamicStorage(*currBmat, rowIdxB, colIdx, m);
      if( !entryExists )
         continue;

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
bool StochPresolverBase::adaptOtherSystemChildB( SystemType system_type, std::vector<COLUMNTOADAPT> & colAdaptBblock, int& newSR )
{
   // Bmat blocks
   bool possFeas = adaptChildBmat( colAdaptBblock, system_type, newSR);
   if( !possFeas ) return false;

   // Blmat blocks
   if( hasLinking(system_type) )
   {
      bool possFeas = adaptChildBlmat( colAdaptBblock, system_type);
      if( !possFeas ) return false;
   }

   return true;
}

/** Adapt the columns for the linking-variable-blocks (the A_i) blocks */
int StochPresolverBase::colAdaptLinkVars(int it, SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int newSingletonRows = 0;

   for( int i=0; i < presData.getNumberColAdParent(); i++)
   {
      int colIdxA = presData.getColAdaptParent(i).colIdx;
      double val = presData.getColAdaptParent(i).val;

      for( int j = currAmatTrans->rowptr[colIdxA].start; j < currAmatTrans->rowptr[colIdxA].end; j++ )
      {
         int rowIdxA = currAmatTrans->jcolM[j];
         double m = 0.0;
         bool entryExists = removeEntryInDynamicStorage(*currAmat, rowIdxA, colIdxA, m);
         if( !entryExists )
            continue;
         cout<<"Removed entry "<<rowIdxA<<", "<<colIdxA<<" with value "<<m<<" in Amat of child "<<it<<" system_type:"<<system_type<<endl;

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

int StochPresolverBase::colAdaptF0(SystemType system_type)
{
   assert(currBlmat != NULL);
   assert( currNnzRow->n == currBlmat->m );

   int newSingletonRows = 0;

   for(int i=0; i<presData.getNumberColAdParent(); i++)
   {
      int colIdx = presData.getColAdaptParent(i).colIdx;
      double val = presData.getColAdaptParent(i).val;

      for( int j = currBlmatTrans->rowptr[colIdx].start; j<currBlmatTrans->rowptr[colIdx].end; j++ )
      {
         int rowIdx = currBlmatTrans->jcolM[j];
         double m = 0.0;
         bool entryExists = removeEntryInDynamicStorage(*currBlmat, rowIdx, colIdx, m);
         if( !entryExists )
            continue;
         cout<<"Removed entry "<<rowIdx<<", "<<colIdx<<" with value "<<m<<" in F0("<<system_type<<")"<<endl;

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
         dynamic_cast<SimpleVector*>(presData.nColElems->vec)->elements()[colIdx] --;
         currNnzRow->elements()[rowIdx] --;

         assert( currNnzRow->elements()[rowIdx] >= 0.0 );
         assert( dynamic_cast<SimpleVector*>(presData.nColElems->vec)->elements()[colIdx] >= 0.0 );

         if(currNnzRow->elements()[rowIdx] == 1.0)
            newSingletonRows++;
      }
      clearRow(*currBlmatTrans, colIdx);
   }
   return newSingletonRows;
}

/** Sum up the individual objective offset on all processes. */
void StochPresolverBase::sumIndivObjOffset()
{
   int myRank;
   bool iAmDistrib = false;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int world_size;
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   if( world_size > 1) iAmDistrib = true;

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &indivObjOffset, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}

