/*
 * StochPresolverTinyEntries.C
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#include "StochPresolverTinyEntries.h"


StochPresolverTinyEntries::StochPresolverTinyEntries(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverTinyEntries::~StochPresolverTinyEntries()
{
}


void StochPresolverTinyEntries::applyPresolving()
{
   int nelims = 0;
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifndef NDEBUG
   if( myRank == 0 )
      cout << "Starting model cleaup..." << endl;
   countRowsCols();
#endif

   nelims += removeTinyEntriesSystemA();
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*presProb->A);
   updateTransposed(A);

   nelims += removeTinyEntriesSystemC();
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*presProb->C);
   updateTransposed(C);

   // todo: removal of redundant constraints
   // set current pointers
   // per row, compute activies and compare to lhs and rhs

#ifndef NDEBUG
   if( myRank == 0)
      cout << "Model cleanup finished. Removed "<< nelims <<" entries in total." << std::endl;
#endif
}

int StochPresolverTinyEntries::removeTinyEntriesSystemA()
{
   int nelims = 0;
   localNelims = 0;

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   presData.resetRedCounters();
   setCurrentPointersToNull();

   int nelimsB0 = 0;
   if( updateCPforTinyEntry( -1, EQUALITY_SYSTEM) )
   {
      nelimsB0 = removeTinyInnerLoop( -1, EQUALITY_SYSTEM, LINKING_VARS_BLOCK );
      assert( nelimsB0 == localNelims );
      updateNnzUsingReductions( (*presData.nRowElemsA).vec, currRedRow);
   }

   if( myRank == 0 )
      nelims += nelimsB0;

   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

   assert( matrix.children.size() == presData.nRowElemsA->children.size() );
   assert( matrix.children.size() == presData.redCol->children.size() );

   // go through the children
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( updateCPforTinyEntry(int(it), EQUALITY_SYSTEM) )
      {
         nelims += removeTinyChild((int)it, EQUALITY_SYSTEM);

         updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nRowElemsA).children[it])->vec, currRedRow);
         updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nColElems).children[it])->vec, currRedColChild);

         assert( dynamic_cast<StochVector*>((*presData.nRowElemsA).children[it])->vecl == NULL );
      }
   }

   // update nColElems.vec via AllReduce
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec)->elements();
      const int message_size = dynamic_cast<SimpleVector*>(presData.redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions((*presData.nColElems).vec, presData.redCol->vec);

   // the linking rows:
   if( hasLinking(EQUALITY_SYSTEM) )
   {
      updateCPforLinkingRows( -1, EQUALITY_SYSTEM );
      int nelimsF0 = removeTinyLinkingRows(-1, EQUALITY_SYSTEM );
      // update nCol and nRowA using reductions:
      updateNnzUsingReductions( (*presData.nColElems).vec, currRedColParent);
      updateNnzUsingReductions( (*presData.nRowElemsA).vecl, currRedRow);

      if( myRank == 0 )
         nelims += nelimsF0;

      // go through the children:
      for( size_t it = 0; it< matrix.children.size(); it++)
      {
         if( updateCPforLinkingRows(int(it), EQUALITY_SYSTEM) )
         {
            nelims += removeTinyLinkingRows((int)it, EQUALITY_SYSTEM);
            // update nCol using reductions:
            updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nColElems).children[it])->vec, currRedColChild);
         }
      }

      // update nRowElems.vecl via AllReduce:
      if( iAmDistrib )
      {
         double* redRowLink = dynamic_cast<SimpleVector*>(presData.redRowA->vecl)->elements();
         const int message_size = dynamic_cast<SimpleVector*>(presData.redRowA->vecl)->length();
         MPI_Allreduce(MPI_IN_PLACE, redRowLink, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
   }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nelims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   return nelims;
}

int StochPresolverTinyEntries::removeTinyEntriesSystemC()
{
   int nelims = 0;
   localNelims = 0;

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   presData.resetRedCounters();
   setCurrentPointersToNull();

   int nelimsB0 = 0;
   if( updateCPforTinyEntry( -1, INEQUALITY_SYSTEM) )
   {
      nelimsB0 = removeTinyInnerLoop( -1, INEQUALITY_SYSTEM, LINKING_VARS_BLOCK );
      assert( nelimsB0 == localNelims );
      updateNnzUsingReductions( (*presData.nRowElemsC).vec, currRedRow);
   }

   if( myRank == 0 )
      nelims += nelimsB0;

   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   assert( matrix.children.size() == presData.nRowElemsC->children.size() );
   assert( matrix.children.size() == presData.redCol->children.size() );

   // go through the children
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( updateCPforTinyEntry((int)it, INEQUALITY_SYSTEM) )
      {
         nelims += removeTinyChild((int)it, INEQUALITY_SYSTEM);

         updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nRowElemsC).children[it])->vec, currRedRow);
         updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nColElems).children[it])->vec, currRedColChild);

         assert( dynamic_cast<StochVector*>((*presData.nRowElemsC).children[it])->vecl == NULL );
      }
   }

   // update nColElems.vec via AllReduce
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec)->elements();
      const int message_size = dynamic_cast<SimpleVector*>(presData.redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions((*presData.nColElems).vec, presData.redCol->vec);

   // the linking rows:
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      updateCPforLinkingRows( -1, INEQUALITY_SYSTEM );
      int nelimsF0 = removeTinyLinkingRows(-1, INEQUALITY_SYSTEM );
      // update nCol and nRowA using reductions:
      updateNnzUsingReductions( (*presData.nColElems).vec, currRedColParent);
      updateNnzUsingReductions( (*presData.nRowElemsC).vecl, currRedRow);

      if( myRank == 0 )
         nelims += nelimsF0;

      // go through the children:
      for( size_t it = 0; it< matrix.children.size(); it++)
      {
         if( updateCPforLinkingRows(int(it), INEQUALITY_SYSTEM) )
         {
            nelims += removeTinyLinkingRows((int)it, INEQUALITY_SYSTEM);
            // update nCol using reductions:
            updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nColElems).children[it])->vec, currRedColChild);
         }
      }

      // update nRowElems.vecl via AllReduce:
      if( iAmDistrib )
      {
         double* redRowLink = dynamic_cast<SimpleVector*>(presData.redRowC->vecl)->elements();
         const int message_size = dynamic_cast<SimpleVector*>(presData.redRowC->vecl)->length();
         MPI_Allreduce(MPI_IN_PLACE, redRowLink, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
   }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nelims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   return nelims;
}

/** Calls removeTinyInnerLoop for a Child on the matrices Amat, Bmat which also adapts the rhs. */
int StochPresolverTinyEntries::removeTinyChild( int it, SystemType system_type )
{
   // for Amat:
   int nelims = removeTinyInnerLoop( it, system_type, LINKING_VARS_BLOCK );

   // for Bmat:
   nelims += removeTinyInnerLoop( it, system_type, CHILD_BLOCK );

   // todo: special treatment for the linking rows

   return nelims;
}

/** Removes tiny entries in storage and adapts the rhs accordingly.
 *  If block_type == LINKING_VARS_BLOCK, then block Amat is considered.
 *  If block_type == CHILD_BLOCK, then block Bmat is considered. */
int StochPresolverTinyEntries::removeTinyInnerLoop( int it, SystemType system_type, BlockType block_type )
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
      int end = storage->rowptr[r].end;
      int k = storage->rowptr[r].start;
      while( k < end )
      {
         const int col = storage->jcolM[k];
         if( fabs(storage->M[k]) < tolerance3 )
         {
            // cout << "Remove entry M ( "<< r << ", " << storage->jcolM[k] << " ) = "<<storage->M[k]<<" (by first test)"<<endl;
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

            // cout << "Remove entry M ( "<< r << ", " << col << " ) = "<<storage->M[k]<<" (by second test)"<<endl;
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

int StochPresolverTinyEntries::removeTinyLinkingRows( int it, SystemType system_type )
{
   // Setting all the pointers correctly
   int nelims = 0;
   double tmp = 0;
   double* xlowElems;
   double* xuppElems;
   double* ixlowElems;
   double* ixuppElems;
   double* redCol;

   if( it == -1 )
   {
      xlowElems = currxlowParent->elements();
      xuppElems = currxuppParent->elements();
      ixlowElems = currIxlowParent->elements();
      ixuppElems = currIxuppParent->elements();
      redCol = currRedColParent->elements();
   }
   else
   {
      xlowElems = currxlowChild->elements();
      xuppElems = currxuppChild->elements();
      ixlowElems = currIxlowChild->elements();
      ixuppElems = currIxuppChild->elements();
      redCol = currRedColChild->elements();
   }

   double* rhsElems;
   double* icuppElems;
   double* iclowElems;
   double* cuppElems;
   double* clowElems;

   if( system_type == EQUALITY_SYSTEM )
      rhsElems = currEqRhsLink->elements();
   else
   {
      assert( system_type == INEQUALITY_SYSTEM );
      icuppElems = currIcuppLink->elements();
      iclowElems = currIclowLink->elements();
      cuppElems = currIneqRhsLink->elements();
      clowElems = currIneqLhsLink->elements();
   }

   // the actual work starts here:
   for( int r = 0; r < currBlmat->m; r++ )   // r: row index
   {
      int end = currBlmat->rowptr[r].end;
      int k = currBlmat->rowptr[r].start;
      while( k < end )
      {
         const int col = currBlmat->jcolM[k];
         if( fabs(currBlmat->M[k]) < tolerance3 )
         {
            // cout << "Remove entry M ( "<< r << ", " <<col<<" (by first test) in Linking Constraint."<<endl;
            updateAndSwap(currBlmat, r, k, end, redCol, nelims);
            removeEntryInDynamicStorage(*currBlmatTrans, col, r, tmp);
         }
         else if( fabs(currBlmat->M[k]) < tolerance1 && ixuppElems[col] != 0.0 && ixlowElems[col] != 0.0
               && fabs(currBlmat->M[k]) * (xuppElems[col] - xlowElems[col]) * currNnzRow->elements()[r] < tolerance2 * feastol )
         {
            if( system_type == EQUALITY_SYSTEM )
               rhsElems[r] -= currBlmat->M[k] * xlowElems[col];
            else
            {
               if( icuppElems[r] != 0.0 )
                  cuppElems[r] -= currBlmat->M[k] * xlowElems[col];
               if( iclowElems[r] != 0.0 )
                  clowElems[r] -= currBlmat->M[k] * xlowElems[col];
            }

            // cout << "Remove entry M ( "<< r << ", " << col << " ) (by second test) in Linking Constraint."<<endl;
            updateAndSwap(currBlmat, r, k, end, redCol, nelims);
            removeEntryInDynamicStorage(*currBlmatTrans, col, r, tmp);
         }
         k++;
      }
   }

   return nelims;
}

/**
 * Set the current pointers to the currently necessary data.
 * If it==-1, case as root.
 */
bool StochPresolverTinyEntries::updateCPforTinyEntry(int it, SystemType system_type)
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

/**
 * Set the current pointers to the linking constraints.
 * If it==-1, case at F_0.
 */
bool StochPresolverTinyEntries::updateCPforLinkingRows(int it, SystemType system_type)
{
   if( it == -1 )
      assert( hasLinking(system_type) );

   setCurrentPointersToNull();
   // set the column pointers:
   if( it == -1 )
      setCPColumnRoot();
   else
      setCPColumnChild(it);

   // set the matrix and row pointers:
   if( system_type == EQUALITY_SYSTEM )
   {
      // dummy child? set currBlmat as F_0:
      if( !setCPLinkConstraint(presProb->A, it, system_type) ) return false;
      // set row pointers:
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      setCPRowLinkEquality();
      currRedRow = currRedRowLink;
   }
   else  // system_type == INEQUALITY_SYSTEM
   {
      // dummy child? set currBlmat as G_0:
      if( !setCPLinkConstraint(presProb->C, it, system_type) ) return false;
      // set row pointers:
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      setCPRhsLinkInequality();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      currRedRow = dynamic_cast<SimpleVector*>(presData.redRowC->vecl);
   }

   return true;
}

