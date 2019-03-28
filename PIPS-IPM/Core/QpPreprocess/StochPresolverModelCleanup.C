/*
 * StochPresolverModelCleanup.C
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#include "StochPresolverModelCleanup.h"
#include <cmath>

StochPresolverModelCleanup::StochPresolverModelCleanup(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverModelCleanup::~StochPresolverModelCleanup()
{
}


void StochPresolverModelCleanup::applyPresolving()
{
   setCurrentPointersToNull();

   int n_elims = 0;

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "Starting model cleanup..." << std::endl;
   // todo this a non const function! watch out
   // countRowsCols();
#endif



   /* remove entries from A and C matrices */
   n_elims += removeTinyEntriesFromSystem(EQUALITY_SYSTEM);
   updateTransposed( dynamic_cast<StochGenMatrix&>(*presProb->A) );

   n_elims += removeTinyEntriesFromSystem(INEQUALITY_SYSTEM);
   updateTransposed( dynamic_cast<StochGenMatrix&>(*presProb->C) );

   // removal of redundant constraints
   int nRemovedRows = 0;
   nRemovedRows += removeRedundantRows(presProb->A, EQUALITY_SYSTEM);
   nRemovedRows += removeRedundantRows(presProb->C, INEQUALITY_SYSTEM);

   // update the nnzColParent counters:
   updateNnzColParent(MPI_COMM_WORLD);
   presData.resetRedCounters();

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nRemovedRows, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#ifndef NDEBUG
   if( myRank == 0)
      std::cout << "Model cleanup finished. Removed " << nRemovedRows << " redundant rows and " << n_elims << " entries in total." << std::endl;
   // countRowsCols();
#endif
}

/** Remove redundant rows in the constraint system. Compares the minimal and maximal row activity
 * with the row bounds (lhs and rhs). If a row is found to be redundant, it is removed.
 * If infeasiblity is detected, then Abort.
 */
int StochPresolverModelCleanup::removeRedundantRows(GenMatrixHandle matrixHandle, SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int nRemovedRows = 0;

   // set current pointers: currAmat, currBmat, lhr, rhs
   // root:
   setCurrentPointersToNull();
   setCPAmatsRoot(matrixHandle);
   setCPColumnRoot();
   if( system_type == EQUALITY_SYSTEM )
      setCPRowRootEquality();
   else
      setCPRowRootInequality();

   int nRemovedRowsRoot = removeRedundantRowsBlockwise(system_type, true);
   // update nnzColParent counters using redColParent:
   updateNnzUsingReductions(presData.nColElems->vec, presData.redCol->vec);
   presData.resetRedCounters();

   if( myRank==0 )
      nRemovedRows += nRemovedRowsRoot;

   // children:
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(matrixHandle));
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( setCPAmatsChild(matrixHandle, (int)it, system_type)) // Amat, AmatTrans
      {
         setCPBmatsChild(matrixHandle, (int)it, system_type);  // Bmat, BmatTrans
         setCPColumnRoot();                                    // redColParent, xlowParent etc
         setCPColumnChild((int)it);
         currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
         if( system_type == EQUALITY_SYSTEM )
            setCPRowChildEquality((int)it);                    // nnzRow, EqRhs
         else
            setCPRowChildInequality((int)it);                  // nnzRow, IneqLhs, IneqRhs

         nRemovedRows += removeRedundantRowsBlockwise(system_type, false);
      }
   }
   // linking rows:
   if( hasLinking(system_type))
   {
      int nRemovedRowsLink = removeRedundantLinkingRows(matrixHandle, system_type);
      if( myRank==0 )
      {
         nRemovedRows += nRemovedRowsLink;
         cout<<"Removed Linking Constraints during model cleanup: "<< nRemovedRowsLink <<endl;
      }
   }

   return nRemovedRows;
}

/** Find and remove redundant rows per child (!atRoot) or for the root block. */
int StochPresolverModelCleanup::removeRedundantRowsBlockwise(SystemType system_type, bool atRoot)
{
   assert( currAmat && currAmatTrans );
   assert( currxlowParent && currIxlowParent && currxuppParent && currIxuppParent );
   assert( currRedColParent && currNnzRow );
   if( atRoot )
   {
      assert( NULL==currBmat && NULL==currBmatTrans && NULL==currNnzColChild );
      assert( NULL==currxlowChild && NULL==currIxlowChild && NULL==currxuppChild && NULL==currIxuppChild );
   }
   else
   {
      assert( currBmat && currBmatTrans );
      assert( currxlowChild && currIxlowChild && currxuppChild && currIxuppChild );
      assert( currNnzColChild );
   }
   if( system_type == EQUALITY_SYSTEM )
      assert( currEqRhs );
   else
      assert( currIclow && currIcupp && currIneqLhs && currIneqRhs );

   int nRemovedRows = 0;

   // per row, compute activies and compare to lhs and rhs:
   for( int r=0; r<currAmat->m; r++)
   {
      if( currNnzRow->elements()[r] == 0.0 ) // empty rows might still have a rhs but should be ignored.
         continue;

      double minAct = 0.0;
      double maxAct = 0.0;
      computeActivityBlockwise(*currAmat, r, -1, minAct, maxAct, *currxlowParent, *currIxlowParent, *currxuppParent, *currIxuppParent);
      if( !atRoot )
         computeActivityBlockwise(*currBmat, r, -1, minAct, maxAct, *currxlowChild, *currIxlowChild, *currxuppChild, *currIxuppChild);

      if( system_type == EQUALITY_SYSTEM )
      {
         if( ( minAct > currEqRhs->elements()[r] + feastol ) || ( maxAct < currEqRhs->elements()[r] - feastol ) )
            abortInfeasible(MPI_COMM_WORLD);
         else if( (minAct >= currEqRhs->elements()[r] - feastol) && (maxAct <= currEqRhs->elements()[r] + feastol) )
         {
            // discard row r
            removeRow(r, *currAmat, *currAmatTrans, currBmat, currBmatTrans, *currNnzRow, *currRedColParent, currNnzColChild);
            nRemovedRows++;
         }
      }
      else  // system_type == INEQUALITY_SYSTEM
      {
         // reformulation of the conditions:
         if( ( currIclow->elements()[r] != 0.0 && maxAct < currIneqLhs->elements()[r] - feastol )
               || ( currIcupp->elements()[r] != 0.0 && minAct > currIneqRhs->elements()[r] + feastol ) )
            abortInfeasible(MPI_COMM_WORLD);
         if( (currIclow->elements()[r] == 0.0 || (currIclow->elements()[r] != 0.0 && currIneqLhs->elements()[r] <= -infinity)) &&
               (currIcupp->elements()[r] == 0.0 || (currIcupp->elements()[r] != 0.0 && currIneqRhs->elements()[r] >= infinity)) )
         {  // discard row r
            removeRow(r, *currAmat, *currAmatTrans, currBmat, currBmatTrans, *currNnzRow, *currRedColParent, currNnzColChild);
            nRemovedRows++;
         }
         else if( (currIclow->elements()[r] == 0.0 || (currIclow->elements()[r] != 0.0 && minAct >= currIneqLhs->elements()[r] - feastol))
               && (currIcupp->elements()[r] == 0.0 || (currIcupp->elements()[r] != 0.0 && maxAct <= currIneqRhs->elements()[r] + feastol)) )
              // ( maxAct <= currIneqRhs->elements()[r] + feastol && minAct >= currIneqLhs->elements()[r] - feastol ) )
         {  // discard row r
            removeRow(r, *currAmat, *currAmatTrans, currBmat, currBmatTrans, *currNnzRow, *currRedColParent, currNnzColChild);
            nRemovedRows++;
         }
      }
   }
   return nRemovedRows;
}

/** Find and remove redundant linking constraints. */
int StochPresolverModelCleanup::removeRedundantLinkingRows(GenMatrixHandle matrixHandle, SystemType system_type)
{
   assert( hasLinking(system_type) );

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   int nRemovedRows = 0;

   int nLinkRows = 0;
   if( system_type== EQUALITY_SYSTEM)
      nLinkRows = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl)->n;
   else
      nLinkRows = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl)->n;

   // todo: use other structure supporting bound checks?
   double* minActivity = new double[nLinkRows];
   double* maxActivity = new double[nLinkRows];
   // bool array of length nLinkRows indicating if row is redundant:
   bool* rowIsRedundant = new bool[nLinkRows];
   for( int i=0; i<nLinkRows; i++)
   {
      minActivity[i] = 0.0;
      maxActivity[i] = 0.0;
      rowIsRedundant[i] = false;
   }

   computeLinkingRowActivity(matrixHandle, system_type, minActivity, maxActivity, nLinkRows);

   // set pointers for checkRedundantLinkingRow():
   if( system_type == EQUALITY_SYSTEM )
   {
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      currEqRhsLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vecl);
   }
   else
   {
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      setCPRhsLinkInequality();
   }
   checkRedundantLinkingRow(matrixHandle, system_type, minActivity, maxActivity, nLinkRows, rowIsRedundant);

   // remove redundant rows indicated in rowIsRedundant:
   // per child, set pointers for F_i blocks
   //    go through all rows and remove the indicated redundant rows
   // set pointers for the F_0 block
   // go again through all rows, remove them in F_0, count number of removed rows and set nnzRow to 0.

   for( int it = 0; it < nChildren; it++)
   {
      if( !nodeIsDummy( it, system_type) )
      {
         setCPBlmatsChild(matrixHandle, (int)it);
         currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);

         for( int r=0; r < nLinkRows; r++)
         {
            if( rowIsRedundant[r] )
               removeRowInBblock(r, currBlmat, currBlmatTrans, currNnzColChild);
         }
      }
   }
   // F_0 block
   setCPBlmatsRoot(matrixHandle);
   currRedColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec);
   if( system_type == EQUALITY_SYSTEM )
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
   else
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);

   for( int r=0; r<nLinkRows; r++)
   {
      if( rowIsRedundant[r] )
      {
         removeRow(r, *currBlmat, *currBlmatTrans, NULL, NULL, *currNnzRow, *currRedColParent, NULL);
         nRemovedRows++;   // only count the removed row once
      }
   }
   // update the nnzColParent counters:
   updateNnzColParent(MPI_COMM_WORLD);
   presData.resetRedCounters();

   delete[] minActivity;
   delete[] maxActivity;
   delete[] rowIsRedundant;

   return nRemovedRows;
}

/** Compute the minimal and maximal row activies of the linking rows. */
void StochPresolverModelCleanup::computeLinkingRowActivity(GenMatrixHandle matrixHandle, SystemType system_type,
      double* minActivity, double* maxActivity, int nLinkRows)
{
   assert( hasLinking(system_type) );

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   // set pointers to F_0 and compute min/max activity of F0 block (only rank==0) :
   if( myRank == 0 )
   {
      setCPBlmatsRoot(matrixHandle);   // set currBlmat, currBlmatTrans
      setCPColumnRoot();               // redColParent, xlowParent etc
      for( int r=0; r<nLinkRows; r++)
      {
         double minAct = 0.0;
         double maxAct = 0.0;
         computeActivityBlockwise(*currBlmat, r, -1, minAct, maxAct, *currxlowParent, *currIxlowParent, *currxuppParent, *currIxuppParent);
         minActivity[r] = minAct;
         maxActivity[r] = maxAct;
      }
   }

   // go through children, set the pointers for F_i blocks. compute activity and sum them up
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(matrixHandle));
   for( size_t it = 0; it< matrix.children.size(); it++)
   {
      if( !nodeIsDummy( (int)it, system_type) )
      {
         setCPBlmatsChild(matrixHandle, (int)it);
         setCPColumnChild((int)it);
         for( int r=0; r<nLinkRows; r++)
         {
            double minAct = 0.0;
            double maxAct = 0.0;
            computeActivityBlockwise(*currBlmat, r, -1, minAct, maxAct, *currxlowChild, *currIxlowChild, *currxuppChild, *currIxuppChild);
            minActivity[r] += minAct;
            maxActivity[r] += maxAct;
         }
      }
   }

   // Allreduce sum over all activites
   if( iAmDistrib )
   {
      MPI_Allreduce(MPI_IN_PLACE, minActivity, nLinkRows, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, maxActivity, nLinkRows, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
}

/** For each linking row, check if it is redundant using the minimal and maximal row activities
 * compared to the lhs and rhs. */
void StochPresolverModelCleanup::checkRedundantLinkingRow(GenMatrixHandle matrixHandle, SystemType system_type,
      double* minActivity, double* maxActivity, int nLinkRows, bool* rowIsRedundant)
{
   // check for redundant rows
   for( int r=0; r<nLinkRows; r++)
   {
      if( currNnzRow->elements()[r] == 0.0 ) // empty rows might still have a rhs but should be ignored.
         continue;
      if( system_type == EQUALITY_SYSTEM )
      {
         if( ( minActivity[r] > currEqRhsLink->elements()[r] + feastol ) || ( maxActivity[r] < currEqRhsLink->elements()[r] - feastol ) )
            abortInfeasible(MPI_COMM_WORLD);
         else if( (minActivity[r] >= currEqRhsLink->elements()[r] - feastol) && (maxActivity[r] <= currEqRhsLink->elements()[r] + feastol) )
            rowIsRedundant[r] = true;  // discard row r
      }
      else  // system_type == INEQUALITY_SYSTEM
      {
         if( (currIclowLink->elements()[r] != 0.0 && maxActivity[r] < currIneqLhsLink->elements()[r] - feastol )
               ||( currIcuppLink->elements()[r] != 0.0 && minActivity[r] > currIneqRhsLink->elements()[r] + feastol ) )
            abortInfeasible(MPI_COMM_WORLD);
         if( (currIclowLink->elements()[r] == 0.0 || (currIclowLink->elements()[r] != 0.0 && currIneqLhsLink->elements()[r] <= -infinity)) &&
               (currIcuppLink->elements()[r] == 0.0 || (currIcuppLink->elements()[r] != 0.0 && currIneqRhsLink->elements()[r] >= infinity)) )
            rowIsRedundant[r] = true;  // discard row r
         else if( (currIclowLink->elements()[r] == 0.0 || (currIclowLink->elements()[r] != 0.0 && minActivity[r] >= currIneqLhsLink->elements()[r] - feastol))
               && (currIcuppLink->elements()[r] == 0.0 || (currIcuppLink->elements()[r] != 0.0 && maxActivity[r] <= currIneqRhsLink->elements()[r] + feastol)) )
            rowIsRedundant[r] = true;  // discard row r
      }
   }
}

int StochPresolverModelCleanup::removeTinyEntriesFromSystem(SystemType system_type)
{
   assert( dynamic_cast<StochGenMatrix&>(*(presProb->A)).children.size() == presData.nRowElemsA->children.size() ); // todo depends on system type
   assert( dynamic_cast<StochGenMatrix&>(*(presProb->A)).children.size() == presData.redCol->children.size() ); // todo depends on system type
   assert( dynamic_cast<StochGenMatrix&>(*(presProb->A)).children.size() == (size_t) nChildren );
   assert( dynamic_cast<StochGenMatrix&>(*(presProb->C)).children.size() == (size_t) nChildren );

   int n_elims = 0;
   localNelims = 0; // todo?

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   presData.resetRedCounters();
   setCurrentPointersToNull();


   /* reductions in root node */
   n_elims += removeEntriesFromRootNode(system_type);
   int nelimsB0 = 0;

   if( myRank == 0 )
      n_elims += nelimsB0;

   // go through the children
   for( int node = 0; node < nChildren; node++)
   {
      // todo if( nodeIsDummy(node, system_type) )
      if( updatePointersForCurrentNode(node, system_type) )
      {
         n_elims += removeTinyChild(node, system_type);

         if(system_type == EQUALITY_SYSTEM)
            updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nRowElemsA).children[node])->vec, currRedRow); // todo system_type wrapper
         else
            updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nRowElemsC).children[node])->vec, currRedRow);

         updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nColElems).children[node])->vec, currRedColChild); // todo system_type wrapper ??

         assert( dynamic_cast<StochVector*>((*presData.nRowElemsA).children[node])->vecl == NULL );
         assert( dynamic_cast<StochVector*>((*presData.nRowElemsC).children[node])->vecl == NULL );
      }
   }

   // update nColElems.vec via AllReduce
   // todo is this needed here? can we do linking first?
   if( iAmDistrib )
   {
      double* redColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec)->elements();
      const int message_size = dynamic_cast<SimpleVector*>(presData.redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   }
   updateNnzUsingReductions((*presData.nColElems).vec, presData.redCol->vec);

   // the linking rows:
   if( hasLinking(system_type) )
   {
      updateCPforLinkingRows( -1, system_type );
      int nelimsF0 = removeTinyLinkingRows(-1, system_type );
      // update nCol and nRowA using reductions:
      updateNnzUsingReductions( (*presData.nColElems).vec, currRedColParent);
      if(system_type == EQUALITY_SYSTEM)
         updateNnzUsingReductions( (*presData.nRowElemsA).vecl, currRedRow); // todo
      else
         updateNnzUsingReductions( (*presData.nRowElemsC).vecl, currRedRow); // todo

      if( myRank == 0 )
         n_elims += nelimsF0;

      // go through the children:
      for( int node = 0; node < nChildren; node++)
      {
         if( updateCPforLinkingRows(node, system_type) )
         {
            n_elims += removeTinyLinkingRows(node, system_type);
            // update nCol using reductions:
            updateNnzUsingReductions( dynamic_cast<StochVector*>((*presData.nColElems).children[node])->vec, currRedColChild);
         }
      }

      // update nRowElems.vecl via AllReduce:
      // todo BUG? update changes of lhs and rhs of linking constraints!!!! remember changes and allreduce
      if( iAmDistrib )
      {
         // TODO wrapper
         if( system_type == EQUALITY_SYSTEM )
         {
            double* redRowLink = dynamic_cast<SimpleVector*>(presData.redRowA->vecl)->elements();
            const int message_size = dynamic_cast<SimpleVector*>(presData.redRowA->vecl)->length();
            MPI_Allreduce(MPI_IN_PLACE, redRowLink, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         }
         else
         {
            double* redRowLink = dynamic_cast<SimpleVector*>(presData.redRowC->vecl)->elements();
            const int message_size = dynamic_cast<SimpleVector*>(presData.redRowC->vecl)->length();
            MPI_Allreduce(MPI_IN_PLACE, redRowLink, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
         }
      }
   }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &n_elims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   return n_elims;
}

int StochPresolverModelCleanup::removeEntriesFromRootNode(SystemType system_type)
{
   int eliminations = 0;

   if( updatePointersForCurrentNode(-1, system_type) )
   {
      eliminations = removeTinyInnerLoop(-1, system_type, LINKING_VARS_BLOCK);

      (system_type == EQUALITY_SYSTEM) ?
            updateNnzUsingReductions((*presData.nRowElemsA).vec, currRedRow) :
            updateNnzUsingReductions((*presData.nRowElemsC).vec, currRedRow);
   }

   assert( eliminations == localNelims);
   return eliminations;
}


/** Calls removeTinyInnerLoop for a Child on the matrices Amat, Bmat which also adapts the rhs. */
int StochPresolverModelCleanup::removeTinyChild( int it, SystemType system_type )
{
   // for Amat:
   int nelims = removeTinyInnerLoop( it, system_type, LINKING_VARS_BLOCK );

   // for Bmat:
   nelims += removeTinyInnerLoop( it, system_type, CHILD_BLOCK );

   return nelims;
}

/** Removes tiny entries in storage and adapts the rhs accordingly.
 *  If block_type == LINKING_VARS_BLOCK, then block Amat is considered.
 *  If block_type == CHILD_BLOCK, then block Bmat is considered. */
/* system type indicates matrix A or C, block_type indicates the block */
int StochPresolverModelCleanup::removeTinyInnerLoop( int it, SystemType system_type, BlockType block_type )
{
   // Setting all the pointers correctly
   int nelims = 0;
   double* nnzRow = currNnzRow->elements();
   SparseStorageDynamic* storage = (block_type == LINKING_VARS_BLOCK) ? currAmat : currBmat;
   double* xlowElems = (block_type == LINKING_VARS_BLOCK) ? currxlowParent->elements() : currxlowChild->elements();
   double* xuppElems = (block_type == LINKING_VARS_BLOCK) ? currxuppParent->elements() : currxuppChild->elements();
   double* ixlowElems = (block_type == LINKING_VARS_BLOCK) ? currIxlowParent->elements() : currIxlowChild->elements();
   double* ixuppElems = (block_type == LINKING_VARS_BLOCK) ? currIxuppParent->elements() : currIxuppChild->elements();
   double* redCol = (block_type == LINKING_VARS_BLOCK) ? currRedColParent->elements() : currRedColChild->elements();

   if( block_type == LINKING_VARS_BLOCK )
      linkVarsBlocks[it+1].start = localNelims;
   else if( block_type == CHILD_BLOCK )
      childBlocks[it+1].start = localNelims;

   /* for every row in row in matrix */
   for( int r = 0; r < storage->m; r++ )
   {
      int start = storage->rowptr[r].start;
      int end = storage->rowptr[r].end;

      /* for every nonzero column in that row */
      for(int k = start; k < end; ++k )
      {
         const int col = storage->jcolM[k];
         const double mat_entry = storage->M[k];

         /* remove all small entries */
         if( fabs( mat_entry ) < tol_matrix_entry )
         {
            storeRemovedEntryIndex(r, storage->jcolM[k], it, block_type);
            updateAndSwap(storage, r, k, end, redCol, nelims);
         }
         else if( fabs( mat_entry ) < tolerance1 && ixuppElems[col] != 0.0 && ixlowElems[col] != 0.0
               && fabs( mat_entry ) * (xuppElems[col] - xlowElems[col]) * nnzRow[r] < tolerance2 * feastol )
         {
            if( system_type == EQUALITY_SYSTEM )
               currEqRhs->elements()[r] -= mat_entry * xlowElems[col];
            else
            {
               if( currIcupp->elements()[r] != 0.0 )
                  currIneqRhs->elements()[r] -= mat_entry * xlowElems[col];
               if( currIclow->elements()[r] != 0.0 )
                  currIneqLhs->elements()[r] -= mat_entry * xlowElems[col];
            }

            // todo third criterion? for linking constraints: call extra function to know whether we have linking cons
            // that link only two blocks (not so urgent for linking)

            storeRemovedEntryIndex(r, storage->jcolM[k], it, block_type);
            updateAndSwap(storage, r, k, end, redCol, nelims);
         }
      }
   }
   if( block_type == LINKING_VARS_BLOCK )
      linkVarsBlocks[it+1].end = localNelims;
   else if( block_type == CHILD_BLOCK )
      childBlocks[it+1].end = localNelims;

   return nelims;
}

int StochPresolverModelCleanup::removeTinyLinkingRows( int it, SystemType system_type )
{
   // Setting all the pointers correctly
   int nelims = 0;
   double tmp = 0;
   double* xlowElems = (it == -1) ? currxlowParent->elements() : currxlowChild->elements();
   double* xuppElems = (it == -1) ? currxuppParent->elements() : currxuppChild->elements();
   double* ixlowElems = (it == -1) ? currIxlowParent->elements() : currIxlowChild->elements();
   double* ixuppElems = (it == -1) ? currIxuppParent->elements() : currIxuppChild->elements();
   double* redCol = (it == -1) ? currRedColParent->elements() : currRedColChild->elements();

   // the actual work starts here:
   for( int r = 0; r < currBlmat->m; r++ )   // r: row index
   {
      int end = currBlmat->rowptr[r].end;
      int k = currBlmat->rowptr[r].start;
      while( k < end )
      {
         const int col = currBlmat->jcolM[k];
         if( fabs(currBlmat->M[k]) < tol_matrix_entry )
         {
            updateAndSwap(currBlmat, r, k, end, redCol, nelims);
            removeEntryInDynamicStorage(*currBlmatTrans, col, r, tmp);
         }
         else if( fabs(currBlmat->M[k]) < tolerance1 && ixuppElems[col] != 0.0 && ixlowElems[col] != 0.0
               && fabs(currBlmat->M[k]) * (xuppElems[col] - xlowElems[col]) * currNnzRow->elements()[r] < tolerance2 * feastol )
         {
            if( system_type == EQUALITY_SYSTEM )
               currEqRhsLink->elements()[r] -= currBlmat->M[k] * xlowElems[col];
            else
            {
               if( currIcuppLink->elements()[r] != 0.0 )
                  currIneqRhsLink->elements()[r] -= currBlmat->M[k] * xlowElems[col];
               if( currIclowLink->elements()[r] != 0.0 )
                  currIneqLhsLink->elements()[r] -= currBlmat->M[k] * xlowElems[col];
            }

            updateAndSwap(currBlmat, r, k, end, redCol, nelims);
            removeEntryInDynamicStorage(*currBlmatTrans, col, r, tmp);
         }
         k++;
      }
   }

   return nelims;
}

/**
 * set all pointers to the currently necessary data
 * If node == -1 we are in the root node
 */
bool StochPresolverModelCleanup::updatePointersForCurrentNode(int node, SystemType system_type)
{
   assert(-1 <= node && node <= nChildren );
   assert(system_type == EQUALITY_SYSTEM || system_type == INEQUALITY_SYSTEM);

   setCurrentPointersToNull(); // todo necessary? maybe set to null inside the set methods?
   setCPColumnRoot(); // todo only one kind for A and C?

   if( !setPointersForAmatBmat(node, system_type) )
      return false;

   if( node == -1 ) // at root node
      setPointersRowRoot(system_type);
   else // at child it
      setPointersRowChildNode(system_type, node);

   if( node > -1 ) // at child node
      setCPColumnChild(node);

   return true;
}

/**
 * Set the current pointers to the linking constraints.
 * If it==-1, case at F_0.
 */
bool StochPresolverModelCleanup::updateCPforLinkingRows(int it, SystemType system_type)
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

