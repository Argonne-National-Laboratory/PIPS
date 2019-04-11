/*
 * StochPresolverModelCleanup.C
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#include "StochPresolverModelCleanup.h"
#include <cmath>
#include <utility>
#include <vector>

StochPresolverModelCleanup::StochPresolverModelCleanup(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverModelCleanup::~StochPresolverModelCleanup()
{
   // todo
}


void StochPresolverModelCleanup::applyPresolving()
{
   /* ideally they are already zeroed */
   presData.resetRedCounters();

   int n_elims = 0;

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

#ifndef NDEBUG
   if( myRank == 0 )
   {
      std::cout << "Starting model cleanup..." << std::endl;
   }
   // todo this a non const function! watch out
//    countRowsCols();
#endif

   /* remove entries from A and C matrices and updates transposed systems */
   n_elims += removeTinyEntriesFromSystem(EQUALITY_SYSTEM);
   assert(verifyNnzcounters());

   n_elims += removeTinyEntriesFromSystem(INEQUALITY_SYSTEM);
   assert(verifyNnzcounters());

   // removal of redundant constraints
   int nRemovedRows = 0;
   nRemovedRows += removeRedundantRows(presProb->A, EQUALITY_SYSTEM);
   nRemovedRows += removeRedundantRows(presProb->C, INEQUALITY_SYSTEM);

   // update all nnzCounters - set reductionStochvecs to zero afterwards
   updateNnzColParent(MPI_COMM_WORLD);
   presData.resetRedCounters();

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &nRemovedRows, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);

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

/** Compute the minimal and maximal row activities of the linking rows. */
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

/* removes all small entries from the specified system
 *
 * While removing the reduction vectors in presData get set - nnzVectors are not updated an that must be
 * done using updateNnzFromReductions if needed.
 * Transposed matrices get updated in a subroutine - so after calling this method, tha matrix should
 * be in a consistent state.
 */

int StochPresolverModelCleanup::removeTinyEntriesFromSystem(SystemType system_type)
{
   assert(dynamic_cast<StochGenMatrix&>(*(presProb->A)).children.size() == presData.nRowElemsA->children.size()); // todo depends on system type
   assert(dynamic_cast<StochGenMatrix&>(*(presProb->A)).children.size() == presData.redCol->children.size()); // todo depends on system type
   assert(dynamic_cast<StochGenMatrix&>(*(presProb->A)).children.size() == (size_t) nChildren);
   assert(dynamic_cast<StochGenMatrix&>(*(presProb->C)).children.size() == (size_t) nChildren);

   /* reset the linking row reduction buffers */
   if(hasLinking(system_type))
      (system_type == EQUALITY_SYSTEM) ? resetEqRhsAdaptionsLink() : resetIneqRhsAdaptionsLink();

   int n_elims = 0;
   presData.resetRedCounters();

   int myRank;
   bool iAmDistrib;
   getRankDistributed(MPI_COMM_WORLD, myRank, iAmDistrib);

   /* reductions in root node */
   if( !nodeIsDummy(-1, system_type) )
   {
      updatePointersForCurrentNode(-1, system_type);

      /* process B0 and Bl0 */
      n_elims += removeTinyInnerLoop(-1, system_type, CHILD_BLOCK);
      if( hasLinking(system_type) )
         n_elims += removeTinyInnerLoop(-1, system_type, LINKING_CONS_BLOCK);
   }

   // todo : traffic can be reduced by only zeroing the linking var col and the linking cons row
   // and not the zero row too
   /* count eliminations in B0 and Bl0 only once */
   if( iAmDistrib && myRank != 0 )
   {
      presData.resetRedCounters();
      /* only rank 0 keeps the changes in Bl0 */
      (system_type == EQUALITY_SYSTEM) ? resetEqRhsAdaptionsLink() : resetIneqRhsAdaptionsLink();
      n_elims = 0;
   }

   // go through the children
   for( int node = 0; node < nChildren; node++ )
   {
      if( !nodeIsDummy(node, system_type) )
      {
         updatePointersForCurrentNode(node, system_type);

         /* Amat */
         n_elims += removeTinyInnerLoop(node, system_type, LINKING_VARS_BLOCK );

         /* Bmat */
         n_elims += removeTinyInnerLoop(node, system_type, CHILD_BLOCK );

         /* this has to be synchronized */
         /* Blmat */
         if( hasLinking(system_type) )
            n_elims += removeTinyInnerLoop(node, system_type, LINKING_CONS_BLOCK);
         }
   }

   //todo make MPI_COMMUNICATOR flexible?
   /* communicate the reductions */
   if( iAmDistrib )
   {
      StochVectorHandle red_row = (system_type == EQUALITY_SYSTEM) ? presData.redRowA : presData.redRowC;

      // allreduce local eliminations
      MPI_Allreduce(MPI_IN_PLACE, &n_elims, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

      // allreduce the linking variables columns
      double* red_col_linking_vars = dynamic_cast<SimpleVector*>(presData.redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(presData.redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, red_col_linking_vars, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      // allreduce B0 row
      double* red_row_b0 = dynamic_cast<SimpleVector*>(red_row->vec)->elements();
      message_size = dynamic_cast<SimpleVector*>(red_row->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, red_row_b0, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

      if( hasLinking(system_type) )
      {
         // allreduce the linking conss rows
         // non-zero counters
         double* red_row_link = dynamic_cast<SimpleVector*>(red_row->vecl)->elements();
         message_size = dynamic_cast<SimpleVector*>(red_row->vecl)->length();
         MPI_Allreduce(MPI_IN_PLACE, red_row_link, message_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );


         // rhs lhs changes (which are actually the same, so it suffices to reduce either rhs or lhs for an
         // INEQUALITY_SYSTEM
         if(system_type == EQUALITY_SYSTEM)
         {
            assert(currEqRhsLink);
            assert(presData.redRowA->vecl->n == currEqRhsLink->n);

            MPI_Allreduce(MPI_IN_PLACE, currEqRhsAdaptionsLink, currEqRhsLink->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

            // apply changes to rhs locally
            for(int i = 0; i < currEqRhsLink->n; ++i)
            {
               currEqRhsLink->elements()[i] += currEqRhsAdaptionsLink[i];
            }

            resetEqRhsAdaptionsLink();
         }
         else
         {
            assert(currIneqRhsLink);
            assert(currIneqLhsLink);
            assert(presData.redRowC->vecl->n == currIneqRhsLink->n);
            assert(presData.redRowC->vecl->n == currIneqLhsLink->n);

            MPI_Allreduce(MPI_IN_PLACE, currInEqRhsAdaptionsLink, currIneqRhsLink->n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

            // apply changes to lhs, rhs locally
            for(int i = 0; i < currIneqRhsLink->n; ++i)
            {
               currIneqRhsLink->elements()[i] += currInEqRhsAdaptionsLink[i];
               currIneqLhsLink->elements()[i] += currInEqLhsAdaptionsLink[i];
            }

            resetIneqRhsAdaptionsLink();
         }
      }
   }

   /* update local nnzCounters */
   updateNnzFromReductions(system_type);

   return n_elims;
}

/** Removes tiny entries in storage and adapts the rhs accordingly.
 *  If block_type == LINKING_VARS_BLOCK, then block Amat is considered.
 *  If block_type == CHILD_BLOCK, then block Bmat is considered. */
/* system type indicates matrix A or C, block_type indicates the block */
// todo what is the proper order for criterion 1 to 3?
// todo for criterion 3 - should ALL eliminations be considered?
int StochPresolverModelCleanup::removeTinyInnerLoop( int node, SystemType system_type, BlockType block_type)
{
   if(nodeIsDummy(node, system_type))
      return 0;

   SparseStorageDynamic* mat = NULL;
   SparseStorageDynamic* mat_transp = NULL;

   SimpleVector* lhs = NULL;
   SimpleVector* lhs_idx = NULL;
   SimpleVector* rhs = NULL;
   SimpleVector* rhs_idx = NULL;
   SimpleVector* x_lower = NULL;
   SimpleVector* x_lower_idx = NULL;
   SimpleVector* x_upper = NULL;
   SimpleVector* x_upper_idx = NULL;
   SimpleVector* redCol = NULL;
   SimpleVector* nnzRow = NULL;

   /* set matrix */
   if( block_type == CHILD_BLOCK )
   {
      mat = currBmat;
      mat_transp = currBmatTrans;
   }
   else if( block_type == LINKING_VARS_BLOCK )
   {
      assert(node != -1);
      mat = currAmat;
      mat_transp = currAmatTrans;
   }
   else if( block_type == LINKING_CONS_BLOCK)
   {
      mat = currBlmat;
      mat_transp = currBlmatTrans;
   }

   /* set lhs rhs */
   /* for linking conss we do not apply the updates directly but rather buffer them
    * Note: the arrays will stay alive after deleting the simpleVecs
    */
   if(system_type == EQUALITY_SYSTEM)
   {
      if(block_type == LINKING_CONS_BLOCK)
      {
         assert(currEqRhsAdaptionsLink);
         //rhs = lhs = lhs_idx = rhs_idx = currEqRhsLink;
         rhs = lhs = lhs_idx = rhs_idx = new SimpleVector(currEqRhsAdaptionsLink, currEqRhsLink->n);
      }
      else
      {
         rhs = lhs = lhs_idx = rhs_idx = currEqRhs;
      }
   }
   else if(system_type == INEQUALITY_SYSTEM)
   {
      if(block_type == LINKING_CONS_BLOCK)
      {
         assert(currInEqRhsAdaptionsLink); assert(currInEqLhsAdaptionsLink);
//         rhs = currIneqRhsLink;
         rhs = new SimpleVector(currInEqRhsAdaptionsLink, currIneqRhsLink->n);
         rhs_idx = currIcuppLink;
//         lhs =  currIneqLhsLink;
         lhs = new SimpleVector(currInEqLhsAdaptionsLink, currIneqLhsLink->n);
         lhs_idx = currIclowLink;
      }
      else
      {
         rhs = currIneqRhs;
         rhs_idx = currIcupp;
         lhs =  currIneqLhs;
         lhs_idx = currIclow;
      }
   }

   /* set variables */
   if(block_type == CHILD_BLOCK || block_type == LINKING_CONS_BLOCK)
   {
      if(node == -1)
      {
         x_lower = currxlowParent;
         x_lower_idx = currIxlowParent;

         x_upper = currxuppParent;
         x_upper_idx = currIxuppParent;
      }
      else
      {
         x_lower = currxlowChild;
         x_lower_idx = currIxlowChild;

         x_upper = currxuppChild;
         x_upper_idx = currIxuppChild;
      }
   }
   else if(block_type == LINKING_VARS_BLOCK)
   {
      assert(node != -1);
      x_lower = currxlowParent;
      x_lower_idx = currIxlowParent;

      x_upper = currxuppParent;
      x_upper_idx = currIxuppParent;

   }

   /* set reduction vectors */
   if( block_type == CHILD_BLOCK || block_type == LINKING_CONS_BLOCK )
   {
      if( node == -1 )
      {
         redCol = currRedColParent;
      }
      else
      {
         redCol = currRedColChild;
      }
   }
   else if( block_type == LINKING_VARS_BLOCK )
   {
      assert(node != -1);
      redCol = currRedColParent;
   }

   /* set non-zero vectors */
   if(node == -1)
   {
      if(block_type == LINKING_CONS_BLOCK)
      {
         assert(hasLinking(system_type));
         nnzRow = currNnzRowLink;
      }
      else
      {
         assert(block_type == CHILD_BLOCK);
         nnzRow = currNnzRow;
      }
   }
   else
   {
      if(block_type == LINKING_CONS_BLOCK)
         nnzRow = currNnzRowLink;
      else
         nnzRow = currNnzRow;
   }

   bool linking_row = (block_type == LINKING_CONS_BLOCK);

   // Setting all the pointers correctly
   int nelims = 0;
   SparseStorageDynamic* storage = mat;
   SparseStorageDynamic* storage_transp = mat_transp;
   assert(storage_transp);

   std::vector<std::pair<int, int> > eliminated_entries;

   /* for every row in row in matrix */
   for( int r = 0; r < storage->m; r++ )
   {
      // todo : is this per row?
      double total_sum_modifications_row = 0.0;

      int start = storage->rowptr[r].start;
      int end = storage->rowptr[r].end;

      /* for every nonzero column in that row */
      for(int k = start; k < end; ++k )
      {
         const int col = storage->jcolM[k];
         const double mat_entry = storage->M[k];

         /* remove all small entries */
         if( fabs( mat_entry ) < tol_matrix_entry ) // todo bugged
         {
            std::pair<int,int> entry(r, col);
            eliminated_entries.push_back(entry);

            updateAndSwap(storage, r, k, end, redCol->elements(), nelims, linking_row);
            continue;
         }

         /* remove entries where their corresponding variables have valid lower and upper bounds, that overall do not have a real influence though */
         if( x_upper_idx->elements()[col] != 0.0 && x_lower_idx->elements()[col] != 0.0 )
         {

            double x = nnzRow->elements()[r];
            x += x_upper->elements()[col];
            x += x_lower->elements()[col];

            if( (fabs( mat_entry ) < tolerance1 && fabs( mat_entry ) * (x_upper->elements()[col] - x_lower->elements()[col]) * nnzRow->elements()[r] < tolerance2 * feastol ))
            {
               if( system_type == EQUALITY_SYSTEM )
               {
                  rhs->elements()[r] -= mat_entry * x_lower->elements()[col];
               }
               else
               {
                  /* if has upper/lower bound */
                  if( rhs_idx->elements()[r] != 0.0 )
                     rhs->elements()[r] -= mat_entry * x_lower->elements()[col];
                  if( lhs_idx->elements()[r] != 0.0 )
                     lhs->elements()[r] -= mat_entry * x_lower->elements()[col];
               }

               std::pair<int,int> entry(r, col);
               eliminated_entries.push_back(entry);

               updateAndSwap(storage, r, k, end, redCol->elements(), nelims, linking_row);
               continue;
            }
         }

         // todo third criterion? for linking constraints: call extra function to know whether we have linking cons
         // that link only two blocks (not so urgent for linking)
         if( false ) //todo if not linking constraints
         {
            /* if valid lower and upper bounds */
            if( x_upper_idx->elements()[col] != 0.0 && x_lower_idx->elements()[col] != 0.0 ){
               if( total_sum_modifications_row + (fabs(mat_entry) * (x_upper->elements()[col] - x_lower->elements()[col])) < 1.0e-1 * feastol)
               {
                  total_sum_modifications_row += fabs(mat_entry) * (x_upper->elements()[col] - x_lower->elements()[col]);

                  std::pair<int,int> entry(r, col);
                  eliminated_entries.push_back(entry);

                  updateAndSwap(storage, r, k, end, redCol->elements(), nelims, linking_row);
                  continue;
               }
            }
         }
         /* not removed */
      }
   }

   /* update the transposed if existent */
   if(storage_transp)
      updateTransposedSubmatrix(storage_transp, eliminated_entries);

   /* free artificial SimpleVectors used for buffering linking constraint changes */
   if( block_type == LINKING_CONS_BLOCK )
   {
      if( system_type == EQUALITY_SYSTEM )
      {
         assert(currEqRhsAdaptionsLink);
         //rhs = lhs = lhs_idx = rhs_idx = currEqRhsLink;
         delete rhs;
      }
      else
      {
         delete rhs;
         delete lhs;
      }
   }

   return nelims;
}

