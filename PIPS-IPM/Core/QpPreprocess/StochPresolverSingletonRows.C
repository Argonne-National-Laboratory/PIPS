/*
 * StochPresolverSingletonRows.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

//#define PIPS_DEBUG
#include "StochPresolverSingletonRows.h"
#include <limits>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

StochPresolverSingletonRows::StochPresolverSingletonRows(PresolveData& presData, const sData& origProb) :
      StochPresolverBase(presData, origProb)
{
   // todo
}

StochPresolverSingletonRows::~StochPresolverSingletonRows()
{
   // todo
}

// todo print singleton rows removed
void StochPresolverSingletonRows::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifndef NDEBUG
   if( myRank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before singleton Row Presolving:" << std::endl;
   }
   countRowsCols();
#endif

   indivObjOffset = 0.0;
   int n_singleton_equality = 0;
   int n_singleton_inequality = 0;

   countSingletonRows(n_singleton_equality, n_singleton_inequality);

#ifndef NDEBUG
   if(myRank == 0)
      std::cout << "Initially found " << n_singleton_equality << " equality and "
         << n_singleton_inequality << " inequality singletons" << std::endl;
#endif

   int iter = 0;

   // main loop:
   while(n_singleton_equality + n_singleton_inequality > 0 && iter < maxIterSR )
   {
      /* eliminate all singleton rows in equality system */
      if( n_singleton_equality > 0 )
      {
         // main method:
         doSingletonRows(n_singleton_equality, n_singleton_inequality,
               EQUALITY_SYSTEM);
      }
      else if( n_singleton_inequality > 0 )
      {
         assert(n_singleton_equality == 0);

         /* main method: */
         doSingletonRows(n_singleton_inequality, n_singleton_equality,
               INEQUALITY_SYSTEM);
      }

      iter++;

      countSingletonRows(n_singleton_equality, n_singleton_inequality);
   }

   assert( (n_singleton_equality == 0 && n_singleton_inequality == 0) || iter >= maxIterSR);

   // Sum up individual objOffset and then add it to the global objOffset:
   sumIndivObjOffset();
   presData.addObjOffset(indivObjOffset);
   indivObjOffset = 0.0;

   if( myRank == 0 )
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- After singleton row presolving:" << std::endl;
   countRowsCols();
   if(myRank == 0)
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);
}

/** Does one round of singleton rows presolving for system A or C
 *
 * the blocks B,D,Fi,Gi (the blocks Bmat and Blmat of both A and C), the fixation and updating
 * of the columns is done. The fixed variables in one of the Amat blocks are stored in the
 * member variable colAdaptParent. Updating the blocks A,C,F0,G0 using colAdaptParent happens
 * in updateLinkingVarsBlocks() which should be called after this method.
 * Returns the number of newly found singleton rows (equality/inequality system) during adaption of B,D,Fi,Gi.
 */
void
StochPresolverSingletonRows::doSingletonRows(int& n_sing_sys, int& n_sing_other_sys, SystemType system_type)
{
   n_sing_sys = 0;

   /* processes root node - finds vars to delete and updates bounds */
   procSingletonRowRoot(system_type);

   /* remove singletons from children */
   for( int node = 0; node < nChildren; node++ )
   {
      procSingletonRowChild(node, n_sing_sys, n_sing_other_sys, system_type);
   }

   /* allreduce found deletions for linking variables */
   if( !presData.combineColAdaptParent() )
   {
      abortInfeasible(MPI_COMM_WORLD );
   }
   int a = 0;
   int b = 0;
   updateLinkingVarsBlocks(a, b);

   /* allreduce rhs lhs changes and rdeuction counters and apply them to non-zero counters */
   allreduceAndApplyNnzReductions(INEQUALITY_SYSTEM);
   allreduceAndApplyNnzReductions(EQUALITY_SYSTEM);

   allreduceAndApplyRhsLhsReductions(system_type);

   /* allreduce new var bounds */
   allreduceAndUpdateVarBounds();
}

/* Finds and stores fixations in colAdapParent, deletes singletons in INEQUALITY_SYSTEM and finds and stores fixations resulting from
 * bound changes. Stores the newly found bounds in newBOundsParent.
 *
 * Later both of the stored changes have to be applied to the whole system.
 */
void StochPresolverSingletonRows::procSingletonRowRoot(SystemType system_type)
{
   /* B0 node */
   processSingletonBlock(system_type, LINKING_VARS_BLOCK, -1);

   /* linking vars Bl */
   processSingletonBlock(system_type, LINKING_CONS_BLOCK, -1);
}

/* Processing the singleton rows in child it, more precisely, goes through all singleton rows in Amat and Bmat and Blmat
 *
 * Those in Amat are stored in colAdaptParent for later processing.
 * Those in Bmat are stored in colAdaptLinkBlock. Furthermore, the corresponding fixed variables (columns)
 * in Bmat and in Blmat are removed.
 * Using this colAdaptLinkBlock, the variables (columns) are removed from the inequalities Bmat, Blmat as well.
 */
void StochPresolverSingletonRows::procSingletonRowChild(int node, int& n_singleton_sys, int& n_singleton_other_sys, SystemType system_type)
{
   if( nodeIsDummy(node, system_type) )
   {
      return;
   }

   /* Amat - store deletions */
   processSingletonBlock(system_type, LINKING_VARS_BLOCK, node);

   /* Bmat */
   processSingletonBlock(system_type, CHILD_BLOCK, node);

   /* Blmat */
   processSingletonBlock(system_type, LINKING_CONS_BLOCK, node);
}

/** Finds singleton rows in the specified block and processes them
 *
 * If the processed block is from the equality system, found singleton entries will get stored for later sync and deletion.
 * If the processed block is from the inequality system, found singleton rows will get deleted, the deletion is stored (for an update of the non-zero
 * counters) and bounds are adapted accordingly.
 *
 * Method works for all types of blocks. Synchronization has to be done accordingly.
 */
void StochPresolverSingletonRows::processSingletonBlock(SystemType system_type, BlockType block_type, int node)
{

   if( block_type == LINKING_CONS_BLOCK )
      if( !hasLinking(system_type) )
         return;
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   updatePointersForCurrentNode(node, system_type);

   double* ixlow = (block_type == LINKING_VARS_BLOCK || node == -1) ? currIxlowParent->elements() : currIxlowChild->elements();
   double* ixupp = (block_type == LINKING_VARS_BLOCK || node == -1) ? currIxuppParent->elements() : currIxuppChild->elements();
   double* xlow = (block_type == LINKING_VARS_BLOCK || node == -1) ? currxlowParent->elements() : currxlowChild->elements();
   double* xupp = (block_type == LINKING_VARS_BLOCK || node == -1) ? currxuppParent->elements() : currxuppChild->elements();

   SimpleVector* nnz_row = (block_type == LINKING_CONS_BLOCK) ? currNnzRowLink : currNnzRow;

   SimpleVector* curr_eq_rhs = (block_type == LINKING_CONS_BLOCK) ? currEqRhsLink : currEqRhs;

   SimpleVector* iclow = (block_type == LINKING_CONS_BLOCK) ? currIclowLink : currIclow;
   SimpleVector* clow = (block_type == LINKING_CONS_BLOCK) ? currIneqLhsLink : currIneqLhs;
   SimpleVector* icupp = (block_type == LINKING_CONS_BLOCK) ? currIcuppLink : currIcupp;
   SimpleVector* cupp = (block_type == LINKING_CONS_BLOCK) ? currIneqRhsLink : currIneqRhs;

   SparseStorageDynamic* matrix;
   SparseStorageDynamic* matrix_transp;

   if( block_type == LINKING_VARS_BLOCK )
   {
      matrix = currAmat;
      matrix_transp = currAmatTrans;
   }
   else if( block_type == LINKING_CONS_BLOCK )
   {
      matrix = currBlmat;
      matrix_transp = currBlmatTrans;
   }
   else
   {
      assert(node != -1);
      matrix = currBmat;
      matrix_transp = currBmatTrans;
   }

   /* set reduction vectors */
   SimpleVector* redCol = (block_type == LINKING_VARS_BLOCK || node == -1) ? currRedColParent : currRedColChild;
   assert(redCol);
   SimpleVector* redRow = (block_type == LINKING_CONS_BLOCK) ? currRedRowLink : currRedRow;

   assert(nnz_row->length() == matrix->m);
   assert(nnz_row->length() == matrix_transp->n);
   /* go through matrix - store variables for deletion - always delete singleton inequality rows right away */
   for( int i = 0; i < nnz_row->length(); ++i )
   {
      /* if singleton row entry is in current block */
      if( nnz_row->elements()[i] - redRow->elements()[i] == 1.0 && matrix->rowptr[i].start + 1 == matrix->rowptr[i].end)
      {

         int colIdx = -1;
         double aik = 0.0;

         getValuesForSR(*matrix, i, colIdx, aik);
         assert( !PIPSisEQ(aik, 0.0) );

         /* if in equality system fix variable */
         if( system_type == EQUALITY_SYSTEM )
         {

            const double rhs = (block_type == LINKING_CONS_BLOCK) ? curr_eq_rhs->elements()[i] + currEqRhsAdaptionsLink[i] : curr_eq_rhs->elements()[i];
            const double fixation_value = rhs / aik;

            if( !variableFixationValid(fixation_value, ixlow[colIdx], xlow[colIdx], ixupp[colIdx], xupp[colIdx], true) )
               abortInfeasible(MPI_COMM_WORLD);

            /* for Amat we store deletions - collect them and apply them later */
            if( node == -1 || block_type == LINKING_VARS_BLOCK )
            {
               /* only 0 process stores fixations in root node B0/Bl0 - this is to reduce MPI communications - it is not necessary */
               if( myRank == 0 && node == -1 )
                  storeColValInColAdaptParent(colIdx, fixation_value);
               else if( block_type == LINKING_VARS_BLOCK )
                  storeColValInColAdaptParent(colIdx, fixation_value);
            }
            else
            {
               /* delete variable */
               deleteNonlinkColumnFromSystem(node, colIdx, fixation_value);

               /* line from matrix removed an rhs set to zero */
               assert(matrix->rowptr[i].start == matrix->rowptr[i].end);

               if( block_type != LINKING_CONS_BLOCK)
                  assert( PIPSisEQ(curr_eq_rhs->elements()[i], 0.0) );
               else
                  assert( PIPSisEQ(curr_eq_rhs->elements()[i] + currEqRhsAdaptionsLink[i], 0.0) );

               assert( nnz_row->elements()[i] - redRow->elements()[i] == 0.0 );
            }
         }
         else if( system_type == INEQUALITY_SYSTEM )
         {
            double new_xlow = -std::numeric_limits<double>::infinity();
            double new_xupp = std::numeric_limits<double>::infinity();
            double fixation_value = 0.0;

            const double lhs = (block_type == LINKING_CONS_BLOCK) ? clow->elements()[i] + currInEqLhsAdaptionsLink[i] : clow->elements()[i];
            const double rhs = (block_type == LINKING_CONS_BLOCK) ? cupp->elements()[i] + currInEqRhsAdaptionsLink[i] : cupp->elements()[i];

            calculateNewBoundsOnVariable(new_xlow, new_xupp, iclow->elements()[i], lhs, icupp->elements()[i], rhs, aik);

            double llow = xlow[colIdx];
            double uuppp = xupp[colIdx];
            if( newBoundsImplyInfeasible(new_xlow, new_xupp, colIdx, ixlow, ixupp, xlow, xupp) )
               abortInfeasible(MPI_COMM_WORLD );
            else if( newBoundsFixVariable(fixation_value, new_xlow, new_xupp, colIdx, ixlow, ixupp, xlow, xupp) )
            {
               if( !variableFixationValid(fixation_value, ixlow[colIdx], xlow[colIdx], ixupp[colIdx], xupp[colIdx], true) )
               {
                  std::cout << new_xlow << "\t" << new_xupp << "\t" << iclow->elements()[i] << "\t" << lhs << "\t" << icupp->elements()[i] << "\t" <<  rhs << "\t" << aik << std::endl;
                  std::cout << xlow[colIdx] << "\t" << xupp[colIdx] << std::endl;
                  std::cout << llow << "\t" << uuppp << std::endl;
                  std::cout << block_type << "\t" << ixlow[colIdx] << "\t" << ixupp[colIdx] << "\t" << lhs << "\t" << rhs << std::endl;
                  abortInfeasible(MPI_COMM_WORLD);
               }

               /* for Amat we store deletions - collect them and apply them later */
               if( node == -1 || block_type == LINKING_VARS_BLOCK)
               {
                  if(myRank == 0 && node == -1)
                     storeColValInColAdaptParent(colIdx, fixation_value);
                  else if(block_type == LINKING_VARS_BLOCK)
                     storeColValInColAdaptParent(colIdx, fixation_value);
               }
               else
               {
                  /* delete variable */
                  deleteNonlinkColumnFromSystem(node, colIdx, fixation_value);

                  /* matrix deleted from storage and rhs lhs set to zero */
                  assert( matrix->rowptr[i].start == matrix->rowptr[i].end );
                  assert( nnz_row->elements()[i] - redRow->elements()[i] == 0.0 );

                  if(iclow->elements()[i] != 0.0 && block_type != LINKING_CONS_BLOCK)
                     assert( PIPSisLE(clow->elements()[i], 0.0));
                  else if(iclow->elements()[i] != 0.0)
                     assert( PIPSisLE(clow->elements()[i] + currInEqLhsAdaptionsLink[i], 0.0) );

                  if(icupp->elements()[i] != 0.0 && block_type != LINKING_CONS_BLOCK)
                     assert( PIPSisLE(0.0, cupp->elements()[i]) );
                  else if(icupp->elements()[i] != 0.0)
                     assert( PIPSisLE(0.0, cupp->elements()[i] + currInEqRhsAdaptionsLink[i]) );
               }
            }
            else
            {
               tightenBounds(new_xlow, new_xupp, ixlow[colIdx], xlow[colIdx], ixupp[colIdx], xupp[colIdx]);

               /* remove entry and update reductions */
               double entry_matrix = 0.0;
               double entry_matrix_transposed = 0.0;
               bool removed_matrix = false;
               bool removed_matrix_transposed = false;

               // todo adjust rhs lhs! // todo check : when removing variable will bounds of var be 0 ? should be i guess.. // todo check when removing row : are rhs lhs bounds still valid?
               removed_matrix = removeEntryInDynamicStorage(*matrix, i, colIdx, entry_matrix);
               removed_matrix_transposed = removeEntryInDynamicStorage(*matrix_transp, colIdx, i, entry_matrix_transposed);

               assert(removed_matrix);
               assert(removed_matrix_transposed);
               assert(entry_matrix == entry_matrix_transposed);
               assert(entry_matrix != 0.0);

               if( node == -1 )
               {
                  if( myRank == 0 )
                  {
                     redCol->elements()[colIdx]++;
                     redRow->elements()[i]++;
                  }
               }
               else
               {
                  redCol->elements()[colIdx]++;
                  redRow->elements()[i]++;
               }

               if(node != -1 || (node == -1 && myRank == 0))
                  assert( nnz_row->elements()[i] - redRow->elements()[i] == 0.0 );
               assert(matrix->rowptr[i].start == matrix->rowptr[i].end);
               }
         }
      }
   }
}

// todo move
void StochPresolverSingletonRows::calculateNewBoundsOnVariable(double& new_xlow, double& new_xupp, const double& iclow, const double& clow,
      const double& icupp, const double& cupp, double aik) const
{
   assert( !PIPSisZero(aik) );
   new_xlow = -std::numeric_limits<double>::infinity();
   new_xupp = std::numeric_limits<double>::infinity();

   if( PIPSisLT(0.0, aik) )
   {
      if( iclow != 0.0 )
         new_xlow = clow / aik;
      if( icupp != 0.0 )
         new_xupp = cupp / aik;
   }
   else
   {
      if( icupp != 0.0 )
         new_xlow = cupp / aik;
      if( iclow != 0.0 )
         new_xupp = clow / aik;
   }
}

void StochPresolverSingletonRows::getValuesForSR( SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const
{
   const int indexK = storage.rowptr[rowIdx].start;
   colIdx = storage.jcolM[indexK];
   aik = storage.M[indexK];

   assert(storage.rowptr[rowIdx].start +1 == storage.rowptr[rowIdx].end);
   assert(aik != 0.0);
}
