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

StochPresolverSingletonRows::StochPresolverSingletonRows(PresolveData& presData) :
      StochPresolverBase(presData)
{
   // todo
}

StochPresolverSingletonRows::~StochPresolverSingletonRows()
{
   // todo
}

void StochPresolverSingletonRows::applyPresolving()
{

   assert(true);
   // todo : assert some stuff :P
   // assert(rootnodedatainsync)
   //assert nnxconters correct
   //   clearNewBoundsParent(); // rather - assert is empty...
   assert(presData.reductionsEmpty());
   presData.resetRedCounters();

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- Before singleton Row Presolving:" << std::endl;
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
         //assert(n_singleton_equality == 0);

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

#ifdef TIMING
   if( myRank == 0 )
   std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;
#endif

   if( myRank == 0 )
      std::cout << "Global objOffset is now: " << presData.getObjOffset() << std::endl;

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- After singleton Row Presolving:" << std::endl;
   countRowsCols();
#endif
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
      abortInfeasible(MPI_COMM_WORLD );
   int a = 0;
   int b = 0;
   updateLinkingVarsBlocks(a, b);

   /* allreduce rhs lhs changes and non-zero counters */
   allreduceAndApplyNnzReductions(INEQUALITY_SYSTEM);
   allreduceAndApplyNnzReductions(EQUALITY_SYSTEM);

   allreduceAndApplyRhsLhsReductions(system_type);

   /* allreduce new var bounds */
   allreduceAndUpdateVarBounds();

   // todo assert all reductions and stuff empty
}

/* Finds and stores fixations in colAdapParent, deletes singletons in INEQUALITY_SYSTEM and finds and stores fixations resulting from
 * bound changes. Stores the newly found bounds in newBOundsParent.
 *
 * Later both of the stored changes have to be applied to the whole system.
 */
void StochPresolverSingletonRows::procSingletonRowRoot(SystemType system_type)
{
   /* B0 node */
   processSingletonBlock(system_type, CHILD_BLOCK, -1);

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
      return;

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

   SimpleVector* nnz_col = (block_type == LINKING_VARS_BLOCK || node == -1) ? currNnzColParent : currNnzColChild;
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
      assert(node != -1);
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
      assert(block_type == CHILD_BLOCK);
      matrix = currBmat;
      matrix_transp = currBmatTrans;
   }

   /* set reduction vectors */
   SimpleVector* redCol = NULL;
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
   assert(redCol);

   SimpleVector* redRow = (block_type == LINKING_CONS_BLOCK) ? currRedRowLink : currRedRow;

   bool linking_row = false;
   if(block_type == LINKING_CONS_BLOCK)
      linking_row = true;


   /* go through matrix - store variables for deletion - always delete singleton inequality rows right away */
   for( int i = 0; i < nnz_row->length(); ++i )
   {
      /* if singleton row entry is in current block */
      if( nnz_row->elements()[i] - redRow->elements()[i] == 1.0 && matrix->rowptr[i].start + 1 == matrix->rowptr[i].end)
      {
         int colIdx = -1;
         double aik = 0.0;

         getValuesForSR(*matrix, i, colIdx, aik);

         assert(!PIPSisEQ(aik, 0.0));

         if( system_type == EQUALITY_SYSTEM )
         {
            const double fixation_value = curr_eq_rhs->elements()[i] / aik;

            if( (ixlow[colIdx] != 0.0 && PIPSisLT(fixation_value, xlow[colIdx]))
                  || (ixupp[colIdx] != 0.0 && PIPSisLT(xupp[colIdx], fixation_value)) )
            {
               std::cout << "Singleton Row Presolving detected infeasibility : fixation of variable to invalid value in EQUALITY_SYSTEM B0" << std::endl;
               std::cout << "variable index: " << colIdx << "\tvalue: " << fixation_value << "\tlower bound: " << xlow[colIdx] << "\tupper bound: " << xupp[colIdx] << std::endl;
               abortInfeasible(MPI_COMM_WORLD );
            }
            else
            {
               /* for Amat we store deletions - collect them and apply them later */
               if( node == -1 || block_type == LINKING_VARS_BLOCK)
               {
                  /* only 0 process stores fixations in root node B0/Bl0 */
                  if(myRank == 0 && node == -1)
                  {
                     assert(block_type != LINKING_VARS_BLOCK);
                     storeColValInColAdaptParent(colIdx, fixation_value);
                  }

                  if(block_type == LINKING_VARS_BLOCK)
                  {
                     assert(node != -1);
                     storeColValInColAdaptParent(colIdx, fixation_value);
                  }
               }
               else
               {
                  /* delete variable */
                  deleteNonlinkColumnFromSystem(node, colIdx, fixation_value);
               }
            }
         }
         else if( system_type == INEQUALITY_SYSTEM )
         {
            double new_xlow = -std::numeric_limits<double>::max();
            double new_xupp = std::numeric_limits<double>::max();
            double val = 0.0;

            calculateNewBoundsOnVariable(new_xlow, new_xupp, iclow->elements()[i], clow->elements()[i], icupp->elements()[i], cupp->elements()[i], aik);

            if( newBoundsImplyInfeasible(new_xlow, new_xupp, colIdx, ixlow, ixupp, xlow, xupp) )
               abortInfeasible(MPI_COMM_WORLD );
            else if( newBoundsFixVariable(val, new_xlow, new_xupp, colIdx, ixlow, ixupp, xlow, xupp) )
            {
               if(myRank == 0)
                  std::cout << val << std::endl;
               /* for Amat we store deletions - collect them and apply them later */
               if( node == -1 || block_type == LINKING_VARS_BLOCK)
               {
                  if(myRank == 0 && node == -1)
                     storeColValInColAdaptParent(colIdx, val);
                  if(block_type == LINKING_VARS_BLOCK)
                     storeColValInColAdaptParent(colIdx, val);
               }
               else
               {
                  /* delete variable */
                  deleteNonlinkColumnFromSystem(node, colIdx, val);
               }
            }
            else
            {
               tightenBounds(new_xlow, new_xupp, ixlow[colIdx], xlow[colIdx], ixupp[colIdx], xupp[colIdx]);

               /* remove entry and update reductions */
               removeEntryInDynamicStorage(*matrix, i, colIdx, val);

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

               // remove entry a_ik in transposed matrix as well
               removeEntryInDynamicStorage(*matrix_transp, colIdx, i, val);
               assert(matrix->rowptr[i].start == matrix->rowptr[i].end);
            }
         }
      }
   }
}

bool StochPresolverSingletonRows::tightenBounds(double new_xlow, double new_xupp, double& ixlow, double& old_xlow, double& ixupp, double& old_xupp) const
{
   assert( !PIPSisEQ(new_xlow, new_xupp) );
   bool tightened = false;

   if( ixlow != 0.0 && PIPSisLT(old_xlow, new_xlow) )
   {
      old_xlow = new_xlow;
      tightened = true;
   }
   else if( ixlow == 0.0 && new_xlow > -std::numeric_limits<double>::max() )
   {
      old_xlow = new_xlow;
      ixlow = 1.0;
      tightened = true;
   }

   if( ixupp != 0.0 && PIPSisLT(new_xupp, old_xupp) )
   {
      old_xupp = new_xupp;
      tightened = true;
   }
   else if( ixupp == 0.0 && new_xupp < std::numeric_limits<double>::max() )
   {
      old_xupp = new_xupp;
      ixupp = 1.0;
      tightened = true;
   }

   assert( !PIPSisEQ(new_xlow, new_xupp) );
   return tightened;
}

void StochPresolverSingletonRows::calculateNewBoundsOnVariable(double& new_xlow, double& new_xupp, const double& iclow, const double& clow,
      const double& icupp, const double& cupp, double aik) const
{

   assert( !PIPSisZero(aik) );

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
//   std::cout << "newxlow: " << new_xlow << "\tnewxupp: " << new_xupp << "\ticlow: " << iclow << "\ticupp: " << icupp <<
//         "\tcupp: " << cupp << "\tclow: " << clow << "\taik: " << aik << std::endl;
}

/** Should be called right after doSingletonRowsC() or another method that stores
 * information to update in newBoundsParent.
 * Updates the bounds on the linking variables.
 */
void StochPresolverSingletonRows::updateLinkingVarsBounds()
{
   setCPColumnRoot();
   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();

   // apply updated newBoundsParent to the variable bounds.
   for( int i = 0; i < getNumberNewBoundsParent(); i++ )
   {
      XBOUNDS newbounds = getNewBoundsParent(i);
      setNewBounds(newbounds.colIdx, newbounds.newxlow, newbounds.newxupp,
            ixlow, xlow, ixupp, xupp);
   }
   clearNewBoundsParent();
}

void StochPresolverSingletonRows::getValuesForSR( SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const
{
   const int indexK = storage.rowptr[rowIdx].start;
   colIdx = storage.jcolM[indexK];
   aik = storage.M[indexK];

   assert(storage.rowptr[rowIdx].start +1 == storage.rowptr[rowIdx].end);
   assert(aik != 0.0);
}
