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

StochPresolverSingletonRows::StochPresolverSingletonRows(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverSingletonRows::~StochPresolverSingletonRows()
{
 // todo
}


void StochPresolverSingletonRows::applyPresolving()
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- Before singleton Row Presolving:" << std::endl;
   countRowsCols();
#endif

   indivObjOffset = 0.0;
   int newSREq = 0;
   int newSRIneq = 0;

   presData.resetRedCounters();
   clearNewBoundsParent();

   newSREq = initSingletonRows(EQUALITY_SYSTEM);
   synchronize(newSREq);
   if( myRank == 0 )
      PIPSdebugMessage("Found %d singleton rows in equality system A. \n", newSREq);

   int iter = 0;
   int globalIter = 0;

   // main loop:
   while( (newSREq + newSRIneq > 0 && iter < maxIterSR) || globalIter == 0 )
   {
      while( newSREq > 0 && iter < maxIterSR)
      {
         if( globalIter > 0 )
            initSingletonRows(EQUALITY_SYSTEM);
         // main method:
         doSingletonRows(newSREq, newSRIneq, EQUALITY_SYSTEM);

         // update the linking variable blocks (A,C,F,G) with the fixations found in doSingletonRowsA:
         updateLinkingVarsBlocks(newSREq, newSRIneq);
         synchronizeSum(newSREq, newSRIneq);

         if( myRank == 0 )
            PIPSdebugMessage("Found new singleton rows that were just created: %d in A and %d in C \n", newSREq, newSRIneq);
         iter++;
      }
      newSREq = 0;
      if( globalIter == 0 )
      {
         newSRIneq = initSingletonRows(INEQUALITY_SYSTEM);
         synchronize(newSRIneq);
         if( myRank == 0 )
            PIPSdebugMessage("Found %d singleton rows in inequality system C. \n", newSRIneq);
      }
      while( newSRIneq > 0 && iter < maxIterSR)
      {
         // if( myRank == 0 ) PIPSdebugMessage("SR(Inequality) loop at iter %d and globalIter %d \n", iter, globalIter);
         if( globalIter > 0 )
         {
            newSRIneq = initSingletonRows(INEQUALITY_SYSTEM);
            // only for debugging:
            synchronize(newSRIneq);
            if( myRank == 0 )
               PIPSdebugMessage("Found %d singleton rows in inequality system C. \n", newSRIneq);
         }
         // main method:
         doSingletonRows(newSRIneq, newSREq, INEQUALITY_SYSTEM);

         // update the variable bounds for the linking variables:
         updateLinkingVarsBounds();
         // update the linking variable blocks (A,C,F,G) with the fixations found in doSingletonRowsC:
         updateLinkingVarsBlocks(newSREq, newSRIneq);
         synchronizeSum(newSREq, newSRIneq);

         if( myRank == 0 )
            PIPSdebugMessage("Found new singleton rows that were just created: %d in A (at most) and %d in C \n", newSREq, newSRIneq);
         iter++;
      }
      newSRIneq = 0;
      globalIter++;

      presData.resetRedCounters();
      int newSREqLink = 0;
      int newSRIneqLink = 0;

      // todo what about inequalities?
      doSingletonLinkRows(newSREqLink, newSRIneqLink);
      synchronizeSum(newSREqLink, newSRIneqLink);
      newSREq += newSREqLink;
      newSRIneq += newSRIneqLink;
   }

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

/** Initializes the singletonRows list (actually a vector<int> singletonRows)
 * and the blocks (int*) pointing to the start and end indices in singletonRows
 * for each block (parent, children, linking rows).
 * Attention: there is no communication over the processes to adapt the number of
 * singleton rows found (for performance reasons). If this is necessary, a simple
 * MPI_Allreduce call should be used right after calling this method.
 * Returns the number of singleton rows found (might be different for each process!).
 *
 * @note: only rank 0 counts singleton rows in B0 and linking rows
 * todo : check linking row
 */
int StochPresolverSingletonRows::initSingletonRows(SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   int nSingletonRows = 0;

   StochVectorHandle nnz_vec =
         (system_type == EQUALITY_SYSTEM) ?
               presData.nRowElemsA : presData.nRowElemsC;

   assert(presData.getNumberSR() == 0);

   /* B0 */
   SimpleVector* nnz_block = dynamic_cast<SimpleVector*>(nnz_vec->vec);
   const int n_singletons_B0 = initSingletonRowsBlock(-1, nnz_block);
   if( myRank == 0 )
      nSingletonRows += n_singletons_B0;

   /* child rows Ai Bi */
   assert( (int)nnz_vec->children.size() == nChildren);
   for( size_t it = 0; it < nnz_vec->children.size(); it++ )
   {
      nnz_block = dynamic_cast<SimpleVector*>(nnz_vec->children[it]->vec);
      nSingletonRows += initSingletonRowsBlock(int(it), nnz_block);
   }
   presData.setBlocks(nChildren + 1, presData.getNumberSR());

   // todo: linking block nnz_vec->vecl these are not added to the overall system.. todo
   //blocks[nChildren+2] = singletonRows.size();
   /* linking rows */
   if( hasLinking(system_type) && myRank == 0 )
   {
      int n_singletons_link = 0;
      nnz_block = dynamic_cast<SimpleVector*>(nnz_vec->vecl);

      for( int i = 0; i < nnz_block->n; i++ )
         if( nnz_block->elements()[i] == 1.0 )
            n_singletons_link++;
      PIPSdebugMessage(
            "There are %d singleton rows among the linking constraints of A. \n",
            n_singletons_link);
   }
   return nSingletonRows;
}

int StochPresolverSingletonRows::initSingletonRowsBlock(int it, SimpleVector const * nnzRowSimple)
{
   int nSingletonRows = 0;

   presData.setBlocks(it+1, presData.getNumberSR());
   double* nnzRow = nnzRowSimple->elements();

   for( int i = 0; i < nnzRowSimple->n; i++ )
      if( nnzRow[i] == 1.0 )
      {
         presData.addSingletonRow(i);
         nSingletonRows++;
      }
   return nSingletonRows;
}

/** Goes through the singleton rows in the equality system A. For those fixing variables in
 * the blocks B,D,Fi,Gi (the blocks Bmat and Blmat of both A and C), the fixation and updating
 * of the columns is done. The fixed variables in one of the Amat blocks are stored in the
 * member variable colAdaptParent. Updating the blocks A,C,F0,G0 using colAdaptParent happens
 * in updateLinkingVarsBlocks() which should be called after this method.
 * Returns the number of newly found singleton rows (equality/inequality system) during adaption of B,D,Fi,Gi.
 */
void StochPresolverSingletonRows::doSingletonRows(int& n_sing_sys,
      int& n_sing_other_sys, SystemType system_type)
{
   n_sing_sys = 0;
   StochGenMatrix& matrix =
         (system_type == EQUALITY_SYSTEM) ?
               dynamic_cast<StochGenMatrix&>(*(presProb->A)) :
               dynamic_cast<StochGenMatrix&>(*(presProb->C));

   updateCPForSingletonRow(-1, system_type);
   procSingletonRowRoot(matrix, system_type);

   assert(nChildren == (int)matrix.children.size());
   for( int it = 0; it < nChildren; it++ )
   {
      // dummy child?
      if( updateCPForSingletonRow(it, system_type) )
      { // main part for each child: go through A and B and adapt F, D and G
         procSingletonRowChild(it, n_sing_sys, n_sing_other_sys, system_type);
      }
   }

   // Update nRowLink and lhs/rhs (Linking part) of both systems:
   updateRhsNRowLink();
   if( system_type == INEQUALITY_SYSTEM ) // todo why?
      combineNewBoundsParent();
   if( !presData.combineColAdaptParent() )
      abortInfeasible(MPI_COMM_WORLD );
}

void StochPresolverSingletonRows::doSingletonLinkRows(int& newSREq, int& newSRIneq)
{
   // todo what if !hasLinking(EQUALITY_SYSTEM), but hasLinking(INEQUALITY_SYSTEM)???
   if( hasLinking(EQUALITY_SYSTEM))
   {
      setCurrentPointersToNull();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      presData.resetRedCounters();
      resetEqRhsAdaptionsLink();
      if( hasLinking(INEQUALITY_SYSTEM) )
         resetIneqRhsAdaptionsLink();

      for( size_t it = 0; it < presData.nRowElemsA->children.size(); it++)
      {
         if( !nodeIsDummy( it, EQUALITY_SYSTEM))
         {
            setCPBlmatsChild(presProb->A, (int)it);
            // set pointers ixlow etc and redColChild
            setCPColumnChild((int) it);
            setCPRowLinkEquality();    // set currEqRhs
            currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec);

            std::vector<COLUMNTOADAPT> varsToBeFixed;

            for( int rowIdx = 0; rowIdx < currNnzRow->n; rowIdx++ )
            {
               if( currNnzRow->elements()[rowIdx] == 1.0 && (currBlmat->rowptr[rowIdx].start != currBlmat->rowptr[rowIdx].end) )
               {  // singleton entry is in this child:
                  assert( currBlmat->rowptr[rowIdx].start +1 == currBlmat->rowptr[rowIdx].end );

                  // collect all variables that can be fixed in varsToBeFixed:
                  int colIdx = -1;
                  double aik = 0.0;
                  getValuesForSR(*currBlmat, rowIdx, colIdx, aik);
                  assert( colIdx >= 0 && colIdx < currgChild->n );
                  assert(!PIPSisZero(aik));

                  // todo
                  //if( aik < tolerance3 )
                  // continue;

                  const double val = currEqRhsLink->elements()[rowIdx] / aik;

                  if( (currIxlowChild->elements()[colIdx] != 0.0 && PIPSisLT(val, currxlowChild->elements()[colIdx]) )
                        || (currIxuppChild->elements()[colIdx] != 0.0 && PIPSisLT(currxuppChild->elements()[colIdx], val) ) )
                  {
                     cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
                     abortInfeasible(MPI_COMM_WORLD);
                  }
                  // check if this variable was already fixed by another row:
                  bool uniqueFixVariable = true;
                  for(int i=0; i<(int)varsToBeFixed.size(); i++)
                  {
                     if( varsToBeFixed[i].colIdx == colIdx )
                     {
                        uniqueFixVariable = false;
                        if( !PIPSisEQ(varsToBeFixed[i].val, val) )
                           abortInfeasible(MPI_COMM_WORLD);
                     }
                  }
                  if( uniqueFixVariable )
                  {
                     // adapting objOffset:
                     indivObjOffset += currgChild->elements()[colIdx] * val;

                     // store in colADaptLinkBlock for G, F, B, D:
                     COLUMNTOADAPT colWithVal = {colIdx, val};
                     varsToBeFixed.push_back(colWithVal);
                  }
               }
            }
            // using colAdaptLinkBlock, go through the columns in this child:

            // go through the columns in Bmat, Blmat of the equality:
            updateCPforAdaptFixationsBChild( it, EQUALITY_SYSTEM );
            adaptOtherSystemChildB( EQUALITY_SYSTEM, varsToBeFixed, newSREq );

            // and go through the columns in Bmat, Blmat of the inequality:
            updateCPforAdaptFixationsBChild( it, INEQUALITY_SYSTEM );
            adaptOtherSystemChildB( INEQUALITY_SYSTEM, varsToBeFixed, newSRIneq );
         }
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      }

      // Update nRowLink and lhs/rhs (Linking part) of both systems:
      updateRhsNRowLink();
      presData.resetRedCounters();

      // F_0 block:
      setCurrentPointersToNull();
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      setCPBlmatsRoot(presProb->A);
      setCPColumnRoot();
      currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vecl);
      for( int rowIdx = 0; rowIdx < currNnzRow->n; rowIdx++ )
      {
         if( currNnzRow->elements()[rowIdx] == 1.0 && (currBlmat->rowptr[rowIdx].start != currBlmat->rowptr[rowIdx].end) )
         {
            assert( currBlmat->rowptr[rowIdx].start +1 == currBlmat->rowptr[rowIdx].end );

            removeSingleRowEntryB0( *currBlmat, rowIdx);
         }
      }
      if( !presData.combineColAdaptParent() )
         abortInfeasible(MPI_COMM_WORLD);

      // todo delete...summed up two times! (after function call)
      updateLinkingVarsBlocks(newSREq, newSRIneq);
      presData.resetRedCounters();
   }
}

void StochPresolverSingletonRows::procSingletonRowRoot(StochGenMatrix& stochMatrix, SystemType system_type)
{
   SparseStorageDynamic& B0_mat = stochMatrix.Bmat->getStorageDynamicRef();
   assert( presData.getNumberColAdParent() == 0 );

   for(int i = presData.getBlocks(0); i<presData.getBlocks(1); i++)
   {
      int rowIdx = presData.getSingletonRow(i);
      presData.setSingletonRow(i, -1);  // for debugging purposes

      if( system_type == EQUALITY_SYSTEM)
         removeSingleRowEntryB0(B0_mat, rowIdx);
      else
      {
         SparseStorageDynamic& B0_trans = stochMatrix.Bmat->getStorageDynamicTransposedRef();
         removeSingleRowEntryB0Inequality(B0_mat, B0_trans, rowIdx);
      }
   }
}

/* Processing the singleton rows in child it, more precisely, goes through all singleton rows in Amat and Bmat.
 * Those in Amat are stored in colAdaptParent for later processing.
 * Those in Bmat are removed and stored in colAdaptLinkBlock. Furthermore, the corresponding fixed variables (columns)
 * in Bmat and in Blmat are removed.
 * Using this colAdaptLinkBlock, the variables (columns) are removed from the inequalities Bmat, Blmat as well.
 */
void StochPresolverSingletonRows::procSingletonRowChild(int it, int& n_singleton_sys, int& n_singleton_other_sys, SystemType system_type)
{
   procSingletonRowChildAmat(it, system_type);

   std::vector<COLUMNTOADAPT> colAdaptLinkBlock;
   procSingletonRowChildBmat(it, colAdaptLinkBlock, n_singleton_sys, system_type);

   // using colAdaptLinkBlock, go through the columns in Blmat
   if( hasLinking(system_type) )
      adaptChildBlmat( colAdaptLinkBlock, system_type);

   updateCPforAdaptFixationsBChild( it, system_type );
   adaptOtherSystemChildB( system_type, colAdaptLinkBlock, n_singleton_other_sys );
}

void StochPresolverSingletonRows::procSingletonRowChildAmat(int it, SystemType system_type)
{
   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();

   for(int i = presData.getBlocks(it+1); i<presData.getBlocks(it+2); i++)
   {
      int rowIdx = presData.getSingletonRow(i);
      if( currAmat->rowptr[rowIdx].start +1 == currAmat->rowptr[rowIdx].end )
      {
         presData.setSingletonRow(i, -1);  // for debugging purposes

         // store the column index with fixed value in colAdaptParent and adapt objOffset:
         int colIdx = -1;
         double aik = 0.0;
         getValuesForSR(*currAmat, rowIdx, colIdx, aik);

         if( system_type == EQUALITY_SYSTEM )
         {
            double val = currEqRhs->elements()[rowIdx] / aik;

            if( (ixlow[colIdx] != 0.0 && PIPSisLT(val, xlow[colIdx]) )
                  || (ixupp[colIdx] != 0.0 && PIPSisLT(xupp[colIdx], val) ) )
            {
               cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<", child="<<it<<endl;
               abortInfeasible(MPI_COMM_WORLD);
            }
            storeColValInColAdaptParent(colIdx, val);
         }
         else  // INEQUALITY_SYSTEM
         {
            // test what the new bounds imply: infeasiblity, fixation, tightening, redundancy
            double newxlow = -std::numeric_limits<double>::max();
            double newxupp = std::numeric_limits<double>::max();
            double val = 0.0;

            // calculate the newly found bounds on variable x_k:
            calculateNewBoundsOnVariable(newxlow, newxupp, rowIdx, aik);

            // test if they imply infeasibility
            if( newBoundsImplyInfeasible(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
               abortInfeasible(MPI_COMM_WORLD);

            // test if they imply fixation
            else if( newBoundsFixVariable(val, newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
            {
               //cout<<"New bounds imply fixation of linking variable "<<colIdx<<" to value: "<<val<<endl;
               // as in SR(equality), store them to remove the column later
               storeColValInColAdaptParent(colIdx, val);

               // nnz/red Counters are not touched yet, they will be set later when colAdaptParent is applied.
            }
            else
            {
               // test if new bounds are tightening: add to newBoundsParent
               if( newBoundsTightenOldBounds(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
               {
                  PIPSdebugMessage("New bounds tighten bounds of variable %d \n", colIdx);
                  // store them to adapt the bounds on all processes later
                  storeNewBoundsParent(colIdx, newxlow, newxupp);
               }
               //else   PIPSdebugMessage("New bounds are redundant for variable %d \n", colIdx);

               // set a_ik=0.0, nRow--, redCol++
               clearRow(*currAmat, rowIdx);
               // remove entry a_ik in transposed matrix as well
               removeEntryInDynamicStorage(*currAmatTrans, colIdx, rowIdx, val);
               currNnzRow->elements()[rowIdx]--;
               assert( currNnzRow->elements()[rowIdx] == 0 );
               currRedColParent->elements()[colIdx]++;
            }
         }
      }
   }
}

void StochPresolverSingletonRows::procSingletonRowChildBmat(int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock,
      int& newSR, SystemType system_type)
{
   assert(currBmat != NULL);
   PIPSdebugMessage("procSingletonRowChildBmat for child %d and system_type %d \n", it, system_type);
   for(int i = presData.getBlocks(it+1); i<presData.getBlocks(it+2); i++)
   {
      const int rowIdx = presData.getSingletonRow(i);
      if( rowIdx == -1 )
         continue;   // entry was already in Amat
      else if( currBmat->rowptr[rowIdx].start == currBmat->rowptr[rowIdx].end)
         presData.setSingletonRow(i, -1); // entry was already in a previous singleton row in Bmat, set singletonRow(i) to -1.
      else
      {
         assert( currBmat->rowptr[rowIdx].start +1 == currBmat->rowptr[rowIdx].end );
         presData.setSingletonRow(i, -1);  // for debugging purposes
         removeSingleRowEntryChildBmat(rowIdx, colAdaptLinkBlock, system_type, newSR);
      }
   }
}

void StochPresolverSingletonRows::removeSingleRowEntryChildBmat(int rowIdx,
      std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR)
{
   double* ixlow = currIxlowChild->elements();
   double* ixupp = currIxuppChild->elements();
   double* xlow = currxlowChild->elements();
   double* xupp = currxuppChild->elements();

   int colIdx = -1;
   double aik = 0.0;
   getValuesForSR(*currBmat, rowIdx, colIdx, aik);

   if( system_type == EQUALITY_SYSTEM )
   {
      const double val = currEqRhs->elements()[rowIdx] / aik;

      if( (ixlow[colIdx] != 0.0 && PIPSisLT(val, xlow[colIdx]) )
            || (ixupp[colIdx] != 0.0 && PIPSisLT(xupp[colIdx], val) ) )
      {
         cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
         abortInfeasible(MPI_COMM_WORLD);
      }
      newSR += fixVarInChildBlockAndStore( colIdx, val, system_type, colAdaptLinkBlock);
   }
   else  // INEQUALITY_SYSTEM
   {
      // test what the new bounds imply: infeasiblity, fixation, tightening, redundancy
      double newxlow = -std::numeric_limits<double>::max();
      double newxupp = std::numeric_limits<double>::max();
      double val = 0.0;

      // calculate the newly found bounds on variable x_k:
      calculateNewBoundsOnVariable(newxlow, newxupp, rowIdx, aik);

      // test if they imply infeasibility
      if( newBoundsImplyInfeasible(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
         abortInfeasible(MPI_COMM_WORLD);

      // test if they imply fixation
      else if( newBoundsFixVariable(val, newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
      {
         newSR += fixVarInChildBlockAndStore( colIdx, val, system_type, colAdaptLinkBlock);
      }
      else
      {
         // test if new bounds are tightening:
         if( newBoundsTightenOldBounds(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
         {
            PIPSdebugMessage("New bounds tighten bounds of variable %d \n", colIdx);
            // adapt immediately the variable bounds
            setNewBounds(colIdx, newxlow, newxupp, ixlow, xlow, ixupp, xupp);
         }
         //else   PIPSdebugMessage("New bounds are redundant for variable %d \n", colIdx);

         // set a_ik=0.0, nRow--, nCol--
         clearRow(*currBmat, rowIdx);
         // remove entry a_ik in transposed matrix as well
         removeEntryInDynamicStorage(*currBmatTrans, colIdx, rowIdx, val);
         currNnzRow->elements()[rowIdx]--;
         assert( currNnzRow->elements()[rowIdx] == 0 );
         currNnzColChild->elements()[colIdx]--;
         assert( currNnzColChild->elements()[colIdx] >= 0 );
      }

   }
}

/** Removes the single entry in row rowIdx in Bmat_0.
 * Stores the corresponding column index in colAdaptParent.
 * Abort false if infeasibility is detected.*/
void StochPresolverSingletonRows::removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();

   int colIdx = -1;
   double aik = 0.0;
   getValuesForSR(storage, rowIdx, colIdx, aik);

   const double val = currEqRhs->elements()[rowIdx] / aik;

   if( (ixlow[colIdx] != 0.0 && PIPSisLT(val, xlow[colIdx]) )
         || (ixupp[colIdx] != 0.0 && PIPSisLT(xupp[colIdx], val) ) )
   {
      cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
      abortInfeasible(MPI_COMM_WORLD);
   }
   else
   {
      if( myRank == 0 )
         storeColValInColAdaptParent(colIdx, val);
   }
}

void StochPresolverSingletonRows::removeSingleRowEntryB0Inequality(SparseStorageDynamic& storage,
      SparseStorageDynamic& storageTransposed, int rowIdx)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();

   int colIdx = -1;
   double aik = 0.0;
   getValuesForSR(storage, rowIdx, colIdx, aik);

   double newxlow = -std::numeric_limits<double>::max();
   double newxupp = std::numeric_limits<double>::max();
   double val = 0.0;

   // calculate the newly found bounds on variable x_k:
   calculateNewBoundsOnVariable(newxlow, newxupp, rowIdx, aik);

   // test if they imply infeasibility
   if( newBoundsImplyInfeasible(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
      abortInfeasible(MPI_COMM_WORLD);

   // test if they imply fixation
   else if( newBoundsFixVariable(val, newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
   {
      PIPSdebugMessage("New bounds imply fixation of variable %d to %f \n", colIdx, val);

      // as in SR(equality), store them to remove the column later
      if( myRank == 0 )
         storeColValInColAdaptParent(colIdx, val);

      // in case of fixation, nnz bzw. red Counters are not touched yet because they will be set
      // correctly later, when colAdaptParent is applied.
   }
   else
   {
      // test if new bounds are tightening: add to newBoundsParent
      if( newBoundsTightenOldBounds(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
      {
         PIPSdebugMessage("New bounds tighten bounds of variable  %d \n", colIdx);

         // store them to adapt the bounds on all processes later
         if( myRank == 0 )
            storeNewBoundsParent(colIdx, newxlow, newxupp);
      }
      //else
         //cout<<"New bounds are redundant for variable "<<colIdx<<endl;

      // set a_ik=0.0, nRow=0, nCol--
      clearRow(storage, rowIdx);
      // remove entry a_ik in transposed matrix as well
      removeEntryInDynamicStorage(storageTransposed, colIdx, rowIdx, val);
      currNnzRow->elements()[rowIdx]--;
      currNnzColParent->elements()[colIdx]--;
      assert( currNnzRow->elements()[rowIdx] == 0 );
   }
}

void StochPresolverSingletonRows::calculateNewBoundsOnVariable(double& newxlow, double& newxupp, int rowIdx, double aik) const
{
   if( PIPSisLT(0.0, aik) )
   {
      if( currIclow->elements()[rowIdx] != 0.0 )
         newxlow = currIneqLhs->elements()[rowIdx] / aik;
      if( currIcupp->elements()[rowIdx] != 0.0 )
         newxupp = currIneqRhs->elements()[rowIdx] / aik;
   }
   else
   {
      if( currIcupp->elements()[rowIdx] != 0.0 )
         newxlow = currIneqRhs->elements()[rowIdx] / aik;
      if( currIclow->elements()[rowIdx] != 0.0 )
         newxupp = currIneqLhs->elements()[rowIdx] / aik;
   }
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
   for(int i=0; i<getNumberNewBoundsParent(); i++)
   {
      XBOUNDS newbounds = getNewBoundsParent(i);
      setNewBounds(newbounds.colIdx, newbounds.newxlow, newbounds.newxupp, ixlow, xlow, ixupp, xupp);
   }
   clearNewBoundsParent();
}

void StochPresolverSingletonRows::getValuesForSR(SparseStorageDynamic const & storage, int rowIdx, int& colIdx, double& aik) const
{
   const int indexK = storage.rowptr[rowIdx].start;
   colIdx = storage.jcolM[indexK];
   aik = storage.M[indexK];

   assert(storage.rowptr[rowIdx].start +1 == storage.rowptr[rowIdx].end);
   assert(aik != 0.0);
}

/** Update the current pointers for the singleton row routine.
 * If it==-1, we are at parent block. Else, et child[it].
 * Return false if child[it] is a dummy child. */
bool StochPresolverSingletonRows::updateCPForSingletonRow(int it, SystemType system_type)
{
   setCurrentPointersToNull();

   setCPColumnRoot();
   currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);
   for(int i=0; i<currgParent->n; i++)
      assert( std::isfinite(currgParent->elements()[i]) );

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
      for(int i=0; i<currgChild->n; i++)
         assert( isfinite(currgChild->elements()[i]) );
      currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
   }
   return true;
}



