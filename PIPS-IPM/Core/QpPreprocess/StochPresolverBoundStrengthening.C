/*
 * StochPresolverBoundStrengthening.C
 *
 *  Created on: 28.05.2018
 *      Author: Svenja Uslu
 */

//#define PIPS_DEBUG
#include "StochPresolverBoundStrengthening.h"

StochPresolverBoundStrengthening::StochPresolverBoundStrengthening(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverBoundStrengthening::~StochPresolverBoundStrengthening()
{
 // todo
}


void StochPresolverBoundStrengthening::applyPresolving()
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- Before Bound Strengthening Presolving:"<<endl;
   countRowsCols();
#endif

   if( myRank == 0 ) cout<<"Start Bound Strengthening Presolving..."<<endl;

   indivObjOffset = 0.0;
   clearNewBoundsParent();

   // root:
   setCPforBounds(presProb->A, -1, EQUALITY_SYSTEM);
   doBoundStrengthParent( EQUALITY_SYSTEM );
   setCPforBounds(presProb->C, -1, INEQUALITY_SYSTEM);
   doBoundStrengthParent( INEQUALITY_SYSTEM );

   // children:
   for( size_t child_it = 0; (int)child_it < nChildren; child_it++)
   {
      // dummy child?
      if( setCPforBounds(presProb->A, (int)child_it, EQUALITY_SYSTEM) )
         doBoundStrengthChild((int)child_it, EQUALITY_SYSTEM);
      if( setCPforBounds(presProb->C, (int)child_it, INEQUALITY_SYSTEM) )
         doBoundStrengthChild((int)child_it, INEQUALITY_SYSTEM);
   }
   // Update nRowLink and lhs/rhs (Linking part) of both systems:
   updateRhsNRowLink();

   // linking rows:
   // todo

   // combine the bounds of linking-variables:
   combineNewBoundsParent();
   // update the the bounds of linking-variables:
   strengthenLinkingVarsBounds();

   // update the linking variable blocks (A,C,F,G) with the fixations found:
   if( !presData.combineColAdaptParent())
      abortInfeasible(MPI_COMM_WORLD);
   int newSREq = 0, newSRIneq = 0;
   updateLinkingVarsBlocks(newSREq, newSRIneq);

   // Sum up individual objOffset and then add it to the global objOffset:
   sumIndivObjOffset();
   presData.addObjOffset(indivObjOffset);
   if( myRank == 0 )
      cout<<"Global objOffset is now: "<<presData.getObjOffset()<<endl;

   if( myRank == 0 ) cout<<"Finished Bound Strengthening Presolving."<<endl;
#ifndef NDEBUG
   if( myRank == 0 )
      cout<<"--- After Bound Strengthening Presolving:"<<endl;
   countRowsCols();
#endif
}

/**
 * Set current pointers to the necessary elements: For all cases: currAmat and currxlowParent &co.
 * For equality system currEqRhs and for inequality system currIneqRhs &co.
 * For a child it: currBmat and currxlowChild &co.
 */
bool StochPresolverBoundStrengthening::setCPforBounds(GenMatrixHandle matrixHandle, int it, SystemType system_type)
{
   assert( it >= -1 && it < nChildren );
   setCurrentPointersToNull();
   currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);
   if(it >= 0)
   {
      if( !setCPAmatsChild( matrixHandle, it, system_type) )  // currAmat(Trans)
         return false;
      setCPBmatsChild( matrixHandle, it, system_type );  // currBmat(Trans)
      setCPColumnChild(it);   //currxlowChild etc.
      if( system_type == INEQUALITY_SYSTEM )
         setCPRowChildInequality(it);
      else
         setCPRowChildEquality(it);
      currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec);
      currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
   }
   else
   {
      setCPAmatsRoot( matrixHandle );  //currAmat(Trans)
      // set rhs/lhs vectors:
      if( system_type == INEQUALITY_SYSTEM )
         setCPRowRootIneqOnlyLhsRhs();
      else
         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
   }
   setCPColumnRoot();   //currxlowParent etc.

   return true;
}

void StochPresolverBoundStrengthening::doBoundStrengthParent(SystemType system_type)
{
   for(int rowIdx=0; rowIdx<currAmat->m; rowIdx++)
      strenghtenBoundsInBlock( *currAmat, false, rowIdx, 0.0, 0.0, system_type, true, NULL);
}

void StochPresolverBoundStrengthening::doBoundStrengthChild(int it, SystemType system_type)
{
   std::vector<COLUMNTOADAPT> colAdaptLinkBlock;
   for(int i=0; i<currAmat->m; i++) // i ~ rowIndex
   {
      // First, compute minActivity and maxActivity of the whole row (per block):
      double partMinActivityA = 0.0, partMaxActivityA = 0.0, partMinActivityB = 0.0, partMaxActivityB = 0.0;
      computeActivityBlockwise(*currAmat, i, -1, partMinActivityA, partMaxActivityA,
            *currxlowParent, *currIxlowParent, *currxuppParent, *currIxuppParent);
      computeActivityBlockwise(*currBmat, i, -1, partMinActivityB, partMaxActivityB,
            *currxlowChild, *currIxlowChild, *currxuppChild, *currIxuppChild);

      strenghtenBoundsInBlock( *currAmat, false, i, partMinActivityB, partMaxActivityB, system_type, false, NULL);
      strenghtenBoundsInBlock( *currBmat, true, i, partMinActivityA, partMaxActivityA, system_type, false,
            &colAdaptLinkBlock);
   }
   // use information in colAdaptLinkBlock to adapt the other blocks:
   // First, go through the columns in Blmat (of the same system)
   if( hasLinking(system_type) )
   {
      updateCPforAdaptFixationsBChild( it, system_type );
      adaptChildBlmat( colAdaptLinkBlock, system_type);
   }
   // Second, go through the columns in Bmat and Blmat of the other system:
   SystemType other_system = EQUALITY_SYSTEM;
   if( system_type == EQUALITY_SYSTEM )
      other_system = INEQUALITY_SYSTEM;
   updateCPforAdaptFixationsBChild( it, other_system );

   int tmp_newSR = 0;
   adaptOtherSystemChildB( other_system, colAdaptLinkBlock, tmp_newSR );
}

/** Computes the new bound (bA - activity)/matrixEntry.
 * If systemType==EQUALITY_SYSTEM, then bA is the rhs currEqRhs.
 * Else, the boolean rhs defines if the lhs or the rhs of the inequality should be used as bA in the computation:
 * If rhs is true, then the right hand side cupp is used. If false, then the lhs clow is used.
 */
double StochPresolverBoundStrengthening::computeNewBound(bool rhs, double activity, double matrixEntry, int rowIdx, SystemType system_type)
{
   assert( matrixEntry != 0.0 );
   assert( rowIdx >= 0 );
   if( system_type == EQUALITY_SYSTEM )
   {
      assert( rowIdx >= 0 && rowIdx<currEqRhs->n);
      return (currEqRhs->elements()[rowIdx] - activity) / matrixEntry;
   }

   // if INEQUALITY_SYSTEM:
   assert( rowIdx >= 0 && rowIdx<currIneqRhs->n);
   if( rhs )
   {
      if( currIcupp->elements()[rowIdx] != 0.0 )
         return ( currIneqRhs->elements()[rowIdx] - activity ) / matrixEntry;
      else
         return std::numeric_limits<double>::max() * matrixEntry; // *matrixEntry necessary for the correct sign
   }
   else  // lhs is used
   {
      if( currIclow->elements()[rowIdx] != 0.0 )
         return ( currIneqLhs->elements()[rowIdx] - activity ) / matrixEntry;
      else
         return -std::numeric_limits<double>::max() * matrixEntry;
   }

}

/**
 * Strengthen the variable bounds in the block matrix. If childBlock==true, then a block B_i or D_i is considered.
 * partMinActivity and partMaxActivity represent the partial row activity of the respective other block.
 */
void StochPresolverBoundStrengthening::strenghtenBoundsInBlock( SparseStorageDynamic& matrix, bool childBlock,
      int rowIdx, double partMinActivity, double partMaxActivity, SystemType system_type, bool atRoot,
      std::vector<COLUMNTOADAPT>* colAdaptLinkBlock)
{
   assert( rowIdx >= 0 && rowIdx < matrix.m );

   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   if( partMinActivity == -std::numeric_limits<double>::max() && partMaxActivity == std::numeric_limits<double>::max())
      return;

   for( int j=matrix.rowptr[rowIdx].start; j<matrix.rowptr[rowIdx].end; j++)
   {
      const int colIdx = matrix.jcolM[j];
      double lis = partMinActivity, uis = partMaxActivity;

      // Compute remaining activity of the row:
      if( childBlock )
         computeActivityBlockwise(matrix, rowIdx, colIdx, lis, uis,
            *currxlowChild, *currIxlowChild, *currxuppChild, *currIxuppChild);
      else
         computeActivityBlockwise(matrix, rowIdx, colIdx, lis, uis,
            *currxlowParent, *currIxlowParent, *currxuppParent, *currIxuppParent);

      if( lis == -std::numeric_limits<double>::max() && uis == std::numeric_limits<double>::max())
         continue;

      // compute the possible new bounds on variable x_colIdx:
      const double a_ik = matrix.M[j];
      double newBoundLow = -std::numeric_limits<double>::max();
      double newBoundUpp = std::numeric_limits<double>::max();
      if( lis > -std::numeric_limits<double>::max() )
      {
         if( PIPSisLT(0.0, a_ik) )
            newBoundUpp = computeNewBound(true, lis, a_ik, rowIdx, system_type);
         else
            newBoundLow = computeNewBound(true, lis, a_ik, rowIdx, system_type);
      }
      if( uis < std::numeric_limits<double>::max() )
      {
         if( PIPSisLT(0.0, a_ik) )
            newBoundLow = computeNewBound(false, uis, a_ik, rowIdx, system_type);
         else
            newBoundUpp = computeNewBound(false, uis, a_ik, rowIdx, system_type);
      }

      //PIPSdebugMessage("Row %d. At variable %d, newly found bounds are [%f, %f]. lis=%f, uis=%f \n", rowIdx, colIdx, newBoundLow, newBoundUpp, lis, uis);
      // tighten the bounds:
      if( childBlock )
      {
         assert( colAdaptLinkBlock );
         if( newBoundsImplyInfeasible(newBoundLow, newBoundUpp, colIdx,
               currIxlowChild->elements(), currIxuppChild->elements(), currxlowChild->elements(), currxuppChild->elements()) )
            abortInfeasible(MPI_COMM_WORLD);
         // test if bounds imply fixation:
         double varvalue = 0.0;
         if( newBoundsFixVariable(varvalue, newBoundLow, newBoundUpp, colIdx,currIxlowChild->elements(),
               currIxuppChild->elements(), currxlowChild->elements(), currxuppChild->elements()) )
         {
            PIPSdebugMessage("New bounds imply fixation of variable %d to value=%f. original bounds:"
                  " [%f, %f] (if existent). Info from row %d \n", colIdx, varvalue, currxlowParent->elements()[colIdx],
                  currxuppParent->elements()[colIdx], rowIdx);
            fixVarInChildBlockAndStore( colIdx, varvalue, system_type, *colAdaptLinkBlock);
         }
         else
            setNewBoundsIfTighter(colIdx, newBoundLow, newBoundUpp,
            *currIxlowChild, *currxlowChild, *currIxuppChild, *currxuppChild);
      }
      else  // Linking Variable Block (Amat)
      {
         // test if new bounds imply infeasible ?
         if( newBoundsImplyInfeasible(newBoundLow, newBoundUpp, colIdx,
               currIxlowParent->elements(), currIxuppParent->elements(), currxlowParent->elements(), currxuppParent->elements()) )
            abortInfeasible(MPI_COMM_WORLD);
         // test if they imply fixation ?
         double varvalue = 0.0;
         if( newBoundsFixVariable(varvalue, newBoundLow, newBoundUpp, colIdx,
               currIxlowParent->elements(), currIxuppParent->elements(), currxlowParent->elements(), currxuppParent->elements()) )
         {
            PIPSdebugMessage("New bounds imply fixation of variable %d to value=%f. original bounds: [%f, %f] (if existent). \n", colIdx, varvalue, currxlowParent->elements()[colIdx], currxuppParent->elements()[colIdx]);
            // store the fixation to remove the column later after communicating
            if( !atRoot || myRank==0 )
              storeColValInColAdaptParent(colIdx, varvalue);
            // tighten the bounds to varvalue on this processor:
            setNewBound(colIdx, varvalue, currxlowParent, currIxlowParent );
            setNewBound(colIdx, varvalue, currxuppParent, currIxuppParent );

            // note : nnz/red Counters are not touched yet, they will be set later when colAdaptParent is applied

         }
         else  // no fixation, but maybe bound strengthening:
         {
            if( atRoot )
               setNewBoundsIfTighter(colIdx, newBoundLow, newBoundUpp,
                     *currIxlowParent, *currxlowParent, *currIxuppParent, *currxuppParent);
            else if( checkNewBoundTightens(true, colIdx, newBoundUpp, *currIxuppParent, *currxuppParent)
                  || checkNewBoundTightens(false, colIdx, newBoundLow, *currIxlowParent, *currxlowParent) )
               storeNewBoundsParent(colIdx, newBoundLow, newBoundUpp);  // store the bounds in newBoundsParent for all processes
         }
      }
   }
}

/**
 * Compares and sets the new lower and upper bounds, provided that:
 * checkNewBoundTightens() returns true.
 */
void StochPresolverBoundStrengthening::setNewBoundsIfTighter(int index, double new_low, double new_upp,
      SimpleVector& ilow, SimpleVector& low, SimpleVector& iupp, SimpleVector& upp)
{
   assert( index >= 0 && index < ilow.n );
   assert( low.n == ilow.n && iupp.n == ilow.n && upp.n == ilow.n );

   if( checkNewBoundTightens(false, index, new_low, ilow, low) )
   {
      ilow.elements()[index] = 1.0;
      low.elements()[index] = new_low;
   }
   if( checkNewBoundTightens(true, index, new_upp, iupp, upp) )
   {
      iupp.elements()[index] = 1.0;
      upp.elements()[index] = new_upp;
   }
}

/** Check if the new bound (upper or lower bound depending on the bool uppLow) tightens the
 * existing bound. Returns true if:
 * - the new bounds tighten the old bounds,
 * - the absolute value of new bounds does not exceed limit2 = 1.0e8,
 * - the change in the bounds is at least limit1*epsilon = 1.0e3*1.0e-6.
 */
bool StochPresolverBoundStrengthening::checkNewBoundTightens(bool uppLow, int colIdx, double newBound,
      SimpleVector& ixbound, SimpleVector& xbound ) const
{
   if( fabs(newBound) >= limit2 )
      return false;
   if(uppLow)   // upper bound
   {
      if( (ixbound.elements()[colIdx] != 0.0 && xbound.elements()[colIdx] - newBound >= limit1 * feastol )
            || (ixbound.elements()[colIdx] == 0.0) )
         return true;
   }
   else  // lower bound
   {
      if( (ixbound.elements()[colIdx] != 0.0 && newBound - xbound.elements()[colIdx] >= limit1 * feastol )
            || (ixbound.elements()[colIdx] == 0.0) )
         return true;
   }
   return false;
}

/** Should be called after combineNewBoundsParent() or another method that stores
 * information to update in newBoundsParent.
 * Updates the bounds on the linking variables, if they satisfy the conditions in checkNewBoundTightens().
 */
void StochPresolverBoundStrengthening::strengthenLinkingVarsBounds()
{
   setCPColumnRoot();

   // apply updated newBoundsParent to the variable bounds.
   for(int i=0; i<getNumberNewBoundsParent(); i++)
   {
      XBOUNDS newbounds = getNewBoundsParent(i);
      setNewBoundsIfTighter(newbounds.colIdx, newbounds.newxlow, newbounds.newxupp, *currIxlowParent, *currxlowParent, *currIxuppParent, *currxuppParent);
   }
   clearNewBoundsParent();
}

