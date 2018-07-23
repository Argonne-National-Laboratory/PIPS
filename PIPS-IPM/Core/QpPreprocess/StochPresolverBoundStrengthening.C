/*
 * StochPresolverBoundStrengthening.C
 *
 *  Created on: 28.05.2018
 *      Author: Svenja Uslu
 */

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


bool StochPresolverBoundStrengthening::applyPresolving(int& nelims)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   if( myRank == 0 ) cout<<"Start Bound Strengthening Presolving..."<<endl;

   clearNewBoundsParent();

   // root:
   setCPforBounds(presProb->A, -1, EQUALITY_SYSTEM);
   if( !doBoundStrengthParent( EQUALITY_SYSTEM ) )
      return false;
   setCPforBounds(presProb->C, -1, INEQUALITY_SYSTEM);
   if( !doBoundStrengthParent( INEQUALITY_SYSTEM ) )
      return false;

   // children:
   for( size_t child_it = 0; (int)child_it < nChildren; child_it++)
   {
      // dummy child?
      if( setCPforBounds(presProb->A, (int)child_it, EQUALITY_SYSTEM) )
         if( !doBoundStrengthChild(EQUALITY_SYSTEM) )
            return false;
      if( setCPforBounds(presProb->C, (int)child_it, INEQUALITY_SYSTEM) )
         if( !doBoundStrengthChild(INEQUALITY_SYSTEM) )
            return false;
   }

   // linking rows:
   // todo

   // combine the bounds of linking-variables:
   if( !combineNewBoundsParent() ) return false;
   // update the the bounds of linking-variables:
   strengthenLinkingVarsBounds();

   if( myRank == 0 ) cout<<"Finished Bound Strengthening Presolving."<<endl;

   return true;
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
   if(it >= 0)
   {
      if( !setCPAmatBmat( matrixHandle, it, system_type) )  //currAmat, currBmat
         return false;
      setCPColumnChild(it);   //currxlowChild etc.
      if( system_type == INEQUALITY_SYSTEM )
         setCPRowChildIneqOnlyLhsRhs(it); //currIneqRhs etc.
      else
         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
   }
   else{
      setCPAmatBmat( matrixHandle, it, system_type);  //currAmat
      // set rhs/lhs vectors:
      if( system_type == INEQUALITY_SYSTEM )
         setCPRowRootIneqOnlyLhsRhs();
      else
         currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
   }
   setCPColumnRoot();   //currxlowParent etc.

   return true;
}

/**
 * Compute the minimum and maximum activity of the row rowIdx in matrix. If colIdx!=-1, then this entry
 * is excluded in the computation of the activities.
 */
void StochPresolverBoundStrengthening::computeActivityBlockwise( SparseStorageDynamic& matrix, int rowIdx, int colIdx,
      double& infRow, double& supRow,
      SimpleVector& xlow, SimpleVector& ixlow, SimpleVector& xupp, SimpleVector& ixupp)
{
   // todo: possibly add two bools: if infty, indicate if at least two bounds were infty
   assert( rowIdx >= 0 && rowIdx < matrix.m );
   const int n = matrix.n;
   assert( colIdx >= -1 && colIdx < n );
   assert( xlow.n == n && ixlow.n == n && xupp.n == n && ixupp.n == n );

   for( int j=matrix.rowptr[rowIdx].start; j<matrix.rowptr[rowIdx].end; j++)
   {
      const int col = matrix.jcolM[j];
      const double entry = matrix.M[j];
      if( col == colIdx )
         continue;
      if( entry > 0)
      {
         // add entry * lower_bound to infRow
         if( ixlow.elements()[col] != 0.0)
            infRow += entry * xlow.elements()[col];
         else
            infRow = -std::numeric_limits<double>::max();
         // add entry * upper_bound to supRow
         if( ixupp.elements()[col] != 0.0 )
            supRow += entry * xupp.elements()[col];
         else
            supRow = std::numeric_limits<double>::max();
      }
      else
      {
         // add entry * upper_bound to infRow
         if( ixupp.elements()[col] != 0.0 )
            infRow += entry * xupp.elements()[col];
         else
            infRow = -std::numeric_limits<double>::max();
         // add entry * lower_bound to supRow
         if( ixlow.elements()[col] != 0.0 )
            supRow += entry * xlow.elements()[col];
         else
            supRow = std::numeric_limits<double>::max();
      }
      if( supRow == std::numeric_limits<double>::max() && infRow == -std::numeric_limits<double>::max() )
         return;
   }
}

bool StochPresolverBoundStrengthening::doBoundStrengthParent(SystemType system_type)
{
   for(int rowIdx=0; rowIdx<currAmat->m; rowIdx++)
      if( !strenghtenBoundsInBlock( *currAmat, false, rowIdx, 0.0, 0.0, system_type))
         return false;

   return true;
}
bool StochPresolverBoundStrengthening::doBoundStrengthChild(SystemType system_type)
{
   for(int i=0; i<currAmat->m; i++) // i ~ rowIndex
   {
      // First, compute minActivity and maxActivity of the whole row (per block):
      double partMinActivityA = 0.0, partMaxActivityA = 0.0, partMinActivityB = 0.0, partMaxActivityB = 0.0;
      computeActivityBlockwise(*currAmat, i, -1, partMinActivityA, partMaxActivityA,
            *currxlowParent, *currIxlowParent, *currxuppParent, *currIxuppParent);
      computeActivityBlockwise(*currBmat, i, -1, partMinActivityB, partMaxActivityB,
            *currxlowChild, *currIxlowChild, *currxuppChild, *currIxuppChild);

      if( !strenghtenBoundsInBlock( *currAmat, false, i, partMinActivityB, partMaxActivityB, system_type))
         return false;
      if( !strenghtenBoundsInBlock( *currBmat, true, i, partMinActivityA, partMaxActivityA, system_type))
         return false;
   }
   return true;
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
bool StochPresolverBoundStrengthening::strenghtenBoundsInBlock( SparseStorageDynamic& matrix, bool childBlock,
      int rowIdx, double partMinActivity, double partMaxActivity, SystemType system_type)
{
   assert( rowIdx >= 0 && rowIdx < matrix.m );
   if( partMinActivity == -std::numeric_limits<double>::max() && partMaxActivity == std::numeric_limits<double>::max())
      return true;

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
         if( a_ik > 0 )
            newBoundUpp = computeNewBound(true, lis, a_ik, rowIdx, system_type);
         else
            newBoundLow = computeNewBound(true, lis, a_ik, rowIdx, system_type);
      }
      if( uis < std::numeric_limits<double>::max() )
      {
         if( a_ik > 0 )
            newBoundLow = computeNewBound(false, uis, a_ik, rowIdx, system_type);
         else
            newBoundUpp = computeNewBound(false, uis, a_ik, rowIdx, system_type);
      }

      // tighten the bounds:
      if( childBlock )
      {
         if( newBoundsImplyInfeasible(newBoundLow, newBoundUpp, colIdx,
               currIxlowChild->elements(), currIxuppChild->elements(), currxlowChild->elements(), currxuppChild->elements()) )
            return false;
         // todo test if bounds imply fixation
         setNewBoundsIfTighter(colIdx, newBoundLow, newBoundUpp,
            *currIxlowChild, *currxlowChild, *currIxuppChild, *currxuppChild);
      }
      else
      {
         double varvalue = 0.0;
         // test if new bounds imply infeasible ?
         if( newBoundsImplyInfeasible(newBoundLow, newBoundUpp, colIdx,
               currIxlowParent->elements(), currIxuppParent->elements(), currxlowParent->elements(), currxuppParent->elements()) )
            return false;
         // todo test if they imply fixation ?
         /*else if( newBoundsFixVariable(varvalue, newBoundLow, newBoundUpp, colIdx,
               currIxlowParent->elements(), currIxuppParent->elements(), currxlowParent->elements(), currxuppParent->elements()) )
         {
            //cout<<"New bounds imply fixation of variable "<<colIdx<<" of child "<<it<<" to value: "<<val<<endl;
            // as in SR(equality), store them to remove the column later
            if( !storeColValInColAdaptParentAndAdaptOffset(colIdx, varvalue, currgParent->elements()) )
               return false;
            // nnz/red Counters are not touched yet, they will be set later when colAdaptParent is applied

         }*/
         // store the bounds in newBoundsParent for all processes:
         if( checkNewBoundTightens(true, colIdx, newBoundUpp, *currIxuppParent, *currxuppParent)
               || checkNewBoundTightens(false, colIdx, newBoundLow, *currIxlowParent, *currxlowParent) )
            storeNewBoundsParent(colIdx, newBoundLow, newBoundUpp);
      }

   }
   return true;
}

/**
 * Compares and sets the new lower and upper bounds, provided that:
 * checkNewBoundTightens() returns true.
 */
void StochPresolverBoundStrengthening::setNewBoundsIfTighter(int index, double new_low, double new_upp,
      SimpleVector& ilow, SimpleVector& low, SimpleVector& iupp, SimpleVector& upp)
{
   const int n = ilow.n;
   assert( index >= 0 && index < n );
   assert( low.n == n && iupp.n == n && upp.n == n );

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

