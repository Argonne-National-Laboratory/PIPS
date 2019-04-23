/*
 * StochPresolverBase.cpp
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

//#define PIPS_DEBUG
#include "StochPresolverBase.h"
#include "DoubleMatrixTypes.h"
#include "SmartPointer.h"
#include "pipsdef.h"

#include <algorithm>
#include <cassert>
#include <limits>
#include <cmath> // std::isfinite


StochPresolverBase::StochPresolverBase(PresolveData& presData) :
      presData(presData)
{
   presProb = presData.presProb;

   setCurrentPointersToNull();

   if( hasLinking(EQUALITY_SYSTEM) )
   {
      currEqRhsAdaptionsLink = new double[presData.redRowA->vecl->n]();
   }
   else
   {
      currEqRhsAdaptionsLink = NULL;
   }

   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      currInEqRhsAdaptionsLink = new double[presData.redRowC->vecl->n]();
      currInEqLhsAdaptionsLink = new double[presData.redRowC->vecl->n]();
   }
   else
   {
      currInEqRhsAdaptionsLink = NULL;
      currInEqLhsAdaptionsLink = NULL;
   }

   localNelims = 0;
   nChildren = presData.getNChildren();

   indivObjOffset = 0.0;
}

StochPresolverBase::~StochPresolverBase()
{
   delete[] currEqRhsAdaptionsLink;
   delete[] currInEqRhsAdaptionsLink;
   delete[] currInEqLhsAdaptionsLink;
}

/**
 * In the dynamic sparse storage, swap entry (rowidx, jcolM[indexK]) with the last entry in this row.
 * Decrement rowptr[rowidx].end by one. Decrement rowend and indexK. Increment currRedRow[rowidx] and redCol[colIdx].
 */
void StochPresolverBase::updateAndSwap(SparseStorageDynamic* storage, int rowidx,
      int& indexK, int& rowEnd, double* redCol, int& nelims, bool linking)
{
   if(linking)
      currRedRowLink->elements()[rowidx]++;
   else
      currRedRow->elements()[rowidx]++;

   redCol[storage->jcolM[indexK]]++;

   std::swap(storage->M[indexK], storage->M[rowEnd - 1]);
   std::swap(storage->jcolM[indexK], storage->jcolM[rowEnd - 1]);
   storage->rowptr[rowidx].end--;
   rowEnd = storage->rowptr[rowidx].end;
   indexK--;

   nelims++;
}

/**
 * Update the linking parts of nnzRow-Vectors using the vectors presData.redRowA->vecl and presData.redRowc->vecl.
 * Update the linking parts of the rhs (and lhs)-vectors using the vectors currEqRhsAdaptionsLink,
 * currInEqRhsAdaptionsLink and currInEqLhsAdaptionsLink.
 */
// todo call me also from clean-up!
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
         const int message_sizeA = dynamic_cast<SimpleVector*>(presData.redRowA->vecl)->length();
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
         const int message_sizeC = dynamic_cast<SimpleVector*>(presData.redRowC->vecl)->length();
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

void StochPresolverBase::updateNnzFromReductions(SystemType system_type)
{

   // todo reductions in root nodes must be synced
   StochVectorHandle red_vector;
   StochVectorHandle nnz_vector;

   /* update column non-zeros */
   red_vector = presData.redCol;
   nnz_vector = presData.nColElems;

   // todo
   if(nnz_vector->vec)
      updateNnzUsingReductions(nnz_vector->vec, red_vector->vec);
   for(size_t node = 0; node < red_vector->children.size(); ++node)
      updateNnzUsingReductions(dynamic_cast<SimpleVector*>(nnz_vector->children[node]->vec),
            dynamic_cast<SimpleVector*>(red_vector->children[node]->vec));

   /* update row non-zeros for equality system */
   if(system_type == EQUALITY_SYSTEM)
   {
      red_vector = presData.redRowA;
      nnz_vector = presData.nRowElemsA;

      updateNnzUsingReductions(nnz_vector, red_vector, EQUALITY_SYSTEM);
   }

   /* update row non-zeros for inequality system */
   if(system_type == INEQUALITY_SYSTEM)
   {
      red_vector = presData.redRowC;
      nnz_vector = presData.nRowElemsC;

      updateNnzUsingReductions(nnz_vector, red_vector, INEQUALITY_SYSTEM);
   }

   presData.resetRedCounters();
}

void StochPresolverBase::updateNnzUsingReductions( StochVectorHandle nnz_vector, StochVectorHandle red_vector, SystemType system_type) const
{
   /* root */
   if(nnz_vector->vec)
      updateNnzUsingReductions(nnz_vector->vec, red_vector->vec);

   if(nnz_vector->vecl)
      updateNnzUsingReductions(nnz_vector->vecl, red_vector->vecl);

   /* children */
   for(int node = 0; node < nChildren; ++node)
   {
      if(!nodeIsDummy(node, system_type)){
         updateNnzUsingReductions(dynamic_cast<SimpleVector*>(nnz_vector->children[node]->vec),
               dynamic_cast<SimpleVector*>(red_vector->children[node]->vec));

         assert(nnz_vector->children[node]->vecl == NULL);
         assert(red_vector->children[node]->vecl == NULL);
      }
   }
}


/** Update the nnzVector by subtracting the reductions vector. */
void StochPresolverBase::updateNnzUsingReductions( OoqpVector* nnzVector, OoqpVector* redVector) const
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

/**
 * Update the vector presData.nColElems->vec (nnzColParent vector).
 * In the distributed case, MPI_allreduce is used to sum up the different
 * redColParent vectors (presData.redCol->vec). Then the reductions are
 * subtracted from the nnz-vector.
 */
void StochPresolverBase::updateNnzColParent(MPI_Comm comm)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);
   if( world_size > 1)
   {
      double* redColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec)->elements();
      int message_size = dynamic_cast<SimpleVector*>(presData.redCol->vec)->length();
      MPI_Allreduce(MPI_IN_PLACE, redColParent, message_size, MPI_DOUBLE, MPI_SUM, comm);
   }
   updateNnzUsingReductions(presData.nColElems->vec, presData.redCol->vec);
}

void StochPresolverBase::allreduceAndUpdate(MPI_Comm comm, SimpleVector& adaptionsVector, SimpleVector& baseVector)
{
   assert(adaptionsVector.n == baseVector.n);

   int world_size;
   MPI_Comm_size(comm, &world_size);
   if( world_size > 1)
   {
      double* adaptionsDouble = adaptionsVector.elements();
      int message_size = adaptionsVector.length();
      MPI_Allreduce(MPI_IN_PLACE, adaptionsDouble, message_size, MPI_DOUBLE, MPI_SUM, comm);
   }
   baseVector.axpy(1.0, adaptionsVector);
}

void StochPresolverBase::updateTransposedSubmatrix(
      SparseStorageDynamic* transStorage,
      std::vector<std::pair<int, int> >& elements) const
{
   for( size_t i = 0; i < elements.size(); ++i)
   {
      std::pair<int,int> entry = elements.at(i);
      const int row_A = entry.first;
      const int row_At = entry.second;

      const int start = transStorage->rowptr[row_At].start;
      const int end = transStorage->rowptr[row_At].end;
      int col_At;

      for( col_At = start; col_At < end; col_At++ )
      {
         if( transStorage->jcolM[col_At] == row_A )
            break;
      }

      std::swap(transStorage->M[col_At], transStorage->M[end - 1]);
      std::swap(transStorage->jcolM[col_At], transStorage->jcolM[end - 1]);
      transStorage->rowptr[row_At].end--;
   }
}


/** Should be called right after doSingletonRowsA() or another method that stores
 * information to update in the member variable colAdaptParent.
 * Updates the blocks A,C,F0,G0 using colAdaptParent.
 * Returns the number of newly found singleton rows (equality/inequality system)
 * during adaption of A,C,F0,G0.
 * Adapts the objective offset g only once for each column (variable).
 */
void StochPresolverBase::updateLinkingVarsBlocks(int& newSREq, int& newSRIneq)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   if( myRank == 0)
   {
      currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);
      for(int i=0; i<currgParent->n; i++)
            assert( std::isfinite(currgParent->elements()[i]) );
      for( int i=0; i<presData.getNumberColAdParent(); i++)
      {
         const int colIdx = presData.getColAdaptParent(i).colIdx;
         const double value = presData.getColAdaptParent(i).val;
         indivObjOffset += currgParent->elements()[colIdx] * value;
      }
   }

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

   updateNnzColParent(MPI_COMM_WORLD);

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

bool StochPresolverBase::newBoundsTightenOldBounds(double new_low, double new_upp, int index,
      double* ilow, double* iupp, double* low, double* upp) const
{
   if( ( ilow[index] != 0.0 && PIPSisLT(low[index], new_low) )
         || ( ilow[index] == 0.0 && new_low > -std::numeric_limits<double>::max() )
         || ( iupp[index] != 0.0 && PIPSisLT(new_upp, upp[index]) )
         || ( iupp[index] == 0.0 && new_upp < std::numeric_limits<double>::max() ) )
      return true;
   return false;
}

/**
 * Compares and sets the new lower and upper bounds, if they tighten the bounds.
 */
void StochPresolverBase::setNewBounds(int index, double new_low, double new_upp,
      double* ilow, double* low, double* iupp, double* upp) const
{
   if( (ilow[index] != 0.0 && PIPSisLT(low[index], new_low) )
      || (ilow[index] == 0.0 && new_low > -std::numeric_limits<double>::max()) )
   {
      ilow[index] = 1.0;
      low[index] = new_low;
   }
   if( (iupp[index] != 0.0 && PIPSisLT(new_upp, upp[index]) )
         || (iupp[index] == 0.0 && new_upp < std::numeric_limits<double>::max()))
   {
      iupp[index] = 1.0;
      upp[index] = new_upp;
   }
}

/**
 * Sets a given new bound at the specified index into the bound_vector
 * and sets 1.0 into the i_bound_vector.
 */
void StochPresolverBase::setNewBound(int index, double new_bound,
      SimpleVector* bound_vector, SimpleVector* i_bound_vector) const
{
   assert( bound_vector->n == i_bound_vector->n );
   assert( index >= 0 && index < bound_vector->n );

   bound_vector->elements()[index] = new_bound;
   if( i_bound_vector->elements()[index] == 0.0 )
      i_bound_vector->elements()[index] = 1.0;
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
   currIxlowParent = NULL;
   currxuppParent = NULL;
   currIxuppParent = NULL;
   currxlowChild = NULL;
   currIxlowChild = NULL;
   currxuppChild = NULL;
   currIxuppChild = NULL;
   currEqRhs = NULL;
   currIneqLhs = NULL;
   currIclow = NULL;
   currIneqRhs = NULL;
   currIcupp = NULL;
   currEqRhsLink = NULL;
   currIneqLhsLink = NULL;
   currIclowLink = NULL;
   currIneqRhsLink = NULL;
   currIcuppLink = NULL;

   currgParent = NULL;
   currgChild = NULL;

   currRedRow = NULL;
   currNnzRow = NULL;
   currRedRowLink = NULL;
   currRedColParent = NULL;
   currRedColChild = NULL;
   currNnzColParent = NULL;
   currNnzColChild = NULL;
}

/** Return false if it is a dummy child. Else, set the current pointers to Amat, AmatTrans, Bmat, BmatTrans,
 * currNnzRow, currRedRow, currNnzColChild, currRedColChild.
 * currxlowChild, currxuppChild, currIxlowChild, currIxuppChild.
 * And depending on the systemType: either currEqRhs or currIneqRhs, currIneqLhs, currIcupp and currIclow. */
bool StochPresolverBase::updateCurrentPointersForColAdapt(int it, SystemType system_type)
{
   assert( it >= -1 && it<nChildren );
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

/** Return false if it has no linking block. */
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


/** Set the current pointers for subsequent application of fixed (non-linking) variables.
 * Return false if it is a dummy child.
 */
bool StochPresolverBase::updateCPforAdaptFixationsBChild( int it, SystemType system_type)
{
   assert( it >= 0);
   setCurrentPointersToNull();

   // dummy child? set currBmat, currBmatTrans
   if( system_type == EQUALITY_SYSTEM )
   {
      if( !setCPBmatsChild(presProb->A, it, EQUALITY_SYSTEM))
         return false;
      setCPRowChildEquality(it);
   }
   else
   {
      if( !setCPBmatsChild(presProb->C, it, INEQUALITY_SYSTEM))
         return false;
      setCPRowChildInequality(it);
   }

   currRedColChild = dynamic_cast<SimpleVector*>(presData.redCol->children[it]->vec);
   currNnzColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);

   if( hasLinking(system_type) )
   {
      if( system_type == EQUALITY_SYSTEM)
      {
         setCPBlmatsChild(presProb->A, it);
         currRedRowLink = dynamic_cast<SimpleVector*>(presData.redRowA->vecl);
      }
      else
      {
         setCPBlmatsChild(presProb->C, it);
         currRedRowLink = dynamic_cast<SimpleVector*>(presData.redRowC->vecl);
      }
      currIcuppLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vecl);
      currIclowLink = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vecl);
   }
   return true;
}



/**
 * set all pointers to the currently necessary data
 * If node == -1 we are in the root node
 */
void StochPresolverBase::updatePointersForCurrentNode(int node, SystemType system_type)
{
   assert( !nodeIsDummy(node, system_type) );
   assert(-1 <= node && node <= nChildren );
   assert(system_type == EQUALITY_SYSTEM || system_type == INEQUALITY_SYSTEM);

   GenMatrixHandle matrix = (system_type == EQUALITY_SYSTEM) ? presProb->A : presProb->C;

   /* set matrix pointers for A B and Bl */
   setPointersMatrices(matrix, node);

   /* set lhs rhs for equations */
   setPointersMatrixBounds(system_type, node);

   /* set x lower upper bounds */
   setPointersVarBounds(node);

   /* set adaptions ? todo */

   /* set objective function pointers */
   setPointersObjective(node);

   /* set reduction pointers columns and rows */
   setReductionPointers(system_type, node);
}

// todo : set pointers NULL if no linking constraints?
void StochPresolverBase::setPointersMatrices(GenMatrixHandle mat, int node)
{
   assert(-1 <= node && node < nChildren);
   StochGenMatrix& smat = dynamic_cast<StochGenMatrix&>(*mat);

   /* in root node only B and Bl are present */
   if( node == -1 )
   {
      currAmat = NULL;
      currAmatTrans = NULL;

      currBmat = dynamic_cast<SparseGenMatrix*>(smat.Bmat)->getStorageDynamic();
      currBmatTrans =
            dynamic_cast<SparseGenMatrix*>(smat.Bmat)->getStorageDynamicTransposed();

      currBlmat =
            dynamic_cast<SparseGenMatrix*>(smat.Blmat)->getStorageDynamic();
      currBlmatTrans =
            dynamic_cast<SparseGenMatrix*>(smat.Blmat)->getStorageDynamicTransposed();
   }
   else
   {
      currAmat =
            dynamic_cast<SparseGenMatrix*>(smat.children[node]->Amat)->getStorageDynamic();
      currAmatTrans =
            dynamic_cast<SparseGenMatrix*>(smat.children[node]->Amat)->getStorageDynamicTransposed();

      currBmat =
            dynamic_cast<SparseGenMatrix*>(smat.children[node]->Bmat)->getStorageDynamic();
      currBmatTrans =
            dynamic_cast<SparseGenMatrix*>(smat.children[node]->Bmat)->getStorageDynamicTransposed();
      currBlmat =
            dynamic_cast<SparseGenMatrix*>(smat.children[node]->Blmat)->getStorageDynamic();
      currBlmatTrans =
            dynamic_cast<SparseGenMatrix*>(smat.children[node]->Blmat)->getStorageDynamicTransposed();
   }
}

void StochPresolverBase::setPointersMatrixBounds(SystemType system_type, int node)
{
   assert(-1 <= node && node < nChildren);

   /* non-linking constraints */
   StochVector& lhs = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochVector&>(*(presProb->bA))
         : dynamic_cast<StochVector&>(*(presProb->bl));
   StochVector& lhs_idx = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochVector&>(*(presProb->bA))
         : dynamic_cast<StochVector&>(*(presProb->iclow));
   StochVector& rhs = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochVector&>(*(presProb->bA))
         : dynamic_cast<StochVector&>(*(presProb->bu));
   StochVector& rhs_idx = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochVector&>(*(presProb->bA))
         : dynamic_cast<StochVector&>(*(presProb->icupp));

   if( system_type == EQUALITY_SYSTEM )
   {
      if( node == -1 )
      {
         currEqRhs = dynamic_cast<SimpleVector*>(rhs.vec);
      }
      else
      {
         currEqRhs = dynamic_cast<SimpleVector*>(rhs.children[node]->vec);
         assert(rhs.children[node]->vecl == NULL);
      }

      currIneqLhs = currIclow = currIneqRhs = currIcupp = currIneqLhsLink =
            currIclowLink = currIneqRhsLink = currIcuppLink = NULL;

      if( hasLinking(system_type) )
         currEqRhsLink = dynamic_cast<SimpleVector*>(rhs.vecl);
      else
         currEqRhsLink = NULL;
   }
   else
   {

      if( node == -1 )
      {
         currIneqLhs = dynamic_cast<SimpleVector*>(lhs.vec);
         currIclow = dynamic_cast<SimpleVector*>(lhs_idx.vec);
         currIneqRhs = dynamic_cast<SimpleVector*>(rhs.vec);
         currIcupp = dynamic_cast<SimpleVector*>(rhs_idx.vec);
      }
      else
      {
         currIneqLhs = dynamic_cast<SimpleVector*>(lhs.children[node]->vec);
         currIclow = dynamic_cast<SimpleVector*>(lhs_idx.children[node]->vec);
         currIneqRhs = dynamic_cast<SimpleVector*>(rhs.children[node]->vec);
         currIcupp = dynamic_cast<SimpleVector*>(rhs_idx.children[node]->vec);

         assert(lhs.children[node]->vecl == NULL);
         assert(lhs_idx.children[node]->vecl == NULL);
         assert(rhs.children[node]->vecl == NULL);
         assert(rhs_idx.children[node]->vecl == NULL);
      }

      currEqRhs = currEqRhsLink = NULL;

      if(hasLinking(system_type))
      {
         currIneqLhsLink = dynamic_cast<SimpleVector*>(lhs.vecl);
         currIclowLink = dynamic_cast<SimpleVector*>(lhs_idx.vecl);
         currIneqRhsLink = dynamic_cast<SimpleVector*>(rhs.vecl);
         currIcuppLink = dynamic_cast<SimpleVector*>(rhs_idx.vecl);
      }
      else
      {
         currIneqLhsLink = currIclowLink = currIneqRhsLink = currIcuppLink = NULL;
      }
   }
}

void StochPresolverBase::setPointersVarBounds(int node)
{
   assert(-1 <= node && node <= nChildren);

   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currIxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);
   currIxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).vec);

   assert(dynamic_cast<StochVector&>(*(presProb->blx)).vecl == NULL);
   assert(dynamic_cast<StochVector&>(*(presProb->ixlow)).vecl == NULL);
   assert(dynamic_cast<StochVector&>(*(presProb->bux)).vecl == NULL);
   assert(dynamic_cast<StochVector&>(*(presProb->ixupp)).vecl == NULL);

   if(node != -1)
   {
      currxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->vec);
      currxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->vec);
      currIxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[node]->vec);
      currIxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[node]->vec);

      assert(dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->vecl == NULL);
      assert(dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->vecl == NULL);
      assert(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[node]->vecl == NULL);
      assert(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[node]->vecl == NULL);
   }
   else
      currxlowChild = currxuppChild = currIxlowChild = currIxuppChild = NULL;
}

void StochPresolverBase::setPointersObjective(int node)
{
   currgParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).vec);
   if(node != -1)
   {
      currgChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[node]->vec);
      assert(dynamic_cast<StochVector&>(*(presProb->g)).children[node]->vecl == NULL);
   }
   else
      currgChild = NULL;

   assert(dynamic_cast<StochVector&>(*(presProb->g)).vecl == NULL);

}


void StochPresolverBase::setReductionPointers(SystemType system_type, int node){
   assert(-1 <= node && node <= nChildren);

   StochVector& row_red = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochVector&>(*(presData.redRowA))
            : dynamic_cast<StochVector&>(*(presData.redRowC));
   StochVector& row_nnz = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochVector&>(*(presData.nRowElemsA))
            : dynamic_cast<StochVector&>(*(presData.nRowElemsC));

   /* rows */
   if( node == -1)
   {
      currRedRow = dynamic_cast<SimpleVector*>(row_red.vec);
      currNnzRow = dynamic_cast<SimpleVector*>(row_nnz.vec);
   }
   else
   {
      assert(row_red.children[node]->vec != NULL);
      assert(row_nnz.children[node]->vec != NULL);

      currNnzRow = dynamic_cast<SimpleVector*>(row_nnz.children[node]->vec);
      currRedRow = dynamic_cast<SimpleVector*>(row_red.children[node]->vec);

      assert(row_red.children[node]->vecl == NULL);
      assert(row_nnz.children[node]->vecl == NULL);
   }

   if(hasLinking(system_type))
   {
      currRedRowLink = dynamic_cast<SimpleVector*>(row_red.vecl);
      currNnzRowLink = dynamic_cast<SimpleVector*>(row_nnz.vecl);;
   }
   else
      currRedRowLink = currNnzRowLink = NULL;

   /* colums */
   currRedColParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presData.redCol)).vec);
   currNnzColParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presData.nColElems)).vec);;

   assert(dynamic_cast<StochVector&>(*(presData.redCol)).vecl == NULL);
   assert(dynamic_cast<StochVector&>(*(presData.nColElems)).vecl == NULL);

   if(node != -1)
   {
      currRedColChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presData.redCol)).children[node]->vec);
      currNnzColChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presData.nColElems)).children[node]->vec);

      assert(dynamic_cast<StochVector&>(*(presData.redCol)).children[node]->vecl == NULL);
      assert(dynamic_cast<StochVector&>(*(presData.nColElems)).children[node]->vecl == NULL);
   }
   else
      currRedColChild = currNnzColChild = NULL;
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
   assert( it >= 0 && it<nChildren );
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   if( nodeIsDummy(it, system_type) )
      return false;
   currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamic();
   currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Amat)->getStorageDynamicTransposed();
   return true;
}

bool StochPresolverBase::setCPBmatsChild(GenMatrixHandle matrixHandle, int it, SystemType system_type)
{
   assert( it >= 0 && it<nChildren );
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   if( nodeIsDummy(it, system_type) )
      return false;
   currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamic();
   currBmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Bmat)->getStorageDynamicTransposed();
   return true;
}

/** Set the current pointers for currAmat and currBmat for child it>=0.
 * If it==-1, set only currAmat to the root Bmat.
 * Return false if it is a dummy child. */
bool StochPresolverBase::setPointersForAmatBmat(int node, SystemType system_type)
{
   assert( node >= -1 && node < nChildren );

   StochGenMatrix& matrix = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochGenMatrix&>(*presProb->A) : dynamic_cast<StochGenMatrix&>(*presProb->C);

   if( node == -1 ) // at root
   {
      currBmat = matrix.Bmat->getStorageDynamic();
      currBmatTrans = matrix.Bmat->getStorageDynamicTransposed();

   }
   else  // at child
   {
      if( nodeIsDummy(node, system_type) )
            return false;

      currAmat = dynamic_cast<SparseGenMatrix*>(matrix.children[node]->Amat)->getStorageDynamic();
      currAmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[node]->Amat)->getStorageDynamicTransposed();
      currBmat = dynamic_cast<SparseGenMatrix*>(matrix.children[node]->Bmat)->getStorageDynamic();
      currBmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[node]->Bmat)->getStorageDynamicTransposed();

   }
   return true;
}

bool StochPresolverBase::setCPLinkConstraint(GenMatrixHandle matrixHandle, int it, SystemType system_type)
{
   assert( it >= -1 && it<nChildren );
   assert( hasLinking(system_type) );

   if( it == -1 ) // F_0
      setCPBlmatsRoot(matrixHandle);
   else  // F_it
   {
      if( nodeIsDummy(it, system_type) )
         return false;
      setCPBlmatsChild( matrixHandle, it);
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
   assert( it >= 0 && it<nChildren );
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*matrixHandle);
   currBlmat = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamic();
   currBlmatTrans = dynamic_cast<SparseGenMatrix*>(matrix.children[it]->Blmat)->getStorageDynamicTransposed();
}

void StochPresolverBase::setCPColumnRoot()
{
   currRedColParent = dynamic_cast<SimpleVector*>(presData.redCol->vec);
   currxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).vec);
   currIxlowParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).vec);
   currxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).vec);
   currIxuppParent = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).vec);
}

void StochPresolverBase::setCPColumnChild(int it)
{
   assert( it >= 0 && it < nChildren );
   currRedColChild = dynamic_cast<SimpleVector*>(presData.redCol->children[it]->vec);
   currxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->blx)).children[it]->vec);
   currxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bux)).children[it]->vec);
   currIxlowChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixlow)).children[it]->vec);
   currIxuppChild = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->ixupp)).children[it]->vec);
}

void StochPresolverBase::setPointersRowRoot(SystemType system_type)
{
   if(system_type == EQUALITY_SYSTEM)
      setCPRowRootEquality();
   else if( system_type == INEQUALITY_SYSTEM)
      setCPRowRootInequality();
}

void StochPresolverBase::setCPRowRootEquality()
{
   currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).vec);
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowA->vec);
}

void StochPresolverBase::setCPRowRootInequality()
{
   setCPRowRootIneqOnlyLhsRhs();
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowC->vec);
}

void StochPresolverBase::setCPRowRootIneqOnlyLhsRhs()
{
   currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).vec);
   currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).vec);
   currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).vec);
   currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).vec);
}

void StochPresolverBase::setPointersRowChildNode(SystemType system_type, int child)
{
   if(system_type == EQUALITY_SYSTEM)
      setCPRowChildEquality(child);
   else if( system_type == INEQUALITY_SYSTEM)
      setCPRowChildInequality(child);
}

void StochPresolverBase::setCPRowChildEquality(int it)
{
   assert( it >= 0 && it < nChildren );

   currEqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bA)).children[it]->vec);
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowA->children[it]->vec);
}

void StochPresolverBase::setCPRowChildInequality(int it)
{
   assert( it >= 0 && it < nChildren );

   setCPRowChildIneqOnlyLhsRhs(it);
   currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
   currRedRow = dynamic_cast<SimpleVector*>(presData.redRowC->children[it]->vec);
}

void StochPresolverBase::setCPRowChildIneqOnlyLhsRhs(int it)
{
   assert( it >= 0 && it<nChildren );
   currIneqRhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bu)).children[it]->vec);
   currIneqLhs = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->bl)).children[it]->vec);
   currIcupp = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->icupp)).children[it]->vec);
   currIclow = dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->iclow)).children[it]->vec);
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
bool StochPresolverBase::removeEntryInDynamicStorage(SparseStorageDynamic& storage, const int rowIdx, const int colIdx, double& m) const
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

void StochPresolverBase::clearRow(SparseStorageDynamic& storage, const int rowIdx) const
{
   storage.rowptr[rowIdx].end = storage.rowptr[rowIdx].start;
}

/**
 * Remove row rowIdx in Ablock and Bblock. Removes the corresponding column in
 * AblockTrans and BblockTrans. Additionally, sets nnzRow[rowIdx] to 0.0.
 * Increments redColParent by one at each column index the row had an entry.
 * Decrements nnzColChild by one at each column index the row had an entry.
 */
// todo : should rhs and lhs be set to zero too?
void StochPresolverBase::removeRow(int rowIdx, SparseStorageDynamic& Ablock, SparseStorageDynamic& AblockTrans,
      SparseStorageDynamic* Bblock, SparseStorageDynamic* BblockTrans, SimpleVector& nnzRow,
      SimpleVector& redColParent, SimpleVector* nnzColChild)
{
   assert( rowIdx>=0 && rowIdx<Ablock.m );
   assert( Ablock.m == nnzRow.n );
   assert( Ablock.n == redColParent.n );

   const int rowStartA = Ablock.rowptr[rowIdx].start;
   const int rowEndA = Ablock.rowptr[rowIdx].end;
   // delete row in AblockTrans:
   for(int k=rowStartA; k<rowEndA; k++)
   {
      const int colIdx = Ablock.jcolM[k];
      double tmp = 0.0;
      removeEntryInDynamicStorage(AblockTrans, colIdx, rowIdx, tmp);
      // increment redColParent[colIdx]:
      redColParent.elements()[colIdx]++;
   }
   // delete row in Ablock:
   clearRow(Ablock, rowIdx);

   if(Bblock)
   {
      assert( Ablock.m == Bblock->m );
      removeRowInBblock( rowIdx, Bblock, BblockTrans, nnzColChild);
   }
   // set nnzRow[rowIdx] to 0.0:
   nnzRow.elements()[rowIdx] = 0.0;
}

/** Remove row rowIdx in Bblock which should not be a linking variable block.
 * Removes the corresponding column in BblockTrans.
 * Decrements nnzColChild by one at each column index the row had an entry.
 */
void StochPresolverBase::removeRowInBblock(int rowIdx, SparseStorageDynamic* Bblock,
      SparseStorageDynamic* BblockTrans, SimpleVector* nnzColChild)
{
   assert( Bblock && BblockTrans );
   assert( nnzColChild );
   assert( Bblock->n == nnzColChild->n );

   const int rowStartB = Bblock->rowptr[rowIdx].start;
   const int rowEndB = Bblock->rowptr[rowIdx].end;

   // delete row in BblockTrans:
   for(int k=rowStartB; k<rowEndB; k++)
   {
      const int colIdx = Bblock->jcolM[k];
      double tmp = 0.0;
      removeEntryInDynamicStorage(*BblockTrans, colIdx, rowIdx, tmp);
      // decrement nnzColChild[colIdx]:
      nnzColChild->elements()[colIdx]--;
   }
   // delete row in Bblock:
   clearRow(*Bblock, rowIdx);
}

bool StochPresolverBase::nodeIsDummy(int node, SystemType system_type) const
{
   assert( node >= -1 && node < nChildren );
   if( node == -1 )
      return false;
   StochGenMatrix& matrix = (system_type == EQUALITY_SYSTEM) ? dynamic_cast<StochGenMatrix&>(*presProb->A) : dynamic_cast<StochGenMatrix&>(*presProb->C);

   if( matrix.children[node]->isKindOf(kStochGenDummyMatrix))
   {
      assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->isKindOf(kStochDummy) );
      assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->isKindOf(kStochDummy) );
      assert( presData.redCol->children[node]->isKindOf(kStochDummy) );

      if( system_type == EQUALITY_SYSTEM)
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bA)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bux)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->blx)).children[node]->isKindOf(kStochDummy) );
         assert( presData.nRowElemsA->children[node]->isKindOf(kStochDummy) );
         assert( presData.redRowA->children[node]->isKindOf(kStochDummy) );
      }
      else
      {
         assert( dynamic_cast<StochVector&>(*(presProb->bu)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->bl)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->icupp)).children[node]->isKindOf(kStochDummy) );
         assert( dynamic_cast<StochVector&>(*(presProb->iclow)).children[node]->isKindOf(kStochDummy) );
         assert( presData.nRowElemsC->children[node]->isKindOf(kStochDummy) );
         assert( presData.redRowC->children[node]->isKindOf(kStochDummy) );
      }
      return true;
   }
   return false;
}

bool StochPresolverBase::hasLinking(SystemType system_type) const
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

void StochPresolverBase::getRankDistributed( MPI_Comm comm, int& myRank, bool& iAmDistrib ) const
{
   MPI_Comm_rank(comm, &myRank);
   int world_size;
   MPI_Comm_size(comm, &world_size);
   if( world_size > 1) iAmDistrib = true;
   else iAmDistrib = false;
}

/** Call MPI_Abort() and print an error message */
void StochPresolverBase::abortInfeasible(MPI_Comm comm) const
{
   cout<<"Infesibility detected in presolving. Aborting now."<<endl;
   MPI_Abort(comm, 1);
}

void StochPresolverBase::synchronize(int& value) const
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );
   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void StochPresolverBase::synchronizeSum(int& first, int& second) const
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );
   if( iAmDistrib )
   {
      int* newSR = new int[2];
      newSR[0] = first;
      newSR[1] = second;
      MPI_Allreduce(MPI_IN_PLACE, newSR, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      first = newSR[0];
      second = newSR[1];
      delete[] newSR;
   }
}

/*
 * Given a vector<COLUMNTOADAPT>, this routine goes through all columns inside and removes them
 * from the current Bmat block. Depending on the system_type, rhs (and lhs) are updated.
 */
void StochPresolverBase::adaptChildBmat( std::vector<COLUMNTOADAPT> const & colAdaptBlock, SystemType system_type, int& newSR )
{
   for(int i=0; i<(int)colAdaptBlock.size(); i++)
   {
      int colIdx = colAdaptBlock[i].colIdx;
      double val = colAdaptBlock[i].val;

      for( int j = currBmatTrans->rowptr[colIdx].start; j<currBmatTrans->rowptr[colIdx].end; j++ )
      {
         int rowIdx = currBmatTrans->jcolM[j];
         double m = 0.0;
         if( !removeEntryInDynamicStorage(*currBmat, rowIdx, colIdx, m) )
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
}

void StochPresolverBase::adaptChildBlmat( std::vector<COLUMNTOADAPT> const & colAdaptBlock, SystemType system_type)
{
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

   for(size_t i = 0; i < colAdaptBlock.size(); i++)
   {
      const int colIdx = colAdaptBlock[i].colIdx;
      const double val = colAdaptBlock[i].val;

      for( int j = currBlmatTrans->rowptr[colIdx].start; j < currBlmatTrans->rowptr[colIdx].end; j++ )
      {
         const int rowIdx = currBlmatTrans->jcolM[j];
         double m = 0.0;

         if( !removeEntryInDynamicStorage(*currBlmat, rowIdx, colIdx, m) )
            continue;

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
}

/** For the given column index and the value to which this variable is to be fixed,
 * go through the current Bmat by going through the corresponding row of the transposed Bmat.
 * Adapt the rhs, the nnzRow, the nnzCol and all concerned entries in Bmat and BmatTransposed.
 * Return the number of newly found singleton rows.
 */
int StochPresolverBase::adaptChildBmatCol(int colIdx, double val, SystemType system_type)
{
   assert(currBmat && currBmatTrans && currNnzColChild && currRedColChild && currNnzRow && currRedRow);
   assert( colIdx >= 0 && colIdx < currNnzColChild->n );
   int newSingletonRows = 0;

   for( int j = currBmatTrans->rowptr[colIdx].start; j<currBmatTrans->rowptr[colIdx].end; j++ )
   {
      const int rowIdxB = currBmatTrans->jcolM[j];
      double m = 0.0;
      if( !removeEntryInDynamicStorage(*currBmat, rowIdxB, colIdx, m) )
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
void StochPresolverBase::adaptOtherSystemChildB(SystemType system_type, std::vector<COLUMNTOADAPT> const & colAdaptBblock, int& newSR )
{
   // Bmat blocks
   adaptChildBmat(colAdaptBblock, system_type, newSR);

   // Blmat blocks
   if( hasLinking(system_type) )
      adaptChildBlmat(colAdaptBblock, system_type);
}

/** Adapt the columns for the linking-variable-blocks (the A_i) blocks */
int StochPresolverBase::colAdaptLinkVars(int it, SystemType system_type)
{
   assert( it >= -1 && it<nChildren );
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int newSingletonRows = 0;

   for( int i=0; i < presData.getNumberColAdParent(); i++)
   {
      const int colIdxA = presData.getColAdaptParent(i).colIdx;
      const double val = presData.getColAdaptParent(i).val;

      for( int j = currAmatTrans->rowptr[colIdxA].start; j < currAmatTrans->rowptr[colIdxA].end; j++ )
      {
         const int rowIdxA = currAmatTrans->jcolM[j];
         double m = 0.0;
         if( !removeEntryInDynamicStorage(*currAmat, rowIdxA, colIdxA, m) )
            continue;
         //PIPSdebugMessage("Removed entry (%d, %d) with value %f in Amat of child %d of system_type %d \n", rowIdxA, colIdxA, m, it, system_type);

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
      const int colIdx = presData.getColAdaptParent(i).colIdx;
      const double val = presData.getColAdaptParent(i).val;

      for( int j = currBlmatTrans->rowptr[colIdx].start; j<currBlmatTrans->rowptr[colIdx].end; j++ )
      {
         const int rowIdx = currBlmatTrans->jcolM[j];
         double m = 0.0;
         if( !removeEntryInDynamicStorage(*currBlmat, rowIdx, colIdx, m) )
            continue;
         //PIPSdebugMessage("Removed entry (%d, %d) with value %f in F_0 of system_type %d \n", rowIdx, colIdx, m, system_type);

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

bool StochPresolverBase::newBoundsImplyInfeasible(double newxlow, double newxupp, int colIdx,
      double* ixlow, double* ixupp, double* xlow, double* xupp) const
{
   assert( colIdx >= 0 );
   if( ( ixlow[colIdx] != 0.0 && PIPSisRelLT(newxupp, xlow[colIdx]) )
         || (ixupp[colIdx] != 0.0 && PIPSisRelLT(xupp[colIdx], newxlow) )
         || (newxlow > newxupp))
   {
      std::cout << "Presolving detected infeasibility: variable: " << colIdx << "\tnew bounds = [" << newxlow << ", " << newxupp << "]" << "\told bounds: [" << xlow <<
    		  ", " << xupp << "]" << endl;
      return true;
   }
   return false;
}

bool StochPresolverBase::newBoundsFixVariable(double& value, double newxlow, double newxupp, int colIdx,
      double* ixlow, double* ixupp, double* xlow, double* xupp) const
{
   assert( colIdx >= 0 );
   if( PIPSisEQ(newxlow, newxupp)
         || ( ixlow[colIdx] != 0.0 && PIPSisEQ( xlow[colIdx], newxupp) ))
   {
      value = newxupp;
      return true;
   }
   else if( ixupp[colIdx] != 0.0 && PIPSisEQ(xupp[colIdx], newxlow) )
   {
      value = newxlow;
      return true;
   }
   // if relative difference between newxlow and newxupp is below a threshold, fix the variable:
   double upperbound = newxupp;
   if( ixupp[colIdx] != 0.0 && xupp[colIdx]<newxupp )
      upperbound = xupp[colIdx];
   double lowerbound = newxlow;
   if( ixlow[colIdx] != 0.0 && xlow[colIdx]>newxlow )
      lowerbound = xlow[colIdx];
   double absmax = std::max(fabs(upperbound), fabs(lowerbound) );
   double absdiff = fabs( upperbound - lowerbound );
   if( absdiff / absmax < tolerance4 )
   {
      // verify if one of the bounds is integer:
      double intpart;
      if( modf(lowerbound, &intpart) == 0.0 )
      {
         value = lowerbound;
         return true;
      }
      else if( modf(upperbound, &intpart) == 0.0 )
      {
         value = upperbound;
         return true;
      }
      else  // set the variable to the arithmetic mean:
      {
         value = (lowerbound + upperbound ) / 2.0;
         return true;
      }
   }
   return false;
}

/** Fix the (non-linking) variable with index colIdx to the value val, by adapting the
 * objective offset accordingly and by removing all variable entries in the currBmat / currBmatTrans
 * (by calling adaptChildBmatCol).
 * Also, adds the pair <colIdx, val> to the vector colAdaptLinkBlock.
 * Return the number of newly found singleton Rows in this matrix block.
 * The following pointers have to be set correctly: currgChild, currBmat, currBmatTrans,
 * currNnzColChild, currRedColChild, currNnzRow, currRedRow and currEqRhs (if EQUALITY_SYSTEM)
 * or the corresponding equivalent (if INEQUALITY_SYSTEM).
 */
int StochPresolverBase::fixVarInChildBlockAndStore( int colIdx, double val, SystemType system_type,
      std::vector<COLUMNTOADAPT> & colAdaptLinkBlock )
{
   assert( currgChild );
   assert( colIdx >= 0 && colIdx < currgChild->n );
   //cout<<"New bounds imply fixation of variable "<<colIdx<<" of a child to value: "<<val<<endl;

   // adapting objOffset:
   indivObjOffset += currgChild->elements()[colIdx] * val;

   // adapt in the currBmat/currBmatTrans by removing column colIdx and store in colADaptLinkBlock for G, F, B:
   int newSR = adaptChildBmatCol(colIdx, val, system_type);

   COLUMNTOADAPT colWithVal = {colIdx, val};
   colAdaptLinkBlock.push_back(colWithVal);

   return newSR;
}

/** Stores the column index colIdx together with the value as a COLUMNTOADAPT in colAdaptParent.
 * Aborts if infeasibility is detected.
 *
 * todo : use std::find and stuff
 */
void StochPresolverBase::storeColValInColAdaptParent(int colIdx, double value)
{
   const COLUMNTOADAPT colWithVal = {colIdx, value};

   bool uniqueAdditionToOffset = true;
   for(int i = 0; i < presData.getNumberColAdParent(); i++)
   {
      if( presData.getColAdaptParent(i).colIdx == colIdx )
      {

         if( !PIPSisEQ(presData.getColAdaptParent(i).val, value) )
         {
            std::cout << "Presolving detected infeasibility : fixation of variable that has previously been fixed to a different value" << std::endl;
            abortInfeasible(MPI_COMM_WORLD);
         }

         uniqueAdditionToOffset = false;
      }
   }
   if( uniqueAdditionToOffset )
      presData.addColToAdaptParent(colWithVal);
}

/** Stores the column index colIdx together with the new bounds as a XBOUNDS in newBoundsParent.
 * Should be called only from Process Zero.
 * Returns false if infeasibility is detected (contradictory bounds).
 */
void StochPresolverBase::storeNewBoundsParent(int colIdx, double newxlow, double newxupp)
{
   assert( colIdx >= 0 );
   XBOUNDS newXbounds = {colIdx, newxlow, newxupp};
   for(size_t i = 0; i < newBoundsParent.size(); i++)
   {
      if( newBoundsParent[i].colIdx == colIdx )
      {
         if( PIPSisLT(newxupp, newBoundsParent[i].newxlow) || PIPSisLT(newBoundsParent[i].newxupp, newxlow) )
         {
        	 std::cout << "Presolving detected infeasibility. Two change of bounds requested to invalid values: bounds_a = [" << newxlow << ", " << newxupp << "]\tbounds_b = ["
        			 << newBoundsParent[i].newxlow << ", " << newBoundsParent[i].newxupp << "]" << std::endl;
            abortInfeasible(MPI_COMM_WORLD);
         }
      }
   }
   newBoundsParent.push_back(newXbounds);
}

/** Method similar to combineColAdaptParent(), that is a method going through newBoundsParent
 * and cleaning it up, removing redundant bounds, checking for infeasibility or more tightening.
 */
// todo why not just MPI_ALLREDUCE with min and max????
void StochPresolverBase::combineNewBoundsParent()
{
   int myRank, world_size;
   bool iAmDistrib = false;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   if( world_size > 1) iAmDistrib = true;

   if( iAmDistrib )
   {
      // allgather the length of each newBoundsParent
      int mylen = getNumberNewBoundsParent();
      int* recvcounts = new int[world_size];

      MPI_Allgather(&mylen, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

      // allgatherv the actual newBoundsParent
      // First, extract the colIdx and val into int* and double* arrays:
      int* colIndicesLocal = new int[mylen];
      double* xlowLocal = new double[mylen];
      double* xuppLocal = new double[mylen];
      for(int i=0; i<mylen; i++)
      {
         colIndicesLocal[i] = getNewBoundsParent(i).colIdx;
         xlowLocal[i] = getNewBoundsParent(i).newxlow;
         xuppLocal[i] = getNewBoundsParent(i).newxupp;
      }
      // Second, prepare the receive buffers:
      int lenghtGlobal = recvcounts[0];
      int* displs = new int[world_size];
      displs[0] = 0;
      for(int i=1; i<world_size; i++)
      {
         lenghtGlobal += recvcounts[i];
         displs[i] = displs[i-1] + recvcounts[i-1];
      }
      int* colIndicesGlobal = new int[lenghtGlobal];
      double* xlowGlobal = new double[lenghtGlobal];
      double* xuppGlobal = new double[lenghtGlobal];
      // Then, do the actual MPI communication:
      MPI_Allgatherv(colIndicesLocal, mylen, MPI_INT, colIndicesGlobal, recvcounts, displs , MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(xlowLocal, mylen, MPI_DOUBLE, xlowGlobal, recvcounts, displs , MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgatherv(xuppLocal, mylen, MPI_DOUBLE, xuppGlobal, recvcounts, displs , MPI_DOUBLE, MPI_COMM_WORLD);

      // Reconstruct a newBoundsParent containing all entries:
      clearNewBoundsParent();
      for(int i=0; i<lenghtGlobal; i++)
      {
         XBOUNDS newXBound = {colIndicesGlobal[i], xlowGlobal[i], xuppGlobal[i]};
         addNewBoundsParent(newXBound);
      }

      delete[] recvcounts;
      delete[] colIndicesLocal;
      delete[] xlowLocal;
      delete[] xuppLocal;
      delete[] displs;
      delete[] colIndicesGlobal;
      delete[] xlowGlobal;
      delete[] xuppGlobal;
   }

   // Sort colIndicesGlobal (and xlowGlobal, xuppGlobal accordingly), remove duplicates,
   // tighten bounds and find infeasibilities
   std::sort(newBoundsParent.begin(), newBoundsParent.end(), xbounds_col_is_smaller());

   if(getNumberNewBoundsParent() > 0)
   {
      int colIdxCurrent = getNewBoundsParent(0).colIdx;
      double xlowCurrent = getNewBoundsParent(0).newxlow;
      double xuppCurrent = getNewBoundsParent(0).newxupp;
      for(int i=1; i<getNumberNewBoundsParent(); i++)
      {
         if( getNewBoundsParent(i).colIdx == colIdxCurrent )
         {
            const double bestLow = max(xlowCurrent, getNewBoundsParent(i).newxlow);
            const double bestUpp = min(xuppCurrent, getNewBoundsParent(i).newxupp);
            if( bestLow > bestUpp )
            {
               cout<<"Detected infeasibility in variable "<<colIdxCurrent<<" of parent. bestLow="<<bestLow<<", bestUpp="<<bestUpp<<endl;
               abortInfeasible(MPI_COMM_WORLD);
            }
            else
            {
               // change the vector element newBoundsParent.begin()+(i-1), also das,
               // welches colIdxCurrent definiert hat:
               setNewBoundsParent(i-1, colIdxCurrent, bestLow, bestUpp);
               newBoundsParent.erase(newBoundsParent.begin()+i);   //todo: implement more efficiently
               i--;  // if i is not decremented, then the next entry in newBoundsParent would be omitted
            }
         }
         else
         {
            colIdxCurrent = getNewBoundsParent(i).colIdx;
            xlowCurrent = getNewBoundsParent(i).newxlow;
            xuppCurrent = getNewBoundsParent(i).newxupp;
         }
      }
   }
   assert( getNumberNewBoundsParent() <= presData.nColElems->vec->n );
}

XBOUNDS StochPresolverBase::getNewBoundsParent(int i) const
{
   assert( i<getNumberNewBoundsParent() );
   return newBoundsParent[i];
}
void StochPresolverBase::setNewBoundsParent(int i, int colIdx, double newxlow, double newxupp)
{
   assert( i<getNumberNewBoundsParent() );
   newBoundsParent[i].colIdx = colIdx;
   newBoundsParent[i].newxlow = newxlow;
   newBoundsParent[i].newxupp = newxupp;
}
int StochPresolverBase::getNumberNewBoundsParent() const
{
   return (int)newBoundsParent.size();
}
void StochPresolverBase::addNewBoundsParent(XBOUNDS newXBounds)
{
   newBoundsParent.push_back(newXBounds);
}
void StochPresolverBase::clearNewBoundsParent()
{
   newBoundsParent.clear();
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

/**
 * Compute the minimum and maximum activity of the row rowIdx in matrix. If colIdx!=-1, then this entry
 * is excluded in the computation of the activities.
 */
void StochPresolverBase::computeActivityBlockwise( SparseStorageDynamic& matrix, int rowIdx, int colIdx,
      double& infRow, double& supRow,
      SimpleVector& xlow, SimpleVector& ixlow, SimpleVector& xupp, SimpleVector& ixupp)
{
   // todo: possibly add two bools: if infty, indicate if at least two bounds were infty
   assert( rowIdx >= 0 && rowIdx < matrix.m );
   assert( colIdx >= -1 && colIdx < matrix.n );
   assert( xlow.n == matrix.n && ixlow.n == matrix.n && xupp.n == matrix.n && ixupp.n == matrix.n );

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

/** Verifies if the nnzCounters are still correct. */
// TODO BUG?
bool StochPresolverBase::verifyNnzcounters()
{
   // todo: distributed case ok?
   bool nnzCorrect = true;
   StochVectorHandle nnzColOrig(dynamic_cast<StochVector*>(presData.nColElems->cloneFull()));
   StochVectorHandle nnzRowAOrig(dynamic_cast<StochVector*>(presData.nRowElemsA->cloneFull()));
   StochVectorHandle nnzRowCOrig(dynamic_cast<StochVector*>(presData.nRowElemsC->cloneFull()));

   presData.nColElems->setToZero();
   presData.nRowElemsA->setToZero();
   presData.nRowElemsC->setToZero();

   // similar to presData.initNnzCounter():
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*(presProb->C));
   StochVectorHandle colClone(dynamic_cast<StochVector*>(presData.nColElems->clone()));
   A.getNnzPerRow(*presData.nRowElemsA);
   C.getNnzPerRow(*presData.nRowElemsC);
   A.getNnzPerCol(*presData.nColElems);
   C.getNnzPerCol(*colClone);
   presData.nColElems->axpy(1.0, *colClone);

   // linking variables:
   SimpleVector* nColOrigSimple = dynamic_cast<SimpleVector*>(nnzColOrig->vec);
   SimpleVector* nColUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nColElems->vec);
   assert( nColUpdatedSimple->n == nColOrigSimple->n );
   for( int i=0; i<nColUpdatedSimple->n; i++)
   {
      if( nColUpdatedSimple->elements()[i] != nColOrigSimple->elements()[i])
      {
         cout<<"Nnz Counter linking column "<<i<<" not correct: "<<nColUpdatedSimple->elements()[i]<<" vs. "<<nColOrigSimple->elements()[i]<<endl;
         nnzCorrect = false;
//         break;
      }
   }
   // non linking variables:
   for( size_t it = 0; it < presData.nColElems->children.size(); it++)
   {
      nColOrigSimple = dynamic_cast<SimpleVector*>(nnzColOrig->children[it]->vec);
      nColUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
      assert( nColUpdatedSimple->n == nColOrigSimple->n );
      for( int i=0; i<nColUpdatedSimple->n; i++)
      {
         if( nColUpdatedSimple->elements()[i] != nColOrigSimple->elements()[i])
         {
            cout<<"Nnz Counter non-linking column "<<i<<" of child "<<(int)it<<" not correct: "<<nColUpdatedSimple->elements()[i]<<" vs. "<<nColOrigSimple->elements()[i]<<endl;
            nnzCorrect = false;
//            break;
         }
      }
   }
   // rows A:
   SimpleVector* nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzRowAOrig->vec);
   SimpleVector* nRowAUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
   assert( nRowAUpdatedSimple->n == nRowAOrigSimple->n );
   for( int i=0; i < nRowAUpdatedSimple->n; i++)
   {
      if( nRowAUpdatedSimple->elements()[i] != nRowAOrigSimple->elements()[i])
      {
         cout<<"Nnz Counter root A row "<<i<<" not correct: "<<nRowAUpdatedSimple->elements()[i]<<" vs. "<<nRowAOrigSimple->elements()[i]<<endl;
         nnzCorrect = false;
//         break;
      }
   }
   // child rows:
   for( size_t it = 0; it < presData.nRowElemsA->children.size(); it++)
   {
      nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzRowAOrig->children[it]->vec);
      nRowAUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
      assert( nRowAUpdatedSimple->n == nRowAOrigSimple->n );
      for( int i = 0; i < nRowAUpdatedSimple->n; i++)
      {
         if( nRowAUpdatedSimple->elements()[i] != nRowAOrigSimple->elements()[i])
         {
            cout<<"Nnz Counter non-linking A row "<<i<<" of child "<<(int)it<<" not correct: "<<nRowAUpdatedSimple->elements()[i]<<" vs. "<<nRowAOrigSimple->elements()[i]<<endl;
            nnzCorrect = false;
//            break;
         }
      }
   }
   if(nnzRowAOrig->vecl) // linking rows:
   {
      nRowAOrigSimple = dynamic_cast<SimpleVector*>(nnzRowAOrig->vecl);
      nRowAUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      assert( nRowAUpdatedSimple->n == nRowAOrigSimple->n );
      for( int i=0; i<nRowAUpdatedSimple->n; i++)
      {
         if( nRowAUpdatedSimple->elements()[i] != nRowAOrigSimple->elements()[i])
         {
            cout<<"Nnz Counter linking row of A "<<i<<" not correct: "<<nRowAUpdatedSimple->elements()[i]<<" vs. "<<nRowAOrigSimple->elements()[i]<<endl;
            nnzCorrect = false;
//            break;
         }
      }
   }
   // rows C:
   SimpleVector* nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzRowCOrig->vec);
   SimpleVector* nRowCUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec);
   assert( nRowCUpdatedSimple->n == nRowCOrigSimple->n );
   for( int i=0; i<nRowCUpdatedSimple->n; i++)
   {
      if( nRowCUpdatedSimple->elements()[i] != nRowCOrigSimple->elements()[i])
      {
         cout<<"Nnz Counter root C row "<<i<<" not correct: "<<nRowCUpdatedSimple->elements()[i]<<" vs. "<<nRowCOrigSimple->elements()[i]<<endl;
         nnzCorrect = false;
//         break;
      }
   }
   // child rows:
   for( size_t it = 0; it < presData.nRowElemsC->children.size(); it++)
   {
      nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzRowCOrig->children[it]->vec);
      nRowCUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
      assert( nRowCUpdatedSimple->n == nRowCOrigSimple->n );
      for( int i=0; i<nRowCUpdatedSimple->n; i++)
      {
         if( nRowCUpdatedSimple->elements()[i] != nRowCOrigSimple->elements()[i])
         {
            cout<<"Nnz Counter non-linking C row "<<i<<" of child "<<(int)it<<" not correct: "<<nRowCUpdatedSimple->elements()[i]<<" vs. "<<nRowCOrigSimple->elements()[i]<<endl;
            nnzCorrect = false;
//            break;
         }
      }
   }
   if(nnzRowCOrig->vecl) // linking rows:
   {
      nRowCOrigSimple = dynamic_cast<SimpleVector*>(nnzRowCOrig->vecl);
      nRowCUpdatedSimple = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      assert( nRowCUpdatedSimple->n == nRowCOrigSimple->n );
      for( int i=0; i<nRowCUpdatedSimple->n; i++)
      {
         if( nRowCUpdatedSimple->elements()[i] != nRowCOrigSimple->elements()[i])
         {
            cout<<"Nnz Counter linking row of C "<<i<<" not correct: "<<nRowCUpdatedSimple->elements()[i]<<" vs. "<<nRowCOrigSimple->elements()[i]<<endl;
            nnzCorrect = false;
//            break;
         }
      }
   }
   return nnzCorrect;
}

void StochPresolverBase::countRowsCols()// method is const but changes pointers
{
   assert( static_cast<int>(presData.nRowElemsC->children.size()) == nChildren);
   assert( presData.nRowElemsA->children.size() == presData.nRowElemsC->children.size());
   assert( static_cast<int>(dynamic_cast<StochVector&>(*(presProb->icupp)).children.size()) == nChildren);

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   int n_rows_eq = 0;
   int n_rows_eq_linking = 0;
   int n_rows_ineq = 0;
   int n_rows_ineq_linking = 0;
   int n_fixed_rows = 0;
   int n_ranged_rows = 0;
   int n_singleton_rows_eq = 0;
   int n_singleton_rows_ineq = 0;
   int n_cols = 0;
   int n_boxed_cols = 0;
   int n_free_cols = 0;

   /* root nodes of equality and inequality system - linking and non linking */
   if( myRank == 0 )
   {
      updatePointersForCurrentNode(-1, EQUALITY_SYSTEM);

      int n_rows_linking_ranged = 0;
      int n_rows_eq_linking_singleton = 0;
      int n_rows_ineq_linking_singleton = 0;

      countRowsBlock(n_rows_eq, n_ranged_rows, n_fixed_rows, n_singleton_rows_eq, EQUALITY_SYSTEM, CHILD_BLOCK);
      countRowsBlock(n_rows_eq_linking, n_rows_linking_ranged, n_fixed_rows, n_rows_eq_linking_singleton, EQUALITY_SYSTEM, LINKING_CONS_BLOCK);
      n_singleton_rows_eq += n_rows_eq_linking_singleton;
      assert(n_rows_linking_ranged == 0);

      updatePointersForCurrentNode(-1, INEQUALITY_SYSTEM);

      countRowsBlock(n_rows_ineq, n_ranged_rows, n_fixed_rows, n_singleton_rows_ineq, INEQUALITY_SYSTEM, CHILD_BLOCK);
      countRowsBlock(n_rows_ineq_linking, n_rows_linking_ranged, n_fixed_rows, n_rows_ineq_linking_singleton, INEQUALITY_SYSTEM, LINKING_CONS_BLOCK);
      n_singleton_rows_ineq += n_rows_ineq_linking_singleton;
      n_ranged_rows += n_rows_linking_ranged;

      countBoxedColumns( n_boxed_cols, n_cols, n_free_cols, LINKING_VARS_BLOCK);

      std::cout << "#Linking_vars:\t" << n_cols << " (#free: " << n_free_cols << ", #boxed: " << n_boxed_cols << ")" << std::endl;

      std::cout << "#rows B0:\t" << n_rows_eq << std::endl;
      std::cout << "#rows Bl_0:\t" << n_rows_eq_linking << " (#singleton: " << n_rows_eq_linking_singleton << ")" << std::endl;
      std::cout << "#rows D0:\t" << n_rows_ineq << std::endl;
      std::cout << "#rows Dl_0:\t" << n_rows_ineq_linking << " (#singleton: " << n_rows_ineq_linking_singleton << ", #ranged: " << n_rows_linking_ranged << ")" << std::endl;
   }

   /* child nodes in both systems */
   for( int node = 0; node < nChildren; node++)
   {
      assert( (nodeIsDummy( node, EQUALITY_SYSTEM) && nodeIsDummy( node, INEQUALITY_SYSTEM)) ||
            (!nodeIsDummy( node, EQUALITY_SYSTEM) && !nodeIsDummy( node, INEQUALITY_SYSTEM) ));

      /* equality system */
      if(!nodeIsDummy( node, EQUALITY_SYSTEM))
      {
         updatePointersForCurrentNode(node, EQUALITY_SYSTEM);
         countRowsBlock(n_rows_eq, n_ranged_rows, n_fixed_rows, n_singleton_rows_eq, EQUALITY_SYSTEM, CHILD_BLOCK);
      }

      /* inequality system */
      if( !nodeIsDummy( node, INEQUALITY_SYSTEM) )
      {
         updatePointersForCurrentNode(node, INEQUALITY_SYSTEM);
         countRowsBlock(n_rows_ineq, n_ranged_rows, n_fixed_rows, n_singleton_rows_ineq, INEQUALITY_SYSTEM, CHILD_BLOCK);

         countBoxedColumns( n_boxed_cols, n_cols, n_free_cols, CHILD_BLOCK);
      }
   }

#if 0//TIMING // TODO
   // count how many linking rows do not really link two blocks:
   if( hasLinking(EQUALITY_SYSTEM) )
   {
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);
      int* rowHasEntryInBlocks = new int[currNnzRow->n];
      for( int i = 0; i < currNnzRow->n; i++ )
         rowHasEntryInBlocks[i] = 0;
      for( size_t it = 0; it < presData.nRowElemsA->children.size(); it++)
      {
         if( !nodeIsDummy( it, EQUALITY_SYSTEM))
         {
            setCPBlmatsChild(presProb->A, (int)it);
            for( int i = 0; i < currNnzRow->n; i++ )
            {
               if( currNnzRow->elements()[i] != 0.0 )
                  if( currBlmat->rowptr[i].start != currBlmat->rowptr[i].end)
                     rowHasEntryInBlocks[i]++;
            }
         }
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vecl);

      }
      MPI_Allreduce(MPI_IN_PLACE, rowHasEntryInBlocks, currNnzRow->n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      int linkRows1Blocks = 0;
      int linkRows2Blocks = 0;
      for( int i = 0; i < currNnzRow->n; i++ )
      {
         if(rowHasEntryInBlocks[i] == 1)
            linkRows1Blocks++;
         if(rowHasEntryInBlocks[i] == 2)
            linkRows2Blocks++;
      }
      if( myRank == 0 )
      {
         cout<<"1-link rows in A: "<<linkRows1Blocks<<endl;
         cout<<"2-link rows in A: "<<linkRows2Blocks<<endl;
      }
   }
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      int* rowHasEntryInBlocks = new int[currNnzRow->n];
      for( int i = 0; i < currNnzRow->n; i++ )
         rowHasEntryInBlocks[i] = 0;
      for( size_t it = 0; it < presData.nRowElemsC->children.size(); it++)
      {
         if( !nodeIsDummy( it, INEQUALITY_SYSTEM))
         {
            setCPBlmatsChild(presProb->C, (int)it);
            for( int i = 0; i < currNnzRow->n; i++ )
            {
               if( currNnzRow->elements()[i] != 0.0 )
                  if( currBlmat->rowptr[i].start != currBlmat->rowptr[i].end)
                     rowHasEntryInBlocks[i]++;
            }
         }
         currNnzRow = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vecl);
      }

      MPI_Allreduce(MPI_IN_PLACE, rowHasEntryInBlocks, currNnzRow->n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      int linkRows1Blocks = 0;
      int linkRows2Blocks = 0;

      for( int i = 0; i < currNnzRow->n; i++ )
      {
         if(rowHasEntryInBlocks[i] == 1)
            linkRows1Blocks++;
         if(rowHasEntryInBlocks[i] == 2)
            linkRows2Blocks++;
      }
      if( myRank == 0 )
      {
         cout<<"1-link rows in C: "<<linkRows1Blocks<<endl;
         cout<<"2-link rows in C: "<<linkRows2Blocks<<endl;
      }
   }
#endif

   /* sync data */
   if( iAmDistrib )
   {
      int* count = new int[9];
      count[0] = n_rows_eq;
      count[1] = n_rows_ineq;
      count[2] = n_fixed_rows;
      count[3] = n_ranged_rows;
      count[4] = n_singleton_rows_eq;
      count[5] = n_singleton_rows_ineq;
      count[6] = n_cols;
      count[7] = n_boxed_cols;
      count[8] = n_free_cols;

      MPI_Allreduce(MPI_IN_PLACE, count, 9, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      n_rows_eq = count[0];
      n_rows_ineq = count[1];
      n_fixed_rows = count[2];
      n_ranged_rows = count[3];
      n_singleton_rows_eq = count[4];
      n_singleton_rows_ineq = count[5];
      n_cols = count[6];
      n_boxed_cols = count[7];
      n_free_cols = count[8];

      delete[] count;
   }

   if( myRank == 0 )
   {
      std::cout << "#rows_total:\t" << n_rows_eq + n_rows_ineq << " (#fixed: " << n_fixed_rows << ", #ranged: " << n_ranged_rows << ", #singleton: " << n_singleton_rows_eq + n_singleton_rows_ineq<< ")" << std::endl;
      std::cout << "#rows A:\t" << n_rows_eq << " (#singleton: " << n_singleton_rows_eq << ")" << std::endl;
      std::cout << "#rows C:\t" << n_rows_ineq << " (#singleton: " << n_singleton_rows_ineq << ")" << std::endl;

      std::cout <<"#vars_total:\t" << n_cols << " (#bounded: " << n_boxed_cols << ", #free: " << n_free_cols << ")" << std::endl;
   }
}

void StochPresolverBase::countRowsBlock(int& n_rows, int& n_ranged_rows, int& n_fixed_rows, int& n_singleton_rows, SystemType system_type, BlockType block_type) const
{
   if(block_type == LINKING_CONS_BLOCK)
      if( !hasLinking(system_type) )
         return;

   SimpleVector* nnz_row = (block_type != LINKING_CONS_BLOCK) ? currNnzRow : currNnzRowLink;
   SimpleVector* iclow = (block_type != LINKING_CONS_BLOCK) ? currIclow : currIclowLink;
   SimpleVector* lhs = (block_type != LINKING_CONS_BLOCK) ? currIneqLhs : currIneqLhsLink;
   SimpleVector* icupp = (block_type != LINKING_CONS_BLOCK) ? currIcupp : currIcuppLink;
   SimpleVector* rhs = (block_type != LINKING_CONS_BLOCK) ? currIneqRhs : currIneqRhsLink;

   if(system_type == EQUALITY_SYSTEM)
      rhs = (block_type != LINKING_CONS_BLOCK) ? currEqRhs : currEqRhsLink;

   assert(nnz_row);
   if(system_type == EQUALITY_SYSTEM)
   {
      assert(rhs); assert(lhs == NULL); assert(iclow == NULL); assert(icupp == NULL);
   }
   else
   {
      assert(lhs); assert(rhs); assert(iclow); assert(icupp); assert( iclow->n == icupp->n );
   }

   for(int i = 0; i < rhs->n; ++i)
   {
      if(nnz_row->elements()[i] != 0.0)
      {
         n_rows++;
         if(nnz_row->elements()[i] == 1.0)
            n_singleton_rows++;

         if(system_type == EQUALITY_SYSTEM)
         {
            n_fixed_rows++;
         }
         else
         {
            if( iclow->elements()[i] != 0.0 && icupp->elements()[i] != 0.0 )
            {
               if( PIPSisEQ(lhs->elements()[i], rhs->elements()[i]))
                  n_fixed_rows++;
               else
                  n_ranged_rows++;
            }
            else
               assert(iclow->elements()[i] != 0.0 || icupp->elements()[i] != 0.0);
         }
      }
   }
}

void StochPresolverBase::countBoxedColumns(int& nBoxCols, int& nColsTotal, int& nFreeVars, BlockType block_type) const
{
   SimpleVector* ixlow = (block_type != LINKING_VARS_BLOCK) ? currIxlowChild : currIxlowParent;
   SimpleVector* ixupp = (block_type != LINKING_VARS_BLOCK) ? currIxuppChild : currIxuppParent;
   SimpleVector* curr_nnz = (block_type != LINKING_VARS_BLOCK) ? currNnzColChild : currNnzColParent;

   assert(curr_nnz); assert(ixlow); assert(ixupp); assert( ixlow->n == ixupp->n );

   for( int i = 0; i < ixlow->n; i++ )
   {
      if( curr_nnz->elements()[i] != 0.0 )
      {
         nColsTotal ++;
         if( ixlow->elements()[i] != 0.0 && ixupp->elements()[i] != 0.0 )
            nBoxCols++;
         else if( ixlow->elements()[i] == 0.0 && ixupp->elements()[i] == 0.0)
            nFreeVars++;
         else
            assert(ixlow->elements()[i] != 0.0 || ixupp->elements()[i] != 0.0);
      }
   }
}

void StochPresolverBase::countSingletonRows(int& n_singletons_equality, int& n_singletons_inequality) const
{
	countSingletonRowsSystem(n_singletons_equality, EQUALITY_SYSTEM);
	countSingletonRowsSystem(n_singletons_inequality, INEQUALITY_SYSTEM);
}

void StochPresolverBase::countSingletonRowsSystem(int& n_singletons, SystemType system_type) const
{
	n_singletons = 0;

	StochVector& nnz_vec = (system_type == EQUALITY_SYSTEM) ? *presData.nRowElemsA : *presData.nRowElemsC;
	assert(nnz_vec.vec);

	/* root node */
	SimpleVector& nnz_b0 = dynamic_cast<SimpleVector&>(*nnz_vec.vec);
	for(long long i = 0; i < nnz_b0.length(); ++i)
		if(nnz_b0[i] == 1.0)
			n_singletons++;

	if(nnz_vec.vecl)
	{
		SimpleVector& nnz_bl = dynamic_cast<SimpleVector&>(*nnz_vec.vecl);
		for(long long i = 0; i < nnz_bl.length(); ++i)
			if(nnz_bl[i] == 1.0)
				n_singletons++;
	}

	for(size_t i = 0; i < nnz_vec.children.size(); ++i)
	{
		SimpleVector& nnz_al = dynamic_cast<SimpleVector&>(*nnz_vec.children[i]->vec);
		for(long long i = 0; i < nnz_al.length(); ++i)
			if(nnz_al[i] == 1.0)
				n_singletons++;
	}
}

// todo verify nonzeros also for dynamic storage
