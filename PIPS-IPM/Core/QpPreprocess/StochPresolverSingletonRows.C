/*
 * StochPresolverSingletonRows.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonRows.h"


StochPresolverSingletonRows::StochPresolverSingletonRows(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverSingletonRows::~StochPresolverSingletonRows()
{
 // todo
}


bool StochPresolverSingletonRows::applyPresolving(int& nelims)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   indivObjOffset = 0.0;
   nelims = 0;
   int newSREq = 0;
   int newSRIneq = 0;

   presData.resetRedCounters();
   newSREq = initSingletonRows(EQUALITY_SYSTEM);
   synchronizeNumberSR(newSREq, newSRIneq);
   if( myRank == 0 ) cout<<"Found "<<newSREq<<" singleton rows in equality system A."<<endl;

   int iter = 0;
   int globalIter = 0;
   bool possibleFeasible = true;

   // main loop:
   while( (newSREq > 0 && iter < maxIterSR) || globalIter == 0 )
   {
      cout<<"Main loop at iter "<<iter<<" and globalIter: "<<globalIter<<endl;
      while( newSREq > 0 && iter < maxIterSR)
      {
         cout<<"SR(Equality) loop at iter "<<iter<<" and globalIter: "<<globalIter<<endl;
         if( globalIter > 0 )
            initSingletonRows(EQUALITY_SYSTEM);
         // main method:
         possibleFeasible = doSingletonRowsA(newSREq, newSRIneq);
         if( !possibleFeasible )
         {
            cout<<"Infeasibility detected: singleton row presolving of A."<<endl;
            return false;
         }
         // update the linking variable blocks (A,C,F,G) with the fixations found in doSingletonRowsA:
         updateLinkingVarsBlocks(newSREq, newSRIneq);
         synchronizeNumberSR(newSREq, newSRIneq);

         if( myRank == 0 )
            cout<<"Found new singleton rows that were just created: "<<newSREq<<" in A, "<<newSRIneq<<" in C."<<endl;
         iter++;
      }
      newSREq = 0;
      if( globalIter == 0 )
      {
         newSRIneq = initSingletonRows(INEQUALITY_SYSTEM);
         synchronizeNumberSR(newSREq, newSRIneq);
         if( myRank == 0 )
            cout<<"Found "<<newSRIneq<<" singleton rows in C."<<endl;
      }
      while( newSRIneq > 0 && iter < maxIterSR)
      {
         cout<<"SR(Inequality) loop at iter "<<iter<<" and globalIter: "<<globalIter<<endl;
         if( globalIter > 0 )
         {
            initSingletonRows(INEQUALITY_SYSTEM);
            // only for debugging:
            synchronizeNumberSR(newSREq, newSRIneq);
            if( myRank == 0 )
               cout<<"Found "<<newSRIneq<<" singleton rows in C."<<endl;
         }
         // main method:
         possibleFeasible = doSingletonRowsC(newSREq, newSRIneq);
         if( !possibleFeasible )
         {
            cout<<"Infeasibility detected: singleton row presolving of C."<<endl;
            return false;
         }
         // update the variable bounds for the linking variables:
         updateLinkingVarsBounds();
         // update the linking variable blocks (A,C,F,G) with the fixations found in doSingletonRowsC:
         updateLinkingVarsBlocks(newSREq, newSRIneq);
         synchronizeNumberSR(newSREq, newSRIneq);

         if( myRank == 0 )
            cout<<"Found new singleton rows that were just created: "<<newSREq<<" in A, "<<newSRIneq<<" in C."<<endl;
         iter++;
      }
      newSRIneq = 0;
      globalIter++;
   }

   // Sum up individual objOffset and then add it to the global objOffset:
   sumIndivObjOffset();
   presData.addObjOffset(indivObjOffset);
   if( myRank == 0 ) cout<<"Global objOffset is now: "<<presData.getObjOffset()<<endl;

   return true;
}

/** Initializes the singletonRows list (acutally a vector<int> singletonRows)
 * and the blocks (int*) pointing to the start and end indices in singletonRows
 * for each block (parent, children, linking rows).
 * Attention: there is no communication over the preocesses to adapt the number of
 * singleton rows found (for performance reasons). If this is necessary, a simple
 * MPI_Allreduce call should be used right after calling this method.
 * Returns the number of singleton rows found (might be different for each process!).
 */
int StochPresolverSingletonRows::initSingletonRows(SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   int nSingletonRows = 0;

   if( system_type == EQUALITY_SYSTEM )
   {
      assert(presData.getNumberSR() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
      const int nSingleRowsA0block =  initSingletonRowsBlock(-1, nRowASimple);
      if( myRank == 0 )
         nSingletonRows += nSingleRowsA0block;

      assert((int)presData.nRowElemsA->children.size() == nChildren);
      for( size_t it = 0; it < presData.nRowElemsA->children.size(); it++)
      {
         SimpleVector* nRowASimpleChild = dynamic_cast<SimpleVector*>(presData.nRowElemsA->children[it]->vec);
         nSingletonRows += initSingletonRowsBlock(int(it), nRowASimpleChild);
      }
      presData.setBlocks(nChildren+1, presData.getNumberSR());

      // todo: linking block nRowElemsA->vecl
      //blocks[nChildren+2] = singletonRows.size();
   }
   else
   {
      assert( system_type == INEQUALITY_SYSTEM );
      //todo: if rows from A and C should be stored, then another variable to store the indices is needed,
      // blocks is not enough.
      assert(presData.getNumberSR() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(presData.nRowElemsC->vec);
      nSingletonRows += initSingletonRowsBlock(-1, nRowASimple);

      assert((int)presData.nRowElemsC->children.size() == nChildren);
      for( size_t it = 0; it < presData.nRowElemsC->children.size(); it++)
      {
         SimpleVector* nRowCSimpleChild = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
         nSingletonRows += initSingletonRowsBlock(int(it), nRowCSimpleChild);
      }
      presData.setBlocks(nChildren+1, presData.getNumberSR());

      // todo: linking block nRowElemsC->vecl
      //blocks[nChildren+2] = singletonRows.size();
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
bool StochPresolverSingletonRows::doSingletonRowsA(int& newSREq, int& newSRIneq)
{
   newSREq = 0;
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

   updateCPForSingletonRow(-1, EQUALITY_SYSTEM);
   bool possFeas = procSingletonRowRoot(matrix, EQUALITY_SYSTEM);
   if( !possFeas ) return false;

   assert(nChildren == (int)matrix.children.size());
   for( int it = 0; it < nChildren; it++ )
   {
      // dummy child?
      if( updateCPForSingletonRow(it, EQUALITY_SYSTEM) )
      {  // main part for each child: go through A and B and adapt F, D and G
         possFeas = procSingletonRowChild( it, newSREq, newSRIneq);
         if( !possFeas ) return false;
      }
   }

   // Update nRowLink and lhs/rhs (Linking part) of both systems:
   updateRhsNRowLink();

   possFeas = presData.combineColAdaptParent();
   if( !possFeas ) return false;

   return true;
}

bool StochPresolverSingletonRows::doSingletonRowsC(int& newSREq, int& newSRIneq)
{
   newSRIneq = 0;
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   updateCPForSingletonRow(-1, INEQUALITY_SYSTEM);
   bool possFeas = procSingletonRowRoot(matrix, INEQUALITY_SYSTEM);
   if( !possFeas ) return false;

   assert(nChildren == (int)matrix.children.size());
   for( int it = 0; it < nChildren; it++ )
   {
      // dummy child?
      if( updateCPForSingletonRow(it, INEQUALITY_SYSTEM) )
      {  // main part for each child: go through A and B and adapt F, D and G ?
         possFeas = procSingletonRowChildInequality(it, newSREq, newSRIneq);
         if( !possFeas ) return false;
      }
   }

   // Update nRowLink and lhs/rhs (Linking part) of both systems:
   updateRhsNRowLink();

   possFeas = combineNewBoundsParent();
   if( !possFeas ) return false;
   possFeas = presData.combineColAdaptParent();
   if( !possFeas ) return false;

   return true;
}

bool StochPresolverSingletonRows::procSingletonRowRoot(StochGenMatrix& stochMatrix, SystemType system_type)
{
   bool possFeas = true;

   SparseStorageDynamic& B0_mat = stochMatrix.Bmat->getStorageDynamicRef();
   assert( presData.getNumberColAdParent() == 0 );

   for(int i = presData.getBlocks(0); i<presData.getBlocks(1); i++)
   {
      int rowIdx = presData.getSingletonRow(i);
      presData.setSingletonRow(i, -1);  // for debugging purposes

      if( system_type == EQUALITY_SYSTEM)
         possFeas = removeSingleRowEntryB0(B0_mat, rowIdx);
      else
      {
         SparseStorageDynamic& B0_trans = stochMatrix.Bmat->getStorageDynamicTransposedRef();
         possFeas = removeSingleRowEntryB0Inequality(B0_mat, B0_trans, rowIdx);
      }
      if( !possFeas ) return false;
   }

   return true;
}

/* Processing the singleton rows in child it, more precisely, goes through all singleton rows in Amat and Bmat.
 * Those in Amat are stored in colAdaptParent for later processing.
 * Those in Bmat are removed and stored in colAdaptLinkBlock. Furthermore, the corresponding fixed variables (columns)
 * in Bmat and in Blmat are removed.
 * Using this colAdaptLinkBlock, the variables (columns) are removed from the inequalities Bmat, Blmat as well.
 */
bool StochPresolverSingletonRows::procSingletonRowChild(int it, int& newSREq, int& newSRIneq)
{
   bool possFeas = procSingletonRowChildAmat( it, EQUALITY_SYSTEM);
   if( !possFeas ) return false;

   std::vector<COLUMNTOADAPT> colAdaptLinkBlock;
   possFeas = procSingletonRowChildBmat(it, colAdaptLinkBlock, newSREq, EQUALITY_SYSTEM);
   if( !possFeas ) return false;

   // using colAdaptLinkBlock, go through the columns in Blmat
   if( hasLinking(EQUALITY_SYSTEM) )
   {
      possFeas = adaptChildBlmat( colAdaptLinkBlock, EQUALITY_SYSTEM);
      if( !possFeas ) return false;
   }
   // and go through the columns in Bmat, Blmat of the inequality
   updateCPForSingletonRowInequalityBChild( it );
   possFeas = adaptOtherSystemChildB( INEQUALITY_SYSTEM, colAdaptLinkBlock, newSRIneq );
   if( !possFeas ) return false;

   return true;
}

bool StochPresolverSingletonRows::procSingletonRowChildAmat(int it, SystemType system_type)
{
   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();
   double* g = currgParent->elements();

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

            if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
                  || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
            {
               cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<", child="<<it<<endl;
               return false;
            }
            if( !storeColValInColAdaptParentAndAdaptOffset(colIdx, val, g) )
               return false;
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
               return false;

            // test if they imply fixation
            else if( newBoundsFixVariable(val, newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
            {
               cout<<"New bounds imply fixation of variable "<<colIdx<<" of child "<<it<<" to value: "<<val<<endl;
               // as in SR(equality), store them to remove the column later
               if( !storeColValInColAdaptParentAndAdaptOffset(colIdx, val, g) )
                  return false;

               // nnz/red Counters are not touched yet, they will be set later when colAdaptParent is applied.
            }
            else
            {
               // test if new bounds are tightening: add to newBoundsParent
               if( newBoundsTightenOldBounds(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
               {
                  cout<<"New bounds tighten bounds of variable "<<colIdx<<endl;
                  // store them to adapt the bounds on all processes later
                  if( !storeNewBoundsParent(colIdx, newxlow, newxupp) )
                     return false;
               }
               else
                  cout<<"New bounds are redundant for variable "<<colIdx<<endl;

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
   return true;
}

bool StochPresolverSingletonRows::procSingletonRowChildBmat(int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock,
      int& newSR, SystemType system_type)
{
   assert(currBmat != NULL);
   cout<<"procSingletonRowChildBmat for child "<<it<<" and system_type "<<system_type<<endl;
   for(int i = presData.getBlocks(it+1); i<presData.getBlocks(it+2); i++)
   {
      const int rowIdx = presData.getSingletonRow(i);
      if( rowIdx == -1 || currBmat->rowptr[rowIdx].start == currBmat->rowptr[rowIdx].end)
         continue;   // entry was already in Amat or in a previous singleton row in Bmat
      else
      {
         assert( currBmat->rowptr[rowIdx].start +1 == currBmat->rowptr[rowIdx].end );
         presData.setSingletonRow(i, -1);  // for debugging purposes
         removeSingleRowEntryChildBmat(rowIdx, colAdaptLinkBlock, system_type, newSR);
      }
   }
   return true;
}

bool StochPresolverSingletonRows::removeSingleRowEntryChildBmat(int rowIdx,
      std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR)
{
   double* ixlow = currIxlowChild->elements();
   double* ixupp = currIxuppChild->elements();
   double* xlow = currxlowChild->elements();
   double* xupp = currxuppChild->elements();
   double* g = currgChild->elements();

   int colIdx = -1;
   double aik = 0.0;
   getValuesForSR(*currBmat, rowIdx, colIdx, aik);

   if( system_type == EQUALITY_SYSTEM )
   {
      const double val = currEqRhs->elements()[rowIdx] / aik;

      if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
            || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
      {
         cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
         return false;
      }
      indivObjOffset += g[colIdx] * val;

      // adapt the col, val immediately in this block B_i and store them for the Blmat
      newSR += adaptChildBmatCol(colIdx, val, system_type);

      const COLUMNTOADAPT colWithVal = {colIdx, val};
      colAdaptLinkBlock.push_back(colWithVal);
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
         return false;

      // test if they imply fixation
      else if( newBoundsFixVariable(val, newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
      {
         cout<<"New bounds imply fixation of variable "<<colIdx<<" of a child to value: "<<val<<endl;
         // adapt immediately in D_i by adapting objOffset, removing column colIdx and store in colADaptLinkBlock for G, F, B
         indivObjOffset += g[colIdx] * val;

         // adapt the col, val immediately in this block B_i and store them for the Blmat
         newSR += adaptChildBmatCol(colIdx, val, INEQUALITY_SYSTEM);

         COLUMNTOADAPT colWithVal = {colIdx, val};
         colAdaptLinkBlock.push_back(colWithVal);
      }
      else
      {
         // test if new bounds are tightening:
         if( newBoundsTightenOldBounds(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
         {
            cout<<"New bounds tighten bounds of variable "<<colIdx<<endl;
            // adapt immediately the variable bounds
            setNewXBounds(colIdx, newxlow, newxupp, ixlow, xlow, ixupp, xupp);
         }
         else
            cout<<"New bounds are redundant for variable "<<colIdx<<endl;

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
   return true;
}

/** Removes the single entry in row rowIdx in Bmat_0.
 * Stores the corresponding column index in colAdaptParent.
 * Return false if infeasibility is detected.*/
bool StochPresolverSingletonRows::removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();
   double* g = currgParent->elements();

   int colIdx = -1;
   double aik = 0.0;
   getValuesForSR(storage, rowIdx, colIdx, aik);

   const double val = currEqRhs->elements()[rowIdx] / aik;

   if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
         || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
   {
      cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
      return false;
   }
   else
   {
      if( myRank == 0 )
      {
         if( !storeColValInColAdaptParentAndAdaptOffset(colIdx, val, g) )
            return false;
      }
   }

   return true;
}

bool StochPresolverSingletonRows::removeSingleRowEntryB0Inequality(SparseStorageDynamic& storage,
      SparseStorageDynamic& storageTransposed, int rowIdx)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();
   double* g = currgParent->elements();

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
      return false;

   // test if they imply fixation
   else if( newBoundsFixVariable(val, newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
   {
      cout<<"New bounds imply fixation of variable "<<colIdx<<" to "<<val<<endl;
      // as in SR(equality), store them to remove the column later
      if( myRank == 0 )
      {
         if( !storeColValInColAdaptParentAndAdaptOffset(colIdx, val, g) )
            return false;
      }
      // in case of fixation, nnz bzw. red Counters are not touched yet because they will be set
      // correctly later, when colAdaptParent is applied.
   }
   else
   {
      // test if new bounds are tightening: add to newBoundsParent
      if( newBoundsTightenOldBounds(newxlow, newxupp, colIdx, ixlow, ixupp, xlow, xupp) )
      {
         cout<<"New bounds tighten bounds of variable "<<colIdx<<endl;
         // store them to adapt the bounds on all processes later
         if( myRank == 0 )
         {
            if( !storeNewBoundsParent(colIdx, newxlow, newxupp) )
               return false;
         }
      }
      else
         cout<<"New bounds are redundant for variable "<<colIdx<<endl;

      // set a_ik=0.0, nRow=0, nCol--
      clearRow(storage, rowIdx);
      // remove entry a_ik in transposed matrix as well
      removeEntryInDynamicStorage(storageTransposed, colIdx, rowIdx, val);
      currNnzRow->elements()[rowIdx]--;
      currNnzColParent->elements()[colIdx]--;
      assert( currNnzRow->elements()[rowIdx] == 0 );
   }

   return true;
}

bool StochPresolverSingletonRows::procSingletonRowChildInequality(int it, int& newSREq, int& newSRIneq)
{
   bool possFeas = true;

   // go through A, storing new bounds in newBoundsParent and possibly fixations in colAdaptParent
   // go through B, adapting new bounds immediately and storing fixations in colAdaptLinkBlock

   possFeas = procSingletonRowChildAmat(it, INEQUALITY_SYSTEM);
   if( !possFeas ) return false;

   std::vector<COLUMNTOADAPT> colAdaptLinkBlock;
   possFeas = procSingletonRowChildBmat( it, colAdaptLinkBlock, newSRIneq, INEQUALITY_SYSTEM);
   if( !possFeas ) return false;

   // using colAdaptLinkBlock, go through the columns in Blmat
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      possFeas = adaptChildBlmat( colAdaptLinkBlock, INEQUALITY_SYSTEM);
      if( !possFeas ) return false;
   }

   // and go through the columns in Bmat, Blmat of the equality system
   updateCPForSingletonRowEqualityBChild( it );
   possFeas = adaptOtherSystemChildB( EQUALITY_SYSTEM, colAdaptLinkBlock, newSREq );
   if( !possFeas ) return false;

   return true;
}

void StochPresolverSingletonRows::calculateNewBoundsOnVariable(double& newxlow, double& newxupp, int rowIdx, double aik) const
{
   if( aik > 0.0 )
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

bool StochPresolverSingletonRows::newBoundsImplyInfeasible(double newxlow, double newxupp, int colIdx,
      double* ixlow, double* ixupp, double* xlow, double* xupp) const
{
   if( ( ixlow[colIdx] != 0.0 && xlow[colIdx] > newxupp)
         || (ixupp[colIdx] != 0.0 && xupp[colIdx] < newxlow )
         || (newxlow > newxupp))
   {
      cout<<"Infeasibility detected at variable "<<colIdx<<", new bounds= ["<<newxlow<<", "<<newxupp<<"]"<<endl;
      return true;
   }
   return false;
}

bool StochPresolverSingletonRows::newBoundsFixVariable(double& value, double newxlow, double newxupp, int colIdx,
      double* ixlow, double* ixupp, double* xlow, double* xupp) const
{
   if( newxlow == newxupp
         || ( ixlow[colIdx] != 0.0 && xlow[colIdx] == newxupp ))
   {
      value = newxupp;
      return true;
   }
   else if( ixupp[colIdx] != 0.0 && xupp[colIdx] == newxlow )
   {
      value = newxlow;
      return true;
   }
   return false;
}

bool StochPresolverSingletonRows::newBoundsTightenOldBounds(double newxlow, double newxupp, int colIdx,
      double* ixlow, double* ixupp, double* xlow, double* xupp) const
{
   if( ( ixlow[colIdx] != 0.0 && newxlow > xlow[colIdx] )
         || ( ixlow[colIdx] == 0.0 && newxlow > -std::numeric_limits<double>::max() )
         || ( ixupp[colIdx] != 0.0 && newxupp < xupp[colIdx] )
         || ( ixupp[colIdx] == 0.0 && newxupp < std::numeric_limits<double>::max() ) )
      return true;
   return false;
}

/** Stores the column index colIdx together with the value as a COLUMNTOADAPT in colAdaptParent.
 * Adapts the objective offset g only once for each column (variable).
 * Returns false if infeasibility is detected.
 */
bool StochPresolverSingletonRows::storeColValInColAdaptParentAndAdaptOffset(int colIdx, double value, double* g)
{
   const COLUMNTOADAPT colWithVal = {colIdx, value};
   bool uniqueAdditionToOffset = true;
   for(int i=0; i<presData.getNumberColAdParent(); i++)
   {
      if( presData.getColAdaptParent(i).colIdx == colIdx )
      {
         if( presData.getColAdaptParent(i).val != value )
         {
            cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<value<<endl;
            return false;
         }
         uniqueAdditionToOffset = false;
      }
   }
   if( uniqueAdditionToOffset )
   {
      presData.addColToAdaptParent(colWithVal);
      indivObjOffset += g[colIdx] * value;
   }
   return true;
}

/** Stores the column index colIdx together with the new bounds as a XBOUNDS in newBoundsParent.
 * Should be called only from Process Zero.
 * Returns false if infeasibility is detected (contradictory bounds).
 */
bool StochPresolverSingletonRows::storeNewBoundsParent(int colIdx, double newxlow, double newxupp)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   assert( myRank == 0 );

   XBOUNDS newXbounds = {colIdx, newxlow, newxupp};
   for(int i=0; i<(int)newBoundsParent.size(); i++)
   {
      if( newBoundsParent[i].colIdx == colIdx )
      {
         if( newBoundsParent[i].newxlow > newxlow || newBoundsParent[i].newxupp < newxlow )
         {
            cout<<"Infeasibility detected at variable "<<colIdx<<" because of tightened bounds."<<endl;
            return false;
         }
      }
   }
   newBoundsParent.push_back(newXbounds);
   return true;
}

/** Method similar to combineColAdaptParent(), that is a method going through newBoundsParent
 * and cleaning it up, removing redundant bounds, checking for infeasibility or more tightening.
 */
bool StochPresolverSingletonRows::combineNewBoundsParent()
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
               cout<<"Detected infeasibility in variable "<<colIdxCurrent<<" of parent."<<endl;
               return false;
            }
            else
            {
               // change the vector element newBoundsParent.begin()+(i-1), also das,
               // welches colIdxCurrent definiert hat:
               setNewBoundsParent(i-1, colIdxCurrent, bestLow, bestUpp);
               newBoundsParent.erase(newBoundsParent.begin()+i);   //todo: implement more efficiently
            }
         }
         else
         {
            colIdxCurrent = getNewBoundsParent(i).colIdx;
            xlowCurrent = getNewBoundsParent(0).newxlow;
            xuppCurrent = getNewBoundsParent(0).newxupp;
         }
      }
   }
   assert( getNumberNewBoundsParent() <= presData.nColElems->vec->n );

   return true;
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
      setNewXBounds(newbounds.colIdx, newbounds.newxlow, newbounds.newxupp, ixlow, xlow, ixupp, xupp);
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
   cout<<"a_ik = "<<aik<<" at ("<<rowIdx<<" ,"<<colIdx<<" )" <<endl;
}

XBOUNDS StochPresolverSingletonRows::getNewBoundsParent(int i) const
{
   assert( i<getNumberNewBoundsParent() );
   return newBoundsParent[i];
}
void StochPresolverSingletonRows::setNewBoundsParent(int i, int colIdx, double newxlow, double newxupp)
{
   assert( i<getNumberNewBoundsParent() );
   newBoundsParent[i].colIdx = colIdx;
   newBoundsParent[i].newxlow = newxlow;
   newBoundsParent[i].newxupp = newxupp;
}
int StochPresolverSingletonRows::getNumberNewBoundsParent() const
{
   return (int)newBoundsParent.size();
}
void StochPresolverSingletonRows::addNewBoundsParent(XBOUNDS newXBounds)
{
   newBoundsParent.push_back(newXBounds);
}
void StochPresolverSingletonRows::clearNewBoundsParent()
{
   newBoundsParent.clear();
}

void StochPresolverSingletonRows::setNewXBounds(int colIdx, double newxlow, double newxupp,
      double* ixlow, double* xlow, double* ixupp, double* xupp) const
{
   if( (ixlow[colIdx] != 0.0 && newxlow > xlow[colIdx])
      || (ixlow[colIdx] == 0.0 && newxlow > -std::numeric_limits<double>::max()) )
   {
      ixlow[colIdx] = 1.0;
      xlow[colIdx] = newxlow;
   }
   if( (ixupp[colIdx] != 0.0 && newxupp < xupp[colIdx])
         || (ixupp[colIdx] == 0.0 && newxupp < std::numeric_limits<double>::max()))
   {
      ixupp[colIdx] = 1.0;
      xupp[colIdx] = newxupp;
   }
}

void StochPresolverSingletonRows::synchronizeNumberSR(int& newSREq, int& newSRIneq) const
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );
   if( iAmDistrib )
   {
      int* newSR = new int[2];
      newSR[0] = newSREq;
      newSR[1] = newSRIneq;
      MPI_Allreduce(MPI_IN_PLACE, newSR, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      newSREq = newSR[0];
      newSRIneq = newSR[1];
      delete[] newSR;
   }
}



