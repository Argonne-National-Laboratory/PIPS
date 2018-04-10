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
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   nelims = 0;
   int newSRIneq = 0;
   int newSREq = initSingletonRows(EQUALITY_SYSTEM);
   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &newSREq, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if( myRank == 0 ) cout<<"Found "<<newSREq<<" singleton rows in equality system A."<<endl;

   int iter = 0;
   bool possibleFeasible = true;

   // main loop:
   while( newSREq > 0 && iter < maxIterSR)
   {
      if( iter > 0 )
         initSingletonRows(EQUALITY_SYSTEM);
      presData.resetRedCounters();
      // main method:
      possibleFeasible = doSingletonRowsA(newSREq, newSRIneq);
      if( !possibleFeasible )
      {
         cout<<"Infeasibility detected during singleton row presolving."<<endl;
         return 0;
      }
      // update the linking variable blocks (A,C,F,G) with the fixations found in doSingletonRowsA:
      updateLinkingVarsBlocks(newSREq, newSRIneq);
      //nelims += doSingletonRowsC();
      if( myRank == 0 )
         cout<<"Found new singleton rows that were just created: "<<newSREq<<" in A, "<<newSRIneq<<" in C."<<endl;
      iter++;
   }

   presData.globalSumObjOffset();

   return nelims;
}

int StochPresolverSingletonRows::initSingletonRows(SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   int nSingletonRows = 0;

   if( system_type == EQUALITY_SYSTEM )
   {
      assert(presData.getNumberSR() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(presData.nRowElemsA->vec);
      int nSingleRowsA0block =  initSingletonRowsBlock(-1, nRowASimple);
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
         SimpleVector* nRowASimpleChild = dynamic_cast<SimpleVector*>(presData.nRowElemsC->children[it]->vec);
         nSingletonRows += initSingletonRowsBlock(int(it), nRowASimpleChild);
      }
      presData.setBlocks(nChildren+1, presData.getNumberSR());

      // todo: linking block nRowElemsC->vecl
      //blocks[nChildren+2] = singletonRows.size();
   }

   return nSingletonRows;
}

int StochPresolverSingletonRows::initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple)
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

bool StochPresolverSingletonRows::doSingletonRowsA(int& newSREq, int& newSRIneq)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   newSREq = 0;
   newSRIneq = 0;
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

   updateCPForSingletonRow(-1, EQUALITY_SYSTEM);
   bool possFeas = procSingletonRowRoot(matrix);
   if( !possFeas ) return false;

   assert(nChildren == (int)matrix.children.size());
   for( int it = 0; it < nChildren; it++ )
   {
      // dummy child?
      if( updateCPForSingletonRow(it, EQUALITY_SYSTEM) )
      {  // main part for each child: go through B and adapt F, D and G
         possFeas = procSingletonRowChild(matrix, it, newSREq, newSRIneq);
         if( !possFeas ) return false;
      }
   }

   // Update nRowLink and lhs/rhs (Linking part) of both systems:
   updateRhsNRowLink();

   possFeas = presData.combineColAdaptParent();
   if( !possFeas ) return false;

   return true;
}

int StochPresolverSingletonRows::doSingletonRowsC()
{
   int nelims = 0;
   return nelims;
}

bool StochPresolverSingletonRows::procSingletonRowRoot(StochGenMatrix& stochMatrix)
{
   bool possFeas = true;

   SparseStorageDynamic& B0_mat = stochMatrix.Bmat->getStorageDynamicRef();
   assert( presData.getNumberColAdParent() == 0 );

   for(int i = presData.getBlocks(0); i<presData.getBlocks(1); i++)
   {
      int rowIdx = presData.getSingletonRow(i);
      presData.setSingletonRow(i, -1);  // for debugging purposes

      possFeas = removeSingleRowEntryB0(B0_mat, rowIdx);
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
bool StochPresolverSingletonRows::procSingletonRowChild(StochGenMatrix& stochMatrix, int it, int& newSREq, int& newSRIneq)
{
   bool possFeas = true;

   SparseStorageDynamic& A_mat = stochMatrix.children[it]->Amat->getStorageDynamicRef();
   SparseStorageDynamic& B_mat = stochMatrix.children[it]->Bmat->getStorageDynamicRef();

   possFeas = procSingletonRowChildAmat( A_mat, it);
   if( !possFeas ) return false;

   std::vector<COLUMNTOADAPT> colAdaptLinkBlock;
   possFeas = procSingletonRowChildBmat( B_mat, it, colAdaptLinkBlock, newSREq);
   if( !possFeas ) return false;

   // using colAdaptLinkBlock, go through the columns in Blmat
   if( hasLinking(EQUALITY_SYSTEM) )
   {
      possFeas = adaptChildBlmat( colAdaptLinkBlock, EQUALITY_SYSTEM);
      if( !possFeas ) return false;
   }
   // and go through the columns in Bmat, Blmat of the inequality
   updateCPForSingletonRowInequalityBChild( it );
   possFeas = adaptInequalityChildB( colAdaptLinkBlock, newSRIneq );
   if( !possFeas ) return false;

   return true;
}

bool StochPresolverSingletonRows::procSingletonRowChildAmat(SparseStorageDynamic& A_mat, int it)
{
   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();
   double* g = currgParent->elements();

   for(int i = presData.getBlocks(it+1); i<presData.getBlocks(it+2); i++)
   {
      int rowIdx = presData.getSingletonRow(i);
      if( A_mat.rowptr[rowIdx].start +1 == A_mat.rowptr[rowIdx].end )
      {
         presData.setSingletonRow(i, -1);  // for debugging purposes

         // store the column index with fixed value in colAdaptParent and adapt objOffset:
         int indexK = A_mat.rowptr[rowIdx].start;
         int colIdx = A_mat.jcolM[indexK];

         assert(A_mat.rowptr[rowIdx].start +1 == A_mat.rowptr[rowIdx].end);

         double aik = A_mat.M[indexK];
         assert(aik != 0.0);
         cout<<"a_ik = "<<aik<<" at ("<<rowIdx<<" ,"<<colIdx<<" )" <<endl;

         double val = currEqRhs->elements()[rowIdx] / aik;

         if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
               || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
         {
            cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
            return false;
         }
         presData.addObjOffset(g[colIdx] * val);

         COLUMNTOADAPT colWithVal = {colIdx, val};
         presData.addColToAdaptParent(colWithVal);
      }
   }
   return true;
}

bool StochPresolverSingletonRows::procSingletonRowChildBmat(SparseStorageDynamic& B_mat, int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR)
{
   for(int i = presData.getBlocks(it+1); i<presData.getBlocks(it+2); i++)
   {
      int rowIdx = presData.getSingletonRow(i);
      if( rowIdx == -1 || B_mat.rowptr[rowIdx].start == B_mat.rowptr[rowIdx].end)
      {
         // entry was already in Amat or in a previous singleton row in Bmat
         continue;
      }
      else
      {
         assert( B_mat.rowptr[rowIdx].start +1 == B_mat.rowptr[rowIdx].end );
         presData.setSingletonRow(i, -1);  // for debugging purposes
         removeSingleRowEntryChildBmat(rowIdx, colAdaptLinkBlock, EQUALITY_SYSTEM, newSR);
      }
   }

   return true;
}

bool StochPresolverSingletonRows::removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR)
{
   double* ixlow = currIxlowChild->elements();
   double* ixupp = currIxuppChild->elements();
   double* xlow = currxlowChild->elements();
   double* xupp = currxuppChild->elements();
   double* g = currgChild->elements();

   assert(currBmat->rowptr[rowIdx].start +1 == currBmat->rowptr[rowIdx].end);

   int indexK = currBmat->rowptr[rowIdx].start;
   int colIdx = currBmat->jcolM[indexK];

   double aik = currBmat->M[indexK];
   assert(aik != 0.0);
   cout<<"a_ik = "<<aik<<" at ("<<rowIdx<<" ,"<<colIdx<<" )" <<endl;

   double val = currEqRhs->elements()[rowIdx] / aik;

   if( (ixlow[colIdx] != 0.0 && xlow[colIdx] > val)
         || (ixupp[colIdx] != 0.0 && xupp[colIdx] < val))
   {
      cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
      return false;
   }
   presData.addObjOffset(g[colIdx] * val);

   // adapt the col, val immediately in this block B_i and store them for the Blmat
   newSR += adaptChildBmatCol(colIdx, val, system_type);

   COLUMNTOADAPT colWithVal = {colIdx, val};
   colAdaptLinkBlock.push_back(colWithVal);

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

   int indexK = storage.rowptr[rowIdx].start;
   int colIdx = storage.jcolM[indexK];

   assert(storage.rowptr[rowIdx].start +1 == storage.rowptr[rowIdx].end);

   double aik = storage.M[indexK];
   assert(aik != 0.0);
   cout<<"a_ik = "<<aik<<" at ("<<rowIdx<<" ,"<<colIdx<<" )" <<endl;

   double val = currEqRhs->elements()[rowIdx] / aik;

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
         COLUMNTOADAPT colWithVal = {colIdx, val};
         bool uniqueAdditionToOffset = true;
         for(int i=0; i<presData.getNumberColAdParent(); i++)
         {
            if( presData.getColAdaptParent(i).colIdx == colIdx )
            {
               if( presData.getColAdaptParent(i).val != val )
               {
                  cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
                  return false;
               }
               uniqueAdditionToOffset = false;
            }

         }
         if( uniqueAdditionToOffset )
         {
            presData.addColToAdaptParent(colWithVal);
            presData.addObjOffset(g[colIdx] * val);
         }
      }
   }

   return true;
}



