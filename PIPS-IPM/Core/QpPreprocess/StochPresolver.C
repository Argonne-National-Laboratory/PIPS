/*
 * StochPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */


#include "StochPresolver.h"
#include <cassert>
#include <iostream>
#include <cmath>
#include <utility>
#include <math.h>
#include <algorithm>

#include "../Abstract/DoubleMatrix.h"
#include "../SparseLinearAlgebra/SparseGenMatrix.h"
#include "../StochLinearAlgebra/StochVectorHandle.h"
#include "../Vector/OoqpVector.h"
#include "StochGenMatrix.h"
#include "sTreeCallbacks.h"
#include "DoubleMatrixTypes.h"
#include "StochPresolverTinyEntries.h"

StochPresolver::StochPresolver(const Data* prob)
 : QpPresolver(prob)
{
   const sData* sorigprob = dynamic_cast<const sData*>(prob);
}


StochPresolver::~StochPresolver()
{
   delete[] blocks;
   delete[] blocksIneq;
}

Data* StochPresolver::presolve()
{

   std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);


   // clone and initialize dynamic storage

   int myRank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   ofstream myfile;
   myfile.open ("before.txt");
   sorigprob->writeToStreamDense(myfile);
   myfile.close();


   int nelims = 0;
// init presolve data
   PresolveData presData(origprob);

   presData.initialize();

   // init presolvers

   StochPresolverTinyEntries presolverTiny(presData);


   // main presolving loop
   // while ( rerun )

   bool infeasible = presolverTiny.applyPresolving(nelims);

  // if( infeasible )
  //    break;


   removedEntries.reserve((int)ceil(0.2 * ((*presProb).my + (*presProb).mz) * (*presProb).nx));
   //int nElimsTinyEntries = removeTinyEntries();
   //if( myRank == 0)
   //   std::cout << "In total, "<<nElimsTinyEntries<<" tiny entries were removed." << std::endl;

   //doSingletonRows();

   myfile.open("middle.txt");
   presProb->writeToStreamDense(myfile);
   myfile.close();

/*   cout<<"nRowElemsA "<<endl;
   nRowElemsA->writeToStreamAll(cout);
   cout<<"nRowElemsC "<<endl;
   nRowElemsC->writeToStreamAll(cout);
   cout<<"nColElems "<<endl;
   nColElems->writeToStreamAll(cout);
*/

   myfile.open("after.txt");
   presProb->writeToStreamDense(myfile);
   myfile.close();

   //presProb->writeToStreamDense(std::cout);

   std::cout << "sorigprob nx, my, mz" << sorigprob->nx << " " << sorigprob->my << " " << sorigprob->mz << std::endl;
   std::cout << "presProb nx, my, mz" << presProb->nx << " " << presProb->my << " " << presProb->mz << std::endl;





   //assert(0);

   return presData.finalize();
}

int StochPresolver::doSingletonRows()
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   int nelims = 0, newSRIneq = 0;
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
      resetRedCounters();
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

   updateObjOffset();

   return nelims;
}

int StochPresolver::initSingletonRows(SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   int nSingletonRows = 0;

   if( system_type == EQUALITY_SYSTEM )
   {
      assert(singletonRows.size() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(nRowElemsA->vec);
      int nSingleRowsA0block =  initSingletonRowsBlock(-1, nRowASimple);
      if( myRank == 0 )
         nSingletonRows += nSingleRowsA0block;

      assert((int)nRowElemsA->children.size() == nChildren);
      for( size_t it = 0; it < nRowElemsA->children.size(); it++)
      {
         SimpleVector* nRowASimpleChild = dynamic_cast<SimpleVector*>(nRowElemsA->children[it]->vec);
         nSingletonRows += initSingletonRowsBlock(int(it), nRowASimpleChild);
      }
      blocks[nChildren+1] = singletonRows.size();

      // todo: linking block nRowElemsA->vecl
      //blocks[nChildren+2] = singletonRows.size();
   }
   else
   {
      assert( system_type == INEQUALITY_SYSTEM );
      //todo: if rows from A and C should be stored, then another variable to store the indices is needed,
      // blocks is not enough.
      assert(singletonRows.size() == 0);

      SimpleVector* nRowASimple = dynamic_cast<SimpleVector*>(nRowElemsC->vec);
      nSingletonRows += initSingletonRowsBlock(-1, nRowASimple);

      assert((int)nRowElemsC->children.size() == nChildren);
      for( size_t it = 0; it < nRowElemsC->children.size(); it++)
      {
         SimpleVector* nRowASimpleChild = dynamic_cast<SimpleVector*>(nRowElemsC->children[it]->vec);
         nSingletonRows += initSingletonRowsBlock(int(it), nRowASimpleChild);
      }
      blocks[nChildren+1] = singletonRows.size();

      // todo: linking block nRowElemsC->vecl
      //blocks[nChildren+2] = singletonRows.size();
   }

   return nSingletonRows;
}

int StochPresolver::initSingletonRowsBlock(int it, SimpleVector* nnzRowSimple)
{
   int nSingletonRows = 0;

   blocks[it + 1] = singletonRows.size();
   double* nnzRow = nnzRowSimple->elements();

   for( int i = 0; i < nnzRowSimple->n; i++ )
      if( nnzRow[i] == 1.0 )
      {
         singletonRows.push_back(i);
         nSingletonRows++;
      }
   return nSingletonRows;
}

bool StochPresolver::doSingletonRowsA(int& newSREq, int& newSRIneq)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   newSREq = 0;
   newSRIneq = 0;
   StochGenMatrix& matrix = dynamic_cast<StochGenMatrix&>(*(presProb->A));

   updateCurrentPointersForSingletonRow(-1, EQUALITY_SYSTEM);
   bool possFeas = procSingletonRowRoot(matrix);
   if( !possFeas ) return false;

   assert(nChildren == (int)matrix.children.size());
   for( int it = 0; it < nChildren; it++ )
   {
      // dummy child?
      if( updateCurrentPointersForSingletonRow(it, EQUALITY_SYSTEM) )
      {  // main part for each child: go through B and adapt F, D and G
         possFeas = procSingletonRowChild(matrix, it, newSREq, newSRIneq);
         if( !possFeas ) return false;
      }
   }

   // Update nRowLink and lhs/rhs (Linking part) of both systems:
   updateRhsNRowLink();

   possFeas = combineColAdaptParent();
   if( !possFeas ) return false;

   return true;
}


/* Processing the singleton rows in child it, more precisely, goes through all singleton rows in Amat and Bmat.
 * Those in Amat are stored in colAdaptParent for later processing.
 * Those in Bmat are removed and stored in colAdaptLinkBlock. Furthermore, the corresponding fixed variables (columns)
 * in Bmat and in Blmat are removed.
 * Using this colAdaptLinkBlock, the variables (columns) are removed from the inequalities Bmat, Blmat as well.
 */
bool StochPresolver::procSingletonRowChild(StochGenMatrix& stochMatrix, int it, int& newSREq, int& newSRIneq)
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

bool StochPresolver::procSingletonRowChildAmat(SparseStorageDynamic& A_mat, int it)
{
   double* ixlow = currIxlowParent->elements();
   double* ixupp = currIxuppParent->elements();
   double* xlow = currxlowParent->elements();
   double* xupp = currxuppParent->elements();
   double* g = currgParent->elements();

   for(int i = blocks[it+1]; i<blocks[it+2]; i++)
   {
      int rowIdx = singletonRows[i];
      if( A_mat.rowptr[rowIdx].start +1 == A_mat.rowptr[rowIdx].end )
      {
         singletonRows[i] = -1;  // for debugging purposes

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
         objOffset += g[colIdx] * val;

         COLUMNTOADAPT colWithVal = {colIdx, val};
         colAdaptParent.push_back(colWithVal);
      }
   }
   return true;
}

bool StochPresolver::procSingletonRowChildBmat(SparseStorageDynamic& B_mat, int it, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, int& newSR)
{
   for(int i = blocks[it+1]; i<blocks[it+2]; i++)
   {
      int rowIdx = singletonRows[i];
      if( rowIdx == -1 || B_mat.rowptr[rowIdx].start == B_mat.rowptr[rowIdx].end)
      {
         // entry was already in Amat or in a previous singleton row in Bmat
         continue;
      }
      else
      {
         assert( B_mat.rowptr[rowIdx].start +1 == B_mat.rowptr[rowIdx].end );
         singletonRows[i] = -1;  // for debugging purposes
         removeSingleRowEntryChildBmat(rowIdx, colAdaptLinkBlock, EQUALITY_SYSTEM, newSR);
      }
   }

   return true;
}

bool StochPresolver::removeSingleRowEntryChildBmat( int rowIdx, std::vector<COLUMNTOADAPT> & colAdaptLinkBlock, SystemType system_type, int& newSR)
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
   objOffset += g[colIdx] * val;

   // adapt the col, val immediately in this block B_i and store them for the Blmat
   newSR += adaptChildBmatCol(colIdx, val, system_type);

   COLUMNTOADAPT colWithVal = {colIdx, val};
   colAdaptLinkBlock.push_back(colWithVal);

   return true;
}

/** For the given column index and the value to which this variable is to be fixed,
 * go through the current Bmat by going through the corresponding row of the transposed Bmat.
 * Adapt the rhs, the nnzRow, the nnzCol and all concerned entries in Bmat and BmatTransposed.
 * Return the number of newly found singleton rows.
 */
int StochPresolver::adaptChildBmatCol(int colIdx, double val, SystemType system_type)
{
   int newSingletonRows = 0;

   for( int j = currBmatTrans->rowptr[colIdx].start; j<currBmatTrans->rowptr[colIdx].end; j++ )
   {
      int rowIdxB = currBmatTrans->jcolM[j];
      double m = removeEntryInDynamicStorage(*currBmat, rowIdxB, colIdx);

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
bool StochPresolver::adaptInequalityChildB( std::vector<COLUMNTOADAPT> & colAdaptBblock, int& newSRIneq )
{
   // Bmat blocks
   bool possFeas = adaptChildBmat( colAdaptBblock, INEQUALITY_SYSTEM, newSRIneq);
   if( !possFeas ) return false;

   // Blmat blocks
   if( hasLinking(INEQUALITY_SYSTEM) )
   {
      bool possFeas = adaptChildBlmat( colAdaptBblock, INEQUALITY_SYSTEM);
      if( !possFeas ) return false;
   }

   return true;
}
/** Removes the single entry in row rowIdx in Bmat_0.
 * Stores the corresponding column index in colAdaptParent.
 * Return false if infeasibility is detected.*/
bool StochPresolver::removeSingleRowEntryB0(SparseStorageDynamic& storage, int rowIdx)
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
         for(int i=0; i<(int)colAdaptParent.size(); i++)
         {
            if( colAdaptParent[i].colIdx == colIdx )
            {
               if( colAdaptParent[i].val != val )
               {
                  cout<<"Infeasibility detected at variable "<<colIdx<<", val= "<<val<<endl;
                  return false;
               }
               uniqueAdditionToOffset = false;
            }

         }
         if( uniqueAdditionToOffset )
         {
            colAdaptParent.push_back(colWithVal);
            objOffset += g[colIdx] * val;
         }
      }
   }

   return true;
}


/** Adapt the columns for the linking-variable-blocks (the A_i) blocks */
int StochPresolver::colAdaptLinkVars(int it, SystemType system_type)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   int newSingletonRows = 0;

   for( int i=0; i < (int)colAdaptParent.size(); i++)
   {
      int colIdxA = colAdaptParent[i].colIdx;
      double val = colAdaptParent[i].val;

      for( int j = currAmatTrans->rowptr[colIdxA].start; j < currAmatTrans->rowptr[colIdxA].end; j++ )
      {
         int rowIdxA = currAmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currAmat, rowIdxA, colIdxA);
         cout<<"Removed entry "<<rowIdxA<<", "<<colIdxA<<" with value "<<m<<" in Amat of child "<<it<<endl;

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

int StochPresolver::colAdaptF0(SystemType system_type)
{
   assert(currBlmat != NULL);
   assert( currNnzRow->n == currBlmat->m );

   int newSingletonRows = 0;

   for(int i=0; i<(int)colAdaptParent.size(); i++)
   {
      int colIdx = colAdaptParent[i].colIdx;
      double val = colAdaptParent[i].val;

      for( int j = currBlmatTrans->rowptr[colIdx].start; j<currBlmatTrans->rowptr[colIdx].end; j++ )
      {
         int rowIdx = currBlmatTrans->jcolM[j];
         double m = removeEntryInDynamicStorage(*currBlmat, rowIdx, colIdx);

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
         dynamic_cast<SimpleVector*>(nColElems->vec)->elements()[colIdx] --;
         currNnzRow->elements()[rowIdx] --;

         assert( currNnzRow->elements()[rowIdx] >= 0.0 );
         assert( dynamic_cast<SimpleVector*>(nColElems->vec)->elements()[colIdx] >= 0.0 );

         if(currNnzRow->elements()[rowIdx] == 1.0)
            newSingletonRows++;
      }
      clearRow(*currBlmatTrans, colIdx);
   }
   return newSingletonRows;
}

bool StochPresolver::combineColAdaptParent()
{
   int myRank, world_size;
   bool iAmDistrib = false;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   if( world_size > 1) iAmDistrib = true;

   if( iAmDistrib )
   {
      // allgather the length of each colAdaptParent
      int mylen = colAdaptParent.size();
      int* recvcounts = new int[world_size];

      MPI_Allgather(&mylen, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

      // allgatherv the actual colAdaptParents
      // First, extract the colIdx and val into int* and double* arrays:
      int* colIndicesLocal = new int[mylen];
      double* valuesLocal = new double[mylen];
      for(int i=0; i<mylen; i++)
      {
         colIndicesLocal[i] = colAdaptParent[i].colIdx;
         valuesLocal[i] = colAdaptParent[i].val;
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
      double* valuesGlobal = new double[lenghtGlobal];
      // Then, do the actual MPI communication:
      MPI_Allgatherv(colIndicesLocal, mylen, MPI_INT, colIndicesGlobal, recvcounts, displs , MPI_INT, MPI_COMM_WORLD);
      MPI_Allgatherv(valuesLocal, mylen, MPI_DOUBLE, valuesGlobal, recvcounts, displs , MPI_DOUBLE, MPI_COMM_WORLD);

      // Reconstruct a colAdaptParent:
      colAdaptParent.clear();
      for(int i=0; i<lenghtGlobal; i++)
      {
         COLUMNTOADAPT colWithVal = {colIndicesGlobal[i], valuesGlobal[i]};
         colAdaptParent.push_back(colWithVal);
      }

      delete[] recvcounts;
      delete[] colIndicesLocal;
      delete[] valuesLocal;
      delete[] colIndicesGlobal;
      delete[] valuesGlobal;
   }

   // Sort colIndicesGlobal (and valuesGlobal accordingly), remove duplicates and find infeasibilities
   std::sort(colAdaptParent.begin(), colAdaptParent.end(), col_is_smaller());

   if(colAdaptParent.size() > 0)
   {
      int colIdxCurrent = colAdaptParent[0].colIdx;
      double valCurrent = colAdaptParent[0].val;
      for(int i=1; i<(int)colAdaptParent.size(); i++)
      {
         if( colAdaptParent[i].colIdx == colIdxCurrent )
         {
            if( colAdaptParent[i].val != valCurrent )
            {
               cout<<"Detected infeasibility (in variable) "<<colIdxCurrent<<endl;
               return false;
            }

            else
               colAdaptParent.erase(colAdaptParent.begin()+i);   //todo: implement more efficiently
         }
         else{
            colIdxCurrent = colAdaptParent[i].colIdx;
            valCurrent = colAdaptParent[i].val;
         }
      }
   }
   assert( (int)colAdaptParent.size() <= nColElems->vec->n );

   return true;
}



int StochPresolver::doSingletonRowsC()
{
   int nelims = 0;
   return nelims;
}


void StochPresolver::updateObjOffset()
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &objOffset, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   if( myRank == 0 )
      cout<<"Objective offset is "<<objOffset<<endl;
}

