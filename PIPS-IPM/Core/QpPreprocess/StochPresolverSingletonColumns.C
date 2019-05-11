/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonColumns.h"


StochPresolverSingletonColumns::StochPresolverSingletonColumns(PresolveData& presData, const sData& origProb, StochPostsolver* postsolver)
   : StochPresolverBase(presData, origProb, postsolver)
{
 // todo
}

StochPresolverSingletonColumns::~StochPresolverSingletonColumns()
{
 // todo
}


void StochPresolverSingletonColumns::applyPresolving()
{
   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);

   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );

#ifndef NDEBUG
   if( myRank == 0 )
   {
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
      std::cout << "--- Before singleton columns presolving:" << std::endl;
   }
   countRowsCols();
#endif


#ifndef NDEBUG
   if( myRank == 0 )
      std::cout << "--- After singleton columns presolving:" << std::endl;
   countRowsCols();
   if( myRank == 0 )
      std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
#endif

   assert(presData.reductionsEmpty());
   assert(presData.presProb->isRootNodeInSync());
   assert(verifyNnzcounters());
   assert(indivObjOffset == 0.0);
   assert(newBoundsParent.size() == 0);
}

void StochPresolverSingletonColumns::countSingletonColumns()
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   int nSC=0;
   int nSColEq=0;
   int nSColIneq=0;
   int nZeroCostSc=0;
   int nLowerBound=0;
   int nUpperBound=0;
   int nBothBounds=0;
   int nNoBounds=0;
   int nSColEqLinkRow=0;
   int nSColIneqLinkrow=0;
   initSingletonColumns(nSColEq, nSColIneq, nZeroCostSc, nLowerBound, nUpperBound, nBothBounds, nNoBounds, nSColEqLinkRow, nSColIneqLinkrow, nSC);

   //synchronizeSum(nSColEq, nSColIneq);
   //synchronize(nZeroCostSc);
   synchronizeSumSeveral(nSColEq, nSColIneq, nZeroCostSc, nLowerBound, nUpperBound, nBothBounds, nNoBounds, nSColEqLinkRow, nSColIneqLinkrow ,nSC);
   assert( nSC == nSColEq + nSColIneq );
   if(myRank == 0)
   {
      cout<<"Number of singleton columns in A: "<<nSColEq<<" (with zero objective: "<<nZeroCostSc<<") and in C: "<<nSColIneq<<
            ". In total: "<<nSColEq+nSColIneq<<endl;
      cout<<"In A, SC have: both bounds, only lower, only upper, no bounds: "<<nBothBounds<<", "<<nLowerBound<<", "<<nUpperBound<<", "<<nNoBounds<<endl;
      cout<<"SCs with entry in Linking Rows in A and C: "<<nSColEqLinkRow<<", "<<nSColIneqLinkrow<<endl;
   }
}

void StochPresolverSingletonColumns::initSingletonColumns(int& nSColEq, int& nSColIneq, int& nZeroCostSc,
      int& nLowerBound, int& nUpperBound, int& nBothBounds, int& nNoBounds, int& nSColEqLinkRow, int& nSColIneqLinkrow, int& nSC)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   // ignore singleton cols in the linking columns for now.
   int nSCLinkVars=0;
   currNnzColParent = dynamic_cast<SimpleVector*>(presData.nColElems->vec);
   for(int i=0; i<currNnzColParent->n; i++)\
   {
      if( currNnzColParent->elements()[i]==1.0 )
         nSCLinkVars++;
   }
   if(myRank == 0)
      cout<<"Number of SCs in the linking variables blocks: "<<nSCLinkVars<<endl;

   // singleton cols in the children: count which of them are in A and C respectively.
   assert((int)presData.nColElems->children.size() == nChildren);
   for( size_t it = 0; it < presData.nColElems->children.size(); it++)
   {
     if( ! presData.nColElems->children[it]->isKindOf(kStochDummy) )
     {
        SimpleVector* nColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
        initSingletonColsBlock((int)it, nColChild, nSColEq, nSColIneq, nZeroCostSc, nLowerBound, nUpperBound, nBothBounds, nNoBounds,
              nSColEqLinkRow, nSColIneqLinkrow, nSC);
     }
   }
}

void StochPresolverSingletonColumns::initSingletonColsBlock(int it, SimpleVector const * nnzColSimple,
      int& nSColEq, int& nSColIneq, int& nZeroCostSC, int& nLowerBound, int& nUpperBound, int& nBothBounds, int& nNoBounds,
      int& nSColEqLinkRow, int& nSColIneqLinkrow, int& nSC)
{
   double* nnzCol = nnzColSimple->elements();
   for( int i = 0; i < nnzColSimple->n; i++ )
      if( nnzCol[i] == 1.0 )
         nSC++;

   if(true)
//   if( setCPBmatsChild(presProb->A, (int)it, EQUALITY_SYSTEM) )
   {
//      setCPBlmatsChild(presProb->A, (int)it);
      for( int i = 0; i < nnzColSimple->n; i++ )
         if( nnzCol[i] == 1.0 )
            if( currBmatTrans->rowptr[i].start +1 == currBmatTrans->rowptr[i].end ||
                  currBlmatTrans->rowptr[i].start +1 == currBlmatTrans->rowptr[i].end )
            {
               nSColEq++;
               if(dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec)->elements()[i] == 0.0 )
                  nZeroCostSC++;
//               setCPColumnChild(it);
               if( currIxlowChild->elements()[i]!=0.0 && currIxuppChild->elements()[i]!=0.0 )
                  nBothBounds++;
               else if( currIxlowChild->elements()[i]!=0.0 )
               {
                  nLowerBound++;
                  if( currxlowChild->elements()[i] < -std::numeric_limits<double>::max()*0.1 )
                     cout<<"Lower x bound close to -inf! "<<endl;
               }
               else if( currIxuppChild->elements()[i]!=0.0 )
                  nUpperBound++;
               else if( currIxlowChild->elements()[i]==0.0 && currIxuppChild->elements()[i]==0.0 )
                  nNoBounds++;

               if( currBlmatTrans->rowptr[i].start +1 == currBlmatTrans->rowptr[i].end )
                  nSColEqLinkRow++;
            }
   }
//   if( setCPBmatsChild(presProb->C, (int)it, INEQUALITY_SYSTEM) )
   {
//      setCPBlmatsChild(presProb->C, (int)it);
      for( int i = 0; i < nnzColSimple->n; i++ )
         if( nnzCol[i] == 1.0 )
            if( currBmatTrans->rowptr[i].start +1 == currBmatTrans->rowptr[i].end ||
                  currBlmatTrans->rowptr[i].start +1 == currBlmatTrans->rowptr[i].end )
            {
               nSColIneq++;
               if( currBlmatTrans->rowptr[i].start +1 == currBlmatTrans->rowptr[i].end )
                  nSColIneqLinkrow++;
            }
   }
}

void StochPresolverSingletonColumns::synchronizeSumSeveral(int& val0, int& val1, int& val2, int& val3, int& val4, int& val5, int& val6,
      int&val7, int& val8, int& val9 )
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );
   if( iAmDistrib )
   {
      int* tmpArray = new int[10];
      tmpArray[0] = val0;
      tmpArray[1] = val1;
      tmpArray[2] = val2;
      tmpArray[3] = val3;
      tmpArray[4] = val4;
      tmpArray[5] = val5;
      tmpArray[6] = val6;
      tmpArray[7] = val7;
      tmpArray[8] = val8;
      tmpArray[9] = val9;
      MPI_Allreduce(MPI_IN_PLACE, tmpArray, 10, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      val0 = tmpArray[0];
      val1 = tmpArray[1];
      val2 = tmpArray[2];
      val3 = tmpArray[3];
      val4 = tmpArray[4];
      val5 = tmpArray[5];
      val6 = tmpArray[6];
      val7 = tmpArray[7];
      val8 = tmpArray[8];
      val9 = tmpArray[9];
      delete[] tmpArray;
   }
}

