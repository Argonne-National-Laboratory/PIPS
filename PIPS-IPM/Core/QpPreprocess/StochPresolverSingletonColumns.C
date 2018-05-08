/*
 * StochPresolverSingletonColumns.C
 *
 *  Created on: 08.05.2018
 *      Author: bzfuslus
 */

#include "StochPresolverSingletonColumns.h"


StochPresolverSingletonColumns::StochPresolverSingletonColumns(PresolveData& presData)
: StochPresolverBase(presData)
{
 // todo
}

StochPresolverSingletonColumns::~StochPresolverSingletonColumns()
{
 // todo
}


bool StochPresolverSingletonColumns::applyPresolving(int& nelims)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   int nSColEq=0;
   int nSColIneq=0;
   int nZeroCostSc=0;
   initSingletonColumns(nSColEq, nSColIneq, nZeroCostSc);

   synchronizeSum(nSColEq, nSColIneq);
   synchronize(nZeroCostSc);
   if(myRank == 0)
      cout<<"Number of singleton columns in A: "<<nSColEq<<" (with zero objective: "<<nZeroCostSc<<") and in C: "<<nSColIneq<<endl;

   return true;
}

void StochPresolverSingletonColumns::initSingletonColumns(int& nSColEq, int& nSColIneq, int& nZeroCostSc)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   // ignore singleton cols in the linking columns for now.
   // singleton cols in the children: count which of them are in A and C respectively.
   assert((int)presData.nColElems->children.size() == nChildren);
   for( size_t it = 0; it < presData.nColElems->children.size(); it++)
   {
     if( ! presData.nColElems->children[it]->isKindOf(kStochDummy) )
     {
        SimpleVector* nColChild = dynamic_cast<SimpleVector*>(presData.nColElems->children[it]->vec);
        initSingletonColsBlock((int)it, nColChild, nSColEq, nSColIneq, nZeroCostSc);
     }
   }
}

void StochPresolverSingletonColumns::initSingletonColsBlock(int it, SimpleVector const * nnzColSimple, int& nSColEq, int& nSColIneq, int& nZeroCostSC)
{
   //todo: check in Blmat as well.
   double* nnzCol = nnzColSimple->elements();

   if( setCPBmatsChild(presProb->A, (int)it, EQUALITY_SYSTEM) )
   {
      setCPBlmatsChild(presProb->A, (int)it);
      for( int i = 0; i < nnzColSimple->n; i++ )
         if( nnzCol[i] == 1.0 )
            if( currBmatTrans->rowptr[i].start +1 == currBmatTrans->rowptr[i].end ||
                  currBlmatTrans->rowptr[i].start +1 == currBlmatTrans->rowptr[i].end )
            {
               nSColEq++;
               if(dynamic_cast<SimpleVector*>(dynamic_cast<StochVector&>(*(presProb->g)).children[it]->vec)->elements()[i] == 0.0 )
                  nZeroCostSC++;
            }
   }
   if( setCPBmatsChild(presProb->C, (int)it, INEQUALITY_SYSTEM) )
   {
      setCPBlmatsChild(presProb->C, (int)it);
      for( int i = 0; i < nnzColSimple->n; i++ )
         if( nnzCol[i] == 1.0 )
            if( currBmatTrans->rowptr[i].start +1 == currBmatTrans->rowptr[i].end ||
                  currBlmatTrans->rowptr[i].start +1 == currBlmatTrans->rowptr[i].end )
               nSColIneq++;
   }
}

