/*
 * PresolveData.C
 *
 *  Created on: 09.04.2018
 *      Author: bzfuslus
 */

#include "PresolveData.h"


PresolveData::PresolveData(const sData* sorigprob)
{
   StochVectorHandle gclone(dynamic_cast<StochVector*>(sorigprob->g->clone()));
   nColElems = gclone;
   StochVectorHandle colClone(dynamic_cast<StochVector*>(sorigprob->g->clone()));
   redCol = colClone;

   StochVectorHandle bAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   nRowElemsA = bAclone;
   StochVectorHandle rowAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   redRowA = rowAclone;

   StochVectorHandle icuppclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   nRowElemsC = icuppclone;
   StochVectorHandle rowCclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   redRowC = rowCclone;

   presProb = sorigprob->cloneFull(true);

   //Apres = dynamic_cast<StochGenMatrix&>(*presProb->A);
   //Cpres = dynamic_cast<StochGenMatrix&>(*presProb->C);
   objOffset = 0.0;

   nChildren = nColElems->children.size();
   blocks = new int[nChildren + 3];
   blocksIneq = new int[nChildren + 3];
   resetBlocks();
}

PresolveData::~PresolveData()
{
   delete[] blocks;
   delete[] blocksIneq;
}

void PresolveData::initialize()
{
   // initialize all dynamic transposed sub matrices
   dynamic_cast<StochGenMatrix&>(*presProb->A).initTransposed(true);
   dynamic_cast<StochGenMatrix&>(*presProb->C).initTransposed(true);

   initNnzCounter();
}

sData* PresolveData::finalize()
{
   presProb->cleanUpPresolvedData(*nRowElemsA, *nRowElemsC, *nColElems);

   dynamic_cast<StochGenMatrix&>(*presProb->A).deleteTransposed();
   dynamic_cast<StochGenMatrix&>(*presProb->C).deleteTransposed();

   return presProb;
}

void PresolveData::initNnzCounter()
{
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   StochVectorHandle colClone(dynamic_cast<StochVector*>(nColElems->clone()));

   A.getNnzPerRow(*nRowElemsA);
   C.getNnzPerRow(*nRowElemsC);

   A.getNnzPerCol(*nColElems);
   C.getNnzPerCol(*colClone);

   nColElems->axpy(1.0, *colClone);


#if 0
   int rank = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if( rank == 0 ) std::cout << "write Cols with all cols" << std::endl;
   nColElems->writeToStreamAll(std::cout);

   if( rank == 0 ) std::cout << "write Rows A" << std::endl;
   nRowElemsA->writeToStreamAll(std::cout);

   if( rank == 0 ) std::cout << "write Rows C " << std::endl;
   nRowElemsC->writeToStreamAll(std::cout);
#endif
}

void PresolveData::resetRedCounters()
{
   redRowA->setToZero();
   redRowC->setToZero();
   redCol->setToZero();
}

void PresolveData::resetBlocks()
{
   for( int i = 0; i < nChildren+3; i++)
   {
      blocks[i] = 0;
      blocksIneq[i] = 0;
   }
}

int PresolveData::getNChildren()
{
   return nChildren;
}

double PresolveData::getObjOffset()
{
   return objOffset;
}

double PresolveData::addObjOffset(double addOffset)
{
   objOffset += addOffset;
   return objOffset;
}

void PresolveData::setObjOffset(double offset)
{
   objOffset = offset;
}


