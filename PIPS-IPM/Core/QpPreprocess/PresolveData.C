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

   //Apresa = dynamic_cast<StochGenMatrix&>(*presProb->A);
   //Cpres = dynamic_cast<StochGenMatrix&>(*presProb->C);
   objOffset = 0.0;

   nChildren = nColElems->children.size();

   /* zero initialized because of () */
   blocks = new int[nChildren + 3]();
   blocksIneq = new int[nChildren + 3]();

   // initialize all dynamic transposed sub matrices
   dynamic_cast<StochGenMatrix&>(*presProb->A).initTransposed(true);
   dynamic_cast<StochGenMatrix&>(*presProb->C).initTransposed(true);

   initNnzCounter();
}

PresolveData::~PresolveData()
{
   delete[] blocks;
   delete[] blocksIneq;
}

sData* PresolveData::finalize()
{
   // this removes all columns and rows that are now empty from the problem
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

bool PresolveData::combineColAdaptParent()
{
   int myRank, world_size;
   bool iAmDistrib = false;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   MPI_Comm_size(MPI_COMM_WORLD, &world_size);
   if( world_size > 1) iAmDistrib = true;

   if( iAmDistrib )
   {
      // allgather the length of each colAdaptParent
      const int mylen = getNumberColAdParent();
      int* recvcounts = new int[world_size];

      MPI_Allgather(&mylen, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

      // allgatherv the actual colAdaptParents
      // First, extract the colIdx and val into int* and double* arrays:
      int* colIndicesLocal = new int[mylen];
      double* valuesLocal = new double[mylen];
      for(int i=0; i<mylen; i++)
      {
         colIndicesLocal[i] = getColAdaptParent(i).colIdx;
         valuesLocal[i] = getColAdaptParent(i).val;
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
      clearColAdaptParent();
      for(int i=0; i<lenghtGlobal; i++)
      {
         COLUMNTOADAPT colWithVal = {colIndicesGlobal[i], valuesGlobal[i]};
         addColToAdaptParent(colWithVal);
      }

      delete[] displs;
      delete[] recvcounts;
      delete[] colIndicesLocal;
      delete[] valuesLocal;
      delete[] colIndicesGlobal;
      delete[] valuesGlobal;
   }

   // Sort colIndicesGlobal (and valuesGlobal accordingly), remove duplicates and find infeasibilities
   std::sort(colAdaptParent.begin(), colAdaptParent.end(), col_is_smaller());
   for(int i=1; i<getNumberColAdParent(); i++)
      assert( getColAdaptParent(i-1).colIdx <= getColAdaptParent(i).colIdx);

   if(getNumberColAdParent() > 0)
   {
      int colIdxCurrent = getColAdaptParent(0).colIdx;
      double valCurrent = getColAdaptParent(0).val;
      for(int i=1; i<getNumberColAdParent(); i++)
      {
         if( getColAdaptParent(i).colIdx == colIdxCurrent )
         {
            if( getColAdaptParent(i).val != valCurrent )
            {
               cout<<"Detected infeasibility (in variable) "<<colIdxCurrent<<endl;
               return false;
            }
            else
            {
               colAdaptParent.erase(colAdaptParent.begin()+i);   //todo: implement more efficiently
               i--;
            }
         }
         else{
            colIdxCurrent = getColAdaptParent(i).colIdx;
            valCurrent = getColAdaptParent(i).val;
         }
      }
   }
   assert( getNumberColAdParent() <= nColElems->vec->n );

   return true;
}

bool PresolveData::reductionsEmpty()
{


	return true;
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

int PresolveData::getNChildren() const
{
   return nChildren;
}

double PresolveData::getObjOffset() const
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

int PresolveData::getSingletonRow(int i) const
{
   assert(i<getNumberSR() && i>=0);
   return singletonRows[i];
}
int PresolveData::getNumberSR() const
{
   return (int)singletonRows.size();
}
void PresolveData::addSingletonRow(int i)
{
   singletonRows.push_back(i);
}
void PresolveData::setSingletonRow(int i, int value)
{
   assert(i>=0 && i<getNumberSR());
   singletonRows[i] = value;
}
/* empties the singletonRows list. Verifies first if all elements were set to -1,
 * then clears the singletonRows list and resets the blocks. */
void PresolveData::clearSingletonRows()
{
   for(int i = 0; i<getNumberSR(); i++)
      assert(getSingletonRow(i) == -1);
   singletonRows.clear();
   resetBlocks();
}
int PresolveData::getSingletonRowIneq(int i) const
{
   assert(i<getNumberSRIneq() && i>=0);
   return singletonRowsIneq[i];
}
int PresolveData::getNumberSRIneq() const
{
   return (int)singletonRowsIneq.size();
}
void PresolveData::addSingletonRowIneq(int i)
{
   singletonRowsIneq.push_back(i);
}
void PresolveData::setSingletonRowIneq(int i, int value)
{
   assert(i>=0 && i<getNumberSRIneq());
   singletonRowsIneq[i] = value;
}
void PresolveData::clearSingletonRowsIneq()
{
   for(int i = 0; i<getNumberSRIneq(); i++)
      assert(getSingletonRowIneq(i) == -1);
   singletonRowsIneq.clear();
   resetBlocks();
}

void PresolveData::setBlocks(int i, double value)
{
   assert(i<nChildren+3 && i>=0);
   blocks[i] = value;
}
double PresolveData::getBlocks(int i) const
{
   assert(i<nChildren+3 && i>=0);
   return blocks[i];
}
void PresolveData::setBlocksIneq(int i, double value)
{
   assert(i<nChildren+3 && i>=0);
   blocksIneq[i] = value;
}
double PresolveData::getBlocksIneq(int i) const
{
   assert(i<nChildren+3 && i>=0);
   return blocksIneq[i];
}

COLUMNTOADAPT PresolveData::getColAdaptParent(int i) const
{
   assert( i<getNumberColAdParent() );
   return colAdaptParent[i];
}
int PresolveData::getNumberColAdParent() const
{
   return (int)colAdaptParent.size();
}
void PresolveData::addColToAdaptParent(COLUMNTOADAPT colToAdapt)
{
   colAdaptParent.push_back(colToAdapt);
}
void PresolveData::clearColAdaptParent()
{
   colAdaptParent.clear();
}
