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

   objOffset = 0.0;

   nChildren = nColElems->children.size();

   // initialize all dynamic transposed sub matrices
   dynamic_cast<StochGenMatrix&>(*presProb->A).initTransposed(true);
   dynamic_cast<StochGenMatrix&>(*presProb->C).initTransposed(true);

   initNnzCounter();
}

PresolveData::~PresolveData()
{
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

   if( world_size > 1)
      iAmDistrib = true;

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
      for(int i = 1; i < world_size; i++)
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
         COLUMNFORDELETION colWithVal = {colIndicesGlobal[i], valuesGlobal[i]};
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
   std::sort(linkingVariablesMarkedForDeletion.begin(), linkingVariablesMarkedForDeletion.end(), col_is_smaller());
   for(int i = 1; i < getNumberColAdParent(); i++)
      assert( getColAdaptParent(i-1).colIdx <= getColAdaptParent(i).colIdx);

   if(getNumberColAdParent() > 0)
   {
      int colIdxCurrent = getColAdaptParent(0).colIdx;
      double valCurrent = getColAdaptParent(0).val;
      for(int i = 1; i < getNumberColAdParent(); i++)
      {
         if( getColAdaptParent(i).colIdx == colIdxCurrent )
         {
            if( getColAdaptParent(i).val != valCurrent )
            {
               std::cout << "Detected infeasibility (in variable) " << colIdxCurrent << std::endl;
               return false;
            }
            else
            {
               linkingVariablesMarkedForDeletion.erase(linkingVariablesMarkedForDeletion.begin()+i);   //todo: implement more efficiently
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

// todo
bool PresolveData::reductionsEmpty()
{
   bool empty = true;

   empty = (empty && redRowA->isZero());
   empty = (empty && redRowC->isZero());
   empty = (empty && redCol->isZero());
	return empty;
}

void PresolveData::resetRedCounters()
{
   redRowA->setToZero();
   redRowC->setToZero();
   redCol->setToZero();
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

COLUMNFORDELETION PresolveData::getColAdaptParent(int i) const
{
   assert( i<getNumberColAdParent() );
   return linkingVariablesMarkedForDeletion[i];
}
int PresolveData::getNumberColAdParent() const
{
   return (int)linkingVariablesMarkedForDeletion.size();
}
void PresolveData::addColToAdaptParent(COLUMNFORDELETION colToAdapt)
{
   linkingVariablesMarkedForDeletion.push_back(colToAdapt);
}
void PresolveData::clearColAdaptParent()
{
   linkingVariablesMarkedForDeletion.clear();
}
