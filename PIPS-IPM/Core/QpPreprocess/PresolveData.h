/*
 * PresolveData.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_
#define PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_


class PresolveData
{
      PresolveData(const sData sorigprob);
      ~PresolveData();

      void initialize();

      sData* finalize();


      /** objective offset created by presolving*/
      double objOffset;



      // number of non-zero elements of each row
      StochVectorHandle nRowElemsA;
      StochVectorHandle nRowElemsC;

      // number of non-zero elements of each column
      StochVectorHandle nColElems;

      // number of removed elements of each row / column
      StochVectorHandle redRowA;
      StochVectorHandle redRowC;
      StochVectorHandle redCol;

      sData* presProb;

      StochGenMatrix& Apres;
      StochGenMatrix& Cpres;

      // todo getter, setter for element access of nnz counter???

   private:
      // initialize row and column nnz counter
      void initNnzCounter();

};



PresolveData::PresolveData(const sData sorigprob)
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

   Apres = dynamic_cast<StochGenMatrix&>(*presProb->A);
   Cpres = dynamic_cast<StochGenMatrix&>(*presProb->C);
   objOffset = 0.0;






}

void
PresolveData::initialize()
{
   // initialized all dynamic transposed sub matrices todo in extra method??? initiliaze?
   Apres.initTransposed(true);
   Cpres.initTransposed(true);

   initNnzCounter();
}

sData*
PresolveData::finalize()
{
   presProb->cleanUpPresolvedData(*nRowElemsA, *nRowElemsC, *nColElems);

   Apres.deleteTransposed();
   Cpres.deleteTransposed();

   return presProb;
}

void
PresolveData::initNnzCounter()
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

#endif /* PIPS_IPM_CORE_QPPREPROCESS_PRESOLVEDATA_H_ */
