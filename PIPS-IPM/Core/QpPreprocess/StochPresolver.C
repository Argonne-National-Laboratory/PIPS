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

#include "../Abstract/DoubleMatrix.h"
#include "../SparseLinearAlgebra/SparseGenMatrix.h"
#include "../StochLinearAlgebra/StochVectorHandle.h"
#include "../Vector/OoqpVector.h"

StochPresolver::StochPresolver(const Data* prob)
 : QpPresolver(prob)
{
   const sData* sorigprob = dynamic_cast<const sData*>(prob);

   StochVectorHandle gclone(dynamic_cast<StochVector*>(sorigprob->g->clone()));
   nColElems = gclone;

   StochVectorHandle bAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   nRowElemsA = bAclone;

   StochVectorHandle icuppclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   nRowElemsC = icuppclone;

   presProb = NULL;
}


StochPresolver::~StochPresolver()
{
}

void
StochPresolver::initNnzCounter()
{
   StochGenMatrix& A = dynamic_cast<StochGenMatrix&>(*(presProb->A));
   StochGenMatrix& C = dynamic_cast<StochGenMatrix&>(*(presProb->C));

   A.addNnzPerRow(*nRowElemsA);
   C.addNnzPerRow(*nRowElemsC);

   C.addNnzPerCol(*nColElems);
   A.addNnzPerCol(*nColElems);


   int rank = 0;

   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   if( rank == 0)
   {
      std::cout << "write Cols" << std::endl;
         nColElems->writeToStreamAll(std::cout);

      std::cout << "writeA" << std::endl;

      nRowElemsA->writeToStreamAll(std::cout);

      std::cout << "writeC " << std::endl;


      nRowElemsC->writeToStreamAll(std::cout);
   }



}

Data* StochPresolver::presolve()
{

   std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

   presProb = sorigprob->cloneFull();

   // initialized all transposed sub matrices
   presProb->A->initTransposed();
   presProb->C->initTransposed();

   initNnzCounter();




  // int cleanup_elims = cleanUp();

   return presProb;
}


int StochPresolver::cleanUp()
{
   int nelims = 0;

   const StochGenMatrix& C = dynamic_cast<const StochGenMatrix&>(*(presProb->C));
   const StochVector& clow = dynamic_cast<const StochVector&>(*(presProb->bl));
   const StochVector& cupp = dynamic_cast<const StochVector&>(*(presProb->bu));
   const StochVector& xlow = dynamic_cast<const StochVector&>(*(presProb->blx));
   const StochVector& xupp = dynamic_cast<const StochVector&>(*(presProb->bux));

   // handle parent
   assert(C.children.size() == cupp.children.size());

   C.writeToStreamDense(cout);
   for( size_t it = 0; it < C.children.size(); it++ )
   {
      // call
      nelims += cleanUpRowsC(C.children[it], clow.children[it], cupp.children[it],
            xlow.children[it], xupp.children[it], nRowElemsC->children[it]);

   }
   C.writeToStreamDense(cout);

   assert(0);


   // todo StochGenDummyMatrix
   return nelims;
}


// new method for Amat, Bmat
int StochPresolver::cleanUpRowsC(StochGenMatrix* matrixC, StochVector* clow, StochVector* cupp,
      StochVector* xlow, StochVector* xupp, StochVector* nnzRowC)
{
   int nelims = 0;
   int m,n;
   matrixC->Bmat->getSize(m, n);

   double* const val = matrixC->Bmat->M();
   const int* const rowStart = matrixC->Bmat->krowM();
   const int* const colIdx = matrixC->Bmat->jcolM();

   double* const xuppElems = (dynamic_cast<SimpleVector*>(xupp->vec))->elements();
   double* const xlowElems = (dynamic_cast<SimpleVector*>(xlow->vec))->elements();
   double* const nnzPerRow = (dynamic_cast<SimpleVector*>(nnzRowC->vec))->elements();

   std::cout << "lower: " << std::endl;
   (dynamic_cast<SimpleVector*>(xlow->vec))->writeToStream(std::cout);

   std::cout << "upper: " << std::endl;
   (dynamic_cast<SimpleVector*>(xupp->vec))->writeToStream(std::cout);

   for( int r = 0; r < m; r++ ) // Row i
   {
      const int rowEnd = rowStart[r + 1];
      std::cout << "length: " << rowEnd - rowStart[r] << std::endl;

      for( int k = rowStart[r]; k < rowEnd; k++ )
      {
         if( val[k] < 1e-3 )
         {
            const int col = colIdx[k]; // Column j

            const double ub = xuppElems[col];
            const double lb = xlowElems[col];

            // todo
            const double supp = nnzPerRow[r];

            if( fabs(val[k]) * (ub - lb) * supp < 1e-2 * feastol)
            {
               // set entry a_ij to 0
               // todo: adapt matrixC and nnz and rhs
               val[k] = 0.0;
               nelims++;
               nnzPerRow[r] -= 1;
               std::cout << "eliminate " << r << "of " << m  << " " << col <<  "of " << n << std::endl;
            }
         }
      }
   }

   return nelims;
}

