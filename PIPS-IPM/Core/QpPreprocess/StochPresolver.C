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

   StochVectorHandle gcloneA(dynamic_cast<StochVector*>(sorigprob->g->clone()));
   nColElemsA = gcloneA;

   StochVectorHandle gcloneC(dynamic_cast<StochVector*>(sorigprob->g->clone()));
   nColElemsC = gcloneC;

   StochVectorHandle bAclone(dynamic_cast<StochVector*>(sorigprob->bA->clone()));
   nRowElemsA = bAclone;

   StochVectorHandle icuppclone(dynamic_cast<StochVector*>(sorigprob->icupp->clone()));
   nRowElemsA = icuppclone;

   presProb = NULL;
}


StochPresolver::~StochPresolver()
{
}


Data* StochPresolver::presolve()
{

   std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);

   presProb = sorigprob->cloneFull();


   // initialized all transposed sub matrices
   presProb->A->initTransposed();
   presProb->C->initTransposed();

   int cleanup_elims = cleanUp();

   return presProb;
}


int StochPresolver::cleanUp()
{
   int nelims = 0;

   GenMatrix * C_p = SpAsPointer(presProb->C);
   StochGenMatrix* C = dynamic_cast<StochGenMatrix*>(C_p);
   OoqpVector * bl_p = SpAsPointer(presProb->bl);
   StochVector* clow = dynamic_cast<StochVector*>(bl_p);
   OoqpVector * bu_p = SpAsPointer(presProb->bu);
   StochVector* cupp = dynamic_cast<StochVector*>(bu_p);
   OoqpVector * blx_p = SpAsPointer(presProb->blx);
   StochVector* xlow = dynamic_cast<StochVector*>(blx_p);
   OoqpVector * bux_p = SpAsPointer(presProb->bux);
   StochVector* xupp = dynamic_cast<StochVector*>(bux_p);

   // handle parent
   assert(C->children.size() == cupp->children.size());

   C->writeToStreamDense(cout);
   for( size_t it = 0; it < C->children.size(); it++ )

   {
      // call
      nelims += cleanUpRowC(C->children[it], clow->children[it], cupp->children[it],
            xlow->children[it], xupp->children[it]);

   }
   C->writeToStreamDense(cout);

   return nelims;
}


// new method for Amat, Bmat
int StochPresolver::cleanUpRowC(StochGenMatrix* matrixC, StochVector* clow, StochVector* cupp,
      StochVector* xlow, StochVector* xupp)
{
   int nelims = 0;
   int m,n;
   matrixC->Bmat->getSize(m, n);

   for( int i = 0; i < m; i++ ) // Row i
   {
      int j = 0; // Column j
      for( int k = matrixC->Bmat->krowM()[i]; k < matrixC->Bmat->krowM()[i + 1];
            k++ )
      {
         if( matrixC->Bmat->M()[k] < 1e-3 )
         {//((SimpleVector*)icupp->vec)->elements();

            double xu_j = ((SimpleVector*)xupp->vec)->elements()[j];
            double xl_j = ((SimpleVector*)xlow->vec)->elements()[j];

            int supp = matrixC->Bmat->krowM()[i+1] - matrixC->Bmat->krowM()[i];

            if( fabs(matrixC->Bmat->M()[k]) * (xu_j - xl_j) * supp < 1e-2 * feastol)
            {
               // set entry a_ij to 0
               // todo: adapt matrixC and nnz and rhs
            matrixC->Bmat->M()[k] = 0;
            nelims++;
            }

         }
         j++;
      }

   }

   return nelims;
}

