/*
 * StochPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */


#include "StochPresolver.h"
#include <cassert>
#include <iostream>
#include "sData.h"

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


int
StochPresolver::cleanUp()
{
   int nelims = 0;
   StochGenMatrix* C =  dynamic_cast<StochGenMatrix*>(presProb);
   StochVector* cupp;

   // handle parent
   assert(C->children.size() == cupp->children());


   for( size_t it = 0; it < C->children.size(); it++ )

   {
      // call
      C->children[it] cupp->children[it];

   }



   return nelims;
}


// new method for Amat, Bmat


