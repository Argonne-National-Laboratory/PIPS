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

}


StochPresolver::~StochPresolver()
{
}


Data* StochPresolver::presolve()
{

   std::cout << "start stoch presolving" << std::endl;

   const sData* sorigprob = dynamic_cast<const sData*>(origprob);



   sData* newprob = sorigprob->cloneFull();

   std::cout << "CLONED" << std::endl;

   return newprob;
}
