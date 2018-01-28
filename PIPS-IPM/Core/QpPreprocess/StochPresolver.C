/*
 * StochPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */


#include "StochPresolver.h"
#include <cassert>
#include <iostream>

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
   assert(0);

   return 0;
}
