/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cstring>
#include <iostream>
#include <cassert>
using namespace std;
#include "Variables.h"
#include "Solver.h"

// default implementation for Variables::print() prints abusive
// message. Since we don't have any knowledge of how the variables are
// stored at this top level, we can't do much except print their
// dimensions.

void Variables::print()
{
  cout << " Complementary Variables = " << nComplementaryVariables << endl;
  cout << "(Cannot tell you more at this level)" << endl;
}

void Variables::saxpy_pd( Variables *b, double alpha_primal, double alpha_dual)
{
  assert(0 && "not implemented here");
}

double Variables::mustep_pd( Variables *step, double alpha_primal, double alpha_dual )
{
  assert(0 && "not implemented here");
  return 0;
}

void Variables::stepbound_pd( Variables *b, double & alpha_primal, double & alpha_dual )
{
  assert(0 && "not implemented here");
}

void Variables::findBlocking_pd( Variables * step,
				double & primalValue,
				double & primalStep,
				double & dualValue,
				double & dualStep,
				double & primalValue_d, double & primalStep_d, double & dualValue_d, double & dualStep_d,
				double& alphaPrimal, double& alphaDual,
				int& primalBlocking, int& dualBlocking )
{
  assert(0 && "not implemented here");
}
