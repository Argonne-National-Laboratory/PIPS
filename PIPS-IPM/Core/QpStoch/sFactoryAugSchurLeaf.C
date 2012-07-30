/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactoryAugSchurLeaf.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"
#include "PardisoSolver.h"
#include "PardisoSchurSolver.h"

#include "sLinsysLeafSchurSlv.h"

sLinsysLeaf* sFactoryAugSchurLeaf::newLinsysLeaf(sData* prob,
						 OoqpVector* dd,OoqpVector* dq,
						 OoqpVector* nomegaInv, OoqpVector* rhs)
{
  //cout << "sFactoryAugSchurLeaf::newLinsysLeaf  returns  a sLinsysLeafSchurSlv" << endl;
  PardisoSchurSolver* linSolver=NULL;
  return new sLinsysLeafSchurSlv(this, prob, dd, dq, nomegaInv, rhs, linSolver);
}
