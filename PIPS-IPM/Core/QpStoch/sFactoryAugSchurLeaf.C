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

#include "pipsport.h"

sLinsysLeaf* sFactoryAugSchurLeaf::newLinsysLeaf(sData* prob,
						 OoqpVector* dd,OoqpVector* dq,
						 OoqpVector* nomegaInv, OoqpVector* rhs)
{
  //cout << "sFactoryAugSchurLeaf::newLinsysLeaf  returns  a sLinsysLeafSchurSlv" << endl;
  PardisoSchurSolver* linSolver=nullptr;
  return new sLinsysLeafSchurSlv(this, prob, dd, dq, nomegaInv, rhs, linSolver);
}

sLinsysLeaf* 
sFactoryAugSchur32Leaf::newLinsysLeaf(sData* prob,
				      OoqpVector* dd,OoqpVector* dq,
				      OoqpVector* nomegaInv, OoqpVector* rhs)
{
  //cout << "sFactoryAugSchurLeaf::newLinsysLeaf  returns  a sLinsysLeafSchurSlv" << endl;
  PardisoSchur32Solver* linSolver=nullptr;
  return new sLinsysLeafSchurSlv(this, prob, dd, dq, nomegaInv, rhs, linSolver);
}

