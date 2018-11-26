/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactoryAugSparseRowMajSC.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAugSparseRowMaj.h"


sLinsysRoot* sFactoryAugSpTripletSC::newLinsysRoot()
{
  sLinsysRoot* l = new sLinsysRootAugSpTriplet(this, data);
  l->initialize(this,data);
  return l;
}

sLinsysRoot* 
sFactoryAugSpTripletSC::newLinsysRoot(sData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag)
{
  sLinsysRoot* l = new sLinsysRootAugSpTriplet(this, data);new sLinsysRootAugSpTriplet(this, prob, dd, dq, nomegaInv, rhs, additiveDiag);
  l->initialize(this,prob);
  return l;
}
