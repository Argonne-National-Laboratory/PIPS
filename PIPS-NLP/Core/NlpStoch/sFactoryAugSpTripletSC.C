/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/   

#include "sFactoryAugSpTripletSC.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAug.h"


sLinsysRoot* sFactoryAugSpTripletSC::newLinsysRoot()
{
  return new sLinsysRootAug(this, data);
}

sLinsysRoot* 
sFactoryAugSpTripletSC::newLinsysRoot(sData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag)
{
  return new sLinsysRootAug(this, prob, dd, dq, nomegaInv, rhs, additiveDiag);
}



