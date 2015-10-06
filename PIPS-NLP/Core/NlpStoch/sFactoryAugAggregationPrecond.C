/* PIPS-NLP                                                           	*
 * Author:  Nai-Yuan Chiang                                       	*
 * (C) 2015 Argonne National Laboratory. 				*/

#include "sFactoryAugAggregationPrecond.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAggregation.h"
#include "Ma57Solver.h"
#include "sLinsysLeaf.h"



sLinsysRoot* 
sFactoryAugAggregationPrecond::newLinsysRoot()
{
  return new sLinsysRootAggregation(this, data);
}

sLinsysRoot* 
sFactoryAugAggregationPrecond::newLinsysRoot(sData* prob,
					   OoqpVector* dd,
					   OoqpVector* dq,
					   OoqpVector* nomegaInv, 
					   OoqpVector* rhs, OoqpVector* additiveDiag)
{
  return new sLinsysRootAggregation(this, prob, dd, dq, nomegaInv, rhs, additiveDiag);
}

sLinsysLeaf* 
sFactoryAugAggregationPrecond::newLinsysLeaf(sData* prob,
			OoqpVector* dd, OoqpVector* dq,
			OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag)
{
  Ma57Solver* sMA57=NULL; 

  return new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs, additiveDiag, sMA57);
}


