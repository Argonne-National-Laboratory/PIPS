/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactoryAug.h"

#include "QpGenStochData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAug.h"

sFactoryAug::sFactoryAug( StochInputTree* inputTree)
  : sFactory(inputTree)
{ };

sFactoryAug::sFactoryAug( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : sFactory(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

sFactoryAug::sFactoryAug()
{ };

sFactoryAug::~sFactoryAug()
{ };


sLinsysRoot* sFactoryAug::newLinsysRoot()
{
  return new sLinsysRootAug(this, data);
}

sLinsysRoot* 
sFactoryAug::newLinsysRoot(QpGenStochData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs)
{
  return new sLinsysRootAug(this, prob, dd, dq, nomegaInv, rhs);
}
