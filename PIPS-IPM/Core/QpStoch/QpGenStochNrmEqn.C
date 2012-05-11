#include "QpGenStochNrmEqn.h"

#include "QpGenStochData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "QpGenStochLinsysRootNrmEqn.h"

QpGenStochNrmEqn::QpGenStochNrmEqn( StochInputTree* inputTree)
  : QpGenStoch(inputTree)
{ };

QpGenStochNrmEqn::QpGenStochNrmEqn( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGenStoch(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

QpGenStochNrmEqn::QpGenStochNrmEqn()
{ };

QpGenStochNrmEqn::~QpGenStochNrmEqn()
{ };


QpGenStochLinsysRoot* QpGenStochNrmEqn::newLinsysRoot()
{
  return new QpGenStochLinsysRootNrmEqn(this, data);
}

QpGenStochLinsysRoot* 
QpGenStochNrmEqn::newLinsysRoot(QpGenStochData* prob,
			  OoqpVector* dd,OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  new QpGenStochLinsysRootNrmEqn(this, prob, dd, dq, nomegaInv, rhs);
}
