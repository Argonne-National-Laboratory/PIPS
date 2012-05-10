#include "QpGenStochAugRed.h"

#include "QpGenStochData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "QpGenStochLinsysRootAugRed.h"

QpGenStochAugRed::QpGenStochAugRed( StochInputTree* inputTree)
  : QpGenStoch(inputTree)
{ };

QpGenStochAugRed::QpGenStochAugRed( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGenStoch(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

QpGenStochAugRed::QpGenStochAugRed()
{ };

QpGenStochAugRed::~QpGenStochAugRed()
{ };


QpGenStochLinsysRoot* QpGenStochAugRed::newLinsysRoot()
{
  return new QpGenStochLinsysRootAugRed(this, data);
}

QpGenStochLinsysRoot* 
QpGenStochAugRed::newLinsysRoot(QpGenStochData* prob,
			  OoqpVector* dd,OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  new QpGenStochLinsysRootAugRed(this, prob, dd, dq, nomegaInv, rhs);
}
