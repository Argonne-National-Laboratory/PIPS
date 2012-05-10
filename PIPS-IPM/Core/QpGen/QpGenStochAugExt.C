#include "QpGenStochAugExt.h"

#include "QpGenStochData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "QpGenStochLinsysRootAugExt.h"

QpGenStochAugExt::QpGenStochAugExt( StochInputTree* inputTree)
  : QpGenStoch(inputTree)
{ };

QpGenStochAugExt::QpGenStochAugExt( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGenStoch(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

QpGenStochAugExt::QpGenStochAugExt()
{ };

QpGenStochAugExt::~QpGenStochAugExt()
{ };


QpGenStochLinsysRoot* QpGenStochAugExt::newLinsysRoot()
{
  return new QpGenStochLinsysRootAugExt(this, data);
}

QpGenStochLinsysRoot* 
QpGenStochAugExt::newLinsysRoot(QpGenStochData* prob,
			  OoqpVector* dd,OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  new QpGenStochLinsysRootAugExt(this, prob, dd, dq, nomegaInv, rhs);
}
