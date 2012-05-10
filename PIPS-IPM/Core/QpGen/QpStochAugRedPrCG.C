#include "QpStochAugRedPrCG.h"

#include "QpGenStochData.h"

#include "StochTreePrecond.h"
#include "StochInputTree.h"

#include "QpStochLinsysRootAugRedPrCG.h"

QpStochAugRedPrCG::QpStochAugRedPrCG( StochInputTree* inputTree)
  : QpGenStochAugRedPr(inputTree)
{};

QpStochAugRedPrCG::QpStochAugRedPrCG( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGenStochAugRedPr(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

QpStochAugRedPrCG::QpStochAugRedPrCG()
{ };

QpStochAugRedPrCG::~QpStochAugRedPrCG()
{ };


QpGenStochLinsysRoot* QpStochAugRedPrCG::newLinsysRoot()
{
  return new QpStochLinsysRootAugRedPrCG(this, data);
}

QpGenStochLinsysRoot* 
QpStochAugRedPrCG::newLinsysRoot(QpGenStochData* prob,
			  OoqpVector* dd,OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  new QpStochLinsysRootAugRedPrCG(this, prob, dd, dq, nomegaInv, rhs);
}
