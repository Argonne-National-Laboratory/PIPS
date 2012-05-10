#include "QpStochAugRedPrPCG.h"

#include "QpGenStochData.h"

#include "StochTreePrecond.h"
#include "StochInputTree.h"

#include "QpStochLinsysRootAugRedPrPCG.h"

QpStochAugRedPrPCG::QpStochAugRedPrPCG( StochInputTree* inputTree)
  : QpGenStochAugRedPr(inputTree)
{};

QpStochAugRedPrPCG::QpStochAugRedPrPCG( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGenStochAugRedPr(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

QpStochAugRedPrPCG::QpStochAugRedPrPCG()
{ };

QpStochAugRedPrPCG::~QpStochAugRedPrPCG()
{ };


QpGenStochLinsysRoot* QpStochAugRedPrPCG::newLinsysRoot()
{
  return new QpStochLinsysRootAugRedPrPCG(this, data);
}

QpGenStochLinsysRoot* 
QpStochAugRedPrPCG::newLinsysRoot(QpGenStochData* prob,
			  OoqpVector* dd,OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  new QpStochLinsysRootAugRedPrPCG(this, prob, dd, dq, nomegaInv, rhs);
}
