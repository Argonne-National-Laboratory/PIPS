#include "QpGenStochAugRedPr.h"

#include "QpGenStochData.h"

#include "StochTreePrecond.h"
#include "StochInputTree.h"

#include "QpGenStochLinsysRootAugRedPrecond.h"

QpGenStochAugRedPr::QpGenStochAugRedPr( StochInputTree* inputTree)
  : QpGenStoch()
{ 
  tree = new StochTreePrecond(inputTree);
  tree->computeGlobalSizes();

  //now the sizes of the problem are available, set them
  tree->GetGlobalSizes(nx, my, mz);

  //decide how the CPUs are assigned 
  tree->assignProcesses(); 

};

QpGenStochAugRedPr::QpGenStochAugRedPr( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGenStoch(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

QpGenStochAugRedPr::QpGenStochAugRedPr()
{ };

QpGenStochAugRedPr::~QpGenStochAugRedPr()
{ };


QpGenStochLinsysRoot* QpGenStochAugRedPr::newLinsysRoot()
{
  return new QpGenStochLinsysRootAugRedPrecond(this, data);
}

QpGenStochLinsysRoot* 
QpGenStochAugRedPr::newLinsysRoot(QpGenStochData* prob,
			  OoqpVector* dd,OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  new QpGenStochLinsysRootAugRedPrecond(this, prob, dd, dq, nomegaInv, rhs);
}
