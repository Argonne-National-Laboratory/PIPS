#include "QpStochNrmEqnPr.h"

#include "QpGenStochData.h"

#include "StochTreePrecond.h"
#include "StochInputTree.h"

#include "QpGenStochLinsysRootNrmEqnPrecond.h"

QpStochNrmEqnPr::QpStochNrmEqnPr( StochInputTree* inputTree)
  : QpGenStoch()
{ 
  tree = new StochTreePrecond(inputTree);
  tree->computeGlobalSizes();

  //now the sizes of the problem are available, set them
  tree->GetGlobalSizes(nx, my, mz);

  //decide how the CPUs are assigned 
  tree->assignProcesses(); 
};

QpStochNrmEqnPr::QpStochNrmEqnPr( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGenStoch(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

QpStochNrmEqnPr::QpStochNrmEqnPr()
{ };

QpStochNrmEqnPr::~QpStochNrmEqnPr()
{ };


QpGenStochLinsysRoot* QpStochNrmEqnPr::newLinsysRoot()
{
  return new QpGenStochLinsysRootNrmEqnPrecond(this, data);
}

QpGenStochLinsysRoot* 
QpStochNrmEqnPr::newLinsysRoot(QpGenStochData* prob,
			       OoqpVector* dd,OoqpVector* dq,
			       OoqpVector* nomegaInv, OoqpVector* rhs)
{
  new QpGenStochLinsysRootNrmEqnPrecond(this, prob, dd, dq, nomegaInv, rhs);
}
