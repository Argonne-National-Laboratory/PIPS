#include "QpGenStoch.h"

#include "QpGenStochData.h"
#include "StochTree.h"
#include "StochInputTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"

#include "QpGenStochVars.h"
#include "QpGenResiduals2.h"
#include "QpGenStochLinsysLeaf.h"


#include "QpGenStochLinsysRootAugRedPrecond.h"

#define ROOTAUG QpGenStochLinsysRootAugRedPrecond

#include "Ma57Solver.h"

QpGenStoch::QpGenStoch( StochInputTree* inputTree)
  :QpGen(0,0,0), m_tmTotal(0.0)
{
  tree = new StochTree(inputTree);
  tree->computeGlobalSizes();

  //now the sizes of the problem are available, set them
  tree->GetGlobalSizes(nx, my, mz);

  //decide how the CPUs are assigned 
  tree->assignProcesses(); 
}

QpGenStoch::QpGenStoch( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGen( nx_, my_, mz_ ),
    data(NULL),tree(NULL), resid(NULL), linsys(NULL),
    nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_), m_tmTotal(0.0)
{ };

QpGenStoch::QpGenStoch()
  : QpGen( 0,0,0 ), m_tmTotal(0.0)
{ };

QpGenStoch::~QpGenStoch()
{
  if(tree) delete tree;
}


Data * QpGenStoch::makeData()
{
  // use the tree to get the data from user and to create OOQP objects
  StochSymMatrixHandle     Q( tree->createQ() );
  StochVectorHandle        c( tree->createc() );

  StochVectorHandle     xlow( tree->createxlow()  );
  StochVectorHandle    ixlow( tree->createixlow() );
  StochVectorHandle     xupp( tree->createxupp()  );
  StochVectorHandle    ixupp( tree->createixupp() );

  StochGenMatrixHandle     A( tree->createA() );
  StochVectorHandle        b( tree->createb() );

  StochGenMatrixHandle     C( tree->createC()     );  
  StochVectorHandle     clow( tree->createclow()  );
  StochVectorHandle    iclow( tree->createiclow() );
  StochVectorHandle     cupp( tree->createcupp()  );
  StochVectorHandle    icupp( tree->createicupp() );

  data = new QpGenStochData( tree, 
			     c, Q, 
			     xlow, ixlow, ixlow->numberOfNonzeros(),
			     xupp, ixupp, ixupp->numberOfNonzeros(),
			     A, b,
			     C, clow, iclow, iclow->numberOfNonzeros(),
			     cupp, icupp, icupp->numberOfNonzeros() );
  
  return data;
}

Variables* QpGenStoch::makeVariables( Data * prob_in )
{
  QpGenStochData* prob = dynamic_cast<QpGenStochData*>(prob_in);


  OoqpVector * x      = tree->newPrimalVector();
  OoqpVector * s      = tree->newDualZVector();
  OoqpVector * y      = tree->newDualYVector();
  OoqpVector * z      = tree->newDualZVector();
  OoqpVector * v      = tree->newPrimalVector(); 
  OoqpVector * gamma  = tree->newPrimalVector();
  OoqpVector * w      = tree->newPrimalVector(); 
  OoqpVector * phi    = tree->newPrimalVector();
  OoqpVector * t      = tree->newDualZVector();
  OoqpVector * lambda = tree->newDualZVector();
  OoqpVector * u      = tree->newDualZVector(); 
  OoqpVector * pi     = tree->newDualZVector();
  OoqpVector * ixlow  = tree->newPrimalVector(); 
  OoqpVector * ixupp  = tree->newPrimalVector();
  OoqpVector * iclow  = tree->newDualZVector(); 
  OoqpVector * icupp  = tree->newDualZVector();

  //vars = new QpGenVars( x, s, y, z,
  QpGenStochVars* vars = new QpGenStochVars( tree, x, s, y, z,
			     v, gamma, w, phi,
			     t, lambda, u, pi, 
			     prob->ixlow, prob->ixlow->numberOfNonzeros(),
			     prob->ixupp, prob->ixupp->numberOfNonzeros(),
			     prob->iclow, prob->iclow->numberOfNonzeros(),
			     prob->icupp, prob->icupp->numberOfNonzeros());
  registeredVars.push_back(vars);

  return vars;
}


Residuals* QpGenStoch::makeResiduals( Data * prob_in )
{
  QpGenStochData* prob = dynamic_cast<QpGenStochData*>(prob_in);

  resid =  new QpGenResiduals2(tree, 
			       prob->ixlow, prob->ixupp,
			       prob->iclow, prob->icupp);
  return resid;
}



LinearSystem* QpGenStoch::makeLinsys( Data * prob_in )
{  
  QpGenStochData* prob = dynamic_cast<QpGenStochData*>(prob_in);
  linsys = newLinsysRoot();
  return linsys; 
}

QpGenStochLinsysLeaf* 
QpGenStoch::newLinsysLeaf()
{
  assert(false && "not supported");
  return NULL;
}
QpGenStochLinsysLeaf* 
QpGenStoch::newLinsysLeaf(QpGenStochData* prob,
			  OoqpVector* dd, OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  return new QpGenStochLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs);
}

void QpGenStoch::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			  OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  assert(0 && "not implemented here");
}

void
QpGenStoch::separateVars( OoqpVector& x_in, OoqpVector& y_in,
			  OoqpVector& z_in, OoqpVector& vars_in )
{
  assert(0 && "not implemented here");
}

void QpGenStoch::iterateStarted()
{
  iterTmMonitor.recIterateTm_start();
  tree->startMonitors();
}


void QpGenStoch::iterateEnded()
{
  tree->stopMonitors();
  if(tree->balanceLoad()) {
    // balance needed
    data->sync();
      
    for(int i=0; i<registeredVars.size(); i++)
      registeredVars[i]->sync();
    
    resid->sync();
    
    linsys->sync();

    printf("Should not get here! OMG OMG OMG\n");
  }

  //logging and monitoring
  iterTmMonitor.recIterateTm_stop();
  m_tmTotal += iterTmMonitor.tmIterate;

  if(tree->rankMe==tree->rankZeroW)
    printf("ITERATION WALLTIME: iter=%g  Total=%g\n", iterTmMonitor.tmIterate, m_tmTotal);
}

