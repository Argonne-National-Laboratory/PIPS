/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactory.h"

#include "QpGenStochData.h"
#include "StochTree.h"
#include "StochInputTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"

#include "sVars.h"
#include "sResiduals.h"

#include "sLinsysRoot.h"
#include "sLinsysLeaf.h"

#include "Ma57Solver.h"

sFactory::sFactory( StochInputTree* inputTree)
  :QpGen(0,0,0), data(NULL), m_tmTotal(0.0)
{
  
  tree = new StochTree(inputTree);
  tree->computeGlobalSizes();
  //now the sizes of the problem are available, set them
  tree->GetGlobalSizes(nx, my, mz);
  //decide how the CPUs are assigned 
  tree->assignProcesses();
}
 
sFactory::sFactory( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : QpGen( nx_, my_, mz_ ),
    nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_),
    tree(NULL), data(NULL), resid(NULL), linsys(NULL),
    m_tmTotal(0.0)
{ };

sFactory::sFactory()
  : QpGen( 0,0,0 ), m_tmTotal(0.0)
{ };

sFactory::~sFactory()
{
  if(tree) delete tree;
}

#define TIM t = MPI_Wtime();
#define REP(s) t = MPI_Wtime()-t;if (p) printf("makeData: %s took %f sec\n",s,t);

Data * sFactory::makeData()
{
  double t,t2=MPI_Wtime();
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  bool p = (mype == 0);

  // use the tree to get the data from user and to create OOQP objects
  TIM;
  StochSymMatrixHandle     Q( tree->createQ() );
  REP("Q");
  TIM;
  StochVectorHandle        c( tree->createc() );
  REP("c");
  TIM;
  StochVectorHandle     xlow( tree->createxlow()  );
  REP("xlow");
  TIM;
  StochVectorHandle    ixlow( tree->createixlow() );
  REP("ixlow");
  TIM;
  StochVectorHandle     xupp( tree->createxupp()  );
  REP("xupp");
  TIM;
  StochVectorHandle    ixupp( tree->createixupp() );
  REP("ixupp");
  TIM;
  StochGenMatrixHandle     A( tree->createA() );
  REP("A");
  TIM;
  StochVectorHandle        b( tree->createb() );
  REP("b");
  TIM;

  StochGenMatrixHandle     C( tree->createC()     );  
  REP("C");
  TIM;
  StochVectorHandle     clow( tree->createclow()  );
  REP("clow");
  TIM;
  StochVectorHandle    iclow( tree->createiclow() );
  REP("iclow");
  TIM;
  StochVectorHandle     cupp( tree->createcupp()  );
  REP("cupp");
  TIM;
  StochVectorHandle    icupp( tree->createicupp() );
  REP("icupp");

  MPI_Barrier(MPI_COMM_WORLD);
  t2 = MPI_Wtime() - t2;
  if (mype == 0) {
    cout << "IO second part took " << t2 << " sec\n";
  }

  data = new QpGenStochData( tree, 
			     c, Q, 
			     xlow, ixlow, ixlow->numberOfNonzeros(),
			     xupp, ixupp, ixupp->numberOfNonzeros(),
			     A, b,
			     C, clow, iclow, iclow->numberOfNonzeros(),
			     cupp, icupp, icupp->numberOfNonzeros() );
  return data;
}

Variables* sFactory::makeVariables( Data * prob_in )
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
  
  sVars* vars = new sVars( tree, x, s, y, z,
			   v, gamma, w, phi,
			   t, lambda, u, pi, 
			   prob->ixlow, prob->ixlow->numberOfNonzeros(),
			   prob->ixupp, prob->ixupp->numberOfNonzeros(),
			   prob->iclow, prob->iclow->numberOfNonzeros(),
			   prob->icupp, prob->icupp->numberOfNonzeros());
  registeredVars.push_back(vars);
  return vars;
}


Residuals* sFactory::makeResiduals( Data * prob_in )
{
  QpGenStochData* prob = dynamic_cast<QpGenStochData*>(prob_in);
  resid =  new sResiduals(tree, 
			  prob->ixlow, prob->ixupp,
			  prob->iclow, prob->icupp);
  return resid; 
}



LinearSystem* sFactory::makeLinsys( Data * prob_in )
{  
  linsys = newLinsysRoot();
  return linsys; 
}

sLinsysLeaf* 
sFactory::newLinsysLeaf()
{
  assert(false && "not supported");
  return NULL;
}
sLinsysLeaf* 
sFactory::newLinsysLeaf(QpGenStochData* prob,
			  OoqpVector* dd, OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  return new sLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs);
}

void sFactory::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			  OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  assert(0 && "not implemented here");
}

void
sFactory::separateVars( OoqpVector& x_in, OoqpVector& y_in,
			  OoqpVector& z_in, OoqpVector& vars_in )
{
  assert(0 && "not implemented here");
}

void sFactory::iterateStarted()
{
  iterTmMonitor.recIterateTm_start();
  tree->startMonitors();
}


void sFactory::iterateEnded()
{
  tree->stopMonitors();
  if(tree->balanceLoad()) {
    // balance needed
    data->sync();
      
    for(size_t i=0; i<registeredVars.size(); i++)
      registeredVars[i]->sync();
    
    resid->sync();
    
    linsys->sync();

    printf("Should not get here! OMG OMG OMG\n");
  }

  //logging and monitoring
  iterTmMonitor.recIterateTm_stop();
  m_tmTotal += iterTmMonitor.tmIterate;

  if(tree->rankMe==tree->rankZeroW) {
#ifdef TIMING
    extern double g_iterNumber;
    printf("TIME %g SOFAR %g ITER %d\n", iterTmMonitor.tmIterate, m_tmTotal, (int)g_iterNumber);
#else
    printf("ITERATION WALLTIME: iter=%g  Total=%g\n", iterTmMonitor.tmIterate, m_tmTotal);
#endif
  }
}

