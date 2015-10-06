#include "NlpGenStoch.h"

#include "NlpGenStochData.h"
/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "StochTree.h"
#include "StochInputTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"

#include "sNlpVars.h"
#include "NlpGenResiduals2.h"
#include "NlpGenStochLinsysLeaf.h"


#include "NlpGenStochLinsysRootAugRedPrecond.h"

#define ROOTAUG NlpGenStochLinsysRootAugRedPrecond

#include "Ma27Solver.h"
#include "Ma57Solver.h"

NlpGenStoch::NlpGenStoch( StochInputTree* inputTree)
  :NlpGen(0,0,0), m_tmTotal(0.0)
{
  tree = new StochTree(inputTree);
  tree->computeGlobalSizes();

  //now the sizes of the problem are available, set them
  tree->GetGlobalSizes(nx, my, mz);

  //decide how the CPUs are assigned 
  tree->assignProcesses(); 
}

NlpGenStoch::NlpGenStoch( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : NlpGen( nx_, my_, mz_ ),
    data(NULL),tree(NULL), resid(NULL), linsys(NULL),
    nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_), m_tmTotal(0.0)
{ };

NlpGenStoch::NlpGenStoch()
  : NlpGen( 0,0,0 ), m_tmTotal(0.0)
{ };

NlpGenStoch::~NlpGenStoch()
{
  if(tree) delete tree;
}


Data * NlpGenStoch::makeData()
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

  data = new NlpGenStochData( tree, 
			     c, Q, 
			     xlow, ixlow, ixlow->numberOfNonzeros(),
			     xupp, ixupp, ixupp->numberOfNonzeros(),
			     A, b,
			     C, clow, iclow, iclow->numberOfNonzeros(),
			     cupp, icupp, icupp->numberOfNonzeros() );
  
  return data;
}

Data * NlpGenStoch::makeData(NlpInfo *updateNlp)
{
  data = dynamic_cast<NlpGenStochData*> (makeData());
	
  data->SetInputNlpPara(updateNlp);
	
  data->inputNlp = updateNlp; 
	
  return data;
}


Variables* NlpGenStoch::makeVariables( Data * prob_in )
{
  NlpGenStochData* prob = dynamic_cast<NlpGenStochData*>(prob_in);


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

  sNlpVars* vars = new sNlpVars( tree, x, s, y, z,
			   v, gamma, w, phi,
			   t, lambda, u, pi, 
			   prob->ixlow, prob->ixlow->numberOfNonzeros(),
			   prob->ixupp, prob->ixupp->numberOfNonzeros(),
			   prob->iclow, prob->iclow->numberOfNonzeros(),
			   prob->icupp, prob->icupp->numberOfNonzeros());
  registeredVars.push_back(vars);

  return vars;
}


Residuals* NlpGenStoch::makeResiduals( Data * prob_in )
{
  NlpGenStochData* prob = dynamic_cast<NlpGenStochData*>(prob_in);

  resid =  new NlpGenResiduals2(tree, 
			       prob->ixlow, prob->ixupp,
			       prob->iclow, prob->icupp);
  return resid;
}



LinearSystem* NlpGenStoch::makeLinsys( Data * prob_in )
{  
  NlpGenStochData* prob = dynamic_cast<NlpGenStochData*>(prob_in);
  linsys = newLinsysRoot();
  return linsys; 
}

NlpGenStochLinsysLeaf* 
NlpGenStoch::newLinsysLeaf()
{
  assert(false && "not supported");
  return NULL;
}
NlpGenStochLinsysLeaf* 
NlpGenStoch::newLinsysLeaf(NlpGenStochData* prob,
			  OoqpVector* dd, OoqpVector* dq,
			  OoqpVector* nomegaInv, OoqpVector* rhs)
{
  return new NlpGenStochLinsysLeaf(this, prob, dd, dq, nomegaInv, rhs);
}

void NlpGenStoch::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			  OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  assert(0 && "not implemented here");
}

void
NlpGenStoch::separateVars( OoqpVector& x_in, OoqpVector& y_in,
			  OoqpVector& z_in, OoqpVector& vars_in )
{
  assert(0 && "not implemented here");
}

void NlpGenStoch::joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			  OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in )
{
  assert(0 && "not implemented here");
}

void
NlpGenStoch::separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in, OoqpVector& y_in,
			  OoqpVector& z_in, OoqpVector& vars_in )
{
  assert(0 && "not implemented here");
}

void NlpGenStoch::iterateStarted()
{
  iterTmMonitor.recIterateTm_start();
  tree->startMonitors();
}


void NlpGenStoch::iterateEnded()
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
#ifdef STOCH_TESTING
  if(tree->rankMe==tree->rankZeroW)
    printf("ITERATION WALLTIME: iter=%g  Total=%g\n", iterTmMonitor.tmIterate, m_tmTotal);
#endif
}


