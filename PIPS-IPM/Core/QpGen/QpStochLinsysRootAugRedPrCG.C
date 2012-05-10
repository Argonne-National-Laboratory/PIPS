#include "QpStochLinsysRootAugRedPrCG.h"
#include "QpGenStoch.h"
#include "QpGenStochData.h"
#include "DeSymIndefSolver.h"

QpStochLinsysRootAugRedPrCG::QpStochLinsysRootAugRedPrCG(QpGenStoch * factory_, 
							 QpGenStochData * prob)
  : QpGenStochLinsysRootAugRedPrecond()
{ 
  factory = factory_;
  kkt = NULL;
  solver = NULL;

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  if( nxupp + nxlow > 0 ) {
    dd      = factory_->tree->newPrimalVector();
    dq      = factory_->tree->newPrimalVector();
    prob->getDiagonalOfQ( *dq );
  }
  nomegaInv   = factory_->tree->newDualZVector();
  rhs         = factory_->tree->newRhs();

  useRefs=0;
  data = prob;
  stochNode = prob->stochNode;

  //QpGenStochLinsysRoot
  iAmDistrib = 0;

  //QpGenStochLinsyRootAug
  prob->getLocalSizes(locnx, locmy, locmz);
  UtV = NULL;

  //QpGenStochLinsyRootAug
  CtDC=NULL;
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);

  //QpGenStochLinsyRootAugRed
  redRhs = new SimpleVector(locnx+locmy+locmz);

  //QpGenStochLinsyRootAugPrecond
  tmpVec1=NULL;

  //intializations related to this class 
  me = whoAmI();
  createChildren(prob);
};

QpStochLinsysRootAugRedPrCG::
QpStochLinsysRootAugRedPrCG(QpGenStoch* factory_,
				  QpGenStochData* prob,
				  OoqpVector* dd_, 
				  OoqpVector* dq_,
				  OoqpVector* nomegaInv_,
				  OoqpVector* rhs_)
  : QpGenStochLinsysRootAugRedPrecond()
{ 
  //QpGenStochLinsy
  factory = factory_;

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = prob->nxlow;
  nxupp = prob->nxupp;
  mclow = prob->mclow;
  mcupp = prob->mcupp;

  if( nxupp + nxlow > 0 ) {
    //dd      = OoqpVectorHandle(dd_);
    dd= dd_;
    //dq      = OoqpVectorHandle(dq_);
    dq = dq_;
  }
  //nomegaInv   = OoqpVectorHandle(nomegaInv_);
  //rhs         = OoqpVectorHandle(rhs_);
  nomegaInv = nomegaInv_;
  rhs = rhs_;

  useRefs=1;
  data = prob;
  stochNode = prob->stochNode;

  //QpGenStochLinsysRoot
  iAmDistrib = 0;

  //QpGenStochLinsyRootAug
  prob->getLocalSizes(locnx, locmy, locmz);
  UtV = NULL;

  //QpGenStochLinsyRootAug
  CtDC=NULL;
  redRhs = new SimpleVector(locnx+locmy+locmz);
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);

  //QpGenStochLinsyRootAugPrecond
  tmpVec1=NULL;

  //intializations related to this class 
  me = whoAmI();
  createChildren(prob);
};

QpStochLinsysRootAugRedPrCG::~QpStochLinsysRootAugRedPrCG()
{ };

SymMatrix*   
QpStochLinsysRootAugRedPrCG::createKKT(QpGenStochData* prob)
{
  assert(locmy==0);
  return new DenseSymMatrix(locnx);
}


DoubleLinearSolver* 
QpStochLinsysRootAugRedPrCG::createSolver(QpGenStochData* prob, 
						SymMatrix* kktmat_)
{
  if(stochNode->rankMe==stochNode->rankPrcnd) {
    /////////////////////////////////////////////////////
    // Preconditioner
    /////////////////////////////////////////////////////

    DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
    return new DeSymIndefSolver(kktmat);
    //return new DeSymPSDSolver(kktmat);

  } else {
    if(stochNode->rankMe==stochNode->rankZeroW) {
      /////////////////////////////////////////////////////
      // Special worker
      /////////////////////////////////////////////////////
      DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
      if(NULL==Amult) Amult = new StoredMatTimesVec(kktmat);

      StochTreePrecond* stochNodePr = dynamic_cast<StochTreePrecond*>(stochNode);
      if(NULL==Pmult) Pmult = new RemoteMatTimesVec(stochNodePr);

      //return new BiCGStabSolver(Amult, Pmult);
      return new CGSolver(Amult, Pmult);
      //return new DeSymIndefSolver(kktmat);
    } else {
      /////////////////////////////////////////////////////
      // Non-special worker
      /////////////////////////////////////////////////////
      return new DummyLinearSolver();
    }
  }
}
