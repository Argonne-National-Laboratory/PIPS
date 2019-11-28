#include "QpGenStochLinsysRootAugRedPrecond.h"
#include "QpGenStoch.h"
#include "QpGenStochData.h"
#include "QpGenStochLinsysLeaf.h"
#include "DeSymIndefSolver.h"
#include "DeSymPSDSolver.h"
#include "pipsport.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

static void addUtVToKKT(double alpha, DenseSymMatrix& UtV, DenseSymMatrix& kkt, int nx); 

QpGenStochLinsysRootAugRedPrecond::
QpGenStochLinsysRootAugRedPrecond(QpGenStoch * factory_, 
				  QpGenStochData * prob)
  : QpGenStochLinsysRootAugRed(), 
    Pmult(nullptr), Amult(nullptr), tmpVec1(nullptr)
{ 
  factory = factory_;
  kkt = nullptr;
  solver = nullptr;

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
  UtV = nullptr;

  //QpGenStochLinsyRootAug
  CtDC=nullptr;
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);

  //QpGenStochLinsyRootAugRed
  redRhs = new SimpleVector(locnx+locmy+locmz);

  //intializations related to this class 
  me = whoAmI();
  createChildren(prob);
};

QpGenStochLinsysRootAugRedPrecond::
QpGenStochLinsysRootAugRedPrecond(QpGenStoch* factory_,
				  QpGenStochData* prob,
				  OoqpVector* dd_, 
				  OoqpVector* dq_,
				  OoqpVector* nomegaInv_,
				  OoqpVector* rhs_)
  : QpGenStochLinsysRootAugRed(),
    Pmult(nullptr), Amult(nullptr), tmpVec1(nullptr)
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
  UtV = nullptr;

  //QpGenStochLinsyRootAug
  CtDC=nullptr;
  redRhs = new SimpleVector(locnx+locmy+locmz);
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);

  //intializations related to this class 
  me = whoAmI();
  createChildren(prob);
};

QpGenStochLinsysRootAugRedPrecond::~QpGenStochLinsysRootAugRedPrecond()
{
  if(Pmult)  delete Pmult;
  if(Amult)  delete Amult;
  if(tmpVec1)delete tmpVec1;
};

void QpGenStochLinsysRootAugRedPrecond::factor2(QpGenStochData *prob, 
						Variables *vars)
{
  assert( children.size() == prob->children.size() );
  double* buffer=nullptr;
  StochTreePrecond* stochNodePrcnd = dynamic_cast<StochTreePrecond*>(stochNode);
  //!!
  DenseSymMatrix * kktd = (DenseSymMatrix*) kkt; 
  myAtPutZeros(kktd, 0, 0, locnx+locmy, locnx+locmy);
  //~~

  // First tell children to factorize.
  for(int it=0; it<children.size(); it++) {
    children[it]->factor2(prob->children[it], vars);
  }

  if(me==ePrecond || me==eSpecialWorker) 
    buffer = new double[locnx*(locnx+locmy)];  

  DenseGenMatrix* U = nullptr;
  DenseGenMatrix* V = nullptr;
  

  //if(me==ePrecond) assert(children.size()==1);
  int commWrkrs = stochNode->commWrkrs;
  ////////////////////////////////////////////////////////
  // DIRECT workers -> all processes in fact
  ////////////////////////////////////////////////////////
  int childrenDone=0;
  for(int it=0; it<children.size(); it++) {

    if(children[it]->mpiComm == MPI_COMM_nullptr)
      continue;
    

    children[it]->stochNode->resMon.recFactTmChildren_start();    
    //-----------------------------------------------------------
    children[it]->allocU(&U, locnx);
    children[it]->allocV(&V, locnx);
    children[it]->computeU_V(prob->children[it], U, V);

    //-----------------------------------------------------------
    children[it]->stochNode->resMon.recSchurMultChildren_start();      
    //-----------------------------------------------------------
    kktd->matMult(-1.0, *U, 1, *V, 0, 1.0); 
    //-----------------------------------------------------------
    children[it]->stochNode->resMon.recSchurMultChildren_stop();
    children[it]->stochNode->resMon.recFactTmChildren_stop();
    childrenDone++;
    ///////////////////////////////////////////////////////////////
    // Stop and engage in communication with preconditioner if
    // enough scenarios were done
    ///////////////////////////////////////////////////////////////
    if(childrenDone==1) {
      int rankPrecond = stochNode->rankPrcnd;
      int rankZeroW   = stochNodePrcnd ->rankZeroW;
      int commP2ZeroW = stochNodePrcnd ->commP2ZeroW;

      if(me!=ePrecond) {
	stochNode->resMon.recFactTmLocal_start();
	///////////////////////////////////////
	// WORKERS  ->   send to precond
	///////////////////////////////////////
	MPI_Reduce(&(kktd->mStorage->M[0][0]), nullptr, locnx*(locnx+locmy), 
		   MPI_DOUBLE, MPI_SUM, rankPrecond, mpiComm);

	//null out Schur complement so the existing info will no more be added
	myAtPutZeros(kktd, 0, 0, locnx, locnx);
	stochNode->resMon.recFactTmLocal_stop();

	if(me==eSpecialWorker) {
	  MPI_Reduce(&(kktd->mStorage->M[0][0]), buffer, locnx*(locnx+locmy), 
		     MPI_DOUBLE, MPI_SUM, rankZeroW, commP2ZeroW);

	  memcpy(&kktd->mStorage->M[0][0], buffer, locnx*(locnx+locmy)*sizeof(double));
	}

      } else {
	////////////////////////////////////////////
	//PRECONDITIONER   ->  receive from workers
	////////////////////////////////////////////
	stochNode->resMon.recFactTmLocal_start();
	stochNode->resMon.recSchurMultLocal_start();
	if(U) delete U; if(V) delete V; 
	//deleteUtV(); reuse this
	stochNode->resMon.recSchurMultLocal_stop();

	MPI_Reduce(&(kktd->mStorage->M[0][0]), buffer, locnx*(locnx+locmy), 
		   MPI_DOUBLE, MPI_SUM, rankPrecond, mpiComm);
    
	memcpy(&kktd->mStorage->M[0][0], buffer, locnx*(locnx+locmy)*sizeof(double));
	delete[] buffer;


	//send the information back to specialWorker
	MPI_Reduce(&(kktd->mStorage->M[0][0]), nullptr, locnx*(locnx+locmy), 
		   MPI_DOUBLE, MPI_SUM, rankZeroW, commP2ZeroW);


	stochNode->resMon.recSchurMultLocal_start();  
	//////////////////////////////////////////////
	// factorize partial schur complement
	//////////////////////////////////////////////

	// update the upper block of the kkt with the UtV block
	int noProcs; MPI_Comm_size(mpiComm, &noProcs);
	double alpha = 1.0*children.size()/(noProcs*childrenDone);

	kktd->scalarMult(alpha);
	updateKKT(prob,vars);
	//addUtVToKKT(alpha, *UtV, *kktd, locnx); 

	//factorize
	double st=MPI_Wtime();
	solver->matrixChanged();
	printf("fact took %g\n", MPI_Wtime()-st);
	stochNode->resMon.recFactTmLocal_stop();
	stochNode->resMon.recSchurMultLocal_stop();  
      }
    }
  }


  if(me!=ePrecond) {
    //printf("Worker finished updates rank=%d\n", stochNode->rankMe);
    stochNode->resMon.recSchurMultLocal_start();
    if(U) delete U; if(V) delete V; 
    //deleteUtV(); reuse this
    stochNode->resMon.recSchurMultLocal_stop();
  }

  /////////////////////////////////////////////////////////
  // Everybody sum the partial Schur complements to
  // special worker who will have the complete matrix
  /////////////////////////////////////////////////////////    
  if(iAmDistrib) {
    int rankZeroW = stochNode->rankZeroW;
    MPI_Comm commWorkers = stochNodePrcnd ->commWorkers;
    if(me==eSpecialWorker) {

      //buffer=new double[locnx*locnx];
      //if(buffer==nullptr) printf("PANIC !!!! not enough memory in doing the reduce !!!!\n");

      MPI_Reduce(&(kktd->mStorage->M[0][0]), buffer, locnx*(locnx+locmy),
		 MPI_DOUBLE, MPI_SUM, 
		 rankZeroW, commWorkers);

      memcpy(&kktd->mStorage->M[0][0], buffer, locnx*(locnx+locmy)*sizeof(double));
      delete[] buffer;

      stochNode->resMon.recFactTmLocal_start();

      updateKKT(prob,vars);
      //addUtVToKKT(1.0, *UtV, *kktd, locnx);
      stochNode->resMon.recFactTmLocal_stop();

    } else {
      //printf("Nonzero worker %d -> reducing...\n", stochNode->rankMe);
      if(me!=ePrecond)
      MPI_Reduce(&(kktd->mStorage->M[0][0]), nullptr, locnx*(locnx+locmy),
		 MPI_DOUBLE, MPI_SUM, rankZeroW, commWorkers);

      //printf("Nonzero worker %d -> finished reducing\n", stochNode->rankMe);
    }
  }
}


void 
QpGenStochLinsysRootAugRedPrecond::Dsolve( QpGenStochData *prob, OoqpVector& b_ )
{
  StochVector& bst = dynamic_cast<StochVector&>(b_);
  for(int it=0; it<children.size(); it++) {
    children[it]->Dsolve(prob->children[it], *bst.children[it]);
  }

  int n = locnx+locmy; int N=n+locmz;
  

  SimpleVector& b = dynamic_cast<SimpleVector&>(*bst.vec);  assert(b.length()==N);

  int rankMe      = stochNode->rankMe;
  int rankPrcnd   = stochNode->rankPrcnd;
  int rankZeroW   = stochNode->rankZeroW;
  int commP2ZeroW = stochNode->commP2ZeroW;

  ///////////////////////////////////////////////////////////////////////
  // precond waits to be signaled
  ///////////////////////////////////////////////////////////////////////
  if(me==ePrecond) {
    if(nullptr==tmpVec1) tmpVec1 = new double[n+PREC_EXTRA_DOUBLES];

    SimpleVector rhs(tmpVec1, n);
    MPI_Status status;
    while(true) {
      MPI_Recv(tmpVec1, n+PREC_EXTRA_DOUBLES, MPI_DOUBLE,
	       rankZeroW, 1, mpiComm,
	       &status);
      if(tmpVec1[n+PREC_EXTRA_DOUBLES-1]==1) break;

      solver->solve(rhs);

      iErr = MPI_Send(&rhs[0], n, MPI_DOUBLE, 
		      rankZeroW, 2, mpiComm);

      assert(iErr==MPI_SUCCESS);
    }
  }  

  ///////////////////////////////////////////////////////////////////////
  // special worker 
  //   - signals the precond and waits for Px matVec via BiCGStabSolver
  //   - broadcasts the solution to the other workers
  ///////////////////////////////////////////////////////////////////////
  if(me==eSpecialWorker) {
    if(nullptr==tmpVec1) tmpVec1 = new double[n+PREC_EXTRA_DOUBLES];
    //allocate the buffer here reuse these buffer in the MatVec
    Pmult->setTmpVec1(tmpVec1);

    solveReduced(prob, b);

    // notify the preconditioner that we converged; 
    // the rest of buffer does not matter.
    tmpVec1[n+PREC_EXTRA_DOUBLES-1] = 1;
    MPI_Send(tmpVec1, n+PREC_EXTRA_DOUBLES, MPI_DOUBLE, rankPrcnd, 1, mpiComm);
  }

  ///////////////////////////////////////////////////////////////////////
  // the rest of the workers wait for special worker to converge
  // and obtain the solution
  ///////////////////////////////////////////////////////////////////////
  MPI_Bcast(&b[0], N, MPI_DOUBLE, rankZeroW, mpiComm);    

  //printf("DSolve rank[%d] finished!!!!!!!\n", rankMe); 
}

void 
QpGenStochLinsysRootAugRedPrecond::createChildren(QpGenStochData* prob)
{
  QpGenStochLinsys* child=nullptr;
  StochVector& ddst = dynamic_cast<StochVector&>(*dd);
  StochVector& dqst = dynamic_cast<StochVector&>(*dq);
  StochVector& nomegaInvst = dynamic_cast<StochVector&>(*nomegaInv);
  StochVector& rhsst = dynamic_cast<StochVector&>(*rhs);
  QpGenStoch* stochFactory = dynamic_cast<QpGenStoch*>(factory);

  ///////////////////////////////////////////////////
  // WORKERS process (all processes)
  ///////////////////////////////////////////////////
  
  //get the communicator from one of the vectors
  this->mpiComm = ddst.mpiComm;
  this->iAmDistrib = ddst.iAmDistrib;
  
  for(int it=0; it<prob->children.size(); it++) {
    
    if(MPI_COMM_nullptr == ddst.children[it]->mpiComm) {
      child = new QpGenStochDummyLinsys(dynamic_cast<QpGenStoch*>(factory), prob->children[it]);
    } else {
      
      if(prob->children[it]->children.size() == 0) {	
	//child = new QpGenStochLinsysLeaf(dynamic_cast<QpGenStoch*>(factory), 
	child = stochFactory->newLinsysLeaf(prob->children[it],
					    ddst.children[it],
					    dqst.children[it],
					    nomegaInvst.children[it],
					    rhsst.children[it]);
      } else {
	//child = new QpGenStochLinsysRoot(dynamic_cast<QpGenStoch*>(factory), 
	child = stochFactory->newLinsysRoot(prob->children[it],
					    ddst.children[it],
					    dqst.children[it],
					    nomegaInvst.children[it],
					    rhsst.children[it]);
      }
    }
    AddChild(child);
  }
}

SymMatrix*   
QpGenStochLinsysRootAugRedPrecond::createKKT(QpGenStochData* prob)
{
  int n = locnx+locmy;
  return new DenseSymMatrix(n);
}


DoubleLinearSolver* 
QpGenStochLinsysRootAugRedPrecond::createSolver(QpGenStochData* prob, 
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
      if(nullptr==Amult) Amult = new StoredMatTimesVec(kktmat);

      StochTreePrecond* stochNodePr = dynamic_cast<StochTreePrecond*>(stochNode);
      if(nullptr==Pmult) Pmult = new RemoteMatTimesVec(stochNodePr);

      return new BiCGStabSolver(Amult, Pmult);
      //return new CGSolver(Amult, Pmult);
      //return new DeSymIndefSolver(kktmat);
    } else {
      /////////////////////////////////////////////////////
      // Non-special worker
      /////////////////////////////////////////////////////
      return new DummyLinearSolver();
    }
  }
}

QpGenStochLinsysRootAugRedPrecond::NodeType QpGenStochLinsysRootAugRedPrecond::whoAmI()
{

  int& me=stochNode->rankMe; int& prcnd=stochNode->rankPrcnd; int& special=stochNode->rankZeroW;

  if(me==prcnd) {
    assert(stochNode->commWrkrs    != MPI_COMM_NULL);
    assert(stochNode->commP2ZeroW  != MPI_COMM_NULL);
    return ePrecond;
  } else {
    if(me==special) {
      assert(stochNode->commWrkrs    != MPI_COMM_NULL);
      assert(stochNode->commP2ZeroW  != MPI_COMM_NULL);
      return eSpecialWorker;
    } else { 
      assert(me>0 && me<prcnd); //non-special workers are 1,2 through P-1
      assert(stochNode->commWrkrs    != MPI_COMM_NULL);
      assert(stochNode->commP2ZeroW  == MPI_COMM_NULL);
      return eWorker;
    }
  }
}

void QpGenStochLinsysRootAugRedPrecond::initializeUtV()
{
  DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt);
  //put zeros
  myAtPutZeros(UtV, 0,0, locnx, locnx);
}


void addUtVToKKT(double alpha, DenseSymMatrix& UtV, DenseSymMatrix& kkt, int n)
{
  if(alpha!=0.0 && alpha!=1.0) {
    for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
	kkt[i][j] += ( alpha*UtV[i][j] );
  } else {
    if(alpha==1.0) {
      for(int i=0; i<n; i++)
	for(int j=0; j<n; j++)
	  kkt[i][j] += UtV[i][j];
    } else {
      //WTF!!! 
      //assert(false);
    }
  }
  

}
