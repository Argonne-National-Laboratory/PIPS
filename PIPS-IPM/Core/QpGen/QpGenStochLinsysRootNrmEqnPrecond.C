#include "QpGenStochLinsysRootNrmEqnPrecond.h"
#include "QpGenStoch.h"
#include "QpGenStochData.h"
#include "QpGenStochLinsysLeaf.h"
#include "DeSymPSDSolver.h"
#include "pipsport.h"


void addUtVToKKT(double alpha, DenseSymMatrix& UtV, DenseSymMatrix& kkt, int nx); 

QpGenStochLinsysRootNrmEqnPrecond::
QpGenStochLinsysRootNrmEqnPrecond(QpGenStoch * factory_, 
				  QpGenStochData * prob)
  : QpGenStochLinsysRootNrmEqn(),
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

  //QpGenStochLinsyRootAugRed
  CtDC=nullptr;
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);

  //QpGenStochLinsysRootNrmEqn parent
  AQinvAt=nullptr; solver2=nullptr;

  //intializations related to this class 
  me = whoAmI();
  createChildren(prob);
};

QpGenStochLinsysRootNrmEqnPrecond::QpGenStochLinsysRootNrmEqnPrecond(QpGenStoch* factory_,
					   QpGenStochData* prob,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : QpGenStochLinsysRootNrmEqn(),
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

    //QpGenStochLinsyRootAugRed
  CtDC=nullptr;
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);

  //QpGenStochLinsysRootNrmEqn parent
  AQinvAt=nullptr; solver2=nullptr;

  //intializations related to this class 
  me = whoAmI();
  createChildren(prob);
};

QpGenStochLinsysRootNrmEqnPrecond::~QpGenStochLinsysRootNrmEqnPrecond()
{
  if(AQinvAt) delete AQinvAt;
  if(solver2) delete solver2;
  delete redRhs;
}


void QpGenStochLinsysRootNrmEqnPrecond::factor2(QpGenStochData *prob, Variables *vars)
{
  assert( children.size() == prob->children.size() );
  double* buffer=nullptr;

  this->updateKKT(prob, vars);

  // First tell children to factorize.
  for(int it=0; it<children.size(); it++) {
    children[it]->factor2(prob->children[it], vars);
  }
  
  if(me==ePrecond) 
    buffer = new double[locnx*locnx]; 
  DenseGenMatrix* U = nullptr;
  DenseGenMatrix* V = nullptr;
  
  stochNode->resMon.recSchurMultLocal_start();
  allocUtV();  initializeUtV();
  stochNode->resMon.recSchurMultLocal_stop();

  int commWrkrs = stochNode->commWrkrs;
  int childrenDone=0;

  for(int it=0; it<children.size(); it++) {

    if(children[it]->mpiComm == MPI_COMM_nullptr) continue;

    children[it]->stochNode->resMon.recFactTmChildren_start();    
    children[it]->allocU(&U, locnx); children[it]->allocV(&V, locnx);
    children[it]->computeU_V(prob->children[it], U, V);
    children[it]->stochNode->resMon.recSchurMultChildren_start();
    UtV->matMult(-1.0, *U, 1, *V, 0, 1.0); 
    children[it]->stochNode->resMon.recSchurMultChildren_stop();
    children[it]->stochNode->resMon.recFactTmChildren_stop();
    childrenDone++;
    ///////////////////////////////////////////////////////////////
    // Stop and engage in communication with preconditioner if
    // enough scenarios were done
    ///////////////////////////////////////////////////////////////   
    if(childrenDone==1) {
      int rankPrecond = stochNode->rankPrcnd;
      if(me!=ePrecond) {
	///////////////////////////////////////
	// WORKERS  ->   send to precond
	///////////////////////////////////////
	MPI_Reduce(&(UtV->mStorage->M[0][0]), nullptr, locnx*locnx, 
		   MPI_DOUBLE, MPI_SUM, 
		   rankPrecond, mpiComm);
	//null out Schur complement so the existing info will no more be added
	myAtPutZeros(UtV, 0, 0, locnx, locnx);
      } else {
	////////////////////////////////////////////
	//PRECONDITIONER   ->  receive from workers
	////////////////////////////////////////////
	stochNode->resMon.recSchurMultLocal_start();
	if(U) delete U; if(V) delete V; 
	stochNode->resMon.recSchurMultLocal_stop();

	MPI_Reduce(&(UtV->mStorage->M[0][0]), buffer, locnx*locnx, 
		   MPI_DOUBLE, MPI_SUM, 
		   rankPrecond, mpiComm);	
	memcpy(&UtV->mStorage->M[0][0], buffer, locnx*locnx*sizeof(double));
	delete[] buffer;

	//////////////////////////////////////////////
	// factorize partial schur complement
	//////////////////////////////////////////////
	stochNode->resMon.recSchurMultLocal_start(); //---1--- 
	int noProcs; MPI_Comm_size(mpiComm, &noProcs);
	double alpha = 1.0*children.size()/(noProcs*childrenDone);
	printf("Regulariz scale in Precond =%g\n", alpha);
	// update the upper block of the kkt with the UtV block
	DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt);
	addUtVToKKT(alpha, *UtV, *kktd, locnx); 


	double st=MPI_Wtime();
	solver->matrixChanged(); printf("fact took %g\n", MPI_Wtime()-st);

	if(locmy) { //perform further Normal Eqns reduction 
	  assert(locmy<=locnx);
	  SparseGenMatrix& A = prob->getLocalB();
	  DenseGenMatrix W(&(UtV->mStorage->M[0][0]), locnx, locmy); //reuse UtV
	  myDenseFromSparseTrans(W,A);
	  //printf("Dense to sparse done  %g sec\n", MPI_Wtime()-st);st=MPI_Wtime();
	  
	  DeSymPSDSolver& cholSolver = dynamic_cast<DeSymPSDSolver&>(*solver);
	  cholSolver.Lsolve(W);
	  //printf("W computed  %g sec\n"   , MPI_Wtime()-st);st=MPI_Wtime(); 
	  
	  //rank-k update
	  if(AQinvAt==nullptr) AQinvAt = new DenseSymMatrix(locmy);
	  AQinvAt->atRankkUpdate(0.0, 1.0, W, 1);
	  //printf("rank-k update done  %g sec\\n" , MPI_Wtime()-st);st=MPI_Wtime();   
	  
	  if(solver2==nullptr) solver2 = new DeSymPSDSolver(AQinvAt);
	  solver2->matrixChanged();

	  stochNode->resMon.recSchurMultLocal_stop();  //~~~1~~~ 
	}//~end locmy non-zero
	
      }//~end preconditioner stuff
    }//~end childrenDone==1
  }

  if(me!=ePrecond) {
    stochNode->resMon.recSchurMultLocal_start();
    delete U; delete V;
    //deleteUtV(); reuse this
    stochNode->resMon.recSchurMultLocal_stop();
  }


  if(iAmDistrib) {
    int rankZeroW = stochNode->rankZeroW;
    if(me==eSpecialWorker) {
      double* buffer=new double[locnx*locnx];
      if(buffer==nullptr) printf("PANIC !!!! not enough memory in doing the reduce !!!!\n");

      MPI_Reduce(&(UtV->mStorage->M[0][0]), buffer, locnx*locnx,
		 MPI_DOUBLE, MPI_SUM, rankZeroW, mpiComm);

      memcpy(&UtV->mStorage->M[0][0], buffer, locnx*locnx*sizeof(double));
      delete[] buffer;
      
      //special worker add UtV to KKT 
      //!optimize
      //   - kkt was set to zero, not used, right? use the space for UtV
      //   - do a memcpy not an addition
      DenseSymMatrix* kktd = dynamic_cast<DenseSymMatrix*>(kkt);
      addUtVToKKT(1.0, *UtV, *kktd, locnx);
    } else {
      // REGULAR worker and preconditioner
      MPI_Reduce(&(UtV->mStorage->M[0][0]), nullptr, locnx*locnx,
		 MPI_DOUBLE, MPI_SUM, 
		 rankZeroW, mpiComm);
    }
  } 



  /////////////////////////////////////////////////////////////////////////
  // NORMAL EQUATIONS stuff
  ////////////////////////////////////////////////////////////////////////
  //--------------------------------------------------------
  // compute A * inverse(Q) * A^T as 
  //  1) Q = LL' so AQinvAt=(L\A')' (L\A')
  //  2) W = L\A'    W is locnx-locmy matrix
  //  3) AQinvAt = W'*W (as a rank-k update)
  //--------------------------------------------------------
  
  //\\ the code is executed only by preconditioner, see above
}

void
QpGenStochLinsysRootNrmEqnPrecond::Dsolve( QpGenStochData *prob, OoqpVector& b_ )
{
  StochVector& bst = dynamic_cast<StochVector&>(b_);
  for(int it=0; it<children.size(); it++) {
    children[it]->Dsolve(prob->children[it], *bst.children[it]);
  }
  
  int n = locnx; int N=locnx+locmy+locmz;
  SimpleVector& b = dynamic_cast<SimpleVector&>(*bst.vec);  assert(b.length()==N);

  int rankMe      = stochNode->rankMe;
  int rankPrcnd   = stochNode->rankPrcnd;
  int rankZeroW   = stochNode->rankZeroW;
  int commP2ZeroW = stochNode->commP2ZeroW;

  ///////////////////////////////////////////////////////////////////////
  // PRECONDITIONER waits to be signaled
  ///////////////////////////////////////////////////////////////////////
  if(me==ePrecond) {
    if(nullptr==tmpVec1) tmpVec1 = new double[n+PREC_EXTRA_DOUBLES];

    SimpleVector rhs(tmpVec1, n); //shortcut
    MPI_Status status;
    while(true) {
      MPI_Recv(tmpVec1, n+PREC_EXTRA_DOUBLES, MPI_DOUBLE,
	       rankZeroW, 1, mpiComm,
	       &status);
      if(tmpVec1[n+PREC_EXTRA_DOUBLES-1]==1) break; //exit the loop on convergence

      solver->solve(rhs);

      iErr = MPI_Send(&rhs[0], n, MPI_DOUBLE, 
		      rankZeroW, 2, mpiComm);

      assert(iErr==MPI_SUCCESS);
    }
  }

  ///////////////////////////////////////////////////////////////////////
  // SPECUAL worker 
  //   - signals the precond and waits for Px matVec via CGSolver
  //   - broadcasts the solution to the other workers
  ///////////////////////////////////////////////////////////////////////
  if(me==eSpecialWorker) {
    if(nullptr==tmpVec1) tmpVec1 = new double[n+PREC_EXTRA_DOUBLES];
    //allocate the buffer here reuse these buffer in the MatVec
    Pmult->setTmpVec1(tmpVec1);

    solveReduced(prob, b);

    //notify the preconditioner that we converged; the rest of buffer does not matter.
    tmpVec1[n+PREC_EXTRA_DOUBLES-1] = 1;
    MPI_Send(tmpVec1, n+PREC_EXTRA_DOUBLES, MPI_DOUBLE, rankPrcnd, 1, mpiComm);
  }

  ///////////////////////////////////////////////////////////////////////
  // SPECIAL worker broadcasts the solution
  // the OTHERs also wait on this call for the convergence to be achieved.
  ///////////////////////////////////////////////////////////////////////
  MPI_Bcast(&b[0], N, MPI_DOUBLE, rankZeroW, mpiComm);    
}

void 
QpGenStochLinsysRootNrmEqnPrecond::createChildren(QpGenStochData* prob)
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
QpGenStochLinsysRootNrmEqnPrecond::createKKT(QpGenStochData* prob)
{
  return new DenseSymMatrix(locnx);
}


DoubleLinearSolver* 
QpGenStochLinsysRootNrmEqnPrecond::createSolver(QpGenStochData* prob, 
						SymMatrix* kktmat_)
{
  if(stochNode->rankMe==stochNode->rankPrcnd) {
    /////////////////////////////////////////////////////
    // Preconditioner
    /////////////////////////////////////////////////////
    DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
    return new DeSymPSDSolver(kktmat);
  } else {
    if(stochNode->rankMe==stochNode->rankZeroW) {
      /////////////////////////////////////////////////////
      // Special worker
      /////////////////////////////////////////////////////
      DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
      if(nullptr==Amult) Amult = new StoredMatTimesVec(kktmat);

      StochTreePrecond* stochNodePr = dynamic_cast<StochTreePrecond*>(stochNode);
      if(nullptr==Pmult) Pmult = new RemoteMatTimesVec(stochNodePr);

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

QpGenStochLinsysRootNrmEqnPrecond::NodeType QpGenStochLinsysRootNrmEqnPrecond::whoAmI()
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

void QpGenStochLinsysRootNrmEqnPrecond::initializeUtV()
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
