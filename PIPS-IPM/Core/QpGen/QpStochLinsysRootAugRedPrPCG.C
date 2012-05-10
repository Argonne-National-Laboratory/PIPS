#include "QpStochLinsysRootAugRedPrPCG.h"
#include "QpGenStoch.h"
#include "QpGenStochData.h"
#include "DeSymIndefSolver.h"
#include "PCGSolver.h"
#include "QpGenStochData.h"
#include "RemoteMatTimesVec.h"
#include "Ma57Solver.h"
#include "Ma27Solver.h"

QpStochLinsysRootAugRedPrPCG::QpStochLinsysRootAugRedPrPCG(QpGenStoch * factory_, QpGenStochData * prob)
  : QpGenStochLinsysRootAugRed(),
    Pmult(NULL), Qmult(NULL), Atmult(NULL), tmpVec1(NULL)
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

  //intializations related to this class 
  AAt = NULL;
  AAtSolver = NULL;
  me = whoAmI();
  createChildren(prob);
};

QpStochLinsysRootAugRedPrPCG::QpStochLinsysRootAugRedPrPCG(QpGenStoch* factory_,
					   QpGenStochData* prob,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : QpGenStochLinsysRootAugRed(),
    Pmult(NULL), Qmult(NULL), Atmult(NULL), tmpVec1(NULL)
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

  //intializations related to this class 
  AAt = NULL;
  AAtSolver = NULL;
  me = whoAmI();
  createChildren(prob);
};

QpStochLinsysRootAugRedPrPCG::~QpStochLinsysRootAugRedPrPCG()
{
  if(AAt) delete AAt;
  if(AAtSolver) delete AAtSolver;
}


SymMatrix* 
QpStochLinsysRootAugRedPrPCG::createKKT(QpGenStochData* prob)
{
  if(stochNode->rankMe!=stochNode->rankPrcnd) 
    return new DenseSymMatrix(locnx);
  else 
    return new DenseSymMatrix(locnx+locmy);
}


void QpStochLinsysRootAugRedPrPCG::factor2(QpGenStochData *prob, Variables *vars)
{
  assert( children.size() == prob->children.size() );
  StochTreePrecond* stochNodePrcnd = dynamic_cast<StochTreePrecond*>(stochNode);
  //!!
  DenseSymMatrix * kktd = (DenseSymMatrix*) kkt;
  if(me==ePrecond) myAtPutZeros(kktd, 0, 0, locnx+locmy, locnx+locmy);
  else             myAtPutZeros(kktd, 0, 0, locnx,       locnx);

  // First tell children to factorize.
  for(int it=0; it<children.size(); it++) {
    children[it]->factor2(prob->children[it], vars);
  }
  
  DenseGenMatrix* U = NULL;
  DenseGenMatrix* V = NULL;

  int childrenDone=0;

  for(int it=0; it<children.size(); it++) {
    // Allocate  U if necessary. If the sizes of 
    // the previous child are good, reuse the mem.

    if(children[it]->mpiComm == MPI_COMM_NULL)
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
      //kktd is locnx-by-locnx on workers, and 
      //(locnx+locmy)-(locnx+locmy) on prcnd. 
      int rankPrecond = stochNodePrcnd ->rankPrcnd;
      int rankZeroW   = stochNodePrcnd ->rankZeroW;
      int commP2ZeroW = stochNodePrcnd ->commP2ZeroW;

      if(me!=ePrecond) {
	stochNode->resMon.recFactTmLocal_start();
	///////////////////////////////////////
	// WORKERS  ->   send to precond
	///////////////////////////////////////
	workerReduce_source(kktd, rankPrecond, mpiComm);

	//null out Schur complement so the existing info will no more be added
	myAtPutZeros(kktd, 0, 0, locnx, locnx);
	//Special worker gets  the info from prcnd
	if(me==eSpecialWorker) workerReduce_target(kktd, rankZeroW, commP2ZeroW);

	stochNode->resMon.recFactTmLocal_stop();
      }
      else { 
	////////////////////////////////////////////
	//PRECONDITIONER   ->  receive from workers
	////////////////////////////////////////////
	stochNode->resMon.recFactTmLocal_start();
	if(U) delete U; if(V) delete V; 
	stochNode->resMon.recSchurMultLocal_stop();

	precndReduce_target(kktd, rankPrecond, mpiComm);

	//send back to the special worker this info before it gets 
	//altered by the factorization phase.
	precndReduce_source(kktd, rankZeroW, commP2ZeroW);

	//"extrapolate" the partial scenarios to full information by
	//multiplying with alpha
	int noProcs; MPI_Comm_size(mpiComm, &noProcs);
	double alpha = 1.0*children.size()/(noProcs*childrenDone);
	kktd->scalarMult(alpha);
	//add Q and insert A
	updateKKT(prob, vars);
	//factorize
	double st=MPI_Wtime(); solver->matrixChanged();
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
  // workers sum the partial Schur complements to
  // special worker who will have the complete matrix
  /////////////////////////////////////////////////////////    
  if(iAmDistrib) {
    MPI_Comm commWorkers = stochNodePrcnd ->commWorkers;
    int rankZeroW = stochNodePrcnd->rankZeroW;
#ifdef DEBUG
    if(me!=ePrecond) assert(commWorkers!=MPI_COMM_NULL);
#endif
    if(me==eSpecialWorker) {
      workerReduce_target(kktd, rankZeroW, commWorkers);
      updateKKT(prob, vars);
      if(AAt==NULL) {
	assert(AAtSolver==NULL);
	SparseGenMatrix& A = prob->getLocalB();
	
	A.matMultTrans(&AAt);
	
	SparseSymMatrix* AAtsp = dynamic_cast<SparseSymMatrix*>(AAt);
	AAtsp->reduceToLower();
	
	AAtSolver = new Ma57Solver(AAtsp);
	AAtSolver->matrixChanged();
      }
    } else {
      if(me!=ePrecond)
	workerReduce_source(kktd, rankZeroW, commWorkers);
    }
  }
}

void QpStochLinsysRootAugRedPrPCG::Dsolve( QpGenStochData *prob, OoqpVector& x )
{
  StochVector& bst = dynamic_cast<StochVector&>(x);
  for(int it=0; it<children.size(); it++) {
    children[it]->Dsolve(prob->children[it], *bst.children[it]);
  }
  int n = locnx+locmy; int N=n+locmz;

  SimpleVector& b = dynamic_cast<SimpleVector&>(*bst.vec); 
  int rankMe      = stochNode->rankMe;
  int rankPrcnd   = stochNode->rankPrcnd;
  int rankZeroW   = stochNode->rankZeroW;
  int commP2ZeroW = stochNode->commP2ZeroW;
  ///////////////////////////////////////////////////////////////////////
  // precond waits to be signaled to apply preconditioner
  ///////////////////////////////////////////////////////////////////////
  if(me==ePrecond) {
    if(NULL==tmpVec1) tmpVec1 = new double[n+PREC_EXTRA_DOUBLES];

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
    if(NULL==tmpVec1) tmpVec1 = new double[n+PREC_EXTRA_DOUBLES];
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


void QpStochLinsysRootAugRedPrPCG::solveReduced( QpGenStochData *prob, SimpleVector& b)
{
  assert(locnx+locmy+locmz==b.length());
  SimpleVector& r = (*redRhs);
  assert(r.length() <= b.length());
  SparseGenMatrix& C = prob->getLocalD();
  SparseGenMatrix& A = prob->getLocalB();

  stochNode->resMon.recDsolveTmLocal_start();

  ///////////////////////////////////////////////////////////////////////
  // LOCAL SOLVE
  ///////////////////////////////////////////////////////////////////////
 
  ///////////////////////////////////////////////////////////////////////
  // b=[b1;b2;b3] is a locnx+locmy+locmz vector 
  // the new rhs should be 
  //           r = [b1-C^T*(zDiag)^{-1}*b3; b2]
  ///////////////////////////////////////////////////////////////////////

  r.copyFromArray(b.elements()); //will copy only as many elems as r has

  // aliases to parts (no mem allocations)
  SimpleVector r3(&r[locnx+locmy], locmz); //r3 is used as a temp
                                           //buffer for b3
  SimpleVector r1(&r[0],           locnx);
  assert(r3.length() == zDiag->length());


  
  ///////////////////////////////////////////////////////////////////////
  // compute r1 = b1 - C^T*(zDiag)^{-1}*b3
  ///////////////////////////////////////////////////////////////////////
  r3.componentDiv(*zDiag);//r3 is a copy of b3
  C.transMult(1.0, r1, -1.0, r3);

  ///////////////////////////////////////////////////////////////////////
  // r contains all the stuff -> solve for it
  ///////////////////////////////////////////////////////////////////////
  SimpleVector rx(locnx); rx.copyFrom(r1);

  //---
  SimpleVector rr1(locnx); rr1.copyFrom(r1);
  SimpleVector rr2(locmy); rr2.copyFromArray(&b[locnx]);
  //---

  SimpleVector realRhs(&r[0], locnx+locmy);
  //cout << "PCG: xin rhs\n"; realRhs.writeToStream(cout);
  //SimpleVector yy(&realRhs[locnx], locmy); cout << "PCG: yin rhs\n"; yy.writeToStream(cout);

  ///////////////////////////////////////////////////////////////////////
  // Solve for x using Proj CG
  ///////////////////////////////////////////////////////////////////////
  solver->Dsolve(realRhs);

  //cout << "PCG:xout rhs\n"; realRhs.writeToStream(cout);
  SimpleVector y(&r[locnx], locmy);
  SimpleVector x(&r[0], locnx);
  //y.negate();
  ///////////////////////////////////////////////////////////////////////
  // solve for y from  AA'y=A(rx-H*x)
  ///////////////////////////////////////////////////////////////////////
  kkt->mult(1.0, rx, -1.0, x);
  A.mult(0.0, y, 1.0, rx);
  AAtSolver->solve(y);  


 
  //cout << "PCG:yout rhs\n"; y.writeToStream(cout);

  //!delete AAt;


  //---
  A.transMult(1.0, rr1, -1.0, y);
  kkt->mult(1.0, rr1, -1.0, x);
  A.mult(1.0, rr2, -1.0, x);
  printf("Actual residuals: %g     %g\n", rr1.twonorm(), rr2.twonorm());
  //---
  
  ///////////////////////////////////////////////////////////////////////
  // r is the sln to the reduced system
  // the sln to the aug system should be 
  //      x = [r1; r2;  (zDiag)^{-1} * (b3-C*r1);
  ///////////////////////////////////////////////////////////////////////
  SimpleVector b1(&b[0],           locnx);
  SimpleVector b2(&b[locnx],       locmy);
  SimpleVector b3(&b[locnx+locmy], locmz);
  SimpleVector r2(&r[locnx],       locmy);
  b1.copyFrom(r1);
  b2.copyFrom(r2);
  C.mult(1.0, b3, -1.0, r1);
  b3.componentDiv(*zDiag);

  //--done
  stochNode->resMon.recDsolveTmLocal_stop();


}


DoubleLinearSolver*
QpStochLinsysRootAugRedPrPCG::createSolver(QpGenStochData* prob, SymMatrix* kktmat_)
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
      if(NULL==Qmult) Qmult = new StoredMatTimesVec(kktmat);

      StochTreePrecond* stochNodePr = dynamic_cast<StochTreePrecond*>(stochNode);
      if(NULL==Pmult) Pmult = new RemoteMatTimesVec(stochNodePr);

      SparseGenMatrix & A = prob->getLocalB();
      if(NULL==Atmult) Atmult = new StoredMatTransTimesVec(&A);

      return new PCGSolver(Qmult, Pmult, Atmult, locnx, locmy);
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


void QpStochLinsysRootAugRedPrPCG::updateKKT(QpGenStochData* prob, Variables* vars)
{
  int j, p, pend; double val;

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  DenseSymMatrix * kktd = (DenseSymMatrix*) kkt; 
  //alias for internal buffer of kkt
  double** dKkt = kktd->Mat();
  
  //////////////////////////////////////////////////////
  // compute Q+diag(xdiag) - C' * diag(zDiag) * C 
  // and update the KKT
  //////////////////////////////////////////////////////
  SparseGenMatrix& C = prob->getLocalD();
  C.matTransDinvMultMat(*zDiag, &CtDC);
  assert(CtDC->size() == locnx);

  /////////////////////////////////////////////////////////////
  // 1. update the KKT with Q (DO NOT PUT DIAG)
  /////////////////////////////////////////////////////////////
  SparseSymMatrix& Q = prob->getLocalQ();
  int* krowQ=Q.krowM(); int* jcolQ=Q.jcolM(); double* dQ=Q.M();
  for(int i=0; i<locnx; i++) {
    pend = krowQ[i+1];
    for(p=krowQ[i]; p<pend; p++) {
      j = jcolQ[p]; 
      if(i==j) continue; 
      val = dQ[p];
      dKkt[i][j] += val;
      dKkt[j][i] += val;
#ifdef DEBUG
      assert(i<j);
#endif
    }
  }
  
  //////////////////////////////////////////////////////////////////
  // 2. update the KKT with the diagonals
  //    xDiag is in fact diag(Q)+X^{-1}S, so overwrite the diagonal 
  //  instead of adding to it
  /////////////////////////////////////////////////////////////////
  SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
  for(int i=0; i<locnx; i++) dKkt[i][i] += sxDiag[i];

  //aliases for internal buffers of CtDC
  SparseSymMatrix* CtDCsp = reinterpret_cast<SparseSymMatrix*>(CtDC);
  int* krowCtDC=CtDCsp->krowM(); int* jcolCtDC=CtDCsp->jcolM(); double* dCtDC=CtDCsp->M();
  /////////////////////////////////////////////////////////////
  // update the KKT with   - C' * diag(zDiag) *C
  /////////////////////////////////////////////////////////////
  for(int i=0; i<locnx; i++) {
    pend = krowCtDC[i+1];
    for(p=krowCtDC[i]; p<pend; p++) {
      j = jcolCtDC[p];
      dKkt[i][j] -= dCtDC[p];
    }
  }

  
  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////

  //only preconditioner does this since the workers does not store A in the kkt matrix
  if(me==ePrecond)
    kktd->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx, 1);

  /////////////////////////////////////////////////////////////
  // update the KKT zeros for the lower right block 
  /////////////////////////////////////////////////////////////
  //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}

QpStochLinsysRootAugRedPrPCG::NodeType QpStochLinsysRootAugRedPrPCG::whoAmI()
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



#define CHUNK_SIZE 1024*1024*16 //(128Mb = 16M of doubles)
void 
QpStochLinsysRootAugRedPrPCG::workerReduce_source(DenseSymMatrix* M, 
						  int target, 
						  MPI_Comm comm)
{
  int n;
#ifdef DEBUG
  int m; M->getSize(m,n);
  assert(m==locnx);
#endif
  n = locnx*locnx;
  int chunksize = min(n, CHUNK_SIZE/n * n);
  int rows_in_chunk = chunksize/locnx;
  //double* buffer = new double[chunksize];
  int irow=0;

  while(true) {
    MPI_Reduce(&(M->mStorage->M[irow][0]), NULL, chunksize, 
	       MPI_DOUBLE, MPI_SUM, target, comm);
    irow += rows_in_chunk;
    if(irow>=locnx) break;
  }
}

void 
QpStochLinsysRootAugRedPrPCG::workerReduce_target(DenseSymMatrix* M, 
						  int target, 
						  MPI_Comm comm)
{
  int n;
#ifdef DEBUG
  int m; M->getSize(m,n);
  assert(m==locnx);
#endif
  n = locnx*locnx;
  int chunksize = min(n, CHUNK_SIZE/n * n);
  int rows_in_chunk = chunksize/locnx;
  double* buffer = new double[chunksize];
  int irow=0;

  while(true) {
    MPI_Reduce(&(M->mStorage->M[irow][0]), buffer, chunksize, 
	       MPI_DOUBLE, MPI_SUM, target, comm);
    memcpy(&(M->mStorage->M[irow][0]), buffer, chunksize*sizeof(double));

    irow += rows_in_chunk;
    if(irow>=locnx) break;
  }
  delete[] buffer;
}

void 
QpStochLinsysRootAugRedPrPCG::precndReduce_source(DenseSymMatrix* M, 
						  int target, 
						  MPI_Comm comm)
{
  int n;
#ifdef DEBUG
  int m; M->getSize(m,n);
  assert(m>=locnx);
#endif
  n = locnx*locnx;
  int chunksize = min(n, CHUNK_SIZE/n * n);
  int rows_in_chunk = chunksize/locnx;
  double* buffer = new double[chunksize];

  int irow=0;
  while(true) {
    for(int i=irow; i<irow+rows_in_chunk; i++)
      memcpy(&buffer[(i-irow)*locnx], &(M->mStorage->M[i][0]), locnx*sizeof(double));
    
    MPI_Reduce(buffer, NULL, chunksize, 
	       MPI_DOUBLE, MPI_SUM, target, comm);
      
    irow += rows_in_chunk;
    if(irow>=locnx) break;
  }

  delete[] buffer;
}
void 
QpStochLinsysRootAugRedPrPCG::precndReduce_target(DenseSymMatrix* M, 
						  int target, 
						  MPI_Comm comm)
{
  int n;
#ifdef DEBUG
  int m; M->getSize(m,n);
  assert(m>=locnx);
#endif
  n = locnx*locnx;
  int chunksize = min(n, CHUNK_SIZE/n * n);
  int rows_in_chunk = chunksize/locnx;
  double* buffer = new double[chunksize];
  double* result = new double[chunksize];
  
  int irow=0, i;
  while(true) {
    for(i=irow; i<irow+rows_in_chunk; i++)
      memcpy(&buffer[(i-irow)*locnx], &(M->mStorage->M[i][0]), locnx*sizeof(double));

    MPI_Reduce(buffer, result, chunksize, 
	       MPI_DOUBLE, MPI_SUM, target, comm);

    for(int i=irow; i<irow+rows_in_chunk; i++)
      memcpy(&(M->mStorage->M[i][0]), &result[(i-irow)*locnx],  locnx*sizeof(double));

    irow += rows_in_chunk;
    if(irow>=locnx) break;
  }
  delete[] buffer; delete[] result;
}
