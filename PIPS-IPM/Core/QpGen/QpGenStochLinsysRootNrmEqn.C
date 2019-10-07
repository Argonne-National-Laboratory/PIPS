#include "QpGenStochLinsysRootNrmEqn.h"
#include "QpGenStoch.h"
#include "DeSymIndefSolver.h"
#include "DeSymPSDSolver.h"
#include "QpGenStochData.h"

QpGenStochLinsysRootNrmEqn::QpGenStochLinsysRootNrmEqn(QpGenStoch * factory_, 
						       QpGenStochData * prob)
  : QpGenStochLinsysRootAugRed(),
    AQinvAt(NULL), solver2(NULL)
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

  //QpGenStochLinsyRootAugRed
  CtDC=NULL;
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);

  createChildren(prob);
};

QpGenStochLinsysRootNrmEqn::QpGenStochLinsysRootNrmEqn(QpGenStoch* factory_,
					   QpGenStochData* prob,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : QpGenStochLinsysRootAugRed(), 
    AQinvAt(NULL), solver2(NULL)
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

    //QpGenStochLinsyRootAugRed
  CtDC=NULL;
  kkt = createKKT(prob);
  solver = createSolver(prob, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);

  createChildren(prob);
};

QpGenStochLinsysRootNrmEqn::~QpGenStochLinsysRootNrmEqn()
{
  if(AQinvAt) delete AQinvAt;
  if(solver2) delete solver2;
  delete redRhs;
}


SymMatrix* 
QpGenStochLinsysRootNrmEqn::createKKT(QpGenStochData* prob)
{
  int n = locnx;
  return new DenseSymMatrix(n);
}

void QpGenStochLinsysRootNrmEqn::updateKKT(QpGenStochData* prob, Variables* vars)
{
  int j, p, pend; double val;

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  DenseSymMatrix * kktd = (DenseSymMatrix*) kkt; 
  myAtPutZeros(kktd, 0, 0, locnx, locnx);


  //////////////////////////////////////////////////////
  // compute Q+diag(xdiag) - C' * diag(zDiag) * C 
  // and update the KKT
  //////////////////////////////////////////////////////
  SparseGenMatrix& C = prob->getLocalD();
  C.matTransDinvMultMat(*zDiag, &CtDC);
  assert(CtDC->size() == locnx);


  //aliases for internal buffers of CtDC
  SparseSymMatrix* CtDCsp = reinterpret_cast<SparseSymMatrix*>(CtDC);
  int* krowCtDC=CtDCsp->krowM(); int* jcolCtDC=CtDCsp->jcolM(); double* dCtDC=CtDCsp->M();

  //alias for internal buffer of kkt
  double** dKkt = kktd->Mat();

  /////////////////////////////////////////////////////////////
  // update the KKT with Q
  /////////////////////////////////////////////////////////////
  SparseSymMatrix& Q = prob->getLocalQ();
  int* krowQ=Q.krowM(); int* jcolQ=Q.jcolM(); double* dQ=Q.M();
  for(int i=0; i<locnx; i++) {
    pend = krowQ[i+1];
    for(p=krowQ[i]; p<pend; p++) {
      j = jcolQ[p];
      val = dQ[p];
      dKkt[i][j] = val;
      dKkt[j][i] = val;
      assert(i<=j);
    }
  }
  /////////////////////////////////////////////////////////////
  // update the KKT with the diagonals
  // xDiag is in fact diag(Q)+X^{-1}S, so overwrite the diagonal 
  // instead of adding to it
  /////////////////////////////////////////////////////////////
  kktd->atPutDiagonal( 0, *xDiag );

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

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}

DoubleLinearSolver*
QpGenStochLinsysRootNrmEqn::createSolver(QpGenStochData* prob, SymMatrix* kktmat_)
{

  DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
  return new DeSymPSDSolver(kktmat);
}


void QpGenStochLinsysRootNrmEqn::factor2(QpGenStochData *prob, Variables *vars)
{
  assert( children.size() == prob->children.size() );

  // Diagonals are already updated but KKT was not formed/updated
  this->updateKKT(prob, vars);

  //dumpMatrix(-1, 0, "M", *kktd);  

  // First tell children to factorize.
  for(int it=0; it<children.size(); it++) {
    children[it]->factor2(prob->children[it], vars);
  }
  
  //!parallel communication here
  
  // each child computes   Gi Li^{-T} Di^{-1} Li^{-1} Gi^T  as
  // BLAS3 mat-mat product of U^T V, where
  //          - U=Li^{-1} Gi^T
  //          - V=Di^{-1} U

  DenseGenMatrix* U = NULL;
  DenseGenMatrix* V = NULL;
  

  stochNode->resMon.recSchurMultLocal_start();

  allocUtV();
  initializeUtV();

  stochNode->resMon.recSchurMultLocal_stop();

  for(int it=0; it<children.size(); it++) {

    // Allocate  U if necessary. If the sizes of 
    // the previous child are good, reuse the mem.

    if(children[it]->mpiComm == MPI_COMM_NULL)
      continue;
    children[it]->stochNode->resMon.recFactTmChildren_start();    

    children[it]->allocU(&U, locnx);
    children[it]->allocV(&V, locnx);

    children[it]->computeU_V(prob->children[it], U, V);

    children[it]->stochNode->resMon.recSchurMultChildren_start();

    //!dumping data
    //DenseSymMatrix* UtV0 = new DenseSymMatrix(locnx);
    //UtV0->scalarMult(0.0);
    //UtV0->matMult(-1.0, *U, 1, *V, 0, 1.0); 
    //dumpMatrix(it, 0, "C", *UtV0);
    //delete UtV0;
    //~dumping data

    UtV->matMult(-1.0, *U, 1, *V, 0, 1.0); 
    children[it]->stochNode->resMon.recSchurMultChildren_stop();

    children[it]->stochNode->resMon.recFactTmChildren_stop();

    //printf("child %d done\n", it);
  }

  stochNode->resMon.recSchurMultLocal_start();
  delete U; delete V;
  //deleteUtV(); reuse this
  stochNode->resMon.recSchurMultLocal_stop();

  //parallel communication
  if(iAmDistrib) {

    double* buffer=new double[locnx*locnx];
    MPI_Allreduce(&(UtV->mStorage->M[0][0]),
		  buffer,
		  locnx*locnx,
		  MPI_DOUBLE,
		  MPI_SUM,
		  MPI_COMM_WORLD);

    memcpy(&UtV->mStorage->M[0][0], buffer, locnx*locnx*sizeof(double));

    delete[] buffer;
  }

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  //update the upper block of the kkt with the UtV block
  kkt->symAtPutSubmatrix(0, 0, *UtV, 0, 0, locnx, locnx);

  stochNode->resMon.recSchurMultLocal_stop();  

  //double st=MPI_Wtime(); 
  solver->matrixChanged(); 
  //printf("fact1 took %g sec\n", MPI_Wtime()-st);

  /////////////////////////////////////////////////////////////////////////
  // NORMAL EQUATIONS stuff
  ////////////////////////////////////////////////////////////////////////
  //--------------------------------------------------------
  // compute A * inverse(Q) * A^T as 
  //  1) Q = LL' so AQinvAt=(L\A')' (L\A')
  //  2) W = L\A'    W is locnx-locmy matrix
  //  3) AQinvAt = W'*W (as a rank-k update)
  //--------------------------------------------------------
  if(locmy) {
    //st=MPI_Wtime();
    
    SparseGenMatrix& A = prob->getLocalB();
    DenseGenMatrix W(&(UtV->mStorage->M[0][0]), locnx, locmy); //reuse UtV
    assert(locmy<=locnx);
    myDenseFromSparseTrans(W,A);
    //printf("Dense to sparse done  %g sec\n", MPI_Wtime()-st);st=MPI_Wtime();

    DeSymPSDSolver& cholSolver = dynamic_cast<DeSymPSDSolver&>(*solver);
    cholSolver.Lsolve(W);
    //printf("W computed  %g sec\n"   , MPI_Wtime()-st);st=MPI_Wtime(); 

    //rank-k update
    if(AQinvAt==NULL) AQinvAt = new DenseSymMatrix(locmy);
    AQinvAt->atRankkUpdate(0.0, 1.0, W, 1);
    //printf("rank-k update done  %g sec\\n" , MPI_Wtime()-st);st=MPI_Wtime();   

    if(solver2==NULL) solver2 = new DeSymPSDSolver(AQinvAt);
    solver2->matrixChanged();
    
    //printf("fact2 took %g\n", MPI_Wtime()-st);st=MPI_Wtime();
  }
  stochNode->resMon.recFactTmLocal_stop();  
}


void QpGenStochLinsysRootNrmEqn::solveReduced( QpGenStochData *prob, SimpleVector& b)
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
  {
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
  // r contains all the stuff for aug-red 
  ///////////////////////////////////////////////////////////////////////
  }

  {
  SimpleVector r1(&r[0], locnx);     //x-part resid 
  //perform further reduction for normal equations
  if(locmy) {
    SimpleVector r2(&r[locnx], locmy); //y-part resid of aug-red resid

    // solve y = (A*Qinv*A) \ (A*Qinv*r1 - r2)
    ///SimpleVector& tmpVec = (*redRhs2);
    SimpleVector tmpVec(&b[0], locnx); //use supplied b as temp buffer
    tmpVec.copyFrom(r1);
    solver->Dsolve(tmpVec); // Qinv*r1

    A.mult(-1.0, r2, 1.0, tmpVec);     //r2=A*Qinv*r1 - r2
    solver2->solve(r2);

    SimpleVector& y = r2; 
    // solve x = Qinv \ (r1-A^T y)
    A.transMult(1.0, r1, -1.0, y);
  }

  solver->Dsolve(r1);
  }
  //SimpleVector x(&r[0],           locnx);
  //x.writeToStream(cout);
  //SimpleVector y(&r[locnx],       locmy);
  //y.writeToStream(cout);


  ///////////////////////////////////////////////////////////////////////
  // r is the sln to the reduced system
  // the sln to the aug system should be 
  //      x = [r1; r2;  (zDiag)^{-1} * (b3-C*r1);
  ///////////////////////////////////////////////////////////////////////
  SimpleVector r3(&r[locnx+locmy], locmz); //r3 is used as a temp
                                           //buffer for b3
  SimpleVector r1(&r[0],           locnx);
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

// D=S'; S is m-ny-n   D is n-by-m
void QpGenStochLinsysRootNrmEqn::
myDenseFromSparseTrans(DenseGenMatrix& D, SparseGenMatrix& S)
{
  int m, n; D.getSize(n,m);
  assert(S.getStorageRef().m==m);   assert(S.getStorageRef().n==n);
  int* krowM = S.krowM(); int* jcolM=S.jcolM(); double *M=S.M();

  int i, j, k, jcurrent;
  
  for ( i=0; i<m; i++ ) {
    // Loop over all rows in range
    jcurrent = -1;
    for( k=krowM[i]; k<krowM[i+1]; k++ ) {
      // Loop over the elements of the sparse row
      j = jcolM[k];
      assert(j<n);//!
      for ( jcurrent++; jcurrent < j; jcurrent++ ) {
	D[jcurrent][i] = 0.0;	//A[(i - row) * lda + jcurrent - col] = 0.0;
      }
      jcurrent = j;

      D[j][i]=M[k];    //A[(i - row) * lda + j - col] = M[k];

    } // End loop over element of the sparse row

    for( jcurrent++; jcurrent < n; jcurrent++ ) {
      D[jcurrent][i] = 0.0;  //A[(i - row) * lda + jcurrent - col] = 0.0;
    }
  }
}
