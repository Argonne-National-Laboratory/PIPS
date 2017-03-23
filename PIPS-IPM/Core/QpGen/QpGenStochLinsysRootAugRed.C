#include "QpGenStochLinsysRootAugRed.h"
#include "DeSymIndefSolver.h"
#include "DeSymPSDSolver.h"
#include "QpGenStochData.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

double randN(double a, double b);

QpGenStochLinsysRootAugRed::QpGenStochLinsysRootAugRed(QpGenStoch * factory_, 
						       QpGenStochData * prob_)
  : QpGenStochLinsysRootAug(factory_, prob_), CtDC(NULL)
{ 
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
};

QpGenStochLinsysRootAugRed::QpGenStochLinsysRootAugRed(QpGenStoch* factory_,
					   QpGenStochData* prob_,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : QpGenStochLinsysRootAug(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), CtDC(NULL)
{ 
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
};

QpGenStochLinsysRootAugRed::~QpGenStochLinsysRootAugRed()
{
  if(CtDC) delete CtDC;
  delete redRhs;
}


SymMatrix* 
QpGenStochLinsysRootAugRed::createKKT(QpGenStochData* prob)
{
  int n = locnx+locmy;
  return new DenseSymMatrix(n);
}


DoubleLinearSolver*
QpGenStochLinsysRootAugRed::createSolver(QpGenStochData* prob, SymMatrix* kktmat_)
{

  DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
  return new DeSymIndefSolver(kktmat);
  //return new DeSymPSDSolver(kktmat);
}

void QpGenStochLinsysRootAugRed::factor2(QpGenStochData *prob, Variables *vars)
{
  assert( children.size() == prob->children.size() );

  //!!
  DenseSymMatrix * kktd = (DenseSymMatrix*) kkt; 
  myAtPutZeros(kktd, 0, 0, locnx+locmy, locnx+locmy);
  //~~

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


    kktd->matMult(-1.0, *U, 1, *V, 0, 1.0);


    children[it]->stochNode->resMon.recSchurMultChildren_stop();

    children[it]->stochNode->resMon.recFactTmChildren_stop();
  }

  stochNode->resMon.recSchurMultLocal_start();
  delete U; delete V;
  stochNode->resMon.recSchurMultLocal_stop();

  //parallel communication
  if(iAmDistrib) submatrixAllReduce(kktd, 0, 0, locnx, locnx, mpiComm);

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  //update the KKT matrix with C'*D*C+Q and put A and A' blocks
  updateKKT(prob, vars);

  stochNode->resMon.recSchurMultLocal_stop();  

  double st=MPI_Wtime();
  solver->matrixChanged();
  printf("fact took %g\n", MPI_Wtime()-st);

  stochNode->resMon.recFactTmLocal_stop();  
}


void QpGenStochLinsysRootAugRed::Dsolve( QpGenStochData *prob, OoqpVector& x )
{
  StochVector& bst = dynamic_cast<StochVector&>(x);

  for(int it=0; it<children.size(); it++) {
    children[it]->Dsolve(prob->children[it], *bst.children[it]);
  }

  SimpleVector& b = dynamic_cast<SimpleVector&>(*bst.vec); 
  solveReduced(prob, b);
}

void QpGenStochLinsysRootAugRed::solveReduced( QpGenStochData *prob, SimpleVector& b)
{
  assert(locnx+locmy+locmz==b.length());
  SimpleVector& r = (*redRhs);
  assert(r.length() <= b.length());
  SparseGenMatrix& C = prob->getLocalD();

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
  //---
  //SimpleVector rr1(locnx); rr1.copyFrom(r1);
  //SimpleVector rr2(locmy); rr2.copyFromArray(&b[locnx]);
  //SimpleVector rr(locnx+locmy); rr.copyFromArray(&b[0]);
  //---
  SimpleVector realRhs(&r[0], locnx+locmy);

  //cout << "AUG: xin rhs\n"; realRhs.writeToStream(cout);
  //SimpleVector yy(&realRhs[locnx], locmy); cout << "AUG: yin rhs\n"; yy.writeToStream(cout);

  solver->Dsolve(realRhs);

  //cout << "AUG:xout rhs\n"; realRhs.writeToStream(cout);
  //cout << "AUG:yout rhs\n"; yy.writeToStream(cout);

  //---
  //SparseGenMatrix& A = prob->getLocalB();
  //SimpleVector y(&realRhs[locnx], locmy);
  //SimpleVector x(&realRhs[0], locnx);

  //A.transMult(1.0, rr1, -1.0, y);
  //K_->mult(1.0, rr1, -1.0, x);double** Kmat = K_->Mat();
  //K_->mult(1.0, rr, -1.0, realRhs);
  //A.mult(1.0, rr2, -1.0, x);
  //printf("Actual residuals: %g    %g\n", rr.twonorm(), rr2.twonorm());

  /*for(int i=7; i<17; i++) {
    for(int j=7; j<17; j++) 
      printf("%18.14f ", Kmat[i][j]);
    printf("\n");
  }
  */
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

void QpGenStochLinsysRootAugRed::updateKKT(QpGenStochData* prob, Variables* vars)
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
 

  /////////////////////////////////////////////////////////////
  // update the KKT with Q (DO NOT PUT DIAG)
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
  
  /////////////////////////////////////////////////////////////
  // update the KKT with the diagonals
  // xDiag is in fact diag(Q)+X^{-1}S
  /////////////////////////////////////////////////////////////
  //kktd->atPutDiagonal( 0, *xDiag );
  SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
  for(int i=0; i<locnx; i++) dKkt[i][i] += sxDiag[i];


  /////////////////////////////////////////////////////////////
  // update the KKT with   - C' * diag(zDiag) *C
  /////////////////////////////////////////////////////////////
  SparseGenMatrix& C = prob->getLocalD();
  C.matTransDinvMultMat(*zDiag, &CtDC);
  assert(CtDC->size() == locnx);

  //aliases for internal buffers of CtDC
  SparseSymMatrix* CtDCsp = reinterpret_cast<SparseSymMatrix*>(CtDC);
  int* krowCtDC=CtDCsp->krowM(); 
  int* jcolCtDC=CtDCsp->jcolM(); 
  double* dCtDC=CtDCsp->M();

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
  kktd->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx, 1);

  /////////////////////////////////////////////////////////////
  // update the KKT zeros for the lower right block 
  /////////////////////////////////////////////////////////////
  //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();

  //~~~~~~~~~~~~~~~~~~~~~~~~` loggging
  //dumpMatrix(-1, 0, "M", *kktd);

  //for(int i=0; i<locnx+locmy; i++)
  //  for(int j=0; j<locnx+locmy; j++) {
  //    dKkt[i][j] += M[i][j];
  //  }
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
/*
void QpGenStochLinsysRootAugRed::updateKKT(QpGenStochData* prob, Variables* vars)
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

  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////
  kktd->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx, 1);

  /////////////////////////////////////////////////////////////
  // update the KKT zeros for the lower right block 
  /////////////////////////////////////////////////////////////
  //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  myAtPutZeros(kktd, locnx, locnx, locmy, locmy);

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}*/

//typedef void (MPI_User_function) ( void * a, 
//               void * b, int * len, MPI_Datatype * ); 

//static int iter;

//void summation(void* a, void* b, int * len, MPI_Datatype* type)
//{
//  printf("len=%d\n", *len);;fflush(stdout);
//  double* da=(double*)a; double* db = (double*)b;
//  for(iter=0;iter<*len; iter++) db[iter] += da[iter]; 
//  printf("done\n");fflush(stdout);
//}

#define CHUNK_SIZE 1024*1024*16 //doubles  = 128 MBytes (maximum)
//todo DR: this function seems to be wrong, should be replaced (as in sLinsysRoot.C); is it ever called?
void QpGenStochLinsysRootAugRed::submatrixAllReduce(DenseSymMatrix* A, 
						    int row, int col, int drow, int dcol,
						    MPI_Comm comm)
{
  double ** M = A->mStorage->M;
  int n = A->mStorage->n;
#ifdef DEBUG 
  assert(n >= row+drow);
  assert(n >= col+dcol);
#endif
  int iErr;
  int chunk_size = CHUNK_SIZE / n * n; 
  chunk_size = min(chunk_size, n*n);
  double* chunk = new double[chunk_size];

  int rows_in_chunk = chunk_size/n;
  int iRow=row;
  do {

    if(iRow+rows_in_chunk > drow)
      rows_in_chunk = drow-iRow;

    //printf("irow=%d rows_in_chunk=%d rows=%d n=%d\n", iRow, rows_in_chunk, drow, n);
    //the REDUCE
    //MPI_Datatype mpiVecType;
    //iErr=MPI_Type_vector(1, dcol, n, MPI_DOUBLE, &mpiVecType); assert(iErr==MPI_SUCCESS);
    //iErr=MPI_Type_commit(&mpiVecType);assert(iErr==MPI_SUCCESS);
    //MPI_Op sumOp; 
    //iErr=MPI_Op_create(&summation, 1, &sumOp);assert(iErr==MPI_SUCCESS);

    //iErr=MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk, mpiVecType, sumOp, comm);assert(iErr==MPI_SUCCESS);
    iErr=MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk*n, MPI_DOUBLE, MPI_SUM, comm);assert(iErr==MPI_SUCCESS);

    //MPI_Op_free(&sumOp);
    //MPI_Type_free(&mpiVecType);
    //printf("MPI_REDUCE done\n");

    //copy data in M
    for(int i=iRow; i<iRow+rows_in_chunk; i++) {
      int shft = (i-iRow)*n;
      for(int j=col; j<col+dcol; j++)
	M[i][j] = chunk[shft+j];
    }
    iRow += rows_in_chunk;
  
  } while(iRow<row+drow);

  delete[] chunk;

}
