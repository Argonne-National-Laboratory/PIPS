/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeaf.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "Ma57Solver.h"
#include "Ma27Solver.h"
#include "PardisoSolver.h"
#include "DeSymIndefSolver.h"

static void mySymAtPutSubmatrix(SymMatrix& kkt, 
			 GenMatrix& B, GenMatrix& D, 
			 int locnx, int locmy, int locmz);

sLinsysLeaf::sLinsysLeaf(sFactory *factory_, sData* prob,
					   OoqpVector* dd_, 
					   OoqpVector* dq_,
					   OoqpVector* nomegaInv_,
					   OoqpVector* rhs_)
  : sLinsys(factory_, prob, dd_, dq_, nomegaInv_, rhs_)
    //,    kkt(NULL), solver(NULL)
{
  int rank; MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //double t = MPI_Wtime();
  
  // create the KKT system matrix
  // size = ?  nnz = ?
  int nnzQ, nnzB, nnzD;

  prob->getLocalSizes(locnx, locmy, locmz);
  int n = locnx+locmy+locmz;

  prob->getLocalNnz(nnzQ, nnzB, nnzD);

  //alocate the matrix and copy the data into
  SparseSymMatrix* kktsp = new SparseSymMatrix(n, n+nnzQ+nnzB+nnzD);
  kkt = kktsp;

  SimpleVectorHandle v( new SimpleVector(n) );
  v->setToZero();
  kkt->setToDiagonal(*v);

  kkt->symAtPutSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);
  if(locmz>0) {
    kkt->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);
    kkt->symAtPutSubmatrix( locnx+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
    
  } else
    mySymAtPutSubmatrix(*kkt, prob->getLocalB(), prob->getLocalD(), locnx, locmy, locmz);

  // create the solver for the linear system
  solver = new Ma57Solver(kktsp);
  //solver = new PardisoSolver(kktsp);
  //solver = new DeSymIndefSolver(kktsp);

  //t = MPI_Wtime() - t;
  //if (rank == 0) printf("new sLinsysLeaf took %f sec\n",t);

  mpiComm = (dynamic_cast<StochVector*>(dd_))->mpiComm;
}

sLinsysLeaf::~sLinsysLeaf()
{

}

void sLinsysLeaf::factor2(sData *prob, Variables *vars)
{
  // Diagonals were already updated, so
  // just trigger a local refactorization (if needed, depends on the type of lin solver).
    
  stochNode->resMon.recFactTmLocal_start();
  solver->matrixChanged();
  stochNode->resMon.recFactTmLocal_stop();
}

void sLinsysLeaf::putXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  kkt->atPutDiagonal( 0, *xdiag.vec );
}

void sLinsysLeaf::putZDiagonal( OoqpVector& zdiag_)
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  kkt->atPutDiagonal( locnx+locmy, *zdiag.vec );
}

void sLinsysLeaf::Lsolve  (  sData *prob, OoqpVector& x_in )
{
  StochVector& x = dynamic_cast<StochVector&>(x_in);
  assert(x.children.size()==0);
 
  stochNode->resMon.recLsolveTmChildren_start();
  solver->Lsolve(*x.vec);
  stochNode->resMon.recLsolveTmChildren_stop();

}

void sLinsysLeaf::Dsolve  (  sData *prob, OoqpVector& x_in )
{
  StochVector& x = dynamic_cast<StochVector&>(x_in);
  assert(x.children.size()==0);
  stochNode->resMon.recDsolveTmChildren_start();
  solver->Dsolve(*x.vec);
  stochNode->resMon.recDsolveTmChildren_stop();
}

void sLinsysLeaf::Ltsolve (  sData *prob, OoqpVector& x_in )
{
  StochVector& x = dynamic_cast<StochVector&>(x_in);
  assert(x.children.size()==0);
  stochNode->resMon.recLtsolveTmChildren_start();
  solver->Ltsolve(*x.vec);
  stochNode->resMon.recLtsolveTmChildren_stop();
}

void sLinsysLeaf::Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp)
{

  StochVector& b   = dynamic_cast<StochVector&>(x);
  SimpleVector& bi = dynamic_cast<SimpleVector&>(*b.vec);

  assert(0==b.children.size());
  stochNode->resMon.recLtsolveTmChildren_start();
  //b_i -= Lni^T x0
  this->LniTransMult(prob, bi, -1.0, xp);
  solver->Ltsolve(bi);
  stochNode->resMon.recLtsolveTmChildren_stop();
}

void sLinsysLeaf::sync()
{ assert(false); }


void sLinsysLeaf::deleteChildren()
{ }

static void mySymAtPutSubmatrix(SymMatrix& kkt_, 
			 GenMatrix& B_, GenMatrix& D_, 
			 int locnx, int locmy, int locmz)
{
  SparseSymMatrix& kkt = reinterpret_cast<SparseSymMatrix&>(kkt_);
  SparseGenMatrix& B   = reinterpret_cast<SparseGenMatrix&>(B_);
  //SparseGenMatrix& D   = reinterpret_cast<SparseGenMatrix&>(D_);

  int* jcolK = kkt.jcolM(); int* jcolB = B.jcolM(); //int* jcolD = D.jcolM(); 
  int* krowK = kkt.krowM(); int* krowB = B.krowM(); //int* krowD =  D.krowM();
  double* MK = kkt.M();     double* MB = B.M();

  for(int i=0; i<locmy; i++) {
    int itK = krowK[i+locnx];
    int j = krowB[i];

    for(; j<krowB[i+1]; j++) { 

      if(jcolB[j]<i+locnx) {
	jcolK[itK]=jcolB[j]; 
	MK[itK]=MB[j]; 
	itK++;
      }
    }
    jcolK[itK]=i+locnx; MK[itK] = 0.0; itK++;

    assert(j==krowB[i+1]);
    //for(; j<krowB[i+1]; j++) { jcolK[itK]=jcolB[j]; MK[itK]=MB[j]; itK++; assert(false);}

    krowK[i+locnx+1]=itK;
  }


  assert(locmz==0);
  //assert(false);
  //12120
  //for(int i=0; i<-240; i++) {
  //  for(int j=krowK[i]; j<krowK[i+1]; j++)
  //    printf("%5d %5d    %8.3f\n", i, jcolK[j], MK[j]);
  //}
}


