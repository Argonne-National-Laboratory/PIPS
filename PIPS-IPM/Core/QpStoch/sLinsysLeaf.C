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

void sLinsysLeaf::mySymAtPutSubmatrix(SymMatrix& kkt_, 
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

    krowK[i+locnx+1]=itK;
  }
  assert(locmz==0);
}


