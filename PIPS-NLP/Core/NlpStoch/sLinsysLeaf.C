/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "sLinsysLeaf.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

extern int gOuterSolve;
extern int separateHandDiag;

extern int gisNLP;

sLinsysLeaf::~sLinsysLeaf()
{

}

int sLinsysLeaf::factor2(sData *prob, Variables *vars)
{  
  int negEValTemp=0;
  // Diagonals were already updated, so
  // just trigger a local refactorization (if needed, depends on the type of lin solver).
  stochNode->resMon.recFactTmLocal_start();
  negEValTemp = solver->matrixChanged();
  stochNode->resMon.recFactTmLocal_stop();
  return negEValTemp;
}

void sLinsysLeaf::putXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  kkt->atPutDiagonal( 0, *xdiag.vec );
}

void sLinsysLeaf::putSDiagonal( OoqpVector& sdiag_ )
{
  StochVector& sdiag = dynamic_cast<StochVector&>(sdiag_);
  kkt->atPutDiagonal( locnx, *sdiag.vec );
}

void sLinsysLeaf::putYDualDiagonal( OoqpVector& ydiag_ )
{
  StochVector& ydiag = dynamic_cast<StochVector&>(ydiag_);

  if(gOuterSolve < 3){
  	kkt->atPutDiagonal( locnx, *ydiag.vec );
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
	kkt->atPutDiagonal( locnx + locmz, *ydiag.vec );
  }else{
	assert(0);
  } 
}

void sLinsysLeaf::putZDiagonal( OoqpVector& zdiag_)
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  
  if(gOuterSolve < 3){
  	kkt->atPutDiagonal( locnx+locmy, *zdiag.vec );
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
    kkt->atPutDiagonal( locnx+locmz+locmy, *zdiag.vec );
  }else{
	assert(0);
  }  
}

void sLinsysLeaf::setXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  kkt->copyDiagonalVal_From( 0, *xdiag.vec, firstXDiagUpdate, xDiagIdxMap);
  firstXDiagUpdate = false;  
}

void sLinsysLeaf::setSDiagonal( OoqpVector& sdiag_ )
{
  StochVector& sdiag = dynamic_cast<StochVector&>(sdiag_);
  kkt->copyDiagonalVal_From( locnx, *sdiag.vec, firstSDiagUpdate, sDiagIdxMap);
  firstSDiagUpdate = false;  
}

void sLinsysLeaf::setYDiagonal( OoqpVector& ydiag_ )
{
  StochVector& ydiag = dynamic_cast<StochVector&>(ydiag_);

  if(gOuterSolve < 3){
  	kkt->copyDiagonalVal_From( locnx, *ydiag.vec, firstYDiagUpdate, yDiagIdxMap);
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
	kkt->copyDiagonalVal_From( locnx + locmz, *ydiag.vec, firstYDiagUpdate, yDiagIdxMap );
  }else{
	assert(0);
  } 
  firstYDiagUpdate = false;
}

void sLinsysLeaf::setZDiagonal( OoqpVector& zdiag_)
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  
  if(gOuterSolve < 3){
  	kkt->copyDiagonalVal_From( locnx+locmy, *zdiag.vec, firstZDiagUpdate, zDiagIdxMap );
  }else if(gOuterSolve >= 3 && separateHandDiag==0){
    kkt->copyDiagonalVal_From( locnx+locmz+locmy, *zdiag.vec, firstZDiagUpdate, zDiagIdxMap );
  }else{
	assert(0);
  }  
  firstZDiagUpdate = false;
}

void sLinsysLeaf::setAdditiveDiagonal()
{
  assert(gOuterSolve >= 3 && separateHandDiag==1);
  StochVector& additiveDiag_ = dynamic_cast<StochVector&>(*additiveDiag);
  kkt->setAdditiveDiagonal(*additiveDiag_.vec);
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

#ifdef TIMING
  stochNode->resMon.eLtsolve.clear();
  stochNode->resMon.recLtsolveTmLocal_start();
#endif

  //b_i -= Lni^T x0
  this->LniTransMult(prob, bi, -1.0, xp);
//  solver->Ltsolve(bi);
#ifdef TIMING
  stochNode->resMon.recLtsolveTmChildren_stop();
#endif
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
  SparseGenMatrix& D   = reinterpret_cast<SparseGenMatrix&>(D_);

  int* jcolK = kkt.jcolM(); 
  int* krowK = kkt.krowM();
  double* MK = kkt.M();

  int* jcolB = B.jcolM(); 
  int* krowB = B.krowM(); 
  double* MB = B.M();

  int* jcolD = D.jcolM();    
  int* krowD =  D.krowM();
  double* MD = D.M();

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

  for(int i=0; i<locmz; i++) {
    int itK = krowK[i+locnx+locmy];
    int j = krowD[i];

    for(; j<krowD[i+1]; j++) { 

      if(jcolD[j]<i+locnx) {
	jcolK[itK]=jcolD[j]; 
	MK[itK]=MD[j]; 
	itK++;
      }
    }
    jcolK[itK]=i+locnx+locmy; MK[itK] = 0.0; itK++;

    assert(j==krowD[i+1]);

    krowK[i+locnx+locmy+1]=itK;
  }

  
//  assert(locmz==0);
}


void
sLinsysLeaf::UpdateMatrices( Data * prob_in,int const updateLevel)
{
  if(!gisNLP) return;

  sData* prob = dynamic_cast<sData*>(prob_in);

  if(updateLevel>=2){
   	if(gOuterSolve < 3){	
      kkt->symAtSetSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx,firstQUpdate, LocQMap);
      kkt->symAtSetSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx,firstBUpdate, LocBMap);
      kkt->symAtSetSubmatrix( locnx+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx,firstDUpdate, LocDMap);	
	}else{
	  kkt->symAtSetSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx,firstQUpdate, LocQMap);
      kkt->symAtSetSubmatrix( locnx + locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx,firstBUpdate, LocBMap);
      kkt->symAtSetSubmatrix( locnx + locmz + locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx,firstDUpdate, LocDMap);
	}
	firstQUpdate = false;
	firstBUpdate = false;
	firstDUpdate = false;
  }

}



