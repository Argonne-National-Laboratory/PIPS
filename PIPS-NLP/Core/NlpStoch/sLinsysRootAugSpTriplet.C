/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAugSpTriplet.h"
#include "DeSymIndefSolver.h"
#include "sData.h"
#include "sTree.h"

#include "MumpsSolver.h"

#include <unistd.h>
#include "math.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif
extern int gInnerSCsolve;
extern int gOuterSolve;
extern int separateHandDiag;

using namespace std;
extern int gBuildSchurComp;

sLinsysRootAugSpTriplet::sLinsysRootAugSpTriplet(sFactory * factory_, sData * prob_)
  : sLinsysRootAug(factory_, prob_)
{ 
  assert(gBuildSchurComp==3);
};

sLinsysRootAugSpTriplet::sLinsysRootAugSpTriplet(sFactory* factory_,
						 sData* prob_,
						 OoqpVector* dd_, 
						 OoqpVector* dq_,
						 OoqpVector* nomegaInv_,
						 OoqpVector* rhs_,
						 OoqpVector* additiveDiag_)
  : sLinsysRootAug(factory_, prob_, dd_, dq_, nomegaInv_, rhs_, additiveDiag_)
{ 
  assert(gOuterSolve>=3);  
  assert(gBuildSchurComp==3);
};

sLinsysRootAugSpTriplet::~sLinsysRootAugSpTriplet()
{

}


SymMatrix* 
sLinsysRootAugSpTriplet::createKKT(sData* prob)
{
  int n;

  if(gOuterSolve < 3){
    n = locnx+locmy;
    assert(locmz==0);
  }else{
    n = locnx+locmy+locmz+locmz;
  }
  //assert(false);
  //return new DenseSymMatrix(n);
  return new SparseSymMatrixRowMajList(n);
}

extern int gBuildSchurComp;

DoubleLinearSolver*
sLinsysRootAugSpTriplet::createSolver(sData* prob, SymMatrix* kktmat_)
{
  assert(gBuildSchurComp==3);

  //1. create the communicator for the subset of processes that MUMPS should use
  int color = MPI_UNDEFINED;
  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  //color only rank 0 and leave the other ones uncolored. MUMPS communicator created below will be MPI_COMM_NULL on
  //the nodes not colored
  if(myRank==0)
    color = 0;
  MPI_Comm mumpsComm;
  MPI_Comm_split(mpiComm, color, 0, &mumpsComm);
  
  //2. create and return the wrapper instance for Mumps
  if(0) {
    DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
    return new MumpsDenseSolver(kktmat, mumpsComm, mpiComm);
  } else {
    SparseSymMatrixRowMajList* kktmat = dynamic_cast<SparseSymMatrixRowMajList*>(kktmat_);
    return new MumpsSolver(kktmat, mumpsComm, mpiComm);
  }

}

#ifdef TIMING
extern double t_start, troot_total, taux, tchild_total, tcomm_total;
#endif


void sLinsysRootAugSpTriplet::finalizeKKT(sData* prob, Variables* vars)
{
//  assert(locmz==0||gOuterSolve<3);
  assert(gOuterSolve>=3);

  int j, p, pend; double val;

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  SparseSymMatrixRowMajList& kktm = dynamic_cast<SparseSymMatrixRowMajList&>(*kkt);
 
  //////////////////////////////////////////////////////
  // compute Q+diag(xdiag)  
  // and update the KKT
  //////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // update the KKT with Q (use diag from xDiag DIAG)
  /////////////////////////////////////////////////////////////
  SparseSymMatrix& Q = prob->getLocalQ();
  SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
  int* krowQ=Q.krowM(); int* jcolQ=Q.jcolM(); double* dQ=Q.M();
  int jstart = krowQ[0], jend;
  for(int i=0; i<locnx; i++) {
    
    jend = krowQ[i+1];
    assert(jcolQ[jstart]==i); 

    //add +1 since we exclude the diagonal entry
    kktm.atAddSpRow(i, jcolQ + jstart+1, dQ+jstart+1, jend-jstart-1);

    jstart = jend;
    // for(p=krowQ[i]; p<pend; p++) {
    //   j = jcolQ[p]; 
    //   assert(j>=i);
    //   if(i==j) {
    // 	colValSrc.push_back(ColVal(j,sxDiag[i]));
    //   } else {
    // 	colValSrc.push_back(ColVal(j,dQ[p]));
    //   }
    // }
    // kktm.atAddSpRow(i, colValSrc);
  }

  kktm.atAddDiagonal(0,sxDiag);

  
  // for(int i=0; i<locnx; i++) {
  //   pend = krowQ[i+1];
  //   for(p=krowQ[i]; p<pend; p++) {
  //     j = jcolQ[p]; 
  //     if(i==j) continue;
  //     val = dQ[p];
  //     dKkt[i][j] += val;
  //     dKkt[j][i] += val;
  //   }
  // }
  

  /////////////////////////////////////////////////////////////
  // update the KKT with the diagonals
  // xDiag is in fact diag(Q)+X^{-1}S
  /////////////////////////////////////////////////////////////

  //SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
  SimpleVector& syDiag = dynamic_cast<SimpleVector&>(*yDiag);

  //for(int i=0; i<locnx; i++) dKkt[i][i] += sxDiag[i];  

  SimpleVector& ssDiag = dynamic_cast<SimpleVector&>(*sDiag);
  SimpleVector& szDiag = dynamic_cast<SimpleVector&>(*zDiag);

  /////////////////////////////////////////////////////////////////////
  // update the KKT with  S part 
  //      and with the left -I -> left upper (locnx,locnx+locmz+locmy)
  //                           ->right lower (locnx+locmz,locnx+locmz+locmy+locmz)
  ////////////////////////////////////////////////////////////////////
  if(locmz>0) {
    int jcol[2]; double M[2];
    for(int i=locnx; i<locnx+locmz; i++) {
      jcol[0] = i; jcol[1] = i+locmz+locmy;
      M[0] = ssDiag[i-locnx]; M[1]=-1.;
      kktm.atAddSpRow(i, jcol, M, 2);
    }
    //kktm.atPutDiagonal(locnx,ssDiag);
    // for(int i=locnx; i<locnx+locmz; i++) {
    //   //dKkt[i][i] += ssDiag[i-locnx];
    // }
  }
  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////
  if(locmy>0){
    //!kktd->symAtAddSubmatrix( locnx+locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx, 1 );`
    kktm.forceSymUpdate(true);
    kktm.symAtAddSubmatrix(locnx+locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx);
    kktm.forceSymUpdate(false);
    //!for(int i=locnx+locmz; i<locnx+locmz+locmy; i++) dKkt[i][i] += syDiag[i-locnx-locmz];
    
  }
  
  int mle = prob->getmle();
  if(mle>0){
    assert(false && "not yet supported");
    //!   kktd->symAtAddSubmatrix( locnx+locmz+locmy-mle, 0, prob->getLocalE(), 0, 0, mle, locnx, 1 );
  }
  // /////////////////////////////////////////////////////////////
  // // update the KKT with C (symmetric update forced)  and dual reg
  // /////////////////////////////////////////////////////////////  
  if(locmz>0){
    kktm.forceSymUpdate(true);
    kktm.symAtAddSubmatrix( locnx+locmz+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
    kktm.forceSymUpdate(false);
    //   kktd->symAtAddSubmatrix( locnx+locmz+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx, 1 );
    kktm.atPutDiagonal(locnx+locmz+locmy, szDiag);
    // 	for(int i=0; i<locmz; i++){
    // 		dKkt[i+locnx+locmz+locmy][i+locnx] -= 1.0;
    // 		dKkt[i+locnx][i+locnx+locmz+locmy] -= 1.0;
    // 		dKkt[i+locnx+locmz+locmy][i+locnx+locmz+locmy] += szDiag[i];
  }
  // }
  int mli = prob->getmli();
  if(mli>0){
    assert(false && "not yet supported");
    //   kktd->symAtAddSubmatrix( locnx+locmz+locmy+locmz-mli, 0, prob->getLocalF(), 0, 0, mli, locnx, 1 );
  }


  // /////////////////////////////////////////////////////////////
  // // update the KKT zeros for the lower right block 
  // /////////////////////////////////////////////////////////////
  // //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  // //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}

//this variable is just reset in this file; children will default to the "safe" linear solver
extern int gLackOfAccuracy;

int sLinsysRootAugSpTriplet::factor2(sData *prob, Variables *vars)
{
  int negEVal=0, tempNegEVal=0;
  int return_NegEval=-1;  
  int matIsSingular=0,matIsSingularAllReduce;
  int mype; MPI_Comm_rank(mpiComm, &mype);
#ifdef TIMING
  double stime=MPI_Wtime();
  double stime1=MPI_Wtime();
  gprof.n_factor2++;
#endif
  	
  SparseSymMatrixRowMajList& kktd = dynamic_cast<SparseSymMatrixRowMajList&>(*kkt);
  initializeKKT(prob, vars);

  // First tell children to factorize. 
  for(size_t c=0; c<children.size(); c++) {
    tempNegEVal = children[c]->factor2(prob->children[c], vars);
	if(tempNegEVal<0){
	  matIsSingular = 1; 
	}else{
	  negEVal += tempNegEVal;
	}
  }

  for(size_t c=0; c<children.size(); c++) {
#ifdef STOCH_TESTING
    g_scenNum=c;
#endif
    if(children[c]->mpiComm == MPI_COMM_NULL)
      continue;

    children[c]->stochNode->resMon.recFactTmChildren_start();    
    //---------------------------------------------
    children[c]->addTermToDenseSchurCompl(prob->children[c], kktd);
    //---------------------------------------------
    children[c]->stochNode->resMon.recFactTmChildren_stop();
  }
#ifdef TIMING
  gprof.t_initializeKKT+=MPI_Wtime()-stime;
  stime=MPI_Wtime();
#endif

#ifdef TIMING
  MPI_Barrier(MPI_COMM_WORLD);
  stochNode->resMon.recReduceTmLocal_start();
#endif 
  reduceKKT();
#ifdef TIMING
  gprof.t_reduceKKT+=MPI_Wtime()-stime;
  stime=MPI_Wtime();
#endif

#ifdef TIMING
  stochNode->resMon.recReduceTmLocal_stop();
#endif  
  finalizeKKT(prob, vars);
#ifdef TIMING
  gprof.t_finalizeKKT+=MPI_Wtime()-stime;
  stime=MPI_Wtime();
#endif
  
  //printf("(%d, %d) --- %f\n", PROW,PCOL, kktd[PROW][PCOL]);

  MPI_Allreduce(&matIsSingular, &matIsSingularAllReduce, 1, MPI_INT, MPI_SUM, mpiComm);

  if(0==matIsSingularAllReduce){
  	// all the diag mat is nonsingular
  	MPI_Allreduce(&negEVal, &return_NegEval, 1, MPI_INT, MPI_SUM, mpiComm);
	negEVal = factorizeKKT();
#ifdef TIMING
  gprof.t_factorizeKKT+=MPI_Wtime()-stime;
  stime=MPI_Wtime();
#endif
	if(negEVal<0){ 
	  return_NegEval = -1;
	}else{
	  return_NegEval += negEVal;
	}
  }
  
#ifdef TIMING
  afterFactor();
  gprof.t_factor2_total+=MPI_Wtime()-stime1;
#endif

  return return_NegEval;

}


void sLinsysRootAugSpTriplet::initializeKKT(sData* prob, Variables* vars)
{
  SparseSymMatrixRowMajList& kktm = dynamic_cast<SparseSymMatrixRowMajList&>(*kkt);
  kktm.setToConstant(0.);
}

void sLinsysRootAugSpTriplet::reduceKKT()
{
  SparseSymMatrixRowMajList& kktm = dynamic_cast<SparseSymMatrixRowMajList&>(*kkt);
  if(iAmDistrib) {
    assert(false);
  }
}

// int sLinsysRootAugSpTriplet::factorizeKKT()
// {
//   assert(false);
//   return -1;
// }

extern int gisNLP;

void
sLinsysRootAugSpTriplet::UpdateMatrices( Data * prob_in, int const updateLevel)
{
  int useUpdate=updateLevel;
  if(!gisNLP) useUpdate=1;

  sData* prob = dynamic_cast<sData*>(prob_in);
  SparseSymMatrixRowMajList* kktm = dynamic_cast<SparseSymMatrixRowMajList*>(kkt);
  assert(kktm!=NULL);

  //kktm stores only upper triangular, but B needs to go in the lower part -> force symmetric update 
  kktm->forceSymUpdate(true);

  if(useUpdate>=2){
    if(gOuterSolve < 3){	
      kktm->symAtSetSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);
      if(locmy>0) {
	kktm->symAtSetSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);

      }
    }else{
      kktm->symAtSetSubmatrix( 0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);
      if(locmy>0)
	kktm->symAtSetSubmatrix( locnx + locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx);
    } 
  }
  kktm->forceSymUpdate(false);

  // propagate it to the subtree
  for(size_t it=0; it<children.size(); it++)
    children[it]->UpdateMatrices(prob->children[it],updateLevel);
  
  if(useUpdate>=1){
    if(gOuterSolve < 3){
      setXDiagonal( *dd );
      setYDiagonal( *temp_diagY);	  
      setZDiagonal( *nomegaInv );	  
    }else if (gOuterSolve >= 3 && separateHandDiag==1){
      joinRHSXSYZ(*additiveDiag,*dd,*temp_diagS,*temp_diagY ,*temp_diagZ);
      setAdditiveDiagonal();
    }else if(gOuterSolve >= 3 && separateHandDiag==0){
      setXDiagonal( *dd );
      setSDiagonal( *temp_diagS );
      setYDiagonal( *temp_diagY);		  
      setZDiagonal( *temp_diagZ );		  
    }else{
      assert(0);
    }
  }
}
