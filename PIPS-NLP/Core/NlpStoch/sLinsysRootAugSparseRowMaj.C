/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAugSparseRowMaj.h"
#include "DeSymIndefSolver.h"
#include "sData.h"
#include "sTree.h"

#ifdef WITH_MUMPS
#include "MumpsSolver.h"
#endif

#include <unistd.h>
#include "math.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
extern double g_scenNum;
#endif

#ifdef TIMING
#include "../../global_var.h"
#include "../PIPS-NLP/Core/Utilities/PerfMetrics.h"
#endif

extern int gInnerSCsolve;
extern int gOuterSolve;
extern int separateHandDiag;

using namespace std;
extern int gBuildSchurComp;

sLinsysRootAugSpTriplet::sLinsysRootAugSpTriplet(sFactory * factory_, sData * prob_)
  : sLinsysRootAug(factory_, prob_), iAmRank0(false)
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
  : sLinsysRootAug(factory_, prob_, dd_, dq_, nomegaInv_, rhs_, additiveDiag_),
    iAmRank0(false)
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
  return new SparseSymMatrixRowMajList(n);
}

extern int gBuildSchurComp;

DoubleLinearSolver*
sLinsysRootAugSpTriplet::createSolver(sData* prob, SymMatrix* kktmat_)
{
  assert(gBuildSchurComp==3);

  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  iAmRank0 = (myRank==0);

  // Build node communicator on each node, with all ranks on the node
  MPI_Comm nodeComm;
  MPI_Comm_split_type(mpiComm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &nodeComm);

  // We do an allreduce on each node to find the node communicator with the rank 0

  int commWithRank0 = 0;
  if(iAmRank0) {
    int one = 1;
    MPI_Allreduce(&one, &commWithRank0, 1, MPI_INT, MPI_SUM, nodeComm);
  }
  else {
    int zero = 0;
    MPI_Allreduce(&zero, &commWithRank0, 1, MPI_INT, MPI_SUM, nodeComm);

  }
  // The process where commWithRank0 = 1 are the ones that have rank 0 on their node.
  // Those processes build the MUMPS communicator
  int color = MPI_UNDEFINED;
  if(commWithRank0 == 1) {
    int nodeRank;
    MPI_Comm_rank(nodeComm, &nodeRank);
  
    //color a couple of ranks, say 0,1,2,4 and leave the other ones uncolored. 
    //MUMPS communicator created below will be MPI_COMM_NULL on the nodes not colored

    // Now we are sure we color ranks that are located on the same node
    if(nodeRank<4)
      color = 0;
    // Make sure myRank 0 is in the MUMPS communicator and does not have nodeRank > 4
    if(myRank==0)
      color = 0;
    
  }
  MPI_Comm mumpsComm;
  int mumpsSize;
  MPI_Comm_split(mpiComm, color, 0, &mumpsComm);
  if (mumpsComm != MPI_COMM_NULL)
  {
    MPI_Comm_size(mumpsComm, &mumpsSize);
    if (myRank==0) std::cout << "MUMPS communicator size: " << mumpsSize << std::endl;
  }

#ifdef WITH_MUMPS
  //2. create and return the wrapper instance for Mumps
  if(0) {
    DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
    return new MumpsDenseSolver(kktmat, mumpsComm, mpiComm);
  } else {
    SparseSymMatrixRowMajList* kktmat = dynamic_cast<SparseSymMatrixRowMajList*>(kktmat_);
    return new MumpsSolver(kktmat, mumpsComm, mpiComm);
  }
#else
  assert(0 && "PIPS was not built with MUMPS");
#endif
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
  
  //int nnz = kktm.numberOfNonZeros();
  //printf("nnz mat------------- %d\n", nnz);

  /////////////////////////////////////////////////////////////
  // update the KKT with Q (use diag from xDiag DIAG)
  /////////////////////////////////////////////////////////////
  SparseSymMatrix& Q = prob->getLocalQ();
  SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
  int* krowQ=Q.krowM(); int* jcolQ=Q.jcolM(); double* dQ=Q.M();
  int jstart = krowQ[0], jend, nelems;
  for(int i=0; i<locnx; i++) {
    
    jend = krowQ[i+1];

    //assert(jend-jstart==0 || jcolQ[jend-1]==i); //safe to remove for the small test examples

    //we only add "jend-jstart-1" since we do not want to add the diagonal entry;
    //this is already in sxDiag; these diag entries are added below by using atAddDiagonal

    nelems = jend-jstart;
    if(nelems<=0) continue;

    //upper or lower triangular Q ?  JuMP sends lower, but the test example are inapropriately coded as upper
    if(jcolQ[jend-1]<=i) {
      //lower!
      assert(jcolQ[jstart]<=i);
      if(jcolQ[jend-1]==i)
	kktm.atAddSpRow(i, jcolQ + jstart, dQ+jstart, nelems-1);
      else 
	kktm.atAddSpRow(i, jcolQ + jstart, dQ+jstart, nelems);
    } else {
      //upper
      assert(jcolQ[jstart]>=i && "we expect the upper triangular part of the matrix");
      printf("Warning: input sparse symmetric matrices are in upper triangular format, which will cause performance issues.\n");
      if(jcolQ[jstart]==i) {
	kktm.atAddSpRow(i, jcolQ+jstart+1, dQ+jstart+1, nelems-1);
      } else {
	kktm.atAddSpRow(i, jcolQ+jstart+1, dQ+jstart+1, nelems);
      }


      //kktm.atAddSpRow(i, jcolQ + jstart, dQ+jstart, jend-jstart-1);
    }
    jstart = jend;
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
  // update the KKT with  S diag part 
  ////////////////////////////////////////////////////////////////////
  if(locmz>0) {
    kktm.atAddDiagonal(locnx, ssDiag);
    // int jcol[2]; double M[2];
    // for(int i=locnx; i<locnx+locmz; i++) {
    //   jcol[0] = i; jcol[1] = i+locmz+locmy;
    //   M[0] = ssDiag[i-locnx]; M[1]=-1.;
    //   kktm.atAddSpRow(i, jcol, M, 2);
    // }
    // //kktm.atPutDiagonal(locnx,ssDiag);
    // // for(int i=locnx; i<locnx+locmz; i++) {
    // //   //dKkt[i][i] += ssDiag[i-locnx];
    // // }
  }
  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////
  if(locmy>0){
    //!kktd->symAtAddSubmatrix( locnx+locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx, 1 );
    kktm.symAtAddSubmatrix(locnx+locmz, 0, prob->getLocalB(), 0, 0, locmy, locnx);
    //!for(int i=locnx+locmz; i<locnx+locmz+locmy; i++) dKkt[i][i] += syDiag[i-locnx-locmz];
    
  }
  
  int mle = prob->getmle();
  if(mle>0){
    assert(false && "not yet supported");
    //!   kktd->symAtAddSubmatrix( locnx+locmz+locmy-mle, 0, prob->getLocalE(), 0, 0, mle, locnx, 1 );
  }
  // ////////////////////////////////////////////////////////////////////////////
  // // update the KKT with C   and (dual reg and -I corresponding to \delta z)
  // ///////////////////////////////////////////////////////////////////////////
  if(locmz>0){
    kktm.symAtAddSubmatrix( locnx+locmz+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
    // //   kktd->symAtAddSubmatrix( locnx+locmz+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx, 1 );

    int jcol[2]; double M[2]; M[0]=-1.;
    for(int i=locnx+locmz+locmy; i<locnx+locmz+locmy + locmz; i++) {
      jcol[0] = i-locmz-locmy; 
      jcol[1] = i; M[1] = szDiag[i - locnx-locmz-locmy];
      kktm.atAddSpRow(i, jcol, M, 2);
    }
    // kktm.atPutDiagonal(locnx+locmz+locmy, szDiag);
    // // 	for(int i=0; i<locmz; i++){
    // // 		dKkt[i+locnx+locmz+locmy][i+locnx] -= 1.0;
    // // 		dKkt[i+locnx][i+locnx+locmz+locmy] -= 1.0;
    // // 		dKkt[i+locnx+locmz+locmy][i+locnx+locmz+locmy] += szDiag[i];
  }
  // }
  int mli = prob->getmli();
  if(mli>0){
    assert(false && "not yet supported");
    //   kktd->symAtAddSubmatrix( locnx+locmz+locmy+locmz-mli, 0, prob->getLocalF(), 0, 0, mli, locnx, 1 );
  }

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

#ifdef TIMING
  gprof.t_initializeKKT+=MPI_Wtime()-stime;
#endif

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
  //MPI_Barrier(MPI_COMM_WORLD);
  stime=MPI_Wtime();
  stochNode->resMon.recReduceTmLocal_start();
#endif 
  ///////////////////////////
  reduceKKT();
  ///////////////////////////
#ifdef TIMING
  gprof.t_reduceKKT+=MPI_Wtime()-stime;
#endif

#ifdef TIMING
  stochNode->resMon.recReduceTmLocal_stop();
  stime=MPI_Wtime();
#endif  
  ///////////////////////////
  finalizeKKT(prob, vars);
  ///////////////////////////
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
  //  if(!iAmDistrib) return;
  SparseSymMatrixRowMajList& kktm = dynamic_cast<SparseSymMatrixRowMajList&>(*kkt);

  //get local contribution in triplet format, so that we can reuse MumpsSolver::mumps_->
  //irn_loc, jcn_loc, and a_loc arrays and save on memory
  //when MumpsSolver is replaced with something else, this code needs reconsideration

  MumpsSolver* mumpsSolver = dynamic_cast<MumpsSolver*>(solver);
  assert(mumpsSolver!=NULL);
  int *irn=NULL, *jcn=NULL; double *M=NULL; 
  bool deleteTriplet=false;
  bool bret = mumpsSolver->getTripletStorageArrays(&irn, &jcn, &M);
  if(false==bret) {
    //this is a rank NOT part of MUMPS and does not have the arrays; 
    //allocate them
    irn = new int[kktm.numberOfNonZeros()];
    jcn = new int[kktm.numberOfNonZeros()];
    M = new double[kktm.numberOfNonZeros()];
    deleteTriplet=true;
  } else {
    assert(irn); 
    assert(jcn);
    assert(M);
  }
  
  //root processor bcasts the indexes
  //root first the nnz
  int nnzRoot = kktm.numberOfNonZeros(); 
  MPI_Bcast(&nnzRoot, 1, MPI_INT, 0, mpiComm);
  assert(nnzRoot==kktm.numberOfNonZeros() && "this case is not handled yet");

  //then the indexes
  if(iAmRank0) {
    assert(irn);
    kktm.atGetSparseTriplet(irn, jcn, M, false);
    assert(irn);
  }
  MPI_Bcast(irn, nnzRoot, MPI_INT, 0, mpiComm);
  MPI_Bcast(jcn, nnzRoot, MPI_INT, 0, mpiComm);

  //all processes check: does its sparsity pattern checks out that of the root processor?
  // - NO: (not yet implemented)
  //    - each process adds the entries common with the root, and saves the other ones
  //    - the common entries are MPI_Reduce-d
  //    - the saved entries are then MPI_Gather-ed to the root, who then add-merge them in 
  //    its matrix
  //
  
  // - YES: MPI_Reduce the entries

  if(!iAmRank0) {
    bool bPatternMatched = kktm.fromGetSparseTriplet_w_patternMatch(irn,jcn,nnzRoot,M);
    assert(bPatternMatched && "this is not yet supported");
  }

  //MPI_Reduce the entries
  double* doublebuffer = NULL;
  if(iAmRank0) {
    doublebuffer = new double[nnzRoot];
  }
  MPI_Reduce(M, doublebuffer, nnzRoot, MPI_DOUBLE, MPI_SUM, 0, mpiComm);
  if(doublebuffer) memcpy(M, doublebuffer, nnzRoot*sizeof(double));
  delete [] doublebuffer;

  if(iAmRank0) {
    bret = kktm.atPutSparseTriplet(irn,jcn,M,nnzRoot);
    assert(bret==true && "something went wrong");
  }

  if(deleteTriplet) {
    delete[] irn;
    delete[] jcn;
    delete[] M;
  }
}

extern int gisNLP;

void
sLinsysRootAugSpTriplet::UpdateMatrices( Data * prob_in, int const updateLevel)
{
  int useUpdate=updateLevel;
  if(!gisNLP) useUpdate=1;

  sData* prob = dynamic_cast<sData*>(prob_in);
  SparseSymMatrixRowMajList* kktm = dynamic_cast<SparseSymMatrixRowMajList*>(kkt);
  assert(kktm!=NULL);

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
