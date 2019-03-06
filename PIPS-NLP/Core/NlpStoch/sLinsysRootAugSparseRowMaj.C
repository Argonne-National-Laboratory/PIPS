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
extern int gMUMPSranks;

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
    if(nodeRank<gMUMPSranks)
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
  SparseSymMatrixRowMajList& kktm = dynamic_cast<SparseSymMatrixRowMajList&>(*kkt);
  int nnzRoot = kktm.numberOfNonZeros(); 

  //get local contribution in triplet format, so that we can reuse MumpsSolver::mumps_->
  //irn_loc, jcn_loc, and a_loc arrays and save on memory
  //when MumpsSolver is replaced with another solver, this code needs reconsideration
  MumpsSolver* mumpsSolver = dynamic_cast<MumpsSolver*>(solver);
  assert(mumpsSolver!=NULL);
  int *irn=NULL, *jcn=NULL; double *M=NULL; 
  bool deleteTriplet=false;

  if(iAmRank0) {
    bool bret = mumpsSolver->getTripletStorageArrays(&irn, &jcn, &M);
    if(bret==false) {
      irn = new int[nnzRoot];
      jcn = new int[nnzRoot];
      M = new double[nnzRoot];
      deleteTriplet=true;
    } else {
      assert(irn); assert(jcn); assert(M);
      deleteTriplet=false;
    }
  } 

  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  //printf("myRank %d ---- nnzRoot %d   localkktm %d \n", myRank, nnzRoot, kktm.numberOfNonZeros());

  //
  //root processor bcasts the its triplet indexes
  //

  //root bcasts first the nnz
  MPI_Bcast(&nnzRoot, 1, MPI_INT, 0, mpiComm);

  if(!iAmRank0) {
    //allocate
    irn = new int[nnzRoot];
    jcn = new int[nnzRoot];
    M = new double[nnzRoot];
    deleteTriplet=true;
  } else {
    //root prepares the sparse triplet
    assert(irn);
    kktm.atGetSparseTriplet(irn, jcn, M, false);
  }
  MPI_Bcast(irn, nnzRoot, MPI_INT, 0, mpiComm);
  MPI_Bcast(jcn, nnzRoot, MPI_INT, 0, mpiComm);

  //all processes check: does its sparsity pattern check that of the root processor?
  //    - each process adds the entries in common with the root, and saves the not-in-the-root/diff 
  //    entries
  //    - the common entries are MPI_Reduce-d
  //    - the not-in-the-root entries are then MPI_Gather-ed to the root; the root rank then 
  //    add-merge them in its matrix
  //
  int *irow_diff=NULL, *jcol_diff=NULL; double* M_diff=NULL; int nnz_diff=0;
  if(!iAmRank0) {

    bool bPatternMatched = kktm.fromGetIntersectionSparseTriplet_w_diff(irn,jcn,nnzRoot,M,
									&irow_diff, &jcol_diff, &M_diff, nnz_diff);
    if(bPatternMatched) { 
      assert(nnz_diff==0); 
      assert(irow_diff==NULL); 
      assert(jcol_diff==NULL); 
      assert(M_diff==NULL); 
    } else {
      assert(nnz_diff!=0); 
      assert(irow_diff!=NULL); 
      assert(jcol_diff!=NULL);
      assert(M_diff!=NULL); 
    }
    //printf("Rank %d  has %d entries that rank0 does not have.\n", myRank, nnz_diff);
  }
  // - the common entries are MPI_Reduce-d
  {
    double* doublebuffer = NULL;
    if(iAmRank0) {
      doublebuffer = new double[nnzRoot];
    }
    MPI_Reduce(M, doublebuffer, nnzRoot, MPI_DOUBLE, MPI_SUM, 0, mpiComm);
    if(doublebuffer) memcpy(M, doublebuffer, nnzRoot*sizeof(double));
    delete [] doublebuffer;

    if(iAmRank0) {
      bool bret = kktm.atPutSparseTriplet(irn,jcn,M,nnzRoot);
      assert(bret==true && "sparsity patterns do not match on rank 0; this should not happen");
    }
  }

  // - the saved entries are then MPI_Gather-ed 
  {
    int mismatch=(nnz_diff>0), auxG;
    MPI_Allreduce(&mismatch, &auxG, 1, MPI_INT, MPI_MAX, mpiComm);
    mismatch=auxG;
    if(mismatch) {

      //first get the sum of nnz of diff over all processes -> needed to allocate recv buffer at root
      int commSize; MPI_Comm_size( mpiComm, &commSize); 
      int* diff_counts = NULL, *irow_diff_dest=NULL, *jcol_diff_dest=NULL; double* M_diff_dest=NULL;
      if(iAmRank0) diff_counts = new int[commSize];

      MPI_Gather(&nnz_diff, 1, MPI_INT, diff_counts, 1, MPI_INT, 0, mpiComm);

      size_t* displs=NULL;
      if(iAmRank0) {
	assert(diff_counts[0]==0); //diff for rank 0 should be empty
	displs = new size_t[commSize];

	displs[0]=0;
	for(int i=1; i<commSize; i++) { 
	  displs[i] = displs[i-1] + diff_counts[i-1];
	  //printf("Rank %d  -> diff_counts=%d\n", i, diff_counts[i]);
	}

	size_t nnz_diff_total = displs[commSize-1]+diff_counts[commSize-1];

	//printf("Rank %d  -> %d entries to be gathered\n", myRank, nnz_diff_total);

	irow_diff_dest = new int[nnz_diff_total];
	jcol_diff_dest = new int[nnz_diff_total];
	M_diff_dest = new double[nnz_diff_total];
      }
      
      if(iAmRank0) {
        MPI_Request request[3][commSize-1];
        for(int i=1; i < commSize; i++) {
          MPI_Irecv(irow_diff_dest+displs[i], diff_counts[i], MPI_INT, i, 0, mpiComm, &request[0][i-1]);
          MPI_Irecv(jcol_diff_dest+displs[i], diff_counts[i], MPI_INT, i, 0, mpiComm, &request[1][i-1]);
          MPI_Irecv(M_diff_dest+displs[i], diff_counts[i], MPI_DOUBLE, i, 0, mpiComm, &request[2][i-1]);
        }
        for(int i=0; i<nnz_diff; i++) irow_diff_dest[i] = irow_diff[i];
        for(int i=0; i<nnz_diff; i++) jcol_diff_dest[i] = jcol_diff[i];
        for(int i=0; i<nnz_diff; i++) M_diff_dest[i] = M_diff[i];
        MPI_Waitall(commSize-1, request[0], MPI_STATUSES_IGNORE);
        MPI_Waitall(commSize-1, request[1], MPI_STATUSES_IGNORE);
        MPI_Waitall(commSize-1, request[2], MPI_STATUSES_IGNORE);
      }
      else {
        MPI_Send(irow_diff, nnz_diff, MPI_INT, 0, 0, mpiComm);
        MPI_Send(jcol_diff, nnz_diff, MPI_INT, 0, 0, mpiComm);
        MPI_Send(M_diff, nnz_diff, MPI_DOUBLE, 0, 0, mpiComm);
      }



      if(iAmRank0) {	
	int nz_start=0, nz_end;
	for(int p=0; p<commSize; p++) {
	  if(diff_counts[p]==0) continue; 

	  //printf("Rank %d  adding %d entries from rank %d\n", myRank, diff_counts[p], p);

	  nz_end=nz_start+diff_counts[p];

	  //for the diff coming from rank p, go over the nnz and add each row to kktm
	  //we assume the row indexes are ordered, and for equal row indexes the col indexes are ordered
  
	  int row_start = nz_start, row_end = nz_start;
	  while(row_end<nz_end) {
	    
	    while(irow_diff_dest[row_start] == irow_diff_dest[row_end] && row_end<nz_end) 
	      row_end++;
	    assert(row_end>=row_start);

	    kktm.atAddSpRow(irow_diff_dest[row_start], jcol_diff_dest+row_start, M_diff_dest+row_start, row_end-row_start);
	    row_start = row_end;
	  }

	  nz_start = nz_end;

	} // end for(int p=0; p<commSize; p++) 
	
      } // end if(iAmRank0) {
      
      delete [] diff_counts;
      delete [] irow_diff_dest;
      delete [] jcol_diff_dest;
      delete [] M_diff_dest;
    } // end if(mismatch)
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
