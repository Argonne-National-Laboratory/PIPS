/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAugSca.h"
#include "QpGenStochData.h"
#include "ScaDenSymMatrix.h"
#include "ScaGenIndefSolver.h"
#include "ScaSymPSDSolver.h"
#include "ScaVector.h"
#include "pipsport.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

sLinsysRootAugSca::sLinsysRootAugSca(sFactory * factory_, QpGenStochData * prob_, COMMINFO &cinfo_)
  : sLinsysRoot(factory_, prob_), CtDC(nullptr)
{ 
  cinfo = cinfo_;
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
  scaRhs = new ScaVector(locnx+locmy, cinfo);
  iAmDistrib = 1;
};

sLinsysRootAugSca::sLinsysRootAugSca(sFactory* factory_,
			       QpGenStochData* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_,
			       COMMINFO &cinfo_)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), CtDC(nullptr)
{ 
  cinfo = cinfo_;
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
  scaRhs = new ScaVector(locnx+locmy, cinfo);
  iAmDistrib = 1;
};

sLinsysRootAugSca::~sLinsysRootAugSca()
{
  if(CtDC) delete CtDC;
  delete redRhs;
  delete scaRhs;
}


SymMatrix* 
sLinsysRootAugSca::createKKT(QpGenStochData* prob)
{
  int n = locnx+locmy;
  return new ScaDenSymMatrix(n, cinfo);
}


DoubleLinearSolver*
sLinsysRootAugSca::createSolver(QpGenStochData* prob, SymMatrix* kktmat_)
{

  ScaDenSymMatrix* kktmat = dynamic_cast<ScaDenSymMatrix*>(kktmat_);
  return new ScaGenIndefSolver(kktmat);
  //return new ScaSymPSDSolver(kktmat);
}

void sLinsysRootAugSca::initializeKKT(QpGenStochData* prob, Variables* vars)
{
  kkt->scalarMult(0.);
  
}



void sLinsysRootAugSca::solveReduced( QpGenStochData *prob, SimpleVector& b)
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

  ///////////////////////////////////////////////////////////////////////
  // compute r1 = b1 - C^T*(zDiag)^{-1}*b3
  ///////////////////////////////////////////////////////////////////////
  if(locmz>0) {
    assert(r3.length() == zDiag->length());
    r3.componentDiv(*zDiag);//r3 is a copy of b3
    C.transMult(1.0, r1, -1.0, r3);
  }
  ///////////////////////////////////////////////////////////////////////
  // r contains all the stuff -> solve for it
  ///////////////////////////////////////////////////////////////////////

  // Plug in scalapack here
  ScaVector& realRhs = dynamic_cast<ScaVector&>(*scaRhs);
  realRhs.copyFromArray(&r[0]);

  solver->Dsolve(realRhs);
  
  realRhs.copyIntoArray(&r[0]);
  // end of scalapack


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
  if(locmz>0) {
    C.mult(1.0, b3, -1.0, r1);
    b3.componentDiv(*zDiag);
  }
  //--done
  stochNode->resMon.recDsolveTmLocal_stop();

  
}

void sLinsysRootAugSca::finalizeKKT(QpGenStochData* prob, Variables* vars)
{
  int j, p, pend; double val;

  stochNode->resMon.recSchurMultLocal_start();

  ScaDenSymMatrix * kktd = dynamic_cast<ScaDenSymMatrix*>(kkt);
 

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
      kktd->symAddAt(i,j,val);
    }
  }

  
  /////////////////////////////////////////////////////////////
  // update the KKT with the diagonals
  // xDiag is in fact diag(Q)+X^{-1}S
  /////////////////////////////////////////////////////////////
  //kktd->atPutDiagonal( 0, *xDiag );
  SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
  for(int i=0; i<locnx; i++) kktd->symAddAt(i,i,sxDiag[i]);


  /////////////////////////////////////////////////////////////
  // update the KKT with   - C' * diag(zDiag) *C
  /////////////////////////////////////////////////////////////
  if(locmz>0) {
    SparseGenMatrix& C = prob->getLocalD();
    C.matTransDinvMultMat(*zDiag, &CtDC);
    assert(CtDC->size() == locnx);
    
    //aliases for internal buffers of CtDC
    SparseSymMatrix* CtDCsp = dynamic_cast<SparseSymMatrix*>(CtDC);
    int* krowCtDC=CtDCsp->krowM(); int* jcolCtDC=CtDCsp->jcolM(); double* dCtDC=CtDCsp->M();
    
    for(int i=0; i<locnx; i++) {
      pend = krowCtDC[i+1];
      for(p=krowCtDC[i]; p<pend; p++) {
	j = jcolCtDC[p];
	kktd->symAddAt(i,j,-dCtDC[p], true);
	//printf("%d %d %f\n", i,j,dCtDC[p]);
      }
    }
  } //~end if locmz>0
  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////
  kktd->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);

  /////////////////////////////////////////////////////////////
  // update the KKT zeros for the lower right block 
  /////////////////////////////////////////////////////////////
  //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);

  stochNode->resMon.recSchurMultLocal_stop();
}

void sLinsysRootAugSca::factor2(QpGenStochData *prob, Variables *vars)
{
  ScaDenSymMatrix& kktd = dynamic_cast<ScaDenSymMatrix&>(*kkt);
  int blocksize = cinfo.nb, nxP = locnx;
  initializeKKT(prob, vars);
  int zero = 0;
  // we're only sending upper block, find out how many rows this processor has of it
  int local_nr = FNAME(numroc)(&nxP, &cinfo.nb, &cinfo.myrow, &zero, &cinfo.nprow);
  local_nr = MAX(local_nr, 0);

  DenseGenMatrix colbuffer(blocksize, nxP);
  double *recvbuffer = new double[local_nr*blocksize];
  double *sendbuffer = new double[nxP*blocksize];
  int *recvcounts = new int[cinfo.nprocs];
  
  // First tell children to factorize.
  for(unsigned int c=0; c<children.size(); c++) {
    children[c]->factor2(prob->children[c], vars);
  }
  
  for (int curcol = 0; curcol < nxP; curcol += blocksize) {
    int proccol = (curcol/blocksize) % cinfo.npcol;
    int endcol = MIN(curcol+blocksize, nxP); // exclusive
    int numcols = endcol-curcol;
    memset(&colbuffer[0][0], 0, blocksize*nxP*sizeof(double));
    
    for (unsigned int c=0; c<children.size(); c++) {
      if(children[c]->mpiComm == MPI_COMM_NULL)
      	continue;
    
      children[c]->stochNode->resMon.recFactTmChildren_start();    
      //---------------------------------------------
      children[c]->addColsToDenseSchurCompl(prob->children[c], colbuffer, curcol, endcol);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();    
    }
    // only to improve timing of reduce 
    //MPI_Barrier(cinfo.mpicomm);
    
    stochNode->resMon.recReduceTmLocal_start();

    // now the fun part, first rearrange the elements so that
    // we can call reducescatter
    // that is, all elements belonging to the first proc go first, etc
    // we send a single column of blocks per iteration
    int desti = 0;
    for (int prow = 0; prow < cinfo.nprow; prow++) {
      for (int j = 0; j < numcols; j++) {
        for (int currow = prow*blocksize; currow < nxP; currow += cinfo.nprow*blocksize) {
          int endrow = MIN(currow+blocksize, nxP);
          int numrows = endrow - currow;
          memcpy(sendbuffer+desti,&colbuffer[j][currow], numrows*sizeof(double));
          desti += numrows;
        }
      }
    }
    assert(desti == nxP*numcols);
    
    memset(recvcounts, 0, cinfo.nprocs*sizeof(int));
    for (int prow = 0; prow < cinfo.nprow; prow++) {
      int destproc = get_pnum(cinfo, prow, proccol);
      int zero = 0, proc_nr;
      proc_nr = FNAME(numroc)(&nxP, &cinfo.nb, &prow, &zero, &cinfo.nprow);
      if (proc_nr > 0) recvcounts[destproc] = proc_nr*numcols;
    }

    
    stochNode->resMon.recReduceScatterTmLocal_start();
    MPI_Reduce_scatter(sendbuffer, recvbuffer, recvcounts, MPI_DOUBLE, 
      MPI_SUM, cinfo.mpicomm);
    stochNode->resMon.recReduceScatterTmLocal_stop();

    // now unpack on the local processor
    // each column is already in continuous memory
    // if zero equality constraints, actually the whole block is continous
    if ( recvcounts[cinfo.mype] > 0 ) {
      desti = 0;
      int lcol, lrow;
      global2local(0, curcol, cinfo, lrow, lcol);
      for (int j = 0; j < numcols; j++) {
        memcpy(kktd.mat->data+(j+lcol)*kktd.getNR(),recvbuffer+j*local_nr,
          local_nr*sizeof(double));
      }
    }

    stochNode->resMon.recReduceTmLocal_stop(); 

    
  }
  
  delete [] recvcounts;
  delete [] recvbuffer;
  delete [] sendbuffer;
  
  finalizeKKT(prob, vars);
  
  //double val = kktd.getVal(PROW,PCOL);
  //if (cinfo.mype == 0) {
  //  printf("(%d,%d) --- %f\n", PROW, PCOL, val);
  //}

  factorizeKKT();


#ifdef TIMING
  afterFactor();
#endif

}



/*
// Old reduce version
void sLinsysRootAugSca::factor2(QpGenStochData *prob, Variables *vars)
{
  ScaDenSymMatrix& kktd = dynamic_cast<ScaDenSymMatrix&>(*kkt);
  int blocksize = cinfo.nb, nxP = locnx;
  initializeKKT(prob, vars);
  
  DenseGenMatrix colbuffer(blocksize, nxP);
  double *recvbuffer = new double[blocksize*blocksize];
  double *sendbuffer = new double[blocksize*blocksize];
  
  // First tell children to factorize.
  for(int c=0; c<children.size(); c++) {
    children[c]->factor2(prob->children[c], vars);
  }
  
  for (int curcol = 0; curcol < nxP; curcol += blocksize) {
   int proccol = (curcol/blocksize) % cinfo.npcol;
   int endcol = MIN(curcol+blocksize, nxP); // exclusive
   int numcols = endcol-curcol;
   memset(&colbuffer[0][0], 0, blocksize*nxP*sizeof(double));
   
    for (int c=0; c<children.size(); c++) {
      if(children[c]->mpiComm == MPI_COMM_NULL)
      	continue;
    
      children[c]->stochNode->resMon.recFactTmChildren_start();    
      //---------------------------------------------
      children[c]->addColsToDenseSchurCompl(prob->children[c], colbuffer, curcol, endcol);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();    
    }
    stochNode->resMon.recReduceTmLocal_start(); 
    // now the fun part, reduce each block to where it needs to be in the scalapack matrix
    for (int currow = 0; currow < nxP; currow += blocksize) {
      int procrow = (currow/blocksize) % cinfo.nprow;
      int destproc = get_pnum(cinfo, procrow, proccol);
      int endrow = MIN(currow+blocksize, nxP);
      int numrows = endrow - currow; // # rows in block to send

      for (int j=0; j < numcols; j++) {
      	memcpy(sendbuffer+j*numrows, &colbuffer[j][currow], numrows*sizeof(double));
      }
      MPI_Reduce(sendbuffer, recvbuffer, numrows*numcols, MPI_DOUBLE,
       MPI_SUM, destproc, cinfo.mpicomm);
      if (cinfo.mype == destproc) {
        int lrow, lcol;
        global2local(currow, curcol, cinfo, lrow, lcol);
        for (int j=0; j < numcols; j++) {
          int doffset = (lcol+j)*kktd.getNR()+lrow;
          memcpy(kktd.mat->data+doffset, recvbuffer+j*numrows, numrows*sizeof(double));
        }
      }     
    }
    stochNode->resMon.recReduceTmLocal_stop(); 

    
  }
  
  delete [] recvbuffer;
  delete [] sendbuffer;
  
  finalizeKKT(prob, vars);
  
  //double val = kktd.getVal(PROW,PCOL);
  //if (cinfo.mype == 0) {
  //  printf("(%d,%d) --- %f\n", PROW, PCOL, val);
  //}

  factorizeKKT();

#ifdef TIMING
  afterFactor();
#endif

}

*/
