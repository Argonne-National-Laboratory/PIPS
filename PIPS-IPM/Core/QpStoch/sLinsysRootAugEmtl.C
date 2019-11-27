/* PIPS
   Authors: Miles Lubin and Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAugEmtl.h"
#include "QpGenStochData.h"
#include "EmtlDenSymMatrix.h"
#include "EmtlGenIndefSolver.h"
#include "EmtlSymPSDSolver.h"
#include "EmtlVector.h"
#include "pipsport.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

sLinsysRootAugEmtl::sLinsysRootAugEmtl(sFactory * factory_, QpGenStochData * prob_, const EmtlContext &ctx_)
  : sLinsysRoot(factory_, prob_), ctx(ctx_), CtDC(nullptr)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = this->createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
  emtlRhs = new EmtlVector(locnx+locmy, ctx);
  iAmDistrib = 1;
};

sLinsysRootAugEmtl::sLinsysRootAugEmtl(sFactory* factory_,
			       QpGenStochData* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_,
			       const EmtlContext &ctx_)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), 
    ctx(ctx_), CtDC(nullptr)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = this->createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
  emtlRhs = new EmtlVector(locnx+locmy, ctx);
  iAmDistrib = 1;
};

sLinsysRootAugEmtl::~sLinsysRootAugEmtl()
{
  if(CtDC) delete CtDC;
  delete redRhs;
  delete emtlRhs;
}


SymMatrix* 
sLinsysRootAugEmtl::createKKT(QpGenStochData* prob)
{
  int n = locnx+locmy;
  return new EmtlDenSymMatrix(n, ctx);
}


DoubleLinearSolver*
sLinsysRootAugEmtl::createSolver(QpGenStochData* prob, SymMatrix* kktmat_)
{

  EmtlDenSymMatrix* kktmat = dynamic_cast<EmtlDenSymMatrix*>(kktmat_);
  return new EmtlGenIndefSolver(kktmat);
}

void sLinsysRootAugEmtl::initializeKKT(QpGenStochData* prob, Variables* vars)
{
  kkt->scalarMult(0.);
  
}



void sLinsysRootAugEmtl::solveReduced( QpGenStochData *prob, SimpleVector& b)
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

  // Plug in elemental here
  EmtlVector& realRhs = dynamic_cast<EmtlVector&>(*emtlRhs);
  realRhs.copyFromArray(&r[0]);

  solver->Dsolve(realRhs);
  
  realRhs.copyIntoArray(&r[0]);
  // end of elemental


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

void sLinsysRootAugEmtl::finalizeKKT(QpGenStochData* prob, Variables* vars)
{
  int j, p, pend; double val;

  stochNode->resMon.recSchurMultLocal_start();

  EmtlDenSymMatrix * kktd = dynamic_cast<EmtlDenSymMatrix*>(kkt);
 

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

static inline int firstcol(const int mycol,const int startcol,const int npcol)
{
  return (mycol - startcol + (startcol/npcol+1)*npcol)%npcol;
}

const double MAX_MB_FOR_COL_BUFFERS = 100;

void sLinsysRootAugEmtl::factor2(QpGenStochData *prob, Variables *vars)
{
  EmtlDenSymMatrix& kktd = dynamic_cast<EmtlDenSymMatrix&>(*kkt);
  int nxP = locnx;
  
  const int BLOCKSIZE = MIN((1048576*MAX_MB_FOR_COL_BUFFERS/
                            (2*sizeof(double)*nxP)),nxP);

  initializeKKT(prob, vars);
  // we're only sending upper Q block,
  // count how many elements each processor has of this block per column
  // indexed by processor row
  int *nr_counts = new int[ctx.nprow()];
  for (int i = 0; i < ctx.nprow(); i++) {
    nr_counts[i] = utilities::LocalLength(nxP, i, ctx.nprow());
  }
  int max_nr = utilities::MaxLocalLength(nxP, ctx.nprow());

  DenseGenMatrix colbuffer(BLOCKSIZE, nxP);
  double *recvbuffer = new double[max_nr*BLOCKSIZE];
  double *sendbuffer = new double[nxP*BLOCKSIZE];
  int *recvcounts = new int[ctx.nprocs()];
  memset(recvcounts, 0, ctx.nprocs()*sizeof(int));

  //printf("got to factorize\n");
  
  // First tell children to factorize.
  for(unsigned int c=0; c<children.size(); c++) {
    children[c]->factor2(prob->children[c], vars);
  }
  
  for (int startcol = 0; startcol < nxP; startcol += BLOCKSIZE) {
    int endcol = MIN(startcol+BLOCKSIZE, nxP); // exclusive
    int numcols = endcol-startcol;
    /*if (ctx.mype() == 0) {
      printf("startcol: %d endcol: %d\n", startcol, endcol);
    }*/
    memset(&colbuffer[0][0], 0, BLOCKSIZE*nxP*sizeof(double));
    
    for (unsigned int c=0; c<children.size(); c++) {
      if(children[c]->mpiComm == MPI_COMM_NULL)
      	continue;
    
      children[c]->stochNode->resMon.recFactTmChildren_start();    
      //---------------------------------------------
      children[c]->addColsToDenseSchurCompl(prob->children[c], colbuffer, startcol, endcol);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();    
    }

    // only to improve timing of reduce 
    //MPI_Barrier(ctx.comm());

    //printf("pe %d got columns\n", ctx.mype());
    stochNode->resMon.recReduceTmLocal_start(); 
    // now the fun part, first rearrange the elements so that
    // we can call reducescatter
    // that is, all elements belonging to the first proc go first, etc
    // TODO: this won't work when using the torus
    assert(!ctx.usingTorus());
    int desti = 0;
    for (int pcol = 0; pcol < ctx.npcol(); pcol++)
    for (int prow = 0; prow < ctx.nprow(); prow++) {
      for(int j = firstcol(pcol,startcol,ctx.npcol()); 
            j < numcols; j+= ctx.npcol()) {
        //printf("pe %d for %d %d at col %d\n",ctx.mype(),prow,pcol,j);
        for (int i = prow; i < nxP; i += ctx.nprow()) {
          sendbuffer[desti++] = colbuffer[j][i];
        }
      }
      //printf("pe %d loaded buffer for %d %d\n",ctx.mype(),prow,pcol); 
    }
    assert(desti == nxP*numcols);
    //printf("pe %d loaded send buffer\n", ctx.mype());
    
    int destproc = 0;
    for (int pcol = 0; pcol < ctx.npcol(); pcol++) {
      for (int prow = 0; prow < ctx.nprow(); prow++) {
        // how many columns are we sending to this processor column
        int ncols = utilities::LocalLength(
                      numcols-firstcol(pcol,startcol,ctx.npcol()), 
                      0, ctx.npcol());
        recvcounts[destproc++] = nr_counts[prow]*ncols;
      }
    }

    stochNode->resMon.recReduceScatterTmLocal_start();
    MPI_Reduce_scatter(sendbuffer, recvbuffer, recvcounts, MPI_DOUBLE, 
      MPI_SUM, ctx.comm());
    stochNode->resMon.recReduceScatterTmLocal_stop();


    // now unpack on the local processor
    // each column is already in continuous memory
    // if zero equality constraints, actually the whole block is continous
    if ( recvcounts[ctx.mype()] > 0 ) {
      desti = 0;
      int lcol = (startcol+firstcol(ctx.mycol(),startcol,ctx.npcol()))
                    /ctx.npcol();
      int local_nr = nr_counts[ctx.myrow()];
      int ncols = recvcounts[ctx.mype()]/local_nr;
      //printf("%d: %d %d %d %d\n", ctx.mype(), lcol, local_nr, numcols,ncols); 
      for (int j = 0; j < ncols; j++) {
        //printf("%d %d\n", ctx.mype(), j);
        memcpy(kktd.mat->data+(j+lcol)*kktd.getNR(),recvbuffer+j*local_nr,
          local_nr*sizeof(double));
      }
    }

    stochNode->resMon.recReduceTmLocal_stop(); 

    
  }
  //printf("done factorizing\n");

  delete [] recvcounts;
  delete [] recvbuffer;
  delete [] sendbuffer;
  delete [] nr_counts;
  
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



