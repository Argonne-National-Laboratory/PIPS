/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAug.h"
#include "DeSymIndefSolver.h"
#include "DeSymIndefSolver2.h"
#include "DeSymPSDSolver.h"
#include "PardisoSolver.h"
#include "sData.h"
#include "sTree.h"

#include <unistd.h>

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

sLinsysRootAug::sLinsysRootAug(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_), CtDC(NULL)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
};

sLinsysRootAug::sLinsysRootAug(sFactory* factory_,
			       sData* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_), CtDC(NULL)
{ 
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  redRhs = new SimpleVector(locnx+locmy+locmz);
};

sLinsysRootAug::~sLinsysRootAug()
{
  if(CtDC) delete CtDC;
  delete redRhs;
}


SymMatrix* 
sLinsysRootAug::createKKT(sData* prob)
{
  int n = locnx+locmy;
  return new DenseSymMatrix(n);
}


DoubleLinearSolver*
sLinsysRootAug::createSolver(sData* prob, SymMatrix* kktmat_)
{

  DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
  //return new PardisoSolver(kktmat);
  return new DeSymIndefSolver(kktmat);
  //return new DeSymIndefSolver2(kktmat, locnx); // saddle point solver
  //return new DeSymPSDSolver(kktmat);
}


extern int gLackOfAccuracy;
void sLinsysRootAug::solveReduced( sData *prob, SimpleVector& b)
{
#ifdef TIMING
  double t_total=MPI_Wtime();
  double troot_total=0.0;
  double taux=MPI_Wtime();
#endif

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
  SimpleVector r2(&r[locnx],       locmy);
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

  //SimpleVector realRhs(&r[0], locnx+locmy);
#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
  double tchild_total=0.0;
  double tcomm_total=0.0;
#endif
  double rhsNorm=r.twonorm(); //r== the initial rhs of the reduced system here

  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  SimpleVector rxy(locnx+locmy); rxy.copyFrom(r);
  SimpleVector   x(locnx+locmy); x.setToZero(); //solution
  SimpleVector  dx(locnx+locmy);                //update from iter refinement
  SimpleVector x_prev(locnx+locmy);
  int refinSteps=0;
  std::vector<double> histResid;
  int maxRefinSteps=(gLackOfAccuracy>0?9:8);
  do {
#ifdef TIMING
    taux=MPI_Wtime();
#endif

    x_prev.copyFrom(x);
    //dx = Ainv * r 
    dx.copyFrom(rxy);
    solver->Dsolve(dx);
    //update x
    x.axpy(1.0,dx);

#ifdef TIMING
    troot_total += (MPI_Wtime()-taux);
#endif  

    if(gLackOfAccuracy<0) break;
    if(refinSteps==maxRefinSteps) break;

    //////////////////////////////////////////////////////////////////////
    //iterative refinement
    //////////////////////////////////////////////////////////////////////
    //compute residual
    
    //if (iAmDistrib) {
    //only one process substracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy ] from r
    //                            [  A*xx                       ]
    if(myRank==0) {
      rxy.copyFrom(r);
      if(locmz>0) {
	SparseSymMatrix* CtDC_sp = dynamic_cast<SparseSymMatrix*>(CtDC);
	CtDC_sp->mult(1.0,&rxy[0],1, 1.0,&x[0],1);
      }
      SparseSymMatrix& Q = prob->getLocalQ();
      Q.mult(1.0,&rxy[0],1, -1.0,&x[0],1);
      
      SimpleVector& xDiagv = dynamic_cast<SimpleVector&>(*xDiag);
      assert(xDiagv.length() == locnx);
      for(int i=0; i<xDiagv.length(); i++)
	rxy[i] -= xDiagv[i]*x[i];
      
      SparseGenMatrix& A=prob->getLocalB();
      A.transMult(1.0,&rxy[0],1, -1.0,&x[locnx],1);
      A.mult(1.0,&rxy[locnx],1, -1.0,&x[0],1);
    } else {
      //other processes set r to zero since they will get this portion from process 0
      rxy.setToZero();
    }

#ifdef TIMING
    taux=MPI_Wtime();
#endif  
    // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    SimpleVector xx(&x[0], locnx);
    for(size_t it=0; it<children.size(); it++) {
      children[it]->addTermToSchurResidual(prob->children[it],rxy,xx);  
    }
#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif
    //~done computing residual 

#ifdef TIMING
    taux=MPI_Wtime();
#endif
    //all-reduce residual
    if(iAmDistrib) {
      dx.setToZero(); //we use dx as the recv buffer
      MPI_Allreduce(rxy.elements(), dx.elements(), locnx+locmy, MPI_DOUBLE, MPI_SUM, mpiComm);
      rxy.copyFrom(dx);
    }
#ifdef TIMING
    tcomm_total += (MPI_Wtime()-taux);
#endif

    double relResNorm=rxy.twonorm()/rhsNorm;
    
    if(relResNorm<1.0e-9) {
      break;
    } else {
      double prevRelResNorm=1.0e10;
      if(histResid.size()) 
	prevRelResNorm=histResid[histResid.size()-1];

      //check for stop, divergence or slow convergence conditions
      if(relResNorm>prevRelResNorm) {
	// diverging; restore iteration
	if(myRank==0) {
	  cout << "1st stg - iter refinement diverges relResNorm=" << relResNorm 
	       << "  before was " << prevRelResNorm << endl;
	  cout << "Restoring iterate." << endl;
	  x.copyFrom(x_prev);
	}
	break;
      }else {
	//check slow convergence for the last xxx iterates.
	// xxx is 1 for now
	//if(relResNorm>0.*prevRelResNorm) {

	//  if(myRank==0) {
	//    cout << "1st stg - iter refinement stuck relResNorm=" << relResNorm 
	//	 << "  before was " << prevRelResNorm << endl;
	//    cout << "exiting refinement." << endl;
	//  }
	//  break;
	//
	//} else {
	//  //really nothing, continue
	//}
      }
      histResid.push_back(relResNorm);
      if(myRank==0)
	cout << "1st stg - sol does NOT  have enough accuracy (" << relResNorm << ") after " 
	     << refinSteps << " refinement steps" << endl;
    }
    refinSteps++;
  }while(refinSteps<=maxRefinSteps);

#ifdef TIMING
  taux = MPI_Wtime();
#endif
  r1.copyFrom(x);
  r2.copyFromArray(&x[locnx]);

  //aaa
  ///////////////////////////////////////////////////////////////////////
  // r is the sln to the reduced system
  // the sln to the aug system should be 
  //      x = [r1; r2;  (zDiag)^{-1} * (b3-C*r1);
  ///////////////////////////////////////////////////////////////////////
  SimpleVector b1(&b[0],           locnx);
  SimpleVector b2(&b[locnx],       locmy);
  SimpleVector b3(&b[locnx+locmy], locmz);
  b1.copyFrom(r1);
  b2.copyFrom(r2);

  if(locmz>0) {
    C.mult(1.0, b3, -1.0, r1);
    b3.componentDiv(*zDiag);
  }
#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
  t_total = (MPI_Wtime()-t_total);
  if(myRank==0)
    cout << "ROOT solve+iter refin: " << t_total 
	 << "  ROOT:" << troot_total << "  CHILD: " << tchild_total << "  COMM:" << tcomm_total << endl;
#endif
  //--done
  stochNode->resMon.recDsolveTmLocal_stop();


}

void sLinsysRootAug::finalizeKKT(sData* prob, Variables* vars)
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
  if(locmz>0) {
    SparseGenMatrix& C = prob->getLocalD();
    C.matTransDinvMultMat(*zDiag, &CtDC);
    assert(CtDC->size() == locnx);

    //aliases for internal buffers of CtDC
    SparseSymMatrix* CtDCsp = reinterpret_cast<SparseSymMatrix*>(CtDC);
    int* krowCtDC=CtDCsp->krowM(); int* jcolCtDC=CtDCsp->jcolM(); double* dCtDC=CtDCsp->M();
    
    for(int i=0; i<locnx; i++) {
      pend = krowCtDC[i+1];
      for(p=krowCtDC[i]; p<pend; p++) {
        j = jcolCtDC[p];
        dKkt[i][j] -= dCtDC[p];
	      //printf("%d %d %f\n", i,j,dCtDC[p]);
      }
    }
  } //~end if locmz>0
  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////
  if(locmy>0)
    kktd->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx, 1);
  //prob->getLocalB().getStorageRef().dump("stage1eqmat2.dump");

  /////////////////////////////////////////////////////////////
  // update the KKT zeros for the lower right block 
  /////////////////////////////////////////////////////////////
  //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}
