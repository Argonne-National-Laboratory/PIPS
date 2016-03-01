/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "sLinsysRootAug.h"
#include "DeSymIndefSolver.h"
#include "sData.h"
#include "sTree.h"

#include <unistd.h>
#include "math.h"

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif
extern int gInnerSCsolve;
extern int gOuterSolve;
extern int separateHandDiag;

using namespace std;

sLinsysRootAug::sLinsysRootAug(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_), CtDC(NULL)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);
  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  assert(gOuterSolve>=3);
  redRhs = new SimpleVector(locnx+locmz+locmy+locmz);
};

sLinsysRootAug::sLinsysRootAug(sFactory* factory_,
			       sData* prob_,
			       OoqpVector* dd_, 
			       OoqpVector* dq_,
			       OoqpVector* nomegaInv_,
			       OoqpVector* rhs_,
			       OoqpVector* additiveDiag_)
  : sLinsysRoot(factory_, prob_, dd_, dq_, nomegaInv_, rhs_, additiveDiag_), CtDC(NULL)
{ 
  prob_->getLocalSizes(locnx, locmy, locmz);

  kkt = createKKT(prob_);
  solver = createSolver(prob_, kkt);
  assert(gOuterSolve>=3);  
  redRhs = new SimpleVector(locnx+locmz+locmy+locmz);
};

sLinsysRootAug::~sLinsysRootAug()
{
  if(CtDC) delete CtDC;
  delete redRhs;
}


SymMatrix* 
sLinsysRootAug::createKKT(sData* prob)
{
  int n;

  if(gOuterSolve < 3){
    n = locnx+locmy;
    assert(locmz==0);
  }else{
    n = locnx+locmy+locmz+locmz;
  }

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

#ifdef TIMING
extern double t_start, troot_total, taux, tchild_total, tcomm_total;
#endif


extern int gLackOfAccuracy;
void sLinsysRootAug::solveReduced( sData *prob, SimpleVector& b)
{
//  assert(locmz==0||gOuterSolve<3);
  assert(gOuterSolve>=3);

  int myRank; MPI_Comm_rank(mpiComm, &myRank);

#ifdef TIMING
  t_start=MPI_Wtime();
  troot_total=tchild_total=tcomm_total=0.0; 
#endif

  assert(locnx+locmz+locmy+locmz==b.length());
  SimpleVector& r = (*redRhs);
  assert(r.length() == b.length());
  SparseGenMatrix& C = prob->getLocalD();

  ///////////////////////////////////////////////////////////////////////
  // LOCAL SOLVE
  ///////////////////////////////////////////////////////////////////////
 
  ///////////////////////////////////////////////////////////////////////
  // b=[b1;b2;b3;b4] is a locnx+locmz+locmy+locmz vector 
  ///////////////////////////////////////////////////////////////////////

  r.copyFromArray(b.elements());
  // aliases to parts (no mem allocations)  
  SimpleVector r4(&r[locnx+locmz+locmy], locmz);
  SimpleVector r3(&r[locnx+locmz], locmy); //r3 is used as a temp buffer for b3
  SimpleVector r2(&r[locnx],       locmz);
  SimpleVector r1(&r[0],           locnx);


  if(gInnerSCsolve==0) {
    // Option 1. - solve with the factors
    solver->Dsolve(r);
  } else if(gInnerSCsolve==1) {
    // Option 2 - solve with the factors and perform iter. ref.
    solveWithIterRef(prob, r);
  } else {
    assert(gInnerSCsolve==2);
    // Option 3 - use the factors as preconditioner and apply BiCGStab
    solveWithBiCGStab(prob, r);
  }
  ///////////////////////////////////////////////////////////////////////
  // r is the sln to the reduced system
  ///////////////////////////////////////////////////////////////////////
  SimpleVector b1(&b[0],           locnx);
  SimpleVector b2(&b[locnx],       locmz);
  SimpleVector b3(&b[locnx+locmz], locmy);
  SimpleVector b4(&b[locnx+locmz+locmy], locmz);  
  b1.copyFrom(r1);
  b2.copyFrom(r2);
  b3.copyFrom(r3);
  b4.copyFrom(r4);

#ifdef TIMING
  if(myRank==0 && gInnerSCsolve>=1)
    cout << "Root - Refin times: child=" << tchild_total << " root=" << troot_total
	 << " comm=" << tcomm_total << " total=" << MPI_Wtime()-t_start << endl;
#endif
}

/** rxy = beta*rxy + alpha * SC * x */
void sLinsysRootAug::SCmult( double beta, SimpleVector& rxy, 
			     double alpha, SimpleVector& x, 
			     sData* prob)
{
  //if (iAmDistrib) {
  //only one process substracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy ] from r
  //                            [  A*xx                       ]
  int myRank; MPI_Comm_rank(mpiComm, &myRank);
  if(myRank==0) {
    //only this proc substracts from rxy
    rxy.scalarMult(beta);
    SparseSymMatrix& Q = prob->getLocalQ();
    Q.mult(1.0,&rxy[0],1, alpha,&x[0],1);

    if(locmz>0) {
      SparseSymMatrix* CtDC_sp = dynamic_cast<SparseSymMatrix*>(CtDC);
      CtDC_sp->mult(1.0,&rxy[0],1, alpha,&x[0],1);
    }
    
    SimpleVector& xDiagv = dynamic_cast<SimpleVector&>(*xDiag);
    assert(xDiagv.length() == locnx);
    for(int i=0; i<xDiagv.length(); i++)
      rxy[i] += alpha*xDiagv[i]*x[i];
    
    SparseGenMatrix& A=prob->getLocalB();
    A.transMult(1.0,&rxy[0],1, alpha,&x[locnx],1);
    A.mult(1.0,&rxy[locnx],1, alpha,&x[0],1);
  } else {
    //other processes set r to zero since they will get this portion from process 0
    rxy.setToZero();
  }
  
#ifdef TIMING
    taux=MPI_Wtime();
#endif  
    // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    SimpleVector xx(locnx);
    xx.copyFromArray(x.elements());
    xx.scalarMult(-alpha);
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
      SimpleVector buf(rxy.length());
      buf.setToZero(); //we use dx as the recv buffer
      MPI_Allreduce(rxy.elements(), buf.elements(), locnx+locmy, MPI_DOUBLE, MPI_SUM, mpiComm);
      rxy.copyFrom(buf);
    }
#ifdef TIMING
    tcomm_total += (MPI_Wtime()-taux);
#endif

}

void sLinsysRootAug::solveWithIterRef( sData *prob, SimpleVector& r)
{
  SimpleVector r2(&r[locnx],       locmy);
  SimpleVector r1(&r[0],           locnx);

  //SimpleVector realRhs(&r[0], locnx+locmy);
#ifdef TIMING
  taux=MPI_Wtime();
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
  do { //iterative refinement
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
    
    if(relResNorm<1.0e-10) {
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
	}
	x.copyFrom(x_prev);
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

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
#endif  
}

void sLinsysRootAug::solveWithBiCGStab( sData *prob, SimpleVector& b)
{
  int n = b.length();

  const int maxit=500;
  const double tol=1e-12, EPS=2e-16;
  double iter=0.0;

  int myRank; MPI_Comm_rank(mpiComm, &myRank);

  SimpleVector r(n);           //residual
  SimpleVector s(n);           //residual associated with half iterate
  SimpleVector rt(n);          //shadow residual
  SimpleVector xmin(n);        //minimal residual iterate
  SimpleVector x(n);           //iterate
  SimpleVector xhalf(n);       // half iterate of BiCG
  SimpleVector p(n),paux(n);
  SimpleVector v(n), t(n);
  int flag; double imin;
  double n2b;                  //norm of b 
  double normr, normrmin;      //norm of the residual and norm of residual at min-resid iterate
  double normr_act;
  double tolb;                 //relative tolerance
  double rho, omega, alpha;
  int stag, maxmsteps, maxstagsteps, moresteps;
  double relres;
  //maxit = n/2+1;

  //////////////////////////////////////////////////////////////////
  //  Problem Setup and intialization
  //////////////////////////////////////////////////////////////////

  n2b = b.twonorm();
  tolb = n2b*tol;


#ifdef TIMING
  taux = MPI_Wtime();
#endif
  //initial guess
  x.copyFrom(b);
  solver->Dsolve(x);
  //initial residual
  r.copyFrom(b);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
  taux = MPI_Wtime();
#endif 
 
  //applyA(1.0, r, -1.0, x);
  SCmult(1.0,r, -1.0,x, prob);

#ifdef TIMING
    tchild_total +=  (MPI_Wtime()-taux);
#endif

  normr=r.twonorm();
#ifdef TIMING
  if(myRank==0)
    cout << "BiCG: initial rel resid: " << normr/n2b << endl;
#endif

  if(normr<tolb) {
    //initial guess is good enough
    b.copyFrom(x); flag=0; return;
  }

  rt.copyFrom(r); //Shadow residual
  double* resvec = new double[2*maxit+1];
  resvec[0] = normr; normrmin=normr;
  rho=1.0; omega=1.0;
  stag=0; maxmsteps=min(min(n/50, 5), n-maxit); 
  maxstagsteps=3; moresteps=0;


  //////////////////////////////////////////////////////////////////
  // loop over maxit iterations
  //////////////////////////////////////////////////////////////////
  int ii=0; while(ii<maxit) {
    //cout << ii << " ";
    flag=-1;
    ///////////////////////////////
    // First half of the iterate
    ///////////////////////////////
    double rho1=rho; double beta;
    rho = rt.dotProductWith(r); 
    //printf("rho=%g\n", rho);
    if(0.0==rho) { flag=4;  break; }

    if(ii==0) p.copyFrom(r);
    else {
      beta = (rho/rho1)*(alpha/omega);
      if(beta==0.0) { flag=4;  break; }

      //-------- p = r + beta*(p - omega*v) --------
      p.axpy(-omega, v); p.scale(beta); p.axpy(1.0, r);
    }

#ifdef TIMING
    taux = MPI_Wtime();
#endif
    //------ v = A*(M2inv*(M1inv*p)) and ph=M2inv*(M1inv*p)
    //first use v as temp storage
    //applyM1(0.0, v,    1.0, p);
    //applyM2(0.0, paux, 1.0, v);
    //applyA (0.0, v,    1.0, paux); 
    paux.copyFrom(p);
    solver->solve(paux);

#ifdef TIMING
  troot_total += (MPI_Wtime()-taux);
#endif 
    
    SCmult(0.0,v, 1.0,paux, prob);

    
    SimpleVector& ph = paux;

    double rtv = rt.dotProductWith(v);
    if(rtv==0.0) { flag=4; break; }

    alpha = rho/rtv;
    if(fabs(alpha)*ph.twonorm()<EPS*x.twonorm()) stag++;
    else                                         stag=0;

    // xhalf = x + alpha*ph and the associated residual
    xhalf.copyFrom(x); xhalf.axpy( alpha, ph);
    s.    copyFrom(r);     s.axpy(-alpha, v);
    normr = s.twonorm(); normr_act = normr;
    resvec[2*ii] = normr;

    //printf("iter %g normr=%g\n", ii+0.5, normr);
    //-------- check for convergence in the middle of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      s.copyFrom(b);
      //applyA(1.0, s, -1.0, xhalf); // s=b-Ax
      SCmult(1.0,s, -1.0,xhalf, prob);
      normr_act = s.twonorm();
      
      if(normr<=tolb) {
	//converged
	x.copyFrom(xhalf);	
	flag = 0; iter = 0.5+ii;
	break;
      } else {
	if(stag>=maxstagsteps && moresteps==0) {
	  stag=0;
	}
	moresteps++;
	if(moresteps>=maxmsteps) {
	  //method stagnated
	  flag=3; x.copyFrom(xhalf);
	  break;
	}
      }
    }
    if(stag>=maxstagsteps) { flag=3; break;} //stagnation

    //update quantities related to minimal norm iterate
    if(normr_act<normrmin) {
      xmin.copyFrom(xhalf); normrmin=normr_act;
      imin=0.5+ii;
    }

#ifdef TIMING
    taux = MPI_Wtime();
#endif
    ///////////////////////////////
    // Second half of the iterate
    //////////////////////////////
    //applyM1(0.0, t,    1.0, s); //applyM1(s,     stemp);
    //applyM2(0.0, paux, 1.0, t); //applyM2(stemp, sh);
    //applyA (0.0, t,    1.0, paux); //applyA (sh, t);
    //kkt->mult(0.0,paux, 1.0,s);
    paux.copyFrom(s);
    solver->solve(paux);
#ifdef TIMING
    troot_total += (MPI_Wtime()-taux);
#endif

    SCmult(0.0,t, 1.0,paux, prob);


    SimpleVector& sh = paux; 
    double tt = t.dotProductWith(t);
    if(tt==0.0) { flag=4; break;}

    omega=t.dotProductWith(s); omega /= tt;

    if(fabs(omega)*sh.twonorm() < EPS*xhalf.twonorm()) stag++;
    else                                               stag=0;

    x.copyFrom(xhalf); x.axpy( omega, sh); // x=xhalf+omega*sh
    r.copyFrom(s);     r.axpy(-omega, t ); // r=s-omega*t

    normr = r.twonorm(); normr_act = normr;
    resvec[2*ii+1] = normr;

    //printf("stag=%d  maxstagsteps=%d moresteps=%d  normr=%g\n",
    //	   stag, maxstagsteps, moresteps, normr);    

    //-------- check for convergence at the end of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      r.copyFrom(b); 
      //applyA(1.0, r, -1.0, x); //r=b-Ax
      SCmult(1.0,r, -1.0,x, prob);
      normr_act=r.twonorm();

      if(normr<=tolb) { flag = 0; iter = 1.0+ii; break; }
      else {
	if(stag>=maxstagsteps && moresteps==0) {
	  stag = 0;
	}
	moresteps++;
	if(moresteps>=maxmsteps) {
	  //method stagnated
	  flag=3; break;
	}
      }
    } // end convergence check
    if(stag>=maxstagsteps) { flag=3; break;} //stagnation

    //update quantities related to minimal norm iterate
    if(normr_act<normrmin) {
      xmin.copyFrom(x); normrmin=normr_act;
      imin=1.5+ii;
    }
    //printf("iter %g normr=%g\n", ii+1.0, normr);
    ///////////////////////////////
    // Next iterate
    ///////////////////////////////
    ii++;
    
  }//end while

  if(ii>=maxit) {
    iter=ii;
    flag=10;
  }
  
  if(flag==0 || flag==-1) {

    relres = normr_act/n2b;
#ifdef TIMING
    if(myRank==0) {
      printf("BiCGStab converged: normResid=%g relResid=%g iter=%g\n", 
	     normr_act, relres, iter);
      }
#endif
  } else {
    if(ii==maxit) flag=10;//aaa
    //FAILURE -> return minimum resid-norm iterate
    r.copyFrom(b); 
    //applyA(1.0, r, -1.0, xmin);
    SCmult(1.0,r, -1.0,xmin, prob);

    normr=r.twonorm();
    if(normr >= normr_act) {
      x.copyFrom(xmin);
      //iter=imin;
      relres=normr/n2b;
    } else {
      iter=1.0+ii;
      relres = normr/n2b;
    }
  
#ifdef TIMING
    if(myRank==0) {
      printf("BiCGStab did not NOT converged after %g[%d] iterations.\n", iter,ii);
      printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g %g\n\n", 
	     flag, normr, relres, normrmin);
    }
#endif
  }

  b.copyFrom(x);
  delete[] resvec;
}



void sLinsysRootAug::finalizeKKT(sData* prob, Variables* vars)
{
//  assert(locmz==0||gOuterSolve<3);
  assert(gOuterSolve>=3);

  int j, p, pend; double val;

  stochNode->resMon.recFactTmLocal_start();
  stochNode->resMon.recSchurMultLocal_start();

  DenseSymMatrix * kktd = (DenseSymMatrix*) kkt;
  //alias for internal buffer of kkt
  double** dKkt = kktd->Mat();
 

  //////////////////////////////////////////////////////
  // compute Q+diag(xdiag)  
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
  SimpleVector& syDiag = dynamic_cast<SimpleVector&>(*yDiag);

  for(int i=0; i<locnx; i++) dKkt[i][i] += sxDiag[i];  
  for(int i=locnx+locmz; i<locnx+locmz+locmy; i++) dKkt[i][i] += syDiag[i-locnx-locmz];

  SimpleVector& ssDiag = dynamic_cast<SimpleVector&>(*sDiag);
  SimpleVector& szDiag = dynamic_cast<SimpleVector&>(*zDiag);

  /////////////////////////////////////////////////////////////
  // update the KKT with  S part
  /////////////////////////////////////////////////////////////
  if(locmz>0) {
    for(int i=locnx; i<locnx+locmz; i++) {
        dKkt[i][i] += ssDiag[i-locnx];
    }
  } //~end if locmz>0
  /////////////////////////////////////////////////////////////
  // update the KKT with A (symmetric update forced)
  /////////////////////////////////////////////////////////////
  if(locmy>0){
    kktd->symAtPutSubmatrix( locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx, 1 );
  }
  /////////////////////////////////////////////////////////////
  // update the KKT with C (symmetric update forced) ,  -I and dual reg
  /////////////////////////////////////////////////////////////  
  if(locmz>0){
    kktd->symAtPutSubmatrix( locnx+locmz+locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx, 1 );
	for(int i=0; i<locmz; i++){
		dKkt[i+locnx+locmz+locmy][i+locnx] -= 1.0;
		dKkt[i+locnx][i+locnx+locmz+locmy] -= 1.0;
		dKkt[i+locnx+locmz+locmy][i+locnx+locmz+locmy] += szDiag[i];
	}
  }


  //prob->getLocalB().getStorageRef().dump("stage1eqmat2.dump");

  /////////////////////////////////////////////////////////////
  // update the KKT zeros for the lower right block 
  /////////////////////////////////////////////////////////////
  //kktd->storage().atPutZeros(locnx, locnx, locmy+locmz, locmy+locmz);
  //myAtPutZeros(kktd, locnx, locnx, locmy, locmy);

  stochNode->resMon.recSchurMultLocal_stop();
  stochNode->resMon.recFactTmLocal_stop();
}
