#include "CGSolver.h"
#include "SimpleVector.h"

#include <math.h>

extern int gOoqpPrintLevel;

#define EPS 2.220e-16

CGSolver::CGSolver( MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2)
  : DoubleIterativeLinearSolver(A, M1, M2)
{ 
  tol = 5.0e-13;
  maxit = 10; iter=-1;
  flag = -1;
  tmpVec1 = tmpVec2 = tmpVec3 = tmpVec4 = tmpVec5 = tmpVec6 =NULL;
  //firstSolve = 1;
};

CGSolver::~CGSolver()
{
  if(tmpVec1) delete[] tmpVec1;
  if(tmpVec2) delete[] tmpVec2;
  if(tmpVec3) delete[] tmpVec3;
  if(tmpVec4) delete[] tmpVec4;
  if(tmpVec5) delete[] tmpVec5;
  if(tmpVec6) delete[] tmpVec6;
}

void CGSolver::solve( OoqpVector& rhs_ )
{
  SimpleVector& b = dynamic_cast<SimpleVector&>(rhs_);
  int n = b.length();

  int flag, imin; int stag, maxmsteps, maxstagsteps, moresteps;
  double normr, normr_act, normrmin; 
  double alpha, beta, rho, rho1, pq; int iter;

  double n2b  = b.twonorm();
  double tolb = n2b*tol;

  if(tmpVec1==NULL) tmpVec1=new double[n];
  if(tmpVec2==NULL) tmpVec2=new double[n];
  if(tmpVec3==NULL) tmpVec3=new double[n];
  if(tmpVec4==NULL) tmpVec4=new double[n];
  if(tmpVec5==NULL) tmpVec5=new double[n];
  if(tmpVec6==NULL) tmpVec6=new double[n];

  SimpleVector x(tmpVec1, n);      //iterate
  SimpleVector r(tmpVec2,n);      //residual
  SimpleVector xmin(tmpVec3,n);   //minimal residual iterate
  SimpleVector y(tmpVec4,n);      //work vectors
  SimpleVector z(tmpVec5,n);      //work vectors
  SimpleVector p(tmpVec6,n);
  //if(firstSolve)
  //  //initial guess is 0, the previous found solution otherwise
  x.setToZero();

  xmin.copyFrom(x); y.setToZero();
  flag=1; imin=0;

  r.copyFrom(b); applyA(1.0, r, -1.0, x);
  normr=r.twonorm();

  maxit=n/2+20;
  if(normr<tolb) {
    printf("lucky!!!!!!!!!!!!\n");
    //initial guess is good enough
    b.copyFrom(x); return;
  }
  normrmin = normr; rho=1.0; stag=0;  
  maxmsteps=min(min(n/50, 5), n-maxit); 
  maxstagsteps=2; moresteps=0; iter=0;

  //////////////////////////////////////////////////////////////////
  // loop over maxit iterations
  //////////////////////////////////////////////////////////////////
  int ii=0; while(ii<maxit) {
    applyM1(0.0, y, 1.0, r);
    applyM2(0.0, z, 1.0, y);
    
    rho1=rho; rho=r.dotProductWith(z);
    if(rho==0.0) {flag=4; break;}

    if(ii==0) p.copyFrom(z);
    else {
      beta = rho/rho1;
      if(beta==0.0) {flag=4; break;}
      p.scale(beta); p.axpy(1.0, z); // p=z + beta*p
    }

    SimpleVector& q = y;
    applyA(0.0, q, 1.0, p); //q=A*p
    pq = p.dotProductWith(q);
    if(pq<=0) {flag=4; break;}
    alpha = rho/pq;
    
    //check for stagnation
    if(p.twonorm()*fabs(alpha) <EPS*x.twonorm()) stag++;
    else                                         stag=0;

    //---- updates ----
    x.axpy( alpha, p); 
    r.axpy(-alpha, q);
    normr=r.twonorm(); normr_act=normr;

    //printf("stag=%d  maxstagsteps=%d moresteps=%d  normr=%g\n",
    //   stag, maxstagsteps, moresteps, normr);
    // check for convergence
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      r.copyFrom(b); applyA(1.0, r, -1.0, x);
      normr_act=r.twonorm();

      if(normr_act<=tolb) { flag=0; iter=1+ii; break; }
      else {
	if(stag>=maxstagsteps && moresteps==0) stag=0;
	moresteps++;
	if(moresteps>=maxmsteps) {flag=3; iter=1+ii; break;}
      }
    }

    if(normr_act<normrmin) {normrmin=normr_act; xmin.copyFrom(x); imin=ii; }

    if(stag>=maxstagsteps) { flag=3; break; }

    ii++;
  }

  //////////////////////////////////////////////////////////
  // status/error output
  /////////////////////////////////////////////////////////
  if(flag==0) {
    double relres = normr_act/n2b;
    printf("CG converged: actual normResid=%g relResid=%g iter=%d\n", 
	   normr_act, relres, iter);
  } else {
    double relres = normr_act/n2b;
    r.copyFrom(b); applyA(1.0, r, -1.0, xmin);
    normr=r.twonorm();
    if(normr<normr_act) { x.copyFrom(xmin); iter=imin; relres=normr/n2b;}
    else                {iter=ii; relres=normr_act/n2b;}

    if(gOoqpPrintLevel>=1) {
      printf("CG did not NOT converged after %d  max of %d iters were made.\n", iter,ii);
      printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g %g\n\n", 
	     flag, normr, relres, normrmin);
    }

  }
  b.copyFrom(x);
}

