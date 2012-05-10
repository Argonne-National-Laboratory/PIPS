#include "BiCGStabSolver.h"
#include "SimpleVector.h"

#include <math.h>

extern int gOoqpPrintLevel;

#define EPS 2.220e-16

BiCGStabSolver::BiCGStabSolver( MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2)
  : DoubleIterativeLinearSolver(A, M1, M2)
{ 
  tol = 2.5e-15;
  maxit = 10; iter=-1;
  flag = -1;
};

void BiCGStabSolver::solve( OoqpVector& rhs_ )
{
  SimpleVector& b = dynamic_cast<SimpleVector&>(rhs_);
  int n = b.length();

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
  maxit = n/2+1;

  //////////////////////////////////////////////////////////////////
  //  Problem Setup and intialization
  //////////////////////////////////////////////////////////////////
  n2b = b.twonorm();
  tolb = n2b*tol;
  x.setToZero(); //initial guess

  //even r=b do the general case, maybe later the class will support user provided initial guesses.
  r.copyFrom(b); applyA(1.0, r, -1.0, x);
  //applyA(x,r); r.axpy(-1.0, b); r.negate(); // residual r=b-Ax
  normr=r.twonorm();

  if(normr<tolb) {
    //initial guess is good enough
    b.copyFrom(x); flag=0; printf("lucky!!!!!!!!!!!!\n");
    return;
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
    ///////////////////////////////
    // First half of the iterate
    ///////////////////////////////
    double rho1=rho; double beta;
    rho = rt.dotProductWith(r); 
    if(0.0==rho) { flag=4;  break; }

    if(ii==0) p.copyFrom(r);
    else {
      beta = (rho/rho1)*(alpha/omega);
      if(beta==0.0) { flag=4;  break; }

      //-------- p = r + beta*(p - omega*v) --------
      p.axpy(-omega, v); p.scale(beta); p.axpy(1.0, r);
    }

    //------ v = A*(M2inv*(M1inv*p)) and ph=M2inv*(M1inv*p)
    //first use v as temp storage
    applyM1(0.0, v,    1.0, p);
    applyM2(0.0, paux, 1.0, v);
    applyA (0.0, v,    1.0, paux); 
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
      applyA(1.0, s, -1.0, xhalf); // s=b-Ax
      //applyA(xhalf, s); s.negate(); s.axpy(1.0, b); // s=b-Ax
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

    ///////////////////////////////
    // Second half of the iterate
    //////////////////////////////

    applyM1(0.0, t,    1.0, s); //applyM1(s,     stemp);
    applyM2(0.0, paux, 1.0, t); //applyM2(stemp, sh);
    applyA (0.0, t,    1.0, paux); //applyA (sh, t);
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
      //applyA(x, r); r.negate(); r.axpy(1.0, b); // r=b-Ax
      r.copyFrom(b); applyA(1.0, r, -1.0, x);
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

  
  if(flag==0) {
    relres = normr_act/n2b;
    printf("BiCGStabConvergence: actual normResid=%g relResid=%g iter=%g\n", 
	   normr_act, relres, iter);
  } else {
    if(ii==maxit) flag=10;//aaa
    //FAILURE -> return minimum resid-norm iterate
    r.copyFrom(b); applyA(1.0, r, -1.0, xmin);
    //applyA(x,r); r.negate(); r.axpy(1.0, b);
    normr=r.twonorm();
    if(normr >= normr_act) {
      x.copyFrom(xmin);
      //iter=imin;
      relres=normr/n2b;
    } else {
      iter=1.0+ii;
      relres = normr/n2b;
    }
  
    if(gOoqpPrintLevel>=1) {
      printf("BiCGStab did not NOT converged after %g[%d] iterations.\n", iter,ii);
      printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g %g\n\n", 
	     flag, normr, relres, normrmin);
    }
  }

  b.copyFrom(x);
  delete[] resvec;
}



/*
void BiCGStabSolver::solve( OoqpVector& rhs_ )
{
  SimpleVector& b = dynamic_cast<SimpleVector&>(rhs_);
  int n = b.length();

  SimpleVector r(n);           //residual
  SimpleVector s(n);           //residual associated with half iterate
  SimpleVector rt(n);          //shadow residual
  SimpleVector xmin(n);        //minimal residual iterate
  SimpleVector x(n);           //iterate
  SimpleVector xhalf(n);       // half iterate of BiCG
  SimpleVector p(n), ph(n), ptemp(n);
  SimpleVector v(n);
  int flag, imin;
  double n2b;                  //norm of b 
  double normr, normrmin;      //norm of the residual and norm of residual at min-resid iterate
  double normr_act;
  double tolb;                 //relative tolerance
  double rho, omega, alpha;
  int stag, maxmsteps, maxstagsteps, moresteps;
  double relres;
  maxit = n/2+1;

  //////////////////////////////////////////////////////////////////
  //  Problem Setup and intialization
  //////////////////////////////////////////////////////////////////
  n2b = b.twonorm();
  tolb = n2b*tol;
  x.setToZero(); //initial guess
  applyA(x,r); r.axpy(-1.0, b); r.negate(); // residual r=b-Ax
  normr=r.twonorm();

  if(normr<tolb) {
    //initial guess is good enough
    b.copyFrom(x);
    return;
  }

  rt.copyFrom(r); //Shadow residual
  double* resvec = new double[2*maxit+1];
  resvec[0] = normr; normrmin=normr;
  rho=1.0; omega=1.0;
  stag=0; maxmsteps=min(min(n/50, 10), n-maxit); 
  maxstagsteps=3; moresteps=0;
  
  //////////////////////////////////////////////////////////////////
  // loop over maxit iterations
  //////////////////////////////////////////////////////////////////
  int ii=0; while(ii<maxit) {
    ///////////////////////////////
    // First half of the iterate
    ///////////////////////////////
    double rho1=rho; double beta;
    rho = rt.dotProductWith(r); 
    if(0.0==rho) { flag=4;  break; }

    if(ii==0) p.copyFrom(r);
    else {
      beta = (rho/rho1)*(alpha/omega);
      if(beta==0.0) { flag=4;  break; }
      //-------- p = r + beta*(p - omega*v) --------
      p.axpy(-omega, v); p.scale(beta); p.axpy(1.0, r);
    }

    ptemp.setToZero(); //! do we need this?
    ph   .setToZero(); //! do we need this?
    applyM1(p,     ptemp);
    applyM2(ptemp, ph);
    applyA(ph,     v);
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

    //-------- check for convergence in the middle of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      applyA(xhalf, s); s.negate(); s.axpy(1.0, b); // s=b-Ax
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

    ///////////////////////////////
    // Second half of the iterate
    ///////////////////////////////
    //aliasses 
    SimpleVector& sh = ph; SimpleVector& stemp = ptemp; SimpleVector& t=v;

    applyM1(s,     stemp);
    applyM2(stemp, sh);
    applyA (sh, t);

    double tt = t.dotProductWith(t);
    if(tt==0.0) { flag=4; break;}

    omega=t.dotProductWith(s); omega /= tt;

    if(fabs(omega)*sh.twonorm() < EPS*xhalf.twonorm()) stag++;
    else                                              stag=0;

    x.copyFrom(xhalf); x.axpy( omega, sh); // x=xhalf+omega*sh
    r.copyFrom(s);     r.axpy(-omega, t ); // r=s-omega*t

    normr = r.twonorm(); normr_act = normr;
    resvec[2*ii+1] = normr;
    
    //-------- check for convergence in the middle of the iterate.  -------- 
    if(normr<=tolb || stag>=maxstagsteps || moresteps) {
      applyA(x, r); r.negate(); r.axpy(1.0, b); // r=b-Ax
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

    ///////////////////////////////
    // Next iterate
    ///////////////////////////////
    ii++;
  }//end while

  
  if(flag==0) {
    relres = normr_act/n2b;
  } else {
    //FAILURE -> return minimum resid-norm iterate
    applyA(x,r); r.negate(); r.axpy(1.0, b);
    normr=r.twonorm();
    if(normr <= normr_act) {
      x.copyFrom(xmin);
      iter=imin;
      relres=normr/n2b;
    } else {
      iter=1.0+ii;
      relres = normr_act/n2b;
    }
  
    if(gOoqpPrintLevel>=1) printf("BiCGStab did NOT converged after %d iterations. Error code %d",
				  iter, flag);
  }

  b.copyFrom(x);
  delete[] resvec;
}

*/
