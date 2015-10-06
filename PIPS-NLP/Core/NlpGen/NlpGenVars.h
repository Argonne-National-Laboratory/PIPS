/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef NLPGENVARS
#define NLPGENVARS

#include "Variables.h"
#include "OoqpVectorHandle.h"

class NlpGen;
class NlpGenData;
class LinearAlgebraPackage;
class MpsReader;

#ifdef TESTING
class NlpGenVarsTester;
#endif

/**
 * Variables for the general NLP formulation.
 * 
 * @ingroup NlpGen
 */
class NlpGenVars : public Variables {
#ifdef TESTING
  friend NlpGenVarsTester;
#endif
//protected:



public:

  long long nx, nxupp, nxlow;
  long long my;
  long long mz, mcupp, mclow;

  OoqpVectorHandle ixlow;
  OoqpVectorHandle ixupp;
  OoqpVectorHandle icupp;
  OoqpVectorHandle iclow;

  int nSlack;
  int nxLOri, nxUOri;
  int nsLOri, nsUOri;
  
  NlpGenVars(){};

//public:  

  OoqpVectorHandle x;
  OoqpVectorHandle s;
  OoqpVectorHandle y;
  OoqpVectorHandle z;

  OoqpVectorHandle v;  
  OoqpVectorHandle gamma;

  OoqpVectorHandle w;
  OoqpVectorHandle phi;


  OoqpVectorHandle t;
  OoqpVectorHandle lambda;
  
  OoqpVectorHandle u;
  OoqpVectorHandle pi;

  OoqpVectorHandle PiLx;
  OoqpVectorHandle PiUx;

  OoqpVectorHandle PiLs;
  OoqpVectorHandle PiUs;


  /** constructor in which the data and variable pointers are set to
      point to the given arguments */
  NlpGenVars( OoqpVector * x_in, OoqpVector * s_in,
	     OoqpVector * y_in, OoqpVector * z_in,
	     OoqpVector * v_in, OoqpVector * gamma_in,
	     OoqpVector * w_in, OoqpVector * phi_in,
	     OoqpVector * t_in, OoqpVector * lambda_in,
	     OoqpVector * u_in, OoqpVector * pi_in,
	     OoqpVector * ixlow_in, OoqpVector * ixupp_in,
	     OoqpVector * iclow_in, OoqpVector * icupp_in );

  /** constructor that creates variables objects of specified
      dimensions. */
  NlpGenVars( LinearAlgebraPackage * la,
	     int nx_, int my_, int mz_,
	     OoqpVector * ixlow, OoqpVector * ixupp,
	     OoqpVector * iclow, OoqpVector * icupp );


  virtual ~NlpGenVars();
  
  /** computes mu = (t'lambda +u'pi + v'gamma + w'phi)/(mclow+mcupp+nxlow+nxupp) */
  virtual double mu();

  virtual double mustep(Variables *step_in, double alpha);

  virtual void saxpy( Variables *b, double alpha );
  virtual void negate();
  
  /** calculate the largest alpha in (0,1] such that the nonnegative
   * variables stay nonnegative in the given search direction. In the
   * general NLP problem formulation, this is the largest value of
   * alpha such that (t,u,v,w,lambda,pi,phi,gamma) + alpha *
   * (b->t,b->u,b->v,b->w,b->lambda,b->pi,b->phi,b->gamma) >= 0.
   *
   * @see findBlocking */
  virtual double stepbound( Variables *b );

  /** Performs the same function as stepbound, and supplies additional
   * information about which component of the nonnegative variables is
   * responsible for restricting alpha. In terms of the abstract
   * formulation, the components have the following meanings.
   *
   * @param primalValue the value of the blocking component of the
   * primal variables (u,t,v,w).
   * 
   * @param primalStep the corresponding value of the blocking
   * component of the primal step variables (b->u,b->t,b->v,b->w).
   * 
   * @param dualValue the value of the blocking component of the dual
   * variables (lambda,pi,phi,gamma).
   *
   * @param dualStep the corresponding value of the blocking component
   * of the dual step variables (b->lambda,b->pi,b->phi,b->gamma).
   * 
   * @param firstOrSecond  1 if the primal step is blocking, 2 if the dual
   * step is block, 0 if no step is blocking.  
   *
   * @see stepbound
   * */
  virtual double findBlocking( Variables * step, 
			       double & primalValue,
			       double & primalStep,
			       double & dualValue,
			       double & dualStep,
			       int& firstOrSecond );


  virtual void findBlockingPriDual( Variables * step, 
				  double & primalValue,
				  double & primalStep,
				  double & dualValue,
				  double & dualStep,
				  int& firstOrSecond, 
				  double tau, double & alphaPri, double & alphaDual );

  

  /** sets components of (u,t,v,w) to alpha and of
      (lambda,pi,phi,gamma) to beta */
  virtual void interiorPoint( double alpha, double beta );
  

  /** add alpha to components of (u,t,v,w) and beta to components of
      (lambda,pi,phi,gamma) */
  virtual void shiftBoundVariables( double alpha, double beta );

  /** check whether this is an interior point. Useful as a sanity check. */
  virtual int isInteriorPoint();

  virtual double violation();

  virtual void print();
  virtual void printSolution( MpsReader * reader, NlpGenData * prob,
			      int& iErr );

  virtual void unscaleSolution( NlpGenData * data);
  virtual void unscaleBounds  ( NlpGenData * data);

  virtual int  validNonZeroPattern();
  
  virtual void copy(Variables *b);

  virtual double onenorm();
  virtual double infnorm();



  /** copy dual variables */
  virtual void copyDualPri0(Variables *b);

  virtual void setPriZero();

  virtual void setZero();








  

  virtual void interiorPointPriX( double alpha);
  virtual void interiorPointPriS( double alpha);
  virtual void interiorBoundSlack( double alpha);  
  virtual void interiorBoundSlackDual( double alpha);  

  
  virtual void interiorPointDualY( double alpha);

  virtual void interiorPointDualZ( double alpha);








/** 
   compute max step size by bound-frac rule, split pri and dual steps
*/
virtual double stepMax_Pri( Variables * b_in, const double tau );
virtual double stepMax_BoundDual( Variables * b_in, const double tau );



virtual void takeStep( Variables *step_in, 
				const double alphaPri, const double alphaDual, const double alphaNu, const int onlyPrimal=0 );

virtual void SetSlackFromMaxXorY( OoqpVector* v_out, 
				OoqpVector* x_in, OoqpVector *y_in);


virtual void takePrimalStep( Variables *step_in, const double alphaPri, double alphaSlack=-1);
virtual void takeDualStep( Variables *step_in, const double alphaPri, double alphaSlackDual=-1);


virtual void updateSlackAndDual( OoqpVector *tempx,OoqpVector *temps,const double k_sigma, const double mu);

virtual double primal_XS_InfNorm(bool XSonly = true);

virtual double dual_YZ_InfNorm(bool YZonly = true);

virtual void push_variables( OoqpVector *vec, OoqpVector *vec_slackLB, OoqpVector *vec_slackUB,
					OoqpVector *vec_Lb,OoqpVector *vec_Ub, 
					OoqpVector *vec_Temp, const double k_1, const double k_2, const int ifX);

virtual void getErrScaling(double s_max, double &scal_comm, double &scal_dualerr);


virtual double computeDD();
virtual double computeXSDD(Variables * step_in);


virtual void mergeNTstep(Variables * Nstep, Variables * Tstep, Variables * Curr_Iter, Residuals * res_in);


  
};

/** Indicates what type is the blocking variable in the step length
 * determination. If tblock, then the blocking variable is one of the
 * slack variables t for a general lower bound, and so on. Special
 * value no_block is for the case in which a step length of 1 can be
 * taken without hitting the bound.  */

enum { no_block = 0,
       t_block = 1,
       lambda_block = 2,
       u_block = 3,
       pi_block = 4,
       v_block = 5,
       gamma_block = 6,
       w_block = 7,
       phi_block = 8
};  


#endif

