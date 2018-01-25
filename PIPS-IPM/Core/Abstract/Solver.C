/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "Solver.h"
#include "OoqpMonitor.h"
#include "Status.h"
#include "SimpleVector.h"
#include "Data.h"
#include "Variables.h"
#include "Residuals.h"
#include "LinearSystem.h"
#include "OoqpStartStrategy.h"

#include <cmath>

//#define BAD_NUMERICS

int gOoqpPrintLevel = 1000;
int gLackOfAccuracy=0;
int onSafeSolver=0;

int gOuterBiCGIter=0; //use for dumping matrices

//controls the type of error absortion at the outer level
// - 0:no error absortion (OOQP works just fine)
// - 1:iterative refinement (used when error absortion is 
// also done at a lower level, for example in the solve with 
// the dense Schur complement
// - 2:BiCGStab with the factorization as preconditioner
int gOuterSolve=0; 

//controls the type of error absortion/correction when solving
//with the dense Schur complement
// - 0: no error correction
// - 1: iter. refin.
// - 2: BiCGStab
int gInnerSCsolve=0;

//number of iterative refinements in the 2nd stage sparse systems
int gInnerStg2solve=3;

Solver::Solver() : itsMonitors(0), status(0), startStrategy(0), dnorm(0.0),
		   mutol(1.e-6), artol(1e-4), phi(0.0), maxit(0), mu_history(0), rnorm_history(0),
		   phi_history(0), phi_min_history(0), iter(0), sys(0)
{
  // define parameters associated with the step length heuristic
#ifdef BAD_NUMERICS
  gamma_f = 0.95;
#else
  gamma_f = 0.99;
#endif
  gamma_a = 1.0 / (1.0 - gamma_f);

}

void Solver::start( ProblemFormulation * formulation,
		    Variables * iterate, Data * prob,
		    Residuals * resid, Variables * step  )
{
  if( startStrategy ) {
    startStrategy->doIt( this, formulation, iterate, prob, resid, step );
  } else {
    this->defaultStart( formulation, iterate, prob, resid, step );
     //this->stevestart( formulation, iterate, prob, resid, step );
    //this->dumbstart( formulation, iterate, prob, resid, step );
  }
}
//#ifdef TIMING
//#include "mpi.h"
//#endif

void Solver::defaultStart( ProblemFormulation * /* formulation */,
			   Variables * iterate, Data * prob,
			   Residuals * resid, Variables * step  )
{
  double sdatanorm = sqrt(dnorm);
  double a  = sdatanorm;
  double b  = sdatanorm;
  iterate->interiorPoint( a, b );

  resid->calcresids( prob, iterate );
  resid->set_r3_xz_alpha( iterate, 0.0 );

  sys->factor( prob, iterate );
  sys->solve( prob, iterate, resid, step );

  step->negate();
 
  // Take the full affine scaling step
  iterate->saxpy( step, 1.0 );
  // resid->calcresids(prob, iterate); // Calc the resids if debugging.
  double shift = 1.e3 + 2*iterate->violation();
  iterate->shiftBoundVariables( shift, shift );
}

void Solver::stevestart(  ProblemFormulation * /* formulation */,
			  Variables * iterate, Data * prob,
			  Residuals * resid, Variables * step  )
{
  double sdatanorm = sqrt(dnorm);
  double a = 0.0, b = 0.0;  

  iterate->interiorPoint( a, b );

  // set the r3 component of the rhs to -(norm of data), and calculate
  // the residuals that are obtained when all values are zero.

  resid->set_r3_xz_alpha( iterate, -sdatanorm );
  resid->calcresids( prob, iterate );

  // next, assign 1 to all the complementary variables, so that there
  // are identities in the coefficient matrix when we do the solve.

  a = 1.0; b = 1.0;
  iterate->interiorPoint( a, b );
  sys->factor(prob, iterate);
  sys->solve (prob, iterate, resid, step);
  step->negate();

  // copy the "step" into the current vector

  iterate->copy(step);

  // calculate the maximum violation of the complementarity
  // conditions, and shift these variables to restore positivity.
  double shift;
  shift = 1.5 * iterate->violation();
  iterate->shiftBoundVariables( shift, shift );

  // do Mehrotra-type adjustment

  double mutemp = iterate->mu();
  double xsnorm = 0.0;
  xsnorm = iterate->onenorm();
  double delta = 0.5 * iterate->nComplementaryVariables * mutemp / xsnorm;
  iterate->shiftBoundVariables( delta, delta );
}

void Solver::dumbstart(  ProblemFormulation * /* formulation */,
			 Variables * iterate, Data * /* prob */,
			 Residuals * /*resid*/, Variables * /*step*/  )
{
  double sdatanorm = dnorm;
  double a, b, bigstart;
  a = 1.e3; b = 1.e5;
  bigstart = a*sdatanorm + b;
  iterate->interiorPoint( bigstart, bigstart );
}

double Solver::finalStepLength( Variables *iterate, Variables *step )
{
	double primalValue, primalStep, dualValue, dualStep,
	  maxAlpha, mufull;
	int firstOrSecond;

	maxAlpha = iterate->findBlocking( step, 
					  primalValue, primalStep,
					  dualValue, dualStep,
					  firstOrSecond );
	mufull = iterate->mustep( step, maxAlpha );
	mufull /= gamma_a;

	double alpha = 1.0;
	switch( firstOrSecond ) {
	case 0:
	  alpha = 1; // No constraints were blocking
	  break;
	case 1:
	  alpha = ( - primalValue +
		    mufull / ( dualValue + maxAlpha * dualStep ) ) /
	    primalStep;
#ifdef TIMING
	  std::cout << "alpha " << primalValue + maxAlpha * primalStep << std::endl;
	  //assert(primalValue + maxAlpha * primalStep >= 0.0);
#endif
	  break;
	case 2:
	  alpha = ( - dualValue +
		    mufull / ( primalValue + maxAlpha * primalStep ) ) /
	    dualStep;
#ifdef TIMING
	  std::cout << "dual alpha " << dualValue + maxAlpha * dualStep << std::endl;
	  //assert(dualValue + maxAlpha * dualStep >= 0.0);
#endif
	  break;
	default:
	  cout << "Can't get here: firstOrSecond=" << firstOrSecond << endl;
	  assert( 0 && "Can't get here" );
          break;
	}

	// make it at least gamma_f * maxStep
	if( alpha < gamma_f * maxAlpha ) alpha = gamma_f * maxAlpha;

	// back off just a touch (or a bit more)
#ifdef BAD_NUMERICS
	alpha *= 0.995;
#else
	alpha *= .99999999;
#endif

	assert(alpha < 1.0);


	return alpha;
}

void Solver::finalStepLength_PD( Variables *iterate, Variables *step,
		  	  	  	  	  	  	  double& alpha_primal, double& alpha_dual )
{
	double primalValue_p, primalStep_p, dualValue_p, dualStep_p,
		  maxAlpha_p;
	double primalValue_d, primalStep_d, dualValue_d, dualStep_d,
		  maxAlpha_d;
	double mufull;
	bool primalBlocking, dualBlocking;

	iterate->findBlocking_pd( step,
					  primalValue_p, primalStep_p, dualValue_p, dualStep_p,
					  primalValue_d, primalStep_d, dualValue_d, dualStep_d,
					  maxAlpha_p, maxAlpha_d,
					  primalBlocking, dualBlocking );
	assert( primalBlocking || !primalBlocking );
	assert( dualBlocking || !dualBlocking );

	mufull = iterate->mustep_pd( step, maxAlpha_p, maxAlpha_d);
	mufull /= gamma_a;

	alpha_primal = 1.0;
	if( !primalBlocking )
	{
		alpha_primal = 1; // No primal constraints were blocking
	}
	else if( primalBlocking )
	{
		alpha_primal = ( - primalValue_p +
						mufull / ( dualValue_p + maxAlpha_d * dualStep_p ) ) /
					primalStep_p;
		#ifdef TIMING
			std::cout << "primal alpha " << primalValue_p + maxAlpha_p * primalStep_p << std::endl;
		#endif
	}

	alpha_dual = 1.0;
	if( !dualBlocking )
	{
		alpha_dual = 1; // No dual constraints were blocking
	}
	else if( dualBlocking )
	{
		alpha_dual = ( - dualValue_d +
						mufull / ( primalValue_d + maxAlpha_p * primalStep_d ) ) /
					dualStep_d;
		#ifdef TIMING
			std::cout << "dual alpha " << dualValue_d + maxAlpha_d * dualStep_d << std::endl;
		#endif
	}

	// make it at least gamma_f * maxAlpha and no bigger than 1
	if( alpha_primal < gamma_f * maxAlpha_p || alpha_primal > 1 ) alpha_primal = gamma_f * maxAlpha_p;
	if( alpha_dual < gamma_f * maxAlpha_d || alpha_dual > 1 ) alpha_dual = gamma_f * maxAlpha_d;

	// back off just a touch (or a bit more)
	#ifdef BAD_NUMERICS
		alpha_primal *= 0.995;
		alpha_dual *= 0.995;
	#else
		alpha_primal *= .99999999;
		alpha_dual *= .99999999;
	#endif

	assert(alpha_primal < 1.0 && alpha_dual < 1.0);
	assert(alpha_primal >= 0 && alpha_dual >= 0 );
}


void Solver::doMonitor( Data * data, Variables * vars,
			Residuals * resids,
			double alpha, double sigma,
			int i, double mu,
                        int stop_code,
			int level )
{
  OoqpMonitor * m = itsMonitors;

  while( m ) {
    m->doIt( this, data, vars, resids, alpha, sigma, i, mu, stop_code, level );
    m = m->nextMonitor;
  }
}

int Solver::doStatus( Data * data, Variables * vars,
		       Residuals * resids,
		       int i, double mu,
		       int level )
{
  if( status ) {
    return status->doIt( this, data, vars, resids, i, mu, level );
  } else {
    return this->defaultStatus( data, vars, resids, i, mu, level );
  }
}


void Solver::monitorSelf()
{
  this->addMonitor( new OoqpSelfMonitor );
}

void Solver::addMonitor( OoqpMonitor * m )
{
  // Push the monitor onto the list
  m->nextMonitor = itsMonitors;
  itsMonitors = m;
}

Solver::~Solver()
{
  OoqpMonitor * m = itsMonitors;
  while( m ) {
    OoqpMonitor * n = m->nextMonitor;
    delete m;
    m = n;
  }
}


int Solver::defaultStatus(Data * /* data */, Variables * /* vars */,
				 Residuals * resids,
				 int iterate, double mu, 
				 int /* level */)
{
  int stop_code = NOT_FINISHED;
  int idx;

  double gap   = fabs( resids->dualityGap() );
  double rnorm = resids->residualNorm();

  idx = iterate-1;
  if(idx <  0     ) idx=0;
  if(idx >= maxit ) idx=maxit-1;

  // store the historical record
  mu_history[idx] = mu;
  rnorm_history[idx] = rnorm;
  phi = (rnorm + gap) / dnorm;
  phi_history[idx] = phi;

#ifdef TIMING
  std::cout << "mu/mutol " << mu << "  " << mutol << "rnorm/limit " << rnorm << " " << artol*dnorm  << std::endl;
#endif

  if(idx > 0) {
    phi_min_history[idx] = phi_min_history[idx-1];
    if(phi < phi_min_history[idx]) phi_min_history[idx] = phi;
  } else
    phi_min_history[idx] = phi;

  if ( iterate >= maxit ) {
    stop_code = MAX_ITS_EXCEEDED;
  } else if ( mu <= mutol && rnorm <= artol*dnorm ) {
    stop_code = SUCCESSFUL_TERMINATION;
  }
  //if(iterate>=50) stop_code = SUCCESSFUL_TERMINATION;
  if(stop_code != NOT_FINISHED)  return stop_code;

  // check infeasibility condition
  if(idx >= 10 && phi >= 1.e-8 && phi >= 1.e4*phi_min_history[idx]) {
#ifdef TIMING
    std::cout << "possible INFEASIBLITY detected, phi: " << phi << std::endl;
#endif
    stop_code = INFEASIBLE;
  }
  if(stop_code != NOT_FINISHED)  return stop_code;

  // check for unknown status: slow convergence first
  if(idx >= 300 && phi_min_history[idx] >= .5 * phi_min_history[idx-30]) {
    stop_code = UNKNOWN;
    printf("hehe dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm, artol);
  }

  if(idx >= 300 && rnorm / dnorm > artol && 
     (rnorm_history[idx]/mu_history[idx]) / (rnorm_history[0]/mu_history[0]) 
     >= 1.e8) {
    stop_code = UNKNOWN;
    printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm, artol);
  }

  //if(mu<50*rnorm/dnorm || mu<1e-5) {
  if(mu<1.0e5*rnorm/dnorm) {
    //if(!onSafeSolver) {
    gLackOfAccuracy=1;
    //cout << "Lack of accuracy detected ---->" << mu << ":" << rnorm/dnorm << endl;
  } else {
      //if(mu>1e7*rnorm/dnorm && mu>1.0e4)
      //gLackOfAccuracy=-1;
      //else
      gLackOfAccuracy=1;
  }
  //onSafeSolver=1;
  //}  
  //gLackOfAccuracy=-1; //disable iter refin in sLinsysRootAug
  
  return stop_code;
}


