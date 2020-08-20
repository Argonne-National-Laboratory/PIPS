/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "MehrotraSolver.h"
#include "Variables.h"
#include "Residuals.h"
#include "LinearSystem.h"
#include "Status.h"
#include "Data.h"
#include "ProblemFormulation.h"

#include "OoqpVector.h"
#include "DoubleMatrix.h"

#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;

#include <cstdio>
#include <cassert>
#include <cmath>


#if NO_BUILT_IN_SUPPORT_FOR_BOOL
typedef int bool;
const bool true = 1;
const bool false = 0;
#endif

// gmu is needed by MA57!
static double gmu;
// double grnorm;
extern int gOoqpPrintLevel;

MehrotraSolver::MehrotraSolver( ProblemFormulation * opt, Data * prob, const Scaler * scaler ) :
      Solver( scaler )
{
  factory = opt;
  sys = 0;
  step = factory->makeVariables( prob );
  
  maxit      = 300;
  printlevel = 0;   // has no meaning right now 
  tsig       = 3.0; // the usual value for the centering exponent (tau)

  // allocate space to track the sequence of complementarity gaps,
  // residual norms, and merit functions.
  mu_history = new double[maxit];
  rnorm_history = new double[maxit];
  phi_history = new double[maxit];
  phi_min_history = new double[maxit];

  // Use the defaultStatus method
  status   = 0; 
}


int MehrotraSolver::solve(Data *prob, Variables *iterate, Residuals * resid )
{
  int done;
  double mu, alpha = 1, sigma = 1, muaff;
  int status_code;

  gmu = 1000;
  //  grnorm = 1000;
  setDnorm( *prob );

  // initialization of (x,y,z) and factorization routine.
  sys = factory->makeLinsys( prob );
  this->start( factory, iterate, prob, resid, step );

  iter = 0;
  done = 0;
  mu = iterate->mu();
  gmu = mu;

  do
    {
      iter ++;

      // evaluate residuals and update algorithm status:
      resid->calcresids(prob, iterate);
      
      // termination test:
      status_code = this->doStatus( prob, iterate, resid, iter, mu, 0 );
      if( status_code != NOT_FINISHED ) break;
      if( gOoqpPrintLevel >= 10 ) {
	this->doMonitor( prob, iterate, resid,
			 alpha, sigma, iter, mu, status_code, 0 );
      }

      // *** Predictor step ***

      resid->set_r3_xz_alpha(iterate, 0.0 );

      sys->factor(prob, iterate);
      sys->solve(prob, iterate, resid, step);

      step->negate();

      alpha = iterate->stepbound(step);

      // calculate centering parameter 
      muaff = iterate->mustep(step, alpha);
      sigma = pow(muaff/mu, tsig);

      assert(sigma<=1);

      // *** Corrector step ***

      // form right hand side of linear system:
      resid->add_r3_xz_alpha( step, -sigma*mu );

      sys->solve(prob, iterate, resid, step);
      step->negate();

 
      // We've finally decided on a step direction, now calculate the
      // length using Mehrotra's heuristic.
      alpha = finalStepLength(iterate, step);

      // alternatively, just use a crude step scaling factor.
      //alpha = 0.995 * iterate->stepbound( step );

      // actually take the step and calculate the new mu
      iterate->saxpy(step, alpha);
      mu = iterate->mu();
      gmu = mu;
    } while(!done);
  
  resid->calcresids(prob,iterate);
  if( gOoqpPrintLevel >= 10 ) {
    this->doMonitor( prob, iterate, resid,
		     alpha, sigma, iter, mu, status_code, 1 );
  }

  return status_code;
}

void MehrotraSolver::defaultMonitor( const Data * /* data */, const Variables * /* vars */,
				     const Residuals * resids,
				     double alpha, double /* sigma */,
				     int i, double mu, 
				     int statusCode,
				     int level ) const
{
  switch( level ) {
  case 0 : case 1: { 

    const Residuals* resids_unscaled = resids;
    if( scaler )
       resids_unscaled = scaler->getResidualsUnscaled(*resids);

    const double resid_norm = resids_unscaled->residualNorm();
    const double resids_gap = resids_unscaled->dualityGap();
    if( scaler )
       delete resids_unscaled;

    std::cout << " *** Iteration " << i << " *** " << std::endl;
    printf(" mu = %16.12e  relative residual norm = %16.12e\n", mu, resid_norm / dnorm_orig);

    std::cout << " Duality Gap: " << resids_gap << std::endl;
    if( i > 1 )
      std::cout << " alpha = " << alpha << std::endl;

    std::cout << endl;

    if( level == 1) { 
      // Termination has been detected by the status check; print
      // appropriate message
      switch( statusCode ) {
      case SUCCESSFUL_TERMINATION:
	std::cout << std::endl << " *** SUCCESSFUL TERMINATION ***" << std::endl;
	break;
      case MAX_ITS_EXCEEDED:
	std::cout << std::endl << " *** MAXIMUM ITERATIONS REACHED *** " << std::endl;
	break;
      case INFEASIBLE:
	std::cout << std::endl << " *** TERMINATION: PROBABLY INFEASIBLE *** " << std::endl;
      case UNKNOWN:
	std::cout << std::endl << " *** TERMINATION: STATUS UNKNOWN *** " << std::endl;
	break;
      } // end switch(statusCode)
    }
  } break; // end case 0: case 1:
  } // end switch(level)
}


MehrotraSolver::~MehrotraSolver()
{
  delete [] mu_history;
  delete [] rnorm_history;
  delete [] phi_history;
  delete [] phi_min_history;

  delete step;
  delete sys;
}


