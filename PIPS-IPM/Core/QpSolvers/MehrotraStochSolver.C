#include "MehrotraStochSolver.h"
#include "Variables.h"
#include "Residuals.h"
#include "LinearSystem.h"
#include "Status.h"
#include "Data.h"
#include "ProblemFormulation.h"

#include "sVars.h"
#include "sData.h"

#include "OoqpVector.h"
#include "DoubleMatrix.h"

#include "StochTree.h"
#include "QpGenStoch.h"
#include "StochResourcesMonitor.h"

#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;

#include <cstdio>
#include <cassert>
#include <cmath>

// gmu is used in EmtlStochSymIndefSolver to decide if the system should be perturbed
double gmu;
// double grnorm;
extern int gOoqpPrintLevel;
extern bool ipStartFound;
extern double g_iterNumber;

int sleepFlag=0;

MehrotraStochSolver::MehrotraStochSolver( ProblemFormulation * opt, Data * prob, const Scaler* scaler )
  : MehrotraSolver(opt, prob, scaler)
{

}
#include "StochVector.h"
#include "mpi.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"

int MehrotraStochSolver::solve(Data *prob, Variables *iterate, Residuals * resid )
{
  int done;
  double mu, alpha = 1, sigma = 1, muaff;
  int status_code;
  QpGenStoch* stochFactory = reinterpret_cast<QpGenStoch*>(factory);
  gmu = 1000;
  //  grnorm = 1000;
  setDnorm( *prob );
  // initialization of (x,y,z) and factorization routine.
  sys = factory->makeLinsys( prob );

  g_iterNumber=0.0;

  stochFactory->iterateStarted();
  this->start( factory, iterate, prob, resid, step );
  stochFactory->iterateEnded();

  assert(!ipStartFound);
  ipStartFound = true;
  iter = 0;
  done = 0;
  mu = iterate->mu();
  gmu = mu;

  do {

    iter++; g_iterNumber=iter;
    stochFactory->iterateStarted();
    
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
    resid->set_r3_xz_alpha(iterate, 0.0 );  //   resid->print();

    sys->factor(prob, iterate);
    sys->solve(prob, iterate, resid, step);
    step->negate();

    alpha = iterate->stepbound(step);

    // calculate centering parameter 
    muaff = iterate->mustep(step, alpha);
    sigma = pow(muaff/mu, tsig);
    
    assert(sigma<=1);
    g_iterNumber+=0.5;
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

    //assert(mu <= gmu);

    gmu = mu;
    
    stochFactory->iterateEnded();

  } while(!done);
  
  resid->calcresids(prob,iterate);
  if( gOoqpPrintLevel >= 10 ) {
    this->doMonitor( prob, iterate, resid,
		     alpha, sigma, iter, mu, status_code, 1 );
  }

  return status_code;
}


MehrotraStochSolver::~MehrotraStochSolver()
{}


