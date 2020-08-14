/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "pipsport.h"
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
#include <limits>

#include "mpi.h"

bool ipStartFound = false;
double g_iterNumber = 0.0;
int gOoqpPrintLevel = 1000;
int gLackOfAccuracy=0;
int onSafeSolver=0;

int gOuterBiCGFails=0;
int gOuterBiCGIter=0;
int gInnerBiCGIter=0;
int gInnerBiCGFails=0;

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
#ifdef STEPLENGTH_CONSERVATIVE
  gamma_f = 0.95;
#else
  gamma_f = 0.99;
#endif
  gamma_a = 1.0 / (1.0 - gamma_f);


  printTimeStamp = false;
  startTime = 0.0;

  // todo proper parameter
  char* var = getenv("PARDISO_PRINT_TIMESTAMP");
  if( var != nullptr )
  {
     int use;
     sscanf(var, "%d", &use);
     if( use == 1 )
     {
        printTimeStamp = true;
        startTime = MPI_Wtime();
     }
  }
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
   double primalValue = -std::numeric_limits<double>::max();
   double primalStep = -std::numeric_limits<double>::max();
   double dualValue = -std::numeric_limits<double>::max();
   double dualStep = -std::numeric_limits<double>::max();
   int firstOrSecond = -1;


#ifdef TIMING
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

	const double maxAlpha = iterate->findBlocking( step,
					  primalValue, primalStep,
					  dualValue, dualStep,
					  firstOrSecond );

	const double mufull = iterate->mustep( step, maxAlpha ) / gamma_a;

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
	  if( myrank == 0 )
	     std::cout << "(primal) original alpha " << alpha << std::endl;
#endif
	  break;
	case 2:
	  alpha = ( - dualValue +
		    mufull / ( primalValue + maxAlpha * primalStep ) ) /
	    dualStep;
#ifdef TIMING
	  if( myrank == 0 )
	     std::cout << "(dual) original alpha " << alpha << std::endl;
#endif
	  break;
	default:
	  std::cout << "Can't get here: firstOrSecond=" << firstOrSecond << std::endl;
	  assert( 0 && "Can't get here" );
          break;
	}
   // safeguard against numerical troubles in the above computations
	alpha = std::min( maxAlpha, alpha );

	// make it at least gamma_f * maxStep
	if( alpha < gamma_f * maxAlpha ) alpha = gamma_f * maxAlpha;

	// back off just a touch (or a bit more)
#ifdef STEPLENGTH_CONSERVATIVE
	alpha *= 0.99;
#else
	alpha *= .99999999;
#endif

	assert(alpha < 1.0);


	return alpha;
}

void Solver::finalStepLength_PD( Variables *iterate, Variables *step,
		  	  	  	  	  	  	  double& alpha_primal, double& alpha_dual )
{
   double primalValue_p = -std::numeric_limits<double>::max();
   double primalStep_p = -std::numeric_limits<double>::max();
   double dualValue_p = -std::numeric_limits<double>::max();
   double dualStep_p = -std::numeric_limits<double>::max();
   double maxAlpha_p;

   double primalValue_d = -std::numeric_limits<double>::max();
   double primalStep_d = -std::numeric_limits<double>::max();
   double dualValue_d = -std::numeric_limits<double>::max();
   double dualStep_d = -std::numeric_limits<double>::max();
	double maxAlpha_d;

	bool primalBlocking, dualBlocking;

	iterate->findBlocking_pd( step,
					  primalValue_p, primalStep_p, dualValue_p, dualStep_p,
					  primalValue_d, primalStep_d, dualValue_d, dualStep_d,
					  maxAlpha_p, maxAlpha_d,
					  primalBlocking, dualBlocking );

	const double mufull = iterate->mustep_pd( step, maxAlpha_p, maxAlpha_d) / gamma_a;

	// No primal constraints were blocking?
	if( !primalBlocking )
	{
		alpha_primal = 1.0;
	}
	else
	{
	   const double dualValueEstim_p = dualValue_p + maxAlpha_d * dualStep_p;

      if( PIPSisEQ(dualValueEstim_p, 0.0 ) )
      {
         alpha_primal = 0.0; // to be corrected below
      }
      else
      {
         alpha_primal = ( - primalValue_p + mufull / ( dualValueEstim_p ) ) /
                  primalStep_p;
      }
	}

	// No dual constraints were blocking?
	if( !dualBlocking )
	{
		alpha_dual = 1.0;
	}
	else
	{
	   const double primValueEstim_d = primalValue_d + maxAlpha_p * primalStep_d;

	   if( PIPSisEQ(primValueEstim_d, 0.0 ) )
	   {
         alpha_dual = 0.0; // to be corrected below
	   }
	   else
	   {
         alpha_dual = ( - dualValue_d + mufull / ( primValueEstim_d ) ) /
               dualStep_d;
	   }
	}

	assert(alpha_primal <= 1.0);
   assert(alpha_dual <= 1.0);

	// safeguard against numerical troubles in the above computations
	alpha_primal = std::min( alpha_primal, maxAlpha_p );
	alpha_dual = std::min( alpha_dual, maxAlpha_d );

	// make it at least gamma_f * maxAlpha and no bigger than 1
	if( alpha_primal < gamma_f * maxAlpha_p ) alpha_primal = gamma_f * maxAlpha_p;
	if( alpha_dual < gamma_f * maxAlpha_d ) alpha_dual = gamma_f * maxAlpha_d;

	// back off just a touch (or a bit more)
	#ifdef STEPLENGTH_CONSERVATIVE
		alpha_primal *= 0.99;
		alpha_dual *= 0.99;
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

void Solver::doMonitorPd( Data * data, Variables * vars,
         Residuals * resids,
         double alpha_primal, double alpha_dual, double sigma,
         int i, double mu,
                        int stop_code,
         int level )
{
  OoqpMonitor * m = itsMonitors;

  while( m ) {
    m->doItPd( this, data, vars, resids, alpha_primal, alpha_dual, sigma, i, mu, stop_code, level );
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

  const double gap   = fabs( resids->dualityGap() );
  const double rnorm = resids->residualNorm();

  const int myrank = PIPS_MPIgetRank();

  idx = iterate - 1;
  if( idx <  0     ) idx = 0;
  if( idx >= maxit ) idx = maxit-1;

  // store the historical record
  mu_history[idx] = mu;
  rnorm_history[idx] = rnorm;
  phi = (rnorm + gap) / dnorm;
  phi_history[idx] = phi;

  if(idx > 0) {
    phi_min_history[idx] = phi_min_history[idx - 1];
    if(phi < phi_min_history[idx]) phi_min_history[idx] = phi;
  } else
    phi_min_history[idx] = phi;

  if ( iterate >= maxit ) {
    stop_code = MAX_ITS_EXCEEDED;
  } else if ( mu <= mutol && rnorm <= artol * dnorm ) {
    stop_code = SUCCESSFUL_TERMINATION;
  }

  if( myrank == 0 )
  {
     std::cout << "mu/mutol: " << mu << "  " << mutol << "  ....   rnorm/limit: " << rnorm << " " << artol * dnorm << std::endl;

     if( printTimeStamp )
     {
        const double timestamp = MPI_Wtime() - startTime;
        std::cout << "time stamp: " << timestamp << std::endl;
     }
  }

  if(stop_code != NOT_FINISHED)  return stop_code;

  // check infeasibility condition
  if(idx >= 10 && phi >= 1.e-8 && phi >= 1.e4 * phi_min_history[idx]) {
#ifdef TIMING
    if( myrank == 0 )
       std::cout << "possible INFEASIBLITY detected, phi: " << phi << std::endl;
#endif
    stop_code = INFEASIBLE;
  }
  if(stop_code != NOT_FINISHED)  return stop_code;

  // check for unknown status: slow convergence first
  if(idx >= 350 && phi_min_history[idx] >= 0.5 * phi_min_history[idx - 30]) {
    stop_code = UNKNOWN;
    printf("hehe dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm, artol);
  }

  if(idx >= 350 && rnorm > artol * dnorm &&
     rnorm_history[idx] * mu_history[0] >= 1.e8 * mu_history[idx] * rnorm_history[0])
  {
    stop_code = UNKNOWN;
    printf("dnorm=%g rnorm=%g artol=%g\n", rnorm, dnorm, artol);
  }

  //if(dnorm * mu < 50 * rnorm || mu < 1e-5) {
  if( mu * dnorm < 1.0e5 * rnorm) {
    //if(!onSafeSolver) {
    gLackOfAccuracy=1;
    //cout << "Lack of accuracy detected ---->" << mu << ":" << rnorm/dnorm << endl;
  } else {
      //if(dnorm * mu > 1e7 * rnorm && mu > 1.0e4)
      //gLackOfAccuracy=-1;
      //else
      gLackOfAccuracy=1;
  }
  //onSafeSolver=1;
  //}  
  //gLackOfAccuracy=-1; //disable iter refin in sLinsysRootAug
  
  return stop_code;
}


