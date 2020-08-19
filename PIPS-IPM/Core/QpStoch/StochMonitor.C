#include "StochMonitor.h"
#include "Status.h"

#include "Residuals.h"
#include "Solver.h"
#include "QpGenData.h"
#include "QpGenVars.h"
#include "sTree.h"
#include <iostream>
#include <cstdio>
#include "pipsport.h"

StochMonitor::StochMonitor(QpGenStoch* qp_, Scaler* scaler)
  : qp(qp_), scaler(scaler)
{
  mpiComm=MPI_COMM_WORLD; //default for old version
  MPI_Comm_rank(mpiComm, &myRank);
  myGlobRank = myRank;
}

StochMonitor::StochMonitor(sFactory* qp_, Scaler* scaler)
  : qp(nullptr), scaler(scaler)
{
  mpiComm = qp_->tree->commWrkrs;
  MPI_Comm_rank(mpiComm, &myRank);
  MPI_Comm_rank(MPI_COMM_WORLD, &myGlobRank);
}

void StochMonitor::doIt( const Solver * solver, const Data * data, const Variables * vars,
			 const Residuals * resids,
			 double alpha, double sigma,
			 int i, double mu,
			 int status_code,
			 int level )
{
   StochMonitor::doItStoch(solver, data, vars, resids, alpha, -1.0, sigma, i, mu, status_code, level);
}

void StochMonitor::doItPd( const Solver * solver, const Data * data, const Variables * vars,
              const Residuals * resids,
              double alpha_primal, double alpha_dual, double sigma,
              int i, double mu,
              int status_code,
              int level )
{
   StochMonitor::doItStoch(solver, data, vars, resids, alpha_primal, alpha_dual, sigma, i, mu, status_code, level);
}

void StochMonitor::doItStoch( const Solver * solver, const Data * data, const Variables * vars,
              const Residuals * resids,
              double alpha_primal, double alpha_dual, double sigma,
              int i, double mu,
              int status_code,
              int level ) const
{
  double objective = dynamic_cast<const QpGenData*>(data)->objectiveValue(dynamic_cast<const QpGenVars*>(vars));

  const Residuals* resids_unscaled = resids;
  if( scaler )
  {
     objective = scaler->getObjUnscaled(objective);
     resids_unscaled = scaler->getResidualsUnscaled(*resids);
  }

  const double dnorm = solver->dataNormOrig();
  const double rnorm = resids_unscaled->residualNorm();
  const double gap = resids_unscaled->dualityGap();

  if( scaler )
     delete resids_unscaled;

  // log only on the first proc
   if( myRank > 0 )
      return;

  switch( level ) {
     case 0 : case 1: {
        std::cout << " --- Iteration " << i << " --- (rank " << myGlobRank << ")" << std::endl;
    if( i == 1 )
      printf(" mu = %16.12e  rel.res.norm=%16.12e  datanorm=%16.12e\n", 
	     mu, rnorm / dnorm, dnorm);
    else
      printf(" mu = %16.12e  rel.res.norm=%16.12e\n",
             mu, rnorm / dnorm);
    //cout << " mu = " << mu << " relative residual norm = " 
    //cout << resids->residualNorm() / dnorm << endl;
    std::cout << " Duality Gap:  " << gap << std::endl;
    if( i > 1 )
    {
       if( alpha_dual != -1.0 )
       {
          std::cout << " alpha primal = " << alpha_primal << std::endl;
          std::cout << " alpha dual = " << alpha_dual << std::endl;
       }
       else
          std::cout << " alpha = " << alpha_primal << std::endl;
    }
    std::cout << " Objective: " << objective << std::endl;
    std::cout << std::endl;
    if( level == 1) { 
      // Termination has been detected by the status check; print
      // appropriate message
      switch( status_code ) {
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
