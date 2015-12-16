/* PIPS-IPM                                                             
 * Author: Cosmin G. Petra
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
 
/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef QPGENSTOCHDRIVER
#define QPGENSTOCHDRIVER

#include <memory>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdlib>

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#include <stdexcept>

#include "StochInputTree.h"
#include "sTree.h"
#include "MehrotraSolver.h"

//#include "NlpGenStoch.h"
#include "StochMonitor.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"
#include "SimpleVector.h"
#include "Status.h"
#include "NlpGenData.h"
#include "OoqpVersion.h"
//#include "StochResourcesMonitor.h"

struct StochRunParams
{
  int scale;
  int printx;
  int printLevel;
};

struct ScaParams
{
  int nscaprocs;
  int nb;
};

inline StochRunParams* defaultStochRunParams()
{
  StochRunParams* params = new StochRunParams();
  params->scale=0; params->printx=0; params->printLevel=0;
  return params;
}

extern int gOoqpPrintLevel;

/** 
 * Solve the stochastic problem specified by argument 'tree'.
 * The function does NOT need the MPI environment to be initialized.
 */
template<class SOLVER, class FORMULATION>
int nlpstoch_solve( int argc, char *argv[], 
		      StochInputTree* tree, 
		      StochRunParams* params,
		      SOLVER* s, 
		      FORMULATION* qp)
{
  MPI_Init(&argc, &argv);

  int iRet = nlpstoch_solve(tree, params, s, qp);

  MPI_Finalize(); 

  return iRet;
}


/**
 * Solve the stochastic problem specified by argument 'tree'.
 * The MPI environment must be intialized and all the processes from 
 * MPI_COMM_WORLD are used in solving the problem.
 */
template<class SOLVER, class FORMULATION>
int nlpstoch_solve( StochInputTree* tree, 
		      StochRunParams* params,
		      SOLVER * , 
		      FORMULATION * )
{
  try {
#ifdef HAVE_GETRUSAGE
    rusage before_read;
    getrusage( RUSAGE_SELF, &before_read );
#endif
    
    FORMULATION* qp = new FORMULATION(tree);
    NlpGenData* prob = (NlpGenData * ) qp->makeData();


#ifdef HAVE_GETRUSAGE
    rusage  after_read;
    getrusage( RUSAGE_SELF, &after_read );
#endif

    NlpGenVars     * vars  = (NlpGenVars * ) qp->makeVariables( prob );
    Residuals     * resid = qp->makeResiduals( prob );
    SOLVER * s            = new SOLVER( qp, prob );
    
    //s->monitorSelf();
    s->addMonitor(new StochMonitor(qp));

#ifdef STOCH_TESTING
    //tree->displayProcessInfo();
#endif

#ifdef HAVE_GETRUSAGE
    rusage before_solve;
    getrusage( RUSAGE_SELF, &before_solve );
#endif
    // this monitor just monitors total time, not iterations
    StochIterateResourcesMonitor execTmMonitor;
    execTmMonitor.recIterateTm_start();

    int result = s->solve(prob, vars, resid);

    execTmMonitor.recIterateTm_stop();
    
    int rank=0, size; MPI_Comm_size(MPI_COMM_WORLD, &size);
    if(size>1) MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if( 0 == result ) {
      if( gOoqpPrintLevel > 0 ) {
#ifdef HAVE_GETRUSAGE
	rusage  after_solve;
	getrusage( RUSAGE_SELF, &after_solve );
	//double read_time =
	//  (after_read.ru_utime.tv_sec - before_read.ru_utime.tv_sec)
	//  + (after_read.ru_utime.tv_usec - before_read.ru_utime.tv_usec)
	//  / 1000000.0;
	//double solve_time =
	//  (after_solve.ru_utime.tv_sec - before_solve.ru_utime.tv_sec)
	//  + (after_solve.ru_utime.tv_usec - before_solve.ru_utime.tv_usec)
	//  / 1000000.0;

#endif
	
	// OUTPUT of the statistics

	double objective = prob->objectiveValue(vars);
	
	if(0==rank) {
	  sTree* stTree = qp->tree;
	  cout << stTree->N << " variables, " 
	       << stTree->MY  << " equality constraints, " 
	       << stTree->MZ  << " inequality constraints.\n";

	    //cout << 11 << " variables, " 
	    //   << 12  << " equality constraints, " 
	    //   << 17  << " inequality constraints.\n";
#ifdef TIMING
    printf("NITER %d\n", s->iter);
    printf("SOL %f\n", objective);
    printf("TOTTIME %f SEC\n", execTmMonitor.tmIterate);
#else
    cout << "Iterates: " << s->iter
	       <<",    Optimal Solution:  " << objective << endl;
	  
	  cout << "Problem solved in " << execTmMonitor.tmIterate << " sec" << endl;
#endif
    if(params->printx) {
	    cout << "The x-solution is:\n";
	    vars->x->writeToStream(cout);
	  } 
	}
      }
    } else {
      if ( gOoqpPrintLevel > 0 && rank==0 ) {
	cout << "Could not solve this QP.\n";
	cout << "Terminated with code " << result;
	if( result > 0 && result <= UNKNOWN ) {
	  cout << " : " << TerminationStrings[result];
	}
	cout << ".\n";
      }
    }
    if (0==rank) {
#ifdef TIMING
      cout << "MAXMEM " << StochIterateResourcesMonitor::getMaxMemUsage() << 
        " KB" << endl;
      //double objective = prob->objectiveValue(vars);
      //cout << "OBJVAL " << objective << endl;
#else
	    cout << "Maximum memory usage: " << 
       StochIterateResourcesMonitor::getMaxMemUsage() << " kb" << endl;
#endif
    }
    delete s;
    delete vars;  
    delete resid;
    delete prob;
    delete qp;

    return result;
  }
  catch( logic_error &e ) {
    cerr << "Elemental error: " << e.what() << endl;
    return -1;
  }
  catch( ... ) {
    cerr << "\nOops, out of memory\n";
    return -1;
  }

}

#endif
