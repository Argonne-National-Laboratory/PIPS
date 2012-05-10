#include "QpGenStochDriver.h"

#include <memory>
#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;
#include <cstdlib>

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#include "QpGenStoch.h"
#include "StochMonitor.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "SimpleVector.h"
#include "Status.h"
#include "QpGenData.h"
#include "OoqpVersion.h"

/*
extern int gOoqpPrintLevel;

//#define SOLVER MehrotraStochSolver
//#define SOLVER MehrotraSolver

template<class SOLVER, class FORMULATION>
int qpgenstoch_solve( int argc, char *argv[], 
		      StochInputTree* tree, 
		      StochRunParams* params,
		      SOLVER* s, 
		      FORMULATION* qp)
{
  MPI_Init(&argc, &argv);

  int iRet = qpgenstoch_solve(tree, params, s, qp);

  MPI_Finalize(); 

  return iRet;
}

template<class SOLVER, class FORMULATION>
int qpgenstoch_solve( StochInputTree* tree, 
		      StochRunParams* params,
		      SOLVER * , 
		      FORMULATION * )
{

  try {
    int iErr;
#ifdef HAVE_GETRUSAGE
    rusage before_read;
    getrusage( RUSAGE_SELF, &before_read );
#endif
    
    //QpGenStoch* qp = new QpGenStochAugExt( tree );
    FORMULATION* qp = new FORMULATION(tree);
    QpGenData* prob = (QpGenData * ) qp->makeData();


#ifdef HAVE_GETRUSAGE
    rusage  after_read;
    getrusage( RUSAGE_SELF, &after_read );
#endif

    QpGenVars     * vars  = (QpGenVars * ) qp->makeVariables( prob );
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
    StochIterateResourcesMonitor execTmMonitor;
    execTmMonitor.recIterateTm_start();

    int result = s->solve(prob, vars, resid);

    execTmMonitor.recIterateTm_stop();

    if( 0 == result ) {
      if( gOoqpPrintLevel > 0 ) {
#ifdef HAVE_GETRUSAGE
	rusage  after_solve;
	getrusage( RUSAGE_SELF, &after_solve );
	double read_time =
	  (after_read.ru_utime.tv_sec - before_read.ru_utime.tv_sec)
	  + (after_read.ru_utime.tv_usec - before_read.ru_utime.tv_usec)
	  / 1000000.0;
	double solve_time =
	  (after_solve.ru_utime.tv_sec - before_solve.ru_utime.tv_sec)
	  + (after_solve.ru_utime.tv_usec - before_solve.ru_utime.tv_usec)
	  / 1000000.0;

#endif
	
	// OUTPUT of the statistics
	int rank=0, size; MPI_Comm_size(MPI_COMM_WORLD, &size);
	if(size>1) MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double objective = prob->objectiveValue(vars);
	
	if(0==rank) {
	  
	  //cout << tree->NX() << " variables, " 
	  //     << tree->MY()  << " equality constraints, " 
	  //     << tree->MZ()  << " inequality constraints.\n";
	  cout << 11 << " variables, " 
	       << 12  << " equality constraints, " 
	       << 17  << " inequality constraints.\n";
	  
	  cout << "Iterates: " << s->iter
	       <<",    Optimal Solution:  " << objective << endl;
	  
	  cout << "Problem solved in " << execTmMonitor.tmIterate << " sec" << endl;
	  if(params->printx) {
	    cout << "The x-solution is:\n";
	    vars->x->writeToStream(cout);
	  } 
	}
      }
    } else {
      if ( gOoqpPrintLevel > 0 ) {
	cout << "Could not solve this QP.\n";
	cout << "Terminated with code " << result;
	if( result > 0 && result <= UNKNOWN ) {
	  cout << " : " << TerminationStrings[result];
	}
	cout << ".\n";
      }
    }

    delete s;
    delete vars;  
    delete resid;
    delete prob;
    delete qp;

    return result;
  } 
  catch( ... ) {
    cerr << "\nOops, out of memory\n";
    return -1;
  }
}

StochRunParams* defaultStochRunParams()
{
  StochRunParams* params = new StochRunParams();
  params->scale=0;
  params->printx=0;
  params->printLevel=0;
  return params;
}

*/
