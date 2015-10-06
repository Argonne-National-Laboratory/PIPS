/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCH_MONITOR
#define STOCH_MONITOR

#include "OoqpMonitor.h"
#include "NlpGenStoch.h"
#include "sFactory.h"
class Solver;
class Data;
class Variables;
class Residuals;

class StochMonitor : public OoqpMonitor 
{
 public:  
  StochMonitor(NlpGenStoch* qp);
  StochMonitor(sFactory* qp);
  virtual void doIt( Solver * solver, Data * data, Variables * vars,
		     Residuals * resids,
		     double alpha, double sigma,
		     int i, double mu,
                     int status_code,
		     int level );
 protected:
  NlpGenStoch* qp;
  MPI_Comm mpiComm;
  int myRank, myGlobRank;
};

#endif
