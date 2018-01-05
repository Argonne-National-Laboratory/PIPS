#ifndef STOCH_MONITOR
#define STOCH_MONITOR

#include "OoqpMonitor.h"
#include "QpGenStoch.h"
#include "sFactory.h"
#include "Scaler.h"

class Solver;
class Data;
class Variables;
class Residuals;

class StochMonitor : public OoqpMonitor 
{
 public:  
  StochMonitor(QpGenStoch* qp, Scaler* scaler = NULL);
  StochMonitor(sFactory* qp, Scaler* scaler = NULL);
  virtual void doIt( Solver * solver, Data * data, Variables * vars,
		     Residuals * resids,
		     double alpha, double sigma,
		     int i, double mu,
                     int status_code,
		     int level );
 protected:
  QpGenStoch* qp;
  Scaler* scaler;

  MPI_Comm mpiComm;
  int myRank, myGlobRank;
};

#endif
