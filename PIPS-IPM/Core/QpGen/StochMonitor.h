#ifndef STOCH_MONITOR
#define STOCH_MONITOR

#include "OoqpMonitor.h"
#include "QpGenStoch.h"
#include "sFactory.h"
class Solver;
class Data;
class Variables;
class Residuals;

class StochMonitor : public OoqpMonitor 
{
 public:  
  StochMonitor(QpGenStoch* qp);
  StochMonitor(sFactory* qp);
  virtual void doIt( Solver * solver, Data * data, Variables * vars,
		     Residuals * resids,
		     double alpha, double sigma,
		     int i, double mu,
                     int status_code,
		     int level );
 protected:
  QpGenStoch* qp;
};

#endif
