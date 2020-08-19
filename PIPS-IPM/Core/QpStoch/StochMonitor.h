#ifndef STOCH_MONITOR
#define STOCH_MONITOR

#include "OoqpMonitor.h"
#include "QpGenStoch.h"
#include "sFactory.h"
#include "Scaler.h"
#include "pipsport.h"

class Solver;
class Data;
class Variables;
class Residuals;

class StochMonitor : public OoqpMonitor 
{
 private:
  void doItStoch( const Solver * solver, const Data * data, const Variables * vars,
        const Residuals * resids,
        double alpha_primal, double alpha_dual, double sigma,
        int i, double mu,
                  int status_code,
        int level ) const;

 public:  
  StochMonitor(QpGenStoch* qp, Scaler* scaler = nullptr);
  StochMonitor(sFactory* qp, Scaler* scaler = nullptr);
  void doIt( const Solver * solver, const Data * data, const Variables * vars,
		     const Residuals * resids,
		     double alpha, double sigma,
		     int i, double mu,
		     int status_code,
		     int level ) override;

  void doItPd( const Solver * solver, const Data * data, const Variables * vars,
                const Residuals * resids,
                double alpha_primal, double alpha_dual, double sigma,
                int i, double mu,
                int status_code,
                int level ) override;

 protected:
  QpGenStoch* qp;
  Scaler* scaler;

  MPI_Comm mpiComm;
  int myRank, myGlobRank;
};

#endif
