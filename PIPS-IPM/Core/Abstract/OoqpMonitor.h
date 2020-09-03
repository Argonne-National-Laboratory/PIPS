/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef OOQPMONITOR
#define OOQPMONITOR
#include <cassert>

class Solver;
class Data;
class Variables;
class Residuals;

/** Represents objects that display progress information for interior
 *  point QP solvers.
 *  @ingroup QpSolvers
 */
class OoqpMonitor {
public:
  OoqpMonitor * nextMonitor;

  OoqpMonitor() { nextMonitor = 0; };

  virtual void doIt( const Solver * solver, const Data * data, const Variables * vars,
					 const Residuals * resids,
					 double alpha, double sigma,
					 int i, double mu, 
                     int status_code,
					 int level ) = 0;

  virtual void doItPd( const Solver * solver, const Data * data, const Variables * vars,
                const Residuals * resids,
                double alpha_primal, double alpha_dual, double sigma,
                int i, double mu,
                     int status_code,
                int level ) { assert(0 && "not implemented here"); };

  virtual ~OoqpMonitor() {};
};  


/** Monitors that simply call the solver's defaultMonitor method.
 * 
 *  Don't create instances of this class. Call the solver's monitorSelf
 *  method instead.
 *
 *  @ingroup QpSolvers
 */
class OoqpSelfMonitor : public OoqpMonitor {
public:
  void doIt( const Solver * solver, const Data * data, const Variables * vars,
					 const Residuals * resids,
					 double alpha, double sigma,
					 int i, double mu,
                     int status_code,
					 int level ) override;
};

#include "OoqpMonitorData.h"

/**
 * Represents monitors that use a C function to print progress information
 * for an algorithm.
 * @ingroup QpSolvers
 */
class COoqpMonitor : public OoqpMonitor {
protected:
  DoItCFunc doItC;
  void * ctx;
public:
  COoqpMonitor( DoItCFunc doItC_, void * ctx_ );
  void doIt( const Solver * solver, const Data * data, const Variables * vars,
					 const Residuals * resids,
					 double alpha, double sigma,
					 int i, double mu,
                     int status_code,
					 int level ) override;
};


#endif
