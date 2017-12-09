/*
 * GondzioStochSolver.h
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_GONDZIOSTOCHSOLVER_H
#define PIPS_IPM_GONDZIOSTOCHSOLVER_H

#include "GondzioSolver.h"


class Data;
class Variables;
class ProblemFormulation;


/**
 * Derived class of Solver implementing Gondzio-correction version of
 * Mehrotra's original predictor-corrector algorithm.
 * @ingroup QpSolvers
 */
class GondzioStochSolver : public GondzioSolver
{
protected:

  /** parameter in range [0,100] determines verbosity. (Higher value
   *  => more verbose.) */
  int        printlevel;

  /** exponent in Mehrotra's centering parameter, which is usually
   *  chosen to me (muaff/mu)^tsig, where muaff is the predicted
   *  complementarity gap obtained from an affine-scaling step, while
   *  mu is the current complementarity gap */
  double     tsig;

  /** maximum number of Gondzio corrector steps */
  int        maximum_correctors;

  /** actual number of Gondzio corrections needed */
  int        NumberGondzioCorrections;

  /** various parameters associated with Gondzio correction */
  double     StepFactor0, StepFactor1, AcceptTol, beta_min, beta_max;

  /**  storage for step vectors */
  Variables *corrector_step, *step;

  /** storage for residual vectors */
  Residuals *corrector_resid;

  ProblemFormulation * factory;

public:

  GondzioStochSolver( ProblemFormulation * of, Data * prob );

  virtual ~GondzioStochSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resid );

  /** reset parameters to their default values */
  virtual void reset_parameters() {};

  virtual void defaultMonitor( Data * data, Variables * vars,
                        Residuals * resids,
                        double alpha, double sigma,
                        int i, double mu,
                        int status_code,
                        int level ) ;

};

#endif /* PIPS_IPM_GONDZIOSTOCHSOLVER_H */
