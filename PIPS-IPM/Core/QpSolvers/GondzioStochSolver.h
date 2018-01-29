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
  unsigned int n_linesearch_points;
  Variables* temp_step;

public:

  GondzioStochSolver( ProblemFormulation * of, Data * prob, unsigned int n_linesearch_points = 10,
        bool adaptive_linesearch = true );

  virtual ~GondzioStochSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resid );

  // returns Gondzio weight for corrector step
  virtual void calculateAlphaWeightCandidate(Variables *iterate, Variables* predictor_step, Variables* corrector_step, double predictor_alpha,
        double& alpha_candidate, double& weight_candidate);

};

#endif /* PIPS_IPM_GONDZIOSTOCHSOLVER_H */
