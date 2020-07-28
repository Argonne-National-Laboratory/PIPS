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

const unsigned int max_linesearch_points = 50;

/**
 * Derived class of Solver implementing Gondzio-correction version of
 * Mehrotra's original predictor-corrector algorithm.
 * @ingroup QpSolvers
 */
class GondzioStochSolver : public GondzioSolver
{
private:
  // returns Gondzio weight for corrector step
  virtual void calculateAlphaWeightCandidate(Variables *iterate, Variables* predictor_step, Variables* corrector_step, double predictor_alpha,
        double& alpha_candidate, double& weight_candidate);

protected:
  unsigned int n_linesearch_points;
  Variables* temp_step;

  void setBiCGStabTol(int iteration) const;
  // controls whether setBiCGTol applies an dynamic schedule for the BiCGStab tolerance or just uses the user defined input (OUTER_BICG_TOL)
  bool dynamic_bicg_tol;

public:

  GondzioStochSolver( ProblemFormulation * of, Data * prob);

  virtual ~GondzioStochSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resid );
};

#endif /* PIPS_IPM_GONDZIOSTOCHSOLVER_H */
