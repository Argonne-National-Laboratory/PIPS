/*
 * GondzioStochSolver.h
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef PIPS_IPM_GONDZIOSTOCHSOLVER_H
#define PIPS_IPM_GONDZIOSTOCHSOLVER_H

#include "GondzioSolver.h"
#include "Observer.h"

class Data;
class Variables;
class ProblemFormulation;

const unsigned int max_linesearch_points = 50;

/**
 * Derived class of Solver implementing Gondzio-correction version of
 * Mehrotra's original predictor-corrector algorithm.
 * @ingroup QpSolvers
 */
class GondzioStochSolver : public GondzioSolver, public Observer
{
private:
  // returns Gondzio weight for corrector step
  virtual void calculateAlphaWeightCandidate(Variables *iterate, Variables* predictor_step, Variables* corrector_step, double predictor_alpha,
        double& alpha_candidate, double& weight_candidate);

protected:
  unsigned int n_linesearch_points;
  Variables* temp_step;

  /** should additional corrector steps for small complementarity pairs be applied */
  const bool additional_correctors_small_comp_pairs;
  /** should additional corrector steps for small complementarity pairs be applied */
  const int max_additional_correctors;
  /** first iteration at which to look for small corrector steps */
  const int first_iter_small_correctors;
  /** alpha must be lower equal to this value for the IPM to try and apply small corrector steps */
  const double max_alpha_small_correctors;

  int NumberSmallCorrectors;

  /* observer stuff for checking convergence of BiCGStab */
  bool bicgstab_converged;
  double norm_bigcstab_res_rel;

  void registerBiCGStabOvserver(LinearSystem* sys);

public:

  GondzioStochSolver( ProblemFormulation * of, Data * prob, unsigned int n_linesearch_points = 10,
        bool adaptive_linesearch = true );

  virtual ~GondzioStochSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resid );

  void notifyFromSubject() override;
};

#endif /* PIPS_IPM_GONDZIOSTOCHSOLVER_H */
