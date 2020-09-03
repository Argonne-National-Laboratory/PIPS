/*
 * GondzioStochLpSolver.h
 *
 *  Created on: 20.12.2017
 *      Author: Svenja Uslu
 */

#ifndef PIPS_IPM_CORE_QPSOLVERS_GONDZIOSTOCHLPSOLVER_H_
#define PIPS_IPM_CORE_QPSOLVERS_GONDZIOSTOCHLPSOLVER_H_

#include "GondzioStochSolver.h"

class Data;
class Variables;
class ProblemFormulation;


/**
 * Derived class of Solver implementing Gondzio-correction version of
 * Mehrotra's original predictor-corrector algorithm
 * allowing different step primal and dual step lengths.
 * @ingroup QpSolvers
 */
class GondzioStochLpSolver : public GondzioStochSolver
{
private:
  // returns Gondzio weight for corrector step for different alpha_primal and alpha_dual
  virtual void calculateAlphaPDWeightCandidate(Variables *iterate, Variables* predictor_step,
         Variables* corrector_step, double alpha_primal, double alpha_dual,
         double& alpha_primal_candidate, double& alpha_dual_candidate,
         double& weight_primal_candidate, double& weight_dual_candidate);

  void computeProbingStep_pd(Variables* probing_step, const Variables* iterate, const Variables* step,
        double alpha_primal, double alpha_dual) const;
public:

  GondzioStochLpSolver( ProblemFormulation * of, Data * prob, const Scaler* scaler = nullptr );

  virtual ~GondzioStochLpSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resid );

};


#endif /* PIPS_IPM_CORE_QPSOLVERS_GONDZIOSTOCHLPSOLVER_H_ */
