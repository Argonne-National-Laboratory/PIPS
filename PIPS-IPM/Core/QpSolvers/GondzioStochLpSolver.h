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
protected:
  //const unsigned int n_linesearch_points;
  //Variables* temp_step;

public:

  GondzioStochLpSolver( ProblemFormulation * of, Data * prob, unsigned int n_linesearch_points = 10 );

  virtual ~GondzioStochLpSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resid );

  // returns Gondzio weight for corrector step
  virtual void calculateAlphaWeightCandidate(Variables *iterate, Variables* predictor_step, Variables* corrector_step, double predictor_alpha,
        double& alpha_candidate, double& weight_candidate);

  // returns Gondzio weight for corrector step for different alpha_primal and alpha_dual
  virtual void calculateAlphaPDWeightCandidate(Variables *iterate, Variables* predictor_step,
	  		Variables* corrector_step, double alpha_primal, double alpha_dual,
	  		double& alpha_primal_candidate, double& alpha_dual_candidate,
	  		double& weight_primal_candidate, double& weight_dual_candidate);

};


#endif /* PIPS_IPM_CORE_QPSOLVERS_GONDZIOSTOCHLPSOLVER_H_ */
