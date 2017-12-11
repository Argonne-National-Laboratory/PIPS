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
public:

  GondzioStochSolver( ProblemFormulation * of, Data * prob );

  ~GondzioStochSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resid );

};

#endif /* PIPS_IPM_GONDZIOSTOCHSOLVER_H */
