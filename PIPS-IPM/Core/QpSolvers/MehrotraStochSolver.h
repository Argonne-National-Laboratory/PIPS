#ifndef MEHALGORITHMSTOCH_H
#define MEHALGORITHMSTOCH_H

#include "MehrotraSolver.h"

class Data;
class Variables;
class ProblemFormulation;

/** Derived class of Solver implementing the original Mehrotra
 *  predictor-corrector algorithm 
 * @ingroup QpSolvers
 */
class MehrotraStochSolver : public MehrotraSolver
{

public:
  MehrotraStochSolver( ProblemFormulation * opt, Data * prob, const Scaler * scaler = nullptr );

  ~MehrotraStochSolver();

  int solve( Data *prob, Variables *iterate, Residuals * resids ) override;

};

#endif

