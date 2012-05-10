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
  MehrotraStochSolver( ProblemFormulation * opt, Data * prob );

  ~MehrotraStochSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resids );

};

#endif

