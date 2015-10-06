/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef FILTERIPMALGORITHM_STOCH_H
#define FILTERIPMALGORITHM_STOCH_H

#include "FilterIPMSolver.h"

class Data;
class Variables;
class ProblemFormulation;

/** 
 * @ingroup NlpSolvers
 */
class FilterIPMStochSolver : public FilterIPMSolver
{

public:
  FilterIPMStochSolver( ProblemFormulation * opt, Data * prob );

  ~FilterIPMStochSolver();

  virtual int solve( Data *prob, Variables *iterate, Residuals * resids );


  virtual int defaultStatus(Data *  data_in, Variables * /* vars */,
			 Residuals * resids_in,
			 int iterate, double mu, 
			 int /* level */);

};

#endif



