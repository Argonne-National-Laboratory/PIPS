/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#ifndef SolverOption_H
#define SolverOption_H


/** 
 * Abstract base class for QP solvers.
 */
class SolverOption
{
public:

  /** parameter in range [0,100] determines verbosity. (Higher value
   *  => more verbose.) */
  int        printlevel;

	
	SolverOption();
	virtual ~SolverOption();
};

#endif
