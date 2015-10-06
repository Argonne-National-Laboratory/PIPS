/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPGENFACTORY_SOLVER
#define NLPGENFACTORY_SOLVER

#include "NlpGenSparse.h"


class NlpGenSparseWithSolver : public NlpGenSparse {

public:
	
  NlpGenSparseWithSolver( int nx, int my, int mz,
		       int nnzQ, int nnzA, int nnzC );
  
  virtual LinearSystem * makeLinsys( Data * prob_in );
  virtual LinearSystem * makeLinsys( Data * prob_in , int setMulti);

  // use full Q or just diagonal part of Q;
  bool fullQ;
  
};


#endif

