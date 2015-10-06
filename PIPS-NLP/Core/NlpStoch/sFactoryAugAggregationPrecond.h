/* PIPS-NLP                                                           	*
 * Author:  Nai-Yuan Chiang                                       	*
 * (C) 2015 Argonne National Laboratory. 				*/

#ifndef STOCHACTORYAUG_AGGREGATION
#define STOCHACTORYAUG_AGGREGATION

#include "sFactory.h"

class sFactoryAugAggregationPrecond : public sFactory {
 public:
  sFactoryAugAggregationPrecond( StochInputTree* in)
    : sFactory(in) {};
  sFactoryAugAggregationPrecond( stochasticInput& in, MPI_Comm comm=MPI_COMM_WORLD)
   : sFactory(in,comm) {};

  virtual sLinsysRoot* newLinsysRoot();
  virtual sLinsysRoot* newLinsysRoot(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag);


  virtual sLinsysLeaf* newLinsysLeaf(sData* prob, OoqpVector* dd, OoqpVector* dq,
			OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag);
};

#endif
