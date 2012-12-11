/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUG_COMM2_SCHURLEAF
#define STOCHACTORYAUG_COMM2_SCHURLEAF

#include "sFactory.h"

class sFactoryAugComm2SchurLeaf : public sFactory {
 public:

  sFactoryAugComm2SchurLeaf( StochInputTree* in)
    : sFactory(in) {};
 sFactoryAugComm2SchurLeaf( stochasticInput& in, MPI_Comm comm=MPI_COMM_WORLD)
   : sFactory(in,comm) {};


  sLinsysLeaf* newLinsysLeaf(sData* prob,
			     OoqpVector* dd,OoqpVector* dq,
			     OoqpVector* nomegaInv, OoqpVector* rhs);

  virtual sLinsysRoot* newLinsysRoot();
  virtual sLinsysRoot* newLinsysRoot(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);
};

#endif
