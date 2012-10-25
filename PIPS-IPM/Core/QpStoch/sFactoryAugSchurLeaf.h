/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUG_SCHURLEAF
#define STOCHACTORYAUG_SCHURLEAF

#include "sFactoryAug.h"

class sFactoryAugSchurLeaf : public sFactoryAug {
 public:

  sFactoryAugSchurLeaf( StochInputTree* in)
    : sFactoryAug(in) {};
 sFactoryAugSchurLeaf( stochasticInput& in, MPI_Comm comm=MPI_COMM_WORLD)
   : sFactoryAug(in,comm) {};


  sLinsysLeaf* newLinsysLeaf(sData* prob,
			     OoqpVector* dd,OoqpVector* dq,
			     OoqpVector* nomegaInv, OoqpVector* rhs);

};
#endif
