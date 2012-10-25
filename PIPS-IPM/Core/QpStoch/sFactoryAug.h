/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUG
#define STOCHACTORYAUG

#include "sFactory.h"

class sFactoryAug : public sFactory {
 public:
  sFactoryAug( StochInputTree* );
  sFactoryAug( stochasticInput&, MPI_Comm comm=MPI_COMM_WORLD );
 private:
  sFactoryAug( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  sFactoryAug();
 public:
  virtual ~sFactoryAug();

  virtual sLinsysRoot* newLinsysRoot();
  virtual sLinsysRoot* newLinsysRoot(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
