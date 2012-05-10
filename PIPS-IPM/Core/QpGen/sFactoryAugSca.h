/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUGSCA
#define STOCHACTORYAUGSCA

#include "sFactory.h"
#include "scalapack.h"

class sFactoryAugSca : public sFactory {
 public:
  sFactoryAugSca( StochInputTree* );
 private:
  sFactoryAugSca( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  sFactoryAugSca();
 public:
  COMMINFO cinfo;
  virtual ~sFactoryAugSca();

  virtual sLinsysRoot* newLinsysRoot();
  virtual sLinsysRoot* newLinsysRoot(QpGenStochData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
