/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUGEMTL
#define STOCHACTORYAUGEMTL

#include "sFactory.h"
#include "EmtlContext.h"

class sFactoryAugEmtl : public sFactory {
 public:
  sFactoryAugEmtl( StochInputTree* );
 private:
  sFactoryAugEmtl( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  sFactoryAugEmtl();
 public:
  EmtlContext *ctx;
  virtual ~sFactoryAugEmtl();

  virtual sLinsysRoot* newLinsysRoot();
  virtual sLinsysRoot* newLinsysRoot(QpGenStochData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
