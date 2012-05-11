/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUGEMTLSYM
#define STOCHACTORYAUGEMTLSYM

#include "sFactory.h"
#include "EmtlContext.h"

class sFactoryAugEmtlSym : public sFactory {
 public:
  sFactoryAugEmtlSym( StochInputTree* );
 private:
  sFactoryAugEmtlSym( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  sFactoryAugEmtlSym();
 public:
  EmtlContext *ctx;
  virtual ~sFactoryAugEmtlSym();

  virtual sLinsysRoot* newLinsysRoot();
  virtual sLinsysRoot* newLinsysRoot(QpGenStochData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);
};
#endif
