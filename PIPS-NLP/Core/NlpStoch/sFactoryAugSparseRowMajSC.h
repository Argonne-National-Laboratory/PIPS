/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHACTORYAUGTRIPSC
#define STOCHACTORYAUGTRIPSC

#include "sFactoryAug.h"

class sFactoryAugSpTripletSC : public sFactoryAug {
public:
  sFactoryAugSpTripletSC( StochInputTree* t) : sFactoryAug(t) {};
  sFactoryAugSpTripletSC( stochasticInput& t, MPI_Comm comm=MPI_COMM_WORLD )
    : sFactoryAug(t,comm) {}
 private:
 public:
  virtual ~sFactoryAugSpTripletSC() {};

  //virtual LinearSystem* makeLinsys( Data * prob_in );

  virtual sLinsysRoot* newLinsysRoot();
  virtual sLinsysRoot* newLinsysRoot(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs, OoqpVector* additiveDiag);


};
#endif
 
