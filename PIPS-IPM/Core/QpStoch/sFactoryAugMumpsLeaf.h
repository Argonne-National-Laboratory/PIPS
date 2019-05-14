/*
 * sFactoryAugMumpsLeaf.h
 *
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPSTOCH_SFACTORYAUGMUMPSLEAF_H_
#define PIPS_IPM_CORE_QPSTOCH_SFACTORYAUGMUMPSLEAF_H_

#include "sFactoryAug.h"

class sFactoryAugMumpsLeaf : public sFactoryAug {
 public:

  sFactoryAugMumpsLeaf( StochInputTree* inputTree, MPI_Comm comm=MPI_COMM_WORLD )
   : sFactoryAug(inputTree, comm) {};

  sFactoryAugMumpsLeaf( stochasticInput& in, MPI_Comm comm=MPI_COMM_WORLD )
   : sFactoryAug(in,comm) {};


  sLinsysLeaf* newLinsysLeaf(sData* prob,
              OoqpVector* dd,OoqpVector* dq,
              OoqpVector* nomegaInv, OoqpVector* rhs);
};



#endif /* PIPS_IPM_CORE_QPSTOCH_SFACTORYAUGMUMPSLEAF_H_ */
