/*
 * sFactoryAugMumpsLeaf.C
 *
 *      Author: bzfrehfe
 */


#include "sFactoryAugMumpsLeaf.h"
#include "sLinsysLeafMumps.h"
#include "MumpsSolver.h"


sLinsysLeaf* sFactoryAugMumpsLeaf::newLinsysLeaf(sData* prob,
                   OoqpVector* dd,OoqpVector* dq,
                   OoqpVector* nomegaInv, OoqpVector* rhs)
{
   MumpsSolver* linSolver = NULL;
   return new sLinsysLeafMumps(this, prob, dd, dq, nomegaInv, rhs, linSolver);
}
