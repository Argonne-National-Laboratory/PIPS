/*
 * sFactoryAugMumpsLeaf.C
 *
 *      Author: bzfrehfe
 */


#include "sFactoryAugMumpsLeaf.h"

#include "../LinearSolvers/MumpsSolver/MumpsSolverLeaf.h"
#include "sLinsysLeafMumps.h"


sLinsysLeaf* sFactoryAugMumpsLeaf::newLinsysLeaf(sData* prob,
                   OoqpVector* dd,OoqpVector* dq,
                   OoqpVector* nomegaInv, OoqpVector* rhs)
{
   MumpsSolverLeaf* linSolver = NULL;
   return new sLinsysLeafMumps(this, prob, dd, dq, nomegaInv, rhs, linSolver);
}
