/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sFactoryAugComm2SchurLeaf.h"

#include "sData.h"

#include "StochTree.h"
#include "StochInputTree.h"
#include "PardisoSolver.h"
#include "PardisoSchurSolver.h"

#include "sLinsysRootAugComm2.h"
#include "sLinsysLeafSchurSlv.h"

sLinsysLeaf* 
sFactoryAugComm2SchurLeaf::newLinsysLeaf(sData* prob,
					 OoqpVector* dd,
					 OoqpVector* dq,
					 OoqpVector* nomegaInv, 
					 OoqpVector* rhs)
{
    //cout << "sFactoryAugComm2SchurLeaf::newLinsysLeaf  returns  a sLinsysLeafSchurSlv" << endl;
  PardisoSchurSolver* linSolver=NULL;
  return new sLinsysLeafSchurSlv(this, prob, dd, dq, nomegaInv, rhs, linSolver);
}

sLinsysRoot* sFactoryAugComm2SchurLeaf::newLinsysRoot()
{
    //cout << "sFactoryAugComm2SchurLeaf::newLinsysRoot()" << endl;
  return new sLinsysRootAugComm2(this, data);
}

sLinsysRoot* 
sFactoryAugComm2SchurLeaf::newLinsysRoot(sData* prob,
					   OoqpVector* dd,
					   OoqpVector* dq,
					   OoqpVector* nomegaInv, 
					   OoqpVector* rhs)
{
    //cout << "sFactoryAugComm2SchurLeaf::newLinsysRoot(arguments)" << endl;
  return new sLinsysRootAugComm2(this, prob, dd, dq, nomegaInv, rhs);
}

/*sFactoryAugComm2Schur32Leaf::sFactoryAugComm2Schur32Leaf( StochInputTree* inputTree)
  : sFactory(inputTree)
{ };

sFactoryAugComm2Schur32Leaf::sFactoryAugComm2Schur32Leaf( stochasticInput& in, MPI_Comm comm)
  : sFactory(in,comm)
{ }
*/
/*
sFactoryAugComm2Schur32Leaf::sFactoryAugComm2Schur32Leaf( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : sFactory(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

sFactoryAugComm2Schur32Leaf::sFactoryAugComm2Schur32Leaf()
{ };
*/
