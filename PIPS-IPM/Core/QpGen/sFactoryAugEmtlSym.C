/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */


#include "sFactoryAugEmtlSym.h"

#include "QpGenStochData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAugEmtlSym.h"
#include "StochResourcePlanner.h"

sFactoryAugEmtlSym::sFactoryAugEmtlSym( StochInputTree* inputTree)
  : sFactory(inputTree)
{ 
  MPI_Comm world = tree->commWrkrs;
  int emtlprocs = StochResourcePlanner::noScaProcesses;
  ctx = new EmtlContext(world, emtlprocs);

  if (ctx->mype() == 0) {
    int nb = Blocksize();
#ifdef TIMING
    printf("EMTL NPROCS %d GRID %d %d BLOCKSIZE %d\n",
      ctx->nprow()*ctx->npcol(), ctx->nprow(), ctx->npcol(), nb);
#else
    printf("Using %d processes for Elemental, %d by %d grid\n, blocksize %d", 
      ctx->nprow()*ctx->npcol(), ctx->nprow(), ctx->npcol(), nb);
#endif
  }
  
  
};

sFactoryAugEmtlSym::sFactoryAugEmtlSym( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : sFactory(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

sFactoryAugEmtlSym::sFactoryAugEmtlSym()
{ 
  assert("Don't call this" && 0);
};

sFactoryAugEmtlSym::~sFactoryAugEmtlSym()
{ 
  delete ctx;
};


sLinsysRoot* sFactoryAugEmtlSym::newLinsysRoot()
{
  return new sLinsysRootAugEmtlSym(this, data, *ctx);
}

sLinsysRoot* 
sFactoryAugEmtlSym::newLinsysRoot(QpGenStochData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs)
{
  return new sLinsysRootAugEmtlSym(this, prob, dd, dq, nomegaInv, rhs, *ctx);
}
