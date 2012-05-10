/* PIPS
   Authors: Miles Lubin
   See license and copyright information in the documentation */

#include "sFactoryAugSca.h"

#include "QpGenStochData.h"

#include "StochTree.h"
#include "StochInputTree.h"

#include "sLinsysRootAugSca.h"
#include "StochResourcePlanner.h"

sFactoryAugSca::sFactoryAugSca( StochInputTree* inputTree)
  : sFactory(inputTree)
{ 
  const int nb = 32; // blocking size
  MPI_Comm world = tree->commWrkrs;
  int scaprocs = StochResourcePlanner::noScaProcesses;
  int bcomm;
  
  int dims[2] = {0,0}, nprocs;  
  
  MPI_Comm_size(world, &nprocs);
  if (scaprocs <= 0) scaprocs = nprocs;
  
  MPI_Dims_create(scaprocs, 2, dims);
  
  FNAME(blacs_pinfo)(&bcomm,&bcomm); // need to call this to link for some reason
  bcomm = world;
  
  // grid is created on first scaprocs processors
  // the rest of the processors get noop
  
  FNAME(blacs_gridinit)(&bcomm, "Row-major", &dims[0], &dims[1]);
  
  
  initcomm(cinfo, world, bcomm, nb, dims[0], dims[1]); 
  if (cinfo.mype == 0) {
#ifdef TIMING
    printf("SCA NPROCS %d GRID %d %d BLOCKSIZE %d\n",
      scaprocs, dims[0], dims[1], nb);
#else
    printf("Using %d processes for ScaLAPACK, %d by %d grid, blocksize %d\n", 
      scaprocs, dims[0], dims[1], nb);
#endif
  }
  
  
};

sFactoryAugSca::sFactoryAugSca( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ )
  : sFactory(nx, my, mz, nnzQ, nnzA, nnzC)
{ };

sFactoryAugSca::sFactoryAugSca()
{ };

sFactoryAugSca::~sFactoryAugSca()
{ };


sLinsysRoot* sFactoryAugSca::newLinsysRoot()
{
  return new sLinsysRootAugSca(this, data, cinfo);
}

sLinsysRoot* 
sFactoryAugSca::newLinsysRoot(QpGenStochData* prob,
			   OoqpVector* dd,OoqpVector* dq,
			   OoqpVector* nomegaInv, OoqpVector* rhs)
{
  return new sLinsysRootAugSca(this, prob, dd, dq, nomegaInv, rhs, cinfo);
}
