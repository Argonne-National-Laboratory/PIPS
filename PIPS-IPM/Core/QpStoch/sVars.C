#include "sVars.h"
#include "QpGenData.h"
#include "OoqpVector.h"
#include "Data.h"
#include "QpGenResiduals.h"
#include "QpGenLinsys.h"
#include "StochVector.h"
#include "sTree.h"

#include "LinearAlgebraPackage.h"

#include <iostream> 
#include <fstream>
using namespace std;

sVars::sVars(sTree* tree, 
	     OoqpVector * ixlow_in, OoqpVector * ixupp_in,
	     OoqpVector * iclow_in, OoqpVector * icupp_in)
  : QpGenVars()
{
  /* not in use anymore */
  assert(false);
  SpReferTo( ixlow, ixlow_in );
  SpReferTo( ixupp, ixupp_in );
  SpReferTo( iclow, iclow_in );
  SpReferTo( icupp, icupp_in );

  assert( nx == ixlow->length() || 0 == ixlow->length() );
  nxlow = ixlow->numberOfNonzeros();
  
  assert( nx == ixlow->length() || 0 == ixlow->length() );
  nxupp = ixupp->numberOfNonzeros();

  assert( mz == iclow->length()  || 0 == iclow->length() );
  mclow = iclow->numberOfNonzeros();
  
  assert( mz == icupp->length() || 0 == icupp->length() );
  mcupp = icupp->numberOfNonzeros();


  s = OoqpVectorHandle( tree->newDualZVector() );
  if ( mclow > 0 ) {
    t      = OoqpVectorHandle( tree->newDualZVector() );
    lambda = OoqpVectorHandle( tree->newDualZVector() );
  } /*else {
    t      = OoqpVectorHandle( la->newVector( 0 ) );
    lambda = OoqpVectorHandle( la->newVector( 0 ) );
    }*/
  if( mcupp > 0 ) {
    u  = OoqpVectorHandle( tree->newDualZVector() );
    pi = OoqpVectorHandle( tree->newDualZVector() );
  } /*else {
    u  = OoqpVectorHandle( la->newVector( 0 ) );
    pi = OoqpVectorHandle( la->newVector( 0 ) );
    }*/
  if( nxlow > 0 ) {
    v     = OoqpVectorHandle( tree->newPrimalVector() );
    gamma = OoqpVectorHandle( tree->newPrimalVector() );
  }
	
  if( nxupp > 0 ) {
    w   = OoqpVectorHandle( tree->newPrimalVector() );
    phi = OoqpVectorHandle( tree->newPrimalVector() );
  } 
  x = OoqpVectorHandle( tree->newPrimalVector() );
  y = OoqpVectorHandle( tree->newDualYVector() );
  z = OoqpVectorHandle( tree->newDualZVector() );

  nComplementaryVariables = mclow + mcupp + nxlow + nxupp;

  createChildren();
}

sVars::sVars( sTree* tree, OoqpVector * x_in, OoqpVector * s_in,
	      OoqpVector * y_in, OoqpVector * z_in,
	      OoqpVector * v_in, OoqpVector * gamma_in,
	      OoqpVector * w_in, OoqpVector * phi_in,
	      OoqpVector * t_in, OoqpVector * lambda_in,
	      OoqpVector * u_in, OoqpVector * pi_in,
	      OoqpVector * ixlow_in, long long nxlowGlobal,
	      OoqpVector * ixupp_in, long long nxuppGlobal,
	      OoqpVector * iclow_in, long long mclowGlobal,
	      OoqpVector * icupp_in, long long mcuppGlobal)
  : QpGenVars()
{

  stochNode = tree;

  SpReferTo( x, x_in );
  SpReferTo( s, s_in );
  SpReferTo( y, y_in );
  SpReferTo( z, z_in );
  SpReferTo( v, v_in );
  SpReferTo( phi, phi_in );
  SpReferTo( w, w_in );
  SpReferTo( gamma, gamma_in );
  SpReferTo( t, t_in );
  SpReferTo( lambda, lambda_in );
  SpReferTo( u, u_in );
  SpReferTo( pi, pi_in );
  SpReferTo( ixlow, ixlow_in );
  SpReferTo( ixupp, ixupp_in );
  SpReferTo( iclow, iclow_in );
  SpReferTo( icupp, icupp_in );

  nx = x->length();
  my = y->length();
  mz = z->length();

  assert( nx == ixlow->length() || 0 == ixlow->length() );
  assert( nx == ixlow->length() || 0 == ixlow->length() );
  assert( mz == iclow->length() || 0 == iclow->length() );
  assert( mz == icupp->length() || 0 == icupp->length() );
  
  nxlow = nxlowGlobal; 
  nxupp = nxuppGlobal; 
  mclow = mclowGlobal;
  mcupp = mcuppGlobal;
  nComplementaryVariables = mclow + mcupp + nxlow + nxupp;

  assert( mz == s->length() );
  assert( nx == v     ->length() || ( 0 == v     ->length() && nxlow == 0 ));
  assert( nx == gamma ->length() || ( 0 == gamma ->length() && nxlow == 0 ));

  assert( nx == w     ->length() || ( 0 == w     ->length() && nxupp == 0 ));
  assert( nx == phi   ->length() || ( 0 == phi   ->length() && nxupp == 0 ));

  assert( mz == t     ->length() || ( 0 == t     ->length() && mclow == 0 ));
  assert( mz == lambda->length() || ( 0 == lambda->length() && mclow == 0 ));

  assert( mz == u     ->length() || ( 0 == u     ->length() && mcupp == 0 ));
  assert( mz == pi    ->length() || ( 0 == pi    ->length() && mcupp == 0 ));

  createChildren();
}

sVars::sVars(const sVars& vars) : QpGenVars(vars)
{
  stochNode = vars.stochNode;
  for(unsigned int i = 0; i < children.size(); ++i)
  {
    children.push_back( new sVars(*vars.children[i]));
  }
}

sVars::sVars(const sVars& vars, OoqpVectorHandle ixlow_, OoqpVectorHandle ixupp_,
      OoqpVectorHandle iclow_, OoqpVectorHandle icupp_) : QpGenVars(vars)
{
   stochNode = vars.stochNode;
   for(unsigned int i = 0; i < children.size(); ++i)
   {
      children.push_back( new sVars(*vars.children[i]));
   }

   ixlow = ixlow_;
   ixupp = ixupp_;
   iclow = iclow_;
   icupp = icupp_;
}

sVars::~sVars()
{ 
  for (size_t c = 0; c < children.size(); c++)
    delete children[c];
}

void sVars::AddChild(sVars* child)
{
  children.push_back(child);
}


void sVars::createChildren()
{
  StochVector& xst     = dynamic_cast<StochVector&>(*x);
  StochVector& sst     = dynamic_cast<StochVector&>(*s);
  StochVector& yst     = dynamic_cast<StochVector&>(*y); 
  StochVector& zst     = dynamic_cast<StochVector&>(*z);
  StochVector& vst     = dynamic_cast<StochVector&>(*v); 
  StochVector& gammast = dynamic_cast<StochVector&>(*gamma);
  StochVector& wst     = dynamic_cast<StochVector&>(*w); 
  StochVector& phist   = dynamic_cast<StochVector&>(*phi);
  StochVector& tst     = dynamic_cast<StochVector&>(*t); 
  StochVector& lambdast= dynamic_cast<StochVector&>(*lambda);
  StochVector& ust     = dynamic_cast<StochVector&>(*u); 
  StochVector& pist    = dynamic_cast<StochVector&>(*pi);
  StochVector& ixlowst = dynamic_cast<StochVector&>(*ixlow); 
  StochVector& ixuppst = dynamic_cast<StochVector&>(*ixupp); 
  StochVector& iclowst = dynamic_cast<StochVector&>(*iclow); 
  StochVector& icuppst = dynamic_cast<StochVector&>(*icupp); 
    

  for (size_t it=0; it<xst.children.size(); it++) {
    AddChild( new sVars( stochNode->children[it],
			 xst.children[it],     sst.children[it],
			 yst.children[it],     zst.children[it],
			 vst.children[it],     gammast.children[it],
			 wst.children[it],     phist.children[it],
			 tst.children[it],     lambdast.children[it],
			 ust.children[it],     pist.children[it],
			 ixlowst.children[it], nxlow,
			 ixuppst.children[it], nxupp,
			 iclowst.children[it], mclow,
			 icuppst.children[it], mcupp));
  }

}

void sVars::sync()
{
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*x));

  stochNode->syncDualYVector(dynamic_cast<StochVector&>(*y));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*z));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*s));

  if ( mclow > 0 ) {
    stochNode->syncDualZVector(dynamic_cast<StochVector&>(*t));
    stochNode->syncDualZVector(dynamic_cast<StochVector&>(*lambda));
  }
  if( mcupp > 0 ) {
    stochNode->syncDualZVector(dynamic_cast<StochVector&>(*u));
    stochNode->syncDualZVector(dynamic_cast<StochVector&>(*pi));
  }
  if( nxlow > 0 ) {
    stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*v));
    stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*gamma));
  }
  if( nxupp > 0 ) {
    stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*w));
    stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*phi));
  }
  //stochNode->syncDualZVector(dynamic_cast<StochVector&>(*));

}

bool sVars::isRootNodeInSync() const
{
   bool in_sync = true;
   const int my_rank = PIPS_MPIgetRank( MPI_COMM_WORLD );
   if( !dynamic_cast<const StochVector&>(*x).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "x not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*s).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "s not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*y).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "y not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*z).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "z not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*v).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "v not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*gamma).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "gamma not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*w).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "w not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*phi).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "phi not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*t).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "t not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*lambda).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "lambda not in sync" << std::endl;
      in_sync = false;
   }

   if( !dynamic_cast<const StochVector&>(*u).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "u not in sync" << std::endl;
      in_sync = false;
   }
   if( !dynamic_cast<const StochVector&>(*pi).isRootNodeInSync() )
   {
      if( my_rank == 0 )
         std::cout << "pi not in sync" << std::endl;
      in_sync = false;
   }

   MPI_Barrier(MPI_COMM_WORLD);
   return in_sync;
}
