#include "sResiduals.h"
#include "sTree.h"
#include "StochVector.h"

sResiduals::sResiduals( sTree* tree, 
			OoqpVector * rQ_,    OoqpVector * rA_, 
			OoqpVector * rC_,    OoqpVector * rz_, 
			OoqpVector * rt_,    OoqpVector * rlambda_, 
			OoqpVector * ru_,    OoqpVector * rpi_, 
			OoqpVector * rv_,    OoqpVector * rgamma_, 
			OoqpVector * rw_,    OoqpVector * rphi_, 
			OoqpVector * ixlow_, double nxlowGlobal,
			OoqpVector * ixupp_, double nxuppGlobal,
			OoqpVector * iclow_, double mclowGlobal, 
			OoqpVector * icupp_, double mcuppGlobal)
  :QpGenResiduals()
{
  assert(false);
  SpReferTo( ixlow, ixlow_ );
  nxlow = nxlowGlobal;

  SpReferTo( ixupp, ixupp_ );
  nxupp = nxuppGlobal;

  SpReferTo( iclow, iclow_ );
  mclow = mclowGlobal;

  SpReferTo( icupp, icupp_ );
  mcupp = mcuppGlobal;

  SpReferTo( rQ, rQ_ );
  SpReferTo( rA, rA_ );
  SpReferTo( rC, rC_ );
  SpReferTo( rz, rz_ );
  SpReferTo( rt     , rt_ );
  SpReferTo( rlambda, rlambda_ );
  SpReferTo( ru    , ru_ );
  SpReferTo( rpi   , rpi_ );
  SpReferTo( rv    , rpi_ );
  SpReferTo( rgamma, rgamma_ );
  SpReferTo( rw  , rw_ );
  SpReferTo( rphi, rphi_ );

  stochNode = tree;
  createChildren();
}


sResiduals::sResiduals( sTree* tree,
			OoqpVector * ixlow_, OoqpVector * ixupp_,
			OoqpVector * iclow_, OoqpVector * icupp_ )
  : QpGenResiduals()
{

  SpReferTo( ixlow, ixlow_ );
  nxlow = ixlow->numberOfNonzeros();

  SpReferTo( ixupp, ixupp_ );
  nxupp = ixupp->numberOfNonzeros();

  SpReferTo( iclow, iclow_ );
  mclow = iclow->numberOfNonzeros();

  SpReferTo( icupp, icupp_ );
  mcupp = icupp->numberOfNonzeros();

  rQ = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
  rA = OoqpVectorHandle( (OoqpVector*) tree->newDualYVector() );
  rC = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );

  rz = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
  if ( mclow > 0 ) {
    rt      = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
    rlambda = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
  } else {
    rt      = OoqpVectorHandle( (OoqpVector*) tree->newDualZVectorEmpty() );
    rlambda = OoqpVectorHandle( (OoqpVector*) tree->newDualZVectorEmpty() );
  }

  if ( mcupp > 0 ) {
    ru     = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
    rpi    = OoqpVectorHandle( (OoqpVector*) tree->newDualZVector() );
  } else {
    ru     = OoqpVectorHandle( (OoqpVector*) tree->newDualZVectorEmpty() );
    rpi    = OoqpVectorHandle( (OoqpVector*) tree->newDualZVectorEmpty() );
  }

  if( nxlow > 0 ) {
    rv     = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
    rgamma = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
  } else {
    rv     = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVectorEmpty() );
    rgamma = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVectorEmpty() );
  }

  if( nxupp > 0 ) {
    rw   = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
    rphi = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVector() );
  } else {
    rw   = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVectorEmpty() );
    rphi = OoqpVectorHandle( (OoqpVector*) tree->newPrimalVectorEmpty() );
  }
  
  stochNode = tree;
}

sResiduals::sResiduals( const sResiduals& res ) : QpGenResiduals( res )
{
   stochNode = res.stochNode;

   for(unsigned int i = 0; i < res.children.size(); ++i)
   {
       children.push_back( new sResiduals(*children[i]) );
   }
}

sResiduals::sResiduals( const sResiduals& res, OoqpVectorHandle ixlow_, OoqpVectorHandle ixupp_,
      OoqpVectorHandle iclow_, OoqpVectorHandle icupp_ ) : QpGenResiduals(res)
{
   stochNode = res.stochNode;

   for(unsigned int i = 0; i < res.children.size(); ++i)
   {
       children.push_back( new sResiduals(*children[i]) );
   }

   ixlow = ixlow_;
   ixupp = ixupp_;
   iclow = iclow_;
   icupp = icupp_;
}


void sResiduals::AddChild(sResiduals* child)
{
  children.push_back(child);
}

void sResiduals::createChildren()
{
  assert(false);
  StochVector& rQSt = dynamic_cast<StochVector&>(*rQ);

  StochVector& rASt = dynamic_cast<StochVector&>(*rA); 
  StochVector& rCSt = dynamic_cast<StochVector&>(*rC);    
  StochVector& rzSt = dynamic_cast<StochVector&>(*rz); 
  StochVector& rtSt = dynamic_cast<StochVector&>(*rt);    
  StochVector& rlambdaSt = dynamic_cast<StochVector&>(*rlambda); 
  StochVector& ruSt = dynamic_cast<StochVector&>(*ru);    
  StochVector& rpiSt = dynamic_cast<StochVector&>(*rpi); 
  StochVector& rvSt = dynamic_cast<StochVector&>(*rv);    
  StochVector& rgammaSt = dynamic_cast<StochVector&>(*rgamma); 
  StochVector& rwSt = dynamic_cast<StochVector&>(*rw);    
  StochVector& rphiSt = dynamic_cast<StochVector&>(*rphi); 
  StochVector& ixlowSt = dynamic_cast<StochVector&>(*ixlow);
  StochVector& ixuppSt = dynamic_cast<StochVector&>(*ixupp);
  StochVector& iclowSt = dynamic_cast<StochVector&>(*iclow);
  StochVector& icuppSt = dynamic_cast<StochVector&>(*icupp);

  //copy the structure of one of the vectors and create children of
  //this
  size_t nChildren=rQSt.children.size();
  for (size_t it=0; it<nChildren; it++) {

    assert(nChildren==stochNode->children.size());
    assert(nChildren==rASt.children.size());
    assert(nChildren==rzSt.children.size());
    assert(nChildren==rCSt.children.size());
    assert(nChildren==rtSt.children.size());
    assert(nChildren==rlambdaSt.children.size());
    assert(nChildren==ruSt.children.size());
    assert(nChildren==rpiSt.children.size());
    assert(nChildren==rvSt.children.size());
    assert(nChildren==rgammaSt.children.size());
    assert(nChildren==rwSt.children.size());
    assert(nChildren==rphiSt.children.size());
    assert(nChildren==ixlowSt.children.size());
    assert(nChildren==ixuppSt.children.size());
    assert(nChildren==iclowSt.children.size());
    assert(nChildren==icuppSt.children.size());
 
    AddChild(new sResiduals(stochNode->children[it], 
			    rQSt.children[it],    rASt.children[it], 
			    rCSt.children[it],    rzSt.children[it], 
			    rtSt.children[it],    rlambdaSt.children[it], 
			    ruSt.children[it],    rpiSt.children[it], 
			    rvSt.children[it],    rgammaSt.children[it], 
			    rwSt.children[it],    rphiSt.children[it], 
			    ixlowSt.children[it], nxlow,
			    ixuppSt.children[it], nxupp,
			    iclowSt.children[it], mclow, 
			    icuppSt.children[it], mcupp));
  }
}


void sResiduals::sync()
{
  //int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*rQ));
  //stochNode->displayVectorVsTreeStructure(dynamic_cast<StochVector&>(*rQ),
  //					  myRank);

  stochNode->syncDualYVector(dynamic_cast<StochVector&>(*rC));
  stochNode->syncDualYVector(dynamic_cast<StochVector&>(*rA));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*rz));

  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*rt));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*rlambda));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*ru));
  stochNode->syncDualZVector(dynamic_cast<StochVector&>(*rpi));

  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*rv));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*rgamma));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*rw));
  stochNode->syncPrimalVector(dynamic_cast<StochVector&>(*rphi));
}

void sResiduals::destroyChildren()
{
  //int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  for (size_t it=0; it<children.size(); it++) {
    //printf("CPU[%d] destroy %d\n", myRank, it);
    children[it]->destroyChildren(); 
  }
  
  for (size_t it=0; it<children.size(); it++) {
    
    delete children[it];
    //printf("deleted %d\n", it);
  }
  children.clear();
}
