/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGenVars.h"
#include "NlpGenData.h"
#include "OoqpVector.h"
#include "Data.h"
#include "NlpGenResiduals.h"
#include "NlpGenLinsys.h"
#include "SimpleVector.h"
#include "MpsReader.h"

#include "LinearAlgebraPackage.h"

#include <iostream> 
#include <fstream>
#include <stdlib.h>
using namespace std;

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif


NlpGenVars::NlpGenVars( OoqpVector * x_in, OoqpVector * s_in,
		      OoqpVector * y_in, OoqpVector * z_in,
		      OoqpVector * v_in, OoqpVector * gamma_in,
		      OoqpVector * w_in, OoqpVector * phi_in,
		      OoqpVector * t_in, OoqpVector * lambda_in,
		      OoqpVector * u_in, OoqpVector * pi_in,
		      OoqpVector * ixlow_in, OoqpVector * ixupp_in,
		      OoqpVector * iclow_in, OoqpVector * icupp_in )
{
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
  assert( mz == iclow->length()  || 0 == iclow->length() );
  assert( mz == icupp->length() || 0 == icupp->length() );
  
  nxlow = ixlow->numberOfNonzeros();
  nxupp = ixupp->numberOfNonzeros();
  mclow = iclow->numberOfNonzeros();
  mcupp = icupp->numberOfNonzeros();
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

  nSlack = nComplementaryVariables;
}



NlpGenVars::NlpGenVars( LinearAlgebraPackage * la,
		      int nx_, int my_, int mz_,
		      OoqpVector * ixlow_in, OoqpVector * ixupp_in,
		      OoqpVector * iclow_in, OoqpVector * icupp_in )
{
  SpReferTo( ixlow, ixlow_in );
  SpReferTo( ixupp, ixupp_in );
  SpReferTo( iclow, iclow_in );
  SpReferTo( icupp, icupp_in );

  nx = nx_;
  my = my_;
  mz = mz_;

  assert( nx == ixlow->length() || 0 == ixlow->length() );
  nxlow = ixlow->numberOfNonzeros();
  
  assert( nx == ixlow->length() || 0 == ixlow->length() );
  nxupp = ixupp->numberOfNonzeros();

  assert( mz == iclow->length()  || 0 == iclow->length() );
  mclow = iclow->numberOfNonzeros();
  
  assert( mz == icupp->length() || 0 == icupp->length() );
  mcupp = icupp->numberOfNonzeros();
  
  s = OoqpVectorHandle( la->newVector( mz ) );
  if ( mclow > 0 ) {
    t      = OoqpVectorHandle( la->newVector( mz ) );
    lambda = OoqpVectorHandle( la->newVector( mz ) );
  } else {
    t      = OoqpVectorHandle( la->newVector( 0 ) );
    lambda = OoqpVectorHandle( la->newVector( 0 ) );
  }
  if( mcupp > 0 ) {
    u  = OoqpVectorHandle( la->newVector( mz ) );
    pi = OoqpVectorHandle( la->newVector( mz ) );
  } else {
    u  = OoqpVectorHandle( la->newVector( 0 ) );
    pi = OoqpVectorHandle( la->newVector( 0 ) );
  }
  if( nxlow > 0 ) {
    v     = OoqpVectorHandle( la->newVector( nx ) );
    gamma = OoqpVectorHandle( la->newVector( nx ) );
  } else {
    v     = OoqpVectorHandle( la->newVector( 0 ) );
    gamma = OoqpVectorHandle( la->newVector( 0 ) );
  }
	
  if( nxupp > 0 ) {
    w   = OoqpVectorHandle( la->newVector( nx ) );
    phi = OoqpVectorHandle( la->newVector( nx ) );
  } else {
    w   = OoqpVectorHandle( la->newVector( 0 ) );
    phi = OoqpVectorHandle( la->newVector( 0 ) );
  }

  x = OoqpVectorHandle( la->newVector( nx ) );
  y = OoqpVectorHandle( la->newVector( my ) );
  z = OoqpVectorHandle( la->newVector( mz ) );
  nComplementaryVariables = mclow + mcupp + nxlow + nxupp;

  nSlack = nComplementaryVariables;

  PiLx = OoqpVectorHandle( la->newVector( nx ) );
  PiUx = OoqpVectorHandle( la->newVector( nx ) );
  PiLs = OoqpVectorHandle( la->newVector( mz ) );
  PiUs = OoqpVectorHandle( la->newVector( mz ) );  

}

double NlpGenVars::mu()
{
  double mu  = 0.0;
  if ( nComplementaryVariables == 0 ) {
    return 0.0;
  } else {
    if( mclow > 0 ) mu += t->dotProductWith( *lambda );
    if( mcupp > 0 ) mu += u->dotProductWith( *pi );
    if( nxlow > 0 ) mu += v->dotProductWith( *gamma );
    if( nxupp > 0 ) mu += w->dotProductWith( *phi );

    mu /= nComplementaryVariables;

    return mu;
  }
}

double NlpGenVars::mustep(Variables * step_in, double alpha)
{
  NlpGenVars * step = (NlpGenVars *) step_in;
  double mu  = 0.0;
  if ( nComplementaryVariables == 0 ) {
    return 0.0;
  } else {
    if( mclow > 0 ) {
      mu += t->shiftedDotProductWith( alpha, *step->t,
				      *lambda,
				      alpha, *step->lambda );
    }
    if( mcupp > 0 ) {
      mu += u->shiftedDotProductWith( alpha, *step->u,
				      *pi,
				      alpha, *step->pi );
    }
    if( nxlow > 0 ) {
      mu += v->shiftedDotProductWith( alpha, *step->v,
				      *gamma,
				      alpha, *step->gamma );
    }
    if( nxupp > 0 ) {
      mu += w->shiftedDotProductWith( alpha, *step->w,
				      *phi,
				      alpha, *step->phi );
    }
    mu /= nComplementaryVariables;
    return mu;
  }
}

void NlpGenVars::saxpy( Variables *b_in, double alpha )
{
  NlpGenVars * b = (NlpGenVars *) b_in;

  x->axpy( alpha, *b->x );
  y->axpy( alpha, *b->y );
  z->axpy( alpha, *b->z );
  s->axpy( alpha, *b->s );
  if( mclow > 0 ) {
    assert( b->t     ->matchesNonZeroPattern( *iclow ) &&
	    b->lambda->matchesNonZeroPattern( *iclow ) );

    t     ->axpy( alpha, *b->t );
    lambda->axpy( alpha, *b->lambda );
  }
  if( mcupp > 0 ) {
    assert( b->u     ->matchesNonZeroPattern( *icupp ) &&
	    b->pi    ->matchesNonZeroPattern( *icupp ) );

    u     ->axpy( alpha, *b->u );
    pi    ->axpy( alpha, *b->pi );
  }
  if( nxlow > 0 ) {
    assert( b->v     ->matchesNonZeroPattern( *ixlow ) &&
	    b->gamma ->matchesNonZeroPattern( *ixlow ) );

    v     ->axpy( alpha, *b->v );
    gamma ->axpy( alpha, *b->gamma );
  }
  if( nxupp > 0 ) {
    assert( b->w     ->matchesNonZeroPattern( *ixupp ) &&
	    b->phi   ->matchesNonZeroPattern( *ixupp ) );

    w     ->axpy( alpha, *b->w );
    phi   ->axpy( alpha, *b->phi );
  }
}

void NlpGenVars::negate()
{
  s->negate();
  x->negate();
  y->negate();
  z->negate();
  if( mclow > 0 ) {
    t     ->negate();
    lambda->negate();
  }
  if( mcupp > 0 ) {
    u     ->negate();
    pi    ->negate();
  }
  if( nxlow > 0 ) {
    v     ->negate();
    gamma ->negate();
  }
  if( nxupp > 0 ) {
    w     ->negate();
    phi   ->negate();
  }
}



double NlpGenVars::stepMax_Pri( Variables * step_in, const double tau )
{
  NlpGenVars * sstep = (NlpGenVars *) step_in;
  double priMaxStep = 1.0/tau;
  
  if( mclow > 0 ) {
    assert( t     ->somePositive( *iclow ) );

    priMaxStep = t     ->stepbound( *sstep->t, priMaxStep );
  }

  if( mcupp > 0 ) {
    assert( u ->somePositive( *icupp ) );

    priMaxStep = u ->stepbound( *sstep->u,  priMaxStep );
  }
  
  if( nxlow > 0 ) {
    assert( v    ->somePositive( *ixlow ) );

    priMaxStep = v    ->stepbound( *sstep->v,     priMaxStep );
  }

  if( nxupp > 0 ) {
    assert( w  ->somePositive( *ixupp ) );

    priMaxStep = w  ->stepbound( *sstep->w,   priMaxStep );
  }

  priMaxStep  *= tau;

  return priMaxStep;
}






double NlpGenVars::stepMax_BoundDual( Variables * step_in, const double tau )
{
  NlpGenVars * sstep = (NlpGenVars *) step_in;
  double dualMaxStep = 1.0/tau;

  if( mclow > 0 ) {
    assert( lambda->somePositive( *iclow ) );

    dualMaxStep = lambda->stepbound( *sstep->lambda, dualMaxStep );
  }

  if( mcupp > 0 ) {
    assert( pi->somePositive( *icupp ) );

    dualMaxStep = pi->stepbound( *sstep->pi, dualMaxStep );
  }
  
  if( nxlow > 0 ) {
    assert( gamma->somePositive( *ixlow ) );

    dualMaxStep = gamma->stepbound( *sstep->gamma, dualMaxStep );
  }

  if( nxupp > 0 ) {
    assert( phi->somePositive( *ixupp ) );

    dualMaxStep = phi->stepbound( *sstep->phi, dualMaxStep );
  }

  dualMaxStep *= tau;

  return dualMaxStep;
}



double NlpGenVars::stepbound( Variables * b_in )
{
  NlpGenVars * b = (NlpGenVars *) b_in;
  double maxStep;

  double tau;

  
  maxStep = 1.0;

  if( mclow > 0 ) {
    assert( t     ->somePositive( *iclow ) );
    assert( lambda->somePositive( *iclow ) );

    maxStep = t     ->stepbound( *b->t, maxStep );
    maxStep = lambda->stepbound( *b->lambda, maxStep );
  }

  if( mcupp > 0 ) {
    assert( u ->somePositive( *icupp ) );
    assert( pi->somePositive( *icupp ) );

    maxStep = u ->stepbound( *b->u,  maxStep );
    maxStep = pi->stepbound( *b->pi, maxStep );
  }
  
  if( nxlow > 0 ) {
    assert( v    ->somePositive( *ixlow ) );
    assert( gamma->somePositive( *ixlow ) );

    maxStep = v    ->stepbound( *b->v,     maxStep );
    maxStep = gamma->stepbound( *b->gamma, maxStep );
  }

  if( nxupp > 0 ) {
    assert( w  ->somePositive( *ixupp ) );
    assert( phi->somePositive( *ixupp ) );

    maxStep = w  ->stepbound( *b->w,   maxStep );
    maxStep = phi->stepbound( *b->phi, maxStep );
  }

  return maxStep;
}
  
int NlpGenVars::isInteriorPoint()
{
  int interior = 1;
  if( mclow > 0 ) {
    interior = interior &&
      t     ->somePositive( *iclow ) &&
      lambda->somePositive( *iclow );
  }

  if( mcupp > 0 ) {
    interior = interior &&
      u ->somePositive( *icupp ) &&
      pi->somePositive( *icupp );
  }
  
  if( nxlow > 0 ) {
    interior = interior &&
      v    ->somePositive( *ixlow ) &&
      gamma->somePositive( *ixlow );
  }

  if( nxupp > 0 ) {
    interior = interior &&
      w  ->somePositive( *ixupp ) &&
      phi->somePositive( *ixupp );
  }

  return interior;
}

double NlpGenVars::findBlocking( Variables * step, 
				double & primalValue,
				double & primalStep,
				double & dualValue,
				double & dualStep,
				int& firstOrSecond )
{
  double alpha = 1.0;
  firstOrSecond = 0;

  NlpGenVars * d = (NlpGenVars *) step;

  if( mclow > 0 ) {
    alpha = t->findBlocking( *d->t, *lambda, *d->lambda, alpha,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,	
			     firstOrSecond );
  }

  if( mcupp > 0 ) {
    alpha = u->findBlocking( *d->u, *pi, *d->pi, alpha,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,
			     firstOrSecond );
  }
  
  if( nxlow > 0 ) {
    alpha = v->findBlocking( *d->v, *gamma, *d->gamma, alpha,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,
			     firstOrSecond );
  }

  if( nxupp > 0 ) {
    alpha = w->findBlocking( *d->w, *phi, *d->phi, alpha,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,
			     firstOrSecond );
  }

  return alpha;
}




void NlpGenVars::findBlockingPriDual( Variables * step, 
				double & primalValue,
				double & primalStep,
				double & dualValue,
				double & dualStep,
				int& firstOrSecond, 
				double tau, double & alphaPri, double & alphaDual )
{

  alphaPri  = 1.0/tau;
  alphaDual = 1.0/tau;
  firstOrSecond = 0;

  NlpGenVars * d = (NlpGenVars *) step;

  if( mclow > 0 ) {
     t->findBlockingPD( *d->t, *lambda, *d->lambda,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,	
			     firstOrSecond, &alphaPri, &alphaDual );
  }

  if( mcupp > 0 ) {
    u->findBlockingPD( *d->t, *lambda, *d->lambda,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,	
			     firstOrSecond, &alphaPri, &alphaDual );
  }
  
  if( nxlow > 0 ) {
    v->findBlockingPD( *d->t, *lambda, *d->lambda,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,	
			     firstOrSecond, &alphaPri, &alphaDual );
  }

  if( nxupp > 0 ) {
    w->findBlockingPD( *d->t, *lambda, *d->lambda,
			     &primalValue, &primalStep,
			     &dualValue, &dualStep,	
			     firstOrSecond, &alphaPri, &alphaDual );
  }

  alphaPri  *= tau;
  alphaDual *= tau;

}


void NlpGenVars::interiorPoint( double alpha, double beta )
{
  s->setToZero();
  x->setToZero();
  y->setToZero();
  z->setToZero();

  if( nxlow > 0 ) {
    v     ->setToConstant ( alpha );
    v     ->selectNonZeros( *ixlow );
    gamma ->setToConstant ( beta  );
    gamma ->selectNonZeros( *ixlow );
  }
  if( nxupp > 0 ) {
    w     ->setToConstant ( alpha );
    w     ->selectNonZeros( *ixupp );
    phi   ->setToConstant ( beta  );
    phi   ->selectNonZeros( *ixupp );
  }

  if( mclow > 0 ) {
    t      ->setToConstant ( alpha );
    t      ->selectNonZeros( *iclow );
    lambda ->setToConstant ( beta  );
    lambda ->selectNonZeros( *iclow );
  }
  if( mcupp > 0 ) {
    u      ->setToConstant ( alpha );
    u      ->selectNonZeros( *icupp );
    pi     ->setToConstant ( beta  );
    pi     ->selectNonZeros( *icupp );
  }
}
    
double NlpGenVars::violation()
{
  double viol = 0.0, cmin = 0.0;
  int iblock;

  if( nxlow > 0 ) {
    v->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
    
    gamma->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
  }
  if( nxupp > 0 ) {
    w->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
    
    phi->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
  }
  if( mclow > 0 ) {
    t->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
    
    lambda->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
  }
  if( mcupp > 0 ) {
    u->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
    
    pi->min( cmin, iblock );
    if( cmin < viol ) viol = cmin;
  }
  return -viol;
}

void NlpGenVars::shiftBoundVariables( double alpha, double beta )
{
  if( nxlow > 0 ) {
    v     ->addSomeConstants( alpha, *ixlow );
    gamma ->addSomeConstants( beta,  *ixlow );
  }
  if( nxupp > 0 ) {
    w     ->addSomeConstants( alpha, *ixupp );
    phi   ->addSomeConstants( beta,  *ixupp );
  }
  if( mclow > 0 ) {
    t     ->addSomeConstants( alpha, *iclow );
    lambda->addSomeConstants( beta,  *iclow );
  }
  if( mcupp > 0 ) {
    u     ->addSomeConstants( alpha, *icupp );
    pi    ->addSomeConstants( beta,  *icupp );
  }
}


void NlpGenVars::print()
{
  x->print();
//  x->writeToStream(cout);
//  x->writefToStream( cout, "x[%{index}] = %{value}" );
}

// this is copyfrom  
void NlpGenVars::copy(Variables *b_in)
{
  NlpGenVars * b = (NlpGenVars *) b_in;
  
  s->copyFrom( *b->s );
  if( nxlow > 0 ) {
    v->copyFrom( *b->v );
    gamma ->copyFrom( *b->gamma );
  }
  if( nxupp > 0 ) {
    w->copyFrom( *b->w );
    phi   ->copyFrom( *b->phi );
  }
  if( mclow > 0 ) {
    t->copyFrom( *b->t );
    lambda->copyFrom( *b->lambda );
  }
  if( mcupp > 0 ) {
    u->copyFrom( *b->u );
    pi    ->copyFrom( *b->pi );
  }
  x->copyFrom( *b->x );
  y->copyFrom( *b->y );
  z->copyFrom( *b->z );
  
}

//copy Dual variable and set primal to zero
void NlpGenVars::setPriZero()
{
  s->setToZero();
  if( nxlow > 0 ) {
    v->setToZero();
  }
  if( nxupp > 0 ) {
    w->setToZero();
  }
  if( mclow > 0 ) {
    t->setToZero();
  }
  if( mcupp > 0 ) {
    u->setToZero();
  }
  x->setToZero();
}

void NlpGenVars::copyDualPri0(Variables *b_in)
{
  NlpGenVars * b = (NlpGenVars *) b_in;
  
  s->setToZero();
  if( nxlow > 0 ) {
    v->setToZero();
    gamma ->copyFrom( *b->gamma );
  }
  if( nxupp > 0 ) {
    w->setToZero();
    phi   ->copyFrom( *b->phi );
  }
  if( mclow > 0 ) {
    t->setToZero();
    lambda->copyFrom( *b->lambda );
  }
  if( mcupp > 0 ) {
    u->setToZero();
    pi    ->copyFrom( *b->pi );
  }
  x->setToZero();
  y->copyFrom( *b->y );
  z->copyFrom( *b->z );
  
}


double NlpGenVars::onenorm()
{
  double norm;
  norm  = x->onenorm();
  norm += s->onenorm();
  norm += y->onenorm();
  norm += z->onenorm();

  norm += v->onenorm();
  norm += phi->onenorm();
  norm += w->onenorm();
  norm += gamma->onenorm();
  norm += t->onenorm();
  norm += lambda->onenorm();
  norm += u->onenorm();
  norm += pi->onenorm();

  return norm;
}


double NlpGenVars::infnorm()
{
  double norm, temp;
  norm = 0.0;

  temp  = x->infnorm();
  if(temp > norm) norm = temp;
  temp = s->infnorm();
  if(temp > norm) norm = temp;
  temp = y->infnorm();
  if(temp > norm) norm = temp;
  temp = z->infnorm();
  if(temp > norm) norm = temp;

  temp = v->infnorm();
  if(temp > norm) norm = temp;
  temp = phi->infnorm();
  if(temp > norm) norm = temp;

  temp = w->infnorm();
  if(temp > norm) norm = temp;
  temp = gamma->infnorm();
  if(temp > norm) norm = temp;

  temp = t->infnorm();
  if(temp > norm) norm = temp;
  temp = lambda->infnorm();
  if(temp > norm) norm = temp;

  temp = u->infnorm();
  if(temp > norm) norm = temp;
  temp = pi->infnorm();
  if(temp > norm) norm = temp;

  return norm;
}

NlpGenVars::~NlpGenVars()
{
}

int NlpGenVars::validNonZeroPattern()
{
  if( nxlow > 0 && 
      ( !v    ->matchesNonZeroPattern( *ixlow ) ||
	!gamma->matchesNonZeroPattern( *ixlow ) ) ) {
    return 0;
  }

  if( nxupp > 0 &&
      ( !w  ->matchesNonZeroPattern( *ixupp ) ||
	!phi->matchesNonZeroPattern( *ixupp ) ) ) {
    return 0;
  }
  if( mclow > 0 &&
      ( !t     ->matchesNonZeroPattern( *iclow ) ||
	!lambda->matchesNonZeroPattern( *iclow ) ) ) {
    return 0;
  }

  if( mcupp > 0 &&
      ( !u ->matchesNonZeroPattern( *icupp ) ||
	!pi->matchesNonZeroPattern( *icupp ) ) ) {
    return 0;
  }
  
  return 1;
}

void NlpGenVars::unscaleSolution(NlpGenData * data)
{

// Modifying sx is equivalent to modifying x
SimpleVector & sx = (SimpleVector &) *this->x;

// x = D * x'
sx.componentMult(data->scale());
}

void NlpGenVars::unscaleBounds(NlpGenData * data)
{

SimpleVector & sxlow = (SimpleVector &) data-> xlowerBound();
SimpleVector & sxupp = (SimpleVector &) data-> xupperBound();

// l = D * l' 
sxlow.componentMult(data->scale());

// u = D * u'
sxupp.componentMult(data->scale());
}

void NlpGenVars::printSolution( MpsReader * reader, NlpGenData * data,
			       int& iErr )
{
  assert( x->isKindOf( kSimpleVector ) ); // Otherwise this routine
  // cannot be used.
  double objective;
  {
    SimpleVectorHandle temp( new SimpleVector(nx) );
    data->getg( *temp );
    data->Qmult( 1.0, *temp, 0.5, *x );
    objective = temp->dotProductWith( *x );
  }

  SimpleVector & sx      = (SimpleVector &) *this->x;
  SimpleVector & sxlow   = (SimpleVector &) data-> xlowerBound();
  SimpleVector & sixlow  = (SimpleVector &) data->ixlowerBound();
  SimpleVector & sxupp   = (SimpleVector &) data-> xupperBound();
  SimpleVector & sixupp  = (SimpleVector &) data->ixupperBound();
  SimpleVector & sgamma  = (SimpleVector &) *this->gamma;
  SimpleVector & sphi    = (SimpleVector &) *this->phi;
  SimpleVector & sy      = (SimpleVector &) *this->y;
  SimpleVector & ss      = (SimpleVector &) *this->s;
  SimpleVector & slambda = (SimpleVector &) *this->lambda;
  SimpleVector & spi     = (SimpleVector &) *this->pi;
  SimpleVector & sz      = (SimpleVector &) *this->z;
  SimpleVector & sclow   = (SimpleVector &) data-> slowerBound();
  SimpleVector & siclow  = (SimpleVector &) data->islowerBound();
  SimpleVector & scupp   = (SimpleVector &) data-> supperBound();
  SimpleVector & sicupp  = (SimpleVector &) data->isupperBound();

  char * cxupp = new char[nx];
  char * cxlow = new char[nx];
  for( int j = 0; j < nx; j++ ) {
    if( nxupp > 0 && sixupp[j] != 0 ) {
      cxupp[j] = 1;
    } else {
      cxupp[j] = 0;
    }
    if( nxlow > 0 && sixlow[j] != 0 ) {
      cxlow[j] = 1;
    } else {
      cxlow[j] = 0;
    }
  }
  char *cclow, *ccupp;
  if( mz <= 0 ) {
    cclow = 0; ccupp = 0;
  } else {
    cclow = new char[mz];
    ccupp = new char[mz];
    for( int i = 0; i < mz; i++ ) {
      if( mclow > 0 && siclow[i] != 0.0 ) {
	cclow[i] = 1;
      } else {
	cclow[i] = 0;
      }
      if( mcupp > 0 && sicupp[i] != 0.0 ) {
	ccupp[i] = 1;
      } else {
	ccupp[i] = 0;
      }
    }
  }

  if( reader->scalingOption == 1){
      // Unscale the solution and bounds before printing
      this->unscaleSolution( data);
      this->unscaleBounds( data);
      }

  reader->printSolution( sx.elements(), nx,
			 sxlow.elements(), cxlow, sxupp.elements(), cxupp,
			 sgamma.elements(), sphi.elements(),
			 sy.elements(), my,
			 ss.elements(), mz,
			 sclow.elements(), cclow,
			 scupp.elements(), ccupp,
			 slambda.elements(), spi.elements(),
			 sz.elements(), 
			 objective,
			 iErr );
  delete [] cclow;
  delete [] ccupp;
  delete [] cxlow;
  delete [] cxupp;
}


void NlpGenVars::interiorPointPriX( double alpha)
{
  x->setToZero();
  x->setToConstant(alpha);
}

void NlpGenVars::interiorPointPriS( double alpha)
{
  s->setToZero();
  s->setToConstant(alpha);
}

void NlpGenVars::interiorPointDualY( double alpha)
{
  y->setToZero();
  y->setToConstant( alpha);
}

void NlpGenVars::interiorPointDualZ( double alpha)
{
  z->setToZero();
  z->setToConstant( alpha);
}

void NlpGenVars::interiorBoundSlack( double alpha)
{
  if( mclow > 0 ) {
	t->setToZero();
	t->setToConstant(alpha);
	t->selectNonZeros(*iclow);	
  }
  if( mcupp > 0 ) {
  	u->setToZero();
	u->setToConstant(alpha);
	u->selectNonZeros(*icupp);	
  }
  if( nxlow > 0 ) {
    v->setToZero();
	v->setToConstant(alpha);
	v->selectNonZeros(*ixlow);	
  }
  if( nxupp > 0 ) {
    w->setToZero();
	w->setToConstant(alpha);
	w->selectNonZeros(*ixupp);	
  }
}


void NlpGenVars::interiorBoundSlackDual( double alpha)
{
  if( mclow > 0 ) {
	lambda->setToZero();
	lambda->setToConstant(alpha);
	lambda->selectNonZeros(*iclow);
  }
  if( mcupp > 0 ) {
  	pi->setToZero();
	pi->setToConstant(alpha);
	pi->selectNonZeros(*icupp);
  }
  if( nxlow > 0 ) {
    gamma->setToZero();
	gamma->setToConstant(alpha);
	gamma->selectNonZeros(*ixlow);
  }
  if( nxupp > 0 ) {
    phi->setToZero();
	phi->setToConstant(alpha);
	phi->selectNonZeros(*ixupp);
  }
}


//copy Dual variable and set primal to zero
void NlpGenVars::setZero()
{
  s->setToZero();
  if( nxlow > 0 ) {
    v->setToZero();
    gamma ->setToZero();
  }
  if( nxupp > 0 ) {
    w->setToZero();
    phi   ->setToZero();
  }
  if( mclow > 0 ) {
    t->setToZero();
    lambda->setToZero();
  }
  if( mcupp > 0 ) {
    u->setToZero();
    pi    ->setToZero();
  }
  x->setToZero();
  y->setToZero();
  z->setToZero();
  
}




void NlpGenVars::takeStep( Variables *step_in, 
				const double alphaPri, const double alphaDualY, const double alphaNu, const int onlyPrimal)
{
  takePrimalStep(step_in,alphaPri,alphaPri);
  if(onlyPrimal == 0)
    takeDualStep(step_in,alphaDualY,alphaNu);
}

void NlpGenVars::SetSlackFromMaxXorY( OoqpVector* v_out, 
				OoqpVector* x_in, OoqpVector *y_in)
{
   v->SetComponentFromMaxXorY(x,y_in,nx-nSlack, nx,nx-nSlack, nx, my-nSlack,my);
}

void NlpGenVars::takePrimalStep( Variables *step_in, const double alphaPri, double alphaSlack)
{
  NlpGenVars * steps = (NlpGenVars *) step_in;

  if(alphaSlack==-1)
  	alphaSlack = alphaPri;
 
  //% primal update	 
  x->axpy( alphaPri, *steps->x );
  s->axpy( alphaSlack, *steps->s );

  //% nu update
  if( mclow > 0 ) {
    assert( steps->t->matchesNonZeroPattern( *iclow ));
    t ->axpy( alphaSlack, *steps->t );
  }
  if( mcupp > 0 ) {
    assert( steps->u->matchesNonZeroPattern( *icupp ));
    u->axpy( alphaSlack, *steps->u );
  }
  if( nxlow > 0 ) {
    assert( steps->v->matchesNonZeroPattern( *ixlow ));
    v->axpy( alphaSlack, *steps->v );
  }
  if( nxupp > 0 ) {
    assert( steps->w->matchesNonZeroPattern( *ixupp ));
    w->axpy( alphaSlack, *steps->w );
  }
}


void NlpGenVars::takeDualStep( Variables *step_in, const double alphaDual, double alphaSlackDual)
{
  NlpGenVars * steps = (NlpGenVars *) step_in;

  if(alphaSlackDual==-1)
  	alphaSlackDual = alphaDual;
 
  //% dual update   
  y->axpy( alphaDual, *steps->y );
  z->axpy( alphaDual, *steps->z );

  //% nu update
  if( mclow > 0 ) {
    assert( steps->lambda->matchesNonZeroPattern( *iclow ));
    lambda->axpy( alphaSlackDual, *steps->lambda );
  }
  if( mcupp > 0 ) {
    assert( steps->pi->matchesNonZeroPattern( *icupp ));
    pi->axpy( alphaSlackDual, *steps->pi );
  }
  if( nxlow > 0 ) {
    assert( steps->gamma->matchesNonZeroPattern( *ixlow ));
    gamma->axpy( alphaSlackDual, *steps->gamma );
  }
  if( nxupp > 0 ) {
    assert( steps->phi->matchesNonZeroPattern( *ixupp ));
    phi->axpy( alphaSlackDual, *steps->phi );
  }
}


void
NlpGenVars::updateSlackAndDual( OoqpVector *tempx,OoqpVector *temps,const double k_sigma, const double mu)
{

  //reset  slack and its corresponding dual var
  if( nxlow > 0 ) {

	tempx->setToConstant(k_sigma*mu);
	tempx->selectNonZeros(*ixlow);
	tempx->divideSome(*v,*ixlow);
	gamma->MinComponentBetween(tempx,ixlow);

	tempx->setToConstant(mu/k_sigma);
	tempx->selectNonZeros(*ixlow);
	tempx->divideSome(*v,*ixlow);
	gamma->SetComponentFromMaxXorY(tempx,gamma,ixlow);
  }
  
  if( nxupp > 0 ) {

	tempx->setToConstant(k_sigma*mu);
	tempx->selectNonZeros(*ixupp);
	tempx->divideSome(*w,*ixupp);
	phi->MinComponentBetween(tempx,ixupp);

	tempx->setToConstant(mu/k_sigma);
	tempx->selectNonZeros(*ixupp);
	tempx->divideSome(*w,*ixupp);
	phi->SetComponentFromMaxXorY(tempx,phi,ixupp);
  }


  if( mclow > 0 ) {

	temps->setToConstant(k_sigma*mu);
	temps->selectNonZeros(*iclow);
	temps->divideSome(*t,*iclow);
	lambda->MinComponentBetween(temps,iclow);

	temps->setToConstant(mu/k_sigma);
	temps->selectNonZeros(*iclow);
	temps->divideSome(*t,*iclow);
	lambda->SetComponentFromMaxXorY(temps,lambda,iclow);
  }
  
  if( mcupp > 0 ) {

	temps->setToConstant(k_sigma*mu);
	temps->selectNonZeros(*icupp);
	temps->divideSome(*u,*icupp);
	pi->MinComponentBetween(temps,icupp);

	temps->setToConstant(mu/k_sigma);
	temps->selectNonZeros(*icupp);
	temps->divideSome(*u,*icupp);
	pi->SetComponentFromMaxXorY(temps,pi,icupp);

  }
  
}


  
double NlpGenVars::primal_XS_InfNorm(bool XSonly)
{
  double result, normWrk;

  result  = x->Norm(PIPS_NORM_INF);
  normWrk = s->Norm(PIPS_NORM_INF);
  if(result<normWrk) result=normWrk;
  
  if(XSonly)
  	return result;
  	
  if( nxlow > 0 ) {
    normWrk = v->Norm(PIPS_NORM_INF);
	if(result<normWrk) result = normWrk;
  }
  if( nxupp > 0 ) {
    normWrk = w->Norm(PIPS_NORM_INF);
    if(result<normWrk) result = normWrk;
  }
  if( mclow > 0 ) {
    normWrk = t->Norm(PIPS_NORM_INF);
    if(result<normWrk) result = normWrk;
  }
  if( mcupp > 0 ) {
    normWrk = u->Norm(PIPS_NORM_INF);
    if(result<normWrk) result = normWrk;
  }

  return result;
}

double NlpGenVars::dual_YZ_InfNorm(bool YZonly)
{
  double result=0, normWrk=0;
  
  result  = y->Norm(PIPS_NORM_INF);
  normWrk = z->Norm(PIPS_NORM_INF);
  if(result<normWrk) result = normWrk; 

  if(YZonly)
  	return result;

  if( nxlow > 0 ) {
	normWrk = gamma->Norm(PIPS_NORM_INF);
	if(result<normWrk) result = normWrk;
  }
  if( nxupp > 0 ) {
	normWrk = phi->Norm(PIPS_NORM_INF);
	if(result<normWrk) result = normWrk;
  }
  if( mclow > 0 ) {
	normWrk = lambda->Norm(PIPS_NORM_INF);
	if(result<normWrk) result = normWrk;
  }
  if( mcupp > 0 ) {
	normWrk = pi->Norm(PIPS_NORM_INF);
	if(result<normWrk) result = normWrk;
  }

  return result;
}




void
NlpGenVars::push_variables( OoqpVector *vec, OoqpVector *vec_slackLB, OoqpVector *vec_slackUB,
					OoqpVector *vec_Lb,OoqpVector *vec_Ub, 
					OoqpVector *vec_Temp, const double k_1, const double k_2, const int ifX)
{
  OoqpVector *PiL, *PiU;
  OoqpVector *idxLow, *idxUp;
  int nUp, nLow;

  // To avoid round-off error, move variables first at the bounds
  if (k_1>0.0 || k_2>0.0) {
	push_variables(vec, vec_slackLB, vec_slackUB, vec_Lb, vec_Ub, 
				vec_Temp,0.0,0.0, ifX);
  }

  if(ifX==1){
  	PiL=&(*PiLx);
	PiU=&(*PiUx);
	idxLow = &(*ixlow);
	idxUp = &(*ixupp);
	nUp = nxupp; 
	nLow = nxlow;
  }
  else{
   	PiL=&(*PiLs);
	PiU=&(*PiUs);
	idxLow = &(*iclow);
	idxUp = &(*icupp);	
	nUp = mcupp; 
	nLow = mclow;	
  }
  
  if(k_1>0.0||k_2>0.0){
    // Calculate p_l = k_1 * max(|lb|,1)
    PiL->absVal(vec_Lb);
    PiL->SetComponentFromMaxXorConstant(PiL,1.0);
    PiL->scale(k_1);

    // Calculate p_u = k_1 * max(|ub|,1)
    PiU->absVal(vec_Ub);
    PiU->SetComponentFromMaxXorConstant(PiU,1.0);
    PiU->scale(k_1);  

    //set vec_Temp=k_2*(ub-lb) with two bounds
    vec_Temp->copyFrom(*vec_Ub);
    vec_Temp->axpy(-1,*vec_Lb);
    vec_Temp->scale(k_2);
    vec_Temp->selectNonZeros(*idxLow);
    vec_Temp->selectNonZeros(*idxUp);
  
    // Calculate  PiL =
    //  min(k_1 * max(|lb|,1), k_2*(ub-lb))  for components 	with two bounds
    //  k_1 * max(|lb|,1) 						   	otherwise
    PiL->MinComponentBetween(vec_Temp,idxUp);
    PiL->selectNonZeros(*idxLow);

    // Calculate PiU =
    //  min(k_1 * max(|ub|,1), k_2*(ub-lb))  for components	with two bounds
    //  k_1 * max(|ub|,1) 						   otherwise
    PiU->MinComponentBetween(vec_Temp,idxLow);
    PiU->selectNonZeros(*idxUp);
  }
  else{
    // set PiL=0 and PiU=0
	PiL->setToZero();
	PiU->setToZero();
  }


  // Calculate  x
  /*
  	only  lower bound:
		x <--  max{ x, lb + PiL} = x + max{ 0, lb + PiL -x}
	only upper bound:
		x <--  min{ x, ub - PiU} = x + min{ 0, ub - PiU -x} = x - max{ 0, x - (ub-PiU)}  
     	with both bounds:
     		piL(i) = min(k_1*max(1,abs(lb(i))),k_2*(ub(i)-lb(i)));
       	piU(i) = min(k_1*max(1,abs(ub(i))),k_2*(ub(i)-lb(i)));

	x	=  x + max{ 0, lb + PiL -x}  - max{ 0, x - (ub-PiU)} 
      	   	=  x + max{ 0, lb + PiL -x} + min{ 0, ub - PiU -x} 
  */

  // Calculate  PiL <- max{ 0, lb + PiL -x}
  PiL->axpy(1,*vec_Lb);
  PiL->axpy(-1,*vec);
  PiL->SetComponentFromMaxXorConstant(PiL,0);
  PiL->selectNonZeros(*idxLow);
  
  // Calculate  PiU <- min{ 0, ub - PiU -x}
  PiU->negate();
  PiU->axpy(1,*vec_Ub);
  PiU->axpy(-1,*vec);
  PiU->MinComponentOrConstant(PiU,0);
  PiU->selectNonZeros(*idxUp);
  
  // Calculate	x  +=  max{ 0, lb + PiL -x} + min{ 0, ub - PiU -x} 
  vec->axpy(1,*PiL);
  vec->axpy(1,*PiU);



  if(k_1>0.0||k_2>0.0){
    // let vec_slackLB = x-lb
    if(nLow>0){
  	  vec_slackLB->copyFrom(*vec);
	  vec_slackLB->axpy(-1,*vec_Lb);
	  vec_slackLB->selectNonZeros(*idxLow);
    }
  
    // let vec_slackUB = ub-x
    if(nUp>0){
  	  vec_slackUB->copyFrom(*vec_Ub);
	  vec_slackUB->axpy(-1,*vec);
	  vec_slackUB->selectNonZeros(*idxUp);
    }
  }
}

void NlpGenVars::getErrScaling(double s_max, double &scal_commerr, double &scal_dualerr)
{
  long long n = nxlow + nxupp + mclow + mcupp;

  assert(n==nComplementaryVariables && "these two should be equal");
  //printf("nx %d nxlow %d nxupp %d    mz %d mclow %d mcupp %d     my %d\n", nx, nxlow, nxupp, mz, mclow, mcupp, my);
  scal_commerr = lambda->onenorm() + pi->onenorm() + gamma->onenorm() + phi->onenorm();

  if (n == 0) {
	scal_commerr = 1.0;
  }
  else {
	scal_commerr = scal_commerr / n;
	scal_commerr = MAX(s_max, scal_commerr)/s_max;
  }

  n  += my + mz;
  scal_dualerr = lambda->onenorm() + pi->onenorm() + gamma->onenorm() + phi->onenorm()+ y->onenorm() + z->onenorm();

  if ( n == 0 ) {
	scal_dualerr = 1.0;
  }
  else {
	scal_dualerr = scal_dualerr / n;
	scal_dualerr = MAX(s_max, scal_dualerr)/s_max;
  }

  assert(scal_commerr>0 && scal_dualerr>0);
}

double NlpGenVars::computeDD()
{
  double result=0.0;

  result += x->dotProductWith( *x );
  if(mz>0) result += s->dotProductWith( *s );
//  if(my>0) result += y->dotProductWith( *y );
//  if(mz>0) result += z->dotProductWith( *z );
  
  if(nxlow){ 
  	result += v->dotProductWith( *v );
//	result += gamma->dotProductWith( *gamma );
  }
  if(nxupp){ 
  	result += w->dotProductWith( *w);
//	result += phi->dotProductWith( *phi );
  }  
  
  if(mclow){ 
  	result += t->dotProductWith( *t );
//	result += lambda->dotProductWith( *lambda );
  }
  if(mcupp){ 
  	result += u->dotProductWith( *u);
//	result += pi->dotProductWith( *pi );
  } 


  return result;
}

double NlpGenVars::computeXSDD(Variables * step_in)
{
  NlpGenVars * steps = (NlpGenVars *) step_in;

  double result=0.0;

  result += steps->x->dotProductWith( *steps->x );
  if(steps->mz>0) result += steps->s->dotProductWith( *steps->s );


/*  
  if(nxlow){ 
  	result += v->dotProductWith( *v );
  }
  if(nxupp){ 
  	result += w->dotProductWith( *w);
  }  
  
  if(mclow){ 
  	result += t->dotProductWith( *t );
  }
  if(mcupp){ 
  	result += u->dotProductWith( *u);
  } 
*/

  return result;
}


void NlpGenVars::mergeNTstep(Variables * Nstep_in, Variables * Tstep_in, Variables * Curr_Iter_in, Residuals * res_in)
{
  NlpGenVars * Nsteps = (NlpGenVars *) Nstep_in;
  NlpGenVars * Tsteps = (NlpGenVars *) Tstep_in;
  NlpGenVars * Cur_Iter = (NlpGenVars *) Curr_Iter_in;
  NlpGenResiduals * res   = (NlpGenResiduals *) res_in;  

  this->copy(Tsteps);

  // for primal step  d = n + t
  x->axpy(1, *Nsteps->x );
  if(mz>0) s->axpy(1, *Nsteps->s );  

/*
  // for dual step d = t - curr_iter  
  if(mz>0) z->axpy(-1, *Cur_Iter->z );
  if(my>0) y->axpy(-1, *Cur_Iter->y );
*/

  if( mclow > 0 ) {
    t->copyFrom(*s );
    t->axpy( -1.0, *res->rt );
    t->selectNonZeros( *iclow );

    lambda->copyFrom( *res->rlambda );
    lambda->axzpy( -1.0, *Cur_Iter->lambda, *t );
    lambda->divideSome( *Cur_Iter->t, *iclow );
  }
  if( mcupp > 0 ) {
    u->copyFrom( *res->ru );
    u->axpy( -1.0, *s );
    u->selectNonZeros( *icupp );

    pi->copyFrom( *res->rpi );
    pi->axzpy( -1.0, *Cur_Iter->pi, *u );
    pi->divideSome( *Cur_Iter->u, *icupp );
  }
  if( nxlow > 0 ) {
    v->copyFrom( *x );
    v->axpy( -1.0, *res->rv );
    v->selectNonZeros( *ixlow );
	
    gamma->copyFrom( *res->rgamma );
    gamma->axzpy( -1.0, *Cur_Iter->gamma, *v );
    gamma->divideSome( *Cur_Iter->v, *ixlow );
  }
  if( nxupp > 0 ) {
    w->copyFrom( *res->rw );
    w->axpy( -1.0, *x );
    w->selectNonZeros( *ixupp );
	
    phi->copyFrom( *res->rphi );
    phi->axzpy( -1.0, *Cur_Iter->phi, *w );
    phi->divideSome( *Cur_Iter->w, *ixupp );
  }

  

}


