/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGenResiduals.h"
#include "NlpGenVars.h"
#include "NlpGenData.h"

#include "OoqpVector.h"
#include "LinearAlgebraPackage.h"

#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#ifndef MAX
#define MAX(a,b) ( (a > b) ? a : b)
#endif
#ifndef MIN
#define MIN(a,b) ( (a > b) ? b : a)
#endif



NlpGenResiduals::NlpGenResiduals( LinearAlgebraPackage * la,
				long long nx_, long long my_, long long mz_,
				OoqpVector * ixlow_in, OoqpVector * ixupp_in,
				OoqpVector * iclow_in, OoqpVector * icupp_in )
{
  nx = nx_;
  my = my_;
  mz = mz_;

  SpReferTo( ixlow, ixlow_in );
  nxlow = ixlow->numberOfNonzeros();

  SpReferTo( ixupp, ixupp_in );
  nxupp = ixupp->numberOfNonzeros();

  SpReferTo( iclow, iclow_in );
  mclow = iclow->numberOfNonzeros();

  SpReferTo( icupp, icupp_in );
  mcupp = icupp->numberOfNonzeros();

  rQ = OoqpVectorHandle( la->newVector( nx ) );
  rA = OoqpVectorHandle( la->newVector( my ) );
  rC = OoqpVectorHandle( la->newVector( mz ) );

  rz = OoqpVectorHandle( la->newVector( mz ) );
  if ( mclow > 0 ) {
    rt      = OoqpVectorHandle( la->newVector( mz ) );
    rlambda = OoqpVectorHandle( la->newVector( mz ) );
  }
  if ( mcupp > 0 ) {
    ru     = OoqpVectorHandle( la->newVector( mz ) );
    rpi    = OoqpVectorHandle( la->newVector( mz ) );
  }
  if( nxlow > 0 ) {
    rv     = OoqpVectorHandle( la->newVector( nx ) );
    rgamma = OoqpVectorHandle( la->newVector( nx ) );
  }
  if( nxupp > 0 ) {
    rw   = OoqpVectorHandle( la->newVector( nx ) );
    rphi = OoqpVectorHandle( la->newVector( nx ) );
  }


  rp_OriSys = OoqpVectorHandle( la->newVector( nx ) );
  rd_OriSys = OoqpVectorHandle( la->newVector( my ) );

  priWrk = OoqpVectorHandle( la->newVector( nx ) );
  dualWrk = OoqpVectorHandle( la->newVector( my ) );
  dualWrk_Z = OoqpVectorHandle( la->newVector( mz ) );
  priWrk_S = OoqpVectorHandle( la->newVector( mz ) );

  Wd  = OoqpVectorHandle( la->newVector( nx ) );

  nrmrho = 0;
  nrmr =0;
  nrmrho0 = 0;
  psi = 0;
  mod0 = 0;
  mod = 0;
  dmod = 0;
  dWd = 0;
  thd = 0;
  dmodu = 0;
  comerr = 0;
}





void NlpGenResiduals::calcresids(Data *prob_in, Variables *vars_in)
{
  NlpGenVars * vars = (NlpGenVars *) vars_in;
  NlpGenData * prob = (NlpGenData *) prob_in;

  eval_kkt_res(prob,vars);
}
  

void NlpGenResiduals::add_r3_xz_alpha(Variables *vars_in, double alpha)
{
  NlpGenVars * vars = (NlpGenVars *) vars_in;

  if( mclow > 0 ) rlambda->axzpy( 1.0, *vars->t, *vars->lambda );
  if( mcupp > 0 ) rpi    ->axzpy( 1.0, *vars->u, *vars->pi );
  if( nxlow > 0 ) rgamma ->axzpy( 1.0, *vars->v, *vars->gamma );
  if( nxupp > 0 ) rphi   ->axzpy( 1.0, *vars->w, *vars->phi );

  if( alpha != 0.0 ) {
    if( mclow > 0 ) rlambda->addSomeConstants( alpha, *iclow );
    if( mcupp > 0 ) rpi    ->addSomeConstants( alpha, *icupp );
    if( nxlow > 0 ) rgamma ->addSomeConstants( alpha, *ixlow );
    if( nxupp > 0 ) rphi   ->addSomeConstants( alpha, *ixupp );
  }
}

void NlpGenResiduals::set_r3_xz_alpha(Variables *vars, double alpha)
{
  this->clear_r3();
  this->add_r3_xz_alpha( vars, alpha );
}
  
void NlpGenResiduals::clear_r3()
{
  if( mclow > 0 ) rlambda->setToZero();
  if( mcupp > 0 ) rpi    ->setToZero();
  if( nxlow > 0 ) rgamma ->setToZero();
  if( nxupp > 0 ) rphi   ->setToZero();
}
  
void NlpGenResiduals::clear_r1r2()
{
  rQ->setToZero();
  rA->setToZero();
  rC->setToZero();
  rz->setToZero();
  if( nxlow > 0 ) rv->setToZero();
  if( nxupp > 0 ) rw->setToZero();
  if( mclow > 0 ) rt->setToZero();
  if( mcupp > 0 ) ru->setToZero();
}

void NlpGenResiduals::project_r3(double rmin, double rmax)
{
  if( mclow > 0 ) {
    rlambda->gondzioProjection( rmin, rmax );
    rlambda->selectNonZeros( *iclow );
  }
  if( mcupp > 0 ) {
    rpi    ->gondzioProjection( rmin, rmax );
    rpi    ->selectNonZeros( *icupp );
  }
  if( nxlow > 0 ) {
    rgamma ->gondzioProjection( rmin, rmax );
    rgamma ->selectNonZeros( *ixlow );
  }
  if( nxupp > 0 ) {
    rphi   ->gondzioProjection( rmin, rmax );
    rphi   ->selectNonZeros( *ixupp );
  }

}
  

int NlpGenResiduals::validNonZeroPattern()
{
  if( nxlow > 0 && 
      ( !rv    ->matchesNonZeroPattern( *ixlow ) ||
	!rgamma->matchesNonZeroPattern( *ixlow ) ) ) {
    return 0;
  }

  if( nxupp > 0 &&
      ( !rw  ->matchesNonZeroPattern( *ixupp ) ||
	!rphi->matchesNonZeroPattern( *ixupp ) ) ) {
    return 0;
  }
  if( mclow > 0 &&
      ( !rt     ->matchesNonZeroPattern( *iclow ) ||
	!rlambda->matchesNonZeroPattern( *iclow ) ) ) {
    return 0;
  }

  if( mcupp > 0 &&
      ( !ru ->matchesNonZeroPattern( *icupp ) ||
	!rpi->matchesNonZeroPattern( *icupp ) ) ) {
    return 0;
  }
  
  return 1;
}

NlpGenResiduals::~NlpGenResiduals()
{
}


void NlpGenResiduals::addDampingTermToOneSidePart(const double DampingTerm)
{
  /////////////////////////    add damping term, see paper section 3.7   /////////////////////////
  if( nxlow > 0 ) {
    //use priWrk as temp vector, set value -1
    priWrk->setToConstant(-1);
    priWrk->selectNonZeros( *ixlow );

    if(nxupp>0){
	  //find index where var have both bounds, val = -1
	  priWrk->selectNonZeros( *ixupp);
	  //find index where var only have lower bound, now val=1 where only lb, otherwise 0
	  priWrk->addSomeConstants(1, *ixlow);
    }else{
	  priWrk->negate();
	}
    rQ->axpy(DampingTerm,*priWrk);
  }
  
  if( nxupp > 0 ) {
    //use priWrk as temp vector, set value -1
    priWrk->setToConstant(-1);
    priWrk->selectNonZeros( *ixupp );

    if(nxlow>0){
	  //find index where var have both bounds, val = -1
	  priWrk->selectNonZeros( *ixlow);
	  //find index where var only have upper bound, now val=1 where only ub, otherwise 0
	  priWrk->addSomeConstants(1, *ixupp);
    }else{
	  priWrk->negate();
	}
    rQ->axpy(-DampingTerm,*priWrk);
  }

  if( mclow > 0 ) {
    //use priWrk_S as temp vector, set value -1
    priWrk_S->setToConstant(-1);
    priWrk_S->selectNonZeros( *iclow );

    if(mcupp>0){
	  //find index where var have both bounds, val = -1
	  priWrk_S->selectNonZeros( *icupp);
	  //find index where var only have lower bound, now val=1 where only lb, otherwise 0
	  priWrk_S->addSomeConstants(1, *iclow);
    }else{
	  priWrk_S->negate();
	}
    rz->axpy(DampingTerm,*priWrk_S);
  }

  if( mcupp > 0 ) {
    //use priWrk_S as temp vector, set value -1
    priWrk_S->setToConstant(-1);
    priWrk_S->selectNonZeros( *icupp );

    if(mclow>0){
	  //find index where var have both bounds, val = -1
	  priWrk_S->selectNonZeros( *iclow);
	  //find index where var only have upper bound, now val=1 where only ub, otherwise 0
	  priWrk_S->addSomeConstants(1, *icupp);
    }else{
	  priWrk_S->negate();
	}
    rz->axpy(-DampingTerm,*priWrk_S);
  }
//////////////////////////////////////////////////////////////////////////////////////////////////

}



void NlpGenResiduals::eval_kkt_res(NlpGenData *prob, NlpGenVars *vars)
{

  double componentNorm, prinorm=0.0, dualnorm=0.0, comnorm=0.0, gap=0.0;
  
  prob->getg( *rQ );  
  // calculate x^T (g+Qx) - contribution to the duality gap
  gap = rQ->dotProductWith(*vars->x);   
  
  // calculate dual residual: g - Ae'*lame - Ai'*lami
  prob->ATransmult( 1.0, *rQ, -1.0, *vars->y ); 
  prob->CTransmult( 1.0, *rQ, -1.0, *vars->z );

  // add dual variables for the var boundary constraints
  if( nxlow > 0 ) 
  	rQ->axpy( -1.0, *vars->gamma );
  if( nxupp > 0 ) 
  	rQ->axpy(  1.0, *vars->phi );

  componentNorm = rQ->infnorm();
  if( componentNorm > dualnorm ) 
  	dualnorm = componentNorm;
  // calculate primal residual for Eq Cons: Ce(x)-b
  rA->copyFrom(*prob->CeqBody);
  rA->axpy( -1.0, *prob->bA);
  // contribution -d^T y to duality gap
  gap -= prob->bA->dotProductWith(*vars->y);

  componentNorm = rA->infnorm();
  if( componentNorm > prinorm ) 
  	prinorm = componentNorm;

  //compute rc = Ci(x)-s
  // additional eq constraints due to the constraint bounds:  C(x)=s
  rC->copyFrom(*prob->CIneqBody);  
  rC->axpy(-1, *vars->s );
  
  componentNorm = rC->infnorm();
  if( componentNorm > prinorm ) 
  	prinorm = componentNorm;
  // calculate primal residual for InEq Cons: clow <=  s <= cupp
  // also calculate dual residual for the slack, as z-lambda, lambda is the dual for the lower bound and pi for the upper
  rz->setToZero();
  if( mclow + mcupp > 0 ) 
  	rz->copyFrom( *vars->z );

  if( mclow > 0 ) {
    rz->axpy( -1.0, *vars->lambda ); // z-lambda

    rt->copyFrom( *vars->s ); 
    rt->axpy( -1.0, prob->slowerBound() );
    rt->selectNonZeros( *iclow );
    rt->axpy( -1.0, *vars->t ); // residual: s-clow-t, t is the slack for the lower bound

	gap -= prob->bl->dotProductWith(*vars->lambda);
	
    componentNorm = rt->infnorm();
    if( componentNorm > prinorm ) 
	  prinorm = componentNorm;
  }
  
  if( mcupp > 0 ) { 
    rz->axpy(  1.0, *vars->pi );   // z-(lambda)+pi

    ru->copyFrom( *vars->s );
    ru->axpy( -1.0, prob->supperBound() );
    ru->selectNonZeros( *icupp );
    ru->axpy( 1.0, *vars->u );  // residual: s-cupp+u, u is the slack for the upper bound
    
    gap += prob->bu->dotProductWith(*vars->pi);

    componentNorm = ru->infnorm();
    if( componentNorm > prinorm ) 
	  prinorm = componentNorm;
  }

  componentNorm = rz->infnorm();
  if( componentNorm > dualnorm ) 
  	dualnorm = componentNorm;
  // additional eq constraints due to the variable bounds:  x-v=xlow, x+w=xupp
  // if nxlow == prob->nxLOri, we do not transfer problem
  if( nxlow > 0 ) {

    rv->copyFrom( *vars->x );
    rv->axpy( -1.0, prob->xlowerBound() );
    rv->selectNonZeros( *ixlow );
    rv->axpy( -1.0, *vars->v );

    gap -= prob->blx->dotProductWith(*vars->gamma);

    componentNorm = rv->infnorm();
    if( componentNorm > prinorm ) 
	  prinorm = componentNorm;
  }
  if( nxupp > 0 ) {
    rw->copyFrom( *vars->x );
    rw->axpy( -1.0, prob->xupperBound() );
    rw->selectNonZeros( *ixupp );
    rw->axpy(  1.0, *vars->w );

    gap += prob->bux->dotProductWith(*vars->phi);

    componentNorm = rw->infnorm();
    if( componentNorm > prinorm ) 
	  prinorm = componentNorm;
  }
  mDualityGap = gap; //no dual gap for nonlinear prob


  perr = prinorm;
  derr = dualnorm;  

  mResidualNorm = perr>derr?perr:derr;

 // evaluate complementariy conditions residual
  this->clear_r3();
  this->add_r3_xz_alpha( vars, 0.0 ); 
  if( mclow > 0 ){ 
  	componentNorm = rlambda->infnorm();
    if( componentNorm > comnorm) 
	  comnorm = componentNorm;
  }
  if( mcupp > 0 ){ 
  	componentNorm = rpi->infnorm();
	if( componentNorm > comnorm) 
	  comnorm = componentNorm;
  }
  if( nxlow > 0 ){
	componentNorm = rgamma->infnorm();
	if( componentNorm > comnorm) 
	  comnorm = componentNorm;
  }
  if( nxupp > 0 ){
	componentNorm = rphi->infnorm();
	if( componentNorm > comnorm) 
	  comnorm = componentNorm;
  }

  comerr_0 = comnorm;

  //!petra
  //valgrind complains that 'comerr' has uninitialised value, which is true
  //mResidualNorm = mResidualNorm>comerr?mResidualNorm:comerr; 
  //I think the above line should be
  mResidualNorm = mResidualNorm>comerr_0?mResidualNorm:comerr_0;  

  //compute XZe-mu  -> comp err here
  comerr = getKKTError_Comp(prob,vars,prob->currMu,PIPS_NORM_INF);

}







double NlpGenResiduals::priErr()
{
  return perr;
}

double NlpGenResiduals::dualErr()
{
  return derr;
}

double NlpGenResiduals::OriPriErr()
{
  return Oriperr;
}

double NlpGenResiduals::OriDualErr()
{
  return Oriderr;
}

double NlpGenResiduals::comp_Err()
{
  return comerr;
}

double NlpGenResiduals::comp_Err_0()
{
  return comerr_0;
}




//compute primal resdual for dual vars y and z
double 
NlpGenResiduals::getKKTRhsNorm_Primal(NlpGenData *prob, NlpGenVars *vars, 
						const int normType, const int isTrialStep)
{
  double prinorm;

  // get Ce(x) and Ci(x)
  if(isTrialStep==0){
    dualWrk->copyFrom(*prob->CeqBody);
	dualWrk_Z->copyFrom(*prob->CIneqBody);
  }else{
  	dualWrk->copyFrom(*prob->trialCeqBody);
	dualWrk_Z->copyFrom(*prob->trialCIneqBody);	
  }

  // calculate primal residual for Eq Cons: Ce(x)-b , do partial diff for y  
  dualWrk->axpy( -1.0, *prob->bA);

  //compute rc = Ci(x)-s
  dualWrk_Z->axpy(-1,*vars->s);
  
  prinorm = dualWrk_Z->Norm(normType,dualWrk);


  // FIXME_NY: do we need to include the additional constraints due to the slacks?
  // prinorm = resid->priErr();


  return prinorm;
}


double 
NlpGenResiduals::getKKTRhsNorm_Dual(NlpGenData *prob, NlpGenVars *vars,
						const int normType,  const int isTrialStep)
{  
  double dualNorm;

  if(isTrialStep==0)
	prob->getg( *priWrk );
  else 
  	assert("not implemented" && 0); 

  prob->ATransmult( 1.0, *priWrk, -1.0, *vars->y );
  prob->CTransmult( 1.0, *priWrk, -1.0, *vars->z );  
  
  // add dual variables for the var boundary constraints
  if( nxlow > 0 ) 
  	priWrk->axpy( -1.0, *vars->gamma );
  if( nxupp > 0 ) 
  	priWrk->axpy(  1.0, *vars->phi );  

  priWrk_S->copyFrom(*vars->z);
  if( mclow > 0 ) {
    priWrk_S->axpy( -1.0, *vars->lambda ); // z-lambda
  }
  if( mcupp > 0 ) { 
    priWrk_S->axpy(  1.0, *vars->pi );   // z-(lambda)+pi
  }

  dualNorm = priWrk_S->Norm(normType,priWrk);

  return dualNorm;
}

double 
NlpGenResiduals::getKKTError_Comp(NlpGenData *prob, NlpGenVars *vars, const double mu, 
						const int normType,  const int isTrialStep)
{
  double componentNorm1=0.0, componentNorm2=0.0;
  
  // evaluate complementariy conditions residual
  this->clear_r3();
  this->add_r3_xz_alpha( vars, -mu ); 
  if( mclow > 0 ){ 
  	componentNorm1 = rlambda->Norm(normType);
  }
  if( mcupp > 0 ){ 
  	componentNorm2 = rpi->Norm(normType);
	componentNorm1 = rpi->Norm(normType,componentNorm1,componentNorm2);
  }
  if( nxlow > 0 ){
	componentNorm2 = rgamma->Norm(normType);
	componentNorm1 = rgamma->Norm(normType,componentNorm1,componentNorm2);
  }
  if( nxupp > 0 ){
	componentNorm2 = rphi->Norm(normType);
	componentNorm1 = rphi->Norm(normType,componentNorm1,componentNorm2);
  }

  return componentNorm1;
}


void
NlpGenResiduals::copyFrom(NlpGenResiduals *residual_in)
{
  rQ->copyFrom(*residual_in->rQ); 
  rA->copyFrom(*residual_in->rA);
  rC->copyFrom(*residual_in->rC);

  rz->copyFrom(*residual_in->rz);

  if ( mclow > 0 ) {
    rt->copyFrom(*residual_in->rt);
    rlambda->copyFrom(*residual_in->rlambda);
  }
  if ( mcupp > 0 ) {
    ru->copyFrom(*residual_in->ru);
    rpi->copyFrom(*residual_in->rpi);
  }
  if ( nxlow > 0 ) {
    rv->copyFrom(*residual_in->rv);
    rgamma->copyFrom(*residual_in->rgamma);
  }
  if ( nxupp > 0 ) {
    rw->copyFrom(*residual_in->rw);
    rphi->copyFrom(*residual_in->rphi);
  }  
}



void
NlpGenResiduals::updateSOCRhs(const double AlphaStep, NlpGenVars *vars_in, NlpGenData *prob_in)
{
  rA->scale(AlphaStep);
  rA->axpy(1,*prob_in->trialCeqBody);
  rA->axpy(-1,*prob_in->bA);

  rC->scale(AlphaStep);
  rC->axpy(1,*prob_in->trialCIneqBody);
  rC->axpy(-1,*vars_in->s);
}



bool
NlpGenResiduals::findSmallStep(NlpGenVars *vars, NlpGenVars *steps, const double tol_mach)
{
  bool testSmall;
  double smallNorm;

  priWrk->absVal(vars->x);
  priWrk->addConstant(1);
  priWrk->componentDiv(*steps->x);
  priWrk->invert();

  priWrk_S->absVal(vars->s);  
  priWrk_S->addConstant(1);
  priWrk_S->componentDiv(*steps->s);
  priWrk_S->invert();

  smallNorm = priWrk->Norm(PIPS_NORM_INF,priWrk_S);
  testSmall = (smallNorm<(10*tol_mach));
  
  return testSmall;
}


