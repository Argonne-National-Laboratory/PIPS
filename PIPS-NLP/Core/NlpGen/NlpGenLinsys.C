/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGenLinsys.h"

#include "NlpGenData.h"
#include "NlpGenResiduals.h"
#include "NlpGenVars.h"

#include "OoqpVector.h"
#include "DoubleMatrix.h"
#include "DoubleLinearSolver.h"
#include "SimpleVector.h"
#include "LinearAlgebraPackage.h"
#include "NlpGen.h"
#include "RegularizationAlg.h"

#include "SparseGenMatrix.h"

#include <stdlib.h>


#ifdef TIMING
#include "mpi.h"
#include <vector>
#endif

#include <fstream>
using namespace std;

#include "constants.h"
extern PreCondInfo *preCond;

extern int gOuterSolve;
extern int gOuterBiCGIter;
extern int separateHandDiag;

extern int gDoIR_Aug;
extern int gDoIR_Full;

extern int gMaxIR;
extern double gIRtol;


extern int gdWd_test;	

extern int gUsePetsc;
extern int gUser_Defined_PC;
extern int gUsePetscOuter;

extern int gUseDualRegAlg;


#ifndef MAX
#define MAX(a,b) ( (a > b) ? a : b)
#endif
#ifndef MIN
#define MIN(a,b) ( (a > b) ? b : a)
#endif

#ifdef WITH_PETSC
#include "petscksp.h"

typedef struct {
  NlpGenLinsys	*NlpGenLinsysPtr;
  NlpGenData	*NlpGenDataPtr;
  OoqpVector 	*solx_;
  OoqpVector 	*sols_;
  OoqpVector 	*soly_;
  OoqpVector 	*solz_;
} Mat_Shell_Petsc;

typedef struct {
  NlpGenLinsys	*NlpGenLinsysPtr;
  NlpGenData	*NlpGenDataPtr;  
} Shell_PC_Petsc;


extern PetscErrorCode _user_MatMultXSYZ(Mat,Vec,Vec);
extern PetscErrorCode PC_Shell_Create(Shell_PC_Petsc **shell);
extern PetscErrorCode PC_Shell_SetUp(PC pc,   NlpGenLinsys	*NlpGenLinsysPtr_in, NlpGenData	*NlpGenDataPtr_in);
extern PetscErrorCode PC_Shell_Destroy(PC pc);
extern PetscErrorCode PC_Shell_Apply_Sparse(PC pc,Vec b_in,Vec x_out);

#endif



NlpGenLinsys::NlpGenLinsys( NlpGen * factory_,
			  NlpGenData * prob,
			  LinearAlgebraPackage * la ) :
  factory( factory_), rhs(NULL), dd(NULL), dq(NULL), useRefs(0),additiveDiag(NULL)
{

  nx = prob->nx; my = prob->my; mz = prob->mz;
  ixlow = prob->ixlow;
  ixupp = prob->ixupp;
  iclow = prob->iclow;
  icupp = prob->icupp;

  nxlow = ixlow->numberOfNonzeros();
  nxupp = ixupp->numberOfNonzeros();
  mclow = iclow->numberOfNonzeros();
  mcupp = icupp->numberOfNonzeros();

  num_slacks = nxlow + nxupp + mclow + mcupp;

  if( nxupp + nxlow >= 0 ) {
    dd      	=  la->newVector( nx ) ;		// diagQ+XZ^{-1}+reg		or XZ^{-1}+reg	if split
    dq      	=  la->newVector( nx ) ;		// diagQ
	temp_diagX 	=  la->newVector( nx ) ;		// diagQ+XZ^{-1}			or XZ^{-1}		if split
    prob->getDiagonalOfQ( *dq );
  }
  nomegaInv   =  la->newVector( mz ) ;
  
  if(gOuterSolve<3) {
  	rhs         =  la->newVector( nx + my + mz ) ;
  }
  else if(gOuterSolve>=3 && separateHandDiag==1){
  	rhs         =  la->newVector( nx + my + mz +mz) ;
	additiveDiag =  la->newVector( nx + my + mz +mz) ;
  }  
  else if(gOuterSolve>=3 && separateHandDiag==0){
  	rhs         =  la->newVector( nx + my + mz +mz) ;
  }else{
	assert(0);
  }
  
  temp_diagS =  la->newVector( mz ) ;
  temp_diagZ =  la->newVector( mz ) ;
  temp_diagY =  la->newVector( my ) ;

  sol2Bicg = res2Bicg = res3Bicg = res4Bicg = res5Bicg = NULL;
  sol2 = res2 = res3 = res4 = res5 = NULL;

  if(gOuterSolve<3) {
    //for iterative refinement or BICGStab
    sol  = la->newVector( nx + my + mz );
    res  = la->newVector( nx + my + mz );
    resx = la->newVector( nx );
    resy = la->newVector( my );
    resz = la->newVector( mz );
    if(gOuterSolve==2) {
      //BiCGStab for compressed; additional vectors needed
      sol2 = la->newVector( nx + my + mz );
      res2 = la->newVector( nx + my + mz );
      res3  = la->newVector( nx + my + mz );
      res4  = la->newVector( nx + my + mz );
      res5  = la->newVector( nx + my + mz );
    }
  } else if(gOuterSolve>=3) {
	sol  = la->newVector( nx + my + mz + mz);
	res  = la->newVector( nx + my + mz + mz);
	resx = la->newVector( nx );
	ress = la->newVector( mz );
	resy = la->newVector( my );
	resz = la->newVector( mz );
	if(gDoIR_Aug ==1 || gUsePetscOuter!=0){
	  sol2 = la->newVector( nx + my + mz + mz );
	  res2 = la->newVector( nx + my + mz + mz );
	  res3 = res4 = res5 = NULL;
	}
    if(gOuterSolve>=4) {
      //BiCGStab; additional vectors needed
      sol2Bicg   = la->newVector( nx + my + mz );
      res2Bicg  = la->newVector( nx + my + mz );
      res3Bicg  = la->newVector( nx + my + mz );
      res4Bicg  = la->newVector( nx + my + mz );
      res5Bicg  = la->newVector( nx + my + mz );
    }	
  }
  else{
    sol = res = resx = ress = resy = resz = NULL;
  }
  
  priReg=0.0;
  dualReg=0.0;

  fullQ = ((gUsePetsc==0) || (gUsePetscOuter==0) || (gUser_Defined_PC!=2));

  KryIter=0;
}


// build the matrix which need to be factorized, in QP form
void NlpGenLinsys::factorNoMatChange(Data *prob_in, Variables *vars_in,RegularizationAlg *RegInfo)
{
  
  int updateDiag=1;
  int updateMatAndDiag=2;
  int updateMethod;

  if(RegInfo==NULL){
    this->NlpGenLinsys::factor( prob_in, vars_in);
	this->UpdateMatrices(prob_in,updateMatAndDiag);
  }
  else{
  	if(RegInfo->newSystem){
  	  updateMethod = updateMatAndDiag;
	  this->NlpGenLinsys::factor( prob_in, vars_in);
  	}else{
	  updateMethod = updateDiag;  	
  	}
	this->NlpGenLinsys::factor( prob_in, vars_in, RegInfo);
	this->UpdateMatrices(prob_in, updateMethod);
  }
  
}

void NlpGenLinsys::factor(Data *  prob_in, Variables *vars_in,RegularizationAlg *RegInfo)
{
  NlpGenVars * vars = (NlpGenVars *) vars_in;
  NlpGenData * prob = dynamic_cast<NlpGenData *>(prob_in);
  priReg=0.0;
  dualReg=0.0;

  assert( vars->validNonZeroPattern() );

   // move to the upper level
  priReg = RegInfo->prim_reg_curr;
  dualReg = RegInfo->dual_reg_curr;
  
  if( nxlow + nxupp >= 0 ){
  	//temp_diagX is computed as the diag val without regularization
	dd->copyFrom(*temp_diagX);
	dd->addConstant(priReg);
  }
  
  if( my > 0 ){
  	temp_diagY->setToConstant(-dualReg);
  }

  if(mclow + mcupp > 0 ){
    if( gOuterSolve<3 ){ 
	  //temp_diagS is computed as the diag val without regularization for slack s. Note that this need to be invert and negate
	  nomegaInv->copyFrom(*temp_diagS);
	  nomegaInv->addConstant(priReg);
	  nomegaInv->invert();
	  nomegaInv->negate();
	  nomegaInv->addConstant(-dualReg);
    }else{
	  //if gOuterSolve>=3, nomegaInv is the diag val without regularization for slack s
	  temp_diagS->copyFrom(*nomegaInv);
	  temp_diagS->addConstant(priReg);
	  temp_diagZ->setToConstant(-dualReg);	  
	}
  }
}



void NlpGenLinsys::factor(Data *  prob_in, Variables *vars_in)
{
  NlpGenVars * vars = (NlpGenVars *) vars_in;
  NlpGenData * prob = dynamic_cast<NlpGenData *>(prob_in);

  priReg=0.0;
  dualReg=0.0;

  assert( vars->validNonZeroPattern() );
  temp_diagS->setToZero();
  temp_diagY->setToZero();
  temp_diagZ->setToZero();

  if( nxlow + nxupp >= 0 ){
    dd->setToZero();
  }

  this->computeDiagonals( *dd, *temp_diagS,
			  *vars->t, *vars->lambda,
			  *vars->u, *vars->pi,
			  *vars->v, *vars->gamma,
			  *vars->w, *vars->phi );

  if( nxlow + nxupp >= 0 ){
  	if(separateHandDiag==0 ){
  	  prob->getDiagonalOfQ(*dq);
	  dd->axpy(1.0,*dq);
  	}
    temp_diagX->copyFrom(*dd);
  }

  nomegaInv->copyFrom(*temp_diagS);

  if(gOuterSolve<3){
    nomegaInv->invert();
    nomegaInv->negate();
  }

 
}

// sompute sigma = x^{-1}z
void NlpGenLinsys::computeDiagonals( OoqpVector& dd_, OoqpVector& omega,
				    OoqpVector& t,  OoqpVector& lambda,
				    OoqpVector& u,  OoqpVector& pi,
				    OoqpVector& v,  OoqpVector& gamma,
				    OoqpVector& w,  OoqpVector& phi )
{
  if( nxupp + nxlow >= 0 ) {
    if( nxlow > 0 ) 
	  dd_.axdzpy( 1.0, gamma, v, *ixlow );
    if( nxupp > 0 ) 
	  dd_.axdzpy( 1.0, phi  , w, *ixupp );
  }
  omega.setToZero();
  if ( mclow > 0 ) 
  	omega.axdzpy( 1.0, lambda, t, *iclow );
  if ( mcupp > 0 ) 
  	omega.axdzpy( 1.0, pi,     u, *icupp );
  //assert( omega->allPositive() );
}




double NlpGenLinsys::computeResidual(NlpGenData * data, OoqpVector & res_, OoqpVector & sol_,
			OoqpVector & solx_,OoqpVector & sols_,OoqpVector & soly_,OoqpVector & solz_)
{
  double result = 0.0;
  
  matXSYZMult(1.0,res_,-1.0, sol_, data,solx_,sols_,soly_,solz_);

  result = res_.infnorm();

  return result;
}





double NlpGenLinsys::computeResidual_FullKKT(NlpGenData * data, NlpGenResiduals *res_, NlpGenVars *sol_, 
											NlpGenVars *var_)
{
  double norm_res, norm_temp;

  sol_->y->negate();
  sol_->z->negate();

  res_->priWrk->setToConstant(priReg);
  res_->priWrk_S->setToConstant(priReg);
  res_->dualWrk->setToConstant(dualReg); res_->dualWrk->negate();
  res_->dualWrk_Z->setToConstant(dualReg); res_->dualWrk_Z->negate();

  // ******************************** 		rQ		**************************************
  // calculate rhs_rQ - (Q+\priReg)x - A'y - C'z;
  data->Qmult(1.0, *res_->rQ, -1.0, *sol_->x);
   
  data->ATransmult(1.0, *res_->rQ, -1.0, *sol_->y);
  if( mclow +mclow> 0 ) data->CTransmult(1.0, *res_->rQ, -1.0, *sol_->z);
  res_->rQ->axzpy(-1.0, *res_->priWrk, *sol_->x);

  // add dual variables for the var boundary constraints
  if( nxlow > 0 ) 
  	res_->rQ->axpy( 1.0, *sol_->gamma );
  if( nxupp > 0 ) 
  	res_->rQ->axpy( -1.0, *sol_->phi );

  norm_res = res_->rQ->infnorm();
  // ******************************** 		rA		**************************************
  data->Amult(1.0, *res_->rA, -1.0, *sol_->x);
  res_->rA->axzpy(-1.0, *res_->dualWrk, *sol_->y);
  
  norm_temp = res_->rA->infnorm();
  norm_res = MAX(norm_res,norm_temp);

  if( mclow +mclow> 0 ){
	// ******************************** 		rC		**************************************
	data->Cmult(1.0, *res_->rC, -1.0, *sol_->x);
	res_->rC->axzpy(-1.0, *res_->dualWrk_Z, *sol_->z);
	res_->rC->axpy(1.0, *sol_->s);

	norm_temp = res_->rC->infnorm();
	norm_res = MAX(norm_res,norm_temp);

	// ******************************** 		rz		**************************************
	res_->rz->axpy(1.0, *sol_->z);
	res_->rz->axzpy(-1.0, *res_->priWrk_S, *sol_->s);
	if( mclow > 0 ) res_->rz->axpy(1.0, *sol_->lambda);
	if( mcupp > 0 ) res_->rz->axpy(-1.0, *sol_->pi);

	norm_temp = res_->rz->infnorm();
	norm_res = MAX(norm_res,norm_temp);
  }
  if( mclow > 0 ) {
    // ******************************** 		rt		**************************************  	
    res_->rt->axpy(-1.0, *sol_->s);
    res_->rt->axpy(1.0, *sol_->t);
	res_->rt->selectNonZeros( *iclow );	

	norm_temp = res_->rt->infnorm();
	norm_res = MAX(norm_res,norm_temp);	
	// ******************************** 	  rlambda	  **************************************
	res_->rlambda->axzpy( -1.0, *var_->t, *sol_->lambda );
	res_->rlambda->axzpy( -1.0, *var_->lambda, *sol_->t );
	res_->rlambda->selectNonZeros( *iclow );	

	norm_temp = res_->rlambda->infnorm();
	norm_res = MAX(norm_res,norm_temp);	
  }
  
  if( mcupp > 0 ) {  	
	// ******************************** 	  ru	  **************************************
	res_->ru->axpy(1.0, *sol_->s);
    res_->ru->axpy(1.0, *sol_->u);
	res_->ru->selectNonZeros( *icupp );	

	norm_temp = res_->ru->infnorm();
	norm_res = MAX(norm_res,norm_temp);		
	// ******************************** 	  rpi	  **************************************
	res_->rpi->axzpy( -1.0, *var_->u, *sol_->pi );
	res_->rpi->axzpy( -1.0, *var_->pi, *sol_->u );
	res_->rpi->selectNonZeros( *icupp );	

	norm_temp = res_->rpi->infnorm();
	norm_res = MAX(norm_res,norm_temp);	
  }
  
  if( nxlow > 0 ) {
	// ********************************		rv		**************************************
    res_->rv->axpy(-1.0, *sol_->x);	
    res_->rv->axpy(1.0, *sol_->v);
	res_->rv->selectNonZeros( *ixlow );	
	
	norm_temp = res_->rv->infnorm();
	norm_res = MAX(norm_res,norm_temp);		
    // ******************************** 		rgamma		**************************************
    res_->rgamma->axzpy( -1.0, *var_->v, *sol_->gamma );
    res_->rgamma->axzpy( -1.0, *var_->gamma, *sol_->v );	
	res_->rgamma->selectNonZeros( *ixlow ); 

	norm_temp = res_->rgamma->infnorm();
	norm_res = MAX(norm_res,norm_temp);	
  }
  
  if( nxupp > 0 ) {
    // ******************************** 		rw		**************************************  	
  	res_->rw->axpy(1.0, *sol_->x);
    res_->rw->axpy(1.0, *sol_->w);
	res_->rw->selectNonZeros( *ixupp );
	
	norm_temp = res_->rw->infnorm();
	norm_res = MAX(norm_res,norm_temp);	
    // ******************************** 		rphi		**************************************  
    res_->rphi->axzpy( -1.0, *var_->w, *sol_->phi );
    res_->rphi->axzpy( -1.0, *var_->phi, *sol_->w );	
	res_->rphi->selectNonZeros( *ixupp );
	
	norm_temp = res_->rphi->infnorm();
	norm_res = MAX(norm_res,norm_temp);		
  }

  sol_->y->negate();
  sol_->z->negate();


  return norm_res;

}


void NlpGenLinsys::solve(Data * prob_in, Variables *vars_in,
			Residuals *res_in, Variables *step_in)
{
  NlpGenData      * prob  = (NlpGenData *) prob_in;
  NlpGenVars      * vars  = (NlpGenVars *) vars_in;
  NlpGenVars      * step  = (NlpGenVars *) step_in;
  NlpGenResiduals * resid   = (NlpGenResiduals *) res_in;

  if(gOuterSolve>=4 && !preCond->setGrad){
  	preCond->setGrad = true;
	preCond->curr_Iter = vars_in;
	preCond->curr_Step = step_in;
  }

  assert( vars->validNonZeroPattern() );
  assert( resid ->validNonZeroPattern() );

  step->x->copyFrom( *resid->rQ );
  
  if( nxlow > 0 && prob->nxlow >0) {
    OoqpVector & vInvGamma = *step->v;
    vInvGamma.copyFrom( *vars->gamma );
    vInvGamma.divideSome( *vars->v, *ixlow );
	
    step->x->axzpy ( 1.0, vInvGamma, *resid->rv );
    step->x->axdzpy( 1.0, *resid->rgamma, *vars->v, *ixlow );
  }
  if( nxupp > 0 ) {
    OoqpVector & wInvPhi   = *step->w;
    wInvPhi.copyFrom( *vars->phi );
    wInvPhi.divideSome( *vars->w, *ixupp );
	  
    step->x->axzpy (  1.0, wInvPhi,   *resid->rw );
    step->x->axdzpy( -1.0, *resid->rphi, *vars->w, *ixupp );
  }
  // start by partially computing step->s
  step->s->copyFrom( *resid->rz );
  if( mclow > 0 ) {
    OoqpVector & tInvLambda = *step->t;
	
    tInvLambda.copyFrom( *vars->lambda );
    tInvLambda.divideSome( *vars->t, *iclow );

    step->s->axzpy ( 1.0, tInvLambda, *resid->rt );
    step->s->axdzpy( 1.0, *resid->rlambda, *vars->t, *iclow );
  }
  if( mcupp > 0 ) {
    OoqpVector & uInvPi = *step->u;
	
    uInvPi.copyFrom( *vars->pi );
    uInvPi.divideSome( *vars->u, *icupp );
	
    step->s->axzpy (  1.0, uInvPi, *resid->ru );
    step->s->axdzpy( -1.0, *resid->rpi, *vars->u, *icupp );
  }
  step->y->copyFrom( *resid->rA );
  step->z->copyFrom( *resid->rC );

  {
    // Unfortunately, we need a temporary  OoqpVector for the solve,
    // Use step->lambda or step->pi
    OoqpVectorHandle ztemp;
    if( mclow > 0 ) {
      ztemp = step->lambda;
    } else {
      ztemp = step->pi;
    }	
    this->solveXYZS( *step->x, *step->y, *step->z, *step->s,
		     *ztemp, prob );
  }

  {
    if( mclow > 0 ) {
      step->t->copyFrom( *step->s );
      step->t->axpy( -1.0, *resid->rt );
      step->t->selectNonZeros( *iclow );

      step->lambda->copyFrom( *resid->rlambda );
      step->lambda->axzpy( -1.0, *vars->lambda, *step->t );
      step->lambda->divideSome( *vars->t, *iclow );
    }
    if( mcupp > 0 ) {
      step->u->copyFrom( *resid->ru );
      step->u->axpy( -1.0, *step->s );
      step->u->selectNonZeros( *icupp );

      step->pi->copyFrom( *resid->rpi );
      step->pi->axzpy( -1.0, *vars->pi, *step->u );
      step->pi->divideSome( *vars->u, *icupp );
    }
    if( nxlow > 0 ) {
      step->v->copyFrom( *step->x );
      step->v->axpy( -1.0, *resid->rv );
      step->v->selectNonZeros( *ixlow );
	
      step->gamma->copyFrom( *resid->rgamma );
      step->gamma->axzpy( -1.0, *vars->gamma, *step->v );
      step->gamma->divideSome( *vars->v, *ixlow );
    }
    if( nxupp > 0 ) {
      step->w->copyFrom( *resid->rw );
      step->w->axpy( -1.0, *step->x );
      step->w->selectNonZeros( *ixupp );
	
      step->phi->copyFrom( *resid->rphi );
      step->phi->axzpy( -1.0, *vars->phi, *step->w );
      step->phi->divideSome( *vars->w, *ixupp );
    }
  }  
  assert( step->validNonZeroPattern() );

}


void NlpGenLinsys::solve_IterRefine(Data * prob_in, Variables *vars_in,
			Residuals *res_in, Variables *step_in, Residuals *KKT_Resid_in ,Variables *KKT_sol_in)
{
  NlpGenData      * prob  	= (NlpGenData *) prob_in;
  NlpGenVars      * vars 	= (NlpGenVars *) vars_in;
  NlpGenVars      * step  	= (NlpGenVars *) step_in;
  NlpGenResiduals * resid   = (NlpGenResiduals *) res_in;
  NlpGenVars      		* kkt_sol  		= (NlpGenVars *) KKT_sol_in;
  NlpGenResiduals      	* kkt_resid  	= (NlpGenResiduals *) KKT_Resid_in;

  solve(prob_in, vars_in,res_in, step_in);

  if (gDoIR_Full==1){
	//set kkt_rhs --- now kkt_resid is the rhs of kkt system
	kkt_resid->copyFrom(resid);
	kkt_sol->copy(step);

	//compute residual of the kkt system, do IR if required;	
	prob->linsysRes_Full = computeResidual_FullKKT(prob, kkt_resid, kkt_sol, vars);
	
    if( prob->linsysRes_Full > gIRtol){
	  double currentRes=0, ir_iter=0; 
	  
  	  while(ir_iter<gMaxIR && prob->linsysRes_Full > gIRtol){

	    solve(prob_in, vars_in,kkt_resid, kkt_sol);
		
	    kkt_sol->saxpy(step,1.0);

		//set resid as kkt rhs
		kkt_resid->copyFrom(resid);
		//get resid 
	    currentRes = computeResidual_FullKKT(prob, kkt_resid, kkt_sol, vars);

	    if(currentRes < prob->linsysRes_Full){
		  prob->linsysRes_Full = currentRes;
		  step->copy(kkt_sol);	  
	    }else{
		  break;
	    }
	    ir_iter++;
	  }
   	}
  }
  
}


void NlpGenLinsys::solve_NTsteps(Data *prob_in, Variables *vars_in, Residuals *resids_in,
		   Variables *Nstep_in, Variables *Tstep_in, Variables *NTstep_in)
{
  NlpGenData      * prob  = (NlpGenData *) prob_in;
  NlpGenVars      * vars  = (NlpGenVars *) vars_in;
  NlpGenVars      * Nstep  = (NlpGenVars *) Nstep_in;
  NlpGenVars      * Tstep  = (NlpGenVars *) Tstep_in;
  NlpGenResiduals * resid   = (NlpGenResiduals *) resids_in;
  NlpGenVars	  * NTstep  = (NlpGenVars *) NTstep_in;  

  assert( vars->validNonZeroPattern() );
  assert( resid ->validNonZeroPattern() );

  // turn off IR for aug sys
  int DoIR_Aug_temp = gDoIR_Aug;
  gDoIR_Aug = 0;

  //******************************    solve N step    ******************************  
  // build rhs for N step
  Nstep->x->setToZero();
  Nstep->s->setToZero();
  Nstep->y->copyFrom( *resid->rA );
  Nstep->z->copyFrom( *resid->rC );
  
  // Unfortunately, we need a temporary  OoqpVector for the solve,
  // Use step->lambda or step->pi
  OoqpVectorHandle ztemp_N;
  if( mclow > 0 ) {
    ztemp_N = Nstep->lambda;
  } else {
    ztemp_N = Nstep->pi;
  }	
  this->solveXYZS( *Nstep->x, *Nstep->y, *Nstep->z, *Nstep->s,
		     *ztemp_N, prob );


  //******************************    solve T step    ******************************  
  
  // build rhs for T step 
  // dq = DiagonalOfQ; dd = complement part + regularizaion + (DiagonalOfQ)   --> (DiagonalOfQ) if separateHandDiag!=1
  if(separateHandDiag != 1 ){
	resid->priWrk->copyFrom(*dd);
	resid->priWrk->axpy(-1,*dq);
  }else{
 	resid->priWrk->copyFrom(*dd);
  }
  resid->priWrk->componentMult(*Nstep->x);
	  
  // compute extra part for t step	--- for x
  prob->Qmult(1.0, *resid->priWrk, 1.0, *Nstep->x ); 
  resid->priWrk->negate();
  
  // compute extra part for t step	--- for s
  resid->priWrk_S->copyFrom(*temp_diagS);
  resid->priWrk_S->componentMult(*Nstep->s);  
  resid->priWrk_S->negate();
	    
  Tstep->x->copyFrom( *resid->rQ );
  if( nxlow > 0 && prob->nxlow >0) {
	OoqpVector & vInvGamma = *Tstep->v;
	vInvGamma.copyFrom( *vars->gamma );
	vInvGamma.divideSome( *vars->v, *ixlow );
	 
	Tstep->x->axzpy ( 1.0, vInvGamma, *resid->rv );
	Tstep->x->axdzpy( 1.0, *resid->rgamma, *vars->v, *ixlow );
  }
  if( nxupp > 0 ) {
	OoqpVector & wInvPhi	= *Tstep->w;
	wInvPhi.copyFrom( *vars->phi );
	wInvPhi.divideSome( *vars->w, *ixupp );
	   
	Tstep->x->axzpy (	1.0, wInvPhi,	*resid->rw );
	Tstep->x->axdzpy( -1.0, *resid->rphi, *vars->w, *ixupp );
  }
   
  Tstep->s->copyFrom( *resid->rz );
  if( mclow > 0 ) {
	OoqpVector & tInvLambda = *Tstep->t; 
	tInvLambda.copyFrom( *vars->lambda );
	tInvLambda.divideSome( *vars->t, *iclow );
  
	Tstep->s->axzpy ( 1.0, tInvLambda, *resid->rt );
	Tstep->s->axdzpy( 1.0, *resid->rlambda, *vars->t, *iclow );
  }
  if( mcupp > 0 ) {
	OoqpVector & uInvPi = *Tstep->u; 
	uInvPi.copyFrom( *vars->pi );
	uInvPi.divideSome( *vars->u, *icupp );
	 
	Tstep->s->axzpy (	1.0, uInvPi, *resid->ru );
	Tstep->s->axdzpy( -1.0, *resid->rpi, *vars->u, *icupp );
  }

  // get primal part of the rhs for  x an s
  Tstep->x->axpy(1,*resid->priWrk); 	
  Tstep->s->axpy(1,*resid->priWrk_S); 

  // get primal part of the rhs for y and z
  Tstep->y->setToZero();
  Tstep->z->setToZero();  	
  
  // Unfortunately, we need a temporary  OoqpVector for the solve,
  // Use step->lambda or step->pi
  OoqpVectorHandle ztemp_T;
  if( mclow > 0 ) {
    ztemp_T = Tstep->lambda;
  } else {
    ztemp_T = Tstep->pi;
  }	
  this->solveXYZS( *Tstep->x, *Tstep->y, *Tstep->z, *Tstep->s,
		     *ztemp_T, prob );


  //******************************    set d_step = n_step+t_step    ******************************  
  if(gdWd_test >= 2 ) 
  	NTstep_in->mergeNTstep(Nstep_in,Tstep_in,vars_in,resids_in); 

  gDoIR_Aug = DoIR_Aug_temp;

}


void NlpGenLinsys::solveXYZS( OoqpVector& stepx, OoqpVector& stepy,
			       OoqpVector& stepz, OoqpVector& steps,
			       OoqpVector& /* ztemp */,
			       NlpGenData *  prob  )
{
  if(gOuterSolve<3)
  	stepz.axzpy( -1.0, *nomegaInv, steps );


  if(gOuterSolve==1) {
    ///////////////////////////////////////////////////////////////
    // Iterative refinement and Schur complement based decomposition, compress s
    ///////////////////////////////////////////////////////////////
    solveCompressedIterRefin(stepx,stepy,stepz,prob);

  } else if(gOuterSolve==0) {
    ///////////////////////////////////////////////////////////////
    // Schur complement based decomposition, compress s
    ///////////////////////////////////////////////////////////////
    this->joinRHS( *rhs, stepx, stepy, stepz );
    this->solveCompressed( *rhs );
    this->separateVars( stepx, stepy, stepz, *rhs );

  } else if(gOuterSolve==3) {
    ///////////////////////////////////////////////////////////////
    // Default solve - Schur complement based decomposition or petsc, do not compress s
    ///////////////////////////////////////////////////////////////
    if(gUsePetscOuter==1&&gUsePetsc==1) 
	  solveCompressedAugXSYZ_PETSC(stepx, steps, stepy, stepz, prob);
	else
	  solveCompressedAugXSYZ(stepx, steps, stepy, stepz, prob);
  } else if(gOuterSolve==2){
    ///////////////////////////////////////////////////////////////
    // BiCGStab, compress s (have never been tested)
    ///////////////////////////////////////////////////////////////
    solveCompressedBiCGStab(stepx,stepy,stepz,prob);
  } else if(gOuterSolve==4 || gOuterSolve==5){
    ///////////////////////////////////////////////////////////////
    // BiCGStab with the aggregation as preconditioner
    ///////////////////////////////////////////////////////////////  
    solveBiCGStab(stepx,steps,stepy,stepz,prob);    
  }

  stepy.negate();
  stepz.negate();

  if(gOuterSolve<3){	
    steps.axpy( -1.0, stepz );
    steps.componentMult( *nomegaInv );
    steps.negate();
  }
}


void NlpGenLinsys::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			     OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  // joinRHS has to be delegated to the factory. This is true because
  // the rhs may be distributed across processors, so the factory is the
  // only object that knows with certainly how to scatter the elements.
  factory->joinRHS( rhs_in, rhs1_in, rhs2_in, rhs3_in );
}

void NlpGenLinsys::joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			     OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in )
{
  factory->joinRHSXSYZ( rhs_in, rhs1_in, rhs2_in, rhs3_in, rhs4_in );
}

void NlpGenLinsys::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				  OoqpVector& z_in, OoqpVector& vars_in )
{
  // separateVars has to be delegated to the factory. This is true because
  // the rhs may be distributed across processors, so the factory is the
  // only object that knows with certainly how to scatter the elements.
  factory->separateVars( x_in, y_in, z_in, vars_in );
}

void NlpGenLinsys::separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in,
				  OoqpVector& y_in,OoqpVector& z_in, OoqpVector& vars_in )
{
  factory->separateVarsXSYZ( x_in, s_in, y_in, z_in, vars_in );
}


void NlpGenLinsys::solveCompressedAugXSYZ(OoqpVector& stepx, OoqpVector& steps,
					   OoqpVector& stepy,
					   OoqpVector& stepz,
					   NlpGenData* prob)
{
    this->joinRHSXSYZ( *rhs, stepx, steps, stepy, stepz );
	sol->copyFrom(*rhs);
    this->solveCompressed( *sol );
	
	//use res as residual
	res->copyFrom(*rhs);	
    //separate var and compute residual,  do IR if required;	
	prob->linsysRes = computeResidual(prob, *res, *sol, stepx, steps, stepy, stepz);
	
	if(gDoIR_Aug==1 && prob->linsysRes > gIRtol){
	  double currentRes=0, ir_iter=0;	
	  int callBackBestSol = 0;
	  
	  while(ir_iter<gMaxIR && prob->linsysRes > gIRtol){
	    sol2->copyFrom(*res);
	    this->solveCompressed( *sol2 );
	  
	    sol2->axpy(1.0,*sol);
	    res->copyFrom(*rhs);	
	    currentRes = computeResidual(prob, *res, *sol2, stepx, steps, stepy, stepz);  

	    if(currentRes < prob->linsysRes){
		  prob->linsysRes = currentRes;
		  sol->copyFrom(*sol2);		
	    }else{
	      callBackBestSol=1;
	      break;
		}
		ir_iter++;
	  }

	  if(callBackBestSol == 1)
	    this->separateVarsXSYZ( stepx, steps, stepy, stepz, *sol);
	}

	prob->KryIter = KryIter;
}


void NlpGenLinsys::solveCompressedAugXSYZ_PETSC(OoqpVector& stepx, OoqpVector& steps,
					   OoqpVector& stepy,
					   OoqpVector& stepz,
					   NlpGenData* prob)
{

#ifdef WITH_PETSC
  int ierr;

  Vec x;  
  Mat LinSysMat_PETSC;
  Mat PCMat_PETSC;
  PC  precond_Method;

  Mat_Shell_Petsc ctx_la;
  ctx_la.NlpGenLinsysPtr = this;
  ctx_la.NlpGenDataPtr   = prob;
  ctx_la.solx_= &stepx;
  ctx_la.sols_= &steps;
  ctx_la.soly_= &stepy;
  ctx_la.solz_= &stepz;

  int nb_col = nx+mz+my+mz;
  
  MatCreateShell(PETSC_COMM_SELF, nb_col, nb_col, PETSC_DETERMINE, PETSC_DETERMINE, (void*)&ctx_la, &LinSysMat_PETSC);
  MatShellSetOperation(LinSysMat_PETSC, MATOP_MULT, (void (*)(void))_user_MatMultXSYZ);


  // ********  Create linear solver context ********
  ierr = KSPCreate(PETSC_COMM_SELF,&mKsp);assert(ierr == 0);
	
  // ********  Set operators. Here the matrix that defines the linear system also serves as the preconditioning matrix. ********
  ierr = KSPSetOperators(mKsp, LinSysMat_PETSC, LinSysMat_PETSC); assert(ierr == 0);

  ierr = KSPGetPC(mKsp,&precond_Method);assert(ierr == 0);

  assert(gUser_Defined_PC!=0);
  if (gUser_Defined_PC != 0) 
  {
    ierr = PCSetType(precond_Method,PCSHELL);assert(ierr == 0);
    Shell_PC_Petsc *shell_la;
    PC_Shell_Create(&shell_la);
    PCShellSetApply(precond_Method,PC_Shell_Apply_Sparse);
    PCShellSetContext(precond_Method,shell_la);
    PCShellSetDestroy(precond_Method,PC_Shell_Destroy);
    PC_Shell_SetUp(precond_Method,this,prob);
  }

  ierr = KSPSetFromOptions(mKsp);assert(ierr == 0);
  ierr = KSPSetTolerances(mKsp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); assert(ierr == 0);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //   Solve the linear system Ax=b by petsc
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  this->joinRHSXSYZ( *rhs, stepx, steps, stepy, stepz );

  double *x_array = new double[nb_col];
  double *b_array = new double[nb_col];
  this->copyXSYZ_toArray( *rhs, b_array, nb_col);

  Vec b_la,x_la;
  VecCreateSeqWithArray(PETSC_COMM_SELF,1,nb_col,b_array,&b_la);
  VecCreateSeqWithArray(PETSC_COMM_SELF,1,nb_col,x_array,&x_la);

  ierr = KSPSolve(mKsp,b_la,x_la); assert( ierr  == 0);

  this->copyXSYZ_fromArray( *sol, x_array, nb_col);  


  VecDestroy(&b_la);
  VecDestroy(&x_la);  


  //use res as residual
  res->copyFrom(*rhs);
  prob->linsysRes = computeResidual(prob, *res, *sol, stepx, steps, stepy, stepz);


  ierr = KSPGetIterationNumber(mKsp, &KryIter); assert( ierr  == 0);
  prob->KryIter = KryIter;
//  std::cout<<"using "<<KryIter<<" GMRES iterations"<<std::endl;
  
  ierr = KSPDestroy( &mKsp );assert(ierr == 0);

  delete [] x_array;
  delete [] b_array;

#else
  assert("require PETSC" &&0);
#endif

}


void NlpGenLinsys::solveCompressedIterRefin(OoqpVector& stepx,
					   OoqpVector& stepy,
					   OoqpVector& stepz,
					   NlpGenData* prob)
{
  printf("Not Finished! STOP! \n" );
  exit(1);
}

void NlpGenLinsys::solveCompressedBiCGStab(OoqpVector& stepx,
					   OoqpVector& stepy,
					   OoqpVector& stepz,
					   NlpGenData* data)
{
#ifdef TIMING
	  vector<double> histRelResid;
	  double tTot=MPI_Wtime(), tSlv=0., tResid=0., tTmp;
#endif
	
	  this->joinRHS( *rhs, stepx, stepy, stepz );
	
	  //aliases
	  OoqpVector &r0=*res2, &dx=*sol2, &v=*res3, &t=*res4, &p=*res5;
	  OoqpVector &x=*sol, &r=*res, &b=*rhs;
	
	  const double tol=1e-10, EPS=2e-16; 
	  const int maxit=500;
	
	  double n2b=b.twonorm(), tolb=n2b*tol;
	  int flag; double iter=0.;
	  double rho=1., omega=1., alpha;  
	
	
#ifdef TIMING
	  gOuterBiCGIter=0;
	  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	  tTmp=MPI_Wtime();
#endif
	  //starting guess/point
	  x.copyFrom(b); solveCompressed(x); //solution to the approx. system
#ifdef TIMING
	  tSlv += (MPI_Wtime()-tTmp);
	  tTmp=MPI_Wtime();
	  gOuterBiCGIter++;
#endif
	  //initial residual: res=res-A*x
	  r.copyFrom(b); 
	  matXYZMult(1.0,r, -1.0,x, data, stepx,stepy,stepz);
#ifdef TIMING
	  tResid += (MPI_Wtime()-tTmp);
#endif  
	  double normr=r.twonorm(), normr_min=normr, normr_act=normr;
#ifdef TIMING
	  histRelResid.push_back(normr/n2b);
#endif
	  
	  //quick return if solve is accurate enough
	  if(normr<tolb) { 
		 this->separateVars( stepx, stepy, stepz, x );
#ifdef TIMING
		 tTot = MPI_Wtime() - tTot;
		 if(0==myRank) {
		   cout << "Outer BiCGStab 0 iterations. Rel.res.nrm:" << normr/n2b  << endl;
		   cout << "solveXYZS w/ BiCGStab times: solve " << tSlv
			<< "  matvec " << tResid
			<< "  total " << tTot << endl; 
		 }
#endif
		 return;
	  }
	
	  //arbitrary vector
	  r0.copyFrom(b);
	
	  //main loop
	  int it=0; while(it<maxit) {
		flag=-1; //reset flag
		double rho1=rho; double beta;
		rho=r0.dotProductWith(r);
		if(0.0==rho) {flag=4; break;}
	
		//first half of the iterate
		{
		  if(it==0) p.copyFrom(r);
		  else {
		beta = (rho/rho1)*(alpha/omega);
		if(beta==0.0) { flag=4; break; }
		
		//-------- p = r + beta*(p - omega*v) --------
		p.axpy(-omega, v); p.scale(beta); p.axpy(1.0, r);
		  }
#ifdef TIMING
		  tTmp = MPI_Wtime();
#endif
		  //precond: ph = \tilde{K}^{-1} p
		  dx.copyFrom(p); solveCompressed(dx); 
#ifdef TIMING
		  tSlv += (MPI_Wtime()-tTmp);
		  tTmp = MPI_Wtime();
		  gOuterBiCGIter++;
#endif
		  //mat-vec: v = K*ph
		  matXYZMult(0.0,v, 1.0,dx, data, stepx,stepy,stepz);
#ifdef TIMING
		  tResid += (MPI_Wtime()-tTmp);
#endif
	
		  double rtv=r0.dotProductWith(v);
		  if(rtv==0.0) { flag=4; break; }
		  
		  alpha=rho/rtv; 
		  // x = x + alpha*dx (x=x+alpha*ph)
		  x.axpy( alpha, dx);
		  // r = r-alpha*v (s=r-alpha*v)
		  r.axpy(-alpha, v);
		  
		  //check for convergence
		  normr=r.twonorm();
	
#ifdef TIMING
		  histRelResid.push_back(normr/n2b);
#endif
		  if(normr<=tolb) {
#ifdef TIMING
		  tTmp=MPI_Wtime();
#endif
		//compute the actual residual
		OoqpVector& res=dx; //use dx 
		res.copyFrom(b); matXYZMult(1.0,res, -1.0,x, data, stepx,stepy,stepz);
#ifdef TIMING
		tResid += (MPI_Wtime()-tTmp);
#endif
	
		normr_act = res.twonorm();
		if(normr_act<=tolb) {
		  //converged
#ifdef TIMING
		  histRelResid[histRelResid.size()-1]=normr_act/n2b;
#endif
		  flag=0; iter=it+0.5; break;
		}
		  } //~end of convergence test
		}
	
		
		//second half of the iterate now
		{
#ifdef TIMING
		  tTmp=MPI_Wtime();
#endif
		  //preconditioner
		  dx.copyFrom(r); solveCompressed(dx);
#ifdef TIMING
		  tSlv += (MPI_Wtime()-tTmp);
		  tTmp=MPI_Wtime();
		  gOuterBiCGIter++;
#endif
		  //mat-vec
		  matXYZMult(0.0,t, 1.0,dx, data, stepx,stepy,stepz);
#ifdef TIMING
		  tResid += (MPI_Wtime()-tTmp);
#endif
		  double tt = t.dotProductWith(t);
		  if(tt==0.0) { flag=4; break;}
		  
		  omega=t.dotProductWith(r); omega /= tt;
		  
		  // x=x+omega*dx  (x=x+omega*sh)
		  x.axpy( omega, dx);
		  // r = r-omega*t (r=s-omega*sh)
		  r.axpy(-omega, t);
		  //check for convergence
		  normr=r.twonorm();
#ifdef TIMING
		  histRelResid.push_back(normr/n2b);
#endif
	
		  if(normr<=tolb) {
#ifdef TIMING
		  tTmp=MPI_Wtime();
#endif
		//compute the actual residual
		OoqpVector& res=dx; //use dx 
		res.copyFrom(b); matXYZMult(1.0,res, -1.0,x, data, stepx,stepy,stepz);
#ifdef TIMING
		  tResid += (MPI_Wtime()-tTmp);
#endif
	
		normr_act = res.twonorm();
		//cout << "Outer BiCG - actual rel.res.nrm: " << normr_act/n2b << endl;
		if(normr_act<=tolb) {
		  //converged
#ifdef TIMING 
		  histRelResid[histRelResid.size()-1]=normr_act/n2b;
#endif
		  flag=0; iter=it+1.; break;
		} // else continue - To Do: detect stagnation (flag==3)
		  } else {
		//To Do: detect stagnation/divergence and rollback to min.norm. iterate
		//for now we print a warning and exit in case residual increases.
		if(normr>normr_min) {
#ifdef TIMING
		  if(0==myRank) 
			cout << "Outer BiCG - Increase in BiCGStab residual. Old=" << normr_min
			 << "  New=" << normr << endl;
#endif
		  flag=5; iter=it+1.; break;
		} else normr_min=normr;
		
		  } //~end of convergence test
		} //~end of scoping
		
		it++;
	  } //~ end of BiCGStab loop
	  
	  //warning/error messaging
	  if(flag!=0) {
#ifdef TIMING
		if(0==myRank) 
		  cout << "Outer BiCG - convergence issues: flag=" << flag << ". "
		   << iter << " iterations" 
		   << " rel.res.norm=" << normr_act/n2b << endl;
#endif
	  }
	  this->separateVars( stepx, stepy, stepz, x );
	
#ifdef TIMING
	  tTot = MPI_Wtime()-tTot;
	  if(0==myRank) {
		cout << "Outer BiCGStab " << iter << " iterations. Rel.res.nrm:";
		for(size_t it=0; it<histRelResid.size(); it++)
		  cout << histRelResid[it] << " | ";	
		cout << endl;
		cout << "solveXYZS w/ BiCGStab times: solve=" << tSlv
		 << "  matvec=" << tResid
		 << "  total=" << tTot << endl; 
	  }
#endif  

}



void NlpGenLinsys::solveBiCGStab(OoqpVector& stepx,
					  OoqpVector& steps,
					  OoqpVector& stepy,
					  OoqpVector& stepz,
					  NlpGenData* data)
{
  this->joinRHSXSYZ( *rhs, stepx, steps, stepy, stepz );

  //aliases
  OoqpVector &r0=*res2Bicg, &dx=*sol2Bicg, &v=*res3Bicg, &t=*res4Bicg, &p=*res5Bicg;
  OoqpVector &x=*sol, &r=*res, &b=*rhs;

  const double tol=1e-10, EPS=2e-16; 
  const int maxit=20;

  double n2b=b.twonorm(), tolb=n2b*tol*1e10;
  int flag; double iter=0.;
  double rho=1., omega=1., alpha=1.;  
  double normr,normr_act, normr_min=1e2; 

  gOuterBiCGIter=0;

  //starting guess/point
  x.setToZero(); 

  //arbitrary vector
  r0.copyFrom(b);

  r.copyFrom(b);
  p.setToZero();
  v.setToZero();

//  bool flagTemp=false;
  //main loop
  int it=1; 
  while(it<maxit) {
    flag=-1; //reset flag
    double rho1=rho; double beta;
	// step 1
    rho=r0.dotProductWith(r);
    if(0.0==rho) {flag=4; break;}

    //first half of the iterate
    {
	  // step 2
	  beta = (rho/rho1)*(alpha/omega);
	  if(beta==0.0) { flag=4; break; }
	
	  // step 3 -------- p = r + beta*(p - omega*v) --------
	  p.axpy(-omega, v); p.scale(beta); p.axpy(1.0, r);
      
      //step 4: precond: ph = \tilde{K}^{-1} p
      dx.copyFrom(p); 
	  solveCompressed(dx); 
	  gOuterBiCGIter++;

      //step 5: mat-vec: v = Mat*ph
      matXSYZMult(0.0,v, 1.0,dx, data, stepx,steps,stepy,stepz);

      //step 6: alpha=rho/(r0.dot.v)
      double rtv=r0.dotProductWith(v);
      if(rtv==0.0) { flag=4; break; }
      alpha = rho/rtv; 

      //step 11.5: add \alpha ph term: x = x + alpha*dx (x=x+alpha*ph)
      x.axpy( alpha, dx);
	  if(gOuterSolve>=5){
		this->separateVarsXSYZ( stepx, steps, stepy, stepz, x );
	  }

      //step 7: r = r-alpha*v (s=r-alpha*v)    
      r.axpy(-alpha, v);
	     
      //check for convergence
      normr=r.twonorm();
	  if(it==1) normr_min = normr;
      if(normr<=tolb ) {
	    //compute the actual residual
	    OoqpVector& res=dx; //use dx 
	    res.copyFrom(b); 
	    matXSYZMult(1.0,res, -1.0,x, data, stepx,steps,stepy,stepz);

	    normr_act = res.twonorm();
	    if(normr_act<=tolb) {
	      //converged
		  flag=0; iter=it+0.5; 
		  break;
	    }
      } //~end of convergence test
    } //~end of first half update

    
    //second half of the iterate now
    {
      //step 8: preconditioner z = \tilde{K}^{-1} s
      dx.copyFrom(r); 
	  solveCompressed(dx);
	  gOuterBiCGIter++;

      //step 9: mat-vec: t = Mat*z 
      matXSYZMult(0.0,t, 1.0,dx, data, stepx,steps,stepy,stepz);

      //step 10: t's/t't  
      double tt = t.dotProductWith(t);
      if(tt==0.0) { flag=4; break;}
      omega = t.dotProductWith(r); 
	  omega /= tt;
      
      // step 11: x=x+omega*dx  (x=x+omega*sh)
      x.axpy( omega, dx);
	  if(gOuterSolve>=5){
		this->separateVarsXSYZ( stepx, steps,stepy, stepz, x );
	  }	  
      // step 13: r = r-omega*t (r=s-omega*sh)
      r.axpy(-omega, t);
	  
      //check for convergence
      normr=r.twonorm();
      if(normr<=tolb) {
		//compute the actual residual
		OoqpVector& res=dx; //use dx 
		res.copyFrom(b); 
		matXSYZMult(1.0,res, -1.0,x, data, stepx,steps,stepy,stepz);

		normr_act = res.twonorm();
		//cout << "Outer BiCG - actual rel.res.nrm: " << normr_act/n2b << endl;
		if(normr_act<=tolb) {
		  //converged
		  flag=0; iter=it+1.; 
		  break;
		} // else continue - To Do: detect stagnation (flag==3)
      } else {
		//To Do: detect stagnation/divergence and rollback to min.norm. iterate
		//for now we print a warning and exit in case residual increases.
//		if(normr>normr_min) {
//		  flag=5; iter=it+1.; break;
//		} else 
//		  normr_min=normr;
      } //~end of convergence test
    } //~end of scoping
    it++;
  } //~ end of BiCGStab loop
  
  //warning/error messaging
  if(flag!=0) {
	  cout << "Outer BiCG - convergence issues: flag=" << flag << ". "
		<< iter << " iterations" 
		<< " rel.res.norm=" << normr_act/n2b << endl;
  }
  this->separateVarsXSYZ( stepx, steps, stepy, stepz, x );

}


/**
 * res = beta*res - alpha*mat*sol
 * stepx, stepy, stepz are used as temporary buffers
 */
void NlpGenLinsys::matXYZMult(double beta,  OoqpVector& res, 
			     double alpha, OoqpVector& sol, 
			     NlpGenData* data,
			     OoqpVector& solx, 
			     OoqpVector& soly, 
			     OoqpVector& solz)
{
  this->separateVars( solx, soly, solz, sol );
  this->separateVars( *resx, *resy, *resz, res);

  data->Qmult(beta, *resx, alpha, solx);
  resx->axzpy(alpha, *dd, solx);
  data->ATransmult(1.0, *resx, alpha, soly);
  data->CTransmult(1.0, *resx, alpha, solz);

  data->Amult(beta, *resy, alpha, solx);
  //cout << "resy norm: " << resy->twonorm() << endl;
  data->Cmult(beta, *resz, alpha, solx);
  resz->axzpy(alpha, *nomegaInv, solz);
  //cout << "resz norm: " << resz->twonorm() << endl;
  this->joinRHS( res, *resx, *resy, *resz );
}



/**
 * res = beta*res + alpha*mat*sol
 * stepx, stepy, stepz are used as temporary buffers
 */
void NlpGenLinsys::matXSYZMult( double beta,  OoqpVector& res_, 
			 double alpha, OoqpVector& sol_, 
			 NlpGenData* data,
			 OoqpVector& solx_, 
			 OoqpVector& sols_,
			 OoqpVector& soly_, 
			 OoqpVector& solz_)
{

  this->separateVarsXSYZ( solx_, sols_, soly_, solz_, sol_ );
  this->separateVarsXSYZ( *resx, *ress,*resy, *resz, res_);

  // if separateHandDiag==1, 	dd = complement part + regularizaion
  // 			 otherwise,	dd = complement part + regularizaion + DiagonalOfQ, 	
  // dq= DiagonalOfQ
  // after this if function,  dd=complement part + regularizaion
  if(separateHandDiag != 1 ){
    dd->axpy(-1,*dq);
  }
  
  data->Qmult(beta, *resx, alpha, solx_);
  resx->axzpy(alpha, *dd, solx_);
  
  data->ATransmult(1.0, *resx, alpha, soly_);
  data->CTransmult(1.0, *resx, alpha, solz_);
  

  ress->scale(beta);
  ress->axpy(-alpha,solz_);
  ress->axzpy(alpha,*temp_diagS,sols_);

  data->Amult(beta, *resy, alpha, solx_);
  resy->axzpy(alpha, *temp_diagY,soly_);

  data->Cmult(beta, *resz, alpha, solx_);
  resz->axpy(-alpha, sols_);
  resz->axzpy(alpha, *temp_diagZ, solz_);

  //cout << "resz norm: " << resz->twonorm() << endl;
  this->joinRHSXSYZ( res_, *resx, *ress, *resy, *resz );

  //FIXME_NY: recover dd, we may need it later? as in eval_xWx. 
  if(separateHandDiag != 1 ){
    dd->axpy(1,*dq);
  }  
}


void NlpGenLinsys::copyXSYZ_fromArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col)
{
  factory->copyXSYZ_fromArray( vec_xsyz, array_in, nb_col );
}
void NlpGenLinsys::copyXSYZ_toArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col)
{
  factory->copyXSYZ_toArray( vec_xsyz, array_in, nb_col );
}






/**
* check curvature xWx
 */
double NlpGenLinsys::eval_xWx( NlpGenData *prob, NlpGenResiduals *resid, NlpGenVars *steps)
{
  double dtWd = 0.;

  //use resid->priWrk as temp vector
  prob->Qmult( 0.0, *resid->priWrk, 1.0, *steps->x ); 
  dtWd += resid->priWrk->dotProductWith(*steps->x);


  // dd - DiagonalOfQ = complement part + regularizaion, otherwise the diagnol term has been computed twice
  // if separateHandDiag==1, dd=complementary part + regularizaion
  // otherwise dd=complement part + regularizaion + DiagonalOfQ
  // dq = DiagonalOfQ;
  if(separateHandDiag != 1 ){  
 	resid->priWrk->copyFrom(*dd);
  	resid->priWrk->axpy(-1,*dq);
  }else{
    resid->priWrk->copyFrom(*dd);
  }
  resid->priWrk->componentMult(*steps->x);
  dtWd += resid->priWrk->dotProductWith(*steps->x);


  if( gOuterSolve>=3 ){ 
  	resid->priWrk_S->setToZero();
	resid->priWrk_S->axzpy(1.,*temp_diagS,*steps->s);
    dtWd += resid->priWrk_S->dotProductWith(*steps->s);
  }


  return dtWd;  
}


/**
* compute Quantities For the usage of Dual Reg
 */
void NlpGenLinsys::computeQuantitiesForDualReg( 	NlpGenData *prob, NlpGenVars *vars,
															NlpGenResiduals *resid, NlpGenVars *steps,
														  	double *dualRegQuantities)
{
  switch ( gUseDualRegAlg ) {
    case 1 : {
	  assert(dualRegQuantities);	
      double nrmp;
	  double nrmc;
	  double dAAd = 0.0;

	  // nrmp = norm(p);
	  // nrmc = norm(c);	
	  // dAAd = d'*(A'*A)*d;

	  nrmp = (steps->y)->Norm(PIPS_NORM_TWO, steps->z);
	  nrmc = resid->getKKTRhsNorm_Primal(prob,vars,PIPS_NORM_TWO);
	  
	  //use resid->priWrk as temp vector
	  prob->Amult( 0.0, *resid->dualWrk, 1.0, *steps->x ); 

	  resid->dualWrk_Z->copyFrom(*steps->s);
	  prob->Cmult( -1.0, *resid->dualWrk_Z, 1.0, *steps->x ); 
	  	  
	  dAAd += resid->dualWrk->dotProductWith(*resid->dualWrk);
	  dAAd += resid->dualWrk_Z->dotProductWith(*resid->dualWrk_Z);

	  dualRegQuantities[0] = nrmp;
	  dualRegQuantities[1] = nrmc;
	  dualRegQuantities[2] = dAAd;
	 
    }; break;
    case 10 : {
	  assert(dualRegQuantities);	
      double nrmp;
	  double nrmc;
	  double dAAd = 0.0;

	  // nrmp = norm(p);
	  // nrmc = norm(c);	
	  // dAAd = d'*(A'*A)*d;

	  nrmp = (steps->y)->Norm(PIPS_NORM_TWO, steps->z);
	  nrmc = resid->getKKTRhsNorm_Primal(prob,vars,PIPS_NORM_TWO);
	  
	  //use resid->priWrk as temp vector
	  prob->Amult( 0.0, *resid->dualWrk, 1.0, *steps->x ); 

	  resid->dualWrk_Z->copyFrom(*steps->s);
	  prob->Cmult( -1.0, *resid->dualWrk_Z, 1.0, *steps->x ); 
	  	  
	  dAAd += resid->dualWrk->dotProductWith(*resid->dualWrk);
	  dAAd += resid->dualWrk_Z->dotProductWith(*resid->dualWrk_Z);

	  dualRegQuantities[0] = nrmp;
	  dualRegQuantities[1] = nrmc;
	  dualRegQuantities[2] = dAAd;
	 
    }; break;	
    default :


		
      break;
  }

  

}





#ifdef WITH_PETSC
// tell petsc how to compute Ax+y,  here we use res and sol2 as temp vec:   y=res2, x=sol2
PetscErrorCode _user_MatMultXSYZ(Mat MatInPetsc, Vec x,  Vec y)
{
  Mat_Shell_Petsc *shell;
  MatShellGetContext(MatInPetsc, &shell);

  NlpGenLinsys  *LinSysMat  = shell->NlpGenLinsysPtr;
  NlpGenData	*data		= shell->NlpGenDataPtr;
  OoqpVector	*solx_ 	  	= shell->solx_;
  OoqpVector	*sols_ 	  	= shell->sols_;
  OoqpVector	*soly_ 	  	= shell->soly_;
  OoqpVector	*solz_ 	  	= shell->solz_;

  
  int nb_row = LinSysMat->nx+LinSysMat->my+LinSysMat->mz+LinSysMat->mz;
  int nb_col = nb_row;

  double  *x_array, *y_array;
  VecGetArray(x,&x_array);
  VecGetArray(y,&y_array);

//  LinSysMat->copyXSYZ_fromArray( *LinSysMat->res, y_array, nb_col);
  LinSysMat->copyXSYZ_fromArray( *LinSysMat->sol2, x_array, nb_col);

  LinSysMat->matXSYZMult(0.0,*LinSysMat->res2,1.0,*LinSysMat->sol2,data,*solx_,*sols_,*soly_,*solz_);

  LinSysMat->copyXSYZ_toArray( *LinSysMat->res2, y_array, nb_col);
  
  VecRestoreArray(x,&x_array);
  VecRestoreArray(y,&y_array);
  	
  return 0;
}


PetscErrorCode PC_Shell_Create(Shell_PC_Petsc **shell)
{
  (*shell)= new Shell_PC_Petsc;
  return 0;
}

PetscErrorCode PC_Shell_SetUp(PC pc,   NlpGenLinsys	*NlpGenLinsysPtr_in, NlpGenData	*NlpGenDataPtr_in)
{
  Shell_PC_Petsc *shell;
  PCShellGetContext(pc,(void**)&shell);
  shell->NlpGenLinsysPtr = NlpGenLinsysPtr_in;
  shell->NlpGenDataPtr   = NlpGenDataPtr_in;
  return 0;
}

PetscErrorCode PC_Shell_Destroy(PC pc)
{
  Shell_PC_Petsc  *shell;
  PCShellGetContext(pc,(void**)&shell);
  delete shell;
  return 0;
}

// solve Px=b, where P is preconditioner
PetscErrorCode PC_Shell_Apply_Sparse(PC pc,Vec b_in,Vec x_out)
{ 
  Shell_PC_Petsc *shell;
  PCShellGetContext(pc,(void**)&shell);

  NlpGenLinsys  *LinSysMat  = shell->NlpGenLinsysPtr;
  NlpGenData	*data		= shell->NlpGenDataPtr;

  int nb_row = LinSysMat->nx+LinSysMat->my+LinSysMat->mz+LinSysMat->mz;
  int nb_col = nb_row;

  double* x_array, *b_array;
  VecGetArray(b_in,&b_array);
  VecGetArray(x_out,&x_array);

  LinSysMat->copyXSYZ_fromArray( *LinSysMat->res, b_array, nb_col);
  LinSysMat->sol->copyFrom(*LinSysMat->res);
  LinSysMat->solveCompressed( *LinSysMat->sol );
  LinSysMat->copyXSYZ_toArray( *LinSysMat->sol, x_array, nb_col);

  VecRestoreArray(b_in,&b_array);  
  VecRestoreArray(x_out,&x_array);
  return 0;

}











#endif


