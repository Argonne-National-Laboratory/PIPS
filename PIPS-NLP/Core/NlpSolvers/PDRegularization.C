/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "PDRegularization.h"
#include "FilterIPMOption.h"
#include "SolverOption.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

extern int gAssumeMatSingular;
extern int gUseDualRegAlg;

PDRegularization::PDRegularization()
	:
	prim_reg_last(0),
	dual_reg_last(0),
	quantitiesForReg(NULL)
{
  if(gUseDualRegAlg>=1)
  	quantitiesForReg = new double[3];
}

PDRegularization::~PDRegularization()
{
  if(gUseDualRegAlg>=1)
	delete []  quantitiesForReg;
}

PDRegularization::PDRegularization(SolverOption *Option)
{
  FilterIPMOption *FilterIPMOpt = dynamic_cast<FilterIPMOption *>(Option);

  prim_reg_last=0;
  dual_reg_last=0;
  num_PrimReg=0;
  num_DualReg=0;
  
  curr_mu = 0;
  
  prim_reg_init = FilterIPMOpt->prim_reg_init;
  prim_reg_min = FilterIPMOpt->prim_reg_min;
  prim_reg_max = FilterIPMOpt->prim_reg_max;
  
  prim_reg_larger_scalar = FilterIPMOpt->prim_reg_larger_scalar;
  prim_reg_increase_scalar = FilterIPMOpt->prim_reg_increase_scalar;
  prim_reg_decrease_scalar = FilterIPMOpt->prim_reg_decrease_scalar;
  
  prim_reg_max = FilterIPMOpt->prim_reg_max;

  dual_reg_scalar = FilterIPMOpt->dual_reg_scalar;
  dual_reg_exp_scalar = FilterIPMOpt->dual_reg_exp_scalar;
  dual_reg_init = FilterIPMOpt->dual_reg_init;

  if(gUseDualRegAlg>=1)
  	quantitiesForReg = new double[3];  
}


int
PDRegularization::newLinearSystem()
{   
  if (prim_reg_curr > 0.) {
	prim_reg_last = prim_reg_curr;
  }
  if (dual_reg_curr > 0.) {
    dual_reg_last = dual_reg_curr;
  }

  prim_reg_curr = 0.;
  dual_reg_curr = 0.;

  ForceReg=false;
  newSystem=true;
  curr_mu = 0;

  if (gUseDualRegAlg!=0){ 
	dual_reg_curr = dual_reg_init;
	ForceReg=true;
  }else{
    if( MatrixSingular==1 && gAssumeMatSingular==1 ){
	  ForceReg=true;
    }else{  
	  MatrixSingular = 0;
	}
  }
  
  return 1;
}


int
PDRegularization::computeRegularization(double &priReg, double &dualReg,  const double mu)
{
  newSystem = false;  
  curr_mu = mu;
  
  int flag = 1;
  
  if(MatrixSingular==1){
	computeReg_Singularity();
	num_DualReg++;
  }
  else{
  	computeReg_WrongInertia();
	num_PrimReg++;
  }

  priReg  = prim_reg_curr;
  dualReg = dual_reg_curr; 

  return flag;
}


void
PDRegularization::computeReg_WrongInertia()
{
  prim_reg_curr = PriRegularization();
}

void
PDRegularization::computeReg_Singularity()
{  
  if (dual_reg_curr > 0.) { 
    prim_reg_curr = PriRegularization();
  }
  else{
    dual_reg_curr = DualRegularization();
  }
}


double
PDRegularization::PriRegularization()
{
  if(prim_reg_curr == 0.) 
  {
	if (prim_reg_last == 0.) 
	  prim_reg_curr = prim_reg_init;
	else 
	  prim_reg_curr = (prim_reg_min > prim_reg_last*prim_reg_decrease_scalar)
	  					?prim_reg_min:prim_reg_last*prim_reg_decrease_scalar;
  }
  else
  {
	if (prim_reg_last == 0. || 1e5*prim_reg_last<prim_reg_curr) 
	  prim_reg_curr = prim_reg_larger_scalar*prim_reg_curr;
	else
	  prim_reg_curr = prim_reg_increase_scalar*prim_reg_curr;
  }
  
  if (prim_reg_curr > prim_reg_max)
  {
	printf("primal regularization is becoming too large: %e\n",prim_reg_curr);
	assert("Wrong Matrix!" && 0);
	return 0;
  }
  
  return prim_reg_curr;
}


double
PDRegularization::DualRegularization()
{ 
  double dual_reg_temp;

  if(gUseDualRegAlg >= 1){
	dual_reg_temp = computeDualRegFromCons();
  }else {
	dual_reg_temp = dual_reg_scalar * pow(curr_mu, dual_reg_exp_scalar);
  }
  dual_reg_temp = 1e-10;
  return dual_reg_temp;
}

double
PDRegularization::computeDualRegFromCons()
{ 
  double dual_reg_temp;

  assert(quantitiesForReg);
  double gamnrmp = quantitiesForReg[0];
  double nrmc = quantitiesForReg[1];
  double dAAd = quantitiesForReg[2];

  if (gamnrmp<=0.99*nrmc || dAAd <=0 ){ // this implies ||Ad+c||<=0.4*||c|| if exact
    dual_reg_temp = dual_reg_curr;
  }else{
    // if system is singular
	dual_reg_temp = 1e-16;
	if(dual_reg_temp<0.1*dual_reg_curr) 
	  dual_reg_temp = 0.1*dual_reg_curr;  
  }

  return dual_reg_temp;
}


