/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef FILTERIPMOPTION_H
#define FILTERIPMOPTION_H

#include "SolverOption.h"

/** 
 * options for Filter line search algorithm 
 */
class FilterIPMOption : public SolverOption
{

  public:

  int maxit;   					//max lPM iterations
  double opt_tol;				// termination tol

  double move_bound_scalar;
  double mu_init;
  double tol_mach;	  
  double s_max;
  double k_ep;			  
  double k_mu;		  
  double Q_mu;
  double tau_Min; 		  
  double k_sigma;
  double r_Q;	
  double r_phi;
  double delta; 
  double r_alpha;
  double powScalar_ConNorm; 
  double powScalar_Obj; 
  double y_phi;
  double k_soc;
  double ConNorm_min;
  double ConNorm_max;
  double k_d;
  int socIter_max;
  
  double tiny_step_dual;

  // for primal regularization:
  double prim_reg_init;
  double prim_reg_min;   
  double prim_reg_max;
  double prim_reg_larger_scalar;
  double prim_reg_increase_scalar;
  double prim_reg_decrease_scalar;

  // for dual regularization:
  double dual_reg_init;
  double dual_reg_scalar;    
  double dual_reg_exp_scalar;

  //for initialization
  double k_1, k_2;
  double dual_var_init;

  // scaling parameter for objective
  double df;				  


  FilterIPMOption();
  ~FilterIPMOption();


};


#endif
