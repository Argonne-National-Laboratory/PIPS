/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "FilterIPMOption.h"

extern double gconv_tol;
extern int gmax_iter;

FilterIPMOption::FilterIPMOption()
 :
  maxit(gmax_iter),
  opt_tol(gconv_tol),
  tol_mach(1e-16),  
  s_max(100),
  k_ep(10),			  
  k_mu(0.2),		  
  Q_mu(1.5),
  tau_Min(0.99), 		  
  k_sigma(1e10),
  r_Q(1e-5 ),	
  r_phi(1e-8),
  delta(1), 
  r_alpha(0.05),
  powScalar_ConNorm(1.1), 
  powScalar_Obj(2.3), 
  y_phi(1e-8), 
  k_soc(0.99),
  socIter_max(4),
  k_d(1e-5),
  prim_reg_init( 1e-4), 
  prim_reg_min(1e-20),	 
  prim_reg_max(1e20),
  prim_reg_larger_scalar(100), 
  prim_reg_increase_scalar(8),  
  prim_reg_decrease_scalar(1./3.),  
  dual_reg_init(1e-10),  
  dual_reg_scalar(1e-8),
  dual_reg_exp_scalar(0.25),
  mu_init(0.1),
  df(1),
  k_1(1e-2), 
  k_2(1e-2),
  dual_var_init(1),
  move_bound_scalar(1e-8),
  tiny_step_dual(0.01)
{}

FilterIPMOption::~FilterIPMOption(){}


