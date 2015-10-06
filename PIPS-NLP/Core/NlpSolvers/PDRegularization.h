/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/
 

#ifndef PDREGALG_H
#define PDREGALG_H

#include "RegularizationAlg.h"

// primal-dual regularization
class PDRegularization: public RegularizationAlg{

protected:
  double prim_reg_min,prim_reg_max,
  		 prim_reg_larger_scalar,prim_reg_increase_scalar, prim_reg_decrease_scalar,
  		 dual_reg_scalar,dual_reg_exp_scalar;

  double curr_mu;

public:

  double *quantitiesForReg;
  
  double prim_reg_init;
  double prim_reg_last;
  double dual_reg_init;
  double dual_reg_last;

  PDRegularization();
  PDRegularization(SolverOption *Option);

  virtual ~PDRegularization();
  
  virtual int
  newLinearSystem();

  virtual int
  computeRegularization(double &priReg, double &dualReg, const double mu);

private:

  virtual double
  PriRegularization();

  virtual double
  DualRegularization();

  void computeReg_WrongInertia();
  void computeReg_Singularity();

  double computeDualRegFromCons();

};

#endif
