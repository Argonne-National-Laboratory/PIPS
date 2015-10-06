/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef REGALG_H
#define REGALG_H

#include <iostream>

class SolverOption;

class RegularizationAlg {

public:

  int DoEvalReg;
  bool ForceReg;
  bool newSystem;
  int MatrixSingular;

  double prim_reg_curr;
  double dual_reg_curr;

  int num_PrimReg;
  int num_DualReg;
  
  RegularizationAlg();

  virtual ~RegularizationAlg();
  
  virtual int
  newLinearSystem()=0;

  virtual int
  computeRegularization(double &priReg, double &dualReg, const double mu=0)=0;

};

#endif

