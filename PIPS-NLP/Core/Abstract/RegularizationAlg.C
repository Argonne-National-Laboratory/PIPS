/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "RegularizationAlg.h"

RegularizationAlg::RegularizationAlg()
  : 
	DoEvalReg(1),
	MatrixSingular(0),
	ForceReg(false),
	newSystem(true),
	prim_reg_curr(0.0),
	dual_reg_curr(0.0),
    num_PrimReg(0),
    num_DualReg(0)	
{}

RegularizationAlg::~RegularizationAlg(){}



