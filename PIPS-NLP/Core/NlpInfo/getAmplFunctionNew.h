/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef AMPLFUNCNEW_H
#define AMPLFUNCNEW_H

#include "asl_pfgh.h"

double Ampl_Eval_Obj(ASL_pfgh *asl_, double *varsX);
void Ampl_Eval_Cons(ASL_pfgh *asl_, double *varsX, double *cons);
void Ampl_Eval_ObjGrad( ASL_pfgh *asl_, double *varsX,  double *Objgrad);

void Ampl_Eval_InitX0(ASL_pfgh *asl_, double *varsX);

void Ampl_Eval_Jac(ASL_pfgh *asl_, double *varsX, double *AmplJacElts);

void Ampl_Eval_Hessian_Tri(ASL_pfgh *asl_, double *varDual, double *Helts, double ObjScale=1);
	

#endif

