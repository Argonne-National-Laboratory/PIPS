/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include <cstdio>
#include <cassert>
#include <climits>

#include "getAmplFunctionNew.h"

extern int nobj;
extern int hes_obj;
extern int hes_con;
extern int hes_tri;

extern double OW[1];

double Ampl_Eval_Obj(ASL_pfgh *asl_, double *varsX)
{
  ASL_pfgh *asl = asl_;
  
  double objv=0;
  fint nerror = -1;

  if(varsX){
    objv = objval(nobj, varsX, &nerror);	
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
    objv = objval(nobj, dxWrk, &nerror);
    free(dxWrk);
  }
  return objv;
}


void Ampl_Eval_Cons(ASL_pfgh *asl_, double *varsX, double *cons)
{
  ASL_pfgh *asl = asl_;

  fint nerror = -1;

  if(varsX){
  	conval(varsX, cons, &nerror);	
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
  	conval(varsX, cons, &nerror);	
    free(dxWrk);
  }
}


void Ampl_Eval_ObjGrad( ASL_pfgh *asl_, double *varsX,  double *Objgrad)
{
  ASL_pfgh *asl = asl_;
  
  fint nerror = -1;

  if(varsX){
	objgrd(nobj, varsX, Objgrad, &nerror);  
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
	objgrd(nobj, dxWrk, Objgrad, &nerror);	
    free(dxWrk);
  }
 
}


void Ampl_Eval_InitX0(ASL_pfgh *asl_, double *varsX)
{
  ASL_pfgh *asl = asl_;

  for (int i=0; i<n_var; i++)
	varsX[i] = havex0[i] ? X0[i] : 0;
}


void Ampl_Eval_Jac(ASL_pfgh *asl_, double *varsX, double *AmplJacElts)
{
  ASL_pfgh *asl = asl_;

  fint nerror = -1;

  if(varsX){
	jacval(varsX, AmplJacElts, &nerror);
  }else{
    double *dxWrk = (double*) malloc (n_var*sizeof(double));
  	jacval(varsX, AmplJacElts, &nerror);	
    free(dxWrk);
  }
}

void Ampl_Eval_Hessian_Tri(ASL_pfgh *asl_, double *varDual, double *Helts, double ObjScale)
{
  ASL_pfgh *asl = asl_;
  
  OW[0] = ObjScale;

  sphes(Helts, -1, OW, varDual);
}


void ampl_write_solution(ASL_pfgh *asl_, double *varsX, double *dual)
{
  ASL_pfgh *asl = asl_;
  write_sol("\pipsnlp_parallel:", varsX, dual,NULL);
}
