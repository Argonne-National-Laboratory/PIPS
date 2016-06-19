/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/


#ifndef GETAMPLFUNC_H
#define GETAMPLFUNC_H

#include "asl_pfgh.h"

enum{
Fixed = -1,
NoBound = 0,
LowBound = 1, 
UppBound = 2,
BothBound = 3,
};

//
void 
doubleLexSort_Init( int first[], int n, int second[], double data[], const int ifEqCon);

void 
doubleLexSort_ValOnly(int nzA, double *dataA, int nzC, double *dataC, double *dataJacFull);

void
ampl_GetHessianInit(double *varsX, double *Helts,double *Yelts=NULL);

int 
ampl_get_nnz_Hessian_Tri();

void
ampl_get_Hessian_Tri(double *varsX, double *Helts, double *Yelts=NULL,double *Zelts=NULL);

void
ampl_get_Jac(double *varsX, const int nzA, double *JacAelts, const int nzC, double *JacCelts);

void 
ampl_get_bounds( 
						  double xlow[], int nx, char ixlow[],
						  double xupp[], char ixupp[],
						  double b[], int /* my */,
						  double clow[], int /* mz */, char iclow[],
						  double cupp[], char icupp[] );

void 
ampl_get_matrices(  fint irow[], fint kcol[], double elts[],
				int /* nx */, int nnzQ,
				int  my , int nnzA,
				int  mz , int nnzC,
				int irowQ[], int jcolQ[], double dQ[],
				int irowA[], int jcolA[], double dA[],
				int irowC[], int jcolC[], double dC[], 
				double *dwrkX, 
				const int full_size );

double 
ampl_get_Obj(double *varsX);

void 
ampl_get_ObjGrad( double *varsX,  double Objgrad[]);

void 
ampl_get_Cons(double *varsX, double *consEqElt, double *consIneqElt);

void
ampl_count_sizes_SplitSlack(  fint irow[], fint kcol[],
			int & nx, int & nnzQ,
			int & my, int & nnzA, int & mz, int & nnzC, 
			const int full_size, 
			int &nnzCL, int & nnzCU,
			int & nxL, int & nxU, int & nsL, int & nsU);


void 
ampl_get_InitX0(double *varsX);

void ampl_write_solution(double *varsX, double *Yelts, double *Zelts);

void
ampl_free_mapinfo();

#endif
