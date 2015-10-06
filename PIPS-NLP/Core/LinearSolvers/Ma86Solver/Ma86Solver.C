/* PIPS-IPM                                                             
 * Author: Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
 
 /* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include <stdlib.h>
#include <iostream>

using namespace std;
#include "Ma86Solver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"

#include <cstdlib>
#include "time.h"
#include "stdio.h" 
#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

extern int gOoqpPrintLevel;


/* Initialise control with default values */
extern "C" void ma86_default_control_d(struct ma86_control_d *control);

extern "C" void ma86_analyse_d(const int n, const int ptr[], const int row[], int order[], void **keep, const struct ma86_control_d *control, struct ma86_info_d *info);

extern "C" void ma86_factor_d(const int n, const int ptr[], const int row[], const ma86pkgtype_d_ val[], const int order[], void **keep, const struct ma86_control_d *control, struct ma86_info_d *info, const ma86realtype_d_ scale[]);

/* To factorize the matrix AND solve AX = B */
extern "C" void ma86_factor_solve_d(const int n, const int ptr[], const int row[], const ma86pkgtype_d_ val[], const int order[], void **keep, const struct ma86_control_d *control, struct ma86_info_d *info, const int nrhs, const int ldx, ma86pkgtype_d_ x[], const ma86realtype_d_ scale[]);

/* To solve AX = B using the computed factors */
extern "C" void ma86_solve_d(const int job, const int nrhs, const int ldx, ma86pkgtype_d_ *x, const int order[], void **keep, const struct ma86_control_d *control, struct ma86_info_d *info, const ma86realtype_d_ scale[]);
/* To clean up memory in keep */
extern "C" void ma86_finalise_d(void **keep, const struct ma86_control_d *control);
                
////////////////////////////////////////////////////////////////////
/* Set default values for control struct */
extern "C" void mc68_default_control(struct mc68_control *control);
/* Perform ordering */
extern "C" void mc68_order(const int ord, const int n, const int ptr[], const int row[], int perm[], const struct mc68_control *control, struct mc68_info *info);

////////////////////////////////////////////////////////////////////

Ma86Solver::Ma86Solver( SparseSymMatrix * sgm )
{

  Msys = sgm;
	// Construct a copy matrix to work on.
  n = sgm->size();
  nnz = sgm->numberOfNonZeros();

	// give some arbitrary pivot ordering, it will be changed later with mc68_order

  order = new int[n];
  for (int i = 0; i < n; ++i) {
    order[i] = i;
  }


  krowM = new int[n+1];
  jcolM = new int[nnz];
  M     = new double[nnz];
  
  first = true;
  second = true;
  ma86_default_control_d(&control);
  mc68_default_control(&control68);
  
  control.scaling=1;

}

void Ma86Solver::firstCall()
{	
	
  Msys->getStorageRef().transpose(krowM, jcolM, M);  
  // sparse matrix storage adjustment for PIPS	
	
  ptr = krowM;
  row = jcolM;
  val = M;

  ////////////////////////////////////////////////////////////////
  mc68_order(3, n, ptr, row, order, &control68, &info68);

 // mc68_order(3,...) uses Metis with parameter 3
  ////////////////////////////////////////////////////////////////
	
  
  ma86_analyse_d(n, ptr, row, order, &keep, &control, &info);
   
  if(info.flag < 0) {
    printf("Failure during analyse with info .flag= %i\n", info.flag);
  }
	
} 

void Ma86Solver::diagonalChanged( int /* idiag */, int /* extent */ )

{
	this->matrixChanged();
}

int Ma86Solver::matrixChanged()
{  
  if (first) {
    firstCall();  
    first=false;
   }
  int matrixSingular = 0;

  Msys->getStorageRef().transpose(krowM, jcolM, M); 


  // Factorize

	//clock_t startTime = clock();
  // if we use clock(), the time is wrong when the   OMP_NUM_THREADS>=2; instead use omp_get_wtime() in the driver file

	ma86_factor_d(n, ptr, row, val, order, &keep, &control, &info, NULL);
  if (info.flag < 0) {
    printf("Failure during factor with info. flag = %i\n", info.flag);
  }

//	clock_t endTime = clock();
 
//	clock_t clockTicksTaken = endTime - startTime;

//	double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;

//	cout<<"factorize time = " << timeInSeconds<< endl ;

  if(info.flag = -3) 
  	matrixSingular=1;
  else
    assert( info.flag >= 0 );
  
  if(matrixSingular==1) 
    negEigVal = -1;
  else
    negEigVal = info.num_neg;

return negEigVal;


}
 
	// Solve
void Ma86Solver::solve( OoqpVector& rhs_in )
{
  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
  double * sol = rhs.elements();


  // Solve

	ma86_solve_d(0, 1, n, sol, order, &keep, &control, &info, NULL);
  if (info.flag < 0) 
	{
    printf("Failure during factor with info. flag = %i\n", info.flag);	
  }

}


void Ma86Solver::solve(GenMatrix& rhs_in)
{ 
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int nrows;
  int ncols;
  rhs.getSize(ncols,nrows);
  double * sol = new double[nrows*ncols];
  for (int i = 0; i < n; ++i) {
    sol[i] = rhs[i][0];
  }
  assert(nrows==n);

  // Solve
  ma86_solve_d(0, 1, n, sol, order, &keep, &control, &info, NULL);
  if (info.flag < 0) {
    printf("Failure during factor with info. flag = %i\n", info.flag);	
	}

	memcpy(&rhs[0][0], sol, nrows*ncols*sizeof(double));

  delete [] sol;

	//clock_t endTime = clock();
	//clock_t clockTicksTaken = endTime - startTime;

	//double timeInSeconds = clockTicksTaken / (double) CLOCKS_PER_SEC;

	//cout<<"SOLVE TIME FOR GenMatrix = " << timeInSeconds<< endl ;

}  


void Ma86Solver::solve(SimpleVector& rhs_in)
{ 

  x = rhs_in.elements();

  // Solve

	ma86_solve_d(0, 1, n, x, order, &keep, &control, &info, NULL);
  if (info.flag < 0) {
    printf("Failure during factor with info. flag = %i\n", info.flag);	
  }

  rhs_in.copyFromArray(x);
}


Ma86Solver::~Ma86Solver()
{

  ma86_finalise_d(&keep,&control);
  delete[] order;
  delete[] jcolM;
  delete[] krowM;
  delete[] M;


}




