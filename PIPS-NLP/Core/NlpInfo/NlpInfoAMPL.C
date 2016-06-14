/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "NlpInfoAMPL.h"
#include "getAmplFunction.h"
#include "NlpGenVars.h"
#include "OoqpVector.h"
#include <cmath>

#include "SimpleVector.h"
#include "DoubleMatrixHandle.h"
#include "DoubleMatrix.h"


#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"

#ifdef TIMING
  #include "mpi.h"
  extern double timeFromAMPL;
#endif


NlpInfoAMPL::NlpInfoAMPL()
{
}

NlpInfoAMPL::~NlpInfoAMPL()
{
}

NlpInfoAMPL::NlpInfoAMPL(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in)
	: NlpInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in)
{
}


NlpInfoAMPL::NlpInfoAMPL(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
					int nxL_in,int nxU_in,int nsL_in,int nsU_in)
	: NlpInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in,nxL_in,nxU_in,nsL_in,nsU_in)
{
}



double 
NlpInfoAMPL::ObjValue( NlpGenVars * vars)
{

#ifdef TIMING
	  double tTot=MPI_Wtime();
#endif


  double *tempDbl = (double*) malloc(nx*sizeof(double));
  double objWrk;
  vars->x->copyIntoArray(tempDbl);

  objWrk = ampl_get_Obj( tempDbl );
  
  free(tempDbl);

#ifdef TIMING
	  timeFromAMPL += MPI_Wtime()-tTot;
#endif  
  
  return objWrk;
}



int 
NlpInfoAMPL::ObjGrad( NlpGenVars * vars, OoqpVector *grad)
{
#ifdef TIMING
		  double tTot=MPI_Wtime();
#endif

  double *tempX = (double*) malloc(nx*sizeof(double));
  double *tempGrad = (double*) malloc(nx*sizeof(double));

  vars->x->copyIntoArray(tempX);

  ampl_get_ObjGrad( tempX, tempGrad);

  grad->copyFromArray(tempGrad);

  free(tempX);
  free(tempGrad);
  
#ifdef TIMING
		timeFromAMPL += MPI_Wtime()-tTot;
#endif  

  return 1;

}


void 
NlpInfoAMPL::ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq)
{
#ifdef TIMING
		  double tTot=MPI_Wtime();
#endif

  int i;
  double *tempX = (double*) malloc(nx*sizeof(double));
  double *tempConEq = (double*) malloc(my*sizeof(double));
  double *tempConInEq = (double*) malloc(mz*sizeof(double));

  vars->x->copyIntoArray(tempX);

  ampl_get_Cons(tempX, tempConEq, tempConInEq);

  conEq->copyFromArray(tempConEq);
  conIneq->copyFromArray(tempConInEq);
 

  free(tempX);
  free(tempConEq);  
  free(tempConInEq);
#ifdef TIMING
		timeFromAMPL += MPI_Wtime()-tTot;
#endif  

}


//make sure we are using Hessian in triangular form!
void NlpInfoAMPL::Hessian( NlpGenVars * vars, SymMatrix* Hess )
{
#ifdef TIMING
		  double tTot=MPI_Wtime();
#endif

   double *tempX = (double*) malloc(nx*sizeof(double));
   double *tempY = (double*) malloc(my*sizeof(double));   
   double *tempZ = (double*) malloc(mz*sizeof(double));    
   double *tempH = (double*) malloc(nzH*sizeof(double));
 
   vars->x->copyIntoArray(tempX);
   vars->y->negate();
   vars->y->copyIntoArray(tempY);
   vars->y->negate();
   vars->z->negate();
   vars->z->copyIntoArray(tempZ);
   vars->z->negate();
   
   ampl_get_Hessian_Tri( tempX, tempH, tempY,tempZ);
 
   Hess->copyMtxFromDouble(nzH,tempH);
   
   free(tempH);
   free(tempZ);
   free(tempY);   
   free(tempX);
#ifdef TIMING
	  timeFromAMPL += MPI_Wtime()-tTot;
#endif     
}


void NlpInfoAMPL::JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC)
{
#ifdef TIMING
		  double tTot=MPI_Wtime();
#endif

  double *tempX = (double*) malloc(nx*sizeof(double));
  double *tempA = (double*) malloc(nzA*sizeof(double));
  double *tempC = (double*) malloc(nzC*sizeof(double));
  
  vars->x->copyIntoArray(tempX);

  ampl_get_Jac( tempX, nzA, tempA, nzC, tempC);

  JacA->copyMtxFromDouble(nzA,tempA);
  JacC->copyMtxFromDouble(nzC,tempC);
  
  free(tempX);
  free(tempA);
  free(tempC);
#ifdef TIMING
		timeFromAMPL += MPI_Wtime()-tTot;
#endif  

}


void
NlpInfoAMPL::get_InitX0(OoqpVector* vX)
{
#ifdef TIMING
		  double tTot=MPI_Wtime();
#endif

  double *tempX = (double*) malloc(nx*sizeof(double));
  
  ampl_get_InitX0(tempX);

  vX->copyFromArray(tempX);

#ifdef TIMING
		timeFromAMPL += MPI_Wtime()-tTot;
#endif   

  free(tempX);
}

