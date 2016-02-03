#include "NlpInfoFIX.h"
#include "getAmplFunction.h"
#include "NlpGenVars.h"
#include "OoqpVector.h"
#include <cmath>

#include "SimpleVector.h"
#include "DoubleMatrixHandle.h"
#include "DoubleMatrix.h"

#include "Data.h"

NlpInfoFIX::~NlpInfoFIX()
{
  if(tempX) delete tempX;
}	  

NlpInfoFIX::NlpInfoFIX(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in)
	 :NlpInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in), tempX(NULL)	
{
  tempX = new SimpleVector(nx_in);
}

NlpInfoFIX::NlpInfoFIX(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
					 int nxL_in,int nxU_in,int nsL_in,int nsU_in)
	 : NlpInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in,nxL_in,nxU_in,nsL_in,nsU_in), tempX(NULL)	
{
  tempX = new SimpleVector(nx_in);
}

NlpInfoFIX::NlpInfoFIX()
	 : tempX(NULL)	
{
}


double NlpInfoFIX::ObjValue( NlpGenVars * vars) 
{
  double temp = 0.0; 
  tempX->copyFrom(*g);
  Q->mult( 1.0, *tempX, 0.5, *vars->x );
  temp = tempX->dotProductWith( *vars->x );
  
  return temp;
}

void NlpInfoFIX::ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq)
{  
  conEq->setToZero();
  A->mult( 0.0, *conEq, 1.0, *vars->x );
	
  conIneq->setToZero();
  C->mult( 0.0, *conIneq, 1.0, *vars->x );
}

int NlpInfoFIX::ObjGrad( NlpGenVars * vars, OoqpVector *grad )
{
  grad->copyFrom(*g);
  Q->mult( 1.0, *grad, 1.0, *vars->x );  
  return 1;
}

void NlpInfoFIX::Hessian( NlpGenVars * vars, SymMatrix *Hess ){}

void NlpInfoFIX::JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC) {}

//void NlpInfoFIX::JacEq( NlpGenVars * vars, GenMatrix* JacA ) {}

//void NlpInfoFIX::JacIneq( NlpGenVars * vars, GenMatrix* JacC ) {}

void NlpInfoFIX::get_InitX0(OoqpVector* vX){}

