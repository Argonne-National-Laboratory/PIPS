#include "sNlpInfoFIX.h"
#include "NlpGenVars.h"
#include "OoqpVector.h"
#include <cmath>

#include "SimpleVector.h"
#include "DoubleMatrixHandle.h"
#include "DoubleMatrix.h"

#include "Data.h"

#include "sData.h"

sNlpInfoFIX::~sNlpInfoFIX()
{
}	  

sNlpInfoFIX::sNlpInfoFIX() 
{
}


sNlpInfoFIX::sNlpInfoFIX(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in)
	 :sInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in)	
{
}

sNlpInfoFIX::sNlpInfoFIX(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
					 int nxL_in,int nxU_in,int nsL_in,int nsU_in)
	 : sInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in,nxL_in,nxU_in,nsL_in,nsU_in)
{
}


sNlpInfoFIX::sNlpInfoFIX(sData *data_in)
	:sInfo(data_in)
{
  data_in->inputNlp = this;

}


double sNlpInfoFIX::ObjValue( NlpGenVars * vars) 
{
  StochVector& x = dynamic_cast<StochVector&>(*vars->x);
  OoqpVectorHandle tempX( x.clone() );

  double temp = 0.0; 
  tempX->copyFrom(*g);
 
  Q->mult( 1.0, *tempX, 0.5, *vars->x );
  temp = tempX->dotProductWith( *vars->x );
  
  return temp;
}

void sNlpInfoFIX::ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq)
{  
  conEq->setToZero();
  A->mult( 0.0, *conEq, 1.0, *vars->x );
	
  conIneq->setToZero();
  C->mult( 0.0, *conIneq, 1.0, *vars->x );
}

int sNlpInfoFIX::ObjGrad( NlpGenVars * vars, OoqpVector *grad )
{
  grad->copyFrom(*g);
  Q->mult( 1.0, *grad, 1.0, *vars->x );  
  return 1;
}


void sNlpInfoFIX::Hessian( NlpGenVars * vars, SymMatrix *Hess ){}

void sNlpInfoFIX::JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC) {}

//void sNlpInfoFIX::JacEq( NlpGenVars * vars, GenMatrix* JacA ) {}
void sNlpInfoFIX::JacIneq( NlpGenVars * vars, GenMatrix* JacC ) {}

void sNlpInfoFIX::get_InitX0(OoqpVector* vX){}

void sNlpInfoFIX::createChildren( sData *data_in,stochasticInput& in) {};

void sNlpInfoFIX::Hessian_FromSon( NlpGenVars * vars, double *tempFromParH ) {};

void sNlpInfoFIX::ObjGrad_FromSon( NlpGenVars * vars, OoqpVector *grad, double *tempFromParH ){};


