#ifndef NLPINFOFIX
#define NLPINFOFIX


#include "NlpInfo.h"
#include "OoqpVectorHandle.h"
#include "OoqpVector.h"
#include "SimpleVector.h"

class NlpGenVars;
	
class NlpInfoFIX : public NlpInfo
{
public:

  NlpInfoFIX();    
  
  NlpInfoFIX(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in);
  NlpInfoFIX( int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
				  int nxL_in,int nxU_in,int nsL_in,int nsU_in);  
  
  virtual ~NlpInfoFIX();

  virtual double ObjValue( NlpGenVars * vars);

  virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq);
 
  virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad);
  
  virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess );

  virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC);
  	
  virtual void get_InitX0(OoqpVector* vX);

  SimpleVector *tempX;

};


#endif

