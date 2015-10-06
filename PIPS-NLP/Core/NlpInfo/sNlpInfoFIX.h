#ifndef NLPINFOFIX_STOCH
#define NLPINFOFIX_STOCH

#include "sInfo.h"
#include "OoqpVectorHandle.h"
#include "OoqpVector.h"

class NlpGenVars;
	
class sNlpInfoFIX : public sInfo
{
public:

  virtual ~sNlpInfoFIX();


  sNlpInfoFIX();	
  
  
  sNlpInfoFIX(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in);
  sNlpInfoFIX( int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
				  int nxL_in,int nxU_in,int nsL_in,int nsU_in);  
  
  sNlpInfoFIX(sData *data_in);   
  sNlpInfoFIX(sData *data_in, stochasticInput &in){assert(0);}

  virtual double ObjValue( NlpGenVars * vars) ;
  
  virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq);

  virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad );

  virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess );

  virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC);

  virtual void JacEq( NlpGenVars * vars, GenMatrix* JacA );

  virtual void JacIneq( NlpGenVars * vars, GenMatrix* JacC );

  virtual void get_InitX0(OoqpVector* vX);

};


#endif


