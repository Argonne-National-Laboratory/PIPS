/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPINFOAMPL
#define NLPINFOAMPL


#include "NlpInfo.h"
#include "OoqpVectorHandle.h"
#include "OoqpVector.h"

class NlpGenVars;


class NlpInfoAMPL : public NlpInfo
{
public:
  
  NlpInfoAMPL();    
  NlpInfoAMPL(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in);
  NlpInfoAMPL( int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
				  int nxL_in,int nxU_in,int nsL_in,int nsU_in);  
  
  virtual ~NlpInfoAMPL();

  virtual double ObjValue( NlpGenVars * vars);

  virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq);
 
  virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad);
  
  virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess );

  virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC);
  	
  virtual void get_InitX0(OoqpVector* vX);



};



#endif
