/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLPINFO
#define NLPINFO


#include "OoqpVectorHandle.h"
#include "OoqpVector.h"
#include "DoubleMatrixHandle.h"

class NlpGenVars;
class Variables;

class NlpInfo
{
public:
  
  NlpInfo();    
  virtual ~NlpInfo();

  long long  nx,my,mz;
  long long nzH,nzA,nzC;
  long long nsL, nsU, nxL, nxU;
  int *rowMap;
  
  NlpInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in);
  NlpInfo( int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
				  int nxL_in,int nxU_in,int nsL_in,int nsU_in);  
  

  

  virtual double ObjValue( NlpGenVars * vars) = 0;
  
  virtual void ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq) = 0;

  virtual int ObjGrad( NlpGenVars * vars, OoqpVector *grad ) = 0;

  

  virtual void Hessian( NlpGenVars * vars, SymMatrix *Hess ) = 0;

  virtual void JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC) = 0;

  virtual void get_InitX0(OoqpVector* vX) = 0;


  SymMatrix *Q;
  GenMatrix *A;
  GenMatrix *C;
  OoqpVector *g;
  OoqpVector *bA;


};

#endif
