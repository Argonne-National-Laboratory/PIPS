/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "NlpInfo.h"

NlpInfo::~NlpInfo()
{
}	  

NlpInfo::NlpInfo()
	 : nx(0),
	   my(0),
	   mz(0),
	   nzA(0),
	   nzC(0),
	   nzH(0),
	   nsU(0),
	   nxL(0),
	   nxU(0)
	   
{
}


 NlpInfo::NlpInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in)
	 : nx(nx_in),
	   my(my_in),
	   mz(mz_in),
	   nzA(nzA_in),
	   nzC(nzC_in),
	   nzH(nzH_in),
	   nsU(0),
	   nxL(0),
	   nxU(0)
 {
 }

 NlpInfo::NlpInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
					 int nxL_in,int nxU_in,int nsL_in,int nsU_in)
	 : nx(nx_in),
	   my(my_in),
	   mz(mz_in),
	   nzA(nzA_in),
	   nzC(nzC_in),
	   nzH(nzH_in),
	   nsL(nsL_in),
	   nsU(nsU_in),
	   nxL(nxL_in),
	   nxU(nxU_in)
 {
 }


 double NlpInfo::ObjValue( NlpGenVars * vars) {assert( "Not supported" && 0 );}

 //note that now ceqbody = c(x)
 void NlpInfo::ConstraintBody( NlpGenVars * vars, OoqpVector *conEq, OoqpVector *conIneq) {assert( "Not supported" && 0 );}

 int NlpInfo::ObjGrad( NlpGenVars * vars, OoqpVector *grad ) {assert( "Not supported" && 0 );}


 void NlpInfo::Hessian( NlpGenVars * vars, SymMatrix *Hess ){assert( "Not supported" && 0 );}

 void NlpInfo::JacFull( NlpGenVars * vars, GenMatrix* JacA, GenMatrix* JacC) {assert( "Not supported" && 0 );}

 void NlpInfo::get_InitX0(OoqpVector* vX)
  {assert( "Not supported" && 0 );}
 
