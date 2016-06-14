/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#include "sInfo.h"
#include "sTree.h"
#include "sData.h"

#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"

sInfo::~sInfo()
{
	destroyChildren();
}	  

sInfo::sInfo()
{

}

sInfo::sInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in)
	 :NlpInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in)	
{
}

sInfo::sInfo(int nx_in,int my_in,int mz_in,int nzH_in,int nzA_in,int nzC_in,
					 int nxL_in,int nxU_in,int nsL_in,int nsU_in)
	 : NlpInfo(nx_in, my_in,mz_in,nzH_in,nzA_in,nzC_in,nxL_in,nxU_in,nsL_in,nsU_in)
{
}


sInfo::sInfo(sData *data_in)
	 : NlpInfo()
	   
{
  stochNode = data_in->stochNode;
  mpiComm = stochNode->commWrkrs;

  locNx = data_in->getLocalnx();
  locMy = data_in->getLocalmy();
  locMz = data_in->getLocalmz();

  A = data_in->Jeq;
  Amat = &(data_in->getLocalA());
  Bmat = &data_in->getLocalB();
  C = data_in->Jineq;
  Cmat = &data_in->getLocalC();
  Dmat = &data_in->getLocalD();  
  Emat = &data_in->getLocalE();
  Fmat = &data_in->getLocalF();
  
  Q = data_in->H;
  Qdiag = &data_in->getLocalQ();
  Qborder = &data_in->getLocalCrossHessian(); 


  g = data_in->grad; 

  int tempH,tempA,tempC;
  data_in->getLocalNnz(tempH,tempA,tempC);

  nzH=tempH;
  nzA=tempA;
  nzC=tempC;

  nx=data_in->nx;
  my=data_in->my;
  mz=data_in->mz;
  
  data_in->inputNlp = this;

}

void sInfo::destroyChildren()
{
  for(size_t it=0; it<children.size(); it++) {
    children[it]->destroyChildren();
    delete children[it];
  }
  children.clear();
}

void sInfo::AddChild(sInfo* child)
{
  children.push_back(child);
}

void sInfo::Emult ( double beta,  OoqpVector& y,
		     double alpha, OoqpVector& x )
{
  Emat->mult(beta, y, alpha, x);
}

void sInfo::Fmult ( double beta,  OoqpVector& y,
		    double alpha, OoqpVector& x )
{
  Fmat->mult(beta, y, alpha, x);
}

//
//void sInfo::createChildren(sData *data_in)
//{
//
//  for (size_t it=0; it<data_in->children.size(); it++) {
//  	if(stochNode->children[it]->commWrkrs != MPI_COMM_NULL)
//      AddChild( new sInfo( data_in->children[it]));
//	else
//	  AddChild( new sInfoDummy());
//  }
//}


