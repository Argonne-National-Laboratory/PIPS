/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "sTreeImpl.h"

#include "StochVector.h"
#include "StochGenMatrix.h"
#include "StochSymMatrix.h"

sTreeImpl::sTreeImpl( stochasticInput &in_)
  : sTree(), m_id(0), in(in_)
{
  if(-1==rankMe) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
  if(-1==numProcs) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  m_nx = in.nFirstStageVars();
  m_my = compute_nFirstStageEq();
  m_mz = in.nFirstStageCons() - m_my;

  for (int scen=0; scen<in.nScenarios(); scen++)
    children.push_back(new sTreeImpl(scen+1,in));

}

sTreeImpl::sTreeImpl(int id, stochasticInput &in_)
  : sTree(), m_id(id), in(in_)
{ 
  m_nx = in.nSecondStageVars(id-1);
  m_my = compute_nSecondStageEq(id-1);
  m_mz = in.nSecondStageCons(id-1) - m_my;
}
  
sTreeImpl::~sTreeImpl()
{ }


StochSymMatrix*   sTreeImpl::createQ() const
{
  return NULL;
}

  StochVector*      sTreeImpl::createc() const
{
  return NULL;
}


StochVector*      sTreeImpl::createxlow()  const
{
  return NULL;
}

StochVector*      sTreeImpl::createixlow() const
{
  return NULL;
}

StochVector*      sTreeImpl::createxupp()  const
{
  return NULL;
}

StochVector*      sTreeImpl::createixupp() const
{
  return NULL;
}



StochGenMatrix*   sTreeImpl::createA() const
{
  return NULL;
}

StochVector*      sTreeImpl::createb() const
{
  return NULL;
}


StochGenMatrix*   sTreeImpl::createC() const
{
  return NULL;
}

StochVector*      sTreeImpl::createclow()  const
{
  return NULL;
}

StochVector*      sTreeImpl::createiclow() const
{
  return NULL;
}

StochVector* sTreeImpl::createcupp()  const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* cupp = new StochVector(mz(), commWrkrs);
  double* vec = ((SimpleVector*)cupp->vec)->elements();  

  //get the data from the stochasticInput
  vector<double> x;
  if(m_id==0)
    x=in.getFirstStageColUB();
  else 
    x=in.getSecondStageColUB(m_id-1);
    
  for(size_t i=0; i<m_mz; i++)
    if(x[i]<1e+20) 
      vec[i]=x[i];
    else 
      vec[i]=0.0;
  
  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->createcupp();
    cupp->AddChild(child);
  }
  return cupp;
}

StochVector* sTreeImpl::createicupp() const
{
  //is this node a dead-end for this process?
  if(commWrkrs==MPI_COMM_NULL)
    return new StochDummyVector();

  StochVector* icupp = new StochVector(mz(), commWrkrs);
  double* vec = ((SimpleVector*)icupp->vec)->elements();  

  //get the data from the stochasticInput
  vector<double> x;
  if(m_id==0)
    x=in.getFirstStageColUB();
  else 
    x=in.getSecondStageColUB(m_id-1);
    
  for(size_t i=0; i<m_mz; i++)
    if(x[i]<1e+20) 
      vec[i]=1.0;
    else 
      vec[i]=0.0;
  
  for(size_t it=0; it<children.size(); it++) {
    StochVector* child = children[it]->createicupp();
    icupp->AddChild(child);
  }
  return icupp;
}


int sTreeImpl::nx() const
{
  return m_nx;
}

int sTreeImpl::my() const
{
  return m_my;
}
 
int sTreeImpl::mz() const
{
  return m_mz;
}
 
int sTreeImpl::id() const
{
  return m_id;
}
 

void sTreeImpl::computeGlobalSizes()
{
}

int sTreeImpl::compute_nFirstStageEq()
{
  int num=0;
  vector<double> lb=in.getFirstStageRowLB();
  vector<double> ub=in.getFirstStageRowUB();

  for (size_t i=0;i<lb.size(); i++)
    if (lb[i]==ub[i]) num++;

  return num;
}

int sTreeImpl::compute_nSecondStageEq(int scen)
{
  int num=0;
  vector<double> lb=in.getSecondStageRowLB(scen);
  vector<double> ub=in.getSecondStageRowUB(scen);

  for (size_t i=0;i<lb.size(); i++)
    if (lb[i]==ub[i]) num++;

  return num;
}
