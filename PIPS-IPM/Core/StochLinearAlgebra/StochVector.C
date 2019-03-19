#include "StochVector.h"
#include "StochTree.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "VectorUtilities.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <math.h>

StochVector::StochVector(int n_, MPI_Comm mpiComm_, int isDistributed/*=-1*/)
  : OoqpVector(n_), vecl(NULL), parent(NULL), mpiComm(mpiComm_),
    iAmDistrib(isDistributed)
{
  vec = new SimpleVector(n_);

  if(-1==iAmDistrib && MPI_COMM_NULL!=mpiComm) {
    int size;
    MPI_Comm_size(mpiComm, &size);
    iAmDistrib = size==1?0:1;
  }
  vecl = NULL;
}

StochVector::StochVector(int n_, int nl_, MPI_Comm mpiComm_, int isDistributed)
  : OoqpVector(n_), parent(NULL), mpiComm(mpiComm_),
    iAmDistrib(isDistributed)
{
  vec = new SimpleVector(n_);

  if( nl_ >= 0 )
    vecl = new SimpleVector(nl_);
  else
	 vecl = NULL;

  if(-1==iAmDistrib && MPI_COMM_NULL!=mpiComm) {
    int size;
    MPI_Comm_size(mpiComm, &size);
    iAmDistrib = size==1?0:1;
  }

}

void StochVector::AddChild(StochVector* child)
{
  child->parent = this;
  children.push_back(child);
  n += child->n;
}

void StochVector::AddChild(OoqpVector* child_)
{
  StochVector* child = reinterpret_cast<StochVector*>(child_);
  AddChild(child);
}

StochVector::~StochVector()
{
  for (size_t it=0; it<children.size(); it++)
    delete children[it];
  
  if( vec )
    delete vec;

  if( vecl )
	 delete vecl;
}

OoqpVector* StochVector::dataClone() const
{
  assert(!vecl);
  OoqpVector* clone = new SimpleVector(vec->length());
  return clone;
}

OoqpVector* StochVector::dataCloneLinkCons() const
{
  assert(vecl);
  OoqpVector* clone = new SimpleVector(vecl->length());
  return clone;
}

StochVector* StochVector::clone() const
{
  StochVector* clone;
  if( vecl )
    clone = new StochVector(vec->length(), vecl->length(), mpiComm, -1);
  else
	 clone = new StochVector(vec->length(), mpiComm);

  for(size_t it=0; it<children.size(); it++) {
    clone->AddChild(children[it]->clone());
  }
  return clone;
}

StochVector* StochVector::cloneFull() const
{
   StochVector* clone = new StochVector(vec->length(), (vecl != NULL) ? vecl->length() : -1, mpiComm, -1);

   clone->vec->copyFrom(*vec);

   if( vecl )
      clone->vecl->copyFrom(*vecl);

   for( size_t it = 0; it < children.size(); it++ )
      clone->AddChild(children[it]->cloneFull());

   return clone;
}


void 
StochVector::jointCopyFrom(StochVector& v1, StochVector& v2, StochVector& v3)
{
  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);
  SimpleVector& sv1 = dynamic_cast<SimpleVector&>(*v1.vec);
  SimpleVector& sv2 = dynamic_cast<SimpleVector&>(*v2.vec);
  SimpleVector& sv3 = dynamic_cast<SimpleVector&>(*v3.vec);

  int n1 = sv1.length();
  int n2 = sv2.length();
  int n3 = sv3.length();

  assert(n1+n2+n3 == sv.length());
  
  if(n1>0)
    memcpy(&sv[0], &sv1[0], n1*sizeof(double));

  if(n2>0)
    memcpy(&sv[n1], &sv2[0], n2*sizeof(double));

  if(n3>0)
    memcpy(&sv[n1+n2], &sv3[0], n3*sizeof(double));

  for(size_t it=0; it<children.size(); it++) {
    children[it]->jointCopyFrom(*v1.children[it], 
				*v2.children[it], 
				*v3.children[it]);
  }

}

void
StochVector::jointCopyFromLinkCons(StochVector& vx, StochVector& vy, StochVector& vz)
{
  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);
  SimpleVector& svx = dynamic_cast<SimpleVector&>(*vx.vec);
  SimpleVector& svy = dynamic_cast<SimpleVector&>(*vy.vec);
  SimpleVector& svz = dynamic_cast<SimpleVector&>(*vz.vec);

  int n1 = svx.length();
  int n2 = svy.length();
  int n3 = svz.length();
  int n4 = 0;
  int n5 = 0;

  assert(n1+n2+n3 <= sv.length());
  assert(sizeof(double) == sizeof(sv[0]));

  if(n1>0)
    memcpy(&sv[0], &svx[0], n1*sizeof(double));

  if(n2>0)
    memcpy(&sv[n1], &svy[0], n2*sizeof(double));

  if(n3>0)
    memcpy(&sv[n1+n2], &svz[0], n3*sizeof(double));

  if( vy.vecl )
  {
    SimpleVector& svyl = dynamic_cast<SimpleVector&>(*vy.vecl);
    n4 = svyl.length();
    assert(n4 >= 0);

    if( n4 > 0 )
      memcpy(&sv[n1+n2+n3], &svyl[0], n4*sizeof(double));
  }

  if( vz.vecl )
  {
    SimpleVector& svzl = dynamic_cast<SimpleVector&>(*vz.vecl);
    n5 = svzl.length();
    assert(n5 >= 0);

    if( n5 > 0 )
      memcpy(&sv[n1+n2+n3+n4], &svzl[0], n5*sizeof(double));
  }

  assert(n1+n2+n3+n4+n5 == sv.length());

  for(size_t it=0; it<children.size(); it++) {
    children[it]->jointCopyFromLinkCons(*vx.children[it],
				*vy.children[it],
				*vz.children[it]);
  }
}


void 
StochVector::jointCopyTo(StochVector& v1, StochVector& v2, StochVector& v3)
{
  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);
  SimpleVector& sv1 = dynamic_cast<SimpleVector&>(*v1.vec);
  SimpleVector& sv2 = dynamic_cast<SimpleVector&>(*v2.vec);
  SimpleVector& sv3 = dynamic_cast<SimpleVector&>(*v3.vec);

  int n1 = sv1.length();
  int n2 = sv2.length();
  int n3 = sv3.length();

  assert(n1+n2+n3 == sv.length());
 
  if(n1>0)
    memcpy(&sv1[0], &sv[0], n1*sizeof(double));

  if(n2>0)
    memcpy(&sv2[0], &sv[n1], n2*sizeof(double));

  if(n3>0)
    memcpy(&sv3[0], &sv[n1+n2], n3*sizeof(double));


  for(size_t it=0; it<children.size(); it++) {
    children[it]->jointCopyTo(*v1.children[it], 
			      *v2.children[it], 
			      *v3.children[it]);
  }
}

void
StochVector::jointCopyToLinkCons(StochVector& vx, StochVector& vy, StochVector& vz)
{
  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);
  SimpleVector& svx = dynamic_cast<SimpleVector&>(*vx.vec);
  SimpleVector& svy = dynamic_cast<SimpleVector&>(*vy.vec);
  SimpleVector& svz = dynamic_cast<SimpleVector&>(*vz.vec);

  int n1 = svx.length();
  int n2 = svy.length();
  int n3 = svz.length();
  int n4 = 0;
  int n5 = 0;

  assert(n1+n2+n3 <= sv.length());
  assert(sizeof(double) == sizeof(sv[0]));

  if(n1>0)
    memcpy(&svx[0], &sv[0], n1*sizeof(double));

  if(n2>0)
    memcpy(&svy[0], &sv[n1], n2*sizeof(double));

  if(n3>0)
    memcpy(&svz[0], &sv[n1+n2], n3*sizeof(double));

  if( vy.vecl )
  {
     SimpleVector& svyl = dynamic_cast<SimpleVector&>(*vy.vecl);
     n4 = svyl.length();
     assert(n4 >= 0);

     if( n4 > 0 )
       memcpy(&svyl[0], &sv[n1+n2+n3], n4*sizeof(double));
  }

  if( vz.vecl )
  {
     SimpleVector& svzl = dynamic_cast<SimpleVector&>(*vz.vecl);
     n5 = svzl.length();
     assert(n5>= 0);

     if( n5 > 0 )
       memcpy(&svzl[0], &sv[n1+n2+n3+n4], n5*sizeof(double));
  }

  assert(n1+n2+n3+n4+n5 == sv.length());

  for(size_t it=0; it<children.size(); it++) {
    children[it]->jointCopyToLinkCons(*vx.children[it],
			      *vy.children[it],
			      *vz.children[it]);
  }
}


int StochVector::isKindOf( int kind )
{
  return kind==kStochVector;
}

void StochVector::scale( double alpha )
{
  vec->scale(alpha);
  
  if( vecl ) vecl->scale(alpha);

  for(size_t it=0; it<children.size(); it++)
    children[it]->scale(alpha);
}

void StochVector::setToZero()
{
  vec->setToZero();

  if( vecl ) vecl->setToZero();
  
  for(size_t it=0; it<children.size(); it++)
    children[it]->setToZero();
}

void StochVector::setToConstant( double c) 
{
  vec->setToConstant(c);
  
  if( vecl ) vecl->setToConstant(c);

  for(size_t it=0; it<children.size(); it++)
    children[it]->setToConstant(c);
}

void StochVector::randomize( double alpha, double beta, double *ix )
{
  assert( "Not implemented" && 0 );
}


void StochVector::copyFrom( OoqpVector& v_ )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);

  this->vec->copyFrom(*v.vec);

  if( this->vecl )
  {
     assert(v.vecl);
     this->vecl->copyFrom(*v.vecl);
  }

  //assert tree compatibility
  assert(children.size() == v.children.size());

  for(size_t it=0; it<children.size(); it++)
    children[it]->copyFrom(*v.children[it]);
}

void StochVector::copyFromAbs(const OoqpVector& v_ )
{
  const StochVector& v = dynamic_cast<const StochVector&>(v_);

  this->vec->copyFromAbs(*v.vec);

  if( this->vecl )
  {
     assert(v.vecl);
     this->vecl->copyFromAbs(*v.vecl);
  }

  //assert tree compatibility
  assert(children.size() == v.children.size());

  for(size_t it=0; it<children.size(); it++)
    children[it]->copyFromAbs(*v.children[it]);
}

double StochVector::infnorm()
{
  double infnrm;

  infnrm = 0.0;

  for(size_t it=0; it<children.size(); it++)
    infnrm = std::max(infnrm, children[it]->infnorm());

  if(iAmDistrib) {
    double infnrmG=0.0;
    MPI_Allreduce(&infnrm, &infnrmG, 1, MPI_DOUBLE, MPI_MAX, mpiComm);
    infnrm = infnrmG;
  }

  infnrm = std::max(vec->infnorm(), infnrm);

  if( vecl ) infnrm = std::max(vecl->infnorm(), infnrm);

  return infnrm; 
}

double StochVector::twonorm()
{
#if 0
  return sqrt(this->dotProductWith(*this));
#else
  const double scale = this->infnorm();
  assert(scale >= 0.0);

  if( PIPSisZero(scale) )
     return 0.0;

  return scale * sqrt(this->dotProductSelf(1 / scale));
#endif
}

double StochVector::onenorm()
{
  double onenrm = 0.0;

  for( size_t it = 0; it < children.size(); it++ )
     onenrm += children[it]->onenorm();

  if( iAmDistrib == 1 )
  {
     double sum;
     MPI_Allreduce(&onenrm, &sum, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
     onenrm = sum;
  }

  onenrm += vec->onenorm();

  if( vecl )
     onenrm += vecl->onenorm();

  return onenrm;
}


void StochVector::min( double& m, int& index )
{
  double lMin; int lInd;

  if(NULL==parent) {
    vec->min(m,index);
    if( vecl )
    {
       vecl->min(lMin,lInd);
       if( lMin < m )
       {
          m = lMin;
          index = lInd + vec->length();
       }
    }
  } else {
    vec->min(lMin,lInd);
    
    if( vecl )
    {
       double lMinlink;
       int lIndlink;
       vecl->min(lMinlink,lIndlink);
       if( lMinlink < lMin )
       {
          lMin = lMinlink;
          lInd = lIndlink + vec->length();
       }
    }

    if(lMin<m) {
      m = lMin;
      index = lInd + parent->n - this->n;
    }
  }

  for(size_t it=0; it<children.size(); it++) {
    children[it]->min(m,index);
  }

  if(iAmDistrib==1) {
    double minG;
    MPI_Allreduce(&m, &minG, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
    m = minG;
  }  
}


void StochVector::max( double& m, int& index )
{
   double lMax;
   int lInd;

   if( NULL == parent )
   {
      vec->max(m, index);
      if( vecl )
      {
         vecl->max(lMax, lInd);
         if( lMax > m )
         {
            m = lMax;
            index = lInd + vec->length();
         }
      }
   }
   else
   {
      vec->max(lMax, lInd);

      if( vecl )
      {
         double lMaxlink;
         int lIndlink;
         vecl->max(lMaxlink, lIndlink);
         if( lMaxlink > lMax )
         {
            lMax = lMaxlink;
            lInd = lIndlink + vec->length();
         }
      }

      if( lMax > m )
      {
         m = lMax;
         index = lInd + parent->n - this->n;
      }
   }

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->max(m, index);

   if( iAmDistrib == 1 )
   {
      double maxG;
      MPI_Allreduce(&m, &maxG, 1, MPI_DOUBLE, MPI_MAX, mpiComm);
      m = maxG;
   }
}

void StochVector::absminVecUpdate(OoqpVector& absminvec)
{
   StochVector& absminvecStoch = dynamic_cast<StochVector&>(absminvec);
   assert(absminvecStoch.children.size() == children.size());

   if( vecl )
   {
      assert(absminvecStoch.vecl);
      vecl->absminVecUpdate(*(absminvecStoch.vecl));
   }

   vec->absminVecUpdate(*(absminvecStoch.vec));

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->absminVecUpdate(*(absminvecStoch.children[it]));
}

void StochVector::absmaxVecUpdate(OoqpVector& absmaxvec)
{
   StochVector& absmaxvecStoch = dynamic_cast<StochVector&>(absmaxvec);
   assert(absmaxvecStoch.children.size() == children.size());

   if( vecl )
   {
      assert(absmaxvecStoch.vecl);
      vecl->absmaxVecUpdate(*(absmaxvecStoch.vecl));
   }

   vec->absmaxVecUpdate(*(absmaxvecStoch.vec));

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->absmaxVecUpdate(*(absmaxvecStoch.children[it]));
}

void StochVector::absmin(double& m)
{
   double lMin;

   if(NULL==parent) {
     vec->absmin(m);
     if( vecl )
     {
        vecl->absmin(lMin);
        if( lMin < m )
           m = lMin;
     }
   } else {
     vec->absmin(lMin);

     if( vecl )
     {
        double lMinlink;
        vecl->absmin(lMinlink);
        if( lMinlink < lMin )
           lMin = lMinlink;
     }

     if(lMin<m)
       m = lMin;
   }

   for(size_t it=0; it<children.size(); it++) {
     children[it]->absmin(m);
   }

   if(iAmDistrib==1) {
     double minG;
     MPI_Allreduce(&m, &minG, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
     m = minG;
   }
   assert( m >= 0.0 );
}

void StochVector::absminNonZero(double& m, double tolerance)
{
   double lMin;
   bool initialized = false;

   if(NULL==parent) {
     vec->absminNonZero(m, tolerance);
     if( m >= tolerance )
        initialized = true;
     if( vecl )
     {
        vecl->absminNonZero(lMin, tolerance);
        if( lMin >= tolerance && (!initialized || lMin < m) )
           m = lMin;
     }
   } else {
     vec->absminNonZero(lMin, tolerance);
     if( lMin >= tolerance )
        initialized = true;

     if( vecl )
     {
        double lMinlink;
        vecl->absminNonZero(lMinlink, tolerance);
        if( lMinlink >= tolerance && (!initialized || lMinlink < lMin) )
        {
           lMin = lMinlink;
           initialized = true;
        }
     }

     if( initialized && lMin<m )
       m = lMin;
   }

   for(size_t it=0; it<children.size(); it++) {
     children[it]->absminNonZero(m, tolerance);
   }

   if(iAmDistrib==1) {
     double minG;
     if( m == 0.0 )
        m = std::numeric_limits<double>::max();
     MPI_Allreduce(&m, &minG, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
     if( minG < std::numeric_limits<double>::max() )
        m = minG;
     else
        m = 0.0;
   }
   assert( m >= tolerance || m == 0.0 );
}


double StochVector::stepbound(OoqpVector & v_, double maxStep )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);

  double step = this->vec->stepbound(*v.vec, maxStep);

  if( vecl )
  {
     assert(v.vecl);
     double stepl = vecl->stepbound(*v.vecl, maxStep);
     if( stepl < step )
        step = stepl;
  }

  //check tree compatibility
  assert(children.size() == v.children.size());

  for(size_t it=0; it<children.size(); it++)
    step = children[it]->stepbound(*v.children[it], step);
  
  if(iAmDistrib==1) {
    double stepG=0.0;
    MPI_Allreduce(&step, &stepG, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
    step = stepG;
  }
  return step;
}

double StochVector::findBlocking(OoqpVector & wstep_vec, 
			      OoqpVector & u_vec, 
			      OoqpVector & ustep_vec, 
			      double maxStep,
			      double *w_elt, 
			      double *wstep_elt,
			      double *u_elt, 
			      double *ustep_elt,
			      int& first_or_second)
{
  StochVector& w = *this;
  StochVector& u = dynamic_cast<StochVector&>(u_vec);

  StochVector& wstep = dynamic_cast<StochVector&>(wstep_vec);
  StochVector& ustep = dynamic_cast<StochVector&>(ustep_vec);
  const double local_eps = 1e-14;

  double step = maxStep;
  
  // todo only if i am special?
  if( w.vecl )
  {
    assert(wstep.vecl);
    assert(u.vecl);
    assert(ustep.vecl);

    step = w.vecl->findBlocking(*wstep.vecl, *u.vecl, *ustep.vecl, step,
                 w_elt, wstep_elt, u_elt, ustep_elt,
                 first_or_second);
  }

  step = w.vec->findBlocking(*wstep.vec, *u.vec, *ustep.vec, step,
			      w_elt, wstep_elt, u_elt, ustep_elt,
			      first_or_second);

  int nChildren=w.children.size();
  //check tree compatibility
  assert( nChildren             - u.children.size() == 0);
  assert( wstep.children.size() == ustep.children.size() );
  assert( nChildren             - ustep.children.size() == 0);

  for(int it=0; it<nChildren; it++) {
    step = w.children[it]->findBlocking(*wstep.children[it],
			       *u.children[it],
			       *ustep.children[it],
			       step,
			       w_elt,
			       wstep_elt,u_elt,ustep_elt, first_or_second);
  }

  if(iAmDistrib==1) {
    double stepG;
    assert(PIPSisLE(step, 1.0));
    assert(PIPSisLE(0.0, step));

    MPI_Allreduce(&step, &stepG, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
    const bool iHaveMinStep = PIPSisEQ(step, stepG, local_eps);

    //we prefer a AllReduce instead of a bcast, since the step==stepG m
    //may occur for two different processes and a deadlock may occur.
    double buffer[5]; //0-primal val, 1-primal step, 2-dual value, 3-step, 4-1st or 2nd

    int count;
    if( iHaveMinStep ) {
      buffer[0]=*w_elt; buffer[1]=*wstep_elt; 
      buffer[2]=*u_elt; buffer[3]=*ustep_elt;
      buffer[4]=first_or_second;

      count = 1;
    } else {

      count = 0;
      buffer[0]=buffer[1]=buffer[2]=buffer[3]=buffer[4]= -std::numeric_limits<double>::max();
    }

    MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, mpiComm);
    assert(count >= 1);

    // is there more than one process with step==stepG?
    if( count > 1 )
    {
       int myrank;
       int mineqrank;

       MPI_Comm_rank(mpiComm, &myrank);

       if( iHaveMinStep )
          mineqrank = myrank;
       else
          mineqrank = std::numeric_limits<int>::max();

       MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm);

       // step==stepG and not smallest rank?
      if( iHaveMinStep && mineqrank != myrank )
         buffer[0]=buffer[1]=buffer[2]=buffer[3]=buffer[4]= -std::numeric_limits<double>::max();
    }

    double bufferOut[5];
    MPI_Allreduce(buffer, bufferOut, 5, MPI_DOUBLE, MPI_MAX, mpiComm);

    *w_elt = bufferOut[0]; *wstep_elt=bufferOut[1];
    *u_elt = bufferOut[2]; *ustep_elt=bufferOut[3];

    // negative or 0 means no blocking, so set first_or_second to 0.
    if( bufferOut[4] <= 0.5 )
       first_or_second = 0;
    else if( bufferOut[4] <= 1.5 )
       first_or_second = 1;
    else
       first_or_second = 2;

    step=stepG;
  }
  return step;
}

void StochVector::findBlocking_pd(const OoqpVector & wstep_vec,
			      const OoqpVector & u_vec,
			      const OoqpVector & ustep_vec,
			      double& maxStepPri, double& maxStepDual,
			      double& w_elt_p, double& wstep_elt_p, double& u_elt_p, double& ustep_elt_p,
				   double& w_elt_d, double& wstep_elt_d, double& u_elt_d, double& ustep_elt_d,
				   bool& primalBlocking, bool& dualBlocking) const
{
  const StochVector& w = *this;
  const StochVector& u = dynamic_cast<const StochVector&>(u_vec);
  const double local_eps = 1e-14;

  const StochVector& wstep = dynamic_cast<const StochVector&>(wstep_vec);
  const StochVector& ustep = dynamic_cast<const StochVector&>(ustep_vec);

  // todo only if i am special?
  if( w.vecl )
  {
    assert(wstep.vecl);
    assert(u.vecl);
    assert(ustep.vecl);

    w.vecl->findBlocking_pd(*wstep.vecl, *u.vecl, *ustep.vecl, maxStepPri, maxStepDual,
                 w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p,
				     w_elt_d, wstep_elt_d, u_elt_d, ustep_elt_d,
                 primalBlocking, dualBlocking);
  }

  w.vec->findBlocking_pd(*wstep.vec, *u.vec, *ustep.vec, maxStepPri, maxStepDual,
		  	  	  w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p,
				  w_elt_d, wstep_elt_d, u_elt_d, ustep_elt_d,
				  primalBlocking, dualBlocking);

  int nChildren=w.children.size();
  //check tree compatibility
  assert( nChildren             - u.children.size() == 0);
  assert( wstep.children.size() == ustep.children.size() );
  assert( nChildren             - ustep.children.size() == 0);

  for(int it=0; it<nChildren; it++) {
    w.children[it]->findBlocking_pd(*wstep.children[it],
			       *u.children[it],
			       *ustep.children[it],
			       maxStepPri, maxStepDual,
				    w_elt_p, wstep_elt_p, u_elt_p, ustep_elt_p,
				    w_elt_d, wstep_elt_d, u_elt_d, ustep_elt_d,
				    primalBlocking, dualBlocking);
   }

   if( iAmDistrib == 1 )
   {
      double maxStepGlobalPri, maxStepGlobalDual;
      assert(PIPSisLE(maxStepPri, 1.0) && PIPSisLE(maxStepDual, 1.0));
      assert(PIPSisLE(0.0, maxStepPri) && PIPSisLE(0.0, maxStepDual));

      MPI_Allreduce(&maxStepPri, &maxStepGlobalPri, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
      MPI_Allreduce(&maxStepDual, &maxStepGlobalDual, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
      const bool iHaveMinStepPri = PIPSisEQ(maxStepPri, maxStepGlobalPri, local_eps);
      const bool iHaveMinStepDual = PIPSisEQ(maxStepDual, maxStepGlobalDual, local_eps);

      //we prefer a AllReduce instead of a bcast, since the step==stepG
      //may occur for two different processes and a deadlock may occur.
      double buffer[10];
      int count[2];
      //values for computation of the primal steplength:
      //0-primal val, 1-primal step, 2-dual value, 3-dual step, 4-primalBlocking
      if( iHaveMinStepPri )
      {
         buffer[0] = w_elt_p;
         buffer[1] = wstep_elt_p;
         buffer[2] = u_elt_p;
         buffer[3] = ustep_elt_p;
         buffer[4] = primalBlocking ? 1.0 : 0.0;

         count[0] = 1;
      }
      else
      {
         buffer[0] = buffer[1] = buffer[2] = buffer[3] = buffer[4] =
               -std::numeric_limits<double>::max();
         count[0] = 0;
      }

      //values for computation of the dual steplength:
      //5-primal val, 6-primal step, 7-dual value, 8-dual step, 9-dualBlocking
      if( iHaveMinStepDual )
      {
         buffer[5] = w_elt_d;
         buffer[6] = wstep_elt_d;
         buffer[7] = u_elt_d;
         buffer[8] = ustep_elt_d;
         buffer[9] = dualBlocking ? 1.0 : 0.0;

         count[1] = 1;
      }
      else
      {
         buffer[5] = buffer[6] = buffer[7] = buffer[8] = buffer[9] =
               -std::numeric_limits<double>::max();
         count[1] = 0;
      }

      MPI_Allreduce(MPI_IN_PLACE, count, 2, MPI_INT, MPI_SUM, mpiComm);

      assert(count[0] >= 1 && count[1] >= 1);

      int myrank;
      MPI_Comm_rank(mpiComm, &myrank);

      // is there more than one process with maxStepPri==stepG?
      if( count[0] > 1 )
      {
         int mineqrank;

         if( iHaveMinStepPri )
            mineqrank = myrank;
         else
            mineqrank = std::numeric_limits<int>::max();

         MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm);

         // step==stepG and not smallest rank?
         if( iHaveMinStepPri && mineqrank != myrank )
            buffer[0] = buffer[1] = buffer[2] = buffer[3]=buffer[4]= -std::numeric_limits<double>::max();
      }

      // is there more than one process with maxStepDual==stepF?
      if( count[1] > 1 )
      {
         int mineqrank;

         if( iHaveMinStepDual )
            mineqrank = myrank;
         else
            mineqrank = std::numeric_limits<int>::max();

         MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm);

         // stepDual==stepF and not smallest rank?
         if( iHaveMinStepDual && mineqrank != myrank )
            buffer[5]=buffer[6]=buffer[7]=buffer[8]=buffer[9]= -std::numeric_limits<double>::max();
      }

      double bufferOut[10];
      MPI_Allreduce(buffer, bufferOut, 10, MPI_DOUBLE, MPI_MAX, mpiComm);

      w_elt_p = bufferOut[0];
      wstep_elt_p = bufferOut[1];
      u_elt_p = bufferOut[2];
      ustep_elt_p = bufferOut[3];

      w_elt_d = bufferOut[5];
      wstep_elt_d = bufferOut[6];
      u_elt_d = bufferOut[7];
      ustep_elt_d = bufferOut[8];

      primalBlocking = bufferOut[4] <= 0.5 ? false : true;
      maxStepPri = maxStepGlobalPri;

      dualBlocking = bufferOut[9] <= 0.5 ? false : true;
      maxStepDual = maxStepGlobalDual;
   }
}

void StochVector::componentMult( OoqpVector& v_ )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);
  assert(v.children.size() == children.size());

  vec->componentMult(*v.vec);
  if( vecl ) vecl->componentMult(*v.vecl);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->componentMult(*v.children[it]);
}

void StochVector::componentDiv ( OoqpVector& v_ )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);
  assert(v.children.size() == children.size());

  vec->componentDiv(*v.vec);
  if( vecl ) vecl->componentDiv(*v.vecl);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->componentDiv(*v.children[it]);
}

void StochVector::scalarMult( double num )
{
  vec->scalarMult(num);
  if( vecl ) vecl->scalarMult(num);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->scalarMult(num);
}

void StochVector::writeToStreamAll( ostream& out ) const
{
   int rank;
   MPI_Comm_rank(mpiComm, &rank);
   int world_size;
   MPI_Comm_size(mpiComm, &world_size);
   MPI_Status status;
   int l;
   stringstream sout;

   if( rank == 0)
   {
      sout << "----" << endl;
      vec->writeToStreamAllStringStream(sout);

      for( size_t it = 0; it < children.size(); it++ )
         children[it]->writeToStreamAllChild(sout);

      out << sout.str();
      sout.str(std::string());

      for( int p = 1; p < world_size; p++ )
      {
         MPI_Probe(p, p, mpiComm, &status);
         MPI_Get_count(&status, MPI_CHAR, &l);
         char *buf = new char[l];
         MPI_Recv(buf, l, MPI_CHAR, p, p, mpiComm, &status);
         string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }
      if( vecl )
      {
         sout << "---" << endl;
         vecl->writeToStreamAllStringStream(sout);
      }
      sout << "----" << endl;
      out << sout.str();
   }
   else if( iAmDistrib==1 )
   { // rank != 0
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->writeToStreamAllChild(sout);

      std::string str = sout.str();
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, rank, mpiComm);

   }

   if( iAmDistrib==1 )
      MPI_Barrier(mpiComm);
}

void StochVector::writeToStreamAllChild( stringstream& sout ) const
{
   sout << "--" << endl;
   vec->writeToStreamAllStringStream(sout);

   for( size_t it = 0; it < children.size(); it++ ){
      sout << "-- " << endl;
      children[it]->writeToStreamAllChild(sout);
   }

   if( vecl )
   {
      sout << "---" << endl;
      vecl->writeToStreamAllStringStream(sout);
   }
}

void StochVector::writeToStream( ostream& out ) const
{
  out << "---" << endl;
  vec->writeToStream(out);
  if( vecl ) vecl->writeToStream(out);
  out << "~~~" << endl;
  //for(size_t it=0; it<children.size(); it++) 
  //  children[it]->writeToStream(out);
}

void StochVector::writefToStream( ostream& out,
				  const char format[] ) const
{
  vec->writefToStream(out, format);
  if( vecl ) vecl->writefToStream(out, format);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->writefToStream(out, format);
}

void StochVector::writeMPSformatRhs(ostream& out, int rowType, OoqpVector* irhs) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   string rt;
   if( rowType == 0 )
      rt = "E";
   else if( rowType == 1 )
      rt = "L";
   else if( rowType == 2 )
      rt = "G";
   else
      assert(0);

   StochVector* ic = NULL;
   if( irhs )
      ic = dynamic_cast<StochVector*>(irhs);

   if( myRank==0 )
   {
      string rowNameStub = " B row_";
      rowNameStub+= rt;
      rowNameStub+="_R_";
      if( irhs && ic )
         vec->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<SimpleVector*>(ic->vec));
      else
         vec->writeMPSformatOnlyRhs( out, rowNameStub, NULL);
      if(vecl)
      {
         rowNameStub = " B row_";
         rowNameStub+= rt;
         rowNameStub+="_L_";
         if( irhs )
            vecl->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<SimpleVector*>(ic->vecl));
         else
            vecl->writeMPSformatOnlyRhs( out, rowNameStub, NULL);
      }
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm;
      sstm << " B row_" << rt << "_" << it << "_";
      string rowNameStub = sstm.str();
      if( irhs )
         children[it]->vec->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<SimpleVector*>(ic->children[it]->vec));
      else
         children[it]->vec->writeMPSformatOnlyRhs( out, rowNameStub, NULL);
   }
}

void StochVector::writeMPSformatBounds(ostream& out, OoqpVector* ix, bool upperBound) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);

   StochVector* ixStoch = dynamic_cast<StochVector*>(ix);

   if( myRank==0 )
   {
      string varNameStub = "var_L_";
      vec->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->vec), upperBound);
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm2;
      sstm2 << "var_" << it << "_";
      string varNameStub = sstm2.str();
      children[it]->vec->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->children[it]->vec), upperBound);
   }
}

/** this += alpha * x */
void StochVector::axpy  ( double alpha, OoqpVector& x_ )
{
  StochVector& x = dynamic_cast<StochVector&>(x_);
  assert(x.children.size() == children.size());

  if( alpha == 0.0)
     return;

  vec->axpy(alpha, *x.vec);

  if( vecl )
  {
     assert(x.vecl);
     vecl->axpy(alpha, *x.vecl);
  }

  for(size_t it=0; it<children.size(); it++)
    children[it]->axpy(alpha, *x.children[it]);
}

/** this += alpha * x * z */
void StochVector::axzpy ( double alpha, OoqpVector& x_, OoqpVector& z_ )
{
  StochVector& x = dynamic_cast<StochVector&>(x_);
  StochVector& z = dynamic_cast<StochVector&>(z_);
  assert(x.children.size() == children.size());
  assert(z.children.size() == children.size());

  vec->axzpy(alpha, *x.vec, *z.vec);

  if( vecl )
  {
     assert(x.vecl);
     assert(z.vecl);
     vecl->axzpy(alpha, *x.vecl, *z.vecl);
  }

  for(size_t it=0; it<children.size(); it++)
    children[it]->axzpy(alpha, *x.children[it], *z.children[it]);
}

/** this += alpha * x / z */
void StochVector::axdzpy( double alpha, OoqpVector& x_, OoqpVector& z_ )
{
  StochVector& x = dynamic_cast<StochVector&>(x_);
  StochVector& z = dynamic_cast<StochVector&>(z_);
  assert(x.children.size() == children.size());
  assert(z.children.size() == children.size());

  vec->axdzpy(alpha, *x.vec, *z.vec);

  if( vecl )
  {
    assert(x.vecl);
    assert(z.vecl);
    vecl->axdzpy(alpha, *x.vecl, *z.vecl);
  }

  for(size_t it=0; it<children.size(); it++)
    children[it]->axdzpy(alpha, *x.children[it], *z.children[it]);
}


void StochVector::addConstant( double c )
{
  vec->addConstant(c);

  if( vecl ) vecl->addConstant(c);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->addConstant(c);
}


void StochVector::gondzioProjection( double rmin, double rmax )
{
  vec->gondzioProjection( rmin, rmax );


  if( vecl ) vecl->gondzioProjection( rmin, rmax );

  for(size_t it=0; it<children.size(); it++)
    children[it]->gondzioProjection( rmin, rmax );
}

double StochVector::dotProductWith( OoqpVector& v_ )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);
  assert(v.children.size() == children.size());

  double dotProd=0.0;

  for(size_t it=0; it<children.size(); it++) 
    dotProd += children[it]->dotProductWith(*v.children[it]);

  assert(!vecl || v.vecl);

  if(iAmDistrib==1) {
    double dotProdG = 0.0;

    MPI_Allreduce(&dotProd, &dotProdG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);

    dotProd = dotProdG;
  }

  dotProd += vec->dotProductWith(*v.vec);

  if( vecl )
    dotProd += vecl->dotProductWith(*v.vecl);

  return dotProd;
}

double StochVector::dotProductSelf(double scaleFactor)
{
  double dotSelf = 0.0;

  for(size_t it=0; it<children.size(); it++)
     dotSelf += children[it]->dotProductSelf(scaleFactor);

  if(iAmDistrib==1) {
    double dotSelfG = 0.0;

    MPI_Allreduce(&dotSelf, &dotSelfG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);

    dotSelf = dotSelfG;
  }

  dotSelf += vec->dotProductSelf(scaleFactor);

  if( vecl )
     dotSelf += vecl->dotProductSelf(scaleFactor);

  return dotSelf;
}

/** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
 */
double StochVector::shiftedDotProductWith( double alpha, OoqpVector& mystep_,
					OoqpVector& yvec_,
					double beta,  OoqpVector& ystep_ )
{
  StochVector& mystep = dynamic_cast<StochVector&>(mystep_);
  StochVector& yvec   = dynamic_cast<StochVector&>(yvec_);
  StochVector& ystep  = dynamic_cast<StochVector&>(ystep_);


  double dotProd = 0.0;
  for(size_t it=0; it<children.size(); it++) 
    dotProd += children[it]->shiftedDotProductWith(alpha, *mystep.children[it], 
						   *yvec.children[it],
						   beta, *ystep.children[it]);
  if(iAmDistrib) {
    double dotProdG=0.0;
    MPI_Allreduce(&dotProd, &dotProdG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
    dotProd = dotProdG;
  }

  dotProd += vec->shiftedDotProductWith(alpha, *mystep.vec, 
					*yvec.vec,
					beta, *ystep.vec);

  if( vecl )
  {
	 assert(mystep.vecl);
	 assert(yvec.vecl);
	 assert(ystep.vecl);
    dotProd += vecl->shiftedDotProductWith(alpha, *mystep.vecl,  *yvec.vecl, beta, *ystep.vecl);
  }

  return dotProd;
}

void StochVector::negate()
{
  vec->negate();
  if( vecl ) vecl->negate();

  for(size_t it=0; it<children.size(); it++) 
    children[it]->negate();
}

void StochVector::invert()
{
  vec->invert();

  if( vecl ) vecl->invert();

  for(size_t it=0; it<children.size(); it++)
    children[it]->invert();
}

void StochVector::invertSave(double zeroReplacementVal)
{
  vec->invertSave(zeroReplacementVal);

  if( vecl ) vecl->invertSave(zeroReplacementVal);

  for( size_t it = 0; it < children.size(); it++ )
    children[it]->invertSave(zeroReplacementVal);
}

void StochVector::applySqrt()
{
   vec->applySqrt();

   if(vecl) vecl->applySqrt();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->applySqrt();
}

void StochVector::roundToPow2()
{
  vec->roundToPow2();

  if( vecl ) vecl->roundToPow2();

  for( size_t it = 0; it < children.size(); it++ )
     children[it]->roundToPow2();
}

int StochVector::allPositive()
{
  //!parallel
  int allPos = vec->allPositive() && ((vecl != NULL) ? vecl->allPositive() : 1);
  if (!allPos) return 0;

  for(size_t it=0; it<children.size() && allPos; it++) 
    allPos = children[it]->allPositive();

  return allPos;
}

int StochVector::matchesNonZeroPattern( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  int match = vec->matchesNonZeroPattern(*select.vec);

  if( vecl )
  {
     assert(select.vecl);
     match = match && vecl->matchesNonZeroPattern(*select.vecl);
  }

  if(!match) return 0;

  for(size_t it=0; it<children.size() && match; it++) 
    match = children[it]->matchesNonZeroPattern(*select.children[it]);

  return match;
}

void StochVector::selectNonZeros( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  vec->selectNonZeros(*select.vec);

  if( vecl )
  {
     assert(select.vecl);
     vecl->selectNonZeros(*select.vecl);
  }

  for(size_t it=0; it<children.size(); it++) 
    children[it]->selectNonZeros(*select.children[it]);
}

long long StochVector::numberOfNonzeros()
{
  //!opt - store the number of nnz to avoid communication
  long long nnz = 0;

  for(size_t it=0; it<children.size(); it++) 
    nnz += children[it]->numberOfNonzeros();

  if(iAmDistrib) {
    long long nnzG = 0;
    MPI_Allreduce(&nnz, &nnzG, 1, MPI_LONG_LONG, MPI_SUM, mpiComm);
    nnz = nnzG;
  }
  nnz += vec->numberOfNonzeros();

  if( vecl ) nnz += vecl->numberOfNonzeros();

  return nnz;
}
void StochVector::addSomeConstants( double c, OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  vec->addSomeConstants(c, *select.vec);

  if( vecl )
  {
     assert(select.vecl);
     vecl->addSomeConstants(c, *select.vecl);
  }

  for(size_t it=0; it<children.size(); it++) 
    children[it]->addSomeConstants(c, *select.children[it]);
}

void StochVector::writefSomeToStream( ostream& out,
			 const char format[],
			 OoqpVector& select_ ) const
{
  assert( "Have not been yet implemented" && 0 );
}

void StochVector::axdzpy( double alpha, OoqpVector& x_,
		       OoqpVector& z_, OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  StochVector&      x = dynamic_cast<StochVector&>(x_);
  StochVector&      z = dynamic_cast<StochVector&>(z_);

  assert(children.size() == select.children.size());
  assert(children.size() == x.     children.size());
  assert(children.size() == z.     children.size());

  vec->axdzpy(alpha, *x.vec, *z.vec, *select.vec);

  if( vecl )
  {
     assert(x.vecl);
     assert(z.vecl);
     assert(select.vecl);
     vecl->axdzpy(alpha, *x.vecl, *z.vecl, *select.vecl);
  }

  for(size_t it=0; it<children.size(); it++)
    children[it]->axdzpy(alpha, *x.children[it], *z.children[it], *select.children[it]);
}

int StochVector::somePositive( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  //!parallel stuff needed

  int somePos = vec->somePositive(*select.vec);

  if( vecl )
  {
    assert(select.vecl);
    somePos = somePos && vecl->somePositive(*select.vecl);
  }

  for(size_t it=0; it<children.size() && somePos; it++)
    somePos = children[it]->somePositive(*select.children[it]);
  
  return somePos;
}

void StochVector::divideSome( OoqpVector& div_, OoqpVector& select_ )
{
  StochVector& div    = dynamic_cast<StochVector&>(div_);
  StochVector& select = dynamic_cast<StochVector&>(select_);

  assert(children.size() == div.   children.size());
  assert(children.size() == select.children.size());

  vec->divideSome(*div.vec, *select.vec);

  if( vecl )
  {
     assert(div.vecl);
     assert(select.vecl);
     vecl->divideSome(*div.vecl, *select.vecl);
  }

  for(size_t it=0; it<children.size(); it++)
    children[it]->divideSome(*div.children[it], *select.children[it]);
}

void StochVector::copyIntoArray( double v[] ) const
{
  assert( "Not supported" && 0 );
}

void StochVector::copyFromArray( double v[] )
{
  assert( "Not supported" && 0 );
}

void StochVector::copyFromArray( char v[] )
{
  assert( "Not supported" && 0 );
}

void StochVector::removeEntries( const OoqpVector& select )
{
   const StochVector& selectStoch = dynamic_cast<const StochVector&>(select);

   assert(children.size() == selectStoch.children.size());

   vec->removeEntries(*selectStoch.vec);
   n = vec->n;

   if( vecl )
   {
      assert(selectStoch.vecl);
      vecl->removeEntries(*selectStoch.vecl);
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->removeEntries(*selectStoch.children[it]);
      n += children[it]->n;
   }
}

void StochVector::permuteVec0Entries(const std::vector<unsigned int>& permvec)
{
   dynamic_cast<SimpleVector*>(vec)->permuteEntries(permvec);
}

void StochVector::permuteLinkingEntries(const std::vector<unsigned int>& permvec)
{
   if( vecl )
      dynamic_cast<SimpleVector*>(vecl)->permuteEntries(permvec);
}

std::vector<double> StochVector::gatherStochVector() const
{
   const SimpleVector& firstvec = dynamic_cast<const SimpleVector&>(*vec);
   const size_t nChildren = children.size();

   int myrank;
   MPI_Comm_rank(mpiComm, &myrank);
   int mysize;
   MPI_Comm_size(mpiComm, &mysize);

   std::vector<double> gatheredVecLocal;

   for( size_t i = 0; i < nChildren; ++i )
   {
      const SimpleVector& vec = dynamic_cast<const SimpleVector&>(*children[i]->vec);

      if( vec.length() > 0 )
         gatheredVecLocal.insert(gatheredVecLocal.end(), &vec[0], &vec[0] + vec.length());
   }

   size_t solLength = firstvec.length();

   // final vector
   std::vector<double> gatheredVec(0);

   if( mysize > 0 )
   {
      // get all lengths
      std::vector<int> recvcounts(mysize);
      std::vector<int> recvoffsets(mysize);

      int mylength = int(gatheredVecLocal.size());

      MPI_Allgather(&mylength, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, mpiComm);

      // all-gather local components
      recvoffsets[0] = 0;
      for( size_t i = 1; i < size_t(mysize); ++i )
         recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

      if( myrank == 0 )
      {
         solLength += recvoffsets[mysize - 1] + recvcounts[mysize - 1];
         gatheredVec = std::vector<double>(solLength);

         MPI_Gatherv(&gatheredVecLocal[0], mylength, MPI_DOUBLE,
               &gatheredVec[0] + firstvec.length(), &recvcounts[0],
               &recvoffsets[0], MPI_DOUBLE, 0, mpiComm);
      }
      else
      {
         MPI_Gatherv(&gatheredVecLocal[0], mylength, MPI_DOUBLE, 0,
               &recvcounts[0], &recvoffsets[0], MPI_DOUBLE, 0, mpiComm);
      }
   }
   else
   {
      solLength += gatheredVecLocal.size();

      gatheredVec = std::vector<double>(solLength);

      std::copy(gatheredVecLocal.begin(), gatheredVecLocal.end(), gatheredVec.begin() + firstvec.length());
   }

   if( myrank == 0 )
   {
      std::copy(&firstvec[0], &firstvec[0] + firstvec.length(), &gatheredVec[0]);

      if( vecl && vecl->length() > 0 )
      {
         const SimpleVector& linkvec = dynamic_cast<const SimpleVector&>(*vecl);
         gatheredVec.insert(gatheredVec.end(), &linkvec[0], &linkvec[0] + linkvec.length());
      }
   }


   return gatheredVec;
}
