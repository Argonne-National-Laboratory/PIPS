#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "VectorUtilities.h"
#include "StochVector.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <math.h>

template<typename T>
StochVectorBase<T>::StochVectorBase(int n_, MPI_Comm mpiComm_, int isDistributed/*=-1*/)
  : OoqpVectorBase<T>(n_), vecl(NULL), parent(NULL), mpiComm(mpiComm_),
    iAmDistrib(isDistributed)
{
  vec = new SimpleVectorBase<T>(n_);

  if( -1 == iAmDistrib && MPI_COMM_NULL != mpiComm) {
    int size;
    MPI_Comm_size(mpiComm, &size);
    iAmDistrib = (size == 1) ? 0 : 1;
  }
  vecl = NULL;
}

template<typename T>
StochVectorBase<T>::StochVectorBase(int n_, int nl_, MPI_Comm mpiComm_, int isDistributed)
  : OoqpVectorBase<T>(n_), parent(NULL), mpiComm(mpiComm_),
    iAmDistrib(isDistributed)
{
  vec = new SimpleVectorBase<T>(n_);

  if( nl_ >= 0 )
    vecl = new SimpleVectorBase<T>(nl_);
  else
	 vecl = NULL;

  if( -1 == iAmDistrib && MPI_COMM_NULL != mpiComm) {
    int size;
    MPI_Comm_size(mpiComm, &size);
    iAmDistrib = (size == 1) ? 0 : 1;
  }

}

template<typename T>
void StochVectorBase<T>::AddChild(StochVectorBase<T>* child)
{
  child->parent = this;
  children.push_back(child);
  this->n += child->n;
}

template<typename T>
void StochVectorBase<T>::AddChild(OoqpVectorBase<T>* child_)
{
  StochVectorBase<T>* child = reinterpret_cast<StochVectorBase<T>*>(child_);
  AddChild(child);
}

template<typename T>
StochVectorBase<T>::~StochVectorBase()
{
  for (size_t it = 0; it < children.size(); it++)
    delete children[it];

  if( vec )
    delete vec;

  if( vecl )
	 delete vecl;
}

template<typename T>
OoqpVectorBase<T>* StochVectorBase<T>::dataClone() const
{
  assert(!vecl);
  OoqpVectorBase<T>* clone = new SimpleVectorBase<T>(vec->length());
  return clone;
}

template<typename T>
OoqpVectorBase<T>* StochVectorBase<T>::dataCloneLinkCons() const
{
  assert(vecl);
  OoqpVectorBase<T>* clone = new SimpleVectorBase<T>(vecl->length());
  return clone;
}

template<typename T>
OoqpVectorBase<T>* StochVectorBase<T>::clone() const
{
  StochVectorBase<T>* clone;
  if( vecl )
    clone = new StochVectorBase<T>(vec->length(), vecl->length(), mpiComm, -1);
  else
	 clone = new StochVectorBase<T>(vec->length(), mpiComm);

  for(size_t it = 0; it < children.size(); it++) {
    clone->AddChild(children[it]->clone());
  }
  return clone;
}

template<typename T>
OoqpVectorBase<T>* StochVectorBase<T>::cloneFull() const
{
   StochVectorBase<T>* clone = new StochVectorBase<T>(vec->length(), (vecl != NULL) ? vecl->length() : -1, mpiComm, -1);

   clone->vec->copyFrom(*vec);

   if( vecl )
      clone->vecl->copyFrom(*vecl);

   for( size_t it = 0; it < children.size(); it++ )
      clone->AddChild(children[it]->cloneFull());

   return clone;
}


template<typename T>
void StochVectorBase<T>::jointCopyFrom(const StochVectorBase<T>& v1, const StochVectorBase<T>& v2, const StochVectorBase<T>& v3)
{
  SimpleVectorBase<T>& sv  = dynamic_cast<SimpleVectorBase<T>&>(*this->vec);
  SimpleVectorBase<T>& sv1 = dynamic_cast<SimpleVectorBase<T>&>(*v1.vec);
  SimpleVectorBase<T>& sv2 = dynamic_cast<SimpleVectorBase<T>&>(*v2.vec);
  SimpleVectorBase<T>& sv3 = dynamic_cast<SimpleVectorBase<T>&>(*v3.vec);

  int n1 = sv1.length();
  int n2 = sv2.length();
  int n3 = sv3.length();

  assert( n1 + n2 + n3 == sv.length() );

  if( n1 > 0 )
    memcpy(&sv[0], &sv1[0], n1 * sizeof(T));

  if( n2 > 0 )
    memcpy(&sv[n1], &sv2[0], n2 * sizeof(T));

  if( n3 > 0 )
    memcpy(&sv[n1 + n2], &sv3[0], n3 * sizeof(T));

  for(size_t it = 0; it < children.size(); it++) {
    children[it]->jointCopyFrom(*v1.children[it],
				*v2.children[it],
				*v3.children[it]);
  }

}

template<typename T>
void StochVectorBase<T>::jointCopyFromLinkCons(const StochVectorBase<T>& vx, const StochVectorBase<T>& vy, const StochVectorBase<T>& vz)
{
  SimpleVectorBase<T>& sv  = dynamic_cast<SimpleVectorBase<T>&>(*this->vec);
  SimpleVectorBase<T>& svx = dynamic_cast<SimpleVectorBase<T>&>(*vx.vec);
  SimpleVectorBase<T>& svy = dynamic_cast<SimpleVectorBase<T>&>(*vy.vec);
  SimpleVectorBase<T>& svz = dynamic_cast<SimpleVectorBase<T>&>(*vz.vec);

  int n1 = svx.length();
  int n2 = svy.length();
  int n3 = svz.length();
  int n4 = 0;
  int n5 = 0;

  assert(n1+n2+n3 <= sv.length());
  assert( sizeof(T) == sizeof(sv[0]) );

  if( n1 > 0 )
    memcpy(&sv[0], &svx[0], n1 * sizeof(T));

  if( n2 > 0 )
    memcpy(&sv[n1], &svy[0], n2 * sizeof(T));

  if( n3 > 0 )
    memcpy(&sv[n1 + n2], &svz[0], n3 * sizeof(T));

  if( vy.vecl )
  {
    SimpleVectorBase<T>& svyl = dynamic_cast<SimpleVectorBase<T>&>(*vy.vecl);
    n4 = svyl.length();
    assert( n4 >= 0 );

    if( n4 > 0 )
      memcpy(&sv[n1 + n2 + n3], &svyl[0], n4 * sizeof(T));
  }

  if( vz.vecl )
  {
    SimpleVectorBase<T>& svzl = dynamic_cast<SimpleVectorBase<T>&>(*vz.vecl);
    n5 = svzl.length();
    assert( n5 >= 0 );

    if( n5 > 0 )
      memcpy(&sv[n1 + n2 + n3 + n4], &svzl[0], n5 * sizeof(T));
  }

  assert(n1+n2+n3+n4+n5 == sv.length());

  for(size_t it = 0; it < children.size(); it++) {
    children[it]->jointCopyFromLinkCons(*vx.children[it],
				*vy.children[it],
				*vz.children[it]);
  }
}


template<typename T>
void StochVectorBase<T>::jointCopyTo(StochVectorBase<T>& v1, StochVectorBase<T>& v2, StochVectorBase<T>& v3) const
{
  const SimpleVectorBase<T>& sv  = dynamic_cast<const SimpleVectorBase<T>&>(*this->vec);
  SimpleVectorBase<T>& sv1 = dynamic_cast<SimpleVectorBase<T>&>(*v1.vec);
  SimpleVectorBase<T>& sv2 = dynamic_cast<SimpleVectorBase<T>&>(*v2.vec);
  SimpleVectorBase<T>& sv3 = dynamic_cast<SimpleVectorBase<T>&>(*v3.vec);

  int n1 = sv1.length();
  int n2 = sv2.length();
  int n3 = sv3.length();

  assert(n1+n2+n3 == sv.length());

  if(n1 > 0)
    memcpy(&sv1[0], &sv[0], n1 * sizeof(T));

  if(n2 > 0)
    memcpy(&sv2[0], &sv[n1], n2 * sizeof(T));

  if(n3 > 0)
    memcpy(&sv3[0], &sv[n1+n2], n3 * sizeof(T));


  for(size_t it = 0; it < children.size(); it++) {
    children[it]->jointCopyTo(*v1.children[it],
			      *v2.children[it],
			      *v3.children[it]);
  }
}

template<typename T>
void StochVectorBase<T>::jointCopyToLinkCons(StochVectorBase<T>& vx, StochVectorBase<T>& vy, StochVectorBase<T>& vz) const
{
  const SimpleVectorBase<T>& sv  = dynamic_cast<const SimpleVectorBase<T>&>(*this->vec);
  SimpleVectorBase<T>& svx = dynamic_cast<SimpleVectorBase<T>&>(*vx.vec);
  SimpleVectorBase<T>& svy = dynamic_cast<SimpleVectorBase<T>&>(*vy.vec);
  SimpleVectorBase<T>& svz = dynamic_cast<SimpleVectorBase<T>&>(*vz.vec);

  int n1 = svx.length();
  int n2 = svy.length();
  int n3 = svz.length();
  int n4 = 0;
  int n5 = 0;

  assert( n1 + n2 + n3 <= sv.length() );
  assert( sizeof(T) == sizeof(sv[0]) );

  if(n1 > 0)
    memcpy(&svx[0], &sv[0], n1 * sizeof(T));

  if(n2 > 0)
    memcpy(&svy[0], &sv[n1], n2 * sizeof(T));

  if(n3 > 0)
    memcpy(&svz[0], &sv[n1 + n2], n3 * sizeof(T));

  if( vy.vecl )
  {
     SimpleVectorBase<T>& svyl = dynamic_cast<SimpleVectorBase<T>&>(*vy.vecl);
     n4 = svyl.length();
     assert(n4 >= 0);

     if( n4 > 0 )
       memcpy(&svyl[0], &sv[n1 + n2 + n3], n4 * sizeof(T));
  }

  if( vz.vecl )
  {
     SimpleVectorBase<T>& svzl = dynamic_cast<SimpleVectorBase<T>&>(*vz.vecl);
     n5 = svzl.length();
     assert(n5>= 0);

     if( n5 > 0 )
       memcpy(&svzl[0], &sv[n1 + n2 + n3 + n4], n5 * sizeof(T));
  }

  assert(n1 + n2 + n3 + n4 + n5 == sv.length());

  for(size_t it = 0; it < children.size(); it++) {
    children[it]->jointCopyToLinkCons(*vx.children[it],
			      *vy.children[it],
			      *vz.children[it]);
  }
}


template<typename T>
bool StochVectorBase<T>::isKindOf( int kind ) const
{
  return kind==kStochVector;
}

template<typename T>
void StochVectorBase<T>::scale( T alpha )
{
  vec->scale(alpha);

  if( vecl ) vecl->scale(alpha);

  for(size_t it = 0; it < children.size(); it++)
    children[it]->scale(alpha);
}

template<typename T>
bool StochVectorBase<T>::isZero() const
{
	bool is_zero = true;

	is_zero = (is_zero && dynamic_cast<SimpleVectorBase<T>&>(*vec).isZero());

	if(vecl)
		is_zero = (is_zero && dynamic_cast<SimpleVectorBase<T>&>(*vecl).isZero());

	for( size_t node = 0; node < children.size(); ++node )
		is_zero = (is_zero && children[node]->isZero());

	return is_zero;
}

template<typename T>
void StochVectorBase<T>::setToZero()
{
  vec->setToZero();

  if( vecl ) vecl->setToZero();

  for(size_t it = 0; it<children.size(); it++)
    children[it]->setToZero();
}

template<typename T>
void StochVectorBase<T>::setToConstant( T c)
{
  vec->setToConstant(c);

  if( vecl ) vecl->setToConstant(c);

  for(size_t it = 0; it<children.size(); it++)
    children[it]->setToConstant(c);
}

template<typename T>
void StochVectorBase<T>::randomize( T alpha, T beta, T *ix )
{
  assert( "Not implemented" && 0 );
}


template<typename T>
void StochVectorBase<T>::copyFrom( const OoqpVectorBase<T>& v_ )
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

  this->vec->copyFrom(*v.vec);

  if( this->vecl )
  {
     assert(v.vecl);
     this->vecl->copyFrom(*v.vecl);
  }

  //assert tree compatibility
  assert(children.size() == v.children.size());

  for(size_t it  =0; it < children.size(); it++)
    children[it]->copyFrom(*v.children[it]);
}

template<typename T>
void StochVectorBase<T>::copyFromAbs(const OoqpVectorBase<T>& v_ )
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

  this->vec->copyFromAbs(*v.vec);

  if( this->vecl )
  {
     assert(v.vecl);
     this->vecl->copyFromAbs(*v.vecl);
  }

  //assert tree compatibility
  assert(children.size() == v.children.size());

  for(size_t it = 0; it < children.size(); it++)
    children[it]->copyFromAbs(*v.children[it]);
}

template<typename T>
T StochVectorBase<T>::infnorm() const
{
  T infnrm;

  infnrm = 0.0;

  for(size_t it = 0; it < children.size(); it++)
    infnrm = std::max(infnrm, children[it]->infnorm());

  if(iAmDistrib) {
    T infnrmG = PIPS_MPIgetMax(infnrm, mpiComm);
    // MPI_Allreduce(&infnrm, &infnrmG, 1, MPI_DOUBLE, MPI_MAX, mpiComm); // not working properly in templated version
    infnrm = infnrmG;
  }

  infnrm = std::max(vec->infnorm(), infnrm);

  if( vecl ) infnrm = std::max(vecl->infnorm(), infnrm);

  return infnrm;
}

template<typename T>
double StochVectorBase<T>::twonorm() const
{
#if 0
  return sqrt(this->dotProductWith(*this));
#else
  const T scale = this->infnorm();
  assert(scale >= 0.0);

  if( PIPSisZero(scale) )
     return 0.0;

  return scale * sqrt(this->dotProductSelf(1 / scale));
#endif
}

template<typename T>
T StochVectorBase<T>::onenorm() const
{
  T onenrm = 0.0;

  for( size_t it = 0; it < children.size(); it++ )
     onenrm += children[it]->onenorm();

  if( iAmDistrib == 1 )
  {
     T sum = PIPS_MPIgetSum(onenrm, mpiComm);
     // MPI_Allreduce(&onenrm, &sum, 1, MPI_DOUBLE, MPI_SUM, mpiComm); // not working properly in templated version
     onenrm = sum;
  }

  onenrm += vec->onenorm();

  if( vecl )
     onenrm += vecl->onenorm();

  return onenrm;
}


template<typename T>
void StochVectorBase<T>::min( T& m, int& index ) const
{
  T lMin; int lInd;

  if(NULL == parent) {
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
       T lMinlink;
       int lIndlink;
       vecl->min(lMinlink,lIndlink);
       if( lMinlink < lMin )
       {
          lMin = lMinlink;
          lInd = lIndlink + vec->length();
       }
    }

    if(lMin < m) {
      m = lMin;
      index = lInd + parent->n - this->n;
    }
  }

  for(size_t it  =0; it < children.size(); it++) {
    children[it]->min(m,index);
  }

  if(iAmDistrib == 1) {
    T minG = PIPS_MPIgetMin(m, mpiComm);
    // MPI_Allreduce(&m, &minG, 1, MPI_DOUBLE, MPI_MIN, mpiComm); // not working properly in templated version
    m = minG;
  }
}


template<typename T>
void StochVectorBase<T>::max( T& m, int& index ) const
{
   T lMax;
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
         T lMaxlink;
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
      T maxG = PIPS_MPIgetMax(m, mpiComm);
      // MPI_Allreduce(&m, &maxG, 1, MPI_DOUBLE, MPI_MAX, mpiComm); // not working properly in templated version
      m = maxG;
   }
}

template<typename T>
void StochVectorBase<T>::absminVecUpdate(OoqpVectorBase<T>& absminvec) const
{
   StochVectorBase<T>& absminvecStoch = dynamic_cast<StochVectorBase<T>&>(absminvec);
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

template<typename T>
void StochVectorBase<T>::absmaxVecUpdate(OoqpVectorBase<T>& absmaxvec) const
{
   StochVectorBase<T>& absmaxvecStoch = dynamic_cast<StochVectorBase<T>&>(absmaxvec);
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

template<typename T>
void StochVectorBase<T>::absmin(T& m) const
{
   T lMin;

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
        T lMinlink;
        vecl->absmin(lMinlink);
        if( lMinlink < lMin )
           lMin = lMinlink;
     }

     if(lMin < m)
       m = lMin;
   }

   for(size_t it=0; it<children.size(); it++) {
     children[it]->absmin(m);
   }

   if(iAmDistrib == 1) {
     T minG = PIPS_MPIgetMin(m, mpiComm);
     // MPI_Allreduce(&m, &minG, 1, MPI_DOUBLE, MPI_MIN, mpiComm); // not working properly in templated version
     m = minG;
   }
   assert( m >= 0.0 );
}

template<typename T>
void StochVectorBase<T>::absminNonZero(T& m, T zero_eps) const
{
   T min;

   assert(zero_eps >= 0.0);

   vec->absminNonZero(m, zero_eps);

   assert(m >= zero_eps || m == -1.0);


   if( vecl )
   {
      vecl->absminNonZero(min, zero_eps);
      if( min >= 0.0 && (min < m  || m < 0.0) )
      {
         m = min;
         assert(m >= zero_eps || m == -1.0);
      }
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->absminNonZero(min, zero_eps);

      if( min >= 0.0 && (min < m  || m < 0.0) )
      {
         m = min;
         assert(m >= zero_eps || m == -1.0);
      }
   }

   if( iAmDistrib == 1 )
   {
      T minG;
      if( m < 0.0 )
      {
         min = std::numeric_limits<T>::max();
      }
      else
      {
         min = m;
         assert(min >= zero_eps);
      }

      minG = PIPS_MPIgetMin(min, mpiComm);
      // MPI_Allreduce(&min, &minG, 1, MPI_DOUBLE, MPI_MIN, mpiComm); // not working properly in templated version

      if( minG < std::numeric_limits<T>::max() )
         m = minG;
      else
         m = -1.0;
   }
   assert(m >= zero_eps || m == -1.0);
}


template<typename T>
T StochVectorBase<T>::stepbound(const OoqpVectorBase<T> & v_, T maxStep ) const
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);

  T step = this->vec->stepbound(*v.vec, maxStep);

  if( vecl )
  {
     assert(v.vecl);
     T stepl = vecl->stepbound(*v.vecl, maxStep);
     if( stepl < step )
        step = stepl;
  }

  //check tree compatibility
  assert(children.size() == v.children.size());

  for(size_t it = 0; it < children.size(); it++)
    step = children[it]->stepbound(*v.children[it], step);

  if(iAmDistrib == 1) {
    T stepG = PIPS_MPIgetMin(step, mpiComm);
    // MPI_Allreduce(&step, &stepG, 1, MPI_DOUBLE, MPI_MIN, mpiComm); // not working properly in templated version
    step = stepG;
  }
  return step;
}

template<typename T>
T StochVectorBase<T>::findBlocking(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T maxStep,
			      T *w_elt,
			      T *wstep_elt,
			      T *u_elt,
			      T *ustep_elt,
			      int& first_or_second) const
{
  const StochVectorBase<T>& w = *this;
  const StochVectorBase<T>& u = dynamic_cast<const StochVectorBase<T>&>(u_vec);

  const StochVectorBase<T>& wstep = dynamic_cast<const StochVectorBase<T>&>(wstep_vec);
  const StochVectorBase<T>& ustep = dynamic_cast<const StochVectorBase<T>&>(ustep_vec);
  const T local_eps = 1e-14;

  T step = maxStep;

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
  assert( nChildren - u.children.size() == 0);
  assert( wstep.children.size() == ustep.children.size() );
  assert( nChildren - ustep.children.size() == 0);

  for(int it = 0; it < nChildren; it++) {
    step = w.children[it]->findBlocking(*wstep.children[it],
			       *u.children[it],
			       *ustep.children[it],
			       step,
			       w_elt,
			       wstep_elt,u_elt,ustep_elt, first_or_second);
  }

  if(iAmDistrib == 1) {
    T stepG;
    assert(PIPSisLE(step, 1.0));
    assert(PIPSisLE(0.0, step));

    stepG = PIPS_MPIgetMin(step, mpiComm);
    // MPI_Allreduce(&step, &stepG, 1, MPI_DOUBLE, MPI_MIN, mpiComm); not working properly in templated version
    const bool iHaveMinStep = PIPSisEQ(step, stepG, local_eps);

    //we prefer a AllReduce instead of a bcast, since the step==stepG m
    //may occur for two different processes and a deadlock may occur.
    T buffer[5]; //0-primal val, 1-primal step, 2-dual value, 3-step, 4-1st or 2nd

    int count;
    if( iHaveMinStep ) {
      buffer[0]=*w_elt; buffer[1]=*wstep_elt;
      buffer[2]=*u_elt; buffer[3]=*ustep_elt;
      buffer[4]=first_or_second;

      count = 1;
    } else {

      count = 0;
      buffer[0]=buffer[1]=buffer[2]=buffer[3]=buffer[4]= -std::numeric_limits<T>::max();
    }

    count = PIPS_MPIgetSum(count, mpiComm);
    // MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, mpiComm); // not working properly in templated version
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

       mineqrank = PIPS_MPIgetMin(mineqrank, mpiComm);
       // MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm); // not working properly in templated version

       // step==stepG and not smallest rank?
      if( iHaveMinStep && mineqrank != myrank )
         buffer[0]=buffer[1]=buffer[2]=buffer[3]=buffer[4]= -std::numeric_limits<T>::max();
    }

    T bufferOut[5];
    PIPS_MPImaxArray(buffer, bufferOut, 5, mpiComm);

    //MPI_Allreduce(buffer, bufferOut, 5, MPI_DOUBLE, MPI_MAX, mpiComm); // not working properly in templated version

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

template<typename T>
void StochVectorBase<T>::findBlocking_pd(const OoqpVectorBase<T> & wstep_vec,
			      const OoqpVectorBase<T> & u_vec,
			      const OoqpVectorBase<T> & ustep_vec,
			      T& maxStepPri, T& maxStepDual,
			      T& w_elt_p, T& wstep_elt_p, T& u_elt_p, T& ustep_elt_p,
				   T& w_elt_d, T& wstep_elt_d, T& u_elt_d, T& ustep_elt_d,
				   bool& primalBlocking, bool& dualBlocking) const
{
  const StochVectorBase<T>& w = *this;
  const StochVectorBase<T>& u = dynamic_cast<const StochVectorBase<T>&>(u_vec);
  const T local_eps = 1e-14;

  const StochVectorBase<T>& wstep = dynamic_cast<const StochVectorBase<T>&>(wstep_vec);
  const StochVectorBase<T>& ustep = dynamic_cast<const StochVectorBase<T>&>(ustep_vec);

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
      T maxStepGlobalPri, maxStepGlobalDual;
      assert(PIPSisLE(maxStepPri, 1.0) && PIPSisLE(maxStepDual, 1.0));
      assert(PIPSisLE(0.0, maxStepPri) && PIPSisLE(0.0, maxStepDual));

      maxStepGlobalPri = PIPS_MPIgetMin(maxStepPri, mpiComm);
      maxStepGlobalDual = PIPS_MPIgetMin(maxStepDual, mpiComm);
      // MPI_Allreduce(&maxStepPri, &maxStepGlobalPri, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
      // MPI_Allreduce(&maxStepDual, &maxStepGlobalDual, 1, MPI_DOUBLE, MPI_MIN, mpiComm);
      const bool iHaveMinStepPri = PIPSisEQ(maxStepPri, maxStepGlobalPri, local_eps);
      const bool iHaveMinStepDual = PIPSisEQ(maxStepDual, maxStepGlobalDual, local_eps);

      //we prefer a AllReduce instead of a bcast, since the step==stepG
      //may occur for two different processes and a deadlock may occur.
      T buffer[10];
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
               -std::numeric_limits<T>::max();
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
               -std::numeric_limits<T>::max();
         count[1] = 0;
      }

      PIPS_MPIsumArrayInPlace(count, 2, mpiComm);
      // MPI_Allreduce(MPI_IN_PLACE, count, 2, MPI_INT, MPI_SUM, mpiComm);

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

          mineqrank = PIPS_MPIgetMin(mineqrank, mpiComm);
         // MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm);

         // step==stepG and not smallest rank?
         if( iHaveMinStepPri && mineqrank != myrank )
            buffer[0] = buffer[1] = buffer[2] = buffer[3]=buffer[4]= -std::numeric_limits<T>::max();
      }

      // is there more than one process with maxStepDual==stepF?
      if( count[1] > 1 )
      {
         int mineqrank;

         if( iHaveMinStepDual )
            mineqrank = myrank;
         else
            mineqrank = std::numeric_limits<int>::max();

          mineqrank = PIPS_MPIgetMin(mineqrank, mpiComm);
         // MPI_Allreduce(MPI_IN_PLACE, &mineqrank, 1, MPI_INT, MPI_MIN, mpiComm);

         // stepDual==stepF and not smallest rank?
         if( iHaveMinStepDual && mineqrank != myrank )
            buffer[5]=buffer[6]=buffer[7]=buffer[8]=buffer[9]= -std::numeric_limits<T>::max();
      }

      T bufferOut[10];
      PIPS_MPImaxArray(buffer, bufferOut, 10, mpiComm);
      // MPI_Allreduce(buffer, bufferOut, 10, MPI_DOUBLE, MPI_MAX, mpiComm); // not working properly in templated version

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

template<typename T>
void StochVectorBase<T>::componentMult( const OoqpVectorBase<T>& v_ )
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);
  assert(v.children.size() == children.size());

  vec->componentMult(*v.vec);
  if( vecl ) vecl->componentMult(*v.vecl);

  for(size_t it = 0; it < children.size(); it++)
    children[it]->componentMult(*v.children[it]);
}

template<typename T>
void StochVectorBase<T>::componentDiv ( const OoqpVectorBase<T>& v_ )
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);
  assert(v.children.size() == children.size());

  vec->componentDiv(*v.vec);
  if( vecl ) vecl->componentDiv(*v.vecl);

  for(size_t it = 0; it < children.size(); it++)
    children[it]->componentDiv(*v.children[it]);
}

template<typename T>
void StochVectorBase<T>::scalarMult( T num )
{
  vec->scalarMult(num);
  if( vecl ) vecl->scalarMult(num);

  for(size_t it = 0; it < children.size(); it++)
    children[it]->scalarMult(num);
}

template<typename T>
void StochVectorBase<T>::writeToStreamAll( std::ostream& out ) const
{
   int rank;
   MPI_Comm_rank(mpiComm, &rank);
   int world_size;
   MPI_Comm_size(mpiComm, &world_size);
   MPI_Status status;
   int l;
   std::stringstream sout;

   if( rank == 0)
   {
      sout << "----" << std::endl;
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
         std::string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }
      if( vecl )
      {
         sout << "---" << std::endl;
         vecl->writeToStreamAllStringStream(sout);
      }
      sout << "----" << std::endl;
      out << sout.str();
   }
   else if( iAmDistrib==1 )
   { // rank != 0
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->writeToStreamAllChild(sout);

      std::string str = sout.str();
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, rank, mpiComm);

   }

   if( iAmDistrib == 1 )
      MPI_Barrier(mpiComm);
}

template<typename T>
void StochVectorBase<T>::writeToStreamAllChild( std::stringstream& sout ) const
{
   sout << "--" << std::endl;
   vec->writeToStreamAllStringStream(sout);

   for( size_t it = 0; it < children.size(); it++ ){
      sout << "-- " << std::endl;
      children[it]->writeToStreamAllChild(sout);
   }

   if( vecl )
   {
      sout << "---" << std::endl;
      vecl->writeToStreamAllStringStream(sout);
   }
}

template<typename T>
void StochVectorBase<T>::writeToStream( std::ostream& out ) const
{
  out << "---" << std::endl;
  vec->writeToStream(out);
  if( vecl ) vecl->writeToStream(out);
  out << "~~~" << std::endl;
  //for(size_t it=0; it<children.size(); it++)
  //  children[it]->writeToStream(out);
}

template<typename T>
void StochVectorBase<T>::writefToStream( std::ostream& out,
				  const char format[] ) const
{
  vec->writefToStream(out, format);
  if( vecl ) vecl->writefToStream(out, format);

  for(size_t it=0; it<children.size(); it++)
    children[it]->writefToStream(out, format);
}

template<typename T>
void StochVectorBase<T>::writeMPSformatRhs( std::ostream& out, int rowType, const OoqpVectorBase<T>* irhs) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   std::string rt;
   if( rowType == 0 )
      rt = "E";
   else if( rowType == 1 )
      rt = "L";
   else if( rowType == 2 )
      rt = "G";
   else
      assert(0);

   const StochVectorBase<T>* ic = NULL;
   if( irhs )
      ic = dynamic_cast<const StochVectorBase<T>*>(irhs);

   if( myRank == 0 )
   {
      std::string rowNameStub = " B row_";
      rowNameStub += rt;
      rowNameStub += "_R_";
      if( irhs && ic )
         vec->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->vec));
      else
         vec->writeMPSformatOnlyRhs( out, rowNameStub, NULL);
      if(vecl)
      {
         rowNameStub = " B row_";
         rowNameStub += rt;
         rowNameStub += "_L_";
         if( irhs )
            vecl->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->vecl));
         else
            vecl->writeMPSformatOnlyRhs( out, rowNameStub, NULL);
      }
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm;
      sstm << " B row_" << rt << "_" << it << "_";
      std::string rowNameStub = sstm.str();
      if( irhs )
         children[it]->vec->writeMPSformatOnlyRhs( out, rowNameStub, dynamic_cast<const SimpleVectorBase<T>*>(ic->children[it]->vec));
      else
         children[it]->vec->writeMPSformatOnlyRhs( out, rowNameStub, NULL);
   }
}

template<typename T>
void StochVectorBase<T>::writeMPSformatBounds(std::ostream& out, const OoqpVectorBase<T>* ix, bool upperBound) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);

   const StochVectorBase<T>* ixStoch = dynamic_cast<const StochVectorBase<T>*>(ix);

   if( myRank==0 )
   {
      std::string varNameStub = "var_L_";
      vec->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->vec), upperBound);
   }
   for(int it=0; it<(int)children.size(); it++)
   {
      std::stringstream sstm2;
      sstm2 << "var_" << it << "_";
      std::string varNameStub = sstm2.str();
      children[it]->vec->writeMPSformatBoundsWithVar(out, varNameStub, (ixStoch->children[it]->vec), upperBound);
   }
}

/** this += alpha * x */
template<typename T>
void StochVectorBase<T>::axpy ( T alpha, const OoqpVectorBase<T>& x_ )
{
  const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);
  assert(x.children.size() == children.size());

  if( alpha == 0.0)
     return;

  vec->axpy(alpha, *x.vec);

  if( vecl )
  {
     assert(x.vecl);
     vecl->axpy(alpha, *x.vecl);
  }

  for(size_t it = 0; it<children.size(); it++)
    children[it]->axpy(alpha, *x.children[it]);
}

/** this += alpha * x * z */
template<typename T>
void StochVectorBase<T>::axzpy ( T alpha, const OoqpVectorBase<T>& x_, const OoqpVectorBase<T>& z_ )
{
  const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);
  const StochVectorBase<T>& z = dynamic_cast<const StochVectorBase<T>&>(z_);
  assert(x.children.size() == children.size());
  assert(z.children.size() == children.size());

  vec->axzpy(alpha, *x.vec, *z.vec);

  if( vecl )
  {
     assert(x.vecl);
     assert(z.vecl);
     vecl->axzpy(alpha, *x.vecl, *z.vecl);
  }

  for(size_t it = 0; it < children.size(); it++)
    children[it]->axzpy(alpha, *x.children[it], *z.children[it]);
}

/** this += alpha * x / z */
template<typename T>
void StochVectorBase<T>::axdzpy( T alpha, const OoqpVectorBase<T>& x_, const OoqpVectorBase<T>& z_ )
{
  const StochVectorBase<T>& x = dynamic_cast<const StochVectorBase<T>&>(x_);
  const StochVectorBase<T>& z = dynamic_cast<const StochVectorBase<T>&>(z_);
  assert(x.children.size() == children.size());
  assert(z.children.size() == children.size());

  vec->axdzpy(alpha, *x.vec, *z.vec);

  if( vecl )
  {
    assert(x.vecl);
    assert(z.vecl);
    vecl->axdzpy(alpha, *x.vecl, *z.vecl);
  }

  for(size_t it = 0; it < children.size(); it++)
    children[it]->axdzpy(alpha, *x.children[it], *z.children[it]);
}


template<typename T>
void StochVectorBase<T>::addConstant( T c )
{
  vec->addConstant(c);

  if( vecl ) vecl->addConstant(c);

  for(size_t it = 0; it < children.size(); it++)
    children[it]->addConstant(c);
}


template<typename T>
void StochVectorBase<T>::gondzioProjection( T rmin, T rmax )
{
  vec->gondzioProjection( rmin, rmax );


  if( vecl ) vecl->gondzioProjection( rmin, rmax );

  for(size_t it = 0; it < children.size(); it++)
    children[it]->gondzioProjection( rmin, rmax );
}

template<typename T>
T StochVectorBase<T>::dotProductWith( const OoqpVectorBase<T>& v_ ) const
{
  const StochVectorBase<T>& v = dynamic_cast<const StochVectorBase<T>&>(v_);
  assert(v.children.size() == children.size());

  T dotProd = 0.0;

  for(size_t it = 0; it < children.size(); it++)
    dotProd += children[it]->dotProductWith(*v.children[it]);

  assert(!vecl || v.vecl);

  if(iAmDistrib == 1) {
    T dotProdG = 0.0;

    dotProdG = PIPS_MPIgetSum(dotProd, mpiComm);
    // MPI_Allreduce(&dotProd, &dotProdG, 1, MPI_DOUBLE, MPI_SUM, mpiComm); // not working properly in templated version

    dotProd = dotProdG;
  }

  dotProd += vec->dotProductWith(*v.vec);

  if( vecl )
    dotProd += vecl->dotProductWith(*v.vecl);

  return dotProd;
}

template<typename T>
T StochVectorBase<T>::dotProductSelf(T scaleFactor) const
{
  T dotSelf = 0.0;

  for(size_t it = 0; it < children.size(); it++)
     dotSelf += children[it]->dotProductSelf(scaleFactor);

  if(iAmDistrib == 1) {
    T dotSelfG = PIPS_MPIgetSum(dotSelf, mpiComm);
    // MPI_Allreduce(&dotSelf, &dotSelfG, 1, MPI_DOUBLE, MPI_SUM, mpiComm); // not working properly in templated version

    dotSelf = dotSelfG;
  }

  dotSelf += vec->dotProductSelf(scaleFactor);

  if( vecl )
     dotSelf += vecl->dotProductSelf(scaleFactor);

  return dotSelf;
}

/** Return the inner product <this + alpha * mystep, yvec + beta * ystep >
 */
template<typename T>
T StochVectorBase<T>::shiftedDotProductWith( T alpha, const OoqpVectorBase<T>& mystep_,
					const OoqpVectorBase<T>& yvec_,
					T beta,  const OoqpVectorBase<T>& ystep_ ) const
{
  const StochVectorBase<T>& mystep = dynamic_cast<const StochVectorBase<T>&>(mystep_);
  const StochVectorBase<T>& yvec   = dynamic_cast<const StochVectorBase<T>&>(yvec_);
  const StochVectorBase<T>& ystep  = dynamic_cast<const StochVectorBase<T>&>(ystep_);


  T dotProd = 0.0;
  for(size_t it=0; it<children.size(); it++)
    dotProd += children[it]->shiftedDotProductWith(alpha, *mystep.children[it],
						   *yvec.children[it],
						   beta, *ystep.children[it]);
  if(iAmDistrib) {
    T dotProdG=0.0;
    dotProdG = PIPS_MPIgetSum(dotProd, mpiComm);
    // MPI_Allreduce(&dotProd, &dotProdG, 1, MPI_DOUBLE, MPI_SUM, mpiComm); // not working properly in templated version
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

template<typename T>
void StochVectorBase<T>::negate()
{
  vec->negate();
  if( vecl ) vecl->negate();

  for(size_t it = 0; it < children.size(); it++)
    children[it]->negate();
}

template<typename T>
void StochVectorBase<T>::invert()
{
  vec->invert();

  if( vecl ) vecl->invert();

  for(size_t it = 0; it < children.size(); it++)
    children[it]->invert();
}

template<typename T>
void StochVectorBase<T>::invertSave(T zeroReplacementVal)
{
  vec->invertSave(zeroReplacementVal);

  if( vecl ) vecl->invertSave(zeroReplacementVal);

  for( size_t it = 0; it < children.size(); it++ )
    children[it]->invertSave(zeroReplacementVal);
}

template<typename T>
void StochVectorBase<T>::applySqrt()
{
   vec->applySqrt();

   if(vecl) vecl->applySqrt();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->applySqrt();
}

template<typename T>
void StochVectorBase<T>::roundToPow2()
{
  vec->roundToPow2();

  if( vecl ) vecl->roundToPow2();

  for( size_t it = 0; it < children.size(); it++ )
     children[it]->roundToPow2();
}

template<typename T>
bool StochVectorBase<T>::allPositive() const
{
  //!parallel
  bool allPos = vec->allPositive() && ((vecl != NULL) ? vecl->allPositive() : true);
  if (!allPos) return false;

  for(size_t it = 0; it < children.size() && allPos; it++)
    allPos = children[it]->allPositive();

  return allPos;
}

template<typename T>
bool StochVectorBase<T>::matchesNonZeroPattern( const OoqpVectorBase<T>& select_ ) const
{
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
  assert(children.size() == select.children.size());

  bool match = vec->matchesNonZeroPattern(*select.vec);

  if( vecl )
  {
     assert(select.vecl);
     match = match && vecl->matchesNonZeroPattern(*select.vecl);
  }

  if(!match) return false;

  for(size_t it = 0; it < children.size() && match; it++)
    match = children[it]->matchesNonZeroPattern(*select.children[it]);

  return match;
}

template<typename T>
void StochVectorBase<T>::selectNonZeros( const OoqpVectorBase<T>& select_ )
{
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
  assert(children.size() == select.children.size());

  vec->selectNonZeros(*select.vec);

  if( vecl )
  {
     assert(select.vecl);
     vecl->selectNonZeros(*select.vecl);
  }

  for(size_t it = 0; it < children.size(); it++)
    children[it]->selectNonZeros(*select.children[it]);
}

template<typename T>
long long StochVectorBase<T>::numberOfNonzeros() const
{
  //!opt - store the number of nnz to avoid communication
  long long nnz = 0;

  for(size_t it = 0; it < children.size(); it++)
    nnz += children[it]->numberOfNonzeros();

  if(iAmDistrib) {
    long long nnzG = 0;
    nnzG = PIPS_MPIgetSum(nnz, mpiComm);
    // MPI_Allreduce(&nnz, &nnzG, 1, MPI_LONG_LONG, MPI_SUM, mpiComm);
    nnz = nnzG;
  }
  nnz += vec->numberOfNonzeros();

  if( vecl ) nnz += vecl->numberOfNonzeros();

  return nnz;
}

template<typename T>
void StochVectorBase<T>::addSomeConstants( T c, const OoqpVectorBase<T>& select_ )
{
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
  assert(children.size() == select.children.size());

  vec->addSomeConstants(c, *select.vec);

  if( vecl )
  {
     assert(select.vecl);
     vecl->addSomeConstants(c, *select.vecl);
  }

  for(size_t it = 0; it < children.size(); it++)
    children[it]->addSomeConstants(c, *select.children[it]);
}

template<typename T>
void StochVectorBase<T>::writefSomeToStream( std::ostream& out,
			 const char format[],
			 const OoqpVectorBase<T>& select_ ) const
{
  assert( "Have not been yet implemented" && 0 );
}

template<typename T>
void StochVectorBase<T>::axdzpy( T alpha, const OoqpVectorBase<T>& x_,
		       const OoqpVectorBase<T>& z_, const OoqpVectorBase<T>& select_ )
{
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
  const StochVectorBase<T>&      x = dynamic_cast<const StochVectorBase<T>&>(x_);
  const StochVectorBase<T>&      z = dynamic_cast<const StochVectorBase<T>&>(z_);

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

  for(size_t it = 0; it < children.size(); it++)
    children[it]->axdzpy(alpha, *x.children[it], *z.children[it], *select.children[it]);
}

template<typename T>
bool StochVectorBase<T>::somePositive( const OoqpVectorBase<T>& select_ ) const
{
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);
  assert(children.size() == select.children.size());

  //!parallel stuff needed

  bool somePos = vec->somePositive(*select.vec);

  if( vecl )
  {
    assert(select.vecl);
    somePos = somePos && vecl->somePositive(*select.vecl);
  }

  for(size_t it = 0; it < children.size() && somePos; it++)
    somePos = children[it]->somePositive(*select.children[it]);

  return somePos;
}

template<typename T>
void StochVectorBase<T>::divideSome( const OoqpVectorBase<T>& div_, const OoqpVectorBase<T>& select_ )
{
  const StochVectorBase<T>& div    = dynamic_cast<const StochVectorBase<T>&>(div_);
  const StochVectorBase<T>& select = dynamic_cast<const StochVectorBase<T>&>(select_);

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

template<typename T>
void StochVectorBase<T>::copyIntoArray( T v[] ) const
{
  assert( "Not supported" && 0 );
}

template<typename T>
void StochVectorBase<T>::copyFromArray( const T v[] )
{
  assert( "Not supported" && 0 );
}

template<typename T>
void StochVectorBase<T>::copyFromArray( const char v[] )
{
  assert( "Not supported" && 0 );
}

template<typename T>
void StochVectorBase<T>::removeEntries( const OoqpVectorBase<T>& select )
{
   const StochVectorBase<T>& selectStoch = dynamic_cast<const StochVectorBase<T>&>(select);

   assert(children.size() == selectStoch.children.size());

   vec->removeEntries(*selectStoch.vec);
   this->n = vec->n;

   if( vecl )
   {
      assert(selectStoch.vecl);
      vecl->removeEntries(*selectStoch.vecl);
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->removeEntries(*selectStoch.children[it]);
      this->n += children[it]->n;
   }
}

template<typename T>
void StochVectorBase<T>::permuteVec0Entries(const std::vector<unsigned int>& permvec)
{
   dynamic_cast<SimpleVectorBase<T>*>(vec)->permuteEntries(permvec);
}

template<typename T>
void StochVectorBase<T>::permuteLinkingEntries(const std::vector<unsigned int>& permvec)
{
   if( vecl )
      dynamic_cast<SimpleVectorBase<T>*>(vecl)->permuteEntries(permvec);
}

template<typename T>
std::vector<T> StochVectorBase<T>::gatherStochVector() const
{
   const SimpleVectorBase<T>& firstvec = dynamic_cast<const SimpleVectorBase<T>&>(*vec);
   const size_t nChildren = children.size();

   int myrank;
   MPI_Comm_rank(mpiComm, &myrank);
   int mysize;
   MPI_Comm_size(mpiComm, &mysize);

   std::vector<T> gatheredVecLocal(0);

   for( size_t i = 0; i < nChildren; ++i )
   {
      const SimpleVectorBase<T>& vec = dynamic_cast<const SimpleVectorBase<T>&>(*children[i]->vec);

      if( vec.length() > 0 )
         gatheredVecLocal.insert(gatheredVecLocal.end(), &vec[0], &vec[0] + vec.length());
   }

   size_t solLength = firstvec.length();

   // final vector
   std::vector<T> gatheredVec(0);

   if( mysize > 0 )
   {
      // get all lengths
      std::vector<int> recvcounts(mysize);
      std::vector<int> recvoffsets(mysize);

      int mylength = int(gatheredVecLocal.size());

      PIPS_MPIallgather(&mylength, 1, &recvcounts[0], 1, mpiComm);
      // MPI_Allgather(&mylength, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, mpiComm);

      // all-gather local components
      recvoffsets[0] = 0;
      for( size_t i = 1; i < size_t(mysize); ++i )
         recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

      if( myrank == 0 )
      {
         solLength += recvoffsets[mysize - 1] + recvcounts[mysize - 1];
         gatheredVec = std::vector<T>(solLength);

         PIPS_MPIgatherv(&gatheredVecLocal[0], mylength, &gatheredVec[0] + firstvec.length(),
            &recvcounts[0], &recvoffsets[0], 0, mpiComm);
         // MPI_Gatherv(&gatheredVecLocal[0], mylength, MPI_DOUBLE,
         //       &gatheredVec[0] + firstvec.length(), &recvcounts[0],
         //       &recvoffsets[0], MPI_DOUBLE, 0, mpiComm);
      }
      else
      {
        T dummy;
        PIPS_MPIgatherv(&gatheredVecLocal[0], mylength, &dummy, &recvcounts[0], &recvoffsets[0], 0, mpiComm);
         // MPI_Gatherv(&gatheredVecLocal[0], mylength, MPI_DOUBLE, 0,
         //       &recvcounts[0], &recvoffsets[0], MPI_DOUBLE, 0, mpiComm);
      }
   }
   else
   {
      solLength += gatheredVecLocal.size();

      gatheredVec = std::vector<T>(solLength);

      std::copy(gatheredVecLocal.begin(), gatheredVecLocal.end(), gatheredVec.begin() + firstvec.length());
   }

   if( myrank == 0 )
   {
      std::copy(&firstvec[0], &firstvec[0] + firstvec.length(), &gatheredVec[0]);

      if( vecl && vecl->length() > 0 )
      {
         const SimpleVectorBase<T>& linkvec = dynamic_cast<const SimpleVectorBase<T>&>(*vecl);
         gatheredVec.insert(gatheredVec.end(), &linkvec[0], &linkvec[0] + linkvec.length());
      }
   }
   return gatheredVec;
}

// is root node data of StochVector same on all procs?
template<typename T>
bool StochVectorBase<T>::isRootNodeInSync() const
{
   assert( vec);
   assert(mpiComm);

   bool in_sync = true;
   const SimpleVectorBase<T>& vec_simple = dynamic_cast<const SimpleVectorBase<T>&>(*vec);

   /* no need to check if not distributed or not at root node */
   if( !iAmDistrib || parent != NULL)
      return in_sync;

   int my_rank, world_size;
   MPI_Comm_rank(mpiComm, &my_rank);
   MPI_Comm_size(mpiComm, &world_size);

   /* if there is a linking part we have to chekc it as well */
   const int vec_length = vec_simple.length();
   const int vecl_length = (vecl) ? dynamic_cast<const SimpleVectorBase<T>&>(*vecl).length() : 0;

   const long long count = vec_length + vecl_length;

   assert( count < std::numeric_limits<int>::max());

   /* mpi reduce on vector */
   std::vector<T> sendbuf(count, 0.0);
   std::vector<T> recvbuf(count, 0.0);
   std::copy(vec_simple.elements(), vec_simple.elements() + vec_simple.length(), sendbuf.begin());

   if( vecl )
   {
      const SimpleVectorBase<T>& vecl_simple = dynamic_cast<const SimpleVectorBase<T>&>(*vecl);
      std::copy(vecl_simple.elements(), vecl_simple.elements() + vecl_simple.length(),
            sendbuf.begin() + vec_simple.length());
   }
   PIPS_MPImaxArray(&sendbuf[0], &recvbuf[0], count, mpiComm);
   // MPI_Allreduce(&sendbuf[0], &recvbuf[0], count, MPI_DOUBLE, MPI_MAX, mpiComm); // not working properly in templated version

   for( int i = 0; i < count; ++i )
   {
      if( !PIPSisEQ(sendbuf[i], recvbuf[i]) )
      {
         /* someone else had a higher value here */
         in_sync = false;
      }
   }

   return in_sync;
}

template class StochVectorBase<int>;
// template class StochVectorBase<bool>;
template class StochVectorBase<double>;
