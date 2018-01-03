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
  return sqrt(this->dotProductWith(*this));
}

double StochVector::onenorm()
{
  double onenrm = vec->onenorm();

  if( vecl ) onenrm += vecl->onenorm();

  for(size_t it=0; it<children.size(); it++)
    onenrm += children[it]->onenorm();

  //!parallel
  assert(false);
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

  double step = maxStep;
  
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
    MPI_Allreduce(&step, &stepG, 1, MPI_DOUBLE, MPI_MIN, mpiComm);

    //we prefer a AllReduce instead of a bcast, since the step==stepG m
    //may occur for two different processes and a deadlock may occur.
    double buffer[5]; //0-primal val, 1-primal step, 2-dual value, 3-step, 4-1st or 2nd
    if(step==stepG) {
      buffer[0]=*w_elt; buffer[1]=*wstep_elt; 
      buffer[2]=*u_elt; buffer[3]=*ustep_elt;
      buffer[4]=first_or_second;
    } else {
      buffer[0]=buffer[1]=buffer[2]=buffer[3]=buffer[4]= -std::numeric_limits<double>::max();
    }
    double bufferOut[5];
    MPI_Allreduce(buffer, bufferOut, 5, MPI_DOUBLE, MPI_MAX, mpiComm);

    *w_elt = bufferOut[0]; *wstep_elt=bufferOut[1];
    *u_elt = bufferOut[2]; *ustep_elt=bufferOut[3];

    //first_or_second  negative means no blocking, so set it to 0.
    first_or_second = bufferOut[4]<0?0:(int)bufferOut[4];
    assert(first_or_second==0 || first_or_second==1 || first_or_second==2);
    step=stepG;
  }
  return step;
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

/** this += alpha * x */
void StochVector::axpy  ( double alpha, OoqpVector& x_ )
{
  StochVector& x = dynamic_cast<StochVector&>(x_);
  assert(x.children.size() == children.size());

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

  if(iAmDistrib==1) {
    double dotProdG=0.0;
    MPI_Allreduce(&dotProd, &dotProdG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);

    dotProd = dotProdG;
  }

  dotProd += vec->dotProductWith(*v.vec);

  if( vecl )
  {
    assert(v.vecl);
    dotProd += vecl->dotProductWith(*v.vecl);
  }

  return dotProd;
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


void StochVector::invertSave(double zeroReplacementVal)
{
  vec->invertSave(zeroReplacementVal);

  if( vecl ) vecl->invertSave(zeroReplacementVal);

  for( size_t it = 0; it < children.size(); it++ )
    children[it]->invertSave(zeroReplacementVal);
}

void StochVector::invert()
{
  vec->invert();

  if( vecl ) vecl->invert();

  for(size_t it=0; it<children.size(); it++) 
    children[it]->invert();
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
