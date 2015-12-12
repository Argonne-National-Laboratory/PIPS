/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "StochVector.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "VectorUtilities.h"

#include <cassert>
#include <cstring>
#include <iostream>
#include <limits>
#include <math.h>

StochVector::StochVector(int n_, MPI_Comm mpiComm_, int isDistributed/*=-1*/)
  : OoqpVector(n_), parent(NULL), mpiComm(mpiComm_), treeIDX(0), firstDoNumOfNonZero(true),
    iAmDistrib(isDistributed)
{
  vec = new SimpleVector(n_);

  if(-1==iAmDistrib && MPI_COMM_NULL!=mpiComm) {
    int size;
    MPI_Comm_size(mpiComm, &size);
    iAmDistrib = size==1?0:1;
  }
}

StochVector::StochVector(int n_, int const treeIDX_in, MPI_Comm mpiComm_, int isDistributed/*=-1*/)
  : OoqpVector(n_), parent(NULL), mpiComm(mpiComm_), treeIDX(treeIDX_in), firstDoNumOfNonZero(true),
    iAmDistrib(isDistributed)
{
  vec = new SimpleVector(n_);

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
  
  if (vec)
    delete vec;
}

OoqpVector* StochVector::dataClone() const
{
  OoqpVector* clone = new SimpleVector(vec->length());
  return clone;
}

StochVector* StochVector::clone() const
{
  StochVector* clone = new StochVector(this->vec->length(), mpiComm);

  for(size_t it=0; it<this->children.size(); it++) {
    clone->AddChild(this->children[it]->clone());
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
StochVector::jointCopyFromXSYZ(StochVector& v1, StochVector& v2, StochVector& v3,StochVector& v4)
{
  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);
  SimpleVector& sv1 = dynamic_cast<SimpleVector&>(*v1.vec);
  SimpleVector& sv2 = dynamic_cast<SimpleVector&>(*v2.vec);
  SimpleVector& sv3 = dynamic_cast<SimpleVector&>(*v3.vec);
  SimpleVector& sv4 = dynamic_cast<SimpleVector&>(*v4.vec);

  int n1 = sv1.length();
  int n2 = sv2.length();
  int n3 = sv3.length();
  int n4 = sv4.length();

  assert(n1+n2+n3+n4 == sv.length());
  
  if(n1>0)
    memcpy(&sv[0], &sv1[0], n1*sizeof(double));

  if(n2>0)
    memcpy(&sv[n1], &sv2[0], n2*sizeof(double));

  if(n3>0)
    memcpy(&sv[n1+n2], &sv3[0], n3*sizeof(double));
  
  if(n4>0)
    memcpy(&sv[n1+n2+n3], &sv4[0], n4*sizeof(double));

  for(size_t it=0; it<children.size(); it++) {
    children[it]->jointCopyFromXSYZ(*v1.children[it], 
				*v2.children[it], *v3.children[it],*v4.children[it]);
  }

}

void 
StochVector::jointCopyToXSYZ(StochVector& v1, StochVector& v2, StochVector& v3,StochVector& v4)
{
  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);
  SimpleVector& sv1 = dynamic_cast<SimpleVector&>(*v1.vec);
  SimpleVector& sv2 = dynamic_cast<SimpleVector&>(*v2.vec);
  SimpleVector& sv3 = dynamic_cast<SimpleVector&>(*v3.vec);
  SimpleVector& sv4 = dynamic_cast<SimpleVector&>(*v4.vec);

  int n1 = sv1.length();
  int n2 = sv2.length();
  int n3 = sv3.length();
  int n4 = sv4.length();

  assert(n1+n2+n3+n4 == sv.length());
 
  if(n1>0)
    memcpy(&sv1[0], &sv[0], n1*sizeof(double));

  if(n2>0)
    memcpy(&sv2[0], &sv[n1], n2*sizeof(double));

  if(n3>0)
    memcpy(&sv3[0], &sv[n1+n2], n3*sizeof(double));

  if(n4>0)
    memcpy(&sv4[0], &sv[n1+n2+n3], n4*sizeof(double));

  for(size_t it=0; it<children.size(); it++) {
    children[it]->jointCopyToXSYZ(*v1.children[it], 
			      *v2.children[it], *v3.children[it], *v4.children[it]);
  }
}


int StochVector::isKindOf( int kind )
{
  return kind==kStochVector;
}

void StochVector::scale( double alpha )
{
  vec->scale(alpha);
  
  for(size_t it=0; it<children.size(); it++)
    children[it]->scale(alpha);
}

void StochVector::setToZero()
{
  vec->setToZero();
  
  for(size_t it=0; it<children.size(); it++)
    children[it]->setToZero();
}

void StochVector::setToConstant( double c) 
{
  vec->setToConstant(c);
  
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

  //check tree compatibility
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

  return infnrm; 
}

double StochVector::twonorm()
{
  return sqrt(this->dotProductWith(*this));
}

double StochVector::onenorm()
{
  double onenrm = 0.0;

  for(size_t it=0; it<children.size(); it++)
    onenrm += children[it]->onenorm();

  if(iAmDistrib) {
    double onenrmG=0.0;
    MPI_Allreduce(&onenrm, &onenrmG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
    onenrm = onenrmG;
  }

  onenrm += vec->onenorm();

  return onenrm; 
}

// FIXME_ME: index is wrong! Do not use it in parallel
void StochVector::min( double& m, int& index )
{
  double lMin; int lInd;

  if(NULL==parent) {
    vec->min(m,index);
  } else {
    vec->min(lMin,lInd);
    
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

double StochVector::stepbound(OoqpVector & v_, double maxStep )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);

  double step = this->vec->stepbound(*v.vec, maxStep);

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
    //ay occur for two different processes and a deadlock may occur.
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

  for(size_t it=0; it<children.size(); it++) 
    children[it]->componentMult(*v.children[it]);
}

void StochVector::componentDiv ( OoqpVector& v_ )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);
  assert(v.children.size() == children.size());

  vec->componentDiv(*v.vec);
  for(size_t it=0; it<children.size(); it++) 
    children[it]->componentDiv(*v.children[it]);
}

void StochVector::scalarMult( double num )
{
  vec->scalarMult(num);
  for(size_t it=0; it<children.size(); it++) 
    children[it]->scalarMult(num);
}

void StochVector::writeToStream( std::ostream& out ) const
{
  out << "---" << std::endl;
  vec->writeToStream(out);
  out << "~~~" << std::endl;
  //for(size_t it=0; it<children.size(); it++) 
  //  children[it]->writeToStream(out);
}

void StochVector::writefToStream( std::ostream& out,
				  const char format[] ) const
{
  vec->writefToStream(out, format);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->writefToStream(out, format);
}

/** this += alpha * x */
void StochVector::axpy  ( double alpha, OoqpVector& x_ )
{
  StochVector& x = dynamic_cast<StochVector&>(x_);
  assert(x.children.size() == children.size());

  vec->axpy(alpha, *x.vec);

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

  for(size_t it=0; it<children.size(); it++)
    children[it]->axdzpy(alpha, *x.children[it], *z.children[it]);
}


void StochVector::addConstant( double c )
{
  vec->addConstant(c);
  for(size_t it=0; it<children.size(); it++) 
    children[it]->addConstant(c);
}


void StochVector::gondzioProjection( double rmin, double rmax )
{
  assert( "Have not been yet implemented" && 0 );

  // FIXME_NY: Is the followings correct? (NY added them)
  vec->gondzioProjection(rmin, rmax);
  for(size_t it=0; it<children.size(); it++) 
    children[it]->gondzioProjection(rmin, rmax);
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
  return dotProd;
}


/** Return the inner product   <this + alpha * mystep, yvec + beta * ystep >
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

  return dotProd;
}

void StochVector::negate()
{
  vec->negate();
  for(size_t it=0; it<children.size(); it++) 
    children[it]->negate();
}

void StochVector::invert()
{
  vec->invert();
  for(size_t it=0; it<children.size(); it++) 
    children[it]->invert();
}

#if 0
int StochVector::allPositive()
{
  //!parallel
  int allPos = vec->allPositive();
  if (!allPos) return 0;

  for(size_t it=0; it<children.size() && allPos; it++) 
    allPos = children[it]->allPositive();

  return allPos;
}
#else
int StochVector::allPositive()
{
  int allPos = vec->allPositive();
  if (!allPos) return 0;

  for(size_t it=0; it<children.size() && allPos; it++) 
    allPos = children[it]->allPositive();

  if(iAmDistrib) {
    int allPosG=0;
    MPI_Allreduce(&allPos, &allPosG, 1, MPI_INT, MPI_MIN, mpiComm);
    allPos = allPosG;
  }

  return allPos;
}
#endif

#if 0
int StochVector::matchesNonZeroPattern( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  int match = vec->matchesNonZeroPattern(*select.vec);
  if(!match) return 0;

  for(size_t it=0; it<children.size() && match; it++) 
    match = children[it]->matchesNonZeroPattern(*select.children[it]);

  return match;
}
#else
int StochVector::matchesNonZeroPattern( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  int match = vec->matchesNonZeroPattern(*select.vec);
  if(!match) return 0;

  for(size_t it=0; it<children.size() && match; it++) 
    match = children[it]->matchesNonZeroPattern(*select.children[it]);

  if(iAmDistrib) {
    int matchG=0;
    MPI_Allreduce(&match, &matchG, 1, MPI_INT, MPI_MIN, mpiComm);
    match = matchG;
  }
  
  return match;
}
#endif

void StochVector::selectNonZeros( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  vec->selectNonZeros(*select.vec);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->selectNonZeros(*select.children[it]);
}


#if 0
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

  return nnz;
}
#else
long long StochVector::numberOfNonzeros()
{
  //!opt - store the number of nnz to avoid communication
  long long nnz = 0;

  if(firstDoNumOfNonZero){

    for(size_t it=0; it<children.size(); it++) 
      nnz += children[it]->numberOfNonzeros();

    if(iAmDistrib) {
      long long nnzG = 0;
      MPI_Allreduce(&nnz, &nnzG, 1, MPI_LONG_LONG, MPI_SUM, mpiComm);
      nnz = nnzG;
    }
    nnz += vec->numberOfNonzeros();
	nnzNonZeros = nnz;
  }
  else{
	nnz = nnzNonZeros;
  }

  return nnz;
}

#endif



void StochVector::addSomeConstants( double c, OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  vec->addSomeConstants(c, *select.vec);

  for(size_t it=0; it<children.size(); it++) 
    children[it]->addSomeConstants(c, *select.children[it]);
}

void StochVector::writefSomeToStream( std::ostream& out,
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

  for(size_t it=0; it<children.size(); it++)
    children[it]->axdzpy(alpha, *x.children[it], *z.children[it], *select.children[it]);
}

#if 0
int StochVector::somePositive( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  //!parallel stuff needed

  int somePos = vec->somePositive(*select.vec);

  for(size_t it=0; it<children.size() && somePos; it++)
    somePos = children[it]->somePositive(*select.children[it]);
  
  return somePos;
}
#else
int StochVector::somePositive( OoqpVector& select_ )
{
  StochVector& select = dynamic_cast<StochVector&>(select_);
  assert(children.size() == select.children.size());

  int somePos = vec->somePositive(*select.vec);
  if(!somePos) return 0;

  for(size_t it=0; it<children.size() && somePos; it++)
    somePos = children[it]->somePositive(*select.children[it]);

  if(iAmDistrib) {
    int somePosG=0;
    MPI_Allreduce(&somePos, &somePosG, 1, MPI_INT, MPI_MIN, mpiComm);
    somePos = somePosG;
  }
  
  return somePos;
}
#endif


void StochVector::divideSome( OoqpVector& div_, OoqpVector& select_ )
{
  StochVector& div    = dynamic_cast<StochVector&>(div_);
  StochVector& select = dynamic_cast<StochVector&>(select_);

  assert(children.size() == div.   children.size());
  assert(children.size() == select.children.size());

  vec->divideSome(*div.vec, *select.vec);

  for(size_t it=0; it<children.size(); it++)
    children[it]->divideSome(*div.children[it], *select.children[it]);
}

void StochVector::copyIntoArray( double v[] )
{
//  assert( "Not supported" && 0 );

  copyIntoArrayFromTo(v,0,n,0,n);
}

void StochVector::copyFromArray( double v[] )
{
//  assert( "Not supported" && 0 );
  this->copyFromArrayFromTo(v,0,n,0,n);
}

void StochVector::copyFromArray( char v[] )
{
  assert( "Not supported" && 0 );
}





/* following routines are added by Naiyuan 2013 */

void StochVector::print()
{
  printf(" *** [%d]: \n", treeIDX);
  vec->print();
  
  for(size_t it=0; it<children.size(); it++) 
    children[it]->print();
}

void StochVector::findBlockingPD(OoqpVector & wstep_vec, 
				      OoqpVector & u_vec, 
				      OoqpVector & ustep_vec, 
				      double *w_elt, 
				      double *wstep_elt,
				      double *u_elt, 
				      double *ustep_elt,
				      int& first_or_second, double * alphaPri, double * alphaDual)
{
  assert( "Not supported" && 0 );
}


// log function
double StochVector::sumLog(OoqpVector* select_)
{
  StochVector* select = NULL;
  double sum = 0;
  
  if(select_ != NULL){
    StochVector* select = dynamic_cast<StochVector*>(select_);
    assert(children.size() == select->children.size());

	for(size_t it=0; it<children.size(); it++)
      sum += children[it]->sumLog(select->children[it]);

  }else{
	for(size_t it=0; it<children.size(); it++)
	  sum += children[it]->sumLog(NULL);
  }

  if(iAmDistrib==1) {
    double sumG=0.0;
    MPI_Allreduce(&sum, &sumG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
    sum = sumG;
  }

  if(select_ != NULL){
  	StochVector* select2 = dynamic_cast<StochVector*>(select_);
    sum += vec->sumLog(select2->vec);
  }else{
	sum += vec->sumLog(NULL);
  }

  return sum;
}


// sum of the power of elts, can be used to evaluate norm
double StochVector::sumPowElt(const double pow_in)
{
  double sum = 0.0;
  
  for(size_t it=0; it<children.size(); it++)
	sum += children[it]->sumPowElt(pow_in);

  if(iAmDistrib==1) {
    double sumG=0.0;
    MPI_Allreduce(&sum, &sumG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
    sum = sumG;
  }

  sum += vec->sumPowElt(pow_in);

  return sum;
}

// sum of elts
double StochVector::sumElt()
{
  double sum = 0.0;
  
  for(size_t it=0; it<children.size(); it++)
	sum += children[it]->sumElt();

  if(iAmDistrib==1) {
    double sumG=0.0;
    MPI_Allreduce(&sum, &sumG, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
    sum = sumG;
  }

  sum += vec->sumElt();
  
  return sum;
}


void StochVector::MinComponentOrConstant( OoqpVector* vec_in, double minVal  )
{
  StochVector* vec_Input = dynamic_cast<StochVector*>(vec_in);
  assert(children.size() == vec_Input->children.size());

  vec->MinComponentOrConstant(vec_Input->vec,minVal);  

  for(size_t it=0; it<children.size(); it++)
    children[it]->MinComponentOrConstant(vec_Input->children[it],minVal);  

}

void StochVector::correctLargeVal( const double testVal, const double corVal, const int absFlag)
{
  vec->correctLargeVal( testVal, corVal,absFlag);

  for(size_t it=0; it<children.size(); it++)
    children[it]->correctLargeVal( testVal, corVal,absFlag);  

}  
  
void StochVector::MinComponentBetween( OoqpVector* vec_in,OoqpVector *select_in)
{
  StochVector* vecSt = dynamic_cast<StochVector*>(vec_in);
  assert(children.size() == vecSt->children.size());

  if(select_in==NULL){
  	vec->MinComponentBetween( vecSt->vec);
    for(size_t it=0; it<children.size(); it++)
      children[it]->MinComponentBetween( vecSt->children[it]);	
  }else{
  	StochVector* selectSt = dynamic_cast<StochVector*>(select_in);
	assert(children.size() == selectSt->children.size());

	vec->MinComponentBetween( vecSt->vec,selectSt->vec);
  	for(size_t it=0; it<children.size(); it++)
      children[it]->MinComponentBetween( vecSt->children[it],selectSt->children[it]);	
  }
}


// this is used to set slack values, no recursive since we need to do it in upper level
void StochVector::setToConstantFromTo( double c, int start, int end )
{
  assert( "Not supported" && 0 );
}  

void StochVector::SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in,OoqpVector *select_in)
{
  StochVector* vecX = dynamic_cast<StochVector*>(x_in);
  StochVector* vecY = dynamic_cast<StochVector*>(y_in);

  assert(children.size() == vecX->children.size());
  assert(children.size() == vecY->children.size());

  if(select_in==NULL){
  	vec->SetComponentFromMaxXorY(vecX->vec,vecY->vec);  
    for(size_t it=0; it<children.size(); it++)
      children[it]->SetComponentFromMaxXorY(vecX->children[it],vecY->children[it]); 
  }else{
  	StochVector* selectSt = dynamic_cast<StochVector*>(select_in);
	assert(children.size() == selectSt->children.size());

	vec->SetComponentFromMaxXorY(vecX->vec,vecY->vec,selectSt->vec);  
  	for(size_t it=0; it<children.size(); it++)
      children[it]->SetComponentFromMaxXorY(vecX->children[it],vecY->children[it],selectSt->children[it]);	
  }
 
}

void StochVector::SetComponentFromMaxXorY( OoqpVector* x_in, OoqpVector *y_in, 
		  int Start, int End,int xStart, int xEnd, int yStart, int yEnd)
{
	assert( "Not supported" && 0 );
}

void StochVector::SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in)
{
  StochVector* vecX = dynamic_cast<StochVector*>(x_in);
  assert(children.size() == vecX->children.size());

  vec->SetComponentFromMaxXorConstant(vecX->vec,y_in);  

  for(size_t it=0; it<children.size(); it++)
    children[it]->SetComponentFromMaxXorConstant(vecX->children[it],y_in);  

}

void StochVector::SetComponentFromMaxXorConstant( OoqpVector* x_in, const double y_in, 
		int Start, int End,int xStart, int xEnd){assert( "Not supported" && 0 );}

void StochVector::copyFromFromTo( OoqpVector* vec_input, int VStart, int VEnd, int VinStart, int VinEnd)
{assert( "Not supported" && 0 );}



void StochVector::copyFromArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend)
{

  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);

  int wStartTemp = Wstart;
  int wEndTemp   = Wend;

#if 1
  if(iAmDistrib) {
	MPI_Bcast(w, n, MPI_DOUBLE, 0, mpiComm);
  }
  for(size_t it=0; it<children.size(); it++){
    int childLength = children[it]->length();
    children[it]->copyFromArrayFromTo(w,0,childLength,wStartTemp,wEndTemp);
	wStartTemp += childLength;
  }

  int vecLength = vec->length();
  sv.copyIntoArrayFromTo(w,0,vecLength,wStartTemp,wEndTemp);

  if(iAmDistrib) {
  	assert(wStartTemp+vecLength==wEndTemp);  
  }
#endif
}


// only support 2 stage problem, no multi stage
void StochVector::copyIntoArrayFromTo( double *w, int VStart, int VEnd, int Wstart, int Wend)
{  
  SimpleVector& sv  = dynamic_cast<SimpleVector&>(*this->vec);

  int wStartTemp = Wstart;
  int wEndTemp   = Wend;

  if(!parent){
    for( int i = 0; i < n; i++ ) w[i] = 0.0;
  }
  
  for(size_t it=0; it<children.size(); it++){
    int childLength = children[it]->length();
    children[it]->copyIntoArrayFromTo(w,0,childLength,wStartTemp,wEndTemp);
	wStartTemp += childLength;
  }

  if(iAmDistrib) {
  	double *buffer = new double[n]; 
	MPI_Reduce(w, buffer, n, MPI_DOUBLE, MPI_SUM, 0, mpiComm);
	memcpy(w,buffer,n*sizeof(double));
	delete[] buffer;
  }

  int vecLength = vec->length();
  sv.copyIntoArrayFromTo(w,0,vecLength,wStartTemp,wEndTemp);

  if(iAmDistrib) {
  	assert(wStartTemp+vecLength==wEndTemp);  
  }

}



void StochVector::absVal( OoqpVector* vec_in )
{
  if(vec_in==NULL){
  	vec->absVal(NULL);  
    for(size_t it=0; it<children.size(); it++)
      children[it]->absVal(NULL); 
  }else{
  	StochVector* vecInSt = dynamic_cast<StochVector*>(vec_in);
	assert(children.size() == vecInSt->children.size());

	vec->absVal(vecInSt->vec);  
  	for(size_t it=0; it<children.size(); it++)
      children[it]->absVal(vecInSt->children[it]);	
  }
}

// only for 2 stage for aggregation preconditioner
#include "constants.h"
extern PreCondInfo *preCond;
void StochVector::copyIntoArrayWithIndex_AggVarCon( double *ResultArray, const int *vecmap, const int _length, bool isVar)
{
  if(children.size()==0)
  	assert(_length>=n);
  	
  for(size_t it=0; it<children.size(); it++){  	
	assert(_length<=n);
    if(isVar) 
   	 children[it]->copyIntoArrayWithIndex_AggVarCon(ResultArray,preCond->varMap[it],_length,true); 
    else 	 
   	 children[it]->copyIntoArrayWithIndex_AggVarCon(ResultArray,preCond->conMap[it],_length,false); 
  }
  
  if(iAmDistrib) {
    double* buffer = new double[_length];
    MPI_Allreduce(ResultArray, buffer, _length, MPI_DOUBLE, MPI_SUM, mpiComm);
    memcpy(buffer,ResultArray, _length*sizeof(double));	
	delete[] buffer;
  }

  if (vecmap){
    vec->copyIntoArrayWithIndex_AggVarCon(ResultArray,vecmap,_length,isVar);
  }else{
	// root: 1st stage
	if(isVar) 
		vec->copyIntoArrayWithIndex_AggVarCon(ResultArray,preCond->firstVarMap,_length,true);
	else   
		vec->copyIntoArrayWithIndex_AggVarCon(ResultArray,preCond->firstConMap,_length,false);
  } 
}

