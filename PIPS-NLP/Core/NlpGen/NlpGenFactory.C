/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGenFactory.h"
#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "myNlpGenSparseLinsys.h"


#include "Ma57Solver.h"
#include "Ma27Solver.h"


#include "SimpleVectorHandle.h"

#include "SimpleVector.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseLinearAlgebraPackage.h"


NlpGenFactory::NlpGenFactory(){}

NlpGenFactory::NlpGenFactory( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_)
  :  nx(nx_), my(my_), mz(mz_),
     nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_)
{}
	


LinearSystem * NlpGenFactory::makeLinsys( Data * prob_in )
{
	NlpGenData * prob = (NlpGenData *) prob_in;
	int n = nx + my + mz;
	
	SparseSymMatrixHandle Mat( new SparseSymMatrix( n,n + nnzQ
							+ nnzA + nnzC ) );
	
	SimpleVectorHandle v( new SimpleVector(n) );
	v->setToZero();
	Mat->setToDiagonal(*v);
	
	prob->putQIntoAt( *Mat, 0, 0 );
	prob->putAIntoAt( *Mat, nx, 0);
	prob->putCIntoAt( *Mat, nx + my, 0 );
	
	Ma57Solver * solver = new Ma57Solver( Mat );


	//QpGenSparseLinsys can be reused for NLP (this was the idea). 
	//However you need to overwrite 'matrixChanged' method to update the nlp linear system 
	//with the new Hessian, Jacobian(s) and all that stuff.
	return new NlpGenSparseLinsys( this, prob, la, Mat, solver );
	

}



Data * NlpGenFactory::makeData()
{
  return new NlpGenData( la, nx, my, mz, nnzQ, nnzA, nnzC );
}




Data         * 
NlpGenFactory::makeData( double    c_[],
			     int    krowQ[],  int  jcolQ[],  double dQ[],
			     double  xlow_[],  char ixlow_[],
			     double  xupp_[],  char ixupp_[],
			     int    krowA[],
			     int    jcolA[],  double dA[],
			     double    b_[],
			     int    krowC[],
			     int    jcolC[],  double dC[],
			     double  clow_[],  char  iclow_[],
			     double  cupp_[],  char  icupp_[] )
{
  // Objective funcition
  SimpleVectorHandle c( new SimpleVector( c_, nx ) );

  nnzQ = krowQ[nx];
  SparseSymMatrixHandle Q( new SparseSymMatrix( nx, nnzQ,
						    krowQ, jcolQ, dQ ) );

  // Bounds on variables
  SimpleVectorHandle xlow( new SimpleVector( xlow_, nx ) );
  SimpleVectorHandle ixlow( new SimpleVector( nx ) );
  ixlow->copyFromArray( ixlow_ );

  SimpleVectorHandle xupp( new SimpleVector( xupp_, nx ) );
  SimpleVectorHandle ixupp( new SimpleVector( nx ) );
  ixupp->copyFromArray( ixupp_ );

  // Equality constraints
  nnzA = 0;
  if( my > 0 )    nnzA = krowA[my];
  SparseGenMatrixHandle A( new SparseGenMatrix( my, nx,
						    nnzA, krowA, jcolA, dA ) );

  SimpleVectorHandle b( new SimpleVector( b_, my ) );

  // Inequality constraints
  nnzC = 0;
  if( mz > 0 )    nnzC = krowC[mz];
  SparseGenMatrixHandle C( new SparseGenMatrix( mz, nx,
						    nnzC, krowC, jcolC, dC ) );

  SimpleVectorHandle clow( new SimpleVector( clow_, mz ) );
  SimpleVectorHandle iclow( new SimpleVector( mz ) );
  iclow->copyFromArray( iclow_ );

  SimpleVectorHandle cupp( new SimpleVector( cupp_, mz ) );
  SimpleVectorHandle icupp( new SimpleVector( mz ) );
  icupp->copyFromArray( icupp_ );

  QpGenData * 
    data = new QpGenData( SparseLinearAlgebraPackage::soleInstance(),
			  c, Q, xlow, ixlow, xupp, ixupp,
			  A, b,
			  C, clow, iclow, cupp, icupp );

  return data;
}


