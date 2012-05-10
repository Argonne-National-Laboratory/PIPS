/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenSca.h"
#include "QpGenData.h"
#include "QpGenScaLinsys.h"
#include "SimpleVector.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseLinearAlgebraPackage.h"
#include "ScaLinearAlgebraPackage.h"
#include "QpGenVars.h"
#include "ScaVector.h"
#include "ScaGenIndefSolver.h"


QpGenSca::QpGenSca(int nx_, int my_, int mz_,
		  int nnzQ_, int nnzA_, int nnzC_,
		  COMMINFO &cinfo)
    : QpGen( nx_, my_, mz_ ),
      nnzQ(nnzQ_), nnzA(nnzA_), nnzC(nnzC_)
{
  la = SparseLinearAlgebraPackage::soleInstance();
  // should delete this somewhere
  sca_la = new ScaLinearAlgebraPackage(cinfo);
}
	

Data * QpGenSca::makeData()
{
  return new QpGenData( la, nx, my, mz, nnzQ, nnzA, nnzC );
}

void QpGenSca::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
				   OoqpVector& rhs2_in, OoqpVector& rhs3_in )
{
  ScaVector & rhs  = dynamic_cast<ScaVector &>(rhs_in);
  SimpleVector & rhs1 = dynamic_cast<SimpleVector &>(rhs1_in);
  SimpleVector & rhs2 = dynamic_cast<SimpleVector &>(rhs2_in);
  SimpleVector & rhs3 = dynamic_cast<SimpleVector &>(rhs3_in);

  rhs.joinRHS(rhs1, nx, rhs2, my, rhs3, mz);
}

void
QpGenSca::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				   OoqpVector& z_in, OoqpVector& vars_in )
{
  ScaVector & vars  = dynamic_cast<ScaVector &>(vars_in);
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  SimpleVector & z = dynamic_cast<SimpleVector &>(z_in);

  vars.separateVars(x, nx, y, my, z, mz);
}


Data         * 
QpGenSca::makeData( double    c_[],
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

LinearSystem  *QpGenSca::makeLinsys( Data * prob_in )
{
  QpGenData * prob = (QpGenData *) prob_in;

  int n = nx + my + mz;
  ScaDenSymMatrixHandle Mat((ScaDenSymMatrix*)sca_la->newSymMatrix(n,0));

  ScaGenIndefSolver * solver = new ScaGenIndefSolver( Mat );

  return new QpGenScaLinsys( this, prob, la, Mat, solver );
}


/*
Data   *
QpGenSca::
copyDataFromSparseTriple( double c[],
			  int irowQ[], int lenQ,  int jcolQ[],  double dQ[],
			  double xlow[],  char ixlow[],
			  double xupp[],  char ixupp[],
			  int irowA[], int lenA,  int jcolA[],  double dA[],
			  double   bA[],
			  int irowC[],  int lenC,  int jcolC[], double dC[],
			  double clow[], char iclow[],
			  double cupp[], char icupp[] )
{
  QpGenData * prob =
    (QpGenData *) new QpGenData( la, nx, my, mz, nnzQ, nnzA, nnzC );  
  int info;

  assert( lenQ <= nnzQ );
  assert( lenA <= nnzA );
  assert( lenC <= nnzC );

  prob->g->copyFromArray( c );
  prob->Q->putSparseTriple( irowQ, lenQ, jcolQ, dQ, info );

  prob-> blx->copyFromArray( xlow );
  prob->ixlow->copyFromArray( ixlow );

  prob-> bux->copyFromArray( xupp );
  prob->ixupp->copyFromArray( ixupp );
  
  prob->A->putSparseTriple( irowA, lenA, jcolA, dA, info );
  prob->bA->copyFromArray( bA );
  
  prob->C->putSparseTriple( irowC, lenC, jcolC, dC, info );

  prob->bl   ->copyFromArray(  clow );
  prob->iclow->copyFromArray( iclow );
  
  prob->bu   ->copyFromArray(  cupp );
  prob->icupp->copyFromArray( icupp );

  return prob;
}*/
