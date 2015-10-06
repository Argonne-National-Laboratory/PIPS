/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
 
 /* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "NlpGenSparse.h"
#include "NlpGenData.h"
#include "SimpleVector.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseLinearAlgebraPackage.h"
#include "NlpGenVars.h"
//#include "NlpInfoAMPL.h"
#include "NlpInfo.h"


Data * NlpGenSparse::makeData()
{
  return new NlpGenData( la, nx, my, mz, nnzQ, nnzA, nnzC );
}

void NlpGenSparse::joinRHS( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
				   OoqpVector& rhs2_in, OoqpVector& rhs3_in)
{
  SimpleVector & rhs  = dynamic_cast<SimpleVector &>(rhs_in);
  SimpleVector & rhs1 = dynamic_cast<SimpleVector &>(rhs1_in);
  SimpleVector & rhs2 = dynamic_cast<SimpleVector &>(rhs2_in);
  SimpleVector & rhs3 = dynamic_cast<SimpleVector &>(rhs3_in);

  memcpy( &rhs[0], &rhs1[0], nx * sizeof( double ) );
  if( my > 0 ) memcpy( &rhs[nx],      &rhs2[0], my * sizeof( double ) );
  if( mz > 0 ) memcpy( &rhs[nx + my], &rhs3[0], mz * sizeof( double ) );
}


void
NlpGenSparse::separateVars( OoqpVector& x_in, OoqpVector& y_in,
				   OoqpVector& z_in, OoqpVector& vars_in)
{
  SimpleVector & vars  = dynamic_cast<SimpleVector &>(vars_in);
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  SimpleVector & z = dynamic_cast<SimpleVector &>(z_in);

  memcpy( &x[0], &vars[0], nx * sizeof( double ) );
  if ( my > 0 ) memcpy( &y[0], &vars[nx],      my * sizeof( double ) );
  if ( mz > 0 ) memcpy( &z[0], &vars[nx + my], mz * sizeof( double ) );
}


Data         * 
NlpGenSparse::makeData( double    c_[],
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
  SimpleVectorHandle grad( new SimpleVector( c_, nx ) );

  nnzQ = krowQ[nx];
  SparseSymMatrixHandle H( new SparseSymMatrix( nx, nnzQ,
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
  SparseGenMatrixHandle Jeq( new SparseGenMatrix( my, nx,
						    nnzA, krowA, jcolA, dA ) );

  SimpleVectorHandle b( new SimpleVector( b_, my ) );

  // Inequality constraints
  nnzC = 0;
  if( mz > 0 )    nnzC = krowC[mz];
  SparseGenMatrixHandle Jineq( new SparseGenMatrix( mz, nx,
						    nnzC, krowC, jcolC, dC ) );

  SimpleVectorHandle clow( new SimpleVector( clow_, mz ) );
  SimpleVectorHandle iclow( new SimpleVector( mz ) );
  iclow->copyFromArray( iclow_ );

  SimpleVectorHandle cupp( new SimpleVector( cupp_, mz ) );
  SimpleVectorHandle icupp( new SimpleVector( mz ) );
  icupp->copyFromArray( icupp_ );

  NlpGenData * 
    data = new NlpGenData( SparseLinearAlgebraPackage::soleInstance(),
			  grad, H, xlow, ixlow, xupp, ixupp,
			  Jeq, b,
			  Jineq, clow, iclow, cupp, icupp );

  return data;
}


Data   *
NlpGenSparse::copyDataFromSparseTriple( double c[],
			  int irowQ[], int lenQ,  int jcolQ[],  double dQ[],
			  double xlow[],  char ixlow[],
			  double xupp[],  char ixupp[],
			  int irowA[], int lenA,  int jcolA[],  double dA[],
			  double   bA[],
			  int irowC[],  int lenC,  int jcolC[], double dC[],
			  double clow[], char iclow[],
			  double cupp[], char icupp[], int *rowMap, NlpInfo * updateNlp)
{
  NlpGenData * prob =
    (NlpGenData *) new NlpGenData( la, nx, my, mz, nnzQ, nnzA, nnzC );  
  int info;

  assert( lenQ <= nnzQ );
  assert( lenA <= nnzA );
  assert( lenC <= nnzC );

  prob->grad->copyFromArray( c );
  prob->H->putSparseTriple( irowQ, lenQ, jcolQ, dQ, info );

  prob->blx->copyFromArray( xlow );
  prob->ixlow->copyFromArray( ixlow );

  prob->bux->copyFromArray( xupp );
  prob->ixupp->copyFromArray( ixupp );
  
  prob->Jeq->putSparseTriple( irowA, lenA, jcolA, dA, info );
  prob->bA->copyFromArray( bA );
  
  prob->Jineq->putSparseTriple( irowC, lenC, jcolC, dC, info );

  prob->bl   ->copyFromArray(  clow );
  prob->iclow->copyFromArray( iclow );
  
  prob->bu   ->copyFromArray(  cupp );
  prob->icupp->copyFromArray( icupp );

  if(rowMap)
  	updateNlp->rowMap=rowMap;

  prob->inputNlp=updateNlp;

  return prob;
}




Data   *
NlpGenSparse::copyDataFromSparseTriple( double c[],
			  int irowQ[], int lenQ,  int jcolQ[],  double dQ[],
			  double xlow[],  char ixlow[],
			  double xupp[],  char ixupp[],
			  int irowA[], int lenA,  int jcolA[],  double dA[],
			  double   bA[],
			  int irowC[],  int lenC,  int jcolC[], double dC[],
			  double clow[], char iclow[],
			  double cupp[], char icupp[], int *rowMap,
			  int nxL,int nxU,int nsL,int nsU, NlpInfo * updateNlp)
{
  NlpGenData * prob =
    (NlpGenData *) new NlpGenData( la, nx, my, mz, nnzQ, nnzA, nnzC, nxL, nxU, nsL, nsU );  
  int info;

  assert( lenQ <= nnzQ );
  assert( lenA <= nnzA );
  assert( lenC <= nnzC );

  prob->grad->copyFromArray( c );
  prob->H->putSparseTriple( irowQ, lenQ, jcolQ, dQ, info );

  prob->blx->copyFromArray( xlow );
  prob->ixlow->copyFromArray( ixlow );

  prob->bux->copyFromArray( xupp );
  prob->ixupp->copyFromArray( ixupp );
  
  prob->Jeq->putSparseTriple( irowA, lenA, jcolA, dA, info );
  prob->bA->copyFromArray( bA );
  
  prob->Jineq->putSparseTriple( irowC, lenC, jcolC, dC, info );

  prob->bl   ->copyFromArray(  clow );
  prob->iclow->copyFromArray( iclow );
  
  prob->bu   ->copyFromArray(  cupp );
  prob->icupp->copyFromArray( icupp );


  

  if(rowMap)
  	updateNlp->rowMap=rowMap;

  prob->inputNlp=updateNlp;


  return prob;
}


void NlpGenSparse::joinRHSXSYZ( OoqpVector& rhs_in,  OoqpVector& rhs1_in,
			     OoqpVector& rhs2_in, OoqpVector& rhs3_in, OoqpVector& rhs4_in )
{
  SimpleVector & rhs  = dynamic_cast<SimpleVector &>(rhs_in);
  SimpleVector & rhs1 = dynamic_cast<SimpleVector &>(rhs1_in);
  SimpleVector & rhs2 = dynamic_cast<SimpleVector &>(rhs2_in);
  SimpleVector & rhs3 = dynamic_cast<SimpleVector &>(rhs3_in);
  SimpleVector & rhs4 = dynamic_cast<SimpleVector &>(rhs4_in);

  memcpy( &rhs[0], &rhs1[0], nx * sizeof( double ) );
  if( mz > 0 ) memcpy( &rhs[nx], &rhs2[0], mz * sizeof( double ) );
  if( my > 0 ) memcpy( &rhs[nx + mz],      &rhs3[0], my * sizeof( double ) );
  if( mz > 0 ) memcpy( &rhs[nx + mz + my], &rhs4[0], mz * sizeof( double ) );


}

void NlpGenSparse::separateVarsXSYZ( OoqpVector& x_in, OoqpVector& s_in,
				  OoqpVector& y_in,OoqpVector& z_in, OoqpVector& vars_in )
{
  SimpleVector & vars  = dynamic_cast<SimpleVector &>(vars_in);
  SimpleVector & x = dynamic_cast<SimpleVector &>(x_in);
  SimpleVector & s = dynamic_cast<SimpleVector &>(s_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  SimpleVector & z = dynamic_cast<SimpleVector &>(z_in);

  memcpy( &x[0], &vars[0], nx * sizeof( double ) );
  if ( mz > 0 ) memcpy( &s[0], &vars[nx], mz * sizeof( double ) );
  if ( my > 0 ) memcpy( &y[0], &vars[nx + mz],      my * sizeof( double ) );
  if ( mz > 0 ) memcpy( &z[0], &vars[nx + mz + my], mz * sizeof( double ) );
}


void NlpGenSparse::copyXSYZ_fromArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col)
{
  SimpleVector & vars  = dynamic_cast<SimpleVector &>(vec_xsyz);
  memcpy( &vars[0], &array_in[0], nb_col * sizeof( double ) );
}

void NlpGenSparse::copyXSYZ_toArray( OoqpVector& vec_xsyz, double* array_in, const int nb_col)
{
  SimpleVector & vars  = dynamic_cast<SimpleVector &>(vec_xsyz);
  memcpy( &array_in[0], &vars[0], nb_col * sizeof( double ) );
}

