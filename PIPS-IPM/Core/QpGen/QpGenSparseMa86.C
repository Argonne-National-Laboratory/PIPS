/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenSparseMa86.h"
#include "QpGenSparseLinsys.h"
#include "QpGenData.h"

#include "SparseLinearAlgebraPackage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "Ma86Solver.h"
#include "SparseLinearAlgebraPackage.h"


QpGenSparseMa86::QpGenSparseMa86( int nx_in, int my_in, int mz_in,
					  int nnzQ_in, int nnzA_in, int nnzC_in ) :
  QpGenSparseSeq( nx_in, my_in, mz_in, nnzQ_in, nnzA_in, nnzC_in )
{
  la = SparseLinearAlgebraPackage::soleInstance();
}

LinearSystem * QpGenSparseMa86::makeLinsys( Data * prob_in )
{
  QpGenData * prob = (QpGenData *) prob_in;
  int n = nx + my + mz;

  SparseSymMatrixHandle Mat( new SparseSymMatrix( n,n + nnzQ
						      + nnzA + nnzC ) );

  SimpleVectorHandle v( new SimpleVector(n) );
  v->setToZero();
  Mat->setToDiagonal(*v);

  prob->putQIntoAt( *Mat, 0, 0 );
  prob->putAIntoAt( *Mat, nx, 0);
  prob->putCIntoAt( *Mat, nx + my, 0 );
  
  Ma86Solver * solver = new Ma86Solver( Mat );
   
  return new QpGenSparseLinsys( this, prob, la, Mat, solver );
}
