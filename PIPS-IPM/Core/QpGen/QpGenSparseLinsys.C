/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenSparseLinsys.h"
#include "DoubleLinearSolver.h"
#include "SparseSymMatrix.h"

QpGenSparseLinsys::QpGenSparseLinsys(  QpGen * factory_in,
				       QpGenData * data,
				       LinearAlgebraPackage * la,
				       SparseSymMatrix * Mat_in,
				       DoubleLinearSolver * solver_in )
  : QpGenLinsys( factory_in, data, la ), solver(solver_in)
{
  SpReferTo( Mat, Mat_in );
}


void QpGenSparseLinsys::putXDiagonal( OoqpVector& xdiag )
{
  Mat->atPutDiagonal( 0, xdiag );
}


void QpGenSparseLinsys::putZDiagonal( OoqpVector& zdiag )
{
  Mat->atPutDiagonal( nx + my, zdiag );
  //zdiag.writeToStream(cout);
  //!assert(false);

  /*//!log
  printf("QpGenSparseLinsys::putZDiagonal\n");
  SparseSymMatrix& M = *Mat;

  int* krowM = M.krowM();
  int* jcolM = M.jcolM();
  double* dM = M.M();

  int nn; M.getSize(nn, nn);
  for(int i=0; i<nn; i++) {
    printf("row %d \n\t", i);

    for(int j=krowM[i];j<krowM[i+1]; j++)
      printf("%9d ", jcolM[j]);
    printf("\n\t");

    for(int j=krowM[i];j<krowM[i+1]; j++)
      printf("%9.6f ", dM[j]);
    printf("\n");    
  }
  */
}


void QpGenSparseLinsys::solveCompressed( OoqpVector & arhs )
{
  //printf("-----\n");arhs.writeToStream(cout);
  solver->solve( arhs );
  //printf("~~~~~\n");arhs.writeToStream(cout);
}  


void QpGenSparseLinsys::factor(Data *prob, Variables *vars)
{
  this->QpGenLinsys::factor( prob, vars );
  solver->matrixChanged();
}

QpGenSparseLinsys::~QpGenSparseLinsys()
{
  delete solver;
}
