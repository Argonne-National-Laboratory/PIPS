/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "QpGenScaLinsys.h"
#include "QpGenData.h"
#include "ScaDenSymMatrix.h"
#include "QpGenSca.h"
#include "ScaGenIndefSolver.h"
#include "LinearAlgebraPackage.h"

QpGenScaLinsys::QpGenScaLinsys( QpGen * qpgen,
			      QpGenData * data,
			      LinearAlgebraPackage * la, ScaDenSymMatrix * Mat,
			      ScaGenIndefSolver * solver_in ) :
  QpGenLinsys( qpgen, data, la ), solver(solver_in)
{
  SpReferTo( kkt, Mat );
  QpGenSca *qp = dynamic_cast<QpGenSca*>(qpgen);
  delete rhs;
  rhs = qp->sca_la->newVector(rhs->length());

}

void QpGenScaLinsys::factor(Data *prob_in, Variables *vars)
{
  QpGenData * prob = (QpGenData *) prob_in;
  

  kkt->scalarMult(0.);
  prob->putQIntoAt( *kkt, 0,       0 );
  prob->putAIntoAt( *kkt, nx,      0 );
  prob->putCIntoAt( *kkt, nx + my, 0 );
  
  //kkt->symAtPutZeros( nx, nx, my+mz, my+mz );
  /*
  cout << nx << " " << my+mz << endl;
  cout << "A\n";
  prob->A->writeToStream(cout);
  cout << "END A\n";
  */
  QpGenLinsys::factor( prob, vars );
  

  //kkt->writeToStream(cout);
  solver->matrixChanged();

}

void QpGenScaLinsys::putXDiagonal( OoqpVector& xdiag )
{
  kkt->atPutDiagonal( 0, xdiag );
}


void QpGenScaLinsys::putZDiagonal( OoqpVector& zdiag )
{
  kkt->atPutDiagonal( nx + my, zdiag );
}


void QpGenScaLinsys::solveCompressed( OoqpVector & compressedRhs )
{
  solver->solve( compressedRhs );
}  


QpGenScaLinsys::~QpGenScaLinsys()
{
  delete solver;
}
