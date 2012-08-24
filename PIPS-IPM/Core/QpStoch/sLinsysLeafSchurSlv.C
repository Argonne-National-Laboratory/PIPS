/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeafSchurSlv.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "PardisoSchurSolver.h"
#include "Ma57Solver.h"

/**
 * Computes U = Gi * inv(H_i) * Gi^T.
 *        [ R 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * A and C are the recourse eq. and ineq. matrices, R is the cross
 * Hessian term.
 */
void sLinsysLeafSchurSlv::addTermToDenseSchurCompl( sData *prob, 
						    DenseSymMatrix& SC)
{
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  PardisoSchurSolver* scSolver=dynamic_cast<PardisoSchurSolver*>(solver);

  scSolver->schur_solve(R,A,C, SC);
}
