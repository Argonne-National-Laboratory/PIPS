/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeafSchurSlv.h"
#include "sTree.h"
#include "sFactory.h"
#include "sData.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "PardisoSolver.h"
#include "PardisoSchurSolver.h"
#include "Ma57Solver.h"

extern int gLackOfAccuracy;

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

  
  if(!gLackOfAccuracy && !switchedToSafeSlv) {
    PardisoSchurSolver* scSolver=dynamic_cast<PardisoSchurSolver*>(solver);
    scSolver->schur_solve(R,A,C, SC);
  } else {
    cout << "\tdefaulting to sLinsysLeaf::addTermToDenseSchurCompl ...";
    sLinsysLeaf::addTermToDenseSchurCompl(prob, SC);
    cout << "done" << endl;
  }
}

void sLinsysLeafSchurSlv::factor2(sData *prob, Variables *vars)
{
  if(gLackOfAccuracy) {
    cout << "sLinsysLeafSchurSlv -> accuracy lost, switching to vanilla PARDISO" << endl;
    delete solver;
    cout << "\tsolver deleted\n";
    SparseSymMatrix* kktsp = dynamic_cast<SparseSymMatrix*>(kkt);
    solver = new PardisoSolver(kktsp);
    cout << "\tnew solver created." << endl;
    switchedToSafeSlv=true;
  }
  
  sLinsysLeaf::factor2(prob, vars);
}
