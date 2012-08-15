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

  //cout  << " sLinsysLeafSchurSlv::addTermToDenseSchurComp is apparently called!" << endl;
  SparseGenMatrix& A = prob->getLocalA();
  SparseGenMatrix& C = prob->getLocalC();
  SparseGenMatrix& R = prob->getLocalCrossHessian();

  PardisoSchurSolver* scSolver=dynamic_cast<PardisoSchurSolver*>(solver);

  //////////////////////////////////////////////////////////
  // int nR, mR; R.getSize(nR,mR);
  // int nA, mA; A.getSize(nA,mA);
  // int nC, mC; C.getSize(nC,mC);

  // int nr=10;  

  // SparseGenMatrix AA(nA,nr,0);
  // SparseGenMatrix CC(nC,nr,nr);
  // SparseGenMatrix RR(nR,nr,0);

  // int* krowC    = CC.getStorageRef().krowM;
  // int* jcolC    = CC.getStorageRef().jcolM;
  // double* eltsC = CC.getStorageRef().M;

  
  // for(int k=0; k<=nr; k++) krowC[k]=k; 
  
  // for(int r=nr+1; r<=nC; r++) 
  //   krowC[r]=krowC[nr];

  // for(int r=0; r<nr; r++) {
  //   jcolC[r]=r;
  //   eltsC[r]=2.5 + r/100.1;
  // } 


  ///////////////////////////////////////////////////////////
  // int nxP=SC.size();
  // //scSolver->schur_solve(RR,AA,CC, SC);
   scSolver->schur_solve(R,A,C, SC);
  // //cout << " sLinsysLeafSchurSlv::addTermToDenseSchurComp done" << endl;
  // cout << "Printing Pardiso-SCUR COMPLEMENT" << endl;
  // int lim=min(8,nxP);
  // for(int i=0; i<lim;i++) {
  //   for(int j=0; j<lim; j++)
  //     printf("%18.8f ", SC[i][j]);
  //   cout<<endl;
  // }

}
