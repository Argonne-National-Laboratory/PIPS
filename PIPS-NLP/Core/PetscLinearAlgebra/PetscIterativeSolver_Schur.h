/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PETSCITERATIVESOLVER_H
#define PETSCITERATIVESOLVER_H

#include "DoubleLinearSolver.h"
#include "DeSymIndefSolver.h"

#include "petscksp.h"
#include "SimpleVector.h"

class PetscIterativeSolver_Schur : public DeSymIndefSolver {
protected:
  KSP mKsp;
  int deleteKSP;
  int total_kry_iter;
  
  Vec x;  
  Mat kspMat;
  PC  precondMat;

  SimpleVectorHandle rhs_back;

public:
	
  PetscIterativeSolver_Schur( DenseSymMatrix * SC_in);
  //, KSPType *ksptype, PCType *pctype=NULL);
  
//  PetscIterativeSolver( KSP ksp, PC pc );

  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve ( OoqpVector&  x );
  
  KSP ksp() { return mKsp; }; 
  virtual ~PetscIterativeSolver_Schur();

};

#endif
