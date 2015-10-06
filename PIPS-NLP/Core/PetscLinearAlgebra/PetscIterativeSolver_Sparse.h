/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef PETSCITERATIVESOLVER_H
#define PETSCITERATIVESOLVER_H

#include "DoubleLinearSolver.h"
//#include "DeSymIndefSolver.h"

#include "petscksp.h"
#include "SimpleVector.h"

class SparseGenMatrix;
class SparseSymMatrix;

class PetscIterativeSolver_Sparse : public DoubleIterativeLinearSolver{
protected:
  SparseGenMatrix *fullSymMat;
  SparseGenMatrix *PCGenMat;
  SparseSymMatrix *PCSymMat;
  SparseSymMatrix *inputMatptr;
  int *goffIDX, *PCgoffIDX;
  int correct_negEigVal;
	
  KSP mKsp;
  int deleteKSP;
  int total_kry_iter;

  DoubleLinearSolver * linear_solver;
  
//  Vec x;  
  Mat LinSysMat_PETSC;
  Mat PCMat_PETSC;
  PC  precond_Method;

  SimpleVectorHandle rhs_back;

public:
	
  PetscIterativeSolver_Sparse( SparseSymMatrix * SC_in, const int numOfNegEigVal_in);
  //, KSPType *ksptype, PCType *pctype=NULL);
  
//  PetscIterativeSolver( KSP ksp, PC pc );

  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve ( OoqpVector&  x );
  
  KSP ksp() { return mKsp; }; 
  virtual ~PetscIterativeSolver_Sparse();

};

#endif

