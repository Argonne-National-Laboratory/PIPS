/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

/*
 * Thanks to Quan H. Nguyen for help porting to newer versions of Petsc.
 */

#include "PetscIterativeSolver_Schur.h"
//#include "PetscSpSymMatrix.h"
//#include "PetscSparseStorage.h"
#include "petscksp.h"
//#include "PetscVector.h"
//#include "PetscVectorHandle.h"

#include "DenseSymMatrix.h"
#include "SparseSymMatrix.h"

#include "SimpleVector.h"

double gKryIterTol=1e-8;

typedef struct {
  DoubleLinearSolver * lin_solver;
  SparseSymMatrix* matS;
  DenseSymMatrix* matD;
  SymMatrix* mat;
} user_PC;


extern PetscErrorCode _user_MatMult(Mat,Vec,Vec);
extern PetscErrorCode _user_PC_Create(user_PC**);
extern PetscErrorCode _user_PC_SetUp(PC pc, DoubleLinearSolver* la, SymMatrix* K);
extern PetscErrorCode _user_PC_Destroy(PC pc);
extern PetscErrorCode _user_PC_Apply_Schur(PC pc,Vec x,Vec y);


PetscIterativeSolver_Schur::PetscIterativeSolver_Schur( DenseSymMatrix * SC_in) :  
	deleteKSP(1), total_kry_iter(0), DeSymIndefSolver(SC_in)
{ 
  int ierr;
  int dummy, n;
  
  int nb_row, nb_col;
  SC_in->getSize(nb_row,nb_col);

  rhs_back = new SimpleVector(nb_row);
  rhs_back->setToZero();

  // ********  Set link to user defined Matrix type ********
  MatCreateShell(PETSC_COMM_SELF, nb_row, nb_col, PETSC_DETERMINE, PETSC_DETERMINE, (void*) SC_in, &kspMat);
  MatShellSetOperation(kspMat, MATOP_MULT, (void (*)(void))_user_MatMult);

  // ********  Create linear solver context ********
  ierr = KSPCreate(PETSC_COMM_SELF,&mKsp);assert(ierr == 0);
	
  // ********  Set operators. Here the matrix that defines the linear system also serves as the preconditioning matrix. ********
  ierr = KSPSetOperators(mKsp, kspMat, kspMat); assert(ierr == 0);

  // ********  Set linear solver defaults for this problem (optional).
  // ********  - By extracting the KSP and PC contexts from the KSP context,
  // ********  we can then directly call any KSP and PC routines to set
  // ********  various options.
  // ********  - The following four statements are optional; all of these
  // ********  parameters could alternatively be specified at runtime via
  // ********  KSPSetFromOptions();
  
//  ierr = KSPSetType(mKsp, ksptype);assert(ierr == 0);

  ierr = KSPGetPC(mKsp,&precondMat);assert(ierr == 0);

  if (guser_defined_pc == 1) 
  {
    PCSetType(precondMat,PCSHELL);
		
//    user_PC* shell_la;
//    _user_PC_Create(&shell_la);
//    PCShellSetApply(precondMat,_user_PC_Apply_Schur);
//    PCShellSetContext(precondMat,shell_la);
//    PCShellSetDestroy(precondMat,_user_PC_Destroy);
//    _user_PC_SetUp(precondMat,this,SC_in);
  }
  else 
  {
	PCSetType(precondMat,PCJACOBI);
  }

  ierr = KSPSetFromOptions(mKsp);assert(ierr == 0);
  ierr = KSPSetTolerances(mKsp,gKryIterTol,gKryIterTol,PETSC_DEFAULT,PETSC_DEFAULT); assert(ierr == 0);

}


PetscIterativeSolver_Schur::~PetscIterativeSolver_Schur()
{
  int ierr;

  if( deleteKSP ) { // We made it, we own it.
    ierr = KSPDestroy( &mKsp ); assert( ierr  == 0);
  }
}

void PetscIterativeSolver_Schur::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

int PetscIterativeSolver_Schur::matrixChanged()
{
  int ierr;

  ierr = KSPSetOperators(mKsp, kspMat, kspMat); assert(ierr == 0);
}



void PetscIterativeSolver_Schur::solve( OoqpVector& rhs_in)
{
  int ierr, its;
  int nb_row, nb_col;
  mStorage->getSize(nb_row,nb_col);


  SimpleVector &  x_sol = dynamic_cast<SimpleVector &>(rhs_in);
  rhs_back->copyFrom(x_sol);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //   Solve the linear system Ax=b by petsc
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  Vec b_la,x_la;
  VecCreateSeqWithArray(PETSC_COMM_WORLD,1,nb_col,rhs_back->elements(),&b_la);
  VecAssemblyBegin(b_la);
  VecAssemblyEnd(b_la);
  VecCreateSeqWithArray(PETSC_COMM_WORLD,1,nb_col,&x_sol[0],&x_la);
  VecAssemblyBegin(x_la);
  VecAssemblyEnd(x_la);
  ierr = KSPSolve(mKsp,b_la,x_la); assert( ierr  == 0);
  VecDestroy(&b_la);
  VecDestroy(&x_la);  
  
  ierr = KSPGetIterationNumber(mKsp, &its); assert( ierr  == 0);

  total_kry_iter += its;  

}





// tell petsc how to compute Ax=b
PetscErrorCode _user_MatMult(Mat A, Vec x,  Vec y)
{
  DenseSymMatrix   *Amat;
  MatShellGetContext(A, &Amat);

  int nb_row, nb_col;
  Amat->getSize(nb_row,nb_col);

  double* x_array, *y_array;
  VecGetArray(x,&x_array);
  VecGetArray(y,&y_array);

  SimpleVector x_temp(x_array,nb_col);
  SimpleVector y_temp(y_array,nb_col);

  Amat->mult(0.0,x_temp,1.0,y_temp);
  
  VecRestoreArray(x,&x_array);
  VecRestoreArray(y,&y_array);

  return 0;
}

PetscErrorCode _user_PC_Create(user_PC **shell)
{
  (*shell)= new user_PC;
  return 0;
}

PetscErrorCode _user_PC_SetUp(PC pc, DoubleLinearSolver* la, SymMatrix* K)
{
  user_PC *shell;
  PCShellGetContext(pc,(void**)&shell);
  shell->lin_solver = la;
  shell->mat = K;
  return 0;
}

PetscErrorCode _user_PC_Destroy(PC pc)
{
  user_PC  *shell;
  PCShellGetContext(pc,(void**)&shell);
  delete shell;
  return 0;
}

// solve Ax=b, where A is preconditioner
PetscErrorCode _user_PC_Apply_Schur(PC pc,Vec x,Vec y)
{ 
  user_PC *shell;
  PCShellGetContext(pc,(void**)&shell);
  
  DenseSymMatrix *Amat = dynamic_cast<DenseSymMatrix*>(shell->mat);

  int nb_row, nb_col;
  Amat->getSize(nb_row,nb_col);

  double* x_array, *y_array;
  VecGetArray(x,&x_array);
  VecGetArray(y,&y_array);

  SimpleVector x_temp(x_array,nb_col);
  SimpleVector y_temp(y_array,nb_col);
  y_temp.copyFrom(x_temp);

  PetscIterativeSolver_Schur *obj_solver = dynamic_cast<PetscIterativeSolver_Schur *>(shell->lin_solver);
  obj_solver->DeSymIndefSolver::solve(y_temp);

  VecRestoreArray(x,&x_array);  
  VecRestoreArray(y,&y_array);
  return 0;
}

