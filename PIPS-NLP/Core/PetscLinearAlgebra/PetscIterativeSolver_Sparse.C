/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "PetscIterativeSolver_Sparse.h"
#include "petscksp.h"

#include "DenseSymMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"


#include "Ma57Solver.h"


#include "SimpleVector.h"

//char gpetsc_preconditioner;

extern int gSymLinearSolver;
extern int gUser_Defined_PC;
extern int gUser_Defined_SymMat;

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
extern PetscErrorCode _user_PC_Apply_Sparse(PC pc,Vec b_in,Vec x_out);


#ifdef TIMING
  static double tTot=0, tMine=0;
#endif


PetscIterativeSolver_Sparse::PetscIterativeSolver_Sparse( SparseSymMatrix * mat_in, const int numOfNegEigVal_in) :  
	deleteKSP(1), total_kry_iter(0), inputMatptr(mat_in), fullSymMat(NULL), goffIDX(NULL), PCgoffIDX(NULL)
{  
  correct_negEigVal = numOfNegEigVal_in;
  linear_solver = NULL;

  int ierr;
  int dummy, n;
  
  int nb_row, nb_col;
  mat_in->getSize(nb_row,nb_col);
  int nnz_in = mat_in->numberOfNonZeros();

  rhs_back = new SimpleVector(nb_row);
  rhs_back->setToZero();

  // ********  Set link to user defined Matrix type ********
  fullSymMat = NULL;
  if (gUser_Defined_SymMat == 1) 
  {
    MatCreateShell(PETSC_COMM_SELF, nb_row, nb_col, PETSC_DETERMINE, PETSC_DETERMINE, (void*) mat_in, &LinSysMat_PETSC);
    MatShellSetOperation(LinSysMat_PETSC, MATOP_MULT, (void (*)(void))_user_MatMult);
  }else{
    int info_;

    fullSymMat = new SparseGenMatrix( nb_row, nb_col, nnz_in*2-nb_row);
    memcpy(&(fullSymMat->krowM()[0]),&(mat_in->krowM()[0]),(nb_row+1)* sizeof( int ));
    memcpy(&(fullSymMat->jcolM()[0]),&(mat_in->jcolM()[0]),(nnz_in)* sizeof( int ));
    memcpy(&(fullSymMat->M()[0]),&(mat_in->M()[0]),(nnz_in)* sizeof( double ));
    assert(goffIDX==NULL);
    goffIDX =  fullSymMat->symmetrize_set(info_);
	
    MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,nb_row,nb_col,fullSymMat->krowM(),fullSymMat->jcolM(),fullSymMat->M(),&LinSysMat_PETSC);
  }

  if(gUser_Defined_PC == 0){
	// use petsc defined pc, PCMat_PETSC must be presented in petsc form. Note that here gUser_Defined_SymMat must be 0
    PCMat_PETSC = LinSysMat_PETSC;
  }else if(gUser_Defined_PC == 1){
    // use user defined PC, here we set pc as the iterate matrix, and latter we do factorization on it, i.e. this is equal to the direct solve
    PCMat_PETSC = LinSysMat_PETSC;
	PCSymMat = mat_in; 
  }else if(gUser_Defined_PC == 2){
    // use user defined PC, here we use Jacek's idea to set preconditioner, i.e. onlt use the diagonals of H 

    int info_;
    int nb_Qrow = nb_row-correct_negEigVal;
    int nnz_Qpart = mat_in->krowM()[nb_Qrow] ;
    int nnz_Apart = nnz_in-nnz_Qpart;
	int nnz_PC    = nnz_Apart+nb_Qrow;

	assert(PCgoffIDX==NULL);
	PCgoffIDX = new int[nnz_PC];
	
    int find_nnz=0;
    PCSymMat = new SparseSymMatrix( nb_row, nnz_PC);

	// copy the diagonal part of Q
    for(int row_i=0; row_i< nb_Qrow; row_i++){
      PCSymMat->krowM()[row_i] = find_nnz;
      PCSymMat->jcolM()[find_nnz] = row_i;
	  PCgoffIDX[find_nnz] = mat_in->krowM()[row_i+1]-1;
	  PCSymMat->M()[find_nnz] = mat_in->M()[PCgoffIDX[find_nnz]];
      find_nnz++;
    }

	// copy the part of A
	int AgoffStart = mat_in->krowM()[nb_Qrow];
    memcpy(&(PCSymMat->krowM()[nb_Qrow]),&(mat_in->krowM()[nb_Qrow]),(nb_row-nb_Qrow+1)* sizeof( int ));
    for(int row_i=nb_Qrow; row_i< nb_row+1; row_i++){
      PCSymMat->krowM()[row_i] -= nnz_Qpart - nb_Qrow;
    }
	
    memcpy(&(PCSymMat->jcolM()[find_nnz]),&(mat_in->jcolM()[AgoffStart]),(nnz_Apart)* sizeof( int ));

	memcpy(&(PCSymMat->M()[find_nnz]),&(mat_in->M()[AgoffStart]),(nnz_Apart)* sizeof( double ));	

    PCgoffIDX[find_nnz] = AgoffStart;
	for(find_nnz=find_nnz+1; find_nnz<nnz_PC; find_nnz++){
      PCgoffIDX[find_nnz] = PCgoffIDX[find_nnz-1]+1;
    }
	
	// no need to buid pc in petsc form, but we may need it in the future
    // MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,nb_row,nb_col,PCSymMat->krowM(),PCSymMat->jcolM(),PCSymMat->M(),&PCMat_PETSC);
  }

  

  // ********  Create linear solver context ********
  ierr = KSPCreate(PETSC_COMM_SELF,&mKsp);assert(ierr == 0);
	
  // ********  Set operators. Here the matrix that defines the linear system also serves as the preconditioning matrix. ********
  ierr = KSPSetOperators(mKsp, LinSysMat_PETSC, LinSysMat_PETSC); assert(ierr == 0);

  // ********  Set linear solver defaults for this problem (optional).
  // ********  - By extracting the KSP and PC contexts from the KSP context,
  // ********  we can then directly call any KSP and PC routines to set
  // ********  various options.
  // ********  - The following four statements are optional; all of these
  // ********  parameters could alternatively be specified at runtime via
  // ********  KSPSetFromOptions();
  
//  ierr = KSPSetType(mKsp, ksptype);assert(ierr == 0);

  ierr = KSPGetPC(mKsp,&precond_Method);assert(ierr == 0);

  if (gUser_Defined_PC != 0) 
  {
    if(1==gSymLinearSolver)
      linear_solver = new Ma57Solver(PCSymMat);

    ierr = PCSetType(precond_Method,PCSHELL);assert(ierr == 0);
    user_PC* shell_la;
    _user_PC_Create(&shell_la);
    PCShellSetApply(precond_Method,_user_PC_Apply_Sparse);
    PCShellSetContext(precond_Method,shell_la);
    PCShellSetDestroy(precond_Method,_user_PC_Destroy);
    _user_PC_SetUp(precond_Method,linear_solver,PCSymMat);
  }else if (gUser_Defined_PC == 0) {
    ierr = PCSetFromOptions(precond_Method);assert(ierr == 0);
  }

  ierr = KSPSetFromOptions(mKsp);assert(ierr == 0);
  ierr = KSPSetTolerances(mKsp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); assert(ierr == 0);

  if (gUser_Defined_PC != 0) 
    assert(linear_solver);


}


PetscIterativeSolver_Sparse::~PetscIterativeSolver_Sparse()
{
  int ierr;

  if( deleteKSP ) { // We made it, we own it.
    ierr = KSPDestroy( &mKsp ); assert( ierr  == 0);
  }
}

void PetscIterativeSolver_Sparse::diagonalChanged( int /* idiag */, int /* extent */ )
{
  linear_solver->matrixChanged();
  this->matrixChanged();
}

int PetscIterativeSolver_Sparse::matrixChanged()
{
  int ierr;
//  negEigVal = linear_solver->matrixChanged();

  if (gUser_Defined_SymMat == 0 && gUser_Defined_PC == 0) 
  {
    fullSymMat->symmetrize_valonly(inputMatptr->M(),goffIDX);
	negEigVal = correct_negEigVal;
  }else if(gUser_Defined_SymMat == 1 && gUser_Defined_PC == 0){
	assert ("this cannot happen"&&0);
  }else if(gUser_Defined_SymMat == 0 && gUser_Defined_PC != 0){
	fullSymMat->symmetrize_valonly(inputMatptr->M(),goffIDX);
	if(PCgoffIDX)
	  PCSymMat->symmetrize_valonly(inputMatptr->M(),PCgoffIDX);
	negEigVal = linear_solver->matrixChanged();
	if(negEigVal>=0) negEigVal = correct_negEigVal;
  }else if(gUser_Defined_SymMat == 1 && gUser_Defined_PC != 0){
	if(PCgoffIDX)
	  PCSymMat->symmetrize_valonly(inputMatptr->M(),PCgoffIDX);
	negEigVal = linear_solver->matrixChanged();
	if(negEigVal>=0) negEigVal = correct_negEigVal;
  } 


  
  ierr = KSPSetOperators(mKsp, LinSysMat_PETSC, LinSysMat_PETSC); assert(ierr == 0);

  return negEigVal;
}



void PetscIterativeSolver_Sparse::solve( OoqpVector& rhs_in)
{
  int ierr;
  int nb_row, nb_col;
  inputMatptr->getSize(nb_row,nb_col);


  SimpleVector &  x_sol = dynamic_cast<SimpleVector &>(rhs_in);
  rhs_back->copyFrom(x_sol);


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  //   Solve the linear system Ax=b by petsc
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  Vec b_la,x_la;
  VecCreateSeqWithArray(PETSC_COMM_SELF,1,nb_col,rhs_back->elements(),&b_la);
  VecCreateSeqWithArray(PETSC_COMM_SELF,1,nb_col,x_sol.elements(),&x_la);

#ifdef TIMING
  double tTemp=MPI_Wtime();
#endif

  ierr = KSPSolve(mKsp,b_la,x_la); assert( ierr  == 0);

#ifdef TIMING
  tTot+=MPI_Wtime()-tTemp;
  std::cout<<"total GMRES time: "<<tTot<<std::endl;
  std::cout<<"total Mine time: "<<tMine<<std::endl;
#endif


  VecDestroy(&b_la);
  VecDestroy(&x_la);  

  
  ierr = KSPGetIterationNumber(mKsp, &KryIter); assert( ierr  == 0);

//  std::cout<<"using "<<KryIter<<" GMRES iterations"<<std::endl;
  
  total_kry_iter += KryIter;  

}





// tell petsc how to compute Ax+y
PetscErrorCode _user_MatMult(Mat A, Vec x,  Vec y)
{
#ifdef TIMING
  double tTemp=MPI_Wtime();
#endif

  SparseSymMatrix   *Amat;
  MatShellGetContext(A, &Amat);

  int nb_row, nb_col;
  Amat->getSize(nb_row,nb_col);

  double* x_array, *y_array;
  VecGetArray(x,&x_array);
  VecGetArray(y,&y_array);

  SimpleVector x_temp(x_array,nb_col);
  SimpleVector y_temp(y_array,nb_col);

  Amat->mult(0.0,y_temp,1.0,x_temp);
  
  VecRestoreArray(x,&x_array);
  VecRestoreArray(y,&y_array);

#ifdef TIMING
  tMine += MPI_Wtime()-tTemp;
#endif

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
PetscErrorCode _user_PC_Apply_Sparse(PC pc,Vec b_in,Vec x_out)
{ 
#ifdef TIMING
  double tTemp=MPI_Wtime();
#endif

  user_PC *shell;
  PCShellGetContext(pc,(void**)&shell);
  
  SparseSymMatrix *Amat = dynamic_cast<SparseSymMatrix*>(shell->mat);

  int nb_row, nb_col;
  Amat->getSize(nb_row,nb_col);

  double* b_array, *x_array;
  VecGetArray(b_in,&b_array);
  VecGetArray(x_out,&x_array);

  SimpleVector b_temp(b_array,nb_col);
  SimpleVector x_temp(x_array,nb_col);
  x_temp.copyFrom(b_temp);

  shell->lin_solver->solve(x_temp);

  VecRestoreArray(b_in,&b_array);  
  VecRestoreArray(x_out,&x_array);

#ifdef TIMING
  tMine += MPI_Wtime()-tTemp;  
#endif
  return 0;
}


