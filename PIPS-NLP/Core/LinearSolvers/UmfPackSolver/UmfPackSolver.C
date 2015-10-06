/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#include "UmfPackSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"

#include <cmath>
#include "stdio.h"


#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

extern int gOoqpPrintLevel;
extern int gOuterSolve;
extern int separateHandDiag;
extern int gRS_MaxIR;

extern int gPipsPrtLV;

extern double gRS_LU_PivotLV;


#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif



UmfPackSolver::UmfPackSolver( SparseGenMatrix * sgm )
{
  mStorage = SparseStorageHandle( sgm->getStorage() );
  assert( mStorage->n ==  mStorage->m);
  n        = mStorage->n;

  assert(gOuterSolve ==3 && separateHandDiag==0);
  nnz = mStorage->krowM[n];

  umfpack_di_defaults (umf_Control) ;

  umf_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
  umf_Control[UMFPACK_PIVOT_TOLERANCE] = gRS_LU_PivotLV;
  umf_Control[UMFPACK_IRSTEP] = gRS_MaxIR;
  firstCallFlag = true;
 
  Numeric=NULL;
  Symbolic=NULL;
  isSymm=0;

}

UmfPackSolver::UmfPackSolver( SparseSymMatrix * sgm )
{
  mStorage = SparseStorageHandle( sgm->getStorage() );
  assert( mStorage->n ==  mStorage->m);
  n        = mStorage->n;

  assert(gOuterSolve ==3 && separateHandDiag==0);
  nnz = 2*mStorage->krowM[n]-n;

  umfpack_di_defaults (umf_Control) ;

  umf_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
  umf_Control[UMFPACK_PIVOT_TOLERANCE] = gRS_LU_PivotLV;
  umf_Control[UMFPACK_IRSTEP] = gRS_MaxIR;
  firstCallFlag = true;
 
  Numeric=NULL;
  Symbolic=NULL;
  isSymm=1;

}

UmfPackSolver::UmfPackSolver( SparseSymMatrix * sgm, const int numOfNegEigVal_in )
{
  mStorage = SparseStorageHandle( sgm->getStorage() );
  assert( mStorage->n ==  mStorage->m);
  n        = mStorage->n;

  assert(gOuterSolve ==3 && separateHandDiag==0);
  nnz = 2*mStorage->krowM[n]-n;

  umfpack_di_defaults (umf_Control) ;

  umf_Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;
  umf_Control[UMFPACK_PIVOT_TOLERANCE] = gRS_LU_PivotLV;
  umf_Control[UMFPACK_IRSTEP] = gRS_MaxIR;
  firstCallFlag = true;
 
  Numeric=NULL;
  Symbolic=NULL;
  isSymm=1;
  negEigVal=numOfNegEigVal_in;
}

UmfPackSolver::~UmfPackSolver()
{

}

void UmfPackSolver::firstCall()
{
  int status;
  
  kcolbegM 	= (int*) malloc((n+1)*sizeof(int));
  irowM 	= (int*) malloc(nnz*sizeof(int));
  eleM 		= (double*) malloc(nnz*sizeof(double));
  eleMap 	= (int*) malloc(nnz*sizeof(int));
  
  int *krowbegM_in = mStorage->krowM;
  int *jcolM_in = mStorage->jcolM;
  double *eleM_in = mStorage->M;
  
  for( int i = 0; i < n+1; i++ ) 
  	kcolbegM[i]=0;

  for( int i = 0; i < n; i++ ) 
	for( int k = krowbegM_in[i]; k < krowbegM_in[i+1]; k++ ) 
	  kcolbegM[jcolM_in[k]+1]++;
	
  for( int i = 1; i < n+1; i++ ) 
  	kcolbegM[i]+=kcolbegM[i-1];

  int *findkincol 	= (int*) malloc(n*sizeof(int));
  memcpy(findkincol,kcolbegM, n*sizeof(int));

  for( int i = 0; i < n; i++ ){ 
	for( int k = krowbegM_in[i]; k < krowbegM_in[i+1]; k++ ){
	  int local_k = findkincol[jcolM_in[k]]++;
	  irowM[local_k]  = i;
	  eleM[local_k]   = eleM_in[k];
	  eleMap[local_k] = k;
    }   
  }
  for( int i = 0; i < n; i++ )
  	assert(findkincol[i]==kcolbegM[i+1]);

  free (findkincol);  

  status = umfpack_di_symbolic (n, n, kcolbegM, irowM, eleM, &Symbolic, umf_Control, umf_Info) ;

  switch( status ) {
      case UMFPACK_OK:		
		break;
	  case UMFPACK_ERROR_n_nonpositive:
		assert("n < 0 !?"&&0);	
	  case UMFPACK_ERROR_invalid_matrix:
		assert("Check input sparse matrix!"&&0);
	  case UMFPACK_ERROR_out_of_memory:
		assert("Out of memory"&&0);	
	  case UMFPACK_ERROR_argument_missing:
		assert("One or more required arguments are missing!"&&0);	
	  case UMFPACK_ERROR_internal_error:
		assert("Call author! :-)"&&0);
	  default:
		assert( "unknowen error!" && 0 );
  }

  firstCallFlag=0;
}  

// input M is lower triangular
void UmfPackSolver::firstCall_Sym()
{
  int status;
  
  kcolbegM 	= (int*) malloc((n+1)*sizeof(int));
  irowM 	= (int*) malloc(nnz*sizeof(int));
  eleM 		= (double*) malloc(nnz*sizeof(double));
  eleMap 	= (int*) malloc(nnz*sizeof(int));
  
  int *krowbegM_in = mStorage->krowM;
  int *jcolM_in = mStorage->jcolM;
  double *eleM_in = mStorage->M;
  
  for( int i = 0; i < n+1; i++ ) 
  	kcolbegM[i]=0;

  for( int i = 0; i < n; i++ ) 
	for( int k = krowbegM_in[i]; k < krowbegM_in[i+1]; k++ ){
	  kcolbegM[jcolM_in[k]+1]++;
	  if(jcolM_in[k]!=i)
	  	kcolbegM[i+1]++;
	} 
	  
  for( int i = 1; i < n+1; i++ ) 
  	kcolbegM[i]+=kcolbegM[i-1];

  int *findkincol 	= (int*) malloc(n*sizeof(int));
  memcpy(findkincol,kcolbegM, n*sizeof(int));

  for( int i = 0; i < n; i++ ){ 
	for( int k = krowbegM_in[i]; k < krowbegM_in[i+1]; k++ ){
	  int local_k = findkincol[jcolM_in[k]]++;
	  irowM[local_k]  = i;
	  eleM[local_k]   = eleM_in[k];
	  eleMap[local_k] = k;
	  if(jcolM_in[k]!=i){
	  	local_k = findkincol[i]++;
	    irowM[local_k]  = jcolM_in[k];
	    eleM[local_k]   = eleM_in[k];
	    eleMap[local_k] = k;
	  }
    }   
  }
  for( int i = 0; i < n; i++ )
  	assert(findkincol[i]==kcolbegM[i+1]);

  free (findkincol);  

  status = umfpack_di_symbolic (n, n, kcolbegM, irowM, eleM, &Symbolic, umf_Control, umf_Info) ;

  switch( status ) {
      case UMFPACK_OK:		
		break;
	  case UMFPACK_ERROR_n_nonpositive:
		assert("n < 0 !?"&&0);	
	  case UMFPACK_ERROR_invalid_matrix:
		assert("Check input sparse matrix!"&&0);
	  case UMFPACK_ERROR_out_of_memory:
		assert("Out of memory"&&0);	
	  case UMFPACK_ERROR_argument_missing:
		assert("One or more required arguments are missing!"&&0);	
	  case UMFPACK_ERROR_internal_error:
		assert("Call author! :-)"&&0);
	  default:
		assert( "unknowen error!" && 0 );
  }

  firstCallFlag=0;
}  


int UmfPackSolver::matrixChanged()
{
  int status;

  matrixSingular=0;
  
  if(firstCallFlag && isSymm==0){
  	this->firstCall();
  }else if(firstCallFlag && isSymm==1){
  	this->firstCall_Sym();
  }else{ 
	for( int i = 0; i < nnz; i++ ){
	  eleM[i]   = mStorage->M[eleMap[i]];
	}
  }
  

  if(Numeric)
  	freeNumFactInfo();
  

  status = umfpack_di_numeric (kcolbegM, irowM, eleM, Symbolic, &Numeric,  umf_Control, umf_Info) ;


  switch( status ) {
      case UMFPACK_OK:		
		break;
	  case UMFPACK_WARNING_singular_matrix:
	  	matrixSingular=1;
		break;
	  case UMFPACK_ERROR_out_of_memory:
		assert("Out of memory"&&0);		
	  case UMFPACK_ERROR_argument_missing:
		assert("One or more required arguments are missing!"&&0);	
	  case UMFPACK_ERROR_invalid_Symbolic_object:
		assert("Symbolic object provided as input is invalid!"&&0);
	  case UMFPACK_ERROR_different_pattern:
		assert("impossible!"&&0);
	  default:
		assert( "unknowen error!" && 0 );
  }

  if(matrixSingular==1) 
  	negEigVal = -1;
//  else
//	negEigVal = 0;

  return negEigVal;  
}



void UmfPackSolver::solve(int solveType, OoqpVector& rhs_in)
{
  int status; 
  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);

  double * drhs   	= rhs.elements();
  double * solX 	= new double[n];

  status = umfpack_di_solve(solveType, kcolbegM, irowM, eleM, solX, drhs, Numeric, umf_Control, umf_Info) ;

  switch( status ) {
      case UMFPACK_OK:		
		break;
	  case UMFPACK_WARNING_singular_matrix:
	  	assert("impossible!"&&0);
	  case UMFPACK_ERROR_out_of_memory:
		assert("Out of memory"&&0);		
	  case UMFPACK_ERROR_argument_missing:
		assert("One or more required arguments are missing!"&&0);	
	  case UMFPACK_ERROR_invalid_Numeric_object:
		assert("Numeric object provided as input is invalid!"&&0);
	  case UMFPACK_ERROR_invalid_system:
		assert("The sys argument is not valid, or the matrix A is not square!"&&0);
	  default:
		assert( "unknowen error!" && 0 );
  }

  memcpy( &drhs[0], &solX[0],  n * sizeof( double ) );



  delete [] solX;
}


void UmfPackSolver::solve( OoqpVector& rhs_in )
{
  solve(UMFPACK_A,rhs_in);
} 

void UmfPackSolver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;
  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert(n==N);
    
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solve(v);
  } 

}

void UmfPackSolver::solveTrans( OoqpVector& rhs_in )
{
  solve(UMFPACK_At,rhs_in);
} 

void UmfPackSolver::solveTrans(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;
  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert(n==N);
    
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solveTrans(v);
  } 
}

void UmfPackSolver::Lsolve( OoqpVector& x )
{
//  solve(UMFPACK_L,x);
}

void UmfPackSolver::Ltsolve( OoqpVector& x )
{
//  solve(UMFPACK_U,x);
}


