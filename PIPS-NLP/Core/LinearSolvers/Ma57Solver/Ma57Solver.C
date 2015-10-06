/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP 
 * Modified by Cosmin Petra to perform solves with the factors.
 */

/* 2015. Rewritten by Nai-Yuan Chiang for NLP*/

#include "Ma57Solver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"

#include <cmath>
#include <cstdio>


#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

extern int gOoqpPrintLevel;
extern int gOuterSolve;
extern int separateHandDiag;


extern double gHSL_PivotLV;
extern double gMA57_Ordering;

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif

#ifdef TIMING
//  #include "mpi.h"
  extern int call_sol_Times_MA57;
  extern int call_fact_Times_MA57;
#endif


#include "mpi.h"

//public Ma57Solver::Ma57Solver( SparseSymMatrix * sgm, const int numOfNegEigVal_in ):this(sgm){}; 
void Ma57Solver::SetUpMa57Solver( SparseSymMatrix * sgm )
{
  irowM = 0;
  jcolM = 0;
  fact  = 0;
  ifact = 0;
  keep  = 0;

  iworkn  = NULL; dworkn = NULL;
  niworkn = ndworkn = 0;
 
  ipessimism = 1.05;
  rpessimism = 1.05;

  FNAME(ma57id)( cntl, icntl );
  icntl[1-1] = 0;	   // don't print warning messages
  icntl[2-1] = 0;	   // no Warning messages
  icntl[4-1] = 1;	   // no statistics messages
  icntl[5-1] = 0;	   // no Print messages.

  icntl[6-1] = gMA57_Ordering; // 4 use Metis; 5 automatic choice(MA47 or Metis); 3 min
	      // degree ordering as in MA27; 2 use MC47;  
  icntl[7-1] = 1;      /* Pivoting strategy. */

  icntl[9-1] = 10; // up to 10 steps of iterative refinement
  icntl[11-1] = 16;
  icntl[12-1] = 16; 

  icntl[15-1] = 1;
  icntl[16-1] = 0;

  setValTreatAsZero(1.e-20);

  // set initial value of Threshold parameter
  setPivotTol(gHSL_PivotLV);


  // set the required precision for each linear system solve
  kPrecision = 1.e-9;

  mStorage = SparseStorageHandle( sgm->getStorage() );
  n        = mStorage->n;

  if (gOuterSolve >=3 && separateHandDiag==1){
  	nnz = mStorage->krowM[n] + n;
	M   = new double[nnz];
	memcpy( &M[0], &mStorage->M[0],  mStorage->krowM[n] * sizeof( double ) );
	memcpy( &M[mStorage->krowM[n]], &mStorage->additiveDiag[0], n * sizeof( double ) );
  }else {
	nnz = mStorage->krowM[n];
	M	= mStorage->M;	
  }
  
}


Ma57Solver::Ma57Solver( SparseSymMatrix * sgm, const int numOfNegEigVal_in )
{
  SetUpMa57Solver(sgm);
}

void Ma57Solver::firstCall()
{
  irowM = new int[nnz];
  jcolM = new int[nnz];

  int * krowM = mStorage->krowM;


  if(separateHandDiag==0){

	for( int i = 0; i < n; i++ ) {
	  for( int k = krowM[i]; k < krowM[i+1]; k++ ) {
		irowM[k] = i + 1;
	  }
	}
	
    for( int k = 0; k < nnz; k++ ) {
      jcolM[k] = mStorage->jcolM[k] + 1;
    }  	
  }else if (gOuterSolve >=3 && separateHandDiag==1){
	int next;
	
	next = 0;
	for( int i = 0; i < n; i++ ) {
	  for( int k = krowM[i]; k < krowM[i+1]; k++ ) {
		irowM[k+next] = i + 1;
	  }
	}
    for( int k = 0; k < nnz - n; k++ ) {
      jcolM[k+next] = mStorage->jcolM[k] + 1;
    }  	
	// add additive diag part
	for( int k = nnz - n; k < nnz; k++ ) {
      irowM[k] = k - (nnz - n) + 1;
	  jcolM[k] = irowM[k];
    }	
  }else{assert(0);}

  lkeep = ( nnz > n ) ? (5 * n + 2 *nnz + 42) : (6 * n + nnz + 42);
  keep = new int[lkeep];

  // Initialize to 0 as otherwise MA57ED can sometimes fail
  for (int k=0; k<lkeep; k++) {
	keep[k] = 0;
  }

  int * iwork = new int[5 * n];
  FNAME(ma57ad)( &n, &nnz, irowM, jcolM, &lkeep, keep, iwork, icntl,
	   info, rinfo );

  delete [] iwork;

  lfact = (int) (rpessimism * info[8]);
  fact  = new double[lfact];

  lifact = (int) (ipessimism * info[9]);
  ifact  = new int[lifact];

}  
void Ma57Solver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

int Ma57Solver::matrixChanged()
{
  if( !keep ) this->firstCall();

#ifdef TIMING
	call_fact_Times_MA57++;
#endif

  int * iwork = new_iworkn(n);

  int done = 0, tries = 0;
  int matrixSingular=0;
  	
  do {
#ifdef HAVE_GETRUSAGE
    rusage before;
    if( gOoqpPrintLevel >= 100 ) {
      getrusage( RUSAGE_SELF, &before );
    }
#endif

	if(gOuterSolve >= 3 && separateHandDiag==1){
  	  memcpy( &M[0], &mStorage->M[0],  mStorage->krowM[n] * sizeof( double ) );
	  memcpy( &M[mStorage->krowM[n]], &mStorage->additiveDiag[0], n * sizeof( double ) );  

  	}
 
    FNAME(ma57bd)( &n,       &nnz,    M,     fact,  &lfact,  ifact,
	     &lifact,  &lkeep,  keep,  iwork,  icntl,  cntl,
	     info,     rinfo );
#ifdef HAVE_GETRUSAGE
    rusage  after;
    if( gOoqpPrintLevel >= 100 ) {
      getrusage( RUSAGE_SELF, &after );
      cout << "For try " << tries + 1 
	   << " the factorization took " 
	   << (double) (after.ru_utime.tv_sec - before.ru_utime.tv_sec)
	+ (after.ru_utime.tv_usec - before.ru_utime.tv_usec) / 1000000.0
	   << " seconds.\n";
    }
#endif

    switch( info[0] ) {
    case 0: done = 1;     
      break;
    case -3: {
      int ic = 0;
      int lnfact = (int) (info[16] * rpessimism);
      double * newfact = new double[lnfact];

	  int intTemp;
	
      FNAME(ma57ed)( &n, &ic, keep, 
	  		fact, &info[1], newfact, &lnfact,
	       ifact, &info[1], &intTemp, &lifact, 
	       info );

    delete [] fact;
      fact = newfact;
      lfact = lnfact;
      rpessimism *= 1.1;
    }; break;
    case -4: {
      int ic = 1;
      int lnifact = (int) (info[17] * ipessimism);
      int * nifact = new int[ lnifact ];
      FNAME(ma57ed)( &n, &ic, keep, fact, &lfact, fact, &lfact,
	       ifact, &lifact, nifact, &lnifact, info );
      delete [] ifact;
      ifact = nifact;
      lifact = lnifact;
      ipessimism *= 1.1;
    }; break;
    case 4: {
 	  matrixSingular=1;
	  done=1;
    }; break;
    default:
      if( info[0] >= 0 ) done = 1;
      assert( "unknowen error!" && 0 );
    } // end switch      
    tries++;
  } while( !done );
  
  if(matrixSingular==1) 
  	negEigVal = -1;
  else
	negEigVal = info[24-1];
  
  return negEigVal;
  
}
 
void Ma57Solver::solve( OoqpVector& rhs_in )
{
#ifdef TIMING
  call_sol_Times_MA57++;
#endif
  int job = 1; // Solve using A
  int one = 1;
  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);

  double * drhs = rhs.elements();

  int * iwork = new int[5*n];
  double * work = new double[n];

  FNAME(ma57cd)( &job,       &n,        
    	   fact,       &lfact,    ifact,  &lifact,  
    	   &one,       drhs,      &n,   
    	   work,      &n,        iwork, 
    	   icntl,      info );
  delete [] iwork;
  delete [] work;
}

Ma57Solver::~Ma57Solver()
{
  delete [] irowM;
  delete [] jcolM;
  delete [] fact;
  delete [] ifact;
  delete [] keep;
  if(iworkn) delete[] iworkn;
  if(dworkn) delete[] dworkn;
}

void Ma57Solver::solve(int solveType, OoqpVector& rhs_in)
{
  if(solveType < 1 || solveType > 4) {
    assert("Unknown JOB assigned for use in MA57CD!" && 0);
  } else if(solveType==1) {
    // we prefer iterative refinement.
    solve(rhs_in);
    return;
  } /*else */

  int job = solveType; // Solve using A
  int one = 1;
  
  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
  //!SimpleVectorHandle work( new SimpleVector(n) );

  double * drhs   = rhs.elements();
  double * dwork  = new_dworkn(n);
  int * iwork     = new_iworkn(n);
  
#ifdef HAVE_GETRUSAGE
  rusage before;
  getrusage( RUSAGE_SELF, &before );
#endif
  
  FNAME(ma57cd)( &job,       &n,        
	   fact,       &lfact,    ifact,  &lifact,  
	   &one,       drhs,      &n,   
	   dwork,      &n,        iwork, 
	   icntl,      info );
  /*  
#ifdef HAVE_GETRUSAGE
  rusage after;
  getrusage( RUSAGE_SELF, &after );
  cout << "Solution with the factored matrix took "
       << (double) (after.ru_utime.tv_sec - before.ru_utime.tv_sec)
    + (after.ru_utime.tv_usec - before.ru_utime.tv_usec) / 1000000.0
       << " seconds.\n";
#endif
  */
  //delete [] iwork; //!using a sort of cache now
}

void Ma57Solver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;
  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert(n==N);
 
  // we need checks on the residuals, can't do that with multiple RHS
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solve(v);
  } 
}

int* Ma57Solver::new_iworkn(int dim)
{
  if(niworkn != dim) {
    if(iworkn) delete[] iworkn;
    iworkn = new int[dim];
    niworkn=dim;
  }else {
    if(NULL == iworkn)
      iworkn = new int[dim];
  }
  return iworkn;
}

double* Ma57Solver::new_dworkn(int dim)
{
  if(ndworkn != dim) {
    if(dworkn) delete[] dworkn;
    dworkn = new double[dim];
    ndworkn = dim;
  }else {
    if(NULL == dworkn)
      dworkn = new double[dim];
  }
  return dworkn;
}
