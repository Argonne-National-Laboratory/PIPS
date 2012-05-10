/* PIPS                                                               *
 * Author: Miles Lubin                                                *
 */

#include "WSMPSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#include <omp.h>

#include <mpi.h>
#include <iostream>
#include <fstream>

extern int gOoqpPrintLevel;

#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif

extern "C" {
	void wssmp_(const int*N,int*IA,int*JA,double*AVALS,
									double*DIAG,int*PERM,int*INVP,
									double*B,const int*LDB,const int*NRHS,double*AUX,
									int*NAUX,int*MRP,int*IPARM,double*DPARM);
  void wsmp_clear();
	void wstoremat_(int*,int*);
	void wrecallmat_(int*,int*);
}

#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

//#ifdef __bg__
#define PIVOT
//#endif


int zero = 0, one = 1;

int WSMPSolver::instances = 0;

WSMPSolver::WSMPSolver( SparseSymMatrix * sgm ) : first(true)
{
	
	nnz = sgm->numberOfNonZeros();
	n = sgm->size();
  mStorage = SparseStorageHandle( sgm->getStorage() );


	krowMt = new int[n+1];
	jcolMt = new int[nnz];
	Mt     = new double[nnz];
	
	//wsmp_initialize();
	perm = new int[n];
	invp = new int[n];
	diagmap = new int[2*n];
  iparm[0] = 0;
	iparm[1] = 0;
	iparm[2] = 0;

	nthreads = omp_get_max_threads();
	// set default values
#ifdef __bg__
  // hack for bgp so that wsmp threads
  // are correctly distributed across cores
  if (nthreads == 4) {
    omp_set_num_threads(5);
  }
  #pragma omp parallel if (nthreads == 4)
  {
    #pragma omp single
    wssmp_(&n,krowMt,jcolMt,Mt,NULL,perm,invp,NULL,&zero,&zero,NULL,&zero,NULL,iparm, dparm);
  }
  omp_set_num_threads(nthreads);
#else
  wssmp_(&n,krowMt,jcolMt,Mt,NULL,perm,invp,NULL,&zero,&zero,NULL,&zero,NULL,iparm, dparm);
#endif

  if (iparm[63] != 0) {
		printf("WSMP initialization error: %d\n", iparm[63]);
	}
	assert(iparm[63] == 0);

	// compressed row format
	iparm[3] = 0;
	// C-style numbering
	iparm[4] = 0;
	// recommended values for interior-point systems
  iparm[15] = 1;
	iparm[16] = 0;
	iparm[17] = 1;
	iparm[18] = 0;
	iparm[19] = 1;
	iparm[26] = 1;
  
#ifdef PIVOT
	iparm[30] = 2; // LDLt with pivoting
#else
	iparm[30] = 1; // NO PIVOTING!!!
#endif

	// max 1 steps of iterative refinement
  iparm[5] = 1;

//#ifndef __bg__
  // turn off scaling, otherwise can't solve against 
  // factors correctly.
  // this took 3 days of debugging to discover
//  iparm[9] = 2;
//#endif

#ifndef PIVOT
  iparm[9] = 1;
#endif

// set the largest value of ThresholdPivoting parameter we are
  // willing to tolerate.
  kThresholdPivotingMax = 0.1;

  // set the increase factor for ThresholdPivoting parameter
  kThresholdPivotingFactor = 10.0;

  // set the required precision for each linear system solve
  kPrecision = 1.e-7;

  //dparm[9] = 1.e-10; // treat as zero parameter

	instance = instances++;

	assert(instances <= 64);
}

void WSMPSolver::firstCall()
{
	// WSMP expects upper-trangle of row-major matrix,
	// but OOQP stores lower triangle
	// so, transpose
  mStorage->transpose(krowMt, jcolMt, Mt);
  // make a map for the diagonal elements
  for (int i = 0; i < n; i++) {
    for (int k = mStorage->krowM[i]; k < mStorage->krowM[i+1];k++) {
      int j = mStorage->jcolM[k];
      if (j==i) {
        diagmap[2*i+1] = k;
        break;
      }
    }
    for (int k = krowMt[i]; k < krowMt[i+1];k++) {
      int j = jcolMt[k];
      if (j==i) {
        diagmap[2*i] = k;
        break;
      }
    }
  }

  // ordering and symbolic factorization
	iparm[1] = 1;
	iparm[2] = 2;
	wssmp_(&n,krowMt,jcolMt,Mt,NULL,perm,invp,NULL,&zero,&zero,NULL,&zero,NULL,iparm, dparm);
	if (iparm[63] != 0) {
		printf("WSMP ordering/symbolic factorization error: %d\n", iparm[63]);
	}
	assert(iparm[63] == 0);
	//printf("symbolic done\n");
  first = false;
}  
void WSMPSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
	printf("diagonal changed\n");
  this->matrixChanged();
}

void WSMPSolver::matrixChanged()
{
	omp_set_num_threads(1);
	
	if (first) {
    firstCall();
  } else {
    // just need to update diagonal
    // should really do this in diagonalChanged,
    // but that isn't called
    for (int i = 0; i < n; i++) {
      Mt[diagmap[2*i]] = mStorage->M[diagmap[2*i+1]];
    }
		wrecallmat_(&instance,&iparm[63]);
  	assert(iparm[63] == 0);
  }
  
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  /*char fname[50];
  int i;
  sprintf(fname, "dump-node-%d-iter-%d.dat",mype,(int)g_iterNumber);
  printf("iter = %d, fname = %s\n", (int)g_iterNumber, fname);
  ofstream fd(fname);
  fd << scientific; 
  fd.precision(16);
  fd << n << endl;
  fd << nnz << endl;
  for (i = 0; i <= n; i++) {
    fd << krowMt[i] << " ";
  }
  fd << endl;
  for (i = 0; i < nnz; i++) {
    fd << jcolMt[i] << " ";
  }
  fd << endl;
  for (i = 0; i < nnz; i++) {
    fd << Mt[i] << " ";
  }
  fd << endl;
*/

	bool redo = false;
  do { 
    // numerical factorization
    iparm[1] = 3;
    iparm[2] = 3;
    if (redo) dparm[14] = 1.0;
    wssmp_(&n,krowMt,jcolMt,Mt,NULL,perm,invp,NULL,&zero,&zero,NULL,&zero,NULL,iparm, dparm);
    if (iparm[63] != 0) {
      printf("WSMP numerical factorization error: %d on rank %d\n", iparm[63],mype);
      if (iparm[63] > 0 && dparm[10] < kThresholdPivotingMax) {
        dparm[10] *= kThresholdPivotingFactor;
        assert(dparm[10] < 1);
        printf("WSMP Retrying with pivot threshold %f\n",dparm[10]);
        redo = true;
      } else{
        // TODO: better error handling
        assert(iparm[63] == 0);
      }
    } else {
      redo = false;
    }
  } while (redo);


	wstoremat_(&instance,&iparm[63]);
	assert(iparm[63] == 0);

	omp_set_num_threads(nthreads);
}
 
void WSMPSolver::solve( OoqpVector& rhs_in )
{
  omp_set_num_threads(1);

  wrecallmat_(&instance,&iparm[63]);
	assert(iparm[63] == 0);

	SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
	double * drhs = rhs.elements();
  double * rhscpy = new double[n];
  memcpy(rhscpy,drhs,n*sizeof(double)); // make copy of rhs in case we need to redo
  bool done = false;

  while (!done) {
    iparm[1] = 4;
    iparm[2] = 5;
    iparm[29] = 0;

    wssmp_(&n,krowMt,jcolMt,Mt,NULL,perm,invp,drhs,&n,&one,NULL,&zero,NULL,iparm, dparm);
    if (iparm[63] != 0) {
      printf("WSMP backsolve error: %d\n", iparm[63]);
    }
    assert(iparm[63] == 0);
    if (dparm[6] > kPrecision) {
      printf("WSMP backsolve large relative error: %f\n",dparm[6]);
      if (dparm[10] < kThresholdPivotingMax) {
        dparm[10] *= kThresholdPivotingFactor;
        printf("WSMP refactoring with pivoting threshold %f\n",dparm[10]);
        matrixChanged();
        memcpy(drhs,rhscpy,n*sizeof(double));
      } else {
        printf("WSMP pivoting threshold already at maximum\n");
        done = true;
      }
    } else {
      done = true;
    }

  }
	wstoremat_(&instance,&iparm[63]);
	assert(iparm[63] == 0);

  delete [] rhscpy;
	
	omp_set_num_threads(nthreads);

}


WSMPSolver::~WSMPSolver()
{
	delete [] krowMt;
	delete [] jcolMt;
	delete [] Mt;
	delete [] perm;
	delete [] invp;
  delete [] diagmap;
  instances--;
}

void WSMPSolver::Lsolve( OoqpVector& x )
{
  solve(1,x);
}

void WSMPSolver::Dsolve( OoqpVector& x )
{
  solve(3,x);
}

void WSMPSolver::Ltsolve( OoqpVector& x )
{
  solve(2,x);
}


void WSMPSolver::solve(int solveType, OoqpVector& rhs_in)
{
  omp_set_num_threads(1);

	wrecallmat_(&instance,&iparm[63]);
	assert(iparm[63] == 0);

	SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
	double * drhs = rhs.elements();

	iparm[1] = 4;
	iparm[2] = 4;
	iparm[29] = solveType;

	wssmp_(&n,krowMt,jcolMt,Mt,NULL,perm,invp,drhs,&n,&one,NULL,&zero,NULL,iparm, dparm);
	if (iparm[63] != 0) {
		printf("WSMP backsolve error: %d\n", iparm[63]);
	}
	assert(iparm[63] == 0);

	wstoremat_(&instance,&iparm[63]);
	assert(iparm[63] == 0);

	omp_set_num_threads(nthreads);
}



void WSMPSolver::solve(GenMatrix& rhs_in)
{
  omp_set_num_threads(1);

	wrecallmat_(&instance,&iparm[63]);
	assert(iparm[63] == 0);

	DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
	int N,NRHS;
	// rhs vectors are on the "rows", for continuous memory
	rhs.getSize(NRHS,N);
	assert(n==N);

  /*
  //static bool first = true;
  int mype; MPI_Comm_rank(MPI_COMM_WORLD,&mype);
  char fname[50];
  int i;
  sprintf(fname, "dump-node-%d-iter-%d.dat",mype,(int)g_iterNumber);
  printf("NRHS = %d, iter = %d, fname = %s\n", NRHS,(int)g_iterNumber, fname);
  ofstream fd(fname);
  fd << scientific; 
  fd.precision(16);
  fd << n << endl;
  fd << nnz << endl;
  for (i = 0; i <= n; i++) {
    fd << krowMt[i] << " ";
  }
  fd << endl;
  for (i = 0; i < nnz; i++) {
    fd << jcolMt[i] << " ";
  }
  fd << endl;
  for (i = 0; i < nnz; i++) {
    fd << Mt[i] << " ";
  }
  fd << endl;
  fd << NRHS << endl;
  for (int k = 0; k < NRHS; k++) {
    for (i = 0; i < n; i++) {
      fd << rhs[k][i] << " ";
    }
    fd << endl;
  }


*/

	iparm[29] = 0;

  const int BLOCKSIZE = 20;


	for (int startcol = 0; startcol < NRHS; startcol += BLOCKSIZE) {
		double *drhs = rhs[startcol];
		int endcol = MIN(startcol+BLOCKSIZE,NRHS);
		int numcols = endcol-startcol;
	
		iparm[1] = 4;
    iparm[2] = 5;
		wssmp_(&n,krowMt,jcolMt,Mt,NULL,perm,invp,drhs,&n,&numcols,NULL,&zero,NULL,iparm, dparm);
		if (iparm[63] != 0) {
			printf("WSMP backsolve error (multiple RHS): %d\n", iparm[63]);
		}
		assert(iparm[63] == 0);

	}
/*
 for (int k = 0; k < NRHS; k++) {
    for (i = 0; i < n; i++) {
      fd << rhs[k][i] << " ";
    }
    fd << endl;
  }*/

	wstoremat_(&instance,&iparm[63]);
	assert(iparm[63] == 0);

	omp_set_num_threads(nthreads);
}



