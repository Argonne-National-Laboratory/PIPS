/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin, Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */

#include <iostream>
using namespace std;

#include "PardisoSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include <cstdlib>

//#include "mpi.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

extern int gOoqpPrintLevel;


extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);


extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);


PardisoSolver::PardisoSolver( SparseSymMatrix * sgm )
{
  Msys = sgm;
  n = sgm->size();
  nnz=sgm->numberOfNonZeros();

  krowM = new int[n+1];
  jcolM = new int[nnz];
  M = new double[nnz];

  first = true;
  nvec=new double[n];
   
  error  = 0;
  solver = 0;  /* use sparse direct solver */
  mtype  = -2; /* real  symmetric with diagonal or Bunch-Kaufman */

  /* Numbers of processors, value of OMP_NUM_THREADS */
  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
    sscanf( var, "%d", &num_threads );
  else {
    printf("Set environment OMP_NUM_THREADS");
    exit(1);
  }
}

void PardisoSolver::firstCall()
{

  pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error); 
  if (error!=0) {
    cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << endl;
    assert(false);
  }
  iparm[2] = num_threads; 
} 

 
void PardisoSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void PardisoSolver::matrixChanged()
{
  if (first) { firstCall(); first = false; }

  
  //double tt=MPI_Wtime();
  //get the matrix in upper triangular
  Msys->getStorageRef().transpose(krowM, jcolM, M);
 
  // pardiso requires diag elems even though they are exactly 0
  for( int i = 0; i < n; i++) {
    bool hasDiag=0;
    for( int j=krowM[i]; j<krowM[i+1] && !hasDiag; j++ ) {
      if( jcolM[j]==i ) hasDiag=true;
    }

    if (!hasDiag)
      assert(false);
      //cout << "NO diag elem in row " << i << endl;
  }
 
  // need Fortran indexes
  for( int i = 0; i < n+1; i++) krowM[i] += 1;
  for( int i = 0; i < nnz; i++) jcolM[i] += 1;


  // compute numerical factorization
  phase = 12; //Analysis, numerical factorization

  int maxfct=1; //max number of fact having same sparsity pattern to keep at the same time
  int mnum=1; //actual matrix (as in index from 1 to maxfct)
  int nrhs=1;
  int msglvl=1; //messaging level


  //iparm[1] = 2; // 2 is for metis, 0 for min degree 
  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 1; // improved accuracy for IPM KKT; used with IPARM(11)=1; 
                 // if needed, use 2 for advanced matchings and higer accuracy.
  
  pardiso (pt , &maxfct , &mnum, &mtype , &phase ,
	   &n, M, krowM, jcolM,
	   NULL, &nrhs,
	   iparm , &msglvl, NULL, NULL, &error, dparm );
  //cout << "factorizing the matrix took:" << MPI_Wtime()-tt << endl;
  if ( error != 0) {
    printf ("PardisoSolver - ERROR during factorization: %d\n", error );
    assert(false);
  }
}
 
void PardisoSolver::solve( OoqpVector& rhs_in )
{
    
  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
  double * drhs = rhs.elements();
  
  double * sol = nvec;

  phase = 33; //solve and iterative refinement
  int maxfct=1; //max number of fact having same sparsity pattern to keep at the same time
  int mnum=1; //actual matrix (as in index from 1 to maxfct)
  int nrhs=1;
  int msglvl=0;

  iparm[7] = 1; /* Max numbers of iterative refinement steps . */
  //iparm[5] = 1; /* replace drhs with the solution */
  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, M, krowM, jcolM, 
	   NULL, &nrhs ,
	   iparm, &msglvl, 
	   drhs, sol,
	   &error, dparm );
  if ( error != 0) {
    printf ("PardisoSolver - ERROR during solve: %d", error ); 
  }
  rhs.copyFromArray(sol);
}


PardisoSolver::~PardisoSolver()
{
  phase = -1; /* Release internal memory . */
  int maxfct=1; //max number of fact having same sparsity pattern to keep at the same time
  int mnum=1; //actual matrix (as in index from 1 to maxfct)
  int nrhs=1;
  int msglvl=0;

  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, NULL, krowM, jcolM, NULL, &nrhs,
	   iparm, &msglvl, NULL, NULL, &error, dparm );
  if ( error != 0) {
    printf ("PardisoSolver - ERROR in pardiso release: %d", error ); 
  }
  delete[] jcolM;
  delete[] krowM;
  delete[] M;
  delete[] nvec;
}



void PardisoSolver::solve(GenMatrix& rhs_in)
{
  //DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  assert(false);
}



 /*void PardisoSolver::Lsolve( OoqpVector& x )
   {
  
   }  
 */

 /*  void PardisoSolver::Dsolve( OoqpVector& x )
  {
  
  }  */

 /* void PardisoSolver::Ltsolve( OoqpVector& x )
    {
  
    } */



