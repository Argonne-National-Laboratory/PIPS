/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <stdlib.h>
#include <iostream>
using namespace std;

#include "PardisoSchurSolver.h"
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
extern "C" void pardiso_schur(void*, int, int, int, double*, int*, int*);

extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);


PardisoSchurSolver::PardisoSchurSolver( SparseSymMatrix * sgm )
{
  Msys = sgm;
  n = sgm->size();
  // - we do not have the augmented system yet
  //nnz=sgm->numberOfNonZeros();

  //krowM = new int[n+1];
  //jcolM = new int[nnz];
  //M = new double[nnz];

  first = true; firstSolve = true;
  //nvec=new double[n];
   
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

void PardisoSchurSolver::firstCall()
{
  pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error); 
  if (error!=0) {
    cout << "PardisoSchurSolver ERROR during pardisoinit:" << error << "." << endl;
    exit(1);
  }
} 

void PardisoSchurSolver::firstSolveCall()
{


} 

 
void PardisoSchurSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void PardisoSchurSolver::matrixChanged()
{
  if(first) { firstCall(); first=false;}

  // we don't have the right hand-size, therefore we can't (re)factorize the augmented system
}
 
void PardisoSchurSolver::schur_solve(SparseGenMatrix& R, 
				     SparseGenMatrix& A,
				     SparseGenMatrix& C,
				     DenseSymMatrix& SC)
{
  if(firstSolve) {
    // create the augmented system

    int nx1,my1,mz1, nxP;
    int nnz=0;
    R.getSize(nx1,nxP); nnz += R.numberOfNonZeros();
    A.getSize(my1,nxP); nnz += A.numberOfNonZeros();
    C.getSize(mz1,nxP); nnz += C.numberOfNonZeros();

    n = nx1+my1+mz1+nxP;
    assert( Msys->size() == nx1+my1+mz1 );
    nnz += Msys->numberOfNonZeros();

    SparseSymMatrix augSys( n, nnz);
      
    //setup the diagonal, as required by Pardiso.
    SimpleVector* v = new SimpleVector(n) ;
    v->setToZero();
    augSys.setToDiagonal(*v);
    delete v;
    
    cout << "Diagonals set" << endl;
    
    //put the blocks in augmented system
    augSys.symAtPutSubmatrix( 0, 0, *Msys, 0, 0, Msys->size(), Msys->size());

    int nR, mR; R.getSize(nR,mR);
    cout << "R: " << nR << " x " << mR << "   nnz=" << R.numberOfNonZeros() << endl;
    if(R.numberOfNonZeros()>0) {
      augSys.symAtPutSubmatrix( Msys->size(), 0, R, 0, 0, nR, mR);
    }
 
   int nA, mA; A.getSize(nA,mA);
    cout << "A: " << nA << " x " << mA << "   nnz=" << A.numberOfNonZeros() << endl;
    if(A.numberOfNonZeros()>0 && nA>0) {
      augSys.symAtPutSubmatrix( Msys->size()+nR, 0, A, 0, 0, nA, mA);
    }

    

    int nC, mC; C.getSize(nC, mC);
    cout << "C: " << nC << " x " << mC << "   nnz=" << C.numberOfNonZeros() << endl;
    if(C.numberOfNonZeros()>0 && nC>0) {
      augSys.symAtPutSubmatrix( Msys->size()+nR+nA, 0, C, 0, 0, nC, mC);
    }

    cout << " AUG SYS nnz=" << augSys.numberOfNonZeros() << endl;

    // we need to transpose before giving it to Pardiso
    krowM = new int[n+1];
    jcolM = new int[nnz];
    M = new double[nnz];

    augSys.getStorageRef().transpose(krowM, jcolM, M);
    //save the indeces for diagonal entries ( first Msys->size() )
    // to do
  
    firstSolve=false; first = false;  // to do
  } //end of if(firstSolve)

  //update diagonal entries in krowM, jcolM and M
  // to do

  // call PARDISO
  int error;
  pardiso_chkmatrix(&mtype,&n, M,krowM,jcolM, &error);
  if(error != 0) {
    cout << "PARDISO check matrix error" << error << endl;
    exit(1);
  }
  phase=12; //Analysis, numerical factorization
  int maxfct=1;
  int mnum=1;
  int nrhs=1;
  iparm[2]=num_threads;
  //iparm[1] = 2; // 2 is for metis, 0 for min degree 
  iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
  iparm[12] = 1; // improved accuracy for IPM KKT; used with IPARM(11)=1; 
                 // if needed, use 2 for advanced matchings and higer accuracy.
  iparm[23] = 1; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

  int msglvl=1;  // with statistical information
  iparm[32] = 1; // compute determinant
  iparm[37] = Msys->size(); //compute Schur-complement
  int nr=n-Msys->size();
  cout << "Factorizing AUG Sys. n=" << n << " nr=" << nr << endl;

  pardiso (pt , &maxfct , &mnum, &mtype , &phase ,
	   &n, M, krowM, jcolM,
	   NULL, &nrhs,
	   iparm , &msglvl, NULL, NULL, &error, dparm );

  cout << "Schur Factorization completed" << endl;
  //cout << "factorizing the matrix took:" << MPI_Wtime()-tt << endl;
  if ( error != 0) {
    printf ("PardisoSolver - ERROR during factorization: %d\n", error );
    assert(false);
  }
  int* colSC=new int[nr+1];
  int* rowSC=new int[iparm[38]];
  double* mSC=new double[iparm[38]];

  pardiso_schur(pt, maxfct, mnum, mtype, mSC, rowSC, colSC);

  //Schur-complement extracted... updating dense SC matrix
  for(int j=0; j<nr; j++) {
    for(int i=colSC[j]; i<colSC[j+1]; i++)
      SC[i][j]=mSC[i];
  }
}

void PardisoSchurSolver::solve( OoqpVector& rhs_in )
{
  assert(false && "Function not supported. Use PardisoSolver for this functionality.");    
}

void PardisoSchurSolver::solve(GenMatrix& rhs_in)
{
  assert(false && "Function not supported. Use PardisoSolver for this functionality.");
}

PardisoSchurSolver::~PardisoSchurSolver()
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
    printf ("PardisoSchurSolver - ERROR in pardiso release: %d", error ); 
  }
  delete[] jcolM;
  delete[] krowM;
  delete[] M;
  //delete[] nvec;
}

