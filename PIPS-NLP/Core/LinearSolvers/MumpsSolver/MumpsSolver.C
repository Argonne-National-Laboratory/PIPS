/* PIPS file - Written by Cosmin Petra, 2018 */

#include "MumpsSolver.h"

#ifndef WITHOUT_PIPS
#include "../../global_var.h"
#include "SimpleVector.h"
#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"
#endif

#include <cassert>
#include <cmath>
#include "stdlib.h"


extern double gHSL_PivotLV;
extern int gSymLinearAlgSolverForDense;



using namespace std;


MumpsSolver::MumpsSolver(const long long& globSize, MPI_Comm mumpsMpiComm, MPI_Comm pipsMpiComm)
{
  mumps_ = new DMUMPS_STRUC_C;
  mumps_->n = globSize;
  mumps_->comm_fortran = getFortranMPIComm(mumpsMpiComm);
  mumps_->a = NULL;
  mumps_->jcn = NULL;
  mumps_->irn = NULL;
  mumps_->sym = 2;//general symetric matrix
  mumps_->job = -1; //initialization
  mumps_->par = 1; //host process involved in parallel computations
  dmumps_c(mumps_);  

  setMumpsVerbosity(3);//3 moderately verbose, 1 only errors, 0 anything is surpressed

  mumps_->icntl[ 5  -1] = 0; //assembled format for matrix (triplet format)
  mumps_->icntl[18  -1] = 3; //distributed assembled format
  dmumps_c(mumps_);  

}

MumpsSolver::~MumpsSolver()
{
  mumps_->job = -2; //terminate mumps
  dmumps_c(mumps_);
  delete [] mumps_->a;
  delete mumps_;

}

// MumpsSolver::MumpsSolver( DenseSymMatrix * dm )
// {
//   mStorage = DenseStorageHandle( dm->getStorage() );

//   DMUMPS_STRUC_C* mumps_ = new DMUMPS_STRUC_C;


//   int size = mStorage->n;
//   ipiv = new int[size];
//   lwork = -1;
//   work = NULL;
//   sparseMat = 0;
  
// }

// MumpsSolver::MumpsSolver( SparseSymMatrix * sm )
// {


//   int size = sm->size();
//   mStorage = DenseStorageHandle( new DenseStorage(size,size) );

//   ipiv = new int[size];
//   lwork = -1;
//   work = NULL;
//   sparseMat = sm;
// }




//#include "mpi.h"
int MumpsSolver::matrixChanged()
{
  int n = mumps_->n;
  /* If matrix is dense use LAPACK, not MA57 */
  char fortranUplo = 'U';
  int info;
  
  // if (sparseMat) {
  //   std::fill(mStorage->M[0],mStorage->M[0]+n*n,0.);

  //   const double *sM = sparseMat->M();
  //   const int *jcolM = sparseMat->jcolM();
  //   const int *krowM = sparseMat->krowM();
  //   for (int i = 0; i < n; i++) {
  //     for (int k = krowM[i]; k < krowM[i+1]; k++) {
  //       int col = jcolM[k];
  //       mStorage->M[i][col] = sM[k];
  //     }
  //   }
  // }
  int negEigVal=0;
  return negEigVal;
}

// Centralized solution (ICNTL(21)=0)
// Distributed solution (ICNTL(21)=1)
void  MumpsSolver::solve ( double* vec )
{

}

#ifndef WITHOUT_PIPS
void MumpsSolver::solve ( OoqpVector& v )
{
  char fortranUplo = 'U';
  int info;
  int one = 1;

  int n = mStorage->n; SimpleVector &  sv = dynamic_cast<SimpleVector &>(v);


  FNAME(dsytrs)( &fortranUplo, &n, &one,	&mStorage->M[0][0],	&n,
	   ipiv, &sv[0],	&n,	&info);


  assert(info==0);
}

void MumpsSolver::solve ( GenMatrix& rhs_in )
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  char fortranUplo = 'U';
  int info;
  int nrows,ncols; rhs.getSize(ncols,nrows);

  int n = mStorage->n;

  FNAME(dsytrs)( &fortranUplo, &n, &ncols,	&mStorage->M[0][0],	&n,
	   ipiv, &rhs[0][0],	&n,	&info);

  assert(info==0);
}
#endif
void MumpsSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}


/* 0  - no output
 * 1  - error messages only
 * 2  - 1+warning messages
 * 3  - 2+diagnostics and statistics in compact output
 * >3 - 2+diagnostics and statistics in detailed output
 *
 * Method to be called only after MUMPS has been initialized (JOB=-1)
 */
void MumpsSolver::setMumpsVerbosity(int level) {

  mumps_->icntl[1 -1] = 0; //surpress error
  mumps_->icntl[2 -1] = 0; //surpress warnings, diagnostics and statistics
  mumps_->icntl[3 -1] = 0; //surpress global information, collected on the host.
  mumps_->icntl[4 -1] = 0; //surpress all of the above 
    
  switch(level) {
  case 0: {
    //everything was setup above
    break;
  }
  case 1: { 
    mumps_->icntl[1 -1] = 6; //standard output stream
    mumps_->icntl[4 -1] = 1; //Only error messages printed.
  }
  case 2: {
    mumps_->icntl[2 -1] = 6; //output stream for diagnostic printing, statistics, and warning messages.
    mumps_->icntl[4 -1] = 2; //Errors, warnings, and main statistics printed.
  }
  case 3: {
    mumps_->icntl[2 -1] = 6; //output stream for diagnostic printing, statistics, and warning messages.
    mumps_->icntl[3 -1] = 6; //output stream for global information, collected on the host.
    mumps_->icntl[4 -1] = 3; //Errors and warnings and terse diagnostics (only first ten entries of arrays) printed.
    break;
  }
  default:
    //level>3
    mumps_->icntl[2 -1] = 6; //output stream for diagnostic printing, statistics, and warning messages.
    mumps_->icntl[3 -1] = 6; //output stream for global information, collected on the host.
    mumps_->icntl[4 -1] = 4; //Errors, warnings and information on input, output parameters printed.
    //assert(false && "verbosity level in MUMPS not recognized");
    
  } //end of switch
};

