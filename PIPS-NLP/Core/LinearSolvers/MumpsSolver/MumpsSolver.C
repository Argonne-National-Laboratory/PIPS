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

static bool gMUMPSStatsOn=true;

using namespace std;


MumpsSolver::MumpsSolver(const long long& globSize, MPI_Comm mumpsMpiComm, MPI_Comm pipsMpiComm)
{
  int error=MPI_Comm_rank(mumpsMpiComm, &my_rank_); assert(error==MPI_SUCCESS);

  n_ = globSize;
  mumps_ = new DMUMPS_STRUC_C;
  mumps_->comm_fortran = getFortranMPIComm(mumpsMpiComm);
  mumps_->sym = 2;//general symetric matrix
  mumps_->job = -1; //initialization
  mumps_->par = 1; //host process involved in parallel computations
  dmumps_c(mumps_);  

  if(mumps_->infog[1-1] != 0)
    if(my_rank_==0) printf("Error occured when initializing MUMPS.\n");

  setMumpsVerbosity(1);//3 moderately verbose, 1 only errors, 0 anything is surpressed

  mumps_->n = (int)globSize;
  mumps_->icntl[ 5  -1] = 0; //assembled format for matrix (triplet format)
  mumps_->icntl[18  -1] = 3; //distributed assembled format
  mumps_->icntl[20  -1] = 0; //the right-hand side is in dense format in the structure component
  mumps_->icntl[21  -1] = 0; // the solution vector is assembled and stored in the structure component RHS
  //ICNTL(29) defines the parallel ordering tool to be used to compute the fill-in reducing permutation.
}

MumpsSolver::~MumpsSolver()
{
  mumps_->job = -2; //terminate mumps
  dmumps_c(mumps_);
  delete mumps_;

}
bool MumpsSolver::setLocalEntries(long long globnnz, long long locnnz, int* locirn, int* locjcn, double* locA)
{
  mumps_->nnz = globnnz;
  mumps_->nnz_loc = locnnz;
  //par%IRN par%JCN par%A
  mumps_->irn_loc = locirn;
  mumps_->jcn_loc = locjcn;
  mumps_->a_loc = locA;

  return true;
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
  mumps_->n = n_;

  //printf("[1] nnz %d nnz_loc %d nelt %d\n", mumps_->nnz, mumps_->nnz_loc, mumps_->nelt);

  mumps_->job = 1; //perform the analysis

  //ICNTL(28) determines whether a sequential (=1) or a parallel analysis (=2) is performed. For the latter case
  //ICNTL(7) is meaningless. Automatic choice would be indicated by 0
  mumps_->icntl[28-1] = 2; 

  //ICNTL(29) defines the parallel ordering tool to be used to compute the fill-in reducing permutation.
  // 1 for PT_SCOTCH, 2 for ParMetis, 0 automatic (default)
  //mumps_->icntl[29-1] =2;

  // 0: parallel factorization of the root node (default)
  // > 0: forces a sequential factorization of the root node (ScaLAPACK will not be used).
  // In this case if the number of working processors is strictly larger than ICNTL(13) then 
  // splitting of the root node is performed, in order to automatically recover part of the 
  // parallelism lost because the root node was processed sequentially (advised value is 1).
  //
  // see also CNTL(1): CNTL(1) is the relative threshold for numerical pivoting.
  //
  //petra: note that parallel factorization of the root nodes results in incorrect inertia calculation
  //petra: because of an issue with ScaLAPACK. 
  //petra: thus, use with > 0, in particular the advised value seem to work fine, but different hardware may
  //petra: give different performance
  mumps_->icntl[13-1] = 1;

  // CNTL(4) determines the threshold for static pivoting. 0.0 will activate it and MUMPS computes the threshold
  // automatically. Negative values will disable it (default). 
  // static pivoting requires iterative refinement.
  //mumps_->cntl[4-1] = 0.0;

  double tm = MPI_Wtime();

  dmumps_c(mumps_);
  int error = mumps_->infog[1-1];
  if(error != 0) {
    if(my_rank_==0) printf("Error INFOG(1)=%d occured in Mumps in the analysis phase. \n", error);
  }

  tm = MPI_Wtime()-tm;
  if(gMUMPSStatsOn)
    if(my_rank_==0) printf("MUMPS analysis phase took %g seconds.\n", tm);

  //printf("[2] nnz %d nnz_loc %d nelt %d\n", mumps_->nnz, mumps_->nnz_loc, mumps_->nelt);

  // INFOG(7) - after analysis: The ordering method actually used. The returned value will depend on
  // the type of analysis performed, e.g. sequential or parallel (see INFOG(32)). Please refer to
  // ICNTL(7) and ICNTL(29) for more details on the ordering methods available in sequential and
  // parallel analysis respectively.
  // 1 scotch, 2 metis
  if(gMUMPSStatsOn)
    if(my_rank_==0) printf("MUMPS ordering actually used %d\n", mumps_->infog[7-1]);
  
  // INFOG(32) - after analysis: The type of analysis actually done (see ICNTL(28)). 
  // 1 sequential, 2 parallel
  if(gMUMPSStatsOn)
    if(my_rank_==0) printf("MUMPS ordering parallelism %d\n", mumps_->infog[32-1]);
  //
  //factorize
  //
  mumps_->job = 2;

  tm = MPI_Wtime();

  dmumps_c(mumps_);
  error = mumps_->infog[1-1];
  if(error != 0) {
    if(my_rank_==0) printf("Error INFOG(1)=%d occured in Mumps in the factorization phase. \n", error);
  }
  //printf("[3] nnz %d nnz_loc %d nelt %d\n", mumps_->nnz, mumps_->nnz_loc, mumps_->nelt);
  saveOrderingPermutation();

  //CNTL(1) is the relative threshold for numerical pivoting.
  //CNTL(4) determines the threshold for static pivoting. See Subsection 3.9
  //INFOG(12) - after factorization: Total number of off-diagonal pivots
  
  //if ScaLAPACK is allowed for the last dense block (default in parallel, ICNTL(13)=0), then the
  //presence of negative pivots in the part of the factorization processed with ScaLAPACK (subroutine
  //P POTRF) will raise an error and the code -40 is then returned in INFOG(1);
  if (error == -8 || error == -9) {
    //not enough memory

    //assert(false && "this code was not fully tested. ");
    const int try_max=10;
    for(int tr=0; tr<try_max; tr++) {
      double mem_percent = mumps_->icntl[13];
      mumps_->icntl[13] = 2 * mumps_->icntl[13];

      if(my_rank_==0) printf("Increased Mumps ICNTL(14) from %g to %d\n", mem_percent, mumps_->icntl[13]);

      dmumps_c(mumps_);
      error = mumps_->infog[1-1];

      saveOrderingPermutation();

      if (error != -8 && error != -9)
	break;
    }

    if (error == -8 || error == -9) {
      if(my_rank_==0) printf("Fatal error in Mumps: not able to obtain more memory (ICNTL(14)-related)\n");
      return -1;
    }
  }

  if (error == -10) {
    //system is singular
    if(my_rank_==0) printf("Warning: Mumps INFO(1) = %d matrix is singular.\n", error);
    return -1;
  }

  tm = MPI_Wtime()-tm;
  if(gMUMPSStatsOn)
    if(my_rank_==0) printf("MUMPS numerical factorization took %g seconds.\n", tm);

  if(gMUMPSStatsOn)
    if(my_rank_==0) printf("Order of the largest frontal matrix: %d\n", mumps_->infog[11-1]);

  int negEigVal = mumps_->infog[12-1];

  if(gMUMPSStatsOn) if(my_rank_==0)  printf("Mumps says matrix has %d negative eigenvalues\n", negEigVal);
  return negEigVal;
}

int MumpsSolver::saveOrderingPermutation() {
  assert(mumps_->job >= 2 && "MUMPS's state must be after analysis.");

  //save ordering computed by MUMPS to disk in ordering.txt. The ordering is on host processor
  //mumps par%SYM PERM
  const char* filename = "ordering.txt";
  if(my_rank_ == 0) {
    FILE* f = fopen(filename, "w");
    if(NULL!=f) {
      for(int it=0; it<n_; it++)
	fprintf(f, "%d\n", mumps_->sym_perm[it]);
      fclose(f);
    } else {
      printf("Could not open files %s to write the ordering permutation.\n", filename);
    }
  }
}

// Centralized solution (ICNTL(21)=0)
// Distributed solution (ICNTL(21)=1)
void  MumpsSolver::solve ( double* vec )
{
  mumps_->rhs = vec;
  //solve
  mumps_->job = 3;

#define MAX_ITER_REFIN 3

  //iterative refinement
  mumps_->icntl[10-1] = MAX_ITER_REFIN; //maximum number of iterative refinements

  // compute main statistics on the error: 0 disabled, 1 full stats (very expensive), 2 main stats
  // in addition to stats computed for 2 (see below), a value of 1 also computes
  // an estimate for the (forward) error in the solution in RINFOG(9), and condition numbers for the linear 
  // system in RINFOG(10) and RINFOG(11) are also returned
  mumps_->icntl[11-1] = 2;

  double t = MPI_Wtime();

  dmumps_c(mumps_);
  int error = mumps_->infog[1-1];
  if(error != 0) {
    if(my_rank_==0) printf("Error INFOG(1)=%d occured in Mumps in the solve phase. \n", error);
    assert(false);
    return;
  }

  t = MPI_Wtime()-t;
  if(gMUMPSStatsOn)
    if(my_rank_==0) printf("MUMPS solve  took %g seconds.\n", t);

  //
  //error analysis ICNTL(11)
  //
  //the infinite norm of the matrix
  double Ainfnrm = mumps_->rinfog[4-1];

  //the infinite norm of the computed solution
  double xinfnrm = mumps_->rinfog[5-1];

  //the scaled residual
  double relresid = mumps_->rinfog[6-1];

  //componentwise backward error estimates
  double omega1=mumps_->rinfog[7-1], omega2=mumps_->rinfog[8-1];

  if(gMUMPSStatsOn)
    if(my_rank_==0) 
      printf("MUMPS solution: backward errors %g %g  scaled resid %g Ainfnorm %g  xinfnorm %g.\n", 
	     omega1, omega2, relresid, Ainfnrm, xinfnrm);

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
    break;
  }
  case 2: {
    mumps_->icntl[2 -1] = 6; //output stream for diagnostic printing, statistics, and warning messages.
    mumps_->icntl[4 -1] = 2; //Errors, warnings, and main statistics printed.
    break;
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

