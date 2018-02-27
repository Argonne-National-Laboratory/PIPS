/* PIPS file - Written by Cosmin Petra, 2018 */

#include "../../global_var.h"
#include <omp.h>
#include "MumpsSolver.h"
#include "SimpleVector.h"
#include <cassert>
#include <cmath>
#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "stdlib.h"


// #ifdef WITH_MA27
// #include "Ma27Solver.h"
// #endif
// #ifdef WITH_MA57
// #include "Ma57Solver.h"
// #endif
// #ifdef WITH_PARDISO
// #include "PardisoSolver.h"
// extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
// extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
//                   double *, int    *,    int *, int *,   int *, int *,
//                      int *, double *, double *, int *, double *);
// #endif

extern double gHSL_PivotLV;
extern int gSymLinearAlgSolverForDense;


#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

// declarations for LAPACK functions used to factor/solve:

// dsytrf_() factors a symmetric indefinite matrix A, see LAPACK 
// documentation for more details.
extern "C" void FNAME(dsytrf)(char *uplo, 
			int *n, 
			double A[], 
			int *lda, 
			int ipiv[], 
			double work[],
			int *lwork, 
			int *info);

// dsytrs_() solves the system Ax = b using the factor obtained by dsytrf_().
extern "C" void FNAME(dsytrs)(char *uplo, 
			int *n, 
			int *nrhs, 
			double A[], 
			int *lda, 
			int ipiv[], 
			double b[], 
			int *ldb,
			int *info);

#ifdef TIMING_FLOPS
extern "C" {
    void HPM_Init(void);
    void HPM_Start(char *);
    void HPM_Stop(char *);
    void HPM_Print(void);
    void HPM_Print_Flops(void);
    void HPM_Print_Flops_Agg(void);
    void HPM_Terminate(char*);
}
#endif

using namespace std;

MumpsSolver::MumpsSolver( DenseSymMatrix * dm )
{
  mStorage = DenseStorageHandle( dm->getStorage() );

  DMUMPS_STRUC_C* mumps_ = new DMUMPS_STRUC_C;


  int size = mStorage->n;
  ipiv = new int[size];
  lwork = -1;
  work = NULL;
  sparseMat = 0;
  
}

MumpsSolver::MumpsSolver( SparseSymMatrix * sm )
{


  int size = sm->size();
  mStorage = DenseStorageHandle( new DenseStorage(size,size) );

  ipiv = new int[size];
  lwork = -1;
  work = NULL;
  sparseMat = sm;
}

MumpsSolver::~MumpsSolver()
{
  DMUMPS_STRUC_C* mumps_ = (DMUMPS_STRUC_C*)mumps_ptr_;
  mumps_->job = -2; //terminate mumps
  dmumps_c(mumps_);
  delete [] mumps_->a;
  delete mumps_;

  delete[] ipiv;
  if(work) delete[] work;
}


//#include "mpi.h"
int MumpsSolver::matrixChanged()
{
  int n = mStorage->n;
  /* If matrix is dense use LAPACK, not MA57 */
  char fortranUplo = 'U';
  int info;
  
  if (sparseMat) {
    std::fill(mStorage->M[0],mStorage->M[0]+n*n,0.);

    const double *sM = sparseMat->M();
    const int *jcolM = sparseMat->jcolM();
    const int *krowM = sparseMat->krowM();
    for (int i = 0; i < n; i++) {
      for (int k = krowM[i]; k < krowM[i+1]; k++) {
        int col = jcolM[k];
        mStorage->M[i][col] = sM[k];
      }
    }
  }

  //query the size of workspace
  lwork=-1;
  double lworkNew;
  FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
	   ipiv, &lworkNew, &lwork, &info );
  
  lwork = (int)lworkNew; 
  if(work) delete[] work;
  work = new double[lwork];


#ifdef TIMING_FLOPS
  HPM_Start("DSYTRFFact");
#endif

  //factorize
  //if(gmyid==0) {
    //printf("OMP_NUM_PROCS: %d\n", omp_get_num_procs());
    //printf("OMP_GET_MAX_THREADS: %d\n", omp_get_max_threads());
    //printf("On node ranks: %d\n", gnprocs_node);
  //}

  /* Create MPI windows used for the scatter of the factorization result. One
   * MPI rank does the factorization and distributes the result through the
   * windows
   */

  if(gwindow==NULL) {
    if(gmyid_node==0) {
      gwindow=new double[n*n];
      gipiv=new int[n+1];
      MPI_Win_create(gwindow, n*n*sizeof(double), sizeof(double), MPI_INFO_NULL, comm_node, &gwin);
      MPI_Win_create(gipiv, (n+1)*sizeof(int), sizeof(int), MPI_INFO_NULL, comm_node, &gwin_ipiv);
    } else {
      MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, comm_node, &gwin); 
      MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, comm_node, &gwin_ipiv); 
      // Allocate something to overwrite NULL
      gwindow=new double[1];
    }
  }
#ifdef DUMP
  static bool dump=true;
  if(gmyid==0 && dump==true) {
    printf("Dumping 1stageM.dmp\n");
    FILE *fp=fopen("1ststageM.dmp","w");
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&mStorage->M[0][0], sizeof(double), n*n, fp);
    fclose(fp);
    dump=false;
  } 
#endif
  if(gmyid_node==0) {
    /* rank 0 does the factorization */
    FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
      ipiv, work, &lwork, &info );
    //  printf("gwindow on 0 before: %lf %d %d\n", gwindow[0],gipiv[0],info);
    for(int i=0;i<n*n;i++) gwindow[i]=(&mStorage->M[0][0])[i];
    for(int i=0;i<n;i++) gipiv[i]=ipiv[i];
    gipiv[n]=info;
    //  printf("gwindow on 0 after: %lf %d %d\n", gwindow[0],gipiv[0],info);
    MPI_Win_fence(0,gwin);
    MPI_Win_fence(0,gwin_ipiv);
    MPI_Win_fence(0,gwin);
    MPI_Win_fence(0,gwin_ipiv);
  }
  else {
    /* All other ranks get the result from rank 0 */
    //printf("gwindow on other before: %lf %d %d\n", mStorage->M[0][0], ipiv[0], info);
    MPI_Win_fence(0,gwin);
    MPI_Win_fence(0,gwin_ipiv);
    MPI_Get(&mStorage->M[0][0], n*n, MPI_DOUBLE, 0, 0, n*n, MPI_DOUBLE, gwin);
    MPI_Get(ipiv, n, MPI_INT, 0, 0, n, MPI_INT, gwin_ipiv);
    MPI_Get(&info, 1, MPI_INT, 0, n, 1, MPI_INT, gwin_ipiv);
    MPI_Win_fence(0,gwin);
    MPI_Win_fence(0,gwin_ipiv);
    //printf("gwindow on other after: %lf %d %d\n", mStorage->M[0][0], ipiv[0], info);
  }

#ifdef TIMING_FLOPS
  HPM_Stop("DSYTRFFact");
#endif
  if(info!=0)
      printf("MumpsSolver::matrixChanged : error - dsytrf returned info=%d\n", info);
/* Compute the inertia. Only negative eigenvalues are returned */
negEigVal=0;
double t=0;
for(int k=0; k<n; k++) {
  double d = mStorage->M[k][k];
  if(ipiv[k] < 0) {
   if(t==0) {
     t=fabs(mStorage->M[k+1][k]);
     d=(d/t)*mStorage->M[k+1][k+1]-t;
   }
   else {
     d=t;
     t=0;
   }
  }
  if(d<0) negEigVal++;
  if(d==0) {
    negEigVal=-1;
    break;
  }
}
return negEigVal;
}

void MumpsSolver::solve ( OoqpVector& v )
{
  char fortranUplo = 'U';
  int info;
  int one = 1;

  int n = mStorage->n; SimpleVector &  sv = dynamic_cast<SimpleVector &>(v);
#ifdef TIMING_FLOPS
  HPM_Start("DSYTRSSolve");
#endif
#ifdef DUMP
  static bool dump=true;
  if(gmyid==0 && dump==true) {
    printf("Dumping 1stageRHS.dmp\n");
    FILE *fp=fopen("1ststageRHS.dmp","w");
    int ione=1;
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&ione, sizeof(int), 1, fp);
    fwrite(&sv[0], sizeof(double), n, fp);
    fclose(fp);
    dump=false;
  } 
#endif

  FNAME(dsytrs)( &fortranUplo, &n, &one,	&mStorage->M[0][0],	&n,
	   ipiv, &sv[0],	&n,	&info);

#ifdef TIMING_FLOPS
  HPM_Stop("DSYTRSSolve");
#endif
#ifdef DUMP
  static bool dump1=true;
  if(gmyid==0 && dump1==true) {
    printf("Dumping 1stageSol.dmp\n");
    FILE *fp=fopen("1ststageSol.dmp","w");
    int ione=1;
    fwrite(&ione, sizeof(int), 1, fp);
    fwrite(&n, sizeof(int), 1, fp);
    fwrite(&sv[0], sizeof(double), n, fp);
    fclose(fp);
    dump1=false;
  } 
#endif
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

void MumpsSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}
