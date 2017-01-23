/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */
/* 2015. Modified by Nai-Yuan Chiang for NLP*/

#include "DeSymIndefSolver.h"
#include "SimpleVector.h"
#include <cassert>

#include "DenseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "stdlib.h"

#ifdef WITH_MA27
#include "Ma27Solver.h"
#endif
#ifdef WITH_MA57
#include "Ma57Solver.h"
#endif
#ifdef WITH_PARDISO
#include "PardisoSolver.h"
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
#endif
#ifdef TIMING
#include <mpi.h>
#include "../../global_var.h"
#include "../PIPS-NLP/Core/Utilities/PerfMetrics.h"
#endif


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
extern "C" void FNAME(dgetrf)(int *n, 
			int *m, 
			double A[], 
			int *lda, 
			int ipiv[], 
			int *info);
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

DeSymIndefSolver::DeSymIndefSolver( DenseSymMatrix * dm )
{
  mStorage = DenseStorageHandle( dm->getStorage() );

  int size = mStorage->n;
  ipiv = new int[size];
  lwork = -1;
  work = NULL;
  sparseMat = 0;
  
}

DeSymIndefSolver::DeSymIndefSolver( SparseSymMatrix * sm )
{
  int size = sm->size();
  mStorage = DenseStorageHandle( new DenseStorage(size,size) );

  ipiv = new int[size];
  lwork = -1;
  work = NULL;
  sparseMat = sm;

  

}


//#include "mpi.h"
int DeSymIndefSolver::matrixChanged()
{

  char fortranUplo = 'U';
  int info;

  int n = mStorage->n;
#ifdef TIMING
  double stime1=MPI_Wtime();
#endif
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
#ifdef TIMING
  std::cout << "dystrf n: " << n << std::endl;
#endif
  FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
	   ipiv, &lworkNew, &lwork, &info );
  
  lwork = (int)lworkNew; 
  if(work) delete[] work;
  work = new double[lwork];


#ifdef TIMING_FLOPS
  HPM_Start("DSYTRFFact");
#endif

//*********************************************************************************
  // try to find inertia information for this dense matrix
  int nnzWrk = (1+n)*n/2;;
  int *rowStartM;  
  int *rowM = (int*) malloc (nnzWrk*sizeof(int));
  int *colM = (int*) malloc (nnzWrk*sizeof(int));
  double *elesM = (double*) malloc (nnzWrk*sizeof(double));
  int findNz=0;

  for(int rowID=0;rowID<n;rowID++){
	for(int colID=0;colID<n;colID++){
	  if(rowID>=colID){
		rowM[findNz] = rowID+1;
		colM[findNz] = colID+1;
		elesM[findNz]= mStorage->M[rowID][colID];
		findNz++;
	  }
	}
	if(gSymLinearAlgSolverForDense>1)
	  rowStartM[rowID+1] = findNz+1;
  }  
  assert(findNz==nnzWrk);
  if(gSymLinearAlgSolverForDense>1) assert(rowStartM[n]-1==nnzWrk);
//*********************************************************************************

  //factorize
  double *temp=new double[n*n];
  for(int i=0; i<n;i++) {
    for(int j=0; j<n;j++) {
      temp[i*n+j]=mStorage->M[i][j];
    }
  }
  FNAME(dgetrf)( &n, &n, temp, &n,
      ipiv, &info );
  FNAME(dsytrf)( &fortranUplo, &n, &mStorage->M[0][0], &n,
      ipiv, work, &lwork, &info );

#ifdef TIMING_FLOPS
  HPM_Stop("DSYTRFFact");
#endif
  if(info!=0)
      printf("DeSymIndefSolver::matrixChanged : error - dsytrf returned info=%d\n", info);
  //assert(info==0);

  //int piv2x2=0;
  //for(int i=0; i<n; i++)
  //  if(ipiv[i]<0) piv2x2++;
  //printf("%d 2x2 pivots were used\n", piv2x2);
#ifdef TIMING
  gprof.t_dsytrf+=MPI_Wtime()-stime1;
#endif
negEigVal=0;
for(int i=0;i<n;i++) {
  if(temp[i*n+i]<0) negEigVal++;
}
std::cout << "negEigVal1: " << negEigVal << std::endl;
//return negEigVal;

//*********************************************************************************
if(gSymLinearAlgSolverForDense<=1){
//*********************************************************************************
  // we use MA57 or MA27 to find inertia
  negEigVal=0;
  
  if(n>1){
	int     icntlTemp[30];
    double  cntlTemp[5];
    int     infoTemp[40];	
    double  rinfoTemp[20];
    double  ipessimism = 1.4;
    double  rpessimism = 1.4;

    int	   lkeepTemp;
    int	   lifactTemp, lfactTemp;
	int   *keepTemp, *iworkTemp;
	int    nsteps;
	double *factTemp;
	int *ifactTemp, *iworkTemp2;

	if(gSymLinearAlgSolverForDense==1){
#ifdef WITH_MA57	
#ifdef TIMING
double stime=MPI_Wtime();
#endif	
	  FNAME(ma57id)( cntlTemp, icntlTemp );
	  
	  icntlTemp[1-1] = 0;	// don't print warning messages
	  icntlTemp[2-1] = 0;	// no Warning messages
	  icntlTemp[4-1] = 0;	// no statistics messages
	  icntlTemp[5-1] = 0;	// no Print messages.
	  icntlTemp[6-1] = 0; 	// no messages output.
	  icntlTemp[7-1] = 1;	// 0 or 2 use MC47;  1 use the one kept from ma57ad; 3 min degree ordering as in MA27; 4 use Metis; 5 automatic choice(MA47 or Metis)
	  icntlTemp[9-1] = 10;  // up to 10 steps of iterative refinement
	  icntlTemp[11-1] = 16;
	  icntlTemp[12-1] = 16; 
	  icntlTemp[15-1] = 0;
      icntlTemp[16-1] = 0;
	  
      cntlTemp[0] = gHSL_PivotLV;
	  cntlTemp[1] = 1.e-20;	  

      lkeepTemp = ( nnzWrk > n ) ? (5 * n + 2 *nnzWrk + 42) : (6 * n + nnzWrk + 42);
      keepTemp  = new int[lkeepTemp];
      iworkTemp = new int[5 * n];
      
#ifdef TIMING
gprof.t_ma57id+=MPI_Wtime()-stime;
stime=MPI_Wtime();
#endif
	  
      FNAME(ma57ad)( &n, &nnzWrk, rowM, colM, &lkeepTemp, keepTemp, iworkTemp, icntlTemp,
	     infoTemp, rinfoTemp );	 
	  
      lfactTemp  = infoTemp[8];
      lfactTemp  = 2*(int) (rpessimism * lfactTemp );
      factTemp   = new double[lfactTemp];
      lifactTemp = infoTemp[9];
      lifactTemp = (int) (ipessimism * lifactTemp );
      ifactTemp  = new int[lifactTemp ];
      iworkTemp2 = new int[n];
	
//      char *name = "matSt.m";
//
//       FILE* outfile = fopen( name, "wr" );
//
//		fprintf(outfile,"n_dim_%s = %d;\n\n", name,n);;
//		fprintf(outfile,"nnz_%s = %d;\n\n", name, nnzWrk);
//
//		int findkk=0;
//		fprintf(outfile,"rowId_%s = [ ",name);
//		for(int kk=0;kk<nnzWrk;kk++)
//			if(findkk%10==0){
//			 fprintf(outfile,"%d, ... \n", rowM[kk]+1);
//					 findkk++;
//			}
//			else{
//			 fprintf(outfile,"%d,", rowM[kk]+1);
//					 findkk++;
//			}
//		fprintf(outfile,"]; \n\n");
//
//		findkk=0;
//		fprintf(outfile,"colId_%s = [ ",name);
//		for(int kk=0;kk<nnzWrk;kk++)
//			if(findkk%10==0){
//			 fprintf(outfile,"%d, ... \n", colM[kk]+1);
//					 findkk++;
//			}
//			else{
//			 fprintf(outfile,"%d,", colM[kk]+1);
//					 findkk++;
//			}
//		fprintf(outfile,"]; \n\n");
//
//		findkk=0;
//		fprintf(outfile,"elts_%s = [ ",name);
//		for(int kk=0;kk<nnzWrk;kk++)
//			if(findkk%10==0){
//			  fprintf(outfile,"%5.17g, ... \n", elesM[kk]);
//					 findkk++;
//			}
//			else{
//			  fprintf(outfile,"%5.17g,", elesM[kk]);
//					 findkk++;
//			}
//		fprintf(outfile,"]; \n\n");
//
//		fclose(outfile);
#ifdef TIMING
gprof.t_ma57ad+=MPI_Wtime()-stime;
stime=MPI_Wtime();
#endif

      FNAME(ma57bd)( &n,	   &nnzWrk,	elesM,	   factTemp,  &lfactTemp,  ifactTemp,
		   &lifactTemp,  &lkeepTemp,  keepTemp,  iworkTemp2,  icntlTemp,  cntlTemp,
		   infoTemp,	 rinfoTemp );
#ifdef TIMING
       gprof.t_ma57bd+=MPI_Wtime()-stime;
#endif

      if( infoTemp[0] == 0 ){
		  negEigVal = infoTemp[24-1]; 
	  }else if (infoTemp[0] == 4){
		  negEigVal = -1; 
	  }else{ 
  	  	cout << "ma57bd: Factorization Fails: info[0]=: " << infoTemp[0] << endl;
//        assert(false);
      }
  
#else
	  assert("ma57 not defined"&&0);
#endif
	}
	else{
#ifdef WITH_MA27		
	  FNAME(ma27id)( icntlTemp,cntlTemp );
	  
	  icntlTemp[1-1] = 0;	// don't print warning messages
	  icntlTemp[2-1] = 0;	// no Warning messages
	  
      cntlTemp[0] = gHSL_PivotLV;
	  cntlTemp[1] = 1.e-20;	  

	  // set array lengths as recommended in ma27 docs
  	  lkeepTemp = (int)(1.3 * (2*nnzWrk + 3*n + 1));
  	  keepTemp = new int[lkeepTemp];

  	  // define ikeep (which stores the pivot sequence)
  	  int *ikeepTemp = new int[3*n];
  	  iworkTemp = new int[2*n];
	  
 	  // set iflag to zero to make ma27ad choose a pivot order.
  	  int iflag = 0, maxfrt;

  	  double ops;
  	  FNAME(ma27ad)( &n, &nnzWrk, rowM, colM, keepTemp, &lkeepTemp, ikeepTemp, iworkTemp, &nsteps, &iflag,
	   		icntlTemp, cntlTemp, infoTemp, &ops);

  	  lfactTemp  = infoTemp[4];
  	  factTemp	 = new double[lfactTemp];
  	  lifactTemp = infoTemp[5];
  	  ifactTemp  = new int[lifactTemp];	  
  	  iworkTemp2 = new int[n];

  	  for (int i=0; i<nnzWrk; i++) factTemp[i] = elesM[i];
  	  for (int i = nnzWrk; i<lfactTemp; i++) factTemp[i] = 0.0;	  
	 
      FNAME(ma27bd)( &n, &nnzWrk, rowM, colM, factTemp, &lfactTemp, ifactTemp, &lifactTemp, ikeepTemp,
	     &nsteps,  &maxfrt, iworkTemp2, icntlTemp, cntlTemp,  infoTemp );  
	  
      if( infoTemp[0] == 0 ){
		  negEigVal = infoTemp[15-1]; 
	  }else if (infoTemp[0] == 4){
		  negEigVal = -1; 
	  }else{ 
  	  	cout << "ma27bd: Factorization Fails: info[0]=: " << infoTemp[0] << endl;
//        assert(false);
      }	
	  
	  delete [] ikeepTemp ;
#else
	  assert("ma27 not defined"&&0);
#endif		
	}

	delete [] iworkTemp2 ;
	delete [] ifactTemp ;
	delete [] factTemp ;
	delete [] iworkTemp ;
	delete [] keepTemp ;
  }
  else{
	assert(n==1);
	if(mStorage->M[0][0]<0) negEigVal=1;	
  }
  
//*********************************************************************************
}
else{
//*********************************************************************************
#ifdef WITH_PARDISO
  // we use Pardiso to find inertia
  void  *pt[64]; 
  int iparm[64];
  int num_threads;
  double dparm[64];

  /* Numbers of processors, value of OMP_NUM_THREADS */
  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
    sscanf( var, "%d", &num_threads );
  else {
    printf("Set environment OMP_NUM_THREADS");
    exit(1);
  }  

  negEigVal=0;

  if(n>1){
    int solverInfo=0, mtype=-2, error;
    pardisoinit (pt,  &mtype, &solverInfo, iparm, dparm, &error); 
    if (error!=0) {
      cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << endl;
      assert(false);
    }
  
    int phase = 12; //Analysis, numerical factorization
    int maxfct=1; //max number of fact having same sparsity pattern to keep at the same time
    int mnum=1; //actual matrix (as in index from 1 to maxfct)
    int nrhs=1;
    int msglvl=0; //messaging level
    int error2,  mtype2=-2;

    iparm[2] = num_threads;
    iparm[1] = 2; // 2 is for metis, 0 for min degree 
    iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
    iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1; 
                 // if needed, use 2 for advanced matchings and higer accuracy.
    iparm[30] = 0; // do not specify sparse rhs at this point

    pardiso (pt , &maxfct , &mnum, &mtype2, &phase,
	   &n, elesM, rowStartM, colM,
	   NULL, &nrhs,
	   iparm, &msglvl, NULL, NULL, &error2, dparm );
    if ( error2 != 0) {
      printf ("PardisoSolver - ERROR during factorization: %d\n", error2 );
     assert(false);
    }

	negEigVal = iparm[22];
  }
  else{
	assert(n==1);
	if(mStorage->M[0][0]<0)
	  negEigVal=1;	
  }

  free(rowStartM);
#else
	  assert("pardiso not defined"&&0);
#endif  
}
//*********************************************************************************

  free (rowM);
  free (colM);
  free (elesM);

std::cout << "negEigVal2: " << negEigVal << std::endl;
  return negEigVal;

}

void DeSymIndefSolver::solve ( OoqpVector& v )
{
  char fortranUplo = 'U';
  int info;
  int one = 1;

  int n = mStorage->n; SimpleVector &  sv = dynamic_cast<SimpleVector &>(v);
#ifdef TIMING_FLOPS
  HPM_Start("DSYTRSSolve");
#endif

  FNAME(dsytrs)( &fortranUplo, &n, &one,	&mStorage->M[0][0],	&n,
	   ipiv, &sv[0],	&n,	&info);

#ifdef TIMING_FLOPS
  HPM_Stop("DSYTRSSolve");
#endif
  assert(info==0);
}

void DeSymIndefSolver::solve ( GenMatrix& rhs_in )
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

void DeSymIndefSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

DeSymIndefSolver::~DeSymIndefSolver()
{
  delete[] ipiv;
  if(work) delete[] work;
}
