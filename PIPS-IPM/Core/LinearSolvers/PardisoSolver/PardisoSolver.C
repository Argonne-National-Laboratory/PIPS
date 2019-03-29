/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin, Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <stdlib.h>
#include <iostream>
using namespace std;

#include "PardisoSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include "pipsdef.h"
#include <cstdlib>

#include "mpi.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#ifdef WITH_MKL_PARDISO
#include "mkl_pardiso.h"
#include "mkl_types.h"
#endif

extern int gOoqpPrintLevel;

#ifndef WITH_MKL_PARDISO
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *, 
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);


extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);
#endif

PardisoSolver::PardisoSolver( SparseSymMatrix * sgm )
{
#ifdef TIMING
    int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if(myRank==0)
	cout << "PardisoSolver::PardisoSolver (sparse input)" << endl;
#endif
  Msys = sgm;
  n = sgm->size();
  nnz=sgm->numberOfNonZeros();

  krowM = new int[n+1];
  jcolM = new int[nnz];
  M = new double[nnz];

  sol = NULL;
  sz_sol = 0;

  first = true;
  nvec=new double[n];

#ifndef WITH_MKL_PARDISO
  num_threads = PIPSgetnOMPthreads();
#endif
  Mdsys=NULL;
}


PardisoSolver::PardisoSolver( DenseSymMatrix * m )
{
#ifdef TIMING
    int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if(myRank==0)
	cout << "PardisoSolver created (dense input)" << endl;
#endif
  Msys = NULL;
  Mdsys = m;
  n = m->size();
  
  nnz=0; //getNumberOfNonZeros(*Mdsys);

  krowM = NULL;//new int[n+1];
  jcolM = NULL;//new int[nnz];
  M = NULL;//new double[nnz];

  sol = NULL;
  sz_sol = 0;

  first = true;
  nvec=new double[n];
   
#ifndef WITH_MKL_PARDISO
  /* Numbers of processors, value of OMP_NUM_THREADS */
  char *var = getenv("OMP_NUM_THREADS");
  if(var != NULL)
    sscanf( var, "%d", &num_threads );
  else {
    printf("Set environment OMP_NUM_THREADS");
    exit(1);
  }
#endif
}
int PardisoSolver::getNumberOfNonZeros(DenseSymMatrix& m)
{
  int nnz=0;
  for(int i=0; i<n; i++) {
    for(int j=i+1; j<n; j++) 
      if(m[i][j]!=0.0) 
	nnz++;
    nnz++; //always have diags
  }
  return nnz;
}

void
PardisoSolver::firstCall()
{

   int mtype = -2;
#ifndef WITH_MKL_PARDISO
   int solver = 0, error;

   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

   if( error != 0 )
   {
      cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << endl;
      assert(false);
   }
#else
   /* enable matrix checker (default disabled) - mkl pardiso does not have chkmatrix */
   //iparm[26] = 1;
   pardisoinit(pt, &mtype, iparm);
#endif

   if( Msys )
   {

      //get the matrix in upper triangular
      Msys->getStorageRef().transpose(krowM, jcolM, M);

      //save the indeces for diagonal entries for a streamlined later update
      int* krowMsys = Msys->getStorageRef().krowM;
      int* jcolMsys = Msys->getStorageRef().jcolM;
      for( int r = 0; r < n; r++ )
      {
         // Msys - find the index in jcol for the diagonal (r,r)
         int idxDiagMsys = -1;
         for( int idx = krowMsys[r]; idx < krowMsys[r + 1]; idx++ )
            if( jcolMsys[idx] == r )
            {
               idxDiagMsys = idx;
               break;
            }
         assert(idxDiagMsys>=0);
         //must have all diagonals entries

         // aug  - find the index in jcol for the diagonal (r,r)
         int idxDiagAug = -1;
         for( int idx = krowM[r]; idx < krowM[r + 1]; idx++ )
            if( jcolM[idx] == r )
            {
               idxDiagAug = idx;
               break;
            }
         assert(idxDiagAug>=0);

         diagMap.insert(pair<int, int>(idxDiagMsys, idxDiagAug));
      }
   }
   else
   {
      // the input is a dense matrix
      // the dense matrix is also processed in matrixChanged everytime the method is called
      if( Mdsys )
      {
         // the input is a dense matrix
         DenseSymMatrix& Md = (*Mdsys);
         nnz = getNumberOfNonZeros(Md);

         delete[] krowM;
         delete[] jcolM;
         delete[] M;
         krowM = new int[n + 1];
         jcolM = new int[nnz];
         M = new double[nnz];

         int nzIt = 0;
         for( int i = 0; i < n; i++ )
         {

            krowM[i] = nzIt;
            //cout << i << " " << krowM[i] << endl;

            jcolM[nzIt] = i; // the diagonal
            M[nzIt] = Md[i][i];
            assert(nzIt<nnz);
            nzIt++;

            for( int j = i + 1; j < n; j++ )
            {
               if( Md[i][j] != 0.0 )
               {
                  jcolM[nzIt] = j;
                  M[nzIt] = Md[i][j];
                  assert(nzIt<nnz);
                  nzIt++;
               }
            }
         }
         //cout << "PardisoSolver::first call nzit=" << nzIt << endl;
         assert(nzIt==nnz);
         krowM[n] = nzIt;

      }
      else
      {
         assert(false);
      }
   }

   // need Fortran indexes
   for( int i = 0; i < n + 1; i++ )
      krowM[i] += 1;
   for( int i = 0; i < nnz; i++ )
      jcolM[i] += 1;


   //
   // symbolic analysis
   //
   /* same for mkl_pardiso and pardiso */
//   mtype = -2;
//   int phase = 11; // analysis
//   int maxfct = 1, mnum = 1, nrhs = 1;
//
//   iparm[1] = 2; // 2 is for metis, 0 for min degree
//   iparm[7] = 3; // # iterative refinements
//   //iparm[9] =10; // pivot perturbation 10^{-xxx}
//   iparm[23] = 1; //Parallel Numerical Factorization (0 = used in the last years/INTEL calls it classical algorithm , 1 = two-level scheduling)
//   int msglvl = 0; // with statistical information
//
//#ifndef WITH_MKL_PARDISO
//   iparm[2] = num_threads;
//   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
//   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
//                  // if needed, use 2 for advanced matchings and higer accuracy.
//#else
//   int error;
//   /* From INTEL (instead of iparm[2] which is not defined there):
//   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
//   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
//   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
//   */
//   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
//   iparm[12] = 1; // 0 disable matching, 1 enable matching, no other settings
//#endif
//
//
//   /* no dparm in MKL PARDISO anymore */
//   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, NULL, &nrhs,
//         iparm, &msglvl, NULL, NULL, &error
//#ifndef WITH_MKL_PARDISO
//         ,dparm
//#endif
//   );
//
//   if( error != 0 )
//   {
//      printf("PardisoSolver - ERROR during symbolic factorization: %d\n",
//            error);
//      assert(false);
//   }

} 

 
void PardisoSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

extern double g_iterNumber;
void PardisoSolver::matrixChanged()
{
   if( first )
   {
      firstCall();
      first = false;
   }
   else
   {
      if( Msys )
      {
         //update diagonal entries in the PARDISO aug sys (if the input is sparse)
         double* eltsMsys = Msys->getStorageRef().M;
         map<int, int>::iterator it;
         for( it = diagMap.begin(); it != diagMap.end(); it++ )
            M[it->second] = eltsMsys[it->first];
      }
   }

   if( Mdsys )
   {
      // the input is a dense matrix
      DenseSymMatrix& Md = (*Mdsys);
      //double tm=MPI_Wtime();
      int nzIt = 0;
      for( int i = 0; i < n; i++ )
      {

         krowM[i] = nzIt;

         jcolM[nzIt] = i; // the diagonal
         M[nzIt] = Md[i][i];
         nzIt++;

         for( int j = i + 1; j < n; j++ )
         {
            if( Md[i][j] != 0.0 )
            {
               jcolM[nzIt] = j;
               M[nzIt] = Md[i][j];
               nzIt++;
            }
         }
      }
      //cout << "PardisoSolver::matrix changed nnzit=" << nzIt << endl;
      krowM[n] = nzIt;
      assert(nzIt==nnz);
      // need Fortran indexes
      for( int i = 0; i < n + 1; i++ )
         krowM[i] += 1;
      for( int i = 0; i < nnz; i++ )
         jcolM[i] += 1;
      //cout << "Forming the matrix took:" << MPI_Wtime()-tm << endl;
   }

   //
   // numerical factorization & symb. analysis
   //
   /* same for mkl_pardiso and pardiso */
   int mtype = -2, error;
   int phase = 12; // Analysis, numerical factorization
   int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
   int mnum = 1; // actual matrix (as in index from 1 to maxfct)
   int nrhs = 1;
   int msglvl = 0; // messaging level (0 = no output, 1 = staticstical info to screen)
   iparm[1] = 2; // 2 is for metis, 0 for min degree
   iparm[30] = 0; // do not specify sparse rhs at this point

#ifndef WITH_MKL_PARDISO
   iparm[2] = num_threads;
   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
                  // if needed, use 2 for advanced matchings and higer accuracy.
#else
   /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
   iparm[12] = 1; // 0 disable matching, 1 enable matching, no other settings
#endif


#ifdef PARDISO_PARALLEL_AGGRESSIVE
   iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

   #ifndef WITH_MKL_PARDISO
   iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   #else
   iparm[24] = 2; // not sure but two seems to be the appropriate equivalent here
                  // one rhs -> parallelization, multiple rhs -> parallel forward backward subst
   #endif
#else
   iparm[23] = 0; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
#endif

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, NULL, &nrhs,
         iparm, &msglvl, NULL, NULL, &error
#ifndef WITH_MKL_PARDISO
         ,dparm
#endif
   );
   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during factorization: %d\n", error);
      assert(false);
   }
}

extern int gLackOfAccuracy; 
void PardisoSolver::solve( OoqpVector& rhs_in )
{
   SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
   double * sol_local = nvec;

   //int maxRefinSteps=(gLackOfAccuracy==0?3:6);

   /* same for mkl_pardiso and pardiso */
   int phase = 33; // solve and iterative refinement
   int mtype = -2, error;
   int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
   int mnum = 1; // actual matrix (as in index from 1 to maxfct)
   int nrhs = 1;
   int msglvl = 0;
   iparm[1] = 2; // 2 is for metis, 0 for min degree
   //iparm[5] = 1; // replace drhs with the solution
   iparm[7] = 8; // max number of iterative refinement steps

#ifndef WITH_MKL_PARDISO
   iparm[2] = num_threads;
   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
                  // if needed, use 2 for advanced matchings and higer accuracy.
#else
   /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
   iparm[12] = 1; // 0 disable matching, 1 enable matching, no other settings
#endif

#ifdef PARDISO_PARALLEL_AGGRESSIVE
   iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)

   #ifndef WITH_MKL_PARDISO
   iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   #else
   iparm[24] = 2; // not sure but two seems to be the appropriate equivalent here
                  // one rhs -> parallelization, multiple rhs -> parallel forward backward subst
   #endif
#else
   iparm[23] = 0; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
#endif


   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM,
         NULL, &nrhs, iparm, &msglvl, rhs.elements(), sol_local, &error
#ifndef WITH_MKL_PARDISO
         ,dparm
#endif
   );

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }
   //iparm[6] //Number of performed iterative refinement steps.
   /*  SimpleVector r(rhs.length());
    r.copyFrom(rhs);
    if(Mdsys) {
    Mdsys->mult(1.0, r.elements(), 1, -1.0, sol,1);
    } else {
    Msys->mult(1.0, r.elements(), 1, -1.0, sol,1);
    }
    if(r.twonorm()/rhs.twonorm()>1e-8) {
    cout << "!!!PardisoSolver - rel resid=" << r.twonorm()/rhs.twonorm() << endl;
    cout << "PardisoSolver - Iter ref step " << iparm[6] << endl;
    }
    */
   /*
    SimpleVector x(n), dx(n), res(n);
    x.setToZero();       res.copyFrom(rhs);    
    int refinSteps=0;
    do {
    iparm[7] = 0; // Max numbers of iterative refinement steps .
    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
    &n, M, krowM, jcolM,
    NULL, &nrhs ,
    iparm, &msglvl,
    res.elements(), dx.elements(),
    &error, dparm );
    x.axpy(1.0,dx);

    res.copyFrom(rhs);
    if(Mdsys) {
    Mdsys->mult(1.0, res.elements(), 1, -1.0, x.elements(),1);
    } else {
    Msys->mult(1.0, res.elements(), 1, -1.0, x.elements(),1);
    }
    cout << "!!!after " << refinSteps << " steps - rel resid=" << res.twonorm()/rhs.twonorm() << endl;
    refinSteps++;
    }while(refinSteps<=8);
    rhs.copyFrom(x);
    } else */
   rhs.copyFromArray(sol_local);
}



void PardisoSolver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  //cout << "Multiple dense rhs " << endl;

  int nrows,ncols; rhs.getSize(ncols,nrows);
  if(sz_sol<nrows*ncols) {
    sz_sol=nrows*ncols;
    if(sol) delete[] sol;
    sol = new double[sz_sol];
  }
  assert(nrows==n);

  /* same for mkl_pardiso and pardiso */
  int phase = 33; // solve and iterative refinement
  int mtype = -2, error;
  int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
  int mnum = 1; // actual matrix (as in index from 1 to maxfct)
  int nrhs = ncols;
  int msglvl = 0;
  iparm[1] = 2; // 2 is for metis, 0 for min degree
  //iparm[5] = 1; // replace drhs with the solution
  iparm[7] = 2; // max number of iterative refinement steps
  iparm[30]=0; // no sparse rhs

#ifndef WITH_MKL_PARDISO
  iparm[2] = num_threads;
#else
  /* From INTEL (instead of iparm[2] which is not defined there):
  *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
  *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
  *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
  */
#endif

  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, M, krowM, jcolM, 
	   NULL, &nrhs ,
	   iparm, &msglvl, 
	   &rhs[0][0], sol,
	   &error
#ifndef WITH_MKL_PARDISO
	   , dparm
#endif
  );

  if ( error != 0) {
    printf ("PardisoSolver - ERROR during solve: %d", error ); 
  }
  memcpy(&rhs[0][0], sol, sz_sol*sizeof(double));
}

void PardisoSolver::solve( GenMatrix& rhs_in, int *colSparsity)
{
  //cout << "PARDISO-multiple sparse rhs" << endl;    
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

  int nrows,ncols; rhs.getSize(ncols,nrows);
  if(sz_sol<nrows*ncols) {
    sz_sol=nrows*ncols;
    if(sol) delete[] sol;
    sol = new double[sz_sol];
  }
  assert(nrows==n);


  /* same for mkl_pardiso and pardiso */
   int phase = 33; // solve and iterative refinement
   int mtype = -2, error;
   int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
   int mnum = 1; // actual matrix (as in index from 1 to maxfct)
   int nrhs = ncols;
   int msglvl = 0;
   iparm[1] = 2; // 2 is for metis, 0 for min degree
   //iparm[5] = 1; // replace drhs with the solution
   iparm[7] = 2; // max number of iterative refinement steps
   iparm[30]= 1; // sparse rhs

 #ifndef WITH_MKL_PARDISO
   iparm[2] = num_threads;
 #else
   /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
 #endif

  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, M, krowM, jcolM, 
	   colSparsity, &nrhs,
	   iparm, &msglvl, 
	   &rhs[0][0], sol,
	   &error
#ifndef WITH_MKL_PARDISO
	   , dparm
#endif
  );
 
  if ( error != 0) {
    printf ("PardisoSolver - ERROR during solve: %d", error ); 
  }
  memcpy(&rhs[0][0], sol, sz_sol*sizeof(double));
}

void PardisoSolver::solve( int nrhss, double* rhss, int* colSparsity )
{
   assert(rhss);
   assert(nrhss >= 1);

   if( sz_sol < nrhss * n )
   {
      sz_sol = nrhss * n;
      delete[] sol;

      sol = new double[sz_sol];
   }

#ifndef NDEBUG
   if( colSparsity )
   {
      for( int nr = 0; nr < nrhss; nr++ )
      {
         for( int i = 0; i < n; i++ )
         {
            const int rhspos = nr * n + i;
            if( rhss[rhspos] != 0.0 )
               assert(colSparsity[i] == 1);
            else if( nrhss == 1 ) // does not work with zeroes in matrix, e.g. callback example
               assert(colSparsity[i] == 0);

         }
      }
   }
#endif

   /* same for mkl_pardiso and pardiso */
   int phase = 33; // solve and iterative refinement
   int mtype = -2, error;
   int maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
   int mnum = 1; // actual matrix (as in index from 1 to maxfct)
   int nrhss_local = nrhss;
   int msglvl = 0;
   iparm[1] = 2; // 2 is for metis, 0 for min degree
   //iparm[5] = 1; // replace drhs with the solution
   iparm[7] = 5; // max number of iterative refinement steps
   if( colSparsity )
      iparm[30] = 1; //sparse rhs
   else
      iparm[30] = 0;

 #ifndef WITH_MKL_PARDISO
   iparm[2] = num_threads;
   //iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
 #else
   /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
   //iparm[24] = 2; // parallelization for the forward and backward solve. 0=sequential, 1=parallel in core solve, 2=parallel solve.
 #endif


   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, colSparsity,
         &nrhss_local, iparm, &msglvl, rhss, sol, &error
#ifndef WITH_MKL_PARDISO
         , dparm
#endif
   );

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }

   memcpy(rhss, sol, n * nrhss * sizeof(double));
}


PardisoSolver::~PardisoSolver()
{

  //cout << "PardisoSolver DESTRUCTOR" << endl;
  int phase = -1; // release internal memory
  int maxfct=1; // max number of fact having same sparsity pattern to keep at the same time
  int mnum=1; // actual matrix (as in index from 1 to maxfct)
  int nrhs=1;
  int msglvl=0;
  int mtype=-2, error;

  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, NULL, krowM, jcolM, NULL, &nrhs,
	   iparm, &msglvl, NULL, NULL, &error
#ifndef WITH_MKL_PARDISO
	   , dparm
#endif
  );


  if ( error != 0) {
    printf ("PardisoSolver - ERROR in pardiso release: %d", error ); 
  }
  delete[] jcolM;
  delete[] krowM;
  delete[] M;
  delete[] nvec;
  if(sol) delete[] sol;
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



