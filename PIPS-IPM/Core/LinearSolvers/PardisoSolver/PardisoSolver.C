/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin, Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

#include "PardisoSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include "pipsdef.h"
#include "pipsport.h"
#include <cstdlib>

#include "mpi.h"

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

extern int gOoqpPrintLevel;

#ifdef WITH_MKL_PARDISO
#include "mkl_pardiso.h"
#include "mkl_types.h"
#else
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
   const int myRank = PIPS_MPIgetRank();
   if( myRank == 0 )
	  std::cout << "PardisoSolver::PardisoSolver (sparse input)" << std::endl;
#endif
  Msys = sgm;
  Mdsys = nullptr;
  n = sgm->size();
  nnz = sgm->numberOfNonZeros();

  krowM = new int[n+1];
  jcolM = new int[nnz];
  M = new double[nnz];

  sol = nullptr;
  sz_sol = 0;

  first = true;
  nvec=new double[n];

#ifndef WITH_MKL_PARDISO
  num_threads = PIPSgetnOMPthreads();
  solver = 0; /* sparse direct solver */
#endif

  mtype = -2;
  error = 0;
  phase = 0;
  nrhs = -1; // do not specify here
  maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
  mnum = 1; // actual matrix (as in index from 1 to maxfct)
  msglvl = 0; // messaging level (0 = no output, 1 = statistical info to screen)
}


PardisoSolver::PardisoSolver( DenseSymMatrix * m )
{
#ifdef TIMING
   const int myRank = PIPS_MPIgetRank();
   if( myRank == 0 )
	  std::cout << "PardisoSolver created (dense input)" << std::endl;
#endif
  Msys = nullptr;
  Mdsys = m;
  n = m->size();
  
  nnz=0; //getNumberOfNonZeros(*Mdsys);

  krowM = nullptr;//new int[n+1];
  jcolM = nullptr;//new int[nnz];
  M = nullptr;//new double[nnz];

  sol = nullptr;
  sz_sol = 0;

  first = true;
  nvec=new double[n];

#ifndef WITH_MKL_PARDISO
  num_threads = PIPSgetnOMPthreads();
  solver = 0; /* sparse direct solver */
#endif

  mtype = -2;
  error = 0;
  nrhs = -1;
  phase = 0;
  maxfct = 1; // max number of fact having same sparsity pattern to keep at the same time
  mnum = 1; // actual matrix (as in index from 1 to maxfct)
  msglvl = 0; // messaging level (0 = no output, 1 = statistical info to screen)
}

void PardisoSolver::firstCall()
{
   iparm[0] = 0;

#ifndef WITH_MKL_PARDISO
   int error = 0;

   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

   if( error != 0 )
   {
      std::cout << "PardisoSolver ERROR during pardisoinit:" << error << "." << std::endl;
      assert(false);
   }
#else
   pardisoinit(pt, &mtype, iparm);
#endif

   setIparm(iparm);

   if( Msys )
   {
      //get the matrix in upper triangular
      Msys->getStorageRef().transpose(krowM, jcolM, M);

      //save the indices for diagonal entries for a streamlined later update
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
         assert(idxDiagMsys >= 0);
         //must have all diagonals entries

         // aug - find the index in jcol for the diagonal (r,r)
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
   else if( Mdsys )
   {
      // the input is a dense matrix
      // the dense matrix is also processed in matrixChanged everytime the method is called

      // the input is a dense matrix
      DenseSymMatrix& Md = (*Mdsys);
      nnz = Md.getNumberOfNonZeros();

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

   // need Fortran indexes
   // extra method?
   for( int i = 0; i < n + 1; i++ )
      krowM[i] += 1;
   for( int i = 0; i < nnz; i++ )
      jcolM[i] += 1;

} 


/*
 * iparm[30] has to be set depending on the circumstances!
 */
void PardisoSolver::setIparm(int* iparm){

   iparm[1] = 2; // 2 is for metis, 0 for min degree

#ifndef WITH_MKL_PARDISO
   iparm[2] = num_threads;
   iparm[7] = 2; // max number of iterative refinements
   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
                  // if needed, use 2 for advanced matchings and higher accuracy.

#else
   /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
   iparm[7] = 0; // mkl runs into numerical problems when setting iparm[7] too high
   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
   iparm[12] = 1; // 0 disable matching, 1 enable matching, no other settings

   #ifndef NDEBUG
   // enable matrix checker - mkl pardiso has no chkmatrix method
   iparm[26] = 1;
   #endif
#endif

}

bool PardisoSolver::iparmUnchanged()
{

   bool same = true;
   bool print = true;

   int iparm_compare[64];
   setIparm(iparm_compare);

   static const int arr[] = { 1, 7, 10, 12 };
   vector<int> to_compare(arr, arr + sizeof(arr) / sizeof(arr[0]) );

   for(int i = 0; i < 64; ++i)
   {
      // if entry should be compared
      if(std::find(to_compare.begin(), to_compare.end(), i) != to_compare.end())
      {
         if(iparm[i] != iparm_compare[i])
         {
            if(print)
               std::cout << "ERROR - PardisoSolver: elements in iparm changed at " << i << ": "
                  << iparm[i] << " != " << iparm_compare[i] << "(new)" << std::endl;
            same = false;
         }
      }

   }
   return same;
}
 
void PardisoSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

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


   phase = 12; // Analysis, numerical factorization
   nrhs = 1;

   iparm[30] = 0; // do not specify sparse rhs
                  // if one wants to use the sparse rhs and partial solves according to
                  // the PARDISO user guide the perm vector has to be present during all
                  // phases of pardiso

   assert(iparmUnchanged());

     pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr, &nrhs,
         iparm, &msglvl, nullptr, nullptr, &error
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
   phase = 33; // solve and iterative refinement
   nrhs = 1;
   iparm[30] = 0; // do not specify sparse rhs at this point

   assert(iparmUnchanged());

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM,
         nullptr, &nrhs, iparm, &msglvl, rhs.elements(), sol_local, &error
#ifndef WITH_MKL_PARDISO
         ,dparm
#endif
   );

   if( error != 0 )
   {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }

   rhs.copyFromArray(sol_local);
}

void PardisoSolver::solve( GenMatrix& rhs_in )
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

  int nrows,ncols; rhs.getSize(ncols,nrows);
  if(sz_sol<nrows*ncols) {
    sz_sol=nrows*ncols;
    if(sol) delete[] sol;
    sol = new double[sz_sol];
  }

  assert(nrows==n);


  phase = 33; // solve and iterative refinement
  nrhs = ncols;
  iparm[30] = 0;
  assert(iparmUnchanged());

  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, M, krowM, jcolM, 
	   nullptr, &nrhs ,
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
  std::cout << "PardisoSolver - using sparse rhs but might lead to numerical troubles .. " << std::endl;
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

  int nrows,ncols; rhs.getSize(ncols,nrows);
  if(sz_sol<nrows*ncols) {
    sz_sol=nrows*ncols;
    if(sol) delete[] sol;
    sol = new double[sz_sol];
  }
  assert(nrows==n);

  /* same for mkl_pardiso and pardiso */
   phase = 33; // solve and iterative refinement
   iparm[30] = 1;

   nrhs = ncols;

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
   phase = 33; // solve and iterative refinement
   int nrhss_local = nrhss;

   assert(iparmUnchanged());

// see notes on [30] earlier - cannot be specified on the go
// seems to detoriate parformance a lot - triggers many BiCGStab steps
//   if( colSparsity )
//   {
//      iparm[30] = 1; //sparse rhs
//   }
//   else
   {
      iparm[30] = 0;
   }

   assert(pt); assert(M); assert(krowM); assert(jcolM); assert(rhss); assert(sol);

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr,
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

  phase = -1; // release internal memory
  nrhs = 1;

  pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	   &n, nullptr, krowM, jcolM, nullptr, &nrhs,
	   iparm, &msglvl, nullptr, nullptr, &error
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
