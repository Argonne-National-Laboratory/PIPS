/*
 * PardisoIndefSolver.C
 *
 *  Created on: 21.03.2018
 *      Author: Daniel Rehfeldt
 */


#include "PardisoIndefSolver.h"


#include "pipschecks.h"
#include "SimpleVector.h"
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <cassert>
#include "mpi.h"
#include "omp.h"
#include "pipsport.h"
#include "pipsdef.h"
#include "StochOptions.h"


#define CHECK_PARDISO

#ifdef WITH_MKL_PARDISO
#include "mkl_pardiso.h"
#include "mkl_types.h"
#endif

/* PARDISO prototype. */
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

PardisoIndefSolver::PardisoIndefSolver( DenseSymMatrix * dm )
{
  mStorage = dm->getStorageHandle();
  mStorageSparse = nullptr;

  assert(mStorage);

  n = mStorage->n;

  initPardiso();
}


PardisoIndefSolver::PardisoIndefSolver( SparseSymMatrix * sm )
{
  mStorage = nullptr;
  mStorageSparse = sm->getStorageHandle();

  assert(mStorageSparse);

  n = mStorageSparse->n;

  initPardiso();
}

void PardisoIndefSolver::setIparm(int* iparm){

   iparm[9] = 13; /* pivot perturbation 10^{-xxx} */

#ifndef WITH_MKL_PARDISO
   /* From INTEL (instead of iparm[2] which is not defined there):
    *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
    *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
    *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
    */
   assert(nThreads >= 1 && pivotPerturbationExp >= 1);
   iparm[9] = pivotPerturbationExp;
   iparm[2] = nThreads;
   iparm[18] = 0; /* don't compute GFLOPS */
   iparm[7] = nIterativeRefins; /* max number of iterative refinement steps. */

   if( highAccuracy )
   {
      iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
      iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
   }
   else
   {
      iparm[10] = 1;
      iparm[12] = 0;
   }

   if( factorizationTwoLevel )
      iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
   else
      iparm[23] = 0;

   if( parallelForwardBackward )
      iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   else
      iparm[24] = 0;

#else

   iparm[7] = 0; /* mkl runs into numerical problems when setting iparm[7] too high */

   /* enable matrix checker (default disabled) - mkl pardiso does not have chkmatrix */
   #ifndef NDEBUG
   iparm[26] = 1;
   #endif

   if( highAccuracy )
   {
      iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
      iparm[12] = 1;// MKL does not have a =2 option here
   }
   else
   {
      iparm[10] = 0;
      iparm[12] = 0;
   }
#endif

}

bool PardisoIndefSolver::iparmUnchanged()
{
   /* put all Parameters that should stay be checked against init into this array */
   static const int check_iparm[] =
      { 7, 10, 12, 23, 24 };

   bool unchanged = true;
   bool print = false;

   int iparm_compare[64];
   setIparm(iparm_compare);

   vector<int> to_compare(check_iparm,
         check_iparm + sizeof(check_iparm) / sizeof(check_iparm[0]));

   for( int i = 0; i < 64; ++i )
   {
      // if entry should be compared
      if( std::find(to_compare.begin(), to_compare.end(), i)
            != to_compare.end() )
      {
         if( iparm[i] != iparm_compare[i] )
         {
            if( print )
               std::cout
                     << "ERROR - PardisoSolver: elements in iparm changed at "
                     << i << ": " << iparm[i] << " != " << iparm_compare[i]
                     << "(new)" << std::endl;
            unchanged = false;
         }
      }
   }
   return unchanged;
}

void PardisoIndefSolver::initPardiso()
{
   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   deleteCSRpointers = false;
   mtype = -2;
   nrhs = 1;
   iparm[0] = 0;

   useSparseRhs = useSparseRhsDefault;

   nThreads = PIPSgetnOMPthreads();
   char* var = getenv("OMP_NUM_THREADS_PIPS_ROOT");
   if( var != nullptr )
   {
      int n = -1;
      sscanf(var, "%d", &n);

      assert(n >= 1);

      nThreads = n;
   }

   pivotPerturbationExp = pips_options::getIntParameter("PARDISO_PIVOT_PERTURBATION_ROOT");
   if( pivotPerturbationExp < 0 )
	   pivotPerturbationExp = pivotPerturbationExpDefault;

   highAccuracy = highAccuracyDefault;
   parallelForwardBackward = parallelForwardBackwardDefault,
   factorizationTwoLevel = factorizationTwoLevelDefault,

   nIterativeRefins = pips_options::getIntParameter("PARDISO_NITERATIVE_REFINS_ROOT");
   if( nIterativeRefins < 0 )
 	  nIterativeRefins = nIterativeRefinsDefault;

   if( myRank == 0 )
   {
      printf("PARDISO root: using pivot perturbation 10^-%d \n", pivotPerturbationExp);

      printf("PARDISO root: using maximum of %d iterative refinements  \n", nIterativeRefins);

      if( parallelForwardBackward )
         printf("PARDISO root: using parallel (forward/backward) solve \n");
      else
         printf("PARDISO root: NOT using parallel (forward/backward) solve \n");

      if( factorizationTwoLevel )
         printf("PARDISO root: using two-level scheduling for numerical factorization \n");
      else
         printf("PARDISO root: NOT using two-level scheduling for numerical factorization \n");

      if( highAccuracy )
         printf("PARDISO root: using high accuracy \n");
      else
         printf("PARDISO root: NOT using high accuracy \n");

      if( useSparseRhs )
         printf("PARDISO root: using sparse rhs \n");
      else
         printf("PARDISO root: NOT using sparse rhs \n");
   }

#ifndef WITH_MKL_PARDISO
   int error = 0;
   solver = 0; /* use sparse direct solver */

   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
   if( error != 0 )
   {
      if( error == -10 )
         printf("No license file found \n");
      if( error == -11 )
         printf("License is expired \n");
      if( error == -12 )
         printf("Wrong username or hostname \n");
      exit(1);
   }
   else if( myRank == 0 )
      printf("[PARDISO]: License check was successful ... \n");

#else
   pardisoinit(pt, &mtype, iparm);
#endif

   setIparm(iparm);

   if( myRank == 0 )
      printf("PARDISO root: using %d threads \n", iparm[2]);

   maxfct = 1; /* Maximum number of numerical factorizations.  */
   mnum = 1; /* Which factorization to use. */

   msglvl = 0; /* Print statistical information?  */
   ia = nullptr;
   ja = nullptr;
   a = nullptr;
   ddum = -1.0;
   idum = -1;
   phase = 11;
   x = nullptr;
}


void PardisoIndefSolver::matrixChanged()
{
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if( myrank == 0 )
   {
      printf("\n Schur complement factorization is starting ...\n ");

      if( mStorageSparse )
         factorizeFromSparse();
      else
         factorizeFromDense();

      printf("\n Schur complement factorization completed \n ");
   }
}


void PardisoIndefSolver::matrixRebuild( DoubleMatrix& matrixNew )
{
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if( myrank == 0 )
   {
      SparseSymMatrix& matrixNewSym = dynamic_cast<SparseSymMatrix&>(matrixNew);

      assert(matrixNewSym.getStorageRef().fortranIndexed());

      printf("\n Schur complement factorization is starting ...\n ");

      assert(mStorageSparse);

      factorizeFromSparse(matrixNewSym);

      printf("\n Schur complement factorization completed \n ");
   }
}


void PardisoIndefSolver::factorizeFromSparse(SparseSymMatrix& matrix_fortran)
{
   assert(n == matrix_fortran.size());
   assert(!deleteCSRpointers);
   assert(matrix_fortran.getStorageRef().fortranIndexed());

   ia = matrix_fortran.krowM();
   ja = matrix_fortran.jcolM();
   a = matrix_fortran.M();

   assert(ia[0] == 1);

   // matrix initialized, now do the actual factorization
   factorize();
}


void PardisoIndefSolver::factorizeFromSparse()
{
   assert(mStorageSparse);

   const int nnz = mStorageSparse->len;
   const int* const iaStorage = mStorageSparse->krowM;
   const int* const jaStorage = mStorageSparse->jcolM;
   const double* const aStorage = mStorageSparse->M;
   const bool usePrecondSparse = pips_options::getBoolParameter("PRECONDITION_SPARSE");

   // first call?
   if( ia == nullptr )
   {
      assert(ja == nullptr && a == nullptr);
      deleteCSRpointers = true;

      ia = new int[n + 1];
      ja = new int[nnz];
      a = new double[nnz];
   }

   assert(n >= 0);

   std::vector<double>diag(n);

   // todo the sparse precond. stuff should be moved out and handled by Sparsifier class
   if( usePrecondSparse )
   {
	  const double t = precondDiagDomBound;

	  for( int r = 0; r < n; r++ )
	  {
		  const int j = iaStorage[r];
		  assert(jaStorage[j] == r);

		  diag[r] = fabs(aStorage[j]) * t;
	  }
   }

   ia[0] = 1;
   int nnznew = 0;

   for( int r = 0; r < n; r++ )
   {
      for( int j = iaStorage[r]; j < iaStorage[r + 1]; j++ )
      {
         if( aStorage[j] != 0.0 || jaStorage[j] == r )
         {
        	if( usePrecondSparse )
        	{
				if( (fabs(aStorage[j]) >= diag[r] || fabs(aStorage[j]) >= diag[jaStorage[j]]) )
				{
				   ja[nnznew] = jaStorage[j] + 1;
				   a[nnznew++] = aStorage[j];
				}
				else
				{
				   assert(jaStorage[j] != r);
				}
        	}
        	else
        	{
               ja[nnznew] = jaStorage[j] + 1;
               a[nnznew++] = aStorage[j];
        	}
         }
      }
      ia[r + 1] = nnznew + 1;
   }

   std::cout << "real nnz in KKT: " << nnznew << " (ratio: " << double(nnznew) / double(iaStorage[n]) << ")" << std::endl;

#if 0
   {
      ofstream myfile;
      int mype;  MPI_Comm_rank(MPI_COMM_WORLD, &mype);

      printf("\n\n ...WRITE OUT! \n\n");

      if( mype == 0 )
      {
         myfile.open("../A.txt");

         for( int i = 0; i < n; i++ )
            for( int k = ia[i]; k < ia[i + 1]; k++ )
               myfile << i << '\t' << ja[k - 1] << '\t' << a[k - 1] << endl;

         myfile.close();
      }

      printf("%d...exiting (pardiso) \n", mype);
      exit(1);
  }
#endif

   // matrix initialized, now do the actual factorization
   factorize();
}


void PardisoIndefSolver::factorizeFromDense()
{
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   assert(mStorage);

#ifndef NDEBUG
  for( int i = 0; i < n; i++ )
  {
     for( int j = 0; j < n; j++ )
        assert(j <= i || PIPSisZero(mStorage->M[i][j]));
  }
#endif
#ifdef TIMING
  if( myrank == 0 )
     std::cout << "from dense, starting factorization" << std::endl;
#endif

   assert(mStorage->n == mStorage->m);
   int nnz = 0;
   for( int i = 0; i < n; i++ )
      for( int j = 0; j <= i; j++ )
         if( mStorage->M[i][j] != 0.0 )
            nnz++;

   if( deleteCSRpointers )
   {
      delete[] ia;
      delete[] ja;
      delete[] a;
   }

   ia = new int[n + 1];
   ja = new int[nnz];
   a = new double[nnz];

   deleteCSRpointers = true;

   nnz = 0;
   for( int j = 0; j < n; j++ )
   {
      ia[j] = nnz;
      for( int i = j; i < n; i++ )
         if( mStorage->M[i][j] != 0.0 )
         {
            ja[nnz] = i;
            a[nnz++] = mStorage->M[i][j];
         }
   }

   ia[n] = nnz;

   for( int i = 0; i < n + 1; i++ )
      ia[i] += 1;
   for( int i = 0; i < nnz; i++ )
      ja[i] += 1;

   factorize();
}


void PardisoIndefSolver::factorize()
{
   int error;
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   assert(ia && ja && a);

#if !defined(NDEBUG) && defined(CHECK_PARDISO) && !defined(WITH_MKL_PARDISO)
   pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
   if( error != 0 )
   {
      printf("\nERROR in consistency of matrix: %d", error);
      exit(1);
   }
#endif

   iparm[17] = -1; /* compute number of nonzeros in factors */

#if 0
   const int nnz = ia[n] - 1;
   double abs_max = 0.0;
   for( int i = 0; i < nnz; i++ )
   {
      const double abs = std::fabs(a[i]);
      if( abs > abs_max)
         abs_max = abs;
   }

   std::cout << "absmax=" << abs_max << " log=" << log10(abs_max) << std::endl;
   if( log10(abs_max) >= 13)
   {
      iparm[9] = min(int(log10(abs_max)), 15);
      std::cout << "new: param " << iparm[9] << std::endl;
   }
else
#endif

   phase = 11;
   assert(iparmUnchanged());

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error
#ifndef WITH_MKL_PARDISO
         , dparm
#endif
   );

   if( error != 0 )
   {
      printf("\nERROR during symbolic factorization: %d", error);
      exit(1);
   }

   if( myrank == 0 )
   {
      printf("\nReordering completed: ");
      printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
   }

   phase = 22;
   assert(iparmUnchanged());

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error
#ifndef WITH_MKL_PARDISO
         , dparm
#endif
   );

   if( error != 0 )
   {
      printf("\nERROR during numerical factorization: %d", error);
      exit(2);
   }
}

void PardisoIndefSolver::solve ( OoqpVector& v )
{
   assert(iparmUnchanged());

   int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   phase = 33;
   SimpleVector& sv = dynamic_cast<SimpleVector&>(v);

   double* b = sv.elements();

   assert(sv.n == n);

#ifdef TIMING_FLOPS
   HPM_Start("DSYTRSSolve");
#endif
   if( myrank == 0 )
   {
      int* rhsSparsity = nullptr;

      int error;

      // first call?
      if( !x )
         x = new double[n];

      if( useSparseRhs )
      {
         iparm[30] = 1; //sparse rhs
         rhsSparsity = new int[n]();

         for( int i = 0; i < n; i++  )
            if( !PIPSisZero(b[i]) )
               rhsSparsity[i] = 1;
      }
      else
      {
         iparm[30] = 0;
      }

      pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, rhsSparsity, &nrhs,
            iparm, &msglvl, b, x, &error
#ifndef WITH_MKL_PARDISO
            ,dparm
#endif
      );
      if( error != 0 )
      {
         printf("\nERROR during solution: %d", error);
         exit(3);
      }

      iparm[30] = 0;
#if 0
      const double b2norm = sv.twonorm();
      const double binfnorm = sv.infnorm();
      double mat_max = 0.0;
      for( int i = 0; i < n; i++ )
      {
         for( int p = ia[i]; p < ia[i + 1]; p++ )
         {
            const int j = ja[p - 1] - 1;

            sv[i] -= a[p - 1] * x[j]; //r[i] = r[i] - M(i,j)*x(j)

            mat_max = std::max(std::fabs(a[p - 1]), mat_max);

            assert(j >= i);

            if( j != i )
               sv[j] -= a[p - 1] * x[i]; //r[j] = r[j] - M(j,i)*x(i)
         }
      }
      const double res2norm = sv.twonorm();
      const double resinfnorm = sv.infnorm();

      std::cout << "GLOBAL SCHUR: res.2norm=" << res2norm << " rel.res2norm=" << res2norm / b2norm  <<
            " res.infnorm=" << resinfnorm << " rel.resinfnorm=" << resinfnorm / binfnorm  << " b2norm=" << b2norm <<
            " abs.mat.elem.=" << mat_max << std::endl;
#endif

      for( int i = 0; i < n; i++ )
         b[i] = x[i];

      if( size > 0 )
         MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      delete[] rhsSparsity;

#ifdef TIMING
      printf("sparse kkt iterative refinement steps: %d \n", iparm[6]);
#endif
   }
   else
   {
      assert(size > 0);
      MPI_Bcast(b, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }
#ifdef TIMING_FLOPS
   HPM_Stop("DSYTRSSolve");
#endif
}

void PardisoIndefSolver::solve ( GenMatrix& rhs_in )
{
   assert(0 && "not supported");
}

void PardisoIndefSolver::diagonalChanged( int /* idiag */, int /* extent */ )
{
   this->matrixChanged();
}

PardisoIndefSolver::~PardisoIndefSolver()
{
   phase = -1; /* Release internal memory. */
   int error;
   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error
#ifndef WITH_MKL_PARDISO
         , dparm
#endif
   );

   if( deleteCSRpointers )
   {
      delete[] ia;
      delete[] ja;
      delete[] a;
   }

   delete[] x;
}
