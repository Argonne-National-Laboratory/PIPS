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
#include <cassert>
#include "mpi.h"
#include "omp.h"

#define CHECK_PARDISO


/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
extern "C" void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
extern "C" void pardiso_chkvec     (int *, int *, double *, int *);
extern "C" void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

PardisoIndefSolver::PardisoIndefSolver( DenseSymMatrix * dm )
{
  mStorage = DenseStorageHandle( dm->getStorage() );
  mStorageSparse = NULL;

  assert(mStorage);

  n = mStorage->n;

  initPardiso();
}


PardisoIndefSolver::PardisoIndefSolver( SparseSymMatrix * sm )
{
  mStorage = NULL;
  mStorageSparse = SparseStorageHandle( sm->getStorage() );

  assert(mStorageSparse);

  n = mStorageSparse->n;

  initPardiso();
}

void PardisoIndefSolver::initPardiso()
{
   int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   mtype = -2;
   nrhs = 1;

   int error = 0;
   solver = 0; /* use sparse direct solver */
   pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

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

   iparm[2] = PIPSgetnOMPthreads();

   maxfct = 1; /* Maximum number of numerical factorizations.  */
   mnum = 1; /* Which factorization to use. */

   msglvl = 0; /* Print statistical information?  */
   ia = NULL;
   ja = NULL;
   a = NULL;
   ddum = -1.0;
   idum = -1;
   phase = 11;
   x = NULL;
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

#define SELECT_NNZS

void PardisoIndefSolver::factorizeFromSparse()
{
   assert(mStorageSparse);

   const int nnz = mStorageSparse->len;
   const int* const iaStorage = mStorageSparse->krowM;
   const int* const jaStorage = mStorageSparse->jcolM;
   const double* const aStorage = mStorageSparse->M;

   // first call?
   if( ia == NULL )
   {
      assert(ja == NULL && a == NULL);

      ia = new int[n + 1];
      ja = new int[nnz];
      a = new double[nnz];

#ifndef SELECT_NNZS

      for( int i = 0; i < n + 1; i++ )
         ia[i] = iaStorage[i] + 1;

      for( int i = 0; i < nnz; i++ )
         ja[i] = jaStorage[i] + 1;
#endif
   }

   assert(n >= 0);

#ifdef SELECT_NNZS
   std::vector<double>diag(n);

   const double t = 0.0001;

   for( int r = 0; r < n; r++ )
   {
      const int j = iaStorage[r];
      assert(jaStorage[j] == r);

      diag[r] = fabs(aStorage[j]) * t;
   }

   ia[0] = 1;
   int kills = 0;
   int nnznew = 0;

   for( int r = 0; r < n; r++ )
   {
      for( int j = iaStorage[r]; j < iaStorage[r + 1]; j++ )
      {
         //if( fabs(aStorage[j]) > 1e-15 || jaStorage[j] == r )
         if( aStorage[j] != 0.0 || jaStorage[j] == r )
         {
#ifdef SPARSE_PRECOND
            if( (fabs(aStorage[j]) >= diag[r] || fabs(aStorage[j]) >= diag[jaStorage[j]]) )
            {
               ja[nnznew] = jaStorage[j] + 1;
               a[nnznew++] = aStorage[j];
            }
            else
            {
               kills++;
               assert(jaStorage[j] != r);
            }
#else
            ja[nnznew] = jaStorage[j] + 1;
            a[nnznew++] = aStorage[j];
#endif

         }
      }
      ia[r + 1] = nnznew + 1;
   }

   std::cout << "real nnz in KKT: " << nnznew << " (kills: " << kills << ")" << std::endl;

#else
   for( int i = 0; i < nnz; i++ )
      a[i] = aStorage[i];
   for( int i = 0; i < n + 1; i++ )
      assert(ia[i] == iaStorage[i] + 1);

   for( int i = 0; i < nnz; i++ )
      assert(ja[i] == jaStorage[i] + 1);
#endif

   // matrix initialized, now do the actual factorization
   factorize();
}


void PardisoIndefSolver::factorizeFromDense()
{
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   assert(mStorage);

#ifdef DENSE_USE_HALF
#ifndef NDEBUG
  for( int i = 0; i < n; i++ )
     for( int j = 0; j < n; j++ )
        assert(j <= i || mStorage->M[i][j] == 0.0);
#endif
#ifdef TIMING
  if( myrank == 0 )
     std::cout << "DENSE_USE_HALF: starting factorization" << std::endl;
#endif
#endif


   assert(mStorage->n == mStorage->m);
   int nnz = 0;
   for( int i = 0; i < n; i++ )
      for( int j = 0; j <= i; j++ )
         if( mStorage->M[i][j] != 0.0 )
            nnz++;

   if( ia )
      delete[] ia;

   if( ja )
      delete[] ja;

   if( a )
      delete[] a;

   ia = new int[n + 1];
   ja = new int[nnz];
   a = new double[nnz];

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

#ifndef NDEBUG
#ifdef CHECK_PARDISO
   pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
   if( error != 0 )
   {
      printf("\nERROR in consistency of matrix: %d", error);
      exit(1);
   }
#endif
#endif

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

   iparm[9] = 13; // pivot perturbation 10^{-xxx}
#ifdef PARDISOINDEF_SCALE
   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
   if( myrank == 0)
      std::cout << "... SCALE PARDISO for global SC" << std::endl;
#endif

#ifdef PARDISO_PARALLEL_AGGRESSIVE
   //iparm[1] = 3; // 3 Metis 5.1 (only for PARDISO >= 6.0)
   iparm[23] = 1; //Parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
   iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   //iparm[27] = 1; // Parallel metis
#endif

   phase = 11;

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error, dparm);

   if( error != 0 )
   {
      printf("\nERROR during symbolic factorization: %d", error);
      exit(1);
   }

   if( myrank == 0 )
   {
      printf("\nReordering completed: ");
      printf("\nNumber of nonzeros in factors  = %d", iparm[17]);
      printf("          Number of factorization MFLOPS = %d", iparm[18]);
   }

   phase = 22;
   // iparm[32] = 1; /* compute determinant */

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error, dparm);

   if( error != 0 )
   {
      printf("\nERROR during numerical factorization: %d", error);
      exit(2);
   }
}

void PardisoIndefSolver::solve ( OoqpVector& v )
{
#ifdef PARDISO_PARALLEL_AGGRESSIVE
   assert(iparm[23] == 1);
   assert(iparm[24] == 1);
#endif

   int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   phase = 33;

   iparm[7] = 8; /* max number of iterative refinement steps. */

   SimpleVector& sv = dynamic_cast<SimpleVector&>(v);

   double* b = sv.elements();

   assert(sv.n == n);

#ifdef TIMING_FLOPS
   HPM_Start("DSYTRSSolve");
#endif
   if( myrank == 0 )
   {
      int error;

      // first call?
      if( !x )
         x = new double[n];

      pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, b, x, &error, dparm);

      if( error != 0 )
      {
         printf("\nERROR during solution: %d", error);
         exit(3);
      }
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
         iparm, &msglvl, &ddum, &ddum, &error, dparm);

   delete[] ia;
   delete[] ja;
   delete[] a;
   delete[] x;
}
