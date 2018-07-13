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

  assert(0);
  // wontworkbeause_NULL;
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

   int num_procs;

   /* Numbers of processors, value of OMP_NUM_THREADS */
   char* var = getenv("OMP_NUM_THREADS");
   if( var != NULL )
      sscanf(var, "%d", &num_procs);
   else
   {
      printf("Set environment OMP_NUM_THREADS to 1");
      exit(1);
   }
   iparm[2] = num_procs;

   maxfct = 1; /* Maximum number of numerical factorizations.  */
   mnum = 1; /* Which factorization to use. */

   msglvl = 0; /* Print statistical information  */
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
   // todo: only on rank 0

   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   if( myrank == 0 )
      printf("\nFactorization starts ...\n ");

   if( mStorageSparse )
      factorizeFromSparse();
   else
      factorizeFromDense();

   if( myrank == 0 )
      printf("\nFactorization completed ...\n ");
}


void PardisoIndefSolver::factorizeFromSparse()
{
   assert(mStorageSparse);

   const int nnz = mStorageSparse->len;

   const double* const aStorage = mStorageSparse->M;

   // first call?
   if( ia == NULL )
   {
      assert(ja == NULL && a == NULL);

      const int* const iaStorage = mStorageSparse->krowM;
      const int* const jaStorage = mStorageSparse->jcolM;

      ia = new int[n + 1];
      ja = new int[nnz];
      a = new double[nnz];

      for( int i = 0; i < n + 1; i++ )
         ia[i] = iaStorage[i] + 1;

      for( int i = 0; i < nnz; i++ )
         ja[i] = jaStorage[i] + 1;
   }

   assert(n >= 0);

   for( int i = 0; i < nnz; i++ )
      a[i] = aStorage[i];

   // matrix initialized, now do the actual factorization
   factorize();
}


void PardisoIndefSolver::factorizeFromDense()
{
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   assert(mStorage);

#if 0
   ofstream myfile;
   myfile.open("xout");

   int all = 0;
   int z = 0;
   for( int i = 0; i < n; i++ )
   {
      for( int j = i + 1; j < n; j++ )
         if( mStorage->M[i][j] != 0.0 )
         {
            std::cout << "FAIL " << std::endl;
            exit(1);
         }
      for( int j = 0; j <= i; j++ )
      {
         all++;
         if( mStorage->M[i][j] != 0.0 )
         {
            z++;

         }
         myfile << mStorage->M[i][j] << ", ";
      }

      for( int j = i + 1; j < n; j++ )
         myfile << "0.0, ";

      myfile << "\n";
   }
   cout << "full " << n * n / 2 << " \n";
   cout << "  all " << all << "\n";
   cout << "non-zeros " << z << "\n";

   myfile.close();
#endif


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

/*
   for( int i = 0; i < n; i++ )
   {
      for( int j = 0; j <= i; j++ )
      {
         cout << mStorage->M[i][j] << "      ";
      }

      for( int j = i + 1; j < n; j++ )
         cout << "0.0, ";

      cout << "\n";
   }

*/

   if( ia ) // todo store and later resize
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

#ifdef CHECK_PARDISO
   pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
   if( error != 0 )
   {
      printf("\nERROR in consistency of matrix: %d", error);
      exit(1);
   }
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
   iparm[32] = 1; /* compute determinant */

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
   int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
   int myrank; MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   int n = mStorage->n;
   phase = 33;

   iparm[7] = 5; /* Max numbers of iterative refinement steps. */

   SimpleVector& sv = dynamic_cast<SimpleVector&>(v);

   double* b = sv.elements();

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

      if( size > 0 )
         MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      for( int i = 0; i < n; i++ )
         b[i] = x[i];

//      printf("number of iterative refinement steps: %d \n", iparm[7]);
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
   int n = mStorage->n;
   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs,
         iparm, &msglvl, &ddum, &ddum, &error, dparm);

   delete[] ia;
   delete[] ja;
   delete[] a;
   delete[] x;
}
