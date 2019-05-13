/*
 * MumpsSolver.C
 */

#include <stdlib.h>

#include "MumpsSolver.h"
#include "SimpleVector.h"


#define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation
#define INFOG(I) infog[(I)-1]
#define RINFOG(I) rinfog[(I)-1]


MumpsSolver::MumpsSolver(long long n, MPI_Comm mpiCommPips, MPI_Comm mpiCommMumps)
 : n(n), mpiCommPips(mpiCommPips), mpiCommMumps(mpiCommMumps)
{
   maxNiterRefinments = defaultMaxNiterRefinments;
   verbosity = defaultVerbosity;

   rankMumps = -1;

   if( mpiCommMumps == MPI_COMM_NULL )
   {
      mumps = NULL;
   }
   else
   {
      setUpMumps();
      MPI_Comm_rank(mpiCommMumps, &rankMumps);
      assert(rankMumps >= 0);
   }

   MPI_Comm_rank(mpiCommPips, &rankPips);

   Msys = NULL;

}

MumpsSolver::~MumpsSolver()
{
   if( mumps != NULL )
   {
      mumps->job = -2;
      dmumps_c (mumps);

      delete mumps;
   }
}

void
MumpsSolver::diagonalChanged(int idiag, int extent)
{
   this->matrixChanged();
}

void
MumpsSolver::matrixChanged()
{

   if( mpiCommMumps == MPI_COMM_NULL )
      return;

   // todo: update only diagonal!

   /*
    if( Msys != NULL )
    sysMatToSpTriplet();
    */

   mumps->n = n;


   // symmetric permutation for factorization, 7: automatic choice; meaningless if mumps->ICNTL(28) == 2
   mumps->ICNTL(7)= 7;

   mumps->ICNTL(28)= 0; // choice of analysis, 0: automatic, 1: sequential, 2: parallel

   mumps->ICNTL(29)= 0; // parallel ordering, 0: automatic, 1: PT-SCOTCH, 2: ParMetis

   // relative threshold for numerical pivoting; 0.01 is default for sym. indef., larger values increase accuracy
   //  mumps->ICNTL(1) = 0.01;

   // relative threshold for static pivoting; -1.0: not used (default), 0.0: use with automatic choice of threshold
   //  mumps->ICNTL(4) = -1.0;


   // analysis phase
   mumps->job = 1;


   double starttime = MPI_Wtime();

   // do analysis
   dmumps_c(mumps);

   processMumpsResultAnalysis(starttime, true);


   // factorization phase
   mumps->job = 2;

   starttime = MPI_Wtime();

   // do factorization
   dmumps_c(mumps);

   processMumpsResultFactor(starttime, true);

   //saveOrderingPermutation();


}


void
MumpsSolver::solve(double* vec)
{
   mumps->rhs = vec;

   // solution phase
   mumps->job = 3;

   mumps->ICNTL(10) = maxNiterRefinments; //maximum number of iterative refinements
   mumps->ICNTL(11) = 2; // error statistics, 0: disabled, 2: main statistics

   const double starttime = MPI_Wtime();

   dmumps_c(mumps);

   processMumpsResultSolve(starttime, true);
}


void
MumpsSolver::solve(OoqpVector& rhs)
{
   SimpleVector& sv = dynamic_cast<SimpleVector &>(rhs);

   if( mpiCommMumps != MPI_COMM_NULL )
      solve(&sv[0]);
}

void
MumpsSolver::processMumpsResultAnalysis(double starttime, bool verbose)
{
   const int errorCode = mumps->INFOG(1);
   if( errorCode != 0 )
   {
      if( rankMumps == 0 )
         printf("Error INFOG(1)=%d in MUMPS analysis phase. \n", errorCode);

      exit(1);
   }

   if( verbose )
   {
      if( starttime >= 0.0 )
      {
         const double timeAnalysis = MPI_Wtime() - starttime;

         if( rankMumps == 0 )
            printf("MUMPS analysis phase took %g seconds.\n", timeAnalysis);
      }

      if( rankMumps != 0 )
         return;

      const int orderingType = mumps->INFOG(7);

      // sequential ordering used?
      if( mumps->INFOG(32) == 1 )
      {
         printf("MUMPS ordering done SEQUENTIALLY \n");
         printf("...Ordering: ");

         switch( orderingType )
         {
            case 0:
               printf("Approximate Minimum Degree (AMD) \n");
               break;
            case 2:
               printf("Approximate Minimum Fill \n");
               break;
            case 3:
               printf("ASCOTCH \n");
               break;
            case 4:
               printf("PORD \n");
               break;
            case 5:
               printf("Metis \n");
               break;
            case 6:
               printf("QAMD \n");
               break;
            default:
               assert(0);
         }
      }
      else
      {
         assert(mumps->INFOG(32) == 2);

         printf("MUMPS ordering done IN PARALLEL \n");
         printf("...Ordering: ");

         assert(orderingType == 1 || orderingType == 2);

         if( orderingType == 1 )
            printf("PT-Scotch \n");

         if( orderingType == 2 )
            printf("ParMetis \n");
      }

   }
}

void
MumpsSolver::processMumpsResultFactor(double starttime, bool verbose)
{
   const int errorCode = mumps->INFOG(1);
   if( errorCode != 0 )
   {
      if( rankMumps == 0 )
         printf("Error INFOG(1)=%d in MUMPS facorization phase. \n", errorCode);

      if( errorCode == -8 || errorCode == -9 )
      {
         // todo allocate more memory
#if 0
         const int try_max = 10;
         for( int tr = 0; tr < try_max; tr++ )
         {
            const double mem_percent = mumps_->icntl[13];
            mumps_->icntl[13] = 2 * mumps_->icntl[13];

            if( my_mumps_rank_ == 0 )
            printf("Increased Mumps ICNTL(14) from %g to %d\n", mem_percent,
                  mumps_->icntl[13]);

            dmumps_c(mumps_);
            error = mumps_->infog[1 - 1];

            saveOrderingPermutation();

            if( error != -8 && error != -9 )
            break;
         }

         if( error == -8 || error == -9 )
         {
            if( my_mumps_rank_ == 0 )
            printf(
                  "Fatal error in Mumps: not able to obtain more memory (ICNTL(14)-related)\n");
            exit(-1);
         }
#endif
      }

      exit(1);
   }

   if( verbose && starttime >= 0.0  )
   {
      const double timeFactor = MPI_Wtime() - starttime;

      if( rankMumps == 0 )
         printf("MUMPS factorization phase took %g seconds.\n", timeFactor);
   }
}

void
MumpsSolver::processMumpsResultSolve(double starttime, bool verbose)
{
   const int errorCode = mumps->INFOG(1);
   if( errorCode != 0 )
   {
      if( rankMumps == 0 )
         printf("Error INFOG(1)=%d in MUMPS solution phase. \n", errorCode);

      exit(1);
   }

   if( verbose  )
   {
      if( starttime >= 0.0 )
      {
         const double timeSolution = MPI_Wtime() - starttime;

         if( rankMumps == 0 )
            printf("MUMPS solution phase took %g seconds.\n", timeSolution);
      }

      assert(mumps->ICNTL(11) == 1 || mumps->ICNTL(11) == 2);

      const double infNormMatrix = mumps->RINFOG(4);
      const double infNormSol = mumps->RINFOG(5);
      const double residualScaled = mumps->RINFOG(6);
      const double omega1 = mumps->RINFOG(7);
      const double omega2 = mumps->RINFOG(8);

      if( rankMumps == 0 )
         printf("backward errors: %f, %f; scaled residual: %f; infNormMatrix: %f infNormSol %f.\n",
                  omega1, omega2, residualScaled, infNormMatrix, infNormSol);
   }
}

void
MumpsSolver::setUpMumps()
{
   mumps = new DMUMPS_STRUC_C;
   mumps->n = 0;
   mumps->nnz = 0;
   mumps->irn = NULL;
   mumps->jcn = NULL;
   mumps->a = NULL;
   memset(mumps->keep, 0, 400 * sizeof(int));
   mumps->comm_fortran = getFortranMPIComm(mpiCommMumps);
   mumps->sym = 2; //general symetric matrix
   mumps->job = -1; //initialization
   mumps->par = 1; //host process involved in parallel computations
   dmumps_c(mumps);

   if( mumps->INFOG(1) != 0 )
   {
      printf("Error occured when initializing MUMPS.\n");
      exit(1);
   }

  // setMumpsVerbosity(1); //3 moderately verbose, 1 only errors, 0 anything is surpressed

   mumps->n = (int) n;

   mumps->ICNTL(5) = 0; // 0: assembled format for matrix (triplet format)
   mumps->ICNTL(18) = 0; // 0: matrix centralized on rank 0

#if 1
   mumps->ICNTL(20) = 3; // 3: exploit sparsity during solve, 1: decide automatically
#else
   mumps->ICNTL(20) = 0; // right-hand side is in dense format
#endif

   mumps->ICNTL(21) = 0; // solution vector is assembled and stored in mumps structure component RHS
}
