/*
 * MumpsSolverBase.C
 *
 *  Created on: 25.06.2019
 *      Author: bzfrehfe
 */

#include <stdlib.h>
#include "MumpsSolverBase.h"
#include "SimpleVector.h"
#include "SparseGenMatrix.h"


MumpsSolverBase::MumpsSolverBase( MPI_Comm mpiCommPips_c, MPI_Comm mpiCommMumps_c, SparseSymMatrix * sgm )
 : verbosity(defaultVerbosity), maxNiterRefinments(defaultMaxNiterRefinments)
{
   PIPSdebugMessage("creating MUMPS solver \n");

   assert(sgm);
   assert(sizeof(MUMPS_INT) == sizeof(int));

   Msys = sgm;
   n = sgm->size();
   tripletIrn = nullptr;
   tripletJcn = nullptr;
   tripletA = nullptr;

   setUpMpiData(mpiCommPips_c, mpiCommMumps_c);
   setUpMumps();
}

MumpsSolverBase::MumpsSolverBase( SparseSymMatrix * sgm )
 : MumpsSolverBase(MPI_COMM_WORLD, MPI_COMM_SELF, sgm)
{

}

MumpsSolverBase::~MumpsSolverBase()
{
   PIPSdebugMessage("deleting MUMPS solver \n");

   if( mumps )
   {
      mumps->job = -2;
      dmumps_c (mumps);

      delete mumps;
   }

   delete[] tripletIrn;
   delete[] tripletJcn;
   delete[] tripletA;
}

void
MumpsSolverBase::diagonalChanged(int idiag, int extent)
{
   PIPSdebugMessage("diagonal changed \n");

   this->matrixChanged();
}


void
MumpsSolverBase::solve(double* vec)
{
   assert(vec);
   assert(mpiCommMumps != MPI_COMM_NULL);

   mumps->rhs = vec;

   // solution phase
   mumps->job = 3;

   mumps->ICNTL(10) = maxNiterRefinments; // maximum number of it. refinements;

   if( verbosity == verb_mute ) // todo only print statistics for high verb
      mumps->ICNTL(11) = 0; // error statistics, 0: disabled, 2: main statistics
   else
      mumps->ICNTL(11) = 2; // error statistics, 0: disabled, 2: main statistics

   if( mumps->nrhs > 1 )
   {
      mumps->ICNTL(10) = 0;
      mumps->ICNTL(11) = 0;
   }

   const double starttime = MPI_Wtime();

   dmumps_c(mumps);

   processMumpsResultSolve(starttime);
}


void
MumpsSolverBase::solve(OoqpVector& rhs)
{
   PIPSdebugMessage("MUMPS solver: solve (single rhs) \n");

   SimpleVector& sv = dynamic_cast<SimpleVector &>(rhs);

   if( mpiCommMumps != MPI_COMM_NULL )
   {
      assert(n == rhs.length());

      mumps->nrhs = 1;
      mumps->lrhs = n;

      mumps->ICNTL(20) = 0; // right-hand side is in dense format
      // todo try sparse also for single rhs?
      solve(sv.elements());
   }
}


void
MumpsSolverBase::processMumpsResultAnalysis(double starttime)
{
   const int errorCode = mumps->INFOG(1);
   if( errorCode != 0 )
   {
      if( errorCode < 0 )
      {
         if( rankMumps == 0 )
         printf("Error INFOG(1)=%d INFOG(2)=%d in MUMPS analysis phase. \n", errorCode, mumps->INFOG(2));

         exit(1);
      }

      if( rankMumps == 0 )
         printf("Warning INFOG(1)=%d INFOG(2)=%d in MUMPS analysis phase. \n", errorCode, mumps->INFOG(2));
   }

   if( verbosity != verb_mute )
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
MumpsSolverBase::processMumpsResultFactor(double starttime)
{
   int errorCode = mumps->INFOG(1);

   if( errorCode != 0 )
   {
      if( errorCode == -8 || errorCode == -9 )
      {
         for( int i = 0; i < maxNreallocs; i++ )
         {
            mumps->ICNTL(14) *= 2;

            if( rankMumps == 0 && verbosity != verb_mute )
               printf("Increased MUMPS memory parameter ICNTL(14) to %d \n", mumps->ICNTL(14));

            dmumps_c(mumps);
            errorCode = mumps->INFOG(1);

            // todo would also need to store permutation here

            if( errorCode != -8 && errorCode != -9 )
               break;
         }

         if( errorCode == -8 || errorCode == -9 )
         {
            if( rankMumps == 0 )
               printf("Fatal error in Mumps: not able to obtain more memory (ICNTL(14)-related)\n");
            exit(1);
         }
      }

      if( errorCode != 0  )
      {
         if( rankMumps == 0 )
            printf("Error INFOG(1)=%d in MUMPS facorization phase. \n", errorCode);
         exit(1);
      }
   }

   if( verbosity != verb_mute  && starttime >= 0.0  )
   {
      const double timeFactor = MPI_Wtime() - starttime;

      if( rankMumps == 0 )
         printf("MUMPS factorization phase took %g seconds.\n", timeFactor);
   }
}

void
MumpsSolverBase::processMumpsResultSolve(double starttime)
{
   const int errorCode = mumps->INFOG(1);
   if( errorCode != 0 )
   {
      if( rankMumps == 0 )
         printf("Error INFOG(1)=%d in MUMPS solution phase. \n", errorCode);

      exit(1);
   }

   if(  verbosity != verb_mute   )
   {
      if( starttime >= 0.0 )
      {
         const double timeSolution = MPI_Wtime() - starttime;

         if( rankMumps == 0 )
            printf("MUMPS solution phase took %g seconds.\n", timeSolution);
      }

      // statistics computed?
      if( mumps->ICNTL(11) != 0 )
      {
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
}

void
MumpsSolverBase::setUpMpiData(MPI_Comm mpiCommPips_c, MPI_Comm mpiCommMumps_c)
{
   rankMumps = -1;
   this->mpiCommPips = mpiCommPips_c;
   this->mpiCommMumps = mpiCommMumps_c;

   if( mpiCommMumps != MPI_COMM_NULL )
   {
      MPI_Comm_rank(mpiCommMumps, &rankMumps);
      assert(rankMumps >= 0);
   }

   MPI_Comm_rank(mpiCommPips, &rankPips);
}


void
MumpsSolverBase::setUpMumps()
{
   mumps = new DMUMPS_STRUC_C;
   mumps->n = 0;
   mumps->nnz = 0;
   mumps->irn = NULL;
   mumps->jcn = NULL;
   mumps->a = NULL;
   memset(mumps->keep, 0, 400 * sizeof(int));
   mumps->comm_fortran = getFortranMPIComm(mpiCommMumps);
   mumps->sym = 2; // general symmetric matrix
   mumps->job = -1; // initialization
   mumps->par = 1; // host process involved in parallel computations
   dmumps_c(mumps);

   if( mumps->INFOG(1) != 0 )
   {
      printf("Error occured when initializing MUMPS.\n");
      exit(1);
   }

   mumps->n = static_cast<int>(n);

   mumps->ICNTL(5) = 0; // 0: assembled format for matrix (triplet format)
   mumps->ICNTL(14) = 50; // percentage increase in the estimated working space; default: 20
   mumps->ICNTL(18) = 0; // 0: matrix centralized on rank 0
   mumps->ICNTL(21) = 0; // solution vector is assembled and stored in MUMPS structure member RHS

   if( verbosity == verb_mute )
      mumps->ICNTL(4) = 0;  // nothing printed
   else if( verbosity == verb_standard )
      mumps->ICNTL(4) = 2; // errors, warnings, and main statistics printed
   else
   {
      assert(verbosity == verb_high);
      mumps->ICNTL(4) = 3; // errors and warnings and terse diagnostics (only first ten entries of arrays) printed.
   }
}
