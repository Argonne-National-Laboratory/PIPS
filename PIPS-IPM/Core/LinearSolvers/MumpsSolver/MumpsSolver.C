/*
 * MumpsSolver.C
 */

//#define PIPS_DEBUG

#include <stdlib.h>

#include "MumpsSolver.h"
#include "SimpleVector.h"
#include "SparseGenMatrix.h"


#define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation
#define INFOG(I) infog[(I)-1]
#define RINFOG(I) rinfog[(I)-1]



MumpsSolver::MumpsSolver( SparseSymMatrix * sgm )
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

   setUpMpiData(MPI_COMM_WORLD, MPI_COMM_SELF);
   setUpMumps();
}

MumpsSolver::~MumpsSolver()
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
MumpsSolver::diagonalChanged(int idiag, int extent)
{
   PIPSdebugMessage("diagonal changed \n");

   this->matrixChanged();
}

void
MumpsSolver::matrixChanged()
{
   PIPSdebugMessage("matrix changed \n");

   if( mpiCommMumps == MPI_COMM_NULL )
      return;

   // todo: update only diagonal!
   assert(Msys);

#ifdef TIME_Triplet_c2fortran
   const double t1 = MPI_Wtime();
#endif

   delete[] tripletIrn;
   delete[] tripletJcn;
   delete[] tripletA;
   tripletIrn = nullptr;
   tripletJcn = nullptr;
   tripletA = nullptr;

   Msys->getSparseTriplet_c2fortran(tripletIrn, tripletJcn, tripletA);

#ifdef TIME_Triplet_c2fortran
   const double t2 = MPI_Wtime();
   std::cout << "Triplet_c2fortran time=" << t2 - t1 << std::endl;
#endif

   mumps->n = n;
   mumps->nnz = Msys->numberOfNonZeros();
   mumps->irn = tripletIrn;
   mumps->jcn = tripletJcn;
   mumps->a = tripletA;

   // todo test whether 7 or 5 is better
   // symmetric permutation for factorization, 5: METIS, 7: automatic choice; meaningless if mumps->ICNTL(28) == 2
   mumps->ICNTL(7) = 5;

   mumps->ICNTL(28) = 0; // choice of analysis, 0: automatic, 1: sequential, 2: parallel

   mumps->ICNTL(29) = 0; // parallel ordering, 0: automatic, 1: PT-SCOTCH, 2: ParMetis

   // relative threshold for numerical pivoting; 0.01 is default for sym. indef., larger values increase accuracy
   //  mumps->ICNTL(1) = 0.01;

   // relative threshold for static pivoting; -1.0: not used (default), 0.0: use with automatic choice of threshold
   //  mumps->ICNTL(4) = -1.0;


   // analysis phase
   mumps->job = 1;

   double starttime = MPI_Wtime();

   // do analysis
   dmumps_c(mumps);

   processMumpsResultAnalysis(starttime);


   // factorization phase
   mumps->job = 2;

   starttime = MPI_Wtime();

   // do factorization
   dmumps_c(mumps);

   processMumpsResultFactor(starttime);

   // todo save permutation for reuse?
}


void
MumpsSolver::solve(double* vec)
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
MumpsSolver::solve(OoqpVector& rhs)
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
MumpsSolver::solve(GenMatrix& rhs_f, int startRow, int range, double* sol)
{
   PIPSdebugMessage("MUMPS solver: solve (multiple rhs) \n");

   assert(sol);
   assert(startRow >= 0 && range >= 1);

   SparseGenMatrix& rhs_matrix = dynamic_cast<SparseGenMatrix &>(rhs_f);

   if( mpiCommMumps == MPI_COMM_NULL )
      return;

   int m_org;
   int n_org;

   rhs_matrix.getSize(m_org, n_org);

   assert(startRow + range <= m_org);

   const int m_sub = range;
   const int n_sub = n_org;

   int* const ia_org = rhs_matrix.krowM();
   int* const ja_org = rhs_matrix.jcolM();
   double* const a_org = rhs_matrix.M();

   // matrix should be in Fortran format
   assert(ia_org[0] == 1 && rhs_matrix.getStorage()->len == ia_org[m_org] - 1);


   int* const ia_sub = new int[m_sub + 1];

   for( int i = 0; i <= m_sub; i++)
   {
      const int newPos = ia_org[i + startRow] - ia_org[startRow] + 1;
      assert(i == 0 || newPos >= ia_sub[i - 1]);

      ia_sub[i] = newPos;
   }

   assert(ia_sub[0] == 1 && ia_sub[m_sub] == (ia_org[startRow + range] - ia_org[startRow] + 1));

   mumps->nrhs = m_sub; // MUMPS expects column major
   mumps->nz_rhs = ia_sub[m_sub] - 1;
   mumps->lrhs = n_sub;
   mumps->irhs_ptr = ia_sub;
   mumps->irhs_sparse = &ja_org[ia_org[startRow] - 1];
   mumps->rhs_sparse = &a_org[ia_org[startRow] - 1];
   mumps->ICNTL(20) = 3; // exploit sparsity during solve

   solve(sol);

   delete[] ia_sub;
}



void
MumpsSolver::processMumpsResultAnalysis(double starttime)
{
   const int errorCode = mumps->INFOG(1);
   if( errorCode != 0 )
   {
      if( rankMumps == 0 )
         printf("Error INFOG(1)=%d in MUMPS analysis phase. \n", errorCode);

      exit(1);
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
MumpsSolver::processMumpsResultFactor(double starttime)
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
MumpsSolver::processMumpsResultSolve(double starttime)
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
MumpsSolver::setUpMpiData(MPI_Comm mpiCommPips_c, MPI_Comm mpiCommMumps_c)
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
