/*
 * MumpsSolverRoot.C
 *
 *  Created on: 25.06.2019
 *      Author: bzfrehfe
 */


#include "MumpsSolverRoot.h"
#include "SimpleVector.h"
#include <stdlib.h>


MumpsSolverRoot::MumpsSolverRoot( MPI_Comm mpiComm, SparseSymMatrix * sgm )
 : MumpsSolverBase(mpiComm, mpiComm, sgm)
{
   assert(rankMumps == rankPips);
}


MumpsSolverRoot::MumpsSolverRoot( SparseSymMatrix * sgm )
 : MumpsSolverBase(sgm)
{
}

MumpsSolverRoot::~MumpsSolverRoot()
{
}

void
MumpsSolverRoot::matrixRebuild( DoubleMatrix& matrixNew )
{
   if( mpiCommMumps == MPI_COMM_NULL )
      return;

   Msys = &dynamic_cast<SparseSymMatrix&>(matrixNew);
   assert(n == Msys->size());
   assert(Msys->getStorageRef().fortranIndexed());

   delete[] tripletIrn;
   delete[] tripletJcn;
   delete[] tripletA;
   tripletIrn = nullptr;
   tripletJcn = nullptr;
   tripletA = nullptr;

   if( rankPips == 0 )
   {
      // todo!
#if 1
      Msys->getStorageRef().fortran2c();
      Msys->getSparseTriplet_c2fortran(tripletIrn, tripletJcn, tripletA);
#else
      Msys->getSparseTriplet_fortran2fortran(tripletIrn, tripletJcn, tripletA);
#endif
   }

   this->factorize();
}


void
MumpsSolverRoot::matrixChanged()
{
   PIPSdebugMessage("matrix changed \n");

   if( mpiCommMumps == MPI_COMM_NULL )
      return;

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

   this->factorize();
}

void
MumpsSolverRoot::factorize()
{
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
}

void
MumpsSolverRoot::solve(OoqpVector& rhs)
{
   PIPSdebugMessage("MUMPS solver: solve (single rhs) \n");

   SimpleVector& sv = dynamic_cast<SimpleVector &>(rhs);

   if( mpiCommMumps != MPI_COMM_NULL )
   {
      assert(n == rhs.length());
      int sizeMPI; MPI_Comm_size(mpiCommPips, &sizeMPI);

      mumps->nrhs = 1;
      mumps->lrhs = n;

      mumps->ICNTL(20) = 0; // right-hand side is in dense format
      // todo try sparse also for single rhs?
      MumpsSolverBase::solve(sv.elements());

      if( sizeMPI > 0 )
         MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, mpiCommPips);
   }
}
