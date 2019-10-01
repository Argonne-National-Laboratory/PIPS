/*
 * MumpsSolverLeaf.C
 */

//#define PIPS_DEBUG

#include "MumpsSolverLeaf.h"
#include <stdlib.h>
#include "SimpleVector.h"
#include "SparseGenMatrix.h"


MumpsSolverLeaf::MumpsSolverLeaf( SparseSymMatrix * sgm )
 : MumpsSolverBase(sgm)
{
}


MumpsSolverLeaf::~MumpsSolverLeaf()
{
}


void
MumpsSolverLeaf::matrixChanged()
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
MumpsSolverLeaf::solve(GenMatrix& rhs_f, int startRow, int range, double* sol)
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
   assert(ia_org[0] == 1 && rhs_matrix.getStorageRef().len == ia_org[m_org] - 1);


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

   MumpsSolverBase::solve(sol);

   delete[] ia_sub;
}


