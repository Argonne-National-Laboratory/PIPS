/*
 * pipschecks.C
 *
 *  Created on: 23.02.2018
 *      Author: bzfrehfe
 */

#include "pipschecks.h"
#include "pipsdef.h"
#include "OoqpVector.h"
#include "DoubleMatrixTypes.h"
#include "sData.h"
#include <cmath>
#include <limits>
#include <cassert>
#include <algorithm>

bool permutationIsValid(const std::vector<unsigned int>& perm)
{
   size_t size = perm.size();
   std::vector<bool> permflag(size, false);

   for( size_t i = 0; i < size; ++i )
   {
      const unsigned int p = perm[i];
      if( p >= size )
         return false;
      permflag[p] = true;
   }

   for( size_t i = 0; i < size; ++i )
      if( !permflag[i] )
         return false;

   return true;
}

bool subMatrixIsOrdered(const int* rowptr, const int* colidx,
      int rowstart, int rowend)
{
   assert(rowptr);
   assert(colidx);
   assert(rowstart >= 0 && rowstart <= rowend);

   for( int r = rowstart; r < rowend; r++ )
      for( int ci = rowptr[r] + 1; ci < rowptr[r + 1]; ci++ )
         if( colidx[ci - 1] >= colidx[ci] )
            return false;

   return true;
}

void computeFortranCSRMatResidualNorms(const int* rowptr, const int* colidx, const double* vals, /*const*/ SimpleVector& rhs,
      /*const*/ SimpleVector& x, double& res_norm2, double& res_nrmInf, double& sol_inf, double& mat_max)
{
   const int dim=rhs.length();
   double* tmp_resid = new double[dim];
   memcpy(tmp_resid, rhs.elements(), dim * sizeof(double));

   mat_max = 0.0, res_norm2 = 0.0, res_nrmInf = 0.0, sol_inf = 0.0;

   for( int i = 0; i < dim; i++ )
   {
      for( int p = rowptr[i]; p < rowptr[i + 1]; p++ )
      {
         const int j = colidx[p - 1] - 1;
         if( j + 1 <= dim )
         {
            //r[i] = r[i] + M(i,j)*x(j)
            tmp_resid[i] -= vals[p - 1] * x[j];

            if( fabs(vals[p - 1]) > mat_max )
               mat_max = fabs(vals[p - 1]);

            if( j != i )
            {
               //r[j] = r[j] + M(j,i)*x(i)
               tmp_resid[j] -= vals[p - 1] * x[i];
            }
         }
      }
   }

   for( int i = 0; i < dim; i++ )
   {
      res_norm2 += tmp_resid[i] * tmp_resid[i];
      if( res_nrmInf < fabs(tmp_resid[i]) )
         res_nrmInf = tmp_resid[i];
      if( fabs(x[i]) > sol_inf )
         sol_inf = fabs(x[i]);
   }
   res_norm2 = sqrt(res_norm2);

   delete[] tmp_resid;
}

// is root node data of sData object same on all procs?
bool rootNodeInSyncSData(const sData& s_data)
{
   int my_rank, world_size;
   MPI_Comm_rank( dynamic_cast<const StochGenMatrix&>(*s_data.A).mpiComm, &my_rank);
   MPI_Comm_size( dynamic_cast<const StochGenMatrix&>(*s_data.A).mpiComm, &world_size);

   if(my_rank == 0)
      std::cout << "checking if root node data over all processes is in sync" << std::endl;

   bool in_sync = true;

   /* matrix Q */
   // todo

   /* matrix A */
   if(!rootNodeInSyncStochGenMatrix(dynamic_cast<const StochGenMatrix&>(*s_data.A)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: matrix A corrupted!" << std::endl;
      in_sync = false;
   }

   /* matrix C */
   if(!rootNodeInSyncStochGenMatrix(dynamic_cast<const StochGenMatrix&>(*s_data.C)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: matrix C corrupted!" << std::endl;
      in_sync = false;
   }

   /* objective g */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.g)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: objective vector corrupted!" << std::endl;
      in_sync = false;
   }

   /* rhs equality bA */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.bA)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: rhs of A corrupted!" << std::endl;
      in_sync = false;
   }

   /* upper bounds x bux */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.bux)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: upper bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* index for upper bounds x ixupp */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.ixupp)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: index upper bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* lower bounds x blx */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.blx)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: lower bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* index for lower bounds x ixlow */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.ixlow)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: index lower bounds x corrupted!" << std::endl;
      in_sync = false;
   }

   /* upper bounds C bu */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.bu)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: rhs C corrupted!" << std::endl;
      in_sync = false;
   }

   /* index upper bounds C icupp */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.icupp)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: index rhs C corrupted!" << std::endl;
      in_sync = false;
   }

   /* lower bounds C bl */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.bl)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: lower bounds C corrupted!" << std::endl;
      in_sync = false;
   }

   /* index for lower bounds C iclow */
   if(!rootNodeInSyncStochVector(dynamic_cast<const StochVector&>(*s_data.iclow)))
   {
      if(my_rank == 0)
         std::cout << "ERROR: index lower bounds C corrupted!" << std::endl;
      in_sync = false;
   }

   /* sacle sc */
   // todo

   if(my_rank == 0 && in_sync)
      std::cout << "root node data over all processes is in sync" << std::endl;

   return in_sync;
}

// is root node data of StochVector same on all procs?
bool rootNodeInSyncStochVector(const StochVector& stoch_vec)
{
   bool in_sync = true;

   /* no need to check not distributed or not root node */
   if( !stoch_vec.iAmDistrib || stoch_vec.parent != NULL)
      return in_sync;

   assert( stoch_vec.vec );

   int my_rank, world_size;
   assert(stoch_vec.mpiComm == MPI_COMM_WORLD);
   MPI_Comm_rank(stoch_vec.mpiComm, &my_rank);
   MPI_Comm_size(stoch_vec.mpiComm, &world_size);

   const int vec_length = dynamic_cast<const SimpleVector&>(*stoch_vec.vec).length();
   const int vecl_length = (stoch_vec.vecl) ? dynamic_cast<const SimpleVector&>(*stoch_vec.vecl).length() : 0;

   const long long count = vec_length + vecl_length;

   assert( count < std::numeric_limits<int>::max());

   /* mpi reduce on vector */
   /* rank 0 recieves max and min over arrays and checks for equality */
   double sendbuf[count];
   double recvbuf_max[count];
   double recvbuf_min[count];

   const SimpleVector& vec = dynamic_cast<const SimpleVector&>(*stoch_vec.vec);
   std::copy(vec.elements(), vec.elements() + vec.length(), sendbuf);

   if(stoch_vec.vecl)
   {
      const SimpleVector& vecl = dynamic_cast<const SimpleVector&>(*stoch_vec.vecl);
      std::copy(vecl.elements(), vecl.elements() + vecl.length(), sendbuf + vec.length());
   }
   MPI_Reduce(sendbuf, recvbuf_max, static_cast<int>(count), MPI_DOUBLE, MPI_MAX, 0, stoch_vec.mpiComm);
   MPI_Reduce(sendbuf, recvbuf_min, static_cast<int>(count), MPI_DOUBLE, MPI_MIN, 0, stoch_vec.mpiComm);


   /* if rank == 0 check for sync */
   if(my_rank == 0)
   {
      for(int i = 0; i < count; ++i){
         if( !PIPSisEQ( recvbuf_max[i], recvbuf_min[i]) )
         {
            in_sync = false;
         }
      }
   }

   /* distribute result */
   MPI_Bcast(&in_sync, 1, MPI_CXX_BOOL, 0, stoch_vec.mpiComm);

   return in_sync;
}

// is root node data of StochMatrix same on all procs?
// not checking dynamic storage !
bool rootNodeInSyncStochGenMatrix(const StochGenMatrix& stoch_mat)
{
   bool in_sync = true;

   /* no need to check not distributed or not root node */
   if( !stoch_mat.iAmDistrib || stoch_mat.children.size() == 0)
      return in_sync;

   assert( stoch_mat.Amat );
   assert( stoch_mat.Bmat );
   assert( stoch_mat.Blmat );
   /* since we are in root node Amat should look as follows */
   assert( stoch_mat.Amat->getStorageRef().len == 0);
   assert( stoch_mat.Amat->getStorageRef().n == -1);

   int my_rank, world_size;
   MPI_Comm_rank(stoch_mat.mpiComm, &my_rank);
   MPI_Comm_size(stoch_mat.mpiComm, &world_size);


   const int lenght_entries_bmat = stoch_mat.Bmat->getStorageRef().len;
   const int length_columns_bmat = stoch_mat.Bmat->getStorageRef().len;
   const int lenght_rowoffest_bmat = stoch_mat.Bmat->getStorageRef().m + 1;

   const int lenght_entries_blmat = stoch_mat.Blmat->getStorageRef().len;
   const int length_columns_blmat = stoch_mat.Blmat->getStorageRef().len;
   const int lenght_rowoffest_blmat = stoch_mat.Blmat->getStorageRef().m + 1;

   const long long count_row_cols = length_columns_bmat + lenght_rowoffest_bmat + length_columns_blmat + lenght_rowoffest_blmat;
   const long long count_entries = lenght_entries_bmat + lenght_entries_blmat;

   assert( count_row_cols < std::numeric_limits<int>::max());
   assert( count_entries < std::numeric_limits<int>::max());

   /* mpi reduce on vector */
   /* rank 0 gets min and max of all arrays and compares them - if not equal -> corrupted*/
   double sendbuf_entries[count_entries];
   double recvbuf_entries_max[count_entries];
   double recvbuf_entries_min[count_entries];

   int sendbuf_row_col[count_row_cols];
   int recvbuf_row_col_max[count_row_cols];
   int recvbuf_row_col_min[count_row_cols];

   /* fill Bmat into send buffers */
   const double * M = stoch_mat.Bmat->getStorageRef().M;
   const int * krowM = stoch_mat.Bmat->getStorageRef().krowM;
   const int * jColM = stoch_mat.Bmat->getStorageRef().jcolM;

   std::copy(M, M + lenght_entries_bmat, sendbuf_entries);

   std::copy(krowM, krowM + lenght_rowoffest_bmat, sendbuf_row_col);
   std::copy(jColM, jColM + lenght_entries_bmat, sendbuf_row_col + lenght_rowoffest_bmat);

   /* fill Blmat into send buffers */
   const double * Ml = stoch_mat.Blmat->getStorageRef().M;
   const int * krowMl = stoch_mat.Blmat->getStorageRef().krowM;
   const int * jColMl = stoch_mat.Blmat->getStorageRef().jcolM;

   std::copy(Ml, Ml + lenght_entries_blmat, sendbuf_entries + lenght_entries_bmat);

   std::copy(krowMl, krowMl + lenght_rowoffest_blmat, sendbuf_row_col + lenght_rowoffest_bmat + lenght_entries_bmat);
   std::copy(jColMl, jColMl + lenght_entries_blmat, sendbuf_row_col + lenght_rowoffest_bmat + lenght_entries_bmat + lenght_rowoffest_blmat);


   /* Reduce Bmat and Blmat buffers */
   MPI_Reduce(sendbuf_entries, recvbuf_entries_max, static_cast<int>(count_entries), MPI_DOUBLE, MPI_MAX, 0, stoch_mat.mpiComm);
   MPI_Reduce(sendbuf_entries, recvbuf_entries_min, static_cast<int>(count_entries), MPI_DOUBLE, MPI_MIN, 0, stoch_mat.mpiComm);

   MPI_Reduce(sendbuf_row_col, recvbuf_row_col_max, static_cast<int>(count_row_cols), MPI_INT, MPI_MAX, 0, stoch_mat.mpiComm);
   MPI_Reduce(sendbuf_row_col, recvbuf_row_col_min, static_cast<int>(count_row_cols), MPI_INT, MPI_MIN, 0, stoch_mat.mpiComm);

   /* if rank == 0 check for sync */
   if(my_rank == 0)
   {
      /* check recvbuf_entries */
      for(int i = 0; i < count_entries; ++i){
         if( !PIPSisEQ(recvbuf_entries_max[i], recvbuf_entries_min[i]) )
         {
            in_sync = false;
         }
      }
      for(int i = 0; i < count_row_cols; ++i)
      {
         if( !PIPSisEQ(recvbuf_row_col_max[i], recvbuf_row_col_min[i]) )
         {
            in_sync = false;
         }
      }

   }

   /* distribute result */
   MPI_Bcast(&in_sync, 1, MPI_CXX_BOOL, 0, stoch_mat.mpiComm);

   return in_sync;
}
