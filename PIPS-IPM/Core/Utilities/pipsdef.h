/*
 * pipsdef.h
 *
 *  Created on: 29.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_UTILITIES_PIPSDEF_H_
#define PIPS_IPM_CORE_UTILITIES_PIPSDEF_H_

#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <mpi.h>
#include <assert.h>

const double pips_eps = 1e-13;
const double pips_eps0 = 1e-40;

static inline double relativeDiff(double val1, double val2)
{
   const double val1Abs = std::fabs(val1);
   const double val2Abs = std::fabs(val2);
   const double div = std::max(1.0, std::max(val1Abs, val2Abs));

   return (val1 - val2) / div;
}

#ifdef PIPS_DEBUG
#define PIPSdebugMessage                printf("[%s:%d] debug: ", __FILE__, __LINE__), printf
#else
#define PIPSdebugMessage      while( 0 ) /*lint -e{530}*/ printf
#endif

inline bool PIPSisEQ(double val1, double val2, double eps = pips_eps)
{
   return (std::fabs(val1 - val2) <= eps);
}

inline bool PIPSisRelEQ(double val1, double val2, double eps = pips_eps)
{
   const double reldiff = relativeDiff(val1, val2);

   return (std::fabs(reldiff) <= eps);
}

inline bool PIPSisLE(double val1, double val2, double eps = pips_eps)
{
   return (val1 <= val2 + eps);
}

inline bool PIPSisLT(double val1, double val2, double eps = pips_eps)
{
   return (val1 < val2 - eps);
}

inline bool PIPSisRelLT(double val1, double val2, double eps = pips_eps)
{
   const double reldiff = relativeDiff(val1, val2);

   return (reldiff < -eps);
}

inline bool PIPSisZero(double val, double eps0 = pips_eps0)
{
   return (std::fabs(val) < eps0);
}

inline int PIPSgetnOMPthreads()
{
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

   return num_procs;
}

inline bool iAmSpecial(int iAmDistrib, MPI_Comm mpiComm)
{
   bool iAmSpecial = true;

   if( iAmDistrib )
   {
      int rank;

      MPI_Comm_rank(mpiComm, &rank);
      if( rank > 0 )
         iAmSpecial = false;
   }

   return iAmSpecial;
}

inline std::vector<int> PIPSallgathervInt(const std::vector<int>& vecLocal, MPI_Comm mpiComm)
{
   int myrank; MPI_Comm_rank(mpiComm, &myrank);
   int mysize; MPI_Comm_size(mpiComm, &mysize);

   std::vector<int> vecGathered;

   if( mysize > 0 )
   {
      // get all lengths
      std::vector<int> recvcounts(mysize);
      std::vector<int> recvoffsets(mysize);

      int lengthLocal = int(vecLocal.size());

      MPI_Allgather(&lengthLocal, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, mpiComm);

      // all-gather local components
      recvoffsets[0] = 0;
      for( int i = 1; i < mysize; ++i )
         recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

      const int lengthGathered = recvoffsets[mysize - 1] + recvcounts[mysize - 1];

      vecGathered = std::vector<int>(lengthGathered);

      MPI_Allgatherv(&vecLocal[0], lengthLocal, MPI_INT, &vecGathered[0],
            &recvcounts[0], &recvoffsets[0], MPI_INT, mpiComm);
   }
   else
   {
      vecGathered = vecLocal;
   }

   return vecGathered;
}

inline std::vector<int> PIPSallgathervInt(const std::vector<int>& vecLocal, MPI_Comm mpiComm, int& startMy, int& endMy)
{
   int myrank; MPI_Comm_rank(mpiComm, &myrank);
   int mysize; MPI_Comm_size(mpiComm, &mysize);

   std::vector<int> vecGathered;

   if( mysize > 0 )
   {
      // get all lengths
      std::vector<int> recvcounts(mysize);
      std::vector<int> recvoffsets(mysize);

      int lengthLocal = int(vecLocal.size());

      MPI_Allgather(&lengthLocal, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, mpiComm);

      // all-gather local components
      recvoffsets[0] = 0;
      for( int i = 1; i < mysize; ++i )
         recvoffsets[i] = recvoffsets[i - 1] + recvcounts[i - 1];

      const int lengthGathered = recvoffsets[mysize - 1] + recvcounts[mysize - 1];

      startMy = recvoffsets[myrank];

      // not the last rank?
      if( myrank < mysize - 1 )
      {
         endMy = recvoffsets[myrank + 1];
      }
      else
      {
         assert(myrank == mysize - 1 );
         endMy = lengthGathered;
      }

      vecGathered = std::vector<int>(lengthGathered);

      MPI_Allgatherv(&vecLocal[0], lengthLocal, MPI_INT, &vecGathered[0],
            &recvcounts[0], &recvoffsets[0], MPI_INT, mpiComm);
   }
   else
   {
      startMy = 0;
      endMy = int(vecLocal.size());
      vecGathered = vecLocal;
   }

   return vecGathered;
}

inline int PIPS_MPIgetRank(MPI_Comm mpiComm)
{
   int myrank;
   MPI_Comm_rank(mpiComm, &myrank);
   return myrank;
}

inline int PIPS_MPIgetSize(MPI_Comm mpiComm)
{
   int mysize;
   MPI_Comm_size(mpiComm, &mysize);
   return mysize;
}

inline double PIPS_MPIgetMin(double localmin, MPI_Comm mpiComm)
{
   double globalmin = 0.0;
   MPI_Allreduce(&localmin, &globalmin, 1, MPI_DOUBLE, MPI_MIN, mpiComm);

   return globalmin;
}

inline double PIPS_MPIgetMax(double localmax, MPI_Comm mpiComm)
{
   double globalmax = 0.0;
   MPI_Allreduce(&localmax, &globalmax, 1, MPI_DOUBLE, MPI_MAX, mpiComm);

   return globalmax;
}

inline double PIPS_MPIgetSum(double localsummand, MPI_Comm mpiComm)
{
   double sum;
   MPI_Allreduce(&localsummand, &sum, 1, MPI_DOUBLE, MPI_SUM, mpiComm);

   return sum;
}

inline void PIPS_MPIsumArray(MPI_Comm mpiComm, int length, double* elements)
{
   assert(length >= 0);

   if( length == 0 )
      return;

   double* buffer = new double[length];

   MPI_Allreduce(elements, buffer, length, MPI_DOUBLE, MPI_SUM, mpiComm);
   memcpy(elements, buffer, length * sizeof( double ));

   delete[] buffer;
}


#endif
