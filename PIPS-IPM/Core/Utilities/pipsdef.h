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

template <typename T>
struct get_mpi_datatype_t;

// specialization for particular types:
template <>
struct get_mpi_datatype_t<int> {
   static constexpr MPI_Datatype value = MPI_INT;
};

template <>
struct get_mpi_datatype_t<long> {
   static constexpr MPI_Datatype value = MPI_LONG;
};

template <>
struct get_mpi_datatype_t<long long> {
   static constexpr MPI_Datatype value = MPI_LONG_LONG;
};

template <>
struct get_mpi_datatype_t<double> {
   static constexpr MPI_Datatype value = MPI_DOUBLE;
};

template <>
struct get_mpi_datatype_t<char> {
   static constexpr MPI_Datatype value = MPI_CHAR;
};

template <>
struct get_mpi_datatype_t<bool> {
   static constexpr MPI_Datatype value = MPI_C_BOOL;
};

template <>
struct get_mpi_datatype_t<unsigned int> {
   static constexpr MPI_Datatype value = MPI_UNSIGNED;
};

template <typename T>
MPI_Datatype get_mpi_datatype(const T& arg) {
   return get_mpi_datatype_t<T>::value;
}

template <typename T>
MPI_Datatype get_mpi_datatype() {
   return get_mpi_datatype(T());
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

template <typename T>
inline T PIPS_MPIgetMin(T localmin, MPI_Comm mpiComm)
{
   T globalmin = 0.0;
   MPI_Allreduce(&localmin, &globalmin, 1, get_mpi_datatype(localmin), MPI_MIN, mpiComm);

   return globalmin;
}

template <typename T>
inline T PIPS_MPIgetMax(T localmax, MPI_Comm mpiComm)
{
   T globalmax = 0.0;
   MPI_Allreduce(&localmax, &globalmax, 1, get_mpi_datatype(localmax), MPI_MAX, mpiComm);

   return globalmax;
}

template <typename T>
inline T PIPS_MPIgetSum(T localsummand, MPI_Comm mpiComm)
{
   T sum;
   MPI_Allreduce(&localsummand, &sum, 1, get_mpi_datatype(localsummand), MPI_SUM, mpiComm);

   return sum;
}

template <typename T>
inline void PIPS_MPIsumArrayInPlace(T* elements, int length, MPI_Comm mpiComm)
{
   assert(length >= 0);

   if( length == 0 )
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements[0]), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPIsumArray(const T* source, T* dest, int length, MPI_Comm mpiComm )
{
   assert(length >= 0);

   if(length == 0)
      return;

   MPI_Allreduce(source, dest, length, get_mpi_datatype(source[0]), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArrayInPlace(T* elements, int length, MPI_Comm mpiComm)
{
   assert(length >= 0);
   if(length == 0)
      return;
   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements[0]), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArray(T* source, T* dest, int length, MPI_Comm mpiComm)
{
   assert(length >= 0);
   if(length == 0)
      return;
   MPI_Allreduce(source, dest, length, get_mpi_datatype(source[0]), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPIgatherv(T* sendbuf, int sendcnt, T* recvbuf, int* recvcnts, int* recvoffsets, int root, MPI_Comm mpiComm)
{
   // todo assert at least one element in any of the buffers.. 
   MPI_Gatherv(sendbuf, sendcnt, get_mpi_datatype(sendbuf[0]), recvbuf, recvcnts, recvoffsets, get_mpi_datatype(recvbuf[0]), root, mpiComm);
}

template <typename T>
inline void PIPS_MPIallgather( T* sendbuf, int sendcnt, T* recvbuf, int recvcnt, MPI_Comm mpiComm)
{
   // todo assert at least one element in any of the buffers.. 

   MPI_Allgather(sendbuf, sendcnt, get_mpi_datatype(recvbuf[0]), recvbuf, recvcnt, get_mpi_datatype(recvbuf[0]), mpiComm);
}

#endif
