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
#include <limits>
#include <assert.h>
#include "pipsport.h"

const double pips_eps = 1e-13;
const double pips_eps0 = 1e-40;

static const double feastol = 1.0e-6; // was 1.0e-6
static const double infinity = 1.0e30;
static const double eps_bounds_nontight = 1.0e-8;

static const double INF_NEG_PRES = -std::numeric_limits<double>::infinity();
static const double INF_POS_PRES = std::numeric_limits<double>::infinity();

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

inline bool PIPSisEQFeas( double val1, double val2 )
{
   return (std::fabs(val1-val2) <= feastol);
}

inline bool PIPSisEQ(double val1, double val2, double eps = pips_eps)
{
   return (std::fabs(val1 - val2) <= eps);
}

inline bool PIPSisRelEQ(double val1, double val2, double eps = pips_eps)
{
   const double reldiff = relativeDiff(val1, val2);

   return (std::fabs(reldiff) <= eps);
}

inline bool PIPSisRelEQFeas(double val1, double val2)
{
   const double reldiff = relativeDiff(val1, val2);

   return std::fabs(reldiff) <= feastol;
}

inline bool PIPSisLE(double val1, double val2, double eps = pips_eps)
{
   return (val1 <= val2 + eps);
}

inline bool PIPSisLEFeas(double val1, double val2)
{
   return (val1 <= val2 + feastol);
}

inline bool PIPSisRelLEFeas(double val1, double val2)
{
   const double reldiff = relativeDiff(val1, val2);

   return reldiff <= feastol;
}

inline bool PIPSisLT(double val1, double val2, double eps = pips_eps)
{
   return (val1 < val2 - eps);
}

inline bool PIPSisLTFeas(double val1, double val2)
{
   return (val1 < val2 - feastol);
}

inline bool PIPSisRelLT(double val1, double val2, double eps = pips_eps)
{
   const double reldiff = relativeDiff(val1, val2);

   return (reldiff < -eps);
}

inline bool PIPSisRelLTFeas(double val1, double val2)
{
   const double reldiff = relativeDiff(val1, val2);

   return (reldiff < -feastol);
}

inline bool PIPSisZero(double val, double eps0 = pips_eps0)
{
   return (std::fabs(val) < eps0);
}

inline bool PIPSisZeroFeas(double val)
{
   return (std::fabs(val) < feastol);
}

inline int PIPSgetnOMPthreads()
{
   int num_procs;

   /* Numbers of processors, value of OMP_NUM_THREADS */
   char* var = getenv("OMP_NUM_THREADS");
   if( var != nullptr )
      sscanf(var, "%d", &num_procs);
   else
   {
      printf("Set environment OMP_NUM_THREADS to 1");
      exit(1);
   }

   return num_procs;
}


inline bool PIPS_MPIiAmSpecial(int iAmDistrib, MPI_Comm mpiComm)
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

inline bool iAmSpecial(int iAmDistrib, MPI_Comm mpiComm)
{
   return PIPS_MPIiAmSpecial(iAmDistrib, mpiComm);
}

inline bool PIPS_MPIgetDistributed(MPI_Comm comm)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);

   if( world_size > 1)
      return true;
   else
      return false;
}

void inline PIPS_MPIabortInfeasible(MPI_Comm comm, std::string message, std::string file, std::string function)
{
   std::cerr << "Infesibility detected in " << file << " function " << function << "!" << std::endl;
   std::cerr << "Message: " << message << std::endl;
   std::cerr << "Aborting now." << std::endl;
   MPI_Abort(comm, 1);
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
   static constexpr MPI_Datatype value = MPI_CXX_BOOL;
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
MPI_Datatype get_mpi_datatype(T* arg) {
   return get_mpi_datatype_t<T>::value;
}

template <typename T>
MPI_Datatype get_mpi_datatype(const T* arg) {
   return get_mpi_datatype_t<T>::value;
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
inline bool PIPS_MPIisValueEqual(const T& val, MPI_Comm mpiComm)
{
   // todo make one vec and + - val
   const int max = PIPS_MPIgetMax(val, MPI_COMM_WORLD);;
   const int min = PIPS_MPIgetMin(val, MPI_COMM_WORLD);
   return (max == min);
}

template <typename T>
inline T PIPS_MPIgetMin(const T& localmin, MPI_Comm mpiComm)
{
   T globalmin = 0.0;
   MPI_Allreduce(&localmin, &globalmin, 1, get_mpi_datatype(localmin), MPI_MIN, mpiComm);

   return globalmin;
}

template <typename T>
inline T PIPS_MPIgetMax(const T& localmax, MPI_Comm mpiComm)
{
   T globalmax = 0.0;
   MPI_Allreduce(&localmax, &globalmax, 1, get_mpi_datatype(localmax), MPI_MAX, mpiComm);

   return globalmax;
}

template <typename T>
void PIPS_MPIgetMaxInPlace(T& max, MPI_Comm mpiComm)
{
   MPI_Allreduce(MPI_IN_PLACE, &max, 1, get_mpi_datatype(max), MPI_MAX, mpiComm);
}

template <typename T>
void PIPS_MPIgetMinInPlace(T& min, MPI_Comm mpiComm)
{
   MPI_Allreduce(MPI_IN_PLACE, &min, 1, get_mpi_datatype(min), MPI_MIN, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArrayInPlace(T* localmax, int length, MPI_Comm mpiComm)
{
   assert(length >= 0);
   if(length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, localmax, length, get_mpi_datatype(localmax), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm)
{
   if(elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPIminArrayInPlace(T* elements, int length, MPI_Comm mpiComm)
{
   assert(length >= 0);
   if(length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_MIN, mpiComm);
}

template <typename T>
inline void PIPS_MPIminArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm)
{
   if(elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_MIN, mpiComm);
}

template <typename T>
inline T PIPS_MPIgetSum(const T& localsummand, MPI_Comm mpiComm)
{
   T sum;
   MPI_Allreduce(&localsummand, &sum, 1, get_mpi_datatype(localsummand), MPI_SUM, mpiComm);

   return sum;
}

template <typename T>
inline void PIPS_MPIgetLogicOrInPlace(T& localval, MPI_Comm mpiComm)
{
   MPI_Allreduce(MPI_IN_PLACE, &localval, 1, get_mpi_datatype(localval), MPI_LOR, mpiComm);
}

template <typename T>
inline T PIPS_MPIgetLogicOr(const T& localval, MPI_Comm mpiComm)
{
   T lor;
   MPI_Allreduce(&localval, &lor, 1, get_mpi_datatype(localval), MPI_LOR, mpiComm);

   return lor;
}

template <typename T>
inline void PIPS_MPIlogicOrArrayInPlace(T* elements, int length, MPI_Comm mpiComm)
{
   assert( length >= 0 );
   if(length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_LOR, mpiComm);
}

template <typename T>
inline void PIPS_MPIlogicOrArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm)
{
   if(elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_LOR, mpiComm);
}

template <typename T>
inline void PIPS_MPIgetSumInPlace(T& sum, MPI_Comm mpiComm)
{
   MPI_Allreduce(MPI_IN_PLACE, &sum, 1, get_mpi_datatype(sum), MPI_SUM, mpiComm);
}

template <typename T>

inline void PIPS_MPIsumArrayInPlace(T* elements, int length, MPI_Comm mpiComm)
{
   assert(length >= 0);

   if( length == 0 )
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPIsumArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm)
{
   if( elements.size() == 0 )
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPIsumArray(const T* source, T* dest, int length, MPI_Comm mpiComm )
{
   assert(length >= 0);

   if(length == 0)
      return;

   MPI_Allreduce(source, dest, length, get_mpi_datatype(source), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArray(const T* source, T* dest, int length, MPI_Comm mpiComm)
{
   assert(length >= 0);
   if(length == 0)
      return;
   MPI_Allreduce(source, dest, length, get_mpi_datatype(source), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPIgatherv(const T* sendbuf, int sendcnt, T* recvbuf, int* recvcnts, const int* recvoffsets, int root, MPI_Comm mpiComm)
{
   assert(sendcnt >= 0);
   MPI_Gatherv(sendbuf, sendcnt, get_mpi_datatype(sendbuf), recvbuf, recvcnts, recvoffsets, get_mpi_datatype(recvbuf), root, mpiComm);
}

template <typename T>
inline void PIPS_MPIallgather( const T* sendbuf, int sendcnt, T* recvbuf, int recvcnt, MPI_Comm mpiComm)
{
   assert(sendcnt >= 0);
   assert(recvcnt >= 0);
   MPI_Allgather(sendbuf, sendcnt, get_mpi_datatype(recvbuf), recvbuf, recvcnt, get_mpi_datatype(recvbuf), mpiComm);
}

#endif
