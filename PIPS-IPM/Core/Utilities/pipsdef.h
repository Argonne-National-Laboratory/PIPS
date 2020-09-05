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


inline bool PIPS_MPIiAmSpecial(int iAmDistrib, MPI_Comm mpiComm = MPI_COMM_WORLD)
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

inline bool iAmSpecial(int iAmDistrib, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   return PIPS_MPIiAmSpecial(iAmDistrib, mpiComm);
}

inline bool PIPS_MPIgetDistributed(MPI_Comm comm = MPI_COMM_WORLD)
{
   int world_size;
   MPI_Comm_size(comm, &world_size);

   if( world_size > 1)
      return true;
   else
      return false;
}

void inline PIPS_MPIabortInfeasible(std::string message, std::string file, std::string function, MPI_Comm comm = MPI_COMM_WORLD)
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
struct get_mpi_datatype_t
{
   static inline MPI_Datatype get();
};

// specialization for particular types:
template <>
inline MPI_Datatype get_mpi_datatype_t<int>::get()
{
   return MPI_INT;
};

template <>
inline MPI_Datatype get_mpi_datatype_t<long>::get()
{
   return MPI_LONG;
};

template <>
inline MPI_Datatype get_mpi_datatype_t<long long>::get()
{
   return MPI_LONG_LONG;
};

template <>
inline MPI_Datatype get_mpi_datatype_t<double>::get()
{
   return MPI_DOUBLE;
};

template <>
inline MPI_Datatype get_mpi_datatype_t<char>::get()
{
   return MPI_CHAR;
};

template <>
inline MPI_Datatype get_mpi_datatype_t<bool>::get()
{
   return MPI_CXX_BOOL;
};

template <>
inline MPI_Datatype get_mpi_datatype_t<unsigned int>::get()
{
   return MPI_UNSIGNED;
};

template <typename T>
inline MPI_Datatype get_mpi_datatype(const T& arg) {
   return get_mpi_datatype_t<T>::get();
};

template <typename T>
inline MPI_Datatype get_mpi_datatype(T* arg) {
   return get_mpi_datatype_t<T>::get();
};

template <typename T>
inline MPI_Datatype get_mpi_datatype(const T* const arg) {
   return get_mpi_datatype_t<T>::get();
};

template <typename T>
struct get_mpi_locdatatype_t
{
   static inline MPI_Datatype get();
};

template <>
inline MPI_Datatype get_mpi_locdatatype_t<float>::get()
{
      return MPI_FLOAT_INT;
};

template <>
inline MPI_Datatype get_mpi_locdatatype_t<double>::get()
{
      return MPI_DOUBLE_INT;
};

template <>
inline MPI_Datatype get_mpi_locdatatype_t<long>::get()
{
      return MPI_LONG_INT;
};

template <>
inline MPI_Datatype get_mpi_locdatatype_t<int>::get()
{
      return MPI_2INT;
};

template <>
inline MPI_Datatype get_mpi_locdatatype_t<short>::get()
{
      return MPI_SHORT_INT;
};

template <>
inline MPI_Datatype get_mpi_locdatatype_t<long double>::get()
{
      return MPI_LONG_DOUBLE_INT;
};

template <typename T>
inline MPI_Datatype get_mpi_locdatatype(const T& arg)
{
   return get_mpi_locdatatype_t<T>::get();
};

template <typename T>
inline MPI_Datatype get_mpi_locdatatype(const T* const arg)
{
   return get_mpi_locdatatype_t<T>::get();
};

template <typename T>
inline MPI_Datatype get_mpi_locdatatype(T* arg)
{
   return get_mpi_locdatatype_t<T>::get();
};


inline int PIPS_MPIgetRank(MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   int myrank;
   MPI_Comm_rank(mpiComm, &myrank);
   return myrank;
}

inline int PIPS_MPIgetSize(MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   int mysize;
   MPI_Comm_size(mpiComm, &mysize);
   return mysize;
}

template <typename T>
inline std::vector<std::pair<T, int> > PIPS_MPIminlocArray(const T* localmin, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if( length <= 0 )
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T,int> > pairs(length);

   for(unsigned int i = 0; i < pairs.size(); ++i)
   {
      pairs[i].first = localmin[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], length, get_mpi_locdatatype(localmin), MPI_MINLOC, mpiComm);

   return pairs;
}

template <typename T>
inline std::vector<std::pair<T, int> > PIPS_MPIminlocArray(const std::vector<T>& localmin, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if( localmin.size() == 0 )
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T,int> > pairs(localmin.size());

   for(unsigned int i = 0; i < pairs.size(); ++i)
   {
      pairs[i].first = localmin[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], localmin.size(), get_mpi_locdatatype(localmin[0]), MPI_MINLOC, mpiComm);

   return pairs;
}

template <typename T>
inline std::vector<std::pair<T, int> > PIPS_MPImaxlocArray(const T* localmax, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if( length <= 0 )
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T,int> > pairs(length);

   for(unsigned int i = 0; i < pairs.size(); ++i)
   {
      pairs[i].first = localmax[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], length, get_mpi_locdatatype(localmax), MPI_MAXLOC, mpiComm);

   return pairs;
}

template <typename T>
inline std::vector<std::pair<T, int> > PIPS_MPImaxlocArray(const std::vector<T>& localmax, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if( localmax.size() == 0 )
      return std::vector<std::pair<T, int> >();

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   std::vector<std::pair<T,int> > pairs(localmax.size());

   for(unsigned int i = 0; i < pairs.size(); ++i)
   {
      pairs[i].first = localmax[i];
      pairs[i].second = my_rank;
   }

   MPI_Allreduce(MPI_IN_PLACE, &pairs[0], localmax.size(), get_mpi_locdatatype(localmax[0]), MPI_MAXLOC, mpiComm);

   return pairs;
}

template <typename T>
inline T PIPS_MPIgetMin(const T& localmin, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   T globalmin = 0.0;
   MPI_Allreduce(&localmin, &globalmin, 1, get_mpi_datatype(localmin), MPI_MIN, mpiComm);

   return globalmin;
}

template <typename T>
inline T PIPS_MPIgetMax(const T& localmax, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   T globalmax = 0.0;
   MPI_Allreduce(&localmax, &globalmax, 1, get_mpi_datatype(localmax), MPI_MAX, mpiComm);

   return globalmax;
}

template <typename T>
void PIPS_MPIgetMaxInPlace(T& max, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   MPI_Allreduce(MPI_IN_PLACE, &max, 1, get_mpi_datatype(max), MPI_MAX, mpiComm);
}

template <typename T>
void PIPS_MPIgetMinInPlace(T& min, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   MPI_Allreduce(MPI_IN_PLACE, &min, 1, get_mpi_datatype(min), MPI_MIN, mpiComm);
}

template <typename T>
inline bool PIPS_MPIisValueEqual(const T& val, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   // todo make one vec and + - val
   const int max = PIPS_MPIgetMax(val, MPI_COMM_WORLD);;
   const int min = PIPS_MPIgetMin(val, MPI_COMM_WORLD);
   return (max == min);
}

template <typename T>
inline void PIPS_MPImaxArrayInPlace(T* localmax, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(length >= 0);
   if(length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, localmax, length, get_mpi_datatype(localmax), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if(elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPIminArrayInPlace(T* elements, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(length >= 0);
   if(length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_MIN, mpiComm);
}

template <typename T>
inline void PIPS_MPIminArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if(elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_MIN, mpiComm);
}

template <typename T>
inline T PIPS_MPIgetSum(const T& localsummand, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   T sum;
   MPI_Allreduce(&localsummand, &sum, 1, get_mpi_datatype(localsummand), MPI_SUM, mpiComm);

   return sum;
}

template <typename T>
inline void PIPS_MPIgetLogicOrInPlace(T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   MPI_Allreduce(MPI_IN_PLACE, &localval, 1, get_mpi_datatype(localval), MPI_LOR, mpiComm);
}

template <typename T>
inline T PIPS_MPIgetLogicOr(const T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   T lor;
   MPI_Allreduce(&localval, &lor, 1, get_mpi_datatype(localval), MPI_LOR, mpiComm);

   return lor;
}

template <typename T>
inline void PIPS_MPIgetLogicAndInPlace(T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   MPI_Allreduce(MPI_IN_PLACE, &localval, 1, get_mpi_datatype(localval), MPI_LAND, mpiComm);
}

template <typename T>
inline T PIPS_MPIgetLogicAnd(const T& localval, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   T land;
   MPI_Allreduce(&localval, &land, 1, get_mpi_datatype(localval), MPI_LAND, mpiComm);

   return land;
}

template <typename T>
inline void PIPS_MPIlogicOrArrayInPlace(T* elements, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert( length >= 0 );
   if(length == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_LOR, mpiComm);
}

template <typename T>
inline void PIPS_MPIlogicOrArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if(elements.size() == 0)
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_LOR, mpiComm);
}

template <typename T>
inline void PIPS_MPIgetSumInPlace(T& sum, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   MPI_Allreduce(MPI_IN_PLACE, &sum, 1, get_mpi_datatype(sum), MPI_SUM, mpiComm);
}

template <typename T>

inline void PIPS_MPIsumArrayInPlace(T* elements, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(length >= 0);

   if( length == 0 )
      return;

   MPI_Allreduce(MPI_IN_PLACE, elements, length, get_mpi_datatype(elements), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPIsumArrayInPlace(std::vector<T>& elements, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   if( elements.size() == 0 )
      return;

   MPI_Allreduce(MPI_IN_PLACE, &elements[0], elements.size(), get_mpi_datatype(&elements[0]), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPIsumArray(const T* source, T* dest, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(length >= 0);

   if(length == 0)
      return;

   MPI_Allreduce(source, dest, length, get_mpi_datatype(source), MPI_SUM, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArray(const T* source, T* dest, int length, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(length >= 0);
   if(length == 0)
      return;
   MPI_Allreduce(source, dest, length, get_mpi_datatype(source), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPImaxArray(const std::vector<T>& source, std::vector<T>& dest, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(source.size() <= dest.size() );
   if(source.size() == 0)
      return;

   MPI_Allreduce(&source[0], &dest[0], source.size(), get_mpi_datatype(&source[0]), MPI_MAX, mpiComm);
}

template <typename T>
inline void PIPS_MPIminArray(const std::vector<T>& source, std::vector<T>& dest, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(source.size() <= dest.size() );
   if(source.size() == 0)
      return;

   MPI_Allreduce(&source[0], &dest[0], source.size(), get_mpi_datatype(&source[0]), MPI_MIN, mpiComm = MPI_COMM_WORLD);
}

template <typename T>
inline void PIPS_MPIgatherv(const T* sendbuf, int sendcnt, T* recvbuf, int* recvcnts, const int* recvoffsets, int root, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(sendcnt >= 0);
   MPI_Gatherv(sendbuf, sendcnt, get_mpi_datatype(sendbuf), recvbuf, recvcnts, recvoffsets, get_mpi_datatype(recvbuf), root, mpiComm = MPI_COMM_WORLD);
}

template <typename T>
inline void PIPS_MPIallgather( const T* sendbuf, int sendcnt, T* recvbuf, int recvcnt, MPI_Comm mpiComm = MPI_COMM_WORLD)
{
   assert(sendcnt >= 0);
   assert(recvcnt >= 0);
   MPI_Allgather(sendbuf, sendcnt, get_mpi_datatype(recvbuf), recvbuf, recvcnt, get_mpi_datatype(recvbuf), mpiComm);
}

#endif
