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
#include <mpi.h>

const double pips_eps = 1e-13;
const double pips_eps0 = 1e-40;

static const double feastol = 1.0e-6; // was 1.0e-6
static const double infinity = 1.0e30;
// todo rename for more clarity
static const double tolerance1 = 1.0e-3;  // for model cleanup // was 1.0e-3
static const double tolerance2 = 1.0e-2;  // for model cleanup // was 1.0e-2
static const double tol_matrix_entry = 1.0e-10;//1.0e-10; // for model cleanup // was 1.0e-10
static const double tolerance4 = 1.0e-12; // for variable fixing
static const double limit1 = 1.0e3;   // for bound strengthening
static const double limit2 = 1.0e8;   // for bound strengthening
static const int maxIterSR = 10;
static const double tol_compare_double = 1.0e-8;

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

inline void getRankDistributed(MPI_Comm comm, int& myRank, bool& iAmDistrib)
{
   MPI_Comm_rank(comm, &myRank);

   int world_size;
   MPI_Comm_size(comm, &world_size);

   if( world_size > 1)
      iAmDistrib = true;
   else
      iAmDistrib = false;
}

inline void synchronize(int& value)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );
   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

inline void synchronize(double& value)
{
   int myRank;
   bool iAmDistrib;
   getRankDistributed( MPI_COMM_WORLD, myRank, iAmDistrib );
   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void inline abortInfeasible(MPI_Comm comm, std::string message, std::string file, std::string function)
{
   std::cerr << "Infesibility detected in " << file << " function " << function << "!" << std::endl;
   std::cerr << "Message: " << message << std::endl;
   std::cerr << "Aborting now." << std::endl;
   MPI_Abort(comm, 1);
}

#endif
