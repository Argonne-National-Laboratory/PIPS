/*
 * MumpsSolverBase.h
 *
 *  Created on: 25.06.2019
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERBASE_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERBASE_H_



#include "dmumps_c.h"
#include "mpi.h"

#include "DoubleLinearSolver.h"
#include "SparseSymMatrix.h"
#include "OoqpVector.h"
#include "pipsport.h"

#define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation
#define INFOG(I) infog[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

enum MumpsVerbosity{verb_mute, verb_standard, verb_high};

/** implements linear solver class that uses the MUMPS solver
 */

class MumpsSolverBase : public DoubleLinearSolver {

 public:
  MumpsSolverBase( SparseSymMatrix * sgm );
  MumpsSolverBase( MPI_Comm mpiCommPips_c, MPI_Comm mpiCommMumps_c, SparseSymMatrix * sgm );

  ~MumpsSolverBase();

  void diagonalChanged( int idiag, int extent ) override;
  void matrixChanged() override = 0;
  void solve( OoqpVector& rhs ) override;

  static constexpr MumpsVerbosity defaultVerbosity = verb_mute;
//  static constexpr MumpsVerbosity defaultVerbosity = verb_standard;

  static constexpr unsigned defaultMaxNiterRefinments = 5;
  static constexpr int maxNreallocs = 5;

 protected:

  static MUMPS_INT getFortranMPIComm(MPI_Comm mpiComm_c)
  {
     return MUMPS_INT(MPI_Comm_c2f(mpiComm_c));
  };

  void setUpMpiData(MPI_Comm mpiCommPips_c, MPI_Comm mpiCommMumps_c);
  void setUpMumps();
  void processMumpsResultAnalysis(double starttime);
  void processMumpsResultFactor(double starttime);
  void processMumpsResultSolve(double starttime);

  void solve(double* vec);

  long long n;
  MumpsVerbosity verbosity;
  unsigned maxNiterRefinments;
  MPI_Comm mpiCommPips, mpiCommMumps;

  int rankMumps, rankPips;

  DMUMPS_STRUC_C* mumps;
  SparseSymMatrix* Msys;

  // matrix pointer in triplet format
  int* tripletIrn;
  int* tripletJcn;
  double* tripletA;
};



#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERBASE_H_ */
