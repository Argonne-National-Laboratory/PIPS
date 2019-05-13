/*
 * MumpsSolver.h
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_

#include "dmumps_c.h"
#include "mpi.h"

#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "OoqpVectorHandle.h"
#include "pipsport.h"

enum MumpsVerbosity{verb_mute, verb_standard, verb_high};

/** implements linear solver class that uses the MUMPS solver
 */

class MumpsSolver : public DoubleLinearSolver {
 private:
  MumpsSolver() {};

 public:
  MumpsSolver(long long n, MPI_Comm mpiCommPips, MPI_Comm mpiCommMumps);

  ~MumpsSolver();

  void diagonalChanged( int idiag, int extent ) override;
  void matrixChanged() override;
  void solve( OoqpVector& rhs ) override;

  static constexpr MumpsVerbosity defaultVerbosity = verb_standard;
  static constexpr unsigned defaultMaxNiterRefinments = 5;

 private:


  static MUMPS_INT getFortranMPIComm(MPI_Comm mpiComm_c)
  {
     return MUMPS_INT(MPI_Comm_c2f(mpiComm_c));
  };

  void solve(double* vec);

  void setUpMumps();
  void processMumpsResultAnalysis(double starttime, bool verbose);
  void processMumpsResultFactor(double starttime, bool verbose);
  void processMumpsResultSolve(double starttime, bool verbose);

  long long n;
  MumpsVerbosity verbosity;
  unsigned maxNiterRefinments;
  MPI_Comm mpiCommPips, mpiCommMumps;

  int rankMumps, rankPips;

  DMUMPS_STRUC_C* mumps;
  SparseSymMatrix* Msys;


};



#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_ */
