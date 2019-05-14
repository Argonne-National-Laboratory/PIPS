/*
 * MumpsSolver.h
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_

#include "dmumps_c.h"
#include "mpi.h"

#include "DoubleLinearSolver.h"
#include "SparseSymMatrix.h"
#include "OoqpVector.h"
#include "pipsport.h"

enum MumpsVerbosity{verb_mute, verb_standard, verb_high};

/** implements linear solver class that uses the MUMPS solver
 */

class MumpsSolver : public DoubleLinearSolver {
 private:
  MumpsSolver() {};

 public:
  MumpsSolver( SparseSymMatrix * sgm );

  ~MumpsSolver();

  void diagonalChanged( int idiag, int extent ) override;
  void matrixChanged() override;
  void solve( OoqpVector& rhs ) override;
  void solve( GenMatrix& rhs ) override;

  static constexpr MumpsVerbosity defaultVerbosity = verb_standard;
  static constexpr unsigned defaultMaxNiterRefinments = 5;

 private:

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


};



#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVER_H_ */
