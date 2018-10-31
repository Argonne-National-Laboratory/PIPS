/* PIPS file - Written by Cosmin Petra, 2018 */

#ifndef MUMPSSOLVER_H
#define MUMPSSOLVER_H

#ifndef WITHOUT_PIPS
#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "SparseSymMatrix.h"
#include "DenseStorageHandle.h"
#include "SparseSymMatrixRowMajList.h"
#else
class DoubleLinearSolver {};
#endif

#include "mpi.h"
#include "dmumps_c.h"

#include <cassert>

/** A linear solver for symmetric indefinite systems using MUMPS
 *
 */
class MumpsSolver : public DoubleLinearSolver {
public:
  MumpsSolver(const long long& globSize, MPI_Comm mumpsMpiComm, MPI_Comm pipsMpiComm);
  virtual ~MumpsSolver();
  /* Set entries that are local to the MPI rank. These arrays are not copied and the
   * the caller of this function needs the keep them for the lifetime of this class. */
  virtual bool setLocalEntries(long long globnnz, long long locnnz, 
			       int* locirn, int* locjcn, double* locA=NULL);
#ifndef WITHOUT_PIPS
  MumpsSolver( SparseSymMatrixRowMajList* storage, MPI_Comm mumpsMpiComm, MPI_Comm pipsMpiComm);
#endif
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual int saveOrderingPermutation();
#ifndef WITHOUT_PIPS
  virtual void solve ( OoqpVector& vec );
  virtual void solve ( GenMatrix& vec );
#else

#endif
  /* solver for vec and return the solution in vec */
  virtual void solve ( double* vec );

  /* Mumps uses Fortran MPI comm, which with some MPI implementations is not compatible
   * with the C MPI comm.
   * For MPI2, and most MPI implementations, you may just do
   *    comm_fortran = (MUMPS_INT) MPI_Comm_c2f(comm_c);
   * (Note that MUMPS INT is defined in [sdcz]mumps c.h and is normally an int.) 
   * For MPI implementations where the Fortran and the C communicators have the same integer representation
   *    comm_fortran = (MUMPS_INT) comm_c;
   * may work.
   * For some MPI implementations, check if 
   *    comm fortran = MPIR FromPointer(comm c) 
   * can be used.
   */
  static MUMPS_INT  getFortranMPIComm(MPI_Comm c_mpiComm)
  {
    return (MUMPS_INT) MPI_Comm_c2f(c_mpiComm);
  };

  /* 0  - no output
   * 1  - error messages only
   * 2  - 1+warning messages
   * 3  - 2+diagnostics and statistics in compact output
   * >3 - 2+diagnostics and statistics in detailed output
   *
   * Method to be called only after MUMPS has been initialized (JOB=-1)
   */
  void setMumpsVerbosity(int level);
protected:
  //mumps data structure
  DMUMPS_STRUC_C* mumps_;
  long long n_;
  int my_mumps_rank_, my_pips_rank_;
  
  MPI_Comm pipsMpiComm, mumpsMpiComm;
#ifndef WITHOUT_PIPS
  SymMatrix * Msys;
#endif
protected:
  void createMumpsStruct();
  void gutsOfconstructor( MPI_Comm mumpsMpiComm_, MPI_Comm pipsMpiComm_ );
  virtual void sysMatToSpTriplet();
#ifndef WITHOUT_PIPS
  MumpsSolver() {};
#endif
}; // end of MumpsSolver class def



class MumpsDenseSolver : public MumpsSolver {
public:
  virtual ~MumpsDenseSolver() {};
#ifndef WITHOUT_PIPS
  MumpsDenseSolver( DenseSymMatrix* storage, MPI_Comm mumpsMpiComm, MPI_Comm pipsMpiComm);
#endif
protected:
  virtual void sysMatToSpTriplet();
}; // end of MumpsSolver class def


#endif
 
