/*
 * PardisoIndefSolver.h
 *
 *  Created on: 21.03.2018
 *      Author: Daniel Rehfeldt
 */

#ifndef _PARDISOINDEFSOLVER_H_
#define _PARDISOINDEFSOLVER_H_

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "DenseSymMatrix.h"
#include "SparseSymMatrix.h"
#include "DenseStorageHandle.h"
#include "pipsport.h"

class PardisoIndefSolver : public DoubleLinearSolver
{
   public:
      DenseStorageHandle mStorage;
      SparseStorageHandle mStorageSparse;
   protected:

      constexpr static double precondDiagDomBound = 0.0001;
      constexpr static int pivotPerturbationExpDefault = 8;
      constexpr static int nIterativeRefinsDefault = 8;
      constexpr static bool highAccuracyDefault = true;
      constexpr static bool useSparseRhsDefault = true;
      constexpr static bool parallelForwardBackwardDefault = true;
      constexpr static bool factorizationTwoLevelDefault = true;


      double* x; /* solution vector */

      int mtype;

      int n; /* size of the matrix */

      int nrhs; /* Number of right hand sides. */

      void *pt[64];  /* Internal solver memory pointer pt */

      /* Pardiso control parameters. */
      int iparm[64];
#ifndef WITH_PARDISO_SOLVER
      double dparm[64];
#endif
      int maxfct, mnum, phase, msglvl, solver;
      int* ia;
      int* ja;
      double* a;
      int idum;
      double ddum;
      int pivotPerturbationExp; // 10^-exp
      int nThreads;
      int nIterativeRefins;
      bool highAccuracy;
      bool parallelForwardBackward;
      bool factorizationTwoLevel;

   public:
      PardisoIndefSolver(DenseSymMatrix * storage);
      PardisoIndefSolver(SparseSymMatrix * storage);
      void diagonalChanged(int idiag, int extent) override;
      void matrixChanged() override;
      void matrixRebuild( DoubleMatrix& matrixNew ) override;
      void solve ( OoqpVector& vec ) override;
      void solve ( GenMatrix& vec ) override;
      virtual ~PardisoIndefSolver();

   private:
      void initPardiso();
      void factorizeFromSparse();
      void factorizeFromSparse(SparseSymMatrix& matrix_fortran);
      void factorizeFromDense();
      void factorize();

      void setIparm(int* iparm);
      bool iparmUnchanged();

      bool useSparseRhs;
      bool deleteCSRpointers;
};

#endif /* _PARDISOINDEFSOLVER_H_ */
