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
#include "SparseSymMatrix.h"
#include "DenseStorageHandle.h"

class PardisoIndefSolver : public DoubleLinearSolver
{
   public:
      DenseStorageHandle mStorage;
   protected:
     // SparseSymMatrix *sparseMat; todo

      int mtype;

      int nrhs; /* Number of right hand sides. */

      void *pt[64];  /* Internal solver memory pointer pt                  */

      /* Pardiso control parameters. */
      int iparm[64];
      double dparm[64];
      int maxfct, mnum, phase, msglvl, solver;
      int* ia;
      int* ja;
      double* a;
      int idum;
      double ddum;

   public:
      PardisoIndefSolver(DenseSymMatrix * storage);
      PardisoIndefSolver(SparseSymMatrix * storage);
      virtual void
      diagonalChanged(int idiag, int extent);
      virtual void matrixChanged();
      virtual void solve ( OoqpVector& vec );
      virtual void solve ( GenMatrix& vec );
      virtual ~PardisoIndefSolver();
};

#endif /* _PARDISOINDEFSOLVER_H_ */
