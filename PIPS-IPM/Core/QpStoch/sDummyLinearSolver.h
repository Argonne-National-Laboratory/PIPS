#ifndef SDUMMYLINEARSOLVER
#define SDUMMYLINEARSOLVER

#include "DoubleLinearSolver.h"

class sDummyLinearSolver : public DoubleLinearSolver {

  void diagonalChanged( int idiag, int extent ) {};
  void matrixChanged() {};
  void solve ( OoqpVector& x ) {};
  void Lsolve  ( OoqpVector& x ) {};
  void Dsolve  ( OoqpVector& x ) {};
  void Ltsolve ( OoqpVector& x ) {};
};


#endif
