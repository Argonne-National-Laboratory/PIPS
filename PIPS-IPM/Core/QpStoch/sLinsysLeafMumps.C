/*
 * sLinsysLeafMumps.C
 *
 *      Author: bzfrehfe
 */


#include "sLinsysLeafMumps.h"
#include "MumpsSolver.h"

void sLinsysLeafMumps::addTermToSparseSchurCompl( sData *prob,
                      SparseSymMatrix& SC)
{
#if 0
  const SparseGenMatrix& A = prob->getLocalA();
  const SparseGenMatrix& C = prob->getLocalC();
  const SparseGenMatrix& F = prob->getLocalF();
  const SparseGenMatrix& G = prob->getLocalG();
  const SparseGenMatrix& R = prob->getLocalCrossHessian();
#endif
  assert(0);

  // todo
}

void sLinsysLeafMumps::addTermToDenseSchurCompl( sData *prob,
                      DenseSymMatrix& SC)
{
#if 0
  const SparseGenMatrix& A = prob->getLocalA();
  const SparseGenMatrix& C = prob->getLocalC();
  const SparseGenMatrix& F = prob->getLocalF();
  const SparseGenMatrix& G = prob->getLocalG();
  const SparseGenMatrix& R = prob->getLocalCrossHessian();
#endif
  assert(0);

  // todo call one function together with sparse...
}
