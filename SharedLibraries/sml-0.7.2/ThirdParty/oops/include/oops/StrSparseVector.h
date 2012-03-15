/* This file is part of OOPS.
 *
 * OOPS is (c) 2003-2009 Jacek Gondzio and Andreas Grothey, 
 *                       University of Edinburgh
 *
 * OOPS is distributed in a restricted form in the hope that it will be a useful
 * example of what can be done with SML, however it is NOT released under a free
 * software license.
 *
 * You may only redistribute this version of OOPS with a version of SML. You
 * may not link OOPS with code which is not part of SML licensed under the
 * LGPL v3.
 *
 * You may NOT modify, disassemble, or otherwise reverse engineer OOPS.
 *
 * OOPS is distributed WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef STRSPARSEVECTOR_H
#define STRSPARSEVECTOR_H

#include "oops/DenseVector.h"
#include "oops/Vector.h"
#include "oops/SparseVector.h"


typedef struct StrSparseVector {

  Tree *node;
  struct StrSparseVector **subvectors;
  SparseVector *sparse;

} StrSparseVector;

#define OuterProdVector(outpd, NonzCC,CompCol,c)  OuterProdDenseVector (GetDenseVectorFromVector (outpd), NonzCC, CompCol, c)
#define NbOfSubVect(v)  ((v->node)->nb_sons)
#define GetSparseFromStrSparseVector(v)  v->subvectors[(v->node)->index]->sparse
#define SubVector(v,i) v->subvectors[v->node->sons[i]->index]

StrSparseVector *
NewStrSparseVector(Tree *T, const char *name);

void
CopySparseToDenseVector (StrSparseVector *x, DenseVector *y);

void
CopyDenseToSparseVector (DenseVector *x, StrSparseVector *y);

void
CopyStrSparseVectorToVector (StrSparseVector *x, Vector *y);

void
CopyVectorToStrSparseVector(Vector *x, StrSparseVector *y);

void
OuterProdStrSparseVector(int n, StrSparseVector **v, Vector *d, double **c);

void 
daxpyStrSparseVectorToVector(StrSparseVector *x, Vector *y, double a);

void
InfNormStrSparseVector (StrSparseVector *x, double *val, int *j);

void
FreeStrSparseVector (StrSparseVector *V);

void
PrintStrSparseVector (StrSparseVector  *y);

StrSparseVector *
NewCopy_StrSparseVector (StrSparseVector *orig);

void
UnpackStrSparseVector (StrSparseVector *sparse, DenseVector* dense);

void
PackStrSparseVector(StrSparseVector *v);

void
CleanStrSparseVector (StrSparseVector *sparse, DenseVector* dense);

void
ZeroStrSparseVector (StrSparseVector *v);

double
ScalPrStrSparseVector (StrSparseVector *v, DenseVector* dense);

int
StrSparseVectorCountElements(StrSparseVector *x);

double
ddotStrSparseVectorPar (StrSparseVector *x, Vector *y);

SparseVector*
GetSparseVectorFromStrSparseVector (StrSparseVector *x);

void
StrSparseVector2SparseVector(StrSparseVector *x, SparseVector *y);

void
PackSparseVectorFromStrSparseVector (StrSparseVector *V, SparseVector *spV);

#ifdef WITH_MPI
void
SharedSparse(StrSparseVector **sparse, int nb_vect);
#endif

#endif /* STRSPARSEVECTOR_H */
