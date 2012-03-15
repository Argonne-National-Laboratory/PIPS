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

#ifndef SPARSEVECTOR_H
#define SPARSEVECTOR_H

#include "oops/DenseVector.h"

typedef struct {

  int dim;
  double *entries;
  int *indices;
  int nb_entries;
  int max_entries;
  char *name;

} SparseVector;

SparseVector *
NewSparseVector(const int dim, const char *name, const int max_nb_el);

void
FreeSparseVector (SparseVector *V);

void 
SparseVectorResize(SparseVector *v, const int new_size);

void
PrintSparseVector(FILE *out, SparseVector *V, const char *format);

void
FromDenseVect (DenseVector *dense, SparseVector *sparse);

void
ToDenseVect (DenseVector *dense, SparseVector *sparse);

void
UnpackSparseVector (SparseVector *sparse, DenseVector* dense);

void
CleanSparseVector (SparseVector *sparse, DenseVector* dense);

void
CopySparseVector (SparseVector *V, SparseVector *CopyV);

double
ScalPrSparseVector (SparseVector *sparse, DenseVector *dense);

void
daxpySparseDense(SparseVector *sparse, DenseVector *dense, const double a);

#endif /* SPARSEVECTOR_H */
