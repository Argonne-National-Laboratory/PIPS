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

#ifndef DENSEVECTOR_H
#define DENSEVECTOR_H
#include <stdio.h>

#ifndef NAME_FCT
#define NAME_FCT
typedef char*   (*name_fct) (int);
typedef void*   (*mem_fct) (size_t, size_t);
typedef void    (*del_fct) (void *);
#endif

class DenseVector {

 private:

  bool mem_alloc;

 public:

  bool is_zero;
  int dim;
  double *elts;
  char *name;
  DenseVector(const int dim, const char *name, double *mem = NULL);
  ~DenseVector();

};

DenseVector *
NewDenseVector(const int dim, const char *name);

DenseVector *
NewDenseVectorMem(const int dim, const char *name, double* mem);

DenseVector *
SubDenseVector(DenseVector *Big, const int begin, const int dim);

DenseVector*
ReadDenseVector(FILE *f);

void
FreeDenseVector (DenseVector *V);

int 
WriteDenseVector(DenseVector *x, FILE *f);

void
PrintDenseVector(FILE *out, DenseVector *V, const char *format, name_fct name);

void
PrintDense2Vector(FILE *out, DenseVector *V, const char *format, name_fct name,
		  const int begin);

void
PrintDense3Vector(DenseVector *V, const char *format);

void
PrintDense4Vector(DenseVector *V, const char *format, FILE *out);

void
OuterProdDenseVector(DenseVector* outpd, int *NonzCC, double* CompCol, 
                   double **c);

void
daxpyDenseVector(DenseVector* x, DenseVector*y, const double a);

#endif
