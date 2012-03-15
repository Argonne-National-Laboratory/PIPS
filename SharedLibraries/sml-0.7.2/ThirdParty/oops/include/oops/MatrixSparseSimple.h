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

#ifndef SPARSESIMPLEMATRIX_H
#define SPARSESIMPLEMATRIX_H

#include <stdlib.h>  
#include <stdio.h>  
#include <stdarg.h>  
#include "oops/Vector.h"
#include "oops/Algebra.h"
#include "oops/CallBack.h"

typedef struct {

  int *links;
  int *headers;
  int *indices;
  int *row_len;

} row_st;

typedef struct {

  int nb_row;
  int nb_col;
  int nb_el;
  int max_nb_el;
  double *element;               
  int *row_nbs;
  int *col_beg;
  int *col_len;
  row_st *rows;
  char *name;
  int row_wise;
  int startnb;
  CallBackFunction cbf;

} SparseSimpleMatrix;

#define forall_elt(A, row_pt, col) \
  for (col = 0; col < A->nb_col; col++)\
      for (row_pt = A->col_beg[col];\
           row_pt < A->col_beg[col] + A->col_len[col];\
           row_pt++)

SparseSimpleMatrix*
NewSparseZeroMatrix(const int row_nbs_sz, const int col_nbs_sz);

SparseSimpleMatrix*
NewSparseMatrixNoMem(const int row_nbs_sz, const int col_nbs_sz,
		     CallBackFunction f, const char *name);

SparseSimpleMatrix*
NewSparseMatrix(const int row_nbs_sz, const int col_nbs_sz,
		const int element_sz, const char *name);

Algebra *
NewAlgebraSparse(int nrow, int ncol, const char *name, 
		 CallBackFunction f, void *id);

void
FreeSparseMatrix (SparseSimpleMatrix* N);

void
row_wise (SparseSimpleMatrix * A);

void
free_row_wise (SparseSimpleMatrix * A);

int
resize (SparseSimpleMatrix * A, int new_size);

Algebra *
NewSparseSimpleAlgebra(SparseSimpleMatrix * M);

void
PrtSparseMtxMatlab(FILE *out, SparseSimpleMatrix * A, const char *name);

void 
AddSparseSimpleMatrix(SparseSimpleMatrix *A, SparseSimpleMatrix *B, double f);

#ifdef WITH_MPI
void
SchurSumSparse(SparseSimpleMatrix *M, MPI_Comm comm);
#endif

Algebra *
MatrixSparseClone(SparseSimpleMatrix *M);

#endif /* SPARSESIMPLEMATRIX_H */
