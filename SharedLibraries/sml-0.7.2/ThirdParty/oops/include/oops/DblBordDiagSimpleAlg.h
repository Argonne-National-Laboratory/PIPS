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

#ifndef DBLBORDDIAGSIMPLEALG_H
#define DBLBORDDIAGSIMPLEALG_H

#include "oops/Algebra.h"

typedef struct {

  Algebra  **D;
  Algebra  **B;
  Algebra  **C;
  int      nb_block;
  int      nb_row;
  int      nb_col;
  int      *col_beg;
  int      *row_beg;
  char     *name;
  
} DblBordDiagSimpleMatrix;

Algebra *
NewAlgebraDblBordDiag(const int nb_block, Algebra **B, Algebra **C,
		      Algebra **D, const char *name);

DblBordDiagSimpleMatrix *
NewDblBordDiagSimpleMatrix(const int nb_block, Algebra **B, Algebra **C,
			   Algebra **D, const char *name);

Algebra *
NewDblBordDiagSimpleAlgebra(DblBordDiagSimpleMatrix *M);

int
DblBordDiagSimpleSetStructure(Algebra *A, int begrow, int begcol,
			      const int level, Tree *Tcol, Tree *Trow,
			      const int copy);

Algebra *
DblBordDiagSimpleMakeAugmentedSystem(Algebra *A, Algebra *Q);

Algebra *
DblBordDiagSimpleMakeAugmentedSystemNoMem(Algebra *A, Algebra *Q);

#endif /* DBLBORDDIAGSIMPLEALG_H */
