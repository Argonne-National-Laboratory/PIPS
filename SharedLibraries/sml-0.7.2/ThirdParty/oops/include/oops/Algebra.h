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

#ifndef ALGEBRA_H
#define ALGEBRA_H

#include "oops/Vector.h"
#include "oops/LogiVector.h"
#include "oops/SparseVector.h"
#include "oops/StrSparseVector.h"
#include "oops/ParAlgList.h"

#ifdef WITH_MPI
#include "oops/Delay.h"
#endif

#define MAXNAME   200


typedef enum {NotAType=0, 
	      BlckAng_type, 
	      SparseSimpleMatrix_type, 
	      SparseAugMatrix_type,
	      FakeAlgebra_type,
	      RnkCor_type,
	      DenseMatrix_type,
              BlockDiagSimple_type,
	      BlockDiag_type,
              DblBordDiag_type,
              BlockRow_type,
              BlockCol_type,
              BlockSparse_type,
              BlockDense_type,
              DblBordDiagSimple_type,
              DenseSimpleMatrix_type, 
              DiagonalSimpleMatrix_type, 
              DenseSimpleAugMatrix_type, 
              LastType
} AlgebraType;

typedef enum {

  UnInit=0,
  Terminal=1,
  Intermediate=2

} AlgProc;

enum {done=1, not_done=0};

typedef struct {

  struct Algebra *A;
  struct Algebra *B;
  struct Algebra *Q;

  double *diag;
  int *mapA;
  int *mapQ;

} AugSystemType;

typedef void* AlgMatrix;

typedef void (*AlgProdFunc)
     (struct Algebra*, Vector*, Vector*, int add, double factor);
typedef void (*AlgPrint)
     (FILE*, AlgMatrix, const char*);
typedef void (*AlgComputeAXAt)
     (AlgMatrix, Vector*, AlgMatrix);
typedef void (*AlgZeroCols)
     (struct Algebra*, LogiVector *cols);
typedef void (*AlgSetDiag)
     (struct Algebra*, double fact, Vector *diag);
typedef void (*AlgDelCols)
     (struct Algebra*, int*);
typedef void (*AlgScale)
     (AlgMatrix, double*, double*, double*);
typedef void (*AlgScaleVec)
     (AlgMatrix, double*, Vector*, Vector*);
typedef void (*AlgComputeCholesky)
     (struct Algebra*, FILE*, Vector*, Vector*, double reginit);
typedef void (*AlgSolveCholesky)
     (struct Algebra*, Vector*, Vector*, FILE*);
typedef void (*AlgSolveSparseCholesky)
     (struct Algebra*, StrSparseVector*, Vector*);
typedef void (*AlgSolveTriangular)
     (struct Algebra*, Vector*, Vector*, FILE*);
typedef void (*AlgSolveTriangularPart)
     (struct Algebra*, Vector*, Vector*, int first_last);
typedef void (*AlgSolveTriangularDPart)
     (struct Algebra*, Vector*, Vector*, int first, int last);
typedef void (*AlgSolveTriangEj)
     (AlgMatrix, int, Vector*, FILE*);
typedef void (*AlgInverseTriangL)
     (AlgMatrix, SparseVector**, FILE*);
typedef void (*AlgGetColumn)
     (AlgMatrix, int, Vector*);
typedef void (*AlgGetDiag)
     (struct Algebra*, Vector*);
typedef void (*AlgFreeAlgebra)
     (AlgMatrix);
typedef struct Algebra* (*AlgBuildAAt)
     (struct Algebra*);
typedef int (*AlgCountCol)
     (AlgMatrix, int);
typedef void (*AlgGetSparseColumn) 
     (AlgMatrix, int, int, int*, int*, double*);

typedef void (*AlgGetStrSparseColumn) 
     (AlgMatrix *, int col, StrSparseVector*);

typedef void (*AlgSolveSparseL)
     (struct Algebra*, StrSparseVector*, StrSparseVector*);

typedef void (*AlgSolveSparsePart)
     (struct Algebra*, StrSparseVector*, StrSparseVector*,
      int first, int last);

typedef int (*AlgSetStructure)
     (struct Algebra*, int, int, int, Tree*, Tree*, int);

typedef void (*AlgCopyParInfo)
     (struct Algebra*);
typedef void (*AlgFakeOutAlgebra)
     (struct Algebra*, int parent_is_faked);
typedef void (*AlgParProcess)
     (struct Algebra*);

#ifdef WITH_MPI
typedef void (*AlgAllocateProcs)
     (struct Algebra*, int first, int last, ParAlgList *alglist,
      host_stat_type type, int set_tree);
#else
typedef void (*AlgMakeFlatList)
     (struct Algebra*, ParAlgList *alglist);
#endif

typedef void (*AlgAddAugmentedSystemDiag)
     (struct Algebra*, double, Vector*, Vector*);
typedef void (*AlgMakeAQMap)
     (struct Algebra*, int *, int *);
typedef struct Algebra* (*AlgMakeAugmentedSystem)
     (struct Algebra*, struct Algebra*); 
typedef struct Algebra* (*AlgMakeAugmentedSystemUnsym)
     (struct Algebra*, struct Algebra*, struct Algebra*); 

typedef struct Algebra {

  AlgMatrix                      Matrix;
  char                           name[MAXNAME];
  void                          *id;
  int                            ix;
  int                            nb_row;
  int                            nb_col;
  int                            first_nz_row;
  int                            last_nz_row;
  int                            first_nz_col;
  int                            last_nz_col;
  int                           *PrimalReg;
  Tree                          *Trow;
  Tree                          *Tcol;
  AlgProc                        complexity;
  int                            level;
  AlgebraType                    Type;
  int                            is_augmented;
  AugSystemType                 *AugSystem;

#ifdef WITH_MPI
  int                             first_proc;
  int                             last_proc;
  MPI_Comm                        comm;
#endif

	AlgPrint                       Print;
        AlgProdFunc                    MatrixTimesVect;
        AlgProdFunc                    MatrixTransTimesVect;
        AlgFreeAlgebra                 FreeAlgebra;
        AlgGetColumn                   GetColumn;
        AlgGetColumn                   GetRow;
        AlgGetSparseColumn             GetSparseColumn;
        AlgGetSparseColumn             GetSparseRow;
        AlgGetDiag                     GetCholDiag;
        AlgCountCol                    CountNzCol;
        AlgCountCol                    CountNzRow;
        AlgGetStrSparseColumn          GetStrSparseColumn;
        AlgGetStrSparseColumn          GetStrSparseRow;
        AlgSetStructure                SetStructure;
        AlgZeroCols                    ZeroColumns;
        AlgZeroCols                    ZeroRows;
        AlgSetDiag                     SetDiag;

#ifdef WITH_MPI
        int                            to_be_faked;
        int                            par_split_node;
        AlgCopyParInfo                 CopyParInfo;
#ifndef NOOLD
        AlgFakeOutAlgebra              FakeOutAlgebra;
#endif
#endif
        AlgParProcess                  ParProcessAlgebra;

#ifdef WITH_MPI
        AlgAllocateProcs               AllocateProcessors;
#else
        AlgMakeFlatList               MakeFlatList;
#endif
        AlgDelCols                     DeleteColumns;
        AlgDelCols                     DeleteRows;
        AlgScale                       ColumnScale;
        AlgScale                       RowScale;
        AlgScaleVec                    ColumnScaleVec;
        AlgScaleVec                    RowScaleVec;

#ifndef NOOLD
        AlgBuildAAt                    BuildAAt;
	AlgComputeAXAt                 ComputeAXAt;
#endif
        AlgAddAugmentedSystemDiag      AddAugmentedSystemDiag;
        AlgComputeCholesky             ComputeCholesky;
	AlgSolveCholesky               SolveCholesky;
	AlgSolveSparseCholesky         SolveSparseCholesky;
        AlgSolveTriangular             SolveL;
        AlgSolveTriangularPart         SolveLPart;
	AlgSolveTriangular             SolveD;
	AlgSolveTriangularDPart     SolveDPart;
	AlgSolveTriangular             SolveLt;
	AlgSolveTriangularPart         SolveLtPart;
        AlgSolveSparseL                SolveSparseL;
        AlgSolveSparsePart             SolveSparseD;
        AlgSolveSparsePart             SolveSparseLt;
	AlgSolveTriangEj               SolveLtej;
        AlgInverseTriangL              InverseTriangL;
        AlgMakeAQMap                   MakeAQMap;

#ifndef NOOLD
        AlgMakeAugmentedSystem         MakeAugmentedSystem;
        AlgMakeAugmentedSystemUnsym    MakeAugmentedSystemUnsym;
#endif
        AlgMakeAugmentedSystem         MakeAugmentedSystemNoMem;
        AlgMakeAugmentedSystemUnsym    MakeAugmentedSystemUnsymNoMem;
} Algebra;

int
MatrixTransTimesVectAlg(Algebra* A, Vector *x, Vector *y,
			const int add, const double fact);

int
MatrixTimesVectAlg(Algebra* A, Vector *x, Vector *y,
		   const int add, const double fact);

int
GetColumnAlg(Algebra *A, const int col, Vector *Col);

int
GetStrSparseColumnAlg(Algebra *A, const int col, StrSparseVector *Col);

int
SolveSparseLAlg(Algebra* A, StrSparseVector *x, StrSparseVector *y);

int 
ComputeCholeskyAlg(Algebra* A, FILE *f, Vector *primal_reg, Vector* dual_reg,
		   const double reginit);

int 
SolveLtejAlg(Algebra* A, const int j, Vector *y, FILE *f);

void
FreeAlgebraAlg(Algebra *);

Algebra*
AlgebraBuildAAt(Algebra* A);

void
PrintAlg(FILE *f, Algebra *A, const char *fmt);

void 
PrintMatrixMatlab(FILE *f, Algebra *A, const char *matrixname);

void 
PrintMatrixMatlabSparse(FILE *f, Algebra *A, const char *matrixname);

Algebra *
FakeSubstitute(Algebra *A, int A_id);

Algebra*
SparseAlgTrans(FILE *out, Algebra* A);

Algebra*
FakeAlgebra(Algebra *Source);

Algebra *
FakeParCopy(Algebra *A, int A_id);

Algebra*
NothingAlgebra(void);

Algebra *
NewAlgebra(void);

AugSystemType *
NewAugSystem(void);

void
FreeAugSystem(AugSystemType *Aug);

void
GetStrSparseColumn(Algebra *Any, int col,  StrSparseVector *v);

void
CommonAllocateProcessors(Algebra *AlgAug, int first, int last,
			 ParAlgList *alglist, host_stat_type type,
			 int set_tree);

#ifndef WITH_MPI


Algebra *
InitAlgebras(FILE *out, Algebra* A, Algebra *Q);
#else

MPI_Comm
CreateCommunicator(const int first, const int last);

void
SetParSplitAlgebra(Algebra *A);

Algebra *
InitParAlgebras(FILE *out, Algebra* A, Algebra *Q, 
		struct delay_st *DelayA, struct delay_st *DelayQ);

#endif

Algebra *
InitAlgebrasNew(Algebra* A, Algebra *Q);

#endif /* ALGEBRA_H */
