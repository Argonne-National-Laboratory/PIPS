#ifndef GMSPIPSIO_H_
#define GMSPIPSIO_H_

#include "inttypes.h"
typedef struct 
{
  int32_t     numBlocks;         /**< total number of blocks */
  int32_t     blockID;           /**< block ID */
/* Counts */
  int32_t     n0;                /**< number of variables in block 0 */
  int32_t     ni;                /**< number of variables in block i */
  int32_t     mA;                /**< number of constraint with equality fA */
  int32_t     mC;                /**< number of constraint with inequality fC */
  int32_t     mBL;               /**< number of linking constraint with equality fBL */
  int32_t     mDL;               /**< number of linking constraint with inequality fDL */
  int64_t     nnzA;              /**< number of nonzeros in fA */
  int64_t     nnzB;              /**< number of nonzeros in recourse matrix fB */
  int64_t     nnzC;              /**< number of nonzeros in fC */
  int64_t     nnzD;              /**< number of nonzeros in recourse matrix fD */
  int64_t     nnzBL;             /**< number of nonzeros in fBL */
  int64_t     nnzDL;             /**< number of nonzeros in fDL */
/* Vectors */
  double*     c;                 /**< objective coefficients c_i, length : blockID==0? n0:ni */
  double*     b;                 /**< right hand side of A_i, length: mA */
  double*     clow;              /**< lower bound of C_i, length: mC */
  double*     cupp;              /**< upper bound of C_i, length: mC */
  int16_t*    iclow;             /**< lower bound indicator of C_i, length: mC */
  int16_t*    icupp;             /**< upper bound indicator of C_i, length: mC */
  double*     xlow;              /**< lower bound of x_i, length: blockID==0? n0:ni */
  double*     xupp;              /**< upper bound of x_i, length: blockID==0? n0:ni */
  int16_t*    ixlow;             /**< lower bound indicator of x_i, length: blockID==0? n0:ni */
  int16_t*    ixupp;             /**< upper bound indicator of x_i, length: blockID==0? n0:ni */
  double*     bL;                /**< right hand side of BL, length: mBL */
  double*     dlow;              /**< lower bound of DL, length: mDL */
  double*     dupp;              /**< upper bound of DL, length: mDL */
  int16_t*    idlow;             /**< lower bound indicator of DL, length: mDL */
  int16_t*    idupp;             /**< upper bound indicator of DL, length: mDL */
  int32_t*    permN;             /**< GDX n vector with pointer to 0..(n0+ni-1) space */
  int32_t*    permM;             /**< GDX m vector with pointer to 0..(mA+mC+mBL+mDL-1) space */
/* Matrices */
  int32_t*    rmA;               /**< row major of A, length: mA+1 */
  int32_t*    ciA;               /**< column index of A, length: nnzA */
  double*     valA;              /**< coefficient value of A, length: nnzA */
  int32_t*    rmB;               /**< row major of B, length: mA+1 */
  int32_t*    ciB;               /**< column index of B, length: nnzB */
  double*     valB;              /**< coefficient value of B, length: nnzB */
  int32_t*    rmC;               /**< row major of C, length: mC+1 */
  int32_t*    ciC;               /**< column index of C, length: nnzC */
  double*     valC;              /**< coefficient value of C, length: nnzC */
  int32_t*    rmD;               /**< row major of D, length: mD+1 */
  int32_t*    ciD;               /**< column index of D, length: nnzD */
  double*     valD;              /**< coefficient value of D, length: nnzD */
  int32_t*    rmBL;              /**< row major of BL, length: (blockID==0? n0:ni)+1 */
  int32_t*    ciBL;              /**< column index of BL, length: nnzBL */
  double*     valBL;             /**< coefficient value of BL, length: nnzBL */
  int32_t*    rmDL;              /**< row major of DL, length: (blockID==0? n0:ni)+1 */
  int32_t*    ciDL;              /**< column index of DL, length: nnzDL */
  double*     valDL;             /**< coefficient value of DL, length: nnzDL */
} GMSPIPSBlockData_t;

#if defined(__cplusplus)
extern "C" {
#endif

int initGMSPIPSIO();

int writeBlock(const char* scrFilename,      /** < scratch file name. If NULL write ASCII to stdout */
               GMSPIPSBlockData_t* blk,      /** < block structure to write */
               int printLevel);
int writeSolution(const char* gdxFileStem,  /** < GDX file stem */
                  const int numcol,         /** < length of varl/varmlo/up array */
                  const int numErow,        /** < length of equEm array */
                  const int numIrow,        /** < length of equIl/equIm array */
                  const double objval,      /** < objective value */
                  double* varl,             /** < variable level (can be NULL) */
                  double* varm,             /** < variable marginals (can be NULL) */
                  double* equEl,            /** < equation =e= level (can be NULL) */
                  double* equIl,            /** < equation =lg= level (can be NULL) */
                  double* equEm,            /** < equation =e= marginals */
                  double* equIm,            /** < equation =lg= marginals */
                  const char* GAMSSysDir);  /** < GAMS system directory to locate shared libraries (can be NULL) */                  
int readBlock(const int numBlocks,           /** < total number of blocks n in problem 0..n */
              const int actBlock,            /** < number of block to read 0..n */
              const int strict,              /** < indicator for clean blocks */
              const int offset,              /** < indicator for clean blocks */
              const char* gdxFilename,       /** < GDX file name with CONVERTD jacobian structure */
              const char* GAMSSysDir,        /** < GAMS system directory to locate shared libraries (can be NULL) */
              GMSPIPSBlockData_t* blk);      /** < block structure to be filled */
int gdxSplitting(const int numBlocks,        /** < total number of blocks n in problem 0..n */
              const int actBlock,            /** < block to split from big GDX file, -1 split all */
              const int offset,              /** < indicator for clean blocks */
              const int skipStrings,         /** < indicator for not registering uels and strings */
              const char* gdxFilename,       /** < GDX file name with CONVERTD jacobian structure */
              const char* GAMSSysDir);       /** < GAMS system directory to locate shared libraries (can be NULL) */
void freeBlock(GMSPIPSBlockData_t* blk);
#if defined(__cplusplus)
}
#endif

#endif /* GMSPIPSIO_H_ */
