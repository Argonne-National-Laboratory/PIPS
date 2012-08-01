/*
 * COPYRIGHT (c) 2011 Science and Technology Facilities Council (STFC)
 * Original date 1 March 2011
 * All rights reserved
 *
 * Written by: Jonathan Hogg
 *
 * THIS FILE ONLY may be redistributed under the below modified BSD licence.
 * All other files distributed as part of the HSL_MC69 package
 * require a licence to be obtained from STFC and may NOT be redistributed
 * without permission. Please refer to your licence for HSL_M69
 * and conditions. STFC may be contacted via hsl(at)stfc.ac.uk.
 *
 * Modified BSD licence (this header file only):
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of STFC nor the names of its contributors may be used
 *    to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL STFC BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef HSL_MC69Z_H
#define HSL_MC69Z_H

/* Following #defines provide backwards compat to older interface version */
#ifndef mc69_verify
#define mc69_verify mc69_verify_z
#define mc69_print mc69_print_z
#define mc69_set_values mc69_set_values_z
#define mc69_cscl_clean mc69_cscl_clean_z
#define mc69_cscl_convert mc69_cscl_convert_z
#define mc69_cscu_convert mc69_cscu_convert_z
#define mc69_csclu_convert mc69_csclu_convert_z
#define mc69_csrl_convert mc69_csrl_convert_z
#define mc69_csru_convert mc69_csru_convert_z
#define mc69_csrlu_convert mc69_csrlu_convert_z
#define mc69_coord_convert mc69_coord_convert_z
#endif

#include <complex.h>

typedef complex double mc69pkgtype_z_;

#ifndef HSL_MATRIX_TYPE
#define HSL_MATRIX_TYPE
typedef enum hsl_matrix_type {
   /* Undefined or Unknown matrix */
   HSL_MATRIX_UNDEFINED       = 0,

   /* Real matrices */
   HSL_MATRIX_REAL_RECT      =  1, /* real rectangular */
   HSL_MATRIX_REAL_UNSYM     =  2, /* real unsymmetric */
   HSL_MATRIX_REAL_SYM_PSDEF =  3, /* real symmetric pos def */
   HSL_MATRIX_REAL_SYM_INDEF =  4, /* real symmetric indef */
   HSL_MATRIX_REAL_SKEW      =  6, /* real skew symmetric */

   /* Complex matrices */
   HSL_MATRIX_CPLX_RECT      = -1, /* complex rectangular */
   HSL_MATRIX_CPLX_UNSYM     = -2, /* complex unsymmetric */
   HSL_MATRIX_CPLX_HERM_PSDEF= -3, /* hermitian pos def */
   HSL_MATRIX_CPLX_HERM_INDEF= -4, /* hermitian indef */
   HSL_MATRIX_CPLX_SYM       = -5, /* complex symmetric */
   HSL_MATRIX_CPLX_SKEW      = -6  /* complex skew symmetric */
} hsl_matrix_type;
#endif

/* Verify is a matrix is in HSL standard form */
int mc69_verify_z(const int unit, const hsl_matrix_type type, const int findex,
   const int m, const int n, const int ptr[], const int row[],
   const mc69pkgtype_z_ val[], int *more);
/* Pretty print a matrix in HSL standard form */
void mc69_print_z(const int unit, const int lines, const hsl_matrix_type type, 
   const int findex, const int m, const int n, const int ptr[], const int row[],
   const mc69pkgtype_z_ val[]);
/* Set values of A following a conversion */
void mc69_set_values_z(const hsl_matrix_type type, const int lmap,
   const int map[], const mc69pkgtype_z_ val_in[], const int ne,
   mc69pkgtype_z_ val_out[]);
/* Convert a matrix from CSC to HSL standard form, in place */
/* Differences from Fortran: noor, ndup are required. map shold be allocated to
 * have size lmap on entry. Error returned if map turns out too short. */
int mc69_cscl_clean_z(const int unit, const hsl_matrix_type type,
   const int findex, const int m, const int n, int ptr[], int row[],
   mc69pkgtype_z_ val[], int *noor, int *ndup, int *lmap, int *map[]);
/* Convert a matrix from CSC to HSL standard form, out of place */
int mc69_cscl_convert_z(const int unit, const hsl_matrix_type type,
   const int findex, const int m, const int n, const int ptr_in[],
   const int row_in[], const mc69pkgtype_z_ val_in[], int ptr_out[],
   const int lrow, int row_out[], mc69pkgtype_z_ val_out[], int *noor,
   int *ndup, int *lmap, int map[]);
/* Convert a symmetric matrix from upper CSC to HSL standard form */
int mc69_cscu_convert_z(const int unit, const hsl_matrix_type type,
   const int findex, const int n, const int ptr_in[],
   const int row_in[], const mc69pkgtype_z_ val_in[], int ptr_out[],
   const int lrow, int row_out[], mc69pkgtype_z_ val_out[], int *noor,
   int *ndup, int *lmap, int map[]);
/* Convert a symmetric matrix from full CSC to HSL standard form */
int mc69_csclu_convert_z(const int unit, const hsl_matrix_type type,
   const int findex, const int n, const int ptr_in[],
   const int row_in[], const mc69pkgtype_z_ val_in[], int ptr_out[],
   const int lrow, int row_out[], mc69pkgtype_z_ val_out[], int *noor,
   int *ndup, int *lmap, int map[]);
/* Convert a matrix from CSR to HSL standard form */
int mc69_csrl_convert_z(const int unit, const hsl_matrix_type type,
   const int findex, const int m, const int n, const int ptr_in[],
   const int col_in[], const mc69pkgtype_z_ val_in[], int ptr_out[],
   const int lrow, int row_out[], mc69pkgtype_z_ val_out[], int *noor,
   int *ndup, int *lmap, int map[]);
/* Convert a symmetric matrix from upper CSR to HSL standard form */
int mc69_csru_convert_z(const int unit, const hsl_matrix_type type,
   const int findex, const int n, const int ptr_in[],
   const int col_in[], const mc69pkgtype_z_ val_in[], int ptr_out[],
   const int lrow, int row_out[], mc69pkgtype_z_ val_out[], int *noor,
   int *ndup, int *lmap, int map[]);
/* Convert a symmetric matrix from full CSR to HSL standard form */
int mc69_csrlu_convert_z(const int unit, const hsl_matrix_type type,
   const int findex, const int n, const int ptr_in[],
   const int col_in[], const mc69pkgtype_z_ val_in[], int ptr_out[],
   const int lrow, int row_out[], mc69pkgtype_z_ val_out[], int *noor,
   int *ndup, int *lmap, int map[]);
/* Convert a matrix from coordinate to HSL standard form */
int mc69_coord_convert_z(const int unit, const hsl_matrix_type type,
   const int findex, const int m, const int n, const int ne, const int row_in[],
   const int col_in[], const mc69pkgtype_z_ val_in[], int ptr_out[],
   const int lrow, int row_out[], mc69pkgtype_z_ val_out[], int *noor,
   int *ndup, int *lmap, int map[]);

#endif
