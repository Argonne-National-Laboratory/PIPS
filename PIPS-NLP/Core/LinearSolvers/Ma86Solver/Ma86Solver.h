/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin, Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
 
 /* 2015. Modified by Nai-Yuan Chiang for NLP*/

#ifndef Ma86LINSYS_H
#define Ma86LINSYS_H


#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "SparseStorageHandle.h"
#include "OoqpVectorHandle.h"
#include "SparseStorage.h"
#include "SimpleVector.h"



#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif


typedef double ma86pkgtype_d_;
typedef double ma86realtype_d_;



struct ma86_control_d {
   /* Note: 0 is false, non-zero is true */

   /* C/Fortran interface related controls */
   int f_arrays; /* Treat arrays as 1-based (Fortran) if true or 0-based (C) if
                    false. */

   /* Printing controls */
   int diagnostics_level; /* Controls diagnostic printing.*/
               /* Possible values are:
                   < 0: no printing.
                     0: error and warning messages only.
                     1: as 0 plus basic diagnostic printing.
                     2: as 1 plus some more detailed diagnostic messages.
                     3: as 2 plus all entries of user-supplied arrays.       */
   int unit_diagnostics;   /* unit for diagnostic messages
                              Printing is suppressed if unit_diagnostics < 0. */
   int unit_error;         /* unit for error messages
                              Printing is suppressed if unit_error  <  0.     */
   int unit_warning;       /* unit for warning messages
                              Printing is suppressed if unit_warning  <  0.   */

   /* Controls used by ma86_analyse */
   int nemin;  /* Node amalgamation parameter. A child node is merged with its
                  parent if they both involve fewer than nemin eliminations.*/
   int nb;     /* Controls the size of the blocks used within each node (used to
                  set nb within node_type)*/

   /* Controls used by ma86_factor and ma86_factor_solve */
   int action; /* Keep going even if matrix is singular if true, or abort
                  if false */
   int nbi;    /* Inner block size for use with ma64*/
   int pool_size; /* Size of task pool arrays*/
   ma86realtype_d_ small_; /* Pivots less than small are treated as zero*/
   ma86realtype_d_ static_;/* Control static pivoting*/
   ma86realtype_d_ u;      /* Pivot tolerance*/
   ma86realtype_d_ umin;   /* Minimum pivot tolerance*/
   int scaling;            /* Scaling algorithm to use */
};

/***************************************************/

/* data type for returning information to user.*/
struct ma86_info_d {
   ma86realtype_d_ detlog;       /* Holds logarithm of abs det A (or 0) */
   int detsign;         /* Holds sign of determinant (+/-1 or 0) */
   int flag;            /* Error return flag (0 on success) */
   int matrix_rank;     /* Rank of matrix */
   int maxdepth;        /* Maximum depth of the tree. */
   int num_delay;       /* Number of delayed pivots */
   long num_factor;     /* Number of entries in the factor. */
   long num_flops;      /* Number of flops for factor. */
   int num_neg;         /* Number of negative pivots */
   int num_nodes;       /* Number of nodes */
   int num_nothresh;    /* Number of pivots not satisfying u */
   int num_perturbed;   /* Number of perturbed pivots */
   int num_two;         /* Number of 2x2 pivots */
   int pool_size;       /* Maximum size of task pool used */
   int stat;            /* STAT value on error return -1. */
   ma86realtype_d_ usmall;       /* smallest threshold parameter used */
};

//////////////////////////////////////////////////////////////////////////

struct mc68_control {
   /* Extra options for C version */
   int f_array_in;      /* 0 for C array indexing, 1 for Fortran indexing */
   int f_array_out;     /* 0 for C array indexing, 1 for Fortran indexing
                         * NOTE: 2x2 pivot information discarded if C indexing
                         * is used for output! */
   int min_l_workspace; /* Initial size of workspace, as argument in Fortran */
   /* Options from Fortran version */
   int lp;              /* stream number for error messages */
   int wp;              /* stream number for warning messages */
   int mp;              /* stream number for diagnostic messages */
   int nemin;           /* stream number for diagnostic messages */
   int print_level;     /* amount of informational output required */
   int row_full_thresh; /* percentage threshold for full row */
   int row_search;      /* Number of rows searched for pivot with ord=6 */
};

struct mc68_info {
   int flag;            /* error/warning flag */
   int iostat;          /* holds Fortran iostat parameter */
   int stat;            /* holds Fortran stat parameter */
   int out_range;       /* holds number of out of range entries ignored */
   int duplicate;       /* holds number of duplicate entries */
   int n_compressions;  /* holds number of compressions in order */
   int n_zero_eigs;     /* holds the number of zero eigs from ma47 */
   int l_workspace;     /* holds length of workspace iw used in order */
   int zb01_info;       /* holds flag from zb01_expand1 call */
   int n_dense_rows;    /* holds number of dense rows from amdd */
};

//////////////////////////////////////////////////////////////////////////


/** implements the linear solver class using the Ma86 solver
 */
 
class Ma86Solver : public DoubleLinearSolver {
private:
  Ma86Solver() {};
  
 public:

  Ma86Solver( SparseSymMatrix * sgm, const int numOfNegEigVal_in = -1 );
 
  virtual void firstCall();  
  virtual void diagonalChanged( int idiag, int extent );
  virtual int matrixChanged();
  virtual void solve( OoqpVector& rhs );
  virtual void solve( SimpleVector& rhs );
  virtual void solve( GenMatrix& rhs);
  
  virtual ~Ma86Solver();

 // virtual void Lsolve( OoqpVector& x );
 // virtual void Dsolve( OoqpVector& x );
 // virtual void Ltsolve( OoqpVector& x );
  
 private:
  SparseSymMatrix* Msys;
  bool first;
  bool second;
  int n; //

  /** storage for the upper triangular (in row-major format) */
  int     *krowM,    *jcolM;
  double  *M;
  double *val;
  
  /** number of nonzeros in the matrix */
  int      nnz;
  
  double* nvec; //temporary vec
  
  /* Derived types */
  void *keep;
  
  struct ma86_control_d  control;
  struct ma86_info_d  info;
  
  ////////////////////////////////////////////////////////
  struct mc68_control control68;
  struct mc68_info info68;
  ////////////////////////////////////////////////////////

  double * x;
  int * order;
  int * ptr;
  int * row;

	

};

#endif

