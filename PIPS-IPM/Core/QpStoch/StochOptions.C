/*
 * StochOptions.C
 *
 *  Created on: 03.04.2020
 *      Author: bzfkempk
 */

#include "StochOptions.h"

#include <limits>
#include "pipsdef.h"

namespace pips_options
{
   StochOptions::StochOptions()
   {
      /* initialize base class options first (QpGenOptions) */
      QpGenOptions::getInstance();

      /* now override with own set of options */
      setDefaults();
   }

   // todo maybe split this up into several submethods?
   void StochOptions::setDefaults()
   {
      /// LINEAR SOLVERS
      bool_options["PARDISO_FOR_GLOBAL_SC"] = true;
      bool_options["PARDISO_SPARSE_RHS_LEAF"] = false;
      /** -1 is choose default */
      int_options["PARDISO_SYMB_INTERVAL"] = -1;
      int_options["PARDISO_PIVOT_PERTURBATION"] = -1;
      int_options["PARDISO_NITERATIVE_REFINS"] = -1;
      int_options["PARDISO_PIVOT_PERTURBATION_ROOT"] = -1;
      int_options["PARDISO_NITERATIVE_REFINS_ROOT"] = -1;

      /// PRECONDITIONERS
      bool_options["PRECONDITION_DISTRIBUTED"] = true;
      bool_options["PRECONDITION_SPARSE"] = true;

      /// INTERIOR-POINT ALGORITHM
      bool_options["IP_ACCURACY_REDUCED"] = false;
      bool_options["IP_PRINT_TIMESTAMP"] = false;
      bool_options["IP_STEPLENGTH_CONSERVATIVE"] = false;

      /// GONDZIO SOLVERS
      bool_options["GONDZIO_ADAPTIVE_LINESEARCH"] = false;

      int_options["GONDZIO_N_LINESEARCH"] = 10;

      /// SOLVER CONTROLS

      int_options["SC_BLOCKWISE_BLOCKSIZE_MAX"] = 64;

      /// ERROR ABSORBTION / ITERATIVE REFINEMENT
      // controls the type of error absorbtion at the outer level of the linear system
      // - 0:no error absortion (OOQP works just fine)
      // - 1:iterative refinement (used when error absortion is
      // also done at a lower level, for example in the solve with
      // the dense Schur complement
      // - 2:BiCGStab with the factorization as preconditioner
      int_options["OUTER_SOLVE"] = 0;
      // controls the type of error absortion/correction at the inner level when solving
      //with the dense Schur complement
      // - 0: no error correction
      // - 1: iter. refin.
      // - 2: BiCGStab
      int_options["INNER_SC_SOLVE"] = 0;

      /// OUTER BIGCSTAB
      double_options["OUTER_BICG_EPSILON"] = 1e-15;

      bool_options["OUTER_BICG_PRINT_STATISTICS"] = false;

      int_options["OUTER_BICG_MAX_ITER"] = 75;
      int_options["OUTER_BICG_MAX_NORMR_DIVERGENCES"] = 4;
      int_options["OUTER_BICG_MAX_STAGNATIONS"] = 4;

      setPresolveDefaults();
   }

   void StochOptions::setPresolveDefaults()
   {
      /** all presolve/postsolve constants and settings */
      // TODO : many of these need adjustments/ have to be thought about
      double_options["PRESOLVE_INFINITY"] = std::numeric_limits<double>::infinity();

      /// STOCH PRESOLVER
      /** limit for max rounds to apply all presolvers */
      int_options["PRESOLVE_MAX_ROUNDS"] = 2;
      /** should the problem be written to std::cout before and after presolve */
      bool_options["PRESOLVE_PRINT_PROBLEM"] = false;
      /** should the presolved problem be written out in MPS format */
      bool_options["PRESOLVE_WRITE_PRESOLVED_PROBLEM_MPS"] = false;
      /** should free variables' bounds be reset after presolve (given the row implying these bounds was not removed */
      bool_options["PRESOLVE_RESET_FREE_VARIABLES"] = false;
      /** verbosity */
      int_options["PRESOLVE_VERBOSITY"] = 1;

      /** turn respective presolvers on/off */
      bool_options["PRESOLVE_BOUND_STRENGTHENING"] = true;
      bool_options["PRESOLVE_PARALLEL_ROWS"] = true;
      bool_options["PRESOLVE_COLUMN_FIXATION"] = true;
      bool_options["PRESOLVE_SINGLETON_ROWS"] = true;
      bool_options["PRESOLVE_SINGLETON_COLUMNS"] = false;

      /// BOUND STRENGTHENING
      /** limit for rounds of bound strengthening per call of presolver */
      int_options["PRESOLVE_BOUND_STR_MAX_ITER"] = 10;
      /** min entry to devide by in order to derive a bound */
      double_options["PRESOLVE_BOUND_STR_NUMERIC_LIMIT_ENTRY"] = 1e-7;
      /** max activity to be devided */
      double_options["PRESOLVE_BOUND_STR_MAX_PARTIAL_ACTIVITY"] = std::numeric_limits<double>::max();
      /** max bounds proposed from bounds strengthening presolver */
      double_options["PRESOLVE_BOUND_STR_NUMERIC_LIMIT_BOUNDS"] = 1e12;

      /// COLUMN FIXATION
      /** limit on the possible impact a column can have on the problem */
      double_options["PRESOLVE_COLUMN_FIXATION_MAX_FIXING_IMPACT"] = 1.0e-12; // for variable fixing

      /// MODEL CLEANUP
      /** limit for the size of a matrix entry below which it will be removed from the problem */
      double_options["PRESOLVE_MODEL_CLEANUP_MIN_MATRIX_ENTRY"] = 1.0e-10;//1.0e-10; // for model cleanup // was 1.0e-10
      /** max for the matrix entry when the impact of entry times (bux-blx) is considered */
      double_options["PRESOLVE_MODEL_CLEANUP_MAX_MATRIX_ENTRY_IMPACT"] = 1.0e-3; // was 1.0e-3
      /** difference in orders between feastol and the impact of entry times (bux-blx) for an entry to get removed */
      double_options["PRESOLVE_MODEL_CLEANUP_MATRIX_ENTRY_IMPACT_FEASDIST"] = 1.0e-2;  // for model cleanup // was 1.0e-2

      /// PARALLEL ROWS
      /** tolerance for comparing two double values in two different rows and for them being considered equal */
      double_options["PRESOLVE_PARALLEL_ROWS_TOL_COMPARE_ENTRIES"] = 1.0e-8;

      /// PRESOLVE DATA
      /** absolute maximum for newly found bounds accepted by presolve */
      double_options["PRESOLVE_MAX_BOUND_ACCEPTED"] = 1e10;

      /** track row/column through presolve */
      bool_options["PRESOLVE_TRACK_ROW"] = false;
      bool_options["PRESOLVE_TRACK_COL"] = false;
      /** the row */
      int_options["PRESOLVE_TRACK_ROW_INDEX"] = 0;
      int_options["PRESOLVE_TRACK_ROW_NODE"] = -1;
      int_options["PRESOLVE_TRACK_ROW_SYSTEM"] = 0; // 0 -> EQ, 1 -> INEQ
      bool_options["PRESOLVE_TRACK_ROW_LINKING"] = false;
      /** the column */
      int_options["PRESOLVE_TRACK_COL_INDEX"] = 0;
      int_options["PRESOLVE_TRACK_COL_NODE"] = -1;

      /// POSTSOLVE
      /** should postsolve be applied */
      bool_options["POSTSOLVE"] = true;

      /** tolerance used for checking residuals after postsolve */
      double_options["POSTSOLVE_TOLERANCE"] = feastol * 1e2;

      /** should the residuals before unscaling, after unscaling before postsolve, after postsolve be printed */
      bool_options["POSTSOLVE_PRINT_RESIDS"] = true;
   }

}
