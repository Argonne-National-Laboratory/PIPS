/*
 * StochPostsolver.h
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_

#include <vector>

#include "QpPostsolver.h"
#include "StochVector.h"
#include "sData.h"
#include "sVars.h"
#include "SystemType.h"
#include "StochRowStorage.h"
#include "StochColumnStorage.h"

class StochPostsolver : public QpPostsolver {
public:

      StochPostsolver( const sData& original_problem );
      virtual ~StochPostsolver();

      void notifyRowModified( const INDEX& row );
      void notifyColModified( const INDEX& col );

      void notifySingletonEqualityRow( int node, int row, BlockType block_type, int col, double coeff, double rhs);
      void notifySingletonIneqalityRow( int node, int row, BlockType block_type, int col, double coeff, double lhs, double rhs );

      void notifySingletonRowBoundsTightened( const INDEX& row, const INDEX& col, double xlow_old, double xupp_old, double xlow_new, double xupp_new, double coeff);
      void notifyRedundantRow( const INDEX& row, int iclow, int icupp, double lhs, double rhs, const StochGenMatrix& matrix_row);
      void notifyFixedColumn( const INDEX& col, double value, double obj_coeff, const StochGenMatrix& eq_mat, const StochGenMatrix& ineq_mat );
      void notifyFixedEmptyColumn( const INDEX& col, double value, double obj_coeff, int ixlow, int ixupp, double lhs, double rhs);

      void putLinkingVarsSyncEvent();
      void notifyFreeColumnSingletonEquality( const INDEX& row, const INDEX& col, double rhs, double obj_coeff, double col_coeff, double xlow, double xupp, const StochGenMatrix& matrix_row );
      void notifyFixedSingletonFromInequalityColumn( const INDEX& col, const INDEX& row, double value, double coeff, double xlow_old, double xupp_old );
      void notifyFreeColumnSingletonInequalityRow( const INDEX& row, const INDEX& col, double rhs, double coeff, double xlow, double xupp, const StochGenMatrix& matrix_row );

      void beginBoundTightening();
      void notifyRowPropagatedBound( const INDEX& row, const INDEX& col, int old_ixlowupp, double old_bound, double new_bound, bool is_upper_bound, const StochGenMatrix& matrix_row);
      void endBoundTightening();
      void notifyDeletedRow( SystemType system_type, int node, int row, bool linking_constraint);
      void notifyParallelColumns();

      void notifyNearlyParallelRowSubstitution( const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double scalar, double translation,
            double obj_col2, double xlow_col2, double xupp_col2, double coeff_col1, double coeff_col2, double parallel_factor );
      void notifyNearlyParallelRowBoundsTightened( const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double xlow_col1,
            double xupp_col1, double xlow_col2, double xupp_col2, double coeff_col1, double coeff_col2, double scalar, double translation, double parallel_factor, double rhs,
            double clow, double cupp);

      void notifyParallelRowsBoundsTightened( const INDEX& row1, const INDEX& row2, double clow_old, double cupp_old, double clow_new, double cupp_new, double factor );

      bool wasColumnRemoved(const INDEX& col) const;
      bool wasRowRemoved(const INDEX& row) const;

private:
      void markColumnRemoved(const INDEX& col);
      void markColumnAdded(const INDEX& col);
      void markRowRemoved(const INDEX& row );
      void markRowAdded(const INDEX& row );

      /// stores row in specified node and returns it's new row index
      int storeRow( const INDEX& row, const StochGenMatrix& matrix_row);
      /// stores col in specified node and returns it's new col index
      int storeColumn( const INDEX& col, const StochGenMatrix& matrix_col_eq, const StochGenMatrix& matrix_col_ineq);

      bool isRowModified(const INDEX& row) const;
      void markRowClean(const INDEX& row);
      void markColClean(const INDEX& col);
      bool isColModified(const INDEX& col) const;

public:
      /// synchronization events

      PostsolveStatus postsolve(const Variables& reduced_solution, Variables& original_solution) const override;
private:

      const int my_rank;
      const bool distributed;

      enum ReductionType
      {
         FIXED_COLUMN = 0,
         SUBSTITUTED_COLUMN = 1,
         BOUND_TIGHTENING_START = 2,
         BOUND_TIGHTENING_END = 3,
         REDUNDANT_ROW = 4,
         BOUNDS_TIGHTENED = 5,
         SINGLETON_EQUALITY_ROW = 6,
         SINGLETON_INEQUALITY_ROW = 7,
         FIXED_EMPTY_COLUMN = 8,
         FREE_COLUMN_SINGLETON_EQUALITY = 9,
         NEARLY_PARALLEL_ROW_SUBSTITUTION = 10,
         DELETED = 11,
         FIXED_COLUMN_SINGLETON_FROM_INEQUALITY = 12,
         FREE_COLUMN_SINGLETON_INEQUALITY_ROW = 13,
         PARALLEL_ROWS_BOUNDS_TIGHTENED = 14,
         NEARLY_PARALLEL_ROW_BOUNDS_TIGHTENED = 15,
         LINKING_VARS_SYNC_EVENT = 16
      };

      const unsigned int n_rows_original;
      const unsigned int n_cols_original;

      /// for now mapping will contain a dummy value for columns that have not been fixed and the value the columns has been fixed to otherwise
      /// 1 indicates that the row / col has not been removed from the problem - -1 indicates the row / col has been removed */
      StochVectorBase<int>* padding_origcol;
      StochVectorBase<int>* padding_origrow_equality;
      StochVectorBase<int>* padding_origrow_inequality;

      /// has a row been modified since last storing it
      /// 1 if yes, -1 if not
      StochVectorBase<int>* eq_row_marked_modified;
      StochVectorBase<int>* ineq_row_marked_modified;
      /// has a column been modified
      StochVectorBase<int>* column_marked_modified;

      /// vectors for storing ints and doubles containting information needed by postsolve
      std::vector<ReductionType> reductions;
      std::vector<INDEX> indices;
      std::vector<unsigned int> start_idx_indices;

      std::vector<double> float_values;
      std::vector<int> int_values;

      std::vector<unsigned int> start_idx_float_values;
      std::vector<unsigned int> start_idx_int_values;

      StochRowStorage row_storage;

      StochColumnStorage col_storage;

      /// stores the index for a row/col indicating where in stored_rows/cols that row/col was stored last
      StochVectorBase<int>* eq_row_stored_last_at;
      StochVectorBase<int>* ineq_row_stored_last_at;
      StochVectorBase<int>* col_stored_last_at;

      /// stores which reduction is last bound-tightening on variable
      StochVectorBase<int>* last_upper_bound_tightened;
      StochVectorBase<int>* last_lower_bound_tightened;

      void finishNotify();

/// postsolve operations
      bool postsolveRedundantRow(sVars& original_vars, int reduction_idx) const;
      bool syncBoundsTightened(sVars& original_vars, int reduction_idx, StochVector& phi, StochVector& gamma) const;
      bool postsolveBoundsTightened(sVars& original_vars, int reduction_idx, StochVector& phi, StochVector& gamma) const;
      bool postsolveFixedColumn(sVars& original_vars, int reduction_idx) const;
      bool postsolveFixedEmptyColumn(sVars& original_vars, int reduction_idx) const;
      bool postsolveFixedColumnSingletonFromInequality(sVars& original_vars, int reduction_idx) const;
      bool postsolveSingletonEqualityRow(sVars& original_vars, int reduction_idx) const;
      bool postsolveSingletonInequalityRow(sVars& original_vars, int reduction_idx) const;
      bool postsolveNearlyParallelRowSubstitution(sVars& original_vars, int reduction_idx) const;
      bool postsolveNearlyParallelRowBoundsTightened(sVars& original_vars, int reduction_idx) const;
      bool postsolveParallelRowsBoundsTightened(sVars& original_vars, int reduction_idx) const;
      bool postsolveFreeColumnSingletonEquality(sVars& original_vars, int reduction_idx) const;
      bool postsolveFreeColumnSingletonInequalityRow(sVars& original_vars, const sData& original_problem, int reduction_idx) const;

      bool syncNewlySetLinkingVars(sVars& original_vars) const;

      void setOriginalVarsFromReduced(const sVars& reduced_vars, sVars& original_vars) const;
      bool allVariablesSet(const sVars& vars) const;
      bool complementarySlackVariablesMet(const sVars& vars, const INDEX& col) const;
      bool complementarySlackRowMet(const sVars& vars, const INDEX& row) const;

      bool sameNonZeroPatternDistributed(const StochVector& svec) const; // TODO: move
      bool sameNonZeroPatternDistributed(const SimpleVector& vec) const; // TODO: move

      template <typename T>
      void setOriginalValuesFromReduced(StochVectorBase<T>& original_vector,
         const StochVectorBase<T>& reduced_vector,
         const StochVectorBase<int>& padding_original) const;

      template <typename T>
      void setOriginalValuesFromReduced(SimpleVectorBase<T>& original_vector,
         const SimpleVectorBase<T>& reduced_vector,
         const SimpleVectorBase<int>& padding_original) const;
};





#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPOSTSOLVER_H_ */
