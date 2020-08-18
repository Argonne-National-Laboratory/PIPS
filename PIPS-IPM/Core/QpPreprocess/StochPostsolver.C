/*
 * StochPostsolver.C
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */

#include "StochPostsolver.h"
#include "OoqpVector.h"
#include "StochOptions.h"
#include "pipsdef.h"
#include "pipsport.h"
#include "StochVectorUtilities.h"

#include <limits>
#include <stdexcept>
#include <iostream>

StochPostsolver::StochPostsolver(const sData& original_problem) :
   QpPostsolver(original_problem),
   my_rank(PIPS_MPIgetRank(MPI_COMM_WORLD)),
   distributed(PIPS_MPIgetDistributed(MPI_COMM_WORLD)),
   postsolve_tol( pips_options::getDoubleParameter("POSTSOLVE_TOLERANCE") ),
   INF_NEG( -pips_options::getDoubleParameter("PRESOLVE_INFINITY") ),
   INF_POS( pips_options::getDoubleParameter("PRESOLVE_INFINITY") ),
   n_rows_original(original_problem.my + original_problem.mz),
   n_cols_original(original_problem.nx),
   padding_origcol( cloneStochVector<double, int>(*original_problem.g) ),
   padding_origrow_equality( cloneStochVector<double, int>(*original_problem.bA) ),
   padding_origrow_inequality( cloneStochVector<double, int>(*original_problem.bu) ),
   eq_row_marked_modified( dynamic_cast<StochVectorBase<int>*>(padding_origrow_equality->clone()) ),
   ineq_row_marked_modified( dynamic_cast<StochVectorBase<int>*>(padding_origrow_inequality->clone()) ),
   column_marked_modified( dynamic_cast<StochVectorBase<int>*>(padding_origcol->clone()) ),
   row_storage( dynamic_cast<const StochGenMatrix&>(*original_problem.A) ),
   col_storage( dynamic_cast<const StochGenMatrix&>(*original_problem.A), dynamic_cast<const StochGenMatrix&>(*original_problem.C) ),
   eq_row_stored_last_at( dynamic_cast<StochVectorBase<int>*>(padding_origrow_equality->clone()) ),
   ineq_row_stored_last_at( dynamic_cast<StochVectorBase<int>*>(padding_origrow_inequality->clone()) ),
   col_stored_last_at( dynamic_cast<StochVectorBase<int>*>(padding_origcol->clone()) ),
   last_upper_bound_tightened( dynamic_cast<StochVectorBase<int>*>(col_stored_last_at->clone()) ),
   last_lower_bound_tightened( dynamic_cast<StochVectorBase<int>*>(col_stored_last_at->clone()) ),
   length_array_outdated_indicators(3),
   array_outdated_indicators(new bool[length_array_outdated_indicators]),
   outdated_linking_vars(array_outdated_indicators[0]),
   outdated_equality_linking_rows(array_outdated_indicators[1]),
   outdated_inequality_linking_rows(array_outdated_indicators[2])
{
   std::memset(array_outdated_indicators, 0, length_array_outdated_indicators * sizeof(bool) );

   const int n_linking_vars = (padding_origcol->vec) ? padding_origcol->vec->n : 0;

   const int n_linking_A = (padding_origrow_equality->vecl) ? padding_origrow_equality->vecl->n : 0;
   const int n_linking_C = (padding_origrow_inequality->vecl) ? padding_origrow_inequality->vecl->n : 0;

   assert( 5.0 * n_linking_vars < std::numeric_limits<int>::max() );

   length_array_linking_var_changes = 5 * n_linking_vars;
   array_linking_var_changes = new double[length_array_linking_var_changes];
   std::memset(array_linking_var_changes, 0, length_array_linking_var_changes * sizeof(double));

   x_changes = SimpleVectorHandle(new SimpleVector(array_linking_var_changes, n_linking_vars));
   v_changes = SimpleVectorHandle(new SimpleVector(array_linking_var_changes + n_linking_vars, n_linking_vars));
   w_changes = SimpleVectorHandle(new SimpleVector(array_linking_var_changes + 2 * n_linking_vars, n_linking_vars));
   gamma_changes = SimpleVectorHandle(new SimpleVector(array_linking_var_changes + 3 * n_linking_vars, n_linking_vars));
   phi_changes = SimpleVectorHandle(new SimpleVector(array_linking_var_changes + 4 * n_linking_vars, n_linking_vars));

   length_array_eq_linking_row_changes = n_linking_A;
   array_eq_linking_row_changes = new double[length_array_eq_linking_row_changes];
   std::memset(array_eq_linking_row_changes, 0, length_array_eq_linking_row_changes * sizeof(double));

   y_changes = SimpleVectorHandle(new SimpleVector(array_eq_linking_row_changes, n_linking_A));

   assert( 4.0 * n_linking_C < std::numeric_limits<int>::max() );

   length_array_ineq_linking_row_changes = 4 * n_linking_C;
   array_ineq_linking_row_changes = new double[length_array_ineq_linking_row_changes];
   std::memset(array_ineq_linking_row_changes, 0, length_array_ineq_linking_row_changes * sizeof(double));

   z_changes = SimpleVectorHandle(new SimpleVector(array_ineq_linking_row_changes, n_linking_C));
   s_changes = SimpleVectorHandle(new SimpleVector(array_ineq_linking_row_changes + n_linking_C, n_linking_C));
   t_changes = SimpleVectorHandle(new SimpleVector(array_ineq_linking_row_changes + 2 * n_linking_C, n_linking_C));
   u_changes = SimpleVectorHandle(new SimpleVector(array_ineq_linking_row_changes + 3 * n_linking_C, n_linking_C));

   padding_origcol->setToConstant(1);
   padding_origrow_equality->setToConstant(1);
   padding_origrow_inequality->setToConstant(1);

   eq_row_marked_modified->setToConstant(1);
   ineq_row_marked_modified->setToConstant(1);
   column_marked_modified->setToConstant(1);

   eq_row_stored_last_at->setToConstant(-1);
   ineq_row_stored_last_at->setToConstant(-1);
   col_stored_last_at->setToConstant(-1);

   last_upper_bound_tightened->setToConstant(-1);
   last_lower_bound_tightened->setToConstant(-1);

   start_idx_float_values.push_back(0);
   start_idx_int_values.push_back(0);
   start_idx_indices.push_back(0);

   /// stuff for synchronization
}

StochPostsolver::~StochPostsolver()
{
   delete last_lower_bound_tightened;
   delete last_upper_bound_tightened;
   delete col_stored_last_at;
   delete ineq_row_stored_last_at;
   delete eq_row_stored_last_at;
   delete column_marked_modified;
   delete ineq_row_marked_modified;
   delete eq_row_marked_modified;
   delete padding_origrow_inequality;
   delete padding_origrow_equality;
   delete padding_origcol;
   delete[] array_outdated_indicators;
   delete[] array_linking_var_changes;
   delete[] array_eq_linking_row_changes;
   delete[] array_ineq_linking_row_changes;
}

void StochPostsolver::notifyRowModified( const INDEX& row )
{
   assert(row.isRow());
   if(row.getSystemType() == EQUALITY_SYSTEM)
      getSimpleVecFromRowStochVec(*eq_row_marked_modified, row) = 1;
   else
      getSimpleVecFromRowStochVec(*ineq_row_marked_modified, row) = 1;
}

void StochPostsolver::notifyColModified( const INDEX& col )
{
   assert(col.isCol());
   getSimpleVecFromColStochVec(*column_marked_modified, col) = 1;
}

bool StochPostsolver::isRowModified( const INDEX& row) const
{
   assert(row.isRow());
   if(row.getSystemType() == EQUALITY_SYSTEM)
      return getSimpleVecFromRowStochVec(*eq_row_marked_modified, row) == 1;
   else
      return getSimpleVecFromRowStochVec(*ineq_row_marked_modified, row) == 1;
}

void StochPostsolver::markRowClean( const INDEX& row )
{
   assert(row.isRow());
   if(row.getSystemType() == EQUALITY_SYSTEM)
      getSimpleVecFromRowStochVec(*eq_row_marked_modified, row) = -1;
   else
      getSimpleVecFromRowStochVec(*ineq_row_marked_modified, row) = -1;
}

void StochPostsolver::markColClean( const INDEX& col )
{
   assert(col.isCol());
   getSimpleVecFromColStochVec(*col_stored_last_at, col) = -1;
}

bool StochPostsolver::isColModified( const INDEX& col) const
{
   assert(col.isCol());
   return getSimpleVecFromColStochVec(*column_marked_modified, col) == 1;
}

int StochPostsolver::storeRow( const INDEX& row, const StochGenMatrix& matrix_row)
{
   assert(row.isRow());

   if( isRowModified(row) )
   {
      markRowClean(row);

      const int stored_at = row_storage.storeRow(row, matrix_row);
      if( row.inEqSys() )
         getSimpleVecFromRowStochVec(*eq_row_stored_last_at, row) = stored_at;
      else
         getSimpleVecFromRowStochVec(*ineq_row_stored_last_at, row) = stored_at;
      return stored_at;
   }
   else
   {
      if( row.inEqSys() )
      {
         assert(getSimpleVecFromRowStochVec(*eq_row_stored_last_at, row) != -1);
         return getSimpleVecFromRowStochVec(*eq_row_stored_last_at, row);
      }
      else
      {
         assert(getSimpleVecFromRowStochVec(*ineq_row_stored_last_at, row) != -1);
         return getSimpleVecFromRowStochVec(*ineq_row_stored_last_at, row);
      }
   }
}

int StochPostsolver::storeColumn( const INDEX& col, const StochGenMatrix& matrix_col_eq, const StochGenMatrix& matrix_col_ineq)
{
   assert(col.isCol());

   if( isColModified(col) )
   {
      markColClean(col);
      const int stored_at = col_storage.storeCol(col, matrix_col_eq, matrix_col_ineq);

      getSimpleVecFromColStochVec(*col_stored_last_at, col) = stored_at;

      return stored_at;
   }
   else
   {
      assert(getSimpleVecFromColStochVec(*col_stored_last_at, col) != -1);
      return getSimpleVecFromColStochVec(*col_stored_last_at, col);
   }
}

void StochPostsolver::notifyFixedSingletonFromInequalityColumn( const INDEX& col, const INDEX& row, double value, double coeff, double xlow_old, double xupp_old )
{
   assert(col.isCol());
   assert(!wasColumnRemoved(col));

   reductions.push_back(FIXED_COLUMN_SINGLETON_FROM_INEQUALITY);

   indices.push_back(col);
   indices.push_back(row);

   float_values.push_back(value);
   float_values.push_back(coeff);
   float_values.push_back(xlow_old);
   float_values.push_back(xupp_old);

   finishNotify();
}

void StochPostsolver::notifyFreeColumnSingletonInequalityRow( const INDEX& row, const INDEX& col, double rhs, double coeff, double xlow, double xupp, const StochGenMatrix& matrix_row )
{
   assert(row.isRow());
   assert( col.isCol() || col.isEmpty() );
   if( col.isEmpty() )
      assert( row.isLinkingRow() );

   markRowRemoved(row);

   reductions.push_back(FREE_COLUMN_SINGLETON_INEQUALITY_ROW);

   indices.push_back(row);
   indices.push_back(col);

   const int index_stored_row = row_storage.storeRow(row, matrix_row);

   int_values.push_back(index_stored_row);

   float_values.push_back(rhs);
   float_values.push_back(coeff);
   float_values.push_back(xlow);
   float_values.push_back(xupp);

   finishNotify();
}

void StochPostsolver::putBoundTighteningLinkingRowSyncEvent()
{
   reductions.push_back(BOUND_TIGHTENING_LINKING_ROW_SYNC_EVENT);

   finishNotify();
}

void StochPostsolver::putLinkingVarsSyncEvent()
{
   reductions.push_back(LINKING_VARS_SYNC_EVENT);
   finishNotify();
}

void StochPostsolver::putLinkingRowIneqSyncEvent()
{
   reductions.push_back(LINKING_INEQ_ROW_SYNC_EVENT);
   finishNotify();
}

void StochPostsolver::putLinkingRowEqSyncEvent()
{
   reductions.push_back(LINKIN_EQ_ROW_SYNC_EVENT);
   finishNotify();
}

void StochPostsolver::notifyFreeColumnSingletonEquality( const INDEX& row, const INDEX& col, double rhs, double obj_coeff, double col_coeff, double xlow, double xupp, const StochGenMatrix& matrix_row )
{
   assert(row.isRow());
   assert(col.isCol() || col.isEmpty() );

   if( col.isCol() )
      assert(!wasColumnRemoved(col));

   assert(!wasRowRemoved(row));

   markRowRemoved(row);
   
   const int stored_row_idx = storeRow(row, matrix_row);

   reductions.push_back( FREE_COLUMN_SINGLETON_EQUALITY );

   indices.push_back(col);
   indices.push_back(row);

   float_values.push_back(rhs);
   float_values.push_back(obj_coeff);
   float_values.push_back(col_coeff);
   float_values.push_back(xlow);
   float_values.push_back(xupp);
   int_values.push_back(stored_row_idx);

   finishNotify();
}

/** substitute col2 with scalar * col1 + translation
 *
 * this happens for:
 *    two equality rows - coeff_col2 != 0
 *    two inequality rows (lot's of prerequisites)
 *
 */
void StochPostsolver::notifyNearlyParallelRowSubstitution(const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double scalar, double translation,
   double obj_col1, double obj_col2, double xlow_col2, double xupp_col2, double coeff_col1, double coeff_col2, double parallel_factor )
{
   assert( row1.isRow() );
   assert( row2.isRow() );
   assert( row1.getNode() == row2.getNode() );

   /* both rows must be in the same system */
   assert( ( row1.inEqSys() && row2.inEqSys() ) ||
      ( row1.inInEqSys() && row2.inInEqSys() ) );

   assert( col2.isCol() );

   if( row1.inInEqSys() )
      assert( col1.isCol() );

   assert( !PIPSisZero(coeff_col2) );

   if( row2.inInEqSys() )
   {
      assert( !PIPSisZero(coeff_col1) );
      assert( PIPSisLT(0.0, coeff_col1 * coeff_col2) );
      assert( PIPSisZero(translation) );
   }

   // todo : linking rows are not possible yet
   assert( !row1.isLinkingRow() && !row2.isLinkingRow() );

   if( col1.isCol() )
      assert( !wasColumnRemoved(col1) );
   else
      assert( PIPSisZero(coeff_col1) );

   assert( !wasColumnRemoved(col2) );

   assert( !wasRowRemoved(row1) );
   assert( !wasRowRemoved(row2) );

   reductions.push_back(NEARLY_PARALLEL_ROW_SUBSTITUTION);

   indices.push_back( col1 );
   indices.push_back( col2 );
   indices.push_back( row1 );
   indices.push_back( row2 );

   float_values.push_back( scalar );
   float_values.push_back( translation );
   float_values.push_back( obj_col1 );
   float_values.push_back( obj_col2 );

   float_values.push_back( xlow_col2 );
   float_values.push_back( xupp_col2 );

   float_values.push_back( coeff_col1 );
   float_values.push_back( coeff_col2 );

   float_values.push_back( parallel_factor );
   finishNotify();
}

/** bounds on col1 get tightened by col2 with col1 = parallelity_factor * col2 */
void StochPostsolver::notifyNearlyParallelRowBoundsTightened( const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2,
      double xlow_col1, double xupp_col1, double xlow_col2, double xupp_col2, double coeff_col1, double coeff_col2, double scalar, double translation, double parallel_factor,
      double rhs, double clow, double cupp)
{
   assert( row1.isRow() );
   assert( row2.isRow() );
   assert( row1.inEqSys() );
   assert( col1.isCol() );

   if( row2.inEqSys() )
      assert( col2.isCol() );
   else
      assert( col2.isEmpty() );

   assert( !wasRowRemoved(row1) );
   assert( !wasRowRemoved(row1) );
   assert( !wasColumnRemoved(col1) );

   if( col2.isCol() )
      assert( !wasColumnRemoved(col2) );
   else
   {
      assert( PIPSisZero(coeff_col2) );
      assert( xlow_col2 == INF_NEG );
      assert( xupp_col2 == INF_POS );
   }

   assert( !PIPSisZero(scalar) );

   reductions.push_back(NEARLY_PARALLEL_ROW_BOUNDS_TIGHTENED);

   indices.push_back(row1);
   indices.push_back(row2);
   indices.push_back(col1);
   indices.push_back(col2);

   float_values.push_back(xlow_col1);
   float_values.push_back(xupp_col1);
   float_values.push_back(xlow_col2);
   float_values.push_back(xupp_col2);

   float_values.push_back(coeff_col1);
   float_values.push_back(coeff_col2);

   float_values.push_back(scalar);
   float_values.push_back(translation);
   float_values.push_back(parallel_factor);

   float_values.push_back(rhs);
   float_values.push_back(clow);
   float_values.push_back(cupp);

   finishNotify();
}

/** tighten row1 bounds with row2 bounds - row1 = s * row2 */
void StochPostsolver::notifyParallelRowsBoundsTightened( const INDEX& row1, const INDEX& row2, double clow_old, double cupp_old, double clow_new, double cupp_new, double factor )
{
   assert(row1.isRow());
   assert(row2.isRow());
   assert(row1.inInEqSys());
   assert(row2.inInEqSys());

   assert(!wasRowRemoved(row1));
   assert(!wasRowRemoved(row2));

   reductions.push_back(PARALLEL_ROWS_BOUNDS_TIGHTENED);

   indices.push_back(row1);
   indices.push_back(row2);

   float_values.push_back(clow_old);
   float_values.push_back(cupp_old);
   float_values.push_back(clow_new);
   float_values.push_back(cupp_new);
   float_values.push_back(factor);

   finishNotify();
}


/** bounds got tightened by propagating a singleton row - not necessary to store whole row */
void StochPostsolver::notifySingletonRowBoundsTightened( const INDEX& row, const INDEX& col, double xlow_old, double xupp_old, double xlow_new, double xupp_new, double coeff )
{
   assert( row.isRow() || row.isEmpty() );
   assert( col.isCol() );

   if( row.isEmpty() )
      assert( col.isLinkingCol() );

   assert( PIPSisLE(xlow_new, xupp_new) );
   assert( xupp_new != INF_POS || xlow_new != INF_NEG );
   assert( !PIPSisZero(coeff) || coeff == NAN );

   if( coeff == NAN )
   {
      assert( col.isCol() );
      assert( row.isEmpty() );
   }
   if( row.isRow() )
      assert( !wasRowRemoved(row) );
   assert( !wasColumnRemoved(col) );

   ReductionType red = row.inEqSys() ? SINGLETON_EQUALITY_ROW : SINGLETON_INEQUALITY_ROW;

   if( red == SINGLETON_EQUALITY_ROW )
      assert(xupp_new == xlow_new);

   reductions.push_back(red);

   indices.push_back(row);
   indices.push_back(col);

   float_values.push_back(xlow_old);
   float_values.push_back(xupp_old);
   float_values.push_back(xlow_new);
   float_values.push_back(xupp_new);
   float_values.push_back(coeff);

   finishNotify();
}

/** postsolve has to compute the optimal dual multipliers here and set the primal value accordingly */
void StochPostsolver::notifyFixedColumn( const INDEX& col, double value, double obj_coeff, const StochGenMatrix& eq_mat, const StochGenMatrix& ineq_mat)
{
   assert(col.isCol());
   assert( !wasColumnRemoved(col) );
   markColumnRemoved(col);
   if( col.isLinkingCol() )
      assert( PIPS_MPIisValueEqual(col.getIndex(), MPI_COMM_WORLD) );

   /* store current upper and lower bounds of x and the local column */
   const int col_index = col_storage.storeCol(col, eq_mat, ineq_mat);

   reductions.push_back(FIXED_COLUMN);

   indices.push_back(col);

   int_values.push_back(col_index);

   float_values.push_back(value);
   float_values.push_back(obj_coeff);

   finishNotify();
}

void StochPostsolver::notifyFixedEmptyColumn( const INDEX& col, double value, double obj_coeff, double xlow, double xupp)
{
   assert(col.isCol());
   assert( !wasColumnRemoved(col) );
   if( col.isLinkingCol() )
      assert( PIPS_MPIisValueEqual(col.getNode(), MPI_COMM_WORLD) );

   markColumnRemoved(col);

   assert(PIPSisLEFeas(xlow, value));
   assert(PIPSisLEFeas(value, xupp));

   reductions.push_back(FIXED_EMPTY_COLUMN);
   indices.push_back(col);
   float_values.push_back(value);
   float_values.push_back(obj_coeff);
   float_values.push_back(xlow);
   float_values.push_back(xupp);

   finishNotify();
}

void StochPostsolver::notifyRedundantRow( const INDEX& row, int iclow, int icupp, double lhs, double rhs, const StochGenMatrix& matrix_row )
{
   assert(row.isRow());
   assert(iclow == 1 || iclow == 0);
   assert(icupp == 1 || icupp == 0);

   if(row.getSystemType() == INEQUALITY_SYSTEM)
      assert(iclow + icupp > 0);

   assert( !wasRowRemoved(row) );
   markRowRemoved(row);

   if( row.isLinkingRow() )
      assert(PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD) );

   /* save row for postsolve */
   reductions.push_back(REDUNDANT_ROW);
   indices.push_back(row);

   int index_stored_row = storeRow(row, matrix_row);

   float_values.push_back(lhs);
   float_values.push_back(rhs);
   int_values.push_back(index_stored_row);
   int_values.push_back(iclow);
   int_values.push_back(icupp);

   finishNotify();
}

void StochPostsolver::beginBoundTightening()
{
   putBoundTighteningLinkingRowSyncEvent();
   putLinkingVarsSyncEvent();
}

void StochPostsolver::endBoundTightening( const std::vector<int>& store_linking_rows_A, const std::vector<int>& store_linking_rows_C,
      const StochGenMatrix& mat_A, const StochGenMatrix& mat_C )
{
   reductions.push_back(STORE_BOUND_TIGHTENING_LINKING_ROWS);
   for( unsigned int i = 0; i < store_linking_rows_A.size(); ++i)
   {
      if( store_linking_rows_A[i] != 0 )
      {
         const INDEX row(ROW, -1, i, true, EQUALITY_SYSTEM);
         const int index = row_storage.storeRow(row, mat_A);

         indices.push_back(row);
         int_values.push_back(index);
      }
   }
   for( unsigned int i = 0; i < store_linking_rows_C.size(); ++i)
   {
      if( store_linking_rows_C[i] != 0 )
      {
         const INDEX row(ROW, -1, i, true, INEQUALITY_SYSTEM);
         const int index = row_storage.storeRow(row, mat_C);

         indices.push_back(row);
         int_values.push_back(index);
      }
   }

   finishNotify();
}

void StochPostsolver::notifyRowPropagatedBound( const INDEX& row, const INDEX& col, double old_bound, double new_bound, bool is_upper_bound, const StochGenMatrix& matrix_row)
{
   assert(!PIPSisEQ(old_bound, new_bound));
   assert(row.isRow());
   assert(col.isCol());

   if(is_upper_bound)
      assert( PIPSisLT(new_bound, old_bound) );
   else
      assert( PIPSisLT(old_bound, new_bound) );

   // TODO : check!
   int& index_last = is_upper_bound ? getSimpleVecFromColStochVec(*last_upper_bound_tightened, col) : getSimpleVecFromColStochVec(*last_lower_bound_tightened, col);
   if( index_last != -1 )
   {
      reductions.at(index_last) = DELETED;
      /* copy old old_bounds */
      assert(int_values.at(start_idx_int_values.at(index_last)) == is_upper_bound);
      old_bound = float_values.at(start_idx_float_values.at(index_last));
      if( is_upper_bound )
         assert( PIPSisLT( new_bound, float_values.at(start_idx_float_values.at(index_last) + 1) ) );
      else
         assert( PIPSisLT( float_values.at(start_idx_float_values.at(index_last) + 1), new_bound ) );
   }
   index_last = reductions.size();

   reductions.push_back( BOUNDS_TIGHTENED );

   indices.push_back( row );
   indices.push_back( col );

   int index_stored_row = storeRow( row, matrix_row );

   int_values.push_back(is_upper_bound);
   int_values.push_back(index_stored_row);

   float_values.push_back(old_bound);
   float_values.push_back(new_bound);

   finishNotify();
}

void StochPostsolver::notifyDeletedRow( SystemType system_type, int node, int row, bool linking_constraint)
{
   throw std::runtime_error("Not yet implemented");
}

void StochPostsolver::notifyParallelColumns()
{
   throw std::runtime_error("Not yet implemented");
}

void StochPostsolver::finishNotify()
{
   assert(reductions.size() == start_idx_float_values.size());
   assert(reductions.size() == start_idx_int_values.size());
   assert(reductions.size() == start_idx_indices.size());

   start_idx_int_values.push_back(int_values.size());
   start_idx_float_values.push_back(float_values.size());
   start_idx_indices.push_back(indices.size());
}

bool StochPostsolver::wasColumnRemoved(const INDEX& col) const
{
   assert(col.isCol());
   return getSimpleVecFromColStochVec(*padding_origcol, col) == -1;
}

void StochPostsolver::markColumnRemoved(const INDEX& col )
{
   assert(col.isCol());
   assert(!wasColumnRemoved(col));
   getSimpleVecFromColStochVec(*padding_origcol, col) = -1;
}

void StochPostsolver::markColumnAdded(const INDEX& col)
{
   assert( col.isCol() );
   assert( wasColumnRemoved(col) );
   getSimpleVecFromColStochVec(*padding_origcol, col) = 1;
}

bool StochPostsolver::wasRowRemoved( const INDEX& row ) const
{
   assert(row.isRow());
   if( row.inEqSys() )
      return getSimpleVecFromRowStochVec(*padding_origrow_equality, row) == -1;
   else
      return getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) == -1;
}

void StochPostsolver::markRowRemoved(const INDEX& row)
{
   assert(row.isRow());
   assert(!wasRowRemoved(row));
   if( row.inEqSys() )
      getSimpleVecFromRowStochVec(*padding_origrow_equality, row) = -1;
   else
      getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) = -1;
}

// todo use somehow?
void StochPostsolver::markRowAdded(const INDEX& row)
{
   assert( row.isRow() );
   assert( wasRowRemoved(row) );
   if( row.inEqSys() )
      getSimpleVecFromRowStochVec(*padding_origrow_equality, row) = 1;
   else
      getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) = 1;
}

// todo : usage and check of padding origrow - can already be done - even without any dual postsolve stuff
// todo : sort reductions by nodes ? and then reverse ?
PostsolveStatus StochPostsolver::postsolve(const Variables& reduced_solution, Variables& original_solution)
{
   if(my_rank == 0)
      std::cout << "start postsolving... " << std::endl;

   const sVars& stoch_reduced_sol = dynamic_cast<const sVars&>(reduced_solution);
   sVars& stoch_original_sol = dynamic_cast<sVars&>(original_solution);

   /* original variables are now reduced vars padded with zeros */
   setOriginalVarsFromReduced(stoch_reduced_sol, stoch_original_sol);

   bool postsolve_success = true;
   /* post-solve the reductions in reverse order */
   /* shift the postsolve of bound tightenings to the very end since they are numerically unstable */
   for( int i = reductions.size() - 1; i >= 0; --i )
   {
      const ReductionType type = reductions.at(i);

      switch( type )
      {
      case DELETED:
      {
         break;
      }
      case REDUNDANT_ROW:
      {
         const bool success = postsolveRedundantRow(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case BOUNDS_TIGHTENED:
      {
         const bool success = postsolveBoundsTightened(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case FIXED_COLUMN:
      {
         const bool success = postsolveFixedColumn(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case FIXED_EMPTY_COLUMN:
      {
         const bool success = postsolveFixedEmptyColumn(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case FIXED_COLUMN_SINGLETON_FROM_INEQUALITY:
      {
         const bool success = postsolveFixedColumnSingletonFromInequality(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case SINGLETON_EQUALITY_ROW:
      {
         const bool success = postsolveSingletonEqualityRow(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case SINGLETON_INEQUALITY_ROW:
      {
         const bool success = postsolveSingletonInequalityRow(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case FREE_COLUMN_SINGLETON_EQUALITY:
      {
         const bool success = postsolveFreeColumnSingletonEquality(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case NEARLY_PARALLEL_ROW_SUBSTITUTION:
      {
         const bool success = postsolveNearlyParallelRowSubstitution(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case NEARLY_PARALLEL_ROW_BOUNDS_TIGHTENED:
      {
         const bool success = postsolveNearlyParallelRowBoundsTightened(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case FREE_COLUMN_SINGLETON_INEQUALITY_ROW:
      {
         const bool success = postsolveFreeColumnSingletonInequalityRow(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case PARALLEL_ROWS_BOUNDS_TIGHTENED:
      {
         const bool success = postsolveParallelRowsBoundsTightened(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case LINKING_VARS_SYNC_EVENT:
      {
         const bool success = syncLinkingVarChanges(stoch_original_sol);
         postsolve_success = postsolve_success && success;
         break;
      }
      case LINKING_INEQ_ROW_SYNC_EVENT:
      {
         const bool success = syncIneqLinkingRowChanges(stoch_original_sol);
         postsolve_success = postsolve_success && success;
         break;
      }
      case LINKIN_EQ_ROW_SYNC_EVENT:
      {
         const bool success = syncEqLinkingRowChanges(stoch_original_sol);
         postsolve_success = postsolve_success && success;
         break;
      }
      case BOUND_TIGHTENING_LINKING_ROW_SYNC_EVENT:
      {
         const bool success = syncLinkingRowsAfterBoundTightening(stoch_original_sol, i);
         postsolve_success = postsolve_success && success;
         break;
      }
      case STORE_BOUND_TIGHTENING_LINKING_ROWS:
      {
         break;
      }
      default:
      {
         throw std::runtime_error("Tried to postsolve not supported reduction type");
         break;
      }
      }
   }

   /* compute all s, t and u that have not yet been computed */

   /* assert that all variables have been set */
   assert( stoch_original_sol.isRootNodeInSync() );
   assert( allVariablesSet(stoch_original_sol) );
   if( my_rank == 0 )
      std::cout << "finished postsolving... " << std::endl;

   if( postsolve_success )
      return PRESOLVE_OK;
   else
      return PRESOLVE_FAIL;
}

/**
 * postsolve for a redundant row is to set all dual variables to zero - the row itself has no primal impact
 * for linking rows:
 *    all processes should have removed them in the same order
 *       -> assert in place to check this
 *    the current activity gets synchronized for slack computation
 */
bool StochPostsolver::postsolveRedundantRow(sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == REDUNDANT_ROW );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 1 == next_first_index);
   assert(first_float_val + 2 == next_first_float_val);
   assert(first_int_val + 3 == next_first_int_val);
#endif

   /* get stored data for postsolve */
   const INDEX& row = indices.at(first_index);
   assert( row.isRow() );
   if( row.isLinkingRow() )
      assert(PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD));

   const double lhs = float_values.at(first_float_val);
   const double rhs = float_values.at(first_float_val + 1);

   const int index_stored_row = int_values.at(first_int_val);
   const int iclow = int_values.at(first_int_val + 1);
   const int icupp = int_values.at(first_int_val + 2);
   assert(iclow + icupp >= 1);

   const INDEX stored_row(ROW, row.getNode(), index_stored_row, row.getLinking(), EQUALITY_SYSTEM);

   // TODO: this could be optimized - all linking rows will be after each other - we could wait with allreduce --- actually not sure whether this would ever happen or not - multiple redundant linking rows...
   // until the current boost of linking rows is processed and then allreduce all activities at once and unpate the slacks then only
   /* get current row activity - redundant linking rows have to lie on the stack in the same order */

   /* primal values */
   double value_row = 0.0;
   StochVector& x_vec = dynamic_cast<StochVector&>(*original_vars.x);

   if( row.isLinkingRow() )
   {
      assert(PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD));

      if(my_rank == 0)
         value_row = row_storage.multRowTimesVec( stored_row, x_vec );
      else
         value_row = row_storage.multLinkingRowTimesVecWithoutBl0(index_stored_row, x_vec);
      /* this might get very expensive if there is many redundant linking rows */
      PIPS_MPIgetSumInPlace(value_row, MPI_COMM_WORLD);
   }
   else
      value_row = row_storage.multRowTimesVec( stored_row, x_vec );

   assert( wasRowRemoved(row) );
   markRowAdded(row);
   if( row.inEqSys() )
   {
      assert(PIPSisEQ(lhs, rhs));
      if( !PIPSisEQFeas(value_row, rhs) )
         PIPSdebugMessage("Postsolve Warning: when reintroducing a redundant equality row it did not meet its rhs with feastol: %f != %f", value_row, rhs);
      assert( PIPSisEQ(value_row, rhs, postsolve_tol) );

      /* set dual multiplier to zero and mark row as added */
      getSimpleVecFromRowStochVec(original_vars.y, row) = 0;
   }
   else
   {
      /* set dual multipliers to zero and mark row as added */
      /* dual of row is zero */
      getSimpleVecFromRowStochVec(original_vars.z, row) = 0;
      getSimpleVecFromRowStochVec(original_vars.lambda, row) = 0;
      getSimpleVecFromRowStochVec(original_vars.pi, row) = 0;

      getSimpleVecFromRowStochVec(original_vars.s, row) = value_row;

      if( iclow == 1 )
      {
         if( !PIPSisLEFeas(lhs, value_row) )
            PIPSdebugMessage("Postsolve Warning: when reintroducing a redundant inequality row it did not meet its lhs with feastol: %f > %f", lhs, value_row);
         assert(PIPSisLE(lhs, value_row, postsolve_tol));
      }

      if( icupp == 1 )
      {
         if( !PIPSisLEFeas(value_row, rhs) )
            PIPSdebugMessage("Postsolve Warning: when reintroducing a redundant inequality row it did not meet its rhs with feastol: %f > %f", value_row, rhs);
         assert(PIPSisLE(value_row, rhs, postsolve_tol));
      }

      /* set correct slacks */
      if( iclow == 1)
         getSimpleVecFromRowStochVec(original_vars.t, row) = value_row - lhs;
      else
         getSimpleVecFromRowStochVec(original_vars.t, row) = 0;

      if( icupp == 1)
         getSimpleVecFromRowStochVec(original_vars.u, row) = rhs - value_row;
      else
         getSimpleVecFromRowStochVec(original_vars.u, row) = 0;

      assert( complementarySlackRowMet(original_vars, row, postsolve_tol) );
   }


   return true;
}

bool StochPostsolver::postsolveBoundsTightened(sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == BOUNDS_TIGHTENED );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 2 == next_first_float_val);
   assert(first_int_val + 2 == next_first_int_val);
#endif

   const INDEX& row = indices.at(first_index);
   const INDEX& col = indices.at(first_index + 1);

   assert( row.isRow() || row.isEmpty() );
   assert( col.isCol() );

   if( row.isRow() )
      assert( !wasRowRemoved(row) );
   assert( !wasColumnRemoved(col) );

   const bool at_root_node = row.isEmpty() ? false : row.getNode() == -1 && col.isLinkingCol();

#ifndef NDEBUG
   if( col.isLinkingCol() && !at_root_node )
   {
      assert( PIPS_MPIisValueEqual( col.getIndex() ) );
      const int my_tightening = row.isEmpty() ? 0 : 1;
      assert( PIPS_MPIgetSum(my_tightening) == 1);
   }
#endif

   const bool is_upper_bound = (int_values[first_int_val] == 1) ? true : false;
   const int index_stored_row = int_values[first_int_val + 1];

   const double old_bound = float_values[first_float_val];
#ifndef NDEBUG
   const double new_bound = float_values[first_float_val + 1];
#endif

   const double curr_x = getSimpleVecFromColStochVec(original_vars.x, col);

   if( is_upper_bound )
   {
      assert( PIPSisLT(new_bound, old_bound) );
      assert( PIPSisLEFeas(curr_x, new_bound) );
   }
   else
   {
      assert( PIPSisLT(old_bound, new_bound) );
      assert( PIPSisLEFeas(new_bound, curr_x) );
   }

   double& slack = is_upper_bound ? getSimpleVecFromColStochVec(original_vars.w, col) :
         getSimpleVecFromColStochVec(original_vars.v, col);
   double dual_bound = is_upper_bound ? getSimpleVecFromColStochVec(original_vars.phi, col) :
         getSimpleVecFromColStochVec(original_vars.gamma, col);
   if( col.isLinkingCol() )
      dual_bound += is_upper_bound ? (*phi_changes)[col.getIndex()] : (*gamma_changes)[col.getIndex()];
   /* If the bound was tight all other variables in that row must have been at their respective upper
    * and lower bounds (depending on sings and orientation and upper/lower).
    * This fact will be used to adjust their duals which can be non-zero then since v/w were zero.
    * TODO : assert that
    */

   // TODO : reintroduce once scaling also scales termination criteria
   //   assert( PIPSisLEFeas(0.0, complementarity) );
   // assert( PIPSisLT( complementarity, feas_tol ) );
   const double old_complementarity = slack * dual_bound;

   /* adjust slack v/w */
   if( std::fabs(old_bound) == INF_POS )
   {
      slack = 0.0;
   }
   else
   {
      if( is_upper_bound )
      {
         assert(PIPSisLE(0, new_bound - curr_x));
         slack = old_bound - curr_x;
      }
      else
      {
         assert(PIPSisLE(0, curr_x - new_bound));
         slack = curr_x - old_bound;
      }
      assert( PIPSisLT(0.0, slack) );

      /* if only adjusting the slack is sufficent and we do not need dual postsolve - return */
      if( PIPSisZeroFeas( dual_bound * slack ) )
         return true;
   }

   /* if this is a local linking variable and we did not do the adjustments we're done */
   if( !row.isRow() )
      return true;

   /* do dual part of postsolve */

   /* adjust column bound dual so that complementary slackness is still as valid as before */
   const double diff_dual_bound = PIPSisZero(slack) ? - dual_bound : std::fabs( old_complementarity ) / slack - dual_bound;
   assert( PIPSisLE(0.0, dual_bound + diff_dual_bound) );
   assert( PIPSisLEFeas(slack * (dual_bound + diff_dual_bound), std::fabs( old_complementarity) ) );

   /* if the change will not introduce a big error into the reduced costs just apply it */
   if( std::fabs(diff_dual_bound) < postsolve_tol * 1e-1 )
   {
      if( col.isLinkingCol() )
      {
         double& dual = is_upper_bound ? (*phi_changes)[col.getIndex()] : (*gamma_changes)[col.getIndex()];
         if( std::fabs(old_bound) == INF_POS )
            dual = 0.0;
         else
            dual += diff_dual_bound;
      }
      else
      {
         double& dual = is_upper_bound ? getSimpleVecFromColStochVec(original_vars.phi, col)
               : getSimpleVecFromColStochVec(original_vars.gamma, col);
         if( std::fabs(old_bound) == INF_POS )
            dual = 0.0;
         else
            dual += diff_dual_bound;
      }
      return true;
   }

   /* at this point we want to adjust the dual_bound by diff_dual_bound introducing an error in the reduced costs which we need to compensate */
   const INDEX stored_row(ROW, row.getNode(), index_stored_row, row.getLinking(), EQUALITY_SYSTEM);
   const double coeff = row_storage.getRowCoefficientAtColumn( stored_row, col );

   /* set z/y of corresponding row such that c_i * delta_z // a_i * delta_y = +- diff_dual_bound */
   assert( !PIPSisZero(coeff) );
   const double change_dual_row = is_upper_bound ? (diff_dual_bound / coeff) : (-diff_dual_bound / coeff);
   assert( PIPSisEQ(-coeff * change_dual_row, is_upper_bound ? -diff_dual_bound : diff_dual_bound) );

   /* add -dz/dy * row to gamma/phi */
   /* store linking variable changes and allreduce them later except when using a row from D/B0 or D/Bl0 */
   StochVector& gamma = dynamic_cast<StochVector&>(*original_vars.gamma);
   StochVector& phi = dynamic_cast<StochVector&>(*original_vars.phi);

   /* linking rows not at root - so far we changed the slacks for the variable - missing is dual for the row and resulting corrections */
   if( !at_root_node && row.isLinkingRow() )
   {
      if( row.inEqSys() )
      {
         (*y_changes)[row.getIndex()] += change_dual_row;
         outdated_equality_linking_rows= true;
      }
      else
      {
         (*z_changes)[row.getIndex()] += change_dual_row;
         outdated_inequality_linking_rows= true;
      }
      return true;
   }

   if( distributed && !at_root_node )
   {
      row_storage.axpyAtRowPosNeg(1.0, &phi, &(*phi_changes), &gamma, &(*gamma_changes), change_dual_row, stored_row );
      outdated_linking_vars = true;
   }
   else
      row_storage.axpyAtRowPosNeg(1.0, &phi, nullptr, &gamma, nullptr, change_dual_row, stored_row );

   /* adjust the row dual */
   if( row.inInEqSys() )
   {
      double& z = getSimpleVecFromRowStochVec(original_vars.z, row);
      double& lambda = getSimpleVecFromRowStochVec(original_vars.lambda, row);
      double& pi = getSimpleVecFromRowStochVec(original_vars.pi, row);

      addIneqRowDual(z, lambda, pi, change_dual_row);
   }
   else
      getSimpleVecFromRowStochVec(original_vars.y, row) += change_dual_row;

   return true;
}

bool StochPostsolver::postsolveFixedColumn(sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == FIXED_COLUMN );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 1 == next_first_index);
   assert(first_float_val + 2 == next_first_float_val);
   assert(first_int_val + 1 == next_first_int_val);
#endif

   const INDEX& col = indices.at(first_index);
   assert( col.isCol() );
   assert( wasColumnRemoved(col) );

   const int index_stored_col = int_values.at(first_int_val);
   const INDEX stored_col(COL, col.getNode(), index_stored_col);

   const double value = float_values.at(first_float_val);
   const double obj_coeff = float_values.at(first_float_val + 1);


   /* mark entry as set and set x value to fixation */
   markColumnAdded(col);

   /* set x value */
   getSimpleVecFromColStochVec(original_vars.x, col) = value;

   /* set slacks for x bounds to zero (bounds were tight) */
   getSimpleVecFromColStochVec(original_vars.v, col) = 0.0;
   getSimpleVecFromColStochVec(original_vars.w, col) = 0.0;

   /* set duals for bounds to satisfy reduced costs of reintroduced column times x */
   double col_times_duals = 0.0;
   if( col.isLinkingCol() )
   {
      assert( PIPS_MPIisValueEqual(col.getIndex(), MPI_COMM_WORLD) );
      /* we need to synchronize the column times duals in this case */
      if( my_rank == 0 )
         col_times_duals = col_storage.multColTimesVec(stored_col, dynamic_cast<const StochVector&>(*original_vars.y),
               dynamic_cast<const StochVector&>(*original_vars.z));
      else
         col_times_duals = col_storage.multColTimesVecWithoutRootNode(stored_col, dynamic_cast<const StochVector&>(*original_vars.y),
               dynamic_cast<const StochVector&>(*original_vars.z));

      PIPS_MPIgetSumInPlace(col_times_duals, MPI_COMM_WORLD);
   }
   else
   {
      col_times_duals = col_storage.multColTimesVec( stored_col, dynamic_cast<const StochVector&>(*original_vars.y),
            dynamic_cast<const StochVector&>(*original_vars.z));
   }
   const double reduced_costs = obj_coeff - col_times_duals;

   /* adjust slacks in inequalities */
   /* add col * value to slacks */
   /* for linking rows we have to store the slack changes if they are not happening in the root node */
   if( col.isLinkingCol() )
      col_storage.axpyAtCol( 1.0, nullptr, &dynamic_cast<StochVector&>(*original_vars.s), nullptr, nullptr, value, stored_col );
   else
   {
      outdated_inequality_linking_rows = true;
      col_storage.axpyAtCol( 1.0, nullptr, &dynamic_cast<StochVector&>(*original_vars.s), nullptr, s_changes, value, stored_col );
   }

   /* set duals of bounds of x */
   double& gamma = getSimpleVecFromColStochVec(original_vars.gamma, col);
   double& phi = getSimpleVecFromColStochVec(original_vars.phi, col);
   gamma = 0.0;
   phi = 0.0;

   if( PIPSisLT(reduced_costs, 0.0) )
      phi = -reduced_costs;
   else if( PIPSisLT(0.0, reduced_costs) )
      gamma = reduced_costs;

   assert( PIPSisZeroFeas(reduced_costs - gamma + phi) );
   return true;
}

/**
 * recover primal value
 * dual multiplies will be set to zero
 * compute slack variables
 * no special treatment for linking variables since all processes should fix them simultaneously and in the same order
 *    -> assert is in place to check this
 */
bool StochPostsolver::postsolveFixedEmptyColumn(sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == FIXED_EMPTY_COLUMN );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);

   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 1 == next_first_index);
   assert(first_float_val + 4 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   const INDEX& col = indices.at(first_index);
   assert( col.isCol() );
   assert( wasColumnRemoved(col) );
   if( col.isLinkingCol() )
      assert( PIPS_MPIisValueEqual(col.getIndex(), MPI_COMM_WORLD) );

   const double value = float_values.at(first_float_val);
   const double obj_coeff = float_values.at(first_float_val + 1);
   const double xlow = float_values.at(first_float_val + 2);
   const double xupp = float_values.at(first_float_val + 3);

   /* primal */
   /* mark entry as set and set x value to fixation */
   markColumnAdded(col);
   getSimpleVecFromColStochVec(original_vars.x, col) = value;

   assert( PIPSisLEFeas(xlow, value) );
   assert( PIPSisLEFeas(value, xupp) );

   /* dual */
   getSimpleVecFromColStochVec(original_vars.gamma, col) = 0.0;
   getSimpleVecFromColStochVec(original_vars.phi, col) = 0.0;

   if( !PIPSisZero(obj_coeff) )
   {
      if( PIPSisLT(obj_coeff, 0.0) )
      {
         assert( xupp != INF_POS );
         assert( PIPSisEQ(value, xupp) );
         getSimpleVecFromColStochVec(original_vars.phi, col) = obj_coeff;
      }
      else if( PIPSisLT(0.0, obj_coeff) )
      {
         assert( xlow != INF_NEG );
         assert( PIPSisEQ(value, xlow) );
         getSimpleVecFromColStochVec(original_vars.gamma, col) = obj_coeff;
      }
   }

   if( xlow != INF_NEG )
      getSimpleVecFromColStochVec(original_vars.v, col) = value - xlow;
   else
      getSimpleVecFromColStochVec(original_vars.v, col) = 0.0;

   if( xupp != INF_POS )
      getSimpleVecFromColStochVec(original_vars.w, col) = xupp - value;
   else
      getSimpleVecFromColStochVec(original_vars.w, col) = 0.0;

   assert( complementarySlackVariablesMet(original_vars, col, postsolve_tol) );

   return true;
}

bool StochPostsolver::postsolveFixedColumnSingletonFromInequality(sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == FIXED_COLUMN_SINGLETON_FROM_INEQUALITY );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 4 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   const INDEX& col = indices.at(first_index);
   const INDEX& row = indices.at(first_index + 1);
   assert( row.isRow() );
   assert( row.inInEqSys() );
   assert( col.isCol() );
   assert( !wasColumnRemoved(col) );

   const double value = float_values.at(first_float_val);
//   const double coeff = float_values.at(first_float_val + 1); // TODO : remove
   const double xlow_old = float_values.at(first_float_val + 2);
   const double xupp_old = float_values.at(first_float_val + 3);

   const bool local_linking_col = col.isLinkingCol() && row.getNode() != -1;

   /* set x value - bound is tight - compute slacks of x - set duals to zero */
   assert( PIPSisEQ( getSimpleVecFromColStochVec(original_vars.x, col), value) );

   /* set slacks for x bounds */
   assert( PIPSisLE(value, xupp_old) || xupp_old == INF_POS);
   assert( PIPSisLE(xlow_old, value) || xlow_old == INF_NEG);

   double& v = local_linking_col ? (*v_changes)[col.getIndex()] : getSimpleVecFromColStochVec(original_vars.v, col);
   double& w = local_linking_col ? (*w_changes)[col.getIndex()] : getSimpleVecFromColStochVec(original_vars.w, col);
   assert( PIPSisZero(v) && PIPSisZero(w) );
   v = (xlow_old == INF_NEG) ? 0 : value - xlow_old;
   w = (xupp_old == INF_POS) ? 0 : xupp_old - value;

   /* one of the bounds has to stay active */
   assert( PIPSisZero(v) || PIPSisZero(w) );

   /* if a bound is no longer tight we do not have to adjust the respective duals since these should be zero anyway */
   if( !PIPSisZero(v) )
      assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.gamma, col) ) );
   if( !PIPSisZero(w) )
      assert( PIPSisZero( getSimpleVecFromColStochVec(original_vars.phi, col) ) );

   assert( complementarySlackVariablesMet( original_vars, col, postsolve_tol) );

   return true;
}

bool StochPostsolver::postsolveSingletonEqualityRow(sVars& original_vars, int reduction_idx) const
{
   assert( reductions.at(reduction_idx) == SINGLETON_EQUALITY_ROW );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 5 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   const INDEX& row = indices.at(first_index);
   const INDEX& col = indices.at(first_index + 1);

   if( row.isRow() )
   {
      assert( row.getSystemType() == EQUALITY_SYSTEM );
      assert( !wasRowRemoved(row) );
   }
   assert( !wasColumnRemoved(col) );

   const double xlow_old = float_values.at(first_float_val);
   const double xupp_old = float_values.at(first_float_val + 1);
   const double coeff = float_values.at(first_float_val + 4);

   assert( !PIPSisZero(coeff) || coeff == NAN );

   const double curr_x = getSimpleVecFromColStochVec(original_vars.x, col);
   double& slack_lower = getSimpleVecFromColStochVec(original_vars.v, col);
   double& slack_upper = getSimpleVecFromColStochVec(original_vars.w, col);

   double& dual_lower = getSimpleVecFromColStochVec(original_vars.gamma, col);
   double& dual_upper = getSimpleVecFromColStochVec(original_vars.phi, col);

   double error_in_reduced_costs = 0.0;
   /* adjust the slacks */
   /* adjust duals in reduced costs and compensate the error with the dual of the newly introduced row */

   /* upper bound */
   if( xupp_old == INF_POS )
   {
      slack_upper = 0.0;
      error_in_reduced_costs -= dual_upper;
      dual_upper = 0.0;
   }
   else
   {
      slack_upper = xupp_old - curr_x;
      assert( PIPSisLE(0.0, slack_upper) );
      assert( std::fabs(slack_upper) != INF_POS );

      if( !PIPSisZero( slack_upper * dual_upper, postsolve_tol ) )
      {
         error_in_reduced_costs -= dual_upper;
         dual_upper = 0.0;
      }
   }

   /* lower bound */
   if(xlow_old == INF_NEG)
   {
      slack_lower = 0.0;
      error_in_reduced_costs += dual_lower;
      dual_lower = 0.0;
   }
   else
   {
      slack_lower = curr_x - xlow_old;
      assert(PIPSisLE(0.0, slack_lower));
      assert(std::fabs(slack_lower) != INF_POS);

      if( !PIPSisZero( slack_lower * dual_lower, postsolve_tol ) )
      {
         error_in_reduced_costs += dual_lower;
         dual_lower = 0.0;
      }
   }

   assert( PIPSisZero( slack_lower * dual_lower, postsolve_tol) );
   assert( PIPSisZero( slack_upper * dual_upper, postsolve_tol) );

   /* adjust duals of bounds so that complementary slackness is still met - iff we are the process that has the row to adjust on it */
   if( row.isRow() )
   {
      if( !PIPSisZero(error_in_reduced_costs) )
      {
         double& dual_singelton_row = getSimpleVecFromRowStochVec(original_vars.y, row);
         assert( coeff != NAN );
         /* we use the coeff * dual_singleton_row to balance the error in the reduced costs */
         const double diff_dual_row = error_in_reduced_costs / coeff;
         dual_singelton_row += diff_dual_row;
         assert(PIPSisEQ(diff_dual_row * coeff, error_in_reduced_costs));
      }
   }
   return true;
}

bool StochPostsolver::postsolveSingletonInequalityRow(sVars& original_vars, int reduction_idx) const
{
   assert( reductions.at(reduction_idx) == SINGLETON_INEQUALITY_ROW );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 5 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   const INDEX& row = indices.at(first_index);
   const INDEX& col = indices.at(first_index + 1);

   if( row.isRow() )
   {
      assert(row.getSystemType() == INEQUALITY_SYSTEM);
      assert( !wasRowRemoved(row) );
   }
   assert( !wasColumnRemoved(col) );

   const double xlow_old = float_values.at(first_float_val);
   const double xupp_old = float_values.at(first_float_val + 1);
   const double xlow_new = float_values.at(first_float_val + 2);
#ifndef NDEBUG
   const double xupp_new = float_values.at(first_float_val + 3);
#endif
   const double coeff = float_values.at(first_float_val + 4);

   assert( !PIPSisZero(coeff) || coeff == NAN );
   assert( xlow_new == INF_NEG || xupp_new == INF_POS );
   assert( xlow_new != INF_NEG || xupp_new != INF_POS );

   bool lower_bound_changed = (xlow_new != INF_NEG);

   const double curr_x = getSimpleVecFromColStochVec(original_vars.x, col);
   double& slack = lower_bound_changed ? getSimpleVecFromColStochVec(original_vars.v, col) :
         getSimpleVecFromColStochVec(original_vars.w, col);
   double& dual_bound = lower_bound_changed ? getSimpleVecFromColStochVec(original_vars.gamma, col) :
         getSimpleVecFromColStochVec(original_vars.phi, col);

   const double old_bound = lower_bound_changed ? xlow_old : xupp_old;
#ifndef NDEBUG
   const double new_bound = lower_bound_changed ? xlow_new : xupp_new; // TODO : remove
   if( lower_bound_changed )
      assert( PIPSisLT( old_bound, new_bound) );
   else
      assert( PIPSisLT( new_bound, old_bound) );
#endif

   double error_in_reduced_costs = 0.0;
   /* adjust bounds slacks and duals so that complementarity condition stays valid*/
   if(std::fabs(old_bound) == INF_POS)
   {
      slack = 0.0;

      error_in_reduced_costs += lower_bound_changed ? dual_bound : -dual_bound;
      dual_bound = 0.0;
   }
   else
   {
      slack = std::fabs(old_bound - curr_x);
      assert(PIPSisLT(0.0, slack));
      assert(std::fabs(slack) != INF_POS);

      if( !PIPSisZero( slack * dual_bound, postsolve_tol ) )
      {
         error_in_reduced_costs += lower_bound_changed ? dual_bound : -dual_bound;
         dual_bound = 0.0;
      }
   }
   assert(PIPSisZero(dual_bound * slack, postsolve_tol));

   /* if we are the corresponding process owning the row and we have an error in the reduced costs */
   if( row.isRow() )
   {
      /* correct optimality conditions and reduced costs */
      if( !PIPSisZero(error_in_reduced_costs) )
      {
         double& dual_singelton_row = getSimpleVecFromRowStochVec(original_vars.z, row);
         assert( PIPSisZero(dual_singelton_row) );
         assert( PIPSisZero(getSimpleVecFromRowStochVec(original_vars.lambda, row)) );
         assert( PIPSisZero(getSimpleVecFromRowStochVec(original_vars.pi, row)) );

         /* we use the coeff * dual_singleton_row to balance the error in the reduced costs */
         const double diff_dual_row = error_in_reduced_costs / coeff;
         dual_singelton_row += diff_dual_row;

         getSimpleVecFromRowStochVec(original_vars.lambda, row) = std::max(0.0, dual_singelton_row);
         getSimpleVecFromRowStochVec(original_vars.pi, row) = -std::min(0.0, dual_singelton_row);

         assert(PIPSisEQFeas(diff_dual_row * coeff, error_in_reduced_costs));
      }
   }
   return true;
}

bool StochPostsolver::postsolveFreeColumnSingletonEquality(sVars& original_vars, int reduction_idx)
{
   /* row can be an equality row but then it must have clow == cupp */
   assert( reductions.at(reduction_idx) == FREE_COLUMN_SINGLETON_EQUALITY );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 5 == next_first_float_val);
   assert(first_int_val + 1 == next_first_int_val);
#endif

   const INDEX& col = indices.at(first_index);
   const INDEX& row = indices.at(first_index + 1);

   assert( row.isRow() );

   if( !row.isLinkingRow() )
      assert( col.isCol() );
   if( !col.isCol() )
      assert( row.isLinkingRow() );

   const double rhs = float_values.at(first_float_val);
   const double obj_coeff = float_values.at(first_float_val + 1);
   const double col_coeff = float_values.at(first_float_val + 2);
   const double xlow = float_values.at(first_float_val + 3);
   const double xupp = float_values.at(first_float_val + 4);

   const int stored_row_idx = int_values.at(first_int_val);
   const INDEX stored_row(ROW, row.getNode(), stored_row_idx, row.getLinking(), row.getSystemType());

   assert( !PIPSisZero(col_coeff) );
   assert( wasRowRemoved(row) );

   /* assert that linking columns and rows are done by all processes simultaneously */
   if( !col.isEmpty() )
      if( col.isLinkingCol() && row.getNode() == -1 )
         assert(PIPS_MPIisValueEqual(col.getIndex(), MPI_COMM_WORLD));

   if( row.isLinkingRow() )
      assert(PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD));

   /* reintroduce row on all processes */
   if( row.inEqSys() )
      getSimpleVecFromRowStochVec(*padding_origrow_equality, row) = 1;
   else
      getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) = 1;

   const double dual_value_row = obj_coeff / col_coeff;

   if( row.isLinkingRow() )
      assert( PIPS_MPIisValueEqual(dual_value_row, MPI_COMM_WORLD) );

   /* set duals of row depending on equality/inequality row */
   if( row.inEqSys() )
   {
      getSimpleVecFromRowStochVec(original_vars.y, row) = dual_value_row;
      assert( PIPSisZero(obj_coeff - getSimpleVecFromRowStochVec(original_vars.y, row) * col_coeff) );
   }
   else
   {
      /* cupp == clow so slacks must be zero close to zero */
      getSimpleVecFromRowStochVec(original_vars.s, row) = 0.0;
      getSimpleVecFromRowStochVec(original_vars.t, row) = 0.0;
      getSimpleVecFromRowStochVec(original_vars.u, row) = 0.0;

      /* set dual */
      getSimpleVecFromRowStochVec(original_vars.z, row) = dual_value_row;
      if( PIPSisLT(0.0, dual_value_row) )
         getSimpleVecFromRowStochVec(original_vars.pi, row) = dual_value_row;
      else
         getSimpleVecFromRowStochVec(original_vars.lambda, row) = -dual_value_row;
   }

   /* synchronize value of row for x_val */
   double value_row = 0.0;
   if( col.isCol() )
      assert( PIPSisZero(getSimpleVecFromColStochVec(*original_vars.x, col)) );

   if( row.isLinkingRow() )
   {
      assert(PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD));

      if(my_rank == 0)
         value_row = row_storage.multRowTimesVec( stored_row, dynamic_cast<const StochVector&>(*original_vars.x));
      else
         value_row = row_storage.multLinkingRowTimesVecWithoutBl0( stored_row_idx, dynamic_cast<const StochVector&>(*original_vars.x));
      PIPS_MPIgetSumInPlace(value_row, MPI_COMM_WORLD);
   }
   else
      value_row = row_storage.multRowTimesVec(stored_row, dynamic_cast<const StochVector&>(*original_vars.x));

   assert(std::abs(value_row) != INF_POS);

   /* reintroduce the removed column on process owning column */
   if( col.isCol() )
   {
      const bool local_col_change = col.isLinkingCol() && row.getNode() != -1;

      assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.x, col)) );
      if( local_col_change )
         outdated_linking_vars = true;
      double& x_val = local_col_change ? (*x_changes)[col.getIndex()] : getSimpleVecFromColStochVec(original_vars.x, col);

      /* mark column as set */
      getSimpleVecFromColStochVec(*padding_origcol, col) = 1;

      /* recover primal value */
      x_val = (rhs - value_row) / col_coeff;
      assert( PIPSisZero(x_val * col_coeff + value_row - rhs) );

      /* compute slacks and set duals for bounds to zero */
      double& slack_lower = local_col_change ? (*v_changes)[col.getIndex()] : getSimpleVecFromColStochVec(original_vars.v, col);
      double& slack_upper = local_col_change ? (*w_changes)[col.getIndex()] : getSimpleVecFromColStochVec(original_vars.w, col);

      if( xlow == INF_NEG )
         slack_lower = 0.0;
      else
      {
         assert( PIPSisLE(xlow, x_val) );
         slack_lower = xlow - x_val;
      }

      if( xupp == INF_POS )
         slack_upper = 0.0;
      else
      {
         assert( PIPSisLE(x_val, xupp) );
         slack_upper = xupp - x_val;
      }

      getSimpleVecFromColStochVec(original_vars.gamma, col) = 0.0;
      getSimpleVecFromColStochVec(original_vars.phi, col) = 0.0;
   }

   return true;
}

bool StochPostsolver::postsolveNearlyParallelRowSubstitution(sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == NEARLY_PARALLEL_ROW_SUBSTITUTION );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 4 == next_first_index);
   assert(first_float_val + 9 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   /* col2 was substituted by col1 via col2 = t * col1 + d */
   const INDEX& col1 = indices.at(first_index);
   const INDEX& col2 = indices.at(first_index + 1);
   const INDEX& row1 = indices.at(first_index + 2);
   const INDEX& row2 = indices.at(first_index + 3);

   assert( row1.isRow() );
   assert( row2.isRow() );
   assert( row1.getNode() == row2.getNode() );
   assert( ( row1.inEqSys() && row2.inEqSys() ) ||
      ( row1.inInEqSys() && row2.inInEqSys() ) );
   assert( col2.isCol() );
   if( row1.inInEqSys() )
      assert( col1.isCol() );

   if( col1.isCol() )
      assert( !wasColumnRemoved(col1) );
   assert( !wasColumnRemoved(col2) );

   assert( !wasRowRemoved(row1) );
   assert( !wasRowRemoved(row2) );

   const bool local_linking_col1 = col1.isCol() ? (col1.isLinkingCol() && row1.getNode() != -1) : false;
   const bool local_linking_col2 = col2.isLinkingCol() && row2.getNode() != -1;

   const double scalar = float_values.at(first_float_val);
   const double translation = float_values.at(first_float_val + 1);
   const double obj_col1 = float_values.at(first_float_val + 2);
   const double obj_col2 = float_values.at(first_float_val + 3);

   const double xlow_col2 = float_values.at(first_float_val + 4);
   const double xupp_col2 = float_values.at(first_float_val + 5);
#ifndef NDEBUG
   const double coeff_col1 = float_values.at(first_float_val + 6);
#endif
   const double coeff_col2 = float_values.at(first_float_val + 7);

   /* row1 = parallel_factor * row2 */
   const double parallel_factor = float_values.at(first_float_val + 8);

   assert( !PIPSisZero(coeff_col2) );
   if( row2.inInEqSys() )
   {
      assert( !PIPSisZero(coeff_col1) );
      assert( PIPSisLT(0.0, coeff_col1 * coeff_col2) );
      assert( PIPSisZero(translation) );
   }

   const double val_col1 = col1.isCol() ? getSimpleVecFromColStochVec(original_vars.x, col1) : 0.0;
   const double val_col2 = scalar * val_col1 + translation;

   /* primal postsolve */

   /* reintroduce the substituted column */
   if( local_linking_col2 )
   {
      (*x_changes)[col2.getIndex()] += val_col2 - getSimpleVecFromColStochVec(original_vars.x, col2);
      outdated_linking_vars = true;
   }
   else
      getSimpleVecFromColStochVec(original_vars.x, col2) = val_col2;

   assert( PIPSisLEFeas( xlow_col2, val_col2 ) );
   assert( PIPSisLEFeas( val_col2, xupp_col2 ) );

   assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.phi, col2)) );
   assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.gamma, col2)) );
   /* slacks for bounds */
   /* bound duals to zero */
   if( local_linking_col2 )
   {
      (*v_changes)[col2.getIndex()] += (xlow_col2 == INF_NEG) ? -getSimpleVecFromColStochVec(original_vars.v, col2) : (*x_changes)[col2.getIndex()];
      (*w_changes)[col2.getIndex()] -= (xupp_col2 == INF_POS) ? getSimpleVecFromColStochVec(original_vars.w, col2) : (*x_changes)[col2.getIndex()];
      outdated_linking_vars = true;
   }
   else
   {
      getSimpleVecFromColStochVec(original_vars.v, col2) = (xlow_col2 == INF_NEG) ? 0 : val_col2 - xlow_col2;
      getSimpleVecFromColStochVec(original_vars.w, col2) = (xupp_col2 == INF_POS) ? 0 : xupp_col2 - val_col2;
   }

   /* dual postsolve substitution */

   /* the substitution itself needs a shift of the row duals to compensate the objective vector change or, if the row was not tight (so only for inequalities)
    * a shift of the column duals
    */

   /* postsolve two equalities */
   if( row1.inEqSys() )
   {
      double& dual_row1 = getSimpleVecFromRowStochVec(original_vars.y, row1);
      double& dual_row2 = getSimpleVecFromRowStochVec(original_vars.y, row2);

      assert( PIPSisZero(dual_row2) );

      /* compensate objective vector change with row2 and row1 */
      dual_row2 = obj_col2 / coeff_col2;
      dual_row1 -= ( dual_row2 / parallel_factor );

      assert( PIPSisZero(obj_col2 - coeff_col2 * dual_row2) );
   }
   /* postsolve two inequalities */
   else
   {
      /* objective coefficients had the same sign */
      assert( row1.inInEqSys() );
      assert( row2.inInEqSys() );
      assert( !PIPSisZero( parallel_factor ) );
      assert( !PIPSisZero( coeff_col2 ) );
      assert( PIPSisEQ( scalar, coeff_col1 / ( parallel_factor * coeff_col2 ) ) );
      assert( PIPSisLE( 0.0, scalar ) );

      /* if the bounds were tight we have to shift their duals - else we have to shift the row duals */
      const double change_objective_row1 = obj_col2 * scalar;

      assert( !PIPSisZero(obj_col1 + change_objective_row1) );

      /* scale row 1 by factor_row1 */
      const double factor_row1 = obj_col1 / (obj_col1 + change_objective_row1);

      /* shift duals from row1 to row2 */
      double& z_row1 = getSimpleVecFromRowStochVec(original_vars.z, row1);
      double& lambda_row1 = getSimpleVecFromRowStochVec(original_vars.lambda, row1);
      double& pi_row1 = getSimpleVecFromRowStochVec(original_vars.pi, row1);

      double& z_row2 = getSimpleVecFromRowStochVec(original_vars.z, row2);
      double& lambda_row2 = getSimpleVecFromRowStochVec(original_vars.lambda, row2);
      double& pi_row2 = getSimpleVecFromRowStochVec(original_vars.pi, row2);

      assert( PIPSisZero( z_row2 ) );
      assert( PIPSisZero( pi_row2 ) );
      assert( PIPSisZero( lambda_row2 ) );

      /* bound duals */
      double& gamma_row1 = getSimpleVecFromRowStochVec(original_vars.gamma, row1);
      double& phi_row1 = getSimpleVecFromRowStochVec(original_vars.phi, row1);

      double& gamma_row2 = getSimpleVecFromRowStochVec(original_vars.gamma, row2);
      double& phi_row2 = getSimpleVecFromRowStochVec(original_vars.phi, row2);

#ifndef NDEBUG
      assert( PIPSisZero(gamma_row2) );
      assert( PIPSisZero(phi_row2) );
      assert( PIPSisZero(z_row1 - lambda_row1 + pi_row1, postsolve_tol) );
      assert( PIPSisZero(z_row2 - lambda_row2 + pi_row2, postsolve_tol) );

      assert( PIPSisZero( (obj_col1 + scalar * obj_col2) - coeff_col1 * z_row1 - gamma_row1 + phi_row1, postsolve_tol ) );
      const double gamma_row1_old = gamma_row1;
      const double phi_row1_old = phi_row1;
#endif

      /* z = lambda - pi */
      const double delta_z_row1 = (factor_row1 - 1) * z_row1;
      const double delta_z_row2 = - parallel_factor * delta_z_row1;
      const double row2_bound_duals = obj_col2 - coeff_col2 * delta_z_row2;

      assert( PIPSisZero( delta_z_row1 + delta_z_row2 / parallel_factor, postsolve_tol ) );

      addIneqRowDual(z_row2, lambda_row2, pi_row2, delta_z_row2);
      addIneqRowDual(z_row1, lambda_row1, pi_row1, delta_z_row1);

      /* row1 bound duals */
      if( local_linking_col2 )
      {
         const int col2_idx = col2.getIndex();
         (*gamma_changes)[col2_idx] = std::max(0.0, row2_bound_duals);
         (*phi_changes)[col2_idx] = std::min(0.0, row2_bound_duals);
         outdated_linking_vars = true;

         assert( PIPSisZero( obj_col2 - coeff_col2 * z_row2 - (gamma_row2 + (*gamma_changes)[col2_idx]) + (phi_row2 + (*phi_changes)[col2_idx]), postsolve_tol ) );
      }
      else
      {
         gamma_row2 = std::max(0.0, row2_bound_duals);
         phi_row2 = std::min(0.0, row2_bound_duals);

         assert( PIPSisZero( obj_col2 - coeff_col2 * z_row2 - gamma_row2 + phi_row2, postsolve_tol ) );
      }

      if( local_linking_col1 )
      {
         const int col1_idx = col1.getIndex();
         (*gamma_changes)[col1_idx] = (factor_row1 - 1) * gamma_row1;
         (*phi_changes)[col1_idx] = (factor_row1 - 1) * phi_row1;
         outdated_linking_vars = true;

         assert( PIPSisZero( obj_col1 - coeff_col1 * z_row1 - (gamma_row1 + (*gamma_changes)[col1_idx]) + (phi_row1 + (*phi_changes)[col1_idx]), postsolve_tol ) );
      }
      else
      {
         gamma_row1 *= factor_row1;
         phi_row1 *= factor_row1;

         assert( PIPSisZero( obj_col1 - coeff_col1 * z_row1 - gamma_row1 + phi_row1, postsolve_tol ) );
      }

      assert( PIPSisZero( obj_col1 - coeff_col1 * z_row1 - gamma_row1_old * factor_row1 + phi_row1_old * factor_row1, postsolve_tol ) );
      assert( PIPSisZero( obj_col2 - coeff_col2 * z_row2 - row2_bound_duals, postsolve_tol ) );

      assert( complementarySlackVariablesMet(original_vars, col1, postsolve_tol) );
      assert( complementarySlackVariablesMet(original_vars, col2, postsolve_tol) );
      assert( complementarySlackRowMet( original_vars, row1, postsolve_tol) );
      assert( complementarySlackRowMet( original_vars, row2, postsolve_tol) );
   }
   return true;
}

bool StochPostsolver::postsolveNearlyParallelRowBoundsTightened(sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == NEARLY_PARALLEL_ROW_BOUNDS_TIGHTENED );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 4 == next_first_index);
   assert(first_float_val + 12 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   /* col2 was substituted by col1 via col2 = t * col1 + d */
   /* dual postsolve implied bounds */
   const INDEX& row1 = indices.at(first_index);
   const INDEX& row2 = indices.at(first_index + 1);
   const INDEX& col1 = indices.at(first_index + 2);
   const INDEX& col2 = indices.at(first_index + 3);

   assert( row1.isRow() );
   assert( row2.isRow() );
   assert( col1.isCol() );
   assert( col2.isCol() || col2.isEmpty() );

   if( col2.isCol() )
      assert( row1.inEqSys() && row2.inEqSys() );

   const bool local_linking_col1 = col1.isLinkingCol() && row1.getNode() != -1;
   const bool local_linking_col2 = col2.isCol() ? (col2.isLinkingCol() && row2.getNode() != -1) : false;

   const double xlow_col1 = float_values.at(first_float_val);
   const double xupp_col1 = float_values.at(first_float_val + 1);
   const double xlow_col2 = float_values.at(first_float_val + 2);
   const double xupp_col2 = float_values.at(first_float_val + 3);

   const double coeff_col1 = float_values.at(first_float_val + 4);
   const double coeff_col2 = float_values.at(first_float_val + 5);

   const double scalar = float_values.at(first_float_val + 6);
   const double translation = float_values.at(first_float_val + 7);
   const double parallel_factor = float_values.at(first_float_val + 8);

   const double rhs = float_values.at(first_float_val + 9);
   const double clow = float_values.at(first_float_val + 10);
   const double cupp = float_values.at(first_float_val + 11);

   assert( !PIPSisZero(scalar) );
#ifndef NDEBUG
   const double val_col1 = getSimpleVecFromColStochVec(original_vars.x, col1);
#endif

   /* if the variable bound of col1 was actually implied via col2 we have to shift it's dual multipliers over via also adjusting the dual of row1 */
   assert( PIPSisLE( xlow_col1, val_col1, postsolve_tol ) );
   assert( PIPSisLE( val_col1, xupp_col1, postsolve_tol ) );

   double xlow_implied = INF_NEG;
   double xupp_implied = INF_POS;

   /* if two nearly parallel equality rows get - get implied bounds */
   if( col2.isCol() )
   {
      assert( xlow_col2 != INF_NEG || xupp_col2 != INF_POS );

      if( PIPSisLT(0.0, scalar) )
      {
         if( xlow_col2 != INF_NEG )
            xlow_implied = std::max(xlow_col1, (xlow_col2 - translation) / scalar);
         else
            xlow_implied = xlow_col1;

         if( xupp_col2 != INF_POS )
            xupp_implied = std::min(xupp_col1, (xupp_col2 - translation) / scalar);
         else
            xupp_implied = xupp_col1;
      }
      else
      {
         if( xlow_col2 != INF_NEG )
            xupp_implied = std::min(xupp_col1, (xlow_col2 - translation) / scalar);
         else
            xupp_implied = xupp_col1;

         if( xupp_col2 != INF_POS )
            xlow_implied = std::max(xlow_col1, (xupp_col2 - translation) / scalar);
         else
            xlow_implied = xlow_col1;
      }
   }
   else
   {
      // TODO : check when example exists
      assert( clow != INF_NEG || cupp != INF_POS );
      assert( !PIPSisZero(coeff_col1) );
      assert( !PIPSisZero(parallel_factor * coeff_col1) );

      const double faq = parallel_factor * coeff_col1;

      assert( !PIPSisZero(faq) );
      if( PIPSisLT( 0.0, faq ) )
      {
         if( cupp != INF_POS )
            xlow_implied = std::max( xlow_col1, (rhs - parallel_factor * cupp ) / coeff_col1 );
         else
            xlow_implied = xlow_col1;

         if( clow != INF_NEG )
            xupp_implied = std::min( xupp_col1, (rhs - parallel_factor * clow) / coeff_col1 );
         else
            xupp_implied = xupp_col1;
      }

      if( PIPSisLT( faq, 0.0 ) )
      {
         if( cupp != INF_POS )
            xupp_implied = std::min( xupp_col1, (rhs - parallel_factor * cupp) / coeff_col1 );
         else
            xupp_implied = xupp_col1;

         if( clow != INF_NEG )
            xlow_implied = std::max( xlow_col1, (rhs - parallel_factor * clow) / coeff_col1 );
         else
            xlow_implied = xlow_col1;
      }
   }

   assert( xlow_implied != INF_NEG || xupp_implied != INF_POS );
   assert( PIPSisLE(xlow_col1, xlow_implied, postsolve_tol) );
   assert( PIPSisLE(xupp_implied, xupp_col1, postsolve_tol) );

#ifndef NDEBUG
   const double old_slack_lower = getSimpleVecFromColStochVec(original_vars.v, col1);
   const double old_slack_upper = getSimpleVecFromColStochVec(original_vars.w, col1);
#endif

   assert( complementarySlackVariablesMet(original_vars, col1, postsolve_tol) );
   if( col2.isCol() )
      assert( complementarySlackVariablesMet(original_vars, col2, postsolve_tol) );

   /* lower bound was implied by substituted column + it's row */
   if( PIPSisLT( xlow_col1, xlow_implied ) )
   {
      assert( PIPSisLT( xlow_col1, val_col1 ) );

      if( local_linking_col1 )
         outdated_linking_vars = true;

      const double gamma_col1_old = local_linking_col1 ? getSimpleVecFromColStochVec(original_vars.gamma, col1) + (*gamma_changes)[col1.getIndex()] :
            getSimpleVecFromColStochVec(original_vars.gamma, col1);
      const double v_col1_old = local_linking_col1 ? getSimpleVecFromColStochVec(original_vars.v, col1) + (*v_changes)[col1.getIndex()] :
            getSimpleVecFromColStochVec(original_vars.v, col1);
      double& gamma_col1 = local_linking_col1 ? (*gamma_changes)[col1.getIndex()] : getSimpleVecFromColStochVec(original_vars.gamma, col1);
      double& v_col1 = local_linking_col1 ? (*v_changes)[col1.getIndex()] : getSimpleVecFromColStochVec(original_vars.v, col1);

      /* reset bounds and adjust slacks */
      if( xlow_col1 != INF_NEG )
         v_col1 += xlow_implied - xlow_col1;
      else
         v_col1 -= v_col1_old;

      /* shift the dual of col1 to col2 */
      if( !PIPSisZero(gamma_col1_old) )
      {
         /* set the dual of col1 to zero, compensate for the error made with the dual of row1, compensate that error with the dual of col2
          * and finally compensate the error in col2's row of the reduced costs with col2's dual
          */

         /* dual to zero and compensate error */
         const double dual_shift_row1 = gamma_col1_old / coeff_col1;
         assert( PIPSisEQ( coeff_col1 * dual_shift_row1, gamma_col1_old ) );
         gamma_col1 -= gamma_col1_old;

         getSimpleVecFromRowStochVec(original_vars.y, row1) += dual_shift_row1;

         /* compensate again and adjust dual of col2 */
         const double dual_shift_row2 = -dual_shift_row1 * parallel_factor;

         if( row2.inEqSys() )
            getSimpleVecFromRowStochVec(original_vars.y, row2) += dual_shift_row2;
         else
         {
            double& z = getSimpleVecFromRowStochVec(original_vars.z, row2);
            double& lambda = getSimpleVecFromRowStochVec(original_vars.lambda, row2);
            double& pi = getSimpleVecFromRowStochVec(original_vars.pi, row2);

            addIneqRowDual(z, lambda, pi, dual_shift_row2);
         }

         /* if we had two equality rows there was a singleton col in row2 an it's dual needs to be adjusted */
         if( col2.isCol() )
         {
            /* assert bounds of substituted variable are similarly tight */
            if( PIPSisLT(0.0, scalar) )
               assert( PIPSisEQFeas( getSimpleVecFromColStochVec(original_vars.x, col2) - xlow_col2, std::fabs(old_slack_lower * scalar) ) );
            else
               assert( PIPSisEQFeas( xupp_col2 - getSimpleVecFromColStochVec(original_vars.x, col2), std::fabs(old_slack_lower * scalar) ) );

            const double dual_shift_col2 = -coeff_col2 * dual_shift_row2;
            PIPSisLT(0.0, scalar) ? assert( PIPSisLE(0.0, dual_shift_col2) ) : assert( PIPSisLE( dual_shift_col2, 0.0) );

            if( PIPSisLT(0.0, scalar) )
            {
               assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.gamma, col2)) );
               double& gamma = local_linking_col2 ? (*gamma_changes)[col2.getIndex()] : getSimpleVecFromColStochVec(original_vars.gamma, col2);

               gamma += std::fabs(dual_shift_col2);
            }
            else
            {
               assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.phi, col2)) );
               double& phi = local_linking_col2 ? (*phi_changes)[col2.getIndex()] : getSimpleVecFromColStochVec(original_vars.phi, col2);

               phi += std::fabs(dual_shift_col2);
            }
         }
      }
   }

   /* upper bound was implied by substituted column + it's row */
   if( PIPSisLT( xupp_implied, xupp_col1 ) )
   {
      assert( PIPSisLT(val_col1, xupp_col1) );

      if( local_linking_col1 )
         outdated_linking_vars = true;

      const double phi_col1_old = local_linking_col1 ? getSimpleVecFromColStochVec(original_vars.phi, col1) + (*phi_changes)[col1.getIndex()] :
            getSimpleVecFromColStochVec(original_vars.phi, col1);
      const double w_col1_old = local_linking_col1 ? getSimpleVecFromColStochVec(original_vars.w, col1) + (*w_changes)[col1.getIndex()] :
            getSimpleVecFromColStochVec(original_vars.w, col1);
      double& phi_col1 = local_linking_col1 ? (*phi_changes)[col1.getIndex()] : getSimpleVecFromColStochVec(original_vars.phi, col1);
      double& w_col1 = local_linking_col1 ? (*w_changes)[col1.getIndex()] : getSimpleVecFromColStochVec(original_vars.w, col1);

      /* reset bounds and adjust slacks */
      if( xupp_col1 != INF_POS )
         w_col1 += xupp_col1 - xupp_implied;
      else
         w_col1 -= w_col1_old;

      if( !PIPSisZero(phi_col1_old) )
      {
         /* set the dual of col1 to zero, compensate for the error made with the dual of row1, compensate that error with the dual of col2
          * and finally compensate the error in col2's row of the reduced costs with col2's dual
          */

         /* dual to zero and compensate error */
         const double dual_shift_row1 = -phi_col1_old / coeff_col1;
         assert( PIPSisEQ( coeff_col1 * dual_shift_row1, -phi_col1_old ) );
         phi_col1 -= phi_col1_old;

         getSimpleVecFromRowStochVec(original_vars.y, row1) += dual_shift_row1;

         /* compensate again and adjust dual of col2 */
         const double dual_shift_row2 = -dual_shift_row1 * parallel_factor;

         if( row2.inEqSys() )
            getSimpleVecFromRowStochVec(original_vars.y, row2) += dual_shift_row2;
         else
         {
            double& z = getSimpleVecFromRowStochVec(original_vars.z, row2);
            double& lambda = getSimpleVecFromRowStochVec(original_vars.lambda, row2);
            double& pi = getSimpleVecFromRowStochVec(original_vars.pi, row2);

            addIneqRowDual(z, lambda, pi, dual_shift_row2);
         }

         if( col2.isCol() )
         {
            /* assert bounds of substituted variable are similarly tight */
            if( PIPSisLT(0.0, scalar) )
               assert( PIPSisEQFeas( xupp_col2 - getSimpleVecFromColStochVec(original_vars.x, col2), std::fabs(old_slack_upper * scalar) ) );
            else
               assert( PIPSisEQFeas( getSimpleVecFromColStochVec(original_vars.x, col2) - xlow_col2, std::fabs(old_slack_upper * scalar) ) );

            const double dual_shift_col2 = -coeff_col2 * dual_shift_row2;
            /* -row^T * dual - gamma + phi */
            PIPSisLT(0.0, scalar) ? assert( PIPSisLE( dual_shift_col2, 0.0 ) ) : assert( PIPSisLE( 0.0, dual_shift_col2 ) );

            if( PIPSisLT(0.0, scalar) )
            {
               assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.phi, col2)) );
               double& phi = local_linking_col2 ? (*phi_changes)[col2.getIndex()] : getSimpleVecFromColStochVec(original_vars.phi, col2);

               phi += std::fabs(dual_shift_col2);
            }
            else
            {
               assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.gamma, col2)) );
               double& gamma = local_linking_col2 ? (*gamma_changes)[col2.getIndex()] : getSimpleVecFromColStochVec(original_vars.gamma, col2);

               gamma += std::fabs(dual_shift_col2);
            }
         }
      }
   }

   if( col2.isCol() )
      assert( complementarySlackVariablesMet(original_vars, col2, postsolve_tol) );

   return true;
}

bool StochPostsolver::postsolveFreeColumnSingletonInequalityRow( sVars& original_vars, int reduction_idx)
{
   assert( reductions.at(reduction_idx) == FREE_COLUMN_SINGLETON_INEQUALITY_ROW );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 4 == next_first_float_val);
   assert(first_int_val + 1 == next_first_int_val);
#endif

   const INDEX& row = indices.at(first_index);
   const INDEX& col = indices.at(first_index + 1);
   assert( row.inInEqSys() );

   const int index_stored_row = int_values.at(first_int_val);
   const INDEX row_stored(ROW, row.getNode(), index_stored_row, row.getLinking(), row.getSystemType());

   const double rhs = float_values.at(first_float_val);
   const double coeff = float_values.at(first_float_val + 1);
   const double xlow = float_values.at(first_float_val + 2);
   const double xupp = float_values.at(first_float_val + 3);

   if( col.isCol() )
      assert( !wasColumnRemoved(col) );
   assert( wasRowRemoved(row) );
   markRowAdded(row);

   double row_value = 0.0;

   if( row.isLinkingRow() )
   {
      assert( PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD) );

      if( col.isCol() && (!col.isLinkingCol() || my_rank == 0) )
         row_value = row_storage.multRowTimesVec( row_stored, dynamic_cast<const StochVector&>(*original_vars.x) );
      else
         row_value = row_storage.multLinkingRowTimesVecWithoutBl0( index_stored_row, dynamic_cast<const StochVector&>(*original_vars.x) );

      PIPS_MPIgetSumInPlace(row_value, MPI_COMM_WORLD);
   }
   else
      row_value = row_storage.multRowTimesVec(row_stored, dynamic_cast<const StochVector&>(*original_vars.x));

   /* duals of row and bounds are zero */
   getSimpleVecFromRowStochVec(original_vars.z, row) = 0;
   getSimpleVecFromRowStochVec(original_vars.lambda, row) = 0;
   getSimpleVecFromRowStochVec(original_vars.pi, row) = 0;

   if( col.isCol() )
   {
      const bool local_linking_col = col.isLinkingCol() && row.getNode() != -1;
      if( local_linking_col )
         outdated_linking_vars = true;

      /* set x value such that row is satisfied */
      assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.x, col)) );
      double& x_val = local_linking_col ? (*x_changes)[col.getIndex()] : getSimpleVecFromColStochVec(original_vars.x, col);
      x_val = (rhs - row_value) / coeff;

      getSimpleVecFromColStochVec(original_vars.gamma, col) = 0;
      getSimpleVecFromColStochVec(original_vars.phi, col) = 0;

      if( xlow == INF_NEG )
         getSimpleVecFromColStochVec(original_vars.v, col) = 0;
      else
      {
         assert( PIPSisLE(xlow, x_val) );
         if( local_linking_col )
            (*v_changes)[col.getIndex()] = x_val - xlow;
         else
            getSimpleVecFromColStochVec(original_vars.v, col) = x_val - xlow;
      }

      if( xupp == INF_POS )
         getSimpleVecFromColStochVec(original_vars.w, col) = 0;
      else
      {
         assert( PIPSisLE(x_val, xupp) );
         if( local_linking_col )
            (*w_changes)[col.getIndex()] = xupp - x_val;
         else
            getSimpleVecFromColStochVec(original_vars.w, col) = xupp - x_val;
      }

      getSimpleVecFromColStochVec(original_vars.gamma, col) = 0.0;
      getSimpleVecFromColStochVec(original_vars.phi, col) = 0.0;

      const double row_slack = row_value + coeff * x_val;

      if( row.isLinkingRow() && !col.isLinkingCol() )
      {
#ifndef NDEBUG
         const double row_slack_recieve = PIPS_MPIgetSum( row_slack, MPI_COMM_WORLD );
         assert( row_slack_recieve == row_slack );
#endif
         /* we just recomputed the slack */
         (*s_changes)[row.getIndex()] = 0;
      }

      getSimpleVecFromRowStochVec(original_vars.s, row) = row_slack;
   }
   else
   {
      assert( row.isLinkingRow() );

      const double row_slack = PIPS_MPIgetSum( 0.0, MPI_COMM_WORLD );
      getSimpleVecFromRowStochVec(original_vars.s, row) = row_slack;

      /* we just recomputed the slack */
      (*s_changes)[row.getIndex()] = 0;
   }

   /* compute slacks for bounds and row */
   /* slack for row is zero by construction */
   getSimpleVecFromRowStochVec(original_vars.t, row) = 0;
   getSimpleVecFromRowStochVec(original_vars.u, row) = 0;

   return true;
}

bool StochPostsolver::postsolveParallelRowsBoundsTightened(sVars& original_vars, int reduction_idx) const
{
   assert( reductions.at(reduction_idx) == PARALLEL_ROWS_BOUNDS_TIGHTENED );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 5 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   const INDEX& row1 = indices.at(first_index);
   const INDEX& row2 = indices.at(first_index + 1);

   assert( row1.isRow() );
   assert( row2.isRow() );
   assert( row1.inInEqSys() );
   assert( row2.inInEqSys() );

   const double clow_old = float_values.at(first_float_val);
   const double cupp_old = float_values.at(first_float_val + 1);
   const double clow_new = float_values.at(first_float_val + 2);
   const double cupp_new = float_values.at(first_float_val + 3);
   const double factor = float_values.at(first_float_val + 4);

   assert (!PIPSisZero(factor) );
   assert( !wasRowRemoved(row1) );
   assert( !wasRowRemoved(row2) );

   bool clow_tightened_by_row2 = PIPSisLT( clow_old, clow_new );
   bool cupp_tightened_by_row2 = PIPSisLT( cupp_new, cupp_old );
   if( clow_old == INF_NEG && clow_new != INF_NEG )
      assert( clow_tightened_by_row2 );
   if( cupp_old == INF_POS && cupp_new != INF_POS )
      assert( cupp_tightened_by_row2 );

   /* recompute duals and slack of both rows - if one bound was tight and thus the dual non-zero we shift it to the row originally implying the bound */
   double& z_row1 = getSimpleVecFromRowStochVec(original_vars.z, row1);
   double& lambda_row1 = getSimpleVecFromRowStochVec(original_vars.lambda, row1);
   double& pi_row1 = getSimpleVecFromRowStochVec(original_vars.pi, row1);

   double& t_row1 = getSimpleVecFromRowStochVec(original_vars.t, row1);
   double& u_row1 = getSimpleVecFromRowStochVec(original_vars.u, row1);

   double& z_row2 = getSimpleVecFromRowStochVec(original_vars.z, row2);
   double& lambda_row2 = getSimpleVecFromRowStochVec(original_vars.lambda, row2);
   double& pi_row2 = getSimpleVecFromRowStochVec(original_vars.pi, row2);

   assert( PIPSisZero(z_row2) );
   assert( PIPSisZero(pi_row2) );
   assert( PIPSisZero(lambda_row2) );

   if( PIPSisLT(factor, 0.0) )
      std::swap(lambda_row2, pi_row2);

   /* if a bound was implied by row2 we have to shift the dual multiplies over to the implying row */
   if( clow_tightened_by_row2 )
   {
      assert( PIPSisLT(clow_old, clow_new) );

      /* adjust slacks */
      if( clow_old == INF_NEG )
         t_row1 = 0.0;
      else
         t_row1 += clow_new - clow_old;

      /* clow is no longer tight on row1 but now on row2 - shift row1's multiplier to row2 */
      z_row1 -= lambda_row1;
      z_row2 += lambda_row1 * factor;

      lambda_row2 = lambda_row1 * factor;
      lambda_row1 = 0.0;
   }

   if( cupp_tightened_by_row2 )
   {
      assert( PIPSisLT(cupp_new, cupp_old) );

      /* adjust slacks */
      if( cupp_old == INF_POS )
         u_row1 = 0.0;
      else
         u_row1 += cupp_old - cupp_new;

      /* cupp is no longer tight on row1 but now on row2 - shift row1's multiplier to row2 */
      z_row1 += pi_row1;
      z_row2 -= pi_row1 * factor;

      pi_row2 = pi_row1 * factor;
      pi_row1 = 0.0;
   }

   /* check slacks - these might have increased by the parallelity factor */
   assert( complementarySlackRowMet(original_vars, row2, postsolve_tol * std::fabs(factor) * 10) );
   assert( complementarySlackRowMet(original_vars, row1, postsolve_tol) );
   return true;
}

/* sync linking variables - either the variables are not set everywhere or on some procs they are set to something non-zero while on others they are zero */
bool StochPostsolver::syncLinkingVarChanges(sVars& original_vars)
{
   assert( dynamic_cast<const StochVector&>(*original_vars.x).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.v).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.w).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.gamma).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.phi).isRootNodeInSync() );

   PIPS_MPIgetLogicOrInPlace(outdated_linking_vars, MPI_COMM_WORLD);

   if( !outdated_linking_vars || x_changes->length() == 0)
      return true;

   PIPS_MPIsumArrayInPlace( array_linking_var_changes, length_array_linking_var_changes, MPI_COMM_WORLD );
   PIPS_MPImaxArrayInPlace( dynamic_cast<SimpleVectorBase<int>*>(padding_origcol->vec)->elements(), dynamic_cast<SimpleVectorBase<int>*>(padding_origcol->vec)->length(), MPI_COMM_WORLD );

   SimpleVector& linking_x = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.x).vec);
   SimpleVector& linking_v = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.v).vec);
   SimpleVector& linking_w = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.w).vec);
   SimpleVector& linking_gamma = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.gamma).vec);
   SimpleVector& linking_phi = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.phi).vec);

   for( int i = 0; i < linking_x.length(); ++i )
   {
      if( wasColumnRemoved( INDEX(COL, -1, i) ) )
         continue;

      if( std::isnan(linking_x[i]) )
      {
         linking_x[i] = (*x_changes)[i];
         linking_v[i] = (*v_changes)[i];
         linking_w[i] = (*w_changes)[i];
         linking_gamma[i] = (*gamma_changes)[i];
         linking_phi[i] = (*phi_changes)[i];
      }
      else
      {
         assert( !std::isnan(linking_v[i]) );
         assert( !std::isnan(linking_w[i]) );
         assert( !std::isnan(linking_gamma[i]) );
         assert( !std::isnan(linking_phi[i]) );

         linking_x[i] += (*x_changes)[i];
         linking_v[i] += (*v_changes)[i];
         linking_w[i] += (*w_changes)[i];
         linking_gamma[i] += (*gamma_changes)[i];
         linking_phi[i] += (*phi_changes)[i];
      }
   }

   std::memset( array_linking_var_changes, 0, length_array_linking_var_changes );
   outdated_linking_vars = false;

   assert( dynamic_cast<const StochVector&>(*original_vars.x).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.v).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.w).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.gamma).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.phi).isRootNodeInSync() );

   return true;
}

bool StochPostsolver::syncEqLinkingRowChanges(sVars& original_vars)
{
   PIPS_MPIgetLogicOrInPlace(outdated_equality_linking_rows, MPI_COMM_WORLD);

   if( !outdated_equality_linking_rows || y_changes->length() == 0)
      return true;

   PIPS_MPIsumArrayInPlace( array_eq_linking_row_changes, length_array_eq_linking_row_changes, MPI_COMM_WORLD );
   PIPS_MPImaxArrayInPlace( dynamic_cast<SimpleVectorBase<int>*>(padding_origrow_equality->vecl)->elements(),
         dynamic_cast<SimpleVectorBase<int>*>(padding_origrow_equality->vecl)->length(), MPI_COMM_WORLD );

   SimpleVector& linking_y = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.y).vecl);

   for( int i = 0; i < linking_y.length(); ++i )
   {
      if( wasRowRemoved( INDEX(ROW, -1, i, true, EQUALITY_SYSTEM) ) )
         continue;

      if( !wasRowRemoved( INDEX(ROW, -1, i, true, EQUALITY_SYSTEM) ) && std::isnan(linking_y[i]) )
         linking_y[i] = (*y_changes)[i];
      else
         linking_y[i] += (*y_changes)[i];
   }

   std::memset( array_eq_linking_row_changes, 0, length_array_eq_linking_row_changes );
   outdated_equality_linking_rows = false;

   assert( dynamic_cast<const StochVector&>(*original_vars.y).isRootNodeInSync() );

   return true;
}

bool StochPostsolver::syncIneqLinkingRowChanges(sVars& original_vars)
{
   assert( dynamic_cast<const StochVector&>(*original_vars.z).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.s).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.t).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.u).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.pi).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.lambda).isRootNodeInSync() );

   PIPS_MPIgetLogicOrInPlace(outdated_inequality_linking_rows, MPI_COMM_WORLD);

   if( !outdated_inequality_linking_rows || z_changes->length() == 0)
      return true;

   PIPS_MPIsumArrayInPlace( array_ineq_linking_row_changes, length_array_ineq_linking_row_changes, MPI_COMM_WORLD );
   PIPS_MPImaxArrayInPlace( dynamic_cast<SimpleVectorBase<int>*>(padding_origrow_inequality->vecl)->elements(),
         dynamic_cast<SimpleVectorBase<int>*>(padding_origrow_inequality->vecl)->length(), MPI_COMM_WORLD );

   SimpleVector& linking_z = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.z).vecl);
   SimpleVector& linking_lambda = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.lambda).vecl);
   SimpleVector& linking_pi = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.pi).vecl);

   SimpleVector& linking_s = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.s).vecl);
   SimpleVector& linking_t = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.t).vecl);
   SimpleVector& linking_u = dynamic_cast<SimpleVector&>(*dynamic_cast<StochVector&>(*original_vars.u).vecl);

   for( int i = 0; i < linking_z.length(); ++i )
   {
      if( wasRowRemoved( INDEX(ROW, -1, i, true, INEQUALITY_SYSTEM) ) )
         continue;

      if( std::isnan(linking_z[i]) )
      {
         linking_z[i] = (*z_changes)[i];
         linking_lambda[i] = std::max(0.0, linking_z[i]);
         linking_pi[i] = -std::min(0.0, linking_z[i]);
         linking_s[i] = (*s_changes)[i];
         linking_t[i] = (*t_changes)[i];
         linking_u[i] = (*u_changes)[i];
      }
      else
      {
         assert( !std::isnan(linking_z[i]) );
         assert( !std::isnan(linking_lambda[i]) );
         assert( !std::isnan(linking_pi[i]) );
         assert( !std::isnan(linking_s[i]) );
         assert( !std::isnan(linking_t[i]) );
         assert( !std::isnan(linking_u[i]) );

         addIneqRowDual(linking_z[i], linking_lambda[i], linking_pi[i], (*z_changes)[i]);
         linking_s[i] += (*s_changes)[i];
         linking_t[i] += (*t_changes)[i];
         linking_u[i] += (*u_changes)[i];
      }
   }

   std::memset( array_ineq_linking_row_changes, 0, length_array_ineq_linking_row_changes );
   outdated_inequality_linking_rows = false;

   assert( dynamic_cast<const StochVector&>(*original_vars.z).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.s).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.t).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.u).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.pi).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.lambda).isRootNodeInSync() );

   return true;
}

bool StochPostsolver::syncLinkingRowsAfterBoundTightening(sVars& original_vars, int i)
{
   assert( dynamic_cast<const StochVector&>(*original_vars.y).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.z).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.lambda).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.pi).isRootNodeInSync() );

   PIPS_MPIgetLogicOrInPlace(outdated_equality_linking_rows);
   PIPS_MPIgetLogicOrInPlace(outdated_inequality_linking_rows);

   if( !outdated_equality_linking_rows && !outdated_inequality_linking_rows )
      return true;

   /* find STORE_BOUND_TIGHTENING_LINKING_ROWS event - there all linking rows have been stored (after bound tightening so it is actually down the stack and has already been processed */
   while( reductions[i] != STORE_BOUND_TIGHTENING_LINKING_ROWS )
   {
      ++i;
      assert( static_cast<unsigned int>(i) < reductions.size() );
   }

   unsigned int current_pos = start_idx_indices.at(i);
   StochVector& gamma = dynamic_cast<StochVector&>(*original_vars.gamma);
   StochVector& phi = dynamic_cast<StochVector&>(*original_vars.phi);

   /* gather all z and y changes if any */
   if( outdated_equality_linking_rows )
   {
      PIPS_MPIsumArrayInPlace(array_eq_linking_row_changes, length_array_eq_linking_row_changes);

      for( int row = 0; row < y_changes->length(); ++row )
      {
         const double change_dual_row = (*y_changes)[row];
         (*y_changes)[row] = 0.0;
         if( !PIPSisZero( change_dual_row ) )
         {
            const INDEX row_INDEX(ROW, -1, row, true, EQUALITY_SYSTEM);
            const int row_stored = findNextRowInStored(i, current_pos, row_INDEX);

            const INDEX stored_row(ROW, -1, row_stored, true, EQUALITY_SYSTEM);
            /// adjust duals of rows and bounds
            row_storage.axpyAtRowPosNeg(1.0, &phi, nullptr, &gamma, nullptr, change_dual_row, stored_row );
            getSimpleVecFromRowStochVec(original_vars.y, row_INDEX) += change_dual_row;
         }
      }
   }

   if( outdated_inequality_linking_rows )
   {
      PIPS_MPIsumArrayInPlace(array_ineq_linking_row_changes, length_array_ineq_linking_row_changes);

      for( int row = 0; row < z_changes->length(); ++row )
      {
         const double change_dual_row = (*z_changes)[row];
         (*z_changes)[row] = 0.0;

         if( !PIPSisZero( change_dual_row ) )
         {
            const INDEX row_INDEX(ROW, -1, row, true, INEQUALITY_SYSTEM);
            const int row_stored = findNextRowInStored(i, current_pos, row_INDEX);

            const INDEX stored_row(ROW, -1, row_stored, true, EQUALITY_SYSTEM);
            /// adjust duals of rows and bounds
            row_storage.axpyAtRowPosNeg(1.0, &phi, nullptr, &gamma, nullptr, change_dual_row, stored_row );

            /* adjust the row dual */
            double& z = getSimpleVecFromRowStochVec(original_vars.z, row_INDEX);
            double& lambda = getSimpleVecFromRowStochVec(original_vars.lambda, row_INDEX);
            double& pi = getSimpleVecFromRowStochVec(original_vars.pi, row_INDEX);

            addIneqRowDual(z, lambda, pi, change_dual_row);
         }
      }
   }
   assert( dynamic_cast<const StochVector&>(*original_vars.y).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.z).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.lambda).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.pi).isRootNodeInSync() );

   outdated_equality_linking_rows = false;
   outdated_inequality_linking_rows = false;

   return true;
}

int StochPostsolver::findNextRowInStored(int pos_reduction, unsigned int& start, const INDEX& row) const
{
   assert( reductions[pos_reduction] == STORE_BOUND_TIGHTENING_LINKING_ROWS );
   assert( start_idx_indices.at(pos_reduction) <= start );
   assert( start_idx_indices.at(pos_reduction + 1) > start );
   assert( start_idx_indices[pos_reduction + 1] - start_idx_indices[pos_reduction] ==
         start_idx_int_values[pos_reduction + 1] - start_idx_int_values[pos_reduction]);

   while( indices[start] != row )
   {
      start++;
      assert( start < start_idx_indices.at(pos_reduction + 1) );
   }

   assert(indices[start] == row);

   const int dist = start - start_idx_indices[pos_reduction];
   const int pos_stored_index = start_idx_int_values.at(pos_reduction) + dist;

   return int_values[pos_stored_index];
}

void StochPostsolver::addIneqRowDual(double& z, double& lambda, double& pi, double value) const
{
   z += value;
   /* z = lambda - pi */
   if( value > 0 )
   {
      if( pi > 0 )
      {
         pi -= value;
         if( pi < 0 )
         {
            value = -pi;
            pi = 0.0;
         }
         else
            value = 0.0;
      }

      lambda += value;
   }
   /* value < 0 */
   else
   {
      if( lambda > 0 )
      {
         lambda += value;
         if( lambda < 0 )
         {
            value = lambda;
            lambda = 0.0;
         }
         else
            value = 0.0;
      }

      pi += -value;
   }

   assert( PIPSisEQ(z, lambda - pi, postsolve_tol) );
}

void StochPostsolver::addIneqRowSlack(double& s, double& t, double& u, double clow, double cupp, double value) const
{
   assert( PIPSisZero( s - clow - t, postsolve_tol) );
   assert( PIPSisZero( s - cupp + u, postsolve_tol) );
   assert( PIPSisLT( 0.0, t, postsolve_tol) );
   assert( PIPSisLT( 0.0, u, postsolve_tol) );

   s += value;
   t += value;
   u += value;

   assert( PIPSisLT( 0.0, t, postsolve_tol) );
   assert( PIPSisLT( 0.0, u, postsolve_tol) );
   assert( PIPSisZero( s - clow - t, postsolve_tol) );
   assert( PIPSisZero( s - cupp + u, postsolve_tol) );
}

bool StochPostsolver::sameNonZeroPatternDistributed(const StochVector& vec) const
{
   assert(vec.vec);

   if( !sameNonZeroPatternDistributed( dynamic_cast<const SimpleVector&>(*vec.vec) ) )
      return false;

   if( vec.vecl )
   {
      if( !sameNonZeroPatternDistributed( dynamic_cast<const SimpleVector&>(*vec.vecl) ) )
         return false;
   }
   return true;
}

bool StochPostsolver::sameNonZeroPatternDistributed(const SimpleVector& vec) const
{
   std::vector<double> v(vec.elements(), vec.elements() + vec.length() );

   bool result = true;

   for( unsigned int i = 0; i < v.size(); ++i )
   {
      if( !PIPSisZero(v[i]) )
         v[i] = 1.0;
   }

   const std::vector<double> ref_vec(v.begin(), v.end());

   PIPS_MPIminArrayInPlace(v, MPI_COMM_WORLD);

   for( unsigned i = 0; i < v.size(); ++i )
   {
      if( v[i] != ref_vec[i] )
      {
         result = false;
      }
   }

   return PIPS_MPIgetLogicAnd( result, MPI_COMM_WORLD );
}


void StochPostsolver::setOriginalVarsFromReduced(const sVars& reduced_vars, sVars& original_vars) const
{
#ifdef ANCIENT_CPP
   const double initial_const = NAN;
#else
   const double initial_const = std::nan("not set");
#endif

   /* x */
   const StochVector &x_reduced = dynamic_cast<const StochVector&>(*reduced_vars.x);
   StochVector &x_orig = dynamic_cast<StochVector&>(*original_vars.x);
   x_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced<>(x_orig, x_reduced, *padding_origcol);

   /* s */
   const StochVector &s_reduced = dynamic_cast<const StochVector&>(*reduced_vars.s);
   StochVector &s_orig = dynamic_cast<StochVector&>(*original_vars.s);
   s_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(s_orig, s_reduced, *padding_origrow_inequality);

   /* y */
   const StochVector &y_reduced = dynamic_cast<const StochVector&>(*reduced_vars.y);
   StochVector &y_orig = dynamic_cast<StochVector&>(*original_vars.y);
   y_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(y_orig, y_reduced, *padding_origrow_equality);

   /* z */
   const StochVector &z_reduced = dynamic_cast<const StochVector&>(*reduced_vars.z);
   StochVector &z_orig = dynamic_cast<StochVector&>(*original_vars.z);
   z_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(z_orig, z_reduced, *padding_origrow_inequality);

   /* v */
   const StochVector &v_reduced = dynamic_cast<const StochVector&>(*reduced_vars.v);
   StochVector &v_orig = dynamic_cast<StochVector&>(*original_vars.v);
   v_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(v_orig, v_reduced, *padding_origcol);

   /* gamma */
   const StochVector &gamma_reduced = dynamic_cast<const StochVector&>(*reduced_vars.gamma);
   StochVector &gamma_orig = dynamic_cast<StochVector&>(*original_vars.gamma);
   gamma_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(gamma_orig, gamma_reduced, *padding_origcol);

   /* w */
   const StochVector &w_reduced = dynamic_cast<const StochVector&>(*reduced_vars.w);
   StochVector &w_orig = dynamic_cast<StochVector&>(*original_vars.w);
   w_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(w_orig, w_reduced, *padding_origcol);

   /* phi */
   const StochVector &phi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.phi);
   StochVector &phi_orig = dynamic_cast<StochVector&>(*original_vars.phi);
   phi_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(phi_orig, phi_reduced, *padding_origcol);

   /* t */
   const StochVector &t_reduced = dynamic_cast<const StochVector&>(*reduced_vars.t);
   StochVector &t_orig = dynamic_cast<StochVector&>(*original_vars.t);
   t_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(t_orig, t_reduced, *padding_origrow_inequality);

   /* lambda */
   const StochVector &lambda_reduced = dynamic_cast<const StochVector&>(*reduced_vars.lambda);
   StochVector &lambda_orig = dynamic_cast<StochVector&>(*original_vars.lambda);
   lambda_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(lambda_orig, lambda_reduced, *padding_origrow_inequality);

   /* u */
   const StochVector &u_reduced = dynamic_cast<const StochVector&>(*reduced_vars.u);
   StochVector &u_orig = dynamic_cast<StochVector&>(*original_vars.u);
   u_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(u_orig, u_reduced, *padding_origrow_inequality);

   /* pi */
   const StochVector &pi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.pi);
   StochVector &pi_orig = dynamic_cast<StochVector&>(*original_vars.pi);
   pi_orig.setToConstant(initial_const);
   setOriginalValuesFromReduced(pi_orig, pi_reduced, *padding_origrow_inequality);
}

bool StochPostsolver::allVariablesSet(const sVars& vars) const
{
   bool all_set = true;
#ifdef ANCIENT_CPP
   const double initial_const = NAN;
#else
   const double initial_const = std::nan("not set"); // TODO not working
#endif

   /* x */
   if( !vars.x->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "x not all set" << std::endl;
      all_set = false;
   }

   /* s */
   if( !vars.s->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "s not all set" << std::endl;
      all_set = false;
   }

   /* y */
   if( !vars.y->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "y not all set" << std::endl;
      all_set = false;
   }

   /* z */
   if( !vars.z->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "z not all set" << std::endl;
      all_set = false;
   }

   /* v */
   if( !vars.v->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "v not all set" << std::endl;
      all_set = false;
   }

   /* gamma */
   if( !vars.gamma->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "gamma not all set" << std::endl;
      all_set = false;
   }

   /* w */
   if( !vars.w->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "w not all set" << std::endl;
      all_set = false;
   }

   /* phi */
   if( !vars.phi->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "phi not all set" << std::endl;
      all_set = false;
   }

   /* t */
   if( !vars.t->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "t not all set" << std::endl;
      all_set = false;
   }

   /* lambda */
   if( !vars.phi->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "phi not all set" << std::endl;
      all_set = false;
   }

   /* u */
   if( !vars.u->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "u not all set" << std::endl;
      all_set = false;
   }

   /* pi */
   if( !vars.pi->componentNotEqual(initial_const) )
   {
      if( my_rank == 0 )
         std::cout << "pi not all set" << std::endl;
      all_set = false;
   }

   return all_set;
}

bool StochPostsolver::complementarySlackVariablesMet(const sVars& vars, const INDEX& col, double tol) const
{
   assert( col.isCol() );
   assert( !wasColumnRemoved(col) );
   assert( tol > 0 );

   const int index = col.getIndex();

   const double v = col.isLinkingCol() ? getSimpleVecFromColStochVec(vars.v, col) + (*v_changes)[index] : getSimpleVecFromColStochVec(vars.v, col);
   const double w = col.isLinkingCol() ? getSimpleVecFromColStochVec(vars.w, col) + (*w_changes)[index] : getSimpleVecFromColStochVec(vars.w, col);
   const double gamma = col.isLinkingCol() ? getSimpleVecFromColStochVec(vars.gamma, col) + (*gamma_changes)[index] : getSimpleVecFromColStochVec(vars.gamma, col);
   const double phi = col.isLinkingCol() ? getSimpleVecFromColStochVec(vars.phi, col) + (*phi_changes)[index] : getSimpleVecFromColStochVec(vars.phi, col);

   assert(!std::isnan(v));
   assert(!std::isnan(w));
   assert(!std::isnan(gamma));
   assert(!std::isnan(phi));

   if( std::fabs(v * gamma) >= tol )
   {
      std::cout << "rv " << v << " " << gamma << " " << v * gamma << " vs " << tol << std::endl;
      return true;
   }
   if( std::fabs(w * phi) >= tol )
   {
      std::cout << "rw " << w << " " << phi << " " << w * phi << " vs " << tol << std::endl;
      return true;
   }
   return std::fabs(v * gamma) < tol && std::fabs(w * phi) < tol ;
}

bool StochPostsolver::complementarySlackRowMet(const sVars& vars, const INDEX& row, double tol) const
{
   assert( row.isRow() );
   assert( !wasRowRemoved(row) );
   assert( tol > 0 );

   const double t = getSimpleVecFromRowStochVec(vars.t, row);
   const double u = getSimpleVecFromRowStochVec(vars.u, row);
   const double lambda = getSimpleVecFromRowStochVec(vars.lambda, row);
   const double pi = getSimpleVecFromRowStochVec(vars.pi, row);

   assert(!std::isnan(t));
   assert(!std::isnan(u));
   assert(!std::isnan(lambda));
   assert(!std::isnan(pi));

   if( std::fabs(t * lambda) >= tol )
   {
      std::cout << "rt " << t << " " << lambda << " " << t * lambda << " vs " << tol << std::endl;
      return true;
   }
   if( !(std::fabs(t * lambda) < tol) )
   {
      std::cout << t * lambda << " " << std::fabs(t * lambda) << " " << tol << std::endl;
      std::cout << "|t * lambda| >= tol => false: " << (std::fabs(t * lambda) >= tol) << std::endl;
      std::cout << "|t * lambda| < tol => false??: " << (std::fabs(t * lambda) < tol) << std::endl;
   }
   assert( std::fabs(t * lambda) < tol );
   if( std::fabs(u * pi) >= tol )
   {
      std::cout << "ru " << u << " " << pi << " " << u * pi << " vs " << tol << std::endl;
      return true;
   }
   assert( std::fabs( u * pi) < tol );

   return std::fabs(t * lambda) < tol && std::fabs(u * pi) < tol;
}

/// fills vars_orig with vars_reduced padded with zeros - padding is done via the padding_map
template <typename T>
void StochPostsolver::setOriginalValuesFromReduced(StochVectorBase<T>& original_vector,
   const StochVectorBase<T>& reduced_vector,
   const StochVectorBase<int>& padding_original) const
{
   assert( reduced_vector.children.size() == original_vector.children.size() );
   assert( padding_original.children.size() == reduced_vector.children.size() );
   assert( reduced_vector.vec != nullptr && original_vector.vec != nullptr && padding_original.vec != nullptr );
   assert( (reduced_vector.vecl != nullptr && original_vector.vecl != nullptr && padding_original.vecl != nullptr)
         || (reduced_vector.vecl == nullptr && original_vector.vecl == nullptr && padding_original.vecl == nullptr) );

   if( reduced_vector.isKindOf(kStochDummy) )
   {
      assert( original_vector.isKindOf(kStochDummy) && padding_original.isKindOf(kStochDummy) );
      return;
   }

   /* root node */
   /* vec */
   setOriginalValuesFromReduced( dynamic_cast<SimpleVectorBase<T>&>(*original_vector.vec), dynamic_cast<const SimpleVectorBase<T>&>(*reduced_vector.vec),
      dynamic_cast<const SimpleVectorBase<int>&>(*padding_original.vec));

   /* vecl */
   if( reduced_vector.vecl )
   {
      setOriginalValuesFromReduced( dynamic_cast<SimpleVectorBase<T>&>(*original_vector.vecl), dynamic_cast<const SimpleVectorBase<T>&>(*reduced_vector.vecl),
         dynamic_cast<const SimpleVectorBase<int>&>(*padding_original.vecl));
   }

   /* child nodes */
   for( int i = 0; i < static_cast<int>(reduced_vector.children.size()); ++i )
   {
      setOriginalValuesFromReduced( *original_vector.children[i],  *reduced_vector.children[i], *padding_original.children[i] );
   }
}

template <typename T>
void StochPostsolver::setOriginalValuesFromReduced(SimpleVectorBase<T>& original_vector,
   const SimpleVectorBase<T>& reduced_vector,
   const SimpleVectorBase<int>& padding_original) const
{
   assert( original_vector.length() == padding_original.length() );

   int col_reduced = 0; 
   for(int i = 0; i < padding_original.length(); ++i)
   {
      if(padding_original[i] == -1)
      {
         continue;
      }
      else
      {
         assert(padding_original[i] == 1);
         original_vector[i] = reduced_vector[col_reduced];
         ++col_reduced;
      }
   }

   /* assert all entries are set */
   assert(col_reduced == reduced_vector.length());
}
