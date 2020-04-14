/*
 * StochPostsolver.C
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */

#include "OoqpVector.h"
#include "StochPostsolver.h"
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
   last_lower_bound_tightened( dynamic_cast<StochVectorBase<int>*>(col_stored_last_at->clone()) )
{
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

   // count !? rows cols...
   start_idx_float_values.push_back(0);
   start_idx_int_values.push_back(0);
   start_idx_indices.push_back(0);
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
}

void StochPostsolver::notifyRowModified( const INDEX& row )
{
   assert(row.isRow());
   if(row.getSystemType() == EQUALITY_SYSTEM)
      getSimpleVecFromRowStochVec(*eq_row_marked_modified, row.getNode(), row.getLinking())[row.getIndex()] = 1;
   else
      getSimpleVecFromRowStochVec(*ineq_row_marked_modified, row.getNode(), row.getLinking())[row.getIndex()] = 1;
}

void StochPostsolver::notifyColModified( const INDEX& col )
{
   assert(col.isCol());
   getSimpleVecFromColStochVec(*column_marked_modified, col.getNode())[col.getIndex()] = 1;
}

bool StochPostsolver::isRowModified( const INDEX& row) const
{
   assert(row.isRow());
   if(row.getSystemType() == EQUALITY_SYSTEM)
      return getSimpleVecFromRowStochVec(*eq_row_marked_modified, row.getNode(), row.getLinking())[row.getIndex()] == 1;
   else
      return getSimpleVecFromRowStochVec(*ineq_row_marked_modified, row.getNode(), row.getLinking())[row.getIndex()] == 1;
}

void StochPostsolver::markRowClean( const INDEX& row )
{
   assert(row.isRow());
   if(row.getSystemType() == EQUALITY_SYSTEM)
      getSimpleVecFromRowStochVec(*eq_row_marked_modified, row.getNode(), row.getLinking())[row.getIndex()] = -1;
   else
      getSimpleVecFromRowStochVec(*ineq_row_marked_modified, row.getNode(), row.getLinking())[row.getIndex()] = -1;
}

void StochPostsolver::markColClean( const INDEX& col )
{
   assert(col.isCol());
   getSimpleVecFromColStochVec(*col_stored_last_at, col.getNode())[col.getIndex()] = -1;
}

bool StochPostsolver::isColModified( const INDEX& col) const
{
   assert(col.isCol());
   return getSimpleVecFromColStochVec(*column_marked_modified, col.getNode())[col.getIndex()] == 1;
}

int StochPostsolver::storeRow( const INDEX& row, const StochGenMatrix& matrix_row)
{
   assert(row.isRow());
   const int node = row.getNode();
   const bool linking_row = row.getLinking();
   const int row_index = row.getIndex();
   const SystemType system_type = row.getSystemType();

   if( isRowModified(row) )
   {
      markRowClean(row);

      const int stored_at = row_storage.storeRow(row, matrix_row);
      if( system_type == EQUALITY_SYSTEM )
         getSimpleVecFromRowStochVec(*eq_row_stored_last_at, node, linking_row)[row_index] = stored_at;
      else
         getSimpleVecFromRowStochVec(*ineq_row_stored_last_at, node, linking_row)[row_index] = stored_at;
      return stored_at;
   }
   else
   {
      if(system_type == EQUALITY_SYSTEM)
      {
         assert(getSimpleVecFromRowStochVec(*eq_row_stored_last_at, node, linking_row)[row_index] != -1);
         return getSimpleVecFromRowStochVec(*eq_row_stored_last_at, node, linking_row)[row_index];
      }
      else
      {
         assert(getSimpleVecFromRowStochVec(*ineq_row_stored_last_at, node, linking_row)[row_index] != -1);
         return getSimpleVecFromRowStochVec(*ineq_row_stored_last_at, node, linking_row)[row_index];
      }
   }
}

int StochPostsolver::storeColumn( const INDEX& col, const StochGenMatrix& matrix_col_eq, const StochGenMatrix& matrix_col_ineq)
{
   assert(col.isCol());
   const int node = col.getNode();
   const int col_index = col.getIndex();

   if( isColModified(col) )
   {
      markColClean(col);
      const int stored_at = col_storage.storeCol(col, matrix_col_eq, matrix_col_ineq);

      getSimpleVecFromColStochVec(*col_stored_last_at, node)[col_index] = stored_at;

      return stored_at;
   }
   else
   {
      assert(getSimpleVecFromColStochVec(*col_stored_last_at, node)[col_index] != -1);
      return getSimpleVecFromColStochVec(*col_stored_last_at, node)[col_index];
   }
}

void StochPostsolver::notifyFixedSingletonFromInequalityColumn( const INDEX& col, double value, double coeff, double xlow_old, double xupp_old )
{
   assert(col.isCol());
   assert(!wasColumnRemoved(col));

   markColumnRemoved(col);

   reductions.push_back(FIXED_COLUMN_SINGLETON_FROM_INEQUALITY);

   indices.push_back(col);

   float_values.push_back(value);
   float_values.push_back(coeff);
   float_values.push_back(xlow_old);
   float_values.push_back(xupp_old);

   finishNotify();
}

void StochPostsolver::notifyFreeColumnSingletonInequalityRow( const INDEX& row, const INDEX& col, double rhs, double coeff, double xlow, double xupp, const StochGenMatrix& matrix_row )
{
   assert(row.isRow());
   assert(col.isCol() || col.isEmpty() );
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

void StochPostsolver::putLinkingVarsSyncEvent()
{
   reductions.push_back(LINKING_VARS_SYNC_EVENT);
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

/** substitute col2 with scalar*col1 + translation
 *
 * this happens for:
 *    two equality rows - coeff_col2 != 0
 *    two inequality rows (lot's of prerequisites)
 *
 */
void StochPostsolver::notifyNearlyParallelRowSubstitution(const INDEX& row1, const INDEX& row2, const INDEX& col1, const INDEX& col2, double scalar, double translation,
   double obj_col2, double xlow_col2, double xupp_col2, double coeff_col1, double coeff_col2, double parallel_factor )
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

   markColumnRemoved(col2);

   float_values.push_back( scalar );
   float_values.push_back( translation );
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
      assert( xlow_col2 == INF_NEG_PRES );
      assert( xupp_col2 == INF_POS_PRES );
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
   assert( xupp_new != INF_POS_PRES || xlow_new != INF_NEG_PRES );
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
   assert(!wasColumnRemoved(col));
   markColumnRemoved(col);

   /* store current upper and lower bounds of x and the local column */
   const int col_index = col_storage.storeCol(col, eq_mat, ineq_mat);

   reductions.push_back(FIXED_COLUMN);

   indices.push_back(col);

   int_values.push_back(col_index);

   float_values.push_back(value);
   float_values.push_back(obj_coeff);

   finishNotify();
}

void StochPostsolver::notifyFixedEmptyColumn( const INDEX& col, double value, double obj_coeff, int ixlow, int ixupp, double lbx, double ubx)
{
   assert(col.isCol());
   assert( !wasColumnRemoved(col) );
   if( col.isLinkingCol() )
      assert( PIPS_MPIisValueEqual(col.getNode(), MPI_COMM_WORLD) );

   markColumnRemoved(col);

   if(ixlow == 1)
      assert(PIPSisLEFeas(lbx, value));
   if(ixupp == 1)
      assert(PIPSisLEFeas(value, ubx));

   reductions.push_back(FIXED_EMPTY_COLUMN);
   indices.push_back(col);
   float_values.push_back(value);
   float_values.push_back(obj_coeff);
   float_values.push_back(lbx);
   float_values.push_back(ubx);
   int_values.push_back(ixlow);
   int_values.push_back(ixupp);

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
   reductions.push_back(BOUND_TIGHTENING_START);
   finishNotify();
}

void StochPostsolver::endBoundTightening()
{
   // TODO : assert bound tightening started
   reductions.push_back(BOUND_TIGHTENING_END);
   finishNotify();
}

void StochPostsolver::notifyRowPropagatedBound( const INDEX& row, const INDEX& col, int old_ixlowupp, double old_bound, double new_bound, bool is_upper_bound, const StochGenMatrix& matrix_row)
{
   assert(!PIPSisEQ(old_bound, new_bound));
   assert(row.isRow());
   assert(col.isCol());

   if(is_upper_bound)
      assert( PIPSisLT(new_bound, old_bound) );
   else
      assert( PIPSisLT(old_bound, new_bound) );

   int& index_last = is_upper_bound ? getSimpleVecFromColStochVec(*last_upper_bound_tightened, col) : getSimpleVecFromColStochVec(*last_lower_bound_tightened, col);
   if( index_last != -1 )
   {
      reductions.at(index_last) = DELETED;
      /* copy old old_bounds */
      old_ixlowupp = int_values.at(start_idx_int_values.at(index_last));
      assert(int_values.at(start_idx_int_values.at(index_last) + 1) == is_upper_bound);
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

   int_values.push_back(old_ixlowupp);
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
   assert(false);
   assert(col.isCol());
   assert(wasColumnRemoved(col));
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
   assert(false);
   assert(row.isRow());
   assert(wasRowRemoved(row));
   if( row.inEqSys() )
      getSimpleVecFromRowStochVec(*padding_origrow_equality, row) = 1;
   else
      getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) = 1;
}

// todo : usage and check of padding origrow - can already be done - even without any dual postsolve stuff
// todo : sort reductions by nodes ? and then reverse ?
PostsolveStatus StochPostsolver::postsolve(const Variables& reduced_solution, Variables& original_solution) const
{
   if(my_rank == 0)
      std::cout << "start postsolving... " << std::endl;

   const sVars& stoch_reduced_sol = dynamic_cast<const sVars&>(reduced_solution);
   sVars& stoch_original_sol = dynamic_cast<sVars&>(original_solution);

   /* original variables are now reduced vars padded with zeros */
   setOriginalVarsFromReduced(stoch_reduced_sol, stoch_original_sol);

   const sData& orig_problem_s = dynamic_cast<const sData&>(original_problem);

   /// stuff for synchronization in-between processes
   StochVectorHandle phi_vec_bound_tightening( dynamic_cast<StochVector*>(orig_problem_s.g->clone()) );
   StochVectorHandle gamma_vec_bound_tightening( dynamic_cast<StochVector*>(phi_vec_bound_tightening->clone()) );

   phi_vec_bound_tightening->setToZero();
   gamma_vec_bound_tightening->setToZero();


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
         postsolve_success = postsolve_success && postsolveRedundantRow(stoch_original_sol, i);
         break;
      }
      case BOUND_TIGHTENING_START:
      {
         // TODO ?
         break;
      }
      case BOUND_TIGHTENING_END:
      {
         postsolve_success = postsolve_success && syncBoundsTightened(stoch_original_sol, i, *phi_vec_bound_tightening, *gamma_vec_bound_tightening);
         break;
      }
      case BOUNDS_TIGHTENED:
      {
         postsolve_success = postsolve_success && postsolveBoundsTightened(stoch_original_sol, i, *phi_vec_bound_tightening, *gamma_vec_bound_tightening);
         break;
      }
      case FIXED_COLUMN:
      {
         postsolve_success = postsolve_success && postsolveFixedColumn(stoch_original_sol, i);
         break;
      }
      case FIXED_EMPTY_COLUMN:
      {
         postsolve_success = postsolve_success && postsolveFixedEmptyColumn(stoch_original_sol, i);
         break;
      }
      case FIXED_COLUMN_SINGLETON_FROM_INEQUALITY:
      {
         postsolve_success = postsolve_success && postsolveFixedColumnSingletonFromInequality(stoch_original_sol, i);
         break;
      }
      case SINGLETON_EQUALITY_ROW:
      {
         postsolve_success = postsolve_success && postsolveSingletonEqualityRow(stoch_original_sol, i);
         break;
      }
      case SINGLETON_INEQUALITY_ROW:
      {
         postsolve_success = postsolve_success && postsolveSingletonInequalityRow(stoch_original_sol, i);
         break;
      }
      case FREE_COLUMN_SINGLETON_EQUALITY:
      {
         postsolve_success = postsolve_success && postsolveFreeColumnSingletonEquality(stoch_original_sol, i);
         break;
      }
      case NEARLY_PARALLEL_ROW_SUBSTITUTION:
      {
         postsolve_success = postsolve_success && postsolveNearlyParallelRowSubstitution(stoch_original_sol, i);
         break;
      }
      case NEARLY_PARALLEL_ROW_BOUNDS_TIGHTENED:
      {
         postsolve_success = postsolve_success && postsolveNearlyParallelRowBoundsTightened(stoch_original_sol, i);
         break;
      }
      case FREE_COLUMN_SINGLETON_INEQUALITY_ROW:
      {
         postsolve_success = postsolve_success && postsolveFreeColumnSingletonInequalityRow(stoch_original_sol, orig_problem_s, i);
         break;
      }
      case PARALLEL_ROWS_BOUNDS_TIGHTENED:
      {
         postsolve_success = postsolve_success && postsolveParallelRowsBoundsTightened(stoch_original_sol, i);
         break;
      }
      case LINKING_VARS_SYNC_EVENT:
      {
         postsolve_success = postsolve_success && syncNewlySetLinkingVars(stoch_original_sol);
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

   return PRESOLVE_OK;
}

/**
 * postsolve for a redundant row is to set all dual variables to zero - the row itself has no primal impact
 * for linking rows:
 *    all processes should have removed them in the same order
 *       -> assert in place to check this
 *    the current activity gets synchronized for slack computation
 */
bool StochPostsolver::postsolveRedundantRow(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == REDUNDANT_ROW );

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

   // TODO: this could be optimized - all linking rows will be after each other - we could wait with allreduce --- acutally not sure whether this would ever happen or not - multiple redundant linking rows...
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
   if( row.inEqSys() )
   {
      assert(PIPSisEQ(lhs, rhs));
      if( !PIPSisEQFeas(value_row, rhs) )
         PIPSdebugMessage("Postsolve Warning: when reintroducing a redundant equality row it did not meet its rhs with feastol: %f != %f", value_row, rhs);
      assert( PIPSisEQ(value_row, rhs, postsolve_tol) );

      /* set dual multiplier to zero and mark row as added */
      getSimpleVecFromRowStochVec(*padding_origrow_equality, row) = 1;
      getSimpleVecFromRowStochVec(original_vars.y, row) = 0;
   }
   else
   {
      /* set dual multipliers to zero and mark row as added */
      getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) = 1;

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
   }

   return true;
}

bool StochPostsolver::syncBoundsTightened(sVars& original_vars, int reduction_idx, StochVector& phi, StochVector& gamma) const
{
   SimpleVector& phi_vec = dynamic_cast<SimpleVector&>(*phi.vec);
   SimpleVector& gamma_vec = dynamic_cast<SimpleVector&>(*gamma.vec);
   PIPS_MPIsumArrayInPlace(phi_vec.elements(), phi_vec.length(), MPI_COMM_WORLD);
   PIPS_MPIsumArrayInPlace(gamma_vec.elements(), gamma_vec.length(), MPI_COMM_WORLD);

   original_vars.gamma->axpy(1.0, gamma);
   original_vars.phi->axpy(1.0, phi);

   assert( dynamic_cast<const StochVector&>(*original_vars.gamma).isRootNodeInSync() );
   assert( dynamic_cast<const StochVector&>(*original_vars.phi).isRootNodeInSync() );

   gamma.setToZero();
   phi.setToZero();

   return true;
}

bool StochPostsolver::postsolveBoundsTightened(sVars& original_vars, int reduction_idx, StochVector& phi, StochVector& gamma) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == BOUNDS_TIGHTENED );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 2 == next_first_float_val);
   assert(first_int_val + 3 == next_first_int_val);
#endif

   const INDEX& row = indices.at(first_index);
   const INDEX& col = indices.at(first_index + 1);
   assert(row.isRow());
   assert(col.isCol());
   assert(!wasRowRemoved(row));
   assert(!wasColumnRemoved(col));

   const int old_ixlowupp = int_values[first_int_val];
   const bool is_upper_bound = (int_values[first_int_val + 1] == 1) ? true : false;
   const int index_stored_row = int_values[first_int_val + 2];
   const INDEX stored_row(ROW, row.getNode(), index_stored_row, row.getLinking(), EQUALITY_SYSTEM);

   const double old_bound = float_values[first_float_val];
   const double new_bound = float_values[first_float_val + 1];

   const double curr_x = getSimpleVecFromColStochVec(original_vars.x, col);

   if(is_upper_bound)
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
   const double dual_bound = is_upper_bound ? getSimpleVecFromColStochVec(original_vars.phi, col) + getSimpleVecFromColStochVec(phi, col) :
         getSimpleVecFromColStochVec(original_vars.gamma, col) + getSimpleVecFromColStochVec(gamma, col);

   assert(!row.isLinkingRow());

   /* If the bound was tight all other variables in that row must have been at their respective upper
    * and lower bounds (depending on sings and orientation and upper/lower).
    * This fact will be used to adjust their duals which can be non-zero then since v/w were zero.
    * TODO : assert that
    */

   const double complementarity = slack * dual_bound;

   // TODO : reintroduce once scaling also scales termination criteria
   //   assert( PIPSisLEFeas(0.0, complementarity) );
   // assert( PIPSisLT( complementarity, feas_tol ) );

   /* adjust slack v/w */
   if( old_ixlowupp == 0 )
   {
      slack = 0;
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
      {
         return true;
      }
   }

   /* do dual part of postsolve */
   const double diff_dual_bound = PIPSisZero(slack) ? - dual_bound : complementarity / slack - dual_bound;

   // TODO : reintroduce
//   assert( PIPSisLEFeas(0.0, dual_bound + diff_dual_bound) );
//   assert( PIPSisLEFeas(slack * (dual_bound + diff_dual_bound), complementarity) );

   /* at this point we want to adjust the dual_bound by diff_dual_bound introducing an error in the reduced costs which we need to compensate */
   const double error_reduced_costs = is_upper_bound ? -diff_dual_bound : diff_dual_bound;

   const double coeff = row_storage.getRowCoefficientAtColumn( stored_row, col );

   /* set z/y of corresponding row such that -c_i * delta_z // -a_i * delta_y = error_reduced_costs */
   assert(!PIPSisZero(coeff));
   const double change_dual_row = -error_reduced_costs / coeff;
   if(coeff < 1e-13)
      std::cout << "Potential numerical issues in postsolve of BoundTightening caused by small coefficient" << std::endl;

   /* add z/y * row to gamma/phi */
   StochVectorHandle tmp_pos(dynamic_cast<StochVector*>(original_vars.gamma->clone()));

   /* calculate multiplier_change times row */
   row_storage.axpyAtRow(0.0, *tmp_pos, change_dual_row, stored_row );
   StochVectorHandle tmp_neg(dynamic_cast<StochVector*>(tmp_pos->cloneFull()));

   tmp_pos->selectPositive();
   tmp_neg->selectNegative();

   /* store linking variable changes and allreduce them later except when using a row from D/B0 or D/Bl0 */
//   if( my_rank == 0 || TODOTODOTDOT)
   phi.axpy(1.0, *tmp_neg);
   gamma.axpy(1.0, *tmp_pos);

   // TODO : reintroduce assert( PIPSisLEFeas(slack * dual_bound changes, complementarity) );

   if( row.inInEqSys() )
   {
      double& z = getSimpleVecFromRowStochVec(original_vars.z, row);
      z += change_dual_row;

      /* z = lambda - pi */
      getSimpleVecFromRowStochVec(original_vars.lambda, row) = std::max(0.0, z);
      getSimpleVecFromRowStochVec(original_vars.pi, row) = -std::min(0.0, z);
   }
   else
   {
      getSimpleVecFromRowStochVec(original_vars.y, row) += change_dual_row;
   }

   return true;
}

bool StochPostsolver::postsolveFixedColumn(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == FIXED_COLUMN );

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
   assert(col.isCol());

   const int index_stored_col = int_values.at(first_int_val);
   const INDEX stored_col(COL, col.getNode(), index_stored_col);

   const double value = float_values.at(first_float_val);
   const double obj_coeff = float_values.at(first_float_val + 1);

   assert(wasColumnRemoved(col));

   /* mark entry as set and set x value to fixation */
   getSimpleVecFromColStochVec(*padding_origcol, col) = 1;

   /* set x value */
   getSimpleVecFromColStochVec(original_vars.x, col) = value;

   /* set slacks for x bounds to zero (bounds were tight) */
   getSimpleVecFromColStochVec(original_vars.v, col) = 0.0;
   getSimpleVecFromColStochVec(original_vars.w, col) = 0.0;

   /* set duals for bounds to satisfy reduced costs of reintroduced column times x */
   double col_times_duals = 0.0;
   if( col.isLinkingCol() )
   {
      /* we need to synchronize the column times duals in this case */
      if(my_rank == 0)
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
bool StochPostsolver::postsolveFixedEmptyColumn(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == FIXED_EMPTY_COLUMN );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 1 == next_first_index);
   assert(first_float_val + 4 == next_first_float_val);
   assert(first_int_val + 2 == next_first_int_val);
#endif

   const INDEX& col = indices.at(first_index);
   assert( col.isCol() );
   assert( wasColumnRemoved(col) );
   if( col.isLinkingCol() )
      assert( PIPS_MPIisValueEqual(col.getIndex(), MPI_COMM_WORLD) );

   const double value = float_values.at(first_float_val);
   const double obj_coeff = float_values.at(first_float_val + 1);
   const double lbx = float_values.at(first_float_val + 2);
   const double ubx = float_values.at(first_float_val + 3);
   const int ixlow = int_values.at(first_int_val);
   const int ixupp = int_values.at(first_int_val + 1);

   /* primal */
   /* mark entry as set and set x value to fixation */
   getSimpleVecFromColStochVec(*padding_origcol, col) = 1;
   getSimpleVecFromColStochVec(original_vars.x, col) = value;

   if( ixlow == 1 )
      assert( PIPSisLEFeas(lbx, value) );
   if( ixupp == 1 )
      assert( PIPSisLEFeas(value, ubx) );

   /* dual */
   getSimpleVecFromColStochVec(original_vars.gamma, col) = 0.0;
   getSimpleVecFromColStochVec(original_vars.phi, col) = 0.0;

   if(!PIPSisZero(obj_coeff))
   {
      if( PIPSisLT(obj_coeff, 0.0) )
      {
         assert( ixupp );
         assert( PIPSisEQ(value, ubx) );
         getSimpleVecFromColStochVec(original_vars.phi, col) = obj_coeff;
      }
      else if( PIPSisLT(0.0, obj_coeff) )
      {
         assert( ixlow );
         assert( PIPSisEQ(value, lbx) );
         getSimpleVecFromColStochVec(original_vars.gamma, col) = obj_coeff;
      }
   }

   if( ixlow == 1 )
      getSimpleVecFromColStochVec(original_vars.v, col) = value - lbx;
   else
      getSimpleVecFromColStochVec(original_vars.v, col)= 0.0;

   if( ixupp == 1)
      getSimpleVecFromColStochVec(original_vars.w, col) = ubx - value;
   else
      getSimpleVecFromColStochVec(original_vars.w, col) = 0.0;

   assert( complementarySlackVariablesMet(original_vars, col) );

   return true;
}

bool StochPostsolver::postsolveFixedColumnSingletonFromInequality(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == FIXED_COLUMN_SINGLETON_FROM_INEQUALITY );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 1 == next_first_index);
   assert(first_float_val + 4 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   const INDEX& col = indices.at(first_index);
   assert(col.isCol());
   assert(wasColumnRemoved(col));

   const double value = float_values.at(first_float_val);
//         const double coeff = float_values.at(first_float_val + 1); TODO: not needed - remove
   const double xlow_old = float_values.at(first_float_val + 2);
   const double xupp_old = float_values.at(first_float_val + 3);

   /* set x value - bound is tight - compute slacks of x - set duals to zero */
   getSimpleVecFromColStochVec(*padding_origcol, col) = 1;
   getSimpleVecFromColStochVec(original_vars.x, col) = value;

   /* set slacks for x bounds */
   assert( PIPSisLE(value, xupp_old) || xupp_old == INF_POS_PRES);
   assert( PIPSisLE(xlow_old, value) || xlow_old == INF_NEG_PRES);

   getSimpleVecFromColStochVec(original_vars.v, col) = (xlow_old == INF_NEG_PRES) ? 0 : value - xlow_old;
   getSimpleVecFromColStochVec(original_vars.w, col) = (xupp_old == INF_POS_PRES) ? 0 : xupp_old - value;

   /* TODO: don't think something needs to be done for the reduced costs */
   getSimpleVecFromColStochVec(original_vars.phi, col) = 0.0;
   getSimpleVecFromColStochVec(original_vars.gamma, col) = 0.0;

   /* set duals of bounds of x */
   getSimpleVecFromColStochVec(original_vars.gamma, col) = 0.0;
   getSimpleVecFromColStochVec(original_vars.phi, col) = 0.0;

   return true;
}

bool StochPostsolver::postsolveSingletonEqualityRow(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == SINGLETON_EQUALITY_ROW );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
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
   if( xupp_old == INF_POS_PRES )
   {
      slack_upper = 0.0;
      error_in_reduced_costs -= dual_upper;
      dual_upper = 0.0;
   }
   else
   {
      slack_upper = xupp_old - curr_x;
      assert( PIPSisLE(0.0, slack_upper) );
      assert( std::fabs(slack_upper) != INF_POS_PRES );

      if( !PIPSisZero( slack_upper * dual_upper, postsolve_tol ) )
      {
         error_in_reduced_costs -= dual_upper;
         dual_upper = 0.0;
      }
   }

   /* lower bound */
   if(xlow_old == INF_NEG_PRES)
   {
      slack_lower = 0.0;
      error_in_reduced_costs += dual_lower;
      dual_lower = 0.0;
   }
   else
   {
      slack_lower = curr_x - xlow_old;
      assert(PIPSisLE(0.0, slack_lower));
      assert(std::fabs(slack_lower) != INF_POS_PRES);

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
   const int type = reductions.at(reduction_idx);
   assert( type == SINGLETON_INEQUALITY_ROW );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
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
   const double xupp_new = float_values.at(first_float_val + 3);
   const double coeff = float_values.at(first_float_val + 4);

   assert( !PIPSisZero(coeff) || coeff == NAN );
   assert(xlow_new == INF_NEG_PRES || xupp_new == INF_POS_PRES);
   assert(xlow_new != INF_NEG_PRES || xupp_new != INF_POS_PRES);

   bool lower_bound_changed = (xlow_new != INF_NEG_PRES);

   const double curr_x = getSimpleVecFromColStochVec(original_vars.x, col);
   double& slack = lower_bound_changed ? getSimpleVecFromColStochVec(original_vars.v, col) :
         getSimpleVecFromColStochVec(original_vars.w, col);
   double& dual_bound = lower_bound_changed ? getSimpleVecFromColStochVec(original_vars.gamma, col) :
         getSimpleVecFromColStochVec(original_vars.phi, col);

   const double old_bound = lower_bound_changed ? xlow_old : xupp_old;
   const double new_bound = lower_bound_changed ? xlow_new : xupp_new;
   assert( std::fabs(new_bound) < std::fabs(old_bound) );

   double error_in_reduced_costs = 0.0;
   /* adjust bounds slacks and duals so that complementarity condition stays valid*/
   if(std::fabs(old_bound) == INF_POS_PRES)
   {
      slack = 0.0;

      error_in_reduced_costs += lower_bound_changed ? dual_bound : -dual_bound;
      dual_bound = 0.0;
   }
   else
   {
      slack = std::fabs(old_bound - curr_x);
      assert(PIPSisLT(0.0, slack));
      assert(std::fabs(slack) != INF_POS_PRES);

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

bool StochPostsolver::postsolveFreeColumnSingletonEquality(sVars& original_vars, int reduction_idx) const
{
   /* row can be an equality row but then it must have clow == cupp */
   const int type = reductions.at(reduction_idx);
   assert( type == FREE_COLUMN_SINGLETON_EQUALITY );

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

   assert(std::abs(value_row) != INF_POS_PRES);

   /* reintroduce the removed column on process owning column */
   if( col.isCol() )
   {
      double& x_val = getSimpleVecFromColStochVec(original_vars.x, col);
      assert(PIPSisZero(x_val));

      /* mark column as set */
      getSimpleVecFromColStochVec(*padding_origcol, col) = 1;

      /* recover primal value */
      x_val = (rhs - value_row) / col_coeff;
      assert( PIPSisZero(x_val * col_coeff + value_row - rhs) );

      /* compute slacks and set duals for bounds to zero */
      double& slack_lower = getSimpleVecFromColStochVec(original_vars.v, col);
      double& slack_upper = getSimpleVecFromColStochVec(original_vars.w, col);

      if( xlow == INF_NEG_PRES )
         slack_lower = 0.0;
      else
      {
         assert( PIPSisLE(xlow, x_val) );
         slack_lower = xlow - x_val;
      }

      if( xupp == INF_POS_PRES )
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

bool StochPostsolver::postsolveNearlyParallelRowSubstitution(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == NEARLY_PARALLEL_ROW_SUBSTITUTION );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 4 == next_first_index);
   assert(first_float_val + 8 == next_first_float_val);
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
   assert( wasColumnRemoved(col2) );

   assert( !wasRowRemoved(row1) );
   assert( !wasRowRemoved(row2) );

   const double scalar = float_values.at(first_float_val);
   const double translation = float_values.at(first_float_val + 1);
   const double obj_col1 = 0.0;
   const double obj_col2 = float_values.at(first_float_val + 2);

   const double xlow_col2 = float_values.at(first_float_val + 3);
   const double xupp_col2 = float_values.at(first_float_val + 4);

   const double coeff_col1 = float_values.at(first_float_val + 5);
   const double coeff_col2 = float_values.at(first_float_val + 6);

   /* row1 = parallel_factor * row2 */
   const double parallel_factor = float_values.at(first_float_val + 7);

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
   getSimpleVecFromColStochVec(*padding_origcol, col2) = 1;
   getSimpleVecFromColStochVec(original_vars.x, col2) = val_col2;
   assert( PIPSisLEFeas( xlow_col2, val_col2 ) );
   assert( PIPSisLEFeas( val_col2, xupp_col2 ) );

   /* slacks for bounds */
   getSimpleVecFromColStochVec(original_vars.v, col2) = (xlow_col2 == INF_NEG_PRES) ? 0 : val_col2 - xlow_col2;
   getSimpleVecFromColStochVec(original_vars.w, col2) = (xupp_col2 == INF_POS_PRES) ? 0 : xupp_col2 - val_col2;

   /* bound duals to zero */
   getSimpleVecFromColStochVec(original_vars.phi, col2) = 0.0;
   getSimpleVecFromColStochVec(original_vars.gamma, col2) = 0.0;

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
      assert( !PIPSisZero(parallel_factor) );
      assert( !PIPSisZero(coeff_col2) );

      const double change_objective_row1 = obj_col2 * coeff_col1 / parallel_factor / coeff_col2;

      assert( !PIPSisZero(obj_col1 + change_objective_row1) );
      const double change_factor_obj_row1 = change_objective_row1 / (obj_col1 + change_objective_row1);

      /* shift duals from row1 to row2 */
      double& z_row1 = getSimpleVecFromRowStochVec(original_vars.z, row1);
      double& lambda_row1 = getSimpleVecFromRowStochVec(original_vars.lambda, row1);
      double& pi_row1 = getSimpleVecFromRowStochVec(original_vars.pi, row1);

      double& z_row2 = getSimpleVecFromRowStochVec(original_vars.z, row2);
      double& lambda_row2 = getSimpleVecFromRowStochVec(original_vars.lambda, row2);
      double& pi_row2 = getSimpleVecFromRowStochVec(original_vars.pi, row2);

      /* z = lambda - pi */
      z_row2 = (1 - change_factor_obj_row1) * parallel_factor * z_row1;
      pi_row2 = (1 - change_factor_obj_row1) * parallel_factor * pi_row1;
      lambda_row2 = (1 - change_factor_obj_row1) * parallel_factor * lambda_row1;

      z_row1 *= change_factor_obj_row1;
      pi_row1 *= change_factor_obj_row1;
      lambda_row1 *= change_factor_obj_row1;

      /* shift bound duals */
      double& gamma_row1 = getSimpleVecFromRowStochVec(original_vars.gamma, row1);
      double& phi_row1 = getSimpleVecFromRowStochVec(original_vars.phi, row1);

      double& gamma_row2 = getSimpleVecFromRowStochVec(original_vars.gamma, row2);
      double& phi_row2 = getSimpleVecFromRowStochVec(original_vars.phi, row2);

      gamma_row2 = (1 - change_factor_obj_row1) * gamma_row1;
      phi_row2 = (1 - change_factor_obj_row1) * phi_row1;

      gamma_row1 *= change_factor_obj_row1;
      phi_row1 *= change_factor_obj_row1;


#ifndef NDEBUG
      const double& v_row1 = getSimpleVecFromColStochVec(original_vars.v, col1);
      const double& v_row2 = getSimpleVecFromColStochVec(original_vars.v, col2);
      const double& w_row1 = getSimpleVecFromColStochVec(original_vars.w, col1);
      const double& w_row2 = getSimpleVecFromColStochVec(original_vars.w, col2);

      const double& t_row1 = getSimpleVecFromColStochVec(original_vars.t, col1);
      const double& t_row2 = getSimpleVecFromColStochVec(original_vars.t, col2);
      const double& u_row1 = getSimpleVecFromColStochVec(original_vars.u, col1);
      const double& u_row2 = getSimpleVecFromColStochVec(original_vars.u, col2);
      /* assert complementary slackness condition */
      assert( PIPSisZeroFeas(v_row1 * gamma_row1) );
      assert( PIPSisZeroFeas(v_row2 * gamma_row2) );
      assert( PIPSisZeroFeas(w_row1 * phi_row1) );
      assert( PIPSisZeroFeas(w_row2 * phi_row2) );
      assert( PIPSisZeroFeas(t_row1 * lambda_row1) );
      assert( PIPSisZeroFeas(t_row2 * lambda_row2) );
      assert( PIPSisZeroFeas(u_row1 * pi_row1) );
      assert( PIPSisZeroFeas(u_row2 * pi_row2) );
#endif
   }

   return true;
}

bool StochPostsolver::postsolveNearlyParallelRowBoundsTightened(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == NEARLY_PARALLEL_ROW_BOUNDS_TIGHTENED );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
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

   const double val_col1 = getSimpleVecFromColStochVec(original_vars.x, col1);

   /* if the variable bound of col1 was actually implied via col2 we have to shift it's dual multipliers over via also adjusting the dual of row1 */
   // TODO : reactivate these once scaled termination criteria exist - found bug elsewhere - might be correct and maybe can be reactivated
//         assert( PIPSisRelLEFeas( xlow_col1, val_col1 ) );
//         assert( PIPSisRelLEFeas( val_col1, xupp_col1 ) );

   double xlow_implied = INF_NEG_PRES;
   double xupp_implied = INF_POS_PRES;

   /* two nearly parallel equality rows */
   if( col2.isCol() )
   {
      assert( xlow_col2 != INF_NEG_PRES || xupp_col2 != INF_POS_PRES );

      if( PIPSisLT(0.0, scalar) )
      {
         if( xlow_col2 != INF_NEG_PRES )
            xlow_implied = std::max(xlow_col1, (xlow_col2 - translation) / scalar);
         else
            xlow_implied = xlow_col1;

         if( xupp_col2 != INF_POS_PRES )
            xupp_implied = std::min(xupp_col1, (xupp_col2 - translation) / scalar);
         else
            xupp_implied = xupp_col1;
      }
      else
      {
         if( xlow_col2 != INF_NEG_PRES )
            xupp_implied = std::min(xupp_col1, (xlow_col2 - translation) / scalar);
         else
            xupp_implied = xupp_col1;

         if( xupp_col2 != INF_POS_PRES )
            xlow_implied = std::max(xlow_col1, (xupp_col2 - translation) / scalar);
         else
            xlow_implied = xlow_col1;
      }
   }
   else
   {
      assert( clow != INF_NEG_PRES || cupp != INF_POS_PRES );
      assert( !PIPSisZero(coeff_col1) );
      assert( !PIPSisZero(parallel_factor * coeff_col1) );

      const double faq = parallel_factor * coeff_col1;

      assert( !PIPSisZero(faq) );
      if( PIPSisLT( 0.0, faq ) )
      {
         if( cupp != INF_POS_PRES )
            xlow_implied = std::max( xlow_col1, (rhs - parallel_factor * cupp ) / coeff_col1 );
         else
            xlow_implied = xlow_col1;

         if( clow != INF_NEG_PRES )
            xupp_implied = std::min( xupp_col1, (rhs - parallel_factor * clow) / coeff_col1 );
         else
            xupp_implied = xupp_col1;
      }

      if( PIPSisLT( faq, 0.0 ) )
      {
         if( cupp != INF_POS_PRES )
            xupp_implied = std::min( xupp_col1, (rhs - parallel_factor * cupp) / coeff_col1 );
         else
            xupp_implied = xupp_col1;

         if( clow != INF_NEG_PRES )
            xlow_implied = std::max( xlow_col1, (rhs - parallel_factor * clow) / coeff_col1 );
         else
            xlow_implied = xlow_col1;
      }
   }

   assert( xlow_implied != INF_NEG_PRES );
   assert( xupp_implied != INF_POS_PRES );
   assert( PIPSisRelLEFeas(xlow_col1, xlow_implied) );
   assert( PIPSisRelLEFeas(xupp_implied, xupp_col1) );

   /* lower bound was implied by substituted column + it's row */
   if( PIPSisEQ( val_col1, xlow_implied ) && !PIPSisEQ( xlow_implied, xlow_col1 ) )
   {
      assert( !PIPSisEQ(val_col1, xlow_col1) );
      double& dual_col1 = getSimpleVecFromColStochVec(original_vars.gamma, col1);

      if( !PIPSisZero(dual_col1) )
      {
         /* set the dual of col1 to zero, compensate for the error made with the dual of row1, compensate that error with the dual of col2
          * and finally compensate the error in col2's row of the reduced costs with col2's dual
          */

         /* only one bound dual variable is non-zero */
         assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.phi, col1)) );

         /* dual to zero and compensate error */
         const double dual_shift_row1 = dual_col1 / coeff_col1;
         dual_col1 = 0;
         assert( PIPSisEQ( -coeff_col1 * dual_shift_row1, -dual_col1 ) );

         getSimpleVecFromRowStochVec(original_vars.y, row1) += dual_shift_row1;

         /* compensate again and adjust dual of col2 */
         const double dual_shift_row2 = -dual_shift_row1 * parallel_factor;

         if( row2.inEqSys() )
            getSimpleVecFromRowStochVec(original_vars.y, row2) += dual_shift_row2;
         else
         {
            getSimpleVecFromRowStochVec(original_vars.z, row2) += dual_shift_row2;
            getSimpleVecFromRowStochVec(original_vars.lambda, row2) = std::max(0.0, getSimpleVecFromRowStochVec(original_vars.z, row2));
            getSimpleVecFromRowStochVec(original_vars.pi, row2) += std::min(0.0, getSimpleVecFromRowStochVec(original_vars.z, row2));
         }

         /* if we had two equality rows there was a singleton col in row2 an it's dual needs to be adjusted */
         if( col2.isCol() )
         {
            /* assert bounds of substituted variable are tight as well */
            assert( PIPSisEQ(getSimpleVecFromColStochVec(original_vars.x, col2) , PIPSisLT(0.0, translation) ? xupp_col2 : xlow_col2) );

            const double dual_shift_col2 = coeff_col2 * dual_shift_row2;
            PIPSisLT(0.0, translation) ? assert( PIPSisLE(0.0, dual_shift_col2) ) : assert( PIPSisLE( dual_shift_col2, 0.0) );

            double& dual_col2 = PIPSisLT(0.0, translation) ? getSimpleVecFromColStochVec(original_vars.phi, col2) :
                  getSimpleVecFromColStochVec(original_vars.gamma, col2);
            assert( PIPSisZero(dual_col2) );
            dual_col2 += std::fabs(dual_shift_col2);
         }
      }
   }

   /* upper bound was implied by substituted column + it's row */
   if( PIPSisEQ( val_col1, xupp_implied ) && !PIPSisEQ( xupp_implied, xupp_col1 ) )
   {
      assert( !PIPSisEQ(val_col1, xupp_col1) );
      double& dual_col1 = getSimpleVecFromColStochVec(original_vars.phi, col1);

      if( !PIPSisZero(dual_col1) )
      {
         /* set the dual of col1 to zero, compensate for the error made with the dual of row1, compensate that error with the dual of col2
          * and finally compensate the error in col2's row of the reduced costs with col2's dual
          */

         /* only one bound dual variable is non-zero */
         assert( PIPSisZero(getSimpleVecFromColStochVec(original_vars.gamma, col1)) );

         /* dual to zero and compensate error */
         const double dual_shift_row1 = dual_col1 / coeff_col1;
         dual_col1 = 0;
         assert( PIPSisEQ( -coeff_col1 * dual_shift_row1, -dual_col1 ) );

         getSimpleVecFromRowStochVec(original_vars.y, row1) += dual_shift_row1;

         /* compensate again and adjust dual of col2 */
         const double dual_shift_row2 = -dual_shift_row1 * parallel_factor;

         if( row2.inEqSys() )
            getSimpleVecFromRowStochVec(original_vars.y, row2) += dual_shift_row2;
         else
         {
            getSimpleVecFromRowStochVec(original_vars.z, row2) += dual_shift_row2;
            getSimpleVecFromRowStochVec(original_vars.lambda, row2) = std::max(0.0, getSimpleVecFromRowStochVec(original_vars.z, row2));
            getSimpleVecFromRowStochVec(original_vars.pi, row2) += std::min(0.0, getSimpleVecFromRowStochVec(original_vars.z, row2));
         }

         if( col2.isCol() )
         {
            /* assert bounds of substituted variable are tight as well */
            assert( PIPSisEQ(getSimpleVecFromColStochVec(original_vars.x, col2), PIPSisLT(0.0, translation) ? xlow_col2 : xupp_col2) );

            double& dual_col2 = PIPSisLT(0.0, translation) ? getSimpleVecFromColStochVec(original_vars.gamma, col2) :
                  getSimpleVecFromColStochVec(original_vars.phi, col2);
            assert( PIPSisZero(dual_col2) );

            const double dual_shift_col2 = coeff_col2 * dual_shift_row2;

            PIPSisLT(0.0, translation) ? assert( PIPSisLE(0.0, dual_shift_col2) ) : assert( PIPSisLE( dual_shift_col2, 0.0) );
            dual_col2 += std::fabs(dual_shift_col2);
         }
      }
   }
   return true;
}

bool StochPostsolver::postsolveFreeColumnSingletonInequalityRow( sVars& original_vars, const sData& original_problem, int reduction_idx ) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == FREE_COLUMN_SINGLETON_INEQUALITY_ROW );

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

   assert( !wasColumnRemoved(col) );
   assert( wasRowRemoved(row) );

   getSimpleVecFromColStochVec(*padding_origcol, col) = 1;

   getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) = 1;

   /* set x value such that row is satisfied */
   double& x_val = getSimpleVecFromColStochVec(original_vars.x, col);
   assert( PIPSisZero(x_val) );

   const double row_value = row_storage.multRowTimesVec(row_stored, dynamic_cast<const StochVector&>(*original_vars.x));
   x_val = (rhs - row_value) / coeff;

   /* duals of row and bounds are zero */
   getSimpleVecFromRowStochVec(original_vars.z, row) = 0;
   getSimpleVecFromRowStochVec(original_vars.lambda, row) = 0;
   getSimpleVecFromRowStochVec(original_vars.pi, row) = 0;

   getSimpleVecFromColStochVec(original_vars.gamma, col) = 0;
   getSimpleVecFromColStochVec(original_vars.phi, col) = 0;

   /* compute slacks for bounds and row */
   /* slack for row is zero by construction */
   getSimpleVecFromRowStochVec(original_vars.t, row) = 0;
   getSimpleVecFromRowStochVec(original_vars.u, row) = 0;
   getSimpleVecFromRowStochVec(original_vars.s, row) = row_value + coeff * x_val;

   if( xlow == INF_NEG_PRES )
      getSimpleVecFromColStochVec(original_vars.v, col) = 0;
   else
   {
      getSimpleVecFromColStochVec(original_vars.v, col) = x_val - xlow;
      assert( PIPSisLE(xlow, x_val) );
   }

   if( xupp == INF_POS_PRES)
      getSimpleVecFromColStochVec(original_vars.w, col) = 0;
   else
   {
      getSimpleVecFromColStochVec(original_vars.w, col) = xupp - x_val;
      assert( PIPSisLE(x_val, xupp) );
   }

   getSimpleVecFromColStochVec(original_vars.gamma, col) = 0.0;
   getSimpleVecFromColStochVec(original_vars.phi, col) = 0.0;

   return true;
}

bool StochPostsolver::postsolveParallelRowsBoundsTightened(sVars& original_vars, int reduction_idx) const
{
   const int type = reductions.at(reduction_idx);
   assert( type == PARALLEL_ROWS_BOUNDS_TIGHTENED );

   const unsigned int first_float_val = start_idx_float_values.at(reduction_idx);
   const unsigned int first_int_val = start_idx_int_values.at(reduction_idx);
   const unsigned int first_index = start_idx_indices.at(reduction_idx);

#ifndef NDEBUG
   const unsigned int next_first_float_val = start_idx_float_values.at(reduction_idx + 1);
   const unsigned int next_first_int_val = start_idx_int_values.at(reduction_idx + 1);
   const unsigned int next_first_index = start_idx_indices.at(reduction_idx + 1);

   assert(first_index + 2 == next_first_index);
   assert(first_float_val + 5 == next_first_float_val);
   assert(first_int_val == next_first_int_val);
#endif

   const INDEX& row1 = indices.at(first_index);
   const INDEX& row2 = indices.at(first_index + 1);

   assert(row1.isRow());
   assert(row2.isRow());
   assert(row1.inInEqSys());
   assert(row2.inInEqSys());

   const double clow_old = float_values.at(first_float_val);
   const double cupp_old = float_values.at(first_float_val + 1);
   const double clow_new = float_values.at(first_float_val + 2);
   const double cupp_new = float_values.at(first_float_val + 3);
   const double factor = float_values.at(first_float_val + 4);

   assert (!PIPSisZero(factor) );
   assert( getSimpleVecFromRowStochVec(*padding_origrow_inequality, row1) == 1 );
   assert( getSimpleVecFromRowStochVec(*padding_origrow_inequality, row2) == 1 );

   /* recompute duals and slack of both rows - if one bound was tight and thus the dual non-zero we shift it to the row originally implying the bound */
   double& z_row1 = getSimpleVecFromRowStochVec(original_vars.z, row1);
   double& lambda_row1 = getSimpleVecFromRowStochVec(original_vars.lambda, row1);
   double& pi_row1 = getSimpleVecFromRowStochVec(original_vars.pi, row1);

   double& t_row1 = getSimpleVecFromRowStochVec(original_vars.t, row1);
   double& u_row1 = getSimpleVecFromRowStochVec(original_vars.u, row1);

   double& z_row2 = getSimpleVecFromRowStochVec(original_vars.z, row2);
   double& lambda_row2 = getSimpleVecFromRowStochVec(original_vars.lambda, row2);
   double& pi_row2 = getSimpleVecFromRowStochVec(original_vars.pi, row2);

#ifndef NDEBUG
   const double& t_row2 = PIPSisLT(factor, 0.0) ? getSimpleVecFromRowStochVec(original_vars.u, row2) :
         getSimpleVecFromRowStochVec(original_vars.t, row2);
   const double& u_row2 = PIPSisLT(factor, 0.0) ? getSimpleVecFromRowStochVec(original_vars.t, row2) :
         getSimpleVecFromRowStochVec(original_vars.u, row2);
#endif
   assert( PIPSisZero(z_row2) );
   assert( PIPSisZero(pi_row2) );
   assert( PIPSisZero(lambda_row2) );

   const bool clow_impied_by_row2 = (clow_old != clow_new);
   const bool cupp_impied_by_row2 = (cupp_old != cupp_new);

   if( PIPSisLT(factor, 0.0) )
      std::swap(lambda_row2, pi_row2);

   /* if parallel row is tight on a bounds implied by another row we have to move the multipliers to the implying row */
   if( clow_impied_by_row2 && PIPSisZero(t_row1) )
   {
      assert( PIPSisLT(clow_old, clow_new) );
      assert( PIPSisZero(t_row2) );

      /* adjust slacks */
      if( clow_old == INF_NEG_PRES )
         t_row1 = 0.0;
      else
         t_row1 += clow_new - clow_old;

      /* clow is no longer tight on row1 but now on row2 - shift row1's multiplier to row2 */
      z_row1 -= lambda_row1;
      z_row2 += lambda_row1 * factor;

      lambda_row2 = lambda_row1 * factor;
      lambda_row1 = 0;

      assert( PIPSisZeroFeas(t_row2 * lambda_row2) );
   }

   if( cupp_impied_by_row2 && PIPSisZero(u_row1) )
   {
      assert( PIPSisLT(cupp_new, cupp_old) );
      assert( PIPSisZero(u_row2) );

      /* adjust slacks */
      if( cupp_old == INF_POS_PRES )
         u_row1 = 0.0;
      else
         u_row1 += cupp_old - cupp_new;

      /* cupp is no longer tight on row1 but now on row2 - shift row1's multiplier to row2 */
      z_row1 += pi_row1;
      z_row2 -= pi_row1 * factor;

      pi_row2 = pi_row1 * factor;
      pi_row1 = 0.0;

      assert( PIPSisZeroFeas(u_row2 * pi_row2) );
   }

   return true;
}

/* sync linking variables - either the variables are not set everywhere or on some procs they are set to something non-zero while on others they are zero */
bool StochPostsolver::syncNewlySetLinkingVars(sVars& original_vars) const
{
   SimpleVector& xl_svec = dynamic_cast<SimpleVector&>( *dynamic_cast<StochVector&>(*original_vars.x).vec );

   std::vector<bool> modified(xl_svec.length(), false);

   for( int i = 0; i < xl_svec.length(); ++i )
   {
      if( xl_svec[i] == 0.0 )
      {
         xl_svec[i] = INF_NEG_PRES;
         modified[i] = true;
      }
   }

   PIPS_MPImaxArrayInPlace(xl_svec.elements(), xl_svec.length(), MPI_COMM_WORLD);

   /* duals and slacks are all component-wise positive thus allreduce max will give the desired result (so either neg inf, 0.0 or some actually set value */
   SimpleVector& v = dynamic_cast<SimpleVector&>( *dynamic_cast<StochVector&>(*original_vars.v).vec );
   SimpleVector& w = dynamic_cast<SimpleVector&>( *dynamic_cast<StochVector&>(*original_vars.w).vec );
   PIPS_MPImaxArrayInPlace(v.elements(), v.length(), MPI_COMM_WORLD);
   PIPS_MPImaxArrayInPlace(w.elements(), w.length(), MPI_COMM_WORLD);

   SimpleVector& gamma = dynamic_cast<SimpleVector&>( *dynamic_cast<StochVector&>(*original_vars.gamma).vec );
   SimpleVector& phi = dynamic_cast<SimpleVector&>( *dynamic_cast<StochVector&>(*original_vars.phi).vec );
   PIPS_MPImaxArrayInPlace(gamma.elements(), gamma.length(), MPI_COMM_WORLD);
   PIPS_MPImaxArrayInPlace(phi.elements(), phi.length(), MPI_COMM_WORLD);

   for( int i = 0; i < xl_svec.length(); ++i)
   {
      if( modified.at(i) && xl_svec[i] == INF_NEG_PRES )
         xl_svec[i] = 0.0;
   }

   return true;
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
   const double initial_const = -std::numeric_limits<double>::infinity();
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
   const double initial_const = -std::numeric_limits<double>::infinity();

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

bool StochPostsolver::complementarySlackVariablesMet(const sVars& vars, const INDEX& col) const
{
   assert( col.isCol() );
   assert( !wasColumnRemoved(col) );
   const double v = getSimpleVecFromColStochVec(vars.v, col);
   const double w = getSimpleVecFromColStochVec(vars.w, col);
   const double gamma = getSimpleVecFromColStochVec(vars.gamma, col);
   const double phi = getSimpleVecFromColStochVec(vars.phi, col);

   return std::fabs(v * gamma) < postsolve_tol && std::fabs(w * phi) < postsolve_tol ;
}

bool StochPostsolver::complementarySlackRowMet(const sVars& vars, const INDEX& row) const
{
   assert( row.isRow() );
   assert( !wasRowRemoved(row) );

   const double t = getSimpleVecFromRowStochVec(vars.t, row);
   const double u = getSimpleVecFromRowStochVec(vars.u, row);
   const double lambda = getSimpleVecFromRowStochVec(vars.lambda, row);
   const double pi = getSimpleVecFromRowStochVec(vars.pi, row);

   return std::fabs(t * lambda) < postsolve_tol && std::fabs(u * pi) < postsolve_tol;
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
