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
   col_stored_last_at( dynamic_cast<StochVectorBase<int>*>(padding_origcol->clone()) )
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

   // count !? rows cols...
   start_idx_float_values.push_back(0);
   start_idx_int_values.push_back(0);
   start_idx_indices.push_back(0);
}

StochPostsolver::~StochPostsolver()
{
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

void StochPostsolver::notifyFreeColumnSingletonInequalityRow( const INDEX& row, const INDEX& col, double lhsrhs, double coeff, const StochGenMatrix& matrix_row )
{
   assert(col.isCol());
   assert(row.isRow());

   markColumnRemoved(col);
   markRowRemoved(row);

   reductions.push_back(FREE_COLUMN_SINGLETON_INEQUALITY_ROW);

   indices.push_back(row);
   indices.push_back(col);

   const int index_stored_row = row_storage.storeRow(row, matrix_row);

   int_values.push_back(index_stored_row);

   float_values.push_back(lhsrhs);
   float_values.push_back(coeff);

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
   assert(row.isRow());
   assert(col.isCol());
   assert(!row.getLinking());
   assert(PIPSisLE(xlow_new, xupp_new));
   assert(xupp_new != INF_POS_PRES || xlow_new != INF_NEG_PRES);
   assert(!PIPSisZero(coeff));

   assert(!wasRowRemoved(row));
   assert(!wasColumnRemoved(col));

   ReductionType red = row.inEqSys() ? SINGLETON_EQUALITY_ROW : SINGLETON_INEQUALITY_ROW;

   if(red == SINGLETON_EQUALITY_ROW)
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
   assert(!wasColumnRemoved(col));
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

   assert(!wasRowRemoved(row));
   markRowRemoved(row);

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

// todo : notify for linking constraints will be a sync event
void StochPostsolver::notifyRowPropagatedBound( const INDEX& row, const INDEX& col, int old_ixlowupp, double old_bound, double new_bound, bool is_upper_bound, const StochGenMatrix& matrix_row)
{
   assert(!PIPSisEQ(old_bound, new_bound));
   assert(row.isRow());
   assert(col.isCol());

   if(is_upper_bound)
      assert( PIPSisLT(new_bound, old_bound) );
   else
      assert( PIPSisLT(old_bound, new_bound) );

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

void StochPostsolver::putLinkingVarsSyncEvent()
{
   reductions.push_back( LINKING_VARS_SYNC_EVENT );
   /// dummy : todo change structure - actually we only need a reduction, nothing else..
   indices.push_back( INDEX(COL, -2, -2 ) );
   finishNotify();
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

   stoch_original_sol.x->setToConstant(std::numeric_limits<double>::infinity());

   /* original variables are now reduced vars padded with zeros */
   setOriginalVarsFromReduced(stoch_reduced_sol, stoch_original_sol);

   /* primal variables */
   StochVector& x_vec = dynamic_cast<StochVector&>(*stoch_original_sol.x);

   /* dual variables */

   StochVector& y_vec = dynamic_cast<StochVector&>(*stoch_original_sol.y);

   StochVector& z_vec = dynamic_cast<StochVector&>(*stoch_original_sol.z);
   StochVector& lambda_vec = dynamic_cast<StochVector&>(*stoch_original_sol.lambda);
   StochVector& pi_vec = dynamic_cast<StochVector&>(*stoch_original_sol.pi);

   StochVector& s_vec = dynamic_cast<StochVector&>(*stoch_original_sol.s);
   StochVector& t_vec = dynamic_cast<StochVector&>(*stoch_original_sol.t);
   StochVector& u_vec = dynamic_cast<StochVector&>(*stoch_original_sol.u);

   StochVector& gamma_vec = dynamic_cast<StochVector&>(*stoch_original_sol.gamma);
   StochVector& phi_vec = dynamic_cast<StochVector&>(*stoch_original_sol.phi);
   StochVector& v_vec = dynamic_cast<StochVector&>(*stoch_original_sol.v);
   StochVector& w_vec = dynamic_cast<StochVector&>(*stoch_original_sol.w);

   const sData& orig_problem_s = dynamic_cast<const sData&>(original_problem);

   /* post-solve the reductions in reverse order */
   for( int i = reductions.size() - 1; i >= 0; --i )
   {
      const int type = reductions.at(i);
      const unsigned int first_float_val = start_idx_float_values.at(i);
      const unsigned int first_int_val = start_idx_int_values.at(i);
      const unsigned int first_index = start_idx_indices.at(i);

#ifndef NDEBUG
      const unsigned int next_first_float_val = start_idx_float_values.at(i + 1);
      const unsigned int next_first_int_val = start_idx_int_values.at(i + 1);
      const unsigned int next_first_index = start_idx_indices.at(i + 1);
#endif

      switch( type )
      {
      case REDUNDANT_ROW:
      {
         /**
          * postsolve for this is simply to set all dual variables to zero - the row itself has no primal impact
          * for linking rows:
          *    all processes should have removed them in the same order
          *       -> assert in place to check this
          *    the current activity gets synchronized for slack computation
          */

         /* only dual postsolve */
         assert(first_index + 1 == next_first_index);
         assert(first_float_val + 2 == next_first_float_val);
         assert(first_int_val + 3 == next_first_int_val);

         const INDEX& row_idx = indices.at(first_index);
         assert(row_idx.isRow());

         const int row = row_idx.getIndex();
         const int node = row_idx.getNode();
         const bool linking_row = row_idx.getLinking();
         const SystemType system_type = row_idx.getSystemType();

         const double lhs = float_values.at(first_float_val);
         const double rhs = float_values.at(first_float_val + 1);

         const int index_stored_row = int_values.at(first_int_val);
         const int iclow = int_values.at(first_int_val + 1);
         const int icupp = int_values.at(first_int_val + 2);
         assert(iclow + icupp >= 1);

         double value_row = 0.0;

         // TODO: this could be optimized - all linking rows will be after each other - we could wait with allreduce
         // until the current boost of linking rows is processed and then allreduce all activities at once and unpate the slacks then only
         /* get current row activity - redundant linking rows have to lie on the stack in the same order */
         if(linking_row)
         {
            assert(node == -1);
            assert(PIPS_MPIisValueEqual(row, MPI_COMM_WORLD));

            if(my_rank == 0)
               value_row = row_storage.multRowTimesVec( INDEX(ROW, node, index_stored_row, linking_row, EQUALITY_SYSTEM), x_vec );
            else
               value_row = row_storage.multLinkingRowTimesVecWithoutBl0(index_stored_row, x_vec);
            /* this might get very expensive if there is many redundant linking rows */
            PIPS_MPIgetSumInPlace(value_row, MPI_COMM_WORLD);
         }
         else
            value_row = row_storage.multRowTimesVec( INDEX(ROW, node, index_stored_row, linking_row, EQUALITY_SYSTEM), x_vec );

         if( system_type == EQUALITY_SYSTEM )
         {
            assert(-1 <= node && node < static_cast<int>(y_vec.children.size()));
            assert(wasRowRemoved(row_idx));
            assert(PIPSisEQ(lhs, rhs));
            assert(PIPSisLEFeas(value_row, rhs));

            /* set dual multiplier to zero and mark row as added */
            getSimpleVecFromRowStochVec(*padding_origrow_equality, node, linking_row)[row] = 1;
            getSimpleVecFromRowStochVec(y_vec, node, linking_row)[row] = 0;
         }
         else
         {
            assert(-1 <= node && node < static_cast<int>(z_vec.children.size()));
            assert(wasRowRemoved(row_idx));

            /* set dual multipliers to zero and mark row as added */
            getSimpleVecFromRowStochVec(*padding_origrow_inequality, node, linking_row)[row] = 1;

            /* dual of row is zero */
            getSimpleVecFromRowStochVec(z_vec, node, linking_row)[row] = 0;
            getSimpleVecFromRowStochVec(lambda_vec, node, linking_row)[row] = 0;
            getSimpleVecFromRowStochVec(pi_vec, node, linking_row)[row] = 0;

            getSimpleVecFromRowStochVec(s_vec, node, linking_row)[row] = value_row;

            if( iclow == 1 )
               assert(PIPSisLEFeas(lhs, value_row));
            if( icupp == 1 )
               assert(PIPSisLEFeas(value_row, rhs));

            /* set correct slacks */
            if( iclow == 1)
               getSimpleVecFromRowStochVec(t_vec, node, linking_row)[row] = value_row - lhs;
            else
               getSimpleVecFromRowStochVec(t_vec, node, linking_row)[row] = 0;

            if( icupp == 1)
               getSimpleVecFromRowStochVec(u_vec, node, linking_row)[row] = rhs - value_row;
            else
               getSimpleVecFromRowStochVec(u_vec, node, linking_row)[row] = 0;
         }
         break;
      }
      case BOUNDS_TIGHTENED:
      {
//         indices.push_back( row );
//         indices.push_back( col );
//
//         int index_stored_row = storeRow( row, matrix_row );
//
//         int_values.push_back(old_ixlowupp);
//         int_values.push_back(is_upper_bound);
//         int_values.push_back(index_stored_row);
//         float_values.push_back(old_bound);
//         float_values.push_back(new_bound);
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 2 == next_first_float_val);
         assert(first_int_val + 3 == next_first_int_val);

         const INDEX& row = indices.at(first_index);
         const INDEX& col = indices.at(first_index + 1);
         assert(row.isRow());
         assert(col.isCol());

//         const int node_row = row_idx.getNode();
//         const int row = row_idx.getIndex();
//         const bool linking_row = row_idx.getLinking();
//         const SystemType system_type = row_idx.getSystemType();

         const int old_ixlowupp = int_values[first_int_val];
         const bool is_upper_bound = (int_values[first_int_val + 1] == 1) ? true : false;
         const int index_stored_row = int_values[first_int_val + 2];

         const double old_bound = float_values[first_float_val];
         const double new_bound = float_values[first_float_val + 1];

         const double curr_x = getSimpleVecFromColStochVec(x_vec, col);

         if(is_upper_bound)
            assert( PIPSisLEFeas(curr_x, new_bound) );
         else
            assert( PIPSisLEFeas(new_bound, curr_x) );

         double& slack = is_upper_bound ? getSimpleVecFromColStochVec(w_vec, col) :
            getSimpleVecFromColStochVec(v_vec, col);

         // TODO : merge slack computations of tight and non-tight == refactoring
         /* if bound is not tight only adjust v/w */
         if( !PIPSisEQ(curr_x, new_bound ) )
         {
            if(old_ixlowupp == 0)
            {
               slack = 0;
            }
            else if(is_upper_bound)
            {
               assert(PIPSisLT(new_bound, old_bound));
               slack += old_bound - new_bound;
            }
            else
            {
               assert(PIPSisLT(old_bound, new_bound));
               slack += new_bound - old_bound;
            }
         }
         /* we always add an epsilon to bounds found so it cannot have been tight - bound was tight */
         else
         {
            assert(false);
            assert(!row.isLinkingRow());
            /* If the bound was tight all other variables in that row must have been at their respective upper
             * and lower bounds (depending on sings and orientation and upper/lower).
             * This fact will be used to adjust their duals which can be non-zero then since v/w were zero.
             * TODO : assert that
             */
            assert(!PIPSisEQ(old_bound, curr_x));

            /* adjust slack v/w */
            if(old_ixlowupp == 0)
            {
               slack = 0;
            }
            else if(is_upper_bound)
            {
               assert(PIPSisLE(0, new_bound - curr_x));
               assert(PIPSisZero(getSimpleVecFromColStochVec(w_vec, col)));
               slack = new_bound - curr_x;
            }
            else
            {
               assert(PIPSisLE(0, curr_x - new_bound));
               assert(PIPSisZero(getSimpleVecFromColStochVec(v_vec, col)));
               slack = curr_x - new_bound;
            }

            /* adjust duals gamma and phi if necessary and use stored col to update */
            double& dual_row = row.inEqSys() ? getSimpleVecFromRowStochVec(y_vec, row) :
               getSimpleVecFromRowStochVec(z_vec, row);

            if(is_upper_bound)
            {
               double& phi = getSimpleVecFromColStochVec(phi_vec, col);
               if(!PIPSisZeroFeas(phi * slack))
               {
                  const double old_phi = phi;

                  const double coeff = row_storage.getRowCoefficientAtColumn( INDEX(ROW, row.getNode(), index_stored_row, row.getLinking(), EQUALITY_SYSTEM), col );
                  /* set z/y of corresponding row such that c_i* deltaz//a_i* deltay = phi */
                  assert(!PIPSisZero(coeff));
                  const double change_dual_row = old_phi/coeff;
                  if(coeff < 1e-13)
                     std::cout << "Potential numerical issues in postsolve of BoundTightening caused by small coefficient" << std::endl;

                  /* add z/y * row to gamma/phi */
                  StochVectorHandle tmp_pos(dynamic_cast<StochVector*>(gamma_vec.clone()));

                  /* calculate multiplier_change times row */
                  row_storage.axpyAtRow(0.0, *tmp_pos, change_dual_row, INDEX(ROW, row.getNode(), index_stored_row, row.getLinking(), EQUALITY_SYSTEM) );
                  StochVectorHandle tmp_neg(dynamic_cast<StochVector*>(tmp_pos->cloneFull()));

                  tmp_pos->selectPositive();
                  tmp_neg->selectNegative();

                  phi_vec.axpy(1.0, *tmp_pos);
                  gamma_vec.axpy(1.0, *tmp_neg);

                  dual_row += change_dual_row;

                  assert(PIPSisZero(phi));
               }
            }
            else
            {
               double& gamma = getSimpleVecFromColStochVec(gamma_vec, col);
               if(!PIPSisZeroFeas(gamma * slack))
               {
                  const double old_gamma = gamma;

                  const double coeff = row_storage.getRowCoefficientAtColumn( INDEX(ROW, row.getNode(), index_stored_row, row.getLinking(), EQUALITY_SYSTEM), col );
                  /* set z/y of corresponding row such that coeff*z//coeff*y = -gamma */
                  assert(!PIPSisZero(coeff));
                  const double change_dual_row = -old_gamma/coeff;

                  /* add z/y * row to gamma/phi */
                  StochVectorHandle tmp_pos(dynamic_cast<StochVector*>(gamma_vec.clone()));

                  /* calculate -multiplier_change times row */
                  row_storage.axpyAtRow(0.0, *tmp_pos, -change_dual_row, INDEX(ROW, row.getNode(), index_stored_row, row.getLinking(), EQUALITY_SYSTEM) );
                  StochVectorHandle tmp_neg(dynamic_cast<StochVector*>(tmp_pos->cloneFull()));

                  tmp_pos->selectPositive();
                  tmp_neg->selectNegative();

                  phi_vec.axpy(1.0, *tmp_pos);
                  gamma_vec.axpy(1.0, *tmp_neg);

                  dual_row += change_dual_row;

                  assert(PIPSisZero(gamma));
               }
            }
         }
         break;
      }
      case FIXED_COLUMN:
      {
         assert(first_index + 1 == next_first_index);
         assert(first_float_val + 2 == next_first_float_val);
         assert(first_int_val + 1 == next_first_int_val);

         const INDEX& idx_col = indices.at(first_index);
         assert(idx_col.isCol());

         const int column = idx_col.getIndex();
         const int node = idx_col.getNode();
         const int index_stored_col = int_values.at(first_int_val);
         const double value = float_values.at(first_float_val);
         const double obj_coeff = float_values.at(first_float_val + 1);

         assert(-1 <= node && node < static_cast<int>(x_vec.children.size()));
         assert(wasColumnRemoved(idx_col));

         /* mark entry as set and set x value to fixation */
         getSimpleVecFromColStochVec(*padding_origcol, node)[column] = 1;

         /* set x value */
         getSimpleVecFromColStochVec(x_vec, node)[column] = value;

         /* set slacks for x bounds to zero (bounds were tight) */
         getSimpleVecFromColStochVec(v_vec, node)[column] = 0.0;
         getSimpleVecFromColStochVec(w_vec, node)[column] = 0.0;

         /* set duals for bounds to satisfy reduced costs of reintroduced column times x */
         double col_times_duals = 0.0;
         if(node == -1)
         {
            /* we need to synchronize the column times duals in this case */
            if(my_rank == 0)
               col_times_duals = col_storage.multColTimesVec(INDEX(COL, node, index_stored_col), y_vec, z_vec);
            else
               col_times_duals = col_storage.multColTimesVecWithoutRootNode(INDEX(COL, node, index_stored_col), y_vec, z_vec);

            PIPS_MPIgetSumInPlace(col_times_duals, MPI_COMM_WORLD);
         }
         else
         {
            col_times_duals = col_storage.multColTimesVec( INDEX(COL, node, index_stored_col), y_vec, z_vec);
         }
         const double reduced_costs = obj_coeff - col_times_duals;

         /* set duals of bounds of x */
         double& gamma = getSimpleVecFromColStochVec(gamma_vec, node)[column];
         double& phi = getSimpleVecFromColStochVec(phi_vec, node)[column];

         gamma = 0.0;
         phi = 0.0;
         if( PIPSisLT(reduced_costs, 0.0) )
            phi = -reduced_costs;
         else if( PIPSisLT(0.0, reduced_costs) )
            gamma = reduced_costs;

         assert( PIPSisZeroFeas(reduced_costs - gamma + phi) );
         break;
      }
      case FIXED_EMPTY_COLUMN:
      {
         /**
          * recover primal value
          * dual multiplies will be set to zero
          * compute slack variables
          * no special treatment for linking variables since all processes should fix them simultaneously and in the same order
          *    -> assert is in place to check this
          */
         assert(first_index + 1 == next_first_index);
         assert(first_float_val + 4 == next_first_float_val);
         assert(first_int_val + 2 == next_first_int_val);

         const INDEX& idx_col = indices.at(first_index);
         assert(idx_col.isCol());

         const int column = idx_col.getIndex();
         const int node = idx_col.getNode();
         const double value = float_values.at(first_float_val);
         const double obj_coeff = float_values.at(first_float_val + 1);
         const double lbx = float_values.at(first_float_val + 2);
         const double ubx = float_values.at(first_float_val + 3);
         const int ixlow = int_values.at(first_int_val);
         const int ixupp = int_values.at(first_int_val + 1);
         assert(-1 <= node && node < static_cast<int>(x_vec.children.size()));
         assert(wasColumnRemoved(idx_col));

         if(node == -1)
            assert(PIPS_MPIisValueEqual(column, MPI_COMM_WORLD));

         /* primal */
         /* mark entry as set and set x value to fixation */
         getSimpleVecFromColStochVec(*padding_origcol, node)[column] = 1;
         getSimpleVecFromColStochVec(x_vec, node)[column] = value;

         if( ixlow == 1 )
            assert(PIPSisLEFeas(lbx, value));
         if( ixupp == 1 )
            assert(PIPSisLEFeas(value, ubx));

         /* dual */
         getSimpleVecFromColStochVec(gamma_vec, node)[column] = 0.0;
         getSimpleVecFromColStochVec(phi_vec, node)[column] = 0.0;

         if(!PIPSisZero(obj_coeff))
         {
            if( PIPSisLT(obj_coeff, 0.0) )
            {
               assert( ixupp );
               assert( PIPSisEQ(value, ubx) );
               getSimpleVecFromColStochVec(phi_vec, node)[column] = obj_coeff;
            }
            else if( PIPSisLT(0.0, obj_coeff) )
            {
               assert( ixlow );
               assert( PIPSisEQ(value, lbx) );
               getSimpleVecFromColStochVec(gamma_vec, node)[column] = obj_coeff;
            }
         }

         if( ixlow == 1 )
            getSimpleVecFromColStochVec(v_vec, node)[column] = value - lbx;
         else
            getSimpleVecFromColStochVec(v_vec, node)[column] = 0.0;

         if( ixupp == 1)
            getSimpleVecFromColStochVec(w_vec, node)[column] = ubx - value;
         else
            getSimpleVecFromColStochVec(w_vec, node)[column] = 0.0;

         assert(PIPSisZeroFeas(getSimpleVecFromColStochVec(v_vec, node)[column] * getSimpleVecFromColStochVec(gamma_vec, node)[column]));
         assert(PIPSisZeroFeas(getSimpleVecFromColStochVec(w_vec, node)[column] * getSimpleVecFromColStochVec(phi_vec, node)[column]));
         break;
      }
      case FIXED_COLUMN_SINGLETON_FROM_INEQUALITY:
      {
         assert(first_index + 1 == next_first_index);
         assert(first_float_val + 4 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

         const INDEX& col = indices.at(first_index);
         assert(col.isCol());

         const int node = col.getNode();
         const int col_index = col.getIndex();

         assert(-1 <= node && node < static_cast<int>(x_vec.children.size()));
         assert(wasColumnRemoved(col));

         const double value = float_values.at(first_float_val);
//         const double coeff = float_values.at(first_float_val + 1);
         const double xlow_old = float_values.at(first_float_val + 2);
         const double xupp_old = float_values.at(first_float_val + 3);

         /* set x value - bound is tight - compute slacks of x - set duals to zero */
         getSimpleVecFromColStochVec(*padding_origcol, node)[col_index] = 1;
         getSimpleVecFromColStochVec(x_vec, node)[col_index] = value;

         /* set slacks for x bounds */
         assert( PIPSisLE(value, xupp_old) || xupp_old == INF_POS_PRES);
         assert( PIPSisLE(xlow_old, value) || xlow_old == INF_NEG_PRES);

         getSimpleVecFromColStochVec(v_vec, node)[col_index] = (xlow_old == INF_NEG_PRES) ? 0 : value - xlow_old;
         getSimpleVecFromColStochVec(w_vec, node)[col_index] = (xupp_old == INF_POS_PRES) ? 0 : xupp_old - value;

         /* TODO: don't think something needs to be done for the reduced costs */
         getSimpleVecFromColStochVec(phi_vec, node)[col_index] = 0.0;
         getSimpleVecFromColStochVec(gamma_vec, node)[col_index] = 0.0;

         /* set duals of bounds of x */
         getSimpleVecFromColStochVec(gamma_vec, node)[col_index] = 0.0;
         getSimpleVecFromColStochVec(phi_vec, node)[col_index] = 0.0;

         break;
      }
      case PARALLEL_COLUMN:
      {
         throw std::runtime_error("PARALLEL_COLUMN not yet implemented");
         break;
      }
      case DELETED_ROW:
      {
         throw std::runtime_error("DELETED_ROW not yet implemented");
         break;
      }
      case SINGLETON_EQUALITY_ROW:
      {
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 5 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

         const INDEX& row = indices.at(first_index);
         const INDEX& col = indices.at(first_index + 1);

         assert(row.isRow());
         assert(row.getSystemType() == EQUALITY_SYSTEM);
         assert(col.isCol());

         /* row should have been re-added by redundant row event by now - same for column by fix column event */
         assert(!wasRowRemoved(row));
         assert(!wasColumnRemoved(col));

         const double xlow_old = float_values.at(first_float_val);
         const double xupp_old = float_values.at(first_float_val + 1);
         const double coeff = float_values.at(first_float_val + 4);

         assert(!PIPSisZero(coeff));

         const double curr_x = getSimpleVecFromColStochVec(x_vec, col.getNode())[col.getIndex()];
         double& slack_lower = getSimpleVecFromColStochVec(v_vec, col.getNode())[col.getIndex()];
         double& slack_upper = getSimpleVecFromColStochVec(w_vec, col.getNode())[col.getIndex()];

         /* if bound is not tight only adjust slack v/w */
         if(xupp_old == INF_POS_PRES)
            slack_upper = 0;
         else
         {
            slack_upper = xupp_old - curr_x;
            assert(PIPSisLT(0.0, slack_upper));
            assert(std::fabs(slack_upper) != INF_POS_PRES);
         }

         if(xlow_old == INF_NEG_PRES)
            slack_lower = 0;
         else
         {
            slack_lower = curr_x - xlow_old;
            assert(PIPSisLE(0.0, slack_lower));
            assert(std::fabs(slack_lower) != INF_POS_PRES);
         }

         double& dual_singelton_row = getSimpleVecFromRowStochVec(y_vec, row.getNode(), row.getLinking())[row.getIndex()];
         double& dual_lower = getSimpleVecFromColStochVec(gamma_vec, col.getNode())[col.getIndex()];
         double& dual_upper = getSimpleVecFromColStochVec(phi_vec, col.getNode())[col.getIndex()];

         double error_in_reduced_costs = 0.0;

         /* adjust duals if necessary */
         if( !PIPSisZero(slack_lower) )
         {
            error_in_reduced_costs -= dual_lower;
            dual_lower = 0.0;
         }

         if( !PIPSisZero(slack_upper) )
         {
            error_in_reduced_costs += dual_upper;
            dual_upper = 0.0;
         }

         if(!PIPSisZero(error_in_reduced_costs))
         {
            /* we use the coeff*dual_singleton_row to balance the error in the reduced costs */
            double diff_dual_row = error_in_reduced_costs / coeff;
            dual_singelton_row += diff_dual_row;
            assert(PIPSisEQ(diff_dual_row * coeff, error_in_reduced_costs));
         }
         break;
      }
      case SINGLETON_INEQUALITY_ROW:
      {
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 5 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

         const INDEX& row = indices.at(first_index);
         const INDEX& col = indices.at(first_index + 1);

         assert(row.isRow());
         assert(row.getSystemType() == INEQUALITY_SYSTEM);
         assert(col.isCol());

         const double xlow_old = float_values.at(first_float_val);
         const double xupp_old = float_values.at(first_float_val + 1);
         const double xlow_new = float_values.at(first_float_val + 2);
         const double xupp_new = float_values.at(first_float_val + 3);
         const double coeff = float_values.at(first_float_val + 4);

         assert(!PIPSisZero(coeff));
         assert(xlow_new == INF_NEG_PRES || xupp_new == INF_POS_PRES);
         assert(xlow_new != INF_NEG_PRES || xupp_new != INF_POS_PRES);

         /* should have been re-added by now */
         assert( !wasColumnRemoved(col) );
         assert( !wasRowRemoved(row) );

         bool lower_bound_changed = (xlow_new != INF_NEG_PRES);

         const double curr_x = getSimpleVecFromColStochVec(x_vec, col.getNode())[col.getIndex()];

         double& slack = lower_bound_changed ? getSimpleVecFromColStochVec(v_vec, col.getNode())[col.getIndex()] : getSimpleVecFromColStochVec(w_vec, col.getNode())[col.getIndex()];
         const double old_bound = (xupp_new != INF_POS_PRES) ? xupp_old : xlow_old;

         /* if bound is not tight only adjust slack v/w */
         if(std::fabs(old_bound) == INF_POS_PRES)
            slack = 0;
         else
         {
            slack = std::fabs(old_bound - curr_x);
            assert(PIPSisLT(0.0, slack));
            assert(std::fabs(slack) != INF_POS_PRES);
         }

         double& dual_singelton_row = getSimpleVecFromRowStochVec(z_vec, row.getNode(), row.getLinking())[row.getIndex()];
         double& dual_bound = lower_bound_changed ? getSimpleVecFromColStochVec(gamma_vec, col.getNode())[col.getIndex()] : getSimpleVecFromColStochVec(phi_vec, col.getNode())[col.getIndex()];

         double error_in_reduced_costs = 0.0;
         if(PIPSisZero(slack))
         {
            error_in_reduced_costs = lower_bound_changed ? -dual_bound : dual_bound;
            dual_bound = 0.0;
         }

         if(!PIPSisZero(error_in_reduced_costs))
         {
            /* we use the coeff*dual_singleton_row to balance the error in th reduced costs */
            double diff_dual_row = error_in_reduced_costs / coeff;
            dual_singelton_row += diff_dual_row;
            assert(PIPSisEQ(diff_dual_row * coeff, error_in_reduced_costs));
         }
         break;
      }
      case FREE_COLUMN_SINGLETON_EQUALITY:
      {
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 5 == next_first_float_val);
         assert(first_int_val + 1 == next_first_int_val);

         const INDEX& col = indices.at(first_index);
         const INDEX& row = indices.at(first_index + 1);

         assert(row.isRow());

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
         /* set duals of row depending on equality/inequality row */
         if( row.inEqSys() )
         {
            getSimpleVecFromRowStochVec(y_vec, row) = dual_value_row;
            assert( PIPSisZero(obj_coeff - getSimpleVecFromRowStochVec(y_vec, row) * col_coeff) );
         }
         else
         {
            /* cupp == clow so slacks must be zero */
            getSimpleVecFromRowStochVec(s_vec, row) = 0.0;
            getSimpleVecFromRowStochVec(t_vec, row) = 0.0;
            getSimpleVecFromRowStochVec(u_vec, row) = 0.0;

            /* set dual */
            getSimpleVecFromRowStochVec(z_vec, row) = dual_value_row;
            if( PIPSisLT(0.0, dual_value_row) )
               getSimpleVecFromRowStochVec(pi_vec, row) = dual_value_row;
            else
               getSimpleVecFromRowStochVec(lambda_vec, row) = -dual_value_row;
         }

         /* synchronize value of row for x_val */
         double value_row = 0.0;

         if( row.isLinkingRow() )
         {
            assert(PIPS_MPIisValueEqual(row.getIndex(), MPI_COMM_WORLD));

            if(my_rank == 0)
               value_row = row_storage.multRowTimesVec( stored_row, x_vec );
            else
               value_row = row_storage.multLinkingRowTimesVecWithoutBl0( stored_row_idx, x_vec);
            PIPS_MPIgetSumInPlace(value_row, MPI_COMM_WORLD);
         }
         else
            value_row = row_storage.multRowTimesVec(stored_row, x_vec);

         assert(std::abs(value_row) != INF_POS_PRES);

         /* reintroduce the removed column on process owning column */
         if( col.isCol() )
         {
            double& x_val = getSimpleVecFromColStochVec(x_vec, col);
            assert(PIPSisZero(x_val));

            /* mark column as set */
            getSimpleVecFromColStochVec(*padding_origcol, col) = 1;

            /* recover primal value */
            x_val = (rhs - value_row) / col_coeff;
            assert( PIPSisZero(x_val * col_coeff + value_row - rhs) );

            /* compute slacks and set duals for bounds to zero */
            double& slack_lower = getSimpleVecFromColStochVec(v_vec, col);
            double& slack_upper = getSimpleVecFromColStochVec(w_vec, col);

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

            getSimpleVecFromColStochVec(gamma_vec, col) = 0.0;
            getSimpleVecFromColStochVec(phi_vec, col) = 0.0;
         }

         break;
      }
      case NEARLY_PARALLEL_ROW_SUBSTITUTION:
      {
         /* col2 was substituted by col1 via col2 = t * col1 + d */
         assert(first_index + 4 == next_first_index);
         assert(first_float_val + 8 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

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

         const double val_col1 = col1.isCol() ? getSimpleVecFromColStochVec(x_vec, col1) : 0.0;
         const double val_col2 = scalar * val_col1 + translation;

         /* primal postsolve */

         /* reintroduce the substituted column */
         getSimpleVecFromColStochVec(*padding_origcol, col2) = 1;
         getSimpleVecFromColStochVec(x_vec, col2) = val_col2;
         assert( PIPSisLEFeas( xlow_col2, val_col2 ) );
         assert( PIPSisLEFeas( val_col2, xupp_col2 ) );

         /* slacks for bounds */
         getSimpleVecFromColStochVec(v_vec, col2) = (xlow_col2 == INF_NEG_PRES) ? 0 : val_col2 - xlow_col2;
         getSimpleVecFromColStochVec(w_vec, col2) = (xupp_col2 == INF_POS_PRES) ? 0 : xupp_col2 - val_col2;

         /* bound duals to zero */
         getSimpleVecFromColStochVec(phi_vec, col2) = 0.0;
         getSimpleVecFromColStochVec(gamma_vec, col2) = 0.0;

         /* dual postsolve substitution */

         /* the substitution itself needs a shift of the row duals to compensate the objective vector change or, if the row was not tight (so only for inequalities)
          * a shift of the column duals
          */

         /* postsolve two equalities */
         if( row1.inEqSys() )
         {
            double& dual_row1 = getSimpleVecFromRowStochVec(y_vec, row1);
            double& dual_row2 = getSimpleVecFromRowStochVec(y_vec, row2);

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
            double& z_row1 = getSimpleVecFromRowStochVec(z_vec, row1);
            double& lambda_row1 = getSimpleVecFromRowStochVec(lambda_vec, row1);
            double& pi_row1 = getSimpleVecFromRowStochVec(pi_vec, row1);

            double& z_row2 = getSimpleVecFromRowStochVec(z_vec, row2);
            double& lambda_row2 = getSimpleVecFromRowStochVec(lambda_vec, row2);
            double& pi_row2 = getSimpleVecFromRowStochVec(pi_vec, row2);

            /* z = lambda - pi */
            z_row2 = (1 - change_factor_obj_row1) * parallel_factor * z_row1;
            pi_row2 = (1 - change_factor_obj_row1) * parallel_factor * pi_row1;
            lambda_row2 = (1 - change_factor_obj_row1) * parallel_factor * lambda_row1;

            z_row1 *= change_factor_obj_row1;
            pi_row1 *= change_factor_obj_row1;
            lambda_row1 *= change_factor_obj_row1;

            /* shift bound duals */
            double& gamma_row1 = getSimpleVecFromRowStochVec(gamma_vec, row1);
            double& phi_row1 = getSimpleVecFromRowStochVec(phi_vec, row1);

            double& gamma_row2 = getSimpleVecFromRowStochVec(gamma_vec, row2);
            double& phi_row2 = getSimpleVecFromRowStochVec(phi_vec, row2);

            gamma_row2 = (1 - change_factor_obj_row1) * gamma_row1;
            phi_row2 = (1 - change_factor_obj_row1) * phi_row1;

            gamma_row1 *= change_factor_obj_row1;
            phi_row1 *= change_factor_obj_row1;


#ifndef NDEBUG
            const double& v_row1 = getSimpleVecFromColStochVec(v_vec, col1);
            const double& v_row2 = getSimpleVecFromColStochVec(v_vec, col2);
            const double& w_row1 = getSimpleVecFromColStochVec(w_vec, col1);
            const double& w_row2 = getSimpleVecFromColStochVec(w_vec, col2);

            const double& t_row1 = getSimpleVecFromColStochVec(t_vec, col1);
            const double& t_row2 = getSimpleVecFromColStochVec(t_vec, col2);
            const double& u_row1 = getSimpleVecFromColStochVec(u_vec, col1);
            const double& u_row2 = getSimpleVecFromColStochVec(u_vec, col2);
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

         break;
      }
      case NEARLY_PARALLEL_ROW_BOUNDS_TIGHTENED:
      {
         /* col2 was substituted by col1 via col2 = t * col1 + d */
         assert(first_index + 4 == next_first_index);
         assert(first_float_val + 12 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

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

         const double val_col1 = getSimpleVecFromColStochVec(x_vec, col1);
         const double val_col2 = col2.isCol() ? getSimpleVecFromColStochVec(x_vec, col2) : 0.0;

         /* if the variable bound of col1 was actually implied via col2 we have to shift it's dual multipliers over via also adjusting the dual of row1 */
         assert( PIPSisRelLEFeas( xlow_col1, val_col1 ) );
         assert( PIPSisRelLEFeas( val_col1, xupp_col1 ) );

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
            double& dual_col1 = getSimpleVecFromColStochVec(gamma_vec, col1);

            if( !PIPSisZero(dual_col1) )
            {
               /* set the dual of col1 to zero, compensate for the error made with the dual of row1, compensate that error with the dual of col2
                * and finally compensate the error in col2's row of the reduced costs with col2's dual
                */

               /* only one bound dual variable is non-zero */
               assert( PIPSisZero(getSimpleVecFromColStochVec(phi_vec, col1)) );

               /* dual to zero and compensate error */
               const double dual_shift_row1 = dual_col1 / coeff_col1;
               dual_col1 = 0;
               assert( PIPSisEQ( -coeff_col1 * dual_shift_row1, -dual_col1 ) );

               getSimpleVecFromRowStochVec(y_vec, row1) += dual_shift_row1;

               /* compensate again and adjust dual of col2 */
               const double dual_shift_row2 = -dual_shift_row1 * parallel_factor;

               if( row2.inEqSys() )
                  getSimpleVecFromRowStochVec(y_vec, row2) += dual_shift_row2;
               else
               {
                  getSimpleVecFromRowStochVec(z_vec, row2) += dual_shift_row2;
                  getSimpleVecFromRowStochVec(lambda_vec, row2) = std::max(0.0, getSimpleVecFromRowStochVec(z_vec, row2));
                  getSimpleVecFromRowStochVec(pi_vec, row2) += std::min(0.0, getSimpleVecFromRowStochVec(z_vec, row2));
               }

               /* if we had two equality rows there was a singleton col in row2 an it's dual needs to be adjusted */
               if( col2.isCol() )
               {
                  /* assert bounds of substituted variable are tight as well */
                  const double implying_bound = PIPSisLT(0.0, translation) ? xupp_col2 : xlow_col2;
                  assert( PIPSisEQ(val_col2, implying_bound) );

                  const double dual_shift_col2 = coeff_col2 * dual_shift_row2;
                  PIPSisLT(0.0, translation) ? assert( PIPSisLE(0.0, dual_shift_col2) ) : assert( PIPSisLE( dual_shift_col2, 0.0) );

                  double& dual_col2 = PIPSisLT(0.0, translation) ? getSimpleVecFromColStochVec(phi_vec, col2) : getSimpleVecFromColStochVec(gamma_vec, col2);
                  assert( PIPSisZero(dual_col2) );
                  dual_col2 += std::fabs(dual_shift_col2);
               }
            }
         }

         /* upper bound was implied by substituted column + it's row */
         if( PIPSisEQ( val_col1, xupp_implied ) && !PIPSisEQ( xupp_implied, xupp_col1 ) )
         {
            assert( !PIPSisEQ(val_col1, xupp_col1) );
            double& dual_col1 = getSimpleVecFromColStochVec(phi_vec, col1);

            if( !PIPSisZero(dual_col1) )
            {
               /* set the dual of col1 to zero, compensate for the error made with the dual of row1, compensate that error with the dual of col2
                * and finally compensate the error in col2's row of the reduced costs with col2's dual
                */

               /* only one bound dual variable is non-zero */
               assert( PIPSisZero(getSimpleVecFromColStochVec(gamma_vec, col1)) );

               /* dual to zero and compensate error */
               const double dual_shift_row1 = dual_col1 / coeff_col1;
               dual_col1 = 0;
               assert( PIPSisEQ( -coeff_col1 * dual_shift_row1, -dual_col1 ) );

               getSimpleVecFromRowStochVec(y_vec, row1) += dual_shift_row1;

               /* compensate again and adjust dual of col2 */
               const double dual_shift_row2 = -dual_shift_row1 * parallel_factor;

               if( row2.inEqSys() )
                  getSimpleVecFromRowStochVec(y_vec, row2) += dual_shift_row2;
               else
               {
                  getSimpleVecFromRowStochVec(z_vec, row2) += dual_shift_row2;
                  getSimpleVecFromRowStochVec(lambda_vec, row2) = std::max(0.0, getSimpleVecFromRowStochVec(z_vec, row2));
                  getSimpleVecFromRowStochVec(pi_vec, row2) += std::min(0.0, getSimpleVecFromRowStochVec(z_vec, row2));
               }

               if( col2.isCol() )
               {
                  /* assert bounds of substituted variable are tight as well */
                  const double implying_bound = PIPSisLT(0.0, translation) ? xlow_col2 : xupp_col2;
                  assert( PIPSisEQ(val_col2, implying_bound) );

                  double& dual_col2 = PIPSisLT(0.0, translation) ? getSimpleVecFromColStochVec(gamma_vec, col2) : getSimpleVecFromColStochVec(phi_vec, col2);
                  assert( PIPSisZero(dual_col2) );

                  const double dual_shift_col2 = coeff_col2 * dual_shift_row2;

                  PIPSisLT(0.0, translation) ? assert( PIPSisLE(0.0, dual_shift_col2) ) : assert( PIPSisLE( dual_shift_col2, 0.0) );
                  dual_col2 += std::fabs(dual_shift_col2);
               }
            }
         }

         break;
      }
      case LINKING_VARS_SYNC_EVENT:
      {
         // todo : dual part of this
         const int length_link_vars = x_vec.vec->length();
         SimpleVector &link_vars = dynamic_cast<SimpleVector&>(*x_vec.vec);

         double *copy_x_link_max = new double[length_link_vars];
         double *copy_x_link_min = new double[length_link_vars];
         std::copy(link_vars.elements(), link_vars.elements() + length_link_vars, copy_x_link_max);
         std::copy(link_vars.elements(), link_vars.elements() + length_link_vars, copy_x_link_min);

         PIPS_MPIminArrayInPlace(copy_x_link_min, length_link_vars, MPI_COMM_WORLD);
         PIPS_MPImaxArrayInPlace(copy_x_link_max, length_link_vars, MPI_COMM_WORLD);

         /* changing vars must have been set to 0 ! */
         for( int j = 0; j < length_link_vars; ++j )
         {
            /* the second check is necessary for things like +- inf used by the postsolver */
            if( !PIPSisEQ(copy_x_link_min[j], copy_x_link_max[j]) && copy_x_link_max[j] != copy_x_link_min[j] )
            {
               assert(PIPSisZero(copy_x_link_max[j]) || PIPSisZero(copy_x_link_min[j]));
               assert(PIPSisEQ(link_vars[j], copy_x_link_min[j]) || PIPSisEQ(link_vars[j], copy_x_link_max[j]));

               if( !PIPSisZero(copy_x_link_min[j]) )
                  link_vars[j] = copy_x_link_min[j];
               else
                  link_vars[j] = copy_x_link_max[j];
            }
         }

         break;
      }
      case FREE_COLUMN_SINGLETON_INEQUALITY_ROW:
      {
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 2 == next_first_float_val);
         assert(first_int_val + 1 == next_first_int_val);

         const INDEX& row = indices.at(first_index);
         const INDEX& col = indices.at(first_index + 1);
         assert( row.inInEqSys() );

         const int index_stored_row = int_values.at(first_int_val);
         const INDEX row_stored(ROW, row.getNode(), index_stored_row, row.getLinking(), row.getSystemType());

         const double lhsrhs = float_values.at(first_float_val);
         const double coeff = float_values.at(first_float_val + 1);

         assert( wasColumnRemoved(col) );
         assert( wasRowRemoved(row) );

         getSimpleVecFromColStochVec(*padding_origcol, col) = 1;

         getSimpleVecFromRowStochVec(*padding_origrow_inequality, row) = 1;

         /* set x value such that row is satisfied */
         double& x_val = getSimpleVecFromColStochVec(x_vec, col);
         assert( PIPSisZero(x_val) );

         x_val = (lhsrhs - row_storage.multRowTimesVec(row_stored, x_vec)) / coeff;

         /* duals of row and bounds are zero */
         getSimpleVecFromRowStochVec(z_vec, row) = 0;
         getSimpleVecFromRowStochVec(lambda_vec, row) = 0;
         getSimpleVecFromRowStochVec(pi_vec, row) = 0;

         getSimpleVecFromColStochVec(gamma_vec, col) = 0;
         getSimpleVecFromColStochVec(phi_vec, col) = 0;

         /* compute slacks for bounds and row */
         /* slack for row is zero by construction */
         getSimpleVecFromRowStochVec(t_vec, row) = 0;
         getSimpleVecFromRowStochVec(u_vec, row) = 0;
         getSimpleVecFromRowStochVec(s_vec, row) = 0;

         const int ixlow = getSimpleVecFromColStochVec( *orig_problem_s.ixlow, col);
         const int ixupp = getSimpleVecFromColStochVec( *orig_problem_s.ixupp, col);
         const double xlow = getSimpleVecFromColStochVec( *orig_problem_s.blx, col);
         const double xupp = getSimpleVecFromColStochVec( *orig_problem_s.bux, col);


         if(PIPSisZero(ixlow))
            getSimpleVecFromColStochVec(v_vec, col) = 0;
         else
         {
            getSimpleVecFromColStochVec(v_vec, col) = x_val - xlow;
            assert( PIPSisLE(xlow, x_val) );
         }

         if(PIPSisZero(ixupp))
            getSimpleVecFromColStochVec(w_vec, col) = 0;
         else
         {
            getSimpleVecFromColStochVec(w_vec, col) = xupp - x_val;
            assert( PIPSisLE(x_val, xupp) );
         }

         getSimpleVecFromColStochVec(gamma_vec, col) = 0.0;
         getSimpleVecFromColStochVec(phi_vec, col) = 0.0;

         break;
      }
      case PARALLEL_ROWS_BOUNDS_TIGHTENED:
      {
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 5 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

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
         double& z_row1 = getSimpleVecFromRowStochVec(z_vec, row1);
         double& lambda_row1 = getSimpleVecFromRowStochVec(lambda_vec, row1);
         double& pi_row1 = getSimpleVecFromRowStochVec(pi_vec, row1);

         double& t_row1 = getSimpleVecFromRowStochVec(t_vec, row1);
         double& u_row1 = getSimpleVecFromRowStochVec(u_vec, row1);

         double& z_row2 = getSimpleVecFromRowStochVec(z_vec, row2);
         double& lambda_row2 = getSimpleVecFromRowStochVec(lambda_vec, row2);
         double& pi_row2 = getSimpleVecFromRowStochVec(pi_vec, row2);

#ifndef NDEBUG
         const double& t_row2 = PIPSisLT(factor, 0.0) ? getSimpleVecFromRowStochVec(u_vec, row2) : getSimpleVecFromRowStochVec(t_vec, row2);
         const double& u_row2 = PIPSisLT(factor, 0.0) ? getSimpleVecFromRowStochVec(t_vec, row2) : getSimpleVecFromRowStochVec(u_vec, row2);
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

         break;
      }
      default:
      {
         throw std::runtime_error("Tried to postsolve not supported reduction type");
         break;
      }

         // todo : an is-set Stochvector? - rather use padding StochVectors -- or check padding vectors for all ones?
      }
   }

   /* compute all s, t and u that have not yet been computed */

#ifndef NDEBUG
   /* assert that all primal variables have been set */
   double absmin_x;
   x_vec.absmin(absmin_x);
   assert( absmin_x < std::numeric_limits<double>::max() ); 
#endif
   assert( x_vec.isRootNodeInSync() );
   assert( s_vec.isRootNodeInSync() );
   assert( y_vec.isRootNodeInSync() );
   assert( z_vec.isRootNodeInSync() );
   assert( v_vec.isRootNodeInSync() );
   assert( gamma_vec.isRootNodeInSync() );
   assert( w_vec.isRootNodeInSync() );
   assert( phi_vec.isRootNodeInSync() );
   assert( t_vec.isRootNodeInSync() );
   assert( lambda_vec.isRootNodeInSync() );
   assert( u_vec.isRootNodeInSync() );
   assert( pi_vec.isRootNodeInSync() );

   if( my_rank == 0 )
      std::cout << "finished postsolving... " << std::endl;

   return PRESOLVE_OK;
}

void StochPostsolver::setOriginalVarsFromReduced(const sVars& reduced_vars, sVars& original_vars) const
{
   /* x */
   const StochVector &x_reduced = dynamic_cast<const StochVector&>(*reduced_vars.x);
   StochVector &x_orig = dynamic_cast<StochVector&>(*original_vars.x);
   setOriginalValuesFromReduced<>(x_orig, x_reduced, *padding_origcol);

   /* s */
   const StochVector &s_reduced = dynamic_cast<const StochVector&>(*reduced_vars.s);
   StochVector &s_orig = dynamic_cast<StochVector&>(*original_vars.s);
   setOriginalValuesFromReduced(s_orig, s_reduced, *padding_origrow_inequality);

   /* y */
   const StochVector &y_reduced = dynamic_cast<const StochVector&>(*reduced_vars.y);
   StochVector &y_orig = dynamic_cast<StochVector&>(*original_vars.y);
   setOriginalValuesFromReduced(y_orig, y_reduced, *padding_origrow_equality);

   /* z */
   const StochVector &z_reduced = dynamic_cast<const StochVector&>(*reduced_vars.z);
   StochVector &z_orig = dynamic_cast<StochVector&>(*original_vars.z);
   setOriginalValuesFromReduced(z_orig, z_reduced, *padding_origrow_inequality);

   /* v */
   const StochVector &v_reduced = dynamic_cast<const StochVector&>(*reduced_vars.v);
   StochVector &v_orig = dynamic_cast<StochVector&>(*original_vars.v);
   setOriginalValuesFromReduced(v_orig, v_reduced, *padding_origcol);

   /* gamma */
   const StochVector &gamma_reduced = dynamic_cast<const StochVector&>(*reduced_vars.gamma);
   StochVector &gamma_orig = dynamic_cast<StochVector&>(*original_vars.gamma);
   setOriginalValuesFromReduced(gamma_orig, gamma_reduced, *padding_origcol);

   /* w */
   const StochVector &w_reduced = dynamic_cast<const StochVector&>(*reduced_vars.w);
   StochVector &w_orig = dynamic_cast<StochVector&>(*original_vars.w);
   setOriginalValuesFromReduced(w_orig, w_reduced, *padding_origcol);

   /* phi */
   const StochVector &phi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.phi);
   StochVector &phi_orig = dynamic_cast<StochVector&>(*original_vars.phi);
   setOriginalValuesFromReduced(phi_orig, phi_reduced, *padding_origcol);

   /* t */
   const StochVector &t_reduced = dynamic_cast<const StochVector&>(*reduced_vars.t);
   StochVector &t_orig = dynamic_cast<StochVector&>(*original_vars.t);
   setOriginalValuesFromReduced(t_orig, t_reduced, *padding_origrow_inequality);

   /* lambda */
   const StochVector &lambda_reduced = dynamic_cast<const StochVector&>(*reduced_vars.lambda);
   StochVector &lambda_orig = dynamic_cast<StochVector&>(*original_vars.lambda);
   setOriginalValuesFromReduced(lambda_orig, lambda_reduced, *padding_origrow_inequality);

   /* u */
   const StochVector &u_reduced = dynamic_cast<const StochVector&>(*reduced_vars.u);
   StochVector &u_orig = dynamic_cast<StochVector&>(*original_vars.u);
   setOriginalValuesFromReduced(u_orig, u_reduced, *padding_origrow_inequality);

   /* pi */
   const StochVector &pi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.pi);
   StochVector &pi_orig = dynamic_cast<StochVector&>(*original_vars.pi);
   setOriginalValuesFromReduced(pi_orig, pi_reduced, *padding_origrow_inequality);
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
