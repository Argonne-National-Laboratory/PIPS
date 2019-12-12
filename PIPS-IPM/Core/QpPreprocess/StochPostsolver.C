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
   stored_rows( dynamic_cast<const StochGenMatrix&>(*original_problem.A).cloneEmptyRows(true) )
{
   assert(stored_rows->children.size() == dynamic_cast<const StochGenMatrix&>(*original_problem.A).children.size() );

   /* 1 indicates that the row / col has not been removed from the problem - -1 indicates the row / col has been removed */
   padding_origcol->setToConstant(1);
   padding_origrow_equality->setToConstant(1);
   padding_origrow_inequality->setToConstant(1);

   // count !? rows cols...
   start_idx_float_values.push_back(0);
   start_idx_int_values.push_back(0);
   start_idx_indices.push_back(0);
}

StochPostsolver::~StochPostsolver()
{
   delete padding_origrow_inequality;
   delete padding_origrow_equality;
   delete padding_origcol;
}

void StochPostsolver::notifySingletonEqualityRow( int node, int row, BlockType block_type, int col, double coeff, double rhs)
{
   throw std::runtime_error("Not yet implemented");
}

void StochPostsolver::notifySingletonIneqalityRow( int node, int row, BlockType block_type, int col, double coeff, double lhs, double rhs )
{
   throw std::runtime_error("Not yet implemented");
}

void StochPostsolver::notifyFreeColumnSingleton( SystemType system_type, int node_row, int row, bool linking_row, double rhs, int node_col, int col, 
   const StochGenMatrix& matrix_row )
{
   //todo INEQUALITY system?  
   assert(system_type == EQUALITY_SYSTEM);
   assert(!wasRowRemoved(system_type, node_row, row, linking_row));
   assert(!wasColumnRemoved(node_col, col));

   markRowRemoved(system_type, node_row, row, linking_row);
   
   // save dual postsolve info
   // save row for primal postsolve info
   int stored_row_idx = stored_rows->appendRow(matrix_row, node_row, row, linking_row);

   indices.push_back(INDEX(COL, node_col, col));
   indices.push_back(INDEX(ROW, node_row, row, linking_row, system_type));

   float_values.push_back(rhs);
   int_values.push_back(stored_row_idx);
   reductions.push_back( FREE_COLUMN_SINGLETON );

   finishNotify();
}

/* substitute var2 with scalar*var1 + translation */
void StochPostsolver::notifyParallelRowSubstitution(SystemType system_type, int node_row, int var1, int row1, int node_var1, int var2, int row2, 
   int node_var2, double scalar, double translation)
{
   // todo : linking rows are not possible yet
   assert( !wasColumnRemoved(node_var1, var1) );
   assert( !wasColumnRemoved(node_var2, var2) );

   assert( !wasRowRemoved(system_type, node_row, row1, false) );
   assert( !wasRowRemoved(system_type, node_row, row2, false) );

   reductions.push_back(PARALLEL_ROW_SUBSTITUTION);

   indices.push_back( INDEX(COL, node_var2, var2) );
   indices.push_back( INDEX(COL, node_var1, var1) );

   markColumnRemoved(node_var2, var2);

   float_values.push_back( scalar );
   float_values.push_back( translation );

   finishNotify();
}

/** postsolve has to compute the optimal dual multipliers here and set the primal value accordingly */
void StochPostsolver::notifyFixedColumn( int node, unsigned int col, double value, const StochGenMatrix& eq_mat, const StochGenMatrix& ineq_mat)
{
   assert( !wasColumnRemoved(node, col) );
   markColumnRemoved(node, col);

   // todo : dual postsolve
   // todo : add matrix to store columns
   // todo : store column

   /* store current upper and lower bounds of x and the local column */
   reductions.push_back(FIXED_COLUMN);
   indices.push_back(INDEX(COL, node, col));
   float_values.push_back(value);

   finishNotify();
}


void StochPostsolver::notifyFixedEmptyColumn(int node, unsigned int col, double value, double obj_value, int ixlow, int ixupp, double lbx, double ubx)
{
   assert(!wasColumnRemoved(node, col));
   markColumnRemoved(node, col);

   // todo : ?
   assert(std::fabs(value) < 1e10);

   reductions.push_back(FIXED_EMPTY_COLUMN);
   indices.push_back(INDEX(COL,node, col));
   float_values.push_back(value);
   float_values.push_back(obj_value);
   float_values.push_back(lbx);
   float_values.push_back(ubx);
   int_values.push_back(ixlow);
   int_values.push_back(ixupp);

   finishNotify();
}

/** postsolve for this is simply to set all dual variables to zero - the row itself has no primal impact */
void StochPostsolver::notifyRedundantRow( SystemType system_type, int node, unsigned int row, bool linking_row,
   int iclow, int icupp, double lhs, double rhs, const StochGenMatrix& matrix_row )
{
   assert(iclow == 1 || iclow == 0);
   assert(icupp == 1 || icupp == 0);
   assert(iclow + icupp > 0);

   assert(!wasRowRemoved(system_type, node, linking_row, row));
   markRowRemoved(system_type, node, row, linking_row);

   /* save row for postsolve */
   reductions.push_back(REDUNDANT_ROW);
   indices.push_back(INDEX(ROW, node, row, linking_row, system_type));

   int index_stored_row = stored_rows->appendRow(matrix_row, node, row, linking_row);

   float_values.push_back(lhs);
   float_values.push_back(rhs);
   int_values.push_back(index_stored_row);
   int_values.push_back(iclow);
   int_values.push_back(icupp);

   finishNotify();
}

// todo : only store each version of each row once!
// todo : store whole row
void StochPostsolver::notifyRowPropagated( SystemType system_type, int node, int row, bool linking_constraint,
      int column, double lb, double ub, double* values_row, int* indices_row, int length)
{
   return;

   // finishNotify();
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

bool StochPostsolver::wasColumnRemoved(int node, int col) const
{
   return getSimpleVecFromColStochVec(*padding_origcol, node)[col] == -1;
}

void StochPostsolver::markColumnRemoved(int node, int col)
{
   assert(!wasColumnRemoved(node, col));
   getSimpleVecFromColStochVec(*padding_origcol, node)[col] = -1;
}

void StochPostsolver::markColumnAdded(int node, int col)
{
   assert(wasColumnRemoved(node, col));
   getSimpleVecFromColStochVec(*padding_origcol, node)[col] = 1;
}

bool StochPostsolver::wasRowRemoved(SystemType system_type, int node, int row, bool linking_row) const
{
   if(system_type == EQUALITY_SYSTEM)
      return getSimpleVecFromRowStochVec(*padding_origrow_equality, node, linking_row)[row] == -1;
   else
      return getSimpleVecFromRowStochVec(*padding_origrow_inequality, node, linking_row)[row] == -1;
}

void StochPostsolver::markRowRemoved(SystemType system_type, int node, int row, bool linking_row)
{
   assert(!wasRowRemoved(system_type, node, row, linking_row));
   if(system_type == EQUALITY_SYSTEM)
      getSimpleVecFromRowStochVec(*padding_origrow_equality, node, linking_row)[row] = -1;
   else
      getSimpleVecFromRowStochVec(*padding_origrow_inequality, node, linking_row)[row] = -1;
}

// todo use somehow?
void StochPostsolver::markRowAdded(SystemType system_type, int node, int row, bool linking_row)
{
   assert(wasRowRemoved(system_type, node, row, linking_row));
   if(system_type == EQUALITY_SYSTEM)
      getSimpleVecFromRowStochVec(*padding_origrow_equality, node, linking_row)[row] = 1;
   else
      getSimpleVecFromRowStochVec(*padding_origrow_inequality, node, linking_row)[row] = 1;
}

// todo : usage and check of padding origrow - can already be done - even without any dual postsolve stuff
// todo : sort reductions by nodes ? and then reverse ?
// todo : at the moment only replaces whatever is given as x with the x solution in the original soution space
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

   // todo
   /* dual solution is now reduced solution padded with zeros */

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
         /* only dual postsolve */
         assert(first_index + 1 == next_first_index);
         assert(first_float_val + 2 == next_first_float_val);
         assert(first_int_val + 3 == next_first_int_val);

         const INDEX& row_idx = indices.at(first_index);
         assert(row_idx.index_type == ROW);

         const int row = row_idx.index;
         const int node = row_idx.node;
         const bool linking_row = row_idx.linking;
         const SystemType system_type = row_idx.system_type;
         assert(!linking_row); // todo

         const double lhs = float_values.at(first_float_val);
         const double rhs = float_values.at(first_float_val + 1);

         const int index_stored_row = int_values.at(first_int_val);
         const int iclow = int_values.at(first_int_val + 1);
         const int icupp = int_values.at(first_int_val + 2);
         assert(iclow + icupp >= 1);

         double value_row = stored_rows->localRowTimesVec(x_vec, node, index_stored_row, linking_row);

         if( system_type == EQUALITY_SYSTEM )
         {
            assert(-1 <= node && node < static_cast<int>(y_vec.children.size()));
            assert(wasRowRemoved(system_type, node, row, linking_row));
            assert(PIPSisEQ(lhs, rhs));
            assert(PIPSisEQ(value_row, rhs));

            /* set dual multiplier to zero and mark row as added */
            getSimpleVecFromRowStochVec(*padding_origrow_equality, node, linking_row)[row] = 1;
            getSimpleVecFromRowStochVec(y_vec, node, linking_row)[row] = 0;
         }
         else
         {
            assert(-1 <= node && node < static_cast<int>(z_vec.children.size()));
            assert(wasRowRemoved(system_type, node, row, linking_row));

            /* set dual multipliers to zero and mark row as added */
            getSimpleVecFromRowStochVec(*padding_origrow_inequality, node, linking_row)[row] = 1;

            getSimpleVecFromRowStochVec(z_vec, node, linking_row)[row] = 0;
            getSimpleVecFromRowStochVec(lambda_vec, node, linking_row)[row] = 0;
            getSimpleVecFromRowStochVec(pi_vec, node, linking_row)[row] = 0;

            getSimpleVecFromRowStochVec(s_vec, node, linking_row)[row] = value_row;

            assert(PIPSisLE(lhs, value_row));
            assert(PIPSisLE(value_row, rhs));
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
         throw std::runtime_error("BOUNDS_TIGHTENED not yet implemented");
         break;
      }
      case FIXED_COLUMN:
      {
         assert(first_index + 1 == next_first_index);
         assert(first_float_val + 1 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

         const INDEX& idx_col = indices.at(first_index);
         assert(idx_col.index_type == COL);

         const int column = idx_col.index;
         const int node = idx_col.node;
         const double value = float_values[first_float_val];

         assert(-1 <= node && node < static_cast<int>(x_vec.children.size()));
         assert(wasColumnRemoved(node, column));

         /* mark entry as set and set x value to fixation */
         getSimpleVecFromColStochVec(*padding_origcol, node)[column] = 1;
         getSimpleVecFromColStochVec(x_vec, node)[column] = value;

         // todo : dual postsolve
         break;
      }
      case FIXED_EMPTY_COLUMN:
      {
         assert(first_index + 1 == next_first_index);
         assert(first_float_val + 4 == next_first_float_val);
         assert(first_int_val + 2 == next_first_int_val);

         const INDEX& idx_col = indices.at(first_index);
         assert(idx_col.index_type == COL);

         const int column = idx_col.index;
         const int node = idx_col.node;
         const double value = float_values.at(first_float_val);
         const double obj_value = float_values.at(first_float_val + 1);
         const double lbx = float_values.at(first_float_val + 2);
         const double ubx = float_values.at(first_float_val + 3);
         const int ixlow = int_values.at(first_int_val);
         const int ixupp = int_values.at(first_int_val + 1);
         assert(-1 <= node && node < static_cast<int>(x_vec.children.size()));
         assert(wasColumnRemoved(node, column));

         /* primal */
         /* mark entry as set and set x value to fixation */
         getSimpleVecFromColStochVec(*padding_origcol, node)[column] = 1;
         getSimpleVecFromColStochVec(x_vec, node)[column] = value;

         if( ixlow == 1 )
            assert(PIPSisLT(lbx, value));
         if( ixupp == 1 )
            assert(PIPSisLT(value, ubx));

         /* dual */
         getSimpleVecFromColStochVec(gamma_vec, node)[column] = 0.0;
         getSimpleVecFromColStochVec(phi_vec, node)[column] = 0.0;

         if(!PIPSisZero(obj_value))
         {
            if( obj_value < 0 )
               getSimpleVecFromColStochVec(gamma_vec, node)[column] = obj_value;
            else
               getSimpleVecFromColStochVec(phi_vec, node)[column] = obj_value;
         }

         if(ixlow == 1)
            getSimpleVecFromColStochVec(v_vec, node)[column] = value - lbx;
         else
            getSimpleVecFromColStochVec(v_vec, node)[column] = 0.0;

         if( ixupp == 1)
            getSimpleVecFromColStochVec(w_vec, node)[column] = ubx - value;
         else
            getSimpleVecFromColStochVec(w_vec, node)[column] = 0.0;

         assert(PIPSisZero(getSimpleVecFromColStochVec(v_vec, node)[column] * getSimpleVecFromColStochVec(gamma_vec, node)[column]));
         assert(PIPSisZero(getSimpleVecFromColStochVec(w_vec, node)[column] * getSimpleVecFromColStochVec(phi_vec, node)[column]));
         break;
      }
      case SUBSTITUTED_COLUMN:
      {
         throw std::runtime_error("SUBSTITUTED_COLUMN not yet implemented");
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
         throw std::runtime_error("SINGLETON_EQUALITY_ROW not yet implemented");
         break;
      }
      case SINGLETON_INEQUALITY_ROW:
      {
         throw std::runtime_error("SINGLETON_INEQUALITY_ROW not yet implemented");
         break;
      }
      case FREE_COLUMN_SINGLETON:
      {
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 1 == next_first_float_val);
         assert(first_int_val + 1 == next_first_int_val);

         const INDEX& idx_col = indices.at(first_index);
         const INDEX& idx_row = indices.at(first_index + 1);
         assert(idx_col.index_type == COL);
         assert(idx_row.index_type == ROW);

         const int column = idx_col.index;
         const int node_column = idx_col.node;

         const bool linking_row = idx_row.linking;
//         const int row_idx = idx_row.index;
         const int node_row = idx_row.node;
//         const SystemType system_type = idx_row.system_type;

         const double rhs = float_values.at(first_float_val);
         const int stored_row_idx = int_values.at(first_int_val);

         assert(!linking_row); // todo

         assert(wasColumnRemoved(node_column, column));
         assert(PIPSisZero(getSimpleVecFromColStochVec(x_vec, node_column)[column]));

         /* mark column as set */
         getSimpleVecFromColStochVec(*padding_origcol, node_column)[column] = 1;
         getSimpleVecFromColStochVec(x_vec, node_column)[column] = 0;

         double value_row = stored_rows->localRowTimesVec(x_vec, node_row, stored_row_idx, linking_row);
         assert(std::abs(value_row) != std::numeric_limits<double>::infinity());

         getSimpleVecFromColStochVec(x_vec, node_column)[column] = rhs - value_row;

         break;
      }
      case PARALLEL_ROW_SUBSTITUTION:
      {
         assert(first_index + 2 == next_first_index);
         assert(first_float_val + 2 == next_first_float_val);
         assert(first_int_val == next_first_int_val);

         const INDEX& var_subst = indices.at(first_index);
         const INDEX& var_used_for_subst = indices.at(first_index + 1);
         assert(var_subst.index_type == COL);
         assert(var_used_for_subst.index_type == COL);

         const int column_subst = var_subst.index;
         const int node_subst = var_subst.node;

         const int column_used_for_subst = var_used_for_subst.index;
         const int node_used_for_subst = var_used_for_subst.node;

         const double scalar = float_values.at(first_float_val);
         const double translation = float_values.at(first_float_val + 1);

         assert(wasColumnRemoved(node_subst, column_subst));

         assert(!wasColumnRemoved(node_used_for_subst, column_used_for_subst));

         const double val_for_subst = getSimpleVecFromColStochVec(x_vec, node_used_for_subst)[column_used_for_subst];

         getSimpleVecFromColStochVec(*padding_origcol, node_subst)[column_subst] = 1;
         getSimpleVecFromColStochVec(x_vec, node_subst)[column_subst] = scalar * val_for_subst + translation;

         break;
      }
      case LINKING_VARS_SYNC_EVENT:
      {
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
      default:
      {
         throw std::runtime_error("Tried to postsolve unknown reduction type"); // todo add what was passed
         break;
      }

         // todo : an is-set Stochvector? - rather use padding StochVectors
      }
   }

   /* compute all s, t and u that have not yet been computed */

#ifdef NDEBUG
   /* assert that all primal variables have been set */
   double absmin_x;
   x_vec.absmin(absmin_x);
   assert( absmin_x < std::numeric_limits<double>::max() ); 
   assert( x_vec->isRootNodeInSync() );
#endif

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
