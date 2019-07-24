/*
 * StochPostsolver.C
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */


#include "StochPostsolver.h"
#include "pipsdef.h"

#include <limits>
#include <stdexcept>
#include <iostream>

StochPostsolver::StochPostsolver(const sData& original_problem) :
   QpPostsolver(original_problem), n_rows_original(original_problem.my + original_problem.mz), n_cols_original(original_problem.nx),
      padding_origcol(dynamic_cast<StochVector*>(original_problem.g->clone())), padding_origrow_equality(dynamic_cast<StochVector*>(original_problem.bA->clone())),
      padding_origrow_inequality(dynamic_cast<StochVector*>(original_problem.bu->clone()))
{
   getRankDistributed(MPI_COMM_WORLD, my_rank, distributed);

   padding_origcol->setToConstant(1);
   padding_origrow_equality->setToConstant(1);
   padding_origrow_inequality->setToConstant(1);

   // count !? rows cols...
   start_idx_values.push_back(0);
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
   if( node == -1 )
      assert(block_type == LINKING_VARS_BLOCK);
   assert( block_type != LINKING_CONS_BLOCK );

   assert( getSimpleVecColFromStochVec(*padding_origcol, node)[col] == 1 );
   getSimpleVecColFromStochVec(*padding_origcol, node)[col] = -1;

   reductions.push_back( SINGLETON_EQUALITY_ROW );
   indices.push_back( INDEX(node, row) );
   values.push_back( coeff );
   values.push_back( rhs );
   val_idx.push_back( col );

   finishNotify();
}





void StochPostsolver::notifySingletonIneqalityRow( int node, int row, BlockType block_type, int col, double coeff, double lhs, double rhs )
{
   if( node == -1 )
      assert(block_type == LINKING_VARS_BLOCK);
   assert( block_type != LINKING_CONS_BLOCK );

   reductions.push_back( SINGLETON_INEQUALITY_ROW);
   indices.push_back( INDEX(node, row) );
   values.push_back( coeff );
   values.push_back( lhs );
   values.push_back( rhs );
   val_idx.push_back( col );

   finishNotify();
}


/** postsolve has to compute the optimal dual multipliers here and set the primal value accordingly 
 * The column is passed in the following format:
 * 
 * node != -1
 * index: idx Bmat, inf, idx Blmat 
 * value: val Bmat, inf, val Blmat
 *
 * node == -1
 * index: idx A0mat, inf, i, Aimat, inf, j, Ajmat, inf,...
 * value: val A0mat, inf, val Aimat, inf, val Ajmat, int,...
 *
 * Only rank 0 gets A0mat and C0mat
 * When a column is fixed one automatically has val = ubx = lbx
 */
void StochPostsolver::notifyFixedColumn( int node, unsigned int col, double value, const std::vector<int>& indices_col,
   const std::vector<double>& values_col)
{
   assert( getSimpleVecColFromStochVec(*padding_origcol, node)[col] == 1 );
   getSimpleVecColFromStochVec(*padding_origcol, node)[col] = -1;

   // todo assert correct format
   /* store current upper and lower bounds of x and the local column */
   reductions.push_back( FIXED_COLUMN );
   indices.push_back( INDEX(node, col) );
   values.push_back( value );
   values.insert( values.end(), values_col.begin(), values_col.end() );
   val_idx.insert( val_idx.end(), indices_col.begin(), indices_col.end() );

   finishNotify();
}

/** postsolve for this is simply to set all dual variables to zero - the row itself has no primal impact */
void StochPostsolver::notifyRedundantRow( SystemType system_type, int node, unsigned int row, bool linking_constraint, const std::vector<int>& indices_row,
   const std::vector<double> values_row )
{
   assert( getSimpleVecRowFromStochVec( (system_type == EQUALITY_SYSTEM) ? *padding_origrow_equality : *padding_origrow_inequality, node, 
      (linking_constraint) ? LINKING_CONS_BLOCK : CHILD_BLOCK)[row] == 1 );
   getSimpleVecRowFromStochVec( (system_type == EQUALITY_SYSTEM) ? *padding_origrow_equality : *padding_origrow_inequality, node,
      (linking_constraint) ? LINKING_CONS_BLOCK : CHILD_BLOCK)[row] = -1;


   reductions.push_back( REDUNDANT_ROW );
   indices.push_back( INDEX(node, row) );
   values.push_back( ( (system_type == EQUALITY_SYSTEM ) ? 0 : 1 ) );
   values.push_back( linking_constraint );

   values.insert( values.end(), values_row.begin(), values_row.end() );
   val_idx.insert( val_idx.end(), indices_row.begin(), indices_row.end() );

   finishNotify();
}

// todo : only store each version of each row once!
// todo : store whole row
void StochPostsolver::notifyRowPropagated( SystemType system_type, int node, int row, bool linking_constraint,
      int column, double lb, double ub, double* values_row, int* indices_row, int length)
{
   return;
   /* store the row with which bound has been tightened */
   reductions.push_back( BOUNDS_TIGHTENED );
   indices.push_back( INDEX(node, row) );

   /* values contains : {system_type (-1 or 1), tightened_column, new ub, new lb, length_row, values_row, col_indices_row */
   values.push_back( (system_type == EQUALITY_SYSTEM ) ? 1 : -1 );
   values.push_back( column );
   values.push_back( ub );
   values.push_back( lb );
   values.push_back( length );
   //values.insert( values.end(), values_row, values_row + length );
   //values.insert( values.end(), indices_row, indices_row + length );

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
   assert( reductions.size() == start_idx_values.size() );
   assert( indices.size() == start_idx_values.size() );
   // assert( values.size() == indices.size() ); // todo why ? probably padding

   start_idx_values.push_back( values.size() );
   start_idx_indices.push_back( val_idx.size() );
}

// todo : sort reductions by nodes ? and then reverse ?
// todo : at the moment only replaces whatever is given as x with the x solution in the original soution space
PostsolveStatus StochPostsolver::postsolve(const Variables& reduced_solution, Variables& original_solution) const
{
   if(my_rank == 0)
      std::cout << "postsolving... " << std::endl;

   const sVars& stoch_reduced_sol = dynamic_cast<const sVars&>(reduced_solution);
   sVars& stoch_original_sol = dynamic_cast<sVars&>(original_solution);

   setOriginalVarsFromReduced(stoch_reduced_sol, stoch_original_sol);

   StochVector& x_vec = dynamic_cast<StochVector&>(*stoch_original_sol.x);
   StochVector& v_vec = dynamic_cast<StochVector&>(*stoch_original_sol.v);
   StochVector& w_vec = dynamic_cast<StochVector&>(*stoch_original_sol.w);

   StochVector& y_vec = dynamic_cast<StochVector&>(*stoch_original_sol.y);
   StochVector& z_vec = dynamic_cast<StochVector&>(*stoch_original_sol.z);
   StochVector& lambda_vec = dynamic_cast<StochVector&>(*stoch_original_sol.lambda);
   StochVector& pi_vec = dynamic_cast<StochVector&>(*stoch_original_sol.pi);

   /* original variables are now reduced vars padded with zeros */


   /* dual variables */

   // todo

   /* dual solution is now reduced solution padded with zeros */

   /* post-solve the reductions in reverse order */
   for( int i = reductions.size() - 1; i >= 0; --i )
   {
      int type = reductions[i];
      unsigned int first_idx = start_idx_values[i];
      unsigned int last_idx = start_idx_values[i + 1];
      unsigned int first_val = start_idx_values[i];
      unsigned int last_val = start_idx_values[i + 1];
      switch( type )
      {
         case REDUNDANT_ROW:
         {
            assert(first_val + 2 == last_val);
            break;
            // todo
            SystemType system_type = (values[first_val] == 0) ? EQUALITY_SYSTEM : INEQUALITY_SYSTEM;
            bool linking = values[first_val + 1];
            int node = indices[i].node;
            int row = indices[i].index;
            BlockType block_type = linking ? LINKING_CONS_BLOCK : CHILD_BLOCK;

            /* all duals set to zero - no primal poststolve necessary */
            if(system_type == EQUALITY_SYSTEM)
            {
               assert( getSimpleVecRowFromStochVec(*padding_origrow_equality, node, block_type)[row] == -1 );
               getSimpleVecRowFromStochVec(y_vec, node, block_type)[row] = 0;
            }
            else
            {
               assert( getSimpleVecRowFromStochVec(*padding_origrow_inequality, node, block_type)[row] == -1 );

               getSimpleVecRowFromStochVec(z_vec, node, block_type)[row] = 0;
               getSimpleVecRowFromStochVec(lambda_vec, node, block_type)[row] = 0;
               getSimpleVecRowFromStochVec(pi_vec, node, block_type)[row] = 0;

               /* we do not compute s t and u here since it is not necessary .. it is though */
            }
            break;
         }
         case BOUNDS_TIGHTENED:
         {
            /* the dual multiplier of the row that was responsible for the bound change gets adjusted */

            
            break;
            throw std::runtime_error("BOUNDS_TIGHTENED not yet implemented");
         }
         case FIXED_COLUMN:
         {
            int column = indices[i].index;
            int node = indices[i].node;
            double value = values[first_val];

            assert( -1 <= node && node < static_cast<int>(x_vec.children.size()) );
            assert( getSimpleVecColFromStochVec(*padding_origcol, node)[column] == -1 );
            getSimpleVecColFromStochVec(x_vec, node)[column] = value;
            getSimpleVecColFromStochVec(v_vec, node)[column] = value;
            getSimpleVecColFromStochVec(w_vec, node)[column] = value;



            break;

   //             reductions.push_back( FIXED_COLUMN );
   // indices.push_back( INDEX(node, col) );
   // values.push_back( value );
   // values.insert( values.end(), values_col.begin(), values_col.end() );
   // val_idx.insert( val_idx.end(), indices_col.begin(), indices_col.end() );
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
            int node = indices[i].node;
            int row = indices[i].index;
            BlockType block_type = (node == -1) ? LINKING_VARS_BLOCK : CHILD_BLOCK;

            break;
         }
         case SINGLETON_INEQUALITY_ROW:
         {
            int node = indices[i].node;
            int row = indices[i].index;
            BlockType block_type = (node == -1) ? LINKING_VARS_BLOCK : CHILD_BLOCK;


            break;
         }
         default:
         {
            throw std::runtime_error("Tried to postsolve unknown reduction type"); // todo add what was passed
            break;
         }

         // todo ? an is-set array
         // todo  check solution for feasibility - kkt checker ?
      }
   }

   /* compute all s, t and u that have not yet been computed */

   return PRESOLVE_OK;
}

void StochPostsolver::setOriginalVarsFromReduced(const sVars& reduced_vars, sVars& original_vars) const
{
   /* x */
   const StochVector& x_reduced = dynamic_cast<const StochVector&>(*reduced_vars.x);
   StochVector& x_orig = dynamic_cast<StochVector&>(*original_vars.x);
   setOriginalValuesFromReduced(x_orig, x_reduced, *padding_origcol);

   /* s */
   const StochVector& s_reduced = dynamic_cast<const StochVector&>(*reduced_vars.s);
   StochVector& s_orig = dynamic_cast<StochVector&>(*original_vars.s);
   setOriginalValuesFromReduced(s_orig, s_reduced, *padding_origrow_inequality);

   /* y */
   const StochVector& y_reduced = dynamic_cast<const StochVector&>(*reduced_vars.y);
   StochVector& y_orig = dynamic_cast<StochVector&>(*original_vars.y);
   setOriginalValuesFromReduced(y_orig, y_reduced, *padding_origrow_equality);

   /* z */
   const StochVector& z_reduced = dynamic_cast<const StochVector&>(*reduced_vars.z);
   StochVector& z_orig = dynamic_cast<StochVector&>(*original_vars.z);
   setOriginalValuesFromReduced(z_orig, z_reduced, *padding_origrow_inequality);
   
   /* v */
   const StochVector& v_reduced = dynamic_cast<const StochVector&>(*reduced_vars.v);
   StochVector& v_orig = dynamic_cast<StochVector&>(*original_vars.v);
   setOriginalValuesFromReduced(v_orig, v_reduced, *padding_origcol);

   /* gamma */
   const StochVector& gamma_reduced = dynamic_cast<const StochVector&>(*reduced_vars.gamma);
   StochVector& gamma_orig = dynamic_cast<StochVector&>(*original_vars.gamma);
   setOriginalValuesFromReduced(gamma_orig, gamma_reduced, *padding_origcol);

   /* w */
   const StochVector& w_reduced = dynamic_cast<const StochVector&>(*reduced_vars.w);
   StochVector& w_orig = dynamic_cast<StochVector&>(*original_vars.w);
   setOriginalValuesFromReduced(w_orig, w_reduced, *padding_origcol);

   /* phi */
   const StochVector& phi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.phi);
   StochVector& phi_orig = dynamic_cast<StochVector&>(*original_vars.phi);
   setOriginalValuesFromReduced(phi_orig, phi_reduced, *padding_origcol);

   /* t */
   const StochVector& t_reduced = dynamic_cast<const StochVector&>(*reduced_vars.t);
   StochVector& t_orig = dynamic_cast<StochVector&>(*original_vars.t);
   setOriginalValuesFromReduced(t_orig, t_reduced, *padding_origrow_inequality);

   /* lambda */
   const StochVector& lambda_reduced = dynamic_cast<const StochVector&>(*reduced_vars.lambda);
   StochVector& lambda_orig = dynamic_cast<StochVector&>(*original_vars.lambda);
   setOriginalValuesFromReduced(lambda_orig, lambda_reduced, *padding_origrow_inequality);

   /* u */
   const StochVector& u_reduced = dynamic_cast<const StochVector&>(*reduced_vars.u);
   StochVector& u_orig = dynamic_cast<StochVector&>(*original_vars.u);
   setOriginalValuesFromReduced(u_orig, u_reduced, *padding_origrow_inequality);

   /* pi */
   const StochVector& pi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.pi);
   StochVector& pi_orig = dynamic_cast<StochVector&>(*original_vars.pi);
   setOriginalValuesFromReduced(pi_orig, pi_reduced, *padding_origrow_inequality);
}


/// fills vars_orig with vars_reduced padded with zeros - padding is done via the padding_map
void StochPostsolver::setOriginalValuesFromReduced(StochVector& original_vector, const StochVector& reduced_vector, const StochVector& padding_original) const
{
   assert( reduced_vector.children.size() == original_vector.children.size() );
   assert( padding_original.children.size() == reduced_vector.children.size() );
   assert( reduced_vector.vec != NULL && original_vector.vec != NULL && padding_original.vec != NULL );
   assert( (reduced_vector.vecl != NULL && original_vector.vecl != NULL && padding_original.vecl != NULL)
         || (reduced_vector.vecl == NULL && original_vector.vecl == NULL && padding_original.vecl == NULL) );

   if( reduced_vector.isKindOf(kStochDummy) )
   {
      assert( original_vector.isKindOf(kStochDummy) && padding_original.isKindOf(kStochDummy) );
      return;
   }

   /* root node */
   /* vec */
   setOriginalValuesFromReduced( dynamic_cast<SimpleVector&>(*original_vector.vec), dynamic_cast<const SimpleVector&>(*reduced_vector.vec),
      dynamic_cast<const SimpleVector&>(*padding_original.vec));

   /* vecl */
   if( reduced_vector.vecl )
   {
      setOriginalValuesFromReduced( dynamic_cast<SimpleVector&>(*original_vector.vecl), dynamic_cast<const SimpleVector&>(*reduced_vector.vecl),
         dynamic_cast<const SimpleVector&>(*padding_original.vecl));
   }

   /* child nodes */
   for( unsigned int i = 0; i < reduced_vector.children.size(); ++i )
   {
      setOriginalValuesFromReduced( *original_vector.children[i],  *reduced_vector.children[i], *padding_original.children[i] );
   }
}

void StochPostsolver::setOriginalValuesFromReduced(SimpleVector& original_vector, const SimpleVector& reduced_vector, const SimpleVector& padding_original) const
{
   assert( original_vector.length() == padding_original.length() );

   unsigned int col_reduced = 0; 
   for(unsigned int i = 0; i < padding_original.length(); ++i)
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
   if( col_reduced != reduced_vector.length())
      std::cout << col_reduced << "\t" << reduced_vector.length() << std::endl;
   assert(col_reduced == reduced_vector.length());
}

/// todo : codu duplication with presolveData.h
SimpleVector& StochPostsolver::getSimpleVecRowFromStochVec(const StochVector& stochvec, int node, BlockType block_type) const
{
   assert(-1 <= node && node < static_cast<int>(stochvec.children.size()));

   if(node == -1)
   {
      if(block_type == LINKING_CONS_BLOCK)
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVector&>(*(stochvec.vecl));
      }
      else
      {
         assert(stochvec.vec);
         return dynamic_cast<SimpleVector&>(*(stochvec.vec));
      }
   }
   else
   {
      if(block_type == CHILD_BLOCK || block_type == LINKING_VARS_BLOCK)
      {
         assert(stochvec.children[node]->vec);
         return dynamic_cast<SimpleVector&>(*(stochvec.children[node]->vec));
      }
      else
      {
         assert(stochvec.vecl);
         return dynamic_cast<SimpleVector&>(*(stochvec.vecl));
      }
   }
}

SimpleVector& StochPostsolver::getSimpleVecColFromStochVec(const StochVector& stochvec, int node) const
{
   assert(-1 <= node && node < static_cast<int>(stochvec.children.size()));

   if(node == -1)
   {
      assert(stochvec.vecl == NULL);
      assert(stochvec.vec);
      return dynamic_cast<SimpleVector&>(*(stochvec.vec));
   }
   else
   {
      assert(stochvec.children[node]->vecl == NULL);
      assert(stochvec.children[node]->vec);
      return dynamic_cast<SimpleVector&>(*(stochvec.children[node]->vec));
   }
}
