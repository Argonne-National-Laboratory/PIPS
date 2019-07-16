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
      padding_origcol(dynamic_cast<StochVector*>(original_problem.g->clone()))
{
   padding_origcol->setToConstant(1);
//   mapping_to_origrow_equality = dynamic_cast<StochVector*>(original_problem.bu->clone());
//   mapping_to_origrow_inequality = dynamic_cast<StochVector*>(original_problem.bu->clone());
//
//   mapping_to_origcol->fillWithIndices();
//   mapping_to_origrow_equality->fillWithIndices();
//   mapping_to_origrow_inequality->fillWithIndices();

   // count !? rows cols...
   start_idx_values.push_back(0);
}

StochPostsolver::~StochPostsolver()
{
   delete padding_origcol;
}

/** postsolve has to compute the optimal dual multipliers here and set the primal value accordingly */
//  todo : is it fine to compute the multiplier for the whole column right away (reduced costs)? // yes, if not yet reintroduced rows have multiplier 0 it does not matter
// (assuming the coefficients don't change - if they do we need to save the state of the row when it got removed)
void StochPostsolver::notifyFixedColumn( int node, unsigned int col, double value)
{
   assert( getSimpleVecColFromStochVec(*padding_origcol, node)[col] == 1 );
   getSimpleVecColFromStochVec(*padding_origcol, node)[col] = -1;

   reductions.push_back( FIXED_COLUMN );
   indices.push_back( INDEX(node, col) );
   values.push_back( value );

   finishNotify();
}

/** postsolve for this is simply to set all dual variables to zero - the row itself has no primal impact */
void StochPostsolver::notifyRedundantRow( SystemType system_type, int node, unsigned int row, bool linking_constraint )
{
   return;
   reductions.push_back( REDUNDANT_ROW );
   indices.push_back( INDEX(node, row) );
   values.push_back( ( (system_type == EQUALITY_SYSTEM ) ? 1 : -1 ) );

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
   values.insert( values.end(), values_row, values_row + length );
   values.insert( values.end(), indices_row, indices_row + length );

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
//   assert( values.size() == indices.size() ); // todo why ? probably padding

   start_idx_values.push_back( values.size() );
}

// todo : sort reductions by nodes ? and then reverse ?
// todo : at the moment only replaces whatever is given as x with the x solution in the original soution space
PostsolveStatus StochPostsolver::postsolve(const Variables& reduced_solution, Variables& original_solution) const
{
   // todo if my_rank == 0
   // if()
   // std::cout << "postsolve" << std::endl;

   const sVars& stoch_reduced_sol = dynamic_cast<const sVars&>(reduced_solution);
   sVars& stoch_original_sol = dynamic_cast<sVars&>(original_solution);

   /* primal variables */
   const StochVector& primal_vars_reduced = dynamic_cast<const StochVector&>(*stoch_reduced_sol.x);
   StochVector& primal_vars_orig = dynamic_cast<StochVector&>(*stoch_original_sol.x);

   setOriginalValuesFromReduced(primal_vars_orig, primal_vars_reduced, *padding_origcol);
   /* original solution is now reduced solution padded with zeros */


   /* dual variables */

   // todo

   /* dual solution is now reduced solution padded with zeros */

   /* post-solve the reductions in reverse order */
   for( int i = reductions.size() - 1; i >= 0; --i )
   {
      int type = reductions[i];
      unsigned int first = start_idx_values[i];
      // unsigned int last = start_idx_values[i + 1];

      switch( type )
      {
         case REDUNDANT_ROW:
         {
            break;
         }
         case BOUNDS_TIGHTENED:
         {
            break;
         }
         case FIXED_COLUMN:
         {
            int column = indices[i].index;
            int node = indices[i].node;

            assert( -1 <= node && node < static_cast<int>(primal_vars_orig.children.size()) );
            assert( getSimpleVecColFromStochVec(*padding_origcol, node)[column] == -1);

            assert( fabs(values[first]) != std::numeric_limits<double>::infinity());
            getSimpleVecColFromStochVec(primal_vars_orig, node)[column] = values[first];
            break;
         }
         case SUBSTITUTED_COLUMN:
         {
            // todo
            break;
         }
         case PARALLEL_COLUMN:
         {
            // todo
            break;
         }
         case DELETED_ROW:
         {
            // todo
            break;
         }
         default:
         {
            throw std::runtime_error("Tried to postsolve unknown reduction type"); // todo add what was passed
            break;
         }

         // todo ? an isset array
         // todo  check solution for feasibility - kkt checker
      }
   }

   return PRESOLVE_OK;
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
