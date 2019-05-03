/*
 * StochPostsolver.C
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */


#include "StochPostsolver.h"
#include "pipsdef.h"

#include <stdexcept>

StochPostsolver::StochPostsolver(const sData& original_problem, int n_rows_original, int n_cols_original) :
   QpPostsolver(original_problem), n_rows_original(n_rows_original), n_cols_original(n_cols_original)
{
   // todo set up stochvectors for mapping
   mapping_to_origcol = dynamic_cast<StochVector*>(original_problem.g->clone());
   mapping_to_origrow_equality = dynamic_cast<StochVector*>(original_problem.bu->clone());
   mapping_to_origrow_inequality = dynamic_cast<StochVector*>(original_problem.bu->clone());

   // count !? rows cols...

   start_indices.push_back(0);
}

StochPostsolver::~StochPostsolver(){}


void StochPostsolver::notifyFixedColumn( int node, unsigned int col, double value)
{
   reductions.push_back(FIXED_COLUMN);
//   indices.push_back( mapping_to_origcol[col]);
   nodes.push_back( node );
   values.push_back( value );

   finishNotify();
}

void StochPostsolver::notifyDeletedRow( SystemType system_type, int node, unsigned int row, bool linking_constraint)
{
   throw std::runtime_error("Not yet implemented");
}

void StochPostsolver::notifyParallelColumns()
{
   throw std::runtime_error("Not yet implemented");
}

void StochPostsolver::finishNotify()
{
   assert( reductions.size() == start_indices.size() );
   assert( nodes.size() == start_indices.size() );
   assert( values.size() == indices.size() ); // todo why ? probably padding

   start_indices.push_back( values.size() );
}

// todo : sort reductions by nodes ? and then reverse ?
// todo : at the moment only replaces whatever is given as x with the x solution in the original soution space
PostsolveStatus StochPostsolver::postsolve(const sVars& reduced_solution, sVars& original_solution) const
{

   /* primal variables */

   const StochVector& primal_vars_reduced = dynamic_cast<const StochVector&>(*reduced_solution.x);
   StochVector& primal_vars_orig = dynamic_cast<StochVector&>(*original_solution.x);

   primal_vars_orig.setToZero(); // todo necessary?
   setReducedValuesInOrigVector( primal_vars_reduced, primal_vars_orig, *mapping_to_origcol);

   /* original solution is now reduced solution padded with zeros */

   /* dual variables */

   // todo

   /* dual solution is now reduced solution padded with zeros */

   /* post-solve the reductions in reverse order */
   for( int i = reductions.size() - 1; i >= 0; --i )
   {
      int type = reductions[i];
      unsigned int first = start_indices[i];
      unsigned int last = start_indices[i + 1];

      switch( type )
      {
         case FIXED_COLUMN:
         {
            int column = indices[first];
            int node = nodes[first];

            assert( -1 <= node && node < static_cast<int>(primal_vars_orig.children.size()) );

            SimpleVector& orig_sol = (node == -1) ? dynamic_cast<SimpleVector&>(*primal_vars_orig.vec)
                  : dynamic_cast<SimpleVector&>(*primal_vars_orig.children[node]->vec);

            orig_sol[column] = values[first];
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

void StochPostsolver::setReducedValuesInOrigVector(const StochVector& reduced_vector, StochVector& original_vector, const StochVector& mapping_to_original) const
{
   assert( reduced_vector.children.size() == original_vector.children.size() );
   assert( mapping_to_original.children.size() == reduced_vector.children.size() );
   assert( reduced_vector.vec != NULL && original_vector.vec != NULL && mapping_to_original.vec != NULL );
   assert( (reduced_vector.vecl != NULL && original_vector.vecl != NULL && mapping_to_original.vecl != NULL)
         || (reduced_vector.vecl == NULL && original_vector.vecl == NULL && mapping_to_original.vecl == NULL) );

   if( reduced_vector.isKindOf(kStochDummy) )
   {
      assert( original_vector.isKindOf(kStochDummy) && mapping_to_original.isKindOf(kStochDummy) );
      return;
   }

   /* root node */
   /* vec */
   setReducedValuesInOrigVector( dynamic_cast<const SimpleVector&>(*reduced_vector.vec),
         dynamic_cast<SimpleVector&>(*original_vector.vec), dynamic_cast<const SimpleVector&>(*mapping_to_original.vec));

   /* vecl */
   if( reduced_vector.vecl )
   {
      setReducedValuesInOrigVector( dynamic_cast<const SimpleVector&>(*reduced_vector.vecl),
            dynamic_cast<SimpleVector&>(*original_vector.vecl), dynamic_cast<const SimpleVector&>(*mapping_to_original.vecl));
   }

   /* child nodes */
   for( unsigned int i = 0; i < reduced_vector.children.size(); ++i )
   {
      setReducedValuesInOrigVector( *reduced_vector.children[i], *original_vector.children[i], *mapping_to_original.children[i] );
   }
}

void StochPostsolver::setReducedValuesInOrigVector(const SimpleVector& reduced_vector, SimpleVector& original_vector, const SimpleVector& mapping_to_original) const
{
   assert( reduced_vector.length() == mapping_to_original.length() );

   for(unsigned int i = 0; i < reduced_vector.length(); ++i)
   {
      unsigned int orig_col = static_cast<unsigned int>(mapping_to_original[i]);

      assert( PIPSisEQ(orig_col, mapping_to_original[i]) );
      assert( orig_col < original_vector.length() );

      original_vector[orig_col] = reduced_vector[i];
   }
}

