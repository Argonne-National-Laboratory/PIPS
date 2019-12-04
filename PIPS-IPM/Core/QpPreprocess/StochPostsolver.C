/*
 * StochPostsolver.C
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */


#include "StochPostsolver.h"
#include "StochVectorUtilities.h"
#include "pipsdef.h"
#include "pipsport.h"

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
   padding_origrow_inequality( cloneStochVector<double, int>(*original_problem.bu) )
{
   padding_origcol->setToConstant(1);
   padding_origrow_equality->setToConstant(1);
   padding_origrow_inequality->setToConstant(1);

   // count !? rows cols...
   start_idx_values.push_back(0);
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

/* substitute var2 with scalar*var1 + translation */
void StochPostsolver::notifyParallelRowSubstitution(SystemType system_type, int node_row, int var1, int row1, int node_var1, int var2, int row2, 
   int node_var2, double scalar, double translation)
{
   assert( getSimpleVecFromColStochVec(*padding_origcol, node_var2)[var2] == 1 );
   assert( getSimpleVecFromColStochVec(*padding_origcol, node_var1)[var1] == 1 );

   reductions.push_back(PARALLEL_ROW_SUBSTITUTION);

   indices.push_back( INDEX(node_var2, var2) );
   getSimpleVecFromColStochVec(*padding_origcol, node_var2)[var2] = -1;

   values.push_back( var1 );
   values.push_back( node_var1 );
   values.push_back( scalar );
   values.push_back( translation );

   finishNotify();
}

/** postsolve has to compute the optimal dual multipliers here and set the primal value accordingly */
void StochPostsolver::notifyFixedColumn( int node, unsigned int col, double value, const std::vector<int>& indices_col,
   const std::vector<double>& values_col)
{
   assert( getSimpleVecFromColStochVec(*padding_origcol, node)[col] == 1 );
   getSimpleVecFromColStochVec(*padding_origcol, node)[col] = -1;
   // todo assert correct format
   // todo add matrix of columns
   /* store current upper and lower bounds of x and the local column */
   reductions.push_back( FIXED_COLUMN );
   indices.push_back( INDEX(node, col) );
   values.push_back( value );

   finishNotify();
}


void StochPostsolver::notifyFixedEmptyColumn(int node, unsigned int col, double value)
{
   if( getSimpleVecFromColStochVec(*padding_origcol, node)[col] != 1 )
   {
      std::cout << "WARNING: attempt to fix col " << col << " on node " 
         << node << " a second time even though it has already been fixed" << std::endl;
      return;
   }

   assert( getSimpleVecFromColStochVec(*padding_origcol, node)[col] == 1);
   getSimpleVecFromColStochVec(*padding_origcol, node)[col] = -1;
   assert( std::fabs(value) < 1e10 );

   reductions.push_back( FIXED_EMPTY_COLUMN );   
   indices.push_back( INDEX(node, col) );
   values.push_back( value );

   finishNotify();
}

/** postsolve for this is simply to set all dual variables to zero - the row itself has no primal impact */
void StochPostsolver::notifyRedundantRow( SystemType system_type, int node, unsigned int row, bool linking_constraint, const std::vector<int>& indices_row,
   const std::vector<double> values_row )
{
   return;

   // finishNotify();
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

void StochPostsolver::notifyParallelColumns()
{
   throw std::runtime_error("Not yet implemented");
}

void StochPostsolver::finishNotify()
{
   assert( reductions.size() == start_idx_values.size() );
   assert( indices.size() == start_idx_values.size() );

   start_idx_values.push_back( values.size() );
}

// todo : sort reductions by nodes ? and then reverse ?
// todo : at the moment only replaces whatever is given as x with the x solution in the original soution space
PostsolveStatus StochPostsolver::postsolve(const Variables& reduced_solution, Variables& original_solution) const
{
   if(my_rank == 0)
      std::cout << "start postsolving... " << std::endl;

   const sVars& stoch_reduced_sol = dynamic_cast<const sVars&>(reduced_solution);
   sVars& stoch_original_sol = dynamic_cast<sVars&>(original_solution);

   stoch_original_sol.x->setToConstant(-std::numeric_limits<double>::infinity());

   /* original variables are now reduced vars padded with zeros */
   setOriginalVarsFromReduced(stoch_reduced_sol, stoch_original_sol);

   StochVector& x_vec = dynamic_cast<StochVector&>(*stoch_original_sol.x);



   /* dual variables */

   // todo

   /* dual solution is now reduced solution padded with zeros */

   /* post-solve the reductions in reverse order */
   for( int i = reductions.size() - 1; i >= 0; --i )
   {
      int type = reductions.at(i);
      unsigned int first_val = start_idx_values.at(i);
      unsigned int last_val = start_idx_values.at(i + 1);
      switch( type )
      {
         case REDUNDANT_ROW:
         {
            throw std::runtime_error("REDUNDANT_ROW not yet implemented");
            break;
         }
         case BOUNDS_TIGHTENED:
         {
            throw std::runtime_error("BOUNDS_TIGHTENED not yet implemented");
            break;
         }
         case FIXED_COLUMN:
         {
            int column = indices.at(i).index;
            int node = indices.at(i).node;
            double value = values[first_val];

            assert( -1 <= node && node < static_cast<int>(x_vec.children.size()) );
            assert( getSimpleVecFromColStochVec(*padding_origcol, node)[column] == -1 );

            getSimpleVecFromColStochVec(*padding_origcol, node)[column] = 1;
            getSimpleVecFromColStochVec(x_vec, node)[column] = value;

            break;
         }
         case FIXED_EMPTY_COLUMN:
         {
            assert( first_val == last_val -1 );

            int column = indices.at(i).index;
            int node = indices.at(i).node;
            double value = values.at(first_val);

            assert( -1 <= node && node < static_cast<int>(x_vec.children.size()) );
            assert( getSimpleVecFromColStochVec(*padding_origcol, node)[column] == -1 );

            getSimpleVecFromColStochVec(*padding_origcol, node)[column] = 1;
            getSimpleVecFromColStochVec(x_vec, node)[column] = value;

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
         case PARALLEL_ROW_SUBSTITUTION:
         {
            int column = indices.at(i).index;
            int node = indices.at(i).node;

            assert( PIPSisEQ( getSimpleVecFromColStochVec(*padding_origcol, node)[column], -1) );
            assert( first_val == last_val - 4 );

            int col_sub = values.at(first_val);
            int node_sub = values.at(first_val + 1);
            double scalar = values.at(first_val + 2);
            double translation = values.at(first_val + 3);

            assert( PIPSisEQ( getSimpleVecFromColStochVec(*padding_origcol, node_sub)[col_sub], 1) );
            const double val_sub = getSimpleVecFromColStochVec(x_vec, node_sub)[col_sub]; 

            getSimpleVecFromColStochVec(*padding_origcol, node)[column] = 1;
            getSimpleVecFromColStochVec(x_vec, node)[column] = scalar * val_sub + translation;


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


#ifdef NDEBUG
   /* assert that all primal variables have been set */
   double absmin_x;
   x_vec.absmin(absmin_x);
   assert( absmin_x < std::numeric_limits<double>::max() ); 
   assert( x_vec->isRootNodeInSync() );
#endif

   if(my_rank == 0)
      std::cout << "finished postsolving... " << std::endl;
   
   return PRESOLVE_OK;
}

void StochPostsolver::setOriginalVarsFromReduced(const sVars& reduced_vars, sVars& original_vars) const
{
   /* x */
   const StochVector& x_reduced = dynamic_cast<const StochVector&>(*reduced_vars.x);
   StochVector& x_orig = dynamic_cast<StochVector&>(*original_vars.x);
   setOriginalValuesFromReduced<>(x_orig, x_reduced, *padding_origcol);

   // /* s */
   // const StochVector& s_reduced = dynamic_cast<const StochVector&>(*reduced_vars.s);
   // StochVector& s_orig = dynamic_cast<StochVector&>(*original_vars.s);
   // setOriginalValuesFromReduced(s_orig, s_reduced, *padding_origrow_inequality);

   // /* y */
   // const StochVector& y_reduced = dynamic_cast<const StochVector&>(*reduced_vars.y);
   // StochVector& y_orig = dynamic_cast<StochVector&>(*original_vars.y);
   // setOriginalValuesFromReduced(y_orig, y_reduced, *padding_origrow_equality);

   // /* z */
   // const StochVector& z_reduced = dynamic_cast<const StochVector&>(*reduced_vars.z);
   // StochVector& z_orig = dynamic_cast<StochVector&>(*original_vars.z);
   // setOriginalValuesFromReduced(z_orig, z_reduced, *padding_origrow_inequality);
   
   // /* v */
   // const StochVector& v_reduced = dynamic_cast<const StochVector&>(*reduced_vars.v);
   // StochVector& v_orig = dynamic_cast<StochVector&>(*original_vars.v);
   // setOriginalValuesFromReduced(v_orig, v_reduced, *padding_origcol);

   // /* gamma */
   // const StochVector& gamma_reduced = dynamic_cast<const StochVector&>(*reduced_vars.gamma);
   // StochVector& gamma_orig = dynamic_cast<StochVector&>(*original_vars.gamma);
   // setOriginalValuesFromReduced(gamma_orig, gamma_reduced, *padding_origcol);

   // /* w */
   // const StochVector& w_reduced = dynamic_cast<const StochVector&>(*reduced_vars.w);
   // StochVector& w_orig = dynamic_cast<StochVector&>(*original_vars.w);
   // setOriginalValuesFromReduced(w_orig, w_reduced, *padding_origcol);

   // /* phi */
   // const StochVector& phi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.phi);
   // StochVector& phi_orig = dynamic_cast<StochVector&>(*original_vars.phi);
   // setOriginalValuesFromReduced(phi_orig, phi_reduced, *padding_origcol);

   // /* t */
   // const StochVector& t_reduced = dynamic_cast<const StochVector&>(*reduced_vars.t);
   // StochVector& t_orig = dynamic_cast<StochVector&>(*original_vars.t);
   // setOriginalValuesFromReduced(t_orig, t_reduced, *padding_origrow_inequality);

   // /* lambda */
   // const StochVector& lambda_reduced = dynamic_cast<const StochVector&>(*reduced_vars.lambda);
   // StochVector& lambda_orig = dynamic_cast<StochVector&>(*original_vars.lambda);
   // setOriginalValuesFromReduced(lambda_orig, lambda_reduced, *padding_origrow_inequality);

   // /* u */
   // const StochVector& u_reduced = dynamic_cast<const StochVector&>(*reduced_vars.u);
   // StochVector& u_orig = dynamic_cast<StochVector&>(*original_vars.u);
   // setOriginalValuesFromReduced(u_orig, u_reduced, *padding_origrow_inequality);

   // /* pi */
   // const StochVector& pi_reduced = dynamic_cast<const StochVector&>(*reduced_vars.pi);
   // StochVector& pi_orig = dynamic_cast<StochVector&>(*original_vars.pi);
   // setOriginalValuesFromReduced(pi_orig, pi_reduced, *padding_origrow_inequality);
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
