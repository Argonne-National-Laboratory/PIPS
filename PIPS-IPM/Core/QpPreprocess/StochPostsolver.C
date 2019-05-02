/*
 * StochPostsolver.C
 *
 *  Created on: 02.05.2019
 *      Author: Nils-Christian Kempke
 */


#include "StochPostsolver.h"


#include <stdexcept>

StochPostsolver::StochPostsolver(const sData& original_problem, int n_rows_original, int n_cols_original) :
      original_problem(original_problem), n_rows_original(n_rows_original), n_cols_original(n_cols_original)
{
   // todo set up stochvectors for mapping
   mapping_to_origcol = dynamic_cast<StochVector*>(original_problem.g->clone());
   mapping_to_origrow = dynamic_cast<StochVector*>(original_problem.g->clone());
   // count !? rows cols...

   start_indices.push_back(0);
}

StochPostsolver::~StochPostsolver(){}


void StochPostsolver::notifyFixedColumn( int node, unsigned int col, double value)
{
   reductions.push_back(FIXED_COLUMN);
   //indices.push_back( mapping_to_origcol[col]);
//   values.push_back( value );
   throw std::runtime_error("Not yet implemented");
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
   assert( values.size() == start_indices.size() );
   start_indices.push_back( values.size() );
}

// todo primal postsolving only (atm)
StochPostsolver::PostsolveStatus StochPostsolver::undo(const sVars& reduced_solution, sVars& original_solution)
{
   const StochVector& primal_vars_orig = dynamic_cast<const StochVector&>(*reduced_solution.x);



   throw std::runtime_error("Not yet implemented");


   return PRESOLVE_OK;
}


