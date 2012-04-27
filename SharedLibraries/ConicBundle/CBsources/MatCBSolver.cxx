/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatCBSolver.cxx

    This file is part of ConciBundle, a C/C++ library for convex optimization.

    ConicBundle is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ConicBundle is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

***************************************************************************** */



#include "MatCBSolver.hxx"
#include "MatFCBSolver.hxx"
#include "MatNBSolver.hxx"

#include "bundle.hxx"
#include "sumproblem.hxx"
#include "funproblem.hxx"
#include "lmaxproblem.hxx"
#include "socproblem.hxx"
#include "coneproblem.hxx"

#include <algorithm>
#include <map>

//------------------------------------------------------------

 
using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {


  //------------------------------------------------------------
  // CBmethod implementation - mostly just wrapped to CBmethodData 
  //------------------------------------------------------------

  //--------------------
  MatrixCBSolver::MatrixCBSolver(bool no_bundle)
  {
    if (no_bundle)
      solver = new MatrixNBSolver;
    else
      solver = new MatrixFCBSolver;
    
    assert( solver );
  }


  //--------------------
  MatrixCBSolver::~MatrixCBSolver()
  {
    assert( solver );
    delete solver;
  }
	

  //--------------------
  void MatrixCBSolver::clear()
  {
    assert( solver );
    solver->clear();
  }
    

  //--------------------
  void MatrixCBSolver::set_defaults()
  {
    assert( solver );
    solver->set_defaults();
  }
    

  //--------------------
  int MatrixCBSolver::init_problem(int dim,
				const Matrix* lbounds,
				const Matrix* ubounds,
				const Matrix* costs)
  {
    assert( solver );
    return solver->init_problem(dim,lbounds,ubounds,costs);
  } 


  //--------------------
  int MatrixCBSolver::add_function( FunctionObject& function )
  {
    assert( solver );
    return solver->add_function(function);
  }
	
  
  //----------------------------------------
  // append new variables (always in last postions in this order)
  int MatrixCBSolver::append_variables(int n_append, 
				    const Matrix* lbounds,
				    const Matrix* ubounds,
				    const Matrix* costs)
  {
    assert( solver );
    return solver->append_variables(n_append,lbounds,ubounds,costs);
  }


  //----------------------------------------
  // delete variables 
  int MatrixCBSolver::delete_variables(const Indexmatrix& del_indices,
				 Indexmatrix& map_to_old)
  {
    assert( solver );
    return solver->delete_variables(del_indices,map_to_old);
  }

    
  //----------------------------------------
  // reassign variables 
  int MatrixCBSolver::reassign_variables(const Indexmatrix& avec)
  {
    assert( solver );
    return solver->reassign_variables(avec);
  }
    

  //----------------------------------------
  int MatrixCBSolver::set_lower_bound(int i,double lb)
  {
    assert( solver );
    return solver->set_lower_bound(i,lb);
  }


  //----------------------------------------
  int MatrixCBSolver::set_upper_bound(int i,double ub)
  {
    assert( solver );
    return solver->set_upper_bound(i,ub);
  }

  
  //--------------------
  int MatrixCBSolver::do_descent_step(int maxsteps)
  {
    assert( solver );
    return solver->do_descent_step(maxsteps);
  }
  

  //--------------------
  int MatrixCBSolver::termination_code() const
  {
    assert( solver );
    return solver->termination_code();  
  }
  

  //--------------------
  std::ostream& MatrixCBSolver::print_termination_code(std::ostream& out) 
  {
    assert( solver );
    return solver->print_termination_code(out);
  }
  

  //--------------------
  double MatrixCBSolver::get_objval() const
  {
    assert( solver );
    return solver->get_objval();
  }
  

  //--------------------
  int MatrixCBSolver::get_center(Matrix& y) const
  {
    assert( solver );
    return solver->get_center(y);
  }


  //--------------------
  double MatrixCBSolver::get_candidate_value() const
  {
    assert( solver );
    return solver->get_candidate_value();
  }
  

  //--------------------
  int MatrixCBSolver::get_candidate(Matrix& y) const
  {
    assert( solver );
    return solver->get_candidate(y);
  }


  //--------------------
  int MatrixCBSolver::get_approximate_slacks(Matrix& eta) const
  {
    assert( solver );
    return solver->get_approximate_slacks(eta);
  }

  //--------------------
  int MatrixCBSolver::get_approximate_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( solver );
    return solver->get_approximate_primal(function,primal);
  }


  //--------------------
  int MatrixCBSolver::get_center_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( solver );
    return solver->get_center_primal(function,primal);
  }

  //--------------------
  int MatrixCBSolver::get_candidate_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( solver );
    return solver->get_candidate_primal(function,primal);
  }

  //--------------------
  double MatrixCBSolver::get_sgnorm() const
  {
    assert( solver );
    return solver->get_sgnorm();
  }
	
  //--------------------
  int MatrixCBSolver::get_subgradient( Matrix& subg) const
  {
    assert( solver );
    return solver->get_subgradient(subg);
  }

  //--------------------
  double MatrixCBSolver::get_cutval() const
  {
    assert( solver );
    return solver->get_cutval();
  }
  
  //--------------------
  int MatrixCBSolver::get_function_status( const FunctionObject& function ) const
  {
    assert( solver );
    return solver->get_function_status(function);
  }
  
  //--------------------
  int MatrixCBSolver::set_max_bundlesize(const FunctionObject& function, 
				   int mb)
  {
    assert( solver );
    return solver->set_max_bundlesize(function,mb);
  }

  //--------------------
  int MatrixCBSolver::set_max_new_subgradients(const FunctionObject& function, 
					 int nnew )
  {
    assert( solver );
    return solver->set_max_new_subgradients(function,nnew);
  }
    
  //--------------------
  int MatrixCBSolver::set_bundle_parameters(const FunctionObject& function,
				      const BundleParameters& bp ) 
  {
    assert( solver );
    return solver->set_bundle_parameters(function,bp);
  }

  //--------------------
  int MatrixCBSolver::get_bundle_parameters(const FunctionObject& function,
				      BundleParameters& bp ) const 
  {
    assert( solver );
    return solver->get_bundle_parameters(function,bp);
  }

  //--------------------
  int MatrixCBSolver::get_bundle_values(const FunctionObject& function,
				      BundleParameters& bp ) const
  {
    assert( solver );
    return solver->get_bundle_values(function,bp);
  }

  //--------------------
  int MatrixCBSolver::reinit_function_model( const FunctionObject& function )
  {
    assert( solver );
    return solver->reinit_function_model(function);
  }
			
  //--------------------
  int MatrixCBSolver::clear_aggregates( const FunctionObject& function )
  {
    assert( solver );
    return solver->clear_aggregates(function);
  }
			
  //--------------------
  int MatrixCBSolver::call_primal_extender(const FunctionObject& function,PrimalExtender& primal_extender)
  {
    assert( solver );
    return solver->call_primal_extender(function,primal_extender);
  }
			
  //--------------------
  int MatrixCBSolver::set_term_relprec( const double term_relprec )
  {
    assert( solver );
    return solver->set_term_relprec(term_relprec);
  }
	
  //--------------------
  double MatrixCBSolver::get_last_weight() const
  {
    assert( solver );
    return solver->get_last_weight();
  }

  //--------------------
  int MatrixCBSolver::set_next_weight( const double weight )
  {
    assert( solver );
    return solver->set_next_weight(weight); 
  }

  //--------------------
  int MatrixCBSolver::set_min_weight( const double weight )
  {
    assert( solver );
    return solver->set_min_weight(weight); 
  }
	
  //--------------------
  int MatrixCBSolver::set_max_weight( const double weight )
  {
    assert( solver );
    return solver->set_max_weight(weight); 
  }

  //--------------------
  int MatrixCBSolver::set_scaling( bool do_scaling )
  {
    assert( solver );
    return solver->set_scaling(int(do_scaling));
  }


  //--------------------
    /** use user defined diagonal scaling of the variables */
  int MatrixCBSolver::set_scaling( const Matrix& scale )
  {
    assert( solver );
    return solver->set_scaling(scale);
  }

  //--------------------
  void MatrixCBSolver::set_active_bounds_fixing( bool allow_fixing )
  {
    assert( solver );
    solver->set_active_bounds_fixing(allow_fixing);
  }
	
  //--------------------
  void MatrixCBSolver::clear_fail_counts(void )
  {
    assert( solver );
    return solver->clear_fail_counts(); 
  }
	
  //--------------------
  void MatrixCBSolver::set_eval_limit(Integer eval_limit)
  {
    assert( solver );
    return solver->set_eval_limit(eval_limit);
  }
	
  //--------------------
  void MatrixCBSolver::set_inner_update_limit(Integer update_limit)
  {
    assert( solver );
    return solver->set_inner_update_limit(update_limit); 
  }
	
  //--------------------
  int MatrixCBSolver::set_new_center_point( const Matrix& center_point )
  {
    assert( solver );
    return solver->set_new_center_point(center_point);
  }
  
  //--------------------
  int MatrixCBSolver::adjust_multiplier( void )
  {
    assert( solver );
    return solver->adjust_multiplier();
  }
  
  //--------------------
  int MatrixCBSolver::get_dim() const
  {
    assert( solver );
    return solver->get_dim();
  }

  //--------------------
  int MatrixCBSolver::get_n_functions() const
  {
    assert( solver );
    return solver->get_n_functions();
  }

  //--------------------
  int MatrixCBSolver::get_n_oracle_calls() const
  {
    assert (solver);
    return solver->get_n_oracle_calls();
  }

  //--------------------
  int MatrixCBSolver::get_n_descent_steps() const
  {
    assert (solver);
    return solver->get_n_descent_steps();
  }

  //--------------------
  int MatrixCBSolver::get_n_inner_iterations() const
  {
    assert (solver);
    return solver->get_n_inner_iterations();
  }

  //--------------------
  int MatrixCBSolver::get_n_inner_updates() const
  {
    assert (solver);
    return solver->get_n_inner_updates();
  }

  //--------------------
  const Matrix& MatrixCBSolver::get_lbounds() const
  {
    assert (solver);
    return solver->get_lbounds();
  }

  //--------------------
  const Matrix& MatrixCBSolver::get_ubounds() const
  {
    assert (solver);
    return solver->get_ubounds();
  }

  //--------------------
  const Indexmatrix& MatrixCBSolver::get_active_bounds_indicator() const
  {
    assert (solver);
    return solver->get_active_bounds_indicator();
  }

  //--------------------
  void MatrixCBSolver::set_out(std::ostream* o,int pril)
  {
    assert( solver );
    solver->set_out(o,pril);
  }


  //--------------------
  std::ostream& MatrixCBSolver::print_line_summary(std::ostream& out) const
  {
    assert( solver );
    return solver->print_line_summary(out);
  }


}

