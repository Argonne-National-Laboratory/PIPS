/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/CBSolver.cxx

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



#include "MatFCBSolver.hxx"
#include "MatNBSolver.hxx"

#include <algorithm>
#include <map>

//------------------------------------------------------------

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  //------------------------------------------------------------
  // CBmethod implementation - mostly passing data to MatrixCBSolver
  //------------------------------------------------------------

  //--------------------
  CBSolver::CBSolver(bool no_bundle)
  {
    if (no_bundle)
      solver = new MatrixNBSolver;
    else 
      solver = new MatrixFCBSolver;
    assert( solver );
  }

  //--------------------
  CBSolver::~CBSolver()
  {
    assert( solver );
    delete solver;
  }
	
  //--------------------
  void CBSolver::clear()
  {
    assert( solver );
    solver->clear();
  }
    
  //--------------------
  void CBSolver::set_defaults()
  {
    assert( solver );
    solver->set_defaults();
  }
    
  //--------------------
  int CBSolver::init_problem(int dim,
			     const DVector* lbounds,
			     const DVector* ubounds)
  {
    assert(solver);
    Matrix lb;
    Matrix *lbp;
    if (lbounds==0) lbp=0;
    else {
      lb.init(*lbounds);
      lbp=&lb;
    }
    Matrix ub;
    Matrix *ubp;
    if (ubounds==0) ubp=0;
    else {
      ub.init(*ubounds);
      ubp=&ub;
    }
    return solver->init_problem(dim,lbp,ubp);
  } 

  //--------------------
  int CBSolver::add_function( FunctionObject& function )
  {
    assert( solver );
    return solver->add_function(function);
  }
	
  
  //----------------------------------------
  // append new variables (always in last postions in this order)
  int CBSolver::append_variables(int n_append, 
				 const DVector* lbounds,
				 const DVector* ubounds)
  {
    assert( solver );
    Matrix lb;
    Matrix *lbp;
    if (lbounds==0) lbp=0;
    else {
      lb.init(*lbounds);
      lbp=&lb;
    }
    Matrix ub;
    Matrix *ubp;
    if (ubounds==0) ubp=0;
    else {
      ub.init(*ubounds);
      ubp=&ub;
    }
    return solver->append_variables(n_append,lbp,ubp);
  }

  //----------------------------------------
  // delete variables 
  int CBSolver::delete_variables(const IVector& del_indices,
				 IVector& map_to_old)
  {
    assert( solver );
    Indexmatrix mapold;
    int ret_code=solver->delete_variables(Indexmatrix(del_indices),mapold);
    assign(map_to_old,mapold);
    return ret_code;
  }

    
  //----------------------------------------
  // reassign variables 
  int CBSolver::reassign_variables(const IVector& aind)
  {
    assert( solver );
    return solver->reassign_variables(Indexmatrix(aind));
  }    

  //----------------------------------------
  int CBSolver::set_lower_bound(int i,double lb)
  {
    assert( solver );
    return solver->set_lower_bound(i,lb);
  }

  //----------------------------------------
  int CBSolver::set_upper_bound(int i,double ub)
  {
    assert( solver );
    return solver->set_upper_bound(i,ub);
  }

  
  //--------------------
  int CBSolver::do_descent_step(int maxsteps)
  {
    assert( solver );
    return solver->do_descent_step(maxsteps);
  }
  
  //--------------------
  int CBSolver::termination_code() const
  {
    assert( solver );
    return solver->termination_code();  
  }
  
  //--------------------
  std::ostream& CBSolver::print_termination_code(std::ostream& out) 
  {
    assert( solver );
    return solver->print_termination_code(out);
  }
  
  //--------------------
  double CBSolver::get_objval() const
  {
    assert( solver );
    return solver->get_objval();
  }
  
  //--------------------
  int CBSolver::get_center( DVector& center ) const
  {
    assert( solver );
    Matrix cen;
    int ret_code=solver->get_center(cen);
    assign(center,cen);
    return ret_code;
  }

  //--------------------
  double CBSolver::get_candidate_value() const
  {
    assert( solver );
    return solver->get_candidate_value();
  }
  
  //--------------------
  int CBSolver::get_candidate( DVector& center ) const
  {
    assert( solver );
    Matrix cen;
    int ret_code=solver->get_candidate(cen);
    assign(center,cen);
    return ret_code;
  }

  //--------------------
  int CBSolver::get_approximate_slacks( DVector& slacks )const
  {
    assert( solver );
    Matrix sl;
    int ret_code=solver->get_approximate_slacks(sl);
    assign(slacks,sl);
    return ret_code;
  }

  //--------------------
  int CBSolver::get_approximate_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( solver );
    return solver->get_approximate_primal(function,primal);
  }
    

  //--------------------
  int CBSolver::get_center_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( solver );
    return solver->get_center_primal(function,primal);
  }

  //--------------------
  int CBSolver::get_candidate_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( solver );
    return solver->get_candidate_primal(function,primal);
  }

  //--------------------
  double CBSolver::get_sgnorm() const
  {
    assert( solver );
    return solver->get_sgnorm();
  }
	
  //--------------------
  int CBSolver::get_subgradient( DVector& subgradient) const
  {
    assert( solver );
    Matrix subg;
    int ret_code=solver->get_subgradient(subg);
    assign(subgradient,subg);
    return ret_code;
  }

  //--------------------
  int CBSolver::get_function_status( const FunctionObject& function ) const
  {
    assert( solver );
    return solver->get_function_status(function);
  }
  
  //--------------------
  int CBSolver::set_max_bundlesize(const FunctionObject& function, 
				   int mb)
  {
    assert( solver );
    return solver->set_max_bundlesize(function,mb);
  }

  //--------------------
  int CBSolver::set_max_new_subgradients(const FunctionObject& function, 
					 int nnew )
  {
    assert( solver );
    return solver->set_max_new_subgradients(function,nnew);
  }
    
  //--------------------
  int CBSolver::set_bundle_parameters(const FunctionObject& function,
				      const BundleParameters& bp ) 
  {
    assert( solver );
    return solver->set_bundle_parameters(function,bp);
  }

  //--------------------
  int CBSolver::get_bundle_parameters(const FunctionObject& function,
				      BundleParameters& bp ) const 
  {
    assert( solver );
    return solver->get_bundle_parameters(function,bp);
  }

  //--------------------
  int CBSolver::get_bundle_values(const FunctionObject& function,
				      BundleParameters& bp ) const
  {
    assert( solver );
    return solver->get_bundle_values(function,bp);
  }

  //--------------------
  int CBSolver::reinit_function_model( const FunctionObject& function )
  {
    assert( solver );
    return solver->reinit_function_model(function);
  }
			
  //--------------------
  int CBSolver::call_primal_extender(const FunctionObject& function,PrimalExtender& primal_extender)
  {
    assert( solver );
    return solver->call_primal_extender(function,primal_extender);
  }
			
  //--------------------
  int CBSolver::set_term_relprec( const double term_relprec )
  {
    assert( solver );
    solver->set_term_relprec(term_relprec);
    return 0;
  }
	
  //--------------------
  double CBSolver::get_last_weight() const
  {
    assert( solver );
    return solver->get_last_weight();
  }

  //--------------------
  int CBSolver::set_next_weight( const double weight )
  {
    assert( solver );
    return solver->set_next_weight(weight); 
  }

  //--------------------
  int CBSolver::set_min_weight( const double weight )
  {
    assert( solver );
    return solver->set_min_weight(weight); 
  }
	
  //--------------------
  int CBSolver::set_max_weight( const double weight )
  {
    assert( solver );
    return solver->set_max_weight(weight); 
  }
	
  //--------------------
  int CBSolver::set_new_center_point( const DVector& center_point )
  {
    assert( solver );
    return solver->set_new_center_point(Matrix(center_point));
  }
  
  //--------------------
  int CBSolver::set_scaling( bool do_scaling )
  {
    assert( solver );
    return solver->set_scaling(int(do_scaling));
  }

  //--------------------
  void CBSolver::set_active_bounds_fixing( bool allow_fixing )
  {
    assert( solver );
    solver->set_active_bounds_fixing(allow_fixing);
  }
    
	
  void CBSolver::clear_fail_counts(void )
  {
    assert( solver );
    return solver->clear_fail_counts(); 
  }
	
  void CBSolver::set_eval_limit(int eval_limit)
  {
    assert( solver );
    return solver->set_eval_limit(eval_limit);
  }
	
  void CBSolver::set_inner_update_limit(int update_limit)
  {
    assert( solver );
    return solver->set_inner_update_limit(update_limit); 
  }
	
  int CBSolver::get_dim()
  {
    assert( solver );
    return solver->get_dim();
  }

  int CBSolver::get_n_functions()
  {
    assert( solver );
    return solver->get_n_functions();
  }

  int CBSolver::get_active_bounds_indicator(IVector& vind) const
  {
    assert( solver );
    const Indexmatrix& ind=solver->get_active_bounds_indicator();
    vind.resize(ind.dim());
    for(Integer i=0;i<ind.dim();i++){
      vind[i]=ind[i];
    }
    return 0;
  }


  void CBSolver::set_out(std::ostream* o,int pril)
  {
    assert( solver );
    solver->set_out(o,pril);
  }

}

