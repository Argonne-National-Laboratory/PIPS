/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatNBSolver.cxx

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



#include "MatNBSolver.hxx"

#include "Nbundle.hxx"
#include "Nsumproblem.hxx"
#include "Nfunproblem.hxx"
//include "Nlmaxproblem.hxx"
//include "Nsocproblem.hxx"
//include "Nconeproblem.hxx"

#include <algorithm>
#include <map>

//------------------------------------------------------------

 
using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  //------------------------------------------------------------
  // Wrapper for transforming FunctionOracle to MatrixFuncitonOracle
  //------------------------------------------------------------

  class NFunctionOracleWrapper : public MatrixFunctionOracle
  {
  private:
    FunctionOracle& oracle;
    std::ostream* out;
    int print_level;
  public:
    NFunctionOracleWrapper(FunctionOracle& o): oracle(o)
    {out=0;print_level=0;}
    ~NFunctionOracleWrapper(){}

    void set_out(std::ostream* o=0,int pril=1)
    {out=o;print_level=pril;}

    /** see MatrixFunctionOracle for explanations */
    int evaluate(const  Matrix& y,double relprec,double&  objective_value,
		 Matrix&  cutvals,   Matrix&  subgs,
		 std::vector<PrimalData*>& primal_data,
		 PrimalExtender*&primal_extender)
    {
      //---- create c++ interface data
      DVector              duals(y.dim());
      DVector              cut_vals;   
      std::vector<DVector> cut_subgs;
      // copy current duals
      for( int i = 0; i < y.dim(); ++i ){
	duals[i] = y(i);
      }

      //---- evaluate			
      int ret_code = oracle.evaluate(duals,relprec,objective_value,
				     cut_vals,cut_subgs,primal_data,primal_extender);
      
      //---- copy solution value
      cutvals.init(cut_vals);
      if (cut_subgs.size()==0){
	subgs.init(0,0,0.);
	return ret_code;
      }
      subgs.init(y.dim(),Integer(cut_subgs.size()),0.);
      for(int j=0;j<int(cut_subgs.size());j++){
	if (y.dim()!=int((cut_subgs[j]).size())){
	  if (out) (*out)<<"**** ERROR: NFunctionOracleWrapper::eval_function(): function returned subgradient "<<j<<" with dimension "<<(cut_subgs[j]).size()<<" instead of "<<y.dim()<<std::endl;
	  ret_code|=1;
	}
	Integer maxi=CH_Matrix_Classes::min(y.dim(),Integer((cut_subgs[j]).size()));
	for( Integer i = 0; i < maxi; ++i){
	  subgs(i,j) = (cut_subgs[j])[i];
	}
      }
      return ret_code;
    }

    int subgradient_extension(const PrimalData* generating_primal,
			      const Indexmatrix& varind,Matrix& newsubgval)
    {
      IVector vind(varind.dim());
      DVector nval;
      for( int i = 0; i < varind.dim(); ++i ){
	vind[i] = varind(i);
      }
      int ret_code=oracle.subgradient_extension(generating_primal,vind,nval);
      newsubgval.init(nval);
      return ret_code;
    }
  };

  //------------------------------------------------------------
  // Data class 
  //------------------------------------------------------------
  /** we do something non-standard here to keep the interface 'clean':
      SBSolver's data are separated into an extra class MatCBSolverData */
  class MatNBSolverData {
    friend class MatrixNBSolver;

    MatNBSolverData            ( const MatNBSolverData& ); // blocked
    MatNBSolverData& operator= ( const MatNBSolverData& ); // blocked
    
    typedef std::map<const FunctionObject*, NConvexProblem*> FunctionMap;
    typedef std::map<const FunctionObject*, NFunctionOracleWrapper*> WrapperMap;
    
  public:
    //----------------------------------------
    // data
    NBundleSolver bundle;
    NSumProblem   problem;
    Integer      m;
    FunctionMap  funprob;		
    WrapperMap   wrapper;
    int          init;
    int          changed_functions;
    Clock        myclock;
    std::ostream*     out;
    int          print_level;

    void set_out(std::ostream* o=0,int pril=1)
    {
      out=o;print_level=pril;
      bundle.set_out(out,print_level-1);
      problem.set_out(out,print_level-1);
      for( WrapperMap::iterator it = wrapper.begin();
	   it != wrapper.end();
	   ++it ){
	it->second->set_out(out,print_level-1);
      }
    }

    void clear(void)
    {
      problem.clear();
      changed_functions=0;
      for( FunctionMap::iterator it = funprob.begin();
	   it != funprob.end();
	   ++it ){
	delete it->second;
      }
      funprob.clear();
      {for( WrapperMap::iterator it = wrapper.begin();
	   it != wrapper.end();
	   ++it ){
	delete it->second;
      }}
      wrapper.clear();
      init=0;
      m=-1;
      
      bundle.clear();
      bundle.set_clock(myclock);
      myclock.start();
    };

    void set_defaults(void)
    {
      bundle.set_defaults(); 
    }
      
    

    ///
    MatNBSolverData(){ 
      print_level=1;
      out=&std::cout;
      clear(); 
      set_defaults();
    }
    
    ///
    ~MatNBSolverData()
    {
      clear();
    };
    
    
  };


  //------------------------------------------------------------
  // CBmethod implementation - mostly just wrapped to CBmethodData 
  //------------------------------------------------------------

  //--------------------
  MatrixNBSolver::MatrixNBSolver()
  {
    data_ = new MatNBSolverData;
    assert( data_ );
  }

  //--------------------
  MatrixNBSolver::~MatrixNBSolver()
  {
    assert( data_ );
    delete data_;
  }
	
  //--------------------
  void MatrixNBSolver::clear()
  {
    assert( data_ );
    data_->clear();
  }
    
  //--------------------
  void MatrixNBSolver::set_defaults()
  {
    assert( data_ );
    data_->set_defaults();
  }
    
  //--------------------
  int MatrixNBSolver::init_problem(int dim,
				const Matrix* lbounds,
				const Matrix* ubounds,
				const Matrix* costs)
  {
    assert( data_ );
    clear();
    if (data_->bundle.init(data_->problem,dim,lbounds,ubounds,costs)){
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::init_problem(...): "
		     << " bundle.init failed" <<std::endl;
      return 1;
    }
    
    data_->m=dim;

    return 0;
  } 

  //--------------------
  int MatrixNBSolver::add_function( FunctionObject& function )
  {
    assert( data_ );

    if (data_->m<0){
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): call set_dimension first" <<std::endl;
      return 1;
    }
 
    if (data_->funprob.find( &function ) != data_->funprob.end() ) {
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): "
	  << "function already added" <<std::endl;
      return 1;
    }
	
    FunctionOracle* fop=dynamic_cast<FunctionOracle*>(&function);
    if (fop) {
      NFunctionOracleWrapper* o = new NFunctionOracleWrapper(*fop);
      if (o==0) {
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): construction of wrapper interface for function failed" <<std::endl;
	return 1;
      }
      NFunctionProblem* p = new NFunctionProblem(*o);
      if (p==0) {
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): construction of subproblem for function failed" <<std::endl;
	return 1;
      }
      
      data_->init=0;
      if (data_->problem.add_problem(p)){
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): could not add function to SumProblem " <<std::endl;
	delete p;
	return 1;
      }

      data_->wrapper[ &function ] = o;
      data_->funprob[ &function ] = p;
      return 0;
    }

    MatrixFunctionOracle* mfop=dynamic_cast<MatrixFunctionOracle*>(&function);
    if (mfop) {
      NFunctionProblem* p = new NFunctionProblem(*mfop);
      if (p==0) {
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): construction of subproblem for function failed" <<std::endl;
	return 1;
      }
      
      data_->init=0;
      if (data_->problem.add_problem(p)){
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): could not add function to SumProblem " <<std::endl;
	delete p;
	return 1;
      }

      data_->funprob[ &function ] = p;
      return 0;
    }

/*
 
    BaseSDPOracle* sdp=dynamic_cast<BaseSDPOracle*>(&function);
    if (sdp) {
      LmaxProblem* lmaxp=new LmaxProblem(sdp);
      if (lmaxp==0) {
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): construction of subproblen for function failed" <<std::endl;
	return 1;
      }
      
      data_->init=0;
      if (data_->problem.add_problem(lmaxp)){
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): could not add function to SumProblem " <<std::endl;
	delete lmaxp;
	return 1;
      }

      data_->funprob[ &function ] = lmaxp;
      return 0;
    }
      
    BaseSOCOracle* soc=dynamic_cast<BaseSOCOracle*>(&function);
    if (soc) {
      SocProblem* socp=new SocProblem(soc);
      if (socp==0) {
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): construction of subproblen for function failed" <<std::endl;
	return 1;
      }
      
      data_->init=0;
      if (data_->problem.add_problem(socp)){
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): could not add function to SumProblem " <<std::endl;
	delete socp;
	return 1;
      }

      data_->funprob[ &function ] = socp;
      return 0;
    }
      
    BaseConeOracle* con=dynamic_cast<BaseConeOracle*>(&function);
    if (con) {
      ConeProblem* conp=new ConeProblem(con);
      if (conp==0) {
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): construction of subproblen for function failed" <<std::endl;
	return 1;
      }
      
      data_->init=0;
      if (data_->problem.add_problem(conp)){
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): could not add function to SumProblem " <<std::endl;
	delete conp;
	return 1;
      }

      data_->funprob[ &function ] = conp;
      return 0;
    }
      
 */

    if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::add_function(...): unknown derivation of FunctionObject"<<std::endl;
     

    return 1;
  }
	
  
  //----------------------------------------
  // append new variables (always in last postions in this order)
  int MatrixNBSolver::append_variables(int n_append, 
				    const Matrix* lbounds,
				    const Matrix* ubounds,
				    const Matrix* costs)
  {
    assert( data_ );

    if (data_->m<0){
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::append_variables(...): call set_dimension first" <<std::endl;
      return 1;
    }
    if (n_append==0) return 0;
    if (n_append<0) {
      if (data_->out) (*data_->out)<<"**** ERROR: MatrixNBSolver::append_variables(...): n_append <0"<<std::endl;
      return 1;
    }
    if ((lbounds!=0)&&(lbounds->dim()!=n_append)) {
      if (data_->out) (*data_->out)<<"**** ERROR: MatrixNBSolver::append_variables(...): lower bounds vector does not match dimension"<<std::endl;
      return 1;
    }
    if ((ubounds!=0)&&(ubounds->dim()!=n_append)) {
      if (data_->out) (*data_->out)<<"**** ERROR: MatrixNBSolver::append_variables(...): lower bounds vector does not match dimension"<<std::endl;
      return 1;
    }

    //check correctnes of values of lbounds and ubounds
    if ((lbounds!=0)||(ubounds!=0)){
      for (Integer i=0;i<n_append;i++){
	if (lbounds){
	  if((*lbounds)[i]>CB_plus_infinity){
	    if (data_->out) (*data_->out)<<"**** ERROR: MatrixNBSolver::append_variables(...): lower bound of coordinate "<<i<<" exceeds plus_infinity: "<<(*lbounds)[i]<<std::endl;
	    return 1;
	  }
	  if ((*lbounds)[i]==CB_plus_infinity){
	    if (data_->out) (*data_->out)<<"**** WARNING: MatrixNBSolver::append_variables(...): lower bound of coordinate "<<i<<" equals plus_infinity: "<<(*lbounds)[i]<<std::endl;
	  }
	  if ((*lbounds)[i]<CB_minus_infinity){
	    if (data_->out) (*data_->out)<<"**** WARNING: MatrixNBSolver::append_variables(...): lower bound of coordinate "<<i<<" is smaller than minus_infinity: "<<(*lbounds)[i]<<std::endl;
	  }
	}
	if (ubounds){
	  if((*ubounds)[i]<CB_minus_infinity){
	    if (data_->out) (*data_->out)<<"**** ERROR: MatrixNBSolver::append_variables(...): upper bound of coordinate "<<i<<" exceeds minus_infinity: "<<(*ubounds)[i]<<std::endl;
	    return 1;
	  }
	  if ((*ubounds)[i]==CB_minus_infinity){
	    if (data_->out) (*data_->out)<<"**** WARNING: MatrixNBSolver::append_variables(...): upper bound of coordinate "<<i<<" equals minus_infinity: "<<(*ubounds)[i]<<std::endl;
	  }
	  if ((*ubounds)[i]>CB_plus_infinity){
	    if (data_->out) (*data_->out)<<"**** WARNING: MatrixNBSolver::append_variables(...): upper bound of coordinate "<<i<<" exceeds plus_infinity: "<<(*ubounds)[i]<<std::endl;
	  }
	  if ((lbounds)&&((*ubounds)[i]<(*lbounds)[i])){
	    if (data_->out) (*data_->out)<<"**** ERROR: MatrixNBSolver::append_variables(...): lower bound "<<(*lbounds)[i]<<" of coordinate "<<i<<" exceeds upper bound "<<(*ubounds)[i]<<std::endl;
	    return 1;
	  }
	} //endif ubounds
      } //endfor
    } //endif lbounds or ubounds
   
    data_->m+=n_append;
    return data_->bundle.append_variables(n_append,lbounds,ubounds,costs);


  }

  //----------------------------------------
  // delete variables 
  int MatrixNBSolver::delete_variables(const Indexmatrix& del_indices,
				 Indexmatrix& map_to_old)
  {
    assert( data_ );

    if (data_->m<0){
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::delete_function(...): call set_dimension first" <<std::endl;
      return 1;
    }
    
    if (del_indices.dim()==0){
      map_to_old.init(Range(0,data_->m-1));
      return 0;
    }

    Indexmatrix del_vec;
    sortindex(del_indices,del_vec);
    del_vec=del_indices(del_vec);
    if ((del_vec(0)<0)||(del_vec(del_vec.dim()-1)>=data_->m)){
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::delete_function(...): variable indices for deletion out of range" <<std::endl;
      return 1;
    }
    for (Integer i=1;i<del_vec.dim();i++){
      if (del_vec(i-1)==del_vec(i)){
	if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::delete_function(...): multiple indices for deletion, e.g. "<<del_vec(i)<<std::endl;
	return 1;
      }
    }

    map_to_old.init(Range(0,data_->m-1));
    map_to_old.delete_rows(del_vec);


    //call problem with changing information
    if (data_->bundle.reassign_variables(map_to_old)){
      return 1;
    }
    data_->m-=del_vec.dim();

    return 0;
  }
    

  //----------------------------------------
  // reassign variables 
  int MatrixNBSolver::reassign_variables(const Indexmatrix& avec)
  {
    assert( data_ );

    if (data_->m<0){
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::reassign_function(...): call set_dimension first" <<std::endl;
      return 1;
    }
    
    if ((avec.dim()>0)&&((min(avec)<0)||(max(avec)>=data_->m))){
      if (data_->out) (*data_->out)<< "**** ERROR: MatrixNBSolver::reassign_function(...): variable indices out of range" <<std::endl;
      return 1;
    }

    //call problem with changing information
    if (data_->bundle.reassign_variables(avec)){
      return 1;
    }
    data_->m=avec.dim();

    return 0;
  }
    

  //----------------------------------------
  int MatrixNBSolver::set_lower_bound(int i,double lb)
  {
    assert( data_ );
    if ((data_->m<0)||(i<0)||(i>=data_->m)) return 1;

    if (lb<=CB_minus_infinity) return data_->bundle.set_lower_bound(i,CB_minus_infinity);
    return data_->bundle.set_lower_bound(i,lb);
  }

  //----------------------------------------
  int MatrixNBSolver::set_upper_bound(int i,double ub)
  {
    assert( data_ );
    if ((data_->m<0)||(i<0)||(i>=data_->m)) return 1;

    if (ub>=CB_plus_infinity) return data_->bundle.set_upper_bound(i,CB_plus_infinity);
    return data_->bundle.set_upper_bound(i,ub);
  }

  
  //--------------------
  int MatrixNBSolver::do_descent_step(int maxsteps)
  {
    assert( data_ );

    //ensure correct initialization		
    if ((data_->funprob.empty())||(data_->m<0)) {
      if (data_->out){
	(*data_->out)<< "**** ERROR: MatrixNBSolver::do_descent_step(...): problem not initialized:"
		     << " number of functions="<<data_->funprob.size()
		     << " dimension="<<data_->m<<std::endl;
      }
      return 1;
    }
    
    int status=0;
    data_->bundle.get_terminator()->clear_terminated();
    status = data_->bundle.inner_loop(maxsteps);
    
    if ((!status)&&(data_->out)&&(data_->print_level>0)) 
      data_->bundle.print_line_summary(*data_->out);
    
    return status;
  }
  
  //--------------------
  int MatrixNBSolver::termination_code() const
  {
    assert( data_ );
    return data_->bundle.get_terminate();  
  }
  
  //--------------------
  std::ostream& MatrixNBSolver::print_termination_code(std::ostream& out) 
  {
    assert( data_ );
    data_->bundle.get_terminator()->print_status(out);  
    return out;
  }
  
  //--------------------
  double MatrixNBSolver::get_objval() const
  {
    assert( data_ );
    
    return data_->bundle.get_oldval();
  }
  
  //--------------------
  int MatrixNBSolver::get_center(Matrix& y) const
  {
    assert( data_ );
    
    y=data_->bundle.get_y();
    return 0;
  }

  //--------------------
  double MatrixNBSolver::get_candidate_value() const
  {
    assert( data_ );
    
    return data_->bundle.get_newval();
  }
  
  //--------------------
  int MatrixNBSolver::get_candidate(Matrix& y) const
  {
    assert( data_ );
    
    y=data_->bundle.get_newy();
    return 0;
  }

  //--------------------
  int MatrixNBSolver::get_approximate_slacks(Matrix& eta) const
  {
    assert( data_ );
    eta=data_->bundle.get_eta();
    return 0;
  }

  //--------------------
  int MatrixNBSolver::get_approximate_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( data_ );
    
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return data_->funprob[ &function ]->get_approximate_primal(primal); 
  }

  //--------------------
  int MatrixNBSolver::get_center_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( data_ );
    
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return data_->funprob[ &function ]->get_center_primal(primal); 
  }

  //--------------------
  int MatrixNBSolver::get_candidate_primal( const FunctionObject& function,
					PrimalData&           primal ) const
  {		
    assert( data_ );
    
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return data_->funprob[ &function ]->get_candidate_primal(primal); 
  }

  //--------------------
  double MatrixNBSolver::get_sgnorm() const
  {
    assert( data_ );
    return sqrt(data_->bundle.get_normsubg2());
  }
	
  //--------------------
  int MatrixNBSolver::get_subgradient( Matrix& subg) const
  {
    assert( data_ );
    subg=data_->bundle.get_aggr();
    return 0;
  }

  //--------------------
  double MatrixNBSolver::get_cutval() const
  {
    assert( data_ );
    
    return data_->bundle.get_modelval();
  }
  
  //--------------------
  int MatrixNBSolver::get_function_status( const FunctionObject& function ) const
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return data_->funprob[&function]->get_ret_code();
  }
  
  //--------------------
  int MatrixNBSolver::set_max_bundlesize(const FunctionObject& function, 
					 int /* mb */)
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return 0;
  }

  //--------------------
  int MatrixNBSolver::set_max_new_subgradients(const FunctionObject& function, 
					       int /* nnew */)
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return 0;
  }
    
  //--------------------
  int MatrixNBSolver::set_bundle_parameters(const FunctionObject& function,
					    const BundleParameters& /* bp */) 
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return 0;
  }

  //--------------------
  int MatrixNBSolver::get_bundle_parameters(const FunctionObject& function,
				      BundleParameters& bp ) const 
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    bp.n_new_subgradients=1;
    bp.n_bundle_size=2;
    return 0;
  }

  //--------------------
  int MatrixNBSolver::get_bundle_values(const FunctionObject& function,
				      BundleParameters& bp ) const
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    bp.n_new_subgradients=1;
    bp.n_bundle_size=2;
    return 0;
  }

  //--------------------
  int MatrixNBSolver::reinit_function_model( const FunctionObject& function )
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return data_->bundle.reinit_function_model();
  }
			
  //--------------------
  int MatrixNBSolver::clear_aggregates( const FunctionObject& function )
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return data_->bundle.clear_aggregate(); 
  }
		
  //--------------------
  int MatrixNBSolver::call_primal_extender(const FunctionObject& function,PrimalExtender& primal_extender)
  {
    assert( data_ );
    if ( data_->funprob.find( &function ) == data_->funprob.end() )
      return 1;
    return data_->funprob[&function]->call_primal_extender(primal_extender);
  }
			
	
  int MatrixNBSolver::set_term_relprec( const double term_relprec )
  {
    assert( data_ );
    data_->bundle.get_terminator()->set_termeps(term_relprec);
    return 0;
  }
	
  double MatrixNBSolver::get_last_weight() const
  {
    assert( data_ );
    return data_->bundle.get_weight();
  }

  int MatrixNBSolver::set_next_weight( const double weight )
  {
    assert( data_ );
    if (weight<=0) return 1;
    data_->bundle.get_bundleweight()->set_next_weight(weight); 
    return 0;
  }

  int MatrixNBSolver::set_min_weight( const double weight )
  {
    assert( data_ );
    data_->bundle.get_bundleweight()->set_minweight(weight); 
    return 0;
  }
	
  int MatrixNBSolver::set_max_weight( const double weight )
  {
    assert( data_ );
    data_->bundle.get_bundleweight()->set_maxweight(weight); 
    return 0;
  }

  int MatrixNBSolver::set_scaling( bool do_scaling )
  {
    assert( data_ );
    data_->bundle.set_do_scaling(int(do_scaling));
    return 0;
  }


    /** use user defined diagonal scaling of the variables */
  int MatrixNBSolver::set_scaling( const Matrix& scale )
  {
    assert( data_ );
    if (data_->m!=scale.dim()) return 1;
    data_->bundle.set_scaling(scale);
    return 0;
  }
	
  void MatrixNBSolver::set_active_bounds_fixing( bool allow_fixing )
  {
    assert( data_ );
    data_->bundle.set_do_yfixing(allow_fixing); 
  }
    
  void MatrixNBSolver::clear_fail_counts(void )
  {
    assert( data_ );
    data_->bundle.clear_fails(); 
  }
	
  void MatrixNBSolver::set_eval_limit(Integer eval_limit)
  {
    assert( data_ );
    data_->bundle.get_terminator()->set_objevallimit(eval_limit); 
  }
	
  void MatrixNBSolver::set_inner_update_limit(Integer update_limit)
  {
    assert( data_ );
    data_->bundle.set_max_updates(update_limit); 
  }
	
  int MatrixNBSolver::set_new_center_point( const Matrix& center_point )
  {
    assert( data_ );
    if (data_->m!=center_point.dim()) return 1;
    return data_->bundle.set_center(center_point);
  }
  
  int MatrixNBSolver::adjust_multiplier( void )
  {
    assert( data_ );
    return data_->problem.adjust_multiplier();
  }
  
  int MatrixNBSolver::get_dim() const
  {
    assert( data_ );
    return data_->m;
  }

  int MatrixNBSolver::get_n_functions() const
  {
    assert( data_ );
    return int(data_->funprob.size());
  }

  int MatrixNBSolver::get_n_oracle_calls() const
  {
    assert (data_);
    return data_->bundle.get_cntobjeval();
  }

  int MatrixNBSolver::get_n_descent_steps() const
  {
    assert (data_);
    return data_->bundle.get_descent_steps();
  }

  int MatrixNBSolver::get_n_inner_iterations() const
  {
    assert (data_);
    return data_->bundle.get_suminnerit();
  }

  int MatrixNBSolver::get_n_inner_updates() const
  {
    assert (data_);
    return data_->bundle.get_sumupdatecnt();
  }

  const Matrix& MatrixNBSolver::get_lbounds() const
  {
    assert (data_);
    return data_->bundle.get_lby();
  }

  const Matrix& MatrixNBSolver::get_ubounds() const
  {
    assert (data_);
    return data_->bundle.get_uby();
  }

  /** returns the indicator vector of variables temporarily fixed to 
      the center value due to significantly positive multipliers
      for the box constraints which indicates that the corresponding 
      variables would like to stay at their bounds.
      If no variables were fixed, the dimension of the vector is zero.
  */
  const Indexmatrix& MatrixNBSolver::get_active_bounds_indicator() const
  {
    assert (data_);
    return data_->bundle.get_yfixed();
  }

  void MatrixNBSolver::set_out(std::ostream* o,int pril)
  {
    assert( data_ );
    data_->set_out(o,pril);
  }


  std::ostream& MatrixNBSolver::print_line_summary(std::ostream& out) const
  {
    assert( data_ );
    return data_->bundle.print_line_summary(out);
  }


}

