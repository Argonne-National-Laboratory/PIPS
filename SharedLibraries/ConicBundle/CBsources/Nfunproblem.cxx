/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/Nfunproblem.cxx

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



#include "mymath.hxx"
#include "Nfunproblem.hxx"

using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                            NFunctionProblem
// *****************************************************************************

NFunctionProblem::NFunctionProblem(MatrixFunctionOracle& fo):
  oracle(fo)
{
  center_x=0;
  new_x=0;
  aggr_x=0;
  ret_code=0;

  out=0;
  print_level=0;

  evaltime=Microseconds(0);
  nr_eval=0;
}

// *****************************************************************************
//                            NFunctionProblem
// *****************************************************************************

NFunctionProblem::~NFunctionProblem()
{
  delete center_x;
  delete new_x;
  delete aggr_x;
}

// *****************************************************************************
//                            eval_function
// *****************************************************************************

//evaluates the objective function in $y$
//if evaluated by an iterative method that provides upper and lower bounds,
//it may stop when the lower bound (lb) is above the nullstep_bound
//evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied
//returns:  0 ... if all is ok, use get_function_sol or do_step
//          1 ... if solution could not be computed to desired precision

int NFunctionProblem::eval_function(const Matrix& iny,
				    Real nullstep_bound,Real& lb, Real& ub, Matrix* subgp,Real relprec)
{
  if (subgp!=0){
    delete new_x; 
    new_x=0;
  }
  std::vector<PrimalData*> cand_x;
  PrimalExtender* pep=0;
  Matrix cand_subg_valvec(0,0,0.);
  Matrix cand_subg(0,0,0.);
			

  //---- evaluate  	
  ub=nullstep_bound;
  ret_code = oracle.evaluate(iny,relprec,ub,
			     cand_subg_valvec,cand_subg,cand_x,pep);
  nr_eval++;

  if (cand_subg_valvec.dim()==0) {
    if (out) (*out)<<"**** ERROR: NFunctionProblem::eval_function(): function returned no subgradient values"<<std::endl;
    if (subgp){
      subgp->init(0,0,0.);
    }
    return 1;
  }

  if (cand_subg_valvec.dim()!=cand_subg.coldim()){
    if (out) (*out)<<"**** ERROR: NFunctionProblem::eval_function(): function returned different number of subgradient values and subgradients"<<std::endl;
    if (subgp){
      subgp->init(0,0,0.);
    }
    return 1;
  }

  if (iny.dim()!=cand_subg.rowdim()){
    if (out) (*out)<<"**** ERROR: NFunctionProblem::eval_function(): function returned subgradients of wrong dimension"<<std::endl;
    if (subgp){
      subgp->init(0,0,0.);
    }
    return 1;
  }

  if (ub<max(cand_subg_valvec)){
    if (out) (*out)<<"**** ERROR: NFunctionProblem::eval_function(): function returned an upper bound="<<ub<<" smaller than the lower bound="<<max(cand_subg_valvec)<<" on the function value"<<std::endl;
    if (subgp){
      subgp->init(0,0,0.);
    }
    return 1;
  }

  if ((cand_x.size()>0)&&(Integer(cand_x.size())!=cand_subg.coldim())){
    if (out) (*out)<<"**** ERROR: NFunctionProblem::eval_function(): function returned different number of primal solutions and subgradients; ignoring primal solutions"<<std::endl;
    for(unsigned int i=0;i<cand_x.size();i++){
      delete cand_x[i];
    }
    cand_x.clear();
  }
  
  Integer cand_subg_maxind;
  lb=max(cand_subg_valvec,&cand_subg_maxind);
  if (subgp){
    *subgp=cand_subg.col(cand_subg_maxind);

    if (cand_x.size()>0){
      new_x=cand_x[cand_subg_maxind];
      cand_x[cand_subg_maxind]=0;
      for(Integer i=Integer(cand_x.size());--i>=0;){
	delete cand_x[i];
      }
      cand_x.clear();
    }
  }

  if (pep){
    if (center_x){
      if (pep->extend(*center_x)){
	if (out) (*out)<<"**** WARNING: NFunctionProblem::eval_function(): PrimalExtender::extend failed"<<std::endl;
      }
    }
    if (aggr_x){
      if (pep->extend(*aggr_x)){
	if (out) (*out)<<"**** WARNING: NFunctionProblem::eval_function(): PrimalExtender::extend failed"<<std::endl;
      }
    }
    delete pep;
  }
    

  if (ret_code){
    if (out) (*out)<<"**** WARNING: NFunctionProblem::eval_function(): function returned code = "<<ret_code<<std::endl;
    return ret_code;
  }

  if ((lb<nullstep_bound)&&
      (ub-lb>relprec*(fabs(ub)+1.))){
    if (out) (*out)<<"**** WARNING: NFunctionProblem::eval_function(): insufficient precision in evaluation routine"<<std::endl;
    return 1;
  }    
  
  return 0;
}


// *****************************************************************************
//                              do_step
// *****************************************************************************


void NFunctionProblem::do_step(void)
{
  delete center_x;
  if (new_x)
    center_x=new_x->clone_primal_data();
  else 
    center_x=0;
}

// *****************************************************************************
//                              init_aggregate
// *****************************************************************************


void NFunctionProblem::init_aggregate(void)
{
  delete aggr_x;
  if (new_x)
    aggr_x=new_x->clone_primal_data();
  else 
    aggr_x=0;
}


// *****************************************************************************
//                              aggregate
// *****************************************************************************


void NFunctionProblem::aggregate(Real alpha)
{
  if ((aggr_x)&&(new_x))
    aggr_x->aggregate_primal_data(1-alpha,alpha,*new_x);
  else{
    delete aggr_x;
    aggr_x=0;
  }
}

// *****************************************************************************
//                              clear_aggregate
// *****************************************************************************


void NFunctionProblem::clear_aggregate(void)
{
  delete aggr_x; aggr_x=0;
}

// *****************************************************************************
//                          
// *****************************************************************************


int NFunctionProblem::subgradient_information(const Indexmatrix& indices,Matrix& missing_aggr_coords,Matrix& missing_new_subg_coords)
{
  if ((aggr_x==0)||(new_x==0)){
    return 1;
  }
  missing_aggr_coords.init(indices.dim(),1,0.);
  if (oracle.subgradient_extension(aggr_x,indices,missing_aggr_coords)){
    if (out) (*out)<<"**** WARNING: NFunctionProblem::subgradient_information(): subgradient extension failed for the aggregate"<<std::endl;
    return 1;
  }
  if (missing_aggr_coords.dim()!=indices.dim()){
    if (out) (*out)<<"**** WARNING: NFunctionProblem::subgradient_information(): oracle.subgradient_extension(...) returns a vector of wrong dimension for the aggregate, subgradient extension failed"<<std::endl;
    return 1;
  }
  missing_new_subg_coords.init(indices.dim(),1,0.);
  if (oracle.subgradient_extension(new_x,indices,missing_new_subg_coords)){
    if (out) (*out)<<"**** WARNING: NFunctionProblem::subgradient_information(): subgradient extension failed for the new subgradient"<<std::endl;
    return 1;
  }
  if (missing_new_subg_coords.dim()!=indices.dim()){
    if (out) (*out)<<"**** WARNING: NFunctionProblem::subgradient_information(): oracle.subgradient_extension(...) returns a vector of wrong dimension for the new subgradient, subgradient extension failed"<<std::endl;
    return 1;
  }
  return 0;
}
  

// *****************************************************************************
//                                 lb_function
// *****************************************************************************

//returns a *quick* lower bound for the function value at y
//(eg by a previous subgradient)

Real NFunctionProblem::lb_function(const Matrix& /* iny */)
{
  return CB_minus_infinity;
}

 
// *****************************************************************************
//                                get_approximate_primal
// *****************************************************************************

int NFunctionProblem::get_approximate_primal(PrimalData& primal) const 
  {
    if (aggr_x==0)
      return 1;
    return primal.assign_primal_data(*aggr_x);
  }

// *****************************************************************************
//                                call_primal_extender
// *****************************************************************************

int NFunctionProblem::call_primal_extender(PrimalExtender& prex) 
{
  bool used=false;   
  if (center_x!=0){
    used=true;
    if (prex.extend(*center_x)){
      if (out) (*out)<<"**** WARNING: NFunctionProblem::call_primal_extender(): PrimalExtender::extend failed"<<std::endl;
    }
  }
    
  if (new_x!=0){
    used=true;
    if (prex.extend(*new_x)){
      if (out) (*out)<<"**** WARNING: NFunctionProblem::call_primal_extender(): PrimalExtender::extend failed"<<std::endl;
    }
  }
    
  if (aggr_x!=0){
    used=true;
    if (prex.extend(*aggr_x)){
      if (out) (*out)<<"**** WARNING: NFunctionProblem::call_primal_extender(): PrimalExtender::extend failed"<<std::endl;
    }
  }
  if (!used) return 2;
  return 0; 
}    


} 
