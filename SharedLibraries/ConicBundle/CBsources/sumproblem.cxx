/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/sumproblem.cxx

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
#include "sumproblem.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              clear()
// *****************************************************************************

void SumProblem::clear(void)
{
  subprobp.clear();  
  b.init(0,0,0.);               
  bounds_index.init(0,0,Integer(0));  
  lby.init(0,0,0.);              
  uby.init(0,0,0.);             

  ncalls=0;

  problem_modified=0;   
  center_available=0;
  y.init(0,0,0.);
  subg.init(0,0,0.);
  subg_val=0.;
  ub_fun_val=0.;
  cand_available=0;
  cand_y.init(0,0,0.);
  cand_subg.init(0,0,0.);
  cand_subg_val=0.;
  cand_ub_fun_val=0.;
  model_available=0;
  model_subg_available=0;
  model_subg.init(0,0,0.);
  model_subg_val=0.;
  model_ub_val=0.;
  aug_available=0;
  aug_ind.init(0,0,Integer(0)); 
  aug_subg.init(0,0,0.);
  aug_linconst=0.; 

  delete solverp;
  solverp=0;
  block.clear_blocks();
}



// *****************************************************************************
//                              set_data
// *****************************************************************************

  int SumProblem::set_data(Integer dim,const Matrix* lbounds, const Matrix* ubounds,const Matrix* costs)
  {
    if (dim<0) {
      if (out) (*out)<<"**** ERROR: SumProblem::set_data(...): dimension <0"<<std::endl;
      return 1;
    }
    if ((costs!=0)&&(costs->dim()!=dim)) {
      if (out) (*out)<<"**** ERROR: SumProblem::set_data(...): costs vector does not match dimension"<<std::endl;
      return 1;
    }
    if ((lbounds!=0)&&(lbounds->dim()!=dim)) {
      if (out) (*out)<<"**** ERROR: SumProblem::set_data(...): lower bounds vector does not match dimension"<<std::endl;
      return 1;
    }
    if ((ubounds!=0)&&(ubounds->dim()!=dim)) {
      if (out) (*out)<<"**** ERROR: SumProblem::set_data(...): upper bounds vector does not match dimension"<<std::endl;
      return 1;
    }
    
    //check correctnes of values of lbounds and ubounds
    if ((lbounds!=0)||(ubounds!=0)){
      for (Integer i=0;i<dim;i++){
	if (lbounds){
	  if((*lbounds)[i]>CB_plus_infinity){
	    if (out) (*out)<<"**** ERROR: SumProblem::set_data(...): lower bound of coordinate "<<i<<" exceeds plus_infinity: "<<(*lbounds)[i]<<std::endl;
	    return 1;
	  }
	  if ((*lbounds)[i]==CB_plus_infinity){
	    if (out) (*out)<<"**** WARNING: SumProblem::set_data(...): lower bound of coordinate "<<i<<" equals plus_infinity: "<<(*lbounds)[i]<<std::endl;
	  }
	  if ((*lbounds)[i]<CB_minus_infinity){
	    if (out) (*out)<<"**** WARNING: SumProblem::set_data(...): lower bound of coordinate "<<i<<" is smaller than minus_infinity: "<<(*lbounds)[i]<<std::endl;
	  }
	}
	if (ubounds){
	  if((*ubounds)[i]<CB_minus_infinity){
	    if (out) (*out)<<"**** ERROR: SumProblem::set_data(...): upper bound of coordinate "<<i<<" exceeds plus_infinity: "<<(*ubounds)[i]<<std::endl;
	    return 1;
	  }
	  if ((*ubounds)[i]==CB_minus_infinity){
	    if (out) (*out)<<"**** WARNING: SumProblem::set_data(...): upper bound of coordinate "<<i<<" equals plus_infinity: "<<(*ubounds)[i]<<std::endl;
	  }
	  if ((*ubounds)[i]>CB_plus_infinity){
	    if (out) (*out)<<"**** WARNING: SumProblem::set_data(...): upper bound of coordinate "<<i<<" exceeds plus_infinity: "<<(*ubounds)[i]<<std::endl;
	  }
	  if ((lbounds)&&((*ubounds)[i]<(*lbounds)[i])){
	    if (out) (*out)<<"**** ERROR: SumProblem::set_data(...): lower bound "<<(*lbounds)[i]<<" of coordinate "<<i<<" exceeds upper bound "<<(*ubounds)[i]<<std::endl;
	    return 1;
	  }
	} //endif ubounds
      } //endfor
    } //endif lbounds or ubounds
	    
    if (costs) b=*costs;
    else b.init(dim,1,0.);

    if ((lbounds==0)&&(ubounds==0)){
      bounds_index.init(0,0,Integer(0));
      lby.init(dim,1,CB_minus_infinity);
      uby.init(dim,1,CB_plus_infinity);
    }
    else if (ubounds==0){
      bounds_index.newsize(dim,1);
      bounds_index.init(0,0,Integer(0));
      lby.newsize(dim,1); chk_set_init(lby,1);
      uby.init(dim,1,CB_plus_infinity);
      for (Integer i=0;i<dim;i++){
	if ((lby(i)=CH_Matrix_Classes::max((*lbounds)[i],CB_minus_infinity))
	    >CB_minus_infinity){
	  bounds_index.concat_below(i);
	}
      }
    }
    else if (lbounds==0){
      bounds_index.newsize(dim,1);
      bounds_index.init(0,0,Integer(0));
      uby.newsize(dim,1); chk_set_init(uby,1);
      lby.init(dim,1,CB_minus_infinity);
      for (Integer i=0;i<dim;i++){
        if ((uby(i)=CH_Matrix_Classes::min((*ubounds)[i],CB_plus_infinity))
	    <CB_plus_infinity){
	  bounds_index.concat_below(i);
	}
      }
    }
    else {
      bounds_index.newsize(dim,1);
      bounds_index.init(0,0,Integer(0));
      lby.newsize(dim,1); chk_set_init(lby,1);
      uby.newsize(dim,1); chk_set_init(uby,1);
      for (Integer i=0;i<dim;i++){
	uby(i)=CH_Matrix_Classes::min((*ubounds)[i],CB_plus_infinity);
	lby(i)=CH_Matrix_Classes::max((*lbounds)[i],CB_minus_infinity);
	if ((lby[i]>CB_minus_infinity)||(uby(i)<CB_plus_infinity)){
	  bounds_index.concat_below(i);
	}
      }
    }

    for(unsigned int i=0;i<subprobp.size();i++){
      if (subprobp[i]->intersect_box(bounds_index,lby,uby)){
	return 1;
      }
    }
    
    y.init(0,0,0.);
    problem_modified=1;   
    center_available=0;
    cand_available=0;
    model_available=0;
    model_subg_available=0;
    aug_available=0;
    
    return 0;
  }

// *****************************************************************************
//                            add_problem
// *****************************************************************************

int SumProblem::add_problem(ConvexProblem* p)
{
  if (p->intersect_box(bounds_index,lby,uby)){
    if (out) (*out)<<"**** WARNING: SumProblem::add_problem(...): intersect_box failed (possibly infeasible)"<<std::endl;
    return 1;
  }
  subprobp.push_back(p);
  p->set_out(out,print_level);
  if (b.dim()==0){
    b.init(lby.dim(),1,0.);
  }

  problem_modified=1;   
  center_available=0;
  cand_available=0;
  model_available=0;
  model_subg_available=0;
  aug_available=0;

  return 0;
}

// *****************************************************************************
//                            init_center
// *****************************************************************************

//returns information about the starting point/current center of stability y; 
//if, on input, init=
//   0 ... initialize all variables regardless of content
//   1 ... initialize only if problem was modified since last call 
//on output init_flag=0 signals "variables not modified", =1 if modified.
//lb gives a lower bound, ub an upper bound on the objective value in y
//if yp is not 0 then the current center of stability y is returned
//if subgp is not 0 then an (arbitrary) eps-subgradient at y is returned,
//the hyperplane corresponding to the eps-subgradient has value lb in y
//(lb and subg do not have to be the same as in eval_function)
//bounds_index holds the indices of y that are bounded (dim=#bounded)
//lby and uby give the lower and upper bounds on y 
//(if bounds exist at all then each dim = y.dim())
//returns: 0 ... if the information is available
//         1 ... if the desired information is not available

int SumProblem::init_center(
	       int& init, Real& lb, Real& ub, Matrix& o_y, Matrix& o_subg,
               Indexmatrix& o_bounds_index, Matrix& o_lby, Matrix& o_uby)
{
  if (!center_available){
    if (out) (*out)<<"**** WARNING: SumProblem::init_center(...): no center available"<<std::endl;
    return 1;
  }
  if ((init==0)||(problem_modified)){
    problem_modified=0;
    lb=subg_val;
    ub=ub_fun_val;
    o_y=y;
    o_subg=subg;
    o_bounds_index=bounds_index;
    o_lby=lby;
    o_uby=uby;
    init=1;
    return 0;
  }
  init=0;
  return 0;
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

int SumProblem::eval_function(const Matrix& iny,
			      Real nullstep_bound,Real relprec)
{
  ncalls++;
  cand_available=0;
  cand_y=iny;
  Real lower_bound=ip(cand_y,b);
  cand_subg_val=lower_bound;
  cand_ub_fun_val=lower_bound;
  
  if (subprobp.size()==1){    //just one subproblem

    //compute value for subproblems 
    int retval1=subprobp[0]->eval_function(cand_y,
			       nullstep_bound-lower_bound,
			       relprec);

    //collect solution informaton
    Real lsg_val;
    Real lub_val;
    int retval2=subprobp[0]->get_function_sol(lsg_val,lub_val,&cand_subg);
    if (retval2){
      return retval1;
    } 
    cand_subg_val+=lsg_val;
    cand_ub_fun_val+=lub_val;
    cand_subg+=b; 
    cand_available=1;

    return retval1;
  }

  //----  nonnegative sum of several subproblems

  //compute a quick lower bound on each subfunction
  Matrix fun_lb(Integer(subprobp.size()),1,0.);   
  for (unsigned int i=0;i<subprobp.size();i++){
    fun_lb(i)=subprobp[i]->lb_function(cand_y);
  }
  lower_bound+=sum(fun_lb);
  
  //compute values for subproblems
  int retval1=0;    
  int retval2=0;
  Matrix tmpvec;
  cand_subg=b;
  {for (unsigned int i=0;i<subprobp.size();i++){
    Real corr=lower_bound-fun_lb(i);
    retval1|=subprobp[i]->eval_function(cand_y,
					(nullstep_bound-corr),
					relprec);
    //collect solution information
    if (retval2) continue;
    Real lsg_val;
    Real lub_val;
    retval2|=subprobp[i]->get_function_sol(lsg_val,lub_val,&tmpvec);
    if (retval2) continue;
    lower_bound+=lsg_val-fun_lb(i); 
    cand_subg_val+=lsg_val;
    cand_ub_fun_val+=lub_val;
    cand_subg+=tmpvec;  //cand_subg +=tmpvec;
  }}
  if (retval2) return retval1;
  
  //all information collected
  cand_available=1;
  return retval1;
}

// *****************************************************************************
//                            get_function_sol
// *****************************************************************************

//returns lower and upper bound on objective value of last eval_function()
//if subgp is not 0 then an eps-subgradient at y is returned,
//the hyperplane corresponding to the subgradient has value lb in y
//returns: 0 ... if the information is available 
//               (if eval_function returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int SumProblem::get_function_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!cand_available) return 1;
  lb=cand_subg_val;
  ub=cand_ub_fun_val;
  if (subgp){
    *subgp=cand_subg;
  }
  return 0;
}

// *****************************************************************************
//                              do_step
// *****************************************************************************

//store y of last eval_function, its subgradient, and objective value
//as new center of stability.

//!!! it might make more sense to use aug_subg and aug_linconst
//    for the information in the current point

int SumProblem::do_step(void)
{
  if (!cand_available) return 1;
  y=cand_y;
  subg=cand_subg;
  subg_val=cand_subg_val;
  ub_fun_val=cand_ub_fun_val;
  center_available=1;
  for (unsigned int i=0;i<subprobp.size();i++){
    subprobp[i]->do_step();
  }

  return 0;
}


// *****************************************************************************
//                              eval_model
// *****************************************************************************

//evaluate the current cutting model in $y$ 
//if evaluated by an iterative method that provides upper and lower bounds
//it may stop when the lower bound (lb) is above the nullstep_bound
//evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied

int SumProblem::eval_model(const Matrix& iny,Real eval_model_bound,Real relprec)
{
  model_available=0;
  model_subg_available=0;
  Real lower_bound=ip(iny,b);
  model_subg_val=lower_bound;
  model_ub_val=lower_bound;
  
  if (subprobp.size()==1){    //just one subproblem

    //compute value for subproblems 
    int retval1=subprobp[0]->eval_model(iny,
			       eval_model_bound-lower_bound,
			       relprec);

    //collect solution informaton
    Real lsg_val;
    Real lub_val;
    int retval2=subprobp[0]->get_model_sol(lsg_val,lub_val,0);
    if (retval2){
      return retval1;
    } 
    model_subg_val+=lsg_val;
    model_ub_val+=lub_val;
    model_available=1;
    return retval1;
  }

  //----  nonnegative sum of several subproblems

  //compute values for subproblems
  int retval1=0;    
  for (unsigned int i=0;i<subprobp.size();i++){
    retval1|=subprobp[i]->eval_model(iny,max_Real,relprec);
  }
  
  //collect solution information
  {for (unsigned int i=0;i<subprobp.size();i++){
    Real lsg_val;
    Real lub_val;
    int retval2 = subprobp[i]->get_model_sol(lsg_val,lub_val,0);
    if (retval2){
      return retval1;
    } 
    model_subg_val+=lsg_val;
    model_ub_val+=lub_val;
  }}
  
  //all information collected
  model_available=1;
  return retval1;
}

// *****************************************************************************
//                              get_model_sol
// *****************************************************************************

  //return lower and upper bound on objective value of last eval_function()
  //if subgp is not 0 then a eps-subgradient at y is returned,
  //the hyperplane corresponding to the subgradient has value lb in y
  //returns: 0 ... if the information is available 
  //               (if eval_model returned 1, the information will not
  //                satisfy the precision requirements)
  //         1 ... if the desired information is not available

int SumProblem::get_model_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!model_available) return 1;
  lb=model_subg_val;
  ub=model_ub_val;
  if (subgp){
    if (!model_subg_available){
      //collect subgradient information
      Matrix tmpvec;
      model_subg=b;
      for (unsigned int i=0;i<subprobp.size();i++){
	Real lsg_val;
	Real lub_val;
	int retval2 =subprobp[i]->get_model_sol(lsg_val,lub_val,&tmpvec);
	if (retval2){
	  return 1;
	} 
	model_subg+=tmpvec;
      }
      model_subg_available=1;
    }
    *subgp=model_subg;
  }
  return 0;
}


// *****************************************************************************
//                              get_augmodel_sol
// *****************************************************************************

int SumProblem::get_augmodel_sol(Real& out_linconst,Matrix& out_subg)
  //return last solution of (re)eval_augmodel in the following sense:
  //the linear function linconst+ip(subg,.) is a minorant of the objective
  //and the minimizer over the corresponding augmented model 
  // newy:= argmin linval+ip(subg-eta,.)+weightu/2*||.-y||^2
  // is an approximate minimizer of the entire augmented model
  //if (re)eval_augmodel returned 1, this quadratic model has minimum value lb
  //returns: 0 ... if the information is available 
  //               (if eval_augmodel  returned 1, the information will not
  //                satisfy the precision requirements)
  //         1 ... if the desired information is not available
{
  if (!aug_available) return 1;
  out_linconst=aug_linconst;
  out_subg=aug_subg;
  return 0;
}

// *****************************************************************************
//                                update_model
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int SumProblem::update_model(bool descent_step)
{
  if ((!cand_available)||(!aug_available)) return 1;
  for (unsigned int i=0;i<subprobp.size();i++){
    subprobp[i]->update_model(descent_step);
  }
  return 0;
}
  
// *****************************************************************************
//                              intersect_box
// *****************************************************************************

int SumProblem::intersect_box(Indexmatrix& inbindex,Matrix& inlb,Matrix& inub)
{
  chk_init(inbindex);
  chk_init(inlb);
  chk_init(inub);
  if ((inbindex.dim()==0)&&(inlb.dim()==0)&&(inub.dim()==0)){
    inbindex=bounds_index;
    inlb=lby;
    inub=uby;
    return 0;
  }
  chk_add(inub,uby);
  chk_add(inlb,lby);
  Indexmatrix tmpind(max(inbindex.dim()+bounds_index.dim(),uby.dim()),1);
  chk_set_init(tmpind,1);
  Integer i=0;
  Integer j=0;
  Integer nind=0;
  while ((i<bounds_index.dim())&&(j<inbindex.dim())){
    if (bounds_index(i)==inbindex(j)){
      tmpind(nind)=bounds_index(i);
      inub(tmpind(nind))=min(uby(tmpind(nind)),inub(tmpind(nind)));
      inlb(tmpind(nind))=max(lby(tmpind(nind)),inlb(tmpind(nind)));
      i++;
      j++;
    }
    else if (bounds_index(i)<=inbindex(j)){
      tmpind(nind)=bounds_index(i);
      inub(tmpind(nind))=uby(tmpind(nind));
      inlb(tmpind(nind))=lby(tmpind(nind));
      i++;
    }
    else {
      tmpind(nind)=inbindex(j);
      j++;
    }
    nind++;
  }
  while (i<bounds_index.dim()){
    tmpind(nind)=bounds_index(i);
    i++;
    inub(tmpind(nind))=uby(tmpind(nind));
    inlb(tmpind(nind))=lby(tmpind(nind));
    nind++;
  }
  while (j<inbindex.dim()){
    tmpind(nind)=inbindex(j);
    j++;
    nind++;
  }
  tmpind.reduce_length(nind);
  inbindex=tmpind;
  return 0;
}

// *****************************************************************************
//                                 lb_function
// *****************************************************************************

//returns a *quick* lower bound for the function value at y
//(eg by a previous subgradient)

Real SumProblem::lb_function(const Matrix& iny)
{
  if (aug_available){
    return aug_linconst+ip(aug_subg,iny);
  }
  if (cand_available){
    return cand_subg_val+ip(cand_subg,iny-cand_y);
  }
  if (center_available){
    return subg_val+ip(subg,iny-y);
  }
  Real lb=0;
  for (unsigned int i=0;i<subprobp.size();i++){
    lb+=subprobp[i]->lb_function(iny);
  }
  return lb;
}

// *****************************************************************************
//                                lb_model
// *****************************************************************************

//returns a *quick* lower bound for the model value at y
//(eg by one model subgradient)

Real SumProblem::lb_model(const Matrix& iny)
{
  if (aug_available){
    return aug_linconst+ip(aug_subg,iny);
  }
  if (cand_available){
    return cand_subg_val+ip(cand_subg,iny-cand_y);
  }
  Real lb=0;
  for (unsigned int i=0;i<subprobp.size();i++){
    lb+=subprobp[i]->lb_model(iny);
  }
  return lb;
}


// *****************************************************************************
//                                start_augmodel
// *****************************************************************************

//return a pointer to the variables/constraints generating the cutting model
//returns 0 on success, 1 on failure

int SumProblem::start_augmodel(QP_Block*& blockp)
{
  //--- set blocks
  //initialize blocks and variable dimensions
  aug_available=0;
  block.clear_blocks();
  aug_ind.init(Integer(subprobp.size()),1,Integer(0));
  Integer dim=0;
  for (unsigned int i=0;i<subprobp.size();i++){
    aug_ind(i)=dim;
    QP_Block* tmpblock;
    if (subprobp[i]->start_augmodel(tmpblock)){
      return 1;
    }
    dim+=tmpblock->xdim();
    block.add_block(tmpblock);
  }
  
  blockp=&block;
  return 0;
}
 
// *****************************************************************************
//                                  get_row
// *****************************************************************************

//store the coefficients corresponding to coordinate index_y of y
//for all model variables in the positions 
//row(startindex,...,startindex+model_dim-1)
//index_y==-1 is used for the cost coefficients

int SumProblem::get_row(Integer index_y,Matrix& row,Real& bi,Integer startindex) const
{
  if (index_y>=0) bi+=b(index_y);
  for (unsigned int i=0;i<subprobp.size();i++){
    if (subprobp[i]->get_row(index_y,row,bi,startindex+aug_ind(i))){
      return 1;
    }
  }
  return 0;
}
 
// *****************************************************************************
//                                make_aug_linmodel
// *****************************************************************************

//form the aggregate from the current solution in the QP_Block
//this yields a minorant  linconst+<subg,.> of the objective
//if the pointers are not nil then
//add linconst to *aug_linconst and subg to *aug_subg
//returns 0 on success, 1 on failure

int SumProblem::make_aug_linmodel(Real* in_aug_linconst,Matrix* in_aug_subg,bool* conemult_increased,Real* function_value)  
{
  aug_linconst=0.;
  aug_subg=b;
  bool new_conemultinc=false;
  bool* ncp=0;
  if (conemult_increased) ncp=&new_conemultinc;
  for (unsigned int i=0;i<subprobp.size();i++){
    if (subprobp[i]->make_aug_linmodel(&aug_linconst,&aug_subg,ncp,function_value)){
      return 1;
    }
  }
  aug_available=1;
  if (new_conemultinc) {
    center_available=0;
    recompute_center();
    *conemult_increased=new_conemultinc;
    if (function_value==0) return 2;
    *function_value=ub_fun_val;
  }

  if (in_aug_linconst) *in_aug_linconst+=aug_linconst;
  if (in_aug_subg) *in_aug_subg+=aug_subg;

  return 0;
} 
 
// *****************************************************************************
//                                adjust_multiplier
// *****************************************************************************

int SumProblem::adjust_multiplier()  
{
  clear_model();
  for (unsigned int i=0;i<subprobp.size();i++){
    if (subprobp[i]->adjust_multiplier()){
      return 1;
    }
  }
  if (recompute_center()){
    return 1;
  }
  return 0;
}
 
// *****************************************************************************
//                                change_variables
// *****************************************************************************

//change variables as described in ChangeVariableInfo
//returns 0 on success, 1 on failure
//(this is not needed in the bundle framework, routine may return 1!)

int SumProblem::change_variables(ChangeVarInfo* cvp)
{
  //--- AppendVars
 
  if (typeid(*cvp)==typeid(AppendVarsSumProblem)){
    AppendVarsSumProblem* cp=dynamic_cast<AppendVarsSumProblem*>(cvp);
    if (cp->subprobp.size()!=subprobp.size()){
      if (out) (*out)<<"*** ERROR: SumProblem::add_variable(): not sufficient information for all subproblems "<<std::endl;
      return 1;
    }
    if (cp->cost) b.concat_below(*(cp->cost));
    else b.concat_below(Matrix(cp->n_append,1,0.));
    if (cp->boundsi) bounds_index.concat_below(*(cp->boundsi));
    if (cp->lb) lby.concat_below(*(cp->lb));
    else lby.concat_below(Matrix(cp->n_append,1,CB_minus_infinity));
    if (cp->ub) uby.concat_below(*(cp->ub));
    else uby.concat_below(Matrix(cp->n_append,1,CB_plus_infinity));
    if (cp->startval) y.concat_below(*(cp->startval));
    else y.concat_below(Matrix(cp->n_append,1,0.));
    problem_modified=1;
    center_available=0;
    cand_available=0;
    model_available=0;
    model_subg_available=0;
    aug_available=0;
    //--- pass information on to subproblems
    for(unsigned int i=0;i<subprobp.size();i++){
      if (subprobp[i]->change_variables(cp->subprobp[i])) {
	if (out) (*out)<<"*** ERROR: SumProblem::add_variable(): failed for subproblem "<<i<<std::endl;
	return 1;
      }
    }
    return 0;
  }

  //--- ReassignVars

  if (typeid(*cvp)==typeid(ReassignVarsSumProblem)){
    ReassignVarsSumProblem* cp=dynamic_cast<ReassignVarsSumProblem*>(cvp);
    if (cp->subprobp.size()!=subprobp.size()){
      if (out) (*out)<<"*** ERROR: SumProblem::change_variables(): not sufficient information for all subproblems "<<std::endl;
      return 1;
    }
    if (center_available){
      //check whether a deleted variable ore one with multiple copies is nonzero
      Indexmatrix tmpind(y.dim(),1,Integer(0));
      for(Integer i=0;i<cp->assign_ind.dim();i++){
	if (((tmpind(cp->assign_ind(i))+=1)>1)&&(y(cp->assign_ind(i))!=0.)){
	  center_available=0;
	  break;
	}
      }
      if (center_available){
	Indexmatrix delind=tmpind.find_number(0);
	for(Integer i=0;i<delind.dim();i++){
	  if (y(delind(i))!=0.){
	    center_available=0;
	    break;
	  }
	}
      }
      if (center_available){ 
	subg=subg(cp->assign_ind);
      }
    }
      
    b=b(cp->assign_ind);
    lby=lby(cp->assign_ind);
    uby=uby(cp->assign_ind);
    if (bounds_index.dim()>0){
      bounds_index.init(0,0,Integer(0));
      for(Integer i=0;i<lby.dim();i++){
	if ((lby(i)>CB_minus_infinity)||(uby(i)<CB_plus_infinity)){
	  bounds_index.concat_below(i);
	}
      }
    }
    y=y(cp->assign_ind);
    
    problem_modified=1;
    cand_available=0;
    aug_available=0;
    model_available=0;
    model_subg_available=0;

    //--- pass information on to subproblems
    for(unsigned int i=0;i<subprobp.size();i++){
      if (subprobp[i]->change_variables(cp->subprobp[i])) {
	if (out) (*out)<<"*** ERROR: SumProblem::del_variable(): failed for subproblem "<<i<<std::endl;
	return 1;
      }
    }
    
    return 0;
  }

  //--- DeleteVars

  if (typeid(*cvp)==typeid(DeleteVarsSumProblem)){
    DeleteVarsSumProblem* cp=dynamic_cast<DeleteVarsSumProblem*>(cvp);
    if (cp->subprobp.size()!=subprobp.size()){
      if (out) (*out)<<"*** ERROR: SumProblem::del_variable(): not sufficient information for all subproblems "<<std::endl;
      return 1;
    }
    if (cp->del_index.dim()==0) return 0;
    if (bounds_index.dim()>0){
      Integer dcnt=0;
      Integer dind=cp->del_index(dcnt);
      Integer bcnt=0;
      Integer bind=bounds_index(bcnt);
      Integer dim=b.dim();
      Indexmatrix tmpind(0,0,Integer(0)); //will be used to delete entries in b_i
      while((dind<dim)||(bind<dim)){
	if (bind<=dind){
	  if (bind==dind) tmpind.concat_below(bcnt);
	  else bounds_index(bcnt)-=dcnt;
	  bcnt++;
	  if (bcnt<bounds_index.dim()) bind=bounds_index(bcnt);
	  else bind=dim;
	}
	else {
	  dcnt++;
	  if (dcnt<cp->del_index.dim()) dind=cp->del_index(dcnt);
	  else dind=dim;
	}
      }
      if (tmpind.dim()>0) bounds_index.delete_rows(tmpind); 
    }
    if (center_available){
      int non_zeros=0;
      for (Integer i=0;i<cp->del_index.dim();i++){
	if (y(cp->del_index(i))!=0.) {
	  non_zeros=1;
	  break;
	}
      }
      if (non_zeros) center_available=0;
      else {
	subg.delete_rows(cp->del_index);
      }
    }
      
    b.delete_rows(cp->del_index);
    lby.delete_rows(cp->del_index);
    uby.delete_rows(cp->del_index);
    y.delete_rows(cp->del_index);
    
    problem_modified=1;
    cand_available=0;
    aug_available=0;
    model_available=0;
    model_subg_available=0;

    //--- pass information on to subproblems
    for(unsigned int i=0;i<subprobp.size();i++){
      if (subprobp[i]->change_variables(cp->subprobp[i])) {
	if (out) (*out)<<"*** ERROR: SumProblem::del_variable(): failed for subproblem "<<i<<std::endl;
	return 1;
      }
    }
    
    return 0;
  }
       
  return 1;
}

  

// *****************************************************************************
//                                recompute_center
// *****************************************************************************

//after modifications of the problem the center information may have
//to be recomputed partially or compeletely. Do the appropriate.
//returns 0 on success, 1 on failure 

int SumProblem::recompute_center()
{
  if (center_available) return 0;
  Matrix tmpy,tmpsubg,tmpylb,tmpyub;
  Indexmatrix tmpind;
  Real lb,ub;
  subg_val=ip(b,y);
  ub_fun_val=subg_val;
  subg=b;
  for(unsigned int i=0;i<subprobp.size();i++){
    if (subprobp[i]->recompute_center()) {
      if (out) (*out)<<"*** ERROR: SumProblem::recompute_center(): failed for subproblem "<<i<<std::endl;
      return 1;
    }

    int init=0;
    if (subprobp[i]->init_center(init,lb,ub,tmpy,tmpsubg,tmpind,tmpylb,tmpyub)){
      if (out) (*out)<<"*** ERROR: SumProblem::recompute_center(): init_center failed for subproblem "<<i<<std::endl;
      return 1;
    }
    subg_val+=lb;
    ub_fun_val+=ub;
    subg+=tmpsubg;
  }
  center_available=1;

  return 0;
}

// *****************************************************************************
//                                set_new_center
// *****************************************************************************

//set a new center point. If the input is NULL then a default
//starting point is constructed.
//returns 0 on success, 1 on failure 

int SumProblem::set_new_center(const Matrix* yp) 
{ 
  center_available=0;
  problem_modified=1;
  if (yp==0) {
    y.init(lby);
    for(Integer i=0;i<y.dim();i++){
      if (y(i)>uby(i)) return 1;
      y(i)=max(y(i),min(0.,uby(i)));
    }
  }
  else {
    y.init(*yp);
    for(Integer i=0;i<y.dim();i++){
      if ((y(i)<lby(i))||(y(i)>uby(i))) return 1;
    }
  }

  int status=eval_function(y,CB_plus_infinity,1e-6);
  status |= do_step();
  for (unsigned int i=0;i<subprobp.size();i++){
    status|= subprobp[i]->update_model(true);
  }
  return status;
}

// *****************************************************************************
//                                set_lower_bound
// *****************************************************************************

int SumProblem::set_lower_bound(Integer i,Real lb)
{
  if ((i<0)||(i>=lby.dim())) {
    if (out) (*out)<<"ERROR: SumProblem::set_lower_bound(): index "<<i<<" out of range"<<std::endl;
    return 1;
  }
  lb=max(lb,CB_minus_infinity);
  if (lb>uby(i)) {
    if (out) (*out)<<"ERROR: SumProblem::set_lower_bound(): new lower bound "<<lb<<" for index "<<i<<" exceeds current upper bound "<<uby(i)<<std::endl;
    return 1;
  }
  if (lb==lby(i)) return 0;

  problem_modified=1;   
  if ((center_available)&&(lb>y(i))) center_available=0;
  cand_available=0;
  model_available=0;
  model_subg_available=0;
  aug_available=0;

  if (lb==CB_minus_infinity){
    lby(i)=lb;
    if (uby(i)==CB_plus_infinity){
      //delete entry in bounds_index.dim();
      Integer lbj=0; 
      Integer ubj=bounds_index.dim()-1;
      while (lbj<ubj){
	Integer j=(lbj+ubj)/2;
	if (i>bounds_index(j)) { lbj=j+1; continue; }
        ubj=j;
      }
      if (bounds_index(lbj)==i){
	bounds_index.delete_rows(Indexmatrix(1,1,lbj));
      }
      else {
	if (out) (*out)<<"ERROR: SumProblem::set_lower_bound(): could not find index "<<i<<" in bounds_index()"<<std::endl;
	return 1;
      }
    }
    return 0;
  }

  if ((lby(i)!=CB_minus_infinity)||(uby(i)!=CB_plus_infinity)){
    lby(i)=lb;
    return 0;
  }

  lby(i)=lb;
  
  Integer lbj=0; 
  Integer ubj=bounds_index.dim()-1;
  while (lbj<ubj){
    Integer j=(lbj+ubj)/2;
    if (i>bounds_index(j)) { lbj=j+1; continue; }
    ubj=j;
  }
  if ((lbj<bounds_index.dim())&&(bounds_index(lbj)==i)){
    if (out) (*out)<<"WARNING: SumProblem::set_lower_bound(): index "<<i<<" already in bounds_index even though it should not"<<std::endl;
    return 0;
  }
  bounds_index.insert_row(lbj,Indexmatrix(1,1,i));
  return 0;
}

// *****************************************************************************
//                                set_upper_bound
// *****************************************************************************

int SumProblem::set_upper_bound(Integer i,Real ub)
{
  if ((i<0)||(i>=uby.dim())) {
    if (out) (*out)<<"ERROR: SumProblem::set_upper_bound(): index "<<i<<" out of range"<<std::endl;
    return 1;
  }
  ub=min(ub,CB_plus_infinity);
  if (ub<lby(i)) {
    if (out) (*out)<<"ERROR: SumProblem::set_upper_bound(): new upper bound "<<ub<<" for index "<<i<<" exceeds current lower bound "<<lby(i)<<std::endl;
    return 1;
  }

  if (ub==uby(i)) return 0;

  problem_modified=1;
  if ((center_available)&&(ub<y(i))) center_available=0;   
  cand_available=0;
  model_available=0;
  model_subg_available=0;
  aug_available=0;

  if (ub==CB_plus_infinity){
    uby(i)=ub;
    if (lby(i)==CB_minus_infinity){
      //delete entry in bounds_index.dim();
      Integer lbj=0; 
      Integer ubj=bounds_index.dim()-1;
      while (lbj<ubj){
	Integer j=(lbj+ubj)/2;
	if (i>bounds_index(j)) { lbj=j+1; continue; }
        ubj=j;
      }
      if (bounds_index(lbj)==i){
	bounds_index.delete_rows(Indexmatrix(1,1,lbj));
      }
      else {
	if (out) (*out)<<"ERROR: SumProblem::set_upper_bound(): could not find index "<<i<<" in bounds_index()"<<std::endl;
	return 1;
      }
    }
    return 0;
  }

  if ((lby(i)!=CB_minus_infinity)||(uby(i)!=CB_plus_infinity)){
    uby(i)=ub;
    return 0;
  }

  uby(i)=ub;
  
  Integer lbj=0; 
  Integer ubj=bounds_index.dim()-1;
  while (lbj<ubj){
    Integer j=(lbj+ubj)/2;
    if (i>bounds_index(j)) { lbj=j+1; continue; }
    ubj=j;
  }
  if ((lbj<bounds_index.dim())&&(bounds_index(lbj)==i)){
    if (out) (*out)<<"WARNING: SumProblem::set_upper_bound(): index "<<i<<" already in bounds_index even though it should not"<<std::endl;
    return 0;
  }
  bounds_index.insert_row(lbj,Indexmatrix(1,1,i));
  return 0;
}

// *****************************************************************************
//                                clear_model
// *****************************************************************************

void SumProblem::clear_model()  
{
  problem_modified=1;   
  center_available=0;
  cand_available=0;
  model_available=0;
  model_subg_available=0;
  aug_available=0;  
}


}

