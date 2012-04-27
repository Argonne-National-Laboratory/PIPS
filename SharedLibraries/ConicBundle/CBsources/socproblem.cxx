/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/socproblem.cxx

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


#include <fstream>
#include "mymath.hxx"
#include "socproblem.hxx"
 

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                            SocProblem
// *****************************************************************************

SocProblem::SocProblem(BaseSOCOracle* oracle_function)
{
  oracle=oracle_function;
  oracle->get_trace_constraint(soc_multiplier,constant_trace);

  problem_modified=0;

  center_available=0;
  subg_val=0.;        
  ub_fun_val=0.;


  cand_available=0;
  cand_subg_val=0.;
  cand_ub_fun_val=0.;

  ret_code=0;

  model_available=0;
  model_subg_available=0;
  model_subg_val=0.;
  model_ub_val=0.;
  model_changed=0;

  model_ret_code=0;
 
  maxkeepvecs=50; 

  aug_available=0;
  aug_linconst=0.;

  last_pos=0;
   
  solverp=0;

  out=0;
  print_level=0;
}

// *****************************************************************************
//                            ~SocProblem
// *****************************************************************************

SocProblem::~SocProblem()
{  
  delete solverp; solverp=0;
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

int SocProblem::init_center(
	       int& init, Real& lb, Real& ub, Matrix& o_y, Matrix& o_subg,
               Indexmatrix& o_bounds_index, Matrix& o_lby, Matrix& o_uby)
{
  if (!center_available) return 1;
  if ((init==0)||(problem_modified)||(o_y.dim()!=y.dim())){
    lb=subg_val;
    ub=ub_fun_val;
    o_y=y;
    o_subg=subg;
    o_bounds_index=bounds_index;
    if (lby.dim()==0) o_lby.init(y.dim(),1,CB_minus_infinity);
    else o_lby=lby;
    if (uby.dim()==0) o_uby.init(y.dim(),1,CB_plus_infinity);
    else o_uby=uby;
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

int SocProblem::eval_function(const Matrix& iny,
			      Real /* nullstep_bound */,Real /* relprec */)
{
  cand_available=0;
  cand_y=iny;
  Real ip_b_iny=ip(oracle->rhs(),cand_y);   //linear term

  Real trace_val;
  SOCtrace con_tr;
  oracle->get_trace_constraint(trace_val,con_tr);
  assert(trace_val>=0.);
  if ((con_tr!=constant_trace)||
      ((con_tr!=SOCtrace_unbounded)&&(trace_val!=soc_multiplier))||
      ((con_tr==SOCtrace_unbounded)&&(trace_val>soc_multiplier))){
    if (out) {
      (*out)<<"**** WARNING: SocProblem::eval_function(...): trace condition changed without reinitialization of model and function values"<<std::endl;
    }
    soc_multiplier=trace_val;
    constant_trace=con_tr;
  }
   
  ret_code=oracle->c_At(cand_y,cand_c);
  if (ret_code) {
    if (out) {
      (*out)<<"**** WARNING: SocProblem::eval_function(...): oracle->c_At returned "<<ret_code<<std::endl;
    }
  }

  if (cand_c.coldim()!=1){
    if (out) {
      (*out)<<"**** WARNING: SocProblem::eval_function(...): oracle->c_At returned not a vector"<<std::endl;
    }

    if (ret_code) return ret_code;
    return 1;
  }
  Real c0=cand_c(0);
  cand_c(0)=0.;
  Real n2=norm2(cand_c);
  if (n2<eps_Real){
    cand_vec.rand(cand_c.dim(),1);
    cand_vec(0)=0.;
    Real nn2=norm2(cand_vec);
    cand_vec/=nn2;
  }
  else {
    cand_vec.init(cand_c,1./n2);
  }
  cand_vec(0)=1.;
  cand_c(0)=c0;

  //--------- compute candidate subgradient and objective values   
  cand_lmax=c0+n2;
  {
    if ((out)&&(print_level>0)){
      out->precision(6);
      (*out)<<" soclmax="<<cand_lmax;
      (*out)<<std::endl;
    }
  }
  if ((constant_trace!=SOCtrace_fixed)&&(cand_lmax<0.)) cand_lmax=0.;
  Real objval=soc_multiplier*cand_lmax;
  cand_subg_val=ip_b_iny+objval;
  cand_ub_fun_val=cand_subg_val;
  //compute cand_subg=b-soc_multiplier*A*cand_vec;
  if ((cand_lmax!=0.)||(constant_trace==SOCtrace_fixed)){
    oracle->ip_A(cand_vec,cand_subg);
    xbpeya(cand_subg,oracle->rhs(),1.,-soc_multiplier);
  }
  else {
    cand_subg=oracle->rhs();
  }
  cand_available=1;

  return ret_code;
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

int SocProblem::get_function_sol(Real& lb, Real& ub, Matrix* subgp)
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

int SocProblem::do_step(void)
{
  if (!cand_available) return 1;
  y=cand_y;
  lmax=cand_lmax;
  subg=cand_subg;
  subg_val=cand_subg_val;
  ub_fun_val=cand_ub_fun_val;
  center_vec=cand_vec;
  center_available=1;
  problem_modified=0;

  return 0;
}


// *****************************************************************************
//                              eval_model
// *****************************************************************************

//evaluate the current cutting model in $y$ 
//if evaluated by an iterative method that provides upper and lower bounds
//it may stop when the lower bound (lb) is above the nullstep_bound
//evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied

int SocProblem::eval_model(const Matrix& iny,Real /* eval_model_bound */,Real /* relprec */)
{
  model_available=0;
  model_subg_available=0;
  model_changed=0;

  if (bundlevecs.dim()==0){
    if (out){
      (*out)<<"**** ERROR: SocProblem::eval_model(...): called without defining a cutting surface model first."<<std::endl;
    }
    return 1;
  }

  ret_code=oracle->c_At(iny,tmpvec);
  if (ret_code) {
    if (out) {
      (*out)<<"**** WARNING: SocProblem::eval_model(...): oracle->c_At returned "<<ret_code<<std::endl;
    }
  }

  if (tmpvec.coldim()!=1){
    if (out) {
      (*out)<<"**** WARNING: SocProblem::eval_modle(...): oracle->c_At returned not a vector"<<std::endl;
    }
    if (ret_code) return ret_code;
    return 1;
  }

  Matrix coords;
  genmult(bundlevecs,tmpvec,coords,1.,0.,1);
  Real c0=coords(0);
  coords(0)=0.;
  Real n2=norm2(coords);
  Real objval=soc_multiplier*(c0+n2);
  if ((constant_trace!=SOCtrace_fixed)&&(objval<0)){
    objval=0.;
    model_vec.init(bundlevecs.coldim(),1,0.);
  }
  else {
    if (n2<eps_Real){
      model_vec.rand(coords.dim(),1);
      model_vec(0)=0.;
      Real nn2=norm2(model_vec);
      model_vec/=nn2;
    }
    else {
      model_vec.init(coords,1./n2);
    }
    model_vec(0)=1;
  }
  Real ip_b_iny=ip(oracle->rhs(),iny);
  
  //--- determine final model
  model_subg_val=ip_b_iny+objval;
  model_ub_val=model_subg_val;

  model_available=1;
  model_subg_available=0;

  return 0;
}

// *****************************************************************************
//                              get_model_sol
// *****************************************************************************

  //return lower and upper bound on objective value of last eval_model()
  //if subgp is not 0 then a eps-subgradient at y is returned,
  //the hyperplane corresponding to the subgradient has value lb in y
  //returns: 0 ... if the information is available 
  //               (if eval_model returned 1, the information will not
  //                satisfy the precision requirements)
  //         1 ... if the desired information is not available

int SocProblem::get_model_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!model_available) return 1;
  lb=model_subg_val;
  ub=model_ub_val;
  if (subgp){
    if (!model_subg_available){
      //collect subgradient information
      if (model_changed) return 1;
      oracle->ip_A(bundlevecs*model_vec,model_subg);
      xbpeya(model_subg,oracle->rhs(),1.,-soc_multiplier);	
      model_subg_available=1;
    }
    *subgp=model_subg;
  }
  return 0;
}


// *****************************************************************************
//                              get_augmodel_sol
// *****************************************************************************

int SocProblem::get_augmodel_sol(Real& out_linconst,Matrix& out_subg)
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

int SocProblem::update_model(bool /*descent_step*/)
{
  //--- if there are no vectors in the bundle intialize to [e_1, candvec] 
  if (bundlevecs.coldim()==0){    
    if (cand_vec.coldim()==0) return 1;
    Integer n=cand_vec.dim();
    model_changed=1;
    bundlevecs.init(n,2,0.);
    bundlevecs(0,0)=1.;
    mat_xey(n-1,bundlevecs.get_store()+n+1,cand_vec.get_store()+1);
    vec_store.newsize(n,max(maxkeepvecs-2,1));
    vec_store.init(cand_vec);
    last_pos=0;
    return 0;
  }

  //---  check that new candidate vector is consistend with existing vectors
  if ((!cand_available)||(primalvec.dim()!=cand_vec.dim())) return 1;

  //--- do nothing if the candidate vector would be the zero vector;
  if ((cand_lmax==0.)&&(constant_trace!=SOCtrace_fixed)) {
    if ((out)&&(print_level>=2)){
      out->precision(8);
      (*out)<<" soc bundle_update: cand_lmax="<<cand_lmax<<", no update needed"<<std::endl;
    }    
    return 0;
  }
  model_changed=1;

  Integer n=bundlevecs.rowdim();
  
  //--- if there are too many bundle vectors extract a few of them
  tmpmat.newsize(n-1,max(2,min(bundlevecs.coldim()+1,maxkeepvecs-1))); 
  chk_set_init(tmpmat,1);
  mat_xey(n-1,tmpmat.get_store(),primalvec.get_store()+1);
  mat_xey(n-1,tmpmat.get_store()+n-1,cand_vec.get_store()+1);

  int strategy=2;

  switch (strategy){

    case 0: 
    case 1:
      tmpind.init(Range(1,bundlevecs.coldim()-1));
      //for strategy 1 sort the vectors by the current cost_function
      if ((strategy==1)&&(bundlevecs.coldim()>=maxkeepvecs)){
	genmult(bundlevecs,cand_c,tmpvec,-1.,0.,1);
	tmpvec(0)=2.*fabs(max(tmpvec));
	sortindex(tmpvec,tmpind);
      }
      {for(Integer i=2;i<tmpmat.coldim();i++){
	Integer ind=tmpind(i-2);
	if (ind==0) continue;
	mat_xey(n-1,tmpmat.get_store()+i*(n-1),bundlevecs.get_store()+ind*n+1);
      }}
      break;
      
    case 2:
      //append the candidates in vec_store
      {for(Integer i=2;i<tmpmat.coldim();i++){
	Integer ind= (last_pos-(i-2))%vec_store.coldim();
	if (ind<0) ind+=vec_store.coldim();
	mat_xey(n-1,tmpmat.get_store()+i*(n-1),vec_store.get_store()+ind*n+1);
      }}
      last_pos++;
      if (last_pos>=max(maxkeepvecs-2,1)){
	last_pos=0;	
      }
      if (last_pos>=vec_store.coldim()){
	vec_store.concat_right(cand_vec);
      }
      else {
	mat_xey(n,vec_store.get_store()+last_pos*n,cand_vec.get_store());
      }
  }
	
  Indexmatrix piv;
  Integer r=tmpmat.QR_factor(piv);
  tmpvec.init(tmpmat.rowdim(),r,0.);
  for(Integer i=0;i<r;i++) tmpvec(i,i)=1.;
  tmpmat.Q_times(tmpvec,r);

  if ((out)&&(print_level>=2)){
    out->precision(8);
    (*out)<<" soc bundle_update: n_old="<<bundlevecs.coldim()-1<<" rank="<<r<<" maxkeep="<<maxkeepvecs<<std::endl;
  }

  bundlevecs.init(n,min(tmpvec.coldim()+1,maxkeepvecs+1),0.);
  bundlevecs(0,0)=1.;
  {for(Integer i=1;i<bundlevecs.coldim();i++){
    mat_xey(n-1,bundlevecs.get_store()+i*n+1,tmpvec.get_store()+(i-1)*(n-1));
  }}

 
  return 0;
}  

  
// *****************************************************************************
//                              intersect_box
// *****************************************************************************

int SocProblem::intersect_box(Indexmatrix& inbindex,Matrix& inlb,Matrix& inub)
{
  chk_init(inbindex);
  chk_init(inlb);
  chk_init(inub);
  Integer dim=oracle->rhs().dim();
  if ((inbindex.dim()==0)&&(inlb.dim()==0)&&(inub.dim()==0)){
    inbindex=bounds_index;
    if (lby.dim()==0) inlb.init(dim,1,CB_minus_infinity);
    else inlb=lby;
    if (uby.dim()==0) inub.init(dim,1,CB_plus_infinity);
    else inub=uby;
    return 0;
  }
  assert(inub.dim()==dim);
  assert(inlb.dim()==dim);
  tmpind.newsize(min(inbindex.dim()+bounds_index.dim(),dim),1);
  chk_set_init(tmpind,1);
  Integer i=0;
  Integer j=0;
  Integer nind=0;
  if ((uby.dim()>0)&&(lby.dim()>0)) {
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
  }
  else if (lby.dim()==0){
    while ((i<bounds_index.dim())&&(j<inbindex.dim())){
      if (bounds_index(i)==inbindex(j)){
	tmpind(nind)=bounds_index(i);
	inub(tmpind(nind))=min(uby(tmpind(nind)),inub(tmpind(nind)));
	i++;
	j++;
      }
      else if (bounds_index(i)<=inbindex(j)){
	tmpind(nind)=bounds_index(i);
	inub(tmpind(nind))=uby(tmpind(nind));
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
      nind++;
    }
    while (j<inbindex.dim()){
      tmpind(nind)=inbindex(j);
      j++;
      nind++;
    }
  }
  else if (uby.dim()==0){
    while ((i<bounds_index.dim())&&(j<inbindex.dim())){
      if (bounds_index(i)==inbindex(j)){
	tmpind(nind)=bounds_index(i);
	inlb(tmpind(nind))=max(lby(tmpind(nind)),inlb(tmpind(nind)));
	i++;
	j++;
      }
      else if (bounds_index(i)<=inbindex(j)){
	tmpind(nind)=bounds_index(i);
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
      inlb(tmpind(nind))=lby(tmpind(nind));
      nind++;
    }
    while (j<inbindex.dim()){
      tmpind(nind)=inbindex(j);
      j++;
      nind++;
    }
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

Real SocProblem::lb_function(const Matrix& iny)
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
  return CB_minus_infinity;
}

// *****************************************************************************
//                                lb_model
// *****************************************************************************

//returns a *quick* lower bound for the model value at y
//(eg by one model subgradient)

//here it is not possible to resort to calling subprobp[i]->lb_model
//since here the model is the task of *this

Real SocProblem::lb_model(const Matrix& iny)
{
  if (aug_available){
    return aug_linconst+ip(aug_subg,iny);
  }
  if (cand_available){
    return cand_subg_val+ip(cand_subg,iny-cand_y);
  }
  return CB_minus_infinity;
}


// *****************************************************************************
//                                start_augmodel
// *****************************************************************************

//return a pointer to the variables/constraints generating the cutting model
//returns 0 on success, 1 on failure

int SocProblem::start_augmodel(QP_Block*& blockp)
{
  //--- if block is not yet available, create one
  if (bundlevecs.coldim()==0){
    if (center_vec.coldim()==0) return 1;
    Integer n=center_vec.dim();
    bundlevecs.init(n,2,0.);
    bundlevecs(0,0)=1.;
    mat_xey(n-1,bundlevecs.get_store()+n+1,center_vec.get_store()+1);
  }
  aug_available=0;
  block.init_block(0,Indexmatrix(1,1,bundlevecs.coldim()),Indexmatrix(0,0,Integer(0)),soc_multiplier,constant_trace!=SOCtrace_fixed);
  
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

int SocProblem::get_row(Integer index_y,Matrix& row,Real& bi,Integer startindex) const
{
  if (index_y==-1){ //cost coefficients
    oracle->project_c(bundlevecs,tmpvec);
    mat_xey(tmpvec.dim(),
	    row.get_store()+startindex,
	    tmpvec.get_store());
  }
  else { //coefficient matrix
    bi+=(oracle->rhs())(index_y);
    oracle->project(index_y,bundlevecs,tmpvec);
    mat_xemy(tmpvec.dim(),
	    row.get_store()+startindex,
	    tmpvec.get_store());
  }  
  return 0;
}
 
// *****************************************************************************
//                               make_aug_linmodel
// *****************************************************************************

//form the aggregate from the current solution in the QP_Block
//this yields a minorant  linconst+<subg,.> of the objective
//if the pointers are not nil then
//add linconst to *aug_linconst and subg to *aug_subg
//returns 0 on success, 1 on failure

int SocProblem::make_aug_linmodel(Real* in_aug_linconst,Matrix* in_aug_subg,
				  bool* conemult_increased,Real* function_value)  
{
  block.get_socx(0,tmpvec);
  genmult(bundlevecs,tmpvec,primalvec);

  oracle->ip_cA(primalvec,aug_linconst,aug_subg);
  xbpeya(aug_subg,oracle->rhs(),1.,-1.);
  
  aug_available=1;

  if (in_aug_linconst) *in_aug_linconst+=aug_linconst;
  if (in_aug_subg) *in_aug_subg+=aug_subg;

#ifdef TRACE_OUTPUT
  if (out) (*out)<<"\n SocProblem: trace="<<tmpvec(0)<<" mult="<<soc_multiplier<<std::endl;
#endif

  //--- check whether it is necessary to increase the cone multiplier
  if ((conemult_increased!=0)&&(constant_trace==SOCtrace_unbounded)&&(tmpvec(0)>.95*soc_multiplier)){
    soc_multiplier=max(1.,soc_multiplier*2.);
#ifdef TRACE_OUTPUT
    if (out) (*out)<<"              multiplier increased to "<<soc_multiplier<<std::endl;
#endif
    *conemult_increased=true;
    block.adjust_trace(soc_multiplier);
    if (lmax>0.){
      Real ip_b_iny=ip(oracle->rhs(),y);
      subg_val=ip_b_iny+lmax*soc_multiplier;
      ub_fun_val=subg_val;
      oracle->ip_A(center_vec,subg);
      xbpeya(subg,oracle->rhs(),1.,-soc_multiplier);
      if (function_value==0) return 2;
      *function_value=ub_fun_val;      
    }
  }

  return 0;
} 
 
// *****************************************************************************
//                                adjust_multiplier
// *****************************************************************************

int SocProblem::adjust_multiplier()  
{
  if (constant_trace!=SOCtrace_unbounded) return 0;
  if (!center_available)  return 1;
  if (block.get_socx(0,tmpvec)) return 1;
#ifdef TRACE_OUTPUT
  if (out) (*out)<<"\n SocProblem: trace="<<tmpvec(0)<<" mult="<<soc_multiplier;
#endif
  soc_multiplier=min(soc_multiplier,max(1.,2.*tmpvec(0)));
  oracle->adjust_multiplier(soc_multiplier);
#ifdef TRACE_OUTPUT
  if (out) (*out)<<" newmult="<<soc_multiplier<<std::endl;
#endif
  block.adjust_trace(soc_multiplier);
  if (lmax>0.){
    Real ip_b_iny=ip(oracle->rhs(),y);
    subg_val=ip_b_iny+lmax*soc_multiplier;
    ub_fun_val=subg_val;
    oracle->ip_A(center_vec,subg);
    xbpeya(subg,oracle->rhs(),1.,-soc_multiplier);
  }
  return 0;
}
 
// *****************************************************************************
//                                change_variables
// *****************************************************************************

//change variables as described in ChangeVariableInfo
//returns 0 on success, 1 on failure
//(this is not needed in the bundle framework, routine may return 1!)

int SocProblem::change_variables(ChangeVarInfo* cvp)
{
  if (typeid(*cvp)==typeid(AppendVars)){
    AppendVars* cp=dynamic_cast<AppendVars*>(cvp);
    Integer n_append=cp->n_append;
    Integer dim=oracle->rhs().dim()-n_append;
    if (cp->boundsi) 
      bounds_index.concat_below(*(cp->boundsi));
    if (cp->lb) {
      if (lby.coldim()==0) 
	lby.init(dim,1,CB_minus_infinity);
      lby.concat_below(*(cp->lb));
    }
    else {
      if (lby.coldim()!=0) 
	lby.concat_below(Matrix(n_append,1,CB_minus_infinity));
    }
    if (cp->ub) {
      if (uby.coldim()==0) 
	uby.init(dim,1,CB_plus_infinity);
      uby.concat_below(*(cp->ub));
    }
    else {
      if (uby.coldim()!=0) 
	uby.concat_below(Matrix(n_append,1,CB_plus_infinity));
    }
    if (cp->startval) {
      if (center_available){
	for(Integer i=0;i<n_append;i++){
	  if ((*cp->startval)(i)!=0.){
	    center_available=0;
	    break;
	  }
	}
      } 
      y.concat_below(*(cp->startval));
    }
    else {
      if (y.coldim()==0)
	y.init(dim,1,0.);
      y.concat_below(Matrix(n_append,1,0.));
    }
    dim=oracle->rhs().dim();
        
    if (center_available){
      //make sure that the subgradient in the center is contained in the model
      if (primalvec.dim()) {
	oracle->ip_cA(primalvec,subg_val,subg);
	xbpeya(subg,oracle->rhs(),1.,-1.);
	subg_val+=ip(subg,y);
      }
      else {
	center_available=0;
      }
    }
    
    problem_modified=1;
    
    //candidate (never mind about the candidate)
    cand_available=0;
    
    //model (never mind about the model value)
    model_available=0;
    model_subg_available=0;
    
    //augmodel (never mind about the augmented model value and subgradient)
    aug_available=0;
    
    return 0;
  }
  
   //--- ReassignVars

  if (typeid(*cvp)==typeid(ReassignVars)) {
    ReassignVars* cp=dynamic_cast<ReassignVars*>(cvp);
    if (center_available){
      //check whether a deleted variable ore one with multiple copies is nonzero
      tmpind.init(y.dim(),1,Integer(0));
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
    }
      
    y=y(cp->assign_ind);
    if (lby.dim()!=0) lby=lby(cp->assign_ind);
    if (uby.dim()!=0) uby=uby(cp->assign_ind);
    if (bounds_index.dim()>0){
      bounds_index.init(0,0,Integer(0));
      for(Integer i=0;i<y.dim();i++){
	if (((lby.dim()!=0)&&(lby(i)>CB_minus_infinity))||
	    ((uby.dim()!=0)&&(uby(i)<CB_plus_infinity))){
	  bounds_index.concat_below(i);
	}
      }
    }
    
    if (center_available) {
      //make sure that the subgradient in the center is contained in the model
      if (primalvec.dim()) {
	oracle->ip_cA(primalvec,subg_val,subg);
	xbpeya(subg,oracle->rhs(),1.,-1.);
	subg_val+=ip(subg,y);
      }
      else {
	center_available=0;
      }
    }

    problem_modified=1;
    cand_available=0;
    aug_available=0;
    model_available=0;
    model_subg_available=0;

    return 0;
  }

  //--- DeleteVars
 
  if (typeid(*cvp)==typeid(DeleteVars)){
    DeleteVars* cp=dynamic_cast<DeleteVars*>(cvp);
    if (cp->del_index.dim()==0) return 0;
    if (bounds_index.dim()>0){
      Integer dim=oracle->rhs().dim();
      Integer dcnt=0;
      Integer dind=cp->del_index(dcnt);
      Integer bcnt=0;
      Integer bind=bounds_index(bcnt);
      tmpind.init(0,0,Integer(0)); //will be used to delete entries in b_i
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
    }
    
    if (lby.dim()!=0) lby.delete_rows(cp->del_index);
    if (uby.dim()!=0) uby.delete_rows(cp->del_index);
    y.delete_rows(cp->del_index);
    
    if (center_available) {
      //make sure that the subgradient in the center is contained in the model
      if (primalvec.dim()) {
	oracle->ip_cA(primalvec,subg_val,subg);
	xbpeya(subg,oracle->rhs(),1.,-1.);
	subg_val+=ip(subg,y);
      }
      else {
	center_available=0;
      }
    }
    
    problem_modified=1;
    cand_available=0;
    model_available=0;
    model_subg_available=0;
    aug_available=0;

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

int SocProblem::recompute_center()
{
  if (!center_available) {
    int status=eval_function(y,max_Real,1e-5);
    if (status) return status;
    int pm=problem_modified;
    status = do_step();
    problem_modified=pm;
    if (status) return status;
  }
  
  center_available=1;
  return 0;  
}

// *****************************************************************************
//                                get_approximate_primal
// *****************************************************************************

int SocProblem::get_approximate_primal(PrimalData& pd) const 
  {
    if (!aug_available) return -1;   
    PrimalDVector* p=dynamic_cast<PrimalDVector*>(&pd);
    if (p!=0) {
      assign(*p,primalvec);
      return 0;
    }
    PrimalMatrix* pm=dynamic_cast<PrimalMatrix*>(&pd);
    if (pm!=0) {
      *pm=primalvec;
      return 0;
    }
    return 1;
  }


// *****************************************************************************
//                                get_center_primal
// *****************************************************************************

int SocProblem::get_center_primal(PrimalData& pd) const 
  {
    PrimalDVector* p=dynamic_cast<PrimalDVector*>(&pd);
    if (p!=0) {
      assign(*p,center_vec);
      return 0;
    }
    PrimalMatrix* pm=dynamic_cast<PrimalMatrix*>(&pd);
    if (pm!=0) {
      *pm=center_vec;
      return 0;
    }
    return 1;
  }


// *****************************************************************************
//                                get_center_primal
// *****************************************************************************

int SocProblem::get_candidate_primal(PrimalData& pd) const 
  {
    PrimalDVector* p=dynamic_cast<PrimalDVector*>(&pd);
    if (p!=0) {
      assign(*p,cand_vec);
      return 0;
    }
    PrimalMatrix* pm=dynamic_cast<PrimalMatrix*>(&pd);
    if (pm!=0) {
      *pm=cand_vec;
      return 0;
    }
    return 1;
  }


// *****************************************************************************
//                                clear_model
// *****************************************************************************

void SocProblem::clear_model()  
{
  problem_modified=1;   
  center_available=0;
  cand_available=0;
  model_available=0;
  model_subg_available=0;
  aug_available=0;
  
  bundlevecs.init(0,0,0.);
  primalvec.init(0,0,0.);
  center_vec.init(0,0,0.);
  
  oracle->get_trace_constraint(soc_multiplier,constant_trace);
}

}

