/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/lmaxproblem.cxx

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
#include "lmaxproblem.hxx"
#include "coeffmat.hxx" 
#include "MatSDPfun.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                            LmaxProblem
// *****************************************************************************

LmaxProblem::LmaxProblem(BaseSDPOracle* oracle_function)
{
  oracle=oracle_function;
  oracle->get_trace_constraint(lmax_multiplier,trace_stat);

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
  model_eig=0.;
  model_ind=-1;

  model_ret_code=0;

  update_rule=0;

  /*
  maxkeepvecs=25; 
  minkeepvecs=5; 
  maxaddvecs=5;  
  maxaggrcols=1; 
  aggregtol=0.01;   
  */
 
  maxkeepvecs=15; 
  minkeepvecs=3; 
  maxaddvecs=7;  
  maxaggrcols=3; 
  aggregtol=0.01;   

  n_new_vecs=0;

  rankdefcnt=0;

  activedim=0;
  keepsize=0;
  curve_delta=0.;

  aug_subg_in_model=0;
  aug_available=0;
  aug_xdim=0;
  aug_linconst=0.;
   
  solverp=0;

  out=0;
  print_level=0;
}

// *****************************************************************************
//                            ~LmaxProblem
// *****************************************************************************

LmaxProblem::~LmaxProblem()
{  
  for(unsigned int i=0;i<primal.size();i++){
    delete primal[i];
  }
  primal.clear();
  delete solverp; solverp=0;
}

// *****************************************************************************
//                                init_subgrad
// *****************************************************************************
 
//Compute the subgradient corresponding to the matrix P*Diag(d)*P'/sum(d)
//and store it in a new column. P has orthonormal columns and d is 
//a nonengative vector of appropriate size.
//If "store_symsubg" is set, then the full matrix itself is stored, as well.
//If "store_sparsesubg" is set, then the matrix is stored on a sparse support
 
void LmaxProblem::init_subgrad(const Matrix& P,const Matrix& d)
{
  Real sumd=sum(d);
  if (sumd<eps_Real) return; //nothing to aggregate
  tmpvec.init(d,1./sumd);
  tmpvec.sqrt();
  tmpmat=P;
  tmpmat.scale_cols(tmpvec);
  
  SDPPrimal* p=oracle->init_primal(tmpmat);
  if ((p!=0)&&(Integer(primal.size())==aggrsubgrads.coldim())){
    primal.push_back(p);
  }
  else {
    delete p;
    for(unsigned int i=0;i<primal.size();i++){
      delete primal[i];
    }
    primal.clear();
  }

  Real cost;
  oracle->gramip(tmpmat,cost,tmpvec);
  if (aggrsubgrads.coldim()==0){
    aggrsubgrads=tmpvec;
    aggrcoeff.init(1,1,sum(d));
    subgCvalues.init(1,1,cost);
  }
  else {
    aggrsubgrads.concat_right(tmpvec);
    aggrcoeff.concat_below(sum(d));
    subgCvalues.concat_below(cost);
  }
}

// ****************************************************************************
//                                update_subgrad
// ****************************************************************************
 
//Compute the subgradient corresponding to the matrix P*Diag(d)*P'/sum(d)
//and add it to the stored subgradient.  P has orthonormal columns and
//d is a nonengative vector of appropriate size.
//If "store_symsubg" is set, then the full supbgradient matrix is updated,
//as well.
 
void LmaxProblem::update_subgrad(const Matrix& P,const Matrix& d,Integer aggr_index)
{
  Real sumd=sum(d);
  if (sumd<eps_Real) return; //nothing to aggregate
  Real alpha=sumd/(sumd+aggrcoeff(aggr_index));
  tmpvec.init(d,alpha/sumd);
  tmpvec.sqrt();
  tmpmat=P;
  tmpmat.scale_cols(tmpvec);
 
  if (Integer(primal.size())==aggrsubgrads.coldim()){
    primal[aggr_index]->aggregate_Gram_matrix(1.-alpha,1.,tmpmat);
  }
  
  Real cost;
  oracle->gramip(tmpmat,cost,tmpvec);
  mat_xbpeya(aggrsubgrads.rowdim(),
	     aggrsubgrads.get_store()+aggrsubgrads.rowdim()*aggr_index,
	     tmpvec.get_store(),1.,1-alpha);
  aggrcoeff(aggr_index)+=sumd;
  subgCvalues(aggr_index)*=(1-alpha);
  subgCvalues(aggr_index)+=cost;
}

 
// ****************************************************************************
//                                delete_subgrads
// ****************************************************************************
 
 
void LmaxProblem::delete_subgrads(const Indexmatrix& delind)
{
  if ((Integer(primal.size())==aggrsubgrads.coldim())&&(primal.size()>0)){
    Indexmatrix sind;
    sortindex(delind,sind);
    Integer cnt=delind(sind(0));
    delete primal[cnt];
    primal[cnt]=0;
    Integer delcnt=0;
    Integer nextdel;
    if (++delcnt<delind.dim()){
      nextdel=delind(sind(delcnt));
    }
    else {
      nextdel=Integer(primal.size());
    }
      for(Integer i=cnt+1;i<Integer(primal.size());i++){
	if (i==nextdel){
	  delete primal[i];
	  primal[i]=0;
	  if (++delcnt<delind.dim()){
	    nextdel=delind(sind(delcnt));
	  }
	  else {
	    nextdel=Integer(primal.size());
	  }
	  continue;
	} 
	primal[cnt]=primal[i];
	cnt++;
      }
      primal.resize(cnt);
  }
  subgCvalues.delete_rows(delind);
  aggrsubgrads.delete_cols(delind);
  aggrcoeff.delete_rows(delind);
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

int LmaxProblem::init_center(
	       int& init, Real& lb, Real& ub, Matrix& o_y, Matrix& o_subg,
               Indexmatrix& o_bounds_index, Matrix& o_lby, Matrix& o_uby)
{
  if (!center_available) return 1;
  if ((init==0)||(problem_modified)||(o_y.dim()!=y.dim())||(!aug_subg_in_model)){
    ub=ub_fun_val;
    lb=subg_val;
    o_y=y;
    o_subg=subg;
    o_bounds_index=bounds_index;
    if (lby.dim()==0) o_lby.init(y.dim(),1,CB_minus_infinity);
    else o_lby=lby;
    if (uby.dim()==0) o_uby.init(y.dim(),1,CB_plus_infinity);
    else o_uby=uby;

    curve_delta=0.;

    if (bundlevecs.coldim()==0){
      if (eigvec.coldim()==0) return 1;
      Integer nc=min(eigvec.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
      bundlevecs=eigvec.cols(Indexmatrix(Range(0,nc)));
    }
	
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

int LmaxProblem::eval_function(const Matrix& iny,
			      Real nullstep_bound,Real relprec)
{
  cand_available=0;
  cand_y=iny;
  Real ip_b_iny=ip(oracle->rhs(),cand_y);   //linear term

  Real trace_val;
  SDPtrace con_tr;
  oracle->get_trace_constraint(trace_val,con_tr);
  assert(trace_val>=0.);
  if ((con_tr!=trace_stat)||
      ((con_tr!=SDPtrace_unbounded)&&(trace_val!=lmax_multiplier))||
      ((con_tr==SDPtrace_unbounded)&&(trace_val>lmax_multiplier))){
    if (out) {
      (*out)<<"**** WARNING: LmaxProblem::eval_function(...): trace condition changed without reinitialization of model and function values"<<std::endl;
    }
    lmax_multiplier=trace_val;
    trace_stat=con_tr;
  }
   
  Real eigval_bound=(nullstep_bound-ip_b_iny)/lmax_multiplier;
  Real eigval_relprec=relprec/lmax_multiplier*(fabs(nullstep_bound)+1.)/(fabs(eigval_bound)+1.);

  //if (topvecs.rowdim()!=bundlevecs.rowdim()){
    ret_code=oracle->evaluate(cand_y,bundlevecs, eigval_relprec,eigval_bound,
			      cand_vecs,cand_eigs);
    /*
  }
  else {
    tmpmat=topvecs;
    tmpmat.concat_right(skippedvecs);
    ret_code=oracle->evaluate(cand_y,tmpmat, eigval_relprec,eigval_bound,
			      cand_vecs,cand_eigs);
			      }
    */
    
  if (ret_code) {
    if (out) {
      (*out)<<"**** WARNING: LmaxProblem::eval_function(...): evaluate returned "<<ret_code<<std::endl;
    }
  }
  if ((cand_vecs.coldim()==0)||(cand_vecs.coldim()!=cand_eigs.rowdim())){
    if (out) {
      (*out)<<"**** WARNING: LmaxProblem::eval_function(...): evaluate returned no eigenvectors or differing number of eigenvalues"<<std::endl;
    }
    if (ret_code) return ret_code;
    return 1;
  }

  //--------- compute candidate subgradient and objective values   
  Integer maxind;
  Real lmax=max(cand_eigs,&maxind);
  cand_ublmax=lmax+eigval_relprec*(fabs(lmax)+1.);
  if ((trace_stat!=SDPtrace_fixed)&&(lmax<0.)) lmax=0.;
  if ((trace_stat!=SDPtrace_fixed)&&(cand_ublmax<0.)) cand_ublmax=0.;
  cand_subg_val=ip_b_iny+lmax*lmax_multiplier;
  cand_ub_fun_val=ip_b_iny+cand_ublmax*lmax_multiplier;


  //compute cand_subg=b-lmax_multiplier*A(vv');
  if ((lmax!=0.)||(trace_stat==SDPtrace_fixed)){
    tmpmat=cand_vecs.col(maxind);
    oracle->gramip_opA(tmpmat,cand_subg);
    xbpeya(cand_subg,oracle->rhs(),1.,-lmax_multiplier);
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

int LmaxProblem::get_function_sol(Real& lb, Real& ub, Matrix* subgp)
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

int LmaxProblem::do_step(void)
{
  if (!cand_available) return 1;
  y=cand_y;
  subg=cand_subg;
  subg_val=cand_subg_val;
  ub_fun_val=cand_ub_fun_val;

  ublmax=cand_ublmax;
  Integer maxind;
  eigval.init(1,1,max(cand_eigs,&maxind));
  eigvec=cand_vecs.col(maxind);
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

int LmaxProblem::eval_model(const Matrix& iny,Real eval_model_bound,Real relprec)
{
  model_available=0;
  model_subg_available=0;
  model_changed=0;

  if (bundlevecs.dim()==0){
    if (out){
      (*out)<<"**** ERROR: LmaxProblem::eval_model(...): called without defining a cutting surface model first."<<std::endl;
    }
    return 1;
  }

  Real ip_b_iny=ip(oracle->rhs(),iny);   //linear term

  //--- determine model values for bundlevectors
  model_eig=CB_minus_infinity;
  model_vec.init(0,0,0.);         //reset 

  Real eigval_bound=(eval_model_bound-ip_b_iny)/lmax_multiplier;
  Real eigval_relprec=relprec/lmax_multiplier*(fabs(eval_model_bound)+1.)/(fabs(eigval_bound)+1.);

  model_ret_code = oracle->evaluate_projection(iny,bundlevecs,eigval_relprec,
				       eigval_bound,tmpmat,tmpvec);
  if (model_ret_code) {
      if (out) {
	(*out)<<"**** WARNING: LmaxProblem::eval_model(...): evaluate_projection returned "<<model_ret_code<<std::endl;
      }
  }
  if ((tmpmat.coldim()==0)||(tmpmat.coldim()!=tmpvec.rowdim())){
    if (out) {
      (*out)<<"**** WARNING: LmaxProblem::eval_model(...): evaluate_projection returned no eigenvectors or differing number of eigenvalues"<<std::endl;
    }
    if (model_ret_code) return model_ret_code;
    return 1;
  }
  
  Integer maxind;
  model_eig=max(tmpvec,&maxind);

  //--- determine model values for aggregates
  Real tmpmax=CB_minus_infinity;
  model_ind=-1;
  if (subgCvalues.dim()>0){
    tmpvec=subgCvalues;
    genmult(aggrsubgrads,iny,tmpvec,-1.,1.,1);
    tmpmax=max(tmpvec,&model_ind);
  }
  
  //--- determine final model
  if ((trace_stat==SDPtrace_fixed)||(model_eig>0.)||(tmpmax>0.)){
    if (model_eig>tmpmax){
      if ((out)&&(print_level>0)){
	(*out)<<" proj,";
      }
      genmult(bundlevecs,tmpmat.col(maxind),model_vec);  //for lazy model_subg eval
      model_ind=-1;
      model_subg_val=lmax_multiplier*model_eig+ip_b_iny;  
      model_ub_val=model_subg_val+lmax_multiplier*eigval_relprec*(fabs(model_eig+1.));
    }
    else {
      if ((out)&&(print_level>0)){
	(*out)<<" aggr,";
      }
      model_eig=tmpmax;
      model_subg_val=lmax_multiplier*model_eig+ip_b_iny;  
      model_ub_val=model_subg_val;      
    }
  }
  else {
    model_eig=0.;
    model_vec.init(bundlevecs.rowdim(),1,0.);
    model_ind=-1;
    model_subg_val=ip_b_iny;  
    model_ub_val=model_subg_val;
  }

  if ((out)&&(print_level>0)){
    (*out)<<" model_eig="<<model_eig;
  }

  model_available=1;
  model_subg_available=0;

  return 0;
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

int LmaxProblem::get_model_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!model_available) return 1;
  lb=model_subg_val;
  ub=model_ub_val;
  if (subgp){
    if (!model_subg_available){
      //collect subgradient information
      if (model_ind>=0){
	if (model_changed) return 1;
	//compute model_subg=b-lmax_multplier*aggrsubgrads.col(model_ind);
	model_subg=oracle->rhs();
	mat_xpeya(model_subg.dim(),model_subg.get_store(),
		  aggrsubgrads.get_store()+model_subg.dim()*model_ind,
		  -lmax_multiplier);
      }
      else {
	//compute model_subg=b-lmax_multiplier*A(vv'); for v=model_vec
	oracle->gramip_opA(model_vec,model_subg);
	xbpeya(model_subg,oracle->rhs(),1.,-lmax_multiplier);	
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

int LmaxProblem::get_augmodel_sol(Real& out_linconst,Matrix& out_subg)
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
//        -1 ... if some major internal error occured

int LmaxProblem::update_model(bool descent_step)
{
  switch(update_rule){
    case 0: return update_model_default(descent_step);
    case 1: return update_model_no_aggregate(descent_step); 
    case 2: return update_model_top_spectrum(descent_step); 
    case 3: return update_model_top_spectrum2(descent_step); 
    case 4: return update_model_top_spectrum3(descent_step); 
    case 5: return update_model_top_spectrum4(descent_step); 
    case 6: return update_model_top_spectrum5(descent_step); 
    case 7: return update_model_top_spectrum6(descent_step); 
    case 8: return update_model_top_spectrum7(descent_step); 
    case 9: return update_model_top_spectrum8(descent_step); 
    case 10: return update_model_new_spectrum(descent_step); 
    case 11: return update_model_tapia_scale(descent_step); 
    case 12: return update_model_tapia_scale1(descent_step); 
    case 13: return update_model_tapia_scale2(descent_step); 
    case 14: return update_model_tapia_scale3(descent_step); 
    case 15: return update_model_tapia_scale4(descent_step); 
    case 16: return update_model_tapia_scale5(descent_step); 
    case 17: return update_model_tapia_scale6(descent_step); 
    case 18: return update_model_tapia_scale7(descent_step); 
    case 19: return update_model_tapia_scale8(descent_step); 
    case 20: return update_model_tapia_scale9(descent_step); 
    case 21: return update_model_tapia_scale10(descent_step); 
    case 22: return update_model_tapia_scale11(descent_step); 
    case 23: return update_model_tapia_scale12(descent_step); 
    default: return update_model_default(descent_step);
  }
  return -1;
}

// *****************************************************************************
//                                update_model_default
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_default(bool)
{
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  if (!cand_available) return 1;

  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- aggregate the aggrsubgrads if coefficients are small
  if (aggrcoeff.dim()>1){
    Indexmatrix ind;
    sortindex(aggrcoeff,ind);
    Integer naggr=0;
    Real aggregval=aggregtol*max(aggrcoeff);
    while((naggr<ind.dim())&&(aggrcoeff(ind(naggr))<aggregval)) naggr++;
    if ((naggr>0)&&(aggrsubgrads.coldim()>1)){
      Integer aggrindex;
      if (naggr<ind.dim()) aggrindex=ind(naggr);
      else aggrindex=0;
      Integer n=aggrsubgrads.rowdim();
      for(Integer i=0;i<naggr;i++){
	Integer indi=ind(i);
	if (indi==aggrindex) continue;
	Real b=aggrcoeff(indi);
	if (b<eps_Real) continue; //nothing to aggregate
	Real a=aggrcoeff(aggrindex);
	Real s=a+b;
	a/=s;
	b/=s;
	subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
	mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		   aggrsubgrads.get_store()+indi*n,b,a);
	if (Integer(primal.size())==aggrsubgrads.coldim()){
	  primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
	}
	aggrcoeff(aggrindex)+=aggrcoeff(indi);  
      }
      if (naggr<aggrsubgrads.coldim()){
	ind.reduce_length(naggr);
	delete_subgrads(ind);
      }
      else { //only column 0 is not deleted
	ind.init(Range(1,aggrsubgrads.coldim()-1));
	delete_subgrads(ind);
      }
    }//endif ((naggr>0)&&(aggrsubgrads.coldim()>1)
  } //endif (aggrcoeff.dim()>1)

  //--- aggregate the bundle vectors
  //primaleigs is sorted non increasingly
  Integer maxadd=maxaddvecs; //min(maxaddvecs,max(primalvecs.rowdim()/4,Integer(1)));
  if (primaleigs.dim()>0){
    Integer maxkeep=min(maxkeepvecs,primalvecs.rowdim()-maxadd);
    Real aggregval=aggregtol*primaleigs(0);
    if  ((out)&&(print_level>2)){
      (*out)<<" ";
      (*out).setf(std::ios::showpoint);
    }
    if ((primaleigs.dim()>maxkeep)||(primaleigs(primaleigs.dim()-1)<aggregval)){
      
      //---- fill ind with column indices that are to be aggregated
      Integer i=min(primaleigs.dim(),maxkeep)-1;
      Integer lbi=max(minkeepvecs,Integer(0));
      while((i>=lbi)&&(primaleigs(i)<aggregval)) --i;
      Indexmatrix ind(Range(i+1,primaleigs.dim()-1));
      
      //---- aggregate
      if  ((out)&&(print_level>2)){(*out)<<" aggrdim="<<ind.dim();out->precision(4);}
      if (ind.dim()>0){
	//find colummn to aggregate to
	Integer aggrindex;   //index of column to aggregate to
	int new_aggrcol=0;   //set to true if a new column should be opened
	if (aggrsubgrads.coldim()<maxaggrcols){
	  aggrindex=aggrsubgrads.coldim();
	  new_aggrcol=1;
	}
	else { //choose the one with minimal contribution
	  min(aggrcoeff,&aggrindex);
	}
	
	//aggregate to the column aggrindex of aggrsubgrads
	if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
	  if (!new_aggrcol){
	    update_subgrad(primalvecs.cols(ind),primaleigs(ind),
			   aggrindex);
	  }
	  else {
	    init_subgrad(primalvecs.cols(ind),primaleigs(ind));
	  }
	  primalvecs.delete_cols(ind);
	  primaleigs.delete_rows(ind);
	}
	else { //aggregate all columns
	  if (!new_aggrcol){
	    update_subgrad(primalvecs,primaleigs,aggrindex);
	  }
	  else {
	    init_subgrad(primalvecs,primaleigs);
	  }
	  primalvecs.init(0,0,0.);
	  primaleigs.init(0,0,0.);
	}
      }
    }
  }

  bundlevecs=primalvecs;
  

  //--- enlarge bundle with new information
  int nr_new=min(maxadd,cand_vecs.coldim());
  n_new_vecs=nr_new;
 
  if  ((out)&&(print_level>2)){
    if (bundlevecs.coldim()>0){
      tmpvec=cand_vecs.col(0);
      tmpvec.transpose();
      tmpvec*=bundlevecs;
      tmpvec.transpose();
      (*out)<<" res="<<norm2(bundlevecs*tmpvec-cand_vecs.col(0))<<std::endl;
    }
    else (*out)<<" res="<<1.<<std::endl;
    (*out)<<" neweigs=";
    for(Integer i=0;i<nr_new;i++){
      out->precision(6);(*out)<<" "<<cand_eigs(i);
    }
    (*out)<<std::endl;
    if (primaleigs.dim()>0){
      (*out)<<" primaleigs=";
      for(Integer i=0;i<primaleigs.dim();i++){
	out->precision(6);(*out)<<" "<<primaleigs(i);
      }
      (*out)<<std::endl;
    }
    if (aggrcoeff.dim()>0){
      (*out)<<" aggrcoeff=";
      for(Integer i=0;i<aggrcoeff.dim();i++){
	out->precision(6);(*out)<<" "<<aggrcoeff(i);
      }
      (*out)<<std::endl;
    }
  }
  
  if (bundlevecs.coldim()==0){ //only new vectors
    bundlevecs.init(cand_vecs.cols(Range(0,nr_new-1)));
    return 0;
  }
  
  //old and new vectors
  bundlevecs.concat_right(cand_vecs.cols(Range(0,nr_new-1)));
  // --- orthogonalize bundlevecs
//     {
//         Integer i,k;
//         Integer n=bundlevecs.rowdim();
//         for(k=0;k<bundlevecs.coldim();k++){
//             Real *xk=bundlevecs.get_store()+k*n;
//             Real t,t1;
//             do{
//                 //compute projection of vector k on vector i and subtract this
//                 t=0.;
//                 for(i=0;i<k;i++){
//                     Real *xi=bundlevecs.get_store()+i*n;
//                     t1=mat_ip(n,xi,xk);
//                     t+=t1*t1;
//                     mat_xpeya(n,xk,xi,-t1);
//                 }
//                 t1=mat_ip(n,xk,xk);
//                 t+=t1;
//             }while(t1<t/100.);
//             t1=sqrt(t1);
//             if (t1>1e-20) {
//                 //vector is not in the span, normalize it
//                 mat_xmultea(n,xk,1./t1);
//             }
//             else {
//                 //vector is in the span of previous vectors, delete it
//                 bundlevecs.delete_cols(Indexmatrix(1,1,k));
//                 k--;  //correct the value of k
//             }
//         }

//         primaleigs.init(bundlevecs.coldim(),1,1.);
//     }
 
 
  //if the old bundlevecs may be changed this is more stable
  tmpmat=bundlevecs;
  Indexmatrix piv;
  Integer r=tmpmat.QR_factor(piv);
  if (r<bundlevecs.coldim()) {
    if (out) {
      (*out)<<"\n*** WARNING: bundle update linearly dependent:";
      (*out)<<" n_columns="<<bundlevecs.coldim()<<" rank="<<r<<std::endl;
    }
    rankdefcnt++;
  }
  
  bundlevecs.init(tmpmat.rowdim(),r,0.);
  for(Integer i=0;i<r;i++) bundlevecs(i,i)=1.;
  tmpmat.Q_times(bundlevecs,r);

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_default(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  

  
// *****************************************************************************
//                                update_model_no_aggregate
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_no_aggregate(bool)
{
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  if (!cand_available) return 1;

  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- aggregate the aggrsubgrads if coefficients are small
  if (aggrcoeff.dim()>1){
    Indexmatrix ind;
    sortindex(aggrcoeff,ind);
    Integer naggr=0;
    Real aggregval=aggregtol*max(aggrcoeff);
    while((naggr<ind.dim())&&(aggrcoeff(ind(naggr))<aggregval)) naggr++;
    if ((naggr>0)&&(aggrsubgrads.coldim()>1)){
      Integer aggrindex;
      if (naggr<ind.dim()) aggrindex=ind(naggr);
      else aggrindex=0;
      Integer n=aggrsubgrads.rowdim();
      for(Integer i=0;i<naggr;i++){
	Integer indi=ind(i);
	if (indi==aggrindex) continue;
	Real b=aggrcoeff(indi);
	if (b<eps_Real) continue; //nothing to aggregate
	Real a=aggrcoeff(aggrindex);
	Real s=a+b;
	a/=s;
	b/=s;
	subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
	mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		   aggrsubgrads.get_store()+indi*n,b,a);
	if (Integer(primal.size())==aggrsubgrads.coldim()){
	  primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
	}
	aggrcoeff(aggrindex)+=aggrcoeff(indi);  
      }
      if (naggr<aggrsubgrads.coldim()){
	ind.reduce_length(naggr);
	delete_subgrads(ind);
      }
      else { //only column 0 is not deleted
	ind.init(Range(1,aggrsubgrads.coldim()-1));
	delete_subgrads(ind);
      }
    }//endif ((naggr>0)&&(aggrsubgrads.coldim()>1)
  } //endif (aggrcoeff.dim()>1)

  //--- aggregate the bundle vectors
  //primaleigs is sorted non increasingly
  Integer maxadd=maxaddvecs; //min(maxaddvecs,max(primalvecs.rowdim()/4,Integer(1)));
  if (primaleigs.dim()>0){
    Integer maxkeep=min(maxkeepvecs,primalvecs.rowdim()-maxadd);
    Real aggregval=aggregtol*primaleigs(0);
    if  ((out)&&(print_level>2)){
      (*out)<<" ";
      (*out).setf(std::ios::showpoint);
    }
    if ((primaleigs.dim()>maxkeep)||(primaleigs(primaleigs.dim()-1)<aggregval)){
      
      //---- fill ind with column indices that are to be aggregated
      Integer i=min(primaleigs.dim(),maxkeep)-1;
      Integer lbi=max(minkeepvecs,Integer(0));
      while((i>=lbi)&&(primaleigs(i)<aggregval)) --i;
      Indexmatrix ind(Range(i+1,primaleigs.dim()-1));
      
      //---- aggregate
      if  ((out)&&(print_level>2)){(*out)<<" aggrdim="<<ind.dim();out->precision(4);}
      if (ind.dim()>0){
	//find colummn to aggregate to
	Integer aggrindex;   //index of column to aggregate to
	int new_aggrcol=0;   //set to true if a new column should be opened
	if (aggrsubgrads.coldim()<maxaggrcols){
	  aggrindex=aggrsubgrads.coldim();
	  new_aggrcol=1;
	}
	else { //choose the one with minimal contribution
	  min(aggrcoeff,&aggrindex);
	}
	
	//aggregate to the column aggrindex of aggrsubgrads
	if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
	  if (!new_aggrcol){
	    update_subgrad(primalvecs.cols(ind),primaleigs(ind),
			   aggrindex);
	  }
	  else {
	    init_subgrad(primalvecs.cols(ind),primaleigs(ind));
	  }
	  primalvecs.delete_cols(ind);
	  primaleigs.delete_rows(ind);
	}
	else { //aggregate all columns
	  if (!new_aggrcol){
	    update_subgrad(primalvecs,primaleigs,aggrindex);
	  }
	  else {
	    init_subgrad(primalvecs,primaleigs);
	  }
	  primalvecs.init(0,0,0.);
	  primaleigs.init(0,0,0.);
	}
      }
    }
  }

  bundlevecs=primalvecs;
  

  //--- enlarge bundle with new information
  int nr_new=min(maxadd,cand_vecs.coldim());
  n_new_vecs=nr_new;
 
  if  ((out)&&(print_level>2)){
    if (bundlevecs.coldim()>0){
      tmpvec=cand_vecs.col(0);
      tmpvec.transpose();
      tmpvec*=bundlevecs;
      tmpvec.transpose();
      (*out)<<" res="<<norm2(bundlevecs*tmpvec-cand_vecs.col(0))<<std::endl;
    }
    else (*out)<<" res="<<1.<<std::endl;
    (*out)<<" neweigs=";
    for(Integer i=0;i<nr_new;i++){
      out->precision(6);(*out)<<" "<<cand_eigs(i);
    }
    (*out)<<std::endl;
    if (primaleigs.dim()>0){
      (*out)<<" primaleigs=";
      for(Integer i=0;i<primaleigs.dim();i++){
	out->precision(6);(*out)<<" "<<primaleigs(i);
      }
      (*out)<<std::endl;
    }
    if (aggrcoeff.dim()>0){
      (*out)<<" aggrcoeff=";
      for(Integer i=0;i<aggrcoeff.dim();i++){
	out->precision(6);(*out)<<" "<<aggrcoeff(i);
      }
      (*out)<<std::endl;
    }
  }
  
  if (bundlevecs.coldim()==0){ //only new vectors
    bundlevecs.init(cand_vecs.cols(Range(0,nr_new-1)));
    return 0;
  }
  
  //old and new vectors
  bundlevecs.concat_right(cand_vecs.cols(Range(0,nr_new-1)));
  // --- orthogonalize bundlevecs
//     {
//         Integer i,k;
//         Integer n=bundlevecs.rowdim();
//         for(k=0;k<bundlevecs.coldim();k++){
//             Real *xk=bundlevecs.get_store()+k*n;
//             Real t,t1;
//             do{
//                 //compute projection of vector k on vector i and subtract this
//                 t=0.;
//                 for(i=0;i<k;i++){
//                     Real *xi=bundlevecs.get_store()+i*n;
//                     t1=mat_ip(n,xi,xk);
//                     t+=t1*t1;
//                     mat_xpeya(n,xk,xi,-t1);
//                 }
//                 t1=mat_ip(n,xk,xk);
//                 t+=t1;
//             }while(t1<t/100.);
//             t1=sqrt(t1);
//             if (t1>1e-20) {
//                 //vector is not in the span, normalize it
//                 mat_xmultea(n,xk,1./t1);
//             }
//             else {
//                 //vector is in the span of previous vectors, delete it
//                 bundlevecs.delete_cols(Indexmatrix(1,1,k));
//                 k--;  //correct the value of k
//             }
//         }

//         primaleigs.init(bundlevecs.coldim(),1,1.);
//     }
 
 
  //if the old bundlevecs may be changed this is more stable
  tmpmat=bundlevecs;
  Indexmatrix piv;
  Integer r=tmpmat.QR_factor(piv);
  if (r<bundlevecs.coldim()) {
    if (out) {
      (*out)<<"\n*** WARNING: bundle update linearly dependent:";
      (*out)<<" n_columns="<<bundlevecs.coldim()<<" rank="<<r<<std::endl;
    }
    rankdefcnt++;
  }
  
  bundlevecs.init(tmpmat.rowdim(),r,0.);
  for(Integer i=0;i<r;i++) bundlevecs(i,i)=1.;
  tmpmat.Q_times(bundlevecs,r);

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_no_aggregate(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}
  
  
// *****************************************************************************
//                                update_model_top_spectrum
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum(bool descent_step)
{
  activedim=0;  
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  if (!cand_available) return 1;

  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- compute the next number of top vectors and the number to keep
  Integer keepsize=0;
  Integer nextsize=0;
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real pmax=primaleigs(0);
    Real dmax=primalZval(0);
    Real lguessmu=min(0.1,0.1*pmax*dmax);
    Real maxRitz=S(0,0);
    Ritzbound=max(1e-4,cand_eigs(0)-maxRitz);
    Real aggregval=min(aggregtol,.1);
    keepsize=1;
    while(keepsize<primaleigs.dim()){
      int primal_dropoff= (primaleigs(keepsize)<pmax*lguessmu);
      int primaldual_ratio=((maxRitz-S(keepsize,keepsize)>2.*Ritzbound)&&(primaleigs(keepsize)*dmax<aggregval*primalZval(keepsize)*pmax));
      int Ritz_gap=(maxRitz-S(keepsize,keepsize)>5.*Ritzbound);
      if (primal_dropoff||primaldual_ratio||Ritz_gap){
	if ((out)&&(print_level>2)){
	  (*out)<<"lmaxproblem::update_model_top_spectrum(...):\n gap located at index "<<keepsize;
	  (*out)<<": pdropoff="<<primal_dropoff;
	  (*out)<<" pdratio="<<primaldual_ratio;
	  (*out)<<" Rgap="<<Ritz_gap;
	}
	break;
      }
      keepsize++;
    }
    activedim=keepsize;
    if (!descent_step)
      keepsize=min(keepsize+3,primaleigs.dim());
    if ((out)&&(print_level>2))
      (*out)<<" lguessmu="<<lguessmu<<" -> keepsize="<<keepsize<<std::endl;
    nextsize=keepsize+min(cand_vecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  if (!descent_step){
    //--- compute an orthonormal basis of bundle and cand_vecs
    //first collect those that have to be included to ensure convergence
    genmult(cand_vecs.col(0),primalvecs,tmpvec,1.,0.,1);
    if ((out)&&(print_level>2))
      (*out)<<" new vec orthog="<<norm2(tmpvec)<<" full orthog="<<tmpvec<<std::endl;
    tmpmat=primalvecs;
    tmpmat.concat_right(cand_vecs.col(0)); 
    //form an orthonormal basis of these
    Indexmatrix piv;
    Integer r1=tmpmat.QR_factor(piv);
    tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
    piv.reduce_length(r1);
    if (r1==keepsize) {
      topvecs=primalvecs;
    }
    else {
      topvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
      for(Integer i=0;i<tmpmat.coldim();i++) topvecs(i,i)=1.;
      tmpmat.Q_times(topvecs,r1);
    }
    //topvecs now holds an orthonormal basis of the essential subspace
    
    //orthogonalize the remaining vectors against this basis
    bundlevecs.concat_right(skippedvecs);
    bundlevecs.concat_right(cand_vecs);
    Integer r2=tmpmat.QR_concat_right(bundlevecs,piv,r1);
    skippedvecs.init(bundlevecs.rowdim(),r2-r1,0.);
    if (r2-r1>0){
      for(Integer i=0;i<r2-r1;i++) skippedvecs(r1+i,i)=1.;
      tmpmat.Q_times(skippedvecs,r2);
    }

    //compute the Ritz_values of topvecs (tmpmat is now free again)
    oracle->evaluate_projection(cand_y,topvecs,1e-10,1e20,tmpmat,Ritz_values);
    bundlevecs=topvecs;
    if ((out)&&(print_level>2))
      (*out)<<" required bundle Ritz_values="<<transpose(Ritz_values);
    //now bundlevecs contains the vectors required for convergence, we enlarge it next
    
    if (skippedvecs.coldim()>0) {
      //compute the eigenvalues/vectors of the matrix projected onto the space spanned by skippedvecs
      oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,tmpvec);
      genmult(skippedvecs,tmpmat,topvecs);
      if ((out)&&(print_level>2))
	(*out)<<"remaining Ritz_values="<<transpose(tmpvec);
      Ritz_values.concat_below(tmpvec);
      
      //fill bundlevecs up to nextsize columns, put some of the rest into skippedvecs
      while ((nextsize>bundlevecs.coldim()+3)
	     &&(Ritz_values(0)-Ritz_values(nextsize-1)>10.*Ritzbound))
	nextsize--;
      if ((out)&&(print_level>2))
	(*out)<<" shortened nextsize="<<nextsize<<std::endl;
      bundlevecs.concat_right(topvecs.cols(Range(0,nextsize-bundlevecs.coldim()-1)));
      skippedvecs=topvecs.cols(Range(nextsize-bundlevecs.coldim(),min(nextsize-bundlevecs.coldim()+minkeepvecs,topvecs.coldim())-1));
      Ritz_values.reduce_length(bundlevecs.coldim()+skippedvecs.coldim());
    }
    
    topvecs=bundlevecs;
  }
  else {
    //---- in the case of a descent step use the primal vectors to set up
    //     the quadratic term but set the bundlevecs to the best subspace
    //compute an orthonormal basis of bundle and cand_vecs
    bundlevecs.concat_right(cand_vecs);
    bundlevecs.concat_right(skippedvecs);
    Indexmatrix piv;
    Integer r=bundlevecs.QR_factor(piv);
    skippedvecs.init(bundlevecs.rowdim(),r,0.);
    for(Integer i=0;i<r;i++) skippedvecs(i,i)=1.;
    bundlevecs.Q_times(skippedvecs,r);
  
    //compute the eigenvalues/vectors of the matrix projected onto this space
    oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
    genmult(skippedvecs,tmpmat,topvecs);
    if ((out)&&(print_level>2))
      (*out)<<" Ritz_values="<<transpose(Ritz_values);
  
    //find the number of vectors to keep for the quadratic term 
    Integer skippedsize=min(minkeepvecs,topvecs.coldim()-keepsize);
    if (skippedsize>3){
      //in "nmat" the factor was 1000. instead of 100.
      Real skip_Ritzbound=min(10.,100.*(Ritz_values(0)-Ritz_values(keepsize)));
      if ((out)&&(print_level>2))
        (*out)<<" mingap="<<Ritz_values(0)-Ritz_values(keepsize)<<" skip_Ritzbound="<<skip_Ritzbound;
      while ((skippedsize>3)&&
	     (Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)>skip_Ritzbound))
	skippedsize--;
    }
    if ((out)&&(print_level>2))
      (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
    
    //check whether the top most vectors are almost the same as primalvecs
    bundlevecs=topvecs.cols(Range(0,keepsize-1));
    genmult(bundlevecs,primalvecs,tmpmat,1.,0.,1);
    if ((out)&&(print_level>2))
      (*out)<<"keepsize="<<keepsize<<" tmpmatnorm2="<<ip(tmpmat,tmpmat)<<" normvec="<<sumrows(tmpmat%tmpmat);
    if ((keepsize-ip(tmpmat,tmpmat)<1e-2)&&(keepsize>3)){
      if ((out)&&(print_level>2))
        (*out)<<"activedim="<<activedim<<std::endl;
    }
    if (keepsize-ip(tmpmat,tmpmat)>1e-8){
      //they are too different, attach some cand_vecs to primalvecs 
      Integer addsize=1;
      while ((addsize<maxaddvecs)&&(addsize<cand_eigs.dim())&&
	     (Ritz_values(0)-cand_eigs(addsize)<10.*Ritzbound))
	addsize++;
      if ((out)&&(print_level>2))
        (*out)<<"enlarging bundle by addsize="<<addsize<<std::endl;
      skippedvecs=primalvecs;
      skippedvecs.concat_right(cand_vecs.cols(Range(0,addsize-1)));
      Indexmatrix piv;
      Integer r=skippedvecs.QR_factor(piv);
      bundlevecs.init(skippedvecs.rowdim(),r,0.);
      for(Integer i=0;i<r;i++) bundlevecs(i,i)=1.;
      skippedvecs.Q_times(bundlevecs,r);
    }
    else {
      //they are almost the same, go for Newton
      aug_subg_in_model=0;
      if ((out)&&(print_level>2))
        (*out)<<"trying Newton without enlarging bundle"<<std::endl;
    }
    skippedvecs=topvecs.cols(Range(keepsize,keepsize+skippedsize-1));
    topvecs=primalvecs;
    Ritz_values.reduce_length(keepsize+skippedsize);
  }

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

  
// *****************************************************************************
//                                update_model_top_spectrum2
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum2(bool /*descent_step*/)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with identical eigenvalues
  activedim=1; //dimension of active subspace
  while((activedim<min(topvecs.dim(),primaleigs.dim()))
	&&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	){
    activedim++;
  }
  if ((out)&&(print_level>2)){
    (*out)<<"lmaxproblem::update_model_top_spectrum2(...): "<<activedim;
    (*out)<<" eigenvalues identical ["<<aggregtol<<"]"<<std::endl;
  }

  //--- decide by primal dual indicators whether to enlarge the active subspace  
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    //Real lguessmu=min(0.001,ip(primaleigs,primalZval)/primaleigs.dim());
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }
    //if (!descent_step)
    //  keepsize=min(activedim+3,primaleigs.dim());
    //else
    //  keepsize=activedim;
    keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<"  lguessmu="<<lguessmu<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
    nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize; //will give the number of columns kept in addition to activedim
  skippedsize=min(max(3,minkeepvecs),topvecs.coldim()-activedim);
  if (skippedsize>3){
    //Real skip_Ritzbound=min(1e-2*(fabs(Ritz_values(0))+1.),1000.*(Ritz_values(0)-Ritz_values(activedim)));
    Real skip_Ritzbound=100.*(Ritz_values(0)-Ritz_values(activedim));
    if ((out)&&(print_level>2))
      (*out)<<" mingap="<<Ritz_values(0)-Ritz_values(activedim)<<" skip_Ritzbound="<<skip_Ritzbound;
    while ((skippedsize>3)&&
	   (Ritz_values(0)-Ritz_values(skippedsize+activedim-1)>skip_Ritzbound))
      skippedsize--;
  }
  if ((out)&&(print_level>2))
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
  Ritz_values.reduce_length(activedim+skippedsize);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum2(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_top_spectrum3
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum3(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with identical eigenvalues
  activedim=1; //dimension of active subspace
  while((activedim<min(topvecs.dim(),primaleigs.dim()))
	&&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	){
    activedim++;
  }
  if ((out)&&(print_level>2)){
    (*out)<<"lmaxproblem::update_model_top_spectrum3(...): "<<activedim;
    (*out)<<" eigenvalues identical ["<<aggregtol<<"]"<<std::endl;
  }

  //--- decide by primal dual indicators whether to enlarge the active subspace  
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    //Real lguessmu=min(0.001,ip(primaleigs,primalZval)/primaleigs.dim());
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }
    if (!descent_step)
      keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    else
      keepsize=activedim;
    //keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<"  lguessmu="<<lguessmu<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
    nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize; //will give the number of columns kept in addition to activedim
  skippedsize=min(max(3,minkeepvecs),topvecs.coldim()-activedim);
  if (skippedsize>3){
    //Real skip_Ritzbound=min(1e-2*(fabs(Ritz_values(0))+1.),1000.*(Ritz_values(0)-Ritz_values(activedim)));
    Real skip_Ritzbound=100.*(Ritz_values(0)-Ritz_values(activedim));
    if ((out)&&(print_level>2))
      (*out)<<" allowed="<<skippedsize<<" mingap="<<Ritz_values(0)-Ritz_values(activedim)<<" skip_Ritzbound="<<skip_Ritzbound;
    while ((skippedsize>3)&&
	   (Ritz_values(0)-Ritz_values(skippedsize+activedim-1)>skip_Ritzbound))
      skippedsize--;
  }
  if ((out)&&(print_level>2))
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
  Ritz_values.reduce_length(activedim+skippedsize);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum3(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_top_spectrum4
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum4(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  activedim=1; //dimension of active subspace
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with identical eigenvalues
  while((activedim<min(topvecs.dim(),primaleigs.dim()))
	&&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	){
    activedim++;
  }
  if ((out)&&(print_level>2)){
    (*out)<<"lmaxproblem::update_model_top_spectrum4(...): "<<activedim;
    (*out)<<" eigenvalues identical ["<<aggregtol<<"]"<<std::endl;
  }

  //--- decide by primal dual indicators whether to enlarge the active subspace  
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    Real pmax=max(primaleigs);
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }
    while((activedim<primaleigs.dim())&&
	  (
	   (primal_tapia(activedim)>.8)
	  ||
	   ((primaleigs(activedim)>pmax*lguessmu)&&(primal_tapia(activedim)>.3)&&(dual_tapia(activedim)<.5))
	   )
	  ){
      activedim++;
    }
    if (!descent_step)
      keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    else
      keepsize=activedim;
    //keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<"  lguessmu="<<lguessmu<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
    nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  }

   //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize; //will give the number of columns kept in addition to activedim
  skippedsize=min(max(3,minkeepvecs),topvecs.coldim()-activedim);
  if (skippedsize>3){
    //Real skip_Ritzbound=min(1e-2*(fabs(Ritz_values(0))+1.),1000.*(Ritz_values(0)-Ritz_values(activedim)));
    Real skip_Ritzbound=100.*(Ritz_values(0)-Ritz_values(activedim));
    if ((out)&&(print_level>2))
      (*out)<<" allowed="<<skippedsize<<" mingap="<<Ritz_values(0)-Ritz_values(activedim)<<" skip_Ritzbound="<<skip_Ritzbound;
    while ((skippedsize>3)&&
	   (Ritz_values(0)-Ritz_values(skippedsize+activedim-1)>skip_Ritzbound))
      skippedsize--;
  }
  if ((out)&&(print_level>2))
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
  Ritz_values.reduce_length(activedim+skippedsize);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum4(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  
// *****************************************************************************
//                                update_model_top_spectrum5
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum5(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with identical eigenvalues
  activedim=1; //dimension of active subspace
  while((activedim<min(topvecs.dim(),primaleigs.dim()))
	&&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	){
    activedim++;
  }
  if ((out)&&(print_level>2)){
    (*out)<<"lmaxproblem::update_model_top_spectrum5(...): "<<activedim;
    (*out)<<" eigenvalues identical ["<<aggregtol<<"]"<<std::endl;
  }

  //--- decide by primal dual indicators whether to enlarge the active subspace  
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    Real pmax=max(primaleigs);
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }
    while((activedim<primaleigs.dim())&&
	  (
	   (primal_tapia(activedim)>.8)
	  ||
	   ((primaleigs(activedim)>pmax*lguessmu)&&(primal_tapia(activedim)>.3)&&(dual_tapia(activedim)<.5))
	   )
	  ){
      activedim++;
    }
    if (!descent_step)
      keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    else
      keepsize=activedim;
    //keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<"  lguessmu="<<lguessmu<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
    nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  Indexmatrix delind;
  if ((descent_step)&&(activedim>1)){    
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>max(keepsize-activedim,3))&&
	  (Ritz_values(0)-Ritz_values(activedim+skippedsize-1)>10*(Ritz_values(0)-Ritz_values(activedim)))
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<Ritz_values(0)-1000*(Ritz_values(0)-Ritz_values(activedim))<<" -> ";
    }
    delind.init(Range(activedim+skippedsize,topvecs.coldim()-1));
    
    bundlevecs=primalvecs;
    bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
    tmpvec=primaleigs;
    tmpvec.reduce_length(activedim);
    tmpvec.sqrt();
    bundlevecs.scale_cols(tmpvec);
    
    skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
    tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
    for(Integer i=0;i<skippedsize;i++){
      tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
    }
    tmpvec.sqrt();
    skippedvecs.scale_cols(tmpvec);
    
    tmpvec.init(1,tmpvec.dim(),0.);
    Matrix tmpmax(1,tmpvec.dim(),0.);
    for(Integer i=0;i<y.dim();i++){
      MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
      if (fun==0) break;
      const Coeffmat* Ak=fun->get_coeffmat(i,0);
      Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
      tmpmat%=tmpmat;
      Matrix sr(sumrows(tmpmat));
      for(Integer j=0;j<tmpmax.dim();j++){
	tmpmax(j)=max(tmpmax(j),sr(j));
      }
      tmpvec+=sr;
    }
    if ((out)&&(print_level>2)){
      (*out)<<" importance="<<tmpvec;
      (*out)<<" importmax="<<tmpmax;
    }
    Indexmatrix sind=sortindex(tmpvec);
    Real sumtmpvec=sum(tmpvec);
    Real tmpsum=0.;
    for(Integer i=0;i<sind.dim();i++){
      Integer ind=sind(i);
      if (ind<3) continue;
      Real tmpval=tmpvec(ind);
      if (tmpsum+tmpval<0.1*sumtmpvec){
	delind.concat_below(activedim+ind);
	tmpsum+=tmpval;
      }
    }
	
    /*
    for(;skippedsize>max(keepsize-activedim,3);skippedsize--){
      if (tmpvec(skippedsize-1)<0.1*y.dim())
	delind.concat_below(activedim+skippedsize-1);
    }
    */
    skippedsize=topvecs.coldim()-delind.dim();
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    delind.init(Range(activedim+skippedsize,topvecs.coldim()-1));
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<") delind="<<transpose(delind);
  }
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(delind);
  Ritz_values.delete_rows(delind);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum5(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  



// *****************************************************************************
//                                update_model_top_spectrum6
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum6(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- estimate the dimension of the active subspace
  activedim=0; //dimension of active subspace
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    //primal-dual solution information coming from the QSP
    Real pmax=primaleigs(0);
    Real dmax=primalZval(0);
    //Real lguessmu=min(0.001,ip(primaleigs,primalZval)/primaleigs.dim());
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }

    Real aggregval=lguessmu; //min(aggregtol,.1);
    activedim=1;
    while(activedim<primaleigs.dim()){
      int no_Ritz_gap=(maxRitz-S(activedim,activedim)<min(1e-3,sqrt(lguessmu))*max(1.,fabs(maxRitz)));
      int primal_dropoff= (primaleigs(activedim)<pmax*lguessmu);
      int Ritz_gap=(maxRitz-S(activedim,activedim)>5.*Ritzbound);
      int primaldual_ratio=((maxRitz-S(activedim,activedim)>2.*Ritzbound)&&(primaleigs(activedim)*dmax<aggregval*primalZval(activedim)*pmax));
      if ((!no_Ritz_gap)&&(primal_dropoff||primaldual_ratio||Ritz_gap)){
	if ((out)&&(print_level>2)){
	  (*out)<<"lmaxproblem::update_model_top_spectrum6(...):\n gap located at index "<<activedim;
	  (*out)<<": pdropoff="<<primal_dropoff;
	  (*out)<<" pdratio="<<primaldual_ratio;
	  (*out)<<" Rgap="<<Ritz_gap;
	}
	break;
      }
      activedim++;
    }
    if (!descent_step)
      keepsize=min(activedim+3,primaleigs.dim());
    else{
      keepsize=activedim;
    }
    //keepsize=min(activedim+1,primaleigs.dim());
    if ((out)&&(print_level>2))
      (*out)<<" lguessmu="<<lguessmu<<" -> keepsize="<<keepsize<<std::endl;
    nextsize=keepsize+min(cand_vecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these
  Integer skippedsize; //will give the number of columns kept in addition to activedim

  //compute a basis of the subspace
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values
  
  //find the number of vectors to keep on top of activedim for the quadratic term 
  skippedsize=min(max(3,minkeepvecs),topvecs.coldim()-activedim);
  if (skippedsize>3){
    //Real skip_Ritzbound=min(1e-2*(fabs(Ritz_values(0))+1.),1000.*(Ritz_values(0)-Ritz_values(activedim)));
    Real skip_Ritzbound=100.*(Ritz_values(0)-Ritz_values(activedim));
    if ((out)&&(print_level>2))
      (*out)<<" mingap="<<Ritz_values(0)-Ritz_values(activedim)<<" skip_Ritzbound="<<skip_Ritzbound;
    while ((skippedsize>3)&&
	   (Ritz_values(0)-Ritz_values(skippedsize+activedim-1)>skip_Ritzbound))
      skippedsize--;
  }
  if ((out)&&(print_level>2))
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
  Ritz_values.reduce_length(activedim+skippedsize);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum6(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_top_spectrum7
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum7(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with identical eigenvalues
  activedim=1; //dimension of active subspace
  while((activedim<min(topvecs.dim(),primaleigs.dim()))
	&&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	){
    activedim++;
  }
  if ((out)&&(print_level>2)){
    (*out)<<"lmaxproblem::update_model_top_spectrum7(...): "<<activedim;
    (*out)<<" eigenvalues identical ["<<aggregtol<<"]"<<std::endl;
  }

  //--- decide by primal dual indicators whether to enlarge the active subspace  
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    //primal-dual solution information coming from the QSP
    Real pmax=primaleigs(0);
    Real dmax=primalZval(0);
    //Real lguessmu=min(0.001,ip(primaleigs,primalZval)/primaleigs.dim());
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }

    Real aggregval=lguessmu; //min(aggregtol,.1);
    while(activedim<primaleigs.dim()){
      int primal_dropoff= (primaleigs(activedim)<pmax*lguessmu);
      int Ritz_gap=(maxRitz-S(activedim,activedim)>5.*Ritzbound);
      int primaldual_ratio=((maxRitz-S(activedim,activedim)>2.*Ritzbound)&&(primaleigs(activedim)*dmax<aggregval*primalZval(activedim)*pmax));
      if (primal_dropoff||primaldual_ratio||Ritz_gap){
	if ((out)&&(print_level>2)){
	  (*out)<<"lmaxproblem::update_model_top_spectrum7(...):\n gap located at index "<<activedim;
	  (*out)<<": pdropoff="<<primal_dropoff;
	  (*out)<<" pdratio="<<primaldual_ratio;
	  (*out)<<" Rgap="<<Ritz_gap;
	}
	break;
      }
      activedim++;
    }
    if (!descent_step)
      keepsize=min(activedim+3,primaleigs.dim());
    else
      keepsize=activedim;
    //keepsize=min(activedim+1,primaleigs.dim());
    if ((out)&&(print_level>2))
      (*out)<<" lguessmu="<<lguessmu<<" -> keepsize="<<keepsize<<std::endl;
    nextsize=keepsize+min(cand_vecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize; //will give the number of columns kept in addition to activedim
  skippedsize=min(max(3,minkeepvecs),topvecs.coldim()-activedim);
  if (skippedsize>3){
    //Real skip_Ritzbound=min(1e-2*(fabs(Ritz_values(0))+1.),1000.*(Ritz_values(0)-Ritz_values(activedim)));
    Real skip_Ritzbound=100.*(Ritz_values(0)-Ritz_values(activedim));
    if ((out)&&(print_level>2))
      (*out)<<" mingap="<<Ritz_values(0)-Ritz_values(activedim)<<" skip_Ritzbound="<<skip_Ritzbound;
    while ((skippedsize>3)&&
	   (Ritz_values(0)-Ritz_values(skippedsize+activedim-1)>skip_Ritzbound))
      skippedsize--;
  }
  if ((out)&&(print_level>2))
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
  Ritz_values.reduce_length(activedim+skippedsize);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum7(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_top_spectrum8
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_top_spectrum8(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with identical eigenvalues
  activedim=1; //dimension of active subspace
  while((activedim<min(topvecs.dim(),primaleigs.dim()))
	&&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	){
    activedim++;
  }
  if ((out)&&(print_level>2)){
    (*out)<<"lmaxproblem::update_model_top_spectrum8(...): "<<activedim;
    (*out)<<" eigenvalues identical ["<<aggregtol<<"]"<<std::endl;
  }

  //--- decide by primal dual indicators whether to enlarge the active subspace  
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    //Real lguessmu=min(0.001,ip(primaleigs,primalZval)/primaleigs.dim());
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }
    if (!descent_step)
      keepsize=min(activedim+3,primaleigs.dim());
    else
      keepsize=activedim;
    //keepsize=min(activedim+1,primaleigs.dim());
    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<"  lguessmu="<<lguessmu<<" -> keepsize="<<keepsize<<std::endl;
    }
    nextsize=keepsize+min(cand_vecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize; //will give the number of columns kept in addition to activedim
  skippedsize=min(max(3,minkeepvecs),topvecs.coldim()-activedim);
  if (skippedsize>3){
    //Real skip_Ritzbound=min(1e-2*(fabs(Ritz_values(0))+1.),1000.*(Ritz_values(0)-Ritz_values(activedim)));
    Real skip_Ritzbound=100.*(Ritz_values(0)-Ritz_values(activedim));
    if ((out)&&(print_level>2))
      (*out)<<" mingap="<<Ritz_values(0)-Ritz_values(activedim)<<" skip_Ritzbound="<<skip_Ritzbound;
    while ((skippedsize>3)&&
	   (Ritz_values(0)-Ritz_values(skippedsize+activedim-1)>skip_Ritzbound))
      skippedsize--;
  }
  if ((out)&&(print_level>2))
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
  Ritz_values.reduce_length(activedim+skippedsize);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_top_spectrum8(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_new_spectrum
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_new_spectrum(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- estimate the dimension of the active subspace
  activedim=0; //dimension of active subspace
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    //primal-dual solution information coming from the QSP
    Real pmax=primaleigs(0);
    Real dmax=primalZval(0);
    //Real lguessmu=min(0.001,ip(primaleigs,primalZval)/primaleigs.dim());
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);

    Real aggregval=lguessmu; //min(aggregtol,.1);
    activedim=1;
    while(activedim<primaleigs.dim()){
      int no_Ritz_gap=(maxRitz-S(activedim,activedim)<min(1e-3,sqrt(lguessmu))*max(1.,fabs(maxRitz)));
      int primal_dropoff= (primaleigs(activedim)<pmax*lguessmu);
      int Ritz_gap=(maxRitz-S(activedim,activedim)>5.*Ritzbound);
      int primaldual_ratio=((maxRitz-S(activedim,activedim)>2.*Ritzbound)&&(primaleigs(activedim)*dmax<aggregval*primalZval(activedim)*pmax));
      if ((!no_Ritz_gap)&&(primal_dropoff||primaldual_ratio||Ritz_gap)){
	if ((out)&&(print_level>2)){
	  (*out)<<"lmaxproblem::update_model_new_spectrum(...):\n gap located at index "<<activedim;
	  (*out)<<": pdropoff="<<primal_dropoff;
	  (*out)<<" pdratio="<<primaldual_ratio;
	  (*out)<<" Rgap="<<Ritz_gap;
	}
	break;
      }
      activedim++;
    }
    if (!descent_step)
      keepsize=min(activedim+3,primaleigs.dim());
    else
      keepsize=activedim;
    //keepsize=min(activedim+1,primaleigs.dim());
    if ((out)&&(print_level>2))
      (*out)<<" lguessmu="<<lguessmu<<" -> keepsize="<<keepsize<<std::endl;
    nextsize=keepsize+min(cand_vecs.coldim(),maxaddvecs);
  }

  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these
  Integer skippedsize; //will give the number of columns kept in addition to activedim

  //compute a basis of the subspace
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values
  
  //find the number of vectors to keep on top of activedim for the quadratic term 
  skippedsize=min(max(3,minkeepvecs),topvecs.coldim()-activedim);
  if (skippedsize>3){
    //Real skip_Ritzbound=min(1e-2*(fabs(Ritz_values(0))+1.),1000.*(Ritz_values(0)-Ritz_values(activedim)));
    Real skip_Ritzbound=100.*(Ritz_values(0)-Ritz_values(activedim));
    if ((out)&&(print_level>2))
      (*out)<<" mingap="<<Ritz_values(0)-Ritz_values(activedim)<<" skip_Ritzbound="<<skip_Ritzbound;
    while ((skippedsize>3)&&
	   (Ritz_values(0)-Ritz_values(skippedsize+activedim-1)>skip_Ritzbound))
      skippedsize--;
  }
  if ((out)&&(print_level>2))
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+keepsize-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  
  //reduce topvecs to hold at most activedim plus skippedsize columns
  topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
  Ritz_values.reduce_length(activedim+skippedsize);
 
  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=1;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<cand_vecs.coldim())&&
	 (Ritz_values(0)-cand_eigs(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,cand_vecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-cand_eigs(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(cand_vecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_new_spectrum(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  
// *****************************************************************************
//                                update_model_tapia_scale
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step)
    activedim=1; //reset dimension of active subspace
  //else we only allow to increase activedim

  //--- decide by primal dual indicators the minimum size of the active subspace  
  Integer keepsize=0;  //number of vectors to keep in primalvecs and in the bundle
  Integer nextsize=0;  //maximum number of columns in the new bundle
  Real Ritzbound=0;
  if (primaleigs.dim()>0){
    Real lguessmu=min(0.1,ip(primaleigs,primalZval)/primaleigs.dim());
    Real pmax=max(primaleigs);
    //actual Ritz values corresponding the primal eigenvectors
    Symmatrix S;
    oracle->compute_projection(cand_y,primalvecs,S);
    Real maxRitz=S(0,0);
    Ritzbound=max(lguessmu*(fabs(cand_eigs(0))+1),cand_eigs(0)-maxRitz);
    if ((out)&&(print_level>2)){
      out->precision(10);
      (*out)<<"Ritzbound="<<Ritzbound<<"  primal Ritz values:"<<transpose(diag(S));
    }
    while((activedim<primaleigs.dim())&&
	  (
	   (primal_tapia(activedim)>.8)
	  ||
	   ((primaleigs(activedim)>pmax*lguessmu)&&(primal_tapia(activedim)>.3)&&(dual_tapia(activedim)<.5))
	   )
	  ){
      activedim++;
    }

    if (!descent_step)
      keepsize=min(activedim+maxkeepvecs,primaleigs.dim());
    else
      keepsize=activedim;

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<"  lguessmu="<<lguessmu<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
    nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(0)-Ritz_values(activedim+skippedsize-1)>0.01*max(1.,Ritz_values(0)))
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<Ritz_values(0)-0.1*max(1.,Ritz_values(0))<<" -> ";
    }
    
    bundlevecs=primalvecs;
    bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
    tmpvec=primaleigs;
    tmpvec.reduce_length(activedim);
    tmpvec.sqrt();
    bundlevecs.scale_cols(tmpvec);
    
    skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
    tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
    for(Integer i=0;i<skippedsize;i++){
      tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
    }
    tmpvec.sqrt();
    skippedvecs.scale_cols(tmpvec);
    
    tmpvec.init(1,tmpvec.dim(),0.);
    Matrix tmpmax(1,tmpvec.dim(),0.);
    for(Integer i=0;i<y.dim();i++){
      MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
      if (fun==0) break;
      const Coeffmat* Ak=fun->get_coeffmat(i,0);
      Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
      tmpmat%=tmpmat;
      Matrix sr(sumrows(tmpmat));
      for(Integer j=0;j<tmpmax.dim();j++){
	tmpmax(j)=max(tmpmax(j),sr(j));
      }
      tmpvec+=sr;
    }
    if ((out)&&(print_level>2)){
      (*out)<<" importance="<<tmpvec;
      (*out)<<" importmax="<<tmpmax;
    }
    Integer i=min(keepsize-activedim,tmpvec.dim());
    bool increase_flag=false;
    while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
      keepsize++;
      i++;
      increase_flag=true;
    }
    if ((out)&&(print_level>2)){
      if (increase_flag){
	(*out)<<" increase keepsize to "<<keepsize;
      }
    }
    i=min(keepsize-activedim+3,tmpvec.dim());
    while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
      i++;
    }
    skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale1
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale1(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim

  Integer keepsize=0;  
  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    /*
    while ((Ritzcnt<3)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    */
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.1*max(1.,fabs(Ritz_values(0))),100.*(Ritz_values(0)-Ritzgapavg));
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	tmpmat%=tmpmat;
	Matrix sr(sumrows(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale1(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale2
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale2(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim

  Integer keepsize=0;  
  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<3)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.01*max(1.,fabs(Ritz_values(0))),10.*(Ritz_values(0)-Ritzgapavg));
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale2(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale3
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale3(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim

  Integer keepsize=0;  
  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    if (!descent_step)
      keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<3)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.01*max(1.,fabs(Ritz_values(0))),10.*(Ritz_values(0)-Ritzgapavg));
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	tmpmat%=tmpmat;
	Matrix sr(sumrows(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      /*
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      */
      Integer i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale3(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  
// *****************************************************************************
//                                update_model_tapia_scale4
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale4(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim

  Integer keepsize=0;  
  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while(
	  (activedim<primaleigs.dim())
	  &&
	  (
	   (primal_tapia(activedim)>.9)
	   ||
	   ((primal_tapia(activedim)>2.*tapia_factor)&&(primal_tapia(activedim)>dual_tapia(activedim)))
	   )
	  ){
      activedim++;
    }

    if ((out)&&(print_level>2)){
      (*out)<<" tapia="<< activedim;
    }
    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real avgRitzgap=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<3)&&(activedim+Ritzcnt<Ritz_values.dim())){
      avgRitzgap+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    avgRitzgap=Ritz_values(0)-avgRitzgap/Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.01*max(1.,fabs(Ritz_values(0))),5.*avgRitzgap);
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" avgRitzgap="<<avgRitzgap<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	tmpmat%=tmpmat;
	Matrix sr(sumrows(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale4(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  


// *****************************************************************************
//                                update_model_tapia_scale5
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale5(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim and keepsize

  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<3)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.01*(1.+fabs(Ritz_values(0))),10.*(Ritz_values(0)-Ritzgapavg));
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=min(max(i,keepsize+3-activedim),topvecs.coldim()-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale5(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale6
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale6(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim

  Integer keepsize=0;  
  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    /*
    while ((Ritzcnt<3)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    */
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.01*max(1.,fabs(Ritz_values(0))),10.*(Ritz_values(0)-Ritzgapavg));
    skippedsize=min(max(2*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      Integer i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale6(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale7
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale7(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim

  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<activedim)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-10.*(Ritz_values(0)-Ritzgapavg);
    skippedsize=min(max(5*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      /*
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      */
      Integer i=0;
      bool increase_flag=false;
      Real trsum=0.;
      while((i<tmpvec.dim())&&(tmpvec(i)>1e-3*trsum)){
	trsum+=tmpvec(i);
	if ((tmpvec(i)/y.dim()>1.)&&(activedim+i==keepsize)&&(keepsize<primaleigs.dim())){
	  keepsize++;
	  increase_flag=true;
	}
	//if (activedim+i<keepsize){
	//  trsum=0.;
	//}
	i++;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" trsum="<<trsum;
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale7(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale8
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale8(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim

  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<activedim)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-10.*(Ritz_values(0)-Ritzgapavg);
    skippedsize=min(max(5*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      skippedsize=min(2*keepsize,topvecs.coldim())-activedim;
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale8(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  
// *****************************************************************************
//                                update_model_tapia_scale9
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale9(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim and keepsize

  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    skippedsize=min(2*primaleigs.dim()-keepsize,topvecs.coldim()-activedim);
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=max(i,min(keepsize+3,topvecs.coldim())-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale9(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale10
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale10(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim and keepsize

  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<3)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.01*(1.+fabs(Ritz_values(0))),10.*(Ritz_values(0)-Ritzgapavg));
    skippedsize=min(max(2*primaleigs.dim(),Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>3*primaleigs.dim()/2-activedim)&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=min(max(i,Integer(3*keepsize/2-activedim)),topvecs.coldim()-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(3*keepsize/2-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale10(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  
// *****************************************************************************
//                                update_model_tapia_scale11
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale11(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim and keepsize

  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Matrix resvec(topvecs.col(0));
  genmult(primalvecs,resvec,tmpvec,1.,0.,1);
  genmult(primalvecs,tmpvec,resvec,-1.,1.);
  if ((out)&&(print_level>2)){
    (*out)<<" n2(resvec)="<<norm2(resvec);
  }
  resvec/=norm2(resvec);
  Real ipC;
  oracle->gramip(resvec,ipC,tmpvec);
  Real cd0=ublmax-(ipC-ip(tmpvec,y));
  curve_delta=max(cd0,curve_delta);
  Integer skippedsize;
  tmpvec.init(y);
  tmpvec-=cand_y;
  Real cd1=(ub_fun_val-(cand_subg_val+ip(cand_subg,tmpvec)))/lmax_multiplier;
  curve_delta=max(cd1,curve_delta);
  if ((out)&&(print_level>2)){
    (*out)<<" cd0="<<cd0<<" cd1="<<cd1<<" curve_delta="<<curve_delta;
  }
 
  if ((activedim>1)&&(descent_step)){    
    Real skipcutoff=Ritz_values(0)-curve_delta;
    skippedsize=topvecs.coldim()-activedim;
    while(
	  (skippedsize>keepsize-activedim+3)&&
	  (Ritz_values(activedim+skippedsize-2)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      /*
      i=min(keepsize-activedim+3,tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      */
      skippedsize=min(max(skippedsize,keepsize+3-activedim),topvecs.coldim()-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,3),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }

  if (descent_step){
    curve_delta=0.;
  }
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale11(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  

// *****************************************************************************
//                                update_model_tapia_scale12
// *****************************************************************************

//generate the next cutting model containing at least the convex combination 
//of the two current subgradients of eval_function and (re)eval_augmodel
//returns: 0 ... if the information is available 
//               (if eval_augmodel  returned 1, the information will not
//                satisfy the precision requirements)
//         1 ... if the desired information is not available

int LmaxProblem::update_model_tapia_scale12(bool descent_step)
{
  //--- initialization if there is no bundle yet
  if (bundlevecs.coldim()==0){
    if (cand_vecs.coldim()==0) return 1;
    model_changed=1;
    Integer nc=min(cand_vecs.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    Indexmatrix ind;
    sortindex(cand_eigs,ind);
    bundlevecs=cand_vecs.cols(ind(Range(ind.dim()-1,ind.dim()-nc,-1)));
    return 0;
  }

  //--- check whether new information is available at all
  if (!cand_available) return 1;

  //--- if the trace is not fixed and there is no active vector, do nothing
  if ((cand_ublmax==0.)&&(trace_stat!=SDPtrace_fixed)){
    if  ((out)&&(print_level>2)){
      (*out)<<"no SDP-update:"; 
      (*out)<<" neweig="<<cand_ublmax;
      (*out)<<std::endl;
      if (primaleigs.dim()>0){
	(*out)<<" primaleigs=";
	for(Integer i=0;i<primaleigs.dim();i++){
	  out->precision(6);(*out)<<" "<<primaleigs(i);
	}
	(*out)<<std::endl;
      }
      if (aggrcoeff.dim()>0){
	(*out)<<" aggrcoeff=";
	for(Integer i=0;i<aggrcoeff.dim();i++){
	  out->precision(6);(*out)<<" "<<aggrcoeff(i);
	}
	(*out)<<std::endl;
      }
    }
    return 0;
  }
  model_changed=1;
  
  //--- here we allow for only one aggregate, aggregate all aggrsubgrads into one
  if (aggrcoeff.dim()>1){
    Integer aggrindex=0;
    Integer n=aggrsubgrads.rowdim();
    for(Integer indi=1;indi<aggrcoeff.dim();indi++){
      Real b=aggrcoeff(indi);
      if (b<eps_Real) continue; //nothing to aggregate
      Real a=aggrcoeff(aggrindex);
      Real s=a+b;
      a/=s;
      b/=s;
      subgCvalues(aggrindex)=subgCvalues(aggrindex)*a+subgCvalues(indi)*b;
      mat_xbpeya(n,aggrsubgrads.get_store()+aggrindex*n,
		 aggrsubgrads.get_store()+indi*n,b,a);
      if (Integer(primal.size())==aggrsubgrads.coldim()){
	primal[aggrindex]->aggregate_primal_data(a,b,*(primal[indi]));
      }
      aggrcoeff(aggrindex)+=aggrcoeff(indi);  
    }
    delete_subgrads(Range(1,aggrsubgrads.coldim()-1));
  } //endif (aggrcoeff.dim()>1)

  //--- form the subspace of all collected vectors, compute the
  //    the projected spectrum and keep the best of these

  //compute a basis of the subspace
  Integer oldskippedsize=skippedvecs.coldim();
  tmpmat=cand_vecs;
  tmpmat.concat_right(topvecs);
  //tmpmat.concat_right(primalvecs);
  //tmpmat.concat_right(bundlevecs);
  tmpmat.concat_right(skippedvecs);
  Indexmatrix piv;
  Integer r1=tmpmat.QR_factor(piv);
  if ((out)&&(print_level>2))
    (*out)<<" subspacedim="<<r1<<"("<<tmpmat.coldim()<<")";
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
  for(Integer i=0;i<tmpmat.coldim();i++) skippedvecs(i,i)=1.;
  tmpmat.Q_times(skippedvecs,r1);
  
  //compute the eigenvalue decomposition of the projected matrix (tmpmat is now free again)
  oracle->evaluate_projection(cand_y,skippedvecs,1e-10,1e20,tmpmat,Ritz_values);
  genmult(skippedvecs,tmpmat,topvecs);
  if ((out)&&(print_level>2))
    (*out)<<" Ritz_values="<<transpose(Ritz_values)<<std::endl;
  
  //topvecs now holds the Ritz-vectors sorted nonincreasingly by Ritz-values


  //--- estimate the dimension of the active subspace, start with tapia indicators
  if (descent_step){
    activedim=1; //reset dimension of active subspace
    keepsize=activedim; //number of vectors to keep in primalvecs and in the bundle
  }
  //else we only allow to increase activedim and keepsize

  if (primaleigs.dim()>0){
    //actual Ritz values corresponding the primal eigenvectors
    while((activedim<primaleigs.dim())&&(primal_tapia(activedim)>.8)){
      activedim++;
    }

    keepsize=min(max(keepsize,activedim+maxkeepvecs),primaleigs.dim());

    //possibly enlarge the dimension of the active subspace for "identical" eigenvalues
    while((activedim<min(topvecs.dim(),primaleigs.dim()))
	  &&(Ritz_values(0)-Ritz_values(activedim)<aggregtol*max(1.,fabs(Ritz_values(0))))
	  ){
      activedim++;
    }
    keepsize=max(activedim,keepsize);

    if ((out)&&(print_level>2)){
      (*out)<<" gap located at index "<< activedim;
      (*out)<<" maxkeepvecs="<<maxkeepvecs<<" -> keepsize="<<keepsize<<std::endl;
    }
  }

  //---- find the number of vectors to keep on top of activedim for the quadratic term 
  Integer skippedsize;
  if ((activedim>1)&&(descent_step)){    
    Real Ritzgapavg=Ritz_values(activedim);
    Integer Ritzcnt=1;
    while ((Ritzcnt<1)&&(activedim+Ritzcnt<Ritz_values.dim())){
      Ritzgapavg+=Ritz_values(activedim+Ritzcnt);
      Ritzcnt++;
    }
    Ritzgapavg/=Real(Ritzcnt);
    Real skipcutoff=Ritz_values(0)-min(0.01*(1.+fabs(Ritz_values(0))),10.*(Ritz_values(0)-Ritzgapavg));
    skippedsize=min(max(3*activedim,Integer(sqrt(topvecs.rowdim()))),topvecs.coldim()-activedim);
    while(
	  (skippedsize>keepsize-activedim+max(3,minkeepvecs))&&
	  (Ritz_values(activedim+skippedsize-1)<skipcutoff)
	  )
      skippedsize--;
    if ((out)&&(print_level>2)){
      (*out)<<" skipcutoff="<<skipcutoff<<" -> "<<skippedsize;
    }
        
    MatrixSDPfunction* fun=dynamic_cast<MatrixSDPfunction*>(oracle);
    if (fun!=0) {
      //estimate importance of each vector for scaling by its contribution to the trace 
      bundlevecs=primalvecs;
      bundlevecs.delete_cols(Range(activedim,bundlevecs.coldim()-1));
      tmpvec=primaleigs;
      tmpvec.reduce_length(activedim);
      tmpvec.sqrt();
      bundlevecs.scale_cols(tmpvec);
    
      skippedvecs.init(topvecs.rowdim(),skippedsize,topvecs.get_store()+topvecs.rowdim()*activedim);
      tmpvec.newsize(skippedsize,1); chk_set_init(tmpvec,1);
      for(Integer i=0;i<skippedsize;i++){
	tmpvec[i]=1./(Ritz_values(0)-Ritz_values(activedim+i));
      }
      tmpvec.sqrt();
      skippedvecs.scale_cols(tmpvec);

      tmpvec.init(1,tmpvec.dim(),0.);
      Matrix tmpmax(1,tmpvec.dim(),0.);
      for(Integer i=0;i<y.dim();i++){
	const Coeffmat* Ak=fun->get_coeffmat(i,0);
	Ak->left_right_prod(bundlevecs,skippedvecs,tmpmat);
	Matrix sr(colsip(tmpmat));
	for(Integer j=0;j<tmpmax.dim();j++){
	  tmpmax(j)=max(tmpmax(j),sr(j));
	}
	tmpvec+=sr;
      }
      if ((out)&&(print_level>2)){
	(*out)<<" importance="<<tmpvec;
	(*out)<<" importmax="<<tmpmax;
      }
      //
      Integer i=min(keepsize-activedim,tmpvec.dim());
      bool increase_flag=false;
      while((keepsize<primaleigs.dim())&&(i<tmpvec.dim())&&(tmpvec(i)/y.dim()>1.)){
	keepsize++;
	i++;
	increase_flag=true;
      }
      if ((out)&&(print_level>2)){
	if (increase_flag){
	  (*out)<<" increase keepsize to "<<keepsize;
	}
      }
      i=min(keepsize-activedim+max(3,minkeepvecs),tmpvec.dim());
      while((i<tmpvec.dim())&&(tmpvec(i)/y.dim()>0.1)){
	i++;
      }
      skippedsize=min(max(i,keepsize+max(3,minkeepvecs)-activedim),topvecs.coldim()-activedim);
    }
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  else {
    skippedsize=min(topvecs.coldim()-activedim,max(max(keepsize-activedim,max(3,minkeepvecs)),oldskippedsize));
    topvecs.delete_cols(Range(activedim+skippedsize,topvecs.coldim()-1));
    Ritz_values.reduce_length(topvecs.coldim());
  }
  if ((out)&&(print_level>2)){
    (*out)<<" skippedsize="<<skippedsize<<" maxgap="<<Ritz_values(0)-Ritz_values(skippedsize+activedim-1)<<" (lmax="<<Ritz_values(0)<<")"<<std::endl;
  }
  
 
  //---- aggregate the remaining primal vectors
  if (keepsize<primaleigs.dim()){
    Indexmatrix ind(Range(keepsize,primaleigs.dim()-1));
    if  ((out)&&(print_level>2)){
      (*out)<<" aggrdim="<<ind.dim();out->precision(4);
    }
    //find colummn to aggregate to
    Integer aggrindex;   //index of column to aggregate to
    int new_aggrcol=0;   //set to true if a new column should be opened
    if (aggrsubgrads.coldim()<maxaggrcols){
      aggrindex=aggrsubgrads.coldim();
      new_aggrcol=1;
    }
    else { //choose the one with minimal contribution
      min(aggrcoeff,&aggrindex);
    }
	
    //aggregate to the column aggrindex of aggrsubgrads
    if (ind.dim()<primalvecs.coldim()){  //only a subset of the columns
      if (!new_aggrcol){
	update_subgrad(primalvecs.cols(ind),primaleigs(ind),aggrindex);
      }
      else {
	init_subgrad(primalvecs.cols(ind),primaleigs(ind));
      }
      primalvecs.delete_cols(ind);
      primaleigs.delete_rows(ind);
    }
    else { //aggregate all columns
      if (!new_aggrcol){
	update_subgrad(primalvecs,primaleigs,aggrindex);
      }
      else {
	init_subgrad(primalvecs,primaleigs);
      }
      primalvecs.init(0,0,0.);
      primaleigs.init(0,0,0.);
    }
  }

  //--- select new bundlevectors

  //--- compute an orthonormal basis of bundle and cand_vecs
  //    first collect those that have to be included to ensure convergence
  tmpmat=primalvecs;
  tmpmat.concat_right(cand_vecs.col(0)); 
  //form an orthonormal basis of these
  r1=tmpmat.QR_factor(piv);
  tmpmat.delete_cols(Range(r1,tmpmat.coldim()-1));
  piv.reduce_length(r1);
  if (r1==primalvecs.dim()) {
    bundlevecs=primalvecs;
    if ((out)&&(print_level>2)){
      (*out)<<"new vec linearly dependent"<<std::endl;
    }
  }
  else {
    bundlevecs.init(tmpmat.rowdim(),tmpmat.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim();i++) bundlevecs(i,i)=1.;
    tmpmat.Q_times(bundlevecs,r1);
  }
  //bundlevecs now holds an orthonormal basis of the essential subspace
    

  //to enlarge this, orthogonalize the remaining vectors against this basis
  Integer nextsize=keepsize+min(topvecs.coldim(),maxaddvecs);
  Integer maxaddcols=max(nextsize-bundlevecs.coldim(),0);
  Integer firstcol=0;
  if ((out)&&(print_level>2))
    (*out)<<" maxaddcols="<<maxaddcols;
  while ((maxaddcols>0)&&(firstcol<topvecs.coldim())&&
	 (Ritz_values(0)-Ritz_values(firstcol)<=10.*(Ritz_values(0)-Ritz_values(activedim)))){
    Integer lastcol=min(firstcol+maxaddcols,topvecs.coldim())-1; 
    //Integer lastcol=firstcol; 
    while ((lastcol>firstcol)
	   &&(Ritz_values(0)-Ritz_values(lastcol)>10.*(Ritz_values(0)-Ritz_values(activedim))))
	lastcol--;
    Integer r2=tmpmat.QR_concat_right(topvecs.cols(Range(firstcol,lastcol)),piv,r1);
    firstcol=lastcol+1;
    tmpmat.delete_cols(Range(r2,tmpmat.coldim()-1));
    piv.reduce_length(r2);
    maxaddcols-=(r2-r1);
    r1=r2;
  }
  if ((out)&&(print_level>2)){
    (*out)<<" reached firstcol="<<firstcol;
  }
  if(tmpmat.coldim()>bundlevecs.coldim()){
    skippedvecs.init(tmpmat.rowdim(),tmpmat.coldim()-bundlevecs.coldim(),0.);
    for(Integer i=0;i<tmpmat.coldim()-bundlevecs.coldim();i++) 
      skippedvecs(bundlevecs.coldim()+i,i)=1.;
    tmpmat.Q_times(skippedvecs,tmpmat.coldim());
    bundlevecs.concat_right(skippedvecs);
    assert(norm2(Diag(Matrix(bundlevecs.coldim(),1,1.))-transpose(bundlevecs)*bundlevecs)<1e-8);
  }
  
  if ((out)&&(print_level>2)){
    oracle->evaluate_projection(cand_y,bundlevecs,1e-10,1e20,tmpmat,tmpvec);
    (*out)<<"\n bundle Ritz_values="<<transpose(tmpvec);
  }

  //split the new Ritzvectors into top block and skipped block
  skippedvecs.init(topvecs.rowdim(),topvecs.coldim()-activedim,topvecs.get_store()+topvecs.rowdim()*activedim);
  topvecs.delete_cols(Range(activedim,topvecs.coldim()-1));

  if ((out)&&(print_level>0))
    (*out)<<"  LmaxProblem::update_model_tapia_scale12(...): bundlesize="<<bundlevecs.coldim()<<std::endl;

  return 0;
}  
  




  
// *****************************************************************************
//                              intersect_box
// *****************************************************************************

int LmaxProblem::intersect_box(Indexmatrix& inbindex,Matrix& inlb,Matrix& inub)
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

Real LmaxProblem::lb_function(const Matrix& iny)
{
  if ((aug_subg_in_model)&&(aug_available)){
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

Real LmaxProblem::lb_model(const Matrix& iny)
{
  if ((aug_subg_in_model)&&(aug_available)){
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

int LmaxProblem::start_augmodel(QP_Block*& blockp)
{
  //--- if block is not yet available, create one
  if (bundlevecs.coldim()==0){
    if (eigvec.coldim()==0) return 1;
    Integer nc=min(eigvec.coldim()-1,max(maxaddvecs,max(minkeepvecs,Integer(3))));
    bundlevecs=eigvec.cols(Indexmatrix(Range(0,nc)));
  }
  aug_available=0;
  block.init_block(subgCvalues.dim(),Indexmatrix(0,0,Integer(0)),Indexmatrix(1,1,bundlevecs.coldim()),lmax_multiplier,trace_stat!=SDPtrace_fixed);
  
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

int LmaxProblem::get_row(Integer index_y,Matrix& row,Real& bi,Integer startindex) const
{
  if (index_y==-1){ //cost coefficients
    mat_xey(subgCvalues.dim(),
	    row.get_store()+startindex,
	    subgCvalues.get_store());
    oracle->project_C(bundlevecs,tmpsym);
    svec(tmpsym,tmpvec);
    mat_xey(tmpvec.dim(),
	    row.get_store()+startindex+subgCvalues.dim(),
	    tmpvec.get_store());
  }
  else { //coefficient matrix
    bi+=(oracle->rhs())(index_y);
    mat_xemy(aggrsubgrads.coldim(),
	    row.get_store()+startindex,1,
	    aggrsubgrads.get_store()+index_y,aggrsubgrads.rowdim());
    oracle->project(index_y,bundlevecs,tmpsym);
    svec(tmpsym,tmpvec);
    mat_xemy(tmpvec.dim(),
	    row.get_store()+startindex+aggrsubgrads.coldim(),
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

int LmaxProblem::make_aug_linmodel(Real* in_aug_linconst,Matrix* in_aug_subg,bool* conemult_increased,Real* function_value)  
{
  block.get_linx(aggrcoeff);
  block.get_X(0,tmpsym);
  tmpsym.eig(tmpmat,primaleigs,false);
  {
    Integer i=primaleigs.dim()-1;
    while((i>=0)&&(primaleigs(i)<0.)){
      if (out) {
        (*out)<<"*** WARNING: LmaxProblem::make_aug_linmodel(): primaleigs("<<i<<")="<<primaleigs(i)<<", setting to 0"<<std::endl;
      }
      primaleigs(i)=0.;
      --i;
    }
  }
  genmult(bundlevecs,tmpmat,primalvecs);

  
  block.get_Z(0,tmpsym);
  primalZval.newsize(primaleigs.dim(),1); chk_set_init(primalZval,1);
  for(Integer i=0;i<primaleigs.dim();i++){
    tmpvec=tmpmat.col(i);
    primalZval(i)=ip(tmpvec,tmpsym*tmpvec);
  }

  block.get_old_X(0,tmpsym);
  Matrix prev_primaleigs;
  tmpsym.eig(tmpmat,prev_primaleigs,false);
  block.get_old_Z(0,tmpsym);
  Matrix prev_Zval(prev_primaleigs.dim(),1);chk_set_init(prev_Zval,1);
  primal_tapia.newsize(primaleigs.dim(),1); chk_set_init(primal_tapia,1);
  dual_tapia.newsize(primaleigs.dim(),1); chk_set_init(dual_tapia,1);
  for(Integer i=0;i<primaleigs.dim();i++){
    tmpvec=tmpmat.col(i);
    prev_Zval(i)=ip(tmpvec,tmpsym*tmpvec);
    primal_tapia(i)=primaleigs(i)/prev_primaleigs(i);
    dual_tapia(i)=primalZval(i)/prev_Zval(i);
  }
  Real mu=ip(primaleigs,primalZval)/primaleigs.dim();
  Real old_mu=ip(prev_primaleigs,prev_Zval)/primaleigs.dim();
  tapia_factor=mu/old_mu;

  if ((out)&&(print_level>2)){
    for(Integer i=0;i<primaleigs.dim();i++){
      (*out)<<" "<<i<<":"<<primaleigs(i)<<","<<primalZval(i);
      (*out)<<" : "<<prev_primaleigs(i)<<","<<prev_Zval(i);
      (*out)<<"   ("<<primal_tapia(i);
      (*out)<<","<<dual_tapia(i)<<","<<tapia_factor<<")"<<std::endl;
    }
  }

  tmpvec=primaleigs;
  tmpvec.sqrt();  
  tmpmat=primalvecs;
  tmpmat.scale_cols(tmpvec);
  oracle->gramip(tmpmat,aug_linconst,aug_subg);
  
  if (aggrcoeff.dim()>0) {
    genmult(aggrsubgrads,aggrcoeff,aug_subg,1.,1.);
    aug_linconst+=ip(aggrcoeff,subgCvalues);
  }
  xbpeya(aug_subg,oracle->rhs(),1.,-1.);
  
  aug_subg_in_model=1;
  aug_available=1;

  if (in_aug_linconst) *in_aug_linconst+=aug_linconst;
  if (in_aug_subg) *in_aug_subg+=aug_subg;

  //--- check whether it is necessary to increase the cone multiplier
  
  if ((out)&&(print_level>1)) 
    (*out)<<"\n  LmaxProblem: trace="<<sum(primaleigs)+sum(aggrcoeff)<<" mult="<<lmax_multiplier<<std::endl;

  if ((conemult_increased)&&((trace_stat==SDPtrace_unbounded)&&(sum(aggrcoeff)+sum(primaleigs)>.95*lmax_multiplier))){
    lmax_multiplier=max(1.,lmax_multiplier*2);
    if ((out)&&(print_level>1)) 
      if (out) (*out)<<"             multiplier increased to "<<lmax_multiplier<<std::endl;
    *conemult_increased=true;
    block.adjust_trace(lmax_multiplier);
    if (ublmax>0.){
      Real lmax=max(0.,eigval(0));
      Real ip_b_iny=ip(oracle->rhs(),y);
      subg_val=ip_b_iny+lmax*lmax_multiplier;
      ub_fun_val=ip_b_iny+ublmax*lmax_multiplier;
      if (lmax!=0.){
	oracle->gramip_opA(eigvec,subg);
	xbpeya(subg,oracle->rhs(),1.,-lmax_multiplier);
      }
      if (function_value==0) return 2;
      *function_value=ub_fun_val;      
    }
  }
  return 0;
} 
 
// *****************************************************************************
//                                adjust_multiplier
// *****************************************************************************

int LmaxProblem::adjust_multiplier()  
{
  if (trace_stat!=SDPtrace_unbounded) return 0;
  if (!center_available)  return 1;
  if (block.get_linx(tmpvec)) return 1;
  if (block.get_X(0,tmpsym)) return 1;
  if ((out)&&(print_level>1)) 
    (*out)<<"\n  LmaxProblem: trace="<<sum(tmpvec)+trace(tmpsym)<<" mult="<<lmax_multiplier;
  lmax_multiplier=min(lmax_multiplier,max(1.,2.*(sum(tmpvec)+trace(tmpsym))));
  oracle->adjust_multiplier(lmax_multiplier);
  if ((out)&&(print_level>1)) 
    (*out)<<" newmult="<<lmax_multiplier<<std::endl;
  block.adjust_trace(lmax_multiplier);
  if (ublmax>0.){
    Real lmax=max(0.,eigval(0));
    Real ip_b_iny=ip(oracle->rhs(),y);
    subg_val=ip_b_iny+lmax*lmax_multiplier;
    ub_fun_val=ip_b_iny+ublmax*lmax_multiplier;
    if (lmax!=0.){
      oracle->gramip_opA(eigvec,subg);
      xbpeya(subg,oracle->rhs(),1.,-lmax_multiplier);
    }
  }
  return 0;
}
 
// *****************************************************************************
//                                change_variables
// *****************************************************************************

//change variables as described in ChangeVariableInfo
//returns 0 on success, 1 on failure
//(this is not needed in the bundle framework, routine may return 1!)

int LmaxProblem::change_variables(ChangeVarInfo* cvp)
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
    
    //try to extend aggregate subgradient vectors
    int no_subg_extension=0;
    tmpind.init(Range(y.dim()-n_append,y.dim()-1));
    tmpmat.init(0,0,0.);
    for(Integer i=0;i<aggrsubgrads.coldim();i++){
      SDPPrimal* pd=0;
      if (i<Integer(primal.size())) { 
	pd=primal[i];
      }
      if (oracle->primalip_opA(pd,tmpind,tmpvec)){
	if (out) (*out)<<"**** WARNING: LmaxProblem::add_variable(): subgradient extension failed"<<std::endl;
	no_subg_extension=1;
	break;
      }
      tmpmat.concat_right(tmpvec);
    }
    if (no_subg_extension) {
      for(unsigned int i=0;i<primal.size();i++){
	delete primal[i];
      }
      primal.clear();
      subgCvalues.init(0,0,0.);
      aggrsubgrads.init(0,0,0.);
      aggrcoeff.init(0,0,0.);
      if ((trace_stat==SDPtrace_fixed)&&(primaleigs.dim()>0)){
	Real sumd=sum(primaleigs);
	if (sumd>eps_Real){
	  primaleigs*=lmax_multiplier/sumd;
	}
	else {	
	  primaleigs.init(primaleigs.dim(),1,lmax_multiplier/primaleigs.dim());
	}
      }
    }
    else {
      if (aggrsubgrads.rowdim()+tmpmat.rowdim()==y.dim()){
	aggrsubgrads.concat_below(tmpmat);
      }
      else {
	subgCvalues.init(0,0,0.);
	aggrsubgrads.init(0,0,0.);
	aggrcoeff.init(0,0,0.);
      }
    }  
    
    if (center_available){
      //find a good subgradient for restarting the method
      if ((primalvecs.dim()>0)&&(aug_subg_in_model)){
	tmpvec=primaleigs;
	if (aggrcoeff.dim()==0){
	  tmpvec*=lmax_multiplier/sum(primaleigs);
	}
	tmpvec.sqrt();
	tmpmat=primalvecs;
	tmpmat.scale_cols(tmpvec);
	oracle->gramip(tmpmat,subg_val,subg);
	if (aggrcoeff.dim()>0){
	  genmult(aggrsubgrads,aggrcoeff,subg,1.,1.,0,1);
	  subg_val+=ip(subgCvalues,aggrcoeff);
	}
	xbpeya(subg,oracle->rhs(),1.,-1.);
	subg_val+=ip(subg,y);
      }
      else if (bundlevecs.dim()) {
	//make sure that the subgradient in the center is contained in the model
	tmpvec=Matrix(bundlevecs.coldim(),1,sqrt(lmax_multiplier/bundlevecs.coldim()));
	tmpmat=bundlevecs;
	tmpmat.scale_cols(tmpvec);
	oracle->gramip(tmpmat,subg_val,subg);
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
    if (aggrsubgrads.dim()!=0){
      aggrsubgrads=aggrsubgrads.rows(cp->assign_ind);
    }
    if (center_available) {
      subg=subg.rows(cp->assign_ind);
    }

    problem_modified=1;
    cand_available=0;
    aug_available=0;
    aug_subg_in_model=0;
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
    if(aggrsubgrads.dim()!=0){
      aggrsubgrads.delete_rows(cp->del_index);
    }
    
    if (center_available) {
      subg.delete_rows(cp->del_index);
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

int LmaxProblem::recompute_center()
{
  if (!center_available) {
    int status=eval_function(y,max_Real,1e-5);
    if (status) return status;
    int pm=problem_modified;
    status = do_step();
    problem_modified=pm;
    curve_delta=0.;
    if (status) return status;
  }
  
  center_available=1;
  return 0;  
}

// *****************************************************************************
//                                get_approximate_primal
// *****************************************************************************

int LmaxProblem::get_approximate_primal(PrimalData& pd) const 
  {
    if (!aug_available) return -1;
    SDPPrimal* p=dynamic_cast<SDPPrimal*>(&pd);
    if (p==0) return 1;
    if (Integer(primal.size())!=aggrcoeff.dim()) return 1;
    if(primaleigs.dim()>0){
      tmpvec=primaleigs;
      tmpvec.sqrt();  
      tmpmat=primalvecs;
      tmpmat.scale_cols(tmpvec);
      int retval=p->assign_Gram_matrix(tmpmat);
      if (retval) return retval; 
      for(Integer i=0;i<aggrcoeff.dim();i++){
	retval=p->aggregate_primal_data(1.,aggrcoeff(i),*primal[i]);
	if (retval) return retval;
      }
    }
    else{ 
      if (aggrcoeff.dim()>0){
	int retval=p->assign_primal_data(*primal[0]);
	if (retval) return retval; 
      }
      if (aggrcoeff.dim()>1){
	int retval=p->aggregate_primal_data(aggrcoeff(0),aggrcoeff(1),*primal[1]);
	if (retval) return retval; 
      }
      else {
	int retval=p->aggregate_primal_data(aggrcoeff(0),0.,*primal[0]);
	return retval; 
      }
      for(Integer i=2;i<aggrcoeff.dim();i++){
	int retval=p->aggregate_primal_data(1.,aggrcoeff(i),*primal[i]);
	if (retval) return retval;
      }
    }
    return 0;
  }


// *****************************************************************************
//                                get_center_primal
// *****************************************************************************

int LmaxProblem::get_center_primal(PrimalData& pd) const 
  {
    SDPPrimal* p=dynamic_cast<SDPPrimal*>(&pd);
    if (p==0) return 1;
    if (eigvec.dim()==0) return 1;
    Integer maxind;
    Real d=max(eigval,&maxind);
    int retval;
    if ((d<0.)&&(trace_stat!=SDPtrace_fixed)){
      tmpmat.init(eigvec.rowdim(),1,0.);
      retval=p->assign_Gram_matrix(tmpmat);
    }
    else {
      if (d<0) d=1.; //constant trace case
      tmpvec.init(1,1,d);
      tmpvec.sqrt();  
      tmpmat=eigvec;
      tmpmat.scale_cols(tmpvec);
      retval=p->assign_Gram_matrix(tmpmat);
    }
    return retval;
  }


// *****************************************************************************
//                                get_candidate_primal
// *****************************************************************************

int LmaxProblem::get_candidate_primal(PrimalData& pd) const 
  {
    SDPPrimal* p=dynamic_cast<SDPPrimal*>(&pd);
    if (p==0) return 1;
    if (cand_eigs.dim()==0) return 1;
    Integer maxind;
    Real d=max(cand_eigs,&maxind);
    int retval;
    if ((d<0.)&&(trace_stat!=SDPtrace_fixed)){
      tmpmat.init(cand_vecs.rowdim(),1,0.);
      retval=p->assign_Gram_matrix(tmpmat);
    }
    else {
      if (d<0) d=1.; //constant trace case
      tmpvec.init(1,1,d);
      tmpvec.sqrt();  
      tmpmat=cand_vecs;
      tmpmat.scale_cols(tmpvec);
      retval=p->assign_Gram_matrix(tmpmat);
    }
    return retval;
  }


// *****************************************************************************
//                                clear_model
// *****************************************************************************

void LmaxProblem::clear_model()  
{
  problem_modified=1;   
  center_available=0;
  cand_available=0;
  model_available=0;
  model_subg_available=0;
  aug_subg_in_model=0;
  aug_available=0;
  
  eigvec.init(0,0,0.);
  eigval.init(0,0,0.);
  bundlevecs.init(0,0,0.);
  primalvecs.init(0,0,0.);
  primaleigs.init(0,0,0.);
  
  subgCvalues.init(0,0,0.);
  aggrsubgrads.init(0,0,0.);
  aggrcoeff.init(0,0,0.);
  for(unsigned int i=0;i<primal.size();i++){
    delete primal[i];
  }
  primal.clear();

  oracle->get_trace_constraint(lmax_multiplier,trace_stat);
}

// *****************************************************************************
//                                clear_aggregates
// *****************************************************************************

//remove all current aggregate cutting planes
//returns 0 on success, 1 on failure 

int LmaxProblem::clear_aggregates()
{
  if (subgCvalues.dim()>0){
    subgCvalues.init(0,0,0.);
    aggrsubgrads.init(0,0,0.);
    aggrcoeff.init(0,0,0.);
    for(unsigned int i=0;i<primal.size();i++){
      delete primal[i];
    }
    primal.clear();
    
    model_available=0;
    model_subg_available=0;
    aug_subg_in_model=0;
    aug_available=0;
    if (bundlevecs.dim()) {
      //make sure that the subgradient in the center is contained in the model
      tmpvec=Matrix(bundlevecs.coldim(),1,sqrt(lmax_multiplier/bundlevecs.coldim()));
      tmpmat=bundlevecs;
      tmpmat.scale_cols(tmpvec);
      oracle->gramip(tmpmat,subg_val,subg);
      xbpeya(subg,oracle->rhs(),1.,-1.);
      subg_val+=ip(subg,y);
    }
  }
  return 0;  
}



}

