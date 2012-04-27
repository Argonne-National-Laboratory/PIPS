/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/coneproblem.cxx

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
#include "coneproblem.hxx"

using namespace CH_Matrix_Classes;
using namespace CH_Tools;

namespace ConicBundle {

// *****************************************************************************
//                            ConeProblem
// *****************************************************************************

ConeProblem::ConeProblem(BaseConeOracle* fo)
{
  oracle=fo;
  oracle->get_trace_constraint(cone_multiplier,trace_stat);

  dim=oracle->rhs().dim();

  ret_code=0;

  center_x=0;

  aggregtol=0.001;
  maxkeepvecs=50;
  max_new_subg=5;
  n_new_subg=0;
  problem_modified =1;
  center_available=0;
  cand_available=0;
  model_subg_available=0;
  model_changed=0;
  aug_available=0;
  solverp=0;
  out=0;
  print_level=0;

  evaltime=Microseconds(0);
  nr_eval=0;
  
}

// *****************************************************************************
//                            ConeProblem
// *****************************************************************************

ConeProblem::~ConeProblem()
{
  for(unsigned int i=0;i<cand_x.size();i++){
    delete cand_x[i];
  }
  cand_x.clear();
  delete center_x;
  {for(unsigned int i=0;i<bundlex.size();i++){
    delete bundlex[i];
  }}
  bundlex.clear();
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

int ConeProblem::init_center(
	       int& init, Real& lb, Real& ub, Matrix& o_y, Matrix& o_subg,
               Indexmatrix& o_bounds_index, Matrix& o_lby, Matrix& o_uby)
{
  if (!center_available) return 1;
  if ((init==0)||(problem_modified)){
    problem_modified=0;
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

int ConeProblem::eval_function(const Matrix& iny,
			      Real nullstep_bound,Real relprec)
{
  cand_available=0;
  cand_y=iny;
  //delete old primal info
  for(unsigned int i=0;i<cand_x.size();i++){
    delete cand_x[i];
  }
  cand_x.clear();
  

  Real ip_b_iny=ip(oracle->rhs(),cand_y);   //linear term

  Real trace_val;
  Conetrace con_tr;
  oracle->get_trace_constraint(trace_val,con_tr);
  assert(trace_val>=0.);
  if ((con_tr!=trace_stat)||
      ((con_tr!=Conetrace_unbounded)&&(trace_val!=cone_multiplier))||
      ((con_tr==Conetrace_unbounded)&&(trace_val>cone_multiplier))){
    if (out) {
      (*out)<<"**** WARNING: ConeProblem::eval_function(...): trace condition changed without reinitialization of model and function values"<<std::endl;
    }
    cone_multiplier=trace_val;
    trace_stat=con_tr;
  }
			

  //---- evaluate  	
  cand_ublmax=(nullstep_bound-ip_b_iny)/cone_multiplier;
  Real ctx_relprec=relprec/cone_multiplier*(fabs(nullstep_bound)+1.)/(fabs(cand_ublmax)+1.);

  cand_ub_fun_val=nullstep_bound;
  ret_code = oracle->evaluate(cand_y,ctx_relprec,cand_ublmax,
			      cand_ctx,cand_Ax,cand_x);
  nr_eval++;

  if (cand_ctx.dim()==0) {
    if (out) (*out)<<"**** ERROR: ConeProblem::eval_function(): function returned no subgradient values"<<std::endl;
    return 1;
  }
  
  if (cand_ctx.dim()!=cand_Ax.coldim()){
    if (out) (*out)<<"**** ERROR: ConeProblem::eval_function(): function returned different number of values and subgradients"<<std::endl;
    return 1;
  }

  if ((cand_x.size()>0)&&(Integer(cand_x.size())!=cand_Ax.coldim())){
    if (out) (*out)<<"**** ERROR: ConeProblem::eval_function(): function returned different number of primal solutions and subgradients; ignoring primal solutions"<<std::endl;
    for(unsigned int i=0;i<cand_x.size();i++){
      delete cand_x[i];
    }
    cand_x.clear();
  }

  if (trace_stat!=Conetrace_fixed){
    if (cand_ublmax<0.) {
      cand_ublmax=0.;
      cand_lmax=0.;
      cand_subg=oracle->rhs();
    }
    else {
      tmpvec=cand_ctx;
      genmult(cand_Ax,cand_y,tmpvec,-1.,1.,1);
      cand_lmax=max(tmpvec,&cand_subg_maxind);
      if (cand_lmax<0.) {
	cand_lmax=0.;
	cand_subg=oracle->rhs();
      }
      else {
	cand_subg.newsize(oracle->rhs().dim(),1); chk_set_init(cand_subg,1);
	mat_xeyapzb(cand_subg.dim(),cand_subg.get_store(),oracle->rhs().get_store(),cand_Ax.get_store()+cand_subg_maxind*cand_Ax.rowdim(),1.,-cone_multiplier);
      }
    }
  }
  cand_subg_val=cand_lmax*cone_multiplier+ip_b_iny;
  cand_ub_fun_val=cand_ublmax*cone_multiplier+ip_b_iny;
  
  cand_available=1;


  if (ret_code){
    if (out) (*out)<<"**** WARNING: ConeProblem::eval_function(): function returned code = "<<ret_code<<std::endl;
    return ret_code;
  }

  if ((cand_subg_val<nullstep_bound)&&
      (cand_ub_fun_val-cand_subg_val>relprec*(fabs(cand_ub_fun_val)+1.))){
    if (out) (*out)<<"**** WARNING: ConeProblem::eval_function(): insufficient precision in evaluation routine";
    return 1;
  }    
  
  return 0;
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

int ConeProblem::get_function_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!cand_available) return 1;
  lb=cand_subg_val;
  ub=cand_ub_fun_val;
  if (subgp){
    if (cand_ctx.dim()==0){
      return 1;
    }
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

int ConeProblem::do_step(void)
{
  if (!cand_available) return 1;
  y=cand_y;
  if(cand_ctx.dim()>0){
    subg=cand_subg;
    Ax=cand_Ax.col(cand_subg_maxind);
    if (cand_x.size()>0){
      delete center_x;
      center_x=cand_x[cand_subg_maxind]->clone_primal_data();
    }
  }
  else {
    subg.init(0,0,0.);
    delete center_x; center_x=0;
  }
  subg_val=cand_subg_val;
  ub_fun_val=cand_ub_fun_val;
  ublmax=cand_ublmax;
  center_available=1;

  return 0;
}


// *****************************************************************************
//                              eval_model
// *****************************************************************************

//evaluate the current cutting model in $y$ 
//if evaluated by an iterative method that provides upper and lower bounds
//it may stop when the lower bound (lb) is above the nullstep_bound
//evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied

int ConeProblem::eval_model(const Matrix& iny,Real /* eval_model_bound */,Real /* relprec */)
{
  model_available=0;
  model_subg_available=0;
  model_changed=0;
  if (bundlevecs.coldim()==0) {
    if ((cand_available)&&(cand_ctx.dim()>0)){
      if (update_model(false)) return 1;
    }
    else if (center_available){
      bundlecosts=Matrix(1,1,(subg_val-ip(subg,y))/cone_multiplier);
      bundlevecs=oracle->rhs()-subg;
      bundlevecs/=cone_multiplier;
      hashsum=sumrows(bundlevecs);
      if (center_x) bundlex.push_back(center_x->clone_primal_data());
      bundlecoeff.init(1,1,cone_multiplier);
    }
    else return 1;
  }
  
  tmpvec=bundlecosts;
  genmult(bundlevecs,iny,tmpvec,-1.,1.,1);
  Real ip_b_iny=ip(oracle->rhs(),iny); 
  model_subg_val=max(tmpvec,&model_ind)*cone_multiplier+ip_b_iny;
  model_ub_val=model_subg_val;
  
  model_available=1;
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

int ConeProblem::get_model_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!model_available) return 1;
  lb=model_subg_val;
  ub=model_ub_val;
  if (subgp){
    if (!model_subg_available){
      if (model_changed) return 1;
      model_subg.init(bundlevecs.col(model_ind),-cone_multiplier);
      model_subg+=oracle->rhs();
      model_subg_available=1;
    }
    *subgp=model_subg;
  }
  return 0;
}


// *****************************************************************************
//                              get_augmodel_sol
// *****************************************************************************

int ConeProblem::get_augmodel_sol(Real& out_linconst,Matrix& out_subg)
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

int ConeProblem::update_model(bool)
{
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));
  if ((!cand_available)||(cand_ctx.dim()==0)) return 1;
  if (bundlevecs.coldim()==0){
    if (maxkeepvecs>=cand_ctx.dim()){
      bundlecosts=cand_ctx;
      bundlevecs=cand_Ax;
      hashsum=sumrows(bundlevecs);
      bundlex.resize(cand_x.size(),0);
      for (unsigned int i=0;i<cand_x.size();i++){
	bundlex[i]=cand_x[i]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(0)=cone_multiplier;
      model_changed=1;
    }
    else {
      tmpvec=cand_ctx;
      genmult(cand_Ax,cand_y,tmpvec,-1,1.,1);
      sortindex(tmpvec,tmpind);
      tmpind=tmpind(Range(tmpind.dim()-maxkeepvecs,tmpind.dim()-1));
      bundlecosts=cand_ctx(tmpind);
      bundlevecs=cand_Ax.cols(tmpind);
      hashsum=sumrows(bundlevecs);
      bundlex.resize(min(Integer(cand_x.size()),maxkeepvecs),0);
      for (unsigned int i=0;i<bundlex.size();i++){
	bundlex[i]=cand_x[tmpind(i)]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(bundlecoeff.dim()-1)=cone_multiplier;
      model_changed=1;
    }
    return 0;
  }

  //compute value of bundle cutting planes at cand_y
  tmpvec=bundlecosts;
  genmult(bundlevecs,cand_y,tmpvec,-1.,1.,1);
  //append values of new epsilon subgradients 
  //(after a null step one of these must have higher value than bundle)
  Matrix tmpmat=cand_ctx;
  genmult(cand_Ax,cand_y,tmpmat,-1.,1.,1);
  tmpvec.concat_below(tmpmat);
  //sort by this value and use this sorting for generating the next bundle
  //(we must keep at least one from the bundle and 
  // try to eliminate identical cuts)
  sortindex(tmpvec,tmpind);
  Integer sind=tmpind.dim();
  Indexmatrix bundledel;
  Indexmatrix bundlekeep;
  Indexmatrix subgkeep;
  Matrix subghashsum=sumrows(cand_Ax);
  while (--sind>=0){
    Integer ind=tmpind(sind);
    int bunvec=(ind<bundlecosts.dim());
    if ((bundlekeep.dim()+subgkeep.dim()>=maxkeepvecs)||
        ((!bunvec)&&(subgkeep.dim()>=maxkeepvecs-1))
	){//there are already enough new vectors
      if (bunvec) { 
	bundledel.concat_below(ind);
      }
      continue;
    }
    //--- is there a cluster of cutting planes with identical value in cand_y?
    Real lbval=tmpvec(ind);
    lbval-=max(fabs(lbval),1.)*(100.*eps_Real);
    
    if ((sind==0)||(tmpvec(tmpind(sind-1))<=lbval)){ 
      //no cluster, only one element
      if (bunvec){
	bundlekeep.concat_below(ind);
      }
      else {
	subgkeep.concat_below(ind-bundlecosts.dim());
      }
      continue;
    }
	
    //--- find the cluster and split it into bundle vectors and new subgradients
    Integer si=sind;
    Indexmatrix bunind;
    Indexmatrix subgind;
    if (bunvec) bunind.concat_below(ind);
    else subgind.concat_below(ind-bundlecosts.dim());
    while ((--si>=0)&&(tmpvec(tmpind(si))>lbval)){
      Integer i=tmpind(si);
      if (i<bundlecosts.dim()){
	bunind.concat_below(i);
      }
      else {
	subgind.concat_below(i-bundlecosts.dim());
      }	   
    }
    sind=si+1;

    //--- first treat cluster of bundle vectors
    Indexmatrix scoeffind;
    sortindex(bundlecoeff(bunind),scoeffind);
    si=scoeffind.dim();
    Integer newbunstart=bundlekeep.dim();
    while(--si>=0){
      Integer ind=bunind(scoeffind(si));
      int deleteind=0;
      if (bundlekeep.dim()+subgkeep.dim()>=maxkeepvecs) {
	deleteind=1;
      }
      else {
	Real sumval=hashsum(ind);
	for(Integer i=newbunstart;i<bundlekeep.dim();i++){
	  Integer j=bundlekeep(i); 
	  if (fabs(sumval-hashsum(j))>max(fabs(sumval),1.)*(100.*eps_Real))
	    continue;
	  if (norm2(bundlevecs.col(ind)-bundlevecs.col(j))<
	       max(max(abs(bundlevecs.col(ind))),1.)*(1000.*eps_Real)){
            deleteind=1;
	    break;
	  }
	} //end for
      } //end else
      if (deleteind){
	bundledel.concat_below(ind);
      }
      else {
        bundlekeep.concat_below(ind);
      }
    }  
    
    //--- now treat cluster of subgradient vectors
    Integer newsubgstart=subgkeep.dim();
    for(si=0;si<subgind.dim();si++){
      Integer ind=subgind(si);
      int deleteind=0;
      if ((bundlekeep.dim()+subgkeep.dim()>=maxkeepvecs)||
          ((bundlekeep.dim()==0)&&(subgkeep.dim()==maxkeepvecs-1))){
	deleteind=1;
      }
      else {
	Real sumval=subghashsum(ind);
	for(Integer i=newbunstart;i<bundlekeep.dim();i++){
	  Integer j=bundlekeep(i); 
	  if (fabs(sumval-hashsum(j))>max(fabs(sumval),1.)*(100.*eps_Real))
	    continue;
	  if (norm2(cand_Ax.col(ind)-bundlevecs.col(j))<
	       max(max(abs(cand_Ax.col(ind))),1.)*(1000.*eps_Real)){
            deleteind=1;
	    break;
	  }
	} //end for
	{for(Integer i=newsubgstart;(i<subgkeep.dim())&&(!deleteind);i++){
	  Integer j=subgkeep(i); 
	  if (fabs(sumval-subghashsum(j))>max(fabs(sumval),1.)*(100.*eps_Real))
	      continue;
	  if (norm2(cand_Ax.col(ind)-cand_Ax.col(j))<
	      max(max(abs(cand_Ax.col(ind))),1.)*(1000.*eps_Real))
	    {
            deleteind=1;
	    break;
	  }
	}} //end for
      } //end else
      if (!deleteind){
	subgkeep.concat_below(ind);
      }
    }//end for cluster of new subgradients

  }//end while(cutting planes)

  assert(bundlekeep.dim()+bundledel.dim()==bundlevecs.coldim());

  if ((out)&&(print_level>=2)){
    out->precision(8);
  (*out)<<" Cone: bundle_update: keep="<<bundlekeep.dim()<<"("<<max(tmpvec(bundlekeep))<<","<<min(tmpvec(bundlekeep))<<")";
  (*out)<<" del="<<bundledel.dim()<<"("<<max(tmpvec(bundledel))<<","<<min(tmpvec(bundledel))<<")";
  (*out)<<" new="<<subgkeep.dim()<<"("<<max(tmpmat(subgkeep))<<","<<min(tmpmat(subgkeep))<<")";
  (*out)<<" ignore="<<cand_ctx.dim()-subgkeep.dim()<<"("<<min(tmpmat)<<")"<<std::endl;
  }
    
  //--- aggregate deleted bundle columns into smallest bundle column kept
  if (bundledel.dim()>0){
    model_changed=1;
    Integer aggregind;
    Real maxcoeff=max(bundlecoeff);
    min(bundlecoeff(bundlekeep),&aggregind);
    aggregind=bundlekeep(aggregind);
    for(Integer i=0;i<bundledel.dim();i++){
      Integer delind=bundledel(i);
      Real a=bundlecoeff(delind);
      Real b=bundlecoeff(aggregind);
      Real s=a+b;
      if (s>max(1.,maxcoeff)*1e-12) {
	mat_xbpeya(bundlevecs.rowdim(),
		   bundlevecs.get_store()+aggregind*bundlevecs.rowdim(),
		   bundlevecs.get_store()+delind*bundlevecs.rowdim(),
		   a/s,b/s);
      }
      if (bundlex.size()>0){
	if (s>max(1.,maxcoeff)*1e-12){
	  bundlex[aggregind]->aggregate_primal_data(b/s,a/s,*bundlex[delind]);
	}
	delete bundlex[delind];
	bundlex[delind]=0;
      }
      if (s>max(1.,maxcoeff)*1e-12){
	bundlecosts(aggregind)=(bundlecosts(delind)*a+bundlecosts(aggregind)*b)/s;
      }
      bundlecoeff(aggregind)=s;
    }
    bundlevecs.delete_cols(bundledel);
    if (bundlex.size()>0){
      sortindex(bundledel,tmpind);
      Integer cnt=bundledel(tmpind(0));
      Integer delcnt=0;
      Integer nextdel;
      if (++delcnt<bundledel.dim()){
	nextdel=bundledel(tmpind(delcnt));
      }
      else {
	nextdel=Integer(bundlex.size());
      }
      for(Integer i=cnt+1;i<Integer(bundlex.size());i++){
	if (i==nextdel){
	  if (++delcnt<bundledel.dim()){
	    nextdel=bundledel(tmpind(delcnt));
	  }
	  else {
	    nextdel=Integer(bundlex.size());
	  }
	  continue;
	} 
	bundlex[cnt]=bundlex[i];
	cnt++;
      }
      bundlex.resize(cnt);
    }
    bundlecosts.delete_rows(bundledel);
    bundlecoeff.delete_rows(bundledel);
    hashsum.delete_cols(bundledel);
  }

  if ((out)&&(print_level>=1)){
    out->precision(8);
    (*out)<<"    max(bundlecoeff)="<<max(bundlecoeff);
    (*out)<<"    min(bundlecoeff)="<<min(bundlecoeff)<<std::endl;
  }

  //--- add new subgradient columns to bundle
  n_new_subg=subgkeep.dim();
  if (subgkeep.dim()>0){
    model_changed=1;
    if (subgkeep.dim()==cand_Ax.coldim()){
      if ((Integer(bundlex.size())==bundlevecs.coldim())&&(cand_x.size()>0)){
	Integer olddim=Integer(bundlex.size());
	bundlex.resize(olddim+subgkeep.dim(),0);
	for(Integer i=0;i<subgkeep.dim();i++){
	  bundlex[olddim+i]=cand_x[i]->clone_primal_data();
	}
      }
      bundlevecs.concat_right(cand_Ax);
      bundlecosts.concat_below(cand_ctx);
      bundlecoeff.concat_below(Matrix(cand_ctx.dim(),1,0.));
      hashsum.concat_right(subghashsum);
    }
    else {
      if ((Integer(bundlex.size())==bundlevecs.coldim())&&(cand_x.size()>0)){
	Integer olddim=Integer(bundlex.size());
	bundlex.resize(olddim+subgkeep.dim(),0);
	for(Integer i=0;i<subgkeep.dim();i++){
	  bundlex[olddim+i]=cand_x[subgkeep(i)]->clone_primal_data();
	}
      }
      for(Integer i=0;i<subgkeep.dim();i++){	
	Integer ind=subgkeep(i);
	tmpvec.init(cand_Ax.col(ind));
	bundlevecs.concat_right(tmpvec);
        bundlecosts.concat_below(cand_ctx(ind));
	bundlecoeff.concat_below(0.);
        hashsum.concat_right(subghashsum(ind));
      }
    }
  }    
     
  /*
  sortindex(bundlecoeff,tmpind);
  Real aggregval=aggregtol*bundlecoeff(tmpind(bundlecoeff.dim()-1));
  if ((bundlevecs.coldim()>=maxkeepvecs)||(bundlecoeff(tmpind(0))<aggregval)){
    Integer ind=max(bundlevecs.coldim()-maxkeepvecs+1,Integer(1));
    while(bundlecoeff(tmpind(ind))<aggregval) ind++;
    Integer aggregind=tmpind(ind);
    tmpind.reduce_length(ind);
    for(Integer i=0;i<ind;i++){
      Real a=bundlecoeff(tmpind(i));
      Real b=bundlecoeff(aggregind);
      Real s=a+b;
      mat_xbpeya(bundlevecs.rowdim(),
		 bundlevecs.get_store()+aggregind*bundlevecs.rowdim(),
		 bundlevecs.get_store()+tmpind(i)*bundlevecs.rowdim(),
		 a/s,b/s);
      if (primdim>0){
	mat_xbpeya(bundlex.rowdim(),
		   bundlex.get_store()+aggregind*bundlex.rowdim(),
		   bundlex.get_store()+tmpind(i)*bundlex.rowdim(),
		   a/s,b/s);
      }
      bundlecosts(aggregind)=(bundlecosts(tmpind(i))*a+bundlecosts(aggregind)*b)/s;
      bundlecoeff(aggregind)=s;
    }
    bundlevecs.delete_cols(tmpind);
    if (primdim>0){
      bundlex.delete_cols(tmpind);
    }
    bundlecosts.delete_rows(tmpind);
    bundlecoeff.delete_rows(tmpind);
  }
  bundlevecs.concat_right(cand_Ax);
  if (primdim>0){
    bundlex.concat_right(cand_x);
  }
  bundlecosts.concat_below(cand_ctx);
  bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
  model_changed=1;
  */
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));

  return 0;
}
  
// *****************************************************************************
//                              intersect_box
// *****************************************************************************

int ConeProblem::intersect_box(Indexmatrix& inbindex,Matrix& inlb,Matrix& inub)
{
  chk_init(inbindex);
  chk_init(inlb);
  chk_init(inub);
  if ((inbindex.dim()==0)&&(inlb.dim()==0)&&(inub.dim()==0)){
    inbindex=bounds_index;
    if (lby.dim()==0) inlb.init(dim,1,CB_minus_infinity);
    else inlb=lby;
    if (uby.dim()==0) inub.init(dim,1,CB_plus_infinity);
    else inub=uby;
    return 0;
  }
  if ((lby.dim()==0)&&(uby.dim()==0)){
    return 0;
  }
  assert(inub.dim()==dim);
  assert(inlb.dim()==dim);
  tmpind(min(inbindex.dim()+bounds_index.dim(),dim),1);
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
  else {
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

Real ConeProblem::lb_function(const Matrix& iny)
{
  if (aug_available){
    return aug_linconst+ip(aug_subg,iny);
  }
  if (cand_available){
    tmpvec=cand_ctx;
    genmult(cand_Ax,iny,tmpvec,-1.,1.,1);
    return max(tmpvec)*cone_multiplier+ip(oracle->rhs(),iny);
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

Real ConeProblem::lb_model(const Matrix& iny)
{
  if (aug_available){
    return aug_linconst+ip(aug_subg,iny);
  }
  if (cand_available){
    tmpvec=cand_ctx;
    genmult(cand_Ax,iny,tmpvec,-1.,1.,1);
    return max(tmpvec)*cone_multiplier+ip(oracle->rhs(),iny);
  }
  return CB_minus_infinity;
}


// *****************************************************************************
//                                start_augmodel
// *****************************************************************************

//return a pointer to the variables/constraints generating the cutting model
//returns 0 on success, 1 on failure

int ConeProblem::start_augmodel(QP_Block*& blockp)
{
  if (bundlevecs.coldim()==0) {
    if ((cand_available)&&(cand_ctx.dim()>0)){
      if (update_model(false)) return 1;
    }
    else if (center_available){
      bundlecosts=Matrix(1,1,(subg_val-ip(subg,y))/cone_multiplier);
      bundlevecs=oracle->rhs()-subg;
      bundlevecs/=cone_multiplier;
      hashsum=sumrows(bundlevecs);
      if (center_x) 
	bundlex.push_back(center_x->clone_primal_data());
      bundlecoeff.init(1,1,cone_multiplier);
    }
    else return 1;
  }
  //--- set blocks
  //initialize blocks and variable dimensions
  tmpind.init(0,0,Integer(0));
  aug_available=0;
  block.init_block(bundlevecs.coldim(),tmpind,tmpind,cone_multiplier,trace_stat!=Conetrace_fixed);

  blockp=&block;

  return 0;
}
 
// *****************************************************************************
//                                  get_row
// *****************************************************************************

//store the coefficients corresponding to coordinate index_y of y
//for all model variables in the positions 
//row(startindex,...,startindex+model_dim-1)
//add rhs to bi
//index_y==-1 is used for the cost coefficients

int ConeProblem::get_row(Integer index_y,Matrix& row,Real& bi,Integer startindex) const
{
  if (index_y==-1){
    mat_xey(bundlecosts.dim(),row.get_store()+startindex,bundlecosts.get_store());
  }
  else {
    bi+=(oracle->rhs())(index_y);
    mat_xemy(bundlevecs.coldim(),
	    row.get_store()+startindex,1,
	    bundlevecs.get_store()+index_y,bundlevecs.rowdim());
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

int ConeProblem::make_aug_linmodel(Real* in_aug_linconst,Matrix* in_aug_subg,bool* conemult_increased, Real* function_value)  
{

  block.get_linx(bundlecoeff);
  aug_linconst=ip(bundlecosts,bundlecoeff);
  aug_subg=oracle->rhs();
  genmult(bundlevecs,bundlecoeff,aug_subg,-1.,1.);
  
  if (in_aug_linconst) *in_aug_linconst+=aug_linconst;
  if (in_aug_subg) *in_aug_subg+=aug_subg;

  //--- check whether it is necessary to increase the cone multiplier
#ifdef TRACE_OUTPUT
  if (out) (*out)<<"\n ConeProblem: trace="<<sum(bundlecoeff)<<" mult="<<cone_multiplier<<std::endl;
#endif

  if ((conemult_increased!=0)&&(trace_stat==Conetrace_unbounded)&&(sum(bundlecoeff)>.95*cone_multiplier)){
    cone_multiplier=max(1.,cone_multiplier*2.);
#ifdef TRACE_OUTPUT
    if (out) (*out)<<"              multiplier increased to "<<cone_multiplier<<std::endl;
#endif
    *conemult_increased=true;
    block.adjust_trace(cone_multiplier);
    if (ublmax>0.){
      Real ip_b_iny=ip(oracle->rhs(),y);
      subg_val=ip_b_iny+lmax*cone_multiplier;
      ub_fun_val=ip_b_iny+ublmax*cone_multiplier;
      if (lmax!=0.){
	subg.init(Ax,-cone_multiplier);
	subg+=oracle->rhs();
      }
      if (function_value==0) return 2;
      *function_value=ub_fun_val;
    }
  }

  aug_available=1;
 
  return 0;
} 

// *****************************************************************************
//                                adjust_multiplier
// *****************************************************************************

int ConeProblem::adjust_multiplier()  
{
  if (trace_stat!=Conetrace_unbounded) return 0;
  if (!center_available)  return 1;
  if (block.get_linx(tmpvec)) return 1;
#ifdef TRACE_OUTPUT
  if (out) (*out)<<"\n ConeProblem: trace="<<sum(tmpvec)<<" mult="<<cone_multiplier;
#endif
  cone_multiplier=min(cone_multiplier,max(1.,2.*sum(tmpvec)));
  oracle->adjust_multiplier(cone_multiplier);
#ifdef TRACE_OUTPUT
  if (out) (*out)<<" newmult="<<cone_multiplier<<std::endl;
#endif
  block.adjust_trace(cone_multiplier);
  if (ublmax>0.){
    Real ip_b_iny=ip(oracle->rhs(),y);
    subg_val=ip_b_iny+lmax*cone_multiplier;
    ub_fun_val=ip_b_iny+ublmax*cone_multiplier;
    if (lmax!=0.){
      subg.init(Ax,-cone_multiplier);
      subg+=oracle->rhs();
    }
  }
  return 0;
}
 
// *****************************************************************************
//                               change_variables
// *****************************************************************************

//change variables as described in ChangeVariableInfo
//returns 0 on success, 1 on failure
//(this is not needed in the bundle framework, routine may return 1!)

int ConeProblem::change_variables(ChangeVarInfo* cvp)
{
  if (typeid(*cvp)==typeid(AppendVars)){
    AppendVars* cp=dynamic_cast<AppendVars*>(cvp);
    Integer n_append=cp->n_append;
    if (cp->boundsi) bounds_index.concat_below(*(cp->boundsi));
    if (cp->lb) {
      if (lby.dim()==0) lby.init(dim,1,CB_minus_infinity);
      lby.concat_below(*(cp->lb));
    }
    else {
      if (lby.dim()!=0) lby.concat_below(Matrix(n_append,1,CB_minus_infinity));
    }
    if (cp->ub) {
      if (uby.dim()==0) uby.init(dim,1,CB_plus_infinity);
      uby.concat_below(*(cp->ub));
    }
    else {
      if (uby.dim()!=0) uby.concat_below(Matrix(n_append,1,CB_plus_infinity));
    }
    if (cp->startval) {
      if (center_available){
	for(Integer i=0;i<n_append;i++){
	  if ((*cp->startval)(i)!=0.){
	    center_available=0;
	    delete center_x; center_x=0;
	    break;
	  }
	}
      } 
      y.concat_below(*(cp->startval));
    }
    else y.concat_below(Matrix(n_append,1,0.));
    dim+=n_append;
    
    
    //try to extend bundlevectors
    tmpvec.newsize(n_append,bundlevecs.coldim());
    tmpvec.init(0,0,0.);
    Indexmatrix ivec(Range(dim-n_append,dim-1));
    Matrix dvec(n_append,1); chk_set_init(dvec,1);
    int no_subg_extension=0;
    for(Integer i=0;i<bundlevecs.coldim();i++){
      PrimalData* pd=0;
      if (i<Integer(bundlex.size())) { 
      pd=bundlex[i];
      }
      if (oracle->subgradient_extension(pd,ivec,dvec)){
	if (out) (*out)<<"**** WARNING: ConeProblem::add_variable(): subgradient extension failed"<<std::endl;
	no_subg_extension=1;
	break;
      }
      if (dvec.dim()!=n_append){
	if (out) (*out)<<"**** WARNING: ConeProblem::add_variable(): oracle.subgradient_extension(...) returns a vector of wrong dimension, subgradient extension failed"<<std::endl;
	no_subg_extension=1;
	break;
      }
      tmpvec.concat_right(dvec);
      for(Integer j=0;j<n_append;j++) {
	tmpvec(j,i)=dvec[j];
      }
    }
    if (no_subg_extension) {
      clear_model();
    }
    else {
      bundlevecs.concat_below(tmpvec);
      hashsum+=sumrows(tmpvec);
      
      if (center_available){
	//make sure that the center subgradient is contained in the model
	if (bundlecoeff.dim()>0) {
	  subg=oracle->rhs();
	  genmult(bundlevecs,bundlecoeff,subg,-cone_multiplier,1.);
	  subg_val=ip(bundlecosts,bundlecoeff)*cone_multiplier+ip(subg,y);
	  delete center_x;center_x=0;
	  if (bundlex.size()>0) {
	    center_x=bundlex[0]->clone_primal_data();
	    if (bundlecoeff.dim()>1){
	      center_x->aggregate_primal_data(bundlecoeff[0],bundlecoeff[1],*bundlex[1]);
	    }
	    for(Integer i=2;i<bundlecoeff.dim();i++){
	      center_x->aggregate_primal_data(1.,bundlecoeff[i],*bundlex[i]);
	    }
	  }
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
    }
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
      if (center_available){ 
	subg=subg(cp->assign_ind);
      }
    }
      
    dim=cp->assign_ind.dim();
    
    if (lby.dim()!=0) lby=lby(cp->assign_ind);
    if (uby.dim()!=0) uby=uby(cp->assign_ind);
    if (bounds_index.dim()>0){
      bounds_index.init(0,0,Integer(0));
      for(Integer i=0;i<dim;i++){
	if (((lby.dim()!=0)&&(lby(i)>CB_minus_infinity))||
	    ((uby.dim()!=0)&&(uby(i)<CB_plus_infinity))){
	  bounds_index.concat_below(i);
	}
      }
    }
    y=y(cp->assign_ind);
    bundlevecs=bundlevecs.rows(cp->assign_ind);
    hashsum=sumrows(bundlevecs);
    
    if (center_available) {
      //make sure that the center subgradient is contained in the model
      if (bundlecoeff.dim()>0){
	subg=oracle->rhs();
	genmult(bundlevecs,bundlecoeff,subg,-cone_multiplier,1.);
	subg_val=ip(bundlecosts,bundlecoeff)*cone_multiplier+ip(subg,y);
	delete center_x;center_x=0;
	if (bundlex.size()>0) {
	  center_x=bundlex[0]->clone_primal_data();
	  if (bundlecoeff.dim()>1){
	    center_x->aggregate_primal_data(bundlecoeff[0],bundlecoeff[1],*bundlex[1]);
	  }
	  for(Integer i=2;i<bundlecoeff.dim();i++){
	    center_x->aggregate_primal_data(1.,bundlecoeff[i],*bundlex[i]);
	  }
	}
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

  //--- delete variables

  if (typeid(*cvp)==typeid(DeleteVars)) {
    DeleteVars* cp=dynamic_cast<DeleteVars*>(cvp);
    if (cp->del_index.dim()==0) return 0;
    if (bounds_index.dim()>0){
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
    
    dim-=cp->del_index.dim();
    if (lby.dim()!=0) lby.delete_rows(cp->del_index);
    if (uby.dim()!=0) uby.delete_rows(cp->del_index);
    y.delete_rows(cp->del_index);
    bundlevecs.delete_rows(cp->del_index);
    hashsum=sumrows(bundlevecs);
    
    if (center_available) {
      //make sure that the center subgradient is contained in the model
      if (bundlecoeff.dim()>0){
	subg=oracle->rhs();
	subg.delete_rows(cp->del_index);
	genmult(bundlevecs,bundlecoeff,subg,-cone_multiplier,1.);
	subg_val=ip(bundlecosts,bundlecoeff)*cone_multiplier+ip(subg,y);
	delete center_x;center_x=0;
	if (bundlex.size()>0) {
	  center_x=bundlex[0]->clone_primal_data();
	  if (bundlecoeff.dim()>1){
	    center_x->aggregate_primal_data(bundlecoeff[0],bundlecoeff[1],*bundlex[1]);
	  }
	  for(Integer i=2;i<bundlecoeff.dim();i++){
	    center_x->aggregate_primal_data(1.,bundlecoeff[i],*bundlex[i]);
	  }
	}
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

int ConeProblem::recompute_center()
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

int ConeProblem::get_approximate_primal(PrimalData& primal) const 
  {
    if (!aug_available) return -1;   
    if (bundlex.size()==0) return 1; 

    int retval = primal.assign_primal_data(*bundlex[0]);
    if (retval) return retval;
    if (bundlecoeff.dim()>1){
      retval=primal.aggregate_primal_data(bundlecoeff[0],bundlecoeff[1],*bundlex[1]);
      if (retval) return retval;
    }
    for(Integer i=2;i<bundlecoeff.dim();i++){
      retval=primal.aggregate_primal_data(1.,bundlecoeff[i],*bundlex[i]);
      if (retval) return retval;
    }
    return 0;
  }


// *****************************************************************************
//                                clear_model
// *****************************************************************************

void ConeProblem::clear_model()  
{
  problem_modified=1;   
  center_available=0;
  cand_available=0;
  model_available=0;
  model_subg_available=0;
  aug_available=0;
  
  bundlecosts.init(0,0,0.);
  bundlevecs.init(0,0,0.); 
  bundlecoeff.init(0,0,0.);
  hashsum.init(0,0,0.);

  for(unsigned int i=0;i<cand_x.size();i++){
    delete cand_x[i];
  }
  cand_x.clear();
  delete center_x;
  center_x=0;
  {for(unsigned int i=0;i<bundlex.size();i++){
    delete bundlex[i];
  }}
  bundlex.clear();
}

} 

