/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/funproblem.cxx

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
#include "funproblem.hxx"

using namespace CH_Tools;
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                            FunctionProblem
// *****************************************************************************

FunctionProblem::FunctionProblem(Integer indim, MatrixFunctionOracle& fo):
  oracle(fo)
{
  dim=indim;

  ret_code=0;

  x=0;

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
  update_rule=0;  
}

// *****************************************************************************
//                            FunctionProblem
// *****************************************************************************

FunctionProblem::~FunctionProblem()
{
  for(unsigned int i=0;i<cand_x.size();i++){
    delete cand_x[i];
  }
  cand_x.clear();
  delete x;
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

int FunctionProblem::init_center(
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
    if (lby.coldim()==0) 
      o_lby.init(y.rowdim(),y.coldim(),CB_minus_infinity);
    else o_lby=lby;
    if (uby.coldim()==0) 
      o_uby.init(y.rowdim(),y.coldim(),CB_plus_infinity);
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

int FunctionProblem::eval_function(const Matrix& iny,
			      Real nullstep_bound,Real relprec)
{
  cand_available=0;
  cand_y=iny;
  //delete old primal info
  for(unsigned int i=0;i<cand_x.size();i++){
    delete cand_x[i];
  }
  cand_x.clear();
  cand_subg_valvec.init(0,0,0.);
  cand_subg.init(0,0,0.);
  PrimalExtender* pep=0;			

  //---- evaluate  	
  cand_ub_fun_val=nullstep_bound;
  ret_code = oracle.evaluate(cand_y,relprec,cand_ub_fun_val,
			     cand_subg_valvec,cand_subg,cand_x,pep);
  nr_eval++;

  if (cand_subg_valvec.dim()==0) {
    if (out) (*out)<<"**** ERROR: FunctionProblem::eval_function(): function returned no subgradients"<<std::endl;
    return 1;
  }

  if (cand_subg_valvec.dim()!=cand_subg.coldim()){
    if (out) (*out)<<"**** ERROR: FunctionProblem::eval_function(): function returned different number of values and subgradients"<<std::endl;
    return 1;
  }

  if (cand_y.dim()!=cand_subg.rowdim()){
    if (out) (*out)<<"**** ERROR: FunctionProblem::eval_function(): function returned subgradients of wrong dimension"<<std::endl;
    return 1;
  }

  if (cand_ub_fun_val<max(cand_subg_valvec)){
    if (out) (*out)<<"**** ERROR: FunctionProblem::eval_function(): function returned an upper bound="<<cand_ub_fun_val<<" smaller than the lower bound="<<max(cand_subg_valvec)<<" on the function value"<<std::endl;
    return 1;
  }

  if ((cand_x.size()>0)&&(Integer(cand_x.size())!=cand_subg.coldim())){
    if (out) (*out)<<"**** ERROR: FunctionProblem::eval_function(): function returned different number of primal solutions and subgradients; ignoring primal solutions"<<std::endl;
    for(unsigned int i=0;i<cand_x.size();i++){
      delete cand_x[i];
    }
    cand_x.clear();
  }
  
  cand_subg_val=max(cand_subg_valvec,&cand_subg_maxind);
  
  cand_available=1;

  if (pep){
    if (x){
      if (pep->extend(*x)){
	if (out) (*out)<<"**** WARNING: FunctionProblem::eval_function(): PrimalExtender::extend failed"<<std::endl;
      }
    }
    for(unsigned int i=0;i<bundlex.size();i++){
      if (pep->extend(*(bundlex[i]))){
	if (out) (*out)<<"**** WARNING: FunctionProblem::eval_function(): PrimalExtender::extend failed"<<std::endl;
      }
    }
    delete pep;
  }

  if (ret_code){
    if (out) (*out)<<"**** WARNING: FunctionProblem::eval_function(): function returned code = "<<ret_code<<std::endl;
    return ret_code;
  }

  if ((cand_subg_val<nullstep_bound)&&
      (cand_ub_fun_val-cand_subg_val>relprec*(fabs(cand_ub_fun_val)+1.))){
    if (out) (*out)<<"**** WARNING: FunctionProblem::eval_function(): insufficient precision in evaluation routine"<<std::endl;
    return 1;
  }    

  if (cand_ub_fun_val>=CB_plus_infinity){
    if (out) (*out)<<"**** WARNING: FunctionProblem::eval_function(): cand_ub_fun_val>=CB_plus_infinity"<<std::endl;
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

int FunctionProblem::get_function_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!cand_available) return 1;
  lb=cand_subg_val;
  ub=cand_ub_fun_val;
  if (subgp){
    if (cand_subg_valvec.dim()==0){
      return 1;
    }
    *subgp=cand_subg.col(cand_subg_maxind);
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

int FunctionProblem::do_step(void)
{
  if (!cand_available) return 1;
  y=cand_y;
  if(cand_subg_valvec.dim()>0){
    subg=cand_subg.col(cand_subg_maxind);
    if (cand_x.size()>0){
      delete x;
      x=cand_x[cand_subg_maxind]->clone_primal_data();
    }
  }
  else {
    subg.init(0,0,0.);
    delete x; x=0;
  }
  subg_val=cand_subg_val;
  ub_fun_val=cand_ub_fun_val;
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

int FunctionProblem::eval_model(const Matrix& iny,Real /* eval_model_bound */,Real /* relprec */)
{
  model_available=0;
  model_subg_available=0;
  model_changed=0;
  if (bundlevecs.coldim()==0) {
    if ((cand_available)&&(cand_subg_valvec.dim()>0)){
      if (update_model(false)) return 1;
    }
    else if (center_available){
      bundlecosts=Matrix(1,1,subg_val-ip(subg,y));
      bundlevecs=subg;
      hashsum=sumrows(bundlevecs);
      if (x) bundlex.push_back(x->clone_primal_data());
      bundlecoeff.init(1,1,1.);
    }
    else return 1;
  }
  
  tmpvec=bundlecosts;
  if (iny.dim()>0){
    genmult(bundlevecs,iny,tmpvec,1.,1.,1);
  }
  model_subg_val=max(tmpvec,&model_ind);
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

int FunctionProblem::get_model_sol(Real& lb, Real& ub, Matrix* subgp)
{
  if (!model_available) return 1;
  lb=model_subg_val;
  ub=model_ub_val;
  if (subgp){
    if (!model_subg_available){
      if (model_changed) return 1;
      model_subg=bundlevecs.col(model_ind);
      model_subg_available=1;
    }
    *subgp=model_subg;
  }
  return 0;
}


// *****************************************************************************
//                              get_augmodel_sol
// *****************************************************************************

int FunctionProblem::get_augmodel_sol(Real& out_linconst,Matrix& out_subg)
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

int FunctionProblem::update_model(bool descent_step)
{
  switch(update_rule){
    case 0: return update_model0(descent_step);
    case 1: return update_model1(descent_step); 
    case 2: return update_model2(descent_step); 
    case 3: return update_model3(descent_step);  
    case 4: return update_model4(descent_step);   
    case 5: return update_model5(descent_step); 
    default: return update_model0(descent_step);
  }
  return -1;
}



int FunctionProblem::update_model0(bool)
{
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));
  if ((!cand_available)||(cand_subg_valvec.dim()==0)) return 1;
  if (bundlevecs.coldim()==0){
    if (maxkeepvecs>=cand_subg_valvec.dim()){
      bundlecosts=cand_subg_valvec;
      if (cand_y.dim()>0){
	genmult(cand_subg,cand_y,bundlecosts,-1.,1.,1);
      }
      bundlevecs=cand_subg;
      hashsum=sumrows(bundlevecs);
      bundlex.resize(cand_x.size(),0);
      for (unsigned int i=0;i<cand_x.size();i++){
	bundlex[i]=cand_x[i]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      Integer maxind;
      max(cand_subg_valvec,&maxind);
      bundlecoeff(maxind)=1.;
      model_changed=1;
    }
    else {
      sortindex(cand_subg_valvec,tmpind);
      tmpind=tmpind(Range(tmpind.dim()-maxkeepvecs,tmpind.dim()-1));
      bundlecosts=cand_subg_valvec(tmpind);
      bundlevecs=cand_subg.cols(tmpind);
      if (cand_y.dim()>0){
	genmult(bundlevecs,cand_y,bundlecosts,-1.,1.,1);
      }
      hashsum=sumrows(bundlevecs);
      bundlex.resize(min(Integer(cand_x.size()),maxkeepvecs),0);
      for (unsigned int i=0;i<bundlex.size();i++){
	bundlex[i]=cand_x[tmpind(i)]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(bundlecoeff.dim()-1)=1.;
      model_changed=1;
    }
    return 0;
  }

  //compute value of bundle cutting planes at cand_x
  tmpvec=bundlecosts;
  if (cand_y.dim()>0){
    genmult(bundlevecs,cand_y,tmpvec,1.,1.,1);
  }
  //append values of new epsilon subgradients 
  //(after a null step one of these must have higher value than bundle)
  tmpvec.concat_below(cand_subg_valvec);
  //sort by this value and use this sorting for generating the next bundle
  //(we must keep at least one from the bundle and 
  // try to eliminate identical cuts)
  sortindex(tmpvec,tmpind);
  Integer sind=tmpind.dim();
  Indexmatrix bundledel;
  Indexmatrix bundlekeep;
  Indexmatrix subgkeep;
  Matrix subghashsum=sumrows(cand_subg);
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
	  if (norm2(cand_subg.col(ind)-bundlevecs.col(j))<
	       max(max(abs(cand_subg.col(ind))),1.)*(1000.*eps_Real)){
            deleteind=1;
	    break;
	  }
	} //end for
	{for(Integer i=newsubgstart;(i<subgkeep.dim())&&(!deleteind);i++){
	  Integer j=subgkeep(i); 
	  if (fabs(sumval-subghashsum(j))>max(fabs(sumval),1.)*(100.*eps_Real))
	      continue;
	  if (norm2(cand_subg.col(ind)-cand_subg.col(j))<
	      max(max(abs(cand_subg.col(ind))),1.)*(1000.*eps_Real))
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
  (*out)<<" bundle_update: keep="<<bundlekeep.dim()<<"("<<max(tmpvec(bundlekeep))<<","<<min(tmpvec(bundlekeep))<<")";
  (*out)<<" del="<<bundledel.dim()<<"("<<max(tmpvec(bundledel))<<","<<min(tmpvec(bundledel))<<")";
  (*out)<<" new="<<subgkeep.dim()<<"("<<max(cand_subg_valvec(subgkeep))<<","<<min(cand_subg_valvec(subgkeep))<<")";
  (*out)<<" ignore="<<cand_subg_valvec.dim()-subgkeep.dim()<<"("<<min(cand_subg_valvec)<<")"<<std::endl;
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
      for(unsigned int i=(unsigned int)(cnt)+1;i<bundlex.size();i++){
	if (Integer(i)==nextdel){
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

  //--- add new subgradient columns to bundle
  n_new_subg=subgkeep.dim();
  if (subgkeep.dim()>0){
    model_changed=1;
    if (subgkeep.dim()==cand_subg.coldim()){
      if ((Integer(bundlex.size())==bundlevecs.coldim())&&(cand_x.size()>0)){
	Integer olddim=Integer(bundlex.size());
	bundlex.resize(olddim+subgkeep.dim(),0);
	for(Integer i=0;i<subgkeep.dim();i++){
	  bundlex[olddim+i]=cand_x[i]->clone_primal_data();
	}
      }
      bundlevecs.concat_right(cand_subg);
      tmpvec=cand_subg_valvec;
      if (cand_y.dim()>0){
	genmult(cand_subg,cand_y,tmpvec,-1.,1.,1);
      }
      bundlecosts.concat_below(tmpvec);
      bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
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
	tmpvec.init(cand_subg.col(ind));
	bundlevecs.concat_right(tmpvec);
        bundlecosts.concat_below(cand_subg_valvec(ind)-ip(tmpvec,cand_y));
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
  bundlevecs.concat_right(cand_subg);
  if (primdim>0){
    bundlex.concat_right(cand_x);
  }
  tmpvec=cand_subg_valvec;
  genmult(cand_subg,cand_y,tmpvec,-1.,1.,1);
  bundlecosts.concat_below(tmpvec);
  bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
  model_changed=1;
  */
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));

  return 0;
}
  
int FunctionProblem::update_model1(bool)
{
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));
  if ((!cand_available)||(cand_subg_valvec.dim()==0)) return 1;
  if (bundlevecs.coldim()==0){
    if (maxkeepvecs>=cand_subg_valvec.dim()){
      bundlecosts=cand_subg_valvec;
      if (cand_y.dim()>0){
	genmult(cand_subg,cand_y,bundlecosts,-1.,1.,1);
      }
      bundlevecs=cand_subg;
      hashsum=sumrows(bundlevecs);
      bundlex.resize(cand_x.size(),0);
      for (unsigned int i=0;i<cand_x.size();i++){
	bundlex[i]=cand_x[i]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      Integer maxind;
      max(cand_subg_valvec,&maxind);
      bundlecoeff(maxind)=1.;
      model_changed=1;
    }
    else {
      sortindex(cand_subg_valvec,tmpind);
      tmpind=tmpind(Range(tmpind.dim()-maxkeepvecs,tmpind.dim()-1));
      bundlecosts=cand_subg_valvec(tmpind);
      bundlevecs=cand_subg.cols(tmpind);
      if (cand_y.dim()>0){
	genmult(bundlevecs,cand_y,bundlecosts,-1.,1.,1);
      }
      hashsum=sumrows(bundlevecs);
      bundlex.resize(min(Integer(cand_x.size()),maxkeepvecs),0);
      for (unsigned int i=0;i<bundlex.size();i++){
	bundlex[i]=cand_x[tmpind(i)]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(bundlecoeff.dim()-1)=1.;
      model_changed=1;
    }
    return 0;
  }

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
      bundlecosts(aggregind)=(bundlecosts(tmpind(i))*a+bundlecosts(aggregind)*b)/s;
      bundlecoeff(aggregind)=s;
      if (bundlex.size()>0){
	bundlex[aggregind]->aggregate_primal_data(b/s,a/s,*bundlex[tmpind(i)]);
      }
    }
    bundlevecs.delete_cols(tmpind);
    bundlecosts.delete_rows(tmpind);
    bundlecoeff.delete_rows(tmpind);
    if ((bundlex.size()>0)&&(tmpind.dim()>0)){
      Indexmatrix sind;
      sortindex(tmpind,sind);
      tmpind=tmpind(sind);
      Integer j=tmpind(0);
      delete bundlex[j];
      Integer ti=1;
      Integer i=j+1;
      for(;(i<Integer(bundlex.size()))&&(ti<tmpind.dim());i++){
        if (tmpind(ti)==i) {
	  delete bundlex[i];
	  ti++;
	  continue;
	}
	bundlex[j]=bundlex[i];
        j++;
      }
      for(;i<Integer(bundlex.size());i++,j++){
	bundlex[j]=bundlex[i];
      }
      bundlex.resize(j);
    }
  }
  bundlevecs.concat_right(cand_subg);
  if (cand_x.size()>0){
    unsigned int sz=(unsigned int)(bundlex.size());
    bundlex.resize(sz+cand_x.size(),0);
    for (unsigned int i=sz;i<sz+cand_x.size();i++){
      bundlex[i]=cand_x[i-sz]->clone_primal_data();
    }
  }
  tmpvec=cand_subg_valvec;
  if (cand_y.dim()>0){
    genmult(cand_subg,cand_y,tmpvec,-1.,1.,1);
  }
  bundlecosts.concat_below(tmpvec);
  bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
  model_changed=1;
  
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));

  return 0;
}
  
int FunctionProblem::update_model2(bool)
{
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));
  if ((!cand_available)||(cand_subg_valvec.dim()==0)) return 1;
  if (bundlevecs.coldim()==0){
    if (maxkeepvecs>=cand_subg_valvec.dim()){
      bundlecosts=cand_subg_valvec;
      if (cand_y.dim()>0){
	genmult(cand_subg,cand_y,bundlecosts,-1.,1.,1);
      }
      bundlevecs=cand_subg;
      hashsum=sumrows(bundlevecs);
      bundlex.resize(Integer(cand_x.size()),0);
      for (unsigned int i=0;i<cand_x.size();i++){
	bundlex[i]=cand_x[i]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      Integer maxind;
      max(cand_subg_valvec,&maxind);
      bundlecoeff(maxind)=1.;
      model_changed=1;
    }
    else {
      sortindex(cand_subg_valvec,tmpind);
      tmpind=tmpind(Range(tmpind.dim()-maxkeepvecs,tmpind.dim()-1));
      bundlecosts=cand_subg_valvec(tmpind);
      bundlevecs=cand_subg.cols(tmpind);
      if (cand_y.dim()>0){
	genmult(bundlevecs,cand_y,bundlecosts,-1.,1.,1);
      }
      hashsum=sumrows(bundlevecs);
      bundlex.resize(min(Integer(cand_x.size()),maxkeepvecs),0);
      for (unsigned int i=0;i<bundlex.size();i++){
	bundlex[i]=cand_x[tmpind(i)]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(bundlecoeff.dim()-1)=1.;
      model_changed=1;
    }
    return 0;
  }

  if (bundlevecs.coldim()>=maxkeepvecs){
    tmpind=Range(1,bundlevecs.coldim()+1-maxkeepvecs);
    Integer aggregind=0;
    Integer ind=min(tmpind.dim(),bundlevecs.coldim()-1);
    for(Integer i=0;i<ind;i++){
      Real a=bundlecoeff(tmpind(i));
      Real b=bundlecoeff(aggregind);
      Real s=a+b;
      mat_xbpeya(bundlevecs.rowdim(),
		 bundlevecs.get_store()+aggregind*bundlevecs.rowdim(),
		 bundlevecs.get_store()+tmpind(i)*bundlevecs.rowdim(),
		 a/s,b/s);
      bundlecosts(aggregind)=(bundlecosts(tmpind(i))*a+bundlecosts(aggregind)*b)/s;
      bundlecoeff(aggregind)=s;
      if (bundlex.size()>0){
	bundlex[aggregind]->aggregate_primal_data(b/s,a/s,*bundlex[tmpind(i)]);
      }
    }
    bundlevecs.delete_cols(tmpind);
    bundlecosts.delete_rows(tmpind);
    bundlecoeff.delete_rows(tmpind);
    if ((bundlex.size()>0)&&(tmpind.dim()>0)){
      Indexmatrix sind;
      sortindex(tmpind,sind);
      tmpind=tmpind(sind);
      Integer j=tmpind(0);
      delete bundlex[j];
      Integer ti=1;
      Integer i=j+1;
      for(;(i<Integer(bundlex.size()))&&(ti<tmpind.dim());i++){
        if (tmpind(ti)==i) {
	  delete bundlex[i];
	  ti++;
	  continue;
	}
	bundlex[j]=bundlex[i];
        j++;
      }
      for(;i<Integer(bundlex.size());i++,j++){
	bundlex[j]=bundlex[i];
      }
      bundlex.resize(j);
    }
  }
  bundlevecs.concat_right(cand_subg);
  if (cand_x.size()>0){
    unsigned int sz=(unsigned int)(bundlex.size());
    bundlex.resize(sz+cand_x.size(),0);
    for (unsigned int i=sz;i<sz+cand_x.size();i++){
      bundlex[i]=cand_x[i-sz]->clone_primal_data();
    }
  }
  tmpvec=cand_subg_valvec;
  if (cand_y.dim()>0){
    genmult(cand_subg,cand_y,tmpvec,-1.,1.,1);
  }
  bundlecosts.concat_below(tmpvec);
  bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
  model_changed=1;
  
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));

  return 0;
}

int FunctionProblem::update_model3(bool)
{
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));
  if ((!cand_available)||(cand_subg_valvec.dim()==0)) return 1;
  if (bundlevecs.coldim()==0){
    if (maxkeepvecs>=cand_subg_valvec.dim()){
      bundlecosts=cand_subg_valvec;
      if (cand_y.dim()>0){
	genmult(cand_subg,cand_y,bundlecosts,-1.,1.,1);
      }
      bundlevecs=cand_subg;
      hashsum=sumrows(bundlevecs);
      bundlex.resize(Integer(cand_x.size()),0);
      for (unsigned int i=0;i<cand_x.size();i++){
	bundlex[i]=cand_x[i]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      Integer maxind;
      max(cand_subg_valvec,&maxind);
      bundlecoeff(maxind)=1.;
      model_changed=1;
    }
    else {
      sortindex(cand_subg_valvec,tmpind);
      tmpind=tmpind(Range(tmpind.dim()-maxkeepvecs,tmpind.dim()-1));
      bundlecosts=cand_subg_valvec(tmpind);
      bundlevecs=cand_subg.cols(tmpind);
      if (cand_y.dim()>0){
	genmult(bundlevecs,cand_y,bundlecosts,-1.,1.,1);
      }
      hashsum=sumrows(bundlevecs);
      bundlex.resize(min(Integer(cand_x.size()),maxkeepvecs),0);
      for (unsigned int i=0;i<bundlex.size();i++){
	bundlex[i]=cand_x[tmpind(i)]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(bundlecoeff.dim()-1)=1.;
      model_changed=1;
    }
    return 0;
  }

  //determine number of new cutting planes to be added
  //(those, that exceed the bundles cut value, but at most max_new_subg)
  tmpvec=bundlecosts;
  if (cand_y.dim()>0){
    genmult(bundlevecs,cand_y,tmpvec,1.,1.,1);
  }
  Real bundle_cutval=max(tmpvec);
  Indexmatrix acceptind;
  sortindex(-cand_subg_valvec,acceptind);
  Integer nnew=1;
  while((nnew<acceptind.dim())&&
	(nnew<max_new_subg)&&
	(cand_subg_valvec(acceptind(nnew))>bundle_cutval))
    nnew++;
  acceptind.reduce_length(nnew);

  //aggregate bundle vectors, if adding nnew exeeds maxkeepvecs
  if ((out)&&(print_level>=2)){
    out->precision(8);
    (*out)<<" bundle_update3: cval="<<bundle_cutval<<" nnew="<<nnew<<std::flush;
  }

  if (bundlecoeff.dim()+nnew>maxkeepvecs){

    //determine the range of bundlevectors, where aggregation will be done
    //(we will remove nnew from 1/3+nnew)
    Integer remove_n=bundlecoeff.dim()+nnew-maxkeepvecs;
    Indexmatrix delindvec;
    assert(bundlecoeff.dim()-remove_n>=1);
    Integer nconsidered=max(2,min(bundlecoeff.dim(),maxkeepvecs/3+remove_n));

    //compute the norm of the nconsidered first vectors in the bundle
    Matrix normvec(nconsidered,1,0.);
    for (Integer i=0;i<nconsidered;i++){
      normvec(i)=norm2(bundlevecs.col(i));
    }


    //first aggregate columns whose coefficient is negligible in comparison to the others
    tmpind=Indexmatrix(Range(0,nconsidered-1));
    Real aggregval=aggregtol*max(bundlecoeff(tmpind));
    if ((out)&&(print_level>=2)){
      out->precision(8);
      (*out)<<" ncons="<<nconsidered<<" norm["<<min(normvec)<<","<<max(normvec)<<"]";
      (*out)<<" coeff["<<min(bundlecoeff(tmpind))<<","<<max(bundlecoeff(tmpind))<<"]"<<std::flush;
    }
    Integer cntvaldel=0;
    Integer minind;
    while ((remove_n>0)&&(min(bundlecoeff(tmpind),&minind)<aggregval)){
      Integer delind=tmpind(minind);
      delindvec.concat_below(delind);
      tmpind.delete_rows(Indexmatrix(1,1,minind));
      remove_n--;
      cntvaldel++;
      //aggregate the column to the bundlevector pointing at most into the same direction
      Integer bestind=-1;
      Real bestipval=-2.;
      Matrix delcol=bundlevecs.col(delind);
      Real delnorm=normvec(delind);
      for (Integer i=0;i<tmpind.dim();i++){
	Integer ind=tmpind(i);
	Real indnorm=normvec(ind);
	Real ipval;
	if ((delnorm<1e-10)||(indnorm<1e-10)){
	  ipval=0.;
	}
	else {
	  ipval=ip(delcol,bundlevecs.col(ind))/indnorm/delnorm;
	}
	if (ipval>bestipval){
	  bestind=ind;
	  bestipval=ipval;
	}
      }
      assert(bestind>=0);
      Real a=bundlecoeff(delind);
      Real b=bundlecoeff(bestind);
      Real s=a+b;
      mat_xbpeya(bundlevecs.rowdim(),
		 bundlevecs.get_store()+bestind*bundlevecs.rowdim(),
		 bundlevecs.get_store()+delind*bundlevecs.rowdim(),
		 a/s,b/s);
      normvec(bestind)=norm2(bundlevecs.col(bestind));
      bundlecosts(bestind)=(bundlecosts(delind)*a+bundlecosts(bestind)*b)/s;
      bundlecoeff(bestind)=s;
      if (bundlex.size()>0){
	bundlex[bestind]->aggregate_primal_data(b/s,a/s,*bundlex[delind]);
      }
      aggregval=max(aggregval,aggregtol*s);
    }

    //next keep aggregating the vectors with largest inner product, till remove_n==0
    int cntmergedel=0;
    if (remove_n>0){
      //first compute all inner products
      tmpvec.init(tmpind.dim(),tmpind.dim(),0.);
      Integer besti=-1;
      Integer bestj=-1;
      Real bestval=-2.;
      for(Integer i=0;i<tmpind.dim();i++){
	Integer indi=tmpind(i);
	Real normi=normvec(indi);
	for(Integer j=i+1;j<tmpind.dim();j++){
	  Integer indj=tmpind(j);
	  Real normj=normvec(indj);
	  Real val;
	  if ((normi<1e-10)||(normj<1e-10)){
	    val=0.;
	  }
	  else {
	    val=mat_ip(bundlevecs.rowdim(),
		       bundlevecs.get_store()+indi*bundlevecs.rowdim(),
		       bundlevecs.get_store()+indj*bundlevecs.rowdim())
	      /normi/normj;
	  }
	  tmpvec(i,j)=val;
	  if (val>bestval){
	    besti=i;
	    bestj=j;
	    bestval=val;
	  }
	}
      }
      while(remove_n>0){
	Integer delind=tmpind(bestj);
	delindvec.concat_below(delind);
	tmpind.delete_rows(Indexmatrix(1,1,delind));
	remove_n--;
	cntmergedel++;
	Integer bestind=tmpind(besti);
	Real a=bundlecoeff(delind);
	Real b=bundlecoeff(bestind);
	Real s=a+b;
	mat_xbpeya(bundlevecs.rowdim(),
		   bundlevecs.get_store()+bestind*bundlevecs.rowdim(),
		   bundlevecs.get_store()+delind*bundlevecs.rowdim(),
		   a/s,b/s);
	Real oldnorm=normvec(bestind);
        Real newnorm=norm2(bundlevecs.col(bestind));
	normvec(bestind)=newnorm;
	bundlecosts(bestind)=(bundlecosts(delind)*a+bundlecosts(bestind)*b)/s;
	bundlecoeff(bestind)=s;
	if (bundlex.size()>0){
	  bundlex[bestind]->aggregate_primal_data(b/s,a/s,*bundlex[delind]);
	}
	if (remove_n>0){
	  //find next best pair
	  //first update tmpvec
	  Real delnorm=normvec(delind);
	  if (newnorm<1e-10){
	    {for(Integer i=0;i<besti;i++){
	      tmpvec(i,besti)=0.;
	    }}
	    {for(Integer i=besti+1;i<tmpvec.coldim();i++){
	      tmpvec(besti,i)=0.;
	    }}
	  }
	  else {
	    {for(Integer i=0;i<besti;i++){
	      tmpvec(i,besti)=(b/s*tmpvec(i,besti)*oldnorm+a/s*tmpvec(i,bestj)*delnorm)/newnorm;
	    }}
	    {for(Integer i=besti+1;i<bestj;i++){
	      tmpvec(besti,i)=(b/s*tmpvec(besti,i)*oldnorm+a/s*tmpvec(i,bestj)*delnorm)/newnorm;
	    }}
	    {for(Integer i=bestj+1;i<tmpvec.coldim();i++){
	      tmpvec(besti,i)=(b/s*tmpvec(besti,i)*oldnorm+a/s*tmpvec(bestj,i)*delnorm)/newnorm;
	    }}	    
	  }
	  tmpvec.delete_rows(Indexmatrix(1,1,bestj));
	  tmpvec.delete_cols(Indexmatrix(1,1,bestj));
	  assert((tmpvec.rowdim()==tmpind.dim())&&(tmpvec.coldim()==tmpind.dim()));
	  //next find the next best pair
	  bestval=-2.;
	  for(Integer i=0;i<tmpind.dim();i++){
	    for(Integer j=i+1;j<tmpind.dim();j++){
	      if (tmpvec(i,j)>bestval){
		besti=i;
		bestj=j;
		bestval=tmpvec(i,j);
	      }
	    }
	  }
	}
      }
    }

    //remove the columns in delindvec
    tmpind=delindvec;
    bundlevecs.delete_cols(tmpind);
    bundlecosts.delete_rows(tmpind);
    bundlecoeff.delete_rows(tmpind);
    if ((bundlex.size()>0)&&(tmpind.dim()>0)){
      Indexmatrix sind;
      sortindex(tmpind,sind);
      tmpind=tmpind(sind);
      Integer j=tmpind(0);
      delete bundlex[j];
      Integer ti=1;
      Integer i=j+1;
      for(;(i<Integer(bundlex.size()))&&(ti<tmpind.dim());i++){
        if (tmpind(ti)==i) {
	  delete bundlex[i];
	  ti++;
	  continue;
	}
	bundlex[j]=bundlex[i];
        j++;
      }
      for(;i<Integer(bundlex.size());i++,j++){
	bundlex[j]=bundlex[i];
      }
      bundlex.resize(j);
    }
    if ((out)&&(print_level>=2)){
      out->precision(8);
      (*out)<<" valdel="<<cntvaldel<<" merge="<<cntmergedel<<std::flush;
    }

  } //done with aggregation

  //append the selected new vectors
  bundlevecs.concat_right(cand_subg.cols(acceptind));
  if (cand_x.size()>0){
    Integer sz=Integer(bundlex.size());
    bundlex.resize(sz+acceptind.dim(),0);
    for (Integer i=sz;i<sz+acceptind.dim();i++){
      bundlex[i]=cand_x[acceptind(i-sz)]->clone_primal_data();
    }
  }
  tmpvec=cand_subg_valvec(acceptind);
  if (cand_y.dim()>0){
    genmult(cand_subg.cols(acceptind),cand_y,tmpvec,-1.,1.,1);
  }
  bundlecosts.concat_below(tmpvec);
  bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
  model_changed=1;

  if ((out)&&(print_level>=2)){
    out->precision(8);
    (*out)<<std::endl;
  }
  
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));

  return 0;
}
  
int FunctionProblem::update_model4(bool)
{
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));
  if ((!cand_available)||(cand_subg_valvec.dim()==0)) return 1;
  if (bundlevecs.coldim()==0){
    if (maxkeepvecs>=cand_subg_valvec.dim()){
      bundlecosts=cand_subg_valvec;
      if (cand_y.dim()>0){
	genmult(cand_subg,cand_y,bundlecosts,-1.,1.,1);
      }
      bundlevecs=cand_subg;
      hashsum=sumrows(bundlevecs);
      bundlex.resize(Integer(cand_x.size()),0);
      for (unsigned int i=0;i<cand_x.size();i++){
	bundlex[i]=cand_x[i]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      Integer maxind;
      max(cand_subg_valvec,&maxind);
      bundlecoeff(maxind)=1.;
      model_changed=1;
    }
    else {
      sortindex(cand_subg_valvec,tmpind);
      tmpind=tmpind(Range(tmpind.dim()-maxkeepvecs,tmpind.dim()-1));
      bundlecosts=cand_subg_valvec(tmpind);
      bundlevecs=cand_subg.cols(tmpind);
      if (cand_y.dim()>0){
	genmult(bundlevecs,cand_y,bundlecosts,-1.,1.,1);
      }
      hashsum=sumrows(bundlevecs);
      bundlex.resize(min(Integer(cand_x.size()),maxkeepvecs),0);
      for (unsigned int i=0;i<bundlex.size();i++){
	bundlex[i]=cand_x[tmpind(i)]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(bundlecoeff.dim()-1)=1.;
      model_changed=1;
    }
    return 0;
  }

  //determine number of new cutting planes to be added
  //(those, that exceed the bundles cut value, but at most max_new_subg)
  tmpvec=bundlecosts;
  if (cand_y.dim()>0){
    genmult(bundlevecs,cand_y,tmpvec,1.,1.,1);
}
  Real bundle_cutval=max(tmpvec);
  Indexmatrix acceptind;
  sortindex(-cand_subg_valvec,acceptind);
  Integer nnew=1;
  while((nnew<acceptind.dim())&&
	(nnew<max_new_subg)&&
	(cand_subg_valvec(acceptind(nnew))>bundle_cutval))
    nnew++;
  acceptind.reduce_length(nnew);

  //aggregate bundle vectors, if adding nnew exeeds maxkeepvecs
  if ((out)&&(print_level>=2)){
    out->precision(8);
    (*out)<<" bundle_update4: cval="<<bundle_cutval<<" nnew="<<nnew<<std::flush;
  }

  if (bundlecoeff.dim()+nnew>maxkeepvecs){

    //determine the range of bundlevectors, where aggregation will be done
    //(we will remove nnew from 1/3+nnew)
    Integer remove_n=bundlecoeff.dim()+nnew-maxkeepvecs;
    Indexmatrix delindvec;
    assert(bundlecoeff.dim()-remove_n>=1);
    Integer nconsidered=max(2,min(bundlecoeff.dim(),maxkeepvecs/3+remove_n));

    //compute the norm of the nconsidered first vectors in the bundle
    Matrix normvec(nconsidered,1,0.);
    for (Integer i=0;i<nconsidered;i++){
      normvec(i)=norm2(bundlevecs.col(i));
    }


    //first aggregate columns whose coefficient is negligible in comparison to the others
    tmpind=Indexmatrix(Range(0,nconsidered-1));
    if ((out)&&(print_level>=2)){
      out->precision(8);
      Indexmatrix sind;
      sortindex(bundlecoeff(tmpind),sind);
      (*out)<<" ncons="<<nconsidered<<" norm["<<min(normvec)<<","<<max(normvec)<<"]";
      (*out)<<" coeff["<<bundlecoeff(sind(sind.dim()-1));
      if (sind.dim()>1)
	(*out)<<","<<bundlecoeff(sind(sind.dim()-2));
      (*out)<<"] del="<<remove_n<<std::flush;
    }
    while (remove_n>0){
      Integer minind;
      min(bundlecoeff(tmpind),&minind);
      Integer delind=tmpind(minind);
      delindvec.concat_below(delind);
      tmpind.delete_rows(Indexmatrix(1,1,minind));
      remove_n--;
      //aggregate the column to the bundlevector pointing at most into the same direction
      Integer bestind=-1;
      Real bestipval=-2.;
      Matrix delcol=bundlevecs.col(delind);
      Real delnorm=normvec(delind);
      for (Integer i=0;i<tmpind.dim();i++){
	Integer ind=tmpind(i);
	Real indnorm=normvec(ind);
	Real ipval;
	if ((delnorm<1e-10)||(indnorm<1e-10)){
	  ipval=0.;
	}
	else {
	  ipval=ip(delcol,bundlevecs.col(ind))/indnorm/delnorm;
	}
	if (ipval>bestipval){
	  bestind=ind;
	  bestipval=ipval;
	}
      }
      assert(bestind>=0);
      Real a=bundlecoeff(delind);
      Real b=bundlecoeff(bestind);
      Real s=a+b;
      mat_xbpeya(bundlevecs.rowdim(),
		 bundlevecs.get_store()+bestind*bundlevecs.rowdim(),
		 bundlevecs.get_store()+delind*bundlevecs.rowdim(),
		 a/s,b/s);
      normvec(bestind)=norm2(bundlevecs.col(bestind));
      bundlecosts(bestind)=(bundlecosts(delind)*a+bundlecosts(bestind)*b)/s;
      bundlecoeff(bestind)=s;
      if (bundlex.size()>0){
	bundlex[bestind]->aggregate_primal_data(b/s,a/s,*bundlex[delind]);
      }
    }

    //remove the columns in delindvec
    tmpind=delindvec;
    bundlevecs.delete_cols(tmpind);
    bundlecosts.delete_rows(tmpind);
    bundlecoeff.delete_rows(tmpind);
    if ((bundlex.size()>0)&&(tmpind.dim()>0)){
      Indexmatrix sind;
      sortindex(tmpind,sind);
      tmpind=tmpind(sind);
      Integer j=tmpind(0);
      delete bundlex[j];
      Integer ti=1;
      Integer i=j+1;
      for(;(i<Integer(bundlex.size()))&&(ti<tmpind.dim());i++){
        if (tmpind(ti)==i) {
	  delete bundlex[i];
	  ti++;
	  continue;
	}
	bundlex[j]=bundlex[i];
        j++;
      }
      for(;i<Integer(bundlex.size());i++,j++){
	bundlex[j]=bundlex[i];
      }
      bundlex.resize(j);
    }

  } //done with aggregation

  //append the selected new vectors
  bundlevecs.concat_right(cand_subg.cols(acceptind));
  if (cand_x.size()>0){
    Integer sz=Integer(bundlex.size());
    bundlex.resize(sz+acceptind.dim(),0);
    for (Integer i=sz;i<sz+acceptind.dim();i++){
      bundlex[i]=cand_x[acceptind(i-sz)]->clone_primal_data();
    }
  }
  tmpvec=cand_subg_valvec(acceptind);
  if (cand_y.dim()>0){
    genmult(cand_subg.cols(acceptind),cand_y,tmpvec,-1.,1.,1);
  }
  bundlecosts.concat_below(tmpvec);
  bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
  model_changed=1;

  if ((out)&&(print_level>=2)){
    out->precision(8);
    (*out)<<std::endl;
  }
  
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));

  return 0;
}
  
int FunctionProblem::update_model5(bool)
{
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));
  if ((!cand_available)||(cand_subg_valvec.dim()==0)) return 1;
  if (bundlevecs.coldim()==0){
    if (maxkeepvecs>=cand_subg_valvec.dim()){
      bundlecosts=cand_subg_valvec;
      if (cand_y.dim()>0){
	genmult(cand_subg,cand_y,bundlecosts,-1.,1.,1);
      }
      bundlevecs=cand_subg;
      hashsum=sumrows(bundlevecs);
      bundlex.resize(Integer(cand_x.size()),0);
      for (unsigned int i=0;i<cand_x.size();i++){
	bundlex[i]=cand_x[i]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      Integer maxind;
      max(cand_subg_valvec,&maxind);
      bundlecoeff(maxind)=1.;
      model_changed=1;
    }
    else {
      sortindex(cand_subg_valvec,tmpind);
      tmpind=tmpind(Range(tmpind.dim()-maxkeepvecs,tmpind.dim()-1));
      bundlecosts=cand_subg_valvec(tmpind);
      bundlevecs=cand_subg.cols(tmpind);
      if (cand_y.dim()>0){
	genmult(bundlevecs,cand_y,bundlecosts,-1.,1.,1);
      }
      hashsum=sumrows(bundlevecs);
      bundlex.resize(min(Integer(cand_x.size()),maxkeepvecs),0);
      for (unsigned int i=0;i<bundlex.size();i++){
	bundlex[i]=cand_x[tmpind(i)]->clone_primal_data();
      }
      bundlecoeff.init(bundlecosts.dim(),1,0.);
      bundlecoeff(bundlecoeff.dim()-1)=1.;
      model_changed=1;
    }
    return 0;
  }

  //determine number of new cutting planes to be added
  //(those, that exceed the bundles cut value, but at most max_new_subg)
  tmpvec=bundlecosts;
  if (cand_y.dim()>0){
    genmult(bundlevecs,cand_y,tmpvec,1.,1.,1);
  }
  Real bundle_cutval=max(tmpvec);
  Indexmatrix acceptind;
  sortindex(-cand_subg_valvec,acceptind);
  Integer nnew=1;
  while((nnew<acceptind.dim())&&
	(nnew<max_new_subg)&&
	(cand_subg_valvec(acceptind(nnew))>bundle_cutval))
    nnew++;
  acceptind.reduce_length(nnew);

  //aggregate bundle vectors, if adding nnew exeeds maxkeepvecs
  if ((out)&&(print_level>=2)){
    out->precision(8);
    (*out)<<" bundle_update5: cval="<<bundle_cutval<<" nnew="<<nnew<<std::flush;
  }

  if (bundlecoeff.dim()+nnew>maxkeepvecs){

    //determine the range of bundlevectors, where aggregation will be done
    //(we will remove nnew from 1/3+nnew)
    Integer remove_n=bundlecoeff.dim()+nnew-maxkeepvecs;
    Indexmatrix delindvec;
    assert(bundlecoeff.dim()-remove_n>=1);
    Integer nconsidered=max(2,min(bundlecoeff.dim(),maxkeepvecs/3+remove_n));

    //compute the norm of the nconsidered first vectors in the bundle
    Matrix normvec(nconsidered,1,0.);
    for (Integer i=0;i<nconsidered;i++){
      normvec(i)=norm2(bundlevecs.col(i));
    }
    tmpind=Indexmatrix(Range(0,nconsidered-1));

    if ((out)&&(print_level>=2)){
      out->precision(8);
      Indexmatrix sind;
      sortindex(bundlecoeff(tmpind),sind);
      (*out)<<" ncons="<<nconsidered<<" norm["<<min(normvec)<<","<<max(normvec)<<"]";
      (*out)<<" coeff["<<bundlecoeff(sind(sind.dim()-1));
      if (sind.dim()>1)
	(*out)<<","<<bundlecoeff(sind(sind.dim()-2));
      (*out)<<"]"<<std::flush;
    }

    //first aggregate the vectors with largest inner product, if almost identical

    int cntmergedel=0;
    if (remove_n>0){
      //first compute all inner products
      tmpvec.init(tmpind.dim(),tmpind.dim(),0.);
      Integer besti=-1;
      Integer bestj=-1;
      Real bestval=-2.;
      for(Integer i=0;i<tmpind.dim();i++){
	Integer indi=tmpind(i);
	Real normi=normvec(indi);
	for(Integer j=i+1;j<tmpind.dim();j++){
	  Integer indj=tmpind(j);
	  Real normj=normvec(indj);
	  Real val;
	  if ((normi<1e-10)||(normj<1e-10)){
	    val=0.;
	  }
	  else {
	    val=mat_ip(bundlevecs.rowdim(),
		       bundlevecs.get_store()+indi*bundlevecs.rowdim(),
		       bundlevecs.get_store()+indj*bundlevecs.rowdim())
	      /normi/normj;
	  }
	  tmpvec(i,j)=val;
	  if (val>bestval){
	    besti=i;
	    bestj=j;
	    bestval=val;
	  }
	}
      }
      if ((out)&&(print_level>=2)){
	out->precision(8);
	(*out)<<" best="<<bestval<<std::flush;
      }
      while((remove_n>0)&&(bestval>1.-1e-4)){
	Integer delind=tmpind(bestj);
	delindvec.concat_below(delind);
	tmpind.delete_rows(Indexmatrix(1,1,delind));
	remove_n--;
	cntmergedel++;
	Integer bestind=tmpind(besti);
	Real a=bundlecoeff(delind);
	Real b=bundlecoeff(bestind);
	Real s=a+b;
	mat_xbpeya(bundlevecs.rowdim(),
		   bundlevecs.get_store()+bestind*bundlevecs.rowdim(),
		   bundlevecs.get_store()+delind*bundlevecs.rowdim(),
		   a/s,b/s);
	Real oldnorm=normvec(bestind);
        Real newnorm=norm2(bundlevecs.col(bestind));
	normvec(bestind)=newnorm;
	bundlecosts(bestind)=(bundlecosts(delind)*a+bundlecosts(bestind)*b)/s;
	bundlecoeff(bestind)=s;
	if (bundlex.size()>0){
	  bundlex[bestind]->aggregate_primal_data(b/s,a/s,*bundlex[delind]);
	}
	if (remove_n>0){
	  //find next best pair
	  //first update tmpvec
	  Real delnorm=normvec(delind);
	  if (newnorm<1e-10){
	    {for(Integer i=0;i<besti;i++){
	      tmpvec(i,besti)=0.;
	    }}
	    {for(Integer i=besti+1;i<tmpvec.coldim();i++){
	      tmpvec(besti,i)=0.;
	    }}
	  }
	  else {
	    {for(Integer i=0;i<besti;i++){
	      tmpvec(i,besti)=(b/s*tmpvec(i,besti)*oldnorm+a/s*tmpvec(i,bestj)*delnorm)/newnorm;
	    }}
	    {for(Integer i=besti+1;i<bestj;i++){
	      tmpvec(besti,i)=(b/s*tmpvec(besti,i)*oldnorm+a/s*tmpvec(i,bestj)*delnorm)/newnorm;
	    }}
	    {for(Integer i=bestj+1;i<tmpvec.coldim();i++){
	      tmpvec(besti,i)=(b/s*tmpvec(besti,i)*oldnorm+a/s*tmpvec(bestj,i)*delnorm)/newnorm;
	    }}	    
	  }
	  tmpvec.delete_rows(Indexmatrix(1,1,bestj));
	  tmpvec.delete_cols(Indexmatrix(1,1,bestj));
	  assert((tmpvec.rowdim()==tmpind.dim())&&(tmpvec.coldim()==tmpind.dim()));
	  //next find the next best pair
	  bestval=-2.;
	  for(Integer i=0;i<tmpind.dim();i++){
	    for(Integer j=i+1;j<tmpind.dim();j++){
	      if (tmpvec(i,j)>bestval){
		besti=i;
		bestj=j;
		bestval=tmpvec(i,j);
	      }
	    }
	  }
	}
      }
    }


    //next aggregate columns whose coefficient is negligible in comparison to the others
    Integer cntvaldel=0;
    while (remove_n>0){
      Integer minind;
      min(bundlecoeff(tmpind),&minind);
      Integer delind=tmpind(minind);
      delindvec.concat_below(delind);
      tmpind.delete_rows(Indexmatrix(1,1,minind));
      remove_n--;
      cntvaldel++;
      //aggregate the column to the bundlevector pointing at most into the same direction
      Integer bestind=-1;
      Real bestipval=-2.;
      Matrix delcol=bundlevecs.col(delind);
      Real delnorm=normvec(delind);
      for (Integer i=0;i<tmpind.dim();i++){
	Integer ind=tmpind(i);
	Real indnorm=normvec(ind);
	Real ipval;
	if ((delnorm<1e-10)||(indnorm<1e-10)){
	  ipval=0.;
	}
	else {
	  ipval=ip(delcol,bundlevecs.col(ind))/indnorm/delnorm;
	}
	if (ipval>bestipval){
	  bestind=ind;
	  bestipval=ipval;
	}
      }
      assert(bestind>=0);
      Real a=bundlecoeff(delind);
      Real b=bundlecoeff(bestind);
      Real s=a+b;
      mat_xbpeya(bundlevecs.rowdim(),
		 bundlevecs.get_store()+bestind*bundlevecs.rowdim(),
		 bundlevecs.get_store()+delind*bundlevecs.rowdim(),
		 a/s,b/s);
      normvec(bestind)=norm2(bundlevecs.col(bestind));
      bundlecosts(bestind)=(bundlecosts(delind)*a+bundlecosts(bestind)*b)/s;
      bundlecoeff(bestind)=s;
      if (bundlex.size()>0){
	bundlex[bestind]->aggregate_primal_data(b/s,a/s,*bundlex[delind]);
      }
    }

    //remove the columns in delindvec
    tmpind=delindvec;
    bundlevecs.delete_cols(tmpind);
    bundlecosts.delete_rows(tmpind);
    bundlecoeff.delete_rows(tmpind);
    if ((bundlex.size()>0)&&(tmpind.dim()>0)){
      Indexmatrix sind;
      sortindex(tmpind,sind);
      tmpind=tmpind(sind);
      Integer j=tmpind(0);
      delete bundlex[j];
      Integer ti=1;
      Integer i=j+1;
      for(;(i<Integer(bundlex.size()))&&(ti<tmpind.dim());i++){
        if (tmpind(ti)==i) {
	  delete bundlex[i];
	  ti++;
	  continue;
	}
	bundlex[j]=bundlex[i];
        j++;
      }
      for(;i<Integer(bundlex.size());i++,j++){
	bundlex[j]=bundlex[i];
      }
      bundlex.resize(j);
    }
    if ((out)&&(print_level>=2)){
      out->precision(8);
      (*out)<<" valdel="<<cntvaldel<<" merge="<<cntmergedel<<std::flush;
    }

  } //done with aggregation

  //append the selected new vectors
  bundlevecs.concat_right(cand_subg.cols(acceptind));
  if (cand_x.size()>0){
    Integer sz=Integer(bundlex.size());
    bundlex.resize(sz+acceptind.dim(),0);
    for (Integer i=sz;i<sz+acceptind.dim();i++){
      bundlex[i]=cand_x[acceptind(i-sz)]->clone_primal_data();
    }
  }
  tmpvec=cand_subg_valvec(acceptind);
  if (cand_y.dim()>0){
    genmult(cand_subg.cols(acceptind),cand_y,tmpvec,-1.,1.,1);
  }
  bundlecosts.concat_below(tmpvec);
  bundlecoeff.concat_below(Matrix(tmpvec.dim(),1,0.));
  model_changed=1;

  if ((out)&&(print_level>=2)){
    out->precision(8);
    (*out)<<std::endl;
  }
  
  assert((bundlex.size()==0)||(Integer(bundlex.size())==bundlevecs.coldim()));
  assert((bundlecoeff.dim()==bundlevecs.coldim()));

  return 0;
}
  
  
// *****************************************************************************
//                              intersect_box
// *****************************************************************************

int FunctionProblem::intersect_box(Indexmatrix& inbindex,Matrix& inlb,Matrix& inub)
{
  chk_init(inbindex);
  chk_init(inlb);
  chk_init(inub);
  if ((inbindex.coldim()==0)&&(inlb.coldim()==0)&&(inub.coldim()==0)){
    inbindex=bounds_index;
    if (lby.coldim()==0) 
      inlb.init(dim,1,CB_minus_infinity);
    else inlb=lby;
    if (uby.coldim()==0) 
      inub.init(dim,1,CB_plus_infinity);
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

Real FunctionProblem::lb_function(const Matrix& iny)
{
  if (aug_available){
    return aug_linconst+ip(aug_subg,iny);
  }
  if (cand_available){
    tmpvec=cand_subg_valvec;
    if (cand_y.dim()>0){
      return max(genmult(cand_subg,iny-cand_y,tmpvec,-1.,1.,1));
    }
    return max(cand_subg_valvec);
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

Real FunctionProblem::lb_model(const Matrix& iny)
{
  if (aug_available){
    return aug_linconst+ip(aug_subg,iny);
  }
  if (cand_available){
    tmpvec=cand_subg_valvec;
    if (cand_y.dim()>0){
      return max(genmult(cand_subg,iny-cand_y,tmpvec,-1.,1.,1));
    }
    return max(cand_subg_valvec);
  }
  return CB_minus_infinity;
}


// *****************************************************************************
//                                start_augmodel
// *****************************************************************************

//return a pointer to the variables/constraints generating the cutting model
//returns 0 on success, 1 on failure

int FunctionProblem::start_augmodel(QP_Block*& blockp)
{
  if (bundlevecs.coldim()==0) {
    if ((cand_available)&&(cand_subg_valvec.dim()>0)){
      if (update_model(false)) return 1;
    }
    else if (center_available){
      bundlecosts=Matrix(1,1,subg_val-ip(subg,y));
      bundlevecs=subg;
      hashsum=sumrows(bundlevecs);
      if (x) 
	bundlex.push_back(x->clone_primal_data());
      bundlecoeff.init(1,1,1.);
    }
    else return 1;
  }
  //--- set blocks
  //initialize blocks and variable dimensions
  aug_available=0;
  tmpind.init(0,0,Integer(0));
  block.init_block(bundlevecs.coldim(),tmpind,tmpind,1.,0);

 
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

int FunctionProblem::get_row(Integer index_y,Matrix& row,Real& /* bi */,Integer startindex) const
{
  if (index_y==-1){
    mat_xey(bundlecosts.dim(),row.get_store()+startindex,bundlecosts.get_store());
  }
  else {
    mat_xey(bundlevecs.coldim(),
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

int FunctionProblem::make_aug_linmodel(Real* in_aug_linconst,Matrix* in_aug_subg,bool* /* conemult_increased */,Real* /* function_value */)  
{

  block.get_linx(bundlecoeff);
  aug_linconst=ip(bundlecosts,bundlecoeff);
  genmult(bundlevecs,bundlecoeff,aug_subg);
  
  if (in_aug_linconst) *in_aug_linconst+=aug_linconst;
  if (in_aug_subg) *in_aug_subg+=aug_subg;

  aug_available=1;
 
  return 0;
} 

// *****************************************************************************
//                               change_variables
// *****************************************************************************

//change variables as described in ChangeVariableInfo
//returns 0 on success, 1 on failure
//(this is not needed in the bundle framework, routine may return 1!)

int FunctionProblem::change_variables(ChangeVarInfo* cvp)
{
  if (typeid(*cvp)==typeid(AppendVars)){
    AppendVars* cp=dynamic_cast<AppendVars*>(cvp);
    Integer n_append=cp->n_append;
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
	    delete x; x=0;
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
      if (oracle.subgradient_extension(pd,ivec,dvec)){
	if (out) (*out)<<"**** WARNING: FunctionProblem::add_variable(): subgradient extension failed"<<std::endl;
	no_subg_extension=1;
	break;
      }
      if (dvec.dim()!=n_append){
	if (out) (*out)<<"**** WARNING: FunctionProblem::add_variable(): oracle.subgradient_extension(...) returns a vector of wrong dimension, subgradient extension failed"<<std::endl;
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
	  genmult(bundlevecs,bundlecoeff,subg);
	  subg_val=ip(bundlecosts,bundlecoeff)+ip(subg,y);
	  delete x;x=0;
	  if (bundlex.size()>0) {
	    x=bundlex[0]->clone_primal_data();
	    if (bundlecoeff.dim()>1){
	      x->aggregate_primal_data(bundlecoeff[0],bundlecoeff[1],*bundlex[1]);
	    }
	    for(Integer i=2;i<bundlecoeff.dim();i++){
	      x->aggregate_primal_data(1.,bundlecoeff[i],*bundlex[i]);
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
      tmpind.init(dim,1,Integer(0));
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
    
    if (lby.coldim()!=0) 
      lby=lby(cp->assign_ind);
    if (uby.coldim()!=0) 
      uby=uby(cp->assign_ind);
    if (bounds_index.coldim()!=0){
      bounds_index.init(0,0,Integer(0));
      for(Integer i=0;i<dim;i++){
	if (((lby.coldim()!=0)&&(lby(i)>CB_minus_infinity))||
	    ((uby.coldim()!=0)&&(uby(i)<CB_plus_infinity))){
	  bounds_index.concat_below(i);
	}
      }
    }
    y=y(cp->assign_ind);
    if (bundlevecs.dim()!=0){
      bundlevecs=bundlevecs.rows(cp->assign_ind);
      hashsum=sumrows(bundlevecs);
    }
    
    if (center_available) {
      //make sure that the center subgradient is contained in the model
      if (bundlecoeff.dim()>0){
	genmult(bundlevecs,bundlecoeff,subg);
	subg_val=ip(bundlecosts,bundlecoeff)+ip(subg,y);
	delete x;x=0;
	if (bundlex.size()>0) {
	  x=bundlex[0]->clone_primal_data();
	  if (bundlecoeff.dim()>1){
	    x->aggregate_primal_data(bundlecoeff[0],bundlecoeff[1],*bundlex[1]);
	  }
	  for(Integer i=2;i<bundlecoeff.dim();i++){
	    x->aggregate_primal_data(1.,bundlecoeff[i],*bundlex[i]);
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
    if (lby.coldim()!=0) lby.delete_rows(cp->del_index);
    if (uby.coldim()!=0) uby.delete_rows(cp->del_index);
    y.delete_rows(cp->del_index);
    if (bundlevecs.dim()!=0){
      bundlevecs.delete_rows(cp->del_index);
      hashsum=sumrows(bundlevecs);
    }
    
    if (center_available) {
      //make sure that the center subgradient is contained in the model
      if (bundlecoeff.dim()>0){
	genmult(bundlevecs,bundlecoeff,subg);
	subg_val=ip(bundlecosts,bundlecoeff)+ip(subg,y);
	delete x;x=0;
	if (bundlex.size()>0) {
	  x=bundlex[0]->clone_primal_data();
	  if (bundlecoeff.dim()>1){
	    x->aggregate_primal_data(bundlecoeff[0],bundlecoeff[1],*bundlex[1]);
	  }
	  for(Integer i=2;i<bundlecoeff.dim();i++){
	    x->aggregate_primal_data(1.,bundlecoeff[i],*bundlex[i]);
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

int FunctionProblem::recompute_center()
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

int FunctionProblem::get_approximate_primal(PrimalData& primal) const 
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

void FunctionProblem::clear_model()  
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
  delete x;
  x=0;
  {for(unsigned int i=0;i<bundlex.size();i++){
    delete bundlex[i];
  }}
  bundlex.clear();
}

// *****************************************************************************
//                                clear_model
// *****************************************************************************

  //if the function is the Lagrangian dual and primal_data of previous calls 
  //has now to be updated due to changes in the primal problem -- e.g., this 
  //may happen in column generation -- the problem updates all its internally 
  //stored primal_data objects by calling PrimalExtender::extend on
  //each of these.
  //returns 0 on success, 1 if not applicable to this function, 2 if it
  //would be applicable but there is no primal data.

int FunctionProblem::call_primal_extender(PrimalExtender& prex)
{
  bool used=false;
  if (x){
    used=true;
    if (prex.extend(*x)){
      if (out) (*out)<<"**** WARNING: FunctionProblem::call_primal_extender(): PrimalExtender::extend failed"<<std::endl;
    }
  }
  for(unsigned int i=0;i<cand_x.size();i++){
    used=true;
    if (prex.extend(*(cand_x[i]))){
      if (out) (*out)<<"**** WARNING: FunctionProblem::call_primal_extender(): PrimalExtender::extend failed"<<std::endl;
    }
  }
 
  for(unsigned int i=0;i<bundlex.size();i++){
    used=true;
    if (prex.extend(*(bundlex[i]))){
      if (out) (*out)<<"**** WARNING: FunctionProblem::call_primal_extender(): PrimalExtender::extend failed"<<std::endl;
    }
  }
  if (!used) return 2;
  return 0; 
}
   

} 
