/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/Nbundle.cxx

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



#include <math.h>
#include <stdlib.h>
#include "Nbundle.hxx"
#include "hkweight.hxx"
#include "mymath.hxx"
#include <algorithm>

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {


// *****************************************************************************
//                                 NBundleSolver
// *****************************************************************************

NBundleSolver::NBundleSolver()
{
  terminator=new BundleTerminator;
  bundleweight=new BundleHKweight(.5);
  problem=0; 
  out=0; 
  print_level=1;
  clockp=0; 
 
  set_defaults();
  clear();
}

// *****************************************************************************
//                                 ~NBundleSolver
// *****************************************************************************

NBundleSolver::~NBundleSolver()
{
  delete terminator; 
  delete bundleweight;
}

// *****************************************************************************
//                                 set_defaults
// *****************************************************************************

// usually called on construction

void NBundleSolver::set_defaults()
{
  modeleps=0.6; 
  mL=0.1;   
  mN=0.1;

  yfixing_allowed=false;
  yfixing_itbound=10;
  yfixing_factor=.5;

  do_scaling=0;
  do_scaling_heuristic=0;
  max_updates=10; //was 100;
  use_linval=1;
  if (terminator) terminator->set_defaults();
  if (bundleweight) bundleweight->set_defaults();
}

// *****************************************************************************
//                                   clear
// *****************************************************************************

// usually called on construction, resets all matrices and calls set_defaults

void NBundleSolver::clear()
{
  problem=0;       
  terminate=0;
  weightu=-1;      
  y.init(0,1,0.);
  lby.init(0,1,0.);
  uby.init(0,1,0.);
  b.init(0,1,0.);
  bound_index.init(0,1,Integer(0)); 
  newy.init(0,1,0.);
  center_subg.init(0,1,0.);
  aggr_subg.init(0,1,0.);
  new_subg.init(0,1,0.);
  eta.init(0,1,0.);
  normsubg2=-1;
  inv_scale.init(0,1,0.); 
  yfixed.init(0,1,Integer(0)); 
  update_index.init(0,1,Integer(0));
  update_value.init(0,1,0.);
  updatecnt=0;
  sumupdatecnt=0;
  retcode=0;
  initialize_center=0;
  recompute_bound_index=false;
  initialize_aggregate=true;
  recompute_normsubg2=true;
  initialize_model=true;
  oldval=0.;
  newval=0.;
  linval=0.;
  cutval=0.;
  modelval=0.;
  augval=0.;
  oldaugval=0.;
  lastaugval=0.;
  innerit=0;
  suminnerit=0; 
  cntobjeval=0; 
  recomp=0;     
  sumrecomp=0;  
  qpfails=0;    
  sumqpfails=0; 
  modelfails=0;   
  summodelfails=0;
  augvalfails=0;                    
  sumaugvalfails=0;
  oraclefails=0;
  sumoraclefails=0;
  cntsolves=0;
  cntresolves=0;
  cntitersolves=0;
  cntiterresolves=0;
  descent_steps=0;        
  bens=0,       
  sens=0,       
  dens=0,       
  shallowcut=0; 
  sensvkappasum=0;
  if (terminator) terminator->clear();
  if (bundleweight) bundleweight->clear();
}

// *****************************************************************************
//                                   clear_fails
// *****************************************************************************

// resets all fail values

void NBundleSolver::clear_fails()
{
  recomp=0;     
  sumrecomp=0;  
  qpfails=0;    
  sumqpfails=0; 
  modelfails=0;   
  summodelfails=0;
  augvalfails=0;                    
  sumaugvalfails=0;
  oraclefails=0;
  sumoraclefails=0;
}

// *****************************************************************************
//                               compute_newy
// *****************************************************************************

// computes also the lagrange multipliers corresponding to the subgradient

int NBundleSolver::compute_newy(const Matrix &subg)
{
  etaval=0.;

  if (do_scaling){

    //------ with scaling
    newy.init(subg,-1./weightu);
    newy%=inv_scale;
    normsubg2=ip(subg,newy)*(-weightu);
    newy+=y;
    if (yfixed.dim()>0){
      for(Integer i=0;i<yfixed.dim();i++){
	if (yfixed(i)) {
	  newy(i)=y(i);
	  eta(i)=0.;
	  normsubg2-=subg(i)*subg(i)*inv_scale(i);
	}
      }
    }

    if (bound_index.dim()<y.dim()/2){ 
      
      //bounds on "few" variables
      for(Integer i=0;i<bound_index.dim();i++){
	Integer ind=bound_index(i);
	Real d=newy(ind);
	Real b=lby(ind);
	if (b>d){
	  Real s=inv_scale(ind);
	  Real e=weightu*(b-d)/s;
	  eta(ind)=e;
	  etaval+=b*e;
	  newy(ind)=b;
	  normsubg2+=(e-2*subg(ind))*e*s;
	  continue;
	}
        b=uby(ind);
	if (b<d){
	  Real s=inv_scale(ind);
	  Real e=weightu*(b-d)/s;
	  eta(ind)=e;
	  etaval+=b*e;
	  newy(ind)=b;
	  normsubg2+=(e-2*subg(ind))*e*s;
	  continue;
	}
	eta(ind)=0.;
      }	  
    }
    
    else {
      
      //bounds on "all" variables
      //newy=y-subg%(weightu*inv_scale)^(-1);
      for(Integer i=0;i<y.dim();i++){
	Real d=newy(i);
	Real b=lby(i);
	if (b>d){
	  Real s=inv_scale(i);
	  Real e=weightu*(b-d)/s;
	  eta(i)=e;
	  etaval+=b*e;
	  newy(i)=b;
	  normsubg2+=(e-2*subg(i))*e*s;
	  continue;
	}
        b=uby(i);
	if (b<d){
	  Real s=inv_scale(i);
	  Real e=weightu*(b-d)/s;
	  eta(i)=e;
	  etaval+=b*e;
	  newy(i)=b;
	  normsubg2+=(e-2*subg(i))*e*s;
	  continue;
	}
	eta(i)=0.;
      }	  
    } //end bounds on "all"

  }
  else {
    //------ without scaling
    normsubg2=ip(subg,subg);
    xeyapzb(newy,y,subg,1.,-1/weightu);  //newy=y-subg/weightu;
    if (yfixed.dim()>0){
      for(Integer i=0;i<yfixed.dim();i++){
	if (yfixed(i)) {
	  newy(i)=y(i);
	  eta(i)=0.;
	  normsubg2-=subg(i)*subg(i);
	}
      }
    }
    
    if (bound_index.dim()<y.dim()/2){ 
      
      //bounds on "few" variables
      for(Integer i=0;i<bound_index.dim();i++){
	Integer ind=bound_index(i);
	Real d=newy(ind);
	Real b=lby(ind);
	if (b>d){
	  Real e=weightu*(b-d);
	  eta(ind)=e;
	  etaval+=b*e;
	  newy(ind)=b;
	  normsubg2+=(e-2.*subg(ind))*e;
	  continue;
	}
        b=uby(ind);
	if (b<d){
	  Real e=weightu*(b-d);
	  eta(ind)=e;
	  etaval+=b*e;
	  newy(ind)=b;
	  normsubg2+=(e-2.*subg(ind))*e;
	  continue;
	}
	eta(ind)=0.;
      }
    }
    
    else {
      
      //bounds on "all" variables
      for(Integer i=0;i<y.dim();i++){
	Real d=newy(i);
	Real b=lby(i);
	if (b>d){
	  Real e=weightu*(b-d);
	  eta(i)=e;
	  etaval+=b*e;
	  newy(i)=b;
	  normsubg2+=(e-2.*subg(i))*e;
	  continue;
	}
        b=uby(i);
	if (b<d){
	  Real e=weightu*(b-d);
	  eta(i)=e;
	  etaval+=b*e;
	  newy(i)=b;
	  normsubg2+=(e-2.*subg(i))*e;
	  continue;
	}
	eta(i)=0.;
      }	  
    } //end bounds on "all"
  } //end of else to do_scaling

  if (normsubg2<0)
    /* may happen because of numerical difficulties */
    normsubg2=0.;

  return 0;
}


// *****************************************************************************
//                                 solve_model
// *****************************************************************************

// loops over model optimization and lagrange updates till model optimizer
// is sufficiently close to a feasible model solution.

int NBundleSolver::solve_model()
{  
  int newy_computed=0;
  updatecnt=0;
  retcode=0;
  
  //--- after a serious step or change of u try to make a good guess for eta
  if ((bound_index.dim()>0)&&
      ((innerit==1)||(bundleweight->weight_changed()))){
    compute_newy(aggr_subg);
    newy_computed=1;
  }
  
  qpfails=0;         //is increased whenever the quadratic semidefinite program
                      //could not be solved to sufficient precision
  modelfails=0;

 //Things that are constant within the eta updating loop
 Matrix subg_diff=new_subg; subg_diff-=aggr_subg;
 for(Integer i=0;i<yfixed.dim();i++){
   if (yfixed(i)) {
     subg_diff(i)=0.;
   }
 }
 Real c2;
 if (do_scaling){
   c2=-ip(subg_diff,inv_scale%subg_diff)/weightu;
 }
 else {
   c2=-ip(subg_diff,subg_diff)/weightu;
 }
 
 //--- iterate till model is sufficiently precise
 do {
   
   updatecnt++;
   sumupdatecnt++;
   if ((out)&&(print_level>=1)){
     (*out)<<"  upd"<<updatecnt<<":";
   }

   //--- solve the quadratic model or update it
   oldaugval=max(oldaugval,augval);

   
   //coefficients of the quadratic c2/2*alpha^2+c1*alpha+c0
   Real c1;
   //Real c0;
   model_subg=aggr_subg;
   if (bound_index.dim()>0){
       model_subg-=eta;
   }   
   if (newy_computed){
     //c0=ip(model_subg,y)+etaval+aggrgamma-normsubg2/2./weightu;
     newy_computed=0;
   }
   else {
     newy.init(model_subg,-1./weightu);
     if (do_scaling) 
       newy%=inv_scale;
     //c0=ip(model_subg,y)+etaval+aggrgamma+ip(model_subg,newy)/2.;
     newy+=y;
   }
   c1=newgamma-aggrgamma+ip(subg_diff,newy);
   //find the alpha maximizing the concave quadratic  
   if (c2>-1e-6) 
     alpha=0.;
   else {
     alpha=-c1/c2;
     if (alpha>1)
       alpha=1.;
     if (alpha<0)
       alpha=0.;
   }
   //augval=c2*alpha*alpha/2.+c1*alpha+c0;

   modelgamma=alpha*newgamma+(1-alpha)*aggrgamma;
   xeyapzb(model_subg,new_subg,aggr_subg,alpha,1-alpha);

   //--- determine step

   compute_newy(model_subg);

   //--- compute the value of the linearized and the augmented model in newy 
   //--- (eta has no influence on linval by complementarity)

   linval=modelgamma+ip(newy,model_subg);
   augval=linval+normsubg2/2./weightu;
   cutval=max(newgamma+ip(newy,new_subg),aggrgamma+ip(newy,aggr_subg));
   modelprec=(cutval-linval)/max(oldval-linval,1e-16);


   //--- output some information on current iteration
   if ((out)&&(print_level>=1)){
     out->precision(10);(*out)<<"\n ";
     (*out)<<" augval="<<augval;out->precision(3);
     out->precision(8);(*out)<<" flin="<<linval;
     if (update_index.dim()) {
       out->precision(8);(*out)<<" fhat="<<cutval;
       (*out)<<" (";out->precision(2);
       (*out)<<modelprec<<")";
     }
     out->precision(2);
     (*out)<<" n2="<<normsubg2;
     (*out)<<" d2="<<norm2(newy-y);
     out->precision(6);(*out)<<std::endl;
   }
   
   //--- check for termination
   if (use_linval){
     modelval=linval;
   }
   else {
     modelval=cutval;
   }

   terminate=terminator->check_termination(this);
   
   if (terminate) break;     
   
 } while ( (modeleps>0) &&
	   (augval>oldaugval+eps_Real*fabs(oldaugval)) &&
           (cutval-linval>modeleps*(oldval-linval)) &&  //modelprecision
           ((max_updates<0)||(updatecnt<max_updates))
	   );
 
 return retcode;
}


// *****************************************************************************
//                                  init
// *****************************************************************************

 int NBundleSolver::init(NBundleProblem& prob,Integer in_dim,const Matrix* in_lb, const Matrix* in_ub,const Matrix* in_b)
{
  set_defaults();
  clear();
  problem=&prob;
  if (in_dim<0){
    if (out) (*out)<<"*** ERROR: NBundleSolver::init(): dimension<0"<<std::endl;
    return 1;
  }    
  if (in_dim==0) {
    if (out) (*out)<<"*** WARNING: NBundleSolver::init(): problem dimension is zero"<<std::endl;
  }
  if ((in_lb)&&(in_lb->dim()!=in_dim)) {
    if (out) (*out)<<"*** ERROR: NBundleSolver::init(): wrong dimension of lower bound vector"<<std::endl;
    return 1;
  }
  if ((in_ub)&&(in_ub->dim()!=in_dim)) {
    if (out) (*out)<<"*** ERROR: NBundleSolver::init(): wrong dimension of upper bound vector"<<std::endl;
    return 1;
  }
  if ((in_b)&&(in_b->dim()!=in_dim)) {
    if (out) (*out)<<"*** ERROR: NBundleSolver::init(): wrong dimension of cost vector"<<std::endl;
    return 1;
  }
  if ((in_ub)&&(in_lb)&&(min(*in_ub-*in_lb)<-1e-12)){
    if (out) (*out)<<"*** ERROR: NBundleSolver::init(): lower bounds exceed upper bounds by more than 1e-12"<<std::endl;
    return 1;
  }
  if (in_lb){
    lby=*in_lb;
    recompute_bound_index=true;
  }
  else {
    lby.init(in_dim,1,CB_minus_infinity);
  }
  if (in_ub){
    uby=*in_ub;
    recompute_bound_index=true;
  }
  else {
    uby.init(in_dim,1,CB_plus_infinity);
  }
  if (in_b){
    b=*in_b;
  }
  else {
    b.init(0,1,0.);
  }
  eta.init(in_dim,1,0.);
  y.init(lby.dim(),1,0.);
  initialize_center=true;
  return 0;
}


// *****************************************************************************
//                                 set_center
// *****************************************************************************


int NBundleSolver::set_center(const Matrix& iny)
{
  if ((iny.coldim()!=1)||(iny.dim()!=lby.dim())) {
    if (out) (*out)<<"*** ERROR: NBundleSolver::set_center(): the dimension of the center does not match the dimension of the boudns"<<std::endl;
    return 1;
  }
  y=iny;
  for(int i=0;i<y.dim();i++){
    if (y(i)>uby(i)){
      if ((out)&&(initialize_center==false))
	(*out)<<"*** WARNING: NBundleSolver::set_center(): coordinate "<<i<<" has value "<<iny(i)<<"and exceeds the upper bound"<< uby(i)<<", truncating ..." <<std::endl;
      y(i)=uby(i);
    }
    if (y(i)<lby(i)){
      if ((out) &&(initialize_center==false))
	(*out)<<"*** WARNING: NBundleSolver::set_center(): coordinate "<<i<<" has value "<<iny(i)<<"smaller than the lower bound"<< lby(i)<<", truncating ..." <<std::endl;
      y(i)=lby(i);
    }
  }
  newy=y;
  new_subg.init(0,1,0.);
  int status1=problem->eval_function(y,CB_plus_infinity,newgamma,oldval,&new_subg,1e-4);
  cntobjeval++;
  if (status1) {//no accurate solution, yet it may suffice to continue
    if (new_subg.dim()==y.dim()) {//no accurate solution, yet it may suffice to continue
      if (out) (*out)<<"*** WARNING: NBundleSolver::set_center(): function evaluation returned: "<<status1<<std::endl;
    }
    else {  //no new solution information at all, return with error
      if (out) (*out)<<"*** ERROR: NBundleSolver::set_center(): function evaluation failed: "<<status1<<std::endl;
      return 1;
    }
  }
  if (b.dim()>0){
    new_subg+=b;
    Real ipby=ip(b,y);
    newgamma+=ipby;
    oldval+=ipby;
  }

  newval=oldval;
  newgamma-=ip(new_subg,y);  
  center_subg=new_subg;
  
  problem->do_step();
  
  initialize_center=false; 
  return 0;
}


// *****************************************************************************
//                                 inner_loop
// *****************************************************************************

// performs null steps till a serious step is encountered, 
// a termination criterion is fulfilled or maxsteps is reached (if maxsteps>0 on input)

int NBundleSolver::inner_loop(int maxsteps)
{
  if (bundleweight) { 
    bundleweight->set_out(out, print_level-1);
  }
  else {
    if (out) (*out)<<"**** ERROR: NBundleSolver::inner_loop(...): no routine for choosing bundleweight specified"<<std::endl;
    return 1;
  }

  if (terminator) { 
    terminator->set_out(out, print_level-1);
  }
  else {
    if (out) (*out)<<"**** ERROR: NBundleSolver::inner_loop(...): no routine for checking termination specified"<<std::endl;
    return 1;
  }
  
 
  terminate=0;
  
  if (recompute_bound_index){
    bound_index.init(0,1,Integer(0));
    for(Integer i=0;i<lby.dim();i++){
      if ((lby(i)>CB_minus_infinity)||(uby(i)<CB_plus_infinity)){
	bound_index.concat_below(i);
      }
    }
    recompute_bound_index=false;
  }
  
  
  if (initialize_center){
    if (y.dim()==lby.dim()){
      if (set_center(y))
	return 1;
    }
    else {
      y.init(lby.dim(),1,0.);
      if (set_center(y))
	return 1;
    }
    initialize_model=true;
  }
  
  //--- check wether the aggregate has to be initialized
  if (initialize_aggregate){
    aggrgamma=newgamma;
    aggr_subg=new_subg;
    problem->init_aggregate();
    initialize_aggregate=false;
    recompute_normsubg2=true;
    initialize_model=true;
  }
  
  //--- do the first steps for initializing the model 
  if (initialize_model){
    
    //free all variables fixed to bounds
    if (yfixed.dim()>0){
      yfixed.init(0,1,Integer(0));
      recompute_normsubg2=true;
    }
    
    //set scaling
    if ((do_scaling)&&(do_scaling_heuristic)){
      Real average=0.;
      Integer nz=0;
      for(Integer i=0;i<y.dim();i++) {
	Real d=fabs(y(i));
	if (d>1e-6){
	  average+=fabs(d);
	  nz++;
	}
      }
      average/=Real(nz);
      inv_scale.init(y.dim(),1,1.);
      for(Integer i=0;i<inv_scale.dim();i++) {
	if (fabs(y(i))>10.*average) 
	  inv_scale(i)=sqr(fabs(y(i))/(10.*average));
      }
      recompute_normsubg2=true;
    }
  }
  
  //--- (re)initialize normsubg2 and bundleweight
  if (recompute_normsubg2){
    if (do_scaling){
      normsubg2=ip(aggr_subg,inv_scale%aggr_subg);
    }
    else {
      normsubg2=ip(aggr_subg,aggr_subg);
    }
    recompute_normsubg2=false;
  }
  if (weightu<0) bundleweight->init(aggr_subg,normsubg2);
  
  //--- complete initialization of the initial model
  if (initialize_model){
    
    //--- g(.)=linval+<subg,.-y> is a minorant of f contained in the model 
    //    compute a lower bound on the model and the augmented model value
    //    for linval and augval, respectively by means of g
    //    (minimizer of g(.)+weightu/2.*||.-y||^2 is  newy = y- subg/weightu)
    
    weightu=bundleweight->get_weight();
    linval=aggrgamma+ip(y,aggr_subg)-normsubg2/weightu; 
    oldaugval=augval=linval+0.999999*normsubg2/weightu/2; // correct value would be /2
    
    innerit=0;
    recomp=0;
    terminate=0;
    augvalfails=0;
    oraclefails=0;
    initialize_model=false;
  }
  
    
  //--- compute null steps till a serious step is achieved
  do{
        
    innerit++;
    suminnerit++;
      
    weightu=bundleweight->get_weight();
    lastaugval=augval;
      
    if ((out)&&(print_level>=1)){
      (*out)<<" ir"<<innerit;
      if (clockp) {(*out)<<"("<<*clockp<<")";}
      (*out)<<": u=";out->precision(4);(*out)<<weightu<<std::endl;
    }
    
    //--- solve quadratic model
    solve_model();
    
    if (augval<=lastaugval+eps_Real*fabs(lastaugval)){
      if ((out)&&(print_level>=1)){
	out->precision(12);
	(*out)<<"*** WARNING: NBundleSolver::inner_loop(): could not increase augval="<<augval<<" lastaugval="<<lastaugval<<"\n";
	(*out)<<"             maybe precision requirements are too high or\n";
	(*out)<<"             the weight was decreased instead of increased during consecutive null steps"<<std::endl;
	out->precision(4);
      }
      augvalfails++;
      sumaugvalfails++;
    }
    
    //--- reduce and aggregate bundle
    aggr_subg=model_subg;
    aggrgamma=modelgamma;
    problem->aggregate(alpha);

    if (terminate) {  //termination checked in solve_model
      break;
    }

    //---- determine nullstep_bound and descent_bound
    Real ipbnewy=0.;
    if (b.dim()>0){
      ipbnewy=ip(b,newy);
    }
    Real nullstep_bound = oldval-mN*(oldval-modelval);
    Real descent_bound  = oldval-mL*(oldval-modelval);
    
    //--- evaluate function at the candidate newy
    
    Real relprec=min(1e-3,.1*(oldval-nullstep_bound)/(fabs(oldval)+1.));
    
    Real newlb;
    int status1=problem->eval_function(newy,nullstep_bound-ipbnewy,newlb,newval,&new_subg,relprec);
    cntobjeval++;
    if (status1) {//no accurate solution, yet it may suffice to continue
      if (new_subg.dim()==y.dim()) {//no accurate solution, yet it may suffice to continue
	if (out) (*out)<<"*** WARNING: NBundleSolver::inner_loop(): function evaluation returned: "<<status1<<std::endl;
	oraclefails++;
	sumoraclefails++;
      }
      else {  //no new solution information at all, return with error
	if (out) (*out)<<"*** ERROR: NBundleSolver::inner_loop(): function evaluation failed: "<<status1<<std::endl;
	return 1;
      }
    }
    if (b.dim()>0){
      new_subg+=b;
      newlb+=ipbnewy;
      newval+=ipbnewy;
    }
    newgamma=newlb-ip(new_subg,newy);

    //--- check for potential problems with relative precision of all kinds 
    if (!status1){
      if (newlb<=descent_bound){ 
	if(newval>descent_bound){  //upper bound will produce a null step
	  Real vkappa=(oldval-newlb)/(oldval-modelval); //ratio of lower bound
	  if ((vkappa>max(.5,mN))||(oldval-newlb>.8*(oldval-cutval))) {
	    if (out){
	      (*out)<<"*** WARNING: NBundleSolver::inner_loop(): relative precision of returned objective interval =["<<newlb<<","<<newval<<"] enforces null step";
	      out->precision(3);(*out)<<" vkappa="<<vkappa;
	      (*out)<<" cutratio="<<(oldval-newval)/(oldval-cutval)<<std::endl;
	    }
	    dens++;
	  }
	  else if (3.*vkappa>4.*mN) {
	    sens++;
	    sensvkappasum+=vkappa;
	  }
	  else bens++;
	}
      }
      else { //lower bound gives already a null step
	if (oldval-newlb>.8*(oldval-cutval)){ 
	  //subgradient won't yield much improvement
	  if ((out)&&(print_level>=1)){
	    (*out)<<"  shallow cut (subgradient won't yield much improvement)";out->precision(3);
	  }
	  shallowcut++;
	}
      }
    }
    
    //--- output solution info
    if ((out)&&(print_level>=1)){
      (*out)<<"  nreval="<<cntobjeval;
      out->precision(8);(*out)<<" oldval="<<oldval;
      out->precision(8);(*out)<<" newval="<<newval;
      if (use_linval) {
	out->precision(8);(*out)<<" lin_model=";
      }
      else {
	out->precision(8);(*out)<<" full_model=";
      }
      (*out)<<modelval;
      out->precision(2);(*out)<<" model_prec="<<(oldval-newlb)/(oldval-modelval);
      out->precision(4);(*out)<<std::endl;
    }
    
    if (modelval>newval){  //minorant cuts above new function value! recompute!
      if (out) {
	out->precision(12);
	(*out)<<"*** WARNING: Bundle::inner_loop(): new objective value "<<newval<<" is smaller than minorizing model value ="<<modelval<<" \n";
	(*out)<<"***          maybe the objective value was not computed to sufficient precision\n";
	(*out)<<"***          or the sugbradient information provided by a previous oracle call was wrong.\n";
 	out->precision(4);
     }
    }

    //--- check for descent step; if so, move to this point and exit the loop
    if ((!status1)&&(newval<=descent_bound)){
      //descent step (serious step)
      bundleweight->descent_update(newval,oldval,modelval,y,newy,
				   normsubg2);
      problem->do_step();
      descent_steps++;
      y=newy;
      oldval=newval;
      center_subg=new_subg;
      initialize_model=true;
      break;
    }

      
    //--- check for validity of the old objective value
    //    by means of the new subgradient = new minorant

    Real lin_approx=oldval-newlb-ip(y-newy,new_subg);

    if (lin_approx<-1e-10*(fabs(oldval)+1.)){  //minorant cuts above old function value! recompute!
      if (out) {
	(*out)<<"*** WARNING: NBundleSolver::inner_loop(): new subgradient yields greater objective value ="<<oldval-lin_approx<<" in center (oldval="<<oldval<<"),\n";
	(*out)<<"***          maybe the old evaluation in the center was not computed to sufficient precision.\n";
        (*out)<<"***          Trying a function reevaluation in the center with higher precision..."<<std::endl;
      }
      recomp++;
      sumrecomp++;
      
      Real try_ubval;
      Real try_lbval;
      status1=problem->eval_function(y,CB_plus_infinity,try_lbval,try_ubval,0,relprec);
      cntobjeval++;
      if (status1) {//no accurate solution, yet it may suffice to continue
	if (out) (*out)<<"*** WARNING: NBundleSolver::inner_loop(): function reevaluation failed returning "<<status1<<", continuing anyway "<<status1<<std::endl;
      }
      if (b.dim()>0){
	model_subg+=b;
	Real ipby=ip(b,y);
	try_lbval+=ipby;
	try_ubval+=ipby;
      }

      if (try_ubval>oldval) {
	oldval=try_ubval;
	if (out) (*out)<<"***          reevaluation in center yields higher value="<<oldval<<std::endl;
	center_subg=model_subg;
      }
      else {
	if (out) (*out)<<"***          reevaluation in center yields no larger value ("<<try_ubval<<")"<<std::endl;
      }
    }
    
    //--- choose a new weight (greater or equal the current one) 
    
    bundleweight->nullstep_update(newval,oldval,modelval,lin_approx,
				  y,newy,normsubg2/weightu);

    //--- fix some box constraint variables if allowed to do so
    if ((yfixing_allowed)&&
	(modelprec>yfixing_factor*modeleps)&&
	(innerit>=yfixing_itbound)){
      Real fixbound=0.;
      Integer cnt=0;
      {for(Integer i=0;i<eta.dim();i++){
	if (eta(i)!=0.){
	  fixbound+=fabs(eta(i));
	  cnt++;
	}
      }}
      fixbound=fixbound/Real(cnt)/2.;
      if (yfixed.dim()==0)
	yfixed.init(y.dim(),1,Integer(0));
      if (out&&(print_level>=2)) (*out)<<" fixing ";
      {for(Integer i=0;i<y.dim();i++){
	if (yfixed(i)){
	  if (out&&(print_level>=2)) (*out)<<" F"<<i;
	}
	else if ((y(i)==lby(i))&&(eta(i)>fixbound)){
	  eta(i)=0.;
	  yfixed(i)=1;
	  if (out&&(print_level>=2)) (*out)<<" L"<<i;
	}
	else if ((y(i)==uby(i))&&(eta(i)<-fixbound)){
	  eta(i)=0.;
	  yfixed(i)=1;
	  if (out&&(print_level>=2)) (*out)<<" U"<<i;
	}
	else if ((y(i)>lby(i)+1e-6)&&(newy(i)==lby(i))&&(eta(i)>fixbound)){
	  if (out&&(print_level>=2)) (*out)<<" #"<<i;
	}
	else if ((y(i)<uby(i)-1e-6)&&(newy(i)==uby(i))&&(eta(i)<-fixbound)){
	  if (out&&(print_level>=2)) (*out)<<" *"<<i;
	}	  
      }}
      if (sum(yfixed)==0){
	yfixed.init(0,1,Integer(0));
      }
      if (out&&(print_level>=2)) (*out)<<" sum="<<sum(yfixed)<<std::endl;
      //yfixed.init(0,1,Integer(0));
    }
			    
    
    
 }while(--maxsteps);

 return 0;
} 

// *****************************************************************************
//                                 set_lower_bound
// *****************************************************************************


int NBundleSolver::set_lower_bound(Integer in_i, Real in_lb)
{
  if ((in_i<0)||(in_i>=lby.dim())){
    if (out) 
      (*out)<<"*** ERROR: NBundleSolver::set_lower_bound: index "<<in_i<<" out of range [0,"<<lby.dim()<<"]"<<std::endl;
    return 1;
  }
  if ((initialize_center==false)&&(in_lb>y(in_i)))
    initialize_center=true;
  lby(in_i)=in_lb;
  recompute_bound_index=true;
  return 0;
}

// *****************************************************************************
//                                 set_upper_bound
// *****************************************************************************


int NBundleSolver::set_upper_bound(Integer in_i, Real in_ub)
{
  if ((in_i<0)||(in_i>=uby.dim())){
    if (out) 
      (*out)<<"*** ERROR: NBundleSolver::set_upper_bound: index "<<in_i<<" out of range [0,"<<uby.dim()<<"]"<<std::endl;
    return 1;
  }
  if ((initialize_center==false)&&(in_ub<y(in_i)))
    initialize_center=true;
  uby(in_i)=in_ub;
  recompute_bound_index=true;
  return 0;
}

// *****************************************************************************
//                                append_variables
// *****************************************************************************


 int NBundleSolver::append_variables(Integer n_append,const Matrix* lbp,const Matrix *ubp,const Matrix *bp)
{
  if (n_append==0) return 0;
#if (CONICBUNDLE_DEBUG>=1)
  if (n_append<0){
    if (out) 
      (*out)<<"*** ERROR: NBundleSolver::append_variables: n_append="<<n_append<<" negative"<<std::endl;
    return 1;
  }
  if (lbp!=0){
    if (lbp->dim()!=n_append){
      if (out) 
	(*out)<<"*** ERROR: NBundleSolver::append_variables: n_append="<<n_append<<" but dimension of lower bounds is "<<lbp->dim()<<std::endl;
      return 1;
    }
    if ((min(*lbp)<CB_minus_infinity)||(max(*lbp)>CB_plus_infinity)){
      if (out) 
	(*out)<<"*** ERROR: NBundleSolver::append_variables: lower bounds not within [CB_minus_infinity,CB_plus_infinity]"<<std::endl;
      return 1;
    }
  }
  if (ubp!=0){
    if (ubp->dim()!=n_append){
      if (out) 
	(*out)<<"*** ERROR: NBundleSolver::append_variables: n_append="<<n_append<<" but dimension of upper bounds is "<<ubp->dim()<<std::endl;
      return 1;
    }
    if ((min(*ubp)<CB_minus_infinity)||(max(*ubp)>CB_plus_infinity)){
      if (out) 
	(*out)<<"*** ERROR: NBundleSolver::append_variables: upper bounds not within [CB_minus_infinity,CB_plus_infinity]"<<std::endl;
      return 1;
    }
  }
  if ((lbp!=0)&&(ubp!=0)){
    for(Integer i=0;i<n_append;i++){
      if ((*lbp)(i)>(*ubp)(i)){
	if (out) 
	  (*out)<<"*** ERROR: NBundleSolver::append_variables: for i="<<i<<" lower_bound(i)="<<(*lbp)(i)<<">"<<(*ubp)(i)<<"=upper_bound(i)"<<std::endl;
	return 1;
      }
    }
  }
  if ((bp!=0)&&(bp->dim()>0)&&(bp->dim()!=n_append)){
    if (out) 
      (*out)<<"*** ERROR: NBundleSolver::append_variables: n_append="<<n_append<<" but dimension of cost vector is "<<bp->dim()<<std::endl;
    return 1;
  }
#endif 
  Real max_newlb=CB_minus_infinity;

  eta.concat_below(Matrix(n_append,1,0.));

  if ((bp==0)||(bp->dim()==0)){
    if (b.dim()>0)
      b.concat_below(Matrix(n_append,1,0.));
  }
  else {
    if (b.dim()==0)
      b.init(lby.dim(),1,0.);
    b.concat_below(*bp);
  }

  if (lbp==0){
    lby.concat_below(Matrix(n_append,1,CB_minus_infinity));
  }
  else {
    lby.concat_below(*lbp);
    recompute_bound_index=true;
    max_newlb=max(*lbp);
  }
  Real min_newub=CB_plus_infinity;

  if (ubp==0){
    uby.concat_below(Matrix(n_append,1,CB_plus_infinity));
  }
  else {
    uby.concat_below(*ubp);
    recompute_bound_index=true;
    min_newub=min(*ubp);
  }

  if (yfixed.dim()>0){
    yfixed.concat_below(Indexmatrix(n_append,1,Integer(0)));
  }

  if (y.dim()>0){
    y.concat_below(Matrix(n_append,1,0.));
    if (do_scaling){
      inv_scale.concat_below(Matrix(n_append,1,1.));
    }
    
    if (initialize_center==false){
      if ((max_newlb<=0)&&(min_newub>=0)){
	if(new_subg.dim()>0){
	  Matrix append_aggr;
	  Matrix append_subg;
	  initialize_center=(problem->subgradient_information(Indexmatrix(Range(y.dim()-n_append,y.dim()-1)),append_aggr,append_subg)!=0);
	  if (!initialize_center){
	    if ((bp!=0)&&(bp->dim()>0)){
	      append_aggr+=*bp;
	      append_subg+=*bp;
	    }
	    aggr_subg.concat_below(append_aggr);
	    new_subg.concat_below(append_subg);
	  }
	}
	else {
	  initialize_center=true;
	}
      }
      else {
	initialize_center=true;
      }
    }
    
    if (initialize_center==false){
      recompute_normsubg2=true;
    }
    else {
      aggr_subg.init(0,1,0.);
      problem->clear_aggregate();
      new_subg.init(0,1,0.);
      center_subg.init(0,1,0.);
      initialize_aggregate=true;
    }
    
  }
  
  return 0;
}

// *****************************************************************************
//                               reassign_variables
// *****************************************************************************


int NBundleSolver::reassign_variables(const Indexmatrix& map_to_old)
{
  if (map_to_old.dim()==0){
    y.init(0,1,0.);
    lby.init(0,1,0.);
    uby.init(0,1,0.);
    b.init(0,1,0.);
    bound_index.init(0,1,Integer(0)); 
    newy.init(0,1,0.);
    center_subg.init(0,1,0.);
    aggr_subg.init(0,1,0.);
    new_subg.init(0,1,0.);
    eta.init(0,1,0.);
    return 0;
  }
  if (max(map_to_old)>=lby.dim()){
    if (out) 
      (*out)<<"*** ERROR: NBundleSolver::reassign_variables: maximum index exceeds current dimension"<<std::endl;
    return 1;
  }
  if (y.dim()==lby.dim()){
    y=y(map_to_old);
  }
  if (center_subg.dim()==lby.dim()){
    center_subg=center_subg(map_to_old);
  }
  if (new_subg.dim()==lby.dim()){
    new_subg=new_subg(map_to_old);
  }
  if ((do_scaling)&&(inv_scale.dim()==lby.dim())){
    inv_scale=inv_scale(map_to_old);
  }
  if (aggr_subg.dim()==lby.dim()){
    aggr_subg=aggr_subg(map_to_old);
    recompute_normsubg2=true;
  }    
  if (eta.dim()==lby.dim()){
    eta=eta(map_to_old);
  }  
  if (yfixed.dim()==lby.dim()){
    yfixed=yfixed(map_to_old);
  }
  if (b.dim()>0)
    b=b(map_to_old);
  lby=lby(map_to_old);
  uby=uby(map_to_old);
  recompute_bound_index=true;
  return 0;
}



// *****************************************************************************
//                             reinit_function_model
// *****************************************************************************


int NBundleSolver::reinit_function_model()
{
  initialize_center=true;
  aggr_subg.init(0,1,0.);
  problem->clear_aggregate();
  initialize_aggregate=true;
  return 0;
}


// *****************************************************************************
//                               clear_aggregates
// *****************************************************************************


int NBundleSolver::clear_aggregate()
{
  aggr_subg.init(0,1,0.);
  problem->clear_aggregate();
  initialize_aggregate=true;
  return 0;
}


// *****************************************************************************
//                                 print_line_summary
// *****************************************************************************

// outputs most important data within one line

std::ostream& NBundleSolver::print_line_summary(std::ostream& o) const
{
 o.setf(std::ios::showpoint);
 o.fill(' ');
 if (clockp) o<<*clockp;
 if (!terminate) o<<" endit ";
 else o<<" _endit ";
 o.width(2);o<<descent_steps;o.precision(3);
 o<<" ";o.width(3);o<<suminnerit;
 o<<" ";o.width(3);o<<sumupdatecnt;
 o<<" ";o.width(6);o<<weightu;
 o<<" ";o.width(7);o<<sqrt(normsubg2);
 o<<" ";o.precision(8);o.width(10);o<<modelval;
 o<<" ";o.width(10);o<<oldval<<std::endl;
 return o;
}

}
