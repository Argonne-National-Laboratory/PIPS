/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/bundle.cxx

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



#include <cmath>
#include <cstdlib>
#include "bundle.hxx"
#include "hkweight.hxx"
#include "hkweight.hxx"
#include "idscaling.hxx"
#include "diagscaling.hxx"
#include "mymath.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {


// *****************************************************************************
//                                 BundleSolver
// *****************************************************************************

BundleSolver::BundleSolver()
{
  terminator=new BundleTerminator;
  bundleweight=new BundleHKweight(.5);
  Hp=0;
  problem=0; 
  out=0; 
  print_level=1;
  clockp=0; 
 
  set_defaults();
  clear();
}

// *****************************************************************************
//                                 ~BundleSolver
// *****************************************************************************

BundleSolver::~BundleSolver()
{
  delete terminator; 
  delete bundleweight;
  delete Hp;
}

// *****************************************************************************
//                                 set_defaults
// *****************************************************************************

// usually called on construction

void BundleSolver::set_defaults()
{
  modeleps=0.6; 
  mL=0.1;   
  mN=0.1;

  yfixing_allowed=false;
  yfixing_itbound=10;
  yfixing_factor=.5;

  do_scaling_heuristic=0;
  delete Hp;
  Hp=new BundleIdScaling;
  Hp->set_out(out,print_level-1);
  
  max_updates=10; //was 100;
  use_linval=1;
  if (terminator) terminator->set_defaults();
  if (bundleweight) bundleweight->set_defaults();
}

// *****************************************************************************
//                                   clear
// *****************************************************************************

// usually called on construction, resets all matrices and calls set_defaults

void BundleSolver::clear()
{
  problem=0;       
  terminate=0;
  weightu=-1;      
  y.init(0,1,0.);
  lby.init(0,1,0.);
  uby.init(0,1,0.);
  bound_index.init(0,1,Integer(0)); 
  newy.init(0,1,0.);
  subg.init(0,1,0.);
  eta.init(0,1,0.);
  normsubg2=-1;
  yfixed.init(0,0,Integer(0)); 
  update_index.init(0,1,Integer(0));
  update_value.init(0,1,0.);
  updatecnt=0;
  sumupdatecnt=0;
  retcode=0;
  problem_changed=0;
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

void BundleSolver::clear_fails()
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
//                                   set_do_scaling
// *****************************************************************************

void BundleSolver::set_do_scaling(int ds)
{
 do_scaling_heuristic=(ds!=0); 
 if (ds==0) { 
   delete Hp; 
   Hp=new BundleIdScaling; 
 }
}

// *****************************************************************************
//                                   get_do_scaling
// *****************************************************************************

int BundleSolver::get_do_scaling(void) const
{
  return (dynamic_cast<BundleIdScaling*>(Hp)==0);
}

// *****************************************************************************
//                                   set_scaling
// *****************************************************************************

void BundleSolver::set_scaling(const CH_Matrix_Classes::Matrix& insc)
{
  delete Hp;
  Matrix tmp(insc);
  tmp.inv();
  BundleDiagonalScaling* dp=new BundleDiagonalScaling(tmp);
  Hp=dp;
  Hp->set_out(out,print_level-1);
  recompute_normsubg2=true;
}

// *****************************************************************************
//                                   set_scaling
// *****************************************************************************

void BundleSolver::set_quadratic_term(BundleScaling* Sp)
{
  if (Sp!=0){
    delete Hp;
    Hp=Sp;
    Hp->set_out(out,print_level-1);
    recompute_normsubg2=true;
  }
  else {
    if (out)
      (*out)<<"*** WARNING BundleSolver::set_scaling(BundleScaling*) was called with 0 pointer"<<std::endl;
  }
}

// *****************************************************************************
//                                 solve_model
// *****************************************************************************

// loops over model optimization and lagrange updates till model optimizer
// is sufficiently close to a feasible model solution.

int BundleSolver::solve_model()
{     
  updatecnt=0;
  retcode=0;
  qpfails=0;         //is increased whenever the quadratic bundle subproblem
                     //could not be solved to sufficient precision
  modelfails=0;
 
  //--- iterate till model is sufficiently precise
  do {
   
    updatecnt++;
    sumupdatecnt++;
    if ((out)&&(print_level>=1)){
      (*out)<<"  upd"<<updatecnt<<":";
    }
    
    //--- solve the quadratic model or update it
    oldaugval=max(oldaugval,augval);   
    int errcode;
    if (updatecnt==1){ 
      if (problem_changed)
	errcode=problem->eval_augmodel(y,lby,uby,eta,oldaugval,oldval,1.e-3,Hp,yfixed);
      else
	errcode=problem->eval_augmodel(y,lby,uby,eta,oldaugval,oldval,0.1,Hp,yfixed);
    }
    else {
      errcode=problem->reeval_augmodel(y,lby,uby,eta,update_index,update_value,oldaugval,oldval,0.1,Hp);
    }
    
    retcode=retcode || errcode;
    
    if (errcode){
      if (out) (*out)<<"*** WARNING BundleSolver::solve_model(): solving quadraticmodel failed "<<errcode<<std::endl;
      qpfails++;
      sumqpfails++;
    }
    
    //--- get model value at y=0 together with the aggregate subgradient
    //    (the subgradient and linval do not include the influence of eta!)
    if (problem->get_augmodel_sol(linval,subg)){
      if (out) (*out)<<"*** ERROR BundleSolver::solve_model(): no augmented model solution available"<<std::endl;
      return 1;
    }
    problem_changed=0;
    
    //--- determine step
    if (Hp->update_eta_step(newy,eta,update_index,update_value,normsubg2,subg,y,lby,uby,bound_index,yfixed)){
      if (out) 
	(*out)<<"*** ERROR BundleSolver::solve_modle(): update_eta_step(...) failed"<<std::endl;
      return 1;
    }

    //--- compute the value of the linearized and the augmented model in newy 
    //--- (eta has no influence on linval by complementarity)
    linval+=ip(newy,subg);
    augval=linval+normsubg2/2.;
    
    //--- if there were changes in newy compute the true model value
    
    if ((update_index.dim())&&(modeleps>0)){
      if ((out)&&(print_level>=1)){ 
	(*out)<<" changed="<<update_index.dim()<<"("<<norm2(update_value)<<")";
      }
      int modelerr=problem->eval_model(newy,linval+modeleps*(oldval-linval),0.1*(oldval-linval)/(fabs(oldval)+1.));
      if (modelerr){
	if (out) (*out)<<"*** WARNING BundleSolver::solve_model(): evaluation of cutting model failed: "<<errcode<<std::endl;
	modelfails++;
	summodelfails++;
      }
      Real dummy;
      problem->get_model_sol(cutval,dummy,0);
    }
    else {
      cutval=linval;
    }
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
//                                 inner_loop
// *****************************************************************************

// performs null steps till a serious step is encountered

int BundleSolver::inner_loop(BundleProblem& prob, int maxsteps)
{
  problem=&prob;

  if (bundleweight) { 
    bundleweight->set_out(out, print_level-1);
  }
  else {
    if (out) (*out)<<"**** ERROR: BundleSolver::inner_loop(...): no routine for choosing bundleweight specified"<<std::endl;
    return 1;
  }

  if (terminator) { 
    terminator->set_out(out, print_level-1);
  }
  else {
    if (out) (*out)<<"**** ERROR: BundleSolver::inner_loop(...): no routine for checking termination specified"<<std::endl;
    return 1;
  }
  
 
  terminate=0;
  
  //--- check whether problem or y needs reinitialization
  
  problem_changed=(weightu>0)||(y.dim()!=0);
  if (problem->init_center(problem_changed,linval,oldval,y,subg,bound_index,lby,uby)){
    if (out) (*out)<<"*** ERROR: BundleSolver::inner_loop(): retrieving information on center of stability failed";
    return 1;
  }
  if (problem_changed){
    eta.init(y.dim(),1,0.); 
    recompute_normsubg2=true;
    initialize_model=true;
  }
    
    
  //--- do the first steps for initializing the model 
  if (initialize_model){

    //free all variables fixed to bounds
    if (yfixed.dim()>0){
      yfixed.init(0,0,Integer(0));
      recompute_normsubg2=true;
    }

    //set scaling
    if (do_scaling_heuristic){
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
      Matrix tmp(y.dim(),1,1.);
      for(Integer i=0;i<y.dim();i++) {
	if (fabs(y(i))>10.*average) 
	  tmp(i)=10.*average/sqr(fabs(y(i)));
      }
      BundleDiagonalScaling* dp=dynamic_cast<BundleDiagonalScaling*>(Hp);
      if (dp==0){
	delete Hp;
	Hp=dp=new BundleDiagonalScaling(tmp);
      }
      else {
	dp->set_D(tmp);
      }
      recompute_normsubg2=true;
    }


  }
  
  //--- (re)initialize normsubg2 and bundleweight
  bool first_eta_computed=false;
  if (weightu<0) {
    if (bound_index.dim()>0){
      Hp->set_weightu(1.);
      if(Hp->compute_first_eta(eta,subg,y,lby,uby,bound_index,yfixed)){
	if (out)
	  (*out)<<"*** ERROR BundleSolver::inner_loop(): compute_first_eta(...) failed"<<std::endl;
	return 1;
      }
      Matrix tmp(subg);
      tmp-=eta;
      bundleweight->init(tmp);
      if (fabs(bundleweight->get_weight()-1.)<1e-10)
	first_eta_computed=true;
    }
    else 
      bundleweight->init(subg);
    initialize_model=true;
    recompute_normsubg2=true;
  }
  if ((fabs(weightu-bundleweight->get_weight())>1e-10*fabs(weightu))||
      (recompute_normsubg2)){
    weightu=bundleweight->get_weight();
    Hp->set_weightu(weightu);
    if (bound_index.dim()>0){
      if (!first_eta_computed){
	if (Hp->compute_first_eta(eta,subg,y,lby,uby,bound_index,yfixed)){
	  if (out)
	    (*out)<<"*** ERROR BundleSolver::inner_loop(): compute_first_eta(...) failed"<<std::endl;
	  return 1;
	}
	first_eta_computed=true;
      }
      Matrix tmp(subg);
      tmp-=eta;
      normsubg2=Hp->dnorm_sqr(tmp);
    }
    else 
      normsubg2=Hp->dnorm_sqr(subg);
    recompute_normsubg2=false;
    initialize_model=true;
  }
  
  //--- complete initialization of the initial model
  if (initialize_model){
    
    //--- g(.)=linval+<subg-eta,.-y> is a minorant of f contained in the model 
    //    compute a lower bound on the model and the augmented model value
    //    for linval and augval, respectively by means of g
    //    (minimizer of g(.)+||.-y||_H^2/2. is  newy = y- H^{-1}subg)
    
    modelval=cutval=linval=linval-normsubg2; 
    oldaugval=augval=linval+0.999999*normsubg2/2.; // correct value would be /2
    
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
    Hp->set_weightu(weightu);
    lastaugval=augval;
    
    if ((out)&&(print_level>=1)){
      (*out)<<" ir"<<innerit;
      if (clockp) {(*out)<<"("<<*clockp<<")";}
      (*out)<<": u=";out->precision(4);(*out)<<weightu<<std::endl;
    }
    
    //--- after a serious step or change of u try to make a good guess for eta
    if ((!first_eta_computed)&&(bound_index.dim()>0)&&
	((innerit==1)||(bundleweight->weight_changed()))){
      if(Hp->compute_first_eta(eta,subg,y,lby,uby,bound_index,yfixed)){
	if (out)
	  (*out)<<"*** ERROR BundleSolver::inner_loop(): compute_first_eta(...) failed"<<std::endl;
	return 1;
      }
    }
    first_eta_computed=false;
  
    //--- solve quadratic model
    solve_model();
    
    if (augval<=lastaugval+eps_Real*fabs(lastaugval)){
      if ((out)&&(print_level>=1)){
	out->precision(12);
	(*out)<<"*** WARNING: Bundle::inner_loop(): could not increase augval="<<augval<<" lastaugval="<<lastaugval<<"\n";
	(*out)<<"             maybe precision requirements are too high or\n";
	(*out)<<"             the weight was decreased instead of increased during consecutive null steps"<<std::endl;
	out->precision(4);
      }
      augvalfails++;
      sumaugvalfails++;
    }
    
    if (terminate) {  //termination checked in solve_model
      break;
    }
    
    //---- determine nullstep_bound and descent_bound
    Real nullstep_bound = oldval-mN*(oldval-modelval);
    Real descent_bound  = oldval-mL*(oldval-modelval);
    
    //--- evaluate function at the candidate newy
    
    Real relprec=min(1e-3,.1*(oldval-nullstep_bound)/(fabs(oldval)+1.));
    
    int status1=problem->eval_function(newy,nullstep_bound,relprec);
    cntobjeval++;
    if (status1) {//no accurate solution, yet it may suffice to continue
      if (out) (*out)<<"*** WARNING: Bundle::inner_loop(): function evaluation returned: "<<status1<<std::endl;
      oraclefails++;
      sumoraclefails++;
    }
    
    Real newlb;
    int status2 = problem->get_function_sol(newlb,newval,0);

    if (status2) {  //no new solution information at all, return with error
      if (out) (*out)<<"*** ERROR: Bundle::inner_loop(): retrieving solution data of function evaluation failed: "<<status2<<std::endl;
      return 1;
    }
    
    //--- check for potential problems with relative precision of all kinds 
    if (!status1){
      if (newlb<=descent_bound){ 
	if(newval>descent_bound){  //upper bound will produce a null step
	  Real vkappa=(oldval-newlb)/(oldval-modelval); //ratio of lower bound
	  if ((vkappa>max(.5,mN))||(oldval-newlb>.8*(oldval-cutval))) {
	    if (out){
	      (*out)<<"*** WARNING: Bundle::inner_loop(): relative precision of returned objective interval =["<<newlb<<","<<newval<<"] enforces null step";
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
      problem->update_model(true); 
      bundleweight->descent_update(newval,oldval,modelval,y,newy,
				   normsubg2);
      problem->do_step();
      descent_steps++;
      y=newy;
      oldval=newval;
      initialize_model=true;
      break;
    }
    else {
      //null step, update the model
      problem->update_model(false); 
    }

      
    //--- check for validity of the old objective value
    //    by means of the new subgradient = new minorant

    Matrix new_subg;
    status2 = problem->get_function_sol(newlb,newval,&new_subg);

    if (status2) {  //no new solution information at all, return with error
      if (out) (*out)<<"*** ERROR: Bundle::inner_loop(): retrieving new subgradient failed: "<<status2<<std::endl;
      return 1;
    }
      
    Real lin_approx=oldval-newlb-ip(y-newy,new_subg);

    if (lin_approx<-1e-10*(fabs(oldval)+1.)){  //minorant cuts above old function value! recompute!
      if (out) {
	(*out)<<"*** WARNING: Bundle::inner_loop(): new subgradient yields greater objective value ="<<oldval-lin_approx<<" in center (oldval="<<oldval<<"),\n";
	(*out)<<"***          maybe the old evaluation in the center was not computed to sufficient precision\n";
	(*out)<<"***          or the sugbradient information provided by the last oracle call was wrong.\n";
        (*out)<<"***          Trying a function reevaluation in the center with higher precision..."<<std::endl;
      }
      recomp++;
      sumrecomp++;
      
      status1=problem->eval_function(y,CB_plus_infinity,.001*relprec);
      cntobjeval++;
      if (status1) {
	if (out) (*out)<<"*** WARNING: Bundle::inner_loop(): reevaluation returned error code: "<<status1<<std::endl;
      }
      
      Real try_ubval;
      Real try_lbval;
      status2 = problem->get_function_sol(try_lbval,try_ubval,0);
      
      if (try_ubval>oldval) {
	oldval=try_ubval;
	problem->do_step();    
	if (out) (*out)<<"***          reevaluation in center yields higher oldval="<<oldval<<std::endl;         
      }
      else {
	if (out) (*out)<<"***          reevaluation in center yields no better value ("<<try_ubval<<")"<<std::endl;
      }
    }
    
    //--- choose a new weight (greater or equal the current one) 
    
    bundleweight->nullstep_update(newval,oldval,modelval,lin_approx,
				  y,newy,normsubg2);

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
    }
			    
    
    
 }while(--maxsteps);

 return 0;
} 


// *****************************************************************************
//                                 print_line_summary
// *****************************************************************************

// outputs most important data within one line

std::ostream& BundleSolver::print_line_summary(std::ostream& o) const
{
  o.setf(std::ios::showpoint);
 o.fill(' ');
 if (clockp) o<<*clockp;
 if (!terminate) o<<" endit ";
 else o<<" _endit ";
 o.width(3);o<<descent_steps;o.precision(3);
 o<<" ";o.width(4);o<<suminnerit;
 o<<" ";o.width(4);o<<sumupdatecnt;
 o<<" ";o.width(6);o<<weightu;
 o<<" ";o.width(7);o<<sqrt(normsubg2);
 o<<" ";o.precision(10);o.width(12);o<<modelval;
 o<<" ";o.width(12);o<<oldval<<std::endl;
 return o;
}


}
