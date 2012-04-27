/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/Nbundle.hxx

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



#ifndef CONICBUNDLE_NBUNDLE_HXX
#define CONICBUNDLE_NBUNDLE_HXX

#include "matrix.hxx"
#include "clock.hxx"
#include "CBSolver.hxx"
#include "bundle.hxx"


namespace ConicBundle {


class NBundleProblem
{

public:
  virtual ~NBundleProblem(){}

  //----------  function evaluation
  virtual int eval_function(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real nullstep_bound,CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp,CH_Matrix_Classes::Real relprec) =0;
  //evaluates the objective function in $y$
  //if evaluated by an iterative method that provides upper and lower bounds,
  //it may stop when the lower bound (lb) is above the nullstep_bound
  //evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied
  //it returns a lower and an upper bound on the objective value.
  //If subgp==0, then no subgradient needs to be returned.
  //if subgp is not 0 then an eps-subgradient at y is returned,
  //the hyperplane corresponding to the subgradient has value lb in y.
  //If subgp!=0 but it returns a zero matrix, then the subgradient is invalid 
  //returns:  0 ... if all is ok, use get_function_sol or do_step
  //          1 ... if solution could not be computed to desired precision

  virtual void do_step(void)=0;
  //this is called whenever the bundle method moves on to the last
  //point evaluated. This may be used to memorize the corresponding primal
  //information in Lagrangean relaxation if needed

  virtual void init_aggregate(void)=0;
  //this is called whenever the bundle method initializes its aggregate
  //to the most recent subgradient computed (the last call to eval_function
  //with subgp!=0).
  //Nothing has to be done by the routine, but it may be used to
  //generate an approximate primal in Lagrangean relaxation if needed

  virtual void aggregate(CH_Matrix_Classes::Real alpha)=0;
  //this is called whenever the bundle method determined a new aggregate as
  //new_aggregate = alpha newsubgradient + (1-alpha) old_aggregate
  //(newsubgradient refers to the last call to eval_function with subgp!=0)
  //nothing has to be done by the routine, but it may be used to
  //generate an approximate primal in Lagrangean relaxation if needed
  
  virtual void clear_aggregate(void)=0;
  //this is called whenever the bundle method removes its aggregate.
  //Nothing has to be done by the routine, but it may be used to
  //generate an approximate primal in Lagrangean relaxation if needed

  virtual int subgradient_information(const CH_Matrix_Classes::Indexmatrix& indices,CH_Matrix_Classes::Matrix& missing_aggr_coords,CH_Matrix_Classes::Matrix& missing_new_subg_coords)=0;
  //this is called whenever the bundle method has to append new variables.
  //If the BundleProblem is able to provide the missing coordinates of
  //the aggregate and the subgradient of the most recent evaluation call 
  //(e.g. because a corresponding primal is known in Lagrangian relaxation),  
  //then it should do so (missing_...[i] has to hold the value of the 
  //subgradient coordinate indices[i]) and return 0,
  // otherwise it should return 1. 

};



  class NBundleSolver: public BundleTerminatorData
{
private:
  //-----------------------------------------------------------------------
  //       Variables and Parameters modified by the algorithm (inner_loop)
  //-----------------------------------------------------------------------
  NBundleProblem *problem;       
  
  int terminate;               //value of last call to terminator
  //see BundleTerminator for values
  
  CH_Matrix_Classes::Real weightu;                //weight u of quadratic term / inverse step size
  
  CH_Matrix_Classes::Matrix y,                     //center of stability
    lby,                        //lower bounds on variables if not free
                                //for free variable use CB_minus_infinity
    uby,                        //upper bounds on variables if not free
                                //for free variable use CB_plus_infinity
    b;                          //additional linear costs on y to be used only if b.dim()!=0;
   
  CH_Matrix_Classes::Indexmatrix bound_index;      //indices of bounded variables
  
  CH_Matrix_Classes::Matrix newy,                  //candidate
    center_subg,                     //subgradient obtained in the center
    new_subg,                        //last subgradient
    aggr_subg,                       //aggregate subgradient 
    model_subg,                       //aggregate subgradient of the model 
    eta;                        //Lagrange multipliers for box constraints

  CH_Matrix_Classes::Real normsubg2;    //=<subg,subg>, or if (do_scaling) =<subg,inv_scale%subg>
                     //wher subg is either the initial subgradient
                     //             or the best subgradient of the aug. model

  CH_Matrix_Classes::Matrix inv_scale;             //diagonal scaling by  1/inv_scale(i);

  CH_Matrix_Classes::Indexmatrix yfixed;           //if some coordinates should be kept fixed for some
                                //time, set the corresponding coordinates to one,
                                //otherwise have yfixed empty or zero

  //for storing updates to lagrange multipliers
  CH_Matrix_Classes::Indexmatrix update_index;
  CH_Matrix_Classes::Matrix update_value;

  //variables in solve_model needed for resuming interrupted optimizations
  CH_Matrix_Classes::Integer updatecnt;
  CH_Matrix_Classes::Integer sumupdatecnt;
  int retcode;
  bool recompute_bound_index;

  bool initialize_center;
  bool initialize_aggregate;
  bool recompute_normsubg2;
  bool initialize_model;

  //objective values
  CH_Matrix_Classes::Real oldval; //function value in the center
  CH_Matrix_Classes::Real newval; //function value in the new point (the most recent candidate)
  CH_Matrix_Classes::Real aggrgamma; //function value of the aggregate supporting hyperplane in 0
  CH_Matrix_Classes::Real newgamma; //function value of the new supporting hyperplane in 0
  CH_Matrix_Classes::Real modelgamma; //function value of the supporting model hyperplane in 0
  CH_Matrix_Classes::Real etaval; //value induced by the bound multipliers
  CH_Matrix_Classes::Real linval; //value of linearized model for feasible (pojected) solution
  CH_Matrix_Classes::Real cutval; //value of cutting surface model for feasible (projected) solution
  CH_Matrix_Classes::Real modelval;
  CH_Matrix_Classes::Real augval;
  CH_Matrix_Classes::Real oldaugval;
  CH_Matrix_Classes::Real lastaugval;
  CH_Matrix_Classes::Real alpha;

  //counters 
  CH_Matrix_Classes::Integer innerit;    //number of nullsteps
  CH_Matrix_Classes::Integer suminnerit; //sum of innerit over all calls
  CH_Matrix_Classes::Integer cntobjeval; //number of calls to eigenvalue computation routine
                      //may differ from suminnerit because of additiona
                      //corrections
  CH_Matrix_Classes::Integer recomp;     //recomputations for the current oldval
  CH_Matrix_Classes::Integer sumrecomp;  //sum over all recomputations
  CH_Matrix_Classes::Integer qpfails;      //counts number of fails in calls to eval_aug_model
  CH_Matrix_Classes::Integer sumqpfails;   //sums up number of fails
  CH_Matrix_Classes::Integer modelfails;   //counts number of fails in calls to eval_model
  CH_Matrix_Classes::Integer summodelfails;//sums up number of fails
  CH_Matrix_Classes::Integer augvalfails;  //counts, within a descent step, null steps
                        //that fail to increase augval 
                        //!(augval>lastaugval+eps_Real*fabs(lastaugval))
  CH_Matrix_Classes::Integer sumaugvalfails;  //sum over all augval failures
  CH_Matrix_Classes::Integer oraclefails;   //counts number of fails in calls to the oracle
  CH_Matrix_Classes::Integer sumoraclefails;//sums up number of fails


  CH_Matrix_Classes::Integer cntsolves;
  CH_Matrix_Classes::Integer cntresolves;
  CH_Matrix_Classes::Integer cntitersolves;
  CH_Matrix_Classes::Integer cntiterresolves;
  CH_Matrix_Classes::Integer descent_steps;        //counts number of calls to inner_loop

  CH_Matrix_Classes::Integer bens,       //number of benign enforced null steps
    sens,             //number of significant enforced null steps
    dens,             //number of dangerous enforced null steps
    shallowcut;       //number of null steps with shallow cut
  CH_Matrix_Classes::Real sensvkappasum;
  CH_Matrix_Classes::Real modelprec;

  //--------------------------------------------------------
  //  parameters that are not modified by the algorithm
  //---------------------------------------------------------

  //----  the following parameters are all reset by set_defaults

  CH_Matrix_Classes::Real modeleps;       //fixed parameter for the precision of the model
  CH_Matrix_Classes::Real mL;             //acceptance factor
  CH_Matrix_Classes::Real mN;             //nullstep factor for inexact evaluation
  int use_linval;      //or use cutval for acceptance crit. default=1

  bool yfixing_allowed;         //sets variables to their bounds if eta updates difficult
  CH_Matrix_Classes::Integer yfixing_itbound;  //do it after so many null steps
  CH_Matrix_Classes::Real yfixing_factor;      //if modelprec is less than factor*modeleps

  int do_scaling;               //scaling is used only if do_scaling=1
  int do_scaling_heuristic;     //set inv_scale=1/sqr(y(i)) for abs(y(i))>1.

  CH_Matrix_Classes::Integer max_updates;  //do at most this number of updates of eta in each step
    

  //----  the following parameters are not influenced by set_defaults

  BundleTerminator *terminator; //pointer to termination criterion routine
                                //(terminator itself may get modified)
  BundleWeight *bundleweight;  //pointer to routine for choosing weight u
                                //(bundleweight itself may get modified)

  //for output
  std::ostream *out;
  int print_level;
  const CH_Tools::Clock* clockp;  

  //--------------------------------------------------------
  //  private subroutines
  //---------------------------------------------------------

  int compute_newy(const CH_Matrix_Classes::Matrix& subg);

  int update_eta(void);

  int solve_model(void);
    
  //--------------------------------------------------------
  //  public functions
  //---------------------------------------------------------
public:
  NBundleSolver();
  ~NBundleSolver();
        
  void set_defaults();  //resets all parameters to default values
                        //and calls set_defaults for terminator and bundleweight
  void clear();         //resets all variables modified by the algorithm 
                        //and calls clear() for terminator and bundleweight

  void clear_fails();   //resets all fail counts to zero

  int init(NBundleProblem& prob,CH_Matrix_Classes::Integer dim,const CH_Matrix_Classes::Matrix* in_lb, const CH_Matrix_Classes::Matrix* in_ub, const CH_Matrix_Classes::Matrix* in_b);
  int set_center(const CH_Matrix_Classes::Matrix& iny);
  int set_lower_bound(CH_Matrix_Classes::Integer in_i, CH_Matrix_Classes::Real in_lb);
  int set_upper_bound(CH_Matrix_Classes::Integer in_i, CH_Matrix_Classes::Real in_ub);
  int append_variables(CH_Matrix_Classes::Integer n_append,const CH_Matrix_Classes::Matrix* lbp,const CH_Matrix_Classes::Matrix*ubp, const CH_Matrix_Classes::Matrix* bp);
  int reassign_variables(const CH_Matrix_Classes::Indexmatrix& map_to_old);
  int reinit_function_model();
  int clear_aggregate();

  inline void set_terminator(BundleTerminator* bt);
  void set_bundleweight(BundleWeight* bw)
    {delete bundleweight; bundleweight=bw;}
  void set_modeleps(CH_Matrix_Classes::Real in_eps){modeleps=in_eps;}
  void set_mL(CH_Matrix_Classes::Real in_mL){mL=in_mL;}
  void set_mN(CH_Matrix_Classes::Real in_mN){mN=in_mN;}
  void set_use_linval(int i){use_linval=i;}
  void set_do_yfixing(bool dofix,CH_Matrix_Classes::Integer nullstepbound=10,
		      CH_Matrix_Classes::Real precfactor=.5)
    {yfixing_allowed=dofix;yfixing_itbound=nullstepbound;yfixing_factor=precfactor;} 
  void set_do_scaling(int ds){do_scaling=ds;do_scaling_heuristic=1;}
  void set_scaling(const CH_Matrix_Classes::Matrix& insc)
    {inv_scale=insc;do_scaling=1;do_scaling_heuristic=0;}
  void set_clock(const CH_Tools::Clock& myclock){clockp=&myclock;}
  void set_out(std::ostream* o=0,int pril=1){out=o;print_level=pril;}
  void set_max_updates(CH_Matrix_Classes::Integer mu){max_updates=mu;}
    
  int get_terminate() const {return terminate;}
  CH_Matrix_Classes::Real get_modeleps() const {return modeleps;}
  int get_do_scaling() const {return do_scaling;}
  int get_do_heuristic_scalig() const {return do_scaling_heuristic;}
    
  CH_Matrix_Classes::Integer get_cntobjeval() const {return cntobjeval;}
  CH_Matrix_Classes::Integer get_descent_steps() const {return descent_steps;}
  CH_Matrix_Classes::Integer get_innerit() const {return innerit;}
  CH_Matrix_Classes::Integer get_suminnerit() const {return suminnerit;}
  CH_Matrix_Classes::Integer get_sumupdatecnt() const {return sumupdatecnt;}
  CH_Matrix_Classes::Integer get_recomp() const {return recomp;}
  CH_Matrix_Classes::Integer get_sumrecomp() const {return sumrecomp;} 
  CH_Matrix_Classes::Integer get_qpfails() const {return qpfails;}
  CH_Matrix_Classes::Integer get_sumqpfails() const {return sumqpfails;}
  CH_Matrix_Classes::Integer get_modelfails() const {return modelfails;}
  CH_Matrix_Classes::Integer get_summodelfails() const {return summodelfails;}
  CH_Matrix_Classes::Integer get_augvalfails() const {return augvalfails;}
  CH_Matrix_Classes::Integer get_sumaugvalfails() const {return sumaugvalfails;}
  CH_Matrix_Classes::Integer get_oraclefails() const {return oraclefails;}
  CH_Matrix_Classes::Integer get_sumoraclefails() const {return sumoraclefails;}
  
  CH_Matrix_Classes::Integer get_cntsolves() const {return cntsolves;}
  CH_Matrix_Classes::Integer get_cntresolves() const {return cntresolves;}
  CH_Matrix_Classes::Integer get_cntitersolves() const {return cntitersolves;}
  CH_Matrix_Classes::Integer get_cntiterresolves() const {return cntiterresolves;}
  CH_Matrix_Classes::Integer get_bens() const {return bens;}
  CH_Matrix_Classes::Integer get_sens() const {return sens;}
  CH_Matrix_Classes::Integer get_dens() const {return dens;}
  CH_Matrix_Classes::Integer get_shallowcut() const {return shallowcut;}
  CH_Matrix_Classes::Real get_sensvkappasum() const{return sensvkappasum;}
  CH_Matrix_Classes::Integer get_calls() const {return descent_steps;}
    
  CH_Matrix_Classes::Real get_oldval() const {return oldval;}
  CH_Matrix_Classes::Real get_newval() const {return newval;}
  CH_Matrix_Classes::Real get_modelval() const {return modelval;}
  int get_use_linval() const {return use_linval;}
  CH_Matrix_Classes::Real get_weight() const {return weightu;}
  CH_Matrix_Classes::Real get_normsubg2() const {return normsubg2;}
  void get_eta(CH_Matrix_Classes::Matrix& outeta) const {outeta=eta;} 
  const CH_Matrix_Classes::Matrix& get_y() const {return y;} 
  const CH_Matrix_Classes::Matrix& get_newy() const {return newy;} 
  const CH_Matrix_Classes::Matrix& get_aggr() const {return aggr_subg;} 
  const CH_Matrix_Classes::Matrix& get_lby() const {return lby;} 
  const CH_Matrix_Classes::Matrix& get_uby() const {return uby;} 
  const CH_Matrix_Classes::Matrix& get_eta() const {return eta;} 
  const CH_Matrix_Classes::Indexmatrix& get_yfixed() const {return yfixed;}
  
  const NBundleProblem* get_problem() const {return problem;}
  BundleTerminator* get_terminator() const {return terminator;}
  BundleWeight*  get_bundleweight() const {return bundleweight;}
  
  int inner_loop(int maxsteps=0);
  
  std::ostream& print_line_summary(std::ostream& out) const;
  
};





//-----------------------------------------------------------------------------
//----                              inline
//-----------------------------------------------------------------------------

inline void NBundleSolver::set_terminator(BundleTerminator* sbt)
{delete terminator; terminator=sbt;}

}

#endif

