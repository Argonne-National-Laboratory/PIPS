/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/bundle.hxx

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



#ifndef CONICBUNDLE_BUNDLE_HXX
#define CONICBUNDLE_BUNDLE_HXX


/**  @file bundle.hxx
    @brief Header declaring the classes ConicBundle::BundleScaling, ConicBundle::BundleProblem, 
         ConicBundle::BundleWeight, ConicBundle::BundleTerminatorData, ConicBundle::BundleSolver, 
         ConicBundle::BundleTerminator
    @version 1.0
    @date 2005-04-21
    @author Christoph Helmberg
*/


#include "matrix.hxx"
#include "clock.hxx"
#include "CBSolver.hxx"

namespace ConicBundle {

/** @defgroup InternalBundleSolver Internal Bundle Solver of ConicBundle
   @brief Solve \f$min_{y\in\mathbf{R}^m}  f(y)\f$ 
    for a convex function f given in the form of a ConicBundle::BundleProblem.
    The y-variables may be bounded or box constrained. 
    The Problem and the Solver are set up and called from
    one of the solver interfaces, e.g., ConicBundle::CBSolver or 
    ConicBundle::MatCBSolver and are typically not accessed directly.

*/
//@{

/** @brief abstract interface for obtaining the bundle data needed in the routine
           BundleScaling::compute_QP_costs
 */
class BundleQPData
{
public:
  virtual ~BundleQPData(){}

  /** @brief *assign* the coefficients corresponding to coordinate index_y of y for all 
    model variables to the positions row(startindex,...,startindex+model_dim-1) 
    and *add* a constant y-coefficient (independent of the model variables) to b,
    index_y==-1 is used for the cost coefficients and constant term.
    Returns 0 on success, 1 on failure
  */

  virtual int get_row(CH_Matrix_Classes::Integer index_y,CH_Matrix_Classes::Matrix& row,CH_Matrix_Classes::Real& b,CH_Matrix_Classes::Integer startindex) const=0;
};

/** @brief abstract interface that allows to use different \f$H\f$-norms \f$\|y-\hat{y}\|_H^2\f$ 
     with a positive definite matrix \f$H\f$ in the augmented model of ConicBundle::BundleSolver.
     ConicBundle::BundleNoScaling is the identitiy \f$H=I\f$. There are further variants 
     for diagonal scaling, low rank scaling or full scaling by a positive definite matrix.      
 */

class BundleScaling
{
protected:
  //for output
  std::ostream *out;
  int print_level;

public:
  BundleScaling(void){out=0;print_level=0;}
  virtual ~BundleScaling(){}

  virtual void set_weightu(CH_Matrix_Classes::Real weightu)=0;
  virtual CH_Matrix_Classes::Real get_weightu() const=0;
   
  /// returns \f$\|B\|^2_H\f$
  virtual CH_Matrix_Classes::Real norm_sqr(const CH_Matrix_Classes::Matrix& B) const =0;
  /// returns \f$\|B\|^2_{H^{-1}}\f$
  virtual CH_Matrix_Classes::Real dnorm_sqr(const CH_Matrix_Classes::Matrix& B) const =0;

  /** @brief computes initial Lagrange multipliers eta for given (aggregate) 
      subgradient subg and center y so that with this eta the next candidate newy 
      would lie within lower bounds lby and upper bounds uby on the support bounds_index 
      and eta and newy would satisfy complementary slackness (i.e., both are optimal).
      
  */
  virtual int compute_first_eta(CH_Matrix_Classes::Matrix& eta,
				const CH_Matrix_Classes::Matrix& subg,
				const CH_Matrix_Classes::Matrix& y,
				const CH_Matrix_Classes::Matrix& lby,
				const CH_Matrix_Classes::Matrix& uby,
				const CH_Matrix_Classes::Indexmatrix& bound_index,
				const CH_Matrix_Classes::Indexmatrix& yfixed) const =0;

  /** @brief computes Lagrange multiplier eta and next candidate newy for 
      given (aggregate) subgradient subg and center y so that newy is within 
      lower bounds lby and upper bounds uby on the support bounds_index and
      eta and newy satisfy complementary slackness (i.e., both are optimal).
      If yfixed.dim()==y.dim() then in all coordinates  i with yfixed(i)!=0
      y cannot change and so these coordinates should be ignored.
      In addition, changes to the previous eta (supplied on input) are stored 
      in sparse form in update_value and update_index and subgnorm2 is
      set to \f$\|newy-y\|_H^2\f$.
  */
  virtual int update_eta_step(CH_Matrix_Classes::Matrix& newy,
			      CH_Matrix_Classes::Matrix& eta,
			      CH_Matrix_Classes::Indexmatrix& update_index,
			      CH_Matrix_Classes::Matrix& update_value,
			      CH_Matrix_Classes::Real& subgnorm2,			      			      const CH_Matrix_Classes::Matrix& subg,
			      const CH_Matrix_Classes::Matrix& y,
			      const CH_Matrix_Classes::Matrix& lby,
			      const CH_Matrix_Classes::Matrix& uby,
			      const CH_Matrix_Classes::Indexmatrix& bound_index,
			      const CH_Matrix_Classes::Indexmatrix& yfixed) const =0;


  /** @brief computes the dual QP costs Q, d, and the constant offset to the bundle subproblem

     The bundle subproblem problem with Lagrangean relaxiation of box constraints on y reads  
     \f[ \min_y \max_{x\in X} (b-\eta-Ax)^\top y + ub^\top \eta^{-}+lb^\top \eta^{+}+c^\top x+\frac{u}{2}\|y-\hat y\|_H^2,\f]
     where \f$ b\f$ is a linear cost term for \f$ y\f$, \f$\eta\f$ are the Lagrange multipliers for the box
     constraints on y and, e.g., \f$ A\f$ contains xdim subgradients \f$s_i\f$ in
     its columns, \f$c_i\f$ give the height of the supporting hyperplane
     corresponding to \f$s_i\f$ for \f$y=0\f$, \f$x\in X\f$ describes convex
     combinations and \f$H\f$ is the quadratic term. Then the dual reads
     \f[ \max_{x\in X}-\frac1{2u}x^\top Qx+d^\top x+offset \f]
     with \f$Q=A^TH^{-1}A\f$, d=... 
  */
  virtual int compute_QP_costs(CH_Matrix_Classes::Symmatrix& Q,
			       CH_Matrix_Classes::Matrix& d,
			       CH_Matrix_Classes::Real& offset,
                               const BundleQPData* datap,
			       const CH_Matrix_Classes::Integer xdim,
			       const CH_Matrix_Classes::Matrix& y,
			       const CH_Matrix_Classes::Matrix& lby,
			       const CH_Matrix_Classes::Matrix& uby,
			       const CH_Matrix_Classes::Matrix& eta,
			       const CH_Matrix_Classes::Indexmatrix& yfixed) const =0;

  

  /** @brief computes the update of the dual QP cost terms d and offset returned by
      compute_QP_costs() for changes in eta given by update_value and update_index,
      more precisely, eta[update_index[i]]=old_eta[update_index[i]]+update_value[i].
  */
  virtual int update_QP_costs(CH_Matrix_Classes::Matrix& delta_d,  
			      CH_Matrix_Classes::Real& delta_offset,
			      const BundleQPData* datap,
			      const CH_Matrix_Classes::Integer xdim,
			      const CH_Matrix_Classes::Matrix& y,
			      const CH_Matrix_Classes::Matrix& lby,
			      const CH_Matrix_Classes::Matrix& uby,
			      const CH_Matrix_Classes::Matrix& eta,
			      const CH_Matrix_Classes::Indexmatrix& update_index,
			      const CH_Matrix_Classes::Matrix& update_value) const =0;

  /** @brief controls output-levle, if out==0 output nothing at all, if out!=0 and 
      print_level<=0 output warnings and errors only, if out!=0 and print_level>=1 
      this may generate more and more log or debug information 
  */
  virtual void set_out(std::ostream* in_out=0,int in_print_level=1)
  {out=in_out;print_level=in_print_level;}

};


  /** @brief abstract interface for BundleSolver to the convex optimization problem and
     the subproblem solvers. It gives access to all problem specific bundle routines. 
     In particular it hides the cutting modle, the quadratic bundle subproblem, and 
     the oracle. So here all the real work is done.
 */
class BundleProblem:public BundleQPData
{

public:
  virtual ~BundleProblem(){} 

  /** @brief returns information about the starting point/current center of stability y; 
    
    @param[in,out] init
       - On input, init=0 asks to initialize all variables regardless of content,
       and init=1 asks to initialize only if the problem was modified since the 
       last call or the last aggregate subgradient is no longer valid due to 
       external modfications in the cutting model.
       - On Output init_flag=0 signals "no variables modified", =1 means "all were modified".

    @param[out] lb
      gives a lower bound on the function value in y; together with subg it defines
      a linear minorant that is part of the bundle (the bundle solver may use this
      to compute a lower bound for the next quadratic bundle subproblem). lb need
      not be the best lower bound known.

    @param[out] ub gives an upper bound on the objective value in y, this should
      be as close as possible to the true objective value in y

    @param[out] y
      is the new center of stability

    @param[out] subg
      subg is an (arbitrary) eps-subgradient in y that is part of the bundle; 
      lb and subg belong together (i.e., lb+<subg,.> is a minorant),
      (lb and subg do not have to be the same as in eval_function!)
   
    @param[out] bounds_index 
       holds the indices of y that are bounded (dim=#bounded) 
   
    @param[out] lby 
       gives lower bounds on y 
       (if lower or upper bounds exist at all then lby.dim() = y.dim())

    @param[out] uby
       gives upper bounds on y 
       (if lower or upper bounds exist at all then uby.dim() = y.dim())

    @return 
        - 0 ... if the information is available
        - 1 ... if the desired information is not available
  */
  virtual int init_center(int& init, 
			  CH_Matrix_Classes::Real& lb, 
			  CH_Matrix_Classes::Real& ub, 
			  CH_Matrix_Classes::Matrix& y, 
			  CH_Matrix_Classes::Matrix& subg,
			  CH_Matrix_Classes::Indexmatrix& bounds_index, 
			  CH_Matrix_Classes::Matrix& lby, 
			  CH_Matrix_Classes::Matrix& uby) =0;

  /** @brief evaluates the objective function in $y$, results are then returned in .

    If evaluated by an iterative method that provides upper (ub) and lower bounds (lb),
    the method may stop when the lower bound is above the nullstep_bound.
    Evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied

    @param[in] y 
        is the point to evaluate at.

    @param[in] nullstep_bound
        if the function value is above this value, a null step will be made
        and above this bound any approximate function value suffices to ensure
        convergence
  
    @param[in] relprec
        if the nullstep_bound is not reached, then evalutation may 
        stop, if (ub-lb)<relprec*(|ub|+1) is satisfied


    @return
        - 0 ... if all is ok, use get_function_sol() or do_step()
        - 1 ... if solution could not be computed to desired precision
  */
  virtual int eval_function(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real nullstep_bound,CH_Matrix_Classes::Real relprec) =0;

  /** @brief returns lower and upper bound on the objective value of the last call to eval_function()
     @param[out] lb
         lower bound on the objective value
    
     @param[out] ub
         upper bound on the objective value
 
     @param[in,out] subgp 
         if not 0 on input, then an  eps-subgradient at y is returned,
         the hyperplane corresponding to this eps-subgradient has value lb in y
     
    @return
        - 0 ... if the information is available 
                (if eval_function() returned 1, the information will not
                satisfy the precision requirements)
        - 1 ... if the desired information is not available
  */
     
  virtual int get_function_sol(CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp) =0;


  /** @brief store y of last eval_function, its subgradient, and objective value 
       as new center of stability.
  */
  virtual int do_step(void)=0;

  /** @brief evaluate the current cutting model in $y$ 
   
    If evaluated by an iterative method that provides upper and lower bounds
    it may stop when the lower bound (lb) is above eval_model_bound .
    Evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied
  */
  virtual int eval_model(
    const CH_Matrix_Classes::Matrix& y,  ///< evaluate the model in this point 
    CH_Matrix_Classes::Real eval_model_bound, ///< no need for exact value above this bound 
    CH_Matrix_Classes::Real relprec ///< otherwise use this relative precision 
  )=0;


  /** @brief return lower and upper bound on objective value of last eval_function()
    
     @param[out] lb
         lower bound on the objective value
    
     @param[out] ub
         upper bound on the objective value
 
     @param[in,out] subgp 
         if not 0 on input, then an  eps-subgradient to the model at y is returned,
         the hyperplane corresponding to this eps-subgradient has value lb in y
     
    @return
        - 0 ... if the information is available 
                (if eval_model() returned 1, the information will not
                satisfy the precision requirements)
        - 1 ... if the desired information is not available
  */
  virtual int get_model_sol(CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp)=0;
	     
  /** @brief Evaluates the augmented model with respect to the center of stability.

  If yfixed is a vector of the same dimension as the center, it 
  must give zeros for indices the may be changed and ones if they should stay fixed.
  If evaluated by an iterative method that provides upper and lower bounds,
  it should stop when  ub-lb <= min(relprec*(fbound-lb),(ub-augbound)/2.)
  If the subproblem needs to adapt some penalty parameter and therefore
  needs to modify the function value it can do so by changing augbound
  and fbound (=the old function value) correspondingly (caution!)
  
     @param[in] center_y   center of stability
     
     @param[in] lby       lower bounds on y
     
     @param[in] uby       upper bounds on y
    
     @param[in] eta       current multipliers for the bound constraints on y (primal slacks)
   
     @param[in,out] augbound the value of the augemented model should be at least this
           (in rare cases of models that adapt penalty parameters, the value may be
            adapted internally and is then returned)

     @param[in,out] fbound the value of the augemented model cannot exceed this
           (in rare cases of models that adapt penalty parameters, the value may be
            adapted internally and is then returned)

     @param[in] relprec relative precision requirements
   
     @param[in] Hp provides the routines for computing the cost matrices dependent on
              the choice of the quadratic term

     @param[in] yfixed allows to specify coordinates of y that may not be changed

    @return
        - 0 ... if all is ok, use get_augmodel_sol() for retrieving the solution
        - -1 ... if setting up the QP failed
        - otherwise it returns the satus returned by the internal QP_Solver

  */

  virtual int eval_augmodel(const CH_Matrix_Classes::Matrix& center_y,
			    const CH_Matrix_Classes::Matrix& lby,
			    const CH_Matrix_Classes::Matrix& uby,
			    const CH_Matrix_Classes::Matrix& eta,
			    CH_Matrix_Classes::Real& augbound,
			    CH_Matrix_Classes::Real& fbound,
			    CH_Matrix_Classes::Real relprec,
			    const BundleScaling *Hp, 
			    const CH_Matrix_Classes::Indexmatrix& yfixed)=0;
	     
  /** @brief reevaluate the augmented model for updated eta w.r.t. the previously called eval_augmodel()
  
    The values of eta-old_eta are listed in update_index and update_value.
    The scaling and the weight u must be the same in both augmodel-calls!
    If yfixed is a vector of the same dimension as the center it 
    must give zeros for indices the may be changed and ones if they should stay fixed.
    If evaluated by an iterative method that provides upper and lower bounds,
    it should stop when  ub-lb <= min(relprec*(fbound-lb),(ub-augbound)/2.).
    If the subproblem needs to adapt some penalty parameter and therefore
    needs to modify the function value it can do so by changing augbound
    and fbound (=the old function value) correspondingly (caution!)
     @param[in] center_y   center of stability
     
     @param[in] lby       lower bounds on y
     
     @param[in] uby       upper bounds on y
    
     @param[in] eta       current multipliers for the bound constraints on y (primal slacks)

     @param[in] update_index lists the indices of eta that were changed, none of them in yfixed of eval_augmodel() 
     
     @param[in] update_value for index i,eta[update_index[i]]=old_eta[update_index[i]]+update_value[i] 

     @param[in,out] augbound the value of the augemented model should be at least this
           (in rare cases of models that adapt penalty parameters, the value may be
            adapted internally and is then returned)

     @param[in,out] fbound the value of the augemented model cannot exceed this
           (in rare cases of models that adapt penalty parameters, the value may be
            adapted internally and is then returned)

     @param[in] relprec relative precision requirements
   
     @param[in] Hp provides the routines for computing the cost matrices dependent on
              the choice of the quadratic term

    @return
        - 0 ... if all is ok, use get_augmodel_sol() for retrieving the solution
        - -1 ... if setting up the QP failed
        - otherwise it returns the satus returned by the internal QP_Solver

  */

  virtual int reeval_augmodel(const CH_Matrix_Classes::Matrix& center_y,
			      const CH_Matrix_Classes::Matrix& lby,
			      const CH_Matrix_Classes::Matrix& uby,
			      const CH_Matrix_Classes::Matrix& eta,
			      const CH_Matrix_Classes::Indexmatrix& update_index,
			      const CH_Matrix_Classes::Matrix& update_value,
			      CH_Matrix_Classes::Real& augbound,
			      CH_Matrix_Classes::Real& fbound,
			      CH_Matrix_Classes::Real relprec,
			      const BundleScaling *Hp)=0;


  /** @brief returns the last solution of (re)eval_augmodel 

    The linear function linconst+ip(subg,.) is a minorant of the objective
    and the minimizer over the corresponding augmented model 
    newy:= argmin linval+ip(subg-eta,.)+weightu/2*||.-y||^2
    is an approximate minimizer of the entire augmented model.

  @param[out] linconst constant of the aggregate linear minorant
 
  @param[out] subg gradient of the aggregate linear minorant, subgradient of the model in newy

  @return 
    - 0 if the information is available 
                  (if eval_augmodel() did not return 0, the information will not
                   satisfy the precision requirements)
    - 1 if the desired information is not available
  */
  virtual int get_augmodel_sol(CH_Matrix_Classes::Real& linconst,CH_Matrix_Classes::Matrix& subg) =0;

  /** @brief generate the next cutting model. 

    If descent_step is false, the next model has to contain at least 
    all convex combinations of the two current subgradients of eval_function 
    and (re)eval_augmodel

    @param[in] descent_step
        If set to true, then the next step is a descent step. In this
        case the model may be restarted from scratch and need not include
        the last aggregate.

    @return
      - 0 ... if the requirements on the model could be met and no errors occured
      - 1 ... otherwise
  */
  virtual int update_model(bool descent_step) =0;
  
};


/** @brief Abstract interface for BundleSolver providing routines that determine 
     the weight of the quadratic term in the augmented model. It also allows the
     user to specify bounds on the weights by setting minweight and maxweight or
     even to choose the weight to be used in the next iteration.
*/
class BundleWeight
{
public:
  BundleWeight(){}
  virtual ~BundleWeight(){}

  /// Sets default values for 'constant' parameters, e.g. minweight and maxweight
  virtual void set_defaults()=0; 

  /// Resets all adaptive variables and parameters
  virtual void clear()=0;        
  
  /** @brief Compute the first weight u and set some parameters, 
     @param[in] subgrad the gradient of a (local) linear minorant (typically the initial or the aggregate subgradient)
     @param[in] normsubg2 the norm squared of the subgradient (dual norm to the quadratic term in the augmented model)
     @return 
     - 0 on success
     - 1 on failure
  */
  virtual int init(const CH_Matrix_Classes::Matrix& subgrad,CH_Matrix_Classes::Real normsubg2)=0;
  
  /** @brief Compute the first weight u and set some parameters, 
     @param[in] subgrad the gradient of a (local) linear minorant (typically the initial or the aggregate subgradient)
     @return 
     - 0 on success
     - 1 on failure
  */
  virtual int init(const CH_Matrix_Classes::Matrix& subgrad)=0;
  
  /** @brief this routine has to be called after variables were added or deleted
     @param[in] subgrad the gradient of a (local) linear minorant (typically the initial or the aggregate subgradient)
     @param[in] normsubg2 the norm squared of the subgradient (dual norm to the quadratic term in the augmented model)
     @return 
     - 0 on success
     - 1 on failure
  */  
  virtual int prob_changed(const CH_Matrix_Classes::Matrix& subgrad,CH_Matrix_Classes::Real normsubg2)=0;

  /** @brief Allows the user to set the weight of the very next iteration of the bundle method

  As it is set by the user, the weight will be accepted as long as it is strictly positive,
  even if it is outside the interval [minweight,maxweight]. Nonpositive values, however,
  will leave everything unchangend.
  */
  virtual void set_next_weight(CH_Matrix_Classes::Real u)=0;

  /** @brief Sets a lower bound for the weight. Nonpositive values may be used to indicate that the weight is allowed to get arbitrarily close to zero.
   */   
  virtual void set_minweight(CH_Matrix_Classes::Real umin) =0;

  /** @brief Sets an upper bound for the weight. Nonpositive values may be used to indicate that the weight is allowed to get arbitrarily large.
   */   
  virtual void set_maxweight(CH_Matrix_Classes::Real umax) =0;
  
  
  /// Returns the current value of the weight 
  virtual CH_Matrix_Classes::Real get_weight() const =0;
  
  /// Returns 0 only if the last call to descent_update() or nullstep_update() did not modify the weight and no other routine influenced the weight since then, otherwise it returns 1.  
  virtual int weight_changed() const =0;
  
  /** @brief The BundleSolver calls this for computing the next weight if if the candidate will result in a descent step

    @param[in] newval objective value in the candidate (next center)

    @param[in] oldval objective value in the current (old) center

    @param[in] modelval objective value of the current model in the candidate

    @param[in] y current (old) center    

    @param[in] newy candidate (next center)    

    @param[in] normsubg2 squared norm of the aggregate subgradient that gave rise to newy

     @return 
     - 0 on success
     - 1 on failure    
  */
  virtual int descent_update(
			     CH_Matrix_Classes::Real newval,CH_Matrix_Classes::Real oldval,CH_Matrix_Classes::Real modelval,
			     const CH_Matrix_Classes::Matrix& y, const CH_Matrix_Classes::Matrix& newy,
			     CH_Matrix_Classes::Real normsubg2)=0;
  
  /** @brief The BundleSolver calls this for computing the next weight if if the candidate will result in a null step

    @param[in] newval objective value in the candidate 

    @param[in] oldval objective value in the current center

    @param[in] modelval objective value of the current model in the candidate

    @param[in] y current center    

    @param[in] newy candidate    

    @param[in] normsubg2 squared norm of the aggregate subgradient that gave rise to newy

    @return 
     - 0 on success
     - 1 on failure    
  */
  virtual int nullstep_update(
			      CH_Matrix_Classes::Real newval,CH_Matrix_Classes::Real oldval,CH_Matrix_Classes::Real modelval,CH_Matrix_Classes::Real lin_approx,
			      const CH_Matrix_Classes::Matrix& y, const CH_Matrix_Classes::Matrix& newy,
			      CH_Matrix_Classes::Real norm2subg)=0;
  
  /** @brief Specifies the output level (out==NULL: no output at all, 
           out!=NULL and level=0: errors and warnings, 
           level>0 increasingly detailed information)

     @param[in] out  (std::ostream*) 
       direct all output to (*out). If out==NULL, there will be no output at all.

     @param[in] print_level (int)
  */
  virtual void set_out(std::ostream* out=0,int print_level=1)=0;
  
};

class BundleTerminator;


  /** @brief abstract interface for BundleTerminator providing the data needed for deciding on termination
   */
class BundleTerminatorData
{
public:
  virtual ~BundleTerminatorData(){}

  /// returns the number of calls to the oracle, i.e., to BundleProblem::eval_function()
  virtual CH_Matrix_Classes::Integer get_cntobjeval() const =0;
  /// returns the number of reevalutions in center points (if their function value violates the most recent subgradient inequality)
  virtual CH_Matrix_Classes::Integer get_sumrecomp() const =0;
  /// returns the number of descent and null steps
  virtual CH_Matrix_Classes::Integer get_suminnerit() const =0;
  /// returns the number of times, the call to BundleProblem::eval_augmodel() or BundleProblem::reeval_augmodel() failed
  virtual CH_Matrix_Classes::Integer get_sumqpfails() const =0;
  /// returns the function value in the current center of stability
  virtual CH_Matrix_Classes::Real get_oldval() const =0;
  /// returns the function value in the current candidate
  virtual CH_Matrix_Classes::Real get_newval() const =0;
  /// returns the model value in the current candidate
  virtual CH_Matrix_Classes::Real get_modelval() const =0;
  /// returns the norm of the aggregate subgradient squared (for the norm dual to quadratic augmented term) 
  virtual CH_Matrix_Classes::Real get_normsubg2() const =0;
  /// returns the number of failed callls to BundleProblem::eval_model() 
  virtual CH_Matrix_Classes::Integer get_summodelfails() const =0;
  /// returns the number of times, the augmented model value could not be increased since the last descent step
  virtual CH_Matrix_Classes::Integer get_augvalfails() const =0;
  /// returns the number of times, the augmented model value could not be increased over all descent and null steps
  virtual CH_Matrix_Classes::Integer get_sumaugvalfails() const =0;
  /// returns the number of times, the oracle returend some error code since the last descent step
  virtual CH_Matrix_Classes::Integer get_oraclefails() const =0;
  /// returns the number of times, the oracle returend some error code over all descent and null steps
  virtual CH_Matrix_Classes::Integer get_sumoraclefails() const =0;

};

/** @brief This is the internal bundle solver managing descent/null steps 
    and inner updates of the multipliers for box constraints, etc.
*/


class BundleSolver: public BundleTerminatorData
{
private:
  //-----------------------------------------------------------------------
  //       Variables and Parameters modified by the algorithm (inner_loop)
  //-----------------------------------------------------------------------
  BundleProblem *problem;       
  
  int terminate;               //value of last call to terminator
  //see BundleTerminator for values
  
  CH_Matrix_Classes::Real weightu;                //weight u of quadratic term / inverse step size
  
  CH_Matrix_Classes::Matrix y,                     //center of stability
    lby,                        //lower bounds on variables if not free
                                //for free variable use CB_minus_infinity
    uby;                        //upper bounds on variables if not free
                                //for free variable use CB_plus_infinity
  CH_Matrix_Classes::Indexmatrix bound_index;      //indices of bounded variables
  
  CH_Matrix_Classes::Matrix newy,                  //candidate
    subg,                       //last subgradient 
    eta;                        //Lagrange multipliers for box constraints

  CH_Matrix_Classes::Real normsubg2;    //=<subg,subg>, or if (do_scaling) =<subg,inv_scale%subg>
                     //wher subg is either the initial subgradient
                     //             or the best subgradient of the aug. model

  /// points to the scaling data, (*Hp) is deleted by BundleSolver before replacement or during destruction 
  BundleScaling* Hp;

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
  int problem_changed;

  bool recompute_normsubg2;
  bool initialize_model;

  //objective values
  CH_Matrix_Classes::Real oldval;
  CH_Matrix_Classes::Real newval;
  CH_Matrix_Classes::Real linval; //value of linearized model for feasible (pojected) solution
  CH_Matrix_Classes::Real cutval; //value of cutting surface model for feasible (projected) solution
  CH_Matrix_Classes::Real modelval;
  CH_Matrix_Classes::Real augval;
  CH_Matrix_Classes::Real oldaugval;
  CH_Matrix_Classes::Real lastaugval;

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
  CH_Matrix_Classes::Integer descent_steps;        //counts number of descent_steps

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

  int solve_model(void);
    
  //--------------------------------------------------------
  //  public functions
  //---------------------------------------------------------
public:
  BundleSolver();
  ~BundleSolver();
        
  void set_defaults();  //resets all parameters to default values
                        //and calls set_defaults for terminator and bundleweight
  void clear();         //resets all variables modified by the algorithm 
                        //and calls clear() for terminator and bundleweight

  void clear_fails();   //resets all fail counts to zero

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
  void set_do_scaling(int ds);  
  void set_scaling(const CH_Matrix_Classes::Matrix& insc);
  void set_quadratic_term(BundleScaling* Sp);
  void set_clock(const CH_Tools::Clock& myclock){clockp=&myclock;}
  void set_out(std::ostream* o=0,int pril=1){out=o;print_level=pril;}
  void set_max_updates(CH_Matrix_Classes::Integer mu){max_updates=mu;}
    
  int get_terminate() const {return terminate;}
  CH_Matrix_Classes::Real get_modeleps() const {return modeleps;}
  int get_do_scaling() const;
  int get_do_heuristic_scaling() const {return do_scaling_heuristic;}
    
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
  const CH_Matrix_Classes::Matrix& get_eta() const {return eta;} 
  const CH_Matrix_Classes::Matrix& get_y() const {return y;} 
  const CH_Matrix_Classes::Matrix& get_newy() const {return newy;} 
  const CH_Matrix_Classes::Indexmatrix& get_yfixed() const {return yfixed;}
  
  const BundleProblem* get_problem() const {return problem;}
  BundleTerminator* get_terminator() const {return terminator;}
  BundleWeight*  get_bundleweight() const {return bundleweight;}
  
  int inner_loop(BundleProblem& problem,int maxsteps=0);
  
  std::ostream& print_line_summary(std::ostream& out) const;
};


class NBundleSolver;

class BundleTerminator
{
  //standard terminator, can be redefined by derivation
  //and using BundleSolver.set_terminator(SBterminator *)
  friend class BundleSolver;
  friend class NBundleSolver;

protected:
  CH_Matrix_Classes::Real termeps;
  CH_Matrix_Classes::Real subgnorm;            //if positive, do not terminate if larger
  const CH_Tools::Clock* clockp;
  //void (*save_function)(const BundleSolver* sb); 
  //if this function is set, it is called on every execution of check_termination
  CH_Tools::Microseconds timelimit;   //upper limit 
  CH_Matrix_Classes::Integer recomplimit;      //upper limit on sumrecomp,  <0 if none 
  CH_Matrix_Classes::Integer qpfailslimit;     //upper limit on sumqpfail, <0 if none
  CH_Matrix_Classes::Integer modelfailslimit;  //upper limit on summodelfail, <0 if none
  CH_Matrix_Classes::Integer augvalfailslimit; //upper limit on augvalfail, <0 if none
  CH_Matrix_Classes::Integer objevallimit;     //upper limit on evaluations, <0 if none
  CH_Matrix_Classes::Integer oraclefailslimit;     //upper limit on evaluations, <0 if none

  int terminated;  //  0    not terminated,
                   //  1    precision achieved,
                   //  2    timelimit exceeded
                   //  4    recomplimit exceeded 
                   //  8    qpfailslimit exceeded
                   // 16    modelfailslimit exceeded
                   // 32    augvalfailslimit exceeded
                   // 64    objevallimit exceeded
                   //128    oraclefailslimit exceeded
    
  int print_level;
  std::ostream* out;
  
  virtual int check_termination(BundleTerminatorData* sb)
  {
    //if (save_function) (*save_function)(sb);
    terminated=0;
    if (clockp)
      terminated|=2*(clockp->time()>=timelimit);
    if (recomplimit>=0)
      terminated|=4*(recomplimit<=sb->get_sumrecomp()); 
    if (qpfailslimit>=0)
      terminated|=8*(qpfailslimit<=sb->get_sumqpfails());
    if (modelfailslimit>=0)
      terminated|=16*(modelfailslimit<=sb->get_summodelfails());
    if (augvalfailslimit>=0)
      terminated|=32*(augvalfailslimit<=sb->get_sumaugvalfails());
    if (objevallimit>=0)
      terminated|=64*(objevallimit<=sb->get_cntobjeval());
    if (oraclefailslimit>=0)
      terminated|=128*(oraclefailslimit<=sb->get_sumoraclefails());
    if ((sb->get_suminnerit()<10)&&(sb->get_normsubg2()>0.1)) return terminated;
    if ((subgnorm>0)&&(sb->get_normsubg2()>subgnorm)) return terminated;
    CH_Matrix_Classes::Real abs_oldval=(sb->get_oldval()>0)?sb->get_oldval():-sb->get_oldval();
    terminated|=(sb->get_oldval()-sb->get_modelval()<=
		termeps*(abs_oldval+1.));
    return terminated;
  }

public:
  virtual void set_defaults()
  {
    termeps=1e-5;clockp=0;timelimit.set_infinity(true);
    recomplimit=100;qpfailslimit=100;modelfailslimit=100;
    augvalfailslimit=10;oraclefailslimit=10;objevallimit=-1; subgnorm=-1.;
  }

  virtual void clear()
  { terminated=0; }
  
  BundleTerminator()
  { /*save_function=0;*/ out=0; print_level=0; set_defaults(); clear(); }
  
  virtual ~BundleTerminator(){}

  virtual void set_termeps(CH_Matrix_Classes::Real teps){termeps=teps;}
  virtual CH_Matrix_Classes::Real get_termeps() const {return termeps;}
  virtual void set_subgnorm(CH_Matrix_Classes::Real sg){subgnorm=sg;}
  virtual CH_Matrix_Classes::Real get_subgnorm() const {return subgnorm;}
  virtual void set_timelimit(const CH_Tools::Clock* cp,CH_Tools::Microseconds tl)
  {clockp=cp;timelimit=tl;}
  virtual CH_Tools::Microseconds get_timelimit() const {return timelimit;}
  virtual void set_recomplimit(CH_Matrix_Classes::Integer rl){recomplimit=rl;}
  virtual CH_Matrix_Classes::Integer get_recomplimit() const {return recomplimit;}
  virtual void set_qpfailslimit(CH_Matrix_Classes::Integer ql){qpfailslimit=ql;}
  virtual CH_Matrix_Classes::Integer get_qpfailslimit() const {return qpfailslimit;}
  virtual void set_modelfailslimit(CH_Matrix_Classes::Integer ml){modelfailslimit=ml;}
  virtual CH_Matrix_Classes::Integer get_modelfailslimit() const {return modelfailslimit;}
  virtual void set_augvalfailslimit(CH_Matrix_Classes::Integer al){augvalfailslimit=al;}
  virtual CH_Matrix_Classes::Integer get_augvalfailslimit() const {return augvalfailslimit;}
  virtual void set_objevallimit(CH_Matrix_Classes::Integer ol){objevallimit=ol;}
  virtual CH_Matrix_Classes::Integer get_objevallimit() const {return objevallimit;}
  virtual void set_oraclefailslimit(CH_Matrix_Classes::Integer ol){oraclefailslimit=ol;}
  virtual CH_Matrix_Classes::Integer get_oraclefailslimit() const {return oraclefailslimit;}
  
  virtual int get_terminated() const {return terminated;}
  virtual void clear_terminated() {terminated=0;}
  
  //virtual void set_save_function(void (*sf)(const BundleSolver*))
  //{ save_function=sf;}    
  
  virtual void print_status(std::ostream& o) const
  {
    o<<"termination status: "<<terminated;
    if (terminated==0) { o<<" (not terminated)"<<std::endl; return;}
    if (terminated & 1){ o<<", relative precision criterion satisfied"; }
    if (terminated & 2){ o<<", timelimit exceeded"; }
    if (terminated & 4){ o<<", function reevaluation limit exceeded";} 
    if (terminated & 8){ o<<", limit of QP failures exceeded";}
    if (terminated & 16){ o<<", limit of model failures exceeded";}
    if (terminated & 32){ o<<", limit of augmented model failures exceeded";} 
    if (terminated & 64){ o<<", limit of calls to evaluation oracle exceeded";} 
    if (terminated & 128){ o<<", limit of failed oracle calls exceeded";} 
    o<<std::endl;
  }

  virtual void set_out(std::ostream* o=0,int pril=1)
  { out=o; print_level=pril; }
  //if out==0 output nothing at all
  //if out!=0 and print_level<=0 output warnings and errors only
  //if out!=0 and print_level>=1 you may output some log or debug information 
  
  virtual std::ostream& save(std::ostream& o) const
  {
    o.precision(20);
    o<<termeps<<"\n"<<timelimit<<"\n"<<recomplimit<<"\n"<<qpfailslimit<<"\n"<<modelfailslimit<<"\n"<<augvalfailslimit<<"\n"<<objevallimit<<"\n"<<oraclefailslimit<<"\n"<<subgnorm<<"\n"<<terminated<<"\n";
    return o;
  }
  
  virtual std::istream& restore(std::istream& in)
  { 
    in>>termeps;
    in>>timelimit;
    in>>recomplimit;
    in>>qpfailslimit;
    in>>modelfailslimit;
    in>>augvalfailslimit;
    in>>objevallimit;
    in>>oraclefailslimit;
    in>>subgnorm;
    in>>terminated; 
    return in;
  }
};



//-----------------------------------------------------------------------------
//----                              inline
//-----------------------------------------------------------------------------

inline void BundleSolver::set_terminator(BundleTerminator* sbt)
{delete terminator; terminator=sbt;}

  //@}

}

#endif

