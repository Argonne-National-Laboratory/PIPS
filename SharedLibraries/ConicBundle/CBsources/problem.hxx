/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/problem.hxx

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



#ifndef CONICBUNDLE_PROBLEM_HXX
#define CONICBUNDLE_PROBLEM_HXX

#include <typeinfo>
#include "qp_solver.hxx"
#include "bundle.hxx"
#include "CBSolver.hxx"


namespace ConicBundle {

// this describes basic routines and classes needed for the class 
// SumProblem, that allows to aplly the bundle code to the sum
// of several convex functions with a separate bundle for each function

//-----------------------------------------------------------------------------
//                               ChangeVarInfo
//-----------------------------------------------------------------------------

// serves for modifying the description of a function on the fly 

class ChangeVarInfo{
public:
  virtual ~ChangeVarInfo(){}
};


//-----------------------------------------------------------------------------
//                                AppendVars
//-----------------------------------------------------------------------------

// typical instance of ChangeVarInfo if new variables have to be 
// introduced but no particular other information has to be passed to the
// function 

// non NULL pointers are assumed to be consistent and have to be
// appended to the corresponding vectors without need for change

class AppendVars:public ChangeVarInfo
{
 public:
  CH_Matrix_Classes::Integer n_append;
  const CH_Matrix_Classes::Indexmatrix* boundsi;  //if NULL then all variables are free
  //else *boundsi gives the final new index (after append) of new bounded vars
  const CH_Matrix_Classes::Matrix* lb;       //if NULL set lower bounds to CB_minus_infinity
  const CH_Matrix_Classes::Matrix* ub;       //if NULL set upper bounds to CB_plus_infinity
  const CH_Matrix_Classes::Matrix* startval; //if NULL set to zero

  AppendVars(CH_Matrix_Classes::Integer n)
  {n_append=n;boundsi=0;lb=ub=startval=0;}

  ~AppendVars(){}
};

//-----------------------------------------------------------------------------
//                               ReassignVars
//-----------------------------------------------------------------------------

// typical instance of ChangeVarInfo if some variables are mapped to new
// positions (possibly copied several times) and unassigned ones have to be 
// deleted but no particular other information has to be passed to the function 

class ReassignVars:public ChangeVarInfo
{
 public:
  const CH_Matrix_Classes::Indexmatrix& assign_ind;

  ReassignVars(const CH_Matrix_Classes::Indexmatrix& aind): assign_ind(aind){}

  ~ReassignVars(){}
};

//-----------------------------------------------------------------------------
//                                DeleteVars
//-----------------------------------------------------------------------------

// typical instance of ChangeVarInfo if some variables have to be 
// deleted but no particular other information has to be passed to the function 

class DeleteVars:public ChangeVarInfo
{
 public:
  const CH_Matrix_Classes::Indexmatrix& del_index; 
  //eliminate all variables indexed by del_index which is sorted increasingly.
  //Remaining indices are successively shifted to the first free position

  DeleteVars(const CH_Matrix_Classes::Indexmatrix& dind): del_index(dind){}

  ~DeleteVars(){}
};


//-----------------------------------------------------------------------------
//                               ConvexProblem
//-----------------------------------------------------------------------------

// In order that a function can be used within SumProblem and CBSolver 
// it is must be derived from this class.

class ConvexProblem:public BundleProblem
{
 protected:
  QP_Solver* solverp;

  //--- output
  std::ostream *out;
  int print_level;

 public:
  virtual ~ConvexProblem(){delete solverp; solverp=0;}

  //----------- from the inherited routines of BundleProblem 
  //            only the following two are implemented here
  //            because they are the same for all derived classes so far

    virtual int eval_augmodel(const CH_Matrix_Classes::Matrix& center_y,
			      const CH_Matrix_Classes::Matrix& lby,
			      const CH_Matrix_Classes::Matrix& uby,
			      const CH_Matrix_Classes::Matrix& eta,
			      CH_Matrix_Classes::Real& augbound,
			      CH_Matrix_Classes::Real& fbound,
			      CH_Matrix_Classes::Real relprec,
			      const BundleScaling* Hp,
			      const CH_Matrix_Classes::Indexmatrix& yfixed);
  //evaluate the augmented model with respect to center of stability
  //if do_scaling!=0, inv_scale gives the inverse of a scaling of y,
  //ie, the quadratic term must be  weightu/2*[sum_i y(i)^2/inv_scale(i)]
  //if yfixed is a vector of the same dimension as the center it 
  //must give zeros for indices the may be changed and ones if they should stay fixed
  //if evaluated by an iterative method that provides upper and lower bounds
  //it should stop when  ub-lb <= min(relprec*(fbound-lb),(ub-augbound)/2.)
  //If the subproblem needs to adapt some penalty parameter and therefore
  //needs to modify the function value it can do so by changing augbound
  //and fbound (=the old function value) correspondingly (caution!)

  virtual int reeval_augmodel(const CH_Matrix_Classes::Matrix& center_y,
			      const CH_Matrix_Classes::Matrix& lby,
			      const CH_Matrix_Classes::Matrix& uby,
			      const CH_Matrix_Classes::Matrix& eta,
			      const CH_Matrix_Classes::Indexmatrix& update_index,
			      const CH_Matrix_Classes::Matrix& update_value,
			      CH_Matrix_Classes::Real& augbound,
			      CH_Matrix_Classes::Real& fbound,
			      CH_Matrix_Classes::Real relprec,
			      const BundleScaling* Hp);
  //reevaluate the augmented model for updated eta
  //the values of eta-old_eta are listed in update_index and update_value
  //if do_scaling!=0, inv_scale gives the inverse of a scaling of y,
  //ie, the quadratic term must be  weightu/2*[sum_i y(i)^2/inv_scale(i)]
  //the scaling and the weight u must be the same in both augmodel-calls!
  //if yfixed is a vector of the same dimension as the center it 
  //must give zeros for indices the may be changed and ones if they should stay fixed
  //if evaluated by an iterative method that provides upper and lower bounds
  //it should stop when  ub-lb <= min(relprec*(fbound-lb),(ub-augbound)/2.)
  //If the subproblem needs to adapt some penalty parameter and therefore
  //needs to modify the function value it can do so by changing augbound
  //and fbound (=the old function value) correspondingly (caution!)


  //----------- information needed by the class SumProblem

  virtual int intersect_box(CH_Matrix_Classes::Indexmatrix& bound_index,CH_Matrix_Classes::Matrix& lb,CH_Matrix_Classes::Matrix& ub)=0;
  //intersect input with local box constraints so that on output
  //the bounds represent the intersection
  //(it is assumed that bound_index is sorted increasingly, on input and output)

  virtual CH_Matrix_Classes::Real lb_function(const CH_Matrix_Classes::Matrix& y)=0;
  //returns a *quick* lower bound for the function value at y
  //(eg by a previous subgradient)

  virtual CH_Matrix_Classes::Real lb_model(const CH_Matrix_Classes::Matrix& y)=0;
  //returns a *quick* lower bound for the model value at y
  //(eg by one model subgradient)

  virtual int start_augmodel(QP_Block*& blockp)=0;
  //return a pointer to the variables/constraints generating the cutting model
  //returns 0 on success, 1 on failure
  
  virtual int make_aug_linmodel(CH_Matrix_Classes::Real* aug_linconst,
				CH_Matrix_Classes::Matrix* aug_subg,
				bool* conemult_increased,
				CH_Matrix_Classes::Real* function_value)=0;  
  //form the aggregate from the current solution in the QP_Block
  //this yields a minorant  linconst+<subg,.> of the objective
  //if the pointers are not nil then
  //add linconst to *aug_linconst and subg to *aug_subg
  //if the solution of QP_Block is at its bounds it may be necessary to
  //increase the bounds. In this case conemult_increased has to be set to true
  //and function_value should be reset to the new value in the center
  //(they must not be modified otherwise) 
  //returns 0 on success, 1 on failure

  
  //-------- messages needed, if online problem modifications are desired
  //         these messages are passed on recursively by sumproblem

  virtual int change_variables(ChangeVarInfo* cvp)=0;
  //delete a variable as described in ChangeVariableInfo
  //returns 0 on success, 1 on failure
  //(this is not needed in the bundle framework, routine may return 1!)

  virtual int recompute_center()=0;
  //after modifications of the problem the center information may have
  //to be recomputed partially or compeletely. Do the appropriate.
  //returns 0 on success, 1 on failure 

  virtual int adjust_multiplier()=0;
  //for conic subproblems with adjustable multipliers, reset the
  //the multiplier to twice the current trace. This message is
  //recursively passed on. 

  //-------- messages for direct get/set requests from CBSolver to problem solvers

  virtual int get_approximate_primal(PrimalData& primal) const =0;
  //if primal data is provided by the oracle then the primal corresponding
  //to the current aggregate is formed in primal (this may differ
  //if bp is a recognized derived class) 

  virtual int get_center_primal(PrimalData& primal) const =0;
  //if primal data is provided by the oracle then the primal corresponding
  //to the best eps-subgradient of the evaluation in the current center is 
  //returned (this may differ if bp is a recognized derived class) 

  virtual int get_candidate_primal(PrimalData& primal) const =0;
  //if primal data is provided by the oracle then the primal corresponding
  //to the best eps-subgradient of the evaluation in the current center is 
  //returned (this may differ if bp is a recognized derived class) 

  virtual int set_bundle_parameters(const BundleParameters& bp)=0;
  //set max_bundle_size and max_new_subgradients (this may differ
  //if bp is a recognized derived class)

  virtual int get_bundle_parameters(BundleParameters& bp) const =0;
  //returns current parameter settings (this may differ if bp is a 
  //recognized derived class)

  virtual int get_bundle_values(BundleParameters& bp) const =0;
  //returns current bundle_size and the number of new subgradients added
  //in the last update (this may differ if bp is a recognized derived class)

  virtual void clear_model()=0;
  //modifications of this specific problem were such that old subgradient
  //data and function values have to be removed completely; do so.

  virtual int clear_aggregates(){ return 1; }
  //remove all current aggregate cutting planes
  //returns 0 on success, 1 on failure 

  virtual int get_ret_code() const =0;
  //for functions given by an oracle: return the return value that 
  //was produced by the last call to the function evaluation routine

  virtual int call_primal_extender(PrimalExtender& ){return 1;}
  //if the function is the Lagrangian dual and primal_data of previous calls 
  //has now to be updated due to changes in the primal problem -- e.g., this 
  //may happen in column generation -- the problem updates all its internally 
  //stored primal_data objects by calling PrimalExtender::extend on
  //each of these.
  //returns 0 on success, 1 if not applicable to this function, 2 if it
  //would be applicable but there is no primal data.
   

  //-------- for setting output options recursively

  virtual void set_out(std::ostream* _out=0,int _print_level=1)
  { out=_out;print_level=_print_level; if (solverp) solverp->set_out(out,_print_level-1);}
   //if out==0 output nothing at all
  //if out!=0 and print_level<=0 output warnings and errors only
  //if out!=0 and print_level>=1 you may output some log or debug information 
};

}

#endif

