/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/Nsumproblem.cxx

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
#include "Nsumproblem.hxx"

using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              clear()
// *****************************************************************************

void NSumProblem::clear(void)
{
  subprobp.clear();  
}



// *****************************************************************************
//                            add_problem
// *****************************************************************************

int NSumProblem::add_problem(NConvexProblem* p)
{
  subprobp.push_back(p);
  p->set_out(out,print_level);

  return 0;
}

// *****************************************************************************
//                            eval_function
// *****************************************************************************

  //evaluates the objective function in $y$
  //if evaluated by an iterative method that provides upper and lower bounds,
  //it may stop when the lower bound (lb) is above the nullstep_bound
  //evalutation may also stop, if (ub-lb)<relprec*(|ub|+1) is satisfied
  //it returns a lower and an upper bound on the objective value 
  //if subgp is not 0 then an eps-subgradient at y is returned,
  //the hyperplane corresponding to the subgradient has value lb in y
  //returns:  0 ... if all is ok, use get_function_sol or do_step
  //          1 ... if solution could not be computed to desired precision


int NSumProblem::eval_function(const Matrix& iny,
			       Real nullstep_bound,Real& lb, Real& ub, Matrix* subgp,Real relprec)
{
  if (subprobp.size()==1){    //just one subproblem
    return subprobp[0]->eval_function(iny,nullstep_bound,lb,ub,subgp,relprec);
  }

  //----  nonnegative sum of several subproblems

  //compute a quick lower bound on each subfunction
  Matrix fun_lb(Integer(subprobp.size()),1,0.);   
  for (unsigned int i=0;i<subprobp.size();i++){
    fun_lb(i)=subprobp[i]->lb_function(iny);
  }
  Real lower_bound=sum(fun_lb);
  lb=0.;
  ub=0.;
  
  //compute values for subproblems
  int retval=0;    
  int subg_valid=1;
  Matrix tmpvec;
  Matrix* tmpvecp=&tmpvec;
  if (subgp==0)
    tmpvecp=0;
  
  {for (unsigned int i=0;i<subprobp.size();i++){
    Real corr=lower_bound-fun_lb(i);
    Real lsg_val;
    Real lub_val;
    retval|=subprobp[i]->eval_function(iny,
				       (nullstep_bound-corr),lsg_val,lub_val,tmpvecp,
					relprec);
    //collect solution information
    lb+=lsg_val;
    ub+=lub_val;
    if ((subg_valid)&&(subgp!=0)){
      if (iny.dim()!=tmpvecp->dim()){
	subgp->init(0,0,0.);
	subg_valid=0;
      }
      else if (i==0){
	subgp->init(*tmpvecp);
      }
      else {
	(*subgp)+=*tmpvecp;
      }
    }
	
  }}

  return retval;
}

// *****************************************************************************
//                              do_step
// *****************************************************************************

void NSumProblem::do_step(void)
{
  for (unsigned int i=0;i<subprobp.size();i++){
    subprobp[i]->do_step();
  }
}

// *****************************************************************************
//                              init_aggregate
// *****************************************************************************

void NSumProblem::init_aggregate(void)
{
  for (unsigned int i=0;i<subprobp.size();i++){
    subprobp[i]->init_aggregate();
  }

}

// *****************************************************************************
//                              aggregate
// *****************************************************************************

void NSumProblem::aggregate(Real alpha)
{
  for (unsigned int i=0;i<subprobp.size();i++){
    subprobp[i]->aggregate(alpha);
  }

}


// *****************************************************************************
//                              clear_aggregate
// *****************************************************************************

void NSumProblem::clear_aggregate(void)
{
  for (unsigned int i=0;i<subprobp.size();i++){
    subprobp[i]->clear_aggregate();
  }
}

// *****************************************************************************
//                        subgradient_information
// *****************************************************************************

int NSumProblem::subgradient_information(const Indexmatrix& indices,CH_Matrix_Classes::Matrix& aggr_subg,CH_Matrix_Classes::Matrix& new_subg)
{
  aggr_subg.init(indices.dim(),1,0.);
  new_subg.init(indices.dim(),1,0.);
  for (unsigned int i=0;i<subprobp.size();i++){
    Matrix tmp_aggr;
    Matrix tmp_subg;
    if (subprobp[i]->subgradient_information(indices,tmp_aggr,tmp_subg)){
      return i+1;
    }
    aggr_subg+=tmp_aggr;
    new_subg+=tmp_subg;
  }
  return 0;
}


// *****************************************************************************
//                                 lb_function
// *****************************************************************************

//returns a *quick* lower bound for the function value at y
//(eg by a previous subgradient)

Real NSumProblem::lb_function(const Matrix& iny)
{
  Real lb=0;
  for (unsigned int i=0;i<subprobp.size();i++){
    lb+=subprobp[i]->lb_function(iny);
  }
  return lb;
}

// *****************************************************************************
//                                adjust_multiplier
// *****************************************************************************

int NSumProblem::adjust_multiplier()  
{
  int retval=0;
  for (unsigned int i=0;i<subprobp.size();i++){
    retval|=subprobp[i]->adjust_multiplier();
  }
  return retval;
}
 

}

