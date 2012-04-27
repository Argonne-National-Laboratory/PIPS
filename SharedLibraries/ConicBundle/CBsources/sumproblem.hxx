/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/sumproblem.hxx

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



#ifndef CONICBUNDLE_SUMPROBLEM_HXX
#define CONICBUNDLE_SUMPROBLEM_HXX

#include "qp_solver.hxx"
#include "qp_sumblock.hxx"
#include "problem.hxx"

namespace ConicBundle {

class AppendVarsSumProblem:public ChangeVarInfo
{
 public:
  CH_Matrix_Classes::Integer n_append;
  const CH_Matrix_Classes::Matrix* cost;     //if NULL set to zero
  const CH_Matrix_Classes::Indexmatrix* boundsi;  //if NULL then all variables are free 
  const CH_Matrix_Classes::Matrix* lb;       //if NULL set to CB_minus_infinity
  const CH_Matrix_Classes::Matrix* ub;       //if NULL set to CB_plus_infinity
  const CH_Matrix_Classes::Matrix* startval; //if NULL set to zero

  std::vector<ChangeVarInfo*> subprobp;

  AppendVarsSumProblem(CH_Matrix_Classes::Integer n)
  {n_append=n;cost=lb=ub=startval=0;boundsi=0;}

  ~AppendVarsSumProblem(){ 
    for(unsigned int i=0;i<subprobp.size();i++) delete subprobp[i];
  }
};

class ReassignVarsSumProblem:public ChangeVarInfo
{
 public:
  const CH_Matrix_Classes::Indexmatrix& assign_ind; 
  //assign to position i the variable with old index assign_ind(i), delete
  //variables that are not assigned.

  std::vector<ChangeVarInfo*> subprobp;

  ReassignVarsSumProblem(const CH_Matrix_Classes::Indexmatrix& aind): assign_ind(aind){}

  ~ReassignVarsSumProblem(){
    for(unsigned int i=0;i<subprobp.size();i++) delete subprobp[i];
  }
};

class DeleteVarsSumProblem:public ChangeVarInfo
{
 public:
  const CH_Matrix_Classes::Indexmatrix& del_index; 
  //eliminate all variables indexed by del_index
  //all later indices are shifted to the next free position

  std::vector<ChangeVarInfo*> subprobp;

  DeleteVarsSumProblem(const CH_Matrix_Classes::Indexmatrix& dind): del_index(dind){}

  ~DeleteVarsSumProblem(){
    for(unsigned int i=0;i<subprobp.size();i++) delete subprobp[i];
  }
};


class SumProblem: public ConvexProblem
{
private:
  //--- problem description
  std::vector<ConvexProblem*> subprobp;  //each describes one convex function
                                         //on destruction these are NOT deleted
  CH_Matrix_Classes::Matrix b;                  //linear cost vector for y
  CH_Matrix_Classes::Indexmatrix bounds_index;  //indices i with y(i) bounded (sorted increasingly)
  CH_Matrix_Classes::Matrix lby;                //lower bounds on y, same dim as y
  CH_Matrix_Classes::Matrix uby;                //upper bounds on y, same dim as y

  //--- flags
  int problem_modified;   //is set to true if the problem description has 
                          //been changed

  //--- data of center point
  int center_available;   //is set to true if the following four items 
                          //   are available
  CH_Matrix_Classes::Matrix y;               //current center of stability
  CH_Matrix_Classes::Real lmax;
  CH_Matrix_Classes::Matrix subg;            //minorant: subg_val+ip(subg,.-y)
  CH_Matrix_Classes::Real subg_val;        
  CH_Matrix_Classes::Real ub_fun_val;      //upper bound on the function value in y

  //--- data of last function evaluation 
  CH_Matrix_Classes::Integer ncalls;

  int cand_available;
  CH_Matrix_Classes::Matrix cand_y;
  CH_Matrix_Classes::Real cand_lmax;
  CH_Matrix_Classes::Matrix cand_subg;     
  CH_Matrix_Classes::Real cand_subg_val;
  CH_Matrix_Classes::Real cand_ub_fun_val;

  //--- data of last model evaluation
  int model_available;
  int model_subg_available;
  CH_Matrix_Classes::Real model_subg_val;
  CH_Matrix_Classes::Real model_ub_val;
  CH_Matrix_Classes::Matrix model_subg;  //this subgradient is usually not needed, lazy evaluation 

  //--- data of last augmented model evaluation
  int aug_available;
  CH_Matrix_Classes::Indexmatrix aug_ind;//starting indices for augmodel variables of each subprob
  CH_Matrix_Classes::Matrix aug_subg;    //this is the aggregate subgradient
  CH_Matrix_Classes::Real aug_linconst;  //aug_linconst+<aug_subg,.> is a minorant of the objective
   
  //--- other
  QP_SumBlock block;

public:
  void clear(void); 

  SumProblem(){solverp=0;out=0;print_level=0;clear();out=0;}
  ~SumProblem(){delete solverp;solverp=0;}

  int set_data(CH_Matrix_Classes::Integer dim,const CH_Matrix_Classes::Matrix* lby=0,const CH_Matrix_Classes::Matrix* uby=0,const CH_Matrix_Classes::Matrix* b=0);

  int add_problem(ConvexProblem* p);

  int subproblem_modified()  //call if a subproblem specification was changed 
  { 
    problem_modified=1; center_available=0; cand_available=0;
    model_available=0; model_subg_available=0; aug_available=0; 
    for(unsigned int i=0;i<subprobp.size();i++){
      if (subprobp[i]->intersect_box(bounds_index,lby,uby)){
	return 1;
      }
    }
    return 0;
  }


  CH_Matrix_Classes::Integer nsubprob() {return CH_Matrix_Classes::Integer(subprobp.size());}
  ConvexProblem* subprob(CH_Matrix_Classes::Integer i){return subprobp[i];}
  CH_Matrix_Classes::Real objval() const {return ub_fun_val;}

  //--- routines  of BundleProblem
  
  int init_center(int& init, CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix& yp, CH_Matrix_Classes::Matrix& subg,
                 CH_Matrix_Classes::Indexmatrix& bounds_index, CH_Matrix_Classes::Matrix& lby, CH_Matrix_Classes::Matrix& uby);

  int eval_function(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real nullstep_bound,CH_Matrix_Classes::Real relprec);

  int get_function_sol(CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp);

  int do_step(void);

  int eval_model(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real eval_model_bound,CH_Matrix_Classes::Real relprec);

  int get_model_sol(CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp);

  int get_augmodel_sol(CH_Matrix_Classes::Real& linconst,CH_Matrix_Classes::Matrix& subg);

  int update_model(bool descent_step);

  //--- routines of ConvexProblem
  
   int intersect_box(CH_Matrix_Classes::Indexmatrix& bound_index,CH_Matrix_Classes::Matrix& lb,CH_Matrix_Classes::Matrix& ub);

   CH_Matrix_Classes::Real lb_function(const CH_Matrix_Classes::Matrix& y);

   CH_Matrix_Classes::Real lb_model(const CH_Matrix_Classes::Matrix& y);

   int start_augmodel(QP_Block*& blockp);
 
   int get_row(CH_Matrix_Classes::Integer index_y,CH_Matrix_Classes::Matrix& row,CH_Matrix_Classes::Real& b,CH_Matrix_Classes::Integer startindex) const;
 
   int make_aug_linmodel(CH_Matrix_Classes::Real* aug_linconst,CH_Matrix_Classes::Matrix* aug_subg, 
			 bool* conemult_increased,
				CH_Matrix_Classes::Real* function_value);  

  
  //-------- messages needed, if online problem modifications are desired

  int adjust_multiplier();
  
  int change_variables(ChangeVarInfo* cvp);

  int recompute_center();
 

  //-------- messages for direct get/set requests from CBSolver to problem solvers

  int get_approximate_primal(PrimalData& /* primal */) const {return 1;}
  //this class has currently no possibility to return primal data
  //(but it could be possible to design a special type of primal that
  //distributes the collection process to its subproblems)

  int get_center_primal(PrimalData& /* primal */) const {return 1;}
  //this class has currently no possibility to return primal data
  //(but it could be possible to design a special type of primal that
  //distributes the collection process to its subproblems)

  int get_candidate_primal(PrimalData& /* primal */) const {return 1;}
  //this class has currently no possibility to return primal data
  //(but it could be possible to design a special type of primal that
  //distributes the collection process to its subproblems)

  int set_bundle_parameters(const BundleParameters& /* bp */) {return 1;}
  //this class has no bundle and therefore no such parameters

  int get_bundle_parameters(BundleParameters& /* bp */) const {return 1;}
  //this class has no bundle and therefore no such parameters

  int get_bundle_values(BundleParameters& /* bp */) const {return 1;}
  //this class has no bundle and therefore no such parameters

  void clear_model();
  //modifications of this specific problem were such that old subgradient
  //data and function values have to be removed completely; do so.

  virtual int get_ret_code() const {return 0;}
  //for functions given by an oracle: return the return value that 
  //was produced by the last call to the function evaluation routine


  //--- other

  void set_out(std::ostream* o=0,int pril=1)
  {
    out=o;print_level=pril; 
    for(unsigned int i=0;i<subprobp.size();i++){
      subprobp[i]->set_out(o,pril);
    }
    if (solverp) solverp->set_out(out,pril-1);
    block.set_out(out,pril-1);
  }

  const CH_Matrix_Classes::Matrix& get_y() const { return y; }
  const CH_Matrix_Classes::Matrix& get_cand_y() const { return cand_y; }
  const CH_Matrix_Classes::Matrix& get_lby() const { return lby; }
  const CH_Matrix_Classes::Matrix& get_uby() const { return uby; }


  CH_Matrix_Classes::Integer get_ncalls() const {return ncalls;}

  int get_center_available() const
  { return center_available; }

  int set_new_center(const CH_Matrix_Classes::Matrix* yp=0);  // yp==0 uses default starting point

  int set_lower_bound(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Real lb);

  int set_upper_bound(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Real ub);
};

}

#endif

