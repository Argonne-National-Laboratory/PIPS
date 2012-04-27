/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/socproblem.hxx

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



#ifndef CONICBUNDLE_SOCPROBLEM_HXX
#define CONICBUNDLE_SOCPROBLEM_HXX

#include "qp_solver.hxx"
#include "qp_sdpblock.hxx"
#include "problem.hxx"
#include "BaseSOCOracle.hxx"

namespace ConicBundle {

class SocProblem: public ConvexProblem
{
private:

  //--- problem description
  BaseSOCOracle* oracle;
  
  CH_Matrix_Classes::Real soc_multiplier;
  SOCtrace constant_trace; 
  
  CH_Matrix_Classes::Indexmatrix bounds_index;  //indices i with y(i) bounded (sorted increasingly)
  CH_Matrix_Classes::Matrix lby;                //lower bounds on y, same dim as y
  CH_Matrix_Classes::Matrix uby;                //upper bounds on y, same dim as y

  int problem_modified;   //is set to true if the problem description has 
                          //been changed

  //--- data of center point
  int center_available;   //is set to true if the following four items 
                          //   are available
  CH_Matrix_Classes::Matrix y;               //current center of stability
  CH_Matrix_Classes::Matrix subg;            //minorant: subg_val+ip(subg,.-y)
  CH_Matrix_Classes::Real subg_val;        
  CH_Matrix_Classes::Real ub_fun_val;        //upper bound on the function value in y
  CH_Matrix_Classes::Real lmax;
  CH_Matrix_Classes::Matrix center_vec;      

  int use_startheur;  //if 1 a heuristic is applied in starting_point to improve
                      //the starting point if the eigenvalues are well separated

  //--- data of last function evaluation 
  int cand_available;
  CH_Matrix_Classes::Matrix cand_y;
  CH_Matrix_Classes::Matrix cand_vec;
  CH_Matrix_Classes::Matrix cand_subg;     
  CH_Matrix_Classes::Real cand_subg_val;
  CH_Matrix_Classes::Real cand_ub_fun_val;
  CH_Matrix_Classes::Real cand_lmax;
  int ret_code; 

  CH_Matrix_Classes::Matrix cand_c;  //primal Lagrange cost vector at cand_y

  //--- data of last model evaluation
  int model_available;
  int model_subg_available;
  CH_Matrix_Classes::Real model_subg_val;
  CH_Matrix_Classes::Real model_ub_val;
  CH_Matrix_Classes::Matrix model_subg;  //this subgradient is usually not needed, lazy evaluation 

  int model_changed;  //if true and model_ind>=0, 
                      //then lazy subg evaluation is impossible
  CH_Matrix_Classes::Matrix  model_vec;   //for lazy evaluation of the modle subgradient

  int model_ret_code; 

  //--- data describing the model
  CH_Matrix_Classes::Integer maxkeepvecs; //max number of vectors kept in bundle from prev iteration
  CH_Matrix_Classes::Matrix bundlevecs;
  CH_Matrix_Classes::Matrix primalvec;    //primal aggregate 


  //--- data of last augmented model evaluation
  int aug_available;
  CH_Matrix_Classes::Matrix aug_subg;    //this is the aggregate subgradient
  CH_Matrix_Classes::Real aug_linconst;  //aug_linconst+<aug_subg,.> is a minorant of the objective
  CH_Matrix_Classes::Integer aug_xdim;

  //--- data for generating the bundle
  CH_Matrix_Classes::Matrix vec_store;
  CH_Matrix_Classes::Integer last_pos;
   
  //--- augmented modle solver
  QP_SDPBlock block;

  //--- temporary matrices
  mutable CH_Matrix_Classes::Matrix tmpmat;
  mutable CH_Matrix_Classes::Matrix tmpvec;
  mutable CH_Matrix_Classes::Indexmatrix tmpind;
  
  std::ostream *out;
  int print_level;
 

public:    
  SocProblem(BaseSOCOracle* oracle_function);
    
  ~SocProblem();

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
  
   
  int adjust_multiplier();
  
  int change_variables(ChangeVarInfo* cvp);
  
  int recompute_center();

  //--- other

  int get_approximate_primal(PrimalData& primal) const;

  int get_center_primal(PrimalData& primal) const;

  int get_candidate_primal(PrimalData& primal) const;

  int set_bundle_parameters(const BundleParameters& bp)
  {
    maxkeepvecs=CH_Matrix_Classes::max(2,bp.n_bundle_size);
    return 0;
  }

  int get_bundle_parameters(BundleParameters& bp) const
  {
    bp.n_bundle_size=maxkeepvecs;
    bp.n_new_subgradients=1;
    return 0;
  }

  int get_bundle_values(BundleParameters& bp) const
  {
    bp.n_bundle_size=bundlevecs.coldim();
    bp.n_new_subgradients=1;
    return 0;
  }

  void clear_model();  

  int get_ret_code() const {return ret_code;}

  //--- other

  void get_y(CH_Matrix_Classes::Matrix& iny) const {iny=y;}
  const CH_Matrix_Classes::Matrix& get_y(void) const {return y;}

  void set_out(std::ostream* o=0,int pril=1)
  { 
    out=o;print_level=pril;     
    if (solverp) solverp->set_out(out,pril-1);
    block.set_out(out,pril-1);
  }

};

}

#endif

