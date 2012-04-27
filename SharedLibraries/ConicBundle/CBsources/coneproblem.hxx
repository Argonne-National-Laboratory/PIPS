/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/coneproblem.hxx

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



#ifndef CONICBUNDLE_CONEPROBLEM_HXX
#define CONICBUNDLE_CONEPROBLEM_HXX

#include "BaseConeOracle.hxx"

#include "qp_solver.hxx"
#include "qp_sdpblock.hxx"
#include "problem.hxx"

namespace ConicBundle {

class ConeProblem: public ConvexProblem
{
private:

  //--- problem description
  CH_Matrix_Classes::Integer dim;               //dimension of design space
  BaseConeOracle* oracle;   //pointer to function for evaluation

  CH_Matrix_Classes::Real cone_multiplier;
  Conetrace trace_stat; 

  CH_Matrix_Classes::Indexmatrix bounds_index;  //indices i with y(i) bounded (sorted increasingly)
  CH_Matrix_Classes::Matrix lby;                //lower bounds on y, same dim as y
  CH_Matrix_Classes::Matrix uby;                //upper bounds on y, same dim as y

  int problem_modified;   //is set to true if the problem description has 
                          //been changed

  //--- data of center point
  int center_available;   //is set to true if the following four items 
                          //   are available
  CH_Matrix_Classes::Matrix y;               // current center of stability
  CH_Matrix_Classes::Real ublmax;            // upper bound on ctx-ytAx
  CH_Matrix_Classes::Real lmax;              // max of ctx-ytAx
  CH_Matrix_Classes::Matrix Ax;
  CH_Matrix_Classes::Matrix subg;            // minorant: subg_val+ip(subg,.-y)
  CH_Matrix_Classes::Real subg_val;        
  CH_Matrix_Classes::Real ub_fun_val;        //upper bound on the function value in y

  PrimalData* center_x;          //may hold primal information

  //--- data of last function evaluation
  int cand_available;
  CH_Matrix_Classes::Matrix cand_y;
  CH_Matrix_Classes::Real cand_ublmax;
  CH_Matrix_Classes::Real cand_lmax;
  std::vector<PrimalData*> cand_x;
  CH_Matrix_Classes::Matrix cand_subg;   
  CH_Matrix_Classes::Matrix cand_ctx;
  CH_Matrix_Classes::Matrix cand_Ax;
  CH_Matrix_Classes::Real cand_subg_val;
  CH_Matrix_Classes::Real cand_ub_fun_val;
  int ret_code; 
  CH_Matrix_Classes::Integer cand_subg_maxind;

  //--- data of last model evaluation
  int model_available;
  int model_subg_available;
  CH_Matrix_Classes::Real model_subg_val;
  CH_Matrix_Classes::Real model_ub_val;
  CH_Matrix_Classes::Matrix model_subg;  //this subgradient is usually not needed, lazy evaluation 

  int model_changed;  //if true and model_ind>=0, 
  CH_Matrix_Classes::Integer model_ind;   //index of the maximizing aggregate subgradient 

  //--- data describing the model
  CH_Matrix_Classes::Integer maxkeepvecs; //max number of vectors kept in the bundle >=2;
  CH_Matrix_Classes::Integer max_new_subg;//maximum number of new subgradients added per evaluation
  CH_Matrix_Classes::Integer n_new_subg;  //actual number of new subgaridents added in last update
  
  CH_Matrix_Classes::Real aggregtol;      //aggregate primalvector if primalvalue smaller than

  CH_Matrix_Classes::Matrix bundlecosts;  //cost corresponding to the subgradients
  CH_Matrix_Classes::Matrix bundlevecs;   //aggregate sugradients of form Ax in the bundle
  CH_Matrix_Classes::Matrix bundlecoeff;  //coefficients determined in last aug_model
  
  CH_Matrix_Classes::Matrix hashsum;      //for recognizing identical sugradients
                       //the sum of the vector elements is stored as fingerprint
 
  std::vector<PrimalData*> bundlex;  //primal vectors corresponding to the bundle

  //--- data of last augmented model evaluation
  int aug_available;
  CH_Matrix_Classes::Integer aug_xdim;
  CH_Matrix_Classes::Matrix aug_subg;    //this is the aggregate subgradient
  CH_Matrix_Classes::Real aug_linconst;  //aug_linconst+<aug_subg,.> is a minorant of the objective

  //--- augmented modle solver
  QP_SDPBlock block;

  //--- temporary matrices
  mutable CH_Matrix_Classes::Matrix tmpvec;
  mutable CH_Matrix_Classes::Indexmatrix tmpind;
 
  //--- accounting information
  CH_Tools::Clock myclock;
  CH_Tools::Microseconds evaltime;    //accumulated time spent in evaluation
  CH_Matrix_Classes::Integer nr_eval;          //counts number of "big" eigenvalue compuations
 
  std::ostream *out;
  int print_level; 
 
public:    
  ConeProblem(BaseConeOracle* fo);
    
  ~ConeProblem();


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

  int get_approximate_primal(PrimalData& primal) const;

  int get_center_primal(PrimalData& primal) const 
  { if ((!center_available)||(center_x==0)) return 1; 
  return primal.assign_primal_data(*center_x); }

  int get_candidate_primal(PrimalData& primal) const 
  { if ((!cand_available)||(cand_x.size()==0)) return 1; 
  return primal.assign_primal_data(*cand_x[cand_subg_maxind]); }

  int set_bundle_parameters(const BundleParameters& bp)
  {
    maxkeepvecs=CH_Matrix_Classes::max(2,bp.n_bundle_size);
    max_new_subg=CH_Matrix_Classes::max(1,bp.n_new_subgradients);
    return 0;
  }

  int get_bundle_parameters(BundleParameters& bp) const
  {
    bp.n_bundle_size=maxkeepvecs;
    bp.n_new_subgradients=max_new_subg;
    return 0;
  }

  int get_bundle_values(BundleParameters& bp) const
  {
    bp.n_bundle_size=bundlecosts.dim();
    bp.n_new_subgradients=n_new_subg;
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

