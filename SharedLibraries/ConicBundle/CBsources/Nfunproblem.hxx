/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/Nfunproblem.hxx

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



#ifndef CONICBUNDLE_NMATFUNPROBLEM_HXX
#define CONICBUNDLE_NMATFUNPROBLEM_HXX

#include "MatNBSolver.hxx"

#include "Nproblem.hxx"

namespace ConicBundle {


class NFunctionProblem: public NConvexProblem
{
private:

  MatrixFunctionOracle& oracle;   //pointer to function for evaluation


  //--- data of center point
  PrimalData* center_x;          //may hold primal information
  PrimalData* new_x;          //may hold primal information
  PrimalData* aggr_x;          //may hold primal information
  
  int ret_code;           //error code returned by the last oracle evaluation call

  //--- accounting information
  CH_Tools::Clock myclock;
  CH_Tools::Microseconds evaltime;    //accumulated time spent in evaluation
  CH_Matrix_Classes::Integer nr_eval;          //counts number of "big" eigenvalue compuations
 
 
public:    
  NFunctionProblem(MatrixFunctionOracle& fo);
    
  ~NFunctionProblem();


   //--- routines  of BundleProblem
  
  virtual int eval_function(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real nullstep_bound,CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp,CH_Matrix_Classes::Real relprec);

  virtual void do_step(void);

  virtual void init_aggregate(void);

  virtual void aggregate(CH_Matrix_Classes::Real alpha);
  
  virtual void clear_aggregate(void);

  virtual int subgradient_information(const CH_Matrix_Classes::Indexmatrix& indices,CH_Matrix_Classes::Matrix& aggr_subg,CH_Matrix_Classes::Matrix& new_subg);


  //--- routines of PConvexProblem
  
   CH_Matrix_Classes::Real lb_function(const CH_Matrix_Classes::Matrix& y);

  //-------- messages needed, if online problem modifications are desired

  int adjust_multiplier(){return 0;}
  
  //-------- messages for direct get/set requests from CBSolver to problem solvers

  int get_approximate_primal(PrimalData& primal) const;

  int get_center_primal(PrimalData& primal) const 
  { if (center_x==0) return 1; return primal.assign_primal_data(*center_x); }

  int get_candidate_primal(PrimalData& primal) const 
  { if (new_x==0) return 1; return primal.assign_primal_data(*new_x); }

  int get_ret_code() const {return ret_code;}

  int call_primal_extender(PrimalExtender& prex); 

  //--- other


  void set_out(std::ostream* o=0,int pril=1)
  { 
    out=o;print_level=pril;     
  }
  
  
};

}

#endif

