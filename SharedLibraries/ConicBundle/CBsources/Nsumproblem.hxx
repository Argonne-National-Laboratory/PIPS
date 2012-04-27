/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/Nsumproblem.hxx

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



#ifndef CONICBUNDLE_NSUMPROBLEM_HXX
#define CONICBUNDLE_NSUMPROBLEM_HXX

#include "Nproblem.hxx"

namespace ConicBundle {


class NSumProblem: public NConvexProblem
{
private:
  //--- problem description
  std::vector<NConvexProblem*> subprobp;  //each describes one convex function
                                         //on destruction these are NOT deleted

public:
  void clear(void); 

  NSumProblem(){}
  ~NSumProblem(){}

  int add_problem(NConvexProblem* p);

  CH_Matrix_Classes::Integer nsubprob() {return CH_Matrix_Classes::Integer(subprobp.size());}
  NConvexProblem* subprob(CH_Matrix_Classes::Integer i){return subprobp[i];}

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

  int adjust_multiplier();
  
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
  }

};

}

#endif

