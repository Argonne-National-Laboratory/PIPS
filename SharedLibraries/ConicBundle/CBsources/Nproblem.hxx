/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/Nproblem.hxx

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



#ifndef CONICBUNDLE_NPROBLEM_HXX
#define CONICBUNDLE_NPROBLEM_HXX

#include <typeinfo>
#include "Nbundle.hxx"
#include "CBSolver.hxx"


namespace ConicBundle {

//-----------------------------------------------------------------------------
//                               ConvexProblem
//-----------------------------------------------------------------------------

// In order that a function can be used within SumProblem and CBSolver 
// it is must be derived from this class.

class NConvexProblem:public NBundleProblem
{
 protected:
  //--- output
  std::ostream *out;
  int print_level;

 public:
  virtual ~NConvexProblem(){}
  
 //-------- messages needed, if online problem modifications are desired
  //         these messages are passed on recursively by sumproblem

  virtual int adjust_multiplier()=0;
  //for conic subproblems with adjustable multipliers, reset the
  //the multiplier to twice the current trace. This message is
  //recursively passed on. 

 //----------- information needed by the class SumProblem

  virtual CH_Matrix_Classes::Real lb_function(const CH_Matrix_Classes::Matrix& y)=0;
  //returns a *quick* lower bound for the function value at y
  //(eg by a previous subgradient)

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

  virtual int get_ret_code() const =0;
  //for functions given by an oracle: return the return value that 
  //was produced by the last call to the function evaluation routine

  virtual int call_primal_extender(PrimalExtender&){return 1;}
  //if the function is the Lagrangian dual and primal_data of previous calls 
  //has now to be updated due to changes in the primal problem -- e.g., this 
  //may happen in column generation -- the problem updates all its internally 
  //stored primal_data objects by calling PrimalExtender::extend on
  //each of these.
  //returns 0 on success, 1 if not applicable to this function, 2 if it
  //would be applicable but there is no primal data.
   

  //-------- for setting output options recursively

  virtual void set_out(std::ostream* _out=0,int _print_level=1)
  { out=_out;print_level=_print_level;}
   //if out==0 output nothing at all
  //if out!=0 and print_level<=0 output warnings and errors only
  //if out!=0 and print_level>=1 you may output some log or debug information 
};

}

#endif

