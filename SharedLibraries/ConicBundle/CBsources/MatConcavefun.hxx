/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatConcavefun.hxx

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

#ifndef __MATCONCAVEFUN_HXX__
#define __MATCONCAVEFUN_HXX__

#include "MatCBSolver.hxx"
#include "sparsmat.hxx"

namespace ConicBundle {

  /**@brief oracle interface for describing the convex
     minimization problem via their dual functions
     (e.g. if the convex problem is the Lagrangian dual
     of a concave maximization problem like maximizing a 
     piecewise linear concave functions over some 
     box constrained set)
     @author  C. Helmberg
  */


  class OnedimConcaveFunction
  {
   public:
    virtual  ~OnedimConcaveFunction(){}
    virtual int evaluate(CH_Matrix_Classes::Real add_linterm,CH_Matrix_Classes::Real relprec,
			 CH_Matrix_Classes::Real& objective_value,CH_Matrix_Classes::Real& solution) const =0;

    virtual CH_Matrix_Classes::Real primal_evaluate(CH_Matrix_Classes::Real x) const =0;
    virtual int primal_evaluate(CH_Matrix_Classes::Real x,CH_Matrix_Classes::Real &fun_value,CH_Matrix_Classes::Real& left_deriv,CH_Matrix_Classes::Real& right_deriv) const=0;

    virtual OnedimConcaveFunction* clone() const=0;
  };
  
  
  class OnedimLinearFunction: public OnedimConcaveFunction
  {
   private: 
    CH_Matrix_Classes::Real cost;   //linear cost function
    CH_Matrix_Classes::Real offset;   //constant term
    CH_Matrix_Classes::Real lb,ub;  //lower and upper bound describing the interval
    

   public:
    OnedimLinearFunction(CH_Matrix_Classes::Real co, CH_Matrix_Classes::Real of, CH_Matrix_Classes::Real l, CH_Matrix_Classes::Real u):
      cost(co),offset(of),lb(l),ub(u) {}
    ~OnedimLinearFunction(){}

    int evaluate(CH_Matrix_Classes::Real add_linterm,CH_Matrix_Classes::Real /* relprec */,
		 CH_Matrix_Classes::Real& objective_value, CH_Matrix_Classes::Real& solution) const
    {
      CH_Matrix_Classes::Real d=add_linterm+cost;
      if (d>=0)
	objective_value=offset+d*(solution=ub);
      else 
	objective_value=offset+d*(solution=lb);
      return 0;
    }		    

    CH_Matrix_Classes::Real primal_evaluate(CH_Matrix_Classes::Real x) const
    {
      if ((x<lb)||(x>ub)) return CB_minus_infinity;
      return offset+cost*x;
    }

    int primal_evaluate(CH_Matrix_Classes::Real x,CH_Matrix_Classes::Real &fun_value,CH_Matrix_Classes::Real& left_deriv,CH_Matrix_Classes::Real& right_deriv) const
    {
      if (x<lb){ 
	fun_value=CB_minus_infinity; 
	left_deriv=right_deriv=CB_plus_infinity; 
	return 0;
      }
      if (x>ub){ 
	fun_value=left_deriv=right_deriv=CB_minus_infinity; 
	return 0;
      }
      fun_value=offset+cost*x;
      left_deriv=right_deriv=cost;
      if (x<=lb) 
	left_deriv=CB_plus_infinity;
      if (x>=ub)
	right_deriv=CB_minus_infinity;
      return 0;
    }
    
    OnedimConcaveFunction* clone() const {
      return new OnedimLinearFunction(cost,offset,lb,ub);
    }
    
  };
  

  class OnedimPiecewiseLinearFunction: public OnedimConcaveFunction
  {
   private: 
    CH_Matrix_Classes::Matrix breakpoints;  // contains breakpoints in increasing order
                         // the interval between first and last entry 
                         // is the domain 
    CH_Matrix_Classes::Matrix fun_value;    // the function value at the breakpoints
    CH_Matrix_Classes::Matrix derivative;   // = (fun_value(i+1)-fun_value(i))/
                         //   (breakpoint(i+1)-breakpoint(i)), 
                         // must be decreasing when computed
    
    int compute_deriv();    
    CH_Matrix_Classes::Integer find_breakpoint(CH_Matrix_Classes::Real deriv) const;
    
   public:
    ~OnedimPiecewiseLinearFunction(){}
    OnedimPiecewiseLinearFunction(const CH_Matrix_Classes::Matrix& bp,const CH_Matrix_Classes::Matrix& fv):
      breakpoints(bp),fun_value(fv) 
    { compute_deriv();}
    
    int evaluate(CH_Matrix_Classes::Real add_linterm,CH_Matrix_Classes::Real relprec,
		 CH_Matrix_Classes::Real& objective_value, CH_Matrix_Classes::Real& solution) const;

    CH_Matrix_Classes::Real primal_evaluate(CH_Matrix_Classes::Real x) const;
    int primal_evaluate(CH_Matrix_Classes::Real x,CH_Matrix_Classes::Real &fun_value,CH_Matrix_Classes::Real& left_deriv,CH_Matrix_Classes::Real& right_deriv) const;

    OnedimConcaveFunction* clone() const 
    { return new OnedimPiecewiseLinearFunction(breakpoints,fun_value); }

    const CH_Matrix_Classes::Matrix& get_breakpoints() const {return breakpoints;}
    const CH_Matrix_Classes::Matrix& get_fun_value() const {return fun_value;}
    const CH_Matrix_Classes::Matrix& get_derivative() const {return derivative;}
    
  };


/**@brief   oracle for LP with matrix classes
	@author  C. Helmberg
*/	
  typedef std::vector<OnedimConcaveFunction*> ConcaveFunVector;

  class ConcaveCostLPFunction: public MatrixFunctionOracle 
  {
  private:

    //Vector of cost function 
    ConcaveFunVector concavefun;

    //linear coupling constraints
    CH_Matrix_Classes::Matrix b;         //column vector of dimension m of right hand sides
    CH_Matrix_Classes::Sparsemat A;      //m x n matrix of constraint coefficients

    bool generate_primals; 
    
    std::ostream* out;
    CH_Matrix_Classes::Integer print_level;
    
  public:

    ConcaveCostLPFunction(){out=0;print_level=0;generate_primals=false;}
    ~ConcaveCostLPFunction();

    void set_generate_primals(bool gp)
    { generate_primals=gp; }

    const ConcaveFunVector& get_concavefun() const {return concavefun;}
    const CH_Matrix_Classes::Sparsemat& get_A() const {return A;}
    const CH_Matrix_Classes::Matrix& get_b() const {return b;}

    const OnedimConcaveFunction* get_function(unsigned int i)
    { assert(i<concavefun.size()); return concavefun[i]; } 
    CH_Matrix_Classes::Real primal_evaluate(const CH_Matrix_Classes::Matrix& x) const;
    int primal_evaluate(const CH_Matrix_Classes::Matrix& x,CH_Matrix_Classes::Matrix &fun_value,CH_Matrix_Classes::Matrix& left_deriv,CH_Matrix_Classes::Matrix& right_deriv) const;
    
   
    //----------- Oracle Implementation of MatrixFunctionOracle ----------

 
    int evaluate(const  CH_Matrix_Classes::Matrix& current_point, CH_Matrix_Classes::Real relprec,
		 CH_Matrix_Classes::Real& objective_value,
		 CH_Matrix_Classes::Matrix& cut_values, CH_Matrix_Classes::Matrix& eps_subgradients,
		 std::vector<PrimalData*>& primal_data, PrimalExtender*& primal_extender);

    int subgradient_extension (const PrimalData* primal,
			       const CH_Matrix_Classes::Indexmatrix& variable_indices, 
			       CH_Matrix_Classes::Matrix& new_subgradient_values);

    //------------------  routines for modifying the problem ----------------

   
    /** appends new variables with their defining data.
        The control over the FunctionObjects pointed to by append_concavefun
        will be take over by this routine. They will be deleted when not
        needed any more.
    */
        
    int append_variables(const ConcaveFunVector& append_convavefun,
			 const CH_Matrix_Classes::Sparsemat& append_A);


    int append_constraints(const CH_Matrix_Classes::Sparsemat& append_A,
			   const CH_Matrix_Classes::Matrix& append_b);

    /** The new assignment of the remaining 
	variables is described by #map_to_old# so that aftwards
	position i is assigned to the variable previously indexed
        by #map_to_old(i)#. Generating copies is allowed. */        
    int reassign_variables(const CH_Matrix_Classes::Indexmatrix& map_to_old);
 
    /** delete variables. If #map_to_old# is not null then
        the correspondence of the new indices to the old indices of
        the remaining variables is described in #map_to_old# so 
        that now position i corresponds to the variable previously 
        indexed by #map_to_old(i)#. */        
    int delete_variables(const CH_Matrix_Classes::Indexmatrix& delete_indices,
			   CH_Matrix_Classes::Indexmatrix* map_to_old=0);
 
    /** The new assignment of the remaining 
	constraints is described by #map_to_old# so that aftwards
	position i is assigned the constraint previously indexed
        by #map_to_old(i)#. Generating copies is allowed. */        
    int reassign_constraints(const CH_Matrix_Classes::Indexmatrix& map_to_old);
 
    /** delete constraints. If #map_to_old# is not null then
        the correspondence of the new indices to the old indices of
        the remaining constraints is described in #map_to_old# so 
        that now position i corresponds to the constraint previously 
        indexed by #map_to_old(i)#. */        
    int delete_constraints(const CH_Matrix_Classes::Indexmatrix& delete_indices,
			   CH_Matrix_Classes::Indexmatrix* map_to_old=0);
 
   
    
    void  set_out(std::ostream* o=0,int pril=1)
    { out=o; print_level=pril; }

   };
  //@}

}

#endif
