/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatConefun.hxx

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



#ifndef CONICBUNDLE_MATCONEFUN_HXX
#define CONICBUNDLE_MATCONEFUN_HXX

//------------------------------------------------------------

#include <map>
#include "BaseConeOracle.hxx"
#include "sparsmat.hxx"

//------------------------------------------------------------

/**@brief   oracle for Cone with matrix classes
	@author  C. Helmberg
*/	
namespace ConicBundle {

  class MatrixConeFunction : public BaseConeOracle
  {
  private:
    CH_Matrix_Classes::Real trace_value;
    Conetrace trace_stat;

    CH_Matrix_Classes::Indexmatrix socdim; //each entry gives dimension of second order cone
    CH_Matrix_Classes::Matrix c;      //column vector of dimension n of costs
    CH_Matrix_Classes::Matrix b;      //column vector of dimension m of right hand sides
    CH_Matrix_Classes::Sparsemat A;   //m x n matrix of constraint coefficients

    CH_Matrix_Classes::Integer maxsol; //the maximum number of solution vectors selected in eval

    bool generate_primals; 

    CH_Matrix_Classes::Indexmatrix socstart;
    CH_Matrix_Classes::Matrix tmpvec;

    std::ostream* out;
    CH_Matrix_Classes::Integer print_level;
    
  public:

    MatrixConeFunction()
    {generate_primals=false;out=0;print_level=0;trace_value=1.;trace_stat=Conetrace_fixed;maxsol=5;}
    ~MatrixConeFunction(){}

    const CH_Matrix_Classes::Indexmatrix& get_socdim() const {return socdim;}
    const CH_Matrix_Classes::Indexmatrix& get_socstart() const {return socstart;}
    const CH_Matrix_Classes::Matrix& get_c() const {return c;}
    const CH_Matrix_Classes::Sparsemat& get_A() const {return A;}
    const CH_Matrix_Classes::Matrix& get_b() const {return b;}

    void set_generate_primals(bool gp)
    { generate_primals=gp; }

    void set_maxsol(CH_Matrix_Classes::Integer ms)
    { maxsol=CH_Matrix_Classes::max(ms,1); }
   
    //----------- Oracle Implementation of MatrixFunctionOracle ----------

    /** this routine is called to set the bound on the trace of the 
	primal Cone vector (=the factor for the lmax-term). 
	If the trace of the primal matrix is supposed to be exactly 
	this value, then set #constant_trace=Conetrace_fixed#. Changing these
	values without reinitializing the bundle algorithm may cause
	malfunction of the algorithm.
    */
    void get_trace_constraint(CH_Matrix_Classes::Real& trace_val,Conetrace& trace_st)
      { 
	trace_val=trace_value; 
	trace_st=trace_stat;
      }

    void adjust_multiplier(CH_Matrix_Classes::Real& mult)
    { 
      if ((trace_stat==Conetrace_unbounded)&&(mult>0.))
	trace_value=mult;
      else mult=trace_value;
    }
    
    /// get right hand side vector
    const CH_Matrix_Classes::Matrix& rhs() const
      { 
	return b;
      } 
        

    int evaluate (const  CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real relprec,
		  CH_Matrix_Classes::Real& ubctx_ytAx, CH_Matrix_Classes::Matrix& ctx, CH_Matrix_Classes::Matrix& Ax,   
		  std::vector<PrimalData*>& primal_data 
		  );

    int subgradient_extension (const PrimalData* generating_primal,
			       const CH_Matrix_Classes::Indexmatrix& variable_indices, 
			       CH_Matrix_Classes::Matrix& Aindx);

    //------------------  routines for modifying the problem ----------------


    /// set the basic constraint <I,X> <=/= trace; #trace# must be nonnegative
    int set_trace(const CH_Matrix_Classes::Real trace, Conetrace constant_tr)
      { 
	if (trace<0.) return 1; 
	trace_value=trace; trace_stat=constant_tr;
	return 0;
      }
        
    int append_variables(const CH_Matrix_Classes::Indexmatrix& append_socdim, 
			 const CH_Matrix_Classes::Matrix& append_c,
			 const CH_Matrix_Classes::Sparsemat& append_A);


    int append_constraints(const CH_Matrix_Classes::Sparsemat& append_A,
			   const CH_Matrix_Classes::Matrix& append_b);

    /** The new assignment of the remaining 
	variables is described by #map_to_old# so that aftwards
	position i is assigned to the variable previously indexed
        by #map_to_old(i)#. Generating copies is allowed. */        
    int reassign_variables(const CH_Matrix_Classes::Indexmatrix& map_to_old_socvars);
 
    /** delete variables. If #map_to_old# is not null then
        the correspondence of the new indices to the old indices of
        the remaining variables is described in #map_to_old# so 
        that now position i corresponds to the variable previously 
        indexed by #map_to_old(i)#. */        
    int delete_variables(const CH_Matrix_Classes::Indexmatrix& delete_socvar_indices,
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

    std::ostream& print_problem_data(std::ostream& out);

  };
  
}

#endif

