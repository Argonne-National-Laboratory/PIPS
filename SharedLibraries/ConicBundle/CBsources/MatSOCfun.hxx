/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatSOCfun.hxx

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



#ifndef CONICBUNDLE_MATSOCFUN_HXX
#define CONICBUNDLE_MATSOCFUN_HXX

//------------------------------------------------------------

#include <map>
#include "BaseSOCOracle.hxx"
#include "sparsmat.hxx"

//------------------------------------------------------------

/**@brief   oracle for SOC with matrix classes
	@author  C. Helmberg
*/	
namespace ConicBundle {

  class MatrixSOCFunction : public BaseSOCOracle
  {
  private:
    CH_Matrix_Classes::Real trace_value;
    SOCtrace constant_trace;

    CH_Matrix_Classes::Matrix c;      //column vector of dimension n of costs
    CH_Matrix_Classes::Matrix b;      //column vector of dimension m of right hand sides
    CH_Matrix_Classes::Sparsemat A;   //m x n matrix of constraint coefficients

    std::ostream* out;
    CH_Matrix_Classes::Integer print_level;
    
  public:

    MatrixSOCFunction(){out=0;print_level=0;trace_value=1.;constant_trace=SOCtrace_fixed;}
    ~MatrixSOCFunction(){}

    const CH_Matrix_Classes::Matrix& get_c() const {return c;}
    const CH_Matrix_Classes::Sparsemat& get_A() const {return A;}
    const CH_Matrix_Classes::Matrix& get_b() const {return b;}
    
   
    //----------- Oracle Implementation of MatrixFunctionOracle ----------

    /** this routine is called to set the bound on the trace of the 
	primal SOC vector (=the factor for the lmax-term). 
	If the trace of the primal matrix is supposed to be exactly 
	this value, then set #constant_trace=SOCtrace_fixed#. Changing these
	values without reinitializing the bundle algorithm may cause
	malfunction of the algorithm.
    */
    void get_trace_constraint(CH_Matrix_Classes::Real& trace_val,SOCtrace& trace_stat)
      { 
	trace_val=trace_value; 
	trace_stat=constant_trace;
      }

    void adjust_multiplier(CH_Matrix_Classes::Real& mult)
    { 
      if ((constant_trace==SOCtrace_unbounded)&&(mult>0.))
	trace_value=mult;
      else mult=trace_value;
    }
    
    /// write inner product with constraint row i to vec(i)=<A_i,vec>
    virtual int ip_A(const CH_Matrix_Classes::Matrix& vec,CH_Matrix_Classes::Matrix& ipvecA)
      { 
	CH_Matrix_Classes::genmult(A,vec,ipvecA); 
	return 0;
      }

    /** write inner product with the cost matrix <c,vec>
	to ipc and with the constraint row i to ipA(i) */
    virtual int ip_cA(const CH_Matrix_Classes::Matrix& vec,double& ipc,CH_Matrix_Classes::Matrix& ipA)
      { 
	ipc=CH_Matrix_Classes::ip(c,vec); 
	return ip_A(vec,ipA); 
      }

    /// get right hand side vector
    virtual const CH_Matrix_Classes::Matrix& rhs() const
      { 
	return b;
      } 
        
    /// compute P^Tc
    virtual int project_c(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& pr)
      { 
	CH_Matrix_Classes::genmult(P,c,pr,1.,0.,1); 
	return 0;  
      }

    /// compute AiP
    virtual int project(const int i, const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& pr)
      { 
	CH_Matrix_Classes::genmult(A.row(i),P,pr);
	pr.transpose();
	return 0; 
      }

    /// compute c-A^Ty
    virtual int c_At(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Matrix &cost)
      { 
	cost=c; 
	CH_Matrix_Classes::genmult(A,y,cost,-1.,1.,1); 
	return 0; 
      }


    //------------------  routines for modifying the problem ----------------


    /// set the basic constraint <I,X> <=/= trace; #trace# must be nonnegative
    int set_trace(const CH_Matrix_Classes::Real trace, SOCtrace constant_tr)
      { 
	if (trace<0.) return 1; 
	trace_value=trace; constant_trace=constant_tr;
	return 0;
      }
        
    int append_variables(const CH_Matrix_Classes::Matrix& append_c,
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

    std::ostream& print_problem_data(std::ostream& out);

  };
  
}

#endif

