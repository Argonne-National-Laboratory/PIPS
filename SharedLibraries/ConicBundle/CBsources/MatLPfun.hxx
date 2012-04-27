/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatLPfun.hxx

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



#ifndef CONICBUNDLE_MATLPFUN_HXX
#define CONICBUNDLE_MATLPFUN_HXX

//------------------------------------------------------------

#include "MatCBSolver.hxx"
#include "sparsmat.hxx"

//------------------------------------------------------------

/**@brief   oracle for LP with matrix classes
	@author  C. Helmberg
*/	
namespace ConicBundle {

  class MatrixLPFunction: public MatrixFunctionOracle 
  {
  private:
    CH_Matrix_Classes::Matrix c;      //column vector of dimension n of costs
    CH_Matrix_Classes::Matrix b;      //column vector of dimension m of right hand sides
    CH_Matrix_Classes::Sparsemat A;   //m x n matrix of constraint coefficients

    CH_Matrix_Classes::Matrix lb;     //column vector of dimension n of (finite!) lower bounds
    CH_Matrix_Classes::Matrix ub;     //column vector of dimension n of (finite!) upper bounds

    bool generate_primals; 
    
    CH_Matrix_Classes::Matrix tmpmat;

    std::ostream* out;
    CH_Matrix_Classes::Integer print_level;
    
  public:

    MatrixLPFunction(){out=0;print_level=0;generate_primals=false;}
    ~MatrixLPFunction(){}

    void set_generate_primals(bool gp)
    { generate_primals=gp; }

    const CH_Matrix_Classes::Matrix& get_c() const {return c;}
    const CH_Matrix_Classes::Sparsemat& get_A() const {return A;}
    const CH_Matrix_Classes::Matrix& get_b() const {return b;}
    
    const CH_Matrix_Classes::Matrix& get_lb() const {return lb;}
    const CH_Matrix_Classes::Matrix& get_ub() const {return ub;}
    
    
   
    //----------- Oracle Implementation of MatrixFunctionOracle ----------

 
    int evaluate(const  CH_Matrix_Classes::Matrix& current_point, double relprec,
		 double& objective_value,
		 CH_Matrix_Classes::Matrix& cut_values, CH_Matrix_Classes::Matrix& eps_subgradients,
		 std::vector<PrimalData*>& primal_data);

    int subgradient_extension (const PrimalData* primal,
			       const CH_Matrix_Classes::Indexmatrix& variable_indices, 
			       CH_Matrix_Classes::Matrix& new_subgradient_values);

    //------------------  routines for modifying the problem ----------------

   
    /** defines in what form primal matrices should be aggregated.
        If the argument is NULL then no primal aggregation will take place.
	The control over the generating primal is
        passed over to this. This will delete an existing generating primal 
        whenever a new generating primal is set or upon destruction. */
        
    int append_variables(const CH_Matrix_Classes::Matrix& append_c,
			 const CH_Matrix_Classes::Matrix& append_lb, const CH_Matrix_Classes::Matrix& append_ub,
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
  //@}
  
}

#endif
