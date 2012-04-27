/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatLPBCfun.hxx

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



#ifndef CONICBUNDLE_MATLPBCFUN_HXX
#define CONICBUNDLE_MATLPBCFUN_HXX

//------------------------------------------------------------

#include <vector>
#include "MatCBSolver.hxx"
#include "sparsmat.hxx"
#include "symmat.hxx"

//------------------------------------------------------------

/**@brief   oracle for LPs over bounded cones with matrix classes
	@author  C. Helmberg
*/	
namespace ConicBundle {

  enum ConeType {CT_sdp,CT_soc,CT_lp};

  class BoundedConeData
  {
  public:
    ConeType conetype;  //default= CT_lp
    CH_Matrix_Classes::Integer dim;        //for SDP give matrix order, for SOC give n, for LP give 1

    //upper and lower bounds on trace 
    CH_Matrix_Classes::Real lb;            //lower bound on trace, default=0. 
                        //>=0. for SDP and SOC, may be negative for LP
    CH_Matrix_Classes::Real ub;            //upper bound on trace, default=1.
                        //>=0. for SDP and SOC, max be negative for LP
    //additionally, bounds must be in [CB_minus_infinity,CB_plus_infinity]

    BoundedConeData():conetype(CT_lp),dim(CH_Matrix_Classes::Integer(1)),lb(0.),ub(1.){}
    BoundedConeData(ConeType ict,CH_Matrix_Classes::Integer idim=1,CH_Matrix_Classes::Real ilb=0.,CH_Matrix_Classes::Real iub=1.)
    { conetype=ict; dim=idim; lb=ilb; ub=iub; }
      
   ~BoundedConeData(){}

    std::ostream& print_data(std::ostream& out);
  };

  typedef std::vector<BoundedConeData> BCDvector;

  class MatrixLPBCFunction : public MatrixFunctionOracle
  {
  private:
    BCDvector conedata;
    CH_Matrix_Classes::Matrix c;      //column vector of dimension n of costs
    CH_Matrix_Classes::Matrix b;      //column vector of dimension m of right hand sides
    CH_Matrix_Classes::Sparsemat A;   //sparse m x n matrix of constraint coefficients
                   //m... number of constraints, n dimension of variables:
                   //each CT_lp takes 1 column
                   //each CT_soc of dim k needs k columns
                   //each CT_sdp of dim k needs k*(k+1)/2 columns
                   //use svec: a11 sqrt(2)*a12 ... sqrt(2)*a1k a22 sqrt(2)*a23 .. akk
                   //they are concatenated in the sequence of conedata

    bool generate_primals;  
    //if true, this generates dense vectors with rows corresponding to cols of A

    CH_Matrix_Classes::Matrix tmpvec;
    CH_Matrix_Classes::Matrix tmpmat;
    CH_Matrix_Classes::Symmatrix tmpsym;

    std::ostream* out;
    CH_Matrix_Classes::Integer print_level;
    
  public:

    MatrixLPBCFunction(){generate_primals=false;out=0;print_level=0;}
    ~MatrixLPBCFunction(){}

    const BoundedConeData& get_conedata(CH_Matrix_Classes::Integer i) const {return conedata[i];}

    const CH_Matrix_Classes::Matrix& get_c() const {return c;}
    const CH_Matrix_Classes::Sparsemat& get_A() const {return A;}
    const CH_Matrix_Classes::Matrix& get_b() const {return b;}

    void set_generate_primals(bool gp)
    { generate_primals=gp; }

    //----------- Oracle Implementation of MatrixFunctionOracle ----------
 
    int evaluate(const  CH_Matrix_Classes::Matrix& current_point, double relprec,
		 double& objective_value,
		 CH_Matrix_Classes::Matrix& cut_values, CH_Matrix_Classes::Matrix& eps_subgradients,
		 std::vector<PrimalData*>& primal_data);

    int subgradient_extension (const PrimalData* primal,
			       const CH_Matrix_Classes::Indexmatrix& variable_indices, 
			       CH_Matrix_Classes::Matrix& new_subgradient_values);

    //------------------  routines for modifying the problem ----------------

    int append_variables(const BCDvector& append_conedata, 
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

