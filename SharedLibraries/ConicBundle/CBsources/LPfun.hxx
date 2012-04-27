/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/LPfun.hxx

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



#ifndef CONICBUNDLE_LPFUN_HXX
#define CONICBUNDLE_LPFUN_HXX

//------------------------------------------------------------

#include "CBSolver.hxx"

//------------------------------------------------------------

/**@brief   oracle for LP with matrix classes
	@author  C. Helmberg
*/	
namespace ConicBundle {

  class MatrixLPFunction;  //internal representation

  class LPMatrix 
  {
  public:
    int n_rows;
    int n_cols;
    IVector row_index;  
    IVector col_index;    
    DVector val;
  }

  class LPFunction: public FunctionObject 
  {
  private:
    MatrixLPFunction* lpf;
    
    
  public:

    LPFunction();
    ~LPFunction();

    MatrixLPFunction* get_function();

    void set_generate_primals(bool gp);
        
    int append_variables(const DVector& append_c,
			 const DVector& append_lb, 
			 const DVector& append_ub,
			 const LPMatrix& append_A);

    int append_constraints(const LPMatrix& append_A,
			   const DVector& append_b);

    /** The new assignment of the remaining 
	variables is described by #map_to_old# so that aftwards
	position i is assigned to the variable previously indexed
        by #map_to_old(i)#. Generating copies is allowed. */        
    int reassign_variables(const IVector& map_to_old);
 
    /** delete variables. If #map_to_old# is not null then
        the correspondence of the new indices to the old indices of
        the remaining variables is described in #map_to_old# so 
        that now position i corresponds to the variable previously 
        indexed by #map_to_old(i)#. */        
    int delete_variables(const IVector& delete_indices,
			   IVector* map_to_old=0);
 
    /** The new assignment of the remaining 
	constraints is described by #map_to_old# so that aftwards
	position i is assigned the constraint previously indexed
        by #map_to_old(i)#. Generating copies is allowed. */        
    int reassign_constraints(const IVector& map_to_old);
 
    /** delete constraints. If #map_to_old# is not null then
        the correspondence of the new indices to the old indices of
        the remaining constraints is described in #map_to_old# so 
        that now position i corresponds to the constraint previously 
        indexed by #map_to_old(i)#. */        
    int delete_constraints(const IVector& delete_indices,
			   IVector* map_to_old=0);
 
    void  set_out(std::ostream* o=0,int pril=1);

  };
  //@}
  
}

#endif
