/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/CFunction.hxx

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



#ifndef CONICBUNDLE_CFUNCTION_HXX
#define CONICBUNDLE_CFUNCTION_HXX

//------------------------------------------------------------

#include "MatCBSolver.hxx"
#include "cb_cinterface.h"

//------------------------------------------------------------

/**@brief   standard function oracle with matrix classes
	@author  C. Helmberg
*/	

namespace ConicBundle {

  class CFunction: public MatrixFunctionOracle 
  {
  private:
    void* function_key;  //identifier for c-code
    cb_functionp oracle; //c-function for evaluate
    cb_subgextp subgext; //c-function for subgradient extension
    CH_Matrix_Classes::Integer primaldim;   //length of primal vectors; uses PrimalMatrix
    CH_Matrix_Classes::Integer max_new;     //maximum number of new vectors per call

  public:
    CFunction(void* fk,cb_functionp fp,cb_subgextp se=0,int prdim=0);
    ~CFunction(){};

    void set_max_new(CH_Matrix_Classes::Integer mn)
    {max_new=CH_Matrix_Classes::max(CH_Matrix_Classes::Integer(1),mn);}

    /** see MatrixFunctionOracle for explanations */
    int evaluate(
		 const  CH_Matrix_Classes::Matrix& current_point,
		 double relprec,
		 double&  objective_value,
		 CH_Matrix_Classes::Matrix&  cut_values,   
		 CH_Matrix_Classes::Matrix&  eps_subgradients,
		 std::vector<PrimalData*>& primal_data,
		 PrimalExtender*& 
     );

    int subgradient_extension(
			      const PrimalData* generating_primal,
			      const CH_Matrix_Classes::Indexmatrix& variable_indices, 
			      CH_Matrix_Classes::Matrix& new_subgradient_values
			      );

       
  };
  //@}

}
#endif

