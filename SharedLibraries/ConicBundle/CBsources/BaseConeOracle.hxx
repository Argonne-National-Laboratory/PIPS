/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/BaseConeOracle.hxx

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



#ifndef CONICBUNDLE_BASECONEORACLE_HXX
#define CONICBUNDLE_BASECONEORACLE_HXX

//------------------------------------------------------------

#include "MatCBSolver.hxx"

//------------------------------------------------------------

/**@brief   oracle for Cone with matrix classes
	@author  C. Helmberg
*/	

namespace CH_Matrix_Classes {
  class Indexmatrix;
  class Matrix;
}



namespace ConicBundle {


  /**@brief Cone oracle interface for 
     max c^Tx s.t. Ax<=>b, e_1^Tx<=a, x>=0 

     It is assumed that all constraints are constant after
     having once been introduced in the problem. It is, however, possible
     to remove old and add new constraints. The constraint  e_1^Tx<=a
     is vital for the approach and may not be modified lateron. 
  */

  enum Conetrace {
    Conetrace_fixed,
    Conetrace_bounded,
    Conetrace_unbounded
  };

  class BaseConeOracle: public FunctionObject 
  {
  public:
    virtual ~BaseConeOracle(){};

    virtual void get_trace_constraint(CH_Matrix_Classes::Real& trace_val,Conetrace& trace_stat)=0;

    /** this is called by the system for automatically updating 
	the multiplier for Conetrace_unbounded. The suggested value
        may be increased but not decreased and must be returned */
    virtual void adjust_multiplier(double& mult)=0;

    virtual int evaluate (const  CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real relprec,
			  CH_Matrix_Classes::Real& ubctx_ytAx, CH_Matrix_Classes::Matrix& ctx, CH_Matrix_Classes::Matrix& Ax,   
			  std::vector<PrimalData*>& primal_data 
			  )= 0;

    virtual int subgradient_extension
    (const PrimalData* /* generating_primal */,
     const CH_Matrix_Classes::Indexmatrix& /* variable_indices */, 
     CH_Matrix_Classes::Matrix& /* Aindx */)
      {return 1;}

    virtual const CH_Matrix_Classes::Matrix& rhs() const=0;

  };
  //@}
  

}
#endif

