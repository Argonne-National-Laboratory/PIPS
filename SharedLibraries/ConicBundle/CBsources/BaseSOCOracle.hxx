/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/BaseSOCOracle.hxx

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



#ifndef CONICBUNDLE_BASESOCORACLE_HXX
#define CONICBUNDLE_BASESOCORACLE_HXX

//------------------------------------------------------------

#include "MatCBSolver.hxx"

//------------------------------------------------------------

/**@brief   oracle for SOC with matrix classes
	@author  C. Helmberg
*/	

namespace CH_Matrix_Classes {
  class Indexmatrix;
  class Matrix;
}

namespace ConicBundle {


  /**@brief SOC oracle interface for 
     max c^Tx s.t. Ax<=>b, e_1^Tx<=a, x>=0 

     It is assumed that all constraints are constant after
     having once been introduced in the problem. It is, however, possible
     to remove old and add new constraints. The constraint  e_1^Tx<=a
     is vital for the approach and may not be modified lateron. 
  */

  enum SOCtrace {
    SOCtrace_fixed,
    SOCtrace_bounded,
    SOCtrace_unbounded
  };

  class BaseSOCOracle: public FunctionObject 
  {
  public:
    virtual ~BaseSOCOracle(){};

    /** this routine is called to set the bound on the trace of the 
	primal SOC vector (=the factor for the lmax-term). 
	If the trace of the primal matrix is supposed to be exactly 
	this value, then set #constant_trace=true#. Changing these
	values without reinitializing the bundle algorithm may cause
	malfunction of the algorithm.
    */
    virtual void get_trace_constraint(double& trace_val,SOCtrace& constant_trace)=0;

    /** this is called by the system for automatically updating 
	the multiplier for SOCtrace_unbounded. The suggested value
        may be increased but not decreased and must be returned */
    virtual void adjust_multiplier(double& mult)=0;

    /** write inner product with the cost matrix <c,vec>
	to ipc and with the constraint row i to ipA(i) */
    virtual int ip_cA(const CH_Matrix_Classes::Matrix& vec,double& ipc,CH_Matrix_Classes::Matrix& ipA)=0;
    
    /// write inner product with constraint row i to vec(i)=<A_i,vec>
    virtual int ip_A(const CH_Matrix_Classes::Matrix& vec,CH_Matrix_Classes::Matrix& ipvecA)=0;

    /// get right hand side vector
    virtual const CH_Matrix_Classes::Matrix& rhs() const=0;
        
    /// compute c*P
    virtual int project_c(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& pr)=0;

    /// compute AiP
    virtual int project(const int i, const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& pr)=0;

    /// compute c-A^Ty
    virtual int c_At(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Matrix &cost)=0;

  };
  //@}
  

}
#endif

