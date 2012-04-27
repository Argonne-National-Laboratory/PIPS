/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/BaseSDPOracle.hxx

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



#ifndef CONICBUNDLE_BASESDPORACLE_HXX
#define CONICBUNDLE_BASESDPORACLE_HXX

//------------------------------------------------------------

#include "CBSolver.hxx"

//------------------------------------------------------------

/**@brief   oracle for SDP with matrix classes
	@author  C. Helmberg
*/	

namespace CH_Matrix_Classes {
  class Indexmatrix;
  class Matrix;
  class Symmatrix;
}


namespace ConicBundle {


  class SDPPrimal: public PrimalData
  {

  public:
    virtual ~SDPPrimal(){}

    /// returns a newly generated identical Object
    virtual PrimalData* clone_primal_data()=0;  

    /// copy the information of pd 
    virtual int assign_primal_data(const PrimalData& pd)=0;

    /// copy the information of pd 
    virtual int assign_Gram_matrix(const CH_Matrix_Classes::Matrix& P)=0;

    /// multiply this with myfactor and add itsfactor*it to this
    virtual int aggregate_primal_data(double myfactor,double itsfactor,const PrimalData& it)=0;
    /// multiply this with myfactor and add itsfactor*P*P^T to this
    virtual int aggregate_Gram_matrix(double myfactor,double itsfactor,const CH_Matrix_Classes::Matrix& P)=0;
	
  };

  /**@brief SDP oracle interface for 
     max <C,X> s.t. <A_i,X><=>b_i, <I,X><=a, X>=0 

     It is assumed that all coefficient matrices are constant after
     having once been introduced in the problem. It is, however, possible
     to remove old and add new constraints. The constraint  <I,X><=a
     is vital for the approach and may not be modified lateron. 
  */

  enum SDPtrace {
    SDPtrace_fixed,
    SDPtrace_bounded,
    SDPtrace_unbounded
  };


  class BaseSDPOracle: public FunctionObject 
  {
  public:
    /** this routine is called to set an upper bound on the trace of the 
	primal SDP matrix (=the factor for the lmax-term). 
	If the trace of the primal matrix is supposed to be exactly 
	this value, then set #SDPtrace=SDPtrace_fixed#. If no
	upper bound is known, use some guess and 
        #SDPtrace=SDPtrace_unbounded#, then the algorithm will 
	increase the value upon need.
        Changing the trace value from outside without reinitializing 
        the bundle algorithm may cause malfunction of the algorithm.
    */
    virtual void get_trace_constraint(double& trace_val,SDPtrace& tr)=0;

    /** this is called by the system for automatically updating 
	the multiplier for SDPtrace_unbounded. The suggested value
        may be increased but not decreased and must be returned */
    virtual void adjust_multiplier(double& mult)=0;

    /** write matrix inner product with the cost matrix =<C,PP^T>
	to ip_C and with the constraint matrices i to ip_opA(i) */
    virtual int gramip(const CH_Matrix_Classes::Matrix& P,double& ip_C,CH_Matrix_Classes::Matrix& ip_opA)=0;
    
    /// write inner product with constraint matrix i to vec(i)=<A_i,PP^T>
    virtual int gramip_opA(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Matrix& opA_vec)=0;

    /// get right hand side vector
    virtual const CH_Matrix_Classes::Matrix& rhs() const=0;
        
    /// compute S=P^TCP
    virtual int project_C(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Symmatrix& S)=0;

    /// compute S=P^TA_iP
    virtual int project(const int i, const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Symmatrix& S)=0;

    /** aggregate the primal matrix if desired, return 0 if not.
        Control over the newly generated SDPPrimal object is passed
        to the caller. The caller has to delete the object. */
    virtual SDPPrimal* init_primal(const CH_Matrix_Classes::Matrix& /* P */){return 0;}
    
    /** the solver calls this routine whenever new variables have been added on
        the fly in order to extend old subgradients to the new coordinates.
        If primal data was supplied for the subgradients then  
	#primal# holds a pointer to this (possibly aggregated) data,
        otherwise it is NULL. 

        In the presence of primal data, the new coordinates correspond to 
        the inner product of the primal matrix with the new constraint matrices.
        These have to be returned in the vector #ip_opA#; more precisely, 
        for i=0 to #ind.dim()# the element #ip_opA(i)#  has to hold the 
        inner product with constraint matrix #ind(i)#; 

        If primal is NULL, then the routine can only successfully
        extend the subgradients, if the constraint matrices in the
        new coordinates are zero; then the new subgradient coordinates 
        are all zero, use ip_opA.init(ind.dim(),1,0.);
        
        You may leave this function unchanged if you do not need this.

        @return 0 on success, 1 if extension is impossible */

    virtual int primalip_opA(const SDPPrimal* /* primal */,
			     const CH_Matrix_Classes::Indexmatrix& /* ind */,
			     CH_Matrix_Classes::Matrix& /* ip_opA */)
    {return 1;}
  
    
    /** compute maximum eigenvalue of the matrix C-opAt(y)
	@return 0 on success */
    virtual
    int
    evaluate
    (
     /// argument = current variables = position where to evaluate the function
     const CH_Matrix_Classes::Matrix& current_point,
     /** The columns of the matrix are orthonormal bundle vectors and 
         may help to construct good startig vectors for the eigenvalue 
         computation by iterative methods but have no other use. 
         bundlevecs may have dimension 0. */
     const CH_Matrix_Classes::Matrix& bundlevecs, 
     /** relative precision requirement for objective values
	 leading to descent steps */
     const double relprec,
     /** gives the threshold for a null step; a vector is good enough to
         yield sufficient improvement if its Ritz value exceeds 
         this Ritz_bound. Otherwise the largest of the returned
         Ritz_values must be guaranteed to lie within 
         relprec*(abs(Ritz_bound)+1.)) of the maximum eigenvalue */
     const double Ritz_bound,
     /** on input: if not zero dimensional, the matrix contains the
         Ritz_vectors returned by the last call; these may help to construct
         good starting vectors for the eigenvalue computation by iterative
         methods.
    
         on output: the matrix has to contain at least one vector who is
         either an eigenvector to the maximum eigenvalue or whose Ritz
         value exceeds the null step bound. The vectors returned must
         form an orthonormal family. */
     CH_Matrix_Classes::Matrix& Ritz_vectors,
     /** Ritz_values corresponding to the Ritz_vectors */
     CH_Matrix_Classes::Matrix&  Ritz_values
     )
     = 0;

    /** compute maximum eigenvalue of P^T(C-opAt(y))P
	@return 0 on success */
    virtual
    int
    evaluate_projection
    (
     /// argument = current variables = position where to evaluate the function
     const CH_Matrix_Classes::Matrix& current_point,
     /** orthogonal matrix defining the projection; */
     const CH_Matrix_Classes::Matrix& P, 
     /** relative precision requirement for objective values
	 leading to descent steps */
     const double relprec,
     /** gives the threshold for a null step; a vector is good enough to
         yield sufficient improvement if its Ritz value exceeds 
         this Ritz_bound. Otherwise the largest of the returned
         Ritz_values must be guaranteed to lie within 
         relprec*(abs(Ritz_bound)+1.)) of the maximum eigenvalue */
     const double Ritz_bound,
     /** on output: orthognal matrix that contains at least one column 
         vector that either is an eigenvector to the maximum eigenvalue 
         of the projected Matrix or has a Ritz value exceeding the null 
         step bound. */
     CH_Matrix_Classes::Matrix& projected_Ritz_vectors,
     /** Ritz_values corresponding to the Ritz_vectors */
     CH_Matrix_Classes::Matrix& projected_Ritz_values  
     )
     = 0;

    /** compute S=P^T(C-opAt(y))P
	@return 0 on success */
    virtual
    int
    compute_projection
    (
     /// argument = current variables = position where to evaluate the function
     const CH_Matrix_Classes::Matrix& current_point,
     /** orthogonal matrix defining the projection; */
     const CH_Matrix_Classes::Matrix& P, 
     /// resulting projected matrix
     CH_Matrix_Classes::Symmatrix& S  
     )
     = 0;

    /** compute a full eigenvector factorizaton of S=C-opAt(y) 
	@return 0 on success */
    virtual
    int
    eig
    (
     /// argument = current variables = position where to evaluate the function
     const CH_Matrix_Classes::Matrix& current_point,
     /// the columns of P give the eigenvectors 
     CH_Matrix_Classes::Matrix& P,
     /// the eigenvalues sorted in nondecreasing order
     CH_Matrix_Classes::Matrix& d  
     )
     = 0;

  };
  //@}
  

  class SDPBundleParameters: public BundleParameters
  {
  public:
    /** #n_new_subgradients# gives an upper bound on the number of new
        Ritz-vectors added in each bundle update
    
       //int n_new_subgradients;  //inherited
        
        #n_bundle_size# gives an upper bound on the number of 
	old Ritz-vectors that may be included in the new bundle
     
       //int n_bundle_size;       //inherited
       */

    /** minimum number of old Ritz-vectors that should be kept in the bundle
	no matter whether they fall below the aggregation tolerance or not. */
    int n_keep;

    /// maximum number of aggregate matrices allowed
    int n_aggregates;

    /** bundle eigenvectors (with respect to the solution matrix of the
        quadratci semidefinite subproblem) whose eigenvalues are below 
        this aggregation tolerance times the maximum eigenvalue may be
        aggregated into an aggregate matrix */
    double aggregation_tolerance;

    /** the current implementation supplies two update rules, 
	the default rule (0) and one without aggregation (1). 
        When switching to (1), all linear aggregates are lost (so the entire
        model should be reinitialized ,otherwise the bundle method may report
        augmented modle failures), the parameters n_bundle_size,n_aggregates 
        are not used, aggregation_tolerance should be set very small and specifies
        the relative size of QSDP-eigenvalues to be considered as zero, 
	n_keep indicates the number of vectors that may and should be kept on top 
	of the active and the new vectors. */
    int update_rule;
    
  };


}
#endif

