/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/qp_block.hxx

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



#ifndef CONICBUNDLE_QP_BLOCK_HXX
#define CONICBUNDLE_QP_BLOCK_HXX

#include "symmat.hxx"

namespace ConicBundle {


//for given psd matrix Q>0, and vector c solve the block quadratic program
//
//   max  -1/2<Qx,x>+<c,x>+offset      min 1/2<Qx,x>+<b,y>+offset
//        A_ix_i + B_is_i = b_i            Qx+ A^Ty - z = c
//        x_i\in K_i1, s_i\in K_i2             B^Ty - t = 0
//        x=[x_i]                          z_i\in K_i1^*, t_i\in K_i2^*       
//
//  the blocks i are completetly independent 
//  the main code only works on variables x,y,z, coefficients A_i, and b_i 
//  the variables s_i and coefficients B_i are local in each block 
//  and have to be treated there (see block definition below).
//  The feasible set of each block has to be full dimensional and compact,
//  each block has to provide a strictly feasible primal-dual starting point,
//  do its line search etc.

// ****************************************************************************
//                                   QP_Block
// ****************************************************************************

// corresponds to one block in the quadratic program
// hides the structure of the constraint set and makes visible only
// the x_i and the dual variables to the constraints.
// The block holds a copy of the variables,
// knows how to find a feasible starting point, computes the step direction,
// the line search, etc.

// an implementation of QP_Block must yield a feasible primal-dual 
// predictor-corrector method with symmetric system matrix, 
// but there are no further restrictions on the choice of the step 
// direction or the step size.
// the QP_solver below will iteratively go through the steps to 
// to generate the next system matrix, compute a predictor step,
// collect suggestions for mu, compute a corrector step,
// find a common step size and move on.

class QP_Block
{
private:
  //the Block has to memorize and update its local x,y,z variables 
  //(and all additional ones that are not visible to the outside)
  //in particular it is assumed to keep a copy of the current variables
  //after the calls to the following functions:
  //  starting_x, starting_yz, restart_x, restart_yz, set_point

public:
  virtual ~QP_Block(){}

  virtual CH_Matrix_Classes::Integer xdim()=0; //dimension of externally visible primal variables
  virtual CH_Matrix_Classes::Integer ydim()=0; //dimension of externally visible dual variables
  
  virtual int set_qp_xstart(CH_Matrix_Classes::Integer x_start_index)=0;
  //the indices of the local variables correspond to the indices 
  //of the qp variables x and z starting with this index
  //returns 0 on success, 1 on failure
  
  virtual int set_qp_ystart(CH_Matrix_Classes::Integer y_start_index)=0;
  //the indices of the local variables correspond to the indices 
  //of the qp variables y starting with this index
  //returns 0 on success, 1 on failure

  virtual int starting_x(CH_Matrix_Classes::Matrix& qp_x)=0;
  //generate a strictly feasible primal starting point
  //store it in the qpx_range of x
  //returns 0 on success, 1 on failure

  virtual int starting_yz(CH_Matrix_Classes::Matrix& qp_y,CH_Matrix_Classes::Matrix& qp_z,
		   const CH_Matrix_Classes::Matrix& qp_Qx, const CH_Matrix_Classes::Matrix& qp_c)=0;
  //generate a strictly feasible dual starting point
  //store it in the qpy_range of y and in the qpx_range of z
  //x is fixed already by a previous call to starting_x and Qx=Q*x
  //returns 0 on success, 1 on failure

  virtual int get_Ab(CH_Matrix_Classes::Matrix& qp_A,CH_Matrix_Classes::Matrix &qp_b) const=0;
  //store the local coefficients of matrices A and b in the positions
  //corresponding to qpy_range (rows) and qpx_range (columns) 
  //returns 0 on success, 1 on failure

  virtual int restart_x(CH_Matrix_Classes::Matrix& qp_x,const CH_Matrix_Classes::Matrix& qp_c,const CH_Matrix_Classes::Matrix& qp_dc)=0;
  //it is assumed that the problem was solved already once and is now
  //resolved for a new linear cost term qp_c that resulted from the old
  //one by adding qp_dc.
  //on input qp_x holds the old optimal solution and on output
  //the coorespoind qpx_range should be replaced by a reasonable 
  //strictly feasible solution for x suitable for restarting
  //(see also restart_yz)
  //returns 0 on success, 1 on failure

  virtual int restart_yz(CH_Matrix_Classes::Matrix& qp_y,CH_Matrix_Classes::Matrix& qp_z,
		  const CH_Matrix_Classes::Matrix& qp_Qx,const CH_Matrix_Classes::Matrix& qp_c,const CH_Matrix_Classes::Matrix& qp_dc)=0;
  //this is called after restart_x (see there)
  //on input qp_y and qp_z hold the old optimal solution and on output
  //the coorespoind qpy/qpx_range should be replaced by a reasonable 
  //strictly feasible solution for y/z suitable for restarting
  //returns 0 on success, 1 on failure
    
  virtual int add_xinv_kron_z(CH_Matrix_Classes::Symmatrix& barQ)=0;
  //add the system term corresponding to (xinv kron z)
  //(that arises from solving the linearized perturbed complementarity system
  // x*z =0 or =mu*I for dx in the preferred search direction)
  //to the diagonal block corresponding to qpx_range x qpx_range

  virtual int add_local_sys(CH_Matrix_Classes::Symmatrix& sysdy,CH_Matrix_Classes::Matrix& rhs)=0;
  //on input: sysdy= A*barQ^{-1}*A^T    (barQ as returned in add_xinv_kron_z)
  //          rhs= A*barQ^{-1}*(c-Q*x-A^T*y)-(b-A*x)
  //if the block uses additional internal variables 
  //(like an additional term + B*s with s>=0 in the primal feasibility constr)
  //then the corresponding block terms have now to be added to sysdy and rhs,
  //eg,
  //   sysdy +=  B*(t^{-1} kron s)*B^T     (if t is the dual variable to s)
  //   rhs   +=  B*s - B*(t^{-1} kron s)*B^T*y

  virtual int suggest_mu(CH_Matrix_Classes::Real& ip_xz,CH_Matrix_Classes::Integer& mu_dim,CH_Matrix_Classes::Real& sigma,
	      const CH_Matrix_Classes::Matrix& qp_dx,const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz)=0;
  //dx, dy, dz is the predictor direction. Based on the predictor 
  //(x+dx,y+dy,z+dz) suggest a value for mu by specifying the
  //inner product of the dual cone variables ip_xz=ip(x,z)+ip(s,t),
  //the dimension of the conic variable space mu_dim= cone_x.dim+cone_s.dim
  //a value for the factor on mu to obtain the new target

  virtual int get_corr(CH_Matrix_Classes::Matrix& xcorr,CH_Matrix_Classes::Matrix& rhs,CH_Matrix_Classes::Real mu,const CH_Matrix_Classes::Matrix& qp_dx,
		       const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz)=0;
  //on input (w.r.t. corresponding positions)
  //    xcorr = 0
  //    rhs as on output of add_local_sys
  //on output the corresponding positions of xcorr should hold the corrector
  //term of the search direction, eg,  xcorr = mu*x^{-1} - x^{-1}*dx*dz,
  //and if the block holds additional local variables as in add_local_sys then
  //   rhs += B*(mu * t^{-1}- t^{-1}*dt*ds) 

  virtual int line_search(CH_Matrix_Classes::Real& alpha,const CH_Matrix_Classes::Matrix& qp_dx,
                          const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz)=0;
  //dx,dy,dz gives the final step direction, alpha is on input
  //an upper bound on the step size. 
  //On output alpha has to be no larger than on input and
  //has to guarantee strict feasibility of the primal/dual step on
  //the local variables.
  //(if the block has additional internal variables it has to compute the step
  // direction for these now and to choose alpha so that
  // strict feasibility is guaranteed for the internal variables as well)

  virtual int set_point(const CH_Matrix_Classes::Matrix& qp_x,const CH_Matrix_Classes::Matrix& qp_y,const CH_Matrix_Classes::Matrix& qp_z,
			CH_Matrix_Classes::Real alpha)=0;
  //x,y,z is the new point and has to be stored
  //alpha is the step size used in the step, it is passed so that 
  //the block can take the same step for internal variables if needed.

  virtual void set_out(std::ostream* out=0, int print_level=1)=0;


  //---------------- for debugging purposes

  virtual int add_Bs(CH_Matrix_Classes::Matrix& qp_vec) const=0;
  //add the local product of matrices B and s in the positions
  //corresponding to qpy_range (rows) 
  //returns 0 on success, 1 on failure



};

}

#endif

