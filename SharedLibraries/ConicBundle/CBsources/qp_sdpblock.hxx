/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/qp_sdpblock.hxx

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



#ifndef CONICBUNDLE_QP_SDPBLOCK_HXX
#define CONICBUNDLE_QP_SDPBLOCK_HXX

#include <vector>
#include "qp_block.hxx"


namespace ConicBundle {

//for given psd matrix Q>0, and vector c solve the block quadratic program
//
//   max  -1/2<Qx,x>+<c,x>+offset      min 1/2<Qx,x>+<b,y>+offset
//        A_ix_i + B_is_i = b_i            Qx+ A^Ty - z = c
//        x_i\in K_i1, s_i\in K_i2             B^Ty - t = 0
//        x=[x_i]                          z_i\in K_i1^*, t_i\in K_i2^*       
//
//  the blocks i are completely independent 
//  the main code only works on variables x,y,z, coefficients A_i, and b_i 
//  the variables s_i and coefficients B_i are local in each block 
//  and have to be treated there (see block definition below).
//  The feasible set of each block has to be full dimensional and compact,
//  each block has to provide a strictly feasible primal-dual starting point,
//  do its line search etc.

// ****************************************************************************
//                                   QP_SDPBlock
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
// find a common step size and move on

class QP_SDPBlock: public QP_Block
{
private:
  CH_Matrix_Classes::Matrix x;                    //joint local primal variables as a vector 
                               //x=(vecx',svec(X1)',svec(X2)',...,svec(Xk)')'
  CH_Matrix_Classes::Matrix z;                    //joint local dual slack variables as a vector
                               //z=(vecz',svec(Z1)',svec(Z2)',...,svec(Zk)')'
  CH_Matrix_Classes::Real y;                      //dual variable
  CH_Matrix_Classes::Matrix A;                    //here, A is only a row vector
                               //A=(ones(vecz)',svec(I1)',....,svec(Ik)')
  CH_Matrix_Classes::Real b;                      //rhs 
  int less_or_equal;           //if !=0, then Ax <= b
  CH_Matrix_Classes::Real s;                      //primal slack variable if inequality
 
  CH_Matrix_Classes::Integer lin_dim;             //dimension of leading linear part

  CH_Matrix_Classes::Indexmatrix soc_dim;         //dimensions of second order cone variables
  CH_Matrix_Classes::Indexmatrix soc_start;       //first element of i-th SOC variable in x

  std::vector<CH_Matrix_Classes::Symmatrix> Xp;   //local primal psd variables X 
  std::vector<CH_Matrix_Classes::Symmatrix> Zp;   //local dual slack psd variables Z
  std::vector<CH_Matrix_Classes::Symmatrix> Xinv; //inverses of local primal psd variables X 
  std::vector<CH_Matrix_Classes::Symmatrix> Xchol; //inverses of local primal psd variables X 
  std::vector<CH_Matrix_Classes::Symmatrix> Zchol; //inverses of local primal psd variables Z


  CH_Matrix_Classes::Integer qp_xstart;            //the first index of the local x in qp_x
  CH_Matrix_Classes::Integer qp_ystart;            //the index of the local y in qp_y

  CH_Matrix_Classes::Integer mu_dim;               //= vecx.dim + sum of the matrix orders [ + 1]
  CH_Matrix_Classes::Real restart_factor;

  //---- copies of old variables for analyzing activity 
  CH_Matrix_Classes::Matrix old_x;                   
  CH_Matrix_Classes::Matrix old_z;                   
  CH_Matrix_Classes::Real old_y;                     
  CH_Matrix_Classes::Real old_s;   
  std::vector<CH_Matrix_Classes::Symmatrix> old_Xp; 
  std::vector<CH_Matrix_Classes::Symmatrix> old_Zp; 
 

  //---- temporary variables
  CH_Matrix_Classes::Matrix tmpvec;
  CH_Matrix_Classes::Matrix tmpsvec;
  CH_Matrix_Classes::Symmatrix tmpsymmat;
  CH_Matrix_Classes::Symmatrix dX;
 
  std::ostream* out;
  int print_level;
  //the Block has to memorize and update its local x,y,z variables 
  //(and all additional ones that are not visible to the outside)
  //in particular it is assumed to keep a copy of the current variables
  //after the calls to the following functions:
  //  starting_x, starting_yz, restart_x, restart_yz, set_point

  CH_Matrix_Classes::Matrix mult_Arw(const CH_Matrix_Classes::Matrix& x,const CH_Matrix_Classes::Matrix& v) const;
  CH_Matrix_Classes::Matrix mult_Arwinv(const CH_Matrix_Classes::Matrix& x,CH_Matrix_Classes::Real gamma,const CH_Matrix_Classes::Matrix& v) const;
  CH_Matrix_Classes::Matrix mult_G(const CH_Matrix_Classes::Matrix& x,CH_Matrix_Classes::Real gamma,const CH_Matrix_Classes::Matrix& v) const;
  CH_Matrix_Classes::Matrix mult_Ginv(const CH_Matrix_Classes::Matrix& x,CH_Matrix_Classes::Real gamma,const CH_Matrix_Classes::Matrix& v) const;



public:
  ~QP_SDPBlock(){}  

  int init_block(CH_Matrix_Classes::Integer lin_dim,const CH_Matrix_Classes::Indexmatrix& soc_dim,const CH_Matrix_Classes::Indexmatrix& sdp_dim,CH_Matrix_Classes::Real b,int less_or_equal); 
  int adjust_trace(CH_Matrix_Classes::Real b); 
  int get_linx(CH_Matrix_Classes::Matrix& vecx){ vecx.init(lin_dim,1,x.get_store());return 0;}
  int get_linz(CH_Matrix_Classes::Matrix& vecz){ vecz.init(lin_dim,1,z.get_store());return 0;}
  int get_socx(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Matrix& vecx)
  { vecx.init(soc_dim(i),1,x.get_store()+soc_start(i));return 0;}
  int get_socz(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Matrix& vecz)
  { vecz.init(soc_dim(i),1,z.get_store()+soc_start(i));return 0;}
  int get_X(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Symmatrix& X){X.init(Xp[i]); return 0;}
  int get_Z(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Symmatrix& Z){Z.init(Zp[i]); return 0;}
  CH_Matrix_Classes::Real get_y(void){return y;}
  CH_Matrix_Classes::Real get_s(void){if (less_or_equal) return s; return 0.;}

  int get_old_linx(CH_Matrix_Classes::Matrix& vecx){ vecx.init(lin_dim,1,old_x.get_store());return 0;}
  int get_old_linz(CH_Matrix_Classes::Matrix& vecz){ vecz.init(lin_dim,1,old_z.get_store());return 0;}
  int get_old_socx(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Matrix& vecx)
  { vecx.init(soc_dim(i),1,old_x.get_store()+soc_start(i));return 0;}
  int get_old_socz(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Matrix& vecz)
  { vecz.init(soc_dim(i),1,old_z.get_store()+soc_start(i));return 0;}
  int get_old_X(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Symmatrix& X){X.init(old_Xp[i]); return 0;}
  int get_old_Z(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Symmatrix& Z){Z.init(old_Zp[i]); return 0;}
  CH_Matrix_Classes::Real get_old_y(void){return old_y;}
  CH_Matrix_Classes::Real get_old_s(void){if (less_or_equal) return old_s; return 0.;}

  int inner_line_search(CH_Matrix_Classes::Real& alpha,
	const CH_Matrix_Classes::Matrix& qp_dx,const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz,CH_Matrix_Classes::Real ds);

  //-----------  QP_Block routines (see there)

  CH_Matrix_Classes::Integer xdim(){return x.dim();}
  CH_Matrix_Classes::Integer ydim(){return 1;}
  
  int set_qp_xstart(CH_Matrix_Classes::Integer x_start_index)
    {qp_xstart=x_start_index; return 0;}
  
  int set_qp_ystart(CH_Matrix_Classes::Integer y_start_index)
    {qp_ystart=y_start_index; return 0;}

  int starting_x(CH_Matrix_Classes::Matrix& qp_x);

  int starting_yz(CH_Matrix_Classes::Matrix& qp_y,CH_Matrix_Classes::Matrix& qp_z,
		   const CH_Matrix_Classes::Matrix& qp_Qx, const CH_Matrix_Classes::Matrix& qp_c);

  int get_Ab(CH_Matrix_Classes::Matrix& qp_A,CH_Matrix_Classes::Matrix &qp_b) const;

  int restart_x(CH_Matrix_Classes::Matrix& qp_x,const CH_Matrix_Classes::Matrix& qp_c,const CH_Matrix_Classes::Matrix& qp_dc);

  int restart_yz(CH_Matrix_Classes::Matrix& qp_y,CH_Matrix_Classes::Matrix& qp_z,
		  const CH_Matrix_Classes::Matrix& qp_Qx,const CH_Matrix_Classes::Matrix& qp_c,const CH_Matrix_Classes::Matrix& qp_dc);
    
  int add_xinv_kron_z(CH_Matrix_Classes::Symmatrix& barQ);

  int add_local_sys(CH_Matrix_Classes::Symmatrix& sysdy,CH_Matrix_Classes::Matrix& rhs);

  int suggest_mu(CH_Matrix_Classes::Real& ip_xz,CH_Matrix_Classes::Integer& mu_dim,CH_Matrix_Classes::Real& sigma,
                 const CH_Matrix_Classes::Matrix& qp_dx,const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz);

  int get_corr(CH_Matrix_Classes::Matrix& xcorr,CH_Matrix_Classes::Matrix& rhs,CH_Matrix_Classes::Real mu,const CH_Matrix_Classes::Matrix& qp_dx,
		       const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz);

  int line_search(CH_Matrix_Classes::Real& alpha,const CH_Matrix_Classes::Matrix& qp_dx,
                          const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz);

  int set_point(const CH_Matrix_Classes::Matrix& qp_x,const CH_Matrix_Classes::Matrix& qp_y,const CH_Matrix_Classes::Matrix& qp_z,
			CH_Matrix_Classes::Real alpha);

  void set_out(std::ostream* o=0,int pril=1)
  {out=o;print_level=pril;}

  //---------------- for debugging purposes

  virtual int add_Bs(CH_Matrix_Classes::Matrix& qp_vec) const;

};

}

#endif

