/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/qp_sumblock.hxx

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



#ifndef CONICBUNDLE_QP_SUMBLOCK_HXX
#define CONICBUNDLE_QP_SUMBLOCK_HXX

#include "qp_block.hxx"

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
//                                   QP_SUMBLOCK
// ****************************************************************************

// corresponds to a block of blocks in the quadratic program

class QP_SumBlock: public QP_Block
{
private:
  std::vector<QP_Block*> blocks;

  CH_Matrix_Classes::Integer xstart;
  CH_Matrix_Classes::Integer xend;
  CH_Matrix_Classes::Integer ystart;
  CH_Matrix_Classes::Integer yend;

  std::ostream* out;
  int print_level;

public:
  ~QP_SumBlock(){out=0;print_level=0;}  

  int clear_blocks(){blocks.clear(); return 0;}
  int add_block(QP_Block* qp_block){blocks.push_back(qp_block); return 0;}

  //-----------  QP_Block routines

  CH_Matrix_Classes::Integer xdim(){
    CH_Matrix_Classes::Integer dim=0;
    for(unsigned int i=0;i<blocks.size();i++) dim+=blocks[i]->xdim();
    return dim;
  }

  CH_Matrix_Classes::Integer ydim(){
    CH_Matrix_Classes::Integer dim=0;
    for(unsigned int i=0;i<blocks.size();i++) dim+=blocks[i]->ydim();
    return dim;
  }
  
  int set_qp_xstart(CH_Matrix_Classes::Integer x_start_index){
    int retval=0;
    xstart=x_start_index;
    xend=xstart;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->set_qp_xstart(xend);
      xend+=blocks[i]->xdim();
    }
    return retval;
  }
  
  int set_qp_ystart(CH_Matrix_Classes::Integer y_start_index){
    int retval=0;
    ystart=y_start_index;
    yend=ystart;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->set_qp_ystart(yend);
      yend+=blocks[i]->ydim();
    }
    return retval;
  }

  int starting_x(CH_Matrix_Classes::Matrix& qp_x){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->starting_x(qp_x);
    }
    return retval;
  }

  int starting_yz(CH_Matrix_Classes::Matrix& qp_y,CH_Matrix_Classes::Matrix& qp_z,
		   const CH_Matrix_Classes::Matrix& qp_Qx, const CH_Matrix_Classes::Matrix& qp_c){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->starting_yz(qp_y,qp_z,qp_Qx,qp_c);
    }
    return retval;
  }
 


  int get_Ab(CH_Matrix_Classes::Matrix& qp_A,CH_Matrix_Classes::Matrix &qp_b) const {
    int retval=0;
    for(CH_Matrix_Classes::Integer j=xstart;j<xend;j++){ //set entire block to zero
      CH_Matrix_Classes::mat_xea(yend-ystart,qp_A.get_store()+qp_A.rowdim()*j+ystart,0.);
    }
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->get_Ab(qp_A,qp_b);
    }
    return retval;
  }
      
    

  int restart_x(CH_Matrix_Classes::Matrix& qp_x,const CH_Matrix_Classes::Matrix& qp_c,const CH_Matrix_Classes::Matrix& qp_dc){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->restart_x(qp_x,qp_c,qp_dc);
    }
    return retval;
  }
    

  int restart_yz(CH_Matrix_Classes::Matrix& qp_y,CH_Matrix_Classes::Matrix& qp_z,
		 const CH_Matrix_Classes::Matrix& qp_Qx,const CH_Matrix_Classes::Matrix& qp_c,const CH_Matrix_Classes::Matrix& qp_dc){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->restart_yz(qp_y,qp_z,qp_Qx,qp_c,qp_dc);
    }
    return retval;
  }
    
  int add_xinv_kron_z(CH_Matrix_Classes::Symmatrix& barQ){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->add_xinv_kron_z(barQ);
    }
    return retval;
  }
    

  int add_local_sys(CH_Matrix_Classes::Symmatrix& sysdy,CH_Matrix_Classes::Matrix& rhs){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->add_local_sys(sysdy,rhs);
    }
    return retval;
  }
    

  int suggest_mu(CH_Matrix_Classes::Real& ip_xz,CH_Matrix_Classes::Integer& mu_dim,CH_Matrix_Classes::Real& sigma,
                 const CH_Matrix_Classes::Matrix& qp_dx,const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz){
    int retval=0;
    ip_xz=0.;
    mu_dim=0;
    sigma=0.;
    CH_Matrix_Classes::Real ipxz,s;
    CH_Matrix_Classes::Integer md;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->suggest_mu(ipxz,md,s,qp_dx,qp_dy,qp_dz);
      ip_xz+=ipxz;
      mu_dim+=md;
      sigma=CH_Matrix_Classes::max(s,sigma);
    }
    return retval;
  }
    

  int get_corr(CH_Matrix_Classes::Matrix& xcorr,CH_Matrix_Classes::Matrix& rhs,CH_Matrix_Classes::Real mu,const CH_Matrix_Classes::Matrix& qp_dx,
	       const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->get_corr(xcorr,rhs,mu,qp_dx,qp_dy,qp_dz);
    }
    return retval;
  }
       

  int line_search(CH_Matrix_Classes::Real& alpha,const CH_Matrix_Classes::Matrix& qp_dx,
		  const CH_Matrix_Classes::Matrix& qp_dy,const CH_Matrix_Classes::Matrix& qp_dz){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->line_search(alpha,qp_dx,qp_dy,qp_dz);
    }
    return retval;
  }
  

  int set_point(const CH_Matrix_Classes::Matrix& qp_x,const CH_Matrix_Classes::Matrix& qp_y,const CH_Matrix_Classes::Matrix& qp_z,
			CH_Matrix_Classes::Real alpha){
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->set_point(qp_x,qp_y,qp_z,alpha);
    }
    return retval;
  }
    
  void set_out(std::ostream* o=0,int pril=1)
  { out=o; print_level=pril; }

  //---------------- for debugging purposes

  virtual int add_Bs(CH_Matrix_Classes::Matrix& qp_vec) const{
    int retval=0;
    for(unsigned int i=0;i<blocks.size();i++) {
      retval|=blocks[i]->add_Bs(qp_vec);
    }
    return retval;
  }


};

}

#endif

