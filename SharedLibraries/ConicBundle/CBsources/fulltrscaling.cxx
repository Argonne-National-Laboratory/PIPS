/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/fulltrscaling.cxx

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



#include <math.h>
#include <stdlib.h>
#include "fulltrscaling.hxx"
#include "mymath.hxx"
#include "sparssym.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                       BundleFullTrustRegionScaling::set_weightu
// *****************************************************************************

void BundleFullTrustRegionScaling::set_weightu(CH_Matrix_Classes::Real in_weightu)
{
  if (fabs(weightu-in_weightu)>1e-10*fabs(weightu)){
    is_factored=false;
    weightu=in_weightu;
  }
}

// *****************************************************************************
//                       BundleFullTrustRegionScaling::norm_sqr
// *****************************************************************************

Real BundleFullTrustRegionScaling::norm_sqr(const Matrix& B) const
{
  return ip(B,H*B)+weightu*ip(B,B);
}


// *****************************************************************************
//                       BundleFullTrustRegionScaling::dnorm_sqr
// *****************************************************************************

Real BundleFullTrustRegionScaling::dnorm_sqr(const Matrix& B) const
{
  if (!is_factored){
    is_factored=true;
    Hchol=H;
    for(Integer i=0;i<Hchol.rowdim();i++)
      Hchol(i,i)+=weightu;
    if (Hchol.Chol_factor()){
      if (out)
	(*out)<<"ERROR in BundleFullTrustRegionScaling::update_eta_step(...): H.Chol_factor() failed"<<std::endl;
      return 1;
    }
  }
  Matrix tmpmat(B);
  Hchol.Chol_Lsolve(tmpmat);
  return ip(tmpmat,tmpmat);
}


// *****************************************************************************
//                       BundleFullTrustRegionScaling::compute_first_eta
// *****************************************************************************


int BundleFullTrustRegionScaling::compute_first_eta(Matrix& ,//eta,
					 const Matrix& ,//subg,
					 const Matrix& ,//y,
					 const Matrix& ,//lby,
					 const Matrix& ,//uby,
					 const Indexmatrix& bound_index,
					 const Indexmatrix& // yfixed
					   ) const
{
  //currently box constraints are not implemented for full scaling
  if (bound_index.dim()>0){
    if (out)
      (*out)<<"ERROR in BundleFullTrustRegionScaling::compute_first_eta(...): box constraints not implemented"<<std::endl;
    return 1;
  }
  
  return 0;
}
  

// *****************************************************************************
//                       BundleFullTrustRegionScaling::update_eta_step
// *****************************************************************************


int BundleFullTrustRegionScaling::update_eta_step(Matrix& newy,
				       Matrix& ,  //eta,
				       Indexmatrix& update_index,
				       Matrix& update_value,
				       Real& normsubg2,
				       const Matrix& subg,
				       const Matrix& y,
				       const Matrix& , //lby,
				       const Matrix& , //uby,
				       const Indexmatrix& bound_index,
				       const Indexmatrix& yfixed
				       ) const
{
  if (!is_factored){
    is_factored=true;
    Hchol=H;
    for(Integer i=0;i<Hchol.rowdim();i++)
      Hchol(i,i)+=weightu;
    if (Hchol.Chol_factor()){
      if (out)
	(*out)<<"ERROR in BundleFullTrustRegionScaling::update_eta_step(...): H.Chol_factor() failed"<<std::endl;
      return 1;
    }
  }
  newy.init(subg,-1.);
  Hchol.Chol_solve(newy);
  if (yfixed.dim()==newy.dim()){
    for(Integer i=0;i<yfixed.dim();i++){
      if (yfixed(i)) 
	newy(i)=0.;
    }
  }
  //as long as there is no eta, it is better to compute normsubg2 right now
  Matrix tmpmat(newy);
  Hchol.Chol_Ltmult(tmpmat);
  assert(norm2(tmpmat-triu(Matrix(Hchol))*newy)<1e-6);
  normsubg2=ip(tmpmat,tmpmat);
  
  newy+=y;
   

  //--- update eta, if necessary 
  //    (so that the step satisfies the box constraints)

  if (bound_index.dim()==0){
    update_value.newsize(0,0);
    update_index.newsize(0,0);
  }
  else {
    //currently box constraints are not implemented for full scaling
    if (out)
      (*out)<<"ERROR in BundleFullTrustRegionScaling::update_eta_step(...): box constraints not implemented"<<std::endl;
    return 1;
  }
    
  /* skip this as long as there is no eta
  Matrix tmpmat;
  tmpmat.xeymz(newy,y);
  Hchol.Chol_Ltmult(tmpmat);
  normsubg2=ip(tmpmat,tmpmat);
  newy+=y;
  */

  return 0;
}    


// *****************************************************************************
//                       BundleFullTrustRegionScaling::compute_QP_costs
// *****************************************************************************


int BundleFullTrustRegionScaling::compute_QP_costs(Symmatrix& Q,
				      Matrix& d,
				      Real& offset,
				      const BundleQPData* datap,
				      const Integer xdim,
				      const Matrix& y,
				      const Matrix& lby,
				      const Matrix& uby,
				      const Matrix& eta,
				      const Indexmatrix& yfixed) const
{
  if (!is_factored){
    is_factored=true;
    Hchol=H;
    for(Integer i=0;i<Hchol.rowdim();i++)
      Hchol(i,i)+=weightu;
    if (Hchol.Chol_factor()){
    if (out)
      (*out)<<"*** ERROR in BundleFullTrustRegionScaling::compute_QP_costs(...): H.Chol_factor failed"<<std::endl;
      return 1;
    }
  }

  //--- initialize linear term and offset
  d.init(xdim,1,0.);
  offset=0;
  if (datap->get_row(-1,d,offset,0)){
    if (out)
      (*out)<<"*** ERROR in BundleFullTrustRegionScaling::compute_QP_costs(...):  datap->get_row(...) failed"<<std::endl;
    return 1;
  }

  //--- collect the bundle information
  bigA.init(y.dim(),xdim,0.);
  rhs.init(y.dim(),1,0.);
  Matrix tmpvec(xdim,1,0.);
  bool yf= (yfixed.dim()==y.dim());

  for(Integer j=0;j<y.dim();j++){
    Real bb=eta(j);
    if (bb<0.) offset+=bb*uby(j);
    else if (bb>0.) offset+=bb*lby(j);
    
    //get cost row corresponding to j
    if (yf && yfixed(j)){
      if (y(j)!=0.){ //(for y==0. there is nothing to do)
	bb=0.;
	if (datap->get_row(j,tmpvec,bb,0)){
	  if (out)
	    (*out)<<"*** ERROR in BundleFullTrustRegionScaling::compute_QP_costs(...):  datap->get_row(...) failed"<<std::endl;
	  return 1;
	}
	
	//add offset and linear term
	offset+=bb*y(j);
	d.xpeya(tmpvec,y(j));
      }
    }
    else{
      bb=-bb;
      if (datap->get_row(j,tmpvec,bb,0)){
	if (out)
	  (*out)<<"*** ERROR in BundleFullTrustRegionScaling::compute_QP_costs(...):  datap->get_row(...) failed"<<std::endl;
	return 1;
      }
      
      offset+=bb*y(j);      
      d.xpeya(tmpvec,y(j));
      
      rhs(j)=bb;
      for(Integer i=0;i<xdim;i++){
	bigA(j,i)=-tmpvec(i);
      }
    }
  }    
  Hchol.Chol_Lsolve(rhs);
  Hchol.Chol_Lsolve(bigA);
  offset-=ip(rhs,rhs)/2.;
  genmult(bigA,rhs,tmpvec,1.,0.,1,0);
  d+=tmpvec;
  rankadd(bigA,Q,1.,0.,1);
 
  //cout<<"offset="<<offset<<std::endl;
  //cout<<"d="<<d<<std::endl;
  //cout<<"Q="<<Q<<std::endl;

  return 0;
}


// *****************************************************************************
//                       BundleFullTrustRegionScaling::update_QP_costs
// *****************************************************************************

int BundleFullTrustRegionScaling::update_QP_costs(
					   Matrix& delta_d,  
					   Real& delta_offset,
					   const BundleQPData* /* datap */,
					   const Integer /* xdim */,
					   const Matrix& y,
					   const Matrix& lby,
					   const Matrix& uby,
					   const Matrix& eta,
					   const Indexmatrix& update_index,
					   const Matrix& update_value) const
{
  if (!is_factored){
    is_factored=true;
    Hchol=H;
    for(Integer i=0;i<Hchol.rowdim();i++)
      Hchol(i,i)+=weightu;
    if (Hchol.Chol_factor()){
      if (out)
	(*out)<<"*** ERROR in BundleFullTrustRegionScaling::compute_QP_costs(...): H.Chol_factor failed"<<std::endl;
      return 1;
    }
  }

  Matrix delta_rhs(y.dim(),1,0.);
  delta_offset=0;

  //--- for each y_index 
  for (Integer j=0;j<update_index.dim();j++){
    Integer ind=update_index(j);
    
    Real uv=update_value(j);
    Real new_eta=eta(ind);
    Real old_eta=new_eta-uv;
    if (new_eta<0.) {
      if (old_eta>0.){
	delta_offset+=new_eta*uby(ind)-old_eta*lby(ind);
      }
      else if (old_eta==0.){
	delta_offset+=new_eta*uby(ind);
      }
    }
    else if (new_eta>0.) {
      if (old_eta<0.){
	delta_offset+=new_eta*lby(ind)-old_eta*uby(ind);
      }
      else if (old_eta==0.){
	delta_offset+=new_eta*lby(ind);
      }
    }

    //get row corresponding to j
    delta_rhs(ind)=uv;
    delta_offset-=uv*y(ind);
      
  }
  delta_offset+=ip(rhs,rhs)/2.;
  Hchol.Chol_solve(delta_rhs);
  rhs+=delta_rhs;
  delta_offset-=ip(rhs,rhs)/2.;
  genmult(bigA,delta_rhs,delta_d,-1.,0.,1,0);
  
  return 0;
}



}  //namespace ConicBundle
