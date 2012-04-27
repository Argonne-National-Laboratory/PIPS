/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/lowranktrscaling.cxx

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
#include "lowranktrscaling.hxx"
#include "idscaling.hxx"
#include "mymath.hxx"
#include "sparssym.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::init
// *****************************************************************************

void BundleLowRankTrustRegionScaling::init(const CH_Matrix_Classes::Matrix& in_vecH, 
	    const CH_Matrix_Classes::Matrix& in_lamH, 
	    const CH_Matrix_Classes::Real& in_regterm)
{
  assert(vecH.coldim()==lamH.dim());
  vecH=in_vecH; 
  lamH=in_lamH; 
  regterm=in_regterm; 
  lamHi.init(Integer(0),Integer(0),0.);
  QRvecH=vecH;
  if (lamH.dim()>0){
    QRvecH.QR_factor();
  }    
}


// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::set_weightu
// *****************************************************************************

void BundleLowRankTrustRegionScaling::set_weightu(CH_Matrix_Classes::Real in_weightu)
{
  if (fabs(weightu-in_weightu)>1e-10*fabs(weightu)){
    weightu=in_weightu;
    lamHi.init(Integer(0),Integer(0),0.);
  }
}

// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::norm_sqr
// *****************************************************************************

Real BundleLowRankTrustRegionScaling::norm_sqr(const Matrix& B) const
{
  if (lamH.dim()==0){
    return (regterm+weightu)*ip(B,B);
  }
  Matrix tmpmat;
  genmult(vecH,B,tmpmat,1.,0.,1);
  return (regterm+weightu)*ip(B,B)+normDsquared(tmpmat,lamH);
}


// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::dnorm_sqr
// *****************************************************************************

Real BundleLowRankTrustRegionScaling::dnorm_sqr(const Matrix& B) const
{
  if (lamH.dim()==0){
    return ip(B,B)/(regterm+weightu);
  }
  if (lamHi.dim()==0){
    lamHi.init(vecH.rowdim(),1,regterm+weightu);
    for(Integer i=0;i<lamH.dim();i++){
      lamHi(i)+=lamH(i);
    }
    lamHi.inv();
  }
  Matrix tmpmat(B);
  QRvecH.Qt_times(tmpmat,QRvecH.coldim());

  return normDsquared(tmpmat,lamHi);
}


// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::compute_first_eta
// *****************************************************************************


int BundleLowRankTrustRegionScaling::compute_first_eta(Matrix& eta,
					 const Matrix& subg,
					 const Matrix& y,
					 const Matrix& lby,
					 const Matrix& uby,
					 const Indexmatrix& bound_index,
					 const Indexmatrix& yfixed
					   ) const
{
  if (lamH.dim()==0){
    BundleIdScaling ids;
    ids.set_weightu(weightu+regterm);
    return ids.compute_first_eta(eta,subg,y,lby,uby,bound_index,yfixed);
  }
  if (lamHi.dim()==0){
    lamHi.init(vecH.rowdim(),1,regterm+weightu);
    for(Integer i=0;i<lamH.dim();i++){
      lamHi(i)+=lamH(i);
    }
    lamHi.inv();
  }

  //currently box constraints are not implemented for full scaling
  if (bound_index.dim()>0){
    if (out)
      (*out)<<"*** ERROR in BundleLowRankTrustRegionScaling::compute_first_eta(...): box constraints not implemented"<<std::endl;
    return 1;
  }
  
  return 0;
}
  

// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::update_eta_step
// *****************************************************************************


int BundleLowRankTrustRegionScaling::update_eta_step(Matrix& newy,
				       Matrix& eta,
				       Indexmatrix& update_index,
				       Matrix& update_value,
				       Real& normsubg2,
				       const Matrix& subg,
				       const Matrix& y,
				       const Matrix& lby,
				       const Matrix& uby,
				       const Indexmatrix& bound_index,
				       const Indexmatrix& yfixed
				       ) const
{
  if (lamH.dim()==0){
    BundleIdScaling ids;
    ids.set_weightu(weightu+regterm);
    return ids.update_eta_step(newy,eta,update_index,update_value,normsubg2,
			       subg,y,lby,uby,bound_index,yfixed);
  }
  if (lamHi.dim()==0){
    lamHi.init(vecH.rowdim(),1,regterm+weightu);
    for(Integer i=0;i<lamH.dim();i++){
      lamHi(i)+=lamH(i);
    }
    lamHi.inv();
  }

  //compute newy=-inv(H)*subg
  newy.init(subg,-1.);
  QRvecH.Qt_times(newy,QRvecH.coldim());
  newy.scale_rows(lamHi);
  QRvecH.Q_times(newy,QRvecH.coldim());
  if (yfixed.dim()==newy.dim()){
    for(Integer i=0;i<yfixed.dim();i++){
      if (yfixed(i)) 
	newy(i)=0.;
    }
  }
  //as long as there is no eta, it is better to compute normsubg2 right now
  Matrix tmpmat;
  genmult(vecH,newy,tmpmat,1.,0.,1);
  normsubg2=((regterm+weightu)*ip(newy,newy)+normDsquared(tmpmat,lamH));

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
      (*out)<<"*** ERROR in BundleLowRankTrustRegionScaling::update_eta_step(...): box constraints not implemented"<<std::endl;
    return 1;
  }
    
  /* skip this as long as there is no eta
  newy-=y;
  Matrix tmpmat;
  genmult(vecH,newy,tmpmat,1.,0.,1);
  normsub2=((regterm+weightu)*ip(B,B)+normDsquared(tmpmat,lamH));
  newy+=y;
  */

  return 0;
}    


// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::compute_QP_costs
// *****************************************************************************


int BundleLowRankTrustRegionScaling::compute_QP_costs(Symmatrix& Q,
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
  if (lamH.dim()==0){
    BundleIdScaling ids;
    ids.set_weightu(weightu+regterm);
    return ids.compute_QP_costs(Q,d,offset,datap,xdim,y,lby,uby,eta,yfixed);
  }
  if (lamHi.dim()==0){
    lamHi.init(vecH.rowdim(),1,regterm+weightu);
    for(Integer i=0;i<lamH.dim();i++){
      lamHi(i)+=lamH(i);
    }
    lamHi.inv();
  }

  //--- initialize coefficient matrices
  Q.init(xdim,0.);
  d.init(xdim,1,0.);
  offset=0;
  if (datap->get_row(-1,d,offset,0)){
    if (out)
      (*out)<<"*** ERROR in BundleLowRankScaing::compute_QP_costs(...): datap->get_row(...) failed"<<std::endl;
    return 1;
  }
  //cout<<"costd="<<transpose(d)<<std::endl;
  //cout<<"offset="<<offset<<std::endl;

  //--- collect the bundle information
  Matrix bigA(y.dim(),xdim,0.);
  Matrix bigb(y.dim(),1,0.);
 
  Matrix rowvecHjt;
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
	    (*out)<<"*** ERROR in BundleLowRankScaing::compute_QP_costs(...): datap->get_row(...) failed"<<std::endl;
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
	  (*out)<<"*** ERROR in BundleLowRankScaing::compute_QP_costs(...): datap->get_row(...) failed"<<std::endl;
	return 1;
      }

      bigb(j)=bb;
      for(Integer i=0;i<xdim;i++){
	bigA(j,i)=tmpvec(i);
      }
      offset+=bb*(y(j));
      d.xpeya(tmpvec,y(j));
    }

  }    

  tmpvec.init(lamHi);
  tmpvec.sqrt();

  vecHtA=bigA;
  QRvecH.Qt_times(vecHtA,QRvecH.coldim());
  vecHtA.scale_rows(tmpvec);


  vecHtb=bigb;
  QRvecH.Qt_times(vecHtb,QRvecH.coldim());
  vecHtb.scale_rows(tmpvec);
  
  offset-=ip(vecHtb,vecHtb)/2.;
  genmult(vecHtA,vecHtb,d,-1.,1.,1,0);  //d-=vecHtA^T*vecHtb    
  rankadd(vecHtA,Q,1.,0.,1);   //Q=tmpmatvecHtA^T*tmpmatvecHtA


  //cout<<"offset="<<offset<<std::endl;
  //cout<<"d="<<d<<std::endl;
  //cout<<"Q="<<Q<<std::endl;
  
  return 0;
}


// *****************************************************************************
//                       BundleLowRankTrustRegionScaling::update_QP_costs
// *****************************************************************************

//not yet tested! Mistakes likely!

int BundleLowRankTrustRegionScaling::update_QP_costs(
					   Matrix& delta_d,  
					   Real& delta_offset,
					   const BundleQPData* datap,
					   const Integer xdim,
					   const Matrix& y,
					   const Matrix& lby,
					   const Matrix& uby,
					   const Matrix& eta,
					   const Indexmatrix& update_index,
					   const Matrix& update_value) const
{
  if (lamH.dim()==0){
    BundleIdScaling ids;
    ids.set_weightu(weightu+regterm);
    return ids.update_QP_costs(delta_d,delta_offset,datap,xdim,
			       y,lby,uby,eta,update_index,update_value);
  }
  if (lamHi.dim()==0){
    lamHi=lamH;
    lamHi+=(regterm+weightu);
    lamHi.inv();
    lamHi*=(regterm+weightu);
    lamHi%=lamH;
  }

  delta_d.init(xdim,1,0.);
  Matrix tmpvec(xdim,1,0.);
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
    Real bb=0.;
    if (datap->get_row(ind,tmpvec,bb,0)){
      if (out)
	(*out)<<"*** ERROR in BundleLowRankTrustRegionScaling::update_QP_costs(...):  datap->get_row(...) failed"<<std::endl;
	return 1;
    }
    delta_offset-=uv*y(ind);
    delta_rhs(ind)=-uv;
  }
  delta_offset-=ip(vecHtb,vecHtb)/2.;
  QRvecH.Qt_times(delta_rhs,QRvecH.coldim());
  tmpvec.init(lamHi);
  tmpvec.sqrt();
  delta_rhs.scale_rows(tmpvec);
  vecHtb+=delta_rhs;
  delta_offset+=ip(vecHtb,vecHtb)/2.;
  genmult(vecHtA,delta_rhs,delta_d,-1.,0.,1,0);  //d=vecHtA^T*delta_rhs
  
  return 0;
}



}  //namespace ConicBundle
