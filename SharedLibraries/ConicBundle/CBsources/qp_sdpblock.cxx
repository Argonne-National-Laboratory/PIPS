/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/qp_sdpblock.cxx

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



#include "mymath.hxx"
#include "qp_sdpblock.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// ****************************************************************************
//                 routines for computing the corrector of SOC
// ****************************************************************************

// x0   = x(0);        (must be >0.)
// barx = x(1..n-1);
//
// Arw(x)=[ x0, barx'; barx, x0*I ]
//
// gamma = sqrt(x0*x0 - barx'*barx);  (must be >0.)
//
// Arwinv(x) = (1/gamma)*[ x0, -barx'; -barx', (barx*barx'+gamma*I)/x0]:
// 
// G(x)= [x0, barx'; barx, gamma*I+barx*barx'/(gamma+x0)];
//
// Ginv(x)= [-1, 0;0, I]/gamma  -  1/( gamma*( barx'*barx / (gamma*x0)-x0)) *
//          [gamma+x0, -barx'; barx, barx*barx'/(gamma+x0)];
//


// Arw(x)*v

  Matrix QP_SDPBlock::mult_Arw(const Matrix& x,const Matrix& v) const
  {
    chk_add(x,v);
    Matrix w(x,v(0));
    mat_xpeya(w.dim()-1,w.get_store()+1,v.get_store()+1,x(0));
    w(0)+=mat_ip(x.dim()-1,x.get_store()+1,v.get_store()+1);
    return w;
  }
  
  
// Arwinv(x)*v    
 
  Matrix QP_SDPBlock::mult_Arwinv(const Matrix& x,Real gamma,const Matrix& v) const
  {
    chk_add(x,v);
    Real x0=x(0);

    //assert(x0>0.);
    //assert(gamma>0.);
    //assert(fabs(gamma-sqrt(x0*x0-mat_ip(x.dim()-1,x.get_store()+1,x.get_store()+1)))<1e-10*gamma);
    
    Real ip_barxv=mat_ip(x.dim()-1,x.get_store()+1,v.get_store()+1);
    Matrix w(x.dim(),1); chk_set_init(w,1);
    w(0)=v(0)*x(0)-ip_barxv;
    mat_xeya(w.dim()-1,w.get_store()+1,x.get_store()+1,(-v(0)+ip_barxv/x0)/gamma);
    mat_xpeya(w.dim()-1,w.get_store()+1,v.get_store()+1,1./x0);
    return w;
  }

  Matrix QP_SDPBlock::mult_G(const Matrix& x,Real gamma,const Matrix& v) const
  {
    chk_add(x,v);
    Real x0=x(0);

    //assert(x0>0.);
    //assert(gamma>0.);
    //assert(fabs(gamma-sqrt(x0*x0-mat_ip(x.dim()-1,x.get_store()+1,x.get_store()+1)))<1e-10*gamma);
    
    
    Real ip_barxv=mat_ip(x.dim()-1,x.get_store()+1,v.get_store()+1);
    Matrix w(x,v(0));
    w(0)+=ip_barxv;
    mat_xpeya(w.dim()-1,w.get_store()+1,x.get_store()+1,ip_barxv/(x0+gamma));
    mat_xpeya(w.dim()-1,w.get_store()+1,v.get_store()+1,gamma);
    return w;
  }
   
  Matrix QP_SDPBlock::mult_Ginv(const Matrix& x,Real gamma,const Matrix& v) const
  {
    chk_add(x,v);
    Real x0=x(0);

    //assert(x0>0.);
    //assert(gamma>0.);
    //assert(fabs(gamma-sqrt(x0*x0-mat_ip(x.dim()-1,x.get_store()+1,x.get_store()+1)))<1e-10*gamma);
    
    Real ip_barxv=mat_ip(x.dim()-1,x.get_store()+1,v.get_store()+1);
    Real d=-1./(gamma*((x0*x0-gamma*gamma)/(gamma+x0)-x0));
    Matrix w(x.dim(),1); chk_set_init(w,1);
    w(0)=v(0)*(-1./gamma+(gamma+x0)*d)-ip_barxv*d;
    mat_xeya(w.dim()-1,w.get_store()+1,v.get_store()+1,1./gamma);
    mat_xpeya(w.dim()-1,w.get_store()+1,x.get_store()+1,(ip_barxv/(gamma+x0)-v(0))*d);
    return w;
  }
    

// ****************************************************************************
//                           init_block
// ****************************************************************************

int QP_SDPBlock::init_block(Integer in_lin_dim,const Indexmatrix& in_soc_dim,const Indexmatrix& sdp_dim,Real in_b,int in_less_or_equal)
{
  lin_dim=in_lin_dim;
  soc_dim=in_soc_dim;
  b=in_b;
  less_or_equal=in_less_or_equal;
  Integer xdim=lin_dim;
  if (less_or_equal) {
    mu_dim=lin_dim+1;
  }
  else {
    mu_dim=lin_dim;
  }
  mu_dim+=soc_dim.dim();
  xdim+=sum(soc_dim);
  Xp.resize(sdp_dim.dim());
  Zp.resize(sdp_dim.dim());
  old_Xp.resize(sdp_dim.dim());
  old_Zp.resize(sdp_dim.dim());
  Xinv.resize(sdp_dim.dim());
  for(unsigned int i=0;i<Xp.size();i++){
    mu_dim+=sdp_dim(i);
    Xp[i].init(sdp_dim(i),0.);
    Xp[i].shift_diag(1.);
    xdim+=(sdp_dim(i)*(sdp_dim(i)+1))/2;
  }
  x.init(xdim,1,0.);
  z.init(xdim,1,0.);
  old_x.init(xdim,1,0.);
  old_z.init(xdim,1,0.);
  
  //--- construct A
  A.init(xdim,1,1.);
  xdim=lin_dim;
  soc_start.newsize(soc_dim.dim(),1); chk_set_init(soc_start,1);
  {for(Integer i=0;i<soc_dim.dim();i++){
    soc_start(i)=xdim;
    mat_xea(soc_dim(i)-1,A.get_store()+xdim+1,0.);
    xdim+=soc_dim(i);
  }}
  {for(unsigned int i=0;i<Xp.size();i++){
    svec(Xp[i],tmpsvec);
    mat_xey(tmpsvec.dim(),A.get_store()+xdim,tmpsvec.get_store());
  }}

  return 0;
}

// ****************************************************************************
//                           adjust_trace
// ****************************************************************************

int QP_SDPBlock::adjust_trace(Real in_b)
{
  b=in_b;
  return 0;
}


// ****************************************************************************
//                          inner_line_search
// ****************************************************************************

// the psd-line-search below should better be done by eigenvalue computation

int QP_SDPBlock::inner_line_search(Real& alpha,
   const Matrix& qp_dx,const Matrix& qp_dy,const Matrix& qp_dz,Real ds)
{
  Real al=alpha;
  Real bdfac=.9;  //boundary factor
  if (less_or_equal){
    if (ds<0.) al=min(al,-bdfac*s/ds);
    if ((y>0.)&&(qp_dy(qp_ystart)<0.)) {
      Real d=-y/qp_dy(qp_ystart);
      if (d<=al){
	al=bdfac*d;
      }
    }
  }

  //---linear part
  Integer xind=qp_xstart;
  {for(Integer i=0;i<lin_dim;i++){
    if (qp_dx(xind+i)<0){
      Real d=-x(i)/qp_dx(xind+i);
      if (d<=al){
	al=bdfac*d;
      }
    } 
    if (qp_dz(xind+i)<0){
      Real d=-z(i)/qp_dz(xind+i);
      if (d<=al){
	al=bdfac*d;
      }
    } 
  }}
  xind+=lin_dim;

  //---second order part
  for(Integer i=0;i<soc_dim.dim();i++){
    //primal variables
    const Real *xp=x.get_store()+xind-qp_xstart;
    const Real *dxp=qp_dx.get_store()+xind;
    Real a=(*dxp)*(*dxp);
    Real b=(*dxp)*(*xp);
    Real c=(*xp)*(*xp);
    xp++;
    dxp++;
    for(Integer j=soc_dim(i)-1;--j>=0;xp++,dxp++){
      a-=(*dxp)*(*dxp);
      b-=(*dxp)*(*xp);
      c-=(*xp)*(*xp);
    }
    if((fabs(a)<=eps_Real)&&(b<-eps_Real)){
      Real d=-c/b/2.;
      if (d<=al){
	al=bdfac*d;
      }
    } 
    else {
      Real d=b*b-a*c; 
      if ((a<-eps_Real)||((b<0.)&&(d>=0.))){
	Real f=(-b-sqrt(d))/a;
	if (f<=al){
	  al=bdfac*f;
	}
      } 
    }
    //dual variables
    xp=z.get_store()+xind-qp_xstart;
    dxp=qp_dz.get_store()+xind;
    a=(*dxp)*(*dxp);
    b=(*dxp)*(*xp);
    c=(*xp)*(*xp);
    xp++;
    dxp++;
    {for(Integer j=soc_dim(i)-1;--j>=0;xp++,dxp++){
      a-=(*dxp)*(*dxp);
      b-=(*dxp)*(*xp);
      c-=(*xp)*(*xp);
    }}
    if((fabs(a)<=eps_Real)&&(b<-eps_Real)){
      Real d=-c/b/2.;
      if (d<=al){
	al=bdfac*d;
      }
    } 
    else {
      Real d=b*b-a*c; 
      if ((a<-eps_Real)||((b<0.)&&(d>=0.))){
	Real f=(-b-sqrt(d))/a;
	if (f<=al){
	  al=bdfac*f;
	}
      } 
    }
    xind+=soc_dim(i);
  }

  //---semidefinite part
  {for(unsigned int i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,qp_dx.get_store()+xind);
    sveci(tmpsvec,dX);
    Real dXtol=norm2(tmpsvec)*1e-10;
    tmpsymmat.init(dX,al); tmpsymmat+=Xp[i];
    while((tmpsymmat.Chol_factor(dXtol*al))&&(al>eps_Real)) {
      al*=.8;
      tmpsymmat.init(dX,al); tmpsymmat+=Xp[i];
    }
    tmpsvec.init(xvdim,1,qp_dz.get_store()+xind);
    sveci(tmpsvec,dX);
    Real dZtol=norm2(tmpsvec)*1e-10;
    tmpsymmat.init(dX,al); tmpsymmat+=Zp[i];
    while((tmpsymmat.Chol_factor(dZtol*al))&&(al>eps_Real)){
      al*=.8;
      tmpsymmat.init(dX,al); tmpsymmat+=Zp[i];
    }
    xind+=xvdim;  
  }}

  if (al<alpha) alpha=al;

  return 0;
}




// ****************************************************************************
//                           starting_x
// ****************************************************************************

//generate a strictly feasible primal starting point
//store it in the qpx_range of x
//returns 0 on success, 1 on failure

int QP_SDPBlock::starting_x(Matrix& qp_x)
{
  Real initval=b/mu_dim;
  if (less_or_equal) {
    s=initval;
  }
  else {
    s=0.;
  }
  x.init(A,initval);
  mat_xey(x.dim(),qp_x.get_store()+qp_xstart,x.get_store());

  for(unsigned int i=0;i<Xp.size();i++){
    Xp[i].init(Xp[i].rowdim(),0.);
    Xp[i].shift_diag(initval);
  }   
  return 0;
}

// ****************************************************************************
//                           starting_yz
// ****************************************************************************

//generate a strictly feasible dual starting point
//store it in the qpy_range of y and in the qpx_range of z
//x is fixed already by a previous call to starting_x and Qx=Q*x
//returns 0 on success, 1 on failure

int QP_SDPBlock::starting_yz(Matrix& qp_y,Matrix& qp_z,
		   const Matrix& qp_Qx, const Matrix& qp_c)
{
  //z = (qp_Qx - qp_c)[qp_xstart .. qp_xstart+x.dim()-1]
  z.init(x.dim(),1,qp_Qx.get_store()+qp_xstart); 
  mat_xmey(z.dim(),z.get_store(),qp_c.get_store()+qp_xstart);

  y=1.;
  {for(Integer i=0;i<lin_dim;i++){
    y=max(y,fabs(z(i))+1.);
  }}
  Integer zind=lin_dim;
  {for(Integer i=0;i<soc_dim.dim();i++){
    const Real *zp=z.get_store()+zind;
    Real a=(*zp++);
    Real b=0.;
    for(Integer j=soc_dim(i)-1;--j>=0;zp++){
      b+=(*zp)*(*zp);
    }
    y=max(y,sqrt(b)-a+1.);
    zind+=soc_dim(i);
  }}

  {for (unsigned int i=0;i<Zp.size();i++){
    Integer Zrdim=Xp[i].rowdim();
    Integer zvdim=(Zrdim*(Zrdim+1))/2;
    tmpsvec.init(zvdim,1,z.get_store()+zind);
    zind+=zvdim;
    sveci(tmpsvec,Zp[i]);
    y=max(y,max(sumrows(abs(Zp[i])))+1.);
  }}
  z.xpeya(A,y);
  {for (unsigned int i=0;i<Zp.size();i++){
    Zp[i].shift_diag(y);
  }}
  
  mat_xey(z.dim(),qp_z.get_store()+qp_xstart,z.get_store());
  qp_y(qp_ystart)=y;

  return 0;
}


// ****************************************************************************
//                            get_AB
// ****************************************************************************
  
//store the local coefficients of matrices A and b in the positions
//corresponding to qpy_range (rows) and qpx_range (columns) 
//returns 0 on success, 1 on failure

int QP_SDPBlock::get_Ab(Matrix& qp_A,Matrix &qp_b) const
{

  mat_xey(A.dim(),qp_A.get_store()+qp_A.rowdim()*qp_xstart+qp_ystart,qp_A.rowdim(),A.get_store(),1);
  qp_b(qp_ystart)=b;

  return 0;
}

// ****************************************************************************
//                           restart_x
// ****************************************************************************
  
//it is assumed that the problem was solved already once and is now
//resolved for a new linear cost term qp_c that resulted from the old
//one by adding qp_dc.
//on input qp_x holds the old optimal solution and on output
//the coorespoind qpx_range should be replaced by a reasonable 
//strictly feasible solution for x suitable for restarting
//(see also restart_yz)
//returns 0 on success, 1 on failure

int QP_SDPBlock::restart_x(Matrix& qp_x,const Matrix& qp_c,const Matrix& qp_dc)
{
  Real normc=mat_ip(x.dim(),qp_c.get_store()+qp_xstart,qp_c.get_store()+qp_xstart);
  Real normdc=mat_ip(x.dim(),qp_dc.get_store()+qp_xstart,qp_dc.get_store()+qp_xstart);
  restart_factor=min(max(0.9,sqrt(1-min(1.,sqrt(normdc/normc)))),0.99999);
  xbpeya(x,A,(b-b*restart_factor)/mu_dim,restart_factor);
  if (less_or_equal){
    s=s*restart_factor+(b-b*restart_factor)/mu_dim;
  }
  mat_xey(x.dim(),qp_x.get_store()+qp_xstart,x.get_store());

  Integer xind=lin_dim;
  for (unsigned int i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,x.get_store()+xind);
    xind+=xvdim;
    sveci(tmpsvec,Xp[i]);
  }
  
  return 0;
}


// ****************************************************************************
//                           restart_yz
// ****************************************************************************
  
//this is called after restart_x (see there)
//on input qp_y and qp_z hold the old optimal solution and on output
//the cooresponding qpy/qpx_range should be replaced by a reasonable 
//strictly feasible solution for y/z suitable for restarting
//returns 0 on success, 1 on failure

int QP_SDPBlock::restart_yz(Matrix& qp_y,Matrix& qp_z,
			    const Matrix& qp_Qx,const Matrix& qp_c,const Matrix& /* qp_dc */)
{
  // tmpvec= A*y-z + qp_Qx - c
  xeyapzb(tmpvec,A,z,y,-1.);
  mat_xpey(tmpvec.dim(),tmpvec.get_store(),qp_Qx.get_store()+qp_xstart);
  mat_xmey(tmpvec.dim(),tmpvec.get_store(),qp_c.get_store()+qp_xstart);

  // make tmpvec positive definit  
  Real dy=0.;
  {for(Integer i=0;i<lin_dim;i++){
    dy=max(dy,fabs(tmpvec(i)));
  }}
  Integer zind=lin_dim;
  {for(Integer i=0;i<soc_dim.dim();i++){
    const Real *zp=tmpvec.get_store()+zind;
    Real a=(*zp++);
    Real b=0.;
    for(Integer j=soc_dim(i)-1;--j>=0;zp++){
      b+=(*zp)*(*zp);
    }
    dy=max(dy,sqrt(b)+fabs(a));
    zind+=soc_dim(i);
  }}
  {for (unsigned int i=0;i<Zp.size();i++){
    Integer Zrdim=Zp[i].rowdim();
    Integer zvdim=(Zrdim*(Zrdim+1))/2;
    tmpsvec.init(zvdim,1,tmpvec.get_store()+zind);
    zind+=zvdim;
    sveci(tmpsvec,Xinv[i]);
    dy=max(dy,max(sumrows(abs(Xinv[i]))));
  }}
  dy+=(1.-restart_factor)*.1;
  dy=max(dy,1e-5);
  z+=tmpvec;
  z.xpeya(A,dy);
  y+=dy;
  {for (unsigned int i=0;i<Zp.size();i++){
    Zp[i].shift_diag(dy);
    Zp[i]+=Xinv[i];
  }}
  
  mat_xey(z.dim(),qp_z.get_store()+qp_xstart,z.get_store());
  qp_y(qp_ystart)=y;
  
  return 0;
}
    
// ****************************************************************************
//                        add_xinv_kron_z
// ****************************************************************************
  
//add the system term corresponding to (xinv kron z)
//(that arises from solving the linearized perturbed complementarity system
// x*z =0 or =mu*I for dx in the preferred search direction)
//to the diagonal block corresponding to qpx_range x qpx_range

int QP_SDPBlock::add_xinv_kron_z(Symmatrix& barQ)
{
  Integer dind=qp_xstart;  //index into the diagonal of barQ
  {for(Integer i=0;i<lin_dim;i++){
    barQ(dind,dind)+=z(i)/x(i);
    dind++;
  }}
  for (Integer i=0;i<soc_dim.dim();i++){
    const Real* xp=x.get_store()+dind-qp_xstart;
    const Real* zp=z.get_store()+dind-qp_xstart;
    Real xs=xp[0]*xp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      xs-=xp[j]*xp[j];
    }
    Real xzip=mat_ip(soc_dim(i),xp,zp);
  
    barQ(dind,dind)+=(-xzip+2*xp[0]*zp[0])/xs;
    {for(Integer j=1;j<soc_dim(i);j++){
      barQ(dind,dind+j)+=(-zp[0]*xp[j]+xp[0]*zp[j])/xs;
      barQ(dind+j,dind+j)+=(xzip-2*zp[j]*xp[j])/xs;
      for(Integer k=j+1;k<soc_dim(i);k++){
	barQ(dind+j,dind+k)+=(-xp[j]*zp[k]-xp[k]*zp[j])/xs;
      }
    }}
    dind+=soc_dim(i);
  }
  {for (unsigned int i=0;i<Xp.size();i++){
    tmpsymmat=Xp[i];
    if (tmpsymmat.Chol_factor(1e-20)){
      if (out) (*out)<<"*** WARNING: QP_SDPBlock::add_xinv_kron_z: factorizing Xp["<<i<<"] failed"<<std::endl;
      return 1;
    }
    tmpsymmat.Chol_inverse(Xinv[i]);
    skron(Xinv[i],Zp[i],tmpsymmat);
    for (Integer j=0;j<tmpsymmat.rowdim();j++){
      mat_xpey(tmpsymmat.rowdim()-j,
	       barQ.get_store()+barQ.rowdim()*dind-(dind*(dind-1))/2,
	       tmpsymmat.get_store()+tmpsymmat.rowdim()*j-(j*(j-1))/2);
      dind++;
    }
  }}
  
  return 0;  
}

// ****************************************************************************
//                           add_local_sys
// ****************************************************************************
  
//on input: sysdy= A*barQ^{-1}*A^T    (barQ as returned in add_xinv_kron_z)
//          rhs= A*barQ^{-1}*(c-Q*x-A^T*y)-(b-A*x)
//if the block uses additional internal variables 
//(like an additional term + B*s with s>=0 in the primal feasibility constr)
//then the corresponding block terms have now to be added to sysdy and rhs,
//eg,
//   sysdy +=  B*(t^{-1} kron s)*B^T     (if t is the dual variable to s)
//   rhs   +=  B*s - B*(t^{-1} kron s)*B^T*y

int QP_SDPBlock::add_local_sys(Symmatrix& sysdy,Matrix& /* rhs */)
{
  if (less_or_equal){
    if (s>=1e10*y) {
      if ((out)&&(print_level>=1)){
	(*out)<<"*** WARNING: QP_SDPBlock::add_local_sys adds large diagonal term: s="<<s<<" y="<<y<<std::endl;
      }
      /*
      y=0;
      for(Integer i=0;i<sysdy.rowdim();i++){
	sysdy(i,qp_ystart)=0.;
      }
      sysdy(qp_ystart,qp_ystart)=1.;
      rhs(qp_ystart)=0;
      return 0;
      */
    }
    sysdy(qp_ystart,qp_ystart)+=s/y;
  }
  return 0;
}

// ****************************************************************************
//                           suggest_mu
// ****************************************************************************
  
//dx, dy, dz is the predictor direction. Based on the predictor 
//(x+dx,y+dy,z+dz) suggest a value for mu by specifying the
//inner product of the dual cone varaibles ip_xz=ip(x,z)+ip(s,t),
//the dimension of the conic variable space mu_dim= cone_x.dim+cone_s.dim
//a value for the factor on mu to obtain the new target

int QP_SDPBlock::suggest_mu(Real& ip_xz,Integer& out_mu_dim,Real& sigma,
	         const Matrix& qp_dx,const Matrix& qp_dy,const Matrix& qp_dz)
{
  out_mu_dim=mu_dim;
  Real ds=0;
  ip_xz=ip(x,z);
  if (less_or_equal){
    if (y>0.) 
      ds=-s*(1+qp_dy(qp_ystart)/y);
    ip_xz+=s*y;
  }
  Real al=1.;
  inner_line_search(al,qp_dx,qp_dy,qp_dz,ds);
  Real ipdxdz=ip_xz;
  if (less_or_equal) {
    ipdxdz+=al*ds*qp_dy(qp_ystart)+al*s*qp_dy(qp_ystart)+al*al*ds*y;
  }
  ipdxdz+=al*mat_ip(x.dim(),x.get_store(),qp_dz.get_store()+qp_xstart);
  ipdxdz+=al*mat_ip(x.dim(),z.get_store(),qp_dx.get_store()+qp_xstart);
  ipdxdz+=al*al*mat_ip(x.dim(),qp_dx.get_store()+qp_xstart,qp_dz.get_store()+qp_xstart);
  sigma=min(1.,pow(ipdxdz/ip_xz,max(2.,3*sqr(al))));
  return 0;
}

// ****************************************************************************
//                           get_corr
// ****************************************************************************
  
//on input (w.r.t. corresponding positions)
//    xcorr = 0
//    rhs as on output of add_local_sys
//on output the corresponding positions of xcorr should hold the corrector
//term of the search direction, eg,  xcorr = mu*x^{-1} - x^{-1}*dx*dz,
//and if the block holds additional local variables as in add_local_sys then
//   rhs += B*(mu * t^{-1}- t^{-1}*dt*ds) 

int QP_SDPBlock::get_corr(Matrix& xcorr,Matrix& rhs,Real mu,const Matrix& qp_dx,
		       const Matrix& qp_dy,const Matrix& qp_dz)
{
  Integer xind=qp_xstart;  
  {for(Integer i=0;i<lin_dim;i++){
    xcorr(i+xind)=(mu-qp_dx(i+xind)*qp_dz(i+xind))/x(i);
  }}
  xind+=lin_dim;

  //------- compute corrector rhs for SOC

  for (Integer i=0;i<soc_dim.dim();i++){
    const Real* xp=x.get_store()+xind-qp_xstart;
    Real x0=xp[0];
    Real xbarsq=mat_ip(soc_dim(i)-1,xp+1,xp+1);
    Real gamma=sqrt(x0*x0-xbarsq);
    Real d=1./(gamma*(xbarsq/(gamma+x0)-x0));
    xcorr(xind)=mu*(-1/gamma-(gamma+x0)*d);
    mat_xeya(soc_dim(i)-1,xcorr.get_store()+xind+1,xp+1,mu*d);
    Matrix socx(soc_dim(i),1,xp);
    Matrix dx(soc_dim(i),1,qp_dx.get_store()+xind);
    Matrix dz(soc_dim(i),1,qp_dz.get_store()+xind);
    Matrix res=mult_Ginv(socx,gamma,mult_Arw(mult_G(socx,gamma,dz),mult_Ginv(socx,gamma,dx)));
    mat_xpeya(soc_dim(i),xcorr.get_store()+xind,res.get_store(),-1.);
    xind+=soc_dim(i);
  }

  //------- compute corrector rhs for SDP

  {for (unsigned int i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,qp_dx.get_store()+xind);  
    sveci(tmpsvec,tmpsymmat);
    genmult(Xinv[i],Matrix(tmpsymmat),tmpvec,-1.);
    tmpsvec.init(xvdim,1,qp_dz.get_store()+xind); 
    sveci(tmpsvec,tmpsymmat);
    genmult(tmpvec,tmpsymmat,tmpsvec);
    tmpsymmat.init(tmpsvec);
    tmpsymmat.xpeya(Xinv[i],mu);
    svec(tmpsymmat,tmpsvec);
    mat_xey(xvdim,xcorr.get_store()+xind,tmpsvec.get_store());
    xind+=xvdim;
  }}
  if (less_or_equal){
    //ds=-s(1+dy/y); (mu-dy*ds)/y=(mu+dy*s*(1+dy/y))/y
    if (y>0.)
      rhs(qp_ystart)+= (mu+qp_dy(qp_ystart)*s*(1+qp_dy(qp_ystart)/y))/y;
    else 
      rhs(qp_ystart)=0.;
  }

  return 0;  
  
}

// ****************************************************************************
//                            line_search
// ****************************************************************************
  
//dx,dy,dz gives the final step direction, alphap and alphad are on input
//upper bounds on the primal and dual step size. 
//On output alpahp and alphad have to be no larger than on input and
//have to guarantee strict feasibility of the primal/dual step on
//the local variables
//(if the block has additional internal variables it has to compute the step
// direction for these now and to choose alphap and alphad so that
// strict feasibility is guaranteed for the internal variables as well)

int QP_SDPBlock::line_search(Real& alpha,const Matrix& qp_dx,
                          const Matrix& qp_dy,const Matrix& qp_dz)
{
  Real ds=0.;
  if (less_or_equal){
    ds=b-s-ip(A,x)-mat_ip(A.dim(),A.get_store(),qp_dx.get_store()+qp_xstart);
  }
  inner_line_search(alpha,qp_dx,qp_dy,qp_dz,ds);

  //in case of doubt or for testing uncomment this
  /*
  for (Integer i=0;i<lin_dim;i++){
    if ((x(i)+alpha*qp_dx(qp_xstart+i)<=0.)||(z(i)+alpha*qp_dz(qp_xstart+i)<=0.)){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::line_search(): liner variables not >=0 for i="<<i<<std::endl;
    }
  }

  Integer xind=lin_dim;
  for (Integer i=0;i<soc_dim.dim();i++){
    Matrix xp(soc_dim(i),1,x.get_store()+xind);
    mat_xpeya(soc_dim(i),xp.get_store(),qp_dx.get_store()+qp_xstart+xind,alpha);
    Real xs=xp[0]*xp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      xs-=xp[j]*xp[j];
    }
    if (xs<0.){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::line_search(): soc_x not soc for i="<<i<<std::endl;
    }
    Matrix zp(soc_dim(i),1,z.get_store()+xind);
    mat_xpeya(soc_dim(i),zp.get_store(),qp_dz.get_store()+qp_xstart+xind,alpha);
    Real zs=zp[0]*zp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      zs-=zp[j]*zp[j];
    }
    if (zs<0.){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::line_search(): soc_z not soc for i="<<i<<std::endl;
    }
    xind+=soc_dim(i);
  }
  for (Integer i=0;i<Xp.size();i++){
    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,x.get_store()+xind);
    mat_xpeya(xvdim,tmpsvec.get_store(),qp_dx.get_store()+qp_xstart+xind,alpha);
    sveci(tmpsvec,tmpsymmat);
    if (tmpsymmat.Chol_factor()){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::line_search(): X not psd for i="<<i<<std::endl;
    }
    tmpsvec.init(xvdim,1,z.get_store()+xind);
    mat_xpeya(xvdim,tmpsvec.get_store(),qp_dz.get_store()+qp_xstart+xind,alpha);
    sveci(tmpsvec,tmpsymmat);
    if (tmpsymmat.Chol_factor()){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::line_search(): X not psd for i="<<i<<std::endl;
    }
    xind+=xvdim;
  }
  */

  return 0;
}

// ****************************************************************************
//                           set_point
// ****************************************************************************
  
  //x,y,z is the new point and has to be stored
  //alphap and alphad are the step size used in the step, they are passed
  //so that the block can take the same step for internal variables if
  //needed.

int QP_SDPBlock::set_point(
		const Matrix& qp_x,const Matrix& qp_y,const Matrix& qp_z,
		Real /* alpha */)
{
  old_x=x;
  old_z=z;
  old_y=y;
  old_s=s;
  
  Integer xind=qp_xstart;
  mat_xey(x.dim(),x.get_store(),qp_x.get_store()+xind);
  mat_xey(z.dim(),z.get_store(),qp_z.get_store()+xind);
  y=qp_y(qp_ystart);
  if (less_or_equal){
    s=b-ip(A,x);
  }

  xind=lin_dim+sum(soc_dim);

  for (unsigned int i=0;i<Xp.size();i++){
    old_Xp[i]=Xp[i];
    old_Zp[i]=Zp[i];

    Integer Xrdim=Xp[i].rowdim();
    Integer xvdim=(Xrdim*(Xrdim+1))/2;
    tmpsvec.init(xvdim,1,x.get_store()+xind);
    sveci(tmpsvec,Xp[i]);
    tmpsvec.init(xvdim,1,z.get_store()+xind);
    sveci(tmpsvec,Zp[i]);
    xind+=xvdim;
  }

  //in case of doubt or for testing uncomment this
  /*
  for (Integer i=0;i<lin_dim;i++){
    if ((x(i)<=0.)||(z(i)<=0.)){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::set_point(): liner variables not >=0 for i="<<i<<std::endl;
    }
  }
  xind=lin_dim;
  for (Integer i=0;i<soc_dim.dim();i++){
    const Real* xp=x.get_store()+xind;
    Real xs=xp[0]*xp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      xs-=xp[j]*xp[j];
    }
    if (xs<0.){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::set_point(): soc_x not soc for i="<<i<<std::endl;
    }
    const Real* zp=z.get_store()+xind;
    Real zs=zp[0]*zp[0];
    for(Integer j=1;j<soc_dim(i);j++){
      zs-=zp[j]*zp[j];
    }
    if (zs<0.){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::set_point(): soc_z not soc for i="<<i<<std::endl;
    }
    xind+=soc_dim(i);
  }
  for (Integer i=0;i<Xp.size();i++){
    tmpsymmat=Xp[i];
    if (tmpsymmat.Chol_factor()){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::set_point(): X not psd for i="<<i<<std::endl;
    }
    tmpsymmat=Zp[i];
    if (tmpsymmat.Chol_factor()){
      if (out) (*out)<<"*** ERROR: QP_SDPBlock::set_point(): X not psd for i="<<i<<std::endl;
    }
  }
  */

  return 0;
}


int QP_SDPBlock::add_Bs(Matrix& qp_vec) const
{
  if (less_or_equal){
    qp_vec(qp_ystart)+=s;
  }
  return 0;
}




}

