/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/qp_solver.cxx

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



#include <iomanip>
#include "qp_solver.hxx"


using namespace CH_Matrix_Classes;

namespace ConicBundle {


void QP_Solver::set_defaults()
{
  termeps=1e-7;
  maxiter=-1;
}

void QP_Solver::clear()
{
  Q.init(0,0.);
  c.init(0,0,0.); 
  offset=0.;  
  blockp.clear();
  A.init(0,0,0.);
  b.init(0,0,0.);
  lowerbound=0;
  upperbound=0;
  x.init(0,0,0.);
  y.init(0,0,0.);
  z.init(0,0,0.);        
  primalval=0.;      
  dualval=0.;        
  mu=0.;             
  Qx.init(0,0,0.);           
  iter=0;        
  status=0;      
  dx.init(0,0,0.);
  dy.init(0,0,0.);
  dz.init(0,0,0.);        
  LDL_Q.init(0,0.); 
  LinvAt.init(0,0,0.);   
  sysdy.init(0,0.); 
  rd.init(0,0,0.);       
  rhs.init(0,0,0.);      
  xcorr.init(0,0,0.);
  tmpvec.init(0,0,0.);   

}


// *************************************************************************
//                            iterate
// *************************************************************************

// loop till convergence to optimal solution

int QP_Solver::iterate()
{

  //cout<<" norm2(A*x-b)="<<norm2(A*x-b);          //TEST
  //cout<<" norm2(Q*x-c+transpose(A)*y-z)="<<norm2(Q*x-c+transpose(A)*y-z)<<std::endl; //TEST
  //cout<<"Q=";Q.display(std::cout);
  //cout<<"c=";c.display(std::cout);
  //cout<<"A=";A.display(std::cout);
  //cout<<"b=";b.display(std::cout);
  //cout<<"x=";x.display(std::cout);
  //cout<<"y=";y.display(std::cout);
  //cout<<"z=";z.display(std::cout);

 // display parameters
 if ((out)&&(print_level>=1)){
     *out<<"    mi="<<maxiter;
     *out<<" te=";out->precision(5);*out<<termeps;
     *out<<" lb=";out->precision(10);out->width(11);*out<<lowerbound;
     *out<<" ub=";out->precision(10);out->width(11);*out<<upperbound;out->precision(2);
     *out<<" Qn="<<Q.rowdim();
     *out<<std::endl;
 }

 // compute primal and dual objective value
 // (x,y,z,and Qx are already computed)

 Real xtQxo2=ip(x,Qx)/2.;
 //x.transpose(); std::cout<<"x= ";x.display(std::cout); x.transpose(); //TEST
 //cout<<"A= ";A.display(std::cout); //TEST
 //cout<<"b= ";b.display(std::cout); //TEST
 //cout<<" xtQxo2="<<xtQxo2; //TEST
 //cout<<" ip(c,x)="<<ip(c,x); //TEST
 //cout<<" offset="<<offset<<std::endl; //TEST
 //Matrix AxBs=A*x; //TEST
 //for(Integer i=0;i<blockp.size();i++) blockp[i]->add_Bs(AxBs); //TEST
 //cout<<" Fx=";out->width(7);*out<<norm2(b-AxBs); //TEST
 //cout<<" Fz=";out->width(7);*out<<norm2(c+z-Qx-transpose(A)*y); //TEST

 dualval=xtQxo2+ip(b,y)+offset;
 primalval=-xtQxo2+ip(c,x)+offset;
 Real ub=upperbound;
 Real lb=lowerbound;
 if (lb<=primalval) lb=primalval;
 
 //output
 if ((out)&&(print_level>=1)){
   out->precision(2);
   *out<<"    ";out->width(2);*out<<iter<<":";
   *out<<" gappd=";out->width(7);*out<<dualval-primalval;
   *out<<" pv=";out->precision(8);out->width(10);*out<<primalval;
   *out<<" dv=";out->precision(8);out->width(10);*out<<dualval;
   out->precision(2);
   *out<<" epsdl=";out->width(7);*out<<.5*(dualval-lb);
   *out<<" epsup=";out->width(7);*out<<termeps*(ub-primalval)<<std::endl;
 }     

 Real alpha=1.;
 Indexmatrix piv;

 while(((maxiter<0)||(iter<maxiter))&&
       (dualval-lb >1e-12*(fabs(dualval)+1.))&&
       (ub-primalval>1e-12*(fabs(primalval)+1.))&&
       (alpha>1e-8)&&
       (dualval-primalval>min(.5*(dualval-lb),termeps*(ub-primalval)))
       ){

  iter++;
   
   //--- build Q-part of system matrix
   LDL_Q=Q;
   for(unsigned int i=0;i<blockp.size();i++) {
     if (blockp[i]->add_xinv_kron_z(LDL_Q)){
       if ((out)&&(print_level>=1)){
	 (*out)<<"*** WARNING: QP_Solver::iterate(): collecting system matrix failed for subprobelm "<<i<<std::endl;
       }
       status=3;
       return status;
     }
   }

   Symmatrix barQ(LDL_Q);                   //TEST
   
   //--- factorize Q-part of system matrix and compute system for y
   if (LDL_Q.Chol_factor(eps_Real)){
     if ((out)&&(print_level>=1)){
       (*out)<<"*** WARNING: QP_Solver::iterate(): factorizing system matrix failed"<<std::endl;
     }
     //cout<<barQ<<std::endl;
     //cout<<Q<<std::endl;
     status=2;
     return status;
   }
   
   //Symmatrix barQinv;                        //TEST
   //Matrix barQinvmat(Diag(Matrix(LDL_Q.dim(),1,1.)));//TEST
   //LDL_Q.Chol_solve(barQinvmat);//TEST
   //barQinv=barQinvmat;//TEST
   //cout<<"  norm2(barQinv*barQ-Diag(Matrix(LDL_Q.dim(),1,1.)))="<<norm2(barQinv*barQ-Diag(Matrix(LDL_Q.dim(),1,1.)))<<std::endl;//TEST

   //LinvAt=L^-1*transpose(A)
   xbpeya(LinvAt,A,1.,0.,1);       
   LDL_Q.Chol_Lsolve(LinvAt);

   //sysdy=LiAt'*LiAt
   rankadd(LinvAt,sysdy,1.,0.,1);

   //cout<<"  norm2(sysdy-A*barQinv*transpose(A))="<<norm2(sysdy-A*barQinv*transpose(A))<<std::endl;//TEST
   
   //rd=c-Qx-At*y;  
   xeyapzb(rd,c,Qx,1.,-1.);
   if (y.dim()>0) {
     genmult(A,y,rd,-1.,1.,1);    

     //cout<<"  norm2(rd-c+Qx+transpose(A)*y)="<<norm2(rd-c+Qx+transpose(A)*y)<<std::endl;//TEST

     //rhs= A*barQinv*rd-(b-Ax)
     tmpvec=rd;    
     LDL_Q.Chol_Lsolve(tmpvec);
     genmult(LinvAt,tmpvec,rhs,1.,0.,1);
     rhs-=b;
     genmult(A,x,rhs,1.,1.);

     //cout<<"  norm2(rhs-A*barQinv*rd+b-A*x)="<<norm2(rhs-A*barQinv*rd+b-A*x)<<std::endl;//TEST
 
   
     //--- build dy-part of system matrix
     for(unsigned int i=0;i<blockp.size();i++) {
       if (blockp[i]->add_local_sys(sysdy,rhs)){
	 if ((out)&&(print_level>=1)){
	   (*out)<<"*** WARNING: QP_Solver::iterate(): adding local system "<<i<<" for predictor failed"<<std::endl;
	 }
	 status=3;
	 return status;
       }
     }

     //Symmatrix testsysdy=sysdy; //TEST

   
     //--- solve predictor direction
     dy=rhs;
     sysdy.Chol_factor(piv,1000.*eps_Real);
     if(piv.dim()<sysdy.rowdim()){
       if ((out)&&(print_level>=1)){
	 (*out)<<"*** WARNING: QP_Solver::iterate(): factorizing for predictor failed"<<std::endl;
       }
       status=2;
       return status;
     }
     sysdy.Chol_solve(dy,piv);

     //cout<<"  norm2(testsysdy*dy-rhs)="<<norm2(testsysdy*dy-rhs)<<std::endl; //TEST

     //dx=barQinv*(c-At*y-Qx-At*dy)
     //dz=-(c-At*y-Qx-At*dy)-z+Q*dx;
     genmult(A,dy,dx,-1.,0.,1);
     dx+=rd;
   }
   else {
     dy.init(0,0,0.);
     rhs.init(0,0,0.);
     dx=rd;
   }
   xeyapzb(dz,dx,z,-1.,-1.);
   LDL_Q.Chol_solve(dx);
   genmult(Q,dx,dz,1.,1.);

   //cout<<" norm2(barQ*dx-c+Qx+transpose(A)*(y+dy))="<<norm2(barQ*dx-c+Qx+transpose(A)*(y+dy))<<std::endl; //TEST
   //cout<<" norm2(Q*(x+dx)-c+transpose(A)*(y+dy)-z-dz)="<<norm2(Q*(x+dx)-c+transpose(A)*(y+dy)-z-dz)<<std::endl; //TEST
   //Symmatrix Xtest,Ztest,dX1test,dZ1test;//TEST 
   //Integer dim1=Integer(::sqrt(Real(8*x.dim()+1))-1+.1)/2;
   //if (x.dim()!=(dim1*(dim1+1))/2){//TEST
   //  sveci(x(Range(1,x.dim()-1)),Xtest);//TEST
   //  sveci(z(Range(1,z.dim()-1)),Ztest);//TEST
   //  sveci(dx(Range(1,x.dim()-1)),dX1test);//TEST
   //  sveci(dz(Range(1,x.dim()-1)),dZ1test);//TEST
   //}                //TEST
   //else {           //TEST
   //  sveci(x,Xtest);//TEST
   //  sveci(z,Ztest);//TEST
   //  sveci(dx,dX1test);//TEST
   //  sveci(dz,dZ1test);//TEST
   //}                  //TEST
   //cout<<" norm2(Symmatrix(Xtest*Ztest*dX1test)+Xtest*(Ztest+dZ1test)*Xtest)="<<norm2(Symmatrix(Xtest*Ztest*dX1test)+Xtest*(Ztest+dZ1test)*Xtest)<<std::endl; //TEST
   
   //--- select new mu
   Real oldmu=mu;
   Integer sum_mu_dim=0;
   Integer mu_dim;
   Real ip_xz;
   mu=0.;
   Real maxsigma=0.01;
   Real sigma;
   {for(unsigned int i=0;i<blockp.size();i++) {
     blockp[i]->suggest_mu(ip_xz,mu_dim,sigma,dx,dy,dz);
     mu+=ip_xz;
     sum_mu_dim+=mu_dim;
     maxsigma=max(maxsigma,sigma);
   }}
   mu/=sum_mu_dim;
   maxsigma=min(1.,maxsigma);
   mu=min(oldmu,maxsigma*mu);

   //--- compute corrector terms (and modify rhs for iternal correctors)
   xcorr.init(x.dim(),1,0.);
   {for(unsigned int i=0;i<blockp.size();i++) {
     if (blockp[i]->get_corr(xcorr,rhs,mu,dx,dy,dz)){
       if ((out)&&(print_level>=1)){
	 (*out)<<"*** WARNING: QP_Solver::iterate(): getting local correctors failed"<<std::endl;
       }
       status=3;
       return status;
     }
   }}
   
   //--- solve corrector direction
   // dy = sysinv*( barQinv*xcorr+rhs) 
   if (y.dim()>0){
     tmpvec=xcorr;    
     LDL_Q.Chol_Lsolve(tmpvec);
     genmult(LinvAt,tmpvec,dy,1.,0.,1);
     dy+=rhs;
     sysdy.Chol_solve(dy,piv);
     //dx=barQinv*(c-At*y-Qx-At*dy+xcorr)
     //dz=-(c-At*y-Qx-At*dy)-z+Q*dx;
     genmult(A,dy,dx,-1.,0.,1);
     dx+=rd;
   }
   else {
     dx=rd;
   }
   xeyapzb(dz,dx,z,-1.,-1.);
   dx+=xcorr;
   LDL_Q.Chol_solve(dx);
   genmult(Q,dx,dz,1.,1.);

   //cout<<" norm2(barQ*dx-c+Qx+transpose(A)*(y+dy)-xcorr)="<<norm2(barQ*dx-c+Qx+transpose(A)*(y+dy)-xcorr)<<std::endl; //TEST
   //cout<<" norm2(Q*(x+dx)-c+transpose(A)*(y+dy)-z-dz)="<<norm2(Q*(x+dx)-c+transpose(A)*(y+dy)-z-dz)<<std::endl; //TEST
   //Symmatrix dX2test,dZ2test;//TEST 
   //Integer dim2=Integer(::sqrt(Real(8*dx.dim()+1))-1+.1)/2;
   //if (x.dim()!=(dim2*(dim2+1))/2){//TEST
   //  sveci(dx(Range(1,x.dim()-1)),dX2test);//TEST
   //  sveci(dz(Range(1,x.dim()-1)),dZ2test);//TEST
   //}                //TEST
   //else {           //TEST
   //  sveci(dx,dX2test);//TEST
   //  sveci(dz,dZ2test);//TEST
   //}                  //TEST
   //cout<<"  norm2(Symmatrix(Xtest*Ztest*dX2test+dX1test*dZ1test*Xtest)+Xtest*(Ztest+dZ2test)*Xtest-mu*Xtest)="<<norm2(Symmatrix(Xtest*Ztest*dX2test+dX1test*dZ1test*Xtest)+Xtest*(Ztest+dZ2test)*Xtest-mu*Xtest)<<std::endl; //TEST

   //--- line search
   alpha=1.;
   {for(unsigned int i=0;i<blockp.size();i++) {
     if (blockp[i]->line_search(alpha,dx,dy,dz)){
       if ((out)&&(print_level>=1)){
	 (*out)<<"*** WARNING: QP_Solver::iterate(): line search failed"<<std::endl;
       }
       status=3;
       return status;
     }
   }}
   
   //--- x+=alpha*dx, y+=alpha*dy, z+=alpha*dz, set new point
   x.xpeya(dx,alpha);
   y.xpeya(dy,alpha);
   z.xpeya(dz,alpha);
   
   {for(unsigned int i=0;i<blockp.size();i++) {
     if (blockp[i]->set_point(x,y,z,alpha)){
       if ((out)&&(print_level>=1)){
	 (*out)<<"*** WARNING: QP_Solver::iterate(): setting new point failed"<<std::endl;
       }
       status=3;
       return status;
     }
   }}
   
   genmult(Q,x,Qx);
   Real xtQxo2=ip(x,Qx)/2.;
   
   dualval=xtQxo2+ip(b,y)+offset;
   primalval=-xtQxo2+ip(c,x)+offset;
   
   //output
   if ((out)&&(print_level>=1)){
     out->precision(2);
     //Matrix AxBs=A*x; //TEST
     //for(Integer i=0;i<blockp.size();i++) blockp[i]->add_Bs(AxBs); //TEST
     //cout<<" Fx=";out->width(7);*out<<norm2(b-AxBs); //TEST
     //cout<<" Fz=";out->width(7);*out<<norm2(c+z-Qx-transpose(A)*y); //TEST

     *out<<"    ";out->width(2);*out<<iter<<":";
     *out<<" gappd=";out->width(7);*out<<dualval-primalval;
     *out<<" pv=";out->precision(8);out->width(10);*out<<primalval;
     *out<<" dv=";out->precision(8);out->width(10);*out<<dualval;
     out->precision(2);
     *out<<" epsdl=";out->width(7);*out<<.5*(dualval-lb);
     *out<<" epsup=";out->width(7);*out<<termeps*(ub-primalval);
     *out<<" mu=";out->width(7);*out<<mu;
     *out<<" alpha=";out->width(4);*out<<alpha<<std::endl;

   }     
   
 } //end while

 if ((out)&&(print_level>=1)){
   (*out)<<"    term: iter=";
   if((maxiter>0)&&(iter<maxiter))
     (*out)<<"T";
   else
     (*out)<<"F";
   (*out)<<" dv-lb=";
   if (dualval-lb <=1e-12*(fabs(dualval)+1.))
     (*out)<<"T";
   else
     (*out)<<"F";
   (*out)<<" ub-pv=";
   if(ub-primalval<=1e-12*(fabs(primalval)+1.))
     (*out)<<"T";
   else
     (*out)<<"F";
   (*out)<<" alpha=";
   if (alpha<=1e-8)
     (*out)<<"T";
   else
     (*out)<<"F";
   (*out)<<" dv-pv=";
   if (dualval-primalval<=min(.5*(dualval-lb),termeps*(ub-primalval)))
     (*out)<<"T";
   else
     (*out)<<"F";
   (*out)<<std::endl;
   
 }
 
 if ((maxiter<0)||(iter<maxiter)){
     status=0;
 }
 else {
     status=1;
 }

 return  status;
}

// *************************************************************************
//                             solve
// *************************************************************************

// call starting_point and loop till convergence to optimal solution

int QP_Solver::solve(const Symmatrix& Qin,const Matrix& cin,Real offsetin)
{
  Q=Qin;
  c=cin;
  offset=offsetin;

  //initialize dimensions and primal starting point
  x.init(Q.rowdim(),1,0.);
  Integer ydim=0;
  Integer xdim=0;
  for(unsigned int i=0;i<blockp.size();i++) {
    blockp[i]->set_qp_xstart(xdim);
    blockp[i]->set_qp_ystart(ydim);
    xdim+=blockp[i]->xdim();
    ydim+=blockp[i]->ydim();
    blockp[i]->starting_x(x);
  }

  //initialize dual starting point and constraints
  z.init(xdim,1,0.);
  y.init(ydim,1,0.);
  A.init(ydim,xdim,0.);
  b.init(ydim,1,0.);
  genmult(Q,x,Qx);   //Qx=Q*x
  {for(unsigned int i=0;i<blockp.size();i++) {
    blockp[i]->starting_yz(y,z,Qx,c);
    blockp[i]->get_Ab(A,b);
  }}
    
  //start solving
  mu=1e100;
  iter=0;

  status=iterate();
  
  return status;
}

// *************************************************************************
//                             resolve
// *************************************************************************

//this routine is called if the right hand sides of the primal constraints 
// were changed
// resolve for the same cost function
// call starting_point and loop till convergence to optimal solution

int QP_Solver::resolve()
{

  //initialize dimensions and primal starting point
  x.init(Q.rowdim(),1,0.);
  Integer ydim=0;
  Integer xdim=0;
  for(unsigned int i=0;i<blockp.size();i++) {
    blockp[i]->set_qp_xstart(xdim);
    blockp[i]->set_qp_ystart(ydim);
    xdim+=blockp[i]->xdim();
    ydim+=blockp[i]->ydim();
    blockp[i]->starting_x(x);
  }

  //initialize dual starting point and constraints
  z.init(xdim,1,0.);
  y.init(ydim,1,0.);
  A.init(ydim,xdim,0.);
  b.init(ydim,1,0.);
  genmult(Q,x,Qx);   //Qx=Q*x
  {for(unsigned int i=0;i<blockp.size();i++) {
    blockp[i]->starting_yz(y,z,Qx,c);
    blockp[i]->get_Ab(A,b);
  }}
    
  //start solving
  mu=1e100;
  iter=0;

  status=iterate();
  
  return status;
}

// *************************************************************************
//                             update
// *************************************************************************

// call starting_point and loop till convergence to optimal solution

int QP_Solver::update(const Matrix& dc,Real doffset)
{
  offset+=doffset;
  c+=dc;
  for(unsigned int i=0;i<blockp.size();i++) {
    blockp[i]->restart_x(x,c,dc);
  }
  genmult(Q,x,Qx);   //Qx=Q*x
  {for(unsigned int i=0;i<blockp.size();i++) {
    blockp[i]->restart_yz(y,z,Qx,c,dc);
  }}
  
  mu=1e100;
  iter=0;

  status=iterate();

  return status;
}

// *************************************************************************
//                             save
// *************************************************************************

// write all variables to out in order to enable resuming interrupted sessions

std::ostream& QP_Solver::save(std::ostream& o) const
{
 o.precision(20);
 o<<Q<<c<<offset<<"\n";
 o<<primalval<<"\n";
 o<<dualval<<"\n";
 o<<termeps<<"\n";
 o<<maxiter<<"\n";
 o<<lowerbound<<"\n";
 o<<upperbound<<"\n";
 o<<iter<<"\n";
 o<<status<<"\n";
 return o;
}

// *************************************************************************
//                             restore
// *************************************************************************

// read all variables from "in" in order to enable resuming interrupted sessions

std::istream& QP_Solver::restore(std::istream& in)
{
 in>>Q>>c>>offset;
 in>>primalval;
 in>>dualval;
 in>>termeps;   
 in>>maxiter;   
 in>>lowerbound;
 in>>upperbound;
 in>>iter;
 in>>status;
 return in;
}

}

