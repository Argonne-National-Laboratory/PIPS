/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/nnls.cxx

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



#include <stdlib.h>
#include "symmat.hxx"
#include "heapsort.hxx"
#include "gb_rand.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

// **************************************************************************
//                                 nnls
// **************************************************************************

// least squares solution of
// min ||Ax-b|| s.t. x >=0;
// kkt: A'*A*x - A'*b - l = 0; x >=0, l>=0, x'*l=0
// solved by interior point method with QR-solution of the extended system
// [P. Matsoms, "Sparse Linear Least Squares Problems in Optimization",
//               Comput. Opt. and Appl., 7, 89-110 (1997)] 

int Matrix::nnls(Matrix& rhs,Matrix *dual,Real tol) const
{
 chk_init(*this);
 chk_init(rhs);
 Integer i,j;

 //--- positive starting point 
 Matrix x(nc,1,1./Real(nr));                   //primal variable
 Matrix y=CH_Matrix_Classes::transpose(*this)*((*this)*x-rhs);  //dual variable
 for(i=0;i<nc;i++) 
     if (y(i)<.1)
         y(i)=.1;
 
 Matrix dx(nc,1),dy(nc,1);                   //step directions
 Matrix xsr(nc,1),ysr(nc,1);                 //sqrt(x),sqrt(y);
 Symmatrix AtA; AtA.newsize(nc);
 chk_set_init(AtA,1);
 // AtA=transpose(*this)*(*this);
 for(i=0;i<nc;i++){       
     for(j=i;j<nc;j++){
         AtA(i,j)=mat_ip(nr,m+i*nr,m+j*nr);
     }
 }
 Matrix Atb(nc,1);
 chk_set_init(Atb,1);
 // Atb=::transpose(*this)*rhs;
 for(i=0;i<nc;i++){
     Atb(i)=mat_ip(nr,m+i*nr,rhs.m);
 }
 Symmatrix AtApXiY;                          //system for search direction
 Matrix sysrhs;                              //rhs for search direction
 Matrix xi;
     
 Real mu=ip(x,y)/Real(nc)/5.;                //barrier parameter
 Integer iter=0;
 Integer failed=0;

 //--- repeat till gap is small enough
 while ((ip(x,y)>tol*max(1.,max(max(x),max(y))))&&(iter<50)){
     iter++;

     //--- compute step direction
     AtApXiY.init(AtA);
     // AtApXiY+=diag(y%::inv(x))
     for(i=0;i<nc;i++){
         AtApXiY(i,i)+=y(i)/x(i);
     }
     // sysrhs=::inv(x)*mu-AtA*x+Atb;
     xi.init(x);
     xi.inv();
     sysrhs.init(xi);
     sysrhs*=mu;
     sysrhs-=AtA*x;
     sysrhs+=Atb;
     // factor and solve
     if (AtApXiY.LDLfactor(1e-14)) {
         MEmessage(MatrixError(ME_warning,"Matrix::nnls(): could not factor AtApXiY",MTmatrix));
         failed=1;
         break;
     }
     AtApXiY.LDLsolve(sysrhs);
     dx=sysrhs;
     dy=mu*xi-y-xi%y%dx;
     
     //--- line search
     Real alpha=1.;                //primal stepsize
     for(i=0;i<nc;i++){
         if (dx(i)<-tol/100.){
             alpha=min(alpha,-.99999*x(i)/dx(i));
         }
     }
     Real beta=1.;                 //dual stepsize
     for(i=0;i<nc;i++){
         if (dy(i)<-tol/100.){
             beta=min(beta,-.99999*y(i)/dy(i));
         }
     }

     //--- update
     
     
     x+=alpha*dx;
     y+=beta*dy;
     mu=ip(x,y)/nc/10.;
     if (min(alpha,beta)>.95) mu/=10.;
 }
 if (iter>=50) {
     MEmessage(MatrixError(ME_warning,"Matrix::nnls(): iter>=50",MTmatrix));
     failed=1;
 }
 
 if (dual) *dual=y;
 rhs=x;
 return failed;
}

/* second version with a 2nc x 2nc system

int Matrix::nnls(Matrix& rhs,Matrix *dual,Real tol)
{
 chk_init(*this);
 chk_init(rhs);
 Integer i,j;

 //--- positive starting point 
 Matrix x(nc,1,1./Real(nr));                   //primal variable
 Matrix y=::transpose(*this)*((*this)*x-rhs);  //dual variable
 for(i=0;i<nc;i++) 
     if (y(i)<.1)
         y(i)=.1;
 
 Matrix dx(nc,1),dy(nc,1);                   //step directions
 Matrix xsr(nc,1),ysr(nc,1);                 //sqrt(x),sqrt(y);
 Matrix AtA=::transpose(*this)*(*this);
 Matrix AtAI(AtA);
 AtAI.concat_right(diag(Matrix(nc,1,-1.)));
 Matrix Atb=::transpose(*this)*rhs;
 Matrix extA;                                //system for search direction
 Matrix extrhs;                              //rhs for search direction
     
 Real mu=ip(x,y)/Real(nc)/5.;                //barrier parameter
 Integer iter=0;

 //--- repeat till gap is small enough
 while ((ip(x,y)>tol*max(1.,max(max(x),max(y))))&&(iter<50)){
     iter++;

     //--- compute step direction
     extA=diag(y);
     extA.concat_right(diag(x));
     extA.concat_below(AtAI);
     extrhs=mu-x%y;
     extrhs.concat_below(-AtA*x+y+Atb);
     if (extA.QR_solve(extrhs)) return 1;
     dx=extrhs(Range(0,nc-1));
     dy=extrhs(Range(nc,2*nc-1));
     //dy=::transpose(*this)*((*this)*dx+((*this)*x-rhs))+y;

     //--- line search
     Real alpha=1.;                //primal stepsize
     for(i=0;i<nc;i++){
         if (dx(i)<-tol/100.){
             alpha=min(alpha,-.99999*x(i)/dx(i));
         }
     }
     Real beta=1.;                 //dual stepsize
     for(i=0;i<nc;i++){
         if (dy(i)<-tol/100.){
             beta=min(beta,-.99999*y(i)/dy(i));
         }
     }

     //--- update
     
     
     x+=alpha*dx;
     y+=beta*dy;
     mu=ip(x,y)/nc/10.;
     if (min(alpha,beta)>.95) mu/=10.;
 }
 if (iter>=50) {
     MEmessage(MatrixError(ME_num,"Matrix::nnls(): iter>=50",MTmatrix));
     return 1;
 }
 
 if (dual) *dual=y;
 rhs=x;
 return 0;
}
*/

/*
// for sparse matrices
// [P. Matsoms, "Sparse Linear Least Squares Problems in Optimization",
//               Comput. Opt. and Appl., 7, 89-110 (1997)] 

int Matrix::nnls(Matrix& rhs,Matrix *dual,Real tol)
{
 chk_init(*this);
 chk_init(rhs);
 Integer i,j;

 //--- positive starting point 
 Matrix x(nc,1,1./Real(nr));                   //primal variable
 Matrix y=::transpose(*this)*((*this)*x-rhs);  //dual variable
 for(i=0;i<nc;i++) 
     if (y(i)<.1)
         y(i)=.1;
 
 Matrix dx(nc,1),dy(nc,1);                   //step directions
 Matrix xsr(nc,1),ysr(nc,1);                 //sqrt(x),sqrt(y);
 Matrix extA;                                //system for search direction
 Matrix extrhs;                              //rhs for search direction
     
 Real mu=ip(x,y)/Real(nc)/5.;                //barrier parameter
 Integer iter=0;

 //--- repeat till gap is small enough
 while ((ip(x,y)>tol*max(1.,max(max(x),max(y))))&&(iter<50)){
     iter++;

     //--- compute step direction
     xsr=::sqrt(x); ysr=::sqrt(y);
     xsr.inv();
     extA=(*this);
     extA.concat_below(diag(ysr%xsr));
     extrhs= rhs-(*this)*x;
     ysr.inv();
     extrhs.concat_below(mu*xsr%ysr);
     //extrhs=transpose(extA)*extrhs;
     if (extA.QR_solve(extrhs)) return 1;
     dx=extrhs;
     dy=-y + (mu - dx % y) % ::inv(x);
     //dy=::transpose(*this)*((*this)*dx+((*this)*x-rhs))+y;

     //--- line search
     Real alpha=1.;                //primal stepsize
     for(i=0;i<nc;i++){
         if (dx(i)<-tol/100.){
             alpha=min(alpha,-.99999*x(i)/dx(i));
         }
     }
     Real beta=1.;                 //dual stepsize
     for(i=0;i<nc;i++){
         if (dy(i)<-tol/100.){
             beta=min(beta,-.99999*y(i)/dy(i));
         }
     }

     //--- update
     
     
     
     x+=alpha*dx;
     y+=beta*dy;
     mu=ip(x,y)/nc/10.;
     if (min(alpha,beta)>.95) mu/=10.;
 }
 if (iter>=50) {
     MEmessage(MatrixError(ME_num,"Matrix::nnls(): iter>=50",MTmatrix));
     return 1;
 }
 
 if (dual) *dual=y;
 rhs=x;
 return 0;
}
*/


}

