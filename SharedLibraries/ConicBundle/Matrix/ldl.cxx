/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/ldl.cxx

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

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {


int Symmatrix::LDLfactor(Real tol)
{
  chk_init(*this);

  Real *v;
  memarray->get(nr,v);
  if (v==0) MEmessage(MEmem((nr*(nr+1))/2,"Symmatrix::LDLfactor not enough Memory for v",MTsymmetric));
  Real *hm;
  memarray->get((nr*(nr+1))/2,hm);
  if (hm==0) MEmessage(MEmem((nr*(nr+1))/2,"Symmatrix::LDLfactor not enough Memory to copy",MTsymmetric));

  const Integer hn=nr;
  {
    //copy to reverse representation
    Real *hmp=hm;
    for(Integer i=0;i<hn;i++){
      const Real *mp=m+i;
      const Integer hhn=hn-i-1;
      for(Integer j=hn;--j>=hhn;){
	(*hmp++)=*mp;
	mp+=j;
      }
    }
  }
  //compute LDL-factorization (Golub/van Loan, Algorithm 4.1.2)
  Integer failure=0;
  for (Integer j=0;j<hn;j++){
    const Real* hmp=hm+(j*(j+1))/2;
    Real* mp=hm;
    Real* vp=v;
    Real d=0.;
    for(Integer i=0;i<j;++i){
      const Real d1=*hmp++;         //A(j,i)
      d+=((*vp++)=(*mp)*d1)*d1;
      mp+=i+2;
    }
    d=*vp=*mp-=d;    
    if (d<tol) {
         failure=1;
         break;
    }
    hmp++;  //points now to A(j+1,j)
    {for(Integer i=j+1;i<hn;){
      Real d1=mat_ip(j,hmp,v);
      mp=const_cast<Real *>(hmp+j);
      *mp=(*mp-d1)/d;
      hmp+=++i;
    }}
  }
  {
    //copy to back to normal representation
    const Real *hmp=hm;
    for(Integer i=0;i<hn;i++){
      Real *mp=m+i;
      const Integer hhn=hn-i-1;
      for(Integer j=hn;--j>=hhn;){
	*mp=(*hmp++);
	mp+=j;
      }
    }
  }
  memarray->free(v);
  memarray->free(hm);
  return failure;
}
  

int Symmatrix::LDLsolve(Matrix& x) const
{
 chk_mult(*this,x);

 for(Integer k=0;k<x.coldim();k++){ //solve for and overwrite column k of x
   Real* xbase=x.m+k*nr;
   //---- solve Lr=xbase
   Integer indi=1;
   for(Integer i=1;i<nr;i++){
     mat_xpeya(nr-i,xbase+i,m+indi,-xbase[i-1]);
     indi+=nr-i+1;
   }
   
   //---- solve Dr=r
   {
     const Real* mp=m;
     Real *rp=xbase;
     for(Integer i=nr;i>0;i--) {
       (*rp++)/=*mp;
       mp+=i;
     }
   }
   
   //---- solve L'r=r
   indi=(nr*(nr+1))/2;
   {for(Integer i=nr-1;--i>=0;){ //compute value for row i
     indi-=nr-i;
     xbase[i]-=mat_ip(nr-i-1,xbase+i+1,m+indi);
   }}
 }
 return 0;
}

int Symmatrix::LDLinverse(Symmatrix &S) const
{
 chk_init(*this);
 S.newsize(nr);

 for(Integer e=0;e<nr;e++){
     Real* sp=S.m+e*nr-((e-1)*e)/2-e;

     //---- solve Lr=v
     sp[e]=1.;
     Integer indi=e*nr-((e-1)*e)/2+1;
     mat_xemy(nr-e-1,sp+e+1,m+indi);
     for(Integer i=e+1;i<nr;i++){
         indi+=nr-i+1;
         mat_xpeya(nr-i-1,sp+i+1,m+indi,-sp[i]);
     }
 
     //---- solve Dr=r
     {
       const Real* mp=m+e*nr-(e*(e-1))/2;
       Real* sp1=sp+e;
       for(Integer i=nr-e;i>0;i--) {
         (*sp1++)/=(*mp);
         mp+=i;
       }
     }
 
     //---- solve L'r=r
     indi=(nr*(nr+1))/2-1;
     {for(Integer i=nr;--i>=e;){
         const Real* mp=m+indi-(nr-i);
         Real *sp1=sp+i;
         const Real f=*sp1;
         for(Integer j=i;--j>=e;) {
             (*(--sp1))-=f*(*mp);
             mp-=nr-j;
         }
         indi-=nr-i+1;
     }} 
 }
 chk_set_init(S,1);
 return 0;
}

}

