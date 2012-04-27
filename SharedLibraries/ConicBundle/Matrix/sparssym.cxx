/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/sparssym.cxx

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
#include <algorithm>
#include "sparssym.hxx"
#include "heapsort.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

Sparsesym& Sparsesym::xeya(const Sparsesym& A,Real d)
//returns (*this)=A*d;
{
 chk_init(A);
 if (d==0.) return init(A.nr);
 nr=A.nr;
 colinfo=A.colinfo;
 colindex=A.colindex;
 colval.xeya(A.colval,d);
 suppind=A.suppind;
 suppcol=A.suppcol;
 chk_set_init(*this,1);
 return *this;
}

Sparsesym& Sparsesym::xeya(const Matrix& A,Real in_d)
//returns (*this)=A*d;
{
 chk_init(A);
 init(A.nr); 
 chk_add(*this,A);
 if (in_d==0.) return *this;
 Integer i,j;

 //--- count nonzeros in each column and nonzero columns of symmetrized matrix 
 Indexmatrix tmpcnt(nr+1,1,Integer(0));
 for(j=0;j<nr;j++){
     Real d=A(j,j);
     if (fabs(d)>tol) tmpcnt(0)++;
     for (i=j+1;i<nr;i++){
         Real d=(A(i,j)+A(j,i))/2.;
         if (fabs(d)>tol) tmpcnt(j+1)++;
     }
 }
 Integer nz=sum(tmpcnt);
 Integer nrcol=sum(tmpcnt>0);

 //--- get memory and copy information
 colinfo.init(nrcol,4,Integer(0));
 colindex.newsize(nz,1); chk_set_init(colindex,1);
 colval.newsize(nz,1); chk_set_init(colval,1);
 nz=0;
 nrcol=0;
 //diagonal part
 if (tmpcnt(nrcol)>0){
     colinfo(nrcol,0)=-1;
     colinfo(nrcol,2)=nz;
     for(i=0;i<nr;i++){
         Real d=A(i,i);
         if (fabs(d)>tol) {
             colindex(nz)=i;
             colval(nz)=in_d*d;
             nz++;
         }
     }
     colinfo(nrcol,1)=nz;
     nrcol++;
 }
 //offdiagonal part
 for(j=0;j<nr;j++){
     if (tmpcnt(j+1)==0) continue;
     colinfo(nrcol,0)=j;
     colinfo(nrcol,2)=nz;
     for (i=j+1;i<nr;i++){         
         Real d=(A(i,j)+A(j,i))/2.;
         if (fabs(d)>tol) {
             colindex(nz)=i-j;  //store index with shift
             colval(nz)=in_d*d;
             nz++;
         }             
     }
     colinfo(nrcol,1)=nz-colinfo(nrcol,2);
     nrcol++;
 }

 update_support();

 return *this;
}

Sparsesym& Sparsesym::xeya(const Indexmatrix& A,Real in_d)
//returns (*this)=A*d;
{
 chk_init(A);
 init(A.nr); 
 chk_add(*this,A);
 if (in_d==0.) return *this;
 Integer i,j;

 //--- count nonzeros in each column and nonzero columns of symmetrized matrix 
 Indexmatrix tmpcnt(nr+1,1,Integer(0));
 for(j=0;j<nr;j++){
     Real d=A(j,j);
     if (fabs(d)>tol) tmpcnt(0)++;
     for (i=j+1;i<nr;i++){
         Real d=(A(i,j)+A(j,i))/2.;
         if (fabs(d)>tol) tmpcnt(j+1)++;
     }
 }
 Integer nz=sum(tmpcnt);
 Integer nrcol=sum(tmpcnt>0);

 //--- get memory and copy information
 colinfo.init(nrcol,4,Integer(0));
 colindex.newsize(nz,1); chk_set_init(colindex,1);
 colval.newsize(nz,1); chk_set_init(colval,1);
 nz=0;
 nrcol=0;
 //diagonal part
 if (tmpcnt(0)>0){
     colinfo(nrcol,0)=-1;
     colinfo(nrcol,2)=nz;
     for(i=0;i<nr;i++){
         Real d=A(i,i);
         if (fabs(d)>tol) {
             colindex(nz)=i;
             colval(nz)=in_d*d;
             nz++;
         }
     }
     colinfo(nrcol,1)=nz;
     nrcol++;
 }
 //offdiagonal part
 for(j=0;j<nr;j++){
     if (tmpcnt(j+1)==0) continue;
     colinfo(nrcol,0)=j;
     colinfo(nrcol,2)=nz;
     for (i=j+1;i<nr;i++){         
         Real d=(A(i,j)+A(j,i))/2.;
         if (fabs(d)>tol) {
             colindex(nz)=i-j; //add shift
             colval(nz)=in_d*d;
             nz++;
         }             
     }
     colinfo(nrcol,1)=nz-colinfo(nrcol,2);
     nrcol++;
 }

 update_support();

 return *this;
}


Sparsesym& xeyapzb(Sparsesym& C,const Sparsesym& A,const Sparsesym& B,Real alpha,Real beta)
  //returns x= alpha*A+beta*B
  //x is initialized to the correct size
{
 chk_add(A,B);
 if (A.colinfo.rowdim()==0) return C.init(B,beta);
 if (B.colinfo.rowdim()==0) return C.init(A,alpha);
 C.init(A.nr);
 Integer nrcols=min(A.nr,A.colinfo.rowdim()+B.colinfo.rowdim());
 Integer nz=B.colindex.dim()+A.colindex.dim();
 C.colinfo.init(nrcols,4,Integer(0));
 C.colindex.newsize(nz,1); chk_set_init(C.colindex,1);
 C.colval.newsize(nz,1); chk_set_init(C.colval,1);
 Integer cind1=0;
 Integer cind2=0;
 Integer cv1=A.colinfo(cind1);
 Integer cv2=B.colinfo(cind2);
 Integer nz1=0;
 Integer nz2=0;
 nrcols=0;
 nz=0;
 do {
     if (cv1<cv2){
         C.colinfo(nrcols,0)=cv1;
         C.colinfo(nrcols,2)=nz;
         for(Integer i=A.colinfo(cind1,1);--i>=0;){
             C.colindex(nz)=A.colindex(nz1);
             C.colval(nz)=alpha*A.colval(nz1);
             nz++;
             nz1++;
         }
         C.colinfo(nrcols,1)=nz-C.colinfo(nrcols,2);
         nrcols++;
         cind1++;
         if (cind1<A.colinfo.rowdim()){
             cv1=A.colinfo(cind1,0);
         }
         else {
             cv1=A.nr;
         }
     }
     else if (cv2<cv1) {
         C.colinfo(nrcols,0)=cv2;
         C.colinfo(nrcols,2)=nz;
         for(Integer i=B.colinfo(cind2,1);--i>=0;){
             C.colindex(nz)=B.colindex(nz2);
             C.colval(nz)=beta*B.colval(nz2);
             nz++;
             nz2++;
         }
         C.colinfo(nrcols,1)=nz-C.colinfo(nrcols,2);
         nrcols++;
         cind2++;
         if (cind2<B.colinfo.rowdim()){
             cv2=B.colinfo(cind2,0);
         }
         else {
             cv2=B.nr;
         }
     }
     else { //cv1==cv2, both work on the same column
         C.colinfo(nrcols,0)=cv1;
         C.colinfo(nrcols,2)=nz;
         Integer uind1=A.colinfo(cind1,2)+A.colinfo(cind1,1);
         Integer uind2=B.colinfo(cind2,2)+B.colinfo(cind2,1);
         Integer rind1,rind2;
         rind1= (nz1<uind1)? A.colindex(nz1):A.nr;
         rind2= (nz2<uind2)? B.colindex(nz2):B.nr;
         while((rind1<A.nr)||(rind2<B.nr)){
             if (rind1<rind2){
                 C.colindex(nz)=rind1;
                 C.colval(nz)=alpha*A.colval(nz1);
                 nz++;
                 nz1++;
                 rind1= (nz1<uind1)? A.colindex(nz1):A.nr;
             }
             else if(rind2<rind1){
                 C.colindex(nz)=rind2;
                 C.colval(nz)=beta*B.colval(nz2);
                 nz++;
                 nz2++;
                 rind2= (nz2<uind2)? B.colindex(nz2):B.nr;
             }
             else {
                 Real d=alpha*A.colval(nz1)+beta*B.colval(nz2);
                 if (fabs(d)>C.tol){
                     C.colindex(nz)=rind1;
                     C.colval(nz)=d;
                     nz++;
                 }
                 nz1++;
                 rind1= (nz1<uind1)? A.colindex(nz1):A.nr;
                 nz2++;
                 rind2= (nz2<uind2)? B.colindex(nz2):B.nr;
             }
         }
         C.colinfo(nrcols,1)=nz-C.colinfo(nrcols,2);
         nrcols++;
         cind1++;
         if (cind1<A.colinfo.rowdim()){
             cv1=A.colinfo(cind1,0);
         }
         else {
             cv1=A.nr;
         }
         cind2++;
         if (cind2<B.colinfo.rowdim()){
             cv2=B.colinfo(cind2,0);
         }
         else {
             cv2=B.nr;
         }                           
     }    
 }while((cv1<A.nr)||(cv2<B.nr));
 C.colindex.reduce_length(nz);
 C.colval.reduce_length(nz);
 Indexmatrix rmind(C.colinfo.rowdim(),1); chk_set_init(rmind,1);
 nz=0;
 for(Integer i=0;i<C.colinfo.rowdim();i++){
     if (C.colinfo(i,1)==0) rmind(nz++)=i;
 }
 rmind.reduce_length(nz);
 C.colinfo.delete_rows(rmind);

 C.update_support();
 
 return C;
}


// ============================================================================
//                              support_rankadd 
// ============================================================================

//to the support of *this the upper triangle of d*A'*B is added

Sparsesym& support_rankadd(const Matrix& A,Sparsesym& C,
                                      Real alpha,Real beta,int trans)
{
 chk_init(A);
 chk_init(C);
 if (trans){ //A^T*A
     chk_colseq(C,A);
     if (beta==0.) {
         C.colval.init(C.colval.rowdim(),1,0.);
     }
     else {
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer n=A.rowdim();
     Integer nz=0;
     for(Integer j=0;j<C.colinfo.rowdim();j++){
         Integer jj=C.colinfo(j,0);
         if (jj<0){ //diagonal case
             for(Integer i=C.colinfo(j,1);--i>=0;){
                 Integer ii=C.colindex(nz);
                 C.colval(nz++)+=alpha*mat_ip(n,A.get_store()+ii*n,A.get_store()+ii*n);
             }
         }
         else { //offdiagonal case
             for(Integer i=C.colinfo(j,1);--i>=0;){
                 Integer ii=C.colindex(nz)+jj;
                 C.colval(nz++)+=alpha*mat_ip(n,A.get_store()+jj*n,A.get_store()+ii*n);
             }
         }
     }
     return C;
 }
 else { //A*A^T
     chk_rowseq(C,A);
     if (beta==0.) {
         C.colval.init(C.colval.rowdim(),1,0.);
     }
     else {
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer n=A.rowdim();
     Integer nc=A.coldim();
     Integer nz=0;
     for(Integer j=0;j<C.colinfo.rowdim();j++){
         Integer jj=C.colinfo(j,0);
         if (jj<0){ //diagonal case
             for(Integer i=C.colinfo(j,1);--i>=0;){
                 Integer ii=C.colindex(nz);
                 C.colval(nz++)+=alpha*mat_ip(nc,A.get_store()+ii,n,A.get_store()+ii,n);
             }
         }
         else { //offdiagonal case
             for(Integer i=C.colinfo(j,1);--i>=0;){
                 Integer ii=C.colindex(nz)+jj;
                 C.colval(nz++)+=alpha*mat_ip(nc,A.get_store()+jj,n,A.get_store()+ii,n);
             }
         }
     }
     return C;
 }
}

// ============================================================================
//                              genmult with Matrix
// ============================================================================

//returns C=beta*C+alpha*A*B, where the Matrix may be transposed
 
Matrix& genmult(const Sparsesym& A,const Matrix& B,Matrix &C,
			 Real alpha,Real beta, int btrans)
{
 chk_init(A);
 chk_init(B);
 Integer nr,nc;
 nr=A.nr;
 if (btrans) {
     nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
     if (A.nr!=B.nc) {
         MEmessage(MatrixError(ME_dim,"genmult(Sparsesym,Matrix,Matrix): dimensions don't match",MTsparsesym));
     }
#endif
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (A.nr!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult(Sparsesym,Matrix,Matrix): dimensions don't match",MTsparsesym));
     }
#endif
 }
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.nr)||(nc!=C.nc)) {
         MEmessage(MatrixError(ME_dim,"genmult(Sparsesym,Matrix,Matrix): dimensions don't match",MTsparsesym));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if (btrans){
   Integer nrcol=A.colinfo.rowdim();
   if ((nrcol==0)||(nc==0)) return C;

   const Integer *jp=A.colinfo.get_store();
   const Integer *cntp=A.colinfo.get_store()+A.colinfo.rowdim();
   register const Integer *ip=A.colindex.get_store();
   register const Real *vp=A.colval.get_store();
   const Real *bbp=B.get_store();
   Integer bnr=B.rowdim();
   Real *cbp=C.get_store();
 
   if (nc>1){ //more than only one column in B
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++,vp++){
	 mat_xpeya(nc,cbp+*ip,nr,bbp+*ip*bnr,1,*vp*alpha);
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       bbp=B.get_store()+*jp*bnr;
       cbp=C.get_store()+*jp++;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 mat_xpeya(nc,cbp,nr,bbp+h*bnr,1,d);
	 mat_xpeya(nc,cbp+h,nr,bbp,1,d);
       }
     }
   }
   else { //B is a vector (only one row, may be treated as vector))
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++){
	 *(cbp+*ip)+=*(bbp+*ip)*(*vp++)*alpha;
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       bbp=B.get_store()+*jp;
       cbp=C.get_store()+*jp++;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 *cbp+=*(bbp+h)*d;
	 *(cbp+h)+=*bbp*d;
       }
     }
   }   
 }
 else { //B is not transposed
   Integer nrcol=A.colinfo.rowdim();
   if ((nrcol==0)||(nc==0)) return C;

   const Integer *jp=A.colinfo.get_store();
   const Integer *cntp=A.colinfo.get_store()+A.colinfo.rowdim();
   register const Integer *ip=A.colindex.get_store();
   register const Real *vp=A.colval.get_store();
   const Real *bbp=B.get_store();
   Real *cbp=C.get_store();
 
   if (nc>1){ //more than only one column in B
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++,vp++){
	 mat_xpeya(nc,cbp+*ip,nr,bbp+*ip,nr,*vp*alpha);
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       bbp=B.get_store()+*jp;
       cbp=C.get_store()+*jp++;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 mat_xpeya(nc,cbp,nr,bbp+h,nr,d);
	 mat_xpeya(nc,cbp+h,nr,bbp,nr,d);
       }
     }
   }
   else { //B is a vector (only one column)
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++){
	 *(cbp+*ip)+=*(bbp+*ip)*(*vp++)*alpha;
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       bbp=B.get_store()+*jp;
       cbp=C.get_store()+*jp++;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 *cbp+=*(bbp+h)*d;
	 *(cbp+h)+=*bbp*d;
       }
     }
   }   
 }
 return C;
}

Matrix& genmult(const Matrix& A,const Sparsesym& B,Matrix& C,
                    Real alpha,Real beta,int atrans)
{
 chk_init(A);
 chk_init(B);
 Integer nr,nc;
 if (atrans){
     nr=A.nc;
 }
 else {
     nr=A.nr;
 }
 nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
 Integer nm;
 if (atrans){
     nm=A.nr;
 }
 else {
     nm=A.nc;
 }
 if (nm!=B.nr) {
     MEmessage(MatrixError(ME_dim,"genmult(Matrix,Sparsesym,Matrix): dimensions don't match",MTsparsesym));
 }
#endif
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.nr)||(nc!=C.nc)) {
         MEmessage(MatrixError(ME_dim,"genmult(Matrix,Sparsesym,Matrix): dimensions don't match",MTsparsesym));
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if (atrans){
   Integer nrcol=B.colinfo.rowdim();
   if ((nrcol==0)||(nc==0)) return C;
   
   const Integer *jp=B.colinfo.get_store();
   const Integer *cntp=B.colinfo.get_store()+B.colinfo.rowdim();
   register const Integer *ip=B.colindex.get_store();
   register const Real *vp=B.colval.get_store();
   const Real *abp=A.m;
   Integer anr=A.nr;
   Real *cbp=C.get_store();
   
   if (nc>1){ //Matrix A consists of several columns
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++,vp++){
	 mat_xpeya(nr,cbp+*ip*nr,1,abp+*ip,anr,*vp*alpha);
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       abp=A.m+*jp;
       cbp=C.get_store()+*(jp++)*nr;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 mat_xpeya(nr,cbp,1,abp+h,anr,d);
	 mat_xpeya(nr,cbp+h*nr,1,abp,anr,d);
       }
     }
   }
   else { //A is a column vector (only one column)
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++){
	 *(cbp+*ip)+=*(abp+*ip)*(*vp++)*alpha;
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       abp=A.m+*jp;
       cbp=C.get_store()+*jp++;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 *cbp+=*(abp+h)*d;
	 *(cbp+h)+=*abp*d;
       }
     }
   }
 }
 else{
   Integer nrcol=B.colinfo.rowdim();
   if ((nrcol==0)||(nc==0)) return C;
   
   const Integer *jp=B.colinfo.get_store();
   const Integer *cntp=B.colinfo.get_store()+B.colinfo.rowdim();
   register const Integer *ip=B.colindex.get_store();
   register const Real *vp=B.colval.get_store();
   const Real *abp=A.m;
   Real *cbp=C.get_store();
   
   if (nc>1){ //Matrix A consists of several rows
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++,vp++){
	 mat_xpeya(nr,cbp+*ip*nr,abp+*ip*nr,*vp*alpha);
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       abp=A.m+*jp*nr;
       cbp=C.get_store()+*(jp++)*nr;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 mat_xpeya(nr,cbp,abp+h*nr,d);
	 mat_xpeya(nr,cbp+h*nr,abp,d);
       }
     }
   }
   else { //A is a column vector (only one column)
     //diagonal part
     if (*jp<0) {
       for(Integer i=*cntp++;--i>=0;ip++){
	 *(cbp+*ip)+=*(abp+*ip)*(*vp++)*alpha;
       }
       --nrcol;
       jp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
       abp=A.m+*jp;
       cbp=C.get_store()+*jp++;
       for(Integer i=*cntp++;--i>=0;){
	 Integer h=*ip++;
	 Real d=*vp++*alpha;
	 *cbp+=*(abp+h)*d;
	 *(cbp+h)+=*abp*d;
       }
     }
   }
 }
 return C;
}

// ============================================================================
//                                  update_support
// ============================================================================

void Sparsesym::update_support()
{
 //remove shift of offdiagonal entries and zero entries in colind and colval
 Integer nrcol=colinfo.rowdim();
 if (nrcol==0) return;
 Indexmatrix sind(colindex.rowdim(),1); chk_set_init(sind,1);
 Integer infind=0;
 Integer oldnz=0;
 Integer newnz=0;
 for(Integer j=0;j<nrcol;j++){
     Integer jj=colinfo(j,0);
     colinfo(j,2)=newnz;
     if (jj<0) {
         for(Integer i=colinfo(j,1);--i>=0;){
             if (fabs(colval(oldnz))>tol){
                 colindex(newnz)=colindex(oldnz);
                 colval(newnz)=colval(oldnz);
                 newnz++;
                 oldnz++;
             }
             else {
                 colinfo(j,1)--;
                 oldnz++;
             }
         }        
     }
     else { 
         for(Integer i=colinfo(j,1);--i>=0;){
             if (fabs(colval(oldnz))>tol){
                 colindex(newnz)=colindex(oldnz)+jj;
                 colval(newnz)=colval(oldnz);
                 newnz++;
                 oldnz++;
             }
             else {
                 colinfo(j,1)--;
                 oldnz++;
             }
         }        
     }
     if (colinfo(j,1)==0) {
         sind(infind)=j;
         infind++;
     }
 }
 colindex.reduce_length(newnz);
 colval.reduce_length(newnz);
 sind.reduce_length(infind);
 colinfo.delete_rows(sind);
 if (colinfo.rowdim()==0) return;
 
 //find the columns that contain support
 suppind.newsize(colindex.rowdim(),1); chk_set_init(suppind,1);
 sortindex(colindex,sind);
 Integer suppcnt=0;
 infind=0;
 Integer infocol=colinfo(infind,0);
 if (infocol<0){
     colinfo(infind,3)=-1;
     infind++;
     infocol=(infind<colinfo.rowdim())?colinfo(infind,0):nr;
 }
 Integer lnz=0;
 Integer elemcol=colindex(sind(lnz));
 Integer suppelemcnt=0;  //counts the first appearances of supports in colindex
                         //the visited entries of sind can then be used to store
                         //the support columns appearing in colindex in increasing order
 Integer lastj;          //holds most recent support column added
 if (infocol<=elemcol){
     lastj=infocol;
     colinfo(infind,3)=suppcnt;
     infind++;
     infocol=(infind<colinfo.rowdim())?colinfo(infind,0):nr;
 }
 else {
     lastj=elemcol;
     suppind(sind(lnz))=suppcnt;    
     sind(suppelemcnt++)=sind(lnz);     //concentrate support information in first fields
     lnz++;
     elemcol=(lnz<colindex.rowdim())?colindex(sind(lnz)):nr;
 }
 while((infocol<nr)||(elemcol<nr)){
      if (infocol<=elemcol) {
          if (lastj<infocol) {
              suppcnt++;
              lastj=infocol;
          }
          colinfo(infind,3)=suppcnt;
          infind++;
          infocol=(infind<colinfo.rowdim())?colinfo(infind,0):nr;
      }
      else {
          if (lastj<elemcol) {
              suppcnt++;
              lastj=elemcol;
              sind(suppelemcnt++)=sind(lnz);  //concentrate support info in first fields
          }
          suppind(sind(lnz))=suppcnt;
          lnz++;
          elemcol=(lnz<colindex.rowdim())?colindex(sind(lnz)):nr;
      }
 }
 suppcnt++;

 //collect support to column information in suppcol, takes only suppcnt iterations now
 suppcol.newsize(suppcnt,1); chk_set_init(suppcol,1);
 infind=0;
 infocol=colinfo(infind,0);
 if (infocol<0){
     infind++;
     infocol=(infind<colinfo.rowdim())?colinfo(infind,0):nr;
 }
 lnz=0;
 elemcol=(lnz<suppelemcnt)?colindex(sind(lnz)):nr;
 suppcnt=0;
 while((infocol<nr)||(elemcol<nr)){
     //(elemcol != infocol) by construction
     if (infocol<elemcol){
         suppcol(suppcnt++)=infocol;
         infind++;
         infocol=(infind<colinfo.rowdim())?colinfo(infind,0):nr;
     }
     else { 
         suppcol(suppcnt++)=elemcol;         //make sure sind is used 
         lnz++;
         elemcol=(lnz<suppelemcnt)?colindex(sind(lnz)):nr;
     }
 }
         
 //shift indices on offdiagonal entries (helps to speed up multiplication)
 lnz=0;
 {for(Integer j=0;j<infind;j++){
     Integer jj=colinfo(j,0);
     Integer ss=colinfo(j,3);
     if (jj<0) {
         lnz+=colinfo(j,1);
         continue;
     }
     for(Integer i=colinfo(j,1);--i>=0;){
         colindex(lnz)-=jj;
         suppind(lnz++)-=ss;
     }
 }}
}


Sparsesym& Sparsesym::init(Integer in_nr,Integer nz,
              const Integer *ini,const Integer *inj,const Real* va)
{
 return init(in_nr,nz,Indexmatrix(nz,1,ini),Indexmatrix(nz,1,inj),Matrix(nz,1,va));
}

// ============================================================================
//                                    init
// ============================================================================

// initialize Sparsesym from an edge representation. Multiple edges are
// added to form a single edge.

Sparsesym& Sparsesym::init(Integer in_nr,Integer nz,
          const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va)
{
 chk_init(ini);
 chk_init(inj);
 chk_init(va);
#if (CONICBUNDLE_DEBUG>=1)
 if (nz<0) {
     MEmessage(MatrixError(ME_unspec,"Sparsesym::init(nr,nc,nz,indi,indj,val): nz<0",MTsparse));
 }
#endif
 if (nz==0) {
     return init(in_nr);
 }
#if (CONICBUNDLE_DEBUG>=1)
 if ((ini.dim()<nz)||(inj.dim()<nz)||(va.dim()<nz)) {  
     MEmessage(MatrixError(ME_unspec,"Sparsesym::init(nr,nz,indi,indj,val): indi, indj, or val has <nz elements",MTsparse));
 }
 if ((min(ini(Range(0,nz-1)))<0)||(max(ini(Range(0,nz-1)))>=in_nr)||
     (min(inj(Range(0,nz-1)))<0)||(max(inj(Range(0,nz-1)))>=in_nr)) {  
     MEmessage(MatrixError(ME_unspec,"Sparsesym::init(nr,nz,indi,indj,val): indices in indi or indj exceed range",MTsparse));
 }
#endif

 nr=in_nr;

 if (nz==1) {
     colinfo.newsize(1,4); chk_set_init(colinfo,1);
     Integer ii=ini(0);
     Integer jj=inj(0);
     if (ii==jj){
         colinfo(0,0)=-1; //signal for diagonal
         colinfo(0,3)=-1;
         colindex.init(1,1,ii);
	 suppind.init(1,1,Integer(0));
         suppcol.init(1,1,ii);
     }
     else {
         if (ii<jj) {Integer h=ii;ii=jj;jj=h;}
         colinfo(0,0)=jj;
         colinfo(0,3)=0;
         colindex.init(1,1,ii-jj); //shift index, beneficial for multiplication
	 suppind.init(1,1,Integer(1));
         suppcol.newsize(2,1); chk_set_init(suppcol,1);
         suppcol(0)=jj;
         suppcol(1)=ii;
     } 
     colinfo(0,1)=1;
     colinfo(0,2)=0;
     colval.init(1,1,va(0));
     return *this;
 }
 else {    
 colindex.newsize(nz,1); chk_set_init(colindex,1);
 suppind.newsize(nz,1); chk_set_init(suppind,1);
 colval.newsize(nz,1); chk_set_init(colval,1);
 
 Integer i,j;

 //--- prepare indices for sorting by guaranteeing i>j
 for(i=0;i<nz;i++){
     Integer ii=ini(i);
     Integer jj=inj(i);
     if (ii<=jj){
         if (ii==jj){ //make negative to move to front
             jj=Integer(-1);
         }
         else {
             Integer h=ii;ii=jj;jj=h;
         }
     }
     colindex(i)=ii;  
     suppind(i)=jj;  
 }

 //--- sort indices by column first
 Indexmatrix sind;
 sortindex(suppind,sind);
 
 //--- run through all elements to count different column indices and to
 //--- sort the corresponding row indices
 Integer nrj=1;
 Integer si=sind(0);
 Integer lastj=suppind(si);
 Integer jcnt=1;
 for(i=1;i<nz;i++){
     si=sind(i);
     if (suppind(si)>lastj){ //a new column starts, sort rows of previous
         Integer *sindp=sind.get_store()+i-jcnt;
         //heapsort(jcnt,sindp,colindex.get_store());
	 std::sort( sindp,sindp+jcnt, mat_less_index<Integer>(colindex.get_store()) );
         jcnt=1;
         nrj++;
         lastj=suppind(si);
         continue;
     }
     jcnt++;
 }
 Integer *sindp=sind.get_store()+i-jcnt;
 //heapsort(jcnt,sindp,colindex.get_store());
 std::sort( sindp,sindp+jcnt, mat_less_index<Integer>(colindex.get_store()) );
         
 //--- resize colinfo and run through the indices to store the 
 //--- sorted column data with corresponding info
 colinfo.newsize(nrj,4); chk_set_init(colinfo,1);
 Integer lnz=0;   //local nonzero counter for contracting identical entries
 nrj=0;
 si=sind(0);
 colinfo(nrj,0)=lastj=suppind(si);
 colinfo(nrj,2)=lnz;
 colindex(lnz)=max(ini(si),inj(si));
 colval(lnz)=va(si);
 lnz++;
 for(i=1;i<nz;i++){
     si=sind(i);
     Integer jj=suppind(si);
     Integer ii=max(ini(si),inj(si));
     if ((jj==lastj)&&(ii==colindex(lnz-1))){
         colval(lnz-1)+=va(si);
         continue;
     }
     if (jj>lastj){ //a new column starts, store info, initialize new
         colinfo(nrj,1)=lnz-colinfo(nrj,2);
         nrj++;
         colinfo(nrj,0)=lastj=suppind(si);
         colinfo(nrj,2)=lnz;
     }
     colindex(lnz)=ii;
     colval(lnz)=va(si);
     lnz++;
 }
 colinfo(nrj,1)=lnz-colinfo(nrj,2);
 nrj++;
 colindex.reduce_length(lnz);
 colval.reduce_length(lnz);

 //shift indices on offdiagonal entries (helps to speed up multiplication)
 lnz=0;
 for(j=0;j<nrj;j++){
     Integer jj=colinfo(j,0);
     if (jj<0) {
         lnz+=colinfo(j,1);
         continue;
     }
     for(i=colinfo(j,1);--i>=0;){
         colindex(lnz++)-=jj;
     }
 }
 }


 update_support();

 return *this;
}

// ============================================================================
//                              get_edge_rep
// ============================================================================

void Sparsesym::get_edge_rep(Indexmatrix& I, Indexmatrix& J, Matrix& val) const
{
 I.newsize(colindex.dim(),1); chk_set_init(I,1);
 J.newsize(colindex.dim(),1); chk_set_init(J,1);
 val.newsize(colindex.dim(),1); chk_set_init(val,1);
 Integer i,j;
 Integer nz=0;
 for(i=0;i<colinfo.rowdim();i++){
     Integer col=colinfo(i,0);
     Integer nrind=colinfo(i,1);
     if (col<0) {
         for(j=0;j<nrind;j++){
             J(nz)=I(nz)=colindex(nz);
             val(nz)=colval(nz);
             nz++;
         }
     }
     else {
         for(j=0;j<nrind;j++){
             I(nz)=col;
             J(nz)=colindex(nz)+col;
             val(nz)=colval(nz);
             nz++;
         }
     }
 }
}

// ============================================================================
//                               contains_support
// ============================================================================

int Sparsesym::contains_support(const Sparsesym& A) const
{
 chk_add(*this,A);
 if (A.colinfo.rowdim()==0) return 1;
 if ((colinfo.rowdim()==0) && (A.colinfo.rowdim()==0)) return 1;
 if (colval.rowdim()<A.colval.rowdim()) return 0;
 if (colinfo.rowdim()<A.colinfo.rowdim()) return 0;
 if (suppcol.rowdim()<A.suppcol.rowdim()) return 0;
 Integer cind1=0;
 Integer cind2=0;
 Integer cv1=colinfo(cind1);
 Integer cv2=A.colinfo(cind2);
 Integer nz1=0;
 Integer nz2=0;
 do {
     if (cv1<cv2){
         nz1+=colinfo(cind1,1);
         cind1++;
         if (cind1<colinfo.rowdim()){
             cv1=colinfo(cind1,0);
         }
         else {
             cv1=nr;
         }
     }
     else if (cv2<cv1) {
         return 0;
     }
     else { //cv1==cv2, both work on the same column
         Integer uind1=colinfo(cind1,2)+colinfo(cind1,1);
         Integer uind2=A.colinfo(cind2,2)+A.colinfo(cind2,1);
         Integer rind1,rind2;
         rind1= (nz1<uind1)? colindex(nz1):nr;
         rind2= (nz2<uind2)? A.colindex(nz2):nr;
         while((rind1<nr)||(rind2<nr)){
             if (rind1<rind2){
                 nz1++;
                 rind1= (nz1<uind1)? colindex(nz1):nr;
             }
             else if(rind2<rind1){
                 return 0;
             }
             else {
                 nz1++;
                 rind1= (nz1<uind1)? colindex(nz1):nr;
                 nz2++;
                 rind2= (nz2<uind2)? A.colindex(nz2):nr;
             }
         }
         cind1++;
         if (cind1<colinfo.rowdim()){
             cv1=colinfo(cind1,0);
         }
         else {
             cv1=nr;
         }
         cind2++;
         if (cind2<A.colinfo.rowdim()){
             cv2=A.colinfo(cind2,0);
         }
         else {
             cv2=nr;
         }
     }
 }while((cv1<nr)||(cv2<nr));
 
 return 1;
}

// ============================================================================
//                                  init_support
// ============================================================================
    
Sparsesym& Sparsesym::init_support(const Sparsesym& A,Real d)
{
 chk_init(A);
 nr=A.nr;
 colinfo=A.colinfo;
 colindex=A.colindex;
 colval.init(A.colval.rowdim(),1,d);
 suppind=A.suppind;
 suppcol=A.suppcol;
 chk_set_init(*this,1);
 return *this;
}

// ============================================================================
//                                  operator(i,j)
// ============================================================================
    
Real Sparsesym::operator()(Integer i,Integer j) const
{
 chk_init(*this);
 chk_range(i,j,nr,nr);
 Integer ncols=colinfo.rowdim();
 if (ncols==0) return 0.;
 if (i<j) {Integer h=i;i=j;j=h;}
 if (i==j) j=-1; //flag for diagonal
 else {
     i-=j;       //add shift
 }
 if ((j<colinfo(0,0))||(j>colinfo(ncols-1,0))) return 0.;
 //binary search for j
 Integer li=0; 
 Integer ui=ncols-1;
 Integer mi=0;
 while(li<=ui) {
     mi=(li+ui)/2;
     if (colinfo(mi,0)==j) break;
     if (colinfo(mi,0)<j){
         li=mi+1;
     }
     else {
         ui=mi-1;
     }
 }
 if (li>ui) return 0.;
 //j found, binary search for i
 li=colinfo(mi,2);
 ui=li+colinfo(mi,1)-1;
 while(li<=ui) {
     mi=(li+ui)/2;
     if (colindex(mi)==i) break;
     if (colindex(mi)<i){
         li=mi+1;
     }
     else {
         ui=mi-1;
     }
 } 
 if (li>ui) return 0.;
 return colval(mi);
}

// ============================================================================
//                                check_support(i,j)
// ============================================================================
    
int Sparsesym::check_support(Integer i,Integer j) const
{
 chk_init(*this);
 chk_range(i,j,nr,nr);
 Integer ncols=colinfo.rowdim();
 if (ncols==0) return 0;
 if (i<j) {Integer h=i;i=j;j=h;}
 if (i==j) j=-1; //flag for diagonal
 else {
     i-=j;       //add shift
 }
 if ((j<colinfo(0,0))||(j>colinfo(ncols-1,0))) return 0;
 //binary search for j
 Integer li=0; 
 Integer ui=ncols-1;
 Integer mi=0;
 while(li<=ui) {
     mi=(li+ui)/2;
     if (colinfo(mi,0)==j) break;
     if (colinfo(mi,0)<j){
         li=mi+1;
     }
     else {
         ui=mi-1;
     }
 }
 if (li>ui) return 0;
 //j found, binary search for i
 li=colinfo(mi,2);
 ui=li+colinfo(mi,1)-1;
 while(li<=ui) {
     mi=(li+ui)/2;
     if (colindex(mi)==i) break;
     if (colindex(mi)<i){
         li=mi+1;
     }
     else {
         ui=mi-1;
     }
 } 
 if (li>ui) return 0;
 return 1;
}

        

void Sparsesym::display(std::ostream& out,
                  int precision,int width,int screenwidth) const
           //for variables of value zero default values are used
           //precision=4,width=precision+6,screenwidth=80
{
 chk_init(*this);
 if (precision==0) precision=4;
 if (width==0) width=precision+6;
 if (screenwidth==0) screenwidth=80;
 out.precision(precision);
 out<<"Sparsesym("<<nr<<") "<<" "<<colindex.dim()<<"\n";
 Integer nz=0;
 for(Integer j=0;j<colinfo.rowdim();j++){
     Integer jj=colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=colinfo(j,1);--i>=0;){
             Integer ii=colindex(nz);
             out.width(5);out<<ii<<" ";
             out.width(5);out<<ii<<" ";
             out.width(width);out<<colval(nz++)<<"\n";
         }
     }
     else { //offdiagonal case
         for(Integer i=colinfo(j,1);--i>=0;){
             Integer ii=colindex(nz)+jj;
             out.width(5);out<<jj<<" ";
             out.width(5);out<<ii<<" ";
             out.width(width);out<<colval(nz++)<<"\n";
        }
     }
 }
}

// *****************************************************************************
//                  Interaction with other classes
// *****************************************************************************


// ******************************    Matrix   ******************************

Matrix& Matrix::xeya(const Sparsesym& A,Real d)
{
 chk_init(A);
 init(A.nr,A.nr,0.);
 Integer nz=0;
 if (d==1.){
   for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz);
	 (*this)(ii,ii)=A.colval(nz++);
       }
     }
     else { //offdiagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz)+jj;
	 Real Ad=A.colval(nz++);
	 (*this)(ii,jj)=Ad;
	 (*this)(jj,ii)=Ad;
       }
     }
   }
 }
 else if (d==-1.){
   for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz);
	 (*this)(ii,ii)=-A.colval(nz++);
       }
     }
     else { //offdiagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz)+jj;
	 Real Ad=A.colval(nz++);
	 (*this)(ii,jj)=-Ad;
	 (*this)(jj,ii)=-Ad;
       }
     }
   }
 }
 else {
   for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz);
	 (*this)(ii,ii)=A.colval(nz++)*d;
       }
     }
     else { //offdiagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz)+jj;
	 Real Ad=d*A.colval(nz++);
	 (*this)(ii,jj)=Ad;
	 (*this)(jj,ii)=Ad;
       }
     }
   }
 }
 return *this;
}

Matrix& Matrix::xpeya(const Sparsesym& A,Real d)
{
 chk_add(*this,A);
 Integer nz=0;
 if (d==1.){
   for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz);
	 (*this)(ii,ii)+=A.colval(nz++);
       }
     }
     else { //offdiagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz)+jj;
	 Real Ad=A.colval(nz++);
	 (*this)(ii,jj)+=Ad;
	 (*this)(jj,ii)+=Ad;
       }
     }
   }
 }
 else if (d==-1.){
   for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz);
	 (*this)(ii,ii)-=A.colval(nz++);
       }
     }
     else { //offdiagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz)+jj;
	 Real Ad=A.colval(nz++);
	 (*this)(ii,jj)-=Ad;
	 (*this)(jj,ii)-=Ad;
       }
     }
   }
 }
 else {
   for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz);
	 (*this)(ii,ii)+=A.colval(nz++)*d;
       }
     }
     else { //offdiagonal case
       for(Integer i=A.colinfo(j,1);--i>=0;){
	 Integer ii=A.colindex(nz)+jj;
	 Real Ad=d*A.colval(nz++);
	 (*this)(ii,jj)+=Ad;
	 (*this)(jj,ii)+=Ad;
       }
     }
   }
 }
 return *this;
}



// ******************************    Symmatrix   **************************

Sparsesym& Sparsesym::xeya(const Symmatrix &A,Real in_d)
{
 chk_init(A);
 init(A.nr); 
 if (in_d==0.) return *this;
 Integer i,j;

 //--- count nonzeros in each column and nonzero columns
 Indexmatrix tmpcnt(nr+1,1,Integer(0));
 for(j=0;j<nr;j++){
     Real d=A(j,j);
     if (fabs(d)>tol) tmpcnt(0)++;
     for (i=j+1;i<nr;i++){
         Real d=A(i,j);
         if (fabs(d)>tol) tmpcnt(j+1)++;
     }
 }
 Integer nz=sum(tmpcnt);
 Integer nrcol=sum(tmpcnt>0);

 //--- get memory and copy information
 colinfo.init(nrcol,4,Integer(0));
 colindex.newsize(nz,1); chk_set_init(colindex,1);
 colval.newsize(nz,1); chk_set_init(colval,1);
 nz=0;
 nrcol=0;
 //diagonal part
 if (tmpcnt(0)>0){
     colinfo(nrcol,0)=-1;
     colinfo(nrcol,2)=nz;
     for(i=0;i<nr;i++){
         Real d=A(i,i);
         if (fabs(d)>tol) {
             colindex(nz)=i;
             colval(nz)=in_d*d;
             nz++;
         }
     }
     colinfo(nrcol,1)=nz;
     nrcol++;
 }
 //offdiagonal part
 for(j=0;j<nr;j++){
     if (tmpcnt(j+1)==0) continue;
     colinfo(nrcol,0)=j;
     colinfo(nrcol,2)=nz;
     for (i=j+1;i<nr;i++){         
         Real d=A(i,j);
         if (fabs(d)>tol) {
             colindex(nz)=i-j;  //add shift
             colval(nz)=in_d*d;
             nz++;
         }             
     }
     colinfo(nrcol,1)=nz-colinfo(nrcol,2);
     nrcol++;
 }

 update_support();
     
 return *this;
}

Symmatrix& Symmatrix::xeya(const Sparsesym& A,Real d)
{
 chk_init(A);
 init(A.nr,0.); 
 if (d==0.) { return *this;}
 if (d==1.) {
     Integer nz=0;
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         if (jj<0){ //diagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz);
                 (*this)(ii,ii)=A.colval(nz++);
             }
         }
         else { //offdiagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz)+jj;
                 (*this)(ii,jj)=A.colval(nz++);
             }
         }
     }
     return *this;
 }
 if (d==-1.){
     Integer nz=0;
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         if (jj<0){ //diagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz);
                 (*this)(ii,ii)=-A.colval(nz++);
             }
         }
         else { //offdiagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz)+jj;
                 (*this)(ii,jj)=-A.colval(nz++);
             }
         }
     }
     return *this;
 }
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz);
             (*this)(ii,ii)=d*A.colval(nz++);
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz)+jj;
             (*this)(ii,jj)=d*A.colval(nz++);
         }
     }
 }
 return *this;
}

Symmatrix& Symmatrix::xpeya(const Sparsesym& A,Real d)
{
 chk_add(*this,A);
 if (d==0.) { return *this;}
 if (d==1.) {
     Integer nz=0;
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         if (jj<0){ //diagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz);
                 (*this)(ii,ii)+=A.colval(nz++);
             }
         }
         else { //offdiagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz)+jj;
                 (*this)(ii,jj)+=A.colval(nz++);
             }
         }
     }
     return *this;
 }
 if (d==-1.){
     Integer nz=0;
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         if (jj<0){ //diagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz);
                 (*this)(ii,ii)-=A.colval(nz++);
             }
         }
         else { //offdiagonal case
             for(Integer i=A.colinfo(j,1);--i>=0;){
                 Integer ii=A.colindex(nz)+jj;
                 (*this)(ii,jj)-=A.colval(nz++);
             }
         }
     }
     return *this;
 }
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz);
             (*this)(ii,ii)+=d*A.colval(nz++);
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz)+jj;
             (*this)(ii,jj)+=d*A.colval(nz++);
         }
     }
 }
 return *this;
}

Real ip(const Symmatrix& B,const Sparsesym& A)
{
 chk_add(B,A);
 Real s=0;
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz);
             s+=B(ii,ii)*A.colval(nz++);
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz)+jj;
             s+=2*B(ii,jj)*A.colval(nz++);
         }
     }
 }
 return s;
}

 
// ******************************    Sparsemat   **************************


//C=beta*C+alpha*A*B;

Matrix& genmult(const Sparsesym& A,const Sparsemat& B,Matrix &C,
                           Real alpha,Real beta,int btrans)
{
 chk_init(A);
 chk_init(B);
 Integer nr,nc;
 const Indexmatrix* ainf;
 const Indexmatrix* aind;
 const Matrix* aval;
 nr=A.nr;
#if (CONICBUNDLE_DEBUG>=1)
 Integer nm=A.nr;
#endif
 ainf=&(A.colinfo);
 aind=&(A.colindex);
 aval=&(A.colval);
 const Indexmatrix* binf;
 const Indexmatrix* bind;
 const Matrix* bval;
 if (btrans) {
     nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nc) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
     binf=&(B.colinfo);
     bind=&(B.colindex);
     bval=&(B.colval);
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
     binf=&(B.rowinfo);
     bind=&(B.rowindex);
     bval=&(B.rowval);
 }
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.rowdim())||(nc!=C.coldim())) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 Integer anz=0;
 Integer brow1=0;
 for(Integer ai=0;ai<ainf->rowdim();ai++){
     Integer ii=(*ainf)(ai,0);
     if (ii<0){ //diagonal elements
         Integer brow2=0;
         for(Integer aj=(*ainf)(ai,1);--aj>=0;){
             ii=(*aind)(anz);
             // fins row ii in second matrix
             while((brow2<binf->rowdim())&&((*binf)(brow2,0)<ii)) brow2++;
             if (brow2==binf->rowdim()){
                 anz+=aj+1;
                 break;
             }
             if ((*binf)(brow2,0)>ii) {
                 anz++;
                 continue;
             }
             // add multiple of second to resulting matrix
             Integer bnz=(*binf)(brow2,2);
             Real d=alpha*(*aval)(anz++);
             for(Integer bj=(*binf)(brow2,1);--bj>=0;){
                 C(ii,(*bind)(bnz))+=d*(*bval)(bnz);
                 bnz++;
             }
         }
         
     }
     else { //offdiagonal elements
         
         //find row ii in second matrix
         while((brow1<binf->rowdim())&&((*binf)(brow1,0)<ii)) brow1++;
         int row_available=0;
         if((brow1<binf->rowdim())&&(ii==(*binf)(brow1,0))){
             row_available=1;
         }
         
         //multiply for each jj with pair (jj,ii) amd (jj,ii)
         Integer brow2=0;
         for(Integer aj=(*ainf)(ai,1);--aj>=0;){
             Integer jj=(*aind)(anz)+ii;
             Real d=alpha*(*aval)(anz++);
             //--- add row ii of second matrix (if it exists)
             //    to row jj of resulting matrix
             if (row_available){
                 Integer bnz=(*binf)(brow1,2);
                 for(Integer bj=(*binf)(brow1,1);--bj>=0;){
                     C(jj,(*bind)(bnz))+=d*(*bval)(bnz);
                     bnz++;
                 }
             }
             //--- add row jj of second matrix (if it exists)
             //    to row ii of resulting matrix
             while((brow2<binf->rowdim())&&((*binf)(brow2,0)<jj)) brow2++;
             if (brow2==binf->rowdim()){
                 if (row_available)
                     continue;
                 anz+=aj;
                 break;
             }
             if ((*binf)(brow2,0)>jj) {
                 continue;
             }
             Integer bnz=(*binf)(brow2,2);
             for(Integer bj=(*binf)(brow2,1);--bj>=0;){
                 C(ii,(*bind)(bnz))+=d*(*bval)(bnz);
                 bnz++;
             }
         }
     }
 }
 return C;
}


Matrix& genmult(const Sparsemat& A,const Sparsesym& B,Matrix &C,
                           Real alpha,Real beta,int atrans)
{
 chk_init(A);
 chk_init(B);
 Integer nr,nc;
 const Indexmatrix* ainf;
 const Indexmatrix* aind;
 const Matrix* aval;
#if (CONICBUNDLE_DEBUG>=1)
 Integer nm;
#endif
 if (atrans) {
     nr=A.nc;
#if (CONICBUNDLE_DEBUG>=1)
     nm=A.nr;
#endif
     ainf=&(A.rowinfo);
     aind=&(A.rowindex);
     aval=&(A.rowval);
 }
 else {
     nr=A.nr;
#if (CONICBUNDLE_DEBUG>=1)
     nm=A.nc;
#endif
     ainf=&(A.colinfo);
     aind=&(A.colindex);
     aval=&(A.colval);
 }

 const Indexmatrix* binf;
 const Indexmatrix* bind;
 const Matrix* bval;
 nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
 if (nm!=B.nr) {
     MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
 }
#endif
 binf=&(B.colinfo);
 bind=&(B.colindex);
 bval=&(B.colval);

 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.rowdim())||(nc!=C.coldim())) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 
 if (alpha==0.) return C;
 
 Integer bnz=0;
 Integer acol1=0;
 for(Integer bi=0;bi<binf->rowdim();bi++){
     Integer ii=(*binf)(bi,0);
     if (ii<0){ //diagonal elements
         Integer acol2=0;
         for(Integer bj=(*binf)(bi,1);--bj>=0;){
             ii=(*bind)(bnz);
             // find column ii in first matrix
             while((acol2<ainf->rowdim())&&((*ainf)(acol2,0)<ii)) acol2++;
             if (acol2==ainf->rowdim()){
                 bnz+=bj+1;
                 break;
             }
             if ((*ainf)(acol2,0)>ii) {
                 bnz++;
                 continue;
             }
             // add multiple of first to resulting matrix
             Integer anz=(*ainf)(acol2,2);
             Real d=alpha*(*bval)(bnz++);
             for(Integer aj=(*ainf)(acol2,1);--aj>=0;){
                 C((*aind)(anz),ii)+=d*(*aval)(anz);
                 anz++;
             }
         }
         
     }
     else { //offdiagonal elements
         
         //find column ii in first matrix
         while((acol1<ainf->rowdim())&&((*ainf)(acol1,0)<ii)) acol1++;
         int col_available=0;
         if((acol1<ainf->rowdim())&&(ii==(*ainf)(acol1,0))){
             col_available=1;
         }
         
         //multiply for each jj with pair (jj,ii) and (jj,ii)
         Integer acol2=0;
         for(Integer bj=(*binf)(bi,1);--bj>=0;){
             Integer jj=(*bind)(bnz)+ii;
             Real d=alpha*(*bval)(bnz++);
             //--- add column ii of first matrix (if it exists)
             //    to column jj of resulting matrix
             if (col_available){
                 Integer anz=(*ainf)(acol1,2);
                 for(Integer aj=(*ainf)(acol1,1);--aj>=0;){
                     C((*aind)(anz),jj)+=d*(*aval)(anz);
                     anz++;
                 }
             }
             //--- add column jj of first matrix (if it exists)
             //    to column ii of resulting matrix
             while((acol2<ainf->rowdim())&&((*ainf)(acol2,0)<jj)) acol2++;
             if (acol2==ainf->rowdim()){
                 if (col_available)
                     continue;
                 bnz+=bj;
                 break;
             }
             if ((*ainf)(acol2,0)>jj) {
                 continue;
             }
             Integer anz=(*ainf)(acol2,2);
             for(Integer aj=(*ainf)(acol2,1);--aj>=0;){
                 C((*aind)(anz),ii)+=d*(*aval)(anz);
                 anz++;
             }
         }
     }
 }
 return C;
}


Sparsesym& Sparsesym::xeya(const Sparsemat &A,Real d)
{
 chk_init(A);
 init(A.nr);
 chk_add(*this,A);
 if (d==0.) return *this;
 Integer nz=A.colval.dim();
 if (nz==0) return *this;
 Indexmatrix I(nz,1); chk_set_init(I,1);
 Indexmatrix J(nz,1); chk_set_init(J,1);
 Matrix val(nz,1); chk_set_init(val,1);
 nz=0;
 for(Integer i=0;i<A.rowinfo.rowdim();i++){
     Integer ii=A.rowinfo(i,0);
     for(Integer j=0;j<A.rowinfo(i,1);j++){
         Integer jj=A.rowindex(nz);
         if (ii==jj) {
             I(nz)=J(nz)=ii;
             val(nz)=d*A.rowval(nz);
         }
         else {
             I(nz)=ii;
             J(nz)=jj;
             val(nz)=d*A.rowval(nz)/2.;
         }
         nz++;
     }
 }

 return init(nr,nz,I,J,val); 
}

Sparsemat Sparsesym::sparsemult(const Matrix& A) const
{
 chk_mult(*this,A);
 Integer nc=A.coldim();
 Sparsemat B(nr,nc);  //B=(*this)*A

 Integer nrcol=colinfo.rowdim();
 if ((nrcol==0)||(nc==0)) return B;

 Integer suppnr=suppcol.rowdim();
 B.rowinfo.newsize(suppnr,3); chk_set_init(B.rowinfo,1);
 B.rowindex.newsize(suppnr*nc,1); chk_set_init(B.rowindex,1);
 B.rowval.newsize(suppnr*nc,1); chk_set_init(B.rowval,1);
 B.colinfo.newsize(nc,3); chk_set_init(B.colinfo,1);
 B.colindex.newsize(suppnr*nc,1); chk_set_init(B.colindex,1);
 B.colval.newsize(suppnr*nc,1); chk_set_init(B.colval,1);


 const Integer *jp=colinfo.get_store();
 const Integer *sjp=colinfo.get_store()+3*colinfo.rowdim();
 const Integer *cntp=colinfo.get_store()+colinfo.rowdim();
 const Integer *ip=colindex.get_store();
 const Integer *sp=suppind.get_store();
 const Real *vp=colval.get_store();
 const Real *abp=A.get_store();
 Real *bbp=B.rowval.get_store();
 
 if (nc>1){
     //initialize B
     Integer *bip=B.rowindex.get_store();
     Real *bvp=bbp;
     for(Integer i=0;i<suppnr;i++){
         B.rowinfo(i,0)=suppcol(i);
         B.rowinfo(i,1)=nc;
         B.rowinfo(i,2)=i*nc;
         for(Integer j=0;j<nc;j++){
             *bip++=j;
             *bvp++=0.;
         }
     }
             
     //diagonal part
     if (*jp<0) {
         for(Integer i=*cntp++;--i>=0;){
             mat_xeya(nc,bbp+(*sp++)*nc,1,abp+*ip++,nr,*vp++);
         }
         --nrcol;
         jp++;
         sjp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
         abp=A.get_store()+*jp++;
         bbp=B.rowval.get_store()+(*sjp++)*nc;
         for(Integer i=*cntp++;--i>=0;){
             Real d=*vp++;
             mat_xpeya(nc,bbp,1,abp+*ip++,nr,d);
             mat_xpeya(nc,bbp+(*sp++)*nc,1,abp,nr,d);
         }
     }

     //initialize colpart of B
     bip=B.colindex.get_store();
     bvp=B.colval.get_store();
     {for(Integer i=0;i<nc;i++){
         B.colinfo(i,0)=i;
         B.colinfo(i,1)=suppnr;
         B.colinfo(i,2)=i*suppnr;
         bbp=B.rowval.get_store()+i;
         sp=suppcol.get_store();
         for(Integer j=0;j<suppnr;j++){
             *bip++=*sp++;
             *bvp++=*bbp;
             bbp+=nc;
         }
     }}
             
 }
 else { //A is a vector (only one column)
     //initialize B
     Integer *bip=B.rowindex.get_store();
     Real *bvp=bbp;
     for(Integer i=0;i<suppnr;i++){
         B.rowinfo(i,0)=suppcol(i);
         B.rowinfo(i,1)=nc;
         B.rowinfo(i,2)=i*nc;
         *bip++=0;
         *bvp++=0.;
     }
             
     //diagonal part
     if (*jp<0) {
         for(Integer i=*cntp++;--i>=0;){
             *(bbp+(*sp++))=*(abp+*ip++)*(*vp++);
         }
         --nrcol;
         jp++;
         sjp++;
     }
     
     //offdiagonal part
     for(;--nrcol>=0;){ //for each column of the Sparsesym Matrix
         abp=A.get_store()+*jp++;
         bbp=B.rowval.get_store()+(*sjp++);
         for(Integer i=*cntp++;--i>=0;){
             Real d=*vp++;
             *bbp+=*(abp+*ip++)*d;
             *(bbp+(*sp++))+=*abp*d;
         }
     }

     //initialize colpart of B
     bip=B.colindex.get_store();
     bvp=B.colval.get_store();
     B.colinfo(0,0)=0;
     B.colinfo(0,1)=suppnr;
     B.colinfo(0,2)=0;
     bbp=B.rowval.get_store();
     sp=suppcol.get_store();
     for(Integer j=0;j<suppnr;j++){
         *bip++=*sp++;
         *bvp++=*bbp++;
     }             
 }

 return B;
}

Sparsemat& Sparsemat::xeya(const Sparsesym& A,Real in_d)
{
 chk_init(A);
 if ((in_d==0.)||(A.colinfo.rowdim()==0)) return init(A.nr,A.nr);
 Integer lnz=2*A.colval.dim();
 if (A.colinfo(0,0)<0) lnz-=A.colinfo(0,1);
 Indexmatrix I(lnz,1); chk_set_init(I,1);
 Indexmatrix J(lnz,1); chk_set_init(J,1);
 Matrix val(lnz,1); chk_set_init(val,1);
 Integer nz=0;
 lnz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz);
             I(lnz)=J(lnz)=ii;
             val(lnz++)=in_d*A.colval(nz++);
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz)+jj;
             Real d=in_d*A.colval(nz++);
             I(lnz)=ii; J(lnz)=jj; val(lnz++)=d;
             I(lnz)=jj; J(lnz)=ii; val(lnz++)=d;
         }
     }
 }
 return init(A.nr,A.nr,lnz,I,J,val);
}


// ******************************************************************************
//                    friends
// ******************************************************************************

void swap(Sparsesym& A, Sparsesym& B)
{
 Integer h=A.nr;A.nr=B.nr;B.nr=h;
 swap(A.colinfo,B.colinfo);
 swap(A.colindex,B.colindex);
 swap(A.colval,B.colval);
 swap(A.suppind,B.suppind);
 swap(A.suppcol,B.suppcol);
 Real d=A.tol;A.tol=B.tol;B.tol=d;
#if (CONICBUNDLE_DEBUG>=1)
 swap(A.is_init,B.is_init);
#endif
}
 
Matrix diag(const Sparsesym& A)
{
 chk_init(A);
 Matrix d(A.nr,1,0.);
 if ((A.colinfo.rowdim()==0)||(A.colinfo(0,0)>=0)) return d;
 const Real *ap=A.colval.get_store();
 const Integer *ip=A.colindex.get_store();
 for (Integer i=A.colinfo(0,1);--i>=0;)
     d(*ip++)=*ap++;
 return d;
}

Sparsesym sparseDiag(const Matrix &A,Real tol)
{
 chk_init(A);
 Sparsesym B;
 B.init(A.dim());

 //--- count nonzeros 
 Integer nz=0;
 for(Integer j=0;j<B.nr;j++){
     Real d=A(j);
     if (fabs(d)>tol) nz++;
 }
 if (nz==0) return B;

 //--- get memory and copy information
 B.colinfo.newsize(1,4); chk_set_init(B.colinfo,1);
 B.colinfo(0,0)=-1; B.colinfo(0,1)=nz; B.colinfo(0,2)=0;B.colinfo(0,3)=-1;
 B.colindex.newsize(nz,1); chk_set_init(B.colindex,1);
 B.colval.newsize(nz,1); chk_set_init(B.colval,1);
 B.suppind.newsize(nz,1); chk_set_init(B.suppind,1);
 B.suppcol.newsize(nz,1); chk_set_init(B.suppcol,1);
 nz=0;
 {for(Integer j=0;j<B.nr;j++){
      Real d=A(j);
      if (fabs(d)>tol) {
          B.colindex(nz)=j;
          B.colval(nz)=d;
          B.suppind(nz)=nz;
          B.suppcol(nz)=j;
          nz++;
      }
 }}
 return B;
}
 
Matrix sumrows(const Sparsesym& A)
{
 chk_init(A);
 Matrix s(1,A.nr,0.);
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;nz++){
             s(A.colindex(nz))+=A.colval(nz);
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;nz++){
             s(A.colindex(nz)+jj)+=A.colval(nz);
             s(jj)+=A.colval(nz);
         }
     }
 }
 return s;
}

Real sum(const Sparsesym& A)
{
 chk_init(A);
 if (A.colinfo.rowdim()==0) return 0.;
 Real s=0.;
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             s+=A.colval(nz++);
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             s+=2.*A.colval(nz++);
         }
     }
 }
 return s;
}

Real trace(const Sparsesym& A)
{
 chk_init(A);
 if ((A.colinfo.rowdim()==0)||(A.colinfo(0,0)>=0)) return 0.;
 Real s=0.;
 const Real *ap=A.colval.get_store();
 for (Integer i=A.colinfo(0,1);--i>=0;)
     s+=*ap++;
 return s;
}
    
Real ip(const Sparsesym& A, const Sparsesym& B)
{
 chk_add(A,B);
 Integer ai=0,bi=0;
 Real s=0.;
 while((ai<A.colinfo.rowdim())&&(bi<B.colinfo.rowdim())){
     if (A.colinfo(ai,0)==B.colinfo(bi,0)){
         Integer aj=A.colinfo(ai,2);
         Integer au=aj+A.colinfo(ai,1);
         Integer bj=B.colinfo(bi,2);
         Integer bu=bj+B.colinfo(bi,1);
         if (A.colinfo(ai,0)<0){
             while((aj<au)&&(bj<bu)){
                 if(A.colindex(aj)==B.colindex(bj)){
                     s+=A.colval(aj)*B.colval(bj);
                     aj++;
                     bj++;
                 }
                 else if (A.colindex(aj)<B.colindex(bj)){
                     aj++;
                 }
                 else {
                     bj++;
                 }
             }
         }
         else {
             while((aj<au)&&(bj<bu)){
                 if(A.colindex(aj)==B.colindex(bj)){
                     s+=2.*A.colval(aj)*B.colval(bj);
                     aj++;
                     bj++;
                 }
                 else if (A.colindex(aj)<B.colindex(bj)){
                     aj++;
                 }
                 else {
                     bj++;
                 }
             }
         }
         ai++;
         bi++;
     }
     else if (A.colinfo(ai,0)<B.colinfo(bi,0)){
         ai++;
     }
     else {
         bi++;
     }
 }
 return s;
}

Real ip(const Matrix& B,const Sparsesym& A)
{
 chk_add(B,A);
 Real s=0;
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz);
             s+=B(ii,ii)*A.colval(nz++);
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz)+jj;
             s+=(B(ii,jj)+B(jj,ii))*A.colval(nz++);
         }
     }
 }
 return s;
}

Real norm2(const Sparsesym& A)
{
 chk_init(A);
 if (A.colinfo.rowdim()==0) return 0.;
 Real s=0.;
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             s+=sqr(A.colval(nz++));
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             s+=2.*sqr(A.colval(nz++));
         }
     }
 }
 return ::sqrt(s);
}

Sparsesym abs(const Sparsesym& A)
{
 chk_init(A);
 Sparsesym B; B.init(A.nr);
 B.colinfo=A.colinfo;
 B.colindex=A.colindex;
 B.suppind=A.suppind;
 B.suppcol=A.suppcol;
 B.colval.newsize(A.colval.dim(),1); chk_set_init(B.colval,1);
 Real *bc=B.colval.get_store();
 const Real *ac=A.colval.get_store();
 Integer i=A.colval.dim();
 for(;--i>=0;){
     *bc++= fabs(*ac++);
 }
 return B;
}

int equal(const Sparsesym& A, const Sparsesym& B,Real eqtol)
{
  chk_init(A);
  chk_init(B);
  if ((A.nr!=B.nr)||
      (A.colinfo.dim()!=B.colinfo.dim())||
      (A.colval.dim()!=B.colval.dim())) return 0;
  const Integer *ai=A.colinfo.get_store();
  const Integer *bi=B.colinfo.get_store();
  for (Integer i=B.colinfo.dim();--i>=0;){
    if (*ai++!=*bi++) return 0;
  }
  ai=A.colindex.get_store();
  bi=B.colindex.get_store();
  const Real *ad=A.colval.get_store();
  const Real *bd=B.colval.get_store();
  {for (Integer i=B.colval.dim();--i>=0;){
    if ((*ai++!=*bi++)||(fabs(*ad++-*bd++))>eqtol) 
      return 0;
  }}
  return 1;
}

std::ostream& operator<<(std::ostream& out,const Sparsesym &A)
{
 chk_init(A);
 out<<A.nr<<" "<<A.colval.dim()<<"\n";
 Integer nz=0;
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     if (jj<0){ //diagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz);
             out<<ii<<" ";
             out<<ii<<" ";
             out<<A.colval(nz++)<<"\n";
         }
     }
     else { //offdiagonal case
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Integer ii=A.colindex(nz)+jj;
             out<<jj<<" ";
             out<<ii<<" ";
             out<<A.colval(nz++)<<"\n";
        }
     }
 }
 return out;
}

std::istream& operator>>(std::istream& in,Sparsesym &A)
{
 const char* format="          format: nr nz i_1 j_1 val_1 ... i_nz j_nz val_nz\n          nr>=0, nz>=0, 0<=i<nr, 0<=j<nr, upper or lower triangle\n          identical elements are added\n"; 
 Integer nr,nz;
 if (!(in>>nr)){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
     if (materrout) (*materrout)<<" failed in reading number of rows"<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if(nr<0){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
     if (materrout) (*materrout)<<" number of rows must be nonnegative but is "<<nr<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if (!(in>>nz)){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
     if (materrout) (*materrout)<<" failed in reading number of nonzero elements"<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if(nz<0){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
     if (materrout) (*materrout)<<" number of nonzeros must be nonnegative but is "<<nz<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if((nr==0)&&(nz>0)){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
     if (materrout) (*materrout)<<" zero row or column dimension but positive number of nonzeros"<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 Integer i;
 Indexmatrix I(nz,1),J(nz,1);
 Matrix val(nz,1);
 chk_set_init(I,1);
 chk_set_init(J,1);
 chk_set_init(val,1);
 for(i=0;i<nz;i++){
     if (!(in>>I(i)>>J(i)>>val(i))){
         if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
         if (materrout) (*materrout)<<" failed in reading nonzero element (i,j,val) #";
         if (materrout) (*materrout)<<i<<std::endl;
         in.clear(std::ios::failbit);
         return in;
     }
     if ((0>I(i))||(I(i)>=nr)){
         if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
         if (materrout) (*materrout)<<" row index of nonzero element #"<<i<<"exceeds range: ";
         if (materrout) (*materrout)<<0<<"<="<<I(i)<<"<"<<nr<<std::endl;
         in.clear(std::ios::failbit);
         return in;
     }
     if ((0>J(i))||(J(i)>=nr)){
         if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsesym&): ";
         if (materrout) (*materrout)<<" column index of nonzero element #"<<i<<"exceeds range: ";
         if (materrout) (*materrout)<<0<<"<="<<J(i)<<"<"<<nr<<std::endl;
         in.clear(std::ios::failbit);
         return in;
     }
 }
 A.init(nr,nz,I,J,val);
 return in;
}

}
