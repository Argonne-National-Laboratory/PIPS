/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/symmat.cxx

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

const Mtype Symmatrix::mtype = MTsymmetric;

// **************************************************************************
//                                Constructors
// **************************************************************************

Symmatrix& Symmatrix::xeya(const Symmatrix& A,Real d)
{
 chk_init(A);
 newsize(A.nr);
 chk_set_init(*this,1);
 if (d==1.) { mat_xey((nr*(nr+1))/2,m,A.get_store()); return *this;}
 if (d==0.) { mat_xea((nr*(nr+1))/2,m,0.); return *this;}
 if (d==-1.) { mat_xemy((nr*(nr+1))/2,m,A.get_store()); return *this;}
 mat_xeya((nr*(nr+1))/2,m,A.get_store(),d);
 return *this;
}

Symmatrix& Symmatrix::xpeya(const Symmatrix& A,Real d)
{
 chk_add(*this,A);
 if (d==1.) { mat_xpey((nr*(nr+1))/2,m,A.get_store()); return *this;}
 if (d==0.) { return *this;}
 if (d==-1.) { mat_xmey((nr*(nr+1))/2,m,A.get_store()); return *this;}
 mat_xpeya((nr*(nr+1))/2,m,A.get_store(),d);
 return *this;
}


// returns C=beta*C+alpha* A*A^T, where A may be transposed
Symmatrix& rankadd(const Matrix& A,Symmatrix& C,Real alpha,Real beta,int trans)
#ifdef WITH_BLAS
{
 chk_init(A);
 Integer nr=(trans?A.coldim():A.rowdim());
 Integer nk=(trans?A.rowdim():A.coldim());
#if (CONICBUNDLE_DEBUG>=1)
 if (beta!=0.){
   chk_init(C);
   if (C.nr!=nr) {
     MEmessage(MatrixError(ME_dim,"rankadd: dimensions don't match",MTmatrix));;
   }
 }
#endif
 if ((nr==0)||(nk==0)) {
   if (beta==0.) return C.init(nr,0.);
   return C*=beta;
 }
 if (beta==0.) {
   C.newsize(nr);
   chk_set_init(C,1);
 }
 Matrix tmp(nr,nr); chk_set_init(tmp,1);
 Real *tp=tmp.get_store();
 Real *cp=C.get_store();
 for(Integer i=nr;i>0;--i){
   mat_xey(i,tp,cp);
   cp+=i;
   tp+=nr+1;
 }
 cblas_dsyrk(CblasColMajor,CblasLower,(trans?CblasTrans:CblasNoTrans),
	     nr,nk,alpha,A.get_store(),A.rowdim(),beta,tmp.get_store(),tmp.rowdim());
 tp=tmp.get_store();
 cp=C.get_store();
 for(Integer i=nr;i>0;--i){
   mat_xey(i,cp,tp);
   cp+=i;
   tp+=nr+1;
 }
 return C;
}
#else
{
 chk_init(A);
 if (trans){ //A^T*A
     if (beta==0.) {
         C.init(A.coldim(),0.);
     }
     else {
         chk_colseq(C,A);
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer nr=C.rowdim();
     Integer n=A.rowdim();
     Real *mp=C.get_store();
     const Real *ap=A.get_store()-n;
     for(Integer i=nr;--i>=0;){ 
         const Real *aap=(ap+=n);
         (*mp++) +=alpha*mat_ip(n,ap,aap);
         for(Integer j=i;--j>=0;){
             (*mp++) +=alpha*mat_ip(n,ap,aap+=n);
         }
     }
     return C;
 }
 else { //A*A^T
     if (beta==0.) {
         C.init(A.rowdim(),0.);
     }
     else {
         chk_rowseq(C,A);
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer nr=C.rowdim();
     Integer nc2=(nr*(nr+1))/2;
     Integer n=A.coldim();
     const Real *ap=A.get_store()-1;
     Real *mp=C.get_store();
     for(Integer i=n;--i>=0;){
         for(Integer j=nr; --j>=0; ){
             const Real *aap;
             Real d1=*(aap=++ap);
             Real d=d1*alpha;
             (*mp++)+=d*d1;
             for(Integer k=j;--k>=0;)
                 (*mp++)+=d*(*(++aap));
         }
         mp-=nc2;
     }
     return C;
 }
}
#endif
     
// returns C=beta*C+alpha* sym(A*B^T) [or sym(A^T*B)]
Symmatrix& rank2add(const Matrix& A, const Matrix& B, Symmatrix& C,
                              Real alpha,Real beta,int trans)
#ifdef WITH_BLAS
{
 chk_add(A,B);
 Integer nr=(trans?A.coldim():A.rowdim());
 Integer nk=(trans?A.rowdim():A.coldim());
#if (CONICBUNDLE_DEBUG>=1)
 if (beta!=0.){
   chk_init(C);
   if (C.nr!=nr) {
     MEmessage(MatrixError(ME_dim,"rank2add: dimensions don't match",MTmatrix));;
   }
 }
#endif
 if ((nr==0)||(nk==0)) {
   if (beta==0.) return C.init(nr,0.);
   return C*=beta;
 }
 if (beta==0.) {
   C.newsize(nr);
   chk_set_init(C,1);
 }
 Matrix tmp(nr,nr); chk_set_init(tmp,1);
 Real *tp=tmp.get_store();
 Real *cp=C.get_store();
 for(Integer i=nr;i>0;--i){
   mat_xey(i,tp,cp);
   cp+=i;
   tp+=nr+1;
 }
 cblas_dsyr2k(CblasColMajor,CblasLower,(trans?CblasTrans:CblasNoTrans),
	      nr,nk,alpha/2.,A.get_store(),A.rowdim(),B.get_store(),B.rowdim(),beta,
	      tmp.get_store(),tmp.rowdim());
 tp=tmp.get_store();
 cp=C.get_store();
 for(Integer i=nr;i>0;--i){
   mat_xey(i,cp,tp);
   cp+=i;
   tp+=nr+1;
 }
 return C;
}
#else
{
 chk_add(A,B);
 if (trans) { //A^T*B
     if (beta==0.) {
         C.init(A.coldim(),0.);
     }
     else {
         chk_colseq(C,A);
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer nr=C.rowdim();
     Integer n=A.rowdim();
     Real *mp=C.get_store();
     const Real *abp=A.get_store()-n;
     const Real *bbp=B.get_store()-n;
     Real d=alpha/2.;
     for(Integer i=nr;--i>=0;){ 
         const Real *bp;
         const Real *ap;
         *mp++ += alpha*mat_ip(n,ap=(abp+=n),bp=(bbp+=n));
         for(Integer j=i;--j>=0;){
             (*mp++) +=d*(mat_ip(n,abp,bp+=n)+mat_ip(n,ap+=n,bbp));
         }       
     }
     return C;
 }
 else {     //A*B^T
     if (beta==0.) {
         C.init(A.rowdim(),0.);
     }
     else {
         chk_rowseq(C,A);
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer nr=C.rowdim();
     Integer nc2=(nr*(nr+1))/2;
     Integer n=A.coldim();
     const Real *abp=A.get_store()-1;
     const Real *bbp=B.get_store()-1;
     Real *mp=C.get_store();
     for(Integer i=n;--i>=0;){
         for(Integer j=nr; --j>=0; ){
             const Real *ap;
             const Real *bp;
             Real ad=*(ap=++abp);
             Real bd=alpha**(bp=++bbp);
             (*mp++)+=ad*bd;
             ad*=alpha/2.;
             bd/=2.;
             for(Integer k=j;--k>=0;)
                 (*mp++)+=bd*(*(++ap))+ad*(*(++bp));
         }
         mp-=nc2;
     }
     return C;
 }
}     
#endif
     

Symmatrix& Symmatrix::xetriu_yza(const Matrix& A,const Matrix& B,Real d)
{
 chk_add(A,B);
 newsize(A.coldim());
 chk_set_init(*this,1);
 if (d==0.) { mat_xea((nr*(nr+1))/2,m,0.); return *this;}
 if (d==1.) {
     Real *mp=m;
     const Real *ap=A.get_store();
     Integer n=A.rowdim();
     const Real *bbp=B.get_store()-n;
     for(Integer i=nr;i>0;i--){ 
         const Real *bp=(bbp+=n);
         for(Integer j=i;--j>=0;){
             *mp++=mat_ip(n,ap,bp);
             bp+=n;
         }
         ap+=n;         
     }
     return *this;
 }
 if (d==-1.) {
     Real *mp=m;
     const Real *ap=A.get_store();
     Integer n=A.rowdim();
     const Real *bbp=B.get_store()-n;
     for(Integer i=nr;i>0;i--){ 
         const Real *bp=(bbp+=n);
         for(Integer j=i;--j>=0;){
             *mp++=-mat_ip(n,ap,bp);
             bp+=n;
         }
         ap+=n;         
     }
     return *this;
 }
 Real *mp=m;
 const Real *ap=A.get_store();
 Integer n=A.rowdim();
 const Real *bbp=B.get_store()-n;
 for(Integer i=nr;i>0;i--){ 
     const Real *bp=(bbp+=n);
     for(Integer j=i;--j>=0;){
         *mp++=d*mat_ip(n,ap,bp);
         bp+=n;
     }
     ap+=n;         
 }
 return *this;
}
         
Symmatrix& Symmatrix::xpetriu_yza(const Matrix& A,const Matrix& B,Real d)
{
 chk_add(A,B);
 chk_colseq(*this,A);
 if (d==0.) { return *this;}
 if (d==1.) {
     Real *mp=m;
     const Real *ap=A.get_store();
     Integer n=A.rowdim();
     for(Integer i=0;i<nr;i++){
         const Real *bp=B.get_store()+i*n;
         for(Integer j=i;j<nr;j++){
             (*mp++) += mat_ip(n,ap,bp);
             bp+=n;
         }
         ap+=n;
     }
     return *this;
 }
 if (d==-1.) {
     Real *mp=m;
     const Real *ap=A.get_store();
     Integer n=A.rowdim();
     for(Integer i=0;i<nr;i++){
         const Real *bp=B.get_store()+i*n;
         for(Integer j=i;j<nr;j++){
             (*mp++)-=mat_ip(n,ap,bp);
             bp+=n;
         }
         ap+=n;
     }
     return *this;
 }
 Real *mp=m;
 const Real *ap=A.get_store();
 Integer n=A.rowdim();
 for(Integer i=0;i<nr;i++){
     const Real *bp=B.get_store()+i*n;
     for(Integer j=i;j<nr;j++){
         (*mp++)+=d*mat_ip(n,ap,bp);
         bp+=n;
     }
     ap+=n;
 }
 return *this;
}

Matrix& genmult(const Symmatrix& A,const Matrix& B,Matrix& C,
                    Real alpha,Real beta,int btrans)
            //returns C=beta*C+alpha*A*B, where A and B may be transposed
#ifdef WITH_BLAS
{
 chk_init(A);
 chk_init(B);
 Integer nr=A.rowdim();
 Integer nc=(btrans?B.rowdim():B.coldim());
#if (CONICBUNDLE_DEBUG>=1)
 if (beta!=0.) chk_init(C);
 if ((nr!=(btrans?B.coldim():B.rowdim()))||
     ((beta!=0.)&&((C.rowdim()!=nr)||(C.coldim()!=nc)))
     ){
   MEmessage(MatrixError(ME_dim,"genmult(Symmatrix,Matrix,Matrix): dimensions don't match",MTmatrix));;
 }
#endif
 if ((nr==0)||(nc==0)) {
   if (beta==0.) return C.init(nr,nc,0.);
   return C*=beta;
 }
 if (beta==0.) {
   C.newsize(nr,nc);
   chk_set_init(C,1);
 }
 Matrix tmp(nr,nr); chk_set_init(tmp,1);
 Real *tp=tmp.get_store();
 const Real *ap=A.get_store();
 for(Integer i=nr;i>0;--i){
   mat_xey(i,tp,ap);
   ap+=i;
   tp+=nr+1;
 }
 if (btrans==0){
   cblas_dsymm(CblasColMajor,CblasLeft,CblasLower,nr,nc,alpha,
	       tmp.get_store(),tmp.rowdim(),B.get_store(),B.rowdim(),beta,
	       C.get_store(),C.rowdim());
 }
 else {
   const Integer br=B.rowdim();
   const Integer bc=B.coldim();
   Matrix tmpB(bc,br); chk_set_init(tmpB,1);
   Real *tp=tmpB.get_store();
   const Real *bp=B.get_store();
   for(Integer i=0;i<br;i++){
     mat_xey(br,tp,bc,bp,1);
     bp+=br;
     tp++;
   }
   cblas_dsymm(CblasColMajor,CblasLeft,CblasLower,nr,nc,alpha,
	       tmp.get_store(),tmp.rowdim(),tmpB.get_store(),tmpB.rowdim(),beta,
	       C.get_store(),C.rowdim());
 }
 return C;
}
#else
{
 chk_init(A);
 chk_init(B);
 Integer nr,nm,nc;
 nr=A.nr;nm=A.nr;
 if (btrans) {
     nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nc) {
         MEmessage(MatrixError(ME_dim,"genmult(Symmatrix,Matrix,Matrix): dimensions don't match",MTsymmetric));;
     }
#endif
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult(Symmatrix,Matrix,Matrix): dimensions don't match",MTsymmetric));;
     }
#endif
 }
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.nr)||(nc!=C.nc)) {
         MEmessage(MatrixError(ME_dim,"genmult(Symmatrix,Matrix,Matrix): dimensions don't match",MTsymmetric));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if ((btrans)&&(nc>1)){
   /*
     const Real *ap=A.get_store();
     for(Integer j=0;j<nm;j++){
         mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+j*B.nc,1,alpha*(*ap++));
         for(Integer i=j+1;i<nr;i++){
             Real d=alpha*(*ap++);
             mat_xpeya(nc,C.get_store()+i,nr,B.get_store()+j*B.nc,1,d);
             mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+i*B.nc,1,d);
         }
     }
   */
   Real* cp=C.m;
   for(Integer col=0;col<nc;col++){
     const Real* bp=B.m+col; 
     const Real* ap=A.m; //the symmetric store
     for(Integer i=nr;--i>=0;){
       Real bi=(*bp);
       Real dsum=alpha*bi*(*ap++);
       Real* cjp=cp+1;
       for(Integer j=i;--j>=0;){
	 Real aij=alpha*(*ap++);
	 (*cjp++) += bi*aij;
	 bp+=nc;
	 dsum +=aij*(*bp);
       }
       (*cp++) += dsum;
       bp-=(i-1)*nc;
     }
   }
 }
 else {
   /*
     const Real *ap=A.get_store();
     for(Integer j=0;j<nm;j++){
         mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+j,B.nr,alpha*(*ap++));
         for(Integer i=j+1;i<nr;i++){
             Real d=alpha*(*ap++);
             mat_xpeya(nc,C.get_store()+i,nr,B.get_store()+j,B.nr,d);
             mat_xpeya(nc,C.get_store()+j,nr,B.get_store()+i,B.nr,d);
         }
     }
   */
   const Real* bp=B.m; 
   Real* cp=C.m;
   for(Integer col=nc;--col>=0;){
     const Real* ap=A.m; //the symmetric store
     for(Integer i=nr;--i>=0;){
       Real bi=(*bp++);
       Real dsum=alpha*bi*(*ap++);
       Real* cjp=cp+1;
       for(Integer j=i;--j>=0;){
	 Real aij=alpha*(*ap++);
	 (*cjp++) += bi*aij;
	 dsum +=aij*(*bp++);
       }
       (*cp++) += dsum;
       bp-=i;
     }
   }
 }
 return C;
}
#endif

Matrix& genmult(const Matrix& A,const Symmatrix& B,Matrix& C,
                    Real alpha,Real beta,int atrans)
            //returns C=beta*C+alpha*A*B, where A and B may be transposed
#ifdef WITH_BLAS
{
 chk_init(B);
 chk_init(A);
 Integer nr=(atrans?A.coldim():A.rowdim());
 Integer nc=B.rowdim();
#if (CONICBUNDLE_DEBUG>=1)
 if (beta!=0.) chk_init(C);
 if ((nc!=(atrans?A.rowdim():A.coldim()))||
     ((beta!=0.)&&((C.rowdim()!=nr)||(C.coldim()!=nc)))
     ){
   MEmessage(MatrixError(ME_dim,"genmult(Symmatrix,Matrix,Matrix): dimensions don't match",MTmatrix));;
 }
#endif
 if ((nr==0)||(nc==0)) {
   if (beta==0.) return C.init(nr,nc,0.);
   return C*=beta;
 }
 if (beta==0.) {
   C.newsize(nr,nc);
   chk_set_init(C,1);
 }
 Matrix tmp(nc,nc); chk_set_init(tmp,1);
 Real *tp=tmp.get_store();
 const Real *bp=B.get_store();
 for(Integer i=nc;i>0;--i){
   mat_xey(i,tp,bp);
   bp+=i;
   tp+=nc+1;
 }
 if (atrans==0){
   cblas_dsymm(CblasColMajor,CblasRight,CblasLower,nr,nc,alpha,
	       tmp.get_store(),tmp.rowdim(),A.get_store(),A.rowdim(),beta,
	       C.get_store(),C.rowdim());
 }
 else {
   const Integer ar=A.rowdim();
   const Integer ac=A.coldim();
   Matrix tmpA(ac,ar); chk_set_init(tmpA,1);
   Real *tp=tmpA.get_store();
   const Real *ap=A.get_store();
   for(Integer i=0;i<ar;i++){
     mat_xey(ar,tp,ac,ap,1);
     ap+=ar;
     tp++;
   }
   cblas_dsymm(CblasColMajor,CblasRight,CblasLower,nr,nc,alpha,
	       tmp.get_store(),tmp.rowdim(),tmpA.get_store(),tmpA.rowdim(),beta,
	       C.get_store(),C.rowdim());
 }
 return C;
}
#else
{
 chk_init(A);
 chk_init(B);
 Integer nr,nm,nc;
 if (atrans){
     nr=A.nc;
     nm=A.nr;
 }
 else {
     nr=A.nr;
     nm=A.nc;
 }
 nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
 if (nm!=B.nr) {
     MEmessage(MatrixError(ME_dim,"genmult(Symmatrix,Matrix,Matrix): dimensions don't match",MTsymmetric));;
 }
#endif
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.nr)||(nc!=C.nc)) {
         MEmessage(MatrixError(ME_dim,"genmult(Symmatrix,Matrix,Matrix): dimensions don't match",MTsymmetric));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if (atrans){
     const Real *bp=B.get_store();
     for(Integer j=0;j<nc;j++){
         mat_xpeya(nr,C.get_store()+j*nr,1,A.get_store()+j,nm,alpha*(*bp++));
         for(Integer i=j+1;i<nm;i++){
             Real d=alpha*(*bp++);
             mat_xpeya(nr,C.get_store()+j*nr,1,A.get_store()+i,nm,d);
             mat_xpeya(nr,C.get_store()+i*nr,1,A.get_store()+j,nm,d);
         }
     }
 }
 else {
     const Real *bp=B.get_store();
     for(Integer j=0;j<nc;j++){
         mat_xpeya(nr,C.get_store()+j*nr,A.get_store()+j*nr,alpha*(*bp++));
         for(Integer i=j+1;i<nm;i++){
             Real d=alpha*(*bp++);
             mat_xpeya(nr,C.get_store()+j*nr,A.get_store()+i*nr,d);
             mat_xpeya(nr,C.get_store()+i*nr,A.get_store()+j*nr,d);
         }
     }
 }
 return C;
}
#endif

Symmatrix& Symmatrix::xeya(const Matrix &A,double d)
{
 chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr!=A.nc)
     MEmessage(MEdim(A.nr,A.nc,0,0,"Symmatrix::Symmatrix(const Matrix& A) A not square",MTsymmetric));
 is_init=0;
#endif
 newsize(A.nr);
 Real *mp=m;
 Real f=d/2.;
 for(Integer i=0;i<nr;i++){
     Real* matcp=A.m+i*nr+i;
     Real* matrp=matcp+nr;
     (*mp++)=(*matcp++)*d;
     for(Integer j=i+1;j<nr;j++,matrp+=nr)
         (*mp++)=((*matcp++)+(*matrp))*f;
 }
 chk_set_init(*this,1);
 return *this;
}

Symmatrix& Symmatrix::xpeya(const Matrix &A,double d)
{
 chk_add(*this,A);
 Real *mp=m;
 Real f=d/2.;
 for(Integer i=0;i<nr;i++){
     Real* matcp=A.m+i*nr+i;
     Real* matrp=matcp+nr;
     (*mp++)+=(*matcp++)*d;
     for(Integer j=i+1;j<nr;j++,matrp+=nr)
         (*mp++)+=((*matcp++)+(*matrp))*f;
 }
 return *this;
}

Symmatrix& Symmatrix::xeya(const Indexmatrix &A,double d)
{
 chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr!=A.nc)
     MEmessage(MEdim(A.nr,A.nc,0,0,"Symmatrix::Symmatrix(const Matrix& A) A not square",MTsymmetric));
 is_init=0;
#endif
 newsize(A.nr);
 Real *mp=m;
 Real f=d/2.;
 for(Integer i=0;i<nr;i++){
     Integer* matcp=A.m+i*nr+i;
     Integer* matrp=matcp+nr;
     (*mp++)=(*matcp++)*d;
     for(Integer j=i+1;j<nr;j++,matrp+=nr)
         (*mp++)=((*matcp++)+(*matrp))*f;
 }
 chk_set_init(*this,1);
 return *this;
}

Symmatrix& Symmatrix::xpeya(const Indexmatrix &A,double d)
{
 chk_add(*this,A);
 Real *mp=m;
 Real f=d/2.;
 for(Integer i=0;i<nr;i++){
     Integer* matcp=A.m+i*nr+i;
     Integer* matrp=matcp+nr;
     (*mp++)+=(*matcp++)*d;
     for(Integer j=i+1;j<nr;j++,matrp+=nr)
         (*mp++)+=((*matcp++)+(*matrp))*f;
 }
 return *this;
}

void Symmatrix::newsize(Integer inr)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (inr<0)
     MEmessage(MEdim(inr,inr,0,0,"Symmatrix::newsize(Integer,Integer) dim<0",MTsymmetric));
 is_init=0;
#endif
 if (inr!=nr) {
     nr=inr;
     Integer n2=(nr*(nr+1))/2;
     if (n2>mem_dim){
         memarray->free(m); m=0;
         mem_dim=Integer(memarray->get(n2,m));
         if (mem_dim<n2)
             MEmessage(MEmem(n2,
                         "Symmatrix::Symmatrix(Integer,Integer,Real) not enough memory",MTsymmetric));
     }
 }
}

Symmatrix& Symmatrix::shift_diag(Real s)
{
  chk_init(*this);
  Real* mp=m;
  for (Integer i=0;i<nr;i++){
    *mp += s;
    mp += nr-i;
  }
  return *this;
}

void Symmatrix::display(std::ostream& out,int precision,int width,
                        int screenwidth) const
{
 chk_init(*this);
 out<<"Symmatrix("<<nr<<")"<<std::endl;
 if (nr==0) return;
 if (precision==0) precision=4;
 out.precision(precision);
 if (width==0) width=precision+6;
 if (screenwidth==0) screenwidth=80;
 Integer colnr=screenwidth/(width+1);
 Integer k,i,j;
 Integer maxk=nr/colnr+((nr%colnr)>0);
 Integer maxj;
 for(k=0;k<maxk;k++){
     out<<"columns "<<k*colnr<<" to "<<min(nr,(k+1)*colnr)-1<<std::endl;
     for(i=0;i<nr;i++){
         maxj=min((k+1)*colnr,nr);
         for(j=k*colnr;j<maxj;j++){
             out<<' ';out.width(width);out<<(*this)(i,j);
         }
         out<<std::endl;
     }     
 }
}

Matrix Symmatrix::col(Integer c) const
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
 if ((c<0)||(c>=nr))
     MEmessage(MErange(nr,0,nr,c,"Symmatrix::col(Integer) index out of range",MTsymmetric));
#endif
 Matrix v(nr,1);
 Real *mp=m+c;
 Real *vp=v.m;
 Integer i=0;
 for(;i<c;i++){
     (*vp++)=*mp;
     mp+=nr-i-1;
 }
 for(;i<nr;i++) (*vp++)=*mp++;
 chk_set_init(v,1);
 return v;
}

Matrix Symmatrix::row(Integer r) const
{
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
 if ((r<0)||(r>=nr))
     MEmessage(MErange(nr,r,nr,0,"Symmatrix::row(Integer) index out of range",MTsymmetric));
#endif
 Matrix v(1,nr);
 Real *mp=m+r;
 Real *vp=v.m;
 Integer i=0;
 for(;i<r;i++){
     (*vp++)=*mp;
     mp+=nr-i-1;
 }
 for(;i<nr;i++) (*vp++)=*mp++;
 chk_set_init(v,1);
 return v;
}

// *****************************************************************************
//                  Interaction with other classes
// *****************************************************************************


// **************************************************************************
//                               friends
// **************************************************************************

void svec(const Symmatrix &A,Matrix& sv)
{
 chk_init(A);
 const Real sqrt2=::sqrt(2.); //1.41421356237310;
 sv.newsize((A.nr*(A.nr+1))/2,1);
 Real *mp=A.m;
 Real *vp=sv.get_store();
 Integer i,j;
 for(i=A.nr;--i>=0;){
     (*vp++)=(*mp++);
     for(j=i;--j>=0;){
         (*vp++)=sqrt2*(*mp++);
     }
 }
 chk_set_init(sv,1);
}

void sveci(const Matrix& sv,Symmatrix &A)
{
 chk_init(sv);
 Integer dim=Integer(::sqrt(Real(8*sv.dim()+1))-1+.1)/2;
#if (CONICBUNDLE_DEBUG>=1)
 if ((dim*(dim+1))/2!=sv.dim())
     MEmessage(MatrixError(ME_unspec,
               "sveci(): dimension of svec does not permit conversion",
               MTsymmetric));
#endif                            
 const Real sqrt2=1./::sqrt(2.); //1./1.41421356237310;
 A.newsize(dim);
 const Real *mp=sv.get_store();
 Real *vp=A.m;
 Integer i,j;
 for(i=A.nr;--i>=0;){
     (*vp++)=(*mp++);
     for(j=i;--j>=0;){
         (*vp++)=sqrt2*(*mp++);
     }
 }
 chk_set_init(A,1);
}

void skron(const Symmatrix& A,const Symmatrix& B,Symmatrix& S)
{
 chk_add(A,B);
 const Real sqrt2=::sqrt(2.); // 1.41421356237310;
 Integer n=A.rowdim(); 
 S.newsize((n*(n+1))/2);
 chk_set_init(S,1);
 Matrix sv((n*(n+1))/2,1);
 chk_set_init(sv,1);
 Real* Sp=S.get_store();
 Integer m=(n*(n+1))/2;
 Matrix Ah(A);
 Matrix Bh(B);
 for(Integer i=0;i<n;i++){

     //row corresponding to diagonal
     {
       const Real* ai=Ah.get_store()+i*n+i; //points now to i-th elem of i-th column of A
       const Real* bi=Bh.get_store()+i*n+i; //points now to i-th elem of i-thcolumn of B
       Real *svp=sv.get_store();
       for(Integer l=i;l<n;l++){
	 const Real ail=*ai++; //ai points afterwards to the element A(l+1,i)
	 const Real bil=*bi++; //bi points afterwards to the element B(l+1,i)
	 const Real* aik=ai;  
	 const Real* bik=bi;  
         (*svp++)=ail * bil;          
         for(Integer k=l+1;k<n;k++){
	   (*svp++)=sqrt2*(ail *(*bik++)+(*aik++)*bil)/2.;
         }
       }
     }
     mat_xey(m,Sp,sv.get_store());
     Sp+=m;
     m--;
     
     for(Integer j=i+1;j<n;j++){
       const Real* ai=Ah.get_store()+i*n+i; //points now to i-th elem in i-th column of A
       const Real* bi=Bh.get_store()+i*n+i; //points now to i-th elem in i-th column of B
       const Real* aj=Ah.get_store()+j*n+i; //points now to i-th elem in j-th column of A
       const Real* bj=Bh.get_store()+j*n+i; //points now to i-th elem in j-th column of B

       //row correpsonding to offdiagonal
       Real *svp=sv.get_store();
       {
	 const Real* aik=ai+j-i; //points now to A(j,i)
	 const Real* bik=bi+j-i; //points now to B(j,i)
	 const Real* ajk=aj+j-i; //points now to A(j,j)
	 const Real* bjk=bj+j-i; //points now to B(j,j)
	 const Real aii=*ai++; //afterwards ai points to A(i+1,i)
	 const Real bii=*bi++; //afterwards bi points to B(i+1,i)
	 const Real aji=*aj++; //afterwards aj points now to A(i+1,j)
	 const Real bji=*bj++; //afterwards bj points now to B(i+1,j)
	 for(Integer k=j;k<n;k++){
	   (*svp++)=(aii*(*bjk++)+(*aik++)*bji+aji*(*bik++)+(*ajk++)*bii)/2.;
	 }
       }
       for(Integer l=i+1;l<n;l++){
	 const Real ail=*ai++; //afterwards ai points to A(l+1,i)
	 const Real bil=*bi++; //afterwards bi points to B(l+1,i)
	 const Real ajl=*aj++; //afterwards aj points to A(l+1,j)
	 const Real bjl=*bj++; //afterwards bj points to B(l+1,j)
	 (*svp++)=sqrt2*(ail*bjl+ajl*bil)/2.;
	 const Real* aik=ai; //aik points now to A(l+1,i)
	 const Real* bik=bi; //bik points now to B(l+1,i)
	 const Real* ajk=aj; //ajk points now to A(l+1,j)
	 const Real* bjk=bj; //bjk points now to B(l+1,j)
	 for(Integer k=l+1;k<n;k++){
	   (*svp++)=(ail*(*bjk++)+(*aik++)*bjl+ajl*(*bik++)+(*ajk++)*bil)/2.;
	 }
       }
       mat_xey(m,Sp,sv.get_store());
       Sp+=m;
       m--;
     }
 }
 
}

Matrix diag(const Symmatrix& A)
{
  chk_init(A);
  Matrix v(A.nr,1);
  for(Integer i=0;i<A.nr;i++)
    v(i)=A(i,i);
  chk_set_init(v,1);
  return v;
}

Symmatrix Diag(const Matrix& A)
{
  chk_init(A);
  Symmatrix S(A.dim(),0.);
  for(Integer i=0;i<A.dim();i++)
    S(i,i)=A(i);
  chk_set_init(S,1);
  return S;
}

Matrix sumrows(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"sumrows(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"sumrows(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Matrix v(1,A.nr);
 Real sum;
 Integer i,j;
 for(j=0;j<A.nr;j++){
     sum=0.;
     for(i=0;i<A.nr;i++)
         sum+=A(i,j);
     v(j)=sum;
 }
#if (CONICBUNDLE_DEBUG>=1)
 v.set_init(1);
#endif 
 return v;
}

Matrix sumcols(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"sumcols(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"sumcols(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Matrix v(A.nr,1);
 Real sum;
 Integer i,j;
 for(i=0;i<A.nr;i++){
     sum=0.;
     for(j=0;j<A.nr;j++)
         sum+=A(i,j);
     v(i)=sum;
 }
#if (CONICBUNDLE_DEBUG>=1)
 v.set_init(1);
#endif 
 return v;
}

Real sum(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"sum(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"sum(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 register Real s=0.;
 register Real s2=0;
 register Integer i,j;
 register const Real *ap=A.m;
 for(i=A.nr;--i>=0;){
     s+=(*ap++);
     s2=0.;
     for(j=i;--j>=0;){
         s2+=(*ap++);
     }
     s+=2.*s2;
 }
 return s;
}

Matrix maxrows(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"maxrows(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"maxrows(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Matrix v(1,A.nr);
 Real maxd;
 Integer i,j;
 for(j=0;j<A.nr;j++){
     maxd=A(0,j);
     for(i=1;i<A.nr;i++)
         maxd=max(maxd,A(i,j));
     v(j)=maxd;
 }
#if (CONICBUNDLE_DEBUG>=1)
 v.set_init(1);
#endif 
 return v;
}

Matrix maxcols(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"maxcols(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"maxcols(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Matrix v(A.nr,1);
 Real maxd;
 Integer i,j;
 for(i=0;i<A.nr;i++){
     maxd=A(i,0);
     for(j=1;j<A.nr;j++)
         maxd=max(maxd,A(i,j));
     v(i)=maxd;
 }
#if (CONICBUNDLE_DEBUG>=1)
 v.set_init(1);
#endif 
 return v;
}

Real max(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"max(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"max(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Integer i,j;
 Real maxd=A(0,0);
 for(i=0;i<A.nr;i++){
     maxd=max(maxd,A(i,i));
     for(j=i+1;j<A.nr;j++){
         maxd=max(maxd,A(i,j));
     }
 }
 return maxd;
}

Matrix minrows(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"minrows(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"minrows(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Matrix v(1,A.nr);
 Real mind;
 Integer i,j;
 for(j=0;j<A.nr;j++){
     mind=A(0,j);
     for(i=1;i<A.nr;i++)
         mind=min(mind,A(i,j));
     v(j)=mind;
 }
#if (CONICBUNDLE_DEBUG>=1)
 v.set_init(1);
#endif 
 return v;
}

Matrix mincols(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"mincols(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"mincols(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Matrix v(A.nr,1);
 Real mind;
 Integer i,j;
 for(i=0;i<A.nr;i++){
     mind=A(i,0);
     for(j=1;j<A.nr;j++)
         mind=min(mind,A(i,j));
     v(i)=mind;
 }
#if (CONICBUNDLE_DEBUG>=1)
 v.set_init(1);
#endif 
 return v;
}

Real min(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (A.nr==0)
     MEmessage(MEdim(A.nr,A.nr,0,0,"min(const Symmatrix&) dimension zero",MTsymmetric));
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"min(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Integer i,j;
 Real mind=A(0,0);
 for(i=0;i<A.nr;i++){
     mind=min(mind,A(i,i));
     for(j=i+1;j<A.nr;j++){
         mind=min(mind,A(i,j));
     }
 }
 return mind;
}

Real trace(const Symmatrix& A)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (!A.is_init)
     MEmessage(MatrixError(ME_unspec,"trace(const Symmatrix& A) A is not initialized",MTsymmetric));
#endif
 Real sum=0.;
 for(Integer i=0;i<A.nr;i++) sum+=A(i,i);
 return sum;
}

Real ip(const Symmatrix& A, const Symmatrix& B)
{
#if (CONICBUNDLE_DEBUG>=1)
 if (B.nr!=A.nr)
     MEmessage(MEdim(A.nr,A.nr,B.nr,B.nr,"ip(const Symmatrix&,const Symmatrix&) wrong dimensions",MTsymmetric));
 if ((!A.is_init)||(!B.is_init))
     MEmessage(MatrixError(ME_unspec,"ip(const Symmatrix&,const Symmatrix&) not initialized",MTsymmetric));
#endif
 register Real sum=0.;
 register Real s2;
 register Integer i,j;
 register const Real *ap=A.m;
 register const Real *bp=B.m;
 for(i=A.nr;--i>=0;){
     sum+=(*ap++)*(*bp++);
     s2=0.;
     for(j=i;--j>=0;){
         s2+=(*ap++)*(*bp++);
     }
     sum+=2.*s2;
 }
 return sum;
}

Real ip(const Symmatrix& A, const Matrix& B)
{
  chk_add(A,B);
  Matrix C(A);
  return ip(C,B);
}

Real ip(const Matrix& A, const Symmatrix& B)
{
  chk_add(A,B);
  Matrix C(B);
  return ip(A,C);
}

Symmatrix abs(const Symmatrix& A)
{
  chk_init(A);
  Symmatrix B; B.newsize(A.nr);
  Real *ap=A.m;
  Real *bp=B.m;
  for(Integer i=(A.nr*(A.nr+1))/2;--i>=0;)
    (*bp++)=fabs(*ap++);
  chk_set_init(B,1);
  return B;
}

std::ostream& operator<<(std::ostream& o,const Symmatrix &A)
{
 chk_init(A);
 o<<A.nr<<'\n';
 Integer i,j;
 for(i=0;i<A.nr;i++){
     for(j=i;j<A.nr;j++) o<<' '<<A(i,j);
     o<<'\n';
 }
 return o;
}

std::istream& operator>>(std::istream& in,Symmatrix &A)
{
 Real d;
 in>>d;
 Integer nr=Integer(d+.5);
 if (nr<0)
     MEmessage(MEdim(nr,nr,0,0,"operator>>(std::istream&,Symmatrix&) dimension negative",MTsymmetric));
 A.newsize(nr);
 for(Integer i=0;i<nr;i++)
     for(Integer j=i;j<nr;j++)
         in>>A(i,j);
 chk_set_init(A,1);
 return in;
}

// **************************************************************************
//                       Matrix:: Symmatrix specific implementations
// **************************************************************************

Matrix& Matrix::xeya(const Symmatrix& A,Real d)
{
  chk_init(A);
  newsize(A.nr,A.nr);
  const Real *matp=A.m;
  if (d==1.){
    for(Integer i=0;i<nr;i++){
      Real *mcp=m+i*nr+i;
      Real *mrp=mcp+nr;
      (*mcp++)=(*matp++);
      for(Integer j=i+1;j<nr;j++,mrp+=nr){
	(*mcp++)=(*mrp)=(*matp++);
      }
    }
  }
  else if (d==-1.){
    for(Integer i=0;i<nr;i++){
      Real *mcp=m+i*nr+i;
      Real *mrp=mcp+nr;
      (*mcp++)=-(*matp++);
      for(Integer j=i+1;j<nr;j++,mrp+=nr){
	(*mcp++)=(*mrp)=-(*matp++);
      }
    }
  }
  else {
    for(Integer i=0;i<nr;i++){
      Real *mcp=m+i*nr+i;
      Real *mrp=mcp+nr;
      (*mcp++)=(*matp++)*d;
      for(Integer j=i+1;j<nr;j++,mrp+=nr){
	(*mcp++)=(*mrp)=(*matp++)*d;
      }
    }
  }
  chk_set_init(*this,1);
  return *this;
}

Matrix& Matrix::xpeya(const Symmatrix& A,Real d)
{
  chk_add(*this,A);
  const Real *matp=A.m;
  if (d==1.){
    for(Integer i=0;i<nr;i++){
      Real *mcp=m+i*nr+i;
      Real *mrp=mcp+nr;
      (*mcp++)+=(*matp++);
      for(Integer j=i+1;j<nr;j++,mrp+=nr){
	(*mcp++)+=(*matp);
        (*mrp)+=(*matp++);
      }
    }
  }
  else if (d==-1.){
    for(Integer i=0;i<nr;i++){
      Real *mcp=m+i*nr+i;
      Real *mrp=mcp+nr;
      (*mcp++)-=(*matp++);
      for(Integer j=i+1;j<nr;j++,mrp+=nr){
	(*mcp++)-=(*matp);
        (*mrp)-=(*matp++);
      }
    }
  }
  else {
    for(Integer i=0;i<nr;i++){
      Real *mcp=m+i*nr+i;
      Real *mrp=mcp+nr;
      (*mcp++)+=(*matp++)*d;
      for(Integer j=i+1;j<nr;j++,mrp+=nr){
        Real f=d*(*matp++);
	(*mcp++)+=f;
        (*mrp)+=f;
      }
    }
  }
  chk_set_init(*this,1);
  return *this;
}

}

