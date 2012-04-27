/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/sparsmat.cxx

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
#include "sparsmat.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

Sparsemat& Sparsemat::xeya(const Sparsemat& A,Real d)
{
 chk_init(A);
 if (d==0.) return init(A.nr,A.nc);
 nr=A.nr;
 nc=A.nc;
 colinfo=A.colinfo;
 rowinfo=A.rowinfo;
 colindex=A.colindex;
 rowindex=A.rowindex;
 colval.xeya(A.colval,d);
 rowval.xeya(A.rowval,d);
 chk_set_init(*this,1);
 return *this;
}

Sparsemat& Sparsemat::xeya(const Matrix& A,Real d)
  //returns (*this)=A*d;
{
 chk_init(A);
 if (d==0.) return init(A.nr,A.nc);
 Integer i,nz=0;
 Integer di=A.nr*A.nc;
 const Real *m=A.get_store();
 for(i=0;i<di;i++){
     if (abs(*m++)>tol) nz++;
 }
 Indexmatrix I(nz,1),J(nz,1);
 Matrix val(nz,1);
 chk_set_init(I,1);
 chk_set_init(J,1);
 chk_set_init(val,1);
 m=A.get_store();
 nz=0;
 for(i=0;i<di;i++,m++){
     if (abs(*m)>tol){
         I(nz)=i%A.nr;
         J(nz)=i/A.nr;
         val(nz++)=*m*d;
     }
 }
 init(A.nr,A.nc,nz,I,J,val);
 return *this;
}

Sparsemat& Sparsemat::xeya(const Indexmatrix& A,Real d)
  //returns (*this)=A*d;
{
 chk_init(A);
 if (d==0.) return init(A.nr,A.nc);
 Integer i,nz=0;
 Integer di=A.nr*A.nc;
 const Integer *m=A.get_store();
 for(i=0;i<di;i++){
     if (*m++!=0) nz++;
 }
 Indexmatrix I(nz,1),J(nz,1);
 Matrix val(nz,1);
 chk_set_init(I,1);
 chk_set_init(J,1);
 chk_set_init(val,1);
 m=A.get_store();
 nz=0;
 for(i=0;i<di;i++,m++){
     if (*m!=0){
         I(nz)=i%A.nr;
         J(nz)=i/A.nr;
         val(nz++)=Real(*m*d);
     }
 }
 init(A.nr,A.nc,nz,I,J,val);
 return *this; 
}


///returns A= beta*B+alpha*A, where y may be transposed (ytrans=1); if alpha==0. then A is initialized to the correct size
Sparsemat& xbpeya(Sparsemat& A,const Sparsemat& B,Real beta,Real alpha,int Btrans)
{
 chk_init(B);
 if (alpha==0.){
   A.init(B,beta);
   if (Btrans)
     A.transpose();
   return A;
 }
 chk_init(A);
 Integer nr=A.nr;
 Integer nc=A.nc;
#if (CONICBUNDLE_DEBUG>=1)
 if (Btrans){
   if ((nr!=B.nc)||(nc!=B.nr)) {
     MEmessage(MatrixError(ME_dim,"xbpeya: dimensions don't match",MTsparse));;
   }
 }
 else {
   if ((nr!=B.nr)||(nc!=B.nc)) {
         MEmessage(MatrixError(ME_dim,"xbpeya: dimensions don't match",MTsparse));;
   }
 }
#endif
 if ((beta==0.)||(B.nonzeros()==0)){
   return A*=alpha;
 }
 if (A.nonzeros()==0){
   if (Btrans) {
     A.rowinfo=B.colinfo;
     A.rowindex=B.colindex;
     A.rowval.init(B.colval,beta);
     A.colinfo=B.rowinfo;
     A.colindex=B.rowindex;
     A.colval.init(B.rowval,beta);
   }
   else
     A.init(B,beta);
   return A;
 }

 //temporary target matrices that are going to be swapped with A's infos 
 Indexmatrix inf;
 Indexmatrix ind;
 Matrix val;
 
 Indexmatrix* Ainf;
 Indexmatrix* Aind;
 Matrix* Aval;
 const Indexmatrix* Binf;
 const Indexmatrix* Bind;
 const Matrix* Bval;

 //treat rows first
 Ainf=&(A.rowinfo);
 Aind=&(A.rowindex);
 Aval=&(A.rowval);
 if (Btrans){
   Binf=&(B.colinfo);
   Bind=&(B.colindex);
   Bval=&(B.colval);
 }
 else {
   Binf=&(B.rowinfo);
   Bind=&(B.rowindex);
   Bval=&(B.rowval);
 }
 Integer maxrows= min(Ainf->rowdim()+Binf->rowdim(),nr);
 Integer maxnz=min(A.nonzeros()+B.nonzeros(),nr*nc);
 Integer maxnr=nr;
 Integer maxnc=nc;
 inf.newsize(maxrows,3); chk_set_init(inf,1);
 ind.newsize(maxnz,1); chk_set_init(ind,1);
 val.newsize(maxnz,1); chk_set_init(val,1);
 Integer ai=0;
 Integer a0=(*Ainf)(ai,0);
 Integer bi=0;
 Integer b0=(*Binf)(bi,0);
 Integer infnz=0;
 Integer nz=0;
 while((a0<maxnr)||(b0<maxnr)){
   while(a0<b0) {
     inf(infnz,0)=a0;
     inf(infnz,1)=(*Ainf)(ai,1);
     inf(infnz,2)=nz;
     Integer* indp=ind.get_store()+nz;
     Real* valp=val.get_store()+nz;
     nz+=inf(infnz,1);
     const Integer* endp=ind.get_store()+nz;
     const Integer* aip=Aind->get_store()+(*Ainf)(ai,2);
     const Real* avp=Aval->get_store()+(*Ainf)(ai,2);
     for(;indp!=endp;){
       (*indp++)=(*aip++);
       (*valp++)=alpha*(*avp++);
     }
     ai++;
     a0=(ai<(*Ainf).rowdim())? (*Ainf)(ai,0) : maxnr;
     infnz++;
     continue;
   }
   while(b0<a0){
     inf(infnz,0)=b0;
     inf(infnz,1)=(*Binf)(bi,1);
     inf(infnz,2)=nz;
     Integer* indp=ind.get_store()+nz;
     Real* valp=val.get_store()+nz;
     nz+=inf(infnz,1);
     const Integer* endp=ind.get_store()+nz;
     const Integer* bip=Bind->get_store()+(*Binf)(bi,2);
     const Real* bvp=Bval->get_store()+(*Binf)(bi,2);
     for(;indp!=endp;){
       (*indp++)=(*bip++);
       (*valp++)=beta*(*bvp++);
     }
     bi++;
     b0=(bi<Binf->rowdim())? (*Binf)(bi,0):maxnr;
     infnz++;
     continue;
   }
   while((a0==b0)&&(a0!=maxnr)){
     inf(infnz,0)=a0;
     inf(infnz,1)=nz;  //stores the old value of nz temporarily
     inf(infnz,2)=nz;
     Integer* indp=ind.get_store()+nz;
     Real* valp=val.get_store()+nz;
     const Integer* bip=Bind->get_store()+(*Binf)(bi,2);
     const Real* bvp=Bval->get_store()+(*Binf)(bi,2);
     const Integer* bend=bip+(*Binf)(bi,1);
     Integer bii=(bip!=bend)?(*bip++):maxnc;
     const Integer* aip=Aind->get_store()+(*Ainf)(ai,2);
     const Real* avp=Aval->get_store()+(*Ainf)(ai,2);
     const Integer* aend=aip+(*Ainf)(ai,1);
     Integer aii=(aip!=aend)?(*aip++):maxnc;
     while((aii<maxnc)||(bii<maxnc)){
       while(aii<bii){
	 (*indp++)=aii;
	 (*valp++)=alpha*(*avp++);
         aii=(aip!=aend)? (*aip++) : maxnc;
	 nz++;
       }
       while(bii<aii){
	 (*indp++)=bii;
	 (*valp++)=beta*(*bvp++);
         bii=(bip!=bend)?(*bip++):maxnc;
	 nz++;
       }
       while((aii==bii)&&(aii!=maxnc)){
	 (*indp++)=aii;
	 (*valp++)=alpha*(*avp++)+beta*(*bvp++);
         aii=(aip!=aend)? (*aip++) : maxnc;
         bii=(bip!=bend)? (*bip++) : maxnc;
	 nz++;
       }
     }
     inf(infnz,1)=nz-inf(infnz,1);
     infnz++;
     bi++;
     b0=(bi<Binf->rowdim())? (*Binf)(bi,0) : maxnr;
     ai++;
     a0=(ai<Ainf->rowdim())? (*Ainf)(ai,0) : maxnr;
   }
 }
 inf.delete_rows(Range(infnz,inf.rowdim()-1));
 A.rowinfo=inf;
 ind.reduce_length(nz);
 A.rowindex=ind;
 val.reduce_length(nz);
 A.rowval=val;

 //next treat columns
 Ainf=&(A.colinfo);
 Aind=&(A.colindex);
 Aval=&(A.colval);
 if (Btrans){
   Binf=&(B.rowinfo);
   Bind=&(B.rowindex);
   Bval=&(B.rowval);
 }
 else {
   Binf=&(B.colinfo);
   Bind=&(B.colindex);
   Bval=&(B.colval);
 }
 maxrows= min(Ainf->rowdim()+Binf->rowdim(),nr);
 maxnr=nc;
 maxnc=nr;
 inf.newsize(maxrows,3); chk_set_init(inf,1);
 ai=0;
 a0=(*Ainf)(ai,0);
 bi=0;
 b0=(*Binf)(bi,0);
 infnz=0;
 nz=0;
 while((a0<maxnr)||(b0<maxnr)){
   while(a0<b0) {
     inf(infnz,0)=a0;
     inf(infnz,1)=(*Ainf)(ai,1);
     inf(infnz,2)=nz;
     Integer* indp=ind.get_store()+nz;
     Real* valp=val.get_store()+nz;
     nz+=inf(infnz,1);
     const Integer* endp=ind.get_store()+nz;
     const Integer* aip=Aind->get_store()+(*Ainf)(ai,2);
     const Real* avp=Aval->get_store()+(*Ainf)(ai,2);
     for(;indp!=endp;){
       (*indp++)=(*aip++);
       (*valp++)=alpha*(*avp++);
     }
     ai++;
     a0=(ai<(*Ainf).rowdim())? (*Ainf)(ai,0) : maxnr;
     infnz++;
     continue;
   }
   while(b0<a0){
     inf(infnz,0)=b0;
     inf(infnz,1)=(*Binf)(bi,1);
     inf(infnz,2)=nz;
     Integer* indp=ind.get_store()+nz;
     Real* valp=val.get_store()+nz;
     nz+=inf(infnz,1);
     const Integer* endp=ind.get_store()+nz;
     const Integer* bip=Bind->get_store()+(*Binf)(bi,2);
     const Real* bvp=Bval->get_store()+(*Binf)(bi,2);
     for(;indp!=endp;){
       (*indp++)=(*bip++);
       (*valp++)=beta*(*bvp++);
     }
     bi++;
     b0=(bi<Binf->rowdim())? (*Binf)(bi,0):maxnr;
     infnz++;
     continue;
   }
   while((a0==b0)&&(a0!=maxnr)){
     inf(infnz,0)=a0;
     inf(infnz,1)=nz;  //stores the old value of nz temporarily
     inf(infnz,2)=nz;
     Integer* indp=ind.get_store()+nz;
     Real* valp=val.get_store()+nz;
     const Integer* bip=Bind->get_store()+(*Binf)(bi,2);
     const Real* bvp=Bval->get_store()+(*Binf)(bi,2);
     const Integer* bend=bip+(*Binf)(bi,1);
     Integer bii=(bip!=bend)?(*bip++):maxnc;
     const Integer* aip=Aind->get_store()+(*Ainf)(ai,2);
     const Real* avp=Aval->get_store()+(*Ainf)(ai,2);
     const Integer* aend=aip+(*Ainf)(ai,1);
     Integer aii=(aip!=aend)?(*aip++):maxnc;
     while((aii<maxnc)||(bii<maxnc)){
       while(aii<bii){
	 (*indp++)=aii;
	 (*valp++)=alpha*(*avp++);
         aii=(aip!=aend)? (*aip++) : maxnc;
	 nz++;
       }
       while(bii<aii){
	 (*indp++)=bii;
	 (*valp++)=beta*(*bvp++);
         bii=(bip!=bend)?(*bip++):maxnc;
	 nz++;
       }
       while((aii==bii)&&(aii!=maxnc)){
	 (*indp++)=aii;
	 (*valp++)=alpha*(*avp++)+beta*(*bvp++);
         aii=(aip!=aend)? (*aip++) : maxnc;
         bii=(bip!=bend)? (*bip++) : maxnc;
	 nz++;
       }
     }
     inf(infnz,1)=nz-inf(infnz,1);
     infnz++;
     bi++;
     b0=(bi<Binf->rowdim())? (*Binf)(bi,0) : maxnr;
     ai++;
     a0=(ai<Ainf->rowdim())? (*Ainf)(ai,0) : maxnr;
   }
 }
 inf.delete_rows(Range(infnz,inf.rowdim()-1));
 swap(A.colinfo,inf);
 ind.reduce_length(nz);
 swap(A.colindex,ind);
 val.reduce_length(nz);
 swap(A.colval,val);
       
 return A;
}
  
 
  
Matrix& genmult(const Sparsemat& A,const Matrix& B,Matrix &C,
                       Real alpha,Real beta, int atrans,int btrans)
{
 chk_init(A);
 chk_init(B);
 Integer nr,nm,nc;
 if (atrans) {nr=A.nc;nm=A.nr;}
 else { nr=A.nr; nm=A.nc;}
 if (btrans) {
     nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nc) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
 }
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.nr)||(nc!=C.nc)) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if (atrans){
     if (btrans){
         const Integer *aip=A.colindex.get_store();
         const Real *avp=A.colval.get_store();
	 Integer i=A.colinfo.rowdim();
	 const Integer *iip=A.colinfo.get_store();
	 const Integer *nzp=iip+i;
         for(;--i>=0;){
             Real *cp=C.get_store()+(*iip++);
             for(Integer j=(*nzp++);--j>=0;){
                 mat_xpeya(nc,cp,nr,B.get_store()+(*aip++)*nc,1,alpha*(*avp++));
             }
         }
     }
     else {
         const Integer *aip=A.colindex.get_store();
         const Real *avp=A.colval.get_store();
	 Integer i=A.colinfo.rowdim();
	 const Integer *iip=A.colinfo.get_store();
	 const Integer *nzp=iip+i;
         for(;--i>=0;){
             Real *cp=C.get_store()+(*iip++);
             for(Integer j=(*nzp++);--j>=0;){
                 mat_xpeya(nc,cp,nr,B.get_store()+*aip++,nm,
                           alpha*(*avp++));
             }
         }
     }
 }
 else {
     if (btrans){
         const Integer *aip=A.rowindex.get_store();
         const Real *avp=A.rowval.get_store();
	 Integer i=A.rowinfo.rowdim();
	 const Integer *iip=A.rowinfo.get_store();
	 const Integer *nzp=iip+i;
         for(;--i>=0;){
             Real *cp=C.get_store()+(*iip++);
             for(Integer j=(*nzp++);--j>=0;){
                 mat_xpeya(nc,cp,nr,B.get_store()+(*aip++)*nc,1,alpha*(*avp++));
             }
         }
     }
     else {
         const Integer *aip=A.rowindex.get_store();
         const Real *avp=A.rowval.get_store();
	 Integer i=A.rowinfo.rowdim();
	 const Integer *iip=A.rowinfo.get_store();
	 const Integer *nzp=iip+i;
         for(;--i>=0;){
             Real *cp=C.get_store()+(*iip++);
             for(Integer j=(*nzp++);--j>=0;){
                 mat_xpeya(nc,cp,nr,B.get_store()+*aip++,nm,
                           alpha*(*avp++));
             }
         }
     }
 }
 return C;
}

Matrix& genmult(const Matrix& A,const Sparsemat& B,Matrix &C,
                       Real alpha,Real beta, int atrans,int btrans)
{
 chk_init(A);
 chk_init(B);
 Integer nr,nm,nc;
 if (atrans) {nr=A.nc;nm=A.nr;}
 else { nr=A.nr; nm=A.nc;}
 if (btrans) {
     nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nc) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
 }
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.nr)||(nc!=C.nc)) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if (atrans){
     if (btrans){ //A^T*B^T
         const Integer *bip=B.rowindex.get_store();
         const Real *bvp=B.rowval.get_store();
         for(Integer i=0;i<B.rowinfo.rowdim();i++){
             Integer ii=B.rowinfo(i,0);
             Real *cp=C.get_store()+ii*nr;
             for(Integer j=B.rowinfo(i,1);--j>=0;){
                 mat_xpeya(nr,cp,1,A.get_store()+(*bip++),nm,alpha*(*bvp++));
             }
         }
     }
     else { //A^T*B
         const Integer *bip=B.colindex.get_store();
         const Real *bvp=B.colval.get_store();
         for(Integer i=0;i<B.colinfo.rowdim();i++){
             Integer ii=B.colinfo(i,0);
             Real *cp=C.get_store()+ii*nr;
             for(Integer j=B.colinfo(i,1);--j>=0;){
                 mat_xpeya(nr,cp,1,A.get_store()+(*bip++),nm,alpha*(*bvp++));
             }
         }
     }
 }
 else {
     if (btrans){ // A*B^T
         const Integer *bip=B.rowindex.get_store();
         const Real *bvp=B.rowval.get_store();
         for(Integer i=0;i<B.rowinfo.rowdim();i++){
             Integer ii=B.rowinfo(i,0);
             Real *cp=C.get_store()+ii*nr;
             for(Integer j=B.rowinfo(i,1);--j>=0;){
                 mat_xpeya(nr,cp,A.get_store()+(*bip++)*nr,alpha*(*bvp++));
             }
         }
     }
     else { // A*B
         const Integer *bip=B.colindex.get_store();
         const Real *bvp=B.colval.get_store();
         for(Integer i=0;i<B.colinfo.rowdim();i++){
             Integer ii=B.colinfo(i,0);
             Real *cp=C.get_store()+ii*nr;
             for(Integer j=B.colinfo(i,1);--j>=0;){
                 mat_xpeya(nr,cp,A.get_store()+(*bip++)*nr,alpha*(*bvp++));
             }
         }
     }
 }
 return C;
}

 
Matrix& genmult(const Sparsemat& A,const Sparsemat& B,Matrix &C,
		Real alpha,Real beta, int atrans,int btrans)
{
  chk_init(A);
  chk_init(B);
  Integer nr,nc;
  const Indexmatrix* ainf;
  const Integer* aind;
  const Real* aval;
#if (CONICBUNDLE_DEBUG>=1)
  Integer nm;
#endif 
  if (atrans) {
    nr=A.nc;
#if (CONICBUNDLE_DEBUG>=1)
    nm=A.nr;
#endif 
    ainf=&(A.rowinfo);
    aind=A.rowindex.get_store();
    aval=A.rowval.get_store();
  }
  else {
    nr=A.nr;
#if (CONICBUNDLE_DEBUG>=1)
    nm=A.nc;
#endif 
    ainf=&(A.colinfo);
    aind=A.colindex.get_store();
    aval=A.colval.get_store();
  }
  
  const Indexmatrix* binf;
  const Integer* bind;
  const Real* bval;
  if (btrans) {
    nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
    if (nm!=B.nc) {
      MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
    }
#endif
    binf=&(B.rowinfo);
    bind=B.rowindex.get_store();
    bval=B.rowval.get_store();
  }
  else {
    nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
    if (nm!=B.nr) {
      MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsparse));;
    }
#endif
    binf=&(B.colinfo);
    bind=B.colindex.get_store();
    bval=B.colval.get_store();
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
  
  Integer Anzrows=ainf->rowdim();
  Integer nzlog2=1;
  while (Anzrows>>=1)
    nzlog2++;
  Anzrows=ainf->rowdim();
  
  for(Integer Bj=0;Bj<binf->rowdim();Bj++){//take column i of matrix B
    Real* cp=C.get_store()+(*binf)(Bj,0)*nr;
    const Integer* bip=bind+(*binf)(Bj,2);
    const Real* bvp=bval+(*binf)(Bj,2);
    const Integer* bend=bip+(*binf)(Bj,1);
    bool binsearch=(Anzrows>(*binf)(Bj,1)*nzlog2);
    Integer Astart=0;
    Integer Aj=0;
    Integer Aend=Anzrows;
    for(;bip!=bend;bip++,bvp++){ //run through all nonzeros of this column of B
      //find the next nonzero column of A with at least this column index
      if (binsearch){
	Aend=Anzrows;
	while(Astart<Aend){
	  Aj=(Astart+Aend)/2;
	   if ((*ainf)(Aj,0)<*bip){
	     Astart=Aj+1;
	   }
	   else if ((*ainf)(Aj,0)>*bip){
	     Aend=Aj;
	   }
	   else break;
	}
	if (Astart==Aend){
	  if (Astart==Anzrows){
	    break;
	  }
	  continue;
	}
      }
      else {
	while((Aj<Anzrows)&&((*ainf)(Aj,0)<*bip))
	  Aj++;
      }
      if (Aj==Anzrows)
	break;
      if ((*ainf)(Aj,0)>*bip)
	continue;
      Real d=alpha*(*bvp);
      const Integer* aip=aind+(*ainf)(Aj,2);
      const Integer* aend=aip+(*ainf)(Aj,1);
      const Real* avp=aval +(*ainf)(Aj,2);
      for (;aip!=aend;)
	*(cp+*aip++)+=d*(*avp++);
    }
  }
  return C;
}


Sparsemat& Sparsemat::init(Integer in_nr,Integer in_nc,Integer nz,
              const Integer *ini,const Integer *inj,const Real* va)
{
 return init(in_nr,in_nc,nz,Indexmatrix(nz,1,ini),Indexmatrix(nz,1,inj),Matrix(nz,1,va));
}


Sparsemat& Sparsemat::init(Integer in_nr,Integer in_nc,Integer nz,
          const Indexmatrix& ini,const Indexmatrix& inj, const Matrix& va)
{
 chk_init(ini);
 chk_init(inj);
 chk_init(va);
#if (CONICBUNDLE_DEBUG>=1)
 if (nz<0) {
     MEmessage(MatrixError(ME_unspec,"Sparsemat::init(nr,nc,nz,indi,indj,val): nz<0",MTsparse));
 }
#endif
 if (nz==0) {
     return init(in_nr,in_nc);
 }
#if (CONICBUNDLE_DEBUG>=1)
 if ((ini.dim()<nz)||(inj.dim()<nz)||(va.dim()<nz)) {  
     MEmessage(MatrixError(ME_unspec,"Sparsemat::init(nr,nc,nz,indi,indj,val): indi, indj, or val has <nz elements",MTsparse));
 }
 if ((min(ini(Range(0,nz-1)))<0)||(max(ini(Range(0,nz-1)))>=in_nr)||
     (min(inj(Range(0,nz-1)))<0)||(max(inj(Range(0,nz-1)))>=in_nc)) {  
     MEmessage(MatrixError(ME_unspec,"Sparsemat::init(nr,nc,nz,indi,indj,val): indices in indi or indj exceed range",MTsparse));
 }
#endif

 nr=in_nr;
 nc=in_nc;

 if (nz==1) {
     colinfo.newsize(1,3); chk_set_init(colinfo,1);
     colinfo(0)=inj(0);
     colinfo(1)=1;
     colinfo(2)=0;
     colindex.init(1,1,ini(0));
     colval.init(1,1,va(0));
     rowinfo.newsize(1,3); chk_set_init(rowinfo,1);
     rowinfo(0)=ini(0);
     rowinfo(1)=1;
     rowinfo(2)=0;
     rowindex.init(1,1,inj(0));
     rowval.init(1,1,va(0));
     return *this;
 }
     
 colindex.newsize(nz,1); chk_set_init(colindex,1);
 colval.newsize(nz,1); chk_set_init(colval,1);
 rowindex.newsize(nz,1); chk_set_init(rowindex,1);
 rowval.newsize(nz,1); chk_set_init(rowval,1);
 
 Integer i;

 //--- generate structure for columns and rows
 Indexmatrix sindc,sindr;

 //--- find number of columns and elements in each column and row by sorting
 colindex.init(nz,1,inj.get_store());   //has to be copied, because length of inj unknown
                                        //values inj*nr+ini might be too big for large nr and nc
 sortindex(colindex,sindc);             
 rowindex.init(nz,1,ini.get_store());   //has to be copied, because length of ini unknown
 sortindex(rowindex,sindr);             
 
 Integer lastcval=inj(sindc(0));  //holds last column index encountered
 Integer lastrval=ini(sindr(0));  //holds last row index encountered
 Integer colcnt=0;                //counts number of columns
 Integer rowcnt=0;                //counts number of rows
 Integer cntcelems=1;             //counts number of elements in current column
                                  //number will be stored in colindex temporarily
 Integer cntrelems=1;             //counts number of elements in current row
                                  //number will be stored in rowindex temporarily
 for(i=1;i<nz;i++){
     if (lastcval!=inj(sindc(i))){  //new column starts
         colindex(colcnt)=cntcelems;
         //sort indices in this column in incereasing rowindex order
         Integer *sindp=sindc.get_store()+i-cntcelems;
         //heapsort(cntcelems,sindp,ini.get_store());
	 std::sort( sindp,sindp+cntcelems, mat_less_index<Integer>(ini.get_store()) );
         //start new column
         lastcval=inj(sindc(i));
         cntcelems=1;
         colcnt++;
     }
     else {
         cntcelems++;
     }
     if (lastrval!=ini(sindr(i))){  //new column starts
         rowindex(rowcnt)=cntrelems;
         //sort indices in this row in increasing colindex order
         Integer *sindp=sindr.get_store()+i-cntrelems;
         //heapsort(cntrelems,sindp,inj.get_store());
	 std::sort( sindp,sindp+cntrelems, mat_less_index<Integer>(inj.get_store()) );
         lastrval=ini(sindr(i));
         cntrelems=1;
         rowcnt++;
     }
     else {
         cntrelems++;
     } 
 }
 colindex(colcnt)=cntcelems;
 //sort indices in this column in incereasing rowindex order
 Integer *sindp=sindc.get_store()+i-cntcelems;
 //heapsort(cntcelems,sindp,ini.get_store());
 std::sort( sindp,sindp+cntcelems, mat_less_index<Integer>(ini.get_store()) );
 colcnt++;
 rowindex(rowcnt)=cntrelems;
 //sort indices in this row in increasing colindex order
 sindp=sindr.get_store()+i-cntrelems;
 //heapsort(cntrelems,sindp,inj.get_store());
 std::sort( sindp,sindp+cntrelems, mat_less_index<Integer>(inj.get_store()) );
 rowcnt++;
 
 //--- initialize colinfo, colindex and colval to final values
 colinfo.newsize(colcnt,3); chk_set_init(colinfo,1);
 colcnt=0;
 lastcval=-1;
 cntcelems=0;
 Integer lastindex=-1;
 Integer new_nz=0;
 for(i=0;i<nz;i++){
   if (lastcval!=inj(sindc(i))){  //new column starts
     if (colcnt==0){	 
       colinfo(0,2)=0;
     }
     else {
       if(cntcelems>0){
	 colinfo(colcnt-1,1)=cntcelems;
	 colinfo(colcnt,2)=colinfo(colcnt-1,2)+cntcelems;
       }
       else {
	 colcnt--;
       }
     }
     lastcval=colinfo(colcnt,0)=inj(sindc(i));
     colcnt++;
     cntcelems=0;
     lastindex=-1;
     if (abs((colval(new_nz)=va(sindc(i))))>=tol){
       lastindex=colindex(new_nz)=ini(sindc(i));
       new_nz++;
       cntcelems++;
     }
   }
   else {
     if (lastindex!=ini(sindc(i))){ 
       if (abs((colval(new_nz)=va(sindc(i))))>=tol){
	 lastindex=colindex(new_nz)=ini(sindc(i));
	 cntcelems++;
	 new_nz++;
       }
     }
     else { //same element as before, add
       if (abs((colval(new_nz-1)+=va(sindc(i))))<tol){
	 cntcelems--;
	 new_nz--;
	 lastindex=-1;
       }
     }
   }
 }
 if(cntcelems>0){
   colinfo(colcnt-1,1)=cntcelems;
 }
 else {
   colcnt--;
 }
 colval.reduce_length(new_nz);
 colindex.reduce_length(new_nz);
 if (colcnt<colinfo.rowdim()){
   colinfo.delete_rows(Range(colcnt,colinfo.rowdim()-1));
 }
 
 //--- initialize rowinfo, rowindex and rowval to final values
 rowinfo.newsize(rowcnt,3); chk_set_init(rowinfo,1);
 rowcnt=0;
 lastrval=-1;
 cntrelems=0;
 lastindex=-1;
 new_nz=0;
 for(i=0;i<nz;i++){
   if (lastrval!=ini(sindr(i))){  //new row starts
     if (rowcnt==0){	 
       rowinfo(0,2)=0;
     }
     else {
       if (cntrelems>0){
	 rowinfo(rowcnt-1,1)=cntrelems;
	 rowinfo(rowcnt,2)=rowinfo(rowcnt-1,2)+cntrelems;
       }
       else {
	 rowcnt--;
       }
     }
     lastrval=rowinfo(rowcnt,0)=ini(sindr(i));
     rowcnt++;
     cntrelems=0;
     lastindex=-1;
     if (abs((rowval(new_nz)=va(sindr(i))))>=tol){
       lastindex=rowindex(new_nz)=inj(sindr(i));
       new_nz++;
       cntrelems++;
     }
   }
   else {
     if (lastindex!=inj(sindr(i))){ 
       if (abs((rowval(new_nz)=va(sindr(i))))>=tol){
	 lastindex=rowindex(new_nz)=inj(sindr(i));
	 cntrelems++;
	 new_nz++;
       }
     }
     else { //same element as before, add
       if (abs((rowval(new_nz-1)+=va(sindr(i))))<tol){
	 cntrelems--;
	 new_nz--;
	 lastindex=-1;
       }
     }
   }
 }
 if (cntrelems>0){
   rowinfo(rowcnt-1,1)=cntrelems;
 }
 else {
   rowcnt--;
 }
 rowval.reduce_length(new_nz);
 rowindex.reduce_length(new_nz);
 if (rowcnt<rowinfo.rowdim()){
   rowinfo.delete_rows(Range(rowcnt,rowinfo.rowdim()-1));
 }
 
 return *this;
}

void Sparsemat::get_edge_rep(Indexmatrix& I, Indexmatrix& J, Matrix& val) const
{
 I.newsize(rowindex.dim(),1); chk_set_init(I,1);
 J.newsize(rowindex.dim(),1); chk_set_init(J,1);
 val.newsize(rowindex.dim(),1); chk_set_init(val,1);
 Integer i,j;
 Integer nz=0;
 for(i=0;i<rowinfo.rowdim();i++){
     Integer row=rowinfo(i,0);
     Integer nrcols=rowinfo(i,1);
     for(j=0;j<nrcols;j++){
         I(nz)=row;
         J(nz)=rowindex(nz);
         val(nz)=rowval(nz);
         nz++;
     }
 }
}

int Sparsemat::get_edge(Integer e,Integer& indi,Integer& indj,Real& val) const
{
  if ((e<0)||(e>=rowinfo.dim())) 
    return 1;
  val=rowval(e);//includes initialization and range check for e
  indj=rowindex(e); 
  //binary search for row of e
  Integer lb=0;
  Integer ub=rowinfo.rowdim()-1;
  Integer ii;
  while (lb<=ub){
    ii=(lb+ub)/2;
    if (rowinfo(ii,2)>e) {
      ub=ii-1;
      continue;
    }
    if (rowinfo(ii,2)<e) { 
      lb=ii+1;
      continue;
    }
    break;
  }
  if (lb>ub){ //then ub holds the correct index
    indi=rowinfo(ub,0);
  }
  else {
    indi=rowinfo(ii);
  }
  return 0;
}

    
// ============================================================================
//                               contains_support
// ============================================================================

int Sparsemat::contains_support(const Sparsemat& A) const
{
 chk_init(*this);
 chk_init(A);
 if ((nr!=A.nr)||(nc!=A.nc)) return 0;
 if (A.colinfo.rowdim()==0) return 1;
 if ((colinfo.rowdim()==0) && (A.colinfo.rowdim()==0)) return 1;
 if (colval.rowdim()<A.colval.rowdim()) return 0;
 if (colinfo.rowdim()<A.colinfo.rowdim()) return 0;
 if (rowinfo.rowdim()<A.rowinfo.rowdim()) return 0;
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
             cv1=nc;
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
             cv1=nc;
         }
         cind2++;
         if (cind2<A.colinfo.rowdim()){
             cv2=A.colinfo(cind2,0);
         }
         else {
             cv2=nc;
         }
     }
 }while((cv1<nc)||(cv2<nc));
 
 return 1;
}


Real Sparsemat::operator()(Integer i,Integer j) const
{
 chk_init(*this);
 chk_range(i,j,nr,nc);
 Integer nrows=rowinfo.rowdim();
 if (nrows==0) return 0.;
 Integer ncols=colinfo.rowdim();
 if ((i<rowinfo(0,0))||(i>rowinfo(nrows-1,0))) return 0.;
 if ((j<colinfo(0,0))||(j>colinfo(ncols-1,0))) return 0.;
 if (nrows<=ncols){
     Integer li=0; 
     Integer ui=nrows-1;
     Integer mi=0;
     while(li<=ui) {
         mi=(li+ui)/2;
         if (rowinfo(mi,0)==i) break;
         if (rowinfo(mi,0)<i){
             li=mi+1;
         }
         else {
             ui=mi-1;
         }
     }
     if (li>ui) return 0.;
     li=rowinfo(mi,2);
     ui=li+rowinfo(mi,1)-1;
     while(li<=ui) {
         mi=(li+ui)/2;
         if (rowindex(mi)==j) break;
         if (rowindex(mi)<j){
             li=mi+1;
         }
         else {
             ui=mi-1;
         }
     }
     if (li>ui) return 0.;
     return rowval(mi);
 }
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

Integer Sparsemat::col_nonzeros(Integer ci,Integer* fi) const
{
  chk_init(*this);
  chk_range(0,ci,nr,nc);
  if (fi!=0) *fi=-1;
  if (colinfo.dim()==0) return 0; //no nonzero columns
  if ((ci<colinfo(0,0))||(ci>colinfo(colinfo.rowdim()-1,0))) return 0;
  //binary search for column ci
  Integer lb=0;
  Integer ub=colinfo.rowdim()-1;
  Integer ii;
  while (lb<=ub){
    ii=(lb+ub)/2;
    if (colinfo(ii,0)<ci) {
      lb=ii+1;
      continue;
    }
    if (colinfo(ii,0)>ci) { 
      ub=ii-1;
      continue;
    }
    break;
  }
  if (lb>ub) return 0; //column ci is zero
  //--- column ci is stored in ii
  if (fi!=0) *fi=colinfo(ii,2); 
  return colinfo(ii,1);
}

Integer Sparsemat::row_nonzeros(Integer ri, Integer *fi) const
{
  chk_init(*this);
  chk_range(ri,0,nr,nc);
  if (fi!=0) *fi=-1;
  if (rowinfo.dim()==0) return 0; //no nonzero rows
  if ((ri<rowinfo(0,0))||(ri>rowinfo(rowinfo.rowdim()-1,0))) return 0;
  //binary search for row ri
  Integer lb=0;
  Integer ub=rowinfo.rowdim()-1;
  Integer ii;
  while (lb<=ub){
    ii=(lb+ub)/2;
    if (rowinfo(ii,0)<ri) {
      lb=ii+1;
      continue;
    }
    if (rowinfo(ii,0)>ri) { 
      ub=ii-1;
      continue;
    }
    break;
  }
  if (lb>ub) return 0; //column ci is zero
  //--- row ri is stored in ii
  if (fi!=0) *fi=rowinfo(ii,2); 
  return rowinfo(ii,1);
}

Sparsemat Sparsemat::col(Integer ci) const
{
  chk_init(*this);
  if ((nr==0)||(nc==0)){ 
     MEmessage(MatrixError(ME_unspec,"Sparsemat::col(Integer ci): n==0 or m==0",MTsparse));
  }
  chk_range(0,ci,nr,nc);
  Sparsemat A;
  A.nr=nr;
  A.nc=1;
  A.tol=tol;
  chk_set_init(A,1);
  if (colinfo.dim()==0) return A; //no nonzero columns
  if ((ci<colinfo(0,0))||(ci>colinfo(colinfo.rowdim()-1,0))) return A;
  //binary search for column ci
  Integer lb=0;
  Integer ub=colinfo.rowdim()-1;
  Integer ii;
  while (lb<=ub){
    ii=(lb+ub)/2;
    if (colinfo(ii,0)<ci) {
      lb=ii+1;
      continue;
    }
    if (colinfo(ii,0)>ci) { 
      ub=ii-1;
      continue;
    }
    break;
  }
  if (lb>ub) return A; //column ci is zero
  //--- column ci is stored in ii
  Integer nz=colinfo(ii,1);
  A.colinfo.newsize(1,3); chk_set_init(A.colinfo,1);
  A.colinfo(0,0)=0;
  A.colinfo(0,1)=nz;
  A.colinfo(0,2)=0;
  A.colindex.init(nz,1,colindex.get_store()+colinfo(ii,2));
  A.colval.init(nz,1,colval.get_store()+colinfo(ii,2));
  A.rowinfo.newsize(nz,3); chk_set_init(A.rowinfo,1);
  A.rowindex.init(nz,1,Integer(0)); //all go into column 0
  A.rowval=A.colval;
  for (Integer i=0;i<nz;i++){
    A.rowinfo(i,0)=A.colindex(i);
    A.rowinfo(i,1)=1;
    A.rowinfo(i,2)=i;
  }
  return A;
}

Sparsemat Sparsemat::row(Integer ri) const
{
  chk_init(*this);
  if ((nr==0)||(nc==0)){ 
     MEmessage(MatrixError(ME_unspec,"Sparsemat::row(Integer ri): n==0 or m==0",MTsparse));
  }
  chk_range(ri,0,nr,nc);
  Sparsemat A;
  A.nr=1;
  A.nc=nc;
  A.tol=tol;
  chk_set_init(A,1);
  if (rowinfo.dim()==0) return A; //no nonzero rows
  if ((ri<rowinfo(0,0))||(ri>rowinfo(rowinfo.rowdim()-1,0))) return A;
  //binary search for row ri
  Integer lb=0;
  Integer ub=rowinfo.rowdim()-1;
  Integer ii;
  while (lb<=ub){
    ii=(lb+ub)/2;
    if (rowinfo(ii,0)<ri) {
      lb=ii+1;
      continue;
    }
    if (rowinfo(ii,0)>ri) { 
      ub=ii-1;
      continue;
    }
    break;
  }
  if (lb>ub) return A; //rowumn ri is zero
  //--- row ri is stored in ii
  Integer nz=rowinfo(ii,1);
  A.rowinfo.newsize(1,3); chk_set_init(A.rowinfo,1);
  A.rowinfo(0,0)=0;
  A.rowinfo(0,1)=nz;
  A.rowinfo(0,2)=0;
  A.rowindex.init(nz,1,rowindex.get_store()+rowinfo(ii,2));
  A.rowval.init(nz,1,rowval.get_store()+rowinfo(ii,2));
  A.colinfo.newsize(nz,3); chk_set_init(A.colinfo,1);
  A.colindex.init(nz,1,Integer(0));
  A.colval=A.rowval;
  for (Integer i=0;i<nz;i++){
    A.colinfo(i,0)=A.rowindex(i);
    A.colinfo(i,1)=1;
    A.colinfo(i,2)=i;
  }
  return A;
}

Sparsemat Sparsemat::cols(const Indexmatrix& ind) const
{
  chk_init(*this);
  chk_init(ind);
#if (CONICBUNDLE_DEBUG>=1)
  if ((ind.dim()!=0)&&((min(ind)<0)||(max(ind)>=nc))){
     MEmessage(MatrixError(ME_unspec,"Sparsemat::cols(const Indexmatrix& ind): ind exceeds range",MTsparse));
  } 
#endif
  if ((ind.dim()==0)||(colinfo.rowdim()==0)) return Sparsemat(nr,ind.dim());
  if (ind.dim()==1) return this->col(ind(0));
  //find nonzero rows and number of nonzeros
  Integer newnz=0;
  Indexmatrix mind(ind.dim(),1,-1);
  for(Integer i=0;i<ind.dim();i++){
    Integer h=ind(i);
    Integer a=0; 
    Integer b=colinfo.rowdim()-1;
    Integer val;
    while(a<b){
      Integer c=(a+b)/2;
      val=colinfo(c,0);
      if (h<val){
	b=c-1;
      }
      else {
	if (h==val) {
	  a=c;
	  break;
	}
	a=c+1;
      }
    }
    if ((a<colinfo.rowdim())&&(h==colinfo(a,0))){
      mind(i)=a;
      newnz+=colinfo(a,1);
    }
  }
  Indexmatrix indi(newnz,1); chk_set_init(indi,1);
  Indexmatrix indj(newnz,1); chk_set_init(indj,1);
  Matrix val(newnz,1); chk_set_init(val,1);
  newnz=0;
  {for(Integer i=0;i<ind.dim();i++){
    Integer mi=mind(i);
    if (mi<0) continue;
    const Integer *ip=colindex.get_store()+colinfo(mi,2);
    const Real *vp=colval.get_store()+colinfo(mi,2);
    for(Integer j=colinfo(mi,1);--j>=0;newnz++){
      indi(newnz)=*ip++;
      indj(newnz)=i;
      val(newnz)=*vp++;
    }
  }}
  return Sparsemat(nr,ind.dim(),newnz,indi,indj,val);
}

Sparsemat Sparsemat::rows(const Indexmatrix& ind) const
{
  chk_init(*this);
  chk_init(ind);
#if (CONICBUNDLE_DEBUG>=1)
  if ((ind.dim()!=0)&&((min(ind)<0)||(max(ind)>=nr))){
     MEmessage(MatrixError(ME_unspec,"Sparsemat::rows(const Indexmatrix& ind): ind exceeds range",MTsparse));
  } 
#endif
  if ((ind.dim()==0)||(rowinfo.rowdim()==0)) return Sparsemat(ind.dim(),nc);
  if (ind.dim()==1) return this->row(ind(0));
  //find nonzero rows and number of nonzeros
  Integer newnz=0;
  Indexmatrix mind(ind.dim(),1,-1);
  for(Integer i=0;i<ind.dim();i++){
    Integer h=ind(i);
    Integer a=0; 
    Integer b=rowinfo.rowdim()-1;
    Integer val;
    while(a<b) {
      Integer c=(a+b)/2;
      val=rowinfo(c,0);
      if (h<val){
	b=c-1;
      }
      else {
	if (h==val) {
	  a=c;
	  break;
	}
	a=c+1;
      }
    }
    if ((a<rowinfo.rowdim())&&(h==rowinfo(a,0))){
      mind(i)=a;
      newnz+=rowinfo(a,1);
    }
  }
  Indexmatrix indi(newnz,1); chk_set_init(indi,1);
  Indexmatrix indj(newnz,1); chk_set_init(indj,1);
  Matrix val(newnz,1); chk_set_init(val,1);
  newnz=0;
  {for(Integer i=0;i<ind.dim();i++){
    Integer mi=mind(i);
    if (mi<0) continue;
    const Integer *ip=rowindex.get_store()+rowinfo(mi,2);
    const Real *vp=rowval.get_store()+rowinfo(mi,2);
    for(Integer j=rowinfo(mi,1);--j>=0;newnz++){
      indi(newnz)=*ip++;
      indj(newnz)=i;
      val(newnz)=*vp++;
    }
  }}
  return Sparsemat(ind.dim(),nc,newnz,indj,indi,val);
}


Sparsemat& Sparsemat::delete_rows(const Indexmatrix& ind)
{
  chk_init(*this);
  chk_init(ind);
  if (ind.dim()==0) return *this;
  Indexmatrix sind;
  sortindex(ind,sind);
  sind=ind(sind);
#if (CONICBUNDLE_DEBUG>=1)
  if ((sind(0)<0)||(sind(sind.dim()-1)>=nr)){
     MEmessage(MatrixError(ME_unspec,"Sparsemat::delete_rows(const Indexmatrix& ind): index exceeds range",MTsparse));
  } 
#endif

  //------------ first change the row representation
  Indexmatrix delinfo;
  Integer sinddim=sind.dim();
  Integer newnz=nonzeros();
  {
    Integer sindi=0;
    Integer sindival=sind(sindi);
    Integer i=0;
    while((i<rowinfo.rowdim())&&(rowinfo(i,0)<sindival))
      i++;
    if (i<rowinfo.rowdim())
      newnz=rowinfo(i,2);
    for(;i<rowinfo.rowdim();i++){
      Integer rowi=rowinfo(i,0);
      while (rowi>sindival) {
	if (++sindi<sinddim){
	  sindival=sind(sindi);
	}
	else {
	  sindival=nr+1;
	}
      }
      if (sindival==rowi){
	delinfo.concat_below(i);
	continue;
      }
      rowinfo(i,0)-=sindi;
      Integer nz=rowinfo(i,1);
      Integer offset=rowinfo(i,2);
      if (newnz!=offset) {//move entries to correct positions
        rowinfo(i,2)=newnz;
        mat_xey(nz,rowindex.get_store()+newnz,rowindex.get_store()+offset);
        mat_xey(nz,rowval.get_store()+newnz,rowval.get_store()+offset);
      }
      newnz+=nz;
    }
  }
  rowinfo.delete_rows(delinfo);
  rowindex.reduce_length(newnz);
  rowval.reduce_length(newnz);
  
  //------------ change the column info
  //compute log 2 of sind.dim()
  Integer log2=1;
  {for (Integer i=sind.dim();(i>>=1)>0;log2++){}}
  newnz=0;
  Integer* ncip=colindex.get_store();  //new column index pointer
  Real* ncvp=colval.get_store();       //new column value pointer
  Integer *ocip=ncip;                  //old ...
  Real* ocvp=ncvp;
  delinfo.init(0,0,Integer(0));
  for(Integer j=0;j<colinfo.rowdim();j++){
    Integer nz=colinfo(j,1);
    colinfo(j,2)=newnz;
    Integer sindi=0;
    Integer* sindp=sind.get_store();
    Integer sindival=(*sindp++);
    if (nz*log2>=sind.dim()){
      //it is more efficient to search for deleted indices going through all indices in sind
      for(Integer i=nz;--i>=0;){
	Integer indi=(*ocip++);
	while (indi>sindival) {
	  if (++sindi<sinddim){
	    sindival=(*sindp++);
	  }
	  else {
	    sindival=nr+1;
	  }
	}
	if (sindival==indi){
	  ocvp++;
	  continue;
	}
	(*ncip++)=indi-sindi;
	(*ncvp++)=(*ocvp++);
	newnz++;
      }
    }
    else {
      //it is more efficient to search in sind, whether the current index appears there 
      Integer lb=0;
      for(Integer i=nz;--i>=0;){
	Integer indi=(*ocip++);
	Integer ub=sind.dim()-1;
	while ((lb<=ub)&&(indi>sind(lb))){
	  Integer ii=(lb+ub)/2;
	  if (sind(ii)<indi) {
	    lb=ii+1;
	    continue;
	  }
	  if (sind(ii)>indi) { 
	    ub=ii-1;
	    continue;
	  }
	  lb=ii;
	  break;
	}
	if ((lb<=ub)&&(indi==sind(lb))){ //entry has to be deleted
	  ocvp++;
	  continue;
	}
	(*ncip++)=indi-lb;
	(*ncvp++)=(*ocvp++);
	newnz++;
      }
    }
    nz=newnz-colinfo(j,2);
    if (nz==0) {
      delinfo.concat_below(j);
    }
    colinfo(j,1)=nz;
  }
  colinfo.delete_rows(delinfo);
  colindex.reduce_length(newnz);
  colval.reduce_length(newnz);

  nr-=ind.dim();
 
  return *this;
}

Sparsemat& Sparsemat::delete_cols(const Indexmatrix& ind)
{
  chk_init(*this);
  chk_init(ind);
  if (ind.dim()==0) return *this;
  Indexmatrix sind;
  sortindex(ind,sind);
  sind=ind(sind);
#if (CONICBUNDLE_DEBUG>=1)
  if ((sind(0)<0)||(sind(sind.dim()-1)>=nc)){
    MEmessage(MatrixError(ME_unspec,"Sparsemat::delete_cols(const Indexmatrix& ind): index exceeds range",MTsparse));
  } 
#endif

  //first change the col representation
  Integer sinddim=sind.dim();
  Indexmatrix delinfo;
  Integer newnz=nonzeros();
  {
    Integer sindi=0;
    Integer sindival=sind(sindi);
    Integer i=0;
    while((i<colinfo.rowdim())&&(colinfo(i,0)<sindival))
      i++;
    if (i<colinfo.rowdim())
      newnz=colinfo(i,2);
    for(;i<colinfo.rowdim();i++){
      Integer coli=colinfo(i,0);
      while (coli>sindival) {
	if (++sindi<sinddim){
	  sindival=sind(sindi);
	}
	else {
	  sindival=nc+1;
	}
      }
      if (sindival==coli){
	delinfo.concat_below(i);
	continue;
      }
      colinfo(i,0)-=sindi;
      Integer nz=colinfo(i,1);
      Integer offset=colinfo(i,2);
      if (newnz!=offset) {//move entries to correct positions
	colinfo(i,2)=newnz;
	mat_xey(nz,colindex.get_store()+newnz,colindex.get_store()+offset);
	mat_xey(nz,colval.get_store()+newnz,colval.get_store()+offset);
      }
      newnz+=nz;
    }
  }
  colinfo.delete_rows(delinfo);
  colindex.reduce_length(newnz);
  colval.reduce_length(newnz);
  
  //------------ change the column info
  //compute log 2 of sind.dim()
  Integer log2=1;
  {for (Integer i=sind.dim();(i>>=1)>0;log2++){}}
  newnz=0;
  Integer* nrip=rowindex.get_store();  //new row index pointer
  Real* nrvp=rowval.get_store();       //new row value pointer
  Integer *orip=nrip;                  //old ...
  Real* orvp=nrvp;
  delinfo.init(0,0,Integer(0));
  for(Integer j=0;j<rowinfo.rowdim();j++){
    Integer nz=rowinfo(j,1);
    rowinfo(j,2)=newnz;
    Integer sindi=0;
    Integer* sindp=sind.get_store();
    Integer sindival=(*sindp++);
    if (nz*log2>=sind.dim()){
      //it is more efficient to search for deleted indices going through all indices in sind
      for(Integer i=nz;--i>=0;){
	Integer indi=(*orip++);
	while (indi>sindival) {
	  if (++sindi<sinddim){
	    sindival=(*sindp++);
	  }
	  else {
	    sindival=nc+1;
	  }
	}
	if (sindival==indi){
	  orvp++;
	  continue;
	}
	(*nrip++)=indi-sindi;
	(*nrvp++)=(*orvp++);
	newnz++;
      }
    }
    else {
      //it is more efficient to search in sind, whether the current index appears there 
      Integer lb=0;
      for(Integer i=nz;--i>=0;){
	Integer indi=(*orip++);
	Integer ub=sind.dim()-1;
	while ((lb<=ub)&&(indi>sind(lb))){
	  Integer ii=(lb+ub)/2;
	  if (sind(ii)<indi) {
	    lb=ii+1;
	    continue;
	  }
	  if (sind(ii)>indi) { 
	    ub=ii-1;
	    continue;
	  }
	  lb=ii;
	  break;
	}
	if ((lb<=ub)&&(indi==sind(lb))){ //entry has to be deleted
	  orvp++;
	  continue;
	}
	(*nrip++)=indi-lb;
	(*nrvp++)=(*orvp++);
	newnz++;
      }
    }
    nz=newnz-rowinfo(j,2);
    if (nz==0) {
      delinfo.concat_below(j);
    }
    rowinfo(j,1)=nz;
  }
  rowinfo.delete_rows(delinfo);
  rowindex.reduce_length(newnz);
  rowval.reduce_length(newnz);

  nc-=ind.dim();
 
  return *this;
}
 

Sparsemat& Sparsemat::insert_row(Integer i,const Sparsemat& v)
{
  chk_init(v);
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (nc!=v.dim()){
     MEmessage(MatrixError(ME_unspec,"Sparsemat::insert_row(Integer i,const Sparsemat& v): column dimensions do not match",MTsparse));
  } 
  if ((v.nr!=1)&&(v.nc!=1)){
     MEmessage(MatrixError(ME_unspec,"Sparsemat::insert_row(Integer i,const Sparsemat& v): v is not a vector",MTsparse));
  } 
#endif
  if ((rowinfo.rowdim()>0)&&(i>rowinfo(rowinfo.rowdim()-1,0))&&(v.nonzeros()==0)){
    nr++;
    return *this;
  }
  Indexmatrix tindi;
  Indexmatrix tindj;
  Matrix tval;
  get_edge_rep(tindi,tindj,tval);
  if ((rowinfo.rowdim()>0)&&(i<=rowinfo(rowinfo.rowdim()-1,0))){
    for(Integer k=0;k<tindi.dim();k++){
      if (tindi(k)>=i) {
	tindi(k)++;
      }
    }
  }
  if (v.nonzeros()){
    Indexmatrix vindi;
    Indexmatrix vindj;
    Matrix vval;
    v.get_edge_rep(vindi,vindj,vval);
    if (v.rowdim()>v.coldim()) swap(vindi,vindj);
    vindi.init(vindj.dim(),1,i);
    
    tindi.concat_below(vindi);
    tindj.concat_below(vindj);
    tval.concat_below(vval);
  }
  init(nr+1,nc,tval.dim(),tindi,tindj,tval);
  return *this;  
}

Sparsemat& Sparsemat::insert_col(Integer i,const Sparsemat& v)
{
  chk_init(v);
  chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
  if (nr!=v.dim()){
    MEmessage(MatrixError(ME_unspec,"Sparsemat::insert_row(Integer i,const Sparsemat& v): column dimensions do not match",MTsparse));
  } 
  if ((v.nr!=1)&&(v.nc!=1)){
    MEmessage(MatrixError(ME_unspec,"Sparsemat::insert_row(Integer i,const Sparsemat& v): v is not a vector",MTsparse));
  } 
#endif
  if ((colinfo.rowdim()>0)&&(i>colinfo(colinfo.rowdim()-1,0))&&(v.nonzeros()==0)){
    nr++;
    return *this;
  }
  Indexmatrix tindi;
  Indexmatrix tindj;
  Matrix tval;
  get_edge_rep(tindi,tindj,tval);
  if((colinfo.rowdim()>0)&&(i>colinfo(colinfo.rowdim()-1,0))){
    for(Integer k=0;k<tindj.dim();k++){
      if (tindj(k)>=i) {
	tindj(k)++;
      }
    }
  }
  if (v.nonzeros()){
    Indexmatrix vindi;
    Indexmatrix vindj;
    Matrix vval;
    v.get_edge_rep(vindi,vindj,vval);
    if (v.coldim()>v.rowdim()) swap(vindi,vindj);
    vindj.init(vindi.dim(),1,i);
    tindi.concat_below(vindi);
    tindj.concat_below(vindj);
    tval.concat_below(vval);
  }
  init(nr,nc+1,tval.dim(),tindi,tindj,tval);
  return *this;  
}    

Sparsemat& Sparsemat::concat_right(const Sparsemat& A)
{
  chk_init(A);
  chk_init(*this);
  if ((A.nr==0)&&(A.nc==0))
    return *this;
  if ((nr==0)&&(nc==0)){
    *this=A;
    return *this;
  } 
#if (CONICBUNDLE_DEBUG>=1)
  if (nr!=A.nr){
     MEmessage(MatrixError(ME_unspec,"Sparsemat::concat_right(Sparsemat& A): number of rows do not match",MTsparse));
  } 
#endif
  if ((nr==0)||(A.nonzeros()==0)){
    nc+=A.nc;
    return *this;
  }
  //column representation
  Integer infonr=colinfo.rowdim();
  Integer oldnz=colindex.dim();
  colinfo.concat_below(A.colinfo);
  colval.concat_below(A.colval);
  colindex.concat_below(A.colindex);
  for(Integer i=infonr;i<colinfo.rowdim();i++){
    colinfo(i,0)+=nc;
    colinfo(i,2)+=oldnz;
  }

  //row representation
  infonr=0;
  Integer thisi=0;
  Integer thisind=nr;
  if (thisi<rowinfo.rowdim()) thisind=rowinfo(thisi,0);
  Integer Ai=0;
  Integer Aind=nr;
  if (Ai<A.rowinfo.rowdim()) Aind=A.rowinfo(Ai,0);
  while((thisind<nr)||(Aind<nr)){
    if (thisind<Aind){
      infonr++;
      thisi++;
      if (thisi<rowinfo.rowdim()){
	thisind=rowinfo(thisi,0);
      }
      else {
	thisind=nr;
      }
      continue;
    }
    if (thisind>Aind){
      infonr++;
      Ai++;
      if (Ai<A.rowinfo.rowdim()){
	Aind=A.rowinfo(Ai,0);
      }
      else {
	Aind=nr;
      }
      continue;
    }
    infonr++;
    thisi++;
    if (thisi<rowinfo.rowdim()){
      thisind=rowinfo(thisi,0);
    }
    else {
      thisind=nr;
    }
    Ai++;
    if (Ai<A.rowinfo.rowdim()){
      Aind=A.rowinfo(Ai,0);
      }
    else {
      Aind=nr;
    }
  }
   
  int free_oldip=0;
  Integer *oldip=rowindex.get_store();
  if (rowindex.mem_dim<rowindex.dim()+A.rowindex.dim()){
    rowindex.mem_dim=Integer(rowindex.memarray->get(rowindex.dim()+A.rowindex.dim(),rowindex.m));
    if (rowindex.mem_dim<rowindex.dim()+A.rowindex.dim()){
      MEmessage(MEmem(rowindex.dim()+A.rowindex.dim(),
		      "Matrix::concat_right(const Matrix&) not enough memory",
		      MTmatrix));
    }
    free_oldip=1;
  } 
  int free_oldvp=0;  
  Real *oldvp=rowval.get_store();
  if (rowval.mem_dim<rowval.dim()+A.rowval.dim()){
    rowval.mem_dim=Integer(rowval.memarray->get(rowval.dim()+A.rowval.dim(),rowval.m));
    if (rowindex.mem_dim<rowindex.dim()+A.rowindex.dim()){
      MEmessage(MEmem(rowindex.dim()+A.rowindex.dim(),
		      "Matrix::concat_right(const Matrix&) not enough memory",
		      MTmatrix));
    }
    free_oldvp=1;
  }   
  Indexmatrix tmpinfo(infonr,3); chk_set_init(tmpinfo,1);
  Integer newnz=rowindex.dim()+A.rowindex.dim();
  rowindex.nc=1;rowindex.nr=newnz;
  rowval.nc=1; rowval.nr=newnz;
  thisi=rowinfo.rowdim();
  thisind=-1;
  if (--thisi>=0) thisind=rowinfo(thisi,0);
  Ai=A.rowinfo.rowdim();
  Aind=-1;
  if (--Ai>=0) Aind=A.rowinfo(Ai,0);
  while((Aind>=0)||(thisind>=0)){
    if (thisind>Aind){
      --infonr;
      tmpinfo(infonr,0)=rowinfo(thisi,0);
      tmpinfo(infonr,1)=rowinfo(thisi,1);
      newnz-=rowinfo(thisi,1);
      mat_xey(rowinfo(thisi,1),rowindex.get_store()+newnz+rowinfo(thisi,1)-1,-1,oldip+rowinfo(thisi,2)+rowinfo(thisi,1)-1,-1);
      mat_xey(rowinfo(thisi,1),rowval.get_store()+newnz+rowinfo(thisi,1)-1,-1,oldvp+rowinfo(thisi,2)+rowinfo(thisi,1)-1,-1);
      tmpinfo(infonr,2)=newnz;
      if (--thisi>=0) {
	thisind=rowinfo(thisi,0);
      }
      else {
	thisind=-1;
      }
      continue;
    }
    if (thisind<Aind){
      --infonr;
      tmpinfo(infonr,0)=A.rowinfo(Ai,0);
      tmpinfo(infonr,1)=A.rowinfo(Ai,1);
      newnz-=A.rowinfo(Ai,1);
      mat_xey(A.rowinfo(Ai,1),rowindex.get_store()+newnz,A.rowindex.get_store()+A.rowinfo(Ai,2));
      mat_xpea(A.rowinfo(Ai,1),rowindex.get_store()+newnz,nc);
      mat_xey(A.rowinfo(Ai,1),rowval.get_store()+newnz,A.rowval.get_store()+A.rowinfo(Ai,2));
      tmpinfo(infonr,2)=newnz;
      if (--Ai>=0){
	Aind=A.rowinfo(Ai,0);
      }
      else {
	Aind=-1;
      }
      continue;
    }
    infonr--;
    tmpinfo(infonr,0)=rowinfo(thisi,0);
    tmpinfo(infonr,1)=A.rowinfo(Ai,1);
    newnz-=A.rowinfo(Ai,1);
    mat_xey(A.rowinfo(Ai,1),rowindex.get_store()+newnz,A.rowindex.get_store()+A.rowinfo(Ai,2));
    mat_xpea(A.rowinfo(Ai,1),rowindex.get_store()+newnz,nc);
    mat_xey(A.rowinfo(Ai,1),rowval.get_store()+newnz,A.rowval.get_store()+A.rowinfo(Ai,2));
    tmpinfo(infonr,1)+=rowinfo(thisi,1);
    newnz-=rowinfo(thisi,1);
    mat_xey(rowinfo(thisi,1),rowindex.get_store()+newnz+rowinfo(thisi,1)-1,-1,oldip+rowinfo(thisi,2)+rowinfo(thisi,1)-1,-1);
    mat_xey(rowinfo(thisi,1),rowval.get_store()+newnz+rowinfo(thisi,1)-1,-1,oldvp+rowinfo(thisi,2)+rowinfo(thisi,1)-1,-1);
    tmpinfo(infonr,2)=newnz;
    if (--thisi>=0) {
      thisind=rowinfo(thisi,0);
    }
    else {
      thisind=-1;
    }
    if (--Ai>=0){
      Aind=A.rowinfo(Ai,0);
    }
    else {
      Aind=-1;
    }
  }
  swap(tmpinfo,rowinfo);

  if (free_oldip) {
    rowindex.memarray->free(oldip);
  }
  if (free_oldvp){
    rowval.memarray->free(oldvp);
  }

  nc+=A.nc;

  return *this;
}

Sparsemat& Sparsemat::concat_below(const Sparsemat& A)
{
  chk_init(A);
  chk_init(*this);
  if ((A.nr==0)&&(A.nc==0))
    return *this;
  if ((nr==0)&&(nc==0)){
    *this=A;
    return *this;
  } 
#if (CONICBUNDLE_DEBUG>=1)
  if ((nr*nc!=0)&&(nc!=A.nc)){
     MEmessage(MatrixError(ME_unspec,"Sparsemat::concat_below(Sparsemat& A): number of rows do not match",MTsparse));
  } 
#endif
  if((nc==0)||(A.nonzeros()==0)){
    nr+=A.nr;
    return *this;
  }
  //row representation
  Integer infonr=rowinfo.rowdim();
  Integer oldnz=rowindex.dim();
  rowinfo.concat_below(A.rowinfo);
  rowval.concat_below(A.rowval);
  rowindex.concat_below(A.rowindex);
  for(Integer i=infonr;i<rowinfo.rowdim();i++){
    rowinfo(i,0)+=nr;
    rowinfo(i,2)+=oldnz;
  }

  //column representation
  //determin joint number of columns in infonr
  infonr=0;
  Integer thisi=0;
  Integer thisind=nc;
  if (thisi<colinfo.rowdim()) thisind=colinfo(thisi,0);
  Integer Ai=0;
  Integer Aind=nc;
  if (Ai<A.colinfo.rowdim()) Aind=A.colinfo(Ai,0);
  while((thisind<nc)||(Aind<nc)){
    if (thisind<Aind){
      infonr++;
      thisi++;
      if (thisi<colinfo.rowdim()){
	thisind=colinfo(thisi,0);
      }
      else {
	thisind=nc;
      }
      continue;
    }
    if (thisind>Aind){
      infonr++;
      Ai++;
      if (Ai<A.colinfo.rowdim()){
	Aind=A.colinfo(Ai,0);
      }
      else {
	Aind=nc;
      }
      continue;
    }
    infonr++;
    thisi++;
    if (thisi<colinfo.rowdim()){
      thisind=colinfo(thisi,0);
    }
    else {
      thisind=nc;
    }
    Ai++;
    if (Ai<A.colinfo.rowdim()){
      Aind=A.colinfo(Ai,0);
      }
    else {
      Aind=nc;
    }
  }
  //infonr now holds cardinality of union of column sets
   
  int free_oldip=0;
  Integer *oldip=colindex.get_store();
  if (colindex.mem_dim<colindex.dim()+A.colindex.dim()){
    colindex.mem_dim=Integer(colindex.memarray->get(colindex.dim()+A.colindex.dim(),colindex.m));
    if (colindex.mem_dim<colindex.dim()+A.colindex.dim()){
      MEmessage(MEmem(colindex.dim()+A.colindex.dim(),
		      "Matrix::concat_right(const Matrix&) not enough memory",
		      MTmatrix));
    }
    free_oldip=1;
  } 
  int free_oldvp=0;  
  Real *oldvp=colval.get_store();
  if (colval.mem_dim<colval.dim()+A.colval.dim()){
    colval.mem_dim=Integer(colval.memarray->get(colval.dim()+A.colval.dim(),colval.m));
    if (colindex.mem_dim<colindex.dim()+A.colindex.dim()){
      MEmessage(MEmem(colindex.dim()+A.colindex.dim(),
		      "Matrix::concat_right(const Matrix&) not enough memory",
		      MTmatrix));
    }
    free_oldvp=1;
  }   
  Indexmatrix tmpinfo(infonr,3); chk_set_init(tmpinfo,1);
  Integer newnz=colindex.dim()+A.colindex.dim();
  colindex.nc=1;colindex.nr=newnz;
  colval.nc=1; colval.nr=newnz;
  thisi=colinfo.rowdim();
  thisind=-1;
  if (--thisi>=0) thisind=colinfo(thisi,0);
  Ai=A.colinfo.rowdim();
  Aind=-1;
  if (--Ai>=0) Aind=A.colinfo(Ai,0);
  while((Aind>=0)||(thisind>=0)){
    if (thisind>Aind){
      --infonr;
      tmpinfo(infonr,0)=colinfo(thisi,0);
      tmpinfo(infonr,1)=colinfo(thisi,1);
      newnz-=colinfo(thisi,1);
      mat_xey(colinfo(thisi,1),colindex.get_store()+newnz+colinfo(thisi,1)-1,-1,oldip+colinfo(thisi,2)+colinfo(thisi,1)-1,-1);
      mat_xey(colinfo(thisi,1),colval.get_store()+newnz+colinfo(thisi,1)-1,-1,oldvp+colinfo(thisi,2)+colinfo(thisi,1)-1,-1);
      tmpinfo(infonr,2)=newnz;
      if (--thisi>=0) {
	thisind=colinfo(thisi,0);
      }
      else {
	thisind=-1;
      }
      continue;
    }
    if (thisind<Aind){
      --infonr;
      tmpinfo(infonr,0)=A.colinfo(Ai,0);
      tmpinfo(infonr,1)=A.colinfo(Ai,1);
      newnz-=A.colinfo(Ai,1);
      mat_xey(A.colinfo(Ai,1),colindex.get_store()+newnz,A.colindex.get_store()+A.colinfo(Ai,2));
      mat_xpea(A.colinfo(Ai,1),colindex.get_store()+newnz,nr);
      mat_xey(A.colinfo(Ai,1),colval.get_store()+newnz,A.colval.get_store()+A.colinfo(Ai,2));
      tmpinfo(infonr,2)=newnz;
      if (--Ai>=0){
	Aind=A.colinfo(Ai,0);
      }
      else {
	Aind=-1;
      }
      continue;
    }
    infonr--;
    tmpinfo(infonr,0)=colinfo(thisi,0);
    tmpinfo(infonr,1)=A.colinfo(Ai,1);
    newnz-=A.colinfo(Ai,1);
    mat_xey(A.colinfo(Ai,1),colindex.get_store()+newnz,A.colindex.get_store()+A.colinfo(Ai,2));
    mat_xpea(A.colinfo(Ai,1),colindex.get_store()+newnz,nr);
    mat_xey(A.colinfo(Ai,1),colval.get_store()+newnz,A.colval.get_store()+A.colinfo(Ai,2));
    tmpinfo(infonr,1)+=colinfo(thisi,1);
    newnz-=colinfo(thisi,1);
    mat_xey(colinfo(thisi,1),colindex.get_store()+newnz+colinfo(thisi,1)-1,-1,oldip+colinfo(thisi,2)+colinfo(thisi,1)-1,-1);
    mat_xey(colinfo(thisi,1),colval.get_store()+newnz+colinfo(thisi,1)-1,-1,oldvp+colinfo(thisi,2)+colinfo(thisi,1)-1,-1);
    tmpinfo(infonr,2)=newnz;
    if (--thisi>=0) {
      thisind=colinfo(thisi,0);
    }
    else {
      thisind=-1;
    }
    if (--Ai>=0){
      Aind=A.colinfo(Ai,0);
    }
    else {
      Aind=-1;
    }
  }
  swap(tmpinfo,colinfo);

  if (free_oldip) {
    colindex.memarray->free(oldip);
  }
  if (free_oldvp){
    colval.memarray->free(oldvp);
  }

  nr+=A.nr;

  return *this;
}


Sparsemat operator*(const Sparsemat &A,const Sparsemat &B) 
{
 chk_mult(A,B);
 Sparsemat C;
 C.nr=A.nr;
 C.nc=B.nc;
 C.rowinfo=A.rowinfo;
 C.colinfo=B.colinfo;

 //first run to find all nonzeros
 Integer i,j,nz=0;
 for(i=0;i<B.colinfo.rowdim();i++){
     C.colinfo(i,1)=0;
 }
 for(i=0;i<C.rowinfo.rowdim();i++){
     C.rowinfo(i,1)=0;
     for(j=0;j<C.colinfo.rowdim();j++){
         Real d=0.;
         Integer rk=A.rowinfo(i,2); 
         Integer ru=rk+A.rowinfo(i,1);
         Integer ck=B.colinfo(j,2);
         Integer cu=ck+B.colinfo(j,1);
         while((rk<ru)&&(ck<cu)){
             if (A.rowindex(rk)==B.colindex(ck)){
                 d+=A.rowval(rk)*B.colval(ck);
                 rk++;
                 ck++;
             }
             else if (A.rowindex(rk)<B.colindex(ck)){
                 rk++;
             }
             else {
                 ck++;
             }
         }
         if (abs(d)>C.tol){
             nz++;
             C.rowinfo(i,1)++;
             C.colinfo(j,1)++;
         }
     }
 }

 //initialize index and val to correct size
 if (nz==0) {
     C.colinfo.init(0,0,Integer(0));
     C.rowinfo.init(0,0,Integer(0));
     return C;
 }
 C.colval.newsize(nz,1); chk_set_init(C.colval,1);
 C.rowval.newsize(nz,1); chk_set_init(C.rowval,1);
 C.colindex.newsize(nz,1); chk_set_init(C.colindex,1);
 C.rowindex.newsize(nz,1); chk_set_init(C.rowindex,1);
 
 //second run to fill with correct contents
 nz=0;
 for(i=0;i<C.colinfo.rowdim();i++){
     //set start index for each column and clear column counter
     //so that both col and row can be stored simultaneously afterwards
     C.colinfo(i,2)=nz;
     nz+=C.colinfo(i,1);
     C.colinfo(i,1)=0;
 }
 nz=0;
 for(i=0;i<C.rowinfo.rowdim();i++){
     C.rowinfo(i,2)=nz;
     for(j=0;j<C.colinfo.rowdim();j++){
         Real d=0.;
         Integer rk=A.rowinfo(i,2); 
         Integer ru=rk+A.rowinfo(i,1);
         Integer ck=B.colinfo(j,2);
         Integer cu=ck+B.colinfo(j,1);
         while((rk<ru)&&(ck<cu)){
             if (A.rowindex(rk)==B.colindex(ck)){
                 d+=A.rowval(rk)*B.colval(ck);
                 rk++;
                 ck++;
             }
             else if (A.rowindex(rk)<B.colindex(ck)){
                 rk++;
             }
             else {
                 ck++;
             }
         }
         if (abs(d)>C.tol){
             C.rowindex(nz)=B.colinfo(j,0);
             C.rowval(nz)=d;
             C.colindex(C.colinfo(j,2)+C.colinfo(j,1))=A.rowinfo(i,0);
             C.colval(C.colinfo(j,2)+C.colinfo(j,1))=d;
             C.colinfo(j,1)++;
             nz++;
         }
     }
 }

 //check for empty columns or rows in info
 Indexmatrix delind=(C.rowinfo.col(1)).find_number(Integer(0));
 C.rowinfo.delete_rows(delind);
 delind=(C.colinfo.col(1)).find_number(Integer(0));
 C.colinfo.delete_rows(delind);

 return C;
}
    
void Sparsemat::display(std::ostream& out,
                  int precision,int width,int screenwidth) const
           //for variables of value zero default values are used
           //precision=4,width=precision+6,screenwidth=80
{
 chk_init(*this);
 if (precision==0) precision=4;
 if (width==0) width=precision+6;
 if (screenwidth==0) screenwidth=80;
 out.precision(precision);
 out<<"Sparsemat:"<<nr<<" "<<nc<<" "<<rowindex.dim()<<"\n";
 Integer i,j,nz=0;
 for(i=0;i<rowinfo.rowdim();i++){
     Integer rind=rowinfo(i,0);
     for(j=0;j<rowinfo(i,1);j++){
         out.width(5);out<<rind<<" ";
         out.width(5);out<<rowindex(nz)<<" ";
         out.width(width);out<<rowval(nz)<<"\n";
         nz++;
     }
 }
}

// *****************************************************************************
//                  Constructors from Sparsemat to other classes
// *****************************************************************************


// ******************************    Matrix   ******************************

Matrix& Matrix::xeya(const Sparsemat &A,Real d)
{
 chk_init(A);
 init(A.nr,A.nc,0.);
 if (d==0.) return *this;
 if (d==1.){
     const Real *avp=A.colval.get_store();
     const Integer *aip=A.colindex.get_store();
     for(Integer i=0;i<A.colinfo.rowdim();i++){
         Real *mp=m+A.colinfo(i,0)*nr;
         for(Integer j=A.colinfo(i,1);--j>=0;){
             *(mp+*aip++)=(*avp++);
         }
     }
     return *this;
 }
 if (d==-1.){
     const Real *avp=A.colval.get_store();
     const Integer *aip=A.colindex.get_store();
     for(Integer i=0;i<A.colinfo.rowdim();i++){
         Real *mp=m+A.colinfo(i,0)*nr;
         for(Integer j=A.colinfo(i,1);--j>=0;){
             *(mp+*aip++)=-(*avp++);
         }
     }
     return *this;
 }
 const Real *avp=A.colval.get_store();
 const Integer *aip=A.colindex.get_store();
 for(Integer i=0;i<A.colinfo.rowdim();i++){
     Real *mp=m+A.colinfo(i,0)*nr;
     for(Integer j=A.colinfo(i,1);--j>=0;){
         *(mp+*aip++)=d*(*avp++);
     }
 }
 return *this;    
}
 
Matrix& Matrix::xpeya(const Sparsemat &A,Real d)
{
 chk_add(*this,A);
 if (d==0.) return *this;
 if (d==1.){
     const Real *avp=A.colval.get_store();
     const Integer *aip=A.colindex.get_store();
     for(Integer i=0;i<A.colinfo.rowdim();i++){
         Real *mp=m+A.colinfo(i,0)*nr;
         for(Integer j=A.colinfo(i,1);--j>=0;){
             *(mp+*aip++)+=(*avp++);
         }
     }
     return *this;
 }
 if (d==-1.){
     const Real *avp=A.colval.get_store();
     const Integer *aip=A.colindex.get_store();
     for(Integer i=0;i<A.colinfo.rowdim();i++){
         Real *mp=m+A.colinfo(i,0)*nr;
         for(Integer j=A.colinfo(i,1);--j>=0;){
             *(mp+*aip++)-=(*avp++);
         }
     }
     return *this;
 }
 const Real *avp=A.colval.get_store();
 const Integer *aip=A.colindex.get_store();
 for(Integer i=0;i<A.colinfo.rowdim();i++){
     Real *mp=m+A.colinfo(i,0)*nr;
     for(Integer j=A.colinfo(i,1);--j>=0;){
         *(mp+*aip++)+=d*(*avp++);
     }
 }
 return *this;    
}
  

Sparsemat& Sparsemat::scale_rows(const Matrix& vec)
{
  chk_init(*this);
  chk_init(vec);
  
  Integer ncnt=0;
  for (Integer i=0;i<rowinfo.rowdim();i++){
    Real d=vec(rowinfo(i,0));
    for(Integer j=rowinfo(i,1);--j>=0;ncnt++){
      rowval(ncnt)*=d;
    }
  }
  {for (Integer ncnt=0;ncnt<colval.dim();ncnt++){
    colval(ncnt)*=vec(colindex(ncnt));
  }}
  return *this;
}
  

Sparsemat& Sparsemat::scale_cols(const Matrix& vec)
{
  chk_init(*this);
  chk_init(vec);
  
  Integer ncnt=0;
  for (Integer i=0;i<colinfo.rowdim();i++){
    Real d=vec(colinfo(i,0));
    for(Integer j=colinfo(i,1);--j>=0;ncnt++){
      colval(ncnt)*=d;
    }
  }
  {for (Integer ncnt=0;ncnt<rowval.dim();ncnt++){
    rowval(ncnt)*=vec(rowindex(ncnt));
  }}
  return *this;
}

// ******************************************************************************
//                    friends
// ******************************************************************************

void swap(Sparsemat& A, Sparsemat& B)
{
 Integer h=A.nr;A.nr=B.nr;B.nr=h;
 h=A.nc;A.nr=B.nc;B.nc=h;
 swap(A.colinfo,B.colinfo);
 swap(A.colindex,B.colindex);
 swap(A.colval,B.colval);
 swap(A.rowinfo,B.rowinfo);
 swap(A.rowindex,B.rowindex);
 swap(A.rowval,B.rowval);
 Real d=A.tol;A.tol=B.tol;B.tol=d;
#if (CONICBUNDLE_DEBUG>=1)
 bool bo=A.is_init;A.is_init=B.is_init;B.is_init=bo;
#endif
}
    

Real trace(const Sparsemat& A)
{
 chk_init(A);
 Integer i;
 Real s=0.;
 if (A.rowinfo.rowdim()<=A.colinfo.rowdim()){
     for(i=0;i<A.rowinfo.rowdim();i++){
         Integer rind=A.rowinfo(i,0);
         if (rind>=A.nc) return s;
         Integer li=A.rowinfo(i,2);
         Integer ui=li+A.rowinfo(i,1)-1;
         Integer mi=0;
         while(li<=ui){
             mi=(li+ui)/2;
             if (A.rowindex(mi)==rind) break;
             if (A.rowindex(mi)<rind){
                 li=mi+1;
             }
             else {
                 ui=mi-1;
             }
         } ;
         if (li<=ui) s+=A.rowval(mi);
     }         
 }
 else { 
     for(i=0;i<A.colinfo.rowdim();i++){
         Integer cind=A.colinfo(i,0);
         if (cind>=A.nr) return s;
         Integer li=A.colinfo(i,2);
         Integer ui=li+A.colinfo(i,1)-1;
         Integer mi=0;
         while(li<=ui){
             mi=(li+ui)/2;
             if (A.colindex(mi)==cind) break;
             if (A.colindex(mi)<cind){
                 li=mi+1;
             }
             else {
                 ui=mi-1;
             }
         }
         if (li<=ui) s+=A.colval(mi);
     }         
 }
 return s;
}
    
Real ip(const Sparsemat& A, const Sparsemat& B)
{
 chk_add(A,B);
 if (&A==&B)
     return mat_ip(A.colval.dim(),A.colval.get_store(),A.colval.get_store());
 Integer ai=0,bi=0;
 Real s=0.;
 while((ai<A.rowinfo.rowdim())&&(bi<B.rowinfo.rowdim())){
     if (A.rowinfo(ai,0)==B.rowinfo(bi,0)){
         Integer aj=A.rowinfo(ai,2);
         Integer au=aj+A.rowinfo(ai,1);
         Integer bj=B.rowinfo(bi,2);
         Integer bu=bj+B.rowinfo(bi,1);
         while((aj<au)&&(bj<bu)){
             if(A.rowindex(aj)==B.rowindex(bj)){
                 s+=A.rowval(aj)*B.rowval(bj);
                 aj++;
                 bj++;
             }
             else if (A.rowindex(aj)<B.rowindex(bj)){
                 aj++;
             }
             else {
                 bj++;
             }
         }
         ai++;
         bi++;
     }
     else if (A.rowinfo(ai,0)<B.rowinfo(bi,0)){
         ai++;
     }
     else {
         bi++;
     }
 }
 return s;
}
             
Real ip(const Sparsemat& A, const Matrix& B)
{
 chk_add(A,B);
 Real s=0.;
 const Real *avp=A.colval.get_store();
 const Integer *aip=A.colindex.get_store();
 for(Integer i=0;i<A.colinfo.rowdim();i++){
     const Real *bbp=B.get_store()+A.colinfo(i,0)*A.rowdim();
     for(Integer j=A.colinfo(i,1);--j>=0;){
         s+=*(bbp+*aip++)*(*avp++);
     }
 }
 return s;
}
             
Matrix sumrows(const Sparsemat& A)   //=(1 1 1 ... 1)*A
{
  chk_init(A);
  Matrix sr(1,A.nc,0.);
  for(Integer j=0;j<A.colinfo.rowdim();j++){
    sr(A.colinfo(j,0))=mat_sum(A.colinfo(j,1),A.colval.get_store()+A.colinfo(j,2));
  }
  return sr;
}

Matrix sumcols(const Sparsemat& A)   //=A*(1 1 ... 1)^t
{
  chk_init(A);
  Matrix sc(A.nr,1,0.);
  for(Integer j=0;j<A.rowinfo.rowdim();j++){
    sc(A.rowinfo(j,0))=mat_sum(A.rowinfo(j,1),A.rowval.get_store()+A.rowinfo(j,2));
  }
  return sc;
}
 
Sparsemat abs(const Sparsemat& A)
{
 chk_init(A);
 Sparsemat B(A.nr,A.nc);
 B.rowinfo=A.rowinfo;
 B.colinfo=A.colinfo;
 B.rowindex=A.rowindex;
 B.colindex=A.colindex;
 B.rowval.newsize(A.rowval.dim(),1); chk_set_init(B.rowval,1);
 B.colval.newsize(A.colval.dim(),1); chk_set_init(B.colval,1);
 Real *br=B.rowval.get_store();
 const Real *ar=A.rowval.get_store();
 Real *bc=B.colval.get_store();
 const Real *ac=A.colval.get_store();
 Integer i=A.rowval.dim();
 for(;--i>=0;){
     *br++= abs(*ar++);
     *bc++= abs(*ac++);
 }
 return B;
}

int equal(const Sparsemat& A, const Sparsemat& B,Real eqtol)
{
  chk_init(A);
  chk_init(B);
  if ((A.nr!=B.nr)||(A.nc!=B.nc)||
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


std::ostream& operator<<(std::ostream& o,const Sparsemat &A)
{
 chk_init(A);
 o<<A.nr<<" "<<A.nc<<" "<<A.rowval.dim()<<"\n";
 Integer i,j,nz=0;
 for(i=0;i<A.rowinfo.rowdim();i++){
     Integer rind=A.rowinfo(i,0);
     for(j=0;j<A.rowinfo(i,1);j++){
         o<<rind<<" ";
         o<<A.rowindex(nz)<<" ";
         o<<A.rowval(nz)<<"\n";
         nz++;
     }
 }
 return o;
}

std::istream& operator>>(std::istream& in,Sparsemat &A)
{
 const char* format="         format: nrows ncols nz i_1 j_1 val_1 ... i_nz j_nz val_nz\n         nrows>=0, ncols>=0, nz>=0, 0<=i<nrows, 0<=j<ncols "; 
 Integer nc,nr,nz;
 if (!(in>>nr)){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
     if (materrout) (*materrout)<<" failed in reading number of rows"<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if(nr<0){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
     if (materrout) (*materrout)<<" number of rows must be nonnegative but is "<<nr<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if (!(in>>nc)){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
     if (materrout) (*materrout)<<" failed in reading number of columns"<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if(nc<0){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
     if (materrout) (*materrout)<<" number of columns must be nonnegative but is "<<nc<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if (!(in>>nz)){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
     if (materrout) (*materrout)<<" failed in reading number of nonzero elements"<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if(nz<0){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
     if (materrout) (*materrout)<<" number of nonzeros must be nonnegative but is "<<nz<<std::endl;
     if (materrout) (*materrout)<<format;
     in.clear(std::ios::failbit);
     return in;
 }
 if(((nr==0)||(nc==0))&&(nz>0)){
     if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
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
         if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
         if (materrout) (*materrout)<<" failed in reading nonzero element (i,j,val) #";
         if (materrout) (*materrout)<<i<<std::endl;
         in.clear(std::ios::failbit);
         return in;
     }
     if ((0>I(i))||(I(i)>=nr)){
         if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
         if (materrout) (*materrout)<<" row index of nonzero element #"<<i<<"exceeds range: ";
         if (materrout) (*materrout)<<0<<"<="<<I(i)<<"<"<<nr<<std::endl;
         in.clear(std::ios::failbit);
         return in;
     }
     if ((0>J(i))||(J(i)>=nc)){
         if (materrout) (*materrout)<<"*** ERROR: operator>>(std::istream&,Sparsemat&): ";
         if (materrout) (*materrout)<<" column index of nonzero element #"<<i<<"exceeds range: ";
         if (materrout) (*materrout)<<0<<"<="<<J(i)<<"<"<<nc<<std::endl;
         in.clear(std::ios::failbit);
         return in;
     }
 }
 A.init(nr,nc,nz,I,J,val);
 return in;
}



// **************************************************************************
//                     specific Indexmatrix 
// **************************************************************************

Indexmatrix::Indexmatrix(const Sparsemat &A)
{
  init_to_zero();
 chk_init(A);
 init(A.nr,A.nc,Integer(0));
 for(Integer i=0;i<A.colinfo.rowdim();i++){
   Integer *mp=m+A.colinfo(i,0)*nr;
   const Real *amp=A.colval.get_store()+A.colinfo(i,2);
   const Integer *aip=A.colindex.get_store()+A.colinfo(i,2);
   for(Integer j=A.colinfo(i,1);--j>=0;++amp){
     *(mp+(*aip++))=Integer((*amp>0)? *amp+.5 : *amp-.5);
   }
 }
}

// ******************************   Symmatrix   ******************************


Symmatrix& rankadd(const Sparsemat& A,Symmatrix& C,
                    Real alpha,Real beta,int trans)
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
     const Integer *ip=A.get_rowindex().get_store()-1;
     const Real *vp=A.get_rowval().get_store()-1;
     for(Integer j=0;j<A.rowinfo.rowdim();j++){
         for(Integer i=A.rowinfo(j,1);--i>=0;){
             const Real *vvp;
             Real d1=*(vvp=++vp);
             Real d=d1*alpha;
             const Integer *iip;
             Integer ii=*(iip=++ip);
             Real *mp=C.get_store()+ii*nr-(ii*(ii+1))/2;
             *(mp+ii)+=d*d1;
             for(Integer k=i;--k>=0;){
                 *(mp+(*(++iip)))+=d*(*(++vvp));
             }
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
     const Integer *ip=A.get_colindex().get_store()-1;
     const Real *vp=A.get_colval().get_store()-1;
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         for(Integer i=A.colinfo(j,1);--i>=0;){
             const Real *vvp;
             Real d1=*(vvp=++vp);
             Real d=d1*alpha;
             const Integer *iip;
             Integer ii=*(iip=++ip);
             Real *mp=C.get_store()+ii*nr-(ii*(ii+1))/2;
             *(mp+ii)+=d*d1;
             for(Integer k=i;--k>=0;){
                 *(mp+(*(++iip)))+=d*(*(++vvp));
             }
         }
     }
     return C;
 }
}
 

Symmatrix& rank2add(const Sparsemat& A,const Matrix& B,Symmatrix& C,
                    Real alpha,Real beta,int trans)
{
 chk_add(A,B);
 if (trans) { //A^TB
     if (beta==0.) {
         C.init(A.coldim(),0.);
     }
     else {
         chk_colseq(C,A);
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer n=A.rowdim();
     Integer nr=C.rowdim();
     const Integer *ip=A.get_colindex().get_store();
     const Real *vp=A.get_colval().get_store();
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         for(Integer i=A.colinfo(j,1);--i>=0;){
             Real d=alpha/2.*(*vp++);
             Real *mp=C.get_store()+jj;
             const Real *bp=B.get_store()+*ip++-n;
             Integer k=0;
             for(;++k<=jj;){
                 *mp += d* *(bp+=n);
                 mp+=nr-k;
             }
             (*mp++) += 2.*d* *(bp+=n);
             for(;++k<=nr;){
                 (*mp++) += d* *(bp+=n);
             }                 
         }
     }
     return C;
 }
 else { // AB^T 
     if (beta==0.) {
         C.init(A.rowdim(),0.);
     }
     else {
         chk_rowseq(C,A);
         if (beta!=1.) C*=beta;
     }
     if (alpha==0.) return C;
     Integer n=A.rowdim();
     Integer nr=C.rowdim();
     const Integer *jp=A.get_rowindex().get_store();
     const Real *vp=A.get_rowval().get_store();
     for(Integer i=0;i<A.rowinfo.rowdim();i++){
         Integer ii=A.rowinfo(i,0);
         for(Integer j=A.rowinfo(i,1);--j>=0;){
             Real d=alpha/2.*(*vp++);
             Real *mp=C.get_store()+ii;
             const Real *bp=B.get_store()+(*jp++)*n;
             Integer k=0;
             for(;++k<=ii;){
                 *mp += d* *(bp++);
                 mp+=nr-k;
             }
             (*mp++) += 2.*d* *(bp++);
             for(;++k<=nr;){
                 (*mp++) += d* *(bp++);
             }                 
         }
     }
     return C;
 }
}     

Symmatrix& Symmatrix::xetriu_yza(const Sparsemat& A,const Matrix& B,Real d)
{
 chk_add(A,B);
 init(A.coldim(),0.);
 chk_set_init(*this,1);
 if (d==0.) { return *this;}
 if (d==1.) {
     Integer n=A.rowdim();
     const Integer *ip=A.get_colindex().get_store();
     const Real *vp=A.get_colval().get_store();
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         Real *mp=m+jj*nr-(jj*(jj-1))/2;
         const Real *bp=B.get_store()+jj*n;
         for(Integer i=A.colinfo(j,1);--i>=0;){
             mat_xpeya(nr-jj,mp,1,bp+*ip++,n,*vp++);
         }
     }
     return *this;
 }
 if (d==-1.) {
     Integer n=A.rowdim();
     const Integer *ip=A.get_colindex().get_store();
     const Real *vp=A.get_colval().get_store();
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         Real *mp=m+jj*nr-(jj*(jj-1))/2;
         const Real *bp=B.get_store()+jj*n;
         for(Integer i=A.colinfo(j,1);--i>=0;){
             mat_xpeya(nr-jj,mp,1,bp+*ip++,n,-*vp++);
         }
     }
     return *this;
 }
 Integer n=A.rowdim();
 const Integer *ip=A.get_colindex().get_store();
 const Real *vp=A.get_colval().get_store();
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     Real *mp=m+jj*nr-(jj*(jj-1))/2;
     const Real *bp=B.get_store()+jj*n;
     for(Integer i=A.colinfo(j,1);--i>=0;){
         mat_xpeya(nr-jj,mp,1,bp+*ip++,n,d*(*vp++));
     }
 }
 return *this;
}
         

// *this is the upper triangle of the matrix product A^T*B
Symmatrix& Symmatrix::xpetriu_yza(const Sparsemat& A,const Matrix& B,Real d)
{
 chk_add(A,B);
 chk_colseq(*this,A);
 if (d==0.) { return *this;}
 if (d==1.) {
     Integer n=A.rowdim();
     const Integer *ip=A.get_colindex().get_store();
     const Real *vp=A.get_colval().get_store();
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         Real *mp=m+jj*nr-(jj*(jj-1))/2;
         const Real *bp=B.get_store()+jj*n;
         for(Integer i=A.colinfo(j,1);--i>=0;){
             mat_xpeya(nr-jj,mp,1,bp+*ip++,n,*vp++);
         }
     }
     return *this;
 }
 if (d==-1.) {
     Integer n=A.rowdim();
     const Integer *ip=A.get_colindex().get_store();
     const Real *vp=A.get_colval().get_store();
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         Real *mp=m+jj*nr-(jj*(jj-1))/2;
         const Real *bp=B.get_store()+jj*n;
         for(Integer i=A.colinfo(j,1);--i>=0;){
             mat_xpeya(nr-jj,mp,1,bp+*ip++,n,-(*vp++));
         }
     }
     return *this;
 }
 Integer n=A.rowdim();
 const Integer *ip=A.get_colindex().get_store();
 const Real *vp=A.get_colval().get_store();
 for(Integer j=0;j<A.colinfo.rowdim();j++){
     Integer jj=A.colinfo(j,0);
     Real *mp=m+jj*nr-(jj*(jj-1))/2;
     const Real *bp=B.get_store()+jj*n;
     for(Integer i=A.colinfo(j,1);--i>=0;){
         mat_xpeya(nr-jj,mp,1,bp+*ip++,n,d*(*vp++));
     }
 }
 return *this;
}

Matrix& genmult(const Symmatrix& A,const Sparsemat& B,Matrix& C,
                    Real alpha,Real beta,int btrans)
            //returns C=beta*C+alpha*A*B, where A and B may be transposed
{
 chk_init(A);
 chk_init(B);
 Integer nr,nc;
 nr=A.nr;
#if (CONICBUNDLE_DEBUG>=1)
 Integer nm=A.nr;
#endif
 if (btrans) {
     nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nc) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsymmetric));;
     }
#endif
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsymmetric));;
     }
#endif
 }
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.rowdim())||(nc!=C.coldim())) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsymmetric));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if (btrans){
     const Integer *bip=B.rowindex.get_store();
     const Real *bvp=B.rowval.get_store();
     for(Integer j=0;j<B.rowinfo.rowdim();j++){
         Integer jj=B.rowinfo(j,0);
         for(Integer i=0;i<B.rowinfo(j,1);i++){
             Real d=alpha*(*bvp++);
             Integer ii=*bip++;
             const Real* ap=A.get_store()+ii;
             Real* cp=C.get_store()+jj*nr;
             Integer k=0;
             for(;k<ii;){
                 *cp++ +=d*(*ap);
                 ap+=nr- ++k;
             }
             for(;k++<nr;){
                 *cp++ +=d*(*ap++);
             }
         }
     }
 }
 else {
     const Integer *bip=B.colindex.get_store();
     const Real *bvp=B.colval.get_store();
     for(Integer j=0;j<B.colinfo.rowdim();j++){
         Integer jj=B.colinfo(j,0);
         for(Integer i=0;i<B.colinfo(j,1);i++){
             Real d=alpha*(*bvp++);
             Integer ii=*bip++;
             const Real* ap=A.get_store()+ii;
             Real* cp=C.get_store()+jj*nr;
             Integer k=0;
             for(;k<ii;){
                 *cp++ +=d*(*ap);
                 ap+=nr- ++k;
             }
             for(;k++<nr;){
                 *cp++ +=d*(*ap++);
             }
         }
     }
 }
 return C;
}
         
Matrix& genmult(const Sparsemat& A,const Symmatrix& B,Matrix& C,
                    Real alpha,Real beta,int atrans)
            //returns C=beta*C+alpha*A*B, where A and B may be transposed
{
 chk_init(A);
 chk_init(B);
 Integer nr,nc;
#if (CONICBUNDLE_DEBUG>=1)
 Integer nm;
#endif
 if (atrans){
     nr=A.nc;
#if (CONICBUNDLE_DEBUG>=1)
     nm=A.nr;
#endif
 }
 else {
     nr=A.nr;
#if (CONICBUNDLE_DEBUG>=1)
     nm=A.nc;
#endif
 }
 nc=B.nr;
#if (CONICBUNDLE_DEBUG>=1)
 if (nm!=B.nr) {
     MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsymmetric));;
 }
#endif
 if (beta!=0.){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.rowdim())||(nc!=C.coldim())) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTsymmetric));;
     }
#endif
     if (beta!=1.) C*=beta;
 }
 else {
     C.init(nr,nc,0.);
 }
 if (alpha==0.) return C;
 if (atrans){
     const Integer *aip=A.colindex.get_store();
     const Real *avp=A.colval.get_store();
     for(Integer j=0;j<A.colinfo.rowdim();j++){
         Integer jj=A.colinfo(j,0);
         for(Integer i=0;i<A.colinfo(j,1);i++){
             Real d=alpha*(*avp++);
             Integer ii=*aip++;
             const Real* bp=B.get_store()+ii;
             Real* cp=C.get_store()+jj;
             Integer k=0;
             for(;k<ii;){
                 *cp +=d*(*bp);
                 cp+=nr;
                 bp+=nr- ++k;
             }
             for(;k++<nr;){
                 *cp +=d*(*bp++);
                 cp+=nr;
             }
         }
     }
 }
 else {
     const Integer *aip=A.rowindex.get_store();
     const Real *avp=A.rowval.get_store();
     for(Integer j=0;j<A.rowinfo.rowdim();j++){
         Integer jj=A.rowinfo(j,0);
         for(Integer i=0;i<A.rowinfo(j,1);i++){
             Real d=alpha*(*avp++);
             Integer ii=*aip++;
             const Real* bp=B.get_store()+ii;
             Real* cp=C.get_store()+jj;
             Integer k=0;
             for(;k<ii;){
                 *cp +=d*(*bp);
                 cp+=nr;
                 bp+=nr- ++k;
             }
             for(;k++<nr;){
                 *cp +=d*(*bp++);
                 cp+=nr;
             }
         }
     }
 }
 return C;
}

}
 
