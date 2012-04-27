/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/indexmat.cxx

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
#include <math.h>
#include <algorithm>
#include "heapsort.hxx"
#include "indexmat.hxx"
#include "mymath.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

const Mtype Indexmatrix::mtype = MTindexmatrix;

// **************************************************************************
//                                BLAS-like routines
// **************************************************************************

Indexmatrix& Indexmatrix::xeya(const Indexmatrix& A,Integer d)
{
 chk_init(A);
 newsize(A.nr,A.nc);
 chk_set_init(*this,1);
 if (d==1) { mat_xey(nr*nc,m,A.m); return *this;}
 if (d==0) { mat_xea(nr*nc,m,Integer(0)); return *this;}
 if (d==-1) { mat_xemy(nr*nc,m,A.m); return *this;}
 mat_xeya(nr*nc,m,A.m,d);
 return *this;
}

Indexmatrix& Indexmatrix::xpeya(const Indexmatrix& A,Integer d)
{
 chk_add(*this,A);
 if (d==1) { mat_xpey(nr*nc,m,A.m); return *this;}
 if (d==0) { return *this;}
 if (d==-1) { mat_xmey(nr*nc,m,A.m); return *this;}
 mat_xpeya(nr*nc,m,A.m,d);
 return *this;
}

Indexmatrix& xbpeya(Indexmatrix& x,const Indexmatrix& y,Integer alpha,Integer beta,int ytrans)
  //returns x= alpha*y+beta*x, where y may be transposed (ytrans=1)
  //if beta==0 then x is initialized to the correct size
{
  chk_init(y);
  if (beta==0){
    if (!ytrans){ //y is not transposed
      x.newsize(y.nr,y.nc);
      if (alpha==0){
	mat_xea(x.dim(),x.m,Integer(0));
      }
      else if (alpha==1){
	mat_xey(y.dim(),x.m,y.m);
      }
      else {
	mat_xeya(y.dim(),x.m,y.m,alpha);
      }
      chk_set_init(x,1);
      return x;
    }
    else{ //y is transposed
      x.newsize(y.nc,y.nr);
      if (alpha==0){
	mat_xea(x.dim(),x.m,Integer(0));
      }
      else if (alpha==1){
	if (x.nr>x.nc){
	  for (Integer i=0;i<x.nc;i++){
	    mat_xey(x.nr,x.m+i*x.nr,1,y.m+i,x.nc);
	  }
	} 
	else {
	  for (Integer i=0;i<x.nr;i++){
	    mat_xey(x.nc,x.m+i,x.nr,y.m+i*x.nc,1);
	  }
	}
      }
      else {
	if (x.nr>x.nc){
	  for (Integer i=0;i<x.nc;i++){
	    mat_xeya(x.nr,x.m+i*x.nr,1,y.m+i,x.nc,alpha);
	  }
	} 
	else {
	  for (Integer i=0;i<x.nr;i++){
	    mat_xeya(x.nc,x.m+i,x.nr,y.m+i*x.nc,1,alpha);
	  }
	}
      }
      chk_set_init(x,1);
      return x;
    }
  }
  //Now beta!=0
  if (!ytrans) { //y is not transposed
    chk_add(x,y);
    if (alpha==0){
      if (beta==1) return x;
      else if (beta==-1) {
	mat_xemx(x.dim(),x.m);
        return x;
      }
      mat_xmultea(x.dim(),x.m,beta);
      return x;   
    }
    else if (beta==1){
      if (alpha==1){
	mat_xpey(x.dim(),x.m,y.m);
        return x;
      }
      else if (alpha==-1){
	mat_xmey(x.dim(),x.m,y.m);
        return x;
      }
      mat_xpeya(x.dim(),x.m,y.m,alpha);
      return x; 
    }
  } 
  //Now y transposed
  chk_init(x);
#if (CONICBUNDLE_DEBUG>=1)
  if ((x.nr!=y.nc)||(x.nc!=y.nr)){
    MEmessage(MatrixError(ME_dim,"xbpeya: dimensions don't match",MTindexmatrix));;
  }
#endif 
  if (alpha==0){
    if (beta==1) return x;
    if (beta==-1) {
      mat_xemx(x.dim(),x.m);
      return x;
    }
    mat_xmultea(x.dim(),x.m,beta);
    return x;   
  }
  if (beta==1){
    if (alpha==1){
      if (x.nr>x.nc){
	for (Integer i=0;i<x.nc;i++){
	  mat_xpey(x.nr,x.m+i*x.nr,1,y.m+i,x.nc);
	}
      } 
      else {
	for (Integer i=0;i<x.nr;i++){
	  mat_xpey(x.nc,x.m+i,x.nr,y.m+i*x.nc,1);
	}
      }
      return x;
    }
    if (alpha==-1){
      if (x.nr>x.nc){
	for (Integer i=0;i<x.nc;i++){
	  mat_xmey(x.nr,x.m+i*x.nr,1,y.m+i,x.nc);
	}
      } 
      else {
	for (Integer i=0;i<x.nr;i++){
	  mat_xmey(x.nc,x.m+i,x.nr,y.m+i*x.nc,1);
	}
      }
      return x;
    }
      
    if (x.nr>x.nc){
      for (Integer i=0;i<x.nc;i++){
	mat_xpeya(x.nr,x.m+i*x.nr,1,y.m+i,x.nc,alpha);
      }
    } 
    else {
      for (Integer i=0;i<x.nr;i++){
	mat_xpeya(x.nc,x.m+i,x.nr,y.m+i*x.nc,1,alpha);
      }
    }
    return x;
  }
  //now beta!=1 
  if (x.nr>x.nc){
    for (Integer i=0;i<x.nc;i++){
      mat_xbpeya(x.nr,x.m+i*x.nr,1,y.m+i,x.nc,alpha,beta);
    }
  } 
  else {
    for (Integer i=0;i<x.nr;i++){
      mat_xbpeya(x.nc,x.m+i,x.nr,y.m+i*x.nc,1,alpha,beta);
    }
  }
  return x;
}
  
Indexmatrix& xeyapzb(Indexmatrix& x,const Indexmatrix& y,const Indexmatrix& z,Integer alpha,Integer beta)
  //returns x= alpha*y+beta*z,
  //x is initialized to the correct size
{
  chk_add(y,z);
  x.newsize(y.nr,y.nc);
  chk_set_init(x,1);
  if (alpha==1){
    if(beta==1){
      mat_xeypz(x.dim(),x.m,y.m,z.m);
      return x;
    }
    if (beta==-1){
      mat_xeymz(x.dim(),x.m,y.m,z.m);
      return x;
    }
  }
  if ((beta==1)&&(alpha==-1)){
    mat_xeymz(x.dim(),x.m,z.m,y.m);
    return x;
  }
  mat_xeyapzb(x.dim(),x.m,y.m,z.m,alpha,beta);
  return x;
}

Indexmatrix& genmult(const Indexmatrix& A,const Indexmatrix& B,Indexmatrix& C,
                    Integer alpha,Integer beta,int atrans,int btrans)
            //returns C=beta*C+alpha*A*B, where A and B may be transposed
            //C may neither be equal to A nor B
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
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTindexmatrix));;
     }
#endif
 }
 else {
     nc=B.nc;
#if (CONICBUNDLE_DEBUG>=1)
     if (nm!=B.nr) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTindexmatrix));;
     }
#endif
 }
 if (beta!=0){
     chk_init(C);
#if (CONICBUNDLE_DEBUG>=1)
     if ((nr!=C.nr)||(nc!=C.nc)) {
         MEmessage(MatrixError(ME_dim,"genmult: dimensions don't match",MTindexmatrix));;
     }
#endif
     if (beta!=1) C*=beta;
 }
 else {
     C.init(nr,nc,Integer(0));
 }
 if (alpha==0) return C;
 if (atrans){
     if (btrans){
         Integer *cp=C.m;
         for(Integer j=0;j<nc;j++){
             for(Integer i=0;i<nr;i++){
                 *cp++ +=alpha*mat_ip(nm,A.m+i*nm,1,B.m+j,nc);
             }
         }
     }
     else {
         Integer *cp=C.m;
         for(Integer j=0;j<nc;j++){
             for(Integer i=0;i<nr;i++){
                 *cp++ +=alpha*mat_ip(nm,A.m+i*nm,B.m+j*nm);
             }
         }
     }
 }
 else {
     if (btrans){
         Integer *cp=C.m;
         for(Integer j=0;j<nc;j++){
             for(Integer i=0;i<nr;i++){
                 *cp++ +=alpha*mat_ip(nm,A.m+i,nr,B.m+j,nc);
             }
         }
     }
     else {
         Integer *cp=C.m;
         for(Integer j=0;j<nc;j++){
             for(Integer i=0;i<nr;i++){
                 *cp++ +=alpha*mat_ip(nm,A.m+i,nr,B.m+j*nm,1);
             }
         }
     }
 }
 return C;
}

// **************************************************************************
//                                Constructors
// **************************************************************************

Indexmatrix& Indexmatrix::init(const Range &r)
{
 Integer i;
 if (r.step>=0) {
     if (r.from>r.to) i=0;
     else i=(r.to-r.from)/r.step+1;
 }
 else {
     if (r.from<r.to) i=0;
     else i=(r.to-r.from)/r.step+1;
 }
 newsize(i,1);
 Integer d=r.from;
 Integer *mp=m;
 for(i=nr;--i>=0;){
     (*mp++)=d;
     d+=r.step;
 }
 chk_set_init(*this,1);
 return *this;
}

void Indexmatrix::newsize(Integer inr,Integer inc)
{
 chk_range(inr,inc,-1,-1);
 chk_set_init(*this,0);
 if ((inr==0)||(inc==0)){
     nr=inr;
     nc=inc;
     return;
 }
 if ((inr!=nr)||(inc!=nc)) {
     nr=inr;
     nc=inc;
     if (nr*nc>mem_dim){
         memarray->free(m); m=0;
         mem_dim=Integer(memarray->get(nr*nc,m));
         if (mem_dim<nr*nc)
             MEmessage(MEmem(nr*nc,
                         "Indexmatrix::Indexmatrix(Integer,Integer,Integer) not enough memory",MTindexmatrix));
     }
 }
}

Indexmatrix Indexmatrix::find(void) const
{
 chk_init(*this);
 Indexmatrix ind(nr*nc,1);
 chk_set_init(ind,1);
 Integer i,k=0;
 Integer *mp=m;
 for(i=0;i<nr*nc;i++){
     if (*mp++!=0) ind(k++)=i;
 }
 ind.nr=k;
 return ind;
}

Indexmatrix Indexmatrix::find_number(Integer num) const
{
 chk_init(*this);
 Indexmatrix ind(nr*nc,1);
 chk_set_init(ind,1);
 Integer i,k=0;
 Integer *mp=m;
 for(i=0;i<nr*nc;i++){
     if (*mp++==num) ind(k++)=i;
 }
 ind.nr=k;
 return ind;
}

Indexmatrix& Indexmatrix::delete_rows(const Indexmatrix& ind)
{
 chk_init(*this);
 chk_init(ind);
 if (ind.dim()==0) return *this;
 Indexmatrix sind=sortindex(ind);
 chk_range(ind(sind(0)),ind(sind.dim()-1),nr,nr);
 Integer j,k;
 Integer* mp=m;
 for(j=0;(j<nc);j++){
     mat_xey(ind(sind(0)),mp,m+j*nr);
     mp+=ind(sind(0));
     for(k=1;k<ind.dim();k++){
#if (CONICBUNDLE_DEBUG>=1)
         if (ind(sind(k))==ind(sind(k-1))) MEmessage(MatrixError(ME_unspec,"Indexmatrix::delete_rows(): a row is to be deleted twice",MTindexmatrix));
#endif
         mat_xey(ind(sind(k))-ind(sind(k-1))-1,mp,m+j*nr+ind(sind(k-1))+1);
         mp+=ind(sind(k))-ind(sind(k-1))-1;
     }
     mat_xey(nr-ind(sind(k-1))-1,mp,m+j*nr+ind(sind(k-1))+1);
     mp+=nr-ind(sind(k-1))-1;
 }
 nr-=ind.dim();
 return *this;
}


Indexmatrix& Indexmatrix::delete_cols(const Indexmatrix& ind)
{
 chk_init(*this);
 chk_init(ind);
 if (ind.dim()==0) return *this;
 Indexmatrix sind=sortindex(ind);
 chk_range(ind(sind(0)),ind(sind.dim()-1),nc,nc);
 Integer j,oldj,k;
 oldj=ind(sind(0));
 k=1;
 for(j=ind(sind(0))+1;(j<nc)&&(k<ind.dim());j++){

#if (CONICBUNDLE_DEBUG>=1)
     if (ind(sind(k))<j) MEmessage(MatrixError(ME_unspec,"Indexmatrix::delete_cols(): a column is to be deleted twice",MTindexmatrix));
#endif

     if (j==ind(sind(k))){
         k++;
         continue;
     }

     mat_xey(nr,m+oldj*nr,m+j*nr);
     oldj++;
 }

 mat_xey(nr*(nc-j),m+oldj*nr,m+j*nr);
 nc-=ind.dim();
 return *this;
}

Indexmatrix&  Indexmatrix::insert_row(Integer ind,const Indexmatrix& v)
{
  chk_init(*this);
  chk_init(v);
#if (CONICBUNDLE_DEBUG>=1)
  if ((nc!=v.dim())&&(!((nc==0)&&(nr==0)))) 
    MEmessage(MatrixError(ME_unspec,"Indexmatrix::insert_row(): column dimensions do not match",MTindexmatrix));
  if ((ind<0)||(ind>nr)) MEmessage(MatrixError(ME_unspec,"Indexmatrix::insert_row(): index out of range",MTindexmatrix));
#endif
  nc=v.dim();
  Integer* mp;
  Integer* oldm=m;
  Integer* op;
  int free_oldm=0;
  if (mem_dim<(nr+1)*nc){
    mem_dim=Integer(memarray->get((nr+1)*nc,m));
    if (mem_dim<(nr+1)*nc){
      MEmessage(MEmem((nr+1)*nc,
		      "Indexmatrix::insert_row(): not enough memory",
		      MTindexmatrix));
    }
    free_oldm=1;
  }
  mp=m+(nr+1)*nc-1;
  op=oldm+nr*nc-1;
  mat_xey(nr-ind,mp,-1,op,-1);
  mp-=nr-ind+1;
  op-=nr-ind;
  for(Integer j=1;j<nc;j++){
    mat_xey(nr,mp,-1,op,-1);
    mp-=nr+1;
    op-=nr;
  }
  if (free_oldm){
    mat_xey(ind,m,oldm);
    memarray->free(oldm);
  }
  mat_xey(nc,m+ind,nr+1,v.get_store(),1);
  nr++;
  return *this;
}

Indexmatrix&  Indexmatrix::insert_col(Integer ind,const Indexmatrix& v)
{
  chk_init(*this);
  chk_init(v);
#if (CONICBUNDLE_DEBUG>=1)
  if ((nr!=v.dim())&&(!((nc==0)&&(nr==0)))) 
    MEmessage(MatrixError(ME_unspec,"Indexmatrix::insert_col(): row dimensions do not match",MTindexmatrix));
  if ((ind<0)||(ind>nc)) MEmessage(MatrixError(ME_unspec,"Indexmatrix::insert_col(): index out of range",MTindexmatrix));
#endif
  nr=v.dim();
  Integer* oldm=m;
  int free_oldm=0;
  if (mem_dim<nr*(nc+1)){
    mem_dim=Integer(memarray->get(nr*(nc+1),m));
    if (mem_dim<nr*(nc+1)){
      MEmessage(MEmem(nr*(nc+1),
		      "Indexmatrix::insert_col(): not enough memory",
		      MTindexmatrix));
    }
    free_oldm=1;
    mat_xey(nr*ind,m,oldm);
  }
  mat_xey(nr*(nc-ind),m+nr*(nc+1)-1,-1,oldm+nr*nc-1,-1);
  mat_xey(nr,m+ind*nr,v.get_store());
  if (free_oldm){
    memarray->free(oldm);
  }
  nc++;
  return *this;
}
   


void Indexmatrix::display(std::ostream& out,int precision,int width,
                        int screenwidth) const
{
 out<<"Indexmatrix("<<nr<<","<<nc<<")"<<std::endl;
 if ((nr==0)||(nc==0)) return;
 chk_init(*this);
 if (precision==0) precision=4;
 out.precision(precision);
 if (width==0) width=precision+6;
 if (screenwidth==0) screenwidth=80;
 Integer colnr=screenwidth/(width+1);
 Integer k,i,j;
 Integer maxk=nc/colnr+((nc%colnr)>0);
 Integer maxj;
 for(k=0;k<maxk;k++){
     out<<"columns "<<k*colnr<<" to "<<min(nc,(k+1)*colnr)-1<<std::endl;
     for(i=0;i<nr;i++){
         maxj=min((k+1)*colnr,nc);
         for(j=k*colnr;j<maxj;j++){
             out<<' ';out.width(width);out<<(*this)(i,j);
         }
         out<<std::endl;
     }     
 }
}

Indexmatrix Indexmatrix::operator()(const Indexmatrix &rv,const Indexmatrix &cv) const
{
 chk_init(rv);
 chk_init(cv);
 chk_init(*this);
 if ((rv.dim()==0)||(cv.dim()==0)) return Indexmatrix(0,0,Integer(0));
 chk_range(min(rv),min(cv),nr,nc);
 chk_range(max(rv),max(cv),nr,nc);
 Indexmatrix A(rv.dim(),cv.dim());
 Integer i,j;
 Integer *ap=A.m;
 Integer *rp;
 Integer *cp=cv.m;
 Integer *mcp;           //points to current column indexed by cv
 for(j=A.nc;--j>=0;){
     mcp=m+nr*(*cp++);
     rp=rv.m;
     for(i=A.nr;--i>=0;)
         (*ap++)=mcp[*rp++];
 }
 chk_set_init(A,1);
 return A;
}

Indexmatrix Indexmatrix::operator()(const Indexmatrix &v) const
{
 chk_init(v);
 chk_init(*this);
 if (v.dim()==0) return Indexmatrix(0,0,Integer(0));
 chk_range(min(v),max(v),nr*nc,nr*nc);
 Indexmatrix A(v.nr,v.nc);
 Integer i;
 Integer *ap=A.m;
 Integer *vp=v.m;
 for(i=A.nr*A.nc;--i>=0;){
         (*ap++)=m[*vp++];
 }
 chk_set_init(A,1);
 return A;
}

Indexmatrix Indexmatrix::col(Integer c) const
{
 chk_init(*this);
 chk_range(c,0,nc,1);
 return Indexmatrix(nr,1,m+c*nr);
}

Indexmatrix Indexmatrix::row(Integer r) const
{
 chk_init(*this);
 chk_range(r,0,nr,1);
 return Indexmatrix(1,nc,m+r,nr);
}

Indexmatrix Indexmatrix::cols(const Indexmatrix& v) const
{
 chk_init(v);
 chk_init(*this);
 if (v.dim()==0) return Indexmatrix(0,0,Integer(0));
 chk_range(min(v),max(v),nc,nc);
 Indexmatrix A(nr,v.dim());
 Integer *mp;
 Integer *ap=A.m;
 Integer *vp=v.m;
 Integer i,j;
 for(i=A.nc;--i>=0;){
     mp=m+nr*(*vp++);
     for(j=nr;--j>=0;)
         *ap++=*mp++;
 }
 chk_set_init(A,1);
 return A;
}

Indexmatrix Indexmatrix::rows(const Indexmatrix& v) const
{
 chk_init(v);
 chk_init(*this);
 if (v.dim()==0) return Indexmatrix(0,0,Integer(0));
 chk_range(min(v),max(v),nr,nr);
 Indexmatrix A(v.dim(),nc);
 Integer *vp=v.m;
 for(Integer i=0;i<A.nr;i++){
     mat_xey(nc,A.m+i,A.nr,m+(*vp++),nr);
 }
 chk_set_init(A,1);
 return A;
}

Indexmatrix& Indexmatrix::transpose()
{
 chk_init(*this);
 if ((nr<=1)||(nc<=1)){
     Integer h=nr;nr=nc;nc=h;
     return *this;
 }
 Integer *nm;
 Integer nmem_dim=Integer(memarray->get(nr*nc,nm));
 if (nmem_dim<nr*nc)
     MEmessage(MEmem(nr*nc,"Indexmatrix::transpose() not enough memory",MTindexmatrix));
 Integer i,j;
 Integer *mp=m;
 Integer *nmp;
 for(j=0;j<nc;j++){
     nmp=nm+j;
     for(i=0;i<nr;i++,nmp+=nc)
         *nmp=*mp++;
 }
 i=nr;nr=nc;nc=i;
 memarray->free(m);
 m=nm;
 mem_dim=nmem_dim;
 return *this;
}

Indexmatrix& Indexmatrix::subassign(const Indexmatrix &rv,const Indexmatrix &cv,
                          const Indexmatrix& A)
{
 chk_init(rv);
 chk_init(cv);
 chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
 if ((rv.dim()!=A.nr)||(cv.dim()!=A.nc))
     MEmessage(MEdim(rv.dim(),cv.dim(),A.nr,A.nc,
                 "Indexmatrix::subassign(const Indexmatrix&,const Indexmatrix&,const Indexmatrix&) dimensions do not match",
                 MTindexmatrix));
#endif
 if ((rv.dim()==0)||(cv.dim()==0)) return *this;
 chk_range(min(rv),min(cv),nr,nc);
 chk_range(max(rv),max(cv),nr,nc);
 Integer i,j;
 Integer *ap=A.m;
 Integer *rp;
 Integer *cp=cv.m;
 Integer *mcp;           //points to current column indexed by cv
 for(j=A.nc;--j>=0;){
     mcp=m+nr*(*cp++);
     rp=rv.m;
     for(i=A.nr;--i>=0;)
         mcp[*rp++]=(*ap++);
 }
 return *this;
}

Indexmatrix& Indexmatrix::subassign(const Indexmatrix &v,const Indexmatrix& A)
{
 chk_init(v);
 chk_init(A);
#if (CONICBUNDLE_DEBUG>=1)
 if (v.dim()!=A.dim())
     MEmessage(MEdim(v.dim(),1,A.dim(),1,
                 "Indexmatrix::subassign(const Indexmatrix&,const Indexmatrix&) dimensions do not match",
                 MTindexmatrix));
#endif
 if (v.dim()==0) return *this;
 chk_range(min(v),max(v),nr*nc,nr*nc);
 Integer i;
 Integer *ap=A.m;
 Integer *vp=v.m;
 for(i=A.nr*A.nc;--i>=0;){
         m[*vp++]=(*ap++);
 }
 return *this;
}

Indexmatrix& Indexmatrix::triu(Integer i)
{
 chk_init(*this);
 Integer j;
 for(j=0;j<nc;j++){
     mat_xea(nr-max(Integer(0),j+1-i),m+j*nr+max(Integer(0),j+1-i),Integer(0));
 }
 return *this;
}

Indexmatrix& Indexmatrix::tril(Integer i)
{
 chk_init(*this);
 Integer j;
 for(j=0;j<nc;j++){
     mat_xea(min(nr,j-i),m+j*nr,Integer(0));
 }
 return *this;
}

Indexmatrix operator<(const Indexmatrix& A,const Indexmatrix& B) 
{
 chk_add(A,B);
 Indexmatrix C(A.nr,A.nc);
 Integer *ap=A.m;
 Integer *bp=B.m;
 Integer *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Integer((*ap++)<(*bp++));
 chk_set_init(C,1);
 return C;
}

Indexmatrix operator<=(const Indexmatrix& A,const Indexmatrix& B) 
{
 chk_add(A,B);
 Indexmatrix C(A.nr,A.nc);
 Integer *ap=A.m;
 Integer *bp=B.m;
 Integer *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Integer((*ap++)<=(*bp++));
 chk_set_init(C,1);
 return C;
}

Indexmatrix operator==(const Indexmatrix& A,const Indexmatrix& B) 
{
 chk_add(A,B);
 Indexmatrix C(A.nr,A.nc);
 Integer *ap=A.m;
 Integer *bp=B.m;
 Integer *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Integer((*ap++)==(*bp++));
 chk_set_init(C,1);
 return C;
}

Indexmatrix operator!=(const Indexmatrix& A,const Indexmatrix& B) 
{
 chk_add(A,B);
 Indexmatrix C(A.nr,A.nc);
 Integer *ap=A.m;
 Integer *bp=B.m;
 Integer *cp=C.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*cp++)=Integer((*ap++)!=(*bp++));
 chk_set_init(C,1);
 return C;
}


Indexmatrix operator<(const Indexmatrix& A,Integer d)
{
 chk_init(A);
 Indexmatrix B(A.nr,A.nc);
 Integer *bp=B.m;
 Integer *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Integer((*ap++)<d);
 chk_set_init(B,1);
 return B;
}

Indexmatrix operator>(const Indexmatrix& A,Integer d)
{
 chk_init(A);
 Indexmatrix B(A.nr,A.nc);
 Integer *bp=B.m;
 Integer *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Integer((*ap++)>d);
 chk_set_init(B,1);
 return B;
}

Indexmatrix operator<=(const Indexmatrix& A,Integer d)
{
 chk_init(A);
 Indexmatrix B(A.nr,A.nc);
 Integer *bp=B.m;
 Integer *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Integer((*ap++)<=d);
 chk_set_init(B,1);
 return B;
}

Indexmatrix operator>=(const Indexmatrix& A,Integer d)
{
 chk_init(A);
 Indexmatrix B(A.nr,A.nc);
 Integer *bp=B.m;
 Integer *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Integer((*ap++)>=d);
 chk_set_init(B,1);
 return B;
}

Indexmatrix operator==(const Indexmatrix& A,Integer d)
{
 chk_init(A);
 Indexmatrix B(A.nr,A.nc);
 Integer *bp=B.m;
 Integer *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Integer((*ap++)==d);
 chk_set_init(B,1);
 return B;
}

Indexmatrix operator!=(const Indexmatrix& A,Integer d)
{
 chk_init(A);
 Indexmatrix B(A.nr,A.nc);
 Integer *bp=B.m;
 Integer *ap=A.m;
 Integer i=A.nr*A.nc;
 while(--i>=0)
     (*bp++)=Integer((*ap++)!=d);
 chk_set_init(B,1);
 return B;
}


Indexmatrix& Indexmatrix::concat_right(const Indexmatrix& A)
{
 chk_init(*this);
 chk_init(A);
 if ((A.nr==0)&&(A.nc==0))
   return *this;
 if ((nr==0)&&(nc==0)){
   *this=A;
   return *this;
 } 
#if (CONICBUNDLE_DEBUG>=1)
 if (nr!=A.nr)
     MEmessage(MEdim(nr,nc,A.nr,A.nc,
                 "Indexmatrix::concat_right(const Indexmatrix&) dimensions do not match",MTindexmatrix));
#endif
 if (nr==0){
   nc+=A.nc;
   return *this;
 }
 if (mem_dim<nr*nc+A.nr*A.nc){
     Integer *mnew;
     mem_dim=Integer(memarray->get(nr*nc+A.nr*A.nc,mnew));
     if (mem_dim<nr*nc+A.nr*A.nc){
         MEmessage(MEmem(nr*nc+A.nr*A.nc,
                     "Indexmatrix::concat_right(const Indexmatrix&) not enough memory",
                     MTindexmatrix));
     }
     mat_xey(nr*nc,mnew,m);
     memarray->free(m);
     m=mnew;
 }
 mat_xey(A.nr*A.nc,m+nr*nc,A.m);
 nc+=A.nc;
 chk_set_init(*this,1);
 return (*this);
}

Indexmatrix& Indexmatrix::concat_below(const Indexmatrix& A)
{
 chk_init(*this);
 chk_init(A);
 if ((A.nr==0)&&(A.nc==0))
   return *this;
 if ((nr==0)&&(nc==0)){
   *this=A;
   return *this;
 } 
#if (CONICBUNDLE_DEBUG>=1)
 if (nc!=A.nc)
     MEmessage(MEdim(nr,nc,A.nr,A.nc,
                 "Indexmatrix::concat_below(const Indexmatrix&) dimensions do not match",MTindexmatrix));
#endif
 if (nc==0){
   nr+=A.nr;
   return *this;
 }
 Integer i,j;
 Integer* mp;
 Integer* ap;
 Integer* oldm=m;
 Integer* op;
 int free_oldm=0;
 if (mem_dim<nr*nc+A.nr*A.nc){
     mem_dim=Integer(memarray->get(nr*nc+A.nr*A.nc,m));
     if (mem_dim<nr*nc+A.nr*A.nc){
         MEmessage(MEmem(nr*nc+A.nr*A.nc,
                     "Indexmatrix::concat_below(const Indexmatrix&) not enough memory",
                     MTindexmatrix));
     }
     free_oldm=1;
 }
 mp=m+nr*nc+A.nr*A.nc-1;
 ap=A.m+A.nr*A.nc-1;
 op=oldm+nr*nc-1;
 for(j=nc;--j>=0;){
     for(i=A.nr;--i>=0;)
         (*mp--)=(*ap--);
     if (mp==op) break;   //only possible for last column on same memory
     for(i=nr;--i>=0;)
         (*mp--)=(*op--);
 }
 nr+=A.nr;
 if (free_oldm) memarray->free(oldm);
 chk_set_init(*this,1);
 return (*this);
}

Indexmatrix& Indexmatrix::concat_below(Integer d)
{
 chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
 if ((nc!=1)&&(!((nc==0)&&(nr==0))))
     MEmessage(MEdim(nr,nc,1,1,
                 "Indexmatrix::concat_below(Integer d) matrix is not a column vector",MTmatrix));
#endif
 Integer n=nr*nc;
 if (mem_dim<n+1){
   Integer* mp=m;
   mem_dim=Integer(memarray->get(n+1,m));
   if (mem_dim<n+1){
     MEmessage(MEmem(n+1,
                     "Indexmatrix::concat_below(Integer d) not enough memory",
                     MTmatrix));
   }
   mat_xey(n,m,mp); 
   memarray->free(mp);
 }
 m[n]=d;
 nr=n+1;nc=1;
 return (*this);
}

Indexmatrix& Indexmatrix::concat_right(Integer d)
{
 chk_init(*this);
#if (CONICBUNDLE_DEBUG>=1)
 if ((nr!=1)&&(!((nr==0)&&(nc==0))))
     MEmessage(MEdim(nr,nc,1,1,
                 "Indexmatrix::concat_right(Integer d) matrix is not a row vector",MTmatrix));
#endif
 Integer n=nr*nc;
 if (mem_dim<n+1){
   Integer* mp=m;
   mem_dim=Integer(memarray->get(n+1,m));
   if (mem_dim<n+1){
     MEmessage(MEmem(n+1,
                     "Indexmatrix::concat_right(Integer d) not enough memory",
                     MTmatrix));
   }
   mat_xey(n,m,mp); 
   memarray->free(mp);
 }
 m[n]=d;
 nr=1;
 nc=n+1;
 return (*this);
}


Indexmatrix diag(const Indexmatrix& A)
{
 chk_init(A);
 if (min(A.nr,A.nc)==1){ //make a diagonal matrix of this vector
     Integer k=max(A.nr,A.nc);
     Indexmatrix B(k,k,Integer(0));
     for(Integer i=0;i<k;i++)
         B(i,i)=A(i);
     return B;
 }
 return Indexmatrix(min(A.nr,A.nc),1,A.m,A.nr+1);
}

Indexmatrix sumrows(const Indexmatrix& A)
{
 chk_init(A);
 Indexmatrix v(1,A.nc,Integer(0));
 for(Integer i=0;i<A.nr;i++)
     mat_xpey(A.nc,v.m,1,A.m+i,A.nr);
 chk_set_init(v,1);
 return v;
}

Indexmatrix sumcols(const Indexmatrix& A)
{
 chk_init(A);
 Indexmatrix v(A.nr,1,Integer(0));
 for(Integer i=0;i<A.nc;i++)
     mat_xpey(A.nr,v.m,A.m+i*A.nr);
 chk_set_init(v,1);
 return v;
}

Integer sum(const Indexmatrix& A)
{
 chk_init(A);
 Integer s=0;
 Integer i;
 Integer *mp=A.m;
 for(i=A.nr*A.nc;--i>=0;){
     s+=(*mp++);
 }
 return s;
}

Indexmatrix maxrows(const Indexmatrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Indexmatrix(0,0,Integer(0));
 Indexmatrix v(1,A.nc);
 Integer maxd;
 Integer i,j;
 for(j=0;j<A.nc;j++){
     maxd=A(0,j);
     for(i=1;i<A.nr;i++)
         maxd=max(maxd,A(i,j));
     v(j)=maxd;
 }
 chk_set_init(v,1);
 return v;
}

Indexmatrix maxcols(const Indexmatrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Indexmatrix(0,0,Integer(0));
 Indexmatrix v(A.nr,1);
 Integer maxd;
 Integer i,j;
 for(i=0;i<A.nr;i++){
     maxd=A(i,0);
     for(j=1;j<A.nc;j++)
         maxd=max(maxd,A(i,j));
     v(i)=maxd;
 }
 chk_set_init(v,1);
 return v;
}

Integer max(const Indexmatrix& A,Integer *iindex,Integer *jindex)
{
 chk_init(A);
 if (A.dim()==0) return min_Integer;
 Integer maxd;
 if (iindex==0){
     Integer i=A.nr*A.nc-1;
     const Integer *ap=A.m;
     maxd=*ap++;
     for(;--i>=0;){
         maxd=max(maxd,*ap++);
     }
 }
 else{
     const Integer *ap=A.m;
     maxd=*ap++;
     Integer besti=0;
     Integer n=A.nr*A.nc;
     for(Integer i=0;++i<n;ap++){
         if (*ap>maxd){
             maxd=*ap;
             besti=i;
         }
     }
     if (jindex!=0){
         *jindex=besti/A.nr;
         *iindex=besti%A.nr;
     }
     else *iindex=besti;
 }
     
 return maxd;
}

Indexmatrix minrows(const Indexmatrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Indexmatrix(0,0,Integer(0));
 Indexmatrix v(1,A.nc);
 Integer mind;
 Integer i,j;
 for(j=0;j<A.nc;j++){
     mind=A(0,j);
     for(i=1;i<A.nr;i++)
         mind=min(mind,A(i,j));
     v(j)=mind;
 }
 chk_set_init(v,1);
 return v;
}

Indexmatrix mincols(const Indexmatrix& A)
{
 chk_init(A);
 if (A.dim()==0) return Indexmatrix(0,0,Integer(0));
 Indexmatrix v(A.nr,1);
 Integer mind;
 Integer i,j;
 for(i=0;i<A.nr;i++){
     mind=A(i,0);
     for(j=1;j<A.nc;j++)
         mind=min(mind,A(i,j));
     v(i)=mind;
 }
#if (CONICBUNDLE_DEBUG>=1)
 chk_set_init(v,1);
#endif 
 return v;
}

Integer min(const Indexmatrix& A,Integer *iindex,Integer *jindex)
{
 chk_init(A);
 if (A.dim()==0) return max_Integer;
 Integer mind;
 if (iindex==0){
     Integer i=A.nr*A.nc-1;
     const Integer *ap=A.m;
     mind=*ap++;
     for(;--i>=0;){
         mind=min(mind,*ap++);
     }
 }
 else{
     Integer n=A.nr*A.nc;
     const Integer *ap=A.m;
     mind=*ap++;
     Integer besti=0;
     for(Integer i=0;++i<n;ap++){
         if (*ap<mind){
             mind=*ap;
             besti=i;
         }
     }
     if (jindex!=0){
         *jindex=besti/A.nr;
         *iindex=besti%A.nr;
     }
     else *iindex=besti;
 }
     
 return mind;
}

void sortindex(const Indexmatrix& vec,Indexmatrix& ind)
{
 ind.init(Range(0,vec.dim()-1));
 Integer* ip=ind.get_store();
 const Integer* vp=vec.get_store();
 std::sort( ip, ip + vec.dim(), mat_less_index<Integer>(vp) );
}

Integer trace(const Indexmatrix& A)
{
 chk_init(A);
 Integer sum=0;
 Integer k=min(A.nr,A.nc);
 for(Integer i=0;i<k;i++) sum+=A(i,i);
 return sum;
}

Indexmatrix& Indexmatrix::rand(Integer innr,Integer innc,Integer lb,Integer ub,GB_rand* rgp)
{
  if (rgp==0)
    rgp=&mat_randgen;
  newsize(innr,innc);
  Integer r=ub-lb+1;
  for(Integer i=0;i<nr*nc;i++)
    m[i]=Integer(rgp->unif_long(r))+lb;
  chk_set_init(*this,1);
  return *this;
}

Indexmatrix& Indexmatrix::shuffle(CH_Tools::GB_rand* rgp)
{
  if (rgp==0)
    rgp=&mat_randgen;
  Integer dim=nr*nc;
  for(Integer i=0;i<dim;i++){
    Integer h=Integer(rgp->unif_long(dim-i));
    Integer d=m[i+h];
    m[i+h]=m[i];
    m[i]=d;
  }
  return *this;
} 

Indexmatrix& Indexmatrix::sign(void)
{
 chk_init(*this);
 for(Integer i=0;i<nr*nc;i++)
   m[i]=CH_Matrix_Classes::sign(m[i]); 
 return *this;
}

Indexmatrix& Indexmatrix::abs(void)
{
 chk_init(*this);
 for(Integer i=0;i<nr*nc;i++) 
   m[i]=CH_Matrix_Classes::abs(m[i]);
 return *this;
}

Indexmatrix abs(const Indexmatrix& A)
{
 chk_init(A);
 Indexmatrix B(A.nr,A.nc);
 for(Integer i=0;i<A.nr*A.nc;i++) B(i)=abs(A(i));
 chk_set_init(B,1);
 return B;
}

Indexmatrix transpose(const Indexmatrix& A)
{
 chk_init(A);
 if ((A.nr<=1)||(A.nc<=1)) return Indexmatrix(A.nc,A.nr,A.m);
 Indexmatrix B(A.nc,A.nr);
 Integer i;
 for(i=0;i<B.nc;i++){
     mat_xey(B.nr,B.m+i*B.nr,1,A.m+i,A.nr);
 }
 chk_set_init(B,1);
 return B;
}

std::ostream& operator<<(std::ostream& o,const Indexmatrix &A)
{
 chk_init(A);
 o<<A.nr<<" "<<A.nc<<'\n';
 Integer i,j;
 for(i=0;i<A.nr;i++){
     for(j=0;j<A.nc;j++) o<<' '<<A(i,j);
     o<<'\n';
 }
 return o;
}

std::istream& operator>>(std::istream& in,Indexmatrix &A)
{
 Integer d;
 Integer nr,nc;
 in>>d;
 nr=Integer(d+.5);
 in>>d;
 nc=Integer(d+.5);
 A.newsize(nr,nc);
 Integer i,j;
 for(i=0;i<nr;i++)
     for(j=0;j<nc;j++)
         in>>A(i,j);
 chk_set_init(A,1);
 return in;
}


}

