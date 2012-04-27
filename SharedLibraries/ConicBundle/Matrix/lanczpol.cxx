/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/lanczpol.cxx

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



#include <cmath>
#include <limits>
#include "mymath.hxx"
#include "lanczpol.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {

// *************************************************************************
//                              constructor(s)
// *************************************************************************

Lanczpol::Lanczpol()
{
 ierr=0;
 eps=1.e-5;
 minval=0.;
 maxop=-1;       //flag for choose default
 maxiter=-1;     //flag for infinite
 guessmult=10;
 maxguessiter=1;
 print_level=0;
 myout=0;
 choicencheb=-1;
 nchebit=20;
 polval=1000;
 iter=0;
 nmult=0;
 blocksz=0;
 neigfound=0;
 choicenbmult=-1;
 nblockmult=20;
 retlanvecs=-1;
 stop_above=0;
 upper_bound=0.;

 nlanczvecs=0;
 mcheps=eps_Real;
 //while(1.+mcheps>1.) mcheps/=2.;
 //mcheps*=2.;

 ncalls=0;

 randgen.init(1);
 
 time_mult=Microseconds(long(0));    
 time_mult_sum=Microseconds(long(0));
 time_iter=Microseconds(long(0));
 time_sum=Microseconds(long(0));
}

Lanczpol::~Lanczpol()
{
  bigmatrix=0;
}

// *************************************************************************
//                            get_lanczosvecs
// *************************************************************************

int Lanczpol::get_lanczosvecs(Matrix& val,Matrix& vecs) const
{
 if (nmult==0) return 1;
 Integer n,m;
 X.dim(n,m);
 m=min(m,nlanczvecs);
 if (retlanvecs>0) m=min(m,retlanvecs);
 val.init(m,1,d.get_store());
 vecs.init(n,m,X.get_store());
 return 0;
}
 

// *************************************************************************
//                              compute
// *************************************************************************

int Lanczpol::compute(const Lanczosmatrix* bigmat,Matrix& eigval,Matrix& eigvec,
                     Integer nreig,Integer in_blocksz,Integer maxj)
{
  ncalls++;
 
 //--- check and set input parameters
 bigmatrix=bigmat;
 Integer n;
 n=bigmatrix->lanczosdim();
 
 if (nreig<=0) return 0;
 else nreig=min(nreig,n);
 
 if (in_blocksz<=0) blocksz=1;
 else blocksz=min(in_blocksz,max(Integer(1),n/10));

 neigfound=0;
 time_mult_sum=Microseconds(long(0));
 time_sum=Microseconds(long(0));
 myclock.start();

 //--- initialize size parameters

 Integer maxnblockmult=min(100,n/blocksz-1);
 if ((choicenbmult<0)&&(choicencheb<0)) {
   if (n<300) maxnblockmult=30; //expect to need this number of Lanczos steps per restart
   else if (n<1000) maxnblockmult=50;
   else maxnblockmult=100;
   maxnblockmult=min(maxnblockmult,n/blocksz-1);
   Integer multflops=bigmat->lanczosflops();
   if (multflops<maxnblockmult*blocksz*n){
     //start with small Chebychev
     nchebit=10;
   }
   else {
     nchebit=0;
   }
   nblockmult=min(10,maxnblockmult);
 }
 else { 
   if (choicenbmult<=0) 
     nblockmult=min(10,maxnblockmult);
   else 
     nblockmult=min(choicenbmult,maxnblockmult);
   if (choicencheb<0) {
     Integer multflops=bigmat->lanczosflops();
     if (multflops<nblockmult*blocksz*n){
       nchebit=10;
     }
     else {
       nchebit=0;
     }
   }
   else 
     nchebit=choicencheb;
 }
 mymaxj=min(n,max(nreig+maxnblockmult*blocksz,nreig+2*guessmult));
 nblockmult=min(nblockmult,(mymaxj-nreig)/blocksz-1);


 maxj=mymaxj;
 Integer maxmult=maxop;
 
 X.newsize(n,maxj); chk_set_init(X,1);
 Xqr.newsize(n,maxj); chk_set_init(Xqr,1);
 C.newsize(maxj,maxj); chk_set_init(C,1);
 d.newsize(maxj,1); chk_set_init(d,1);
 e.newsize(maxj,1); chk_set_init(e,1);
 u.newsize(n,max(blocksz,Integer(2))); chk_set_init(u,1);
 v.newsize(n,max(blocksz,Integer(2))); chk_set_init(v,1);
 w.newsize(n,max(blocksz,Integer(2))); chk_set_init(w,1);
 /*
 X.init(n,maxj,std::numeric_limits<double>::signaling_NaN()); chk_set_init(X,1);
 Xqr.init(n,maxj,std::numeric_limits<double>::signaling_NaN()); chk_set_init(Xqr,1);
 C.init(maxj,maxj,std::numeric_limits<double>::signaling_NaN()); chk_set_init(C,1);
 d.init(maxj,1,std::numeric_limits<double>::signaling_NaN()); chk_set_init(d,1);
 e.init(maxj,1,std::numeric_limits<double>::signaling_NaN()); chk_set_init(e,1);
 u.init(n,max(blocksz,Integer(2)),std::numeric_limits<double>::signaling_NaN()); chk_set_init(u,1);
 v.init(n,max(blocksz,Integer(2)),std::numeric_limits<double>::signaling_NaN()); chk_set_init(v,1);
 w.init(n,max(blocksz,Integer(2)),std::numeric_limits<double>::signaling_NaN()); chk_set_init(w,1);
 */
 

 nmult=0;
 iter=0;
 errc=0.;
 maxval=-1e40;
 minval=0.;

 //minvec.init(0,0,0.);
 //randgen.init(1);  
 //minvec and mat_randgen may influence repeated runs with the same 
 //calling parameters to remove this influence, uncomment the two
 //lines above. 
          

 //--- generate starting vectors from proposed vectors
 Integer m,i,nproposed;
 eigvec.dim(m,nproposed);
 if (m==n){
     nproposed-=neigfound;
     nproposed=min(nproposed,maxj-1-neigfound);
     //- copy found and proposed vectors that fit into storage
     mat_xey(n*(neigfound+nproposed),X.get_store(),eigvec.get_store());

     //- if there are still too many, take the eigenvectors
     //- to the largest eigenvalues of the projected matrix
     if (nproposed>blocksz){
         orthog(neigfound,nproposed,X,C);
	 time_mult=Microseconds(long(0)); 
         sectn(X,neigfound,nproposed,C,d,u,v,minval);
         time_mult_sum+=time_mult;
         nmult+=nproposed;
         if (nproposed>=3){
             maxval=sum(d(Range(neigfound,neigfound+nproposed-1)))/nproposed;
         }
         else {
             maxval=d(neigfound+1);
         }
         rotate(neigfound,nproposed,blocksz,C,X,v);
     }
     nproposed=blocksz;
     //- perturb the vectors a little bit
     Real randfac=0.05/::sqrt(double(n));
     for(i=n*neigfound;i<n*(neigfound+blocksz);i++)
         X(i)+=randfac*(2.*randgen.next()-1.);
 }
 else {
     neigfound=0;
     nproposed=0;
 }
     
 //--- compute guess for interval of Chebychev polynomial
 if (nchebit>0){
     if ((minvec.rowdim()==n)&&(minvec.coldim()==1)){
         nproposed=min(nproposed,maxj-2-neigfound);
         for(i=0;i<n;i++)
             X(n*(neigfound+nproposed)+i)=minvec(i);
         nproposed+=1;
     }
     ierr=guess_extremes(nproposed);
     minval= d(neigfound+1);
     minvec= X.col(neigfound+1);
     //fill remaining vectors again with starting vectors
     if (m==n){
         eigvec.dim(m,nproposed);
         nproposed-=neigfound;
         Real randfac=0.05/::sqrt(double(n));
         nproposed=min(nproposed,blocksz-1);
         for(i=n*(neigfound+1);i<n*(neigfound+1+nproposed);i++)
             X(i)=eigvec(i)+randfac*(2.*randgen.next()-1.);
         nproposed+=1;
     }
 }
 
 //--- call lanczos
 if ((myout)&&(print_level>0)) {
     myout->precision(8);
     (*myout)<<" maxval="<<maxval<<" minval="<<minval<<" ncheb="<<nchebit;
     (*myout)<<" nblock="<<nblockmult<<" relprec="<<eps<<std::endl;
 }
 ierr|=dhua(nreig,nproposed,maxmult);
 if (ierr>0){
   time_sum=myclock.time();     
   return 1;
 }

 //--- store the parameters now known
 eigvec.init(n,neigfound,X.get_store());
 eigval.init(neigfound,1,d.get_store());
 time_sum=myclock.time();  
 return (ierr!=0);
}

// *****************************************************************************
//                                guess_extremes
// *****************************************************************************

//iterates a few time in order to obtain a good guess for minimal and
//maximal eigenvalue. These will later be used to determine the
//interval of the Chebychef polynomial.
//The routine performs a "guessmult" steps Lanczosapproximation and
//returns the largest value in d(neigfound) and the smallest in d(neigfound+1)

int Lanczpol::guess_extremes(Integer nproposed)
{

 Integer n,q;
 X.dim(n,q);
 Integer lblocksz=2;
 Integer usemult=min(guessmult,(mymaxj-neigfound)/lblocksz-1);
 if (usemult<=0){
   if (myout){
     (*myout)<<"**** ERROR in Lanczpol::guess_extremes(...): no space left to compute the guesses"<<std::endl;
     (*myout)<<"columns available="<<mymaxj<<" neigfound="<<neigfound<<" lblocksz="<<lblocksz<<std::endl;
   }
   return 1;
 }
 time_mult=Microseconds(long(0)); 

 //fill empty starting vectors with random vectors
 for(Integer i=neigfound+nproposed;i<neigfound+lblocksz;i++){ 
     random(X,i);
 } 
  
 orthog(neigfound,max(lblocksz,nproposed),X,C);
 sectn(X,neigfound,max(lblocksz,nproposed),C,d,u,v,minval);
 nmult+=max(nproposed,lblocksz);
 rotate_extremes(neigfound,max(lblocksz,nproposed),d,C,X,v);

 time_mult_sum+=time_mult;
 
 do {
     if ((myout)&&(print_level>0)){
         (*myout)<<"G: ";myout->width(2);(*myout)<<iter<<":";
         (*myout)<<"  blocksz=";myout->width(2);(*myout)<<lblocksz;
	 (*myout)<<"  s=";myout->width(2);(*myout)<<guessmult;
         (*myout)<<"  neigfound=";myout->width(2);(*myout)<<neigfound;
	 (*myout)<<"  nmult=";myout->width(4);(*myout)<<nmult;
         (*myout)<<"  maxval="<<maxval<<std::endl;
         (*myout)<<"  minval="<<minval<<std::endl;
     }

     iter++;
     Microseconds secs=myclock.time();
     time_mult=Microseconds(long(0)); 
   
     Integer sbs=guessmult*lblocksz;
     bklanc(neigfound,lblocksz,guessmult,d,C,X,e,u,v);
     nmult+=sbs;
     eigen(neigfound,lblocksz,sbs,C,d,u,v,minval);
     maxval=max(maxval,(d(1)+d(2))/2.);
     if ((myout)&&(print_level>0)){
         (*myout)<<"GE: lblocksz="<<lblocksz;
         (*myout)<<"  nmult=";myout->width(4);(*myout)<<nmult<<"\n  ";
         (*myout)<<"  maxval="<<maxval<<std::endl;
         (*myout)<<neigfound<<":";
         for(Integer i=0;i<neigfound;i++) (*myout)<<" "<<d(i);
         (*myout)<<"\n  "<<sbs<<":";
         {for(Integer i=0;i<sbs;i++) (*myout)<<" "<<d(neigfound+i);}
         (*myout)<<std::endl;
     }   
     Integer nconv;
     cnvtst(neigfound,lblocksz,errc,eps,d,e,nconv);
     rotate_extremes(neigfound,sbs,d,C,X,v);
     
     time_mult_sum+=time_mult;
     time_iter=myclock.time()-secs;
 }while(iter<maxguessiter);

 if ((myout)&&(print_level>0)){
         (*myout)<<"G__";myout->width(2);(*myout)<<iter<<":";
         (*myout)<<"  blocksz=";myout->width(2);(*myout)<<lblocksz;
	 (*myout)<<"  s=";myout->width(2);(*myout)<<guessmult;
         (*myout)<<"  neigfound=";myout->width(2);(*myout)<<neigfound;
	 (*myout)<<"  nmult=";myout->width(4);(*myout)<<nmult;
     (*myout)<<"  guessvals:";
     for(Integer i=0;i<neigfound+lblocksz;i++) (*myout)<<" "<<d(i);
     (*myout)<<"  minval="<<minval<<std::endl;
 }

 return 0;
}

// *****************************************************************************
//                                dhua
// *****************************************************************************

// X       Matrix(n,q), n dimension der sym Matrix.
// pinit      initial block size
// nreig      number of eigenvalues to be computed
// maxmult    maximal number of matrix vector computations 
// neigfound  number of eigenvalues found (initially zero unless some known)
// nmult      counts number of matrix vector computations
// iter       number of program iterations
// nchebit    number of chebychev iterations
// randseed   
// eps        relative precision for eigenvalues
// errc       accumulated error, set to 0. initally
// minval         block Cheb parameter (best least eigenval of A), 0. ok

/*
c on input: x (of size (n,q) contains eigenvector estimates
c on output: vec ... eigenvector
c            lambda ... eigenvalue
c	     ie (integer): 0 if everything ok, else 1
c
C       This program implements the Block Chebychev-Lanczos method
C     for computing the largest R eigenvalues and eigenvectors of
C     an NxN symmetric matrix A. If the least R eigenvalues and
C     corresponding eigenvectors of A are found, then the Block
C     Chebychev-Lanczos method is applied to symmetric matrix -A. 
C
      INTEGER n,Q,PINIT,R,S,PS,P,ZZ, nmax, ie
	parameter( nmax = 2001, q = 20)
	double precision x(nmax,q), c(q,q),u(nmax),v(nmax),w(nmax)
	double precision lambda, vec(n), vecin( n)
	double precision EPS,ERRC,AF,D(q),E(q)     
      EXTERNAL OP 
C 
C     Description of parameters:
C
C       N: Integer variable. The order of the symmetric matrix A.
C
C       Q: Integer variable. The number of column vectors of length N 
C          contained in array X.     
C
C       PINIT: Integer variable. The initial block size to be used in
C              the Block Chebychev-Lanczos method.
C
C       R: Integer variable. The number of eigenvalues and eigenvectors
C          being computed.
C
C       MMAX: Integer variable. The maximum number of matrix-vector 
C             products A*X where X is a vector. MMAX should be given
C             a very large value.
C
C       M: Integer variable. M gives the number of eigenvalues and
C          eigenvectors already computed. Thus, initially, M should be
C          zero. If M is greater than zero, then columns one through M
C          of the array X are assumed to contain the computed 
C          approximations to the eigenvectors corresponding to the
C          largest M eigenvalues. At exit, M contains a value equal to 
C          the total number of eigenvalues computed including any
C          already computed when the program was entered. Thus, at exit,
C          the first M elements of D and the first m columns of X will
C          contain approximations to the largest M eigenvalues and 
C          corresponding eigenvectors of A, respectively.
C
C       IE: Integer variable. The value of IE indicates whether the
C           program terminated successfully. If IE=1, the value of
C           MMAX was exceeded before R eigenvalues were computed.
C
C       IMM: Integer variable. IMM counts the number of matrix-vector
C            products computed which is the number of times the 
C            Subroutine named by OP is called. Thus, initially, IMM
C            must be set to zero.
C
C       ITER: Integer variable. ITER counts the number of the program
C             iterations. Thus, initially, ITER must be set to zero.
C
C       IM: Integer variable. IM is the number of internal cycle in
C           one Block Chebychev iteration. Commonly, IM is taken as
C           the value between 20 and 50.
C
C       ZZ: Integer variable. ZZ is used to produce random number.
C           ZZ may be taken as any value. Commonly set ZZ=0 or 1.
C
C       EPS: Double precision variable. Initially, EPS should contain
C            a value indicating the relative precision to which the
C            program will attempt to compute the eigenvalues and
C            eigenvectors of matrix A. For eigenvalues less in modulus
C            than 1, EPS will be an absolute tolerance.
C
C       ERRC: Double precision variable. ERRC measures the accumuted
C             error in the eigenvalues and eigenvectors. Initially,
C             ERRC must be set to zero.
C
C       AF: Double precision variable. AF is a parameter in the Block
C           Chebychev iteration. It is the best to take AF as the
C           least eigenvalue of matrix A. Initially, AF may be taken
C           as zero.
C
C       X: Double precision variable. X contains the computed eigenvectors.
C          X should be an array containing at least N*Q elements. X is
C          used not only to store the eigenvectors computed by the
C          program, but also as working storage for the Block Chebychev-
C          Lanczos method. At exit, the first M columns of X contain
C          the corresponding eigenvector approximations.
C 
C      OP: Subroutine name. The actual argument corresponding to OP
C          should be the name of a subroutine used to define the matrix A.
C          This subroutine should have three arguments N, U and V, where
C          N is the order of the matrix A, U and V are two one dimensional
C          arrays of length N. Then the statement 
C                          CALL OP(N,U,V)
C          should result in vector A*U being computed and stored in V.
C
*/

int Lanczpol::dhua(Integer nreig,Integer nproposed, Integer maxmult)
{
 Integer n,q;
 X.dim(n,q);
 time_mult=Microseconds(long(0)); 

 nproposed=min(blocksz,nproposed);
 for(Integer i=neigfound+nproposed;i<neigfound+blocksz;i++){ 
     random(X,i);
 } 
 
 orthog(neigfound,blocksz,X,C);
 
 if(nchebit==0){
   sectn(X,neigfound,blocksz,C,d,u,v,minval);
   nmult+=nproposed;
   if (nproposed>=3){
     maxval=sum(d(Range(neigfound,neigfound+nproposed-1)))/nproposed;
   }
   else {
     maxval=d(neigfound+1);
   }
   rotate(neigfound,blocksz,blocksz,C,X,v);
 }

 time_mult_sum+=Microseconds(long(0))-time_mult;
 
 //for detecting when the maximum numerical precision is reached
 Integer stalldim=5;
 Matrix eigval_stall(stalldim,1,min_Real);
 Matrix eigvec_stall(stalldim,1,max_Real);
 Integer stall_pos=0;
 bool stalled=false;

 Integer nconv;
 Integer sbs;
 do {

     iter++;
     Microseconds secs=myclock.time();
     time_mult=Microseconds(long(0));

     if ((myout)&&(print_level>0)){
         (*myout)<<"L: ";myout->width(2);(*myout)<<iter<<":";
         (*myout)<<"  blocksz=";myout->width(2);(*myout)<<blocksz;
	 (*myout)<<"  s=";myout->width(2);(*myout)<<nblockmult;
         (*myout)<<"  neigfound=";myout->width(2);(*myout)<<neigfound;
	 (*myout)<<"  nmult=";myout->width(4);(*myout)<<nmult;
         if (nchebit) (*myout)<<"  nchebit="<<nchebit<<" maxval="<<maxval;
         if (stop_above) (*myout)<<" stop>"<<upper_bound;
         (*myout)<<std::endl;
     }

     sbs=nblockmult*blocksz;
     if (nchebit>0){
       if (iter<=20){
	 if (bklanccheb(neigfound,blocksz,nblockmult,C,X,e,v)){
	   if (myout){
	     (*myout)<<"**** WARNING: Lanczpol::dhua(): bklanccheb failed in iteration "<<iter<<std::endl;
	     (*myout)<<" nchebit="<<nchebit<<" neigfound="<<neigfound<<" blocksz="<<blocksz<<" nblockmult="<<nblockmult<<std::endl;
	   }
	 }
       }
       else {
	 if (bkqrlanccheb(neigfound,blocksz,nblockmult,C,X,e,v)){
	   if (myout){
	     (*myout)<<"**** WARNING: Lanczpol::dhua(): bkqrlanccheb failed in iteration "<<iter<<std::endl;
	     (*myout)<<" nchebit="<<nchebit<<" neigfound="<<neigfound<<" blocksz="<<blocksz<<" nblockmult="<<nblockmult<<std::endl;
	   }
	 }	
       }   
     }
     else {
       if (iter<=20) {
	 if (bklanc(neigfound,blocksz,nblockmult,d,C,X,e,u,v)){
	   if (myout){
	     (*myout)<<"**** WARNING: Lanczpol::dhua(): bklanc failed in iteration "<<iter<<std::endl;
	     (*myout)<<" nchebit="<<nchebit<<" neigfound="<<neigfound<<" blocksz="<<blocksz<<" nblockmult="<<nblockmult<<std::endl;
	   }
	 }	   
       }
       else {
	 if (bkqrlanc(neigfound,blocksz,nblockmult,d,C,X,e,u,v)){
	   if (myout){
	     (*myout)<<"**** WARNING: Lanczpol::dhua(): bkqrlanccheb failed in iteration "<<iter<<std::endl;
	     (*myout)<<" nchebit="<<nchebit<<" neigfound="<<neigfound<<" blocksz="<<blocksz<<" nblockmult="<<nblockmult<<std::endl;
	   }
	 }
       }
     }
     nmult+=sbs*max(Integer(1),nchebit);
     Real dummy=-10.;
     eigen(neigfound,blocksz,sbs,C,d,u,v,dummy);
     if ((myout)&&(print_level>0)){
         (*myout)<<"LE: blocksz="<<blocksz<<"\n  "<<neigfound<<":";
         for(Integer i=0;i<neigfound;i++) (*myout)<<" "<<d(i);
         (*myout)<<"\n  "<<sbs<<":";
         {for(Integer i=0;i<sbs;i++) (*myout)<<" "<<d(neigfound+i);}
         (*myout)<<std::endl;
     }
     
     if (nchebit>0){
         bool maxval_too_close=(d(neigfound)<10.);
         rotate(neigfound,sbs,sbs,C,X,v);
         sectn(X,neigfound,sbs,C,d,u,v,minval);
	 if (maxval_too_close){
	   Integer i=1;
	   while((i<sbs)
		 &&(abs(d(neigfound)-d(neigfound+i))<1e-3*max(abs(d(neigfound)),1.))
		 ){
	     i++;
	   }
	   if (i<sbs){
	     maxval=d(neigfound+i);
	   }
	   else {
	     maxval=maxval-(maxval-minval)*1e-3;
	   }
	 }
         //maxval=min(d(neigfound+min(Integer(4),sbs)),
	 //	    minval+(d(neigfound)-minval)/(1.+1./Real(nchebit)));
         nmult+=sbs;
	 if ((myout)&&(print_level>0)){
	   (*myout)<<"  {"<<sbs<<":";
	   {for(Integer i=0;i<sbs;i++) (*myout)<<" "<<d(neigfound+i);}
	   (*myout)<<"}\n  [";
	   {for(Integer i=0;i<sbs;i++) (*myout)<<" "<<e(neigfound+i);}
	   (*myout)<<"]"<<std::endl;
	 }
	 cnvtst(neigfound,sbs,errc,eps,d,e,nconv);
     }
     else {
         cnvtst(neigfound,blocksz,errc,eps,d,e,nconv);
     }
     nlanczvecs=neigfound+sbs;
     neigfound+=nconv;

     //collect stall information
     if (nconv>0){
       //reinitialize stall
       eigval_stall.init(stalldim,1,min_Real);
       eigvec_stall.init(stalldim,1,max_Real);
       stall_pos=0;
       stalled=false;
     }
     else{
       //if there was no improvement in both eigenvalue and eigenvector for the last stalldim iterations
       Real t=fabs(d(neigfound));
       Real mev=max(eigvec_stall);
       if (
	   (d(neigfound)-min(eigval_stall)<t*10.*mcheps)
	   && (
	       (e(neigfound)>mev)
	       ||
	       (100*(log(e(neigfound))-log(mev))>log(t*(eps+10.*Real(n)*mcheps)+errc)-log(mev))
	       )
	   ){
	 stalled=true;
       }
     }
     eigval_stall(stall_pos)=d(neigfound);
     eigvec_stall(stall_pos)=e(neigfound);
     stall_pos=(stall_pos+1)%stalldim;

     //check termination
     if (neigfound>=nreig) {
       time_iter=myclock.time()-secs;
       time_mult_sum+=time_mult;
       break;
     } 
              
     if ((nchebit>0)&&(nconv>0)) maxval=d(min(nreig+1,nlanczvecs-1));

     if ((stop_above)&&
         (nlanczvecs>=nreig)&&
         (upper_bound<d(nreig-1))){
         time_iter=myclock.time()-secs;
	 time_mult_sum+=time_mult;
         break;
     }

     if ((maxmult>=0)&&(nmult>maxmult)){
         if ((myout)&&(print_level>0)){
             (*myout)<<"\nLanczpol: nmult>maxmult: "<<nmult<<">"<<maxmult<<std::endl;
             (*myout)<<"         neigfound="<<neigfound<<"  nreig="<<nreig<<std::endl;
             (*myout)<<"         iter="<<iter<<"  s="<<nblockmult<<"  blocksz="<<blocksz<<std::endl;
             //d.display(*myout);
         }
         rotate(neigfound-nconv,sbs,sbs,C,X,v);
	 time_iter=myclock.time()-secs;
	 time_mult_sum+=time_mult;
         return 1;
     }
     if ((maxiter>=0)&&(iter>=maxiter)){
         if ((myout)&&(print_level>0)){
             (*myout)<<"\nLanczpol: iter>maxiter: "<<iter<<">"<<maxiter<<std::endl;
             (*myout)<<"         neigfound="<<neigfound<<"  nreig="<<nreig<<std::endl;
             (*myout)<<"         iter="<<iter<<"  s="<<nblockmult<<"  blocksz="<<blocksz<<std::endl;
             //d.display(*myout);
         }
         rotate(neigfound-nconv,sbs,sbs,C,X,v);
	 time_iter=myclock.time()-secs;
	 time_mult_sum+=time_mult;
         return 2;
     }
     if (stalled){
         if ((myout)&&(print_level>0)){
             (*myout)<<"\nLanczpol: stalled: "<<iter<<std::endl;
             (*myout)<<"         neigfound="<<neigfound<<"  nreig="<<nreig<<std::endl;
             (*myout)<<"         iter="<<iter<<"  s="<<nblockmult<<"  blocksz="<<blocksz<<std::endl;
             //d.display(*myout);
         }
         rotate(neigfound-nconv,sbs,sbs,C,X,v);
	 time_iter=myclock.time()-secs;
	 time_mult_sum+=time_mult;
         return 3;
     }

     //check whether it is necessary to reduce the block size
     if (neigfound+2*blocksz>mymaxj){
       blocksz=max(Integer(1),(mymaxj-neigfound)/2);
     }
     rotate(neigfound-nconv,sbs,nconv+blocksz,C,X,v);

    //--- determine new nblockmult and ncheb if default
    time_iter=myclock.time()-secs;
    time_mult_sum+=time_mult;
    if ((choicenbmult<0)&&(choicencheb<0)){
      if (nconv>0) {
	if (nchebit>0){
	  nchebit=10;
	}
	nblockmult=min(Integer(20),(mymaxj-neigfound)/blocksz-1);
	continue;
      }
      if (nchebit>0){
	nchebit=min(nchebit+10,Integer(100));
	nblockmult=min(Integer(20),(mymaxj-neigfound)/blocksz-1);
	continue;
      }
      if (nblockmult==min(Integer(50),(mymaxj-neigfound)/blocksz-1)){
	nchebit=10;
	if (sbs>=3)
	  maxval=d(neigfound+2);
	else
	  maxval=d(neigfound+1);
	minval=min(minval,d(neigfound+sbs-1));
	nblockmult=min(Integer(20),(mymaxj-neigfound)/blocksz-1);
	continue;
      }
      nblockmult=min(nblockmult+10,(mymaxj-neigfound)/blocksz-1);
      continue;
    }
    

    /*
    if ((choicenbmult<0)&&(choicencheb<0)){
      if (nchebit>0){
	if ((iter==2)||(nconv>0)){ 
	  nchebit=30;
	  nblockmult=min(Integer(25),(mymaxj-neigfound)/blocksz-1);
	  continue;
	}
	if (iter%2==0){ 
	  nchebit=min(nchebit+10,Integer(100));
	}  
	continue;
      }
      if (nconv>0) {
	nblockmult=min(Integer(50),(mymaxj-neigfound)/blocksz-1);
      }
      else if ((neigfound==0)&&(iter%5==0)&&(nblockmult==(mymaxj-neigfound)/blocksz-1)){
	nchebit=10;
	if (sbs>=3)
	  maxval=d(neigfound+2);
	else
	  maxval=d(neigfound+1);
	minval=min(minval,d(neigfound+sbs-1));
      }
      else if (iter%2==0){
	  nblockmult=min(nblockmult+50,(mymaxj-neigfound)/blocksz-1);
      }
      continue;
    }
    */
      
 }while(1);

 rotate(neigfound-nconv,sbs,sbs,C,X,v);

 if ((myout)&&(print_level>0)){
   (*myout)<<"L__";myout->width(2);(*myout)<<iter<<":";
   (*myout)<<"  blocksz=";myout->width(2);(*myout)<<blocksz;
   (*myout)<<"  s=";myout->width(2);(*myout)<<nblockmult;
   (*myout)<<"  neigfound=";myout->width(2);(*myout)<<neigfound;
   (*myout)<<"  nmult=";myout->width(4);(*myout)<<nmult;
   (*myout)<<"  minval="<<minval<<std::endl;
   (*myout)<<"LE: blocksz="<<blocksz<<"\n  "<<neigfound<<":";
   for(Integer i=0;i<neigfound;i++) (*myout)<<" "<<d(i);
   (*myout)<<"\n other vecs("<<nlanczvecs-neigfound<<"):";
   {for(Integer i=neigfound;i<nlanczvecs;i++) (*myout)<<" "<<d(i);}
   (*myout)<<std::endl;
 }

 return 0;
}

// *****************************************************************************
//                                bklanc
// *****************************************************************************

//Input: neigfound,blocksz,s,d,X
//output: X,e,C
//usus: u,v

int Lanczpol::bklanc(Integer l_neigfound,Integer l_blocksz,Integer s,
                    const Matrix &l_d,Matrix &l_C,Matrix &l_X,
                    Matrix &l_e,Matrix &l_u,Matrix &l_v)
{
  if ((myout)&&(print_level>1)){
    (*myout)<<"start Lanczpol::bklanc"<<std::endl;
  }

 Integer n=l_X.rowdim();

 for(Integer l=0;l<s;l++){

     Integer ll=l_neigfound+l*l_blocksz;
     Integer lu=l_neigfound+(l+1)*l_blocksz;

     //---  l_v = A * l_X(:,block l)
     l_u.newsize(n,l_blocksz); chk_set_init(l_u,1);
     l_v.newsize(n,l_blocksz); chk_set_init(l_v,1);
     mat_xey(n*l_blocksz,l_u.get_store(),l_X.get_store()+ll*n);
     Microseconds secs=myclock.time();
     bigmatrix->lanczosmult(l_u,l_v);
     time_mult+=myclock.time()-secs; 


     //- compute l_C(block l)=l_v'*l_X(block l)
     if (l==0){
         for(Integer k=ll;k<lu;k++){    
             mat_xea(lu-k+1,l_C.get_store()+k*l_C.rowdim()+k+1,0.);
             l_C(k,k)=l_d(k);
         }
     }
     else {
         for(Integer k=ll;k<lu;k++){                  
             for(Integer i=k;i<lu;i++){
	         l_C(i,k)=mat_ip(n,l_v.get_store()+(k-ll)*n,l_X.get_store()+n*i);
             }
         }
     }
   
     if ((l>0)&&(l==s-1))
         break;

     //- compute residual l into l_X(block l+1)
     if (l==0){
         for(Integer k=ll;k<lu;k++){
             mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+n*k,-l_d(k));
         }
     }
     else {
         for (Integer k=ll;k<lu;k++){
             for (Integer i=ll;i<k;i++){
	          mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+i*n,-l_C(k,i));
             }
             for (Integer i=k;i<lu;i++){
	          mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+i*n,-l_C(i,k));
             }
         }
         for (Integer k=ll;k<lu;k++){
             for (Integer i=k;i<lu;i++){
	          mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+(i-l_blocksz)*n,-l_C(k,i-l_blocksz));
             }
         }
     }

     if ((l==0)&&(l==s-1)){
         for (Integer k=0;k<l_blocksz;k++){
              l_e(neigfound+k)= ::sqrt(mat_ip(n,l_v.get_store()+n*k,l_v.get_store()+n*k));
         }
         break;
     }
     else {
         mat_xey((lu-ll)*n,l_X.get_store()+n*lu,l_v.get_store());
     }

     if (l==0) err(l_neigfound,l_blocksz,l_X,l_e);
     if (l==s-1) continue;
     orthog(lu,l_blocksz,l_X,l_C);  //orthogonalize new vectors
     //--- copy R (l_C) of QR factorization to correct position
     Integer it=lu-1;
     for(Integer j=0;j<l_blocksz;j++){
         it++;
         //for(i=lu;i<=it;i++) l_C(i,it-l_blocksz)=l_C(i,it);
         mat_xey(it-lu+1,l_C.get_store()+(it-l_blocksz)*l_C.rowdim()+lu,
                 l_C.get_store()+it*l_C.rowdim()+lu);
     }
 }
 if ((myout)&&(print_level>1)){
   (*myout)<<"end Lanczpol::bklanc"<<std::endl;
 }

 return 0;
}


// *****************************************************************************
//                                bkqrlanc
// *****************************************************************************

//Input: neigfound,blocksz,s,d,X
//output: X,e,C
//usus: u,v

int Lanczpol::bkqrlanc(Integer l_neigfound,Integer l_blocksz,Integer s,
                    const Matrix &l_d,Matrix &l_C,Matrix &l_X,
                    Matrix &l_e,Matrix &l_u,Matrix &l_v)
{
 if ((myout)&&(print_level>1)){
   (*myout)<<"start Lanczpol::bkqrlanc"<<std::endl;
 }

 Integer n=l_X.rowdim();

 //--- start with a QR-factorization of the current columns of l_X

 //Xqr(:,0:l_neigfound+l_blocksz)
 Integer ncols=l_neigfound+l_blocksz;
 mat_xey(n*ncols,Xqr.get_store(),l_X.get_store());
 Matrix hfactor(l_neigfound+(s+1)*l_blocksz,1,0.);  //householder factor of each column

 for(Integer j=0;j<ncols;j++){

   //determine the householder vector for column j
   Real* xbase=Xqr.get_store()+j*n+j;
   Real Xjj=*(xbase);
   Real normj=mat_ip(n-j,xbase,xbase);   //will be used for householder vector norm
   Real mu=::sqrt(normj);
   Real beta=1.;
   if (mu>1e-10){
     beta=Xjj;
     if (beta<0) beta-=mu;
     else beta+=mu;
   }
   else {
     if (myout){
       (*myout)<<"*** ERROR in Lanczpol::bkqrlanc(): rank deficiency in first block"<<std::endl;
       (*myout)<<"n="<<n<<" ncols="<<ncols<<" l_blocksz="<<l_blocksz<<" l_neigfound="<<l_neigfound<<" j="<<j<<std::endl;
     }
     return 1;
   }
   // compute the squared norm of the householder vector
   Real hf=(normj-Xjj*Xjj)/(beta*beta)+1.;
   hf=-2./hf;
   hfactor(j)=hf;
   *xbase = Xjj+hf*(Xjj+(normj-Xjj*Xjj)/beta);
   mat_xmultea(n-j-1,xbase+1,1./beta);

   //rotate the remaining columns 
   //w=beta*transpose(A(i:m-1,j:n-1))*v(i:m-1);
   //A(i:m-1,j:n-1)+=v(i:m-1)*transpose(w);
   Real *colbase=xbase+n;
   for(Integer k=j+1;k<ncols;k++){
     Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
     *colbase+=d;
     mat_xpeya(n-j-1,colbase+1,xbase+1,d);
     colbase+=n;
   }
 }

 /*
 { // **** TEST
   Matrix QR(n,ncols,Xqr.get_store()); 
   Matrix Q1(n,ncols,l_X.get_store());
   Q1.QR_factor();
   cout<<"\n\n QR-diff="<<norm2(QR-Q1)<<std::endl;
 }
 */

 //--- compute s Lanczos steps
 for(Integer l=0;l<s;l++){

   //- call the (Chebycheff) matrix block multiplication  l_v=cheb(A)*l_X(block l) 
   Integer ll=l_neigfound+l*l_blocksz;
   Integer lu=l_neigfound+(l+1)*l_blocksz;
   l_u.newsize(n,l_blocksz); chk_set_init(l_u,1);
   l_v.newsize(n,l_blocksz); chk_set_init(l_v,1);
   mat_xey(n*l_blocksz,l_u.get_store(),l_X.get_store()+ll*n);
   Microseconds secs=myclock.time();
   bigmatrix->lanczosmult(l_u,l_v);
   time_mult+=myclock.time()-secs; 

   // **** TEST
   /*
   Matrix V(n,l_blocksz,l_v.get_store());
   */

   //- compute l_C(block l)=l_v'*l_X(block l)
   if (l==0){
     for(Integer k=ll;k<lu;k++){    
       mat_xea(lu-k+1,l_C.get_store()+k*l_C.rowdim()+k+1,0.);
       l_C(k,k)=l_d(k);
     }
   }
   else {
     for(Integer k=ll;k<lu;k++){                  
       for(Integer i=k;i<lu;i++){
	 l_C(i,k)=mat_ip(n,l_v.get_store()+(k-ll)*n,l_X.get_store()+n*i);
       }
     }
   }
   

   if ((l>0)&&(l==s-1))
     break;

   //- compute residual l into Xqr(block l+1)
   if (l==0){
     for(Integer k=ll;k<lu;k++){
       mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+n*k,-l_d(k));
     }
   }
   else {
     for(Integer k=ll;k<lu;k++){
       for (Integer i=ll;i<k;i++){
	 mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+i*n,-l_C(k,i));
       }
       for (Integer i=k;i<lu;i++){
	 mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+i*n,-l_C(i,k));
       }
     }
     for (Integer k=ll;k<lu;k++){
       for (Integer i=k;i<lu;i++){
	 mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+(i-l_blocksz)*n,-l_C(k,i-l_blocksz));
       }
     }
   }

   if ((l==0)&&(l==s-1)){
     for (Integer k=0;k<l_blocksz;k++){
       l_e(neigfound+k)= ::sqrt(mat_ip(n,l_v.get_store()+n*k,l_v.get_store()+n*k));
     }
     break;
   }
   else {
     mat_xey((lu-ll)*n,Xqr.get_store()+n*lu,l_v.get_store());
   }



   // **** TEST
   /*
   Matrix Ql(n,l_blocksz,l_X.get_store()+n*ll);
   Matrix C=tril(l_C(Range(ll,lu-1),Range(ll,lu-1)));
   if (l==0)
     C=triu(C);
   else if (l_blocksz>1) 
     C=C+transpose(tril(C,-1));
   cout<<"  C-diff="<<norm2(tril(C)-tril(transpose(V)*Ql))<<std::endl;
   Matrix R(n,blocksz,Xqr.get_store()+n*lu);
   if (l==0){
     cout<<" R-diff="<<norm2(R-V+Ql*C)<<std::endl;
   }
   else {
     Matrix B=triu(Xqr(Range(ll,lu-1),Range(ll,lu-1)));
     Matrix Qlm1(n,blocksz,l_X.get_store()+n*(ll-l_blocksz));
     cout<<" R-diff="<<norm2(R-V+Ql*C+Qlm1*transpose(B))<<std::endl;
   }
   if (l>0){
     C=triu(tril(l_C(Range(l_neigfound,lu-1),Range(l_neigfound,lu-1))),-l_blocksz);
     C=C+transpose(tril(C,-1));
   }
   Matrix Qall(n,lu,l_X.get_store());
   Matrix AQ(n,lu,0.);
   bigmatrix->lanczosmult(Qall,AQ);
   Matrix El(l_blocksz,lu,0.);
   for(Integer j=0;j<l_blocksz;j++)
     El(j,lu-l_blocksz+j)=1.;
   cout<<" Lanequation = "<<norm2(AQ-Qall*C-R*El)<<std::endl;
   */
   
   //- extend QR-factorization to new vectors in Xqr
   //transform the new columns by the Householder vectors of the r first columns
   for(Integer j=0;j<lu;j++){
     const Real* xbase=Xqr.get_store()+j*n+j;
     Real hf=hfactor(j);
     Real *colbase=Xqr.get_store()+lu*n+j;
     for(Integer k=lu;k<lu+l_blocksz;k++){
       Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
      *colbase+=d;
       mat_xpeya(n-j-1,colbase+1,xbase+1,d);
       colbase+=n;
     }
   }

   // **** TEST
   /*
   Matrix QR(n,lu,Xqr.get_store());
   Matrix R1(R);
   QR.Qt_times(R1,lu);
   cout<<" QR-extension1 ="<<norm2(R1-Matrix(n,l_blocksz,Xqr.get_store()+lu*n))<<std::endl;
   */

   //transform the new vectors by themselves
   for(Integer j=lu;j<lu+l_blocksz;j++){
     Real* xbase=Xqr.get_store()+j*n+j;
     Real Xjj=*(xbase);
     Real normj=mat_ip(n-j,xbase,xbase);   //will be used for householder vector norm
     Real mu=::sqrt(normj);
     if (l==0){
       l_e(j-l_blocksz)=mu;
     }
     Real beta=1.;
     if (mu>eps_Real){
       beta=Xjj;
       if (beta<0) beta-=mu;
       else beta+=mu;
     }
     // compute the squared norm of the householder vector
     Real hf=(normj-Xjj*Xjj)/(beta*beta)+1.;
     hf=-2./hf;
     hfactor(j)=hf;
     *xbase = Xjj+hf*(Xjj+(normj-Xjj*Xjj)/beta);
     mat_xmultea(n-j-1,xbase+1,1./beta);
     
     //rotate the remaining columns 
     //w=beta*transpose(A(i:m-1,j:n-1))*v(i:m-1);
     //A(i:m-1,j:n-1)+=v(i:m-1)*transpose(w);
     Real *colbase=xbase+n;
     for(Integer k=j+1;k<lu+l_blocksz;k++){
       Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
       *colbase+=d;
       mat_xpeya(n-j-1,colbase+1,xbase+1,d);
       colbase+=n;
     } 
   }

   /*
   // **** TEST
   //Indexmatrix piv(Range(0,lu-1));
   //QR.QR_concat_right(R,piv,lu);
   //cout<<" ext QR diff ="<<norm2(QR-Matrix(n,lu+l_blocksz,Xqr.get_store()))<<std::endl;
   QR.init(n,lu+l_blocksz,Xqr.get_store());
   Matrix Qt(n,n,0.);
   for(Integer j=0;j<n;j++) 
     Qt(j,j)=1.;
   QR.Qt_times(Qt,lu+l_blocksz);
   //cout<<" Qt = "<<norm2(Diag(Matrix(n,1,1.))-Qt*transpose(Qt))<<std::endl;
   cout<<" Qt*R diff ="<<norm2(Qt*R-triu(Matrix(n,l_blocksz,Xqr.get_store()+lu*n),-lu))<<std::endl;
   if (l==0) cout<<"l_e="<<l_e(Range(l_neigfound,l_neigfound+l_blocksz-1))<<std::endl;
   //cout<<" Qt*R ="<<Qt*R<<std::endl;
   */

   //- generate the new vectors in l_X
   mat_xea(n*l_blocksz,l_X.get_store()+lu*n,0.);
   for(Integer j=lu;j<lu+l_blocksz;j++)
     l_X(j,j)=1.;
   for(Integer j=lu+l_blocksz;--j>=0;){
     const Real* xbase=Xqr.get_store()+j*n+j;
     const Real hf=hfactor(j);
     Real *colbase=l_X.get_store()+lu*n+j;
     for(Integer k=lu;k<lu+l_blocksz;k++){
       const Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
       *colbase+=d;
       mat_xpeya(n-j-1,colbase+1,xbase+1,d);
       colbase+=n;
     }
   }
  
   /*
   // *** TEST
   Matrix Qnew(n,l_blocksz,0.);
   for (Integer j=lu;j<lu+l_blocksz;j++)
     Qnew(j,j-lu)=1.;
   QR.Q_times(Qnew,lu+l_blocksz);
   cout<<" Qnew diff = "<<norm2(Qnew-Matrix(n,l_blocksz,l_X.get_store()+lu*n))<<std::endl;
   Qnew.init(n,lu+l_blocksz,l_X.get_store());
   cout<<" orthog = "<<norm2(Diag(Matrix(lu+l_blocksz,1,1.))-transpose(Qnew)*Qnew)<<std::endl;
   */
     
   //- copy R (l_C) of QR factorization to correct position
   Integer it=lu-1;
   for(Integer j=0;j<l_blocksz;j++){
     it++;
     //for(i=lu;i<=it;i++) l_C(i,it-l_blocksz)=l_C(i,it);
     mat_xey(it-lu+1,l_C.get_store()+(it-l_blocksz)*l_C.rowdim()+lu,
	     Xqr.get_store()+it*n+lu);
   }

   // *** TEST
   //cout<<" CB = "<<triu(l_C(Range(lu,lu+l_blocksz-1),Range(lu-l_blocksz,lu-1)))<<std::endl;

 }

 if ((myout)&&(print_level>1)){
   (*myout)<<"end Lanczpol::bkqrlanc"<<std::endl;
 }

 return 0;
}


// *****************************************************************************
//                                bklanccheb
// *****************************************************************************

/** @brief block lanczos with Gram Schmidt orthogonalization 
 
  @param[in] l_neigfound number of eigenvalues found (first columns of l_X)
  
  @param[in] l_blocksz the block size of the block lanczos mutliplication

  @parma[in] s the number of block multiplications to be performed

  @param[in,out] l_X the matrix storing the Ritz vectors

  @param[out] l_e error vector consisting of    
  Input: l_neigfound,l_blocksz,s,l_X
  output: l_X,e,l_C
  uses: l_u,l_v
*/

int Lanczpol::bklanccheb(Integer l_neigfound,Integer l_blocksz,Integer s,
                    Matrix &l_C,Matrix &l_X,Matrix &l_e,Matrix &l_v)
{
 if ((myout)&&(print_level>1)){
   (*myout)<<"start Lanczpol::bklanccheb"<<std::endl;
 }

 Integer n=l_X.rowdim();

 for(Integer l=0;l<s;l++){

     Integer ll=l_neigfound+l*l_blocksz;
     Integer lu=l_neigfound+(l+1)*l_blocksz;
     blockcheby(ll,l_X,l_v);              // l_v <- cheb_polynom(Matrix) * l_X(block l)

     //- compute l_C(block l)=l_v'*l_X(block l)
     for(Integer k=ll;k<lu;k++){                  
         for(Integer i=k;i<lu;i++){
             l_C(i,k)=mat_ip(n,l_v.get_store()+(k-ll)*n,l_X.get_store()+n*i);
         }  
     }

     if ((l>0)&&(l==s-1))
       break;

     //- compute residual l into X(block l+1)
     mat_xey((lu-ll)*n,l_X.get_store()+n*lu,l_v.get_store());
     for(Integer k=ll;k<lu;k++){
         for (Integer i=ll;i<k;i++){
            mat_xpeya(n,l_X.get_store()+n*(k+l_blocksz),l_X.get_store()+i*n,-l_C(k,i));
         }
         for (Integer i=k;i<lu;i++){
            mat_xpeya(n,l_X.get_store()+n*(k+l_blocksz),l_X.get_store()+i*n,-l_C(i,k));
         }
     }
     if (l>0){
         for (Integer k=ll;k<lu;k++){
             for (Integer i=k;i<lu;i++){
	          mat_xpeya(n,l_X.get_store()+n*(k+l_blocksz),l_X.get_store()+(i-l_blocksz)*n,-l_C(i,k-l_blocksz));
             }
         }
     }

     // norm of residual vectors is used for convergence check later 
     if (l==0) err(l_neigfound,l_blocksz,l_X,l_e);
     if (l==s-1) continue;

     //orthogonalize new vectors to previous ones
     orthog(lu,l_blocksz,l_X,l_C);  
     //--- copy R (l_C) of QR factorization to correct position
     Integer it=lu-1;
     for(Integer j=0;j<l_blocksz;j++){
         it++;
         //for(i=lu;i<=it;i++) l_C(i,it-l_blocksz)=l_C(i,it);
         mat_xey(it-lu+1,l_C.get_store()+(it-l_blocksz)*l_C.rowdim()+lu,
                 l_C.get_store()+it*l_C.rowdim()+lu);
     }
 }
 if ((myout)&&(print_level>1)){
   (*myout)<<"end Lanczpol::bklanccheb"<<std::endl;
 }

 return 0;
}


// *****************************************************************************
//                                bkqrlanccheb
// *****************************************************************************

/** @brief block lanczos with Gram Schmidt orthogonalization 
 
  @param[in] l_neigfound number of eigenvalues found (first columns of l_X)
  
  @param[in] l_blocksz the block size of the block lanczos mutliplication

  @parma[in] s the number of block multiplications to be performed

  @param[in,out] l_X the matrix storing the Ritz vectors

  @param[out] l_e error vector consisting of    
  Input: l_neigfound,l_blocksz,s,l_X
  output: l_X,e,l_C
  uses: l_u,l_v
*/

int Lanczpol::bkqrlanccheb(Integer l_neigfound,Integer l_blocksz,Integer s,
                    Matrix &l_C,Matrix &l_X,Matrix &l_e,Matrix &l_v)
{
 if ((myout)&&(print_level>1)){
   (*myout)<<"start Lanczpol::bkqrlanccheb"<<std::endl;
 }

 Integer n=l_X.rowdim();

 //--- start with a QR-factorization of the current columns of l_X

 //Xqr(:,0:l_neigfound+l_blocksz)
 Integer ncols=l_neigfound+l_blocksz;
 mat_xey(n*ncols,Xqr.get_store(),l_X.get_store());
 Matrix hfactor(l_neigfound+(s+1)*l_blocksz,1,0.);  //householder factor of each column

 for(Integer j=0;j<ncols;j++){

   //determine the householder vector for column j
   Real* xbase=Xqr.get_store()+j*n+j;
   Real Xjj=*(xbase);
   Real normj=mat_ip(n-j,xbase,xbase);   //will be used for householder vector norm
   Real mu=::sqrt(normj);
   Real beta=1.;
   if (mu>1e-10){
     beta=Xjj;
     if (beta<0) beta-=mu;
     else beta+=mu;
   }
   else {
     if (myout){
       (*myout)<<"*** ERROR in Lanczpol::bkqrlanccheb(): rank deficiency in first block"<<std::endl;
       (*myout)<<"n="<<n<<" ncols="<<ncols<<" l_blocksz="<<l_blocksz<<" l_neigfound="<<l_neigfound<<" j="<<j<<std::endl;
     }
     return 1;
   }
   // compute the squared norm of the householder vector
   Real hf=(normj-Xjj*Xjj)/(beta*beta)+1.;
   hf=-2./hf;
   hfactor(j)=hf;
   *xbase = Xjj+hf*(Xjj+(normj-Xjj*Xjj)/beta);
   mat_xmultea(n-j-1,xbase+1,1./beta);

   //rotate the remaining columns 
   //w=beta*transpose(A(i:m-1,j:n-1))*v(i:m-1);
   //A(i:m-1,j:n-1)+=v(i:m-1)*transpose(w);
   Real *colbase=xbase+n;
   for(Integer k=j+1;k<ncols;k++){
     Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
     *colbase+=d;
     mat_xpeya(n-j-1,colbase+1,xbase+1,d);
     colbase+=n;
   }
 }

 /*
 { // **** TEST
   Matrix QR(n,ncols,Xqr.get_store()); 
   Matrix Q1(n,ncols,l_X.get_store());
   Q1.QR_factor();
   cout<<"\n\n QR-diff="<<norm2(QR-Q1)<<std::endl;
 }
 */
   
 //--- compute s Lanczos steps
 for(Integer l=0;l<s;l++){

   //- call the (Chebycheff) matrix block multiplication  l_v=cheb(A)*l_X(block l) 
   Integer ll=l_neigfound+l*l_blocksz;
   Integer lu=l_neigfound+(l+1)*l_blocksz;
   blockcheby(ll,l_X,l_v);              

   /*
   // **** TEST
   Matrix V(n,l_blocksz,l_v.get_store());
   */

   //- compute l_C(block l)=l_v'*l_X(block l)
   for(Integer k=ll;k<lu;k++){                  
     for(Integer i=k;i<lu;i++){
       l_C(i,k)=mat_ip(n,l_v.get_store()+(k-ll)*n,l_X.get_store()+n*i);
     }
   }

   if ((l>0)&&(l==s-1))
     break;

   //- compute residual l into Xqr(block l+1)
   for(Integer k=ll;k<lu;k++){
     for (Integer i=ll;i<k;i++){
       mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+i*n,-l_C(k,i));
     }
     for (Integer i=k;i<lu;i++){
       mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+i*n,-l_C(i,k));
     }
   }
   if (l>0){
     for (Integer k=ll;k<lu;k++){
       for (Integer i=k;i<lu;i++){
	 mat_xpeya(n,l_v.get_store()+n*(k-ll),l_X.get_store()+(i-l_blocksz)*n,-Xqr(k,i));
       }
     }
   }
   
   if ((l==0)&&(l==s-1)){
     for (Integer k=0;k<l_blocksz;k++){
       l_e(neigfound+k)= ::sqrt(mat_ip(n,l_v.get_store()+n*k,l_v.get_store()+n*k));
     }
     break;
   }
   else {
     mat_xey((lu-ll)*n,Xqr.get_store()+n*lu,l_v.get_store());
   }

   /*
   // **** TEST
   Matrix Ql(n,l_blocksz,l_X.get_store()+n*ll);
   Matrix C=tril(l_C(Range(ll,lu-1),Range(ll,lu-1)));
   if (l_blocksz>1) 
     C=C+transpose(tril(C,-1));
   cout<<"  C-diff="<<norm2(tril(C)-tril(transpose(V)*Ql))<<std::endl;
   Matrix R(n,blocksz,Xqr.get_store()+n*lu);
   if (l==0){
     cout<<" R-diff="<<norm2(R-V+Ql*C)<<std::endl;
   }
   else {
     Matrix B=triu(Xqr(Range(ll,lu-1),Range(ll,lu-1)));
     Matrix Qlm1(n,blocksz,l_X.get_store()+n*(ll-l_blocksz));
     cout<<" R-diff="<<norm2(R-V+Ql*C+Qlm1*transpose(B))<<std::endl;
   }
   if (l>0){
     C=triu(tril(l_C(Range(l_neigfound,lu-1),Range(l_neigfound,lu-1))),-l_blocksz);
     C=C+transpose(tril(C,-1));
   }
   Matrix Qall(n,lu,l_X.get_store());
   Matrix AQ(n,lu,0.);
   bigmatrix->lanczosmult(Qall,AQ);
   Matrix El(l_blocksz,lu,0.);
   for(Integer j=0;j<l_blocksz;j++)
     El(j,lu-l_blocksz+j)=1.;
   cout<<" Lanequation = "<<norm2(AQ-Qall*C-R*El)<<std::endl;
   */
   
   //- extend QR-factorization to new vectors in Xqr
   //transform the new columns by the Householder vectors of the r first columns
   for(Integer j=0;j<lu;j++){
     const Real* xbase=Xqr.get_store()+j*n+j;
     Real hf=hfactor(j);
     Real *colbase=Xqr.get_store()+lu*n+j;
     for(Integer k=lu;k<lu+l_blocksz;k++){
       Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
      *colbase+=d;
       mat_xpeya(n-j-1,colbase+1,xbase+1,d);
       colbase+=n;
     }
   }

   /*
   // **** TEST
   Matrix QR(n,lu,Xqr.get_store());
   Matrix R1(R);
   QR.Qt_times(R1,lu);
   cout<<" QR-extension1 ="<<norm2(R1-Matrix(n,l_blocksz,Xqr.get_store()+lu*n))<<std::endl;
   */

   //transform the new vectors by themselves
   for(Integer j=lu;j<lu+l_blocksz;j++){
     Real* xbase=Xqr.get_store()+j*n+j;
     Real Xjj=*(xbase);
     Real normj=mat_ip(n-j,xbase,xbase);   //will be used for householder vector norm
     Real mu=::sqrt(normj);
     if (l==0){
       l_e(j-l_blocksz)=mu;
     }
     Real beta=1.;
     if (mu>eps_Real){
       beta=Xjj;
       if (beta<0) beta-=mu;
       else beta+=mu;
     }
     // compute the squared norm of the householder vector
     Real hf=(normj-Xjj*Xjj)/(beta*beta)+1.;
     hf=-2./hf;
     hfactor(j)=hf;
     *xbase = Xjj+hf*(Xjj+(normj-Xjj*Xjj)/beta);
     mat_xmultea(n-j-1,xbase+1,1./beta);
     
     //rotate the remaining columns 
     //w=beta*transpose(A(i:m-1,j:n-1))*v(i:m-1);
     //A(i:m-1,j:n-1)+=v(i:m-1)*transpose(w);
     Real *colbase=xbase+n;
     for(Integer k=j+1;k<lu+l_blocksz;k++){
       Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
       *colbase+=d;
       mat_xpeya(n-j-1,colbase+1,xbase+1,d);
       colbase+=n;
     } 
   }

   /*
   // **** TEST
   //Indexmatrix piv(Range(0,lu-1));
   //QR.QR_concat_right(R,piv,lu);
   //cout<<" ext QR diff ="<<norm2(QR-Matrix(n,lu+l_blocksz,Xqr.get_store()))<<std::endl;
   QR.init(n,lu+l_blocksz,Xqr.get_store());
   Matrix Qt(n,n,0.);
   for(Integer j=0;j<n;j++) 
     Qt(j,j)=1.;
   QR.Qt_times(Qt,lu+l_blocksz);
   //cout<<" Qt = "<<norm2(Diag(Matrix(n,1,1.))-Qt*transpose(Qt))<<std::endl;
   cout<<" Qt*R diff ="<<norm2(Qt*R-triu(Matrix(n,l_blocksz,Xqr.get_store()+lu*n),-lu))<<std::endl;
   if (l==0) cout<<"l_e="<<l_e(Range(l_neigfound,l_neigfound+l_blocksz-1))<<std::endl;
   //cout<<" Qt*R ="<<Qt*R<<std::endl;
   */

   //- generate the new vectors in l_X
   mat_xea(n*l_blocksz,l_X.get_store()+lu*n,0.);
   for(Integer j=lu;j<lu+l_blocksz;j++)
     l_X(j,j)=1.;
   for(Integer j=lu+l_blocksz;--j>=0;){
     const Real* xbase=Xqr.get_store()+j*n+j;
     const Real hf=hfactor(j);
     Real *colbase=l_X.get_store()+lu*n+j;
     for(Integer k=lu;k<lu+l_blocksz;k++){
       const Real d=hf*(*colbase+mat_ip(n-j-1,colbase+1,xbase+1));
       *colbase+=d;
       mat_xpeya(n-j-1,colbase+1,xbase+1,d);
       colbase+=n;
     }
   }

   /*
   // *** TEST
   Matrix Qnew(n,l_blocksz,0.);
   for (Integer j=lu;j<lu+l_blocksz;j++)
     Qnew(j,j-lu)=1.;
   QR.Q_times(Qnew,lu+l_blocksz);
   cout<<" Qnew diff = "<<norm2(Qnew-Matrix(n,l_blocksz,l_X.get_store()+lu*n))<<std::endl;
   Qnew.init(n,lu+l_blocksz,l_X.get_store());
   cout<<" orthog = "<<norm2(Diag(Matrix(lu+l_blocksz,1,1.))-transpose(Qnew)*Qnew)<<std::endl;
   */
     
   //- copy R (l_C) of QR factorization to correct position
   Integer it=lu-1;
   for(Integer j=0;j<l_blocksz;j++){
     it++;
     //for(i=lu;i<=it;i++) l_C(i,it-l_blocksz)=l_C(i,it);
     mat_xey(it-lu+1,l_C.get_store()+(it-l_blocksz)*l_C.rowdim()+lu,
	     Xqr.get_store()+it*n+lu);
   }

   // *** TEST
   //cout<<" CB = "<<triu(l_C(Range(lu,lu+l_blocksz-1),Range(lu-l_blocksz,lu-1)))<<std::endl;

 }

 if ((myout)&&(print_level>1)){
   (*myout)<<"end Lanczpol::bkqrlanccheb"<<std::endl;
 }

 return 0;
}


// *****************************************************************************
//                                err
// *****************************************************************************

// e(k)=norm_2(l_X(:,k))   k=l_neigfound+l_blocksz:l_neigfound+2*l_blocksz-1;
//Input:  l_neigfound,l_blocksz,l_X
//Output: e

int Lanczpol::err(Integer l_neigfound,Integer l_blocksz,const Matrix& l_X,Matrix &l_e)
{
 Integer n=l_X.rowdim();
 for(Integer k=l_neigfound+l_blocksz;k<l_neigfound+2*l_blocksz;k++){
     //t=0.; for(i=0;i<n;i++) t+=l_X(i,k)*l_X(i,k);
     Real t=mat_ip(n,l_X.get_store()+n*k,l_X.get_store()+n*k);
     l_e(k-l_blocksz)=::sqrt(t);
 }
 return 0;
}

// *****************************************************************************
//                                cnvtst
// *****************************************************************************

//Input: l_neigfound, l_blocksz,errc,eps,d,e
//output: nconv,errc

int Lanczpol::cnvtst(Integer l_neigfound,Integer l_blocksz,Real& l_errc,Real l_eps,
                    const Matrix& l_d,const Matrix &l_e,Integer &nconv)
{
 //const Real mcheps=2.22e-16;
 Integer n=l_d.dim();
 Integer i,k=0;
 Real t=1.;
 if ((myout)&&(print_level>0))
   (*myout)<<" cnvtst:";
 for(i=0;i<l_neigfound;i++){t=max(t,fabs(l_d(i)));}
 for(i=0;i<l_blocksz;i++){
     t=max(t,fabs(l_d(i+l_neigfound)));
     if ((myout)&&(print_level>0))
       (*myout)<<" ("<<i+l_neigfound<<","<<l_e(i+l_neigfound)<<"<"<<t*(l_eps+10.*Real(n)*mcheps)+l_errc<<")";
     if (l_e(i+l_neigfound)>t*(l_eps+10.*Real(n)*mcheps)+l_errc) break;
     k=i+1;
 }
 if ((myout)&&(print_level>0))
   (*myout)<<std::endl;
 nconv=k;
 if (k==0) return 0;
 t=0.;
 for (i=0;i<k;i++){
     t+=l_e(i+l_neigfound)*l_e(i+l_neigfound);
 }
 l_errc=::sqrt(l_errc*l_errc+t);
 return 0;
}
               
 
// *****************************************************************************
//                                eigen
// *****************************************************************************

//Input: l_neigfound,l_blocksz,sbsz,l_C,minval
//Output: l_C,d,minval

int Lanczpol::eigen(Integer l_neigfound,Integer l_blocksz,Integer sbsz,Matrix& l_C,
                   Matrix& l_d,Matrix& l_u,Matrix& l_v,Real& l_minval)
{
 Integer i,j;
 Integer lim;
 for(i=0;i<sbsz;i++){
     lim=i-l_blocksz;
     if (i<l_blocksz) lim=0;
     for(j=0;j<lim;j++)             
         l_C(i,j)=0.;               //l_C(i,0:lim)=0
     for(j=lim;j<=i;j++)
         l_C(i,j)=l_C(i+l_neigfound,j+l_neigfound);
 }
 tred2(sbsz,l_C,l_u,l_v,l_C);
 if ((i=tql2(sbsz,l_u,l_v,l_C))){
     (*myout)<<"\n"<<i<<std::endl;
     error("eigen failed on computing tql2");
 }
 for(i=0;i<sbsz;i++){
     l_d(i+l_neigfound)=u(i);
 }
 l_minval=min(l_minval,l_d(l_neigfound+sbsz-1));
 return 0;
}

// *****************************************************************************
//                                sectn
// *****************************************************************************



//Input: l_X,l_neigfound,l_blocksz,minval,
//Output: l_X,minval,d

int Lanczpol::sectn(Matrix& l_X,Integer l_neigfound,Integer l_blocksz,
                   Matrix& l_C,Matrix& l_d,Matrix& l_u,Matrix& l_v,
                   Real& l_minval)
{
 Integer n=l_X.rowdim();

 //---  l_C=l_X'*A*l_X
 //    lower triangle of l_C(l_neigfound+1:l_neigfound+l_blocksz,l_neigfound:nef+l_blocksz)
 //    l_X is only l_X(:,l_neigfound+1,l_neigfound:l_blocksz)
 Integer col1,col2;
 col1=l_neigfound-1;
 u.newsize(n,l_blocksz);
 l_v.newsize(n,l_blocksz);
 w.newsize(n,1);
 chk_set_init(l_u,1);
 chk_set_init(l_v,1);
 chk_set_init(w,1);
 mat_xey(n*l_blocksz,l_u.get_store(),l_X.get_store()+l_neigfound*n);
 Microseconds secs=myclock.time();
 bigmatrix->lanczosmult(l_u,l_v);
 time_mult+=myclock.time()-secs;
 for(Integer j=0;j<l_blocksz;j++){
     col1++;
     col2=l_neigfound-1;
     for(Integer i=0;i<=j;i++){
         col2++;
         l_C(col1,col2)=mat_ip(n,l_v.get_store()+j*n,l_X.get_store()+col2*n);
         if (i==j){
             mat_xeya(n,w.get_store(),l_X.get_store()+col2*n,l_C(col1,col2));
             mat_xmey(n,w.get_store(),l_v.get_store()+j*n);
             e(l_neigfound+i)=norm2(w);
         }
             
     }
 }

 //--- diagonalize shifted l_C
 //==> diag d, l_C trafo matrix, minval updated small eig est 
 eigen(l_neigfound,l_blocksz,l_blocksz,l_C,l_d,l_u,l_v,l_minval);
 return 0;
}


// *****************************************************************************
//                                rotate_extremes
// *****************************************************************************

//rotates the eigenvectors of the largest and smallest eigenvalues of
//the tridiagonal matrix.
// l_X(:,l_neigfound[+0,+sbs-1])=l_X(:,l_neigfound:l_neigfound+sbs-1)*l_C(:,[0,sbs-1])

// Input: l_neigfound,sbs,l_C,l_X
// Output: l_X
// uses: l_v as storage

int Lanczpol::rotate_extremes(Integer l_neigfound,Integer sbs,Matrix& l_d,const Matrix& l_C,Matrix &l_X,Matrix& l_v)
{
 Integer q=l_C.rowdim();
 Integer n=l_X.rowdim();
 for(Integer i=0;i<n;i++){ //for each row
     Real *xi=l_X.get_store()+l_neigfound*n+i;
     l_v(0)=mat_ip(sbs,xi,n,l_C.get_store(),1);            //rotate largest EV
     l_v(1)=mat_ip(sbs,xi,n,l_C.get_store()+(sbs-1)*q,1);  //rotate smallest EV
     mat_xey(2,xi,n,l_v.get_store(),1);
 }
 l_d(l_neigfound+1)=l_d(l_neigfound+sbs-1);
 return 0;
}
 
           
// *****************************************************************************
//                                rotate
// *****************************************************************************

// l_X(:,l_neigfound:l_neigfound+l-1)=l_X(:,l_neigfound:l_neigfound+sbs-1)*l_C

// Input: l_neigfound,sbs,l,l_C,l_X
// Output: l_X
// uses: v as storage

int Lanczpol::rotate(Integer l_neigfound,Integer sbs,Integer l,
                    const Matrix& l_C,Matrix &l_X,Matrix& l_v)
{
 Integer q=l_C.rowdim();
 Integer n=l_X.rowdim();
 for(Integer i=0;i<n;i++){
     Real *xi=l_X.get_store()+l_neigfound*n+i; Real *vp=l_v.get_store();const Real *cp=l_C.get_store();
     for(Integer k=0;k<l;k++){
         //t=0.; for(j=0;j<sbs;j++) t+=l_X(i,j+l_neigfound)*l_C(j,k);
         *vp++=mat_ip(sbs,xi,n,cp,1); cp+=q;
     }
     //for(k=0;k<l;k++) l_X(i,k+l_neigfound)=l_v(k);
     mat_xey(l,xi,n,l_v.get_store(),1);
 }
 return 0;
}
 
           
// *****************************************************************************
//                                random
// *****************************************************************************

//assigns a random vector to column j of l_X

//Input: l_X,j,randomseed
//Output: l_X,randomseed

int Lanczpol::random(Matrix& l_X,Integer j)
{
 Integer n=l_X.rowdim();
 for(Integer i=0;i<n;i++){
     l_X(i,j)=2.*randgen.next()-1.;
 }
 return 0;
}

// *****************************************************************************
//                                orthog
// *****************************************************************************

// orthonormalize columns l_X(:,offset:offset+l_blocksz-1) also with respect to
// the first offset columns of l_X (these are already assumed to be orthonormal). 
// On output the upper triangle of B(:,offset:offset+l_blocksz-1) 
// contains the coefficients of the 
// components (B is upper triangular with diagonal but only for the new vectors)
// uses Gram-Schmidt for the QR factorization
// returns rank deficiency

//Input: l_X,l_blocksz,offset
//Output: l_X,B

Integer Lanczpol::orthog(Integer offset,Integer l_blocksz,Matrix &l_X,Matrix& B)
{
 Integer n=l_X.rowdim();
 int orig;
 Integer  rankdef=0;
 for(Integer k=offset;k<offset+l_blocksz;k++){
     Real *xk=l_X.get_store()+k*l_X.rowdim();
     orig=1;
     Real t,t1;
     do{
         // compute projection on l_X(:,i) and subtract this 
         t=0.;
	 Real *xi=l_X.get_store();
         for(Integer i=0;i<k-rankdef;i++){
             //t1=0.; for(j=0;j<n;j++) t1+=l_X(j,i)*l_X(j,k);
             t1=mat_ip(n,xi,xk);
             if ((orig)&&(i>=offset)) B(i,k)=t1;
             t+=t1*t1;
             //for(j=0;j<n;j++) l_X(j,k)-=t1*l_X(j,i);
             mat_xpeya(n,xk,xi,-t1);
             xi+=n;
         }
         //t1=0.;for(j=0;j<n;j++) t1+=l_X(j,k)*l_X(j,k);
         t1=mat_ip(n,xk,xk);
         t+=t1;
         orig=0;
     }while(t1<t/100.);
     t1=::sqrt(t1);
     for(Integer i=k-rankdef+1;i<=k;i++) B(i,k)=0.;
     if (t1>eps_Real) {
         B(k-rankdef,k)=t1;
         t1=1./t1;
         //for(j=0;j<n;j++) l_X(j,k-rankdef)=l_X(j,k)*t1;
         if (rankdef==0) {
	     mat_xmultea(n,xk,t1);
         }
         else {
	     mat_xeya(n,l_X.get_store()+(k-rankdef)*l_X.rowdim(),xk,t1);
         }
     }
     else { //residual is zero -> rank deficiency, eliminate column
         B(k-rankdef,k)=0.;
         rankdef++;
     }
 }  
 //fill up rank deficiency by random columns and orthonormalize them
 {for(Integer k=offset+l_blocksz-rankdef;k<min(n,offset+l_blocksz);k++){
     random(l_X,k);
     Real *xk=l_X.get_store()+k*l_X.rowdim();
     orig=1;
     Real t,t1;
     do{
         // compute projection on l_X(:,i) and subtract this 
         t=0.;
	 Real *xi=l_X.get_store();
         for(Integer i=0;i<k;i++){
             //t1=0.; for(j=0;j<n;j++) t1+=l_X(j,i)*l_X(j,k);
             t1=mat_ip(n,xi,xk);
             t+=t1*t1;
             //for(j=0;j<n;j++) l_X(j,k)-=t1*l_X(j,i);
             mat_xpeya(n,xk,xi,-t1);
             xi+=n;
         }
         //t1=0.;for(j=0;j<n;j++) t1+=l_X(j,k)*l_X(j,k);
         t1=mat_ip(n,xk,xk);
         t+=t1;
         orig=0;
     }while(t1<t/100.);
     t1=::sqrt(t1);
     if (t1>eps_Real) {
         B(k,k)=0.;
         t1=1./t1;
         //for(j=0;j<n;j++) l_X(j,k)*=t1;
	 mat_xmultea(n,xk,t1);
     }
     else { //random column was again rank deficient, try again
       //(k<n) is guaranteed in the intilialization of the parameters
       k--;
     }
 }}  
 return rankdef;
}      

// *****************************************************************************
//                                blockcheby
// *****************************************************************************

//The routine takes the l_blocksz first vectors of l_X after colum col_offset
//and applies the matrix cheb polynomial of order nchebit with respect to the
//interval [minval, maxval] to these vectors. The result is returned in v.
//Input: l_X, col_offset
//global params: maxval,minval,nchebit, l_blocksz, polval
//output: v
//uses: u,w

int Lanczpol::blockcheby(Integer col_offset,const Matrix& l_X,Matrix& l_v)
{
 //T0: u=l_X(:,col_offset:col_offset+blocksz-1)
 Integer n=l_X.rowdim();
 u.newsize(n,blocksz); chk_set_init(u,1);
 l_v.newsize(n,blocksz); chk_set_init(l_v,1);
 mat_xey(n*blocksz,u.get_store(),l_X.get_store()+n*col_offset);
 Microseconds secs=myclock.time();
 bigmatrix->lanczosmult(u,l_v);
 time_mult+=myclock.time()-secs;
 if ((maxval-minval<1e-8)||(nchebit<=0)){
     return 0;
 }

 //determine interval limits [a,b] for polynomial
 if ((nchebit%2)==0) nchebit++;
 Real a=minval-1e-3*(maxval-minval);
 Real b=maxval;
 //Real b=(2*maxval+(cosh(acosh(polval)/nchebit)-1)*a)/(cosh(acosh(polval)/nchebit)+1);
 
 //shift away known eigenvalues
 for(Integer i=0;i<neigfound;i++){
   for(Integer k=0;k<blocksz;k++){
     Real dip=mat_ip(n,u.get_store()+k*n,l_X.get_store()+i*n);
     mat_xpeya(n,l_v.get_store()+k*n,l_X.get_store()+i*n,(-d(i)+(a+b)/2.)*dip);
   }
 }
 
 //T1: l_v=(a+b)/(a-b)*u-2./(a-b)*l_v
 mat_xmultea(n*blocksz,l_v.get_store(),-2./(a-b));
 mat_xpeya(n*blocksz,l_v.get_store(),u.get_store(),(a+b)/(a-b));
 Real c1=2.*(a+b)/(a-b);
 Real c2=4./(a-b);

 //T_n = c1*T_{n-1} - c2 * bigmatrix * T_{n-1} - T_{n-2}
 w.newsize(n,blocksz); chk_set_init(w,1);
 Integer j;
 for(j=2;j<nchebit;j+=2){
     secs=myclock.time();
     bigmatrix->lanczosmult(l_v,w);
     time_mult+=myclock.time()-secs;
     for(Integer i=0;i<neigfound;i++){
       for(Integer k=0;k<blocksz;k++){
	 Real dip=mat_ip(n,l_v.get_store()+k*n,l_X.get_store()+i*n);
	 mat_xpeya(n,w.get_store()+k*n,l_X.get_store()+i*n,(-d(i)+(a+b)/2.)*dip);
       }
     }
     //T_{2n}= u = c1*v-c2*w-u;
     mat_xemx(n*blocksz,u.get_store());
     mat_xpeya(n*blocksz,u.get_store(),w.get_store(),-c2);
     mat_xpeya(n*blocksz,u.get_store(),l_v.get_store(),c1);
     secs=myclock.time();
     bigmatrix->lanczosmult(u,w);
     time_mult+=myclock.time()-secs;
     {for(Integer i=0;i<neigfound;i++){
       for(Integer k=0;k<blocksz;k++){
	 Real dip=mat_ip(n,u.get_store()+k*n,l_X.get_store()+i*n);
	 mat_xpeya(n,w.get_store()+k*n,l_X.get_store()+i*n,(-d(i)+(a+b)/2.)*dip);
       }
     }}
     //T_{2n+1}= l_v = c1*u-c2*w-l_v;
     mat_xemx(n*blocksz,l_v.get_store());
     mat_xpeya(n*blocksz,l_v.get_store(),w.get_store(),-c2);
     mat_xpeya(n*blocksz,l_v.get_store(),u.get_store(),c1);
 }
 return 0;
}

// *****************************************************************************
//                                scalarcheby
// *****************************************************************************

//computes the same l_Chebycheff polynomial as blockcheby but for scalar values
//Input: xval, col_offset
//global params: maxval,minval,nchebit, blocksz, polval
//output: v
//uses: u,w

Real Lanczpol::scalarcheby(Real xval)
{
 if ((maxval-minval<1e-8)||(nchebit==0)){
     return xval;
 }
 //determine interval limits [a,b] for polynomial
 if ((nchebit%2)==0) nchebit++;
 Real a=minval-1e-3*(maxval-minval);
#ifndef __unix
 Real coshpolval=cosh(::log(double(1./(polval-::sqrt(double(polval*polval-1)))))/nchebit);
#else
 Real coshpolval=cosh(acosh(polval)/nchebit);
#endif
 Real b=(2*maxval+(coshpolval-1.)*a)/(coshpolval/nchebit+1.);

 Real ru,rv;
 //T0= ru = 1;
 ru=1;

 //T1: v=(a+b)/(a-b)-2./(a-b)*xval
 rv=(a+b)/(a-b)-2./(a-b)*xval;

 Real c1=2.*(a+b)/(a-b);
 Real c2=4./(a-b);
 //T_n = c1*T_{n-1} - c2 * xval * T_{n-1} - T_{n-2}
 Integer j;
 for(j=2;j<nchebit;j+=2){
     //T_{2n}= u = v*(c1-c2*xval)-u;
     ru=rv*(c1-c2*xval)-ru;

     //T_{2n+1}= v = c1*u-c2*w-v;
     rv=ru*(c1-c2*xval)-rv;
 }
 return rv;
}

// *****************************************************************************
//                                tred2
// *****************************************************************************

//transform symmetric input matrix l_C (stored in the lower triagonal) into
//tridiagonal form (diagonal u and subdiagonal v), the orthogonal 
//transformation matrix is stored in Z; l_C and Z may be the same  
//the subdiagonal elements v(i) are stored in positions 1 to n,
//v(0) is set to zero.

//uses householder transformation

int Lanczpol::tred2(Integer n,const Matrix& l_C,Matrix& l_u,Matrix& l_v,Matrix& Z)
{
 Integer i,j,k,l;
 l_u.init(n,1,0.);
 l_v.init(n,1,0.);
 
 //---- copy lower triangular part from l_C to Z (may be the same)
 if (Z.get_store()!=l_C.get_store()){
     for(i=0;i<n;i++)
         for(j=0;j<=i;j++)
             Z(i,j)=l_C(i,j);
 }
 
 //---- compute householder recursively starting with the last column 
 for(i=n-1;i>0;i--){
     l=i-1;             
     Real h=0.;
     Real scale=0.;
     if (l<1){
         l_v(i)=Z(i,l);
         l_u(i)=h;
         continue;
     }
     for(k=0;k<=l;k++)
         scale+=fabs(Z(i,k));
     if(scale==0.){
         l_v(i)=Z(i,l);
         l_u(i)=h;
         continue;
     }                                //scale=sum(abs(Z(i,0:l)) != 0.

     //---- compute Z(i,0:l) /= scale
     //---- compute h = Z(i,0:l)*Z(i,0:l)'
     for(k=0;k<=l;k++){               //scaling Z(i,0:l)
         Z(i,k)/=scale;
         h+=Z(i,k)*Z(i,k);
     }

     //---- compute householder vector into Z(i,0:l)
     //       Z(i,0:l) = Z(i,0:l) + sign(Z(i,l))*norm(Z(i,0:l))*e_l;
     //       h = 2*Z(i,0:l)*Z(i,0:l)';
     Real f=Z(i,l);
     Real g=-d_sign(::sqrt(h),f);       //g=-sign(f)*norm(Z(i,0:l))
     l_v(i)=scale*g;                    //offdiagonal element is now fixed
     h-=f*g;                         
     Z(i,l)=f-g;

     //---- compute l_v(0:l) = Z(0:l,0:l) * Z(i,0:l)' /h;
     //                  f = l_v(0:l)' * Z(i,0:l)';
     f=0.;
     Real *zi=Z.get_store()+i;
     Real *zj=Z.get_store();
     Integer Zrowdim=Z.rowdim();
     for(j=0;j<=l;j++,zj++){
         //Z(j,i)=Z(i,j)/(scale*h);      // stores Z(i,0:l) in upper triangle
         *(zj+i*Zrowdim)=*(zi+j*Zrowdim)/(scale*h);
         //g=0.;
         //for(k=0;k<=j;k++) g+=Z(j,k)*Z(i,k);
         //for(;k<=l;k++)    g+=Z(k,j)*Z(i,k);
         g= mat_ip(j+1,zj,Zrowdim,zi,Zrowdim);
         g+=mat_ip(l-j,zj+j*Zrowdim+1,1,
                   zi+(j+1)*Zrowdim,Zrowdim);
         //l_v(j)=g/h;
         //f+=l_v(j)*Z(i,j);
         f+=(l_v(j)=g/h)*(*(zi+j*Zrowdim));
     }

     //---- compute   l_v = l_v - Z(i,0:l)'*f/2h;   (second vector for rank2-update)
     //---- compute   Z(0:l,0:l) -= Z(i,0:l)'*l_v'+ l_v*Z(i,0:l);     (rank2-update)
     Real hh=f/(h+h);
     for(j=0;j<=l;j++){
         zj=Z.get_store()+j;
         f=Z(i,j);
         g=(l_v(j)-=hh*f);
         //for(k=0;k<=j;k++) Z(j,k)-=f*v(k)+g*Z(i,k);
         mat_xpeya(j+1,zj,Z.rowdim(),zi,Z.rowdim(),-g);
         mat_xpeya(j+1,zj,Z.rowdim(),l_v.get_store(),1,-f);
     } 

     //---- unscale Z(i,0:l);
     //---- store squared norm of householdervec in u(i)
     //for(k=0;k<=l;k++) Z(i,k)*=scale;
     mat_xmultea(l+1,zi,Z.rowdim(),scale);
     l_u(i)=h;
 }
 l_u(0)=0.;
 l_v(0)=0.;

 //---- reconstruct orthogonal transformation by succesively computing
 //----  Z= Z*(I-2l_vl_v'/l_v'l_v) for the householder vectors l_v (starting with Z=I)

 Real g;
 
 for(i=0;i<n;i++){
     l=i-1;
     if (l_u(i)!=0.){ //if the householder vector is not the zero vector
         //--- compute Z(0:l,0:l) -= Z(i,0:l)*Z(0:l,0:l)*Z(0:l,j)
         //---                     (== z'*Z*z/scale/scale/(2*z'*z)) 
         for(j=0;j<=l;j++){
             //--- g= Z(0:l,j)*Z(i,0:l)
             //g=0.; for(k=0;k<=l;k++) g+= Z(i,k)*Z(k,j);
             g=mat_ip(l+1,Z.get_store()+i,Z.rowdim(),
                      Z.get_store()+j*Z.rowdim(),1);
             //--- Z(0:l,j) -= g*Z(0:l,i)
             //for(k=0;k<=l;k++) Z(k,j)-=g*Z(k,i);
             mat_xpeya(l+1,Z.get_store()+j*Z.rowdim(),
                       Z.get_store()+i*Z.rowdim(),-g);
         }
     }
    
     //---- store the diagonal of tridiagonal matrix in u before overwriting 
     l_u(i)=Z(i,i);

     //--- replace used householder vector by identity matrix part
     Z(i,i)=1.;
     for(j=0;j<=l;j++){ 
         Z(i,j)=0.;
         Z(j,i)=0.;
     }
 }
 return 0;
}

// *****************************************************************************
//                                tql2
// *****************************************************************************

//computes diagonal representation of the tridiagnal symmetric matrix given by
//the diagonal u and the subdiagonal v (v(0) is empty) by the Symmetric
//QR Algorithm (Golub/van Loan 2nd p423). The transformations are added to
//Z which is assumed to be the transoformation matrix bringing the
//original matrix into tridiagonal form.
//on output eigenvalues appear in u, eigenvectors are the columns of Z
//Input: u,v,Z
//Output: u,Z (v is destroyed)

int Lanczpol::tql2(Integer n,Matrix &l_u,Matrix &l_v, Matrix& Z)
{

 if (n<=1) return 0;

 const Real machep=2.e-47;
 Integer i,j,l,m;
 Integer k;
 
 //shift subdiagonal elements to position 0 to n-2
 //for (i=1;i<n;i++) l_v(i-1)=l_v(i);
 mat_xey(n-1,l_v.get_store(),l_v.get_store()+1);

 Real f=0.;
 Real subdiagzero=machep;
 l_v(n-1)=0;
 Integer l1;
 Real g;
 Real h;
 Real t2;
 Real t5;
 Real t3,t4;
 
 for(l=0;l<n;l++){  
     Integer l_iter=0;
     h=machep*(fabs(l_u(l))+fabs(l_v(l)));
     if (subdiagzero<h) subdiagzero=h;
     for(m=l;m<n;m++){
         if (fabs(l_v(m))<subdiagzero) break;
     }
     //for(m=l;m<n-1;m++){
     //    Real tst1=fabs(l_u(m))+fabs(l_u(m+1));
     //    Real tst2=tst1+fabs(l_v(m));
     //    if (tst1==tst2) break;
     //}
     if (m==l) {
         l_u(l)+=f;
         continue;
     }
     do{
         if (l_iter==30) return l+1;   //at most 30 iterations
         l_iter++;
         l1=l+1;

         //---- comput shift mu (Golub/van Loan p423)
         g=l_u(l);
         t3=(l_u(l1)-g)/(2.*l_v(l));      //t3= l_u(l+1)-l_u(l)/(2.*l_v(l))
         t4=::sqrt(t3*t3+1.);           //=::sqrt((l_u(l+1)-l_u(l))^2/4+l_v(l)^2)/|l_v(l)|
         l_u(l)=l_v(l)/(t3+d_sign(t4,t3)); 
         h=g-l_u(l);                    //=mu

         //---- subtract shift from remaining diagonal elements
         for(i=l1;i<n;i++) l_u(i)-=h;
         //mat_xpea(n-l1,l_u.get_store()+l1,-h);
         f+=h;

         
         t3=l_u(m);
         t2=1.;
         t5=0.;
         Real *vi=l_v.get_store()+m-1;
         Real *vi1=vi+1;
         Real *ui=l_u.get_store()+m-1;
         Real *ui1=ui+1;
         for(i=m-1;i>=l;i--,vi--,vi1--,ui--,ui1--){

             //---- compute Givens coefficients (Golub/van Loan p202) 
             g=t2* (*vi);               //g=t2*v(i);     
             h=t2*t3;
             if (fabs(t3)>=fabs(*vi)){    //if (fabs(t3)>=fabs(v(i))){
                 t2=(*vi)/t3;           //t2=v(i)/t3;
                 t4=::sqrt(t2*t2+1.);
                 *vi1=t5*t3*t4;         //v(i+1)=t5*t3*t4;
                 t5=t2/t4;
                 t2=1./t4;
             }
             else {
                 t2=t3/(*vi);           //t2=t3/v(i);
                 t4=::sqrt(t2*t2+1.);
                 *vi1=t5*(*vi)*t4;      //v(i+1)=t5*v(i)*t4;
                 t5=1./t4;
                 t2=t2*t5;
             }
             t3=t2*(*ui)-t5*g;          //t3=t2*u(i)-t5*g;
             *ui1=h+t5*(t2*g+t5*(*ui)); //l_u(i+1)=h+t5*(t2*g+t5*l_u(i));

             //---- change vectors
             Real *zi=Z.get_store()+i*Z.rowdim();
             Real *zi1=zi+Z.rowdim();
             for(k=n;--k>=0;){
                 //g=Z(k,i); h=Z(k,i+1);
                 g=*zi;
                 h=*zi1;
                 (*zi1++) = t5*g+t2*h;      //Z(k,i+1)=t5*Z(k,i)+t2*h;
                 (*zi++)  = t2*g-t5*h;      //Z(k,i)=t2*Z(k,i)-t5*h;
             }
         }
         l_v(l)=t5*t3;
         l_u(l)=t2*t3;
     } while(fabs(l_v(l))>subdiagzero); 
     l_u(l)+=f;
 }

 //--- sort eigenvalues nonincreasingly and rearrange eigenvectors
 //    (quite stupid, should probably be replaced) ?does it work at all? 
 for(i=0;i<n-1;i++){
     k=i;
     t3=l_u(i);
     //--- find remaining maximal eigenvalue 
     for(j=i+1;j<n;j++){
         if(l_u(j)<=t3) continue;
         k=j;
         t3=l_u(j);
     }
     if (k==i) continue;   //i is already maximal, skip rest

     //--- swap maximal eigenvalue and eigenvector k with i
     l_u(k)=l_u(i);            
     l_u(i)=t3;
     //for(j=0;j<n;j++){
     //    t3=Z(j,i); Z(j,i)=Z(j,k); Z(j,k)=t3;
     //}
     mat_swap(n,Z.get_store()+i*Z.rowdim(),Z.get_store()+k*Z.rowdim());
 }
 return 0;
}


// *****************************************************************************
//                                save
// *****************************************************************************

std::ostream& Lanczpol::save(std::ostream& out) const
{
 out<<ierr<<"\n"<<maxop<<"\n"<<maxiter<<"\n";
 out<<maxguessiter<<"\n"<<guessmult<<"\n"<<choicencheb<<"\n"; 
 out<<nchebit<<"\n"<<choicenbmult<<"\n"<<nblockmult<<"\n";
 out<<nlanczvecs<<"\n"<<neigfound<<"\n"<<blocksz<<"\n";
 out<<iter<<"\n"<<nmult<<"\n"<<errc<<"\n"<<eps<<"\n";
 out<<mcheps<<"\n"<<maxval<<"\n"<<minval<<"\n"<<polval<<"\n"; 
 out<<X<<C<<d<<e<<u<<v<<w<<minvec;      
 out<<stop_above<<"\n"<<upper_bound<<"\n"<<print_level<<"\n"<<mymaxj<<"\n";
 out<<time_mult<<"\n"<<time_mult_sum<<"\n"<<time_iter<<"\n"<<time_sum<<"\n";
 return out;
}

// *****************************************************************************
//                                restore
// *****************************************************************************

std::istream& Lanczpol::restore(std::istream& in)
{
 in>>ierr>>maxop>>maxiter;
 in>>maxguessiter>>guessmult>>choicencheb; 
 in>>nchebit>>choicenbmult>>nblockmult;
 in>>nlanczvecs>>neigfound>>blocksz;
 in>>iter>>nmult>>errc>>eps;
 in>>mcheps>>maxval>>minval>>polval; 
 in>>X>>C>>d>>e>>u>>v>>w>>minvec;      
 in>>stop_above>>upper_bound>>print_level>>mymaxj;
 in>>time_mult>>time_mult_sum>>time_iter>>time_sum;
 return in;
}

}

