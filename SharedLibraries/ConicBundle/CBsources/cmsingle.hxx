/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/cmsingle.hxx

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



#ifndef CONICBUNDLE_CMSINGLE_HXX
#define CONICBUNDLE_CMSINGLE_HXX

//defines a base class for constraint matrices.
//the idea is to help exploiting special structures for
//the computation of tr(AiXAjZi) without having to know
//the special structure of both, Ai and Aj. 

#include "cmgramde.hxx"
#include "cmrankdd.hxx"
#include "sparssym.hxx"

namespace ConicBundle {

  class CMsingleton: public Coeffmat
  {
  private:
    CH_Matrix_Classes::Integer nr;
    CH_Matrix_Classes::Integer ii;
    CH_Matrix_Classes::Integer jj;
    CH_Matrix_Classes::Real    val;
  public:
    CMsingleton(CH_Matrix_Classes::Integer innr, CH_Matrix_Classes::Integer ini, CH_Matrix_Classes::Integer inj, CH_Matrix_Classes::Real inval,CH_Matrix_Classes::Integer k=0)
    {nr=innr;ii=ini;jj=inj;val=inval;CM_type=CM_singleton;userkey=k;}
    virtual ~CMsingleton(){}

    //virtual CM_type get_type() const {return CM_type;} //no need to overload
    //virtual CH_Matrix_Classes::Integer get_userkey() const {return userkey;}
    //virtual void set_userkey(CH_Matrix_Classes::Integer k) const {userkey=k;}
    
    virtual Coeffmat* clone() const 
    {return new CMsingleton(nr,ii,jj,val,userkey);}
    //makes an explicit copy of itself and returns a pointer to it 

    virtual CH_Matrix_Classes::Integer dim() const {return nr;}
    //returns the order of the represented symmetric matrix
    
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j) const 
    { if (((i==ii)&&(j==jj))||((i==jj)&&(j==ii))) return val; return 0.; }

    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const
    {S.init(nr,0.); S(ii,jj)=val;}
    //return dense symmetric constraint matrix

    virtual CH_Matrix_Classes::Real norm(void) const
    {if (ii==jj) return fabs(val); return sqrt(2.)*fabs(val);}
    //compute Frobenius norm of matrix

    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const
    {
     if (ii==jj){ return new CMgramdense(sqrt(fabs(val))*(P.row(ii)).transpose(),val>0); }
     return new CMlowrankdd((P.row(ii)).transpose()*val,(P.row(jj)).transpose());
    }
    //delivers a new object on the heap corresponding
    //to the matrix P^TAP, the caller is responsible for deleting the object
    
    virtual void multiply(CH_Matrix_Classes::Real d) { val*=d; }
    //multiply constraint permanentely by d
    //this is to allow scaling or sign changes in the constraints

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const
    {if (ii==jj) return val*S(ii,ii); return 2.*val*S(ii,jj);}
    //=ip(*this,S)=trace(*this*S) trace inner product
    
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const
    {
     if (ii==jj)
        return val*CH_Matrix_Classes::mat_ip(P.coldim(),P.get_store()+ii,P.rowdim(),P.get_store()+ii,P.rowdim());
     return 2.*val*CH_Matrix_Classes::mat_ip(P.coldim(),P.get_store()+ii,P.rowdim(),P.get_store()+jj,P.rowdim());  
    }
    //=ip(*this,PP^T)=trace P^T(*this)P

    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    {
     if (ii==jj)
        return val*CH_Matrix_Classes::mat_ip(P.coldim(),P.get_store()+ii+start_row,P.rowdim(),P.get_store()+ii+start_row,P.rowdim());
     return 2.*val*CH_Matrix_Classes::mat_ip(P.coldim(),P.get_store()+ii+start_row,P.rowdim(),P.get_store()+jj+start_row,P.rowdim());  
    }

    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S,CH_Matrix_Classes::Real d=1.) const {S(ii,jj)+=d*val;}
    //S+=d*(*this);

    virtual void addprodto(CH_Matrix_Classes::Matrix& B,const CH_Matrix_Classes::Matrix&C ,CH_Matrix_Classes::Real d=1.) const
    {
     CH_Matrix_Classes::mat_xpeya(C.coldim(),B.get_store()+ii,B.rowdim(),C.get_store()+jj,C.rowdim(),d*val);
     if (ii==jj) return;
     CH_Matrix_Classes::mat_xpeya(C.coldim(),B.get_store()+jj,B.rowdim(),C.get_store()+ii,C.rowdim(),d*val);
    }
    //B+=d*(*this)*C

    virtual void addprodto(CH_Matrix_Classes::Matrix& B,const CH_Matrix_Classes::Sparsemat&C ,CH_Matrix_Classes::Real d=1.) const
    {
      if (C.coldim()==1){
	B(ii)+=d*val*C(jj);
	if (ii!=jj) B(jj)+=d*val*C(ii);
        return;
      }
      CH_Matrix_Classes::Sparsemat tmp(C.row(jj));
      for(CH_Matrix_Classes::Integer i=0;i<tmp.nonzeros();i++){
	CH_Matrix_Classes::Integer indi, indj; CH_Matrix_Classes::Real v;
	tmp.get_edge(i,indi,indj,v);
	B(ii,indj)+=d*val*v;
      }
      if (ii==jj) return;
      tmp.init(C.row(ii));
      {for(CH_Matrix_Classes::Integer i=0;i<tmp.nonzeros();i++){
	CH_Matrix_Classes::Integer indi, indj; CH_Matrix_Classes::Real v;
	tmp.get_edge(i,indi,indj,v);
	B(jj,indj)+=d*val*v;
      }}
    }
    //B+=d*(*this)*C

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& Q,CH_Matrix_Classes::Matrix& R) const
    {
      if (ii==jj) CH_Matrix_Classes::genmult(P.row(ii),Q.row(ii),R,val,0.,1,0);
      else {
	CH_Matrix_Classes::genmult(P.row(ii),Q.row(jj),R,val,0.,1,0);
	CH_Matrix_Classes::genmult(P.row(jj),Q.row(ii),R,val,1.,1,0);
      }
    }
        
    virtual CH_Matrix_Classes::Integer prodvec_flops() const 
    { return (ii==jj)?2:4; }
    //return estimate of number of flops to compute addprodto for a vector

    virtual int dense() const
    {return 0;}
    //returns 1 if its structure as bad as its dense symmetric representation,
    //otherwise 1
    
    virtual int sparse() const
    { return 1;}
    //returns 0 if not sparse otherwise 1
        
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& I,CH_Matrix_Classes::Indexmatrix& J,CH_Matrix_Classes::Matrix& v,CH_Matrix_Classes::Real d=1.)const
    { I.init(1,1,ii); J.init(1,1,jj); v.init(1,1,val*d); return 1;}
    //returns 0 if not sparse. If it is spars it returns 1 and
    //the nonzero structure in I,J and v, where v is multiplied by d.
    //Only the upper triangle (including diagonal) is delivered
    
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& S) const
    {return S.check_support(ii,jj);}
    //returns 0 if the support of the costraint matrix is not contained in the
    //support of the sparse symmetric matrix A, 1 if it is contained.

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& S) const
    {if (ii==jj) return val*S(ii,jj); return 2.*val*S(ii,jj);}
    //returns the inner product of the constraint matrix with A
    
    virtual void project(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P) const // S=P^t*A*P
    {
     S.newsize(P.coldim());chk_set_init(S,1);
     const CH_Matrix_Classes::Real *ap=P.get_store()+ii-nr;
     CH_Matrix_Classes::Real *sp=S.get_store();
     if (ii==jj){
         for(CH_Matrix_Classes::Integer i=P.coldim();--i>=0;){
             const CH_Matrix_Classes::Real *aap=ap; CH_Matrix_Classes::Real a=*(aap=(ap+=nr));
             *sp++=a*a*val; a*=val;
             for(CH_Matrix_Classes::Integer j=i;--j>=0;) *sp++=a*(*(aap+=nr));
         }
         return;
     }
     const CH_Matrix_Classes::Real *bp=P.get_store()+jj-nr;
     for(CH_Matrix_Classes::Integer i=P.coldim();--i>=0;){
         const CH_Matrix_Classes::Real *aap; CH_Matrix_Classes::Real a=*(aap=(ap+=nr))*val;
         const CH_Matrix_Classes::Real *bbp; CH_Matrix_Classes::Real b=*(bbp=(bp+=nr));
         *sp++=2.*a*b; b*=val;
         for(CH_Matrix_Classes::Integer j=i;--j>=0;) *sp++=a*(*(bbp+=nr))+b*(*(aap+=nr));
     }
    }

    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    {
     CH_Matrix_Classes::Integer pnr=P.rowdim();
     const CH_Matrix_Classes::Real *ap=P.get_store()+start_row+ii-pnr;
     CH_Matrix_Classes::Real *sp=S.get_store();
     if (ii==jj){
         for(CH_Matrix_Classes::Integer i=P.coldim();--i>=0;){
             const CH_Matrix_Classes::Real *aap=ap; CH_Matrix_Classes::Real a=*(aap=(ap+=pnr));
             (*sp++)+=a*a*val; a*=val;
             for(CH_Matrix_Classes::Integer j=i;--j>=0;) (*sp++)+=a*(*(aap+=pnr));
         }
         return;
     }
     const CH_Matrix_Classes::Real *bp=P.get_store()+start_row+jj-pnr;
     for(CH_Matrix_Classes::Integer i=P.coldim();--i>=0;){
         const CH_Matrix_Classes::Real *aap; CH_Matrix_Classes::Real a=*(aap=(ap+=pnr))*val;
         const CH_Matrix_Classes::Real *bbp; CH_Matrix_Classes::Real b=*(bbp=(bp+=pnr));
         (*sp++)+=2.*a*b; b*=val;
         for(CH_Matrix_Classes::Integer j=i;--j>=0;) (*sp++)+=a*(*(bbp+=pnr))+b*(*(aap+=pnr));
     }
    }

    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const
    { 
      CH_Matrix_Classes::Sparsesym S(nr,1,&ii,&jj,&val);
      return CH_Matrix_Classes::genmult(S,B,C,alpha,beta,btrans);
    }
    // C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned

    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const
    { 
      CH_Matrix_Classes::Sparsesym S(nr,1,&ii,&jj,&val);
      CH_Matrix_Classes::Matrix D; 
      return CH_Matrix_Classes::genmult(B,S,C,alpha,beta,btrans);
    }
    // C= alpha*B^(T if btrans)*(*this) + beta*C 

    virtual int equal(const Coeffmat* p,double tol=1e-6) const
    {
      const CMsingleton *pp=dynamic_cast<const CMsingleton *>(p);
      if (pp==0) 
	return 0;
      return ((nr==pp->nr)&&(ii==pp->ii)&&(jj==pp->jj)
	      &&(fabs(val-pp->val)<tol));
    }
    //return true, if p is the same derived class and 
    //entries differ by less than tol

    virtual std::ostream& display(std::ostream& o) const 
    {o<<"CMsingleton\n"; o<<nr<<" "<<ii<<" "<<jj<<" "<<val<<"\n"; return o;}
    //display constraint information
    
    virtual std::ostream& out(std::ostream& o) const
    {return o<<"SINGLETON\n"<<nr<<" "<<ii<<" "<<jj<<" "<<val<<"\n";}
    //put entire contents onto outstream with the class type in the beginning so
    //that the derived class can be recognized.

    virtual std::istream& in(std::istream& is)
    {
     if (!(is>>nr>>ii>>jj>>val)) { 
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: CMsingleton::in(): reading from input failed";
         is.clear(std::ios::failbit);
         return is;
     }
     if (nr<0) {
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: CMsingleton::in(): dimension of matrix must positive";
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"but is "<<nr<<std::endl;
         is.clear(std::ios::failbit);
         return is;
     }
     if ((ii<0)||(ii>nr)){
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: CMsingleton::in(): row index outside range, ";
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<0<<"<="<<ii<<"<"<<nr<<std::endl;
         is.clear(std::ios::failbit);
         return is;
     }
     if ((jj<0)||(jj>nr)){
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<"*** ERROR: CMsingleton::in(): column index outside range, ";
         if (CH_Matrix_Classes::materrout) (*CH_Matrix_Classes::materrout)<<0<<"<="<<jj<<"<"<<nr<<std::endl;
         is.clear(std::ios::failbit);
         return is;
     }
     return is;
    }
    //counterpart to out, does not read the class type, though.
    //This is assumed to have been read in order to generate the correct class

    CMsingleton(std::istream& is,CH_Matrix_Classes::Integer k=0)
    {CM_type=CM_singleton;userkey=k;in(is);}

    //--- specific routines
    int get_ijval(CH_Matrix_Classes::Integer& i,CH_Matrix_Classes::Integer &j, CH_Matrix_Classes::Real& v) const
    {i=ii; j=jj; v=val; return 0;}

  };
}

#endif

