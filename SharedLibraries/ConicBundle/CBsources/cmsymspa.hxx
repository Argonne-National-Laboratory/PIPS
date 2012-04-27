/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/cmsymspa.hxx

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


#ifndef CONICBUNDLE_CMSYMSPA_HXX
#define CONICBUNDLE_CMSYMSPA_HXX

//defines a base class for constraint matrices.
//the idea is to help exploiting special structures for
//the computation of tr(AiXAjZi) without having to know
//the special structure of both, Ai and Aj. 

#include "cmsymden.hxx"
#include "sparssym.hxx"

namespace ConicBundle {

  class CMsymsparse: public Coeffmat
  {
  private:
    CH_Matrix_Classes::Sparsesym A;
    bool use_sparsemult;
  public:
    CMsymsparse(const CH_Matrix_Classes::Sparsesym& Ain,CH_Matrix_Classes::Integer k=0)
    {A=Ain;CM_type=CM_symsparse;userkey=k; use_sparsemult=false;
    if (A.get_suppcol().rowdim()<A.rowdim()/2) use_sparsemult=true; }
    virtual ~CMsymsparse(){}

    //virtual CM_type get_type() const {return CM_type;} //no need to overload
    //virtual CH_Matrix_Classes::Integer get_userkey() const {return userkey;}
    //virtual void set_userkey(CH_Matrix_Classes::Integer k) const {userkey=k;}

    virtual Coeffmat* clone() const 
    {return new CMsymsparse(A,userkey);}
    //makes an explicit copy of itself and returns a pointer to it 
    
    virtual CH_Matrix_Classes::Integer dim() const {return A.rowdim();}
    //returns the order of the represented symmetric matrix
    
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j) const {return A(i,j);}
    //returns the value of the matrix element (i,j)
    
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const {S=A;}
    //return dense symmetric constraint matrix

    virtual CH_Matrix_Classes::Real norm(void) const {return CH_Matrix_Classes::norm2(A);}
    //compute Frobenius norm of matrix

    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const
    {
     CH_Matrix_Classes::Symmatrix S; 
     if (use_sparsemult) S.xetriu_yza(A.sparsemult(P),P);
     else S.xetriu_yza(A*P,P);
     return new CMsymdense(S,userkey);
    }
    //delivers a new object on the heap corresponding
    //to the matrix P^TAP, the caller is responsible for deleting the object
    
    virtual void multiply(CH_Matrix_Classes::Real d) {A*=d;}
    //multiply constraint permanentely by d
    //this is to allow scaling or sign changes in the constraints

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const {return CH_Matrix_Classes::ip(S,A);}
    //=ip(*this,S)=trace(*this*S) trace inner product
    
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const
    {if (use_sparsemult) return CH_Matrix_Classes::ip(A.sparsemult(P),P); return CH_Matrix_Classes::ip(A*P,P);}
    //=ip(*this,PP^T)=trace P^T(*this)P

    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    {
      const CH_Matrix_Classes::Indexmatrix &colinfo=A.get_colinfo();
      const CH_Matrix_Classes::Indexmatrix &suppcol=A.get_suppcol();
      CH_Matrix_Classes::Matrix B(suppcol.dim(),P.coldim(),0.);
      CH_Matrix_Classes::Real *bp=B.get_store();
      const CH_Matrix_Classes::Integer brd=B.rowdim();
      const CH_Matrix_Classes::Integer pcd=P.coldim();
      const CH_Matrix_Classes::Integer prd=P.rowdim();
      const CH_Matrix_Classes::Real *pp=P.get_store()+start_row;
      const CH_Matrix_Classes::Real *vp=(A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer *ip=(A.get_colindex()).get_store();
      const CH_Matrix_Classes::Integer *sip=(A.get_suppind()).get_store();
      for(CH_Matrix_Classes::Integer j=0;j<colinfo.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo(j,0);
	if (indj<0){ //"column" is the diagonal 
	  for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	    CH_Matrix_Classes::mat_xpeya(pcd,bp+(*sip++),brd,pp+(*ip++),prd,(*vp++));
	  }
	  continue;
	}
	//offdiagonal
	CH_Matrix_Classes::Integer sindj=colinfo(j,3);
	for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+sindj,brd,pp+(*ip++)+indj,prd,(*vp));
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+(*sip++)+sindj,brd,pp+indj,prd,(*vp++));
	}
      } 
      //make S
      CH_Matrix_Classes::Real sumval=0.;
      for(CH_Matrix_Classes::Integer i=0;i<brd;i++){
	const CH_Matrix_Classes::Real* pp2=pp+suppcol(i);
	const CH_Matrix_Classes::Real* bp=B.get_store()+i;
	sumval+=CH_Matrix_Classes::mat_ip(pcd,bp,brd,pp2,prd);
      }
      return sumval;
    }

    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S,CH_Matrix_Classes::Real d=1.) const {S.xpeya(A,d);}
    //S+=d*(*this);

    virtual void addprodto(CH_Matrix_Classes::Matrix& B,const CH_Matrix_Classes::Matrix&C ,CH_Matrix_Classes::Real d=1.) const
    {if (use_sparsemult) B.xpeya(A.sparsemult(C),d); else B.xpeya(A*C,d); }
    //B+=d*(*this)*C

    virtual void addprodto(CH_Matrix_Classes::Matrix& B,const CH_Matrix_Classes::Sparsemat&C ,CH_Matrix_Classes::Real d=1.) const
    { CH_Matrix_Classes::genmult(A,C,B,d,1.); }
    //B+=d*(*this)*C

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& Q,CH_Matrix_Classes::Matrix& R) const
    {
      if (A.get_suppcol().dim()>A.rowdim()/2){
	if (P.coldim()<Q.coldim()){
	  CH_Matrix_Classes::Matrix tmp1; 
          CH_Matrix_Classes::genmult(A,P,tmp1,1.,0.,0);
	  CH_Matrix_Classes::genmult(tmp1,Q,R,1.,0.,1,0);
	}
	else {
	  CH_Matrix_Classes::Matrix tmp1; 
	  CH_Matrix_Classes::genmult(A,Q,tmp1,1.,0.,0);
	  CH_Matrix_Classes::genmult(P,tmp1,R,1.,0.,1,0);
	}
      }
      else {
	if (P.coldim()<Q.coldim()) CH_Matrix_Classes::genmult(A.sparsemult(P),Q,R,1.,0.,1,0);
        else CH_Matrix_Classes::genmult(P,A.sparsemult(Q),R,1.,0.,1,0);
      }
    }

    virtual CH_Matrix_Classes::Integer prodvec_flops() const 
    { return 2*A.nonzeros(); }
    //return estimate of number of flops to compute addprodto for a vector

    virtual int dense() const
    {return 0;}
    //returns 1 if its structure as bad as its dense symmetric representation,
    //otherwise 1
    
    virtual int sparse() const
    { return 1;}
    //returns 0 if not sparse otherwise 1
        
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& I,CH_Matrix_Classes::Indexmatrix& J,CH_Matrix_Classes::Matrix& val,CH_Matrix_Classes::Real d=1.)const
    { A.get_edge_rep(I,J,val); val*=d; return 1;}
    //returns 0 if not sparse. If it is spars it returns 1 and
    //the nonzero structure in I,J and val, where val is multiplied by d.
    //Only the upper triangle (including diagonal) is delivered
    
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& S) const
    {return S.contains_support(A);}
    //returns 0 if the support of the costraint matrix is not contained in the
    //support of the sparse symmetric matrix A, 1 if it is contained.

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& S) const
    {return CH_Matrix_Classes::ip(A,S);}
    //returns the inner product of the constraint matrix with A
    
    virtual void project(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P) const
    {
     if (use_sparsemult) S.xetriu_yza(A.sparsemult(P),P);
     else S.xetriu_yza(A*P,P);
    }
    // S=P^t*A*P

    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    {
      const CH_Matrix_Classes::Indexmatrix &colinfo=A.get_colinfo();
      const CH_Matrix_Classes::Indexmatrix &suppcol=A.get_suppcol();
      CH_Matrix_Classes::Matrix B(suppcol.dim(),P.coldim(),0.);
      CH_Matrix_Classes::Real *bp=B.get_store();
      const CH_Matrix_Classes::Integer brd=B.rowdim();
      const CH_Matrix_Classes::Integer pcd=P.coldim();
      const CH_Matrix_Classes::Integer prd=P.rowdim();
      const CH_Matrix_Classes::Real *pp=P.get_store()+start_row;
      const CH_Matrix_Classes::Real *vp=(A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer *ip=(A.get_colindex()).get_store();
      const CH_Matrix_Classes::Integer *sip=(A.get_suppind()).get_store();
      for(CH_Matrix_Classes::Integer j=0;j<colinfo.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo(j,0);
	if (indj<0){ //"column" is the diagonal 
	  for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	    CH_Matrix_Classes::mat_xpeya(pcd,bp+(*sip++),brd,pp+(*ip++),prd,(*vp++));
	  }
	  continue;
	}
	//offdiagonal
	CH_Matrix_Classes::Integer sindj=colinfo(j,3);
	for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+sindj,brd,pp+(*ip++)+indj,prd,(*vp));
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+(*sip++)+sindj,brd,pp+indj,prd,(*vp++));
	}
      } 
      //add to S
      for(CH_Matrix_Classes::Integer i=0;i<brd;i++){
        CH_Matrix_Classes::Real* sp=S.get_store();
	const CH_Matrix_Classes::Real* pp2=pp+suppcol(i);
	const CH_Matrix_Classes::Real* bp=B.get_store()+i;
	for(CH_Matrix_Classes::Integer j=pcd;j>0;--j){
	  CH_Matrix_Classes::mat_xpeya(j,sp,1,bp,brd,(*pp2));
	  sp+=j;
	  bp+=brd;
	  pp2+=prd;
	}
      }
    }
    // S=P^t*A*P

    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const
    { 
      return CH_Matrix_Classes::genmult(A,B,C,alpha,beta,btrans);
    }
    // C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned

    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const
    { 
      CH_Matrix_Classes::Matrix D; 
      return CH_Matrix_Classes::genmult(B,A,C,alpha,beta,btrans);
    }
    // C= alpha*B^(T if btrans)*(*this) + beta*C 

    virtual int equal(const Coeffmat* p,double tol=1e-6) const
    {
      const CMsymsparse *pp=dynamic_cast<const CMsymsparse *>(p);
      if (pp==0) 
	return 0;
      return CH_Matrix_Classes::equal(A,pp->A,tol);
    }
    //return true, if p is the same derived class and 
    //entries differ by less than tol

    virtual std::ostream& display(std::ostream& o) const 
    {o<<"CMsymsparse\n";A.display(o); return o;}
    //display constraint information
    
    virtual std::ostream& out(std::ostream& o) const
    {return o<<"SYMMETRIC_SPARSE\n"<<A;}
    //put entire contents onto outstream with the class type in the beginning so
    //that the derived class can be recognized.

    virtual std::istream& in(std::istream& i)
    {return i>>A;}
    //counterpart to out, does not read the class type, though.
    //This is assumed to have been read in order to generate the correct class

    CMsymsparse(std::istream& is,CH_Matrix_Classes::Integer k=0)
    {CM_type=CM_symsparse;userkey=k;in(is);use_sparsemult=false;
     if (A.get_suppcol().rowdim()<A.rowdim()/2) use_sparsemult=true;
    }

    //--- specific routines
    const CH_Matrix_Classes::Sparsesym& get_A() const{return A;}

  };

}

#endif

