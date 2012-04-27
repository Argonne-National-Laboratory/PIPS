/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/cmrankss.hxx

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


#ifndef CONICBUNDLE_CMRANKSS_HXX
#define CONICBUNDLE_CMRANKSS_HXX

//defines a base class for constraint matrices.
//the idea is to help exploiting special structures for
//the computation of tr(AiXAjZi) without having to know
//the special structure of both, Ai and Aj. 

#include "coeffmat.hxx"
#include "cmrankdd.hxx"

namespace ConicBundle {

  class CMlowrankss: public Coeffmat
  {
  private:
    CH_Matrix_Classes::Sparsemat A,B;  //CH_Matrix_Classes::Symmatrix=A*B^T+B*A^T
  public:
    CMlowrankss(const CH_Matrix_Classes::Sparsemat& Ain,const CH_Matrix_Classes::Sparsemat& Bin,CH_Matrix_Classes::Integer k=0)
    {A=Ain;B=Bin;CM_type=CM_lowrankss;userkey=k;}
    virtual ~CMlowrankss(){}

    //virtual CM_type get_type() const {return CM_type;} //no need to overload
    //virtual CH_Matrix_Classes::Integer get_userkey() const {return userkey;}
    //virtual void set_userkey(CH_Matrix_Classes::Integer k) const {userkey=k;}
    
    virtual Coeffmat* clone() const 
    {return new CMlowrankss(A,B,userkey);}
    //makes an explicit copy of itself and returns a pointer to it 

    virtual CH_Matrix_Classes::Integer dim() const { return A.rowdim(); }
    //returns the order of the represented symmetric matrix
    
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j) const 
    { return CH_Matrix_Classes::ip(A.row(i),B.row(j))+ CH_Matrix_Classes::ip(B.row(i),A.row(j));}

    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const
    { CH_Matrix_Classes::rank2add(A,CH_Matrix_Classes::Matrix(B),S,2.); }
    //return dense symmetric constraint matrix

    virtual CH_Matrix_Classes::Real norm(void) const
    { CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::Matrix(B),C,1.,0.,1); CH_Matrix_Classes::genmult(C,C,D);
      CH_Matrix_Classes::Real d=2.*CH_Matrix_Classes::trace(D); CH_Matrix_Classes::genmult(A,A,C,1.,0.,1); CH_Matrix_Classes::genmult(B,B,D,1.,0.,1);
      return sqrt(2.*CH_Matrix_Classes::ip(C,D)+d);}
    //compute Frobenius norm of matrix

    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const
    {CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(P,A,C,1.,0.,1); CH_Matrix_Classes::genmult(P,B,D,1.,0.,1);
     return new CMlowrankdd(C,D,userkey); }
    //delivers a new object on the heap corresponding
    //to the matrix P^TAP, the caller is responsible for deleting the object
    
    virtual void multiply(CH_Matrix_Classes::Real d)
    { A*=d; }
    //multiply constraint permanentely by d
    //this is to allow scaling or sign changes in the constraints

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const
    { CH_Matrix_Classes::Matrix C; return 2.*CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S,A,C),B); }
    //=ip(*this,S)=trace(*this*S) trace inner product
    
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const
    { CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(P,A,C,1.,0.,1); CH_Matrix_Classes::genmult(P,B,D,1.,0.,1);
      return 2.*CH_Matrix_Classes::ip(C,D); }
    //=ip(*this,PP^T)=trace P^T(*this)P

    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    { 
      CH_Matrix_Classes::Matrix C(A.coldim(),P.coldim(),0.);
      CH_Matrix_Classes::Real *bp=C.get_store();
      const CH_Matrix_Classes::Integer brd=C.rowdim();
      const CH_Matrix_Classes::Integer pcd=P.coldim();
      const CH_Matrix_Classes::Integer prd=P.rowdim();
      const CH_Matrix_Classes::Real *pp=P.get_store()+start_row;
      const CH_Matrix_Classes::Indexmatrix &colinfo=A.get_colinfo();
      const CH_Matrix_Classes::Real *vp=(A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer *ip=(A.get_colindex()).get_store();
      for(CH_Matrix_Classes::Integer j=0;j<colinfo.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }
      CH_Matrix_Classes::Matrix D(B.coldim(),P.coldim(),0.);
      bp=D.get_store();
      const CH_Matrix_Classes::Indexmatrix &colinfo2=B.get_colinfo();
      vp=(B.get_colval()).get_store();
      ip=(B.get_colindex()).get_store();
      {for(CH_Matrix_Classes::Integer j=0;j<colinfo2.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo2(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo2(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }}
      return 2.*CH_Matrix_Classes::ip(C,D); 
    }
    //=ip(*this,PP^T)=trace P^T(*this)P

    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S,CH_Matrix_Classes::Real d=1.) const
    { CH_Matrix_Classes::rank2add(A,CH_Matrix_Classes::Matrix(B),S,2.*d,1.);}
    //S+=d*(*this);

    virtual void addprodto(CH_Matrix_Classes::Matrix& D,const CH_Matrix_Classes::Matrix&C ,CH_Matrix_Classes::Real d=1.) const
    {CH_Matrix_Classes::Matrix E; CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(B,C,E,1.,0.,1),D,d,1.);
     CH_Matrix_Classes::genmult(B,CH_Matrix_Classes::genmult(A,C,E,1.,0.,1),D,d,1.);}
    //B+=d*(*this)*C

    virtual void addprodto(CH_Matrix_Classes::Matrix& D,const CH_Matrix_Classes::Sparsemat&C ,CH_Matrix_Classes::Real d=1.) const
    {CH_Matrix_Classes::Matrix E; CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(B,C,E,1.,0.,1),D,d,1.);
     CH_Matrix_Classes::genmult(B,CH_Matrix_Classes::genmult(A,C,E,1.,0.,1),D,d,1.);}
    //B+=d*(*this)*C

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& Q,CH_Matrix_Classes::Matrix& R) const
    {
      CH_Matrix_Classes::Matrix tmp1; CH_Matrix_Classes::genmult(P,A,tmp1,1.,0.,1,0);
      CH_Matrix_Classes::Matrix tmp2; CH_Matrix_Classes::genmult(B,Q,tmp2,1.,0.,1,0);
      CH_Matrix_Classes::genmult(tmp1,tmp2,R,1.,0.,0,0);
      CH_Matrix_Classes::genmult(P,B,tmp1,1.,0.,1,0);
      CH_Matrix_Classes::genmult(A,Q,tmp2,1.,0.,1,0);
      CH_Matrix_Classes::genmult(tmp1,tmp2,R,1.,1.,0,0);
    }

    virtual CH_Matrix_Classes::Integer prodvec_flops() const 
    { return 4*A.nonzeros()+4*B.nonzeros(); }
    //return estimate of number of flops to compute addprodto for a vector

    virtual int dense() const
    {return 0;}
    //returns 1 if its structure as bad as its dense symmetric representation,
    //otherwise 1
    
    virtual int sparse() const
    { return 0;}
    //returns 0 if not sparse otherwise 1
    
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& /* I */,
		       CH_Matrix_Classes::Indexmatrix& /* J */,
		       CH_Matrix_Classes::Matrix& /* val */,
		       CH_Matrix_Classes::Real /* d=1. */)const
    {return 0;}
    //returns 0 if not sparse. If it is spars it returns 1 and
    //the nonzero structure in I,J and val, where val is multiplied by d.
    //Only the upper triangle (including diagonal) is delivered
    
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& /* S */) const
    {return 0;}
    //returns 0 if the support of the costraint matrix is not contained in the
    //support of the sparse symmetric matrix A, 1 if it is contained.

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& S) const
    {return 2.*CH_Matrix_Classes::ip(A,S*B);}
    //returns the inner product of the constraint matrix with A
    
    virtual void project(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P) const
    {CH_Matrix_Classes::Matrix C,D; CH_Matrix_Classes::genmult(P,A,C,1.,0.,1); CH_Matrix_Classes::genmult(P,B,D,1.,0.,1);
     CH_Matrix_Classes::rank2add(C,D,S,2.);}
    // S=P^t*(*this)*P

    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    {
      CH_Matrix_Classes::Matrix C(A.coldim(),P.coldim(),0.);
      CH_Matrix_Classes::Real *bp=C.get_store();
      const CH_Matrix_Classes::Integer brd=C.rowdim();
      const CH_Matrix_Classes::Integer pcd=P.coldim();
      const CH_Matrix_Classes::Integer prd=P.rowdim();
      const CH_Matrix_Classes::Real *pp=P.get_store()+start_row;
      const CH_Matrix_Classes::Indexmatrix &colinfo=A.get_colinfo();
      const CH_Matrix_Classes::Real *vp=(A.get_colval()).get_store();
      const CH_Matrix_Classes::Integer *ip=(A.get_colindex()).get_store();
      for(CH_Matrix_Classes::Integer j=0;j<colinfo.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }
      CH_Matrix_Classes::Matrix D(B.coldim(),P.coldim(),0.);
      bp=D.get_store();
      const CH_Matrix_Classes::Indexmatrix &colinfo2=B.get_colinfo();
      vp=(B.get_colval()).get_store();
      ip=(B.get_colindex()).get_store();
      {for(CH_Matrix_Classes::Integer j=0;j<colinfo2.rowdim();j++){
	CH_Matrix_Classes::Integer indj=colinfo2(j,0);
	for(CH_Matrix_Classes::Integer i=colinfo2(j,1);--i>=0;){
	  CH_Matrix_Classes::mat_xpeya(pcd,bp+indj,brd,pp+(*ip++),prd,(*vp++));
	}
      }}
      CH_Matrix_Classes::rank2add(C,D,S,2.,1.,1); 
    }

    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& D,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int dtrans=0) const
    { 
      CH_Matrix_Classes::Matrix E; 
      CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(B,D,E,1.,0.,1,dtrans),C,alpha,beta);
      return CH_Matrix_Classes::genmult(B,CH_Matrix_Classes::genmult(A,D,E,1.,0.,1,dtrans),C,alpha,1.);
    }
    // C= alpha*(*this)*D^(T if btrans) + beta*C, C is also returned

    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& D,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int dtrans=0) const
    { 
      CH_Matrix_Classes::Matrix E; 
      CH_Matrix_Classes::genmult(CH_Matrix_Classes::genmult(D,A,E,1.,0.,dtrans),B,C,alpha,beta,0,1);
      return CH_Matrix_Classes::genmult(CH_Matrix_Classes::genmult(D,B,E,1.,0.,dtrans),A,C,alpha,1.,0,1);
    }
    // C= alpha*D^(T if btrans)*(*this) + beta*C 

    virtual int equal(const Coeffmat* p,double tol=1e-6) const
    {
      const CMlowrankss *pp=dynamic_cast<const CMlowrankss *>(p);
      if (pp==0) 
	return 0;
      if (CH_Matrix_Classes::equal(A,pp->A,tol) && CH_Matrix_Classes::equal(B,pp->B,tol))
	return 1;
      if (CH_Matrix_Classes::equal(A,pp->B,tol) && CH_Matrix_Classes::equal(B,pp->A,tol))
	return 1;
      return 0;
    }
    //return true, if p is the same derived class and 
    //entries differ by less than tol

    virtual std::ostream& display(std::ostream& o) const 
    {o<<"CMlowrankss\n";A.display(o);B.display(o);return o;}
    //display constraint information
    
    virtual std::ostream& out(std::ostream& o) const
    {return o<<"LOWRANK_SPARSE_SPARSE\n"<<A<<B;}
    //put entire contents onto outstream with the class type in the beginning so
    //that the derived class can be recognized.

    virtual std::istream& in(std::istream& i)
    {return i>>A>>B;}
    //counterpart to out, does not read the class type, though.
    //This is assumed to have been read in order to generate the correct class

    CMlowrankss(std::istream& is,CH_Matrix_Classes::Integer k=0)
    {CM_type=CM_lowrankss;userkey=k;in(is);}

  };

}

#endif

