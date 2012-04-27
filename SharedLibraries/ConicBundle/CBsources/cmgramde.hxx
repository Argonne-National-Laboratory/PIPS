/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/cmgramde.hxx

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


#ifndef CONICBUNDLE_CMGRAMDE_HXX
#define CONICBUNDLE_CMGRAMDE_HXX

//defines a base class for constraint matrices.
//the idea is to help exploiting special structures for
//the computation of tr(AiXAjZi) without having to know
//the special structure of both, Ai and Aj. 

#include "coeffmat.hxx"


namespace ConicBundle{
 
  class CMgramdense: public Coeffmat
  {
  private:
    CH_Matrix_Classes::Matrix A;  //Symmatrix=A*A^T
    bool positive;
  public:
    CMgramdense(const CH_Matrix_Classes::Matrix& Ain,bool pos=true,CH_Matrix_Classes::Integer k=0)
    {A=Ain;positive=pos;CM_type=CM_gramdense;userkey=k;}
    virtual ~CMgramdense(){}
    
    //virtual CM_type get_type() const {return CM_type;} //no need to overload
    //virtual CH_Matrix_Classes::Integer get_userkey() const {return userkey;}
    //virtual void set_userkey(CH_Matrix_Classes::Integer k) const {userkey=k;}

    virtual Coeffmat* clone() const 
    {return new CMgramdense(A,positive,userkey);}
    //makes an explicit copy of itself and returns a pointer to it 
    
    virtual CH_Matrix_Classes::Integer dim() const  { return A.rowdim(); }
    //returns the order of the represented symmetric matrix
    
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j) const 
    {return CH_Matrix_Classes::mat_ip(A.coldim(),A.get_store()+i,A.rowdim(),A.get_store()+j,A.rowdim());}
    //returns the value of the matrix element (i,j)
    
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const
    {if (positive) { CH_Matrix_Classes::rankadd(A,S);} else { CH_Matrix_Classes::rankadd(A,S,-1.);}}
    //return dense symmetric constraint matrix

    virtual CH_Matrix_Classes::Real norm(void) const
    {CH_Matrix_Classes::Symmatrix B; return CH_Matrix_Classes::norm2(CH_Matrix_Classes::rankadd(A,B,1.,0.,1));}
    //compute Frobenius norm of matrix

    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const
    { CH_Matrix_Classes::Matrix B; CH_Matrix_Classes::genmult(P,A,B,1.,0.,1); return new CMgramdense(B,positive,userkey); }
    //delivers a new object on the heap corresponding
    //to the matrix P^TAP, the caller is responsible for deleting the object
    
    virtual void multiply(CH_Matrix_Classes::Real d)
    { if (d<0.) {A*=sqrt(-d); positive=!positive;} else {A*=sqrt(d);} }
    //multiply constraint permanentely by d
    //this is to allow scaling or sign changes in the constraints

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const
    {CH_Matrix_Classes::Matrix B; if (positive) return CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S,A,B),A);
     else return -CH_Matrix_Classes::ip(CH_Matrix_Classes::genmult(S,A,B),A); }
    //=ip(*this,S)=trace(*this*S) trace inner product
    
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const
    {CH_Matrix_Classes::Matrix B; CH_Matrix_Classes::genmult(P,A,B,1.,0.,1);
    if (positive) return CH_Matrix_Classes::ip(B,B); else return -CH_Matrix_Classes::ip(B,B);}
    //=ip(*this,PP^T)=trace P^T(*this)P

    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    {
      CH_Matrix_Classes::Matrix B(A.coldim(),P.coldim()); //B=A^T*P
      CH_Matrix_Classes::Real *bp=B.get_store();
      for(CH_Matrix_Classes::Integer j=0;j<B.coldim();j++){
	const CH_Matrix_Classes::Real *pp=P.get_store()+start_row+j*P.rowdim();
	const CH_Matrix_Classes::Real *ap=A.get_store();
	CH_Matrix_Classes::Integer ard=A.rowdim();
	for(CH_Matrix_Classes::Integer i=B.rowdim();--i>=0;){
	  (*bp++)=CH_Matrix_Classes::mat_ip(ard,ap,pp);
	  ap+=ard;
	}
      }	
      chk_set_init(B,1);
      if (positive) return CH_Matrix_Classes::ip(B,B); else return -CH_Matrix_Classes::ip(B,B);
    }
    //=ip(*this,PP^T)=trace P^T(*this)P

    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S,CH_Matrix_Classes::Real d=1.) const
    {if (positive) CH_Matrix_Classes::rankadd(A,S,d,1.); else CH_Matrix_Classes::rankadd(A,S,-d,1.);}
    //S+=d*(*this);

    virtual void addprodto(CH_Matrix_Classes::Matrix& B,const CH_Matrix_Classes::Matrix&C ,CH_Matrix_Classes::Real d=1.) const
    {CH_Matrix_Classes::Matrix D;
     if (positive) CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(A,C,D,1.,0.,1),B,d,1.);
     else CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(A,C,D,1.,0.,1),B,-d,1.);}
    //B+=d*(*this)*C

    virtual void addprodto(CH_Matrix_Classes::Matrix& B,const CH_Matrix_Classes::Sparsemat&C ,CH_Matrix_Classes::Real d=1.) const
    {CH_Matrix_Classes::Matrix D;
     if (positive) CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(A,C,D,1.,0.,1),B,d,1.);
     else CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(A,C,D,1.,0.,1),B,-d,1.);}
    //B+=d*(*this)*C

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& Q,CH_Matrix_Classes::Matrix& R) const
    {
      CH_Matrix_Classes::Matrix tmp1; CH_Matrix_Classes::genmult(P,A,tmp1,1.,0.,1,0);
      CH_Matrix_Classes::Matrix tmp2; CH_Matrix_Classes::genmult(A,Q,tmp2,1.,0.,1,0);
      if (positive) CH_Matrix_Classes::genmult(tmp1,tmp2,R,1.,0.,0,0);
      else CH_Matrix_Classes::genmult(tmp1,tmp2,R,-1.,0.,0,0);
    }

    virtual CH_Matrix_Classes::Integer prodvec_flops() const 
    { return 4*A.rowdim()*A.coldim(); }
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
		       CH_Matrix_Classes::Real /* d=1. */ )const
    {return 0;}
    //returns 0 if not sparse. If it is spars it returns 1 and
    //the nonzero structure in I,J and val, where val is multiplied by d.
    //Only the upper triangle (including diagonal) is delivered
    
    virtual int support_in(const CH_Matrix_Classes::Sparsesym& /* S */) const
    {return 0;}
    //returns 0 if the support of the costraint matrix is not contained in the
    //support of the sparse symmetric matrix A, 1 if it is contained.

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& S) const
    {if (positive) return CH_Matrix_Classes::ip(A,S*A); else return -CH_Matrix_Classes::ip(A,S*A);}
    //returns the inner product of the constraint matrix with A
    
    virtual void project(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P) const
    {CH_Matrix_Classes::Matrix B; CH_Matrix_Classes::genmult(P,A,B,1.,0.,1);
     if (positive) CH_Matrix_Classes::rankadd(B,S); else CH_Matrix_Classes::rankadd(B,S,-1.);
    }
    // S=P^t*A*P

    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const
    {
      CH_Matrix_Classes::Matrix B(A.coldim(),P.coldim()); //B=A^T*P
      CH_Matrix_Classes::Real *bp=B.get_store();
      for(CH_Matrix_Classes::Integer j=0;j<B.coldim();j++){
	const CH_Matrix_Classes::Real *pp=P.get_store()+start_row+j*P.rowdim();
	const CH_Matrix_Classes::Real *ap=A.get_store();
	CH_Matrix_Classes::Integer ard=A.rowdim();
	for(CH_Matrix_Classes::Integer i=B.rowdim();--i>=0;){
	  (*bp++)=CH_Matrix_Classes::mat_ip(ard,ap,pp);
	  ap+=ard;
	}
      }	
      chk_set_init(B,1);
      if (positive) CH_Matrix_Classes::rankadd(B,S,-1.,1.,1); else CH_Matrix_Classes::rankadd(B,S,-1.,1.,1);
    }
    // S=Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 

    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const
    { 
      CH_Matrix_Classes::Matrix D; 
      return CH_Matrix_Classes::genmult(A,CH_Matrix_Classes::genmult(A,B,D,1.,0.,1,btrans),C,alpha,beta);
    }
    // C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned

    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const
    { 
      CH_Matrix_Classes::Matrix D; 
      return CH_Matrix_Classes::genmult(CH_Matrix_Classes::genmult(B,A,D,1.,0.,btrans),A,C,alpha,beta,0,1);
    }
    // C= alpha*B^(T if btrans)*(*this) + beta*C 

    virtual int equal(const Coeffmat* p,double tol=1e-6) const
    {
      const CMgramdense *pp=dynamic_cast<const CMgramdense *>(p);
      if (pp==0) 
	return 0;
      if ((A.rowdim()!=(pp->A).rowdim())||(A.coldim()!=(pp->A).coldim()))
	return 0;
      return (CH_Matrix_Classes::norm2(A-pp->A)<tol);
    }
    //return true, if p is the same derived class and 
    //entries differ by less than tol

    virtual std::ostream& display(std::ostream& o) const 
    {o<<"CMgramdense\n"; A.display(o); return o;}
    //display constraint information
    
    virtual std::ostream& out(std::ostream& o) const
    {return o<<"GRAM_DENSE\n"<<positive<<"\n"<<A;}
    //put entire contents onto outstream with the class type in the beginning so
    //that the derived class can be recognized.

    virtual std::istream& in(std::istream& i)
    {return i>>positive>>A;}
    //counterpart to out, does not read the class type, though.
    //This is assumed to have been read in order to generate the correct class

    CMgramdense(std::istream& is,CH_Matrix_Classes::Integer k=0)
    {CM_type=CM_gramdense;userkey=k;in(is);}

    //--- specific routines
    const CH_Matrix_Classes::Matrix& get_A() const {return A;}
    bool get_positive() const {return positive;}
  };

}
#endif

