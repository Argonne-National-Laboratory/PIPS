/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/coeffmat.hxx

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



#ifndef CONICBUNDLE_COEFFMAT_HXX
#define CONICBUNDLE_COEFFMAT_HXX

//defines a base class for constraint matrices.
//the idea is to help exploiting special structures for
//the computation of tr(AiXAjZi) without having to know
//the special structure of both, Ai and Aj. 

#include "memarray.hxx"
#include "symmat.hxx"
#include "sparssym.hxx"


namespace ConicBundle {

#define GEQ  '>'
#define LEQ  '<'
#define EQU  '='

enum Coeffmattype {
                  CM_unspec=0,    //any user defined constraint may use this
                  CM_symdense=1,
                  CM_symsparse=2,
                  CM_lowrankdd=3,
                  CM_lowranksd=4,
                  CM_lowrankss=5,
                  CM_gramdense=6,
                  CM_gramsparse=7,
                  CM_singleton=8
                 };

class Coeffmat: protected CH_Matrix_Classes::Memarrayuser
{
protected:
    Coeffmattype CM_type; //in order to enable type identification
    CH_Matrix_Classes::Integer userkey;      //for the user to specify additional information
 public:
    Coeffmat(){CM_type=CM_unspec;userkey=0;}
    virtual ~Coeffmat(){}

    virtual Coeffmattype get_type() const {return CM_type;}
    virtual CH_Matrix_Classes::Integer get_userkey() const {return userkey;}
    virtual void set_userkey(CH_Matrix_Classes::Integer k) {userkey=k;}

    virtual Coeffmat* clone() const =0;
    //makes an explicit copy of itself and returns a pointer to it 

    virtual CH_Matrix_Classes::Integer dim() const =0;
    //returns the order of the represented symmetric matrix
    
    virtual CH_Matrix_Classes::Real operator()(CH_Matrix_Classes::Integer i,CH_Matrix_Classes::Integer j) const =0;
    //returns the value of the matrix element (i,j)
    
    virtual void make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const =0;
    //return dense symmetric constraint matrix

    virtual CH_Matrix_Classes::Real norm(void) const =0;
    //compute Frobenius norm of matrix

    virtual Coeffmat* subspace(const CH_Matrix_Classes::Matrix& P) const =0;
    //delivers a new object on the heap corresponding
    //to the matrix P^TAP, the caller is responsible for deleting the object
    
    virtual void multiply(CH_Matrix_Classes::Real d) =0;
    //multiply constraint permanentely by d
    //this is to allow scaling or sign changes in the constraints

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Symmatrix& S) const =0;
    //=ip(*this,S)=trace(*this*S) trace inner product
    
    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P) const =0;
    //=ip(*this,PP^T)=trace P^T(*this)P

    virtual CH_Matrix_Classes::Real gramip(const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const =0;
    //=ip(*this,QQ^T)=trace Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 

    virtual void addmeto(CH_Matrix_Classes::Symmatrix& S,CH_Matrix_Classes::Real d=1.) const =0;
    //S+=d*(*this);

    virtual void addprodto(CH_Matrix_Classes::Matrix& A,const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Real d=1.) const =0;
    //A+=d*(*this)*B

    virtual void addprodto(CH_Matrix_Classes::Matrix& A,const CH_Matrix_Classes::Sparsemat& B,CH_Matrix_Classes::Real d=1.) const =0;
    //A+=d*(*this)*B

    /// computes R=P^T*(*this)*Q
    virtual void left_right_prod(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& Q,CH_Matrix_Classes::Matrix& R) const =0;

    virtual CH_Matrix_Classes::Integer prodvec_flops() const =0;
    //return estimate of number of flops to compute addprodto for a vector

    virtual int dense() const =0;
    //returns 1 if its structure is as bad as its dense symmetric representation,
    //otherwise 0
    
    virtual int sparse() const =0;
    //returns 0 if not sparse, otherwise 1
    
    virtual int sparse(CH_Matrix_Classes::Indexmatrix& I,CH_Matrix_Classes::Indexmatrix& J,CH_Matrix_Classes::Matrix& val,CH_Matrix_Classes::Real d=1.)const=0;
    //returns 0 if not sparse. If it is sparse it returns 1 and
    //the nonzero structure in I,J and val, where val is multiplied by d.
    //Only the upper triangle (including diagonal) is delivered

    virtual int support_in(const CH_Matrix_Classes::Sparsesym& A) const =0;
    //returns 0 if the support of the costraint matrix is not contained in the
    //support of the sparse symmetric matrix A, 1 if it is contained.

    virtual CH_Matrix_Classes::Real ip(const CH_Matrix_Classes::Sparsesym& A) const =0;
    //returns the inner product of the constraint matrix with A
    
    virtual void project(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P) const=0;
    // S=P^t*A*P
    
    virtual void add_projection(CH_Matrix_Classes::Symmatrix& S,const CH_Matrix_Classes::Matrix& P,CH_Matrix_Classes::Integer start_row) const =0;
    // S+=Q^T(*this)Q for Q=P.rows(start_row,start_row+dim-1) 

    virtual const CH_Matrix_Classes::Matrix& postgenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const =0;
    // C= alpha*(*this)*B^(T if btrans) + beta*C, C is also returned

    virtual const CH_Matrix_Classes::Matrix& pregenmult(const CH_Matrix_Classes::Matrix& B,CH_Matrix_Classes::Matrix& C,
			     CH_Matrix_Classes::Real alpha=1.,CH_Matrix_Classes::Real beta=0.,int btrans=0) const =0;
    // C= alpha*B^(T if btrans)*(*this) + beta*C, C is also returned

    virtual int equal(const Coeffmat* p,double tol=1e-6) const =0;
    //returns 1, if p is the same derived class and 
    //entries differ by less than tol, otherwise zero

    virtual std::ostream& display(std::ostream& o) const =0;
    //display constraint information

    virtual std::ostream& out(std::ostream& o) const =0;
    //put entire contents onto outstream with the class type in the beginning so
    //that the derived class can be recognized.
    
    virtual std::istream& in(std::istream& i) =0;
    //counterpart to out, does not read the class type, though.
    //This is assumed to have been read in order to generate the correct class

};

Coeffmat* coeffmat_read(std::istream& in);
//reads the next Coeffmat from in into an object on the heap and
//returns a pointer to it. The caller has to destruct the object.

}

#endif

