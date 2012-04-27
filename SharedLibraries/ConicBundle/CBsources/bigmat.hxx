/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/bigmat.hxx

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


#ifndef CONICBUNDLE_BIGMAT_HXX
#define CONICBUNDLE_BIGMAT_HXX

#include <vector>
#include <map>
#include "lanczos.hxx"
#include "coeffmat.hxx"

namespace ConicBundle {

  class Bigmatrix:
    public CH_Matrix_Classes::Lanczosmatrix,
    protected CH_Matrix_Classes::Memarrayuser
  {
  private:
    typedef std::map<CH_Matrix_Classes::Integer,CH_Matrix_Classes::Real> Indexhash;

    CH_Matrix_Classes::Real    tol;        //tolerance for considering an element to be zero

    CH_Matrix_Classes::Integer dim;      //size of the matrix

    CH_Matrix_Classes::Integer nz;
    CH_Matrix_Classes::Indexmatrix colnz;
    std::vector<CH_Matrix_Classes::Indexmatrix> rowind;
    std::vector<CH_Matrix_Classes::Matrix> rowval;

    CH_Matrix_Classes::Matrix di;   //diagonal elements
    std::vector<Indexhash> rowhash;
    
    std::vector<const Coeffmat*> mcp;  //pointer to these constraints
                                       //(matrix constraints pointer)

    CH_Matrix_Classes::Matrix mcv;        //multiplier values for these constraints
    
    mutable CH_Matrix_Classes::Integer nmult;    //counts number of Mat*vec operations
    
    //---
    bool use_dense; //usually false, set to true if matrix is sufficiently dense
    mutable CH_Matrix_Classes::Symmatrix symrep;  //in this case the matrix is stored here
    bool symrep_init;
    
    //--- temporary variables for lanczosmult
    mutable CH_Matrix_Classes::Matrix At;
    mutable CH_Matrix_Classes::Matrix Bt;
    CH_Matrix_Classes::Indexmatrix collecti;
    CH_Matrix_Classes::Indexmatrix collectj;
    CH_Matrix_Classes::Matrix collectval;
    
 public:
    Bigmatrix();
    ~Bigmatrix();

    void clear();
    CH_Matrix_Classes::Integer get_nmult() const {return nmult;}
    void reset_nmult() {nmult=0;}
    
    int init(const CH_Matrix_Classes::Matrix& yin,CH_Matrix_Classes::Integer indim,const Coeffmat* C,
	     const std::map<CH_Matrix_Classes::Integer,Coeffmat*>* opA,const bool dense=false);
    //dense should be set to true if a dense
    //representation is desired

    void set_tol(CH_Matrix_Classes::Real t){tol=t;}

    int get_dense() const {return use_dense;}

       
    virtual CH_Matrix_Classes::Integer lanczosdim() const;

    virtual CH_Matrix_Classes::Integer lanczosflops() const;

    virtual int lanczosmult(const CH_Matrix_Classes::Matrix& A,CH_Matrix_Classes::Matrix& B) const;
        //computes  B = bigMatrix * A
        //A and B must not be the same object!

    int make_symmatrix(CH_Matrix_Classes::Symmatrix& S) const;

    const CH_Matrix_Classes::Symmatrix& get_symrep() const
    { make_symmatrix(symrep); return symrep; }
    
    friend std::ostream& operator<<(std::ostream& out,const Bigmatrix&);
};

} //end namespace ConicBundle

#endif

