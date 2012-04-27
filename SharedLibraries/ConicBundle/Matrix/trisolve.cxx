/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  Matrix/trisolve.cxx

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
#include "heapsort.hxx"
#include "gb_rand.hxx"

 
using namespace CH_Tools;

namespace CH_Matrix_Classes {


// **************************************************************************
//                                 triu_solve
// **************************************************************************

//On input matrix X is the right hand side of an equation (*this)*X=rhs where
//only the upper triangle of the first n lines of (*this) forming a square matrix
//is used. On output X is the solution.
//Back substitution.
//returns:
//  0      if all diagonal elements were larger than tol in absolut value ->ok
//  index+1  of the first diag elem violating this condition -> error

int Matrix::triu_solve(Matrix& X,Real tol)
{
 chk_init(*this);
 chk_init(X);
 Integer n=min(nr,nc); 
#if (CONICBUNDLE_DEBUG>=1)
 if ((X.nr!=n)||(X.nc<1)) MEmessage(MEdim(nr,nc,X.nr,X.nc,"Matrix::triu_solve(): rhs does no match min(nr,nc)",MTmatrix));
#endif
 Integer i,j;
 Real d;
 Real *xp;
 for(i=n;--i>=0;){  //solve for variable row i
     d=(*this)(i,i);
     if (::fabs(d)<tol) {
         xp=X.m+i;
         j=0;
         while((j<X.nc)&&(::fabs(*xp)<tol)) {*xp=0.; j++; xp+=X.nr;}
         if (j<X.nc) return i+1;
         continue;
     }
     xp=X.m+i;
     for(j=0;j<X.nc;j++,xp+=X.nr){
         (*xp)/=d;
         mat_xpeya(i,X.m+j*X.nr,m+i*nr,-(*xp));
     }
 }
 return 0;
}

// **************************************************************************
//                                 tril_solve
// **************************************************************************

//On input matrix X is the right hand side of an equation (*this)*X=rhs where
//only the lower triangle of the first n lines of (*this) forming a square matrix
//is used. On output X is the solution.
//Forward substitution.
//returns:
//  0      if all diagonal elements were larger than tol in absolut vlaue ->ok
//  index+1  of the first diag elem violating this condition -> error

int Matrix::tril_solve(Matrix& X,Real tol)
{
 chk_init(*this);
 chk_init(X);
 Integer n=min(nr,nc); 
#if (CONICBUNDLE_DEBUG>=1)
 if ((X.nr!=n)||(X.nc<1)) MEmessage(MEdim(nr,nc,X.nr,X.nc,"Matrix::tril_solve(): rhs does not match min(nr,nc)"));
#endif
 chk_rowseq(*this,X);
 Integer i,j;
 Real d;
 Real *xp;
 for(i=0;i<n;i++){  //solve for variable row i
     d=(*this)(i,i);
     if (::fabs(d)<tol)  {
         xp=X.m+i;
         j=0;
         while((j<X.nc)&&(::fabs(*xp)<tol)) {*xp=0.; j++; xp+=X.nr;}
         if (j<X.nc) return i+1;
         continue;
     }
     xp=X.m+i;
     for(j=0;j<X.nc;j++,xp+=X.nr){
         (*xp)/=d;
         mat_xpeya(nr-i-1,X.m+j*X.nr+i+1,m+i*nr+i+1,-(*xp));
     }
 }
 return 0;
}

}

