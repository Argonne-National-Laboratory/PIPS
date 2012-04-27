/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatSOCfun.cxx

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



#include <math.h>
#include <stdlib.h>
#include "MatSOCfun.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {
        
  int MatrixSOCFunction::append_variables(const Matrix& append_c,
					 const Sparsemat& append_A)
  {
    if ((append_c.coldim()!=1)||
	((append_A.rowdim()!=0)&&(append_c.rowdim()!=append_A.coldim()))||
	((append_A.rowdim()!=0)&&(append_A.rowdim()!=b.rowdim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixSOCFunction::append_variables(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    c.concat_below(append_c);
    A.concat_right(append_A);
    return 0;
  }

  int MatrixSOCFunction::append_constraints(const Sparsemat& append_A,
					   const Matrix& append_b)
  {
    if ((append_b.coldim()!=1)||
	((append_A.coldim()!=0)&&(append_A.rowdim()!=append_b.rowdim())) ||
	((append_A.coldim()!=0)&&(c.rowdim()!=append_A.coldim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixSOCFunction::append_constraints(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    b.concat_below(append_b);
    A.concat_below(append_A);
    return 0;
  }
   
  int MatrixSOCFunction::reassign_variables(const Indexmatrix& map_to_old)
  {
    A=A.cols(map_to_old);
    c=c(map_to_old);
    return 0;
  }
 
  int MatrixSOCFunction::delete_variables(const Indexmatrix& delete_indices,
					 Indexmatrix* map_to_old)
  {
    if (map_to_old) {
      map_to_old->init(Range(0,c.dim()-1));
      map_to_old->delete_rows(delete_indices);
    }
    A.delete_cols(delete_indices);
    c.delete_rows(delete_indices);
    return 0;
  }
 
  int MatrixSOCFunction::reassign_constraints(const Indexmatrix& map_to_old)
  {
    A=A.rows(map_to_old);
    b=b.rows(map_to_old);
    return 0;
  }
  
  
 
  int MatrixSOCFunction::delete_constraints(const Indexmatrix& delete_indices,
					   Indexmatrix* map_to_old)
  {
    if (map_to_old) {
      map_to_old->init(Range(0,b.dim()-1));
      map_to_old->delete_rows(delete_indices);
    }
    A.delete_rows(delete_indices);
    b.delete_rows(delete_indices);
    return 0;
  }

  std::ostream& MatrixSOCFunction::print_problem_data(std::ostream& o)
  {
    if (constant_trace) o<<" = ";
    else o<<" < ";
    o<<trace_value<<"\n";
    o<<c;
    o<<b;
    o<<A;
    return o;
  }

}

