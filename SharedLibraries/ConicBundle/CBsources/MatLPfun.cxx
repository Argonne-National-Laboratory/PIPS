/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatLPfun.cxx

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
#include "MatLPfun.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  int MatrixLPFunction::evaluate(const  Matrix& y, double /* relprec */,
				 double& objective_value,
				 Matrix& cut_values, Matrix& eps_subgradients,
				 std::vector<PrimalData*>& primal_data)
  {
    tmpmat=c;
    genmult(A,y,tmpmat,-1.,1.,1);

    //--- find the best vertex of the box
    PrimalMatrix pr(ub);
    for(Integer i=0;i<tmpmat.dim();i++){
      if (tmpmat(i)<0.){
	pr(i)=lb(i);
      }
    }

    objective_value=ip(b,y)+ip(tmpmat,pr);
    cut_values.init(1,1,objective_value);
    eps_subgradients=b;
    genmult(A,pr,eps_subgradients,-1.,1.);
    primal_data.clear();
    if (generate_primals) {
      primal_data.push_back(pr.clone_primal_data());
    }
    return 0;
  }

  int MatrixLPFunction::subgradient_extension(const PrimalData* primal,
					      const Indexmatrix& variable_indices, 
					      Matrix& new_subgradient_values)
  {
    new_subgradient_values=b(variable_indices);
    Sparsemat subA=A.rows(variable_indices);
    if ((primal==0)&&(subA.nonzeros()==0)) return 0;
    const PrimalMatrix* pr=dynamic_cast<const PrimalMatrix*>(primal);
    if (pr==0) return 1;
    genmult(subA,*pr,new_subgradient_values,-1.,1.);
    return 0;
  }

        
  int MatrixLPFunction::append_variables(const Matrix& append_c,
					 const Matrix& append_lb, 
					 const Matrix& append_ub,
					 const Sparsemat& append_A)
  {
    if ((append_c.coldim()!=1)||
	(append_lb.coldim()!=1)||
	(append_ub.coldim()!=1)||
	((append_A.rowdim()!=0)&&(append_c.rowdim()!=append_A.coldim()))||
	(append_c.rowdim()!=append_lb.rowdim())||
	(append_c.rowdim()!=append_ub.rowdim())||
	((append_A.rowdim()!=0)&&(append_A.rowdim()!=b.rowdim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixLPFunction::append_variables(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    c.concat_below(append_c);
    lb.concat_below(append_lb);
    ub.concat_below(append_ub);
    A.concat_right(append_A);
    return 0;
  }

  int MatrixLPFunction::append_constraints(const Sparsemat& append_A,
					   const Matrix& append_b)
  {
    if ((append_b.coldim()!=1)||
	((append_A.coldim()!=0)&&(append_A.rowdim()!=append_b.rowdim())) ||
	((append_A.coldim()!=0)&&(c.rowdim()!=append_A.coldim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixLPFunction::append_constraints(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    b.concat_below(append_b);
    A.concat_below(append_A);
    return 0;
  }
   
  int MatrixLPFunction::reassign_variables(const Indexmatrix& map_to_old)
  {
    A=A.cols(map_to_old);
    c=c(map_to_old);
    lb=lb(map_to_old);
    ub=ub(map_to_old);
    return 0;
  }
 
  int MatrixLPFunction::delete_variables(const Indexmatrix& delete_indices,
					 Indexmatrix* map_to_old)
  {
    if (map_to_old) {
      map_to_old->init(Range(0,c.dim()-1));
      map_to_old->delete_rows(delete_indices);
    }
    A.delete_cols(delete_indices);
    c.delete_rows(delete_indices);
    lb.delete_rows(delete_indices);
    ub.delete_rows(delete_indices);
    return 0;
  }
 
  int MatrixLPFunction::reassign_constraints(const Indexmatrix& map_to_old)
  {
    A=A.rows(map_to_old);
    b=b.rows(map_to_old);
    return 0;
  }
  
  
 
  int MatrixLPFunction::delete_constraints(const Indexmatrix& delete_indices,
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

  std::ostream& MatrixLPFunction::print_problem_data(std::ostream& o)
  {
    o<<c;
    o<<b;
    o<<A;
    o<<lb;
    o<<ub;
    return o;
  }

}

