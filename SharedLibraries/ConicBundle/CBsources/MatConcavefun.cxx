/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatConcavefun.cxx

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
#include "MatConcavefun.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  int OnedimPiecewiseLinearFunction::compute_deriv()
  {
    if ((breakpoints.dim()<=1)||(fun_value.dim()!=breakpoints.dim())){
      //exit(1);
      abort();
    }
    derivative.newsize(breakpoints.dim()-1,1);
    chk_set_init(derivative,1);
    Integer cnt=0;
    for (Integer i=0;i<derivative.dim();i++){
      if (breakpoints(i+1)-breakpoints(cnt)<=eps_Real) return 2;
      Real deriv=0.;
      ++cnt;
      do {  //concavify if necessary
	--cnt;
	deriv=(fun_value(i+1)-fun_value(cnt))/(breakpoints(i+1)-breakpoints(cnt));
      } while ((cnt>=1)&&(deriv>=derivative(cnt-1)));
      derivative(cnt)=deriv;
      cnt++;
      breakpoints(cnt)=breakpoints(i+1);
      fun_value(cnt)=fun_value(i+1);
    }
    breakpoints.reduce_length(cnt+1);
    fun_value.reduce_length(cnt+1); 
    derivative.reduce_length(cnt);  
    return 0;
  }
    
  Integer OnedimPiecewiseLinearFunction::find_breakpoint(Real deriv) const
  {
    Integer lb=0;
    Integer ub=breakpoints.dim()-1;
    while (lb<ub){
      Integer i=(lb+ub)/2;            //therefore i<ub
      if (derivative(i)+deriv>=0){
	lb=i+1;                       //and i<lb<=ub
      }
	else {
	  ub=i;                         //or lb<=i=ub
	}
    }
    return lb;
  }
      
    
    
  int OnedimPiecewiseLinearFunction::evaluate(Real add_linterm,Real /* relprec */,
					      Real& objective_value, Real& solution) const
  {
    Integer bpi=find_breakpoint(add_linterm);
    solution = breakpoints(bpi);
    objective_value=fun_value(bpi)+solution*add_linterm;
    return 0;
  }

  Real OnedimPiecewiseLinearFunction::primal_evaluate(Real x) const
  {
    if (breakpoints.dim()==0) 
       return CB_minus_infinity;
    Integer ilb=0;
    Real lb=breakpoints(ilb);
    if (x<lb) 
      return CB_minus_infinity;
    Integer iub=breakpoints.dim()-1;
    Real ub=breakpoints(iub);
    if (x>ub) 
      return CB_minus_infinity;
    //now lb<=x<=ub
    while (ilb<iub-1){
      Integer i=(ilb+iub)/2;            //therefore i<ub
      Real cb=breakpoints(i);
      if (x>cb-1e-10){
	ilb=i;                       //and i<lb<=ub
	lb=cb;
      }
      else {
	iub=i;                         //or lb<=i=ub
	ub=cb;
      }
      //now lb<=x<=ub
    }
    if (derivative.dim()>0){
      return fun_value(ilb)+derivative(ilb)*(x-lb);
    }
    return fun_value(ilb);
  }

  int OnedimPiecewiseLinearFunction::primal_evaluate(Real x, Real& fv,
						     Real& ld,Real& rd) const
  {
    if (breakpoints.dim()==0) {
      fv=rd=CB_minus_infinity;
      ld=CB_plus_infinity;
      return 0;
    }
    Integer ilb=0;
    Real lb=breakpoints(ilb);
    if (x<lb) {
      fv=CB_minus_infinity;
      ld=rd=CB_plus_infinity;
      return 0;
    }
    Integer iub=breakpoints.dim()-1;
    Real ub=breakpoints(iub);
    if (x>ub) {
      fv=CB_minus_infinity;
      ld=rd=CB_minus_infinity;
      return 0;
    }
    //now lb<=x<=ub
    while (ilb<iub-1){
      Integer i=(ilb+iub)/2;            //therefore i<ub
      Real cb=breakpoints(i);
      if (x>cb-1e-10){
	ilb=i;                       //and i<lb<=ub
	lb=cb;
      }
      else {
	iub=i;                         //or lb<=i=ub
	ub=cb;
      }
      //now lb<=x<=ub
    }
    if (derivative.dim()>0){
      ld=rd=derivative(ilb);
      fv=fun_value(ilb)+ld*(x-lb);
    }
    else 
      fv=fun_value(ilb);
    if (x<=lb){
      if (ilb==0)
	ld=CB_plus_infinity;
      else
	ld=derivative(ilb-1);
    }
    if (x>=ub){
      if (iub==derivative.dim())
	rd=CB_plus_infinity;
      else
	rd=derivative(iub);
    }
    return 0;
  }


  ConcaveCostLPFunction::~ConcaveCostLPFunction()
  {
    for(unsigned int i=0;i<concavefun.size();i++){
      delete concavefun[i];
    }
    concavefun.clear();
  }


  Real ConcaveCostLPFunction::primal_evaluate(const Matrix& x) const
  {
    assert(x.dim()==Integer(concavefun.size()));
    Real val=0.;
    for(unsigned int i=0;i<concavefun.size();i++){
      Real lval=concavefun[i]->primal_evaluate(x(i));
      if (lval<=CB_minus_infinity) 
	return lval;
      val+= lval;
    }
    return val;
  }

  int ConcaveCostLPFunction::primal_evaluate(const Matrix& x,Matrix &fun_value,Matrix& left_deriv,Matrix& right_deriv) const
  {
    assert(x.dim()==Integer(concavefun.size()));
    fun_value.newsize(x.dim(),1); chk_set_init(fun_value,1);
    left_deriv.newsize(x.dim(),1); chk_set_init(left_deriv,1);
    right_deriv.newsize(x.dim(),1); chk_set_init(right_deriv,1);
    for(unsigned int i=0;i<concavefun.size();i++){
      if (concavefun[i]->primal_evaluate(x(i),fun_value(i),
					 left_deriv(i),right_deriv(i)))
	return i+1;
    }
    return 0;
  }


  int ConcaveCostLPFunction::evaluate(const  Matrix& y, Real relprec,
				      Real& objective_value,
				      Matrix& cut_values, Matrix& eps_subgradients,
				      std::vector<PrimalData*>& primal_data,
				      PrimalExtender*&)
  {
    Matrix tmpmat;
    genmult(A,y,tmpmat,-1.,0.,1);

    //--- find the best vertex of the box
    PrimalMatrix pr(tmpmat.dim(),1,0.);
    Real objsum=0;
    for(Integer i=0;i<tmpmat.dim();i++){
      Real d;
      concavefun[i]->evaluate(tmpmat(i),relprec,d,pr(i));
      objsum+=d;
    }

    objective_value=objsum+ip(b,y);
    cut_values.init(1,1,objective_value);
    eps_subgradients=b;
    genmult(A,pr,eps_subgradients,-1.,1.);
    primal_data.clear();
    if (generate_primals) {
      primal_data.push_back(pr.clone_primal_data());
    }
    return 0;
  }

  int ConcaveCostLPFunction::subgradient_extension(const PrimalData* primal,
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

        
  int ConcaveCostLPFunction::append_variables(const ConcaveFunVector& append_f,
					      const Sparsemat& append_A)
  {
    if (((append_A.rowdim()!=0)&&(Integer(append_f.size())!=append_A.coldim()))||
	(append_A.rowdim()!=b.rowdim())){
      if (out){
	(*out)<<"**** ERROR: ConcaveCostLPFunction::append_variables(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    concavefun.insert(concavefun.end(),append_f.begin(),append_f.end());
    A.concat_right(append_A);
    return 0;
  }

  int ConcaveCostLPFunction::append_constraints(const Sparsemat& append_A,
					   const Matrix& append_b)
  {
    if ((append_b.coldim()!=1)||
	((append_A.rowdim()!=0)&&(append_A.rowdim()!=append_b.rowdim())) ||
	((append_A.coldim()!=0)&&(Integer(concavefun.size())!=append_A.coldim()))){
      if (out){
	(*out)<<"**** ERROR: ConcaveCostLPFunction::append_constraints(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    b.concat_below(append_b);
    A.concat_below(append_A);
    return 0;
  }
   
  int ConcaveCostLPFunction::reassign_variables(const Indexmatrix& map_to_old)
  {
    ConcaveFunVector tmp=concavefun;
    concavefun.resize(map_to_old.dim());
    for(unsigned int i=0;i<unsigned(map_to_old.dim());i++){
      concavefun[i]=tmp[i]->clone();
    }
    {for(unsigned int i=0;i<tmp.size();i++){
      delete tmp[i];
    }}
    tmp.clear();
    A=A.cols(map_to_old);
    return 0;
  }
 
  int ConcaveCostLPFunction::delete_variables(const Indexmatrix& delete_indices,
					 Indexmatrix* map_to_old)
  {
    Indexmatrix mto(Range(0,A.coldim()-1));
    mto.delete_rows(delete_indices);
    if (map_to_old) {
      *map_to_old=mto;
    }
    reassign_variables(mto);
    return 0;
  }
 
  int ConcaveCostLPFunction::reassign_constraints(const Indexmatrix& map_to_old)
  {
    A=A.rows(map_to_old);
    b=b.rows(map_to_old);
    return 0;
  }
  
  
 
  int ConcaveCostLPFunction::delete_constraints(const Indexmatrix& delete_indices,
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

}

/* ---------------------------------------------------------------------------
 *    Change log $Log: MatLPfun.cc,v $
 *    Change log Revision 1.1.1.1.2.1  2000/02/28 13:55:44  bzfhelmb
 *    Change log Reording inline functions and correction of minor details
 *    Change log
 *    Change log Revision 1.1.1.1  1999/12/15 13:47:50  bzfhelmb
 *    Change log Imported sources
 *    Change log
 *
 *    End of $Source: /kombadon/cvs/cvsroot/bzfhelmb/SBmethod/spectral/MatLPfun.cc,v $
 --------------------------------------------------------------------------- */
