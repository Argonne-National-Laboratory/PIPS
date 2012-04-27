/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatLPBCfun.cxx

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
#include "MatLPBCfun.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  std::ostream& BoundedConeData::print_data(std::ostream& out)
  {
    switch(conetype){
    case CT_lp: out<<" 0 "; break;
    case CT_soc: out<<" 1 "; break;
    case CT_sdp: out<<" 2 ";
    }
    int oldprec=int(out.precision(12));
    out<<dim<<" "<<lb<<" "<<ub;
    out.precision(oldprec);
    return out;
  }

  int MatrixLPBCFunction::evaluate(const  Matrix& y, double /* relprec */,
				 double& objective_value,
				 Matrix& cut_values, Matrix& eps_subgradients,
				 std::vector<PrimalData*>& primal_data)
  {
    tmpmat=c;
    genmult(A,y,tmpmat,-1.,1.,1);

    //--- find the best vertex of the box
    PrimalMatrix pr(c.dim(),1); chk_set_init(pr,1);
    Integer ind=0;
    for(BCDvector::iterator bci=conedata.begin();bci!=conedata.end();bci++){
      switch(bci->conetype){
      case CT_lp:
	{
	  pr(ind)=(tmpmat(ind)<0)?bci->lb:bci->ub; 
	  ind++; 
	  break;
	}
      case CT_soc: 
	{
	  tmpvec.init(bci->dim,1,tmpmat.get_store()+ind);
	  tmpvec(0)=0.;
	  Real n2=norm2(tmpvec);
	  while(n2<1e4*eps_Real){
	    tmpvec.rand(bci->dim,1);
	    tmpvec(0)=0.;
	    n2=norm2(tmpvec);
	  }
	  tmpvec/=n2;
	  tmpvec(0)=1.;
	  if (mat_ip(bci->dim,tmpmat.get_store()+ind,tmpvec.get_store())<0.){
	    mat_xeya(bci->dim,pr.get_store()+ind,tmpvec.get_store(),bci->lb);
	  }
	  else {
	    mat_xeya(bci->dim,pr.get_store()+ind,tmpvec.get_store(),bci->ub);
	  }
	  ind+=bci->dim;
	  break;
	}
      case CT_sdp:
	{
	  Integer sdim=(bci->dim*(bci->dim+1))/2;
	  tmpvec.init(sdim,1,tmpmat.get_store()+ind);
	  sveci(tmpvec,tmpsym);
	  Matrix P;
	  tmpsym.eig(P,tmpvec,false); //sort non_increasingly
	  if (tmpvec(0)<0.){
	    rankadd(P.col(0),tmpsym,bci->lb);
	  }
	  else {
	    rankadd(P.col(0),tmpsym,bci->ub);
	  }
	  svec(tmpsym,tmpvec);
	  mat_xey(sdim,pr.get_store()+ind,tmpvec.get_store());
	  ind+=sdim;
	}
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

  int MatrixLPBCFunction::subgradient_extension(const PrimalData* primal,
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

        
  int MatrixLPBCFunction::append_variables(const BCDvector& append_conedata,
					   const Matrix& append_c,
					   const Sparsemat& append_A)
  {
    if ((append_c.coldim()!=1)||
	((append_A.rowdim()!=0)&&(append_c.rowdim()!=append_A.coldim()))||
	((append_A.rowdim()!=0)&&(append_A.rowdim()!=b.rowdim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixLPBCFunction::append_variables(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    Integer adim=0;
    for(unsigned int i=0;i<append_conedata.size();i++){
      if (append_conedata[i].dim<1) {
	if (out){
	  (*out)<<"**** ERROR: MatrixLPBCFunction::append_variables(...): append_conedata["<<i<<"].dim<1"<<std::endl;
	}
	return 1;
      }

      switch(append_conedata[i].conetype){
      case CT_lp: 
	if (append_conedata[i].dim>1){
	  if (out){
	    (*out)<<"**** ERROR: MatrixLPBCFunction::append_variables(...): append_conedata["<<i<<"].dim>1 for LP-variable"<<std::endl;
	  }
	  return 1;
	}  
	adim++; 
	break;
      case CT_soc:
	adim+=append_conedata[i].dim;
	break;
      case CT_sdp:
	adim+=(append_conedata[i].dim*(append_conedata[i].dim+1))/2;
      }
	
      if ((append_conedata[i].ub<0.)&&(append_conedata[i].conetype!=CT_lp)){
	if (out){
	  (*out)<<"**** ERROR: MatrixLPBCFunction::append_variables(...): append_conedata["<<i<<"].ub<0. but not LP-variable"<<std::endl;
	}
	return 1;
      }

      if ((append_conedata[i].lb<0.)&&(append_conedata[i].conetype!=CT_lp)){
 	if (out){
 	  (*out)<<"**** ERROR: MatrixLPBCFunction::append_variables(...): append_conedata["<<i<<"].lb<0. but not LP-variable"<<std::endl;
 	}
 	return 1;
       }
      
      if (append_conedata[i].lb>append_conedata[i].ub){
	if (out){
	  (*out)<<"**** ERROR: MatrixLPBCFunction::append_variables(...): append_conedata["<<i<<"]: ub<lb, trivially infeasible"<<std::endl;
	}
	return 1;
      }
      
    }

    if (adim!=append_c.rowdim()){
      if (out){
	(*out)<<"**** ERROR: MatrixLPBCFunction::append_variables(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }

    {unsigned int i=(unsigned int)(conedata.size());
    conedata.insert(conedata.end(),append_conedata.begin(),append_conedata.end());
    for(;i<conedata.size();i++){
      conedata[i].lb=max(conedata[i].lb,CB_minus_infinity);
      conedata[i].lb=min(conedata[i].lb,CB_plus_infinity);
      conedata[i].ub=max(conedata[i].ub,CB_minus_infinity);
      conedata[i].ub=min(conedata[i].ub,CB_plus_infinity);
    }}
    c.concat_below(append_c);
    A.concat_right(append_A);
    return 0;
  }

  int MatrixLPBCFunction::append_constraints(const Sparsemat& append_A,
					     const Matrix& append_b)
  {
    if ((append_b.coldim()!=1)||
	((append_A.coldim()!=0)&&(append_A.rowdim()!=append_b.rowdim())) ||
	((append_A.coldim()!=0)&&(c.rowdim()!=append_A.coldim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixLPBCFunction::append_constraints(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    b.concat_below(append_b);
    A.concat_below(append_A);
    return 0;
  }
   
  int MatrixLPBCFunction::reassign_variables(const Indexmatrix& map_to_old)
  {
    BCDvector nconedata;
    Indexmatrix col_to_old;
    Indexmatrix col_ind(Integer(conedata.size())+1,1,Integer(0));
    Integer ind=0;
    for(unsigned int i=0;i<conedata.size();i++){
      switch(conedata[i].conetype){
      case CT_lp: ind++; break;
      case CT_soc: ind+=conedata[i].dim; break;
      case CT_sdp: ind+=(conedata[i].dim*(conedata[i].dim+1)/2);
      }
      col_ind(i+1)=ind;
    }
      
    {for(Integer i=0;i<map_to_old.dim();i++){
      nconedata.push_back(conedata[map_to_old(i)]);
      col_to_old.concat_below(Range(col_ind(map_to_old(i)),col_ind(map_to_old(i+1))-1));
    }}
    conedata.swap(nconedata);
    A=A.cols(col_to_old);
    c=c(col_to_old);
    return 0;
  }
 
  int MatrixLPBCFunction::delete_variables(const Indexmatrix& delete_indices,
					 Indexmatrix* map_to_old)
  {
    Indexmatrix tmpind(Range(0,Integer(conedata.size())-1));
    tmpind.delete_rows(delete_indices);
    if (map_to_old) {
      *map_to_old=tmpind;
    }
    return reassign_variables(tmpind);
  }
 
  int MatrixLPBCFunction::reassign_constraints(const Indexmatrix& map_to_old)
  {
    A=A.rows(map_to_old);
    b=b.rows(map_to_old);
    return 0;
  }
  
  
 
  int MatrixLPBCFunction::delete_constraints(const Indexmatrix& delete_indices,
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

  std::ostream& MatrixLPBCFunction::print_problem_data(std::ostream& o)
  {
    o<<"\nBEGIN_LPBC_PROBLEM\n";
    o<<"\nBLOCKS\n";
    for(unsigned int i=0;i<conedata.size();i++){
      conedata[i].print_data(o)<<"\n";
    }
    o<<"\nCOST_MATRIX\n";
    o<<c;
    o<<"\nCONSTRAINT_MATRIX\n";
    o<<A;
    o<<"\nRIGHT_HAND_SIDE\n";
    o<<b;
    o<<"\nEND_LPBC_PROBLEM"<<std::endl;
    return o;
  }

}

