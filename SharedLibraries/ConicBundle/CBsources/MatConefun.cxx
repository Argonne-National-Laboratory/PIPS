/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatConefun.cxx

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
#include "MatConefun.hxx"
#include <queue>
#include <utility>

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {


  int MatrixConeFunction::evaluate(const  Matrix& y, double /* relprec */,
				   Real& ubctx_ytAx, Matrix& ctx, Matrix& Ax,   
				   std::vector<PrimalData*>& primal_data 
				   )
  {
    Matrix tmpmat(c);
    genmult(A,y,tmpmat,-1.,1.,1);

    //--- find the best solution vectors
    typedef std::pair<Real,Integer> Sol;
    std::priority_queue<Sol> priq;
    for (Integer i=0;i<socdim.dim();i++){
      const Real* p=tmpmat.get_store()+socstart(i);
      Real firstval=*p++;
      Real objval=0.;
      for(Integer j=1;j<socdim(i);j++,p++){
	objval+=(*p)*(*p);
      }
      objval=sqrt(objval)+firstval;
      if (Integer(priq.size())>maxsol) {
	if (-priq.top().first<objval){
	  priq.pop();
	}
	else {
	  continue;
	}
      }
      priq.push(Sol(-objval,i));
    }

    //--- generate output information
    ctx.newsize(Integer(priq.size()),1); chk_set_init(ctx,1);
    Ax.newsize(b.dim(),Integer(priq.size())); chk_set_init(Ax,1);
    if (generate_primals){
      primal_data.resize(priq.size());
    }
    else {
      primal_data.clear();
    }      
    while(priq.size()>0){
      ubctx_ytAx=-priq.top().first;
      Integer ind=priq.top().second;
      Integer dim=socdim(ind);
      if (dim==1){
	ctx(Integer(priq.size())-1)=c(ind);
	Sparsemat xvec(c.dim(),1,1,Indexmatrix(1,1,ind),Indexmatrix(1,1,Integer(0)),Matrix(1,1,1.));	
	genmult(A,xvec,tmpvec);
	if (generate_primals){
	  primal_data[priq.size()-1]=new PrimalMatrix(xvec);
	}
      }
      else {
	//compute x (into tmpvec)
	tmpvec.newsize(dim,1); chk_set_init(tmpvec,1);
	const Real* p=tmpmat.get_store()+socstart(ind)+1;
	Real norm=0.;
	for(Integer j=1;j<dim;j++,p++){
	  tmpvec(j)=*p;
	  norm+=(*p)*(*p);
	}
	norm=sqrt(norm);
	if (norm<eps_Real){
	  tmpvec.rand(dim,1);
	  tmpvec(0)=0.;
	  Real nn2=norm2(tmpvec);
	  tmpvec/=nn2;
	}
	else {
	  tmpvec/=norm;
	}
	tmpvec(0)=1.;
	//compute ctx
	ctx(Integer(priq.size())-1)=mat_ip(dim,c.get_store()+socstart(ind),tmpvec.get_store());
	//compute Ax (into tmpvec and then into the appropriate row
	if (dim<c.dim()/3){
	  //use sparse vector
	  Sparsemat xvec(c.dim(),1,dim,Indexmatrix(Range(socstart(ind),socstart(ind)+dim-1)),Indexmatrix(dim,1,Integer(0)),tmpvec);	
	  genmult(A,xvec,tmpvec);
	  if (generate_primals){
	    primal_data[priq.size()-1]=new PrimalMatrix(xvec);
	  }
	}
	else {
	  Matrix xvec(c.dim(),1,0.);
	  mat_xey(dim,xvec.get_store()+socstart(ind),tmpvec.get_store());
	  genmult(A,xvec,tmpvec);
	  if (generate_primals){
	    primal_data[priq.size()-1]=new PrimalMatrix(xvec);
	  }
	}
      }
      mat_xey(b.dim(),Ax.get_store()+(priq.size()-1)*b.dim(),tmpvec.get_store());

      if (priq.size()<6){ 
	std::ostream* out=&std::cout; print_level=1;
	if ((out)&&(print_level>0)){
	   out->precision(6);
	   (*out)<<" socmax("<<priq.size()<<":"<<ind<<")="<<ubctx_ytAx;
	  (*out)<<std::endl;
	}
      }

      priq.pop();
    }
    return 0;
  }

  int MatrixConeFunction::subgradient_extension(const PrimalData* primal,
					      const Indexmatrix& constraint_indices, 
					      Matrix& new_subgradient_values)
  {
    new_subgradient_values=b(constraint_indices);
    Sparsemat subA=A.rows(constraint_indices);
    if ((primal==0)&&(subA.nonzeros()==0)) return 0;
    const PrimalMatrix* pr=dynamic_cast<const PrimalMatrix*>(primal);
    if (pr==0) return 1;
    genmult(subA,*pr,new_subgradient_values,-1.,1.);
    return 0;
  }

        
  int MatrixConeFunction::append_variables(
					   const Indexmatrix& append_socdim,
					   const Matrix& append_c,
					   const Sparsemat& append_A)
  {
    if ((append_socdim.coldim()!=1)||
	(sum(append_socdim)!=append_c.rowdim())||
	(append_c.coldim()!=1)||
	((append_A.rowdim()!=0)&&(append_c.rowdim()!=append_A.coldim()))||
	((append_A.rowdim()!=0)&&(append_A.rowdim()!=b.rowdim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixConeFunction::append_variables(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    if (min(append_socdim)<1){
      if (out){
	(*out)<<"**** ERROR: MatrixConeFunction::append_variables(...): dimension vector of SOC variables has entries <1"<<std::endl;
      }
      return 1;
    }
    socdim.concat_below(append_socdim);
    Integer olddim=socstart.dim();
    socstart.concat_below(append_socdim);
    Integer oldstart=(olddim>0)?socstart(olddim-1):0;
    Integer offset=(olddim>0)?socdim(olddim-1):0;
    for(Integer i=olddim;i<socstart.dim();i++){
      oldstart+=offset;
      offset=socstart(i);
      socstart(i)=oldstart;
    }
    c.concat_below(append_c);
    A.concat_right(append_A);
    return 0;
  }

  int MatrixConeFunction::append_constraints(const Sparsemat& append_A,
					   const Matrix& append_b)
  {
    if ((append_b.coldim()!=1)||
	((append_A.coldim()!=0)&&(append_A.rowdim()!=append_b.rowdim())) ||
	((append_A.coldim()!=0)&&(c.rowdim()!=append_A.coldim()))){
      if (out){
	(*out)<<"**** ERROR: MatrixConeFunction::append_constraints(...): dimensions do not match"<<std::endl;
      }
      return 1;
    }
    b.concat_below(append_b);
    A.concat_below(append_A);
    return 0;
  }
   
  int MatrixConeFunction::reassign_variables(const Indexmatrix& map_to_old)
  {
    Indexmatrix mto;
    for(Integer i=0;i<map_to_old.dim();i++){
      mto.concat_below(Range(socstart(map_to_old(i)),socstart(map_to_old(i))+socdim(map_to_old(i))-1));
    }
    socdim=socdim(map_to_old);
    socstart=socstart(map_to_old);
    A=A.cols(mto);
    c=c(mto);
    return 0;
  }
 
  int MatrixConeFunction::delete_variables(const Indexmatrix& delete_indices,
					 Indexmatrix* map_to_old)
  {
    Indexmatrix mto(Range(0,socdim.dim()-1));
    mto.delete_rows(delete_indices);
    reassign_variables(mto);
    
    if (map_to_old) {
      *map_to_old=mto;
    }
    return 0;
  }
 
  int MatrixConeFunction::reassign_constraints(const Indexmatrix& map_to_old)
  {
    A=A.rows(map_to_old);
    b=b.rows(map_to_old);
    return 0;
  }
  
  
 
  int MatrixConeFunction::delete_constraints(const Indexmatrix& delete_indices,
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

  std::ostream& MatrixConeFunction::print_problem_data(std::ostream& o)
  {
    o<<trace_value;
    switch(trace_stat){
    case Conetrace_fixed: o<<" fixed\n"; break;
    case Conetrace_bounded: o<<" bounded\n"; break;
    case Conetrace_unbounded: o<<" unbounded\n"; break;
    default: o<<" error_undefined\n";
    }
    o<<socdim;
    o<<c;
    o<<b;
    o<<A;
    return o;
  }

}

