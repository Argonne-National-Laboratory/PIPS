/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/CFunction.cxx

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



#include "CFunction.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

  CFunction::CFunction(void* fk,cb_functionp fp,cb_subgextp se,int prdim)
  {
    function_key=fk;
    assert(fp!=0);
    oracle=fp;
    subgext=se;
    primaldim=prdim;
    max_new=1;
  }
       
    /** evaluation function
        for each i: cut_values[i] and epssubgradients[i] describe a
        linear minorant of the convex function, at least one must be returned
	@return 0 on success */
  int CFunction::evaluate(
			  const  Matrix& current_point,
			  double relprec,
			  double&  obval,
			  Matrix&  val,   
			  Matrix&  subg,
			  std::vector<PrimalData*>&    pdata,
			  PrimalExtender*& 
			)
  {
    Integer dim=current_point.dim();
    subg.newsize(dim,Integer(max_new)); chk_set_init(subg,1); 
    val.newsize(max_new,1);             chk_set_init(val,1);
    
    int n_new;
    int ret_code;
    if (primaldim>0){
      Matrix x(primaldim,max_new);  chk_set_init(x,1);
      ret_code =(*oracle)(
			  function_key,
			  const_cast<double *>(current_point.get_store()),
			  relprec,
			  max_new,
			  &obval,
			  &n_new,
			  val.get_store(),
			  subg.get_store(),
			  x.get_store()
			  );
      pdata.resize(n_new);
      for (Integer i=0;i<n_new;i++){
	pdata[i]=new PrimalMatrix(x.col(i));
      }
    }
    else {
      ret_code =(*oracle)(
			  function_key,
			  const_cast<double *>(current_point.get_store()),
			  relprec,
			  max_new,
			  &obval,
			  &n_new,
			  val.get_store(),
			  subg.get_store(),
			  0
			  );
      pdata.clear();
    }
    if (n_new<max_new){
      val.delete_rows(Range(Integer(n_new),Integer(max_new-1)));
      subg.delete_cols(Range(Integer(n_new),Integer(max_new-1)));
    }
    
    return ret_code;
  }



  int CFunction::subgradient_extension(
				       const PrimalData* generating_primal,
				       const Indexmatrix& variable_indices, 
				       Matrix& new_subgradient_values
				       )
  { 
    if (subgext==0) return 1;
    int ret_code;
    new_subgradient_values.init(variable_indices.dim(),1,0.);
    if (generating_primal==0){
      ret_code = subgext(
			 function_key,
			 0,
			 variable_indices.dim(),
			 const_cast<Integer *>(variable_indices.get_store()),
			 new_subgradient_values.get_store()
			 );
      return ret_code;
    }
      
    const PrimalMatrix* xp=dynamic_cast<const PrimalMatrix *>(generating_primal);
    assert(xp!=0);
    assert(xp->dim()==primaldim);
    ret_code = subgext(
		       function_key,
		       const_cast<double *>(xp->get_store()),
		       variable_indices.dim(),
		       const_cast<Integer *>(variable_indices.get_store()),
		       new_subgradient_values.get_store()
		       );
    return ret_code;

       
  }

}

