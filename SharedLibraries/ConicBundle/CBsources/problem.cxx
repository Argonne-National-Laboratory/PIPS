/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/problem.cxx

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



#include "mymath.hxx"
#include "problem.hxx"

 
using namespace CH_Matrix_Classes;

namespace ConicBundle {

// *****************************************************************************
//                              eval_augmodel
// *****************************************************************************

//evaluate the augmented model with respect to center of stability
//if do_scaling!=0, inv_scale gives the inverse of a scaling of y,
//ie, the quadratic term must be  weightu/2*[sum_i y(i)^2/inv_scale(i)]
//if evaluated by an iterative method that provides upper and lower bounds
//it should stop when  ub-lb <= min(relprec*(fbound-lb),(ub-augbound)/2.)

int ConvexProblem::eval_augmodel(const Matrix& y,
				 const Matrix& lby,
				 const Matrix& uby,
				 const Matrix& eta,
				 Real& augbound,
				 Real& fbound,Real relprec,
				 const BundleScaling* Hp,
				 const Indexmatrix& yfixed)
{
  //initialize blocks and variable dimensions
  if (solverp==0) {
    solverp=new QP_Solver;
    solverp->set_out(out,print_level-1);    
    solverp->set_maxiter(100);
  }

  solverp->clear_blocks();
  QP_Block* blockp;
  if (start_augmodel(blockp)){
    return -1;
  }
  Integer aug_xdim=blockp->xdim();
  solverp->add_block(blockp);

  //--- get cost matrices  
  Symmatrix Q;
  Matrix c;
  Real offset;

  if (Hp->compute_QP_costs(Q,c,offset,this,aug_xdim,y,lby,uby,eta,yfixed)){
    if (out) 
      (*out)<<"*** ERROR ConvexProblem::reeval_augmodel(): update_QP_costs failed"<<std::endl;
    return -1;
  }

  //--- call the quadratic program solver

  solverp->set_termbounds(augbound,fbound);
  solverp->set_termeps(relprec);
  //std::cout<<"Q="<<Q<<" c="<<c<<" offset="<<offset<<std::endl;
  int status=solverp->solve(Q,c,offset);
  if (status){
    if (out) (*out)<<"*** WARNING: ConvexProblem::eval_augmodel(): solverp->solve() failed and returned "<<status<<std::endl;
  }

  //--- get the objective values and the aggregate subgradient
  //    if the cone the cone multiplier was increased recompute the solution
  Integer it=0;
  do {
    bool conemult_increased=false;
    bool* cp=0;
    if (it++<1) cp=&conemult_increased;
    int linstat=make_aug_linmodel(0,0,cp,&fbound);
    if (linstat){
      if (out) (*out)<<"*** WARNING: ConvexProblem::eval_augmodel(): make_aug_linmodel failed and returned "<<linstat<<std::endl;
      status|=linstat;
      break;
    }
    if (!conemult_increased) break;
    solverp->set_termbounds(augbound,fbound);
    int status=solverp->resolve();
    if (status){
      if (out) (*out)<<"*** WARNING: ConvexProblem::eval_augmodel(): solverp->solve() failed and returned "<<status<<std::endl;
    }
      
  }while (true);
 
  return status;
}


// *****************************************************************************
//                              reeval_augmodel
// *****************************************************************************

//reevaluate the augmented model for updated eta
//the values of eta-old_eta are listed in update_index and update_value
//if inv_scale.dim()==y.dim(), it gives the inverse of a scaling of y,
//ie, the quadratic term must be  weightu/2*[sum_i y(i)^2/inv_scale(i)]
//the scaling and the weight u must be the same in both augmodel-calls!
//if evaluated by an iterative method that provides upper and lower bounds
//it should stop when  ub-lb <= min(relprec*(fbound-lb),(ub-augbound)/2.)

int ConvexProblem::reeval_augmodel(const Matrix& y,
				   const Matrix& lby,
				   const Matrix& uby,
				   const Matrix& eta,
				   const Indexmatrix& update_index,
				   const Matrix& update_value,
				   Real& augbound,
				   Real& fbound,
				   Real relprec, 
				   const BundleScaling *Hp)
{

  //--- get linear cost matrix
  if (solverp==0) {
    return -1;
  }
  Integer aug_xdim=solverp->get_c().dim();

  Matrix dc;
  Real doffset;

  if (Hp->update_QP_costs(dc,doffset,this,aug_xdim,y,lby,uby,eta,update_index,update_value)){
    if (out) 
      (*out)<<"*** ERROR ConvexProblem::reeval_augmodel(): update_QP_costs failed"<<std::endl;
    return 1;
  }

  //--- call the quadratic program solver

  solverp->set_termbounds(augbound,fbound);
  solverp->set_termeps(relprec);
  int status=solverp->update(dc,doffset);
  if (status){
    if (out) (*out)<<"ConvexProblem::reeval_augmodel(): solverp->update() failed ..."<<std::endl;
  }

  //--- get the objective values and the aggregate subgradient
  //    if the cone the cone multiplier was increased recompute the solution
  Integer it=0;
  do {
    bool conemult_increased=false;
    bool* cp=0;
    if (it++<1) cp=&conemult_increased;
    int linstat=make_aug_linmodel(0,0,cp,&fbound);
    if (linstat){
      if (out) (*out)<<"*** WARNING: ConvexProblem::reeval_augmodel(): make_aug_linmodel failed and returned "<<linstat<<std::endl;
      status|=linstat;
      break;
    }
    if (!conemult_increased) break;
    solverp->set_termbounds(augbound,fbound);
    int status=solverp->resolve();
    if (status){
      if (out) (*out)<<"*** WARNING: ConvexProblem::reeval_augmodel(): solverp->solve() failed and returned "<<status<<std::endl;
    }
      
  }while (true);
  
  return status;
}
  

}

