/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/MatSDPfun.cxx

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



#include "MatSDPfun.hxx"
#include "cmsymspa.hxx"
#include "lanczpol.hxx"
#include "bigmat.hxx"
#include <queue>
#include <utility>
#include <string.h>
 

 
using namespace CH_Matrix_Classes;


namespace ConicBundle{


  BlockSDPPrimal::BlockSDPPrimal(const BlockSDPPrimal& pr): SDPPrimal() 
  {
    Xdim=pr.Xdim; 
    for(std::map<Integer,SDPPrimal*>::const_iterator i=pr.primal.begin();i!=pr.primal.end();++i){
      primal[i->first]=dynamic_cast<SDPPrimal *>(i->second->clone_primal_data());
      assert(primal[i->first]);
    }
    
  }

  BlockSDPPrimal::BlockSDPPrimal(const Indexmatrix& Xd,const std::map<Integer,SDPPrimal*>& pr): SDPPrimal()
  {
    Xdim=Xd; 
    for(std::map<Integer,SDPPrimal*>::const_iterator i=pr.begin();i!=pr.end();++i){
      if (i->second==0) continue;
      assert(i->first<Xdim.dim());
      primal[i->first]=i->second;
    }
  }
 
  BlockSDPPrimal::~BlockSDPPrimal()
  {
    for(std::map<Integer,SDPPrimal*>::iterator i=primal.begin();i!=primal.end();++i){
      delete i->second;
    }
    primal.clear();
  }
 
  /// returns a newly generated identical Object
  PrimalData* BlockSDPPrimal::clone_primal_data()
  {
    return new BlockSDPPrimal(*this);
  }

  /// this assign is only feasible for selected derivations of PrimalData 
  int BlockSDPPrimal::assign_primal_data(const PrimalData& pd)
  {
    const BlockSDPPrimal* p=dynamic_cast<const BlockSDPPrimal*>(&pd);
    if (p==0) return 1;
    Xdim=p->Xdim; 
    for(std::map<Integer,SDPPrimal*>::const_iterator i=primal.begin();i!=primal.end();++i){
      delete i->second;
    } 
    primal.clear();
    {for(std::map<Integer,SDPPrimal*>::const_iterator i=p->primal.begin();i!=p->primal.end();++i){
      primal[i->first]=dynamic_cast<SDPPrimal *>(i->second->clone_primal_data());
      assert(primal[i->first]);
    }}       
    return 0;
  }

  /// for each element aij in the support set aij=<P.row(i),P.row(j)> 
  int BlockSDPPrimal::assign_Gram_matrix(const Matrix& P)
  {
    assert(sum(Xdim)==P.rowdim());
    Integer start_row=0;
    Integer cnt=0;
    for(std::map<Integer,SDPPrimal*>::iterator i=primal.begin();i!=primal.end();++i){
      while(cnt<i->first){
	start_row+=Xdim(cnt);
	cnt++;
      }
      i->second->assign_Gram_matrix(P.rows(Range(start_row,start_row+Xdim(cnt)-1)));
    }
    return 0;
  }
    
    /** multiply this with myfactor and add itsfactor*it to this
        (it must also be a SparseSDPPrimal and on the same support) */
  int BlockSDPPrimal::aggregate_primal_data(double myfactor,double itsfactor,const PrimalData& it)
  {
    const BlockSDPPrimal* p=dynamic_cast<const BlockSDPPrimal*>(&it);
    if (p==0) return 1;
    assert(norm2(Xdim-p->Xdim)<0.1); 
    for(std::map<Integer,SDPPrimal*>::iterator i=primal.begin();i!=primal.end();++i){
      std::map<Integer,SDPPrimal*>::const_iterator j=p->primal.find(i->first);
      if ((j==p->primal.end())||(j->second==0)) return 1;
      i->second->aggregate_primal_data(myfactor,itsfactor,*j->second);
    }       
    return 0;
  }
  
  /// multiply this with myfactor and add itsfactor*P*P^T to this
  int BlockSDPPrimal::aggregate_Gram_matrix(double myfactor,double itsfactor,const Matrix& P)
  {
    assert(sum(Xdim)==P.rowdim());
    Integer start_row=0;
    Integer cnt=0;
    for(std::map<Integer,SDPPrimal*>::iterator i=primal.begin();i!=primal.end();++i){
      while(cnt<i->first){
	start_row+=Xdim(cnt);
	cnt++;
      }
      i->second->aggregate_Gram_matrix(myfactor,itsfactor,P.rows(Range(start_row,start_row+Xdim(cnt)-1)));
    }
    return 0;
  }
  
  SDPPrimal* BlockSDPPrimal::block(Integer i) const
  {
    std::map<Integer,SDPPrimal*>::const_iterator j=primal.find(i);
    if (j==primal.end()) return 0;
    return j->second;
  }
     


  
  class MaxEigSolver
  {
  private:
    
    Integer dim;
    Bigmatrix bigmat;
    bool bigmat_init;
    
    Lanczos* lanczos;
    
    /** if zero, then imprecise eigenvalue computation is allowed;
        for nonegative values the eigenvalue solver has to deliver
        at least this number of eigenvalues and eigenvectors */
    Integer exacteigs;
    
    Integer dense_limit;
    
    std::ostream* out;
    Integer print_level;
    
    
    
  public:
    void clear()
    { 
      dim=-1;
      bigmat.clear();
      bigmat_init=false;
    }
    
    MaxEigSolver()
    {

      dense_limit=50;
      out=0;
      print_level=0;
      lanczos=new Lanczpol;
      assert(lanczos);
      exacteigs=0;
      clear();
    }
    
    ~MaxEigSolver()
    { 
      delete lanczos;
    }
    
    bool is_init(){return bigmat_init;}
    
    void set_exacteigs(Integer i){ exacteigs=max(i,Integer(0));}
    
    void set_dense_limit(Integer lim){dense_limit=max(lim,Integer(0));}

    const Bigmatrix& get_bigmat() const {return bigmat;}
    
    int init(const Matrix& y,const Integer indim, const Coeffmat* C,
	     const std::map<Integer,Coeffmat*>* opA, bool dense=false)
    {
      clear();
      dim=indim;
      int status=bigmat.init(y,indim,C,opA,dense);
      if (status==0) {
	bigmat_init=true;
      }
      return status;
    }
    
    
    int evaluate_projection (const Matrix& P, 
			     const double /* relprec */,const double /* Ritz_bound */,
			     Matrix& projected_Ritz_vectors,
			     Matrix& projected_Ritz_values)
    {
      assert(bigmat_init);
      //----  compute tmpsym=P'*bigm*P   for P=bundlevecs (same span as bundle)
      bigmat.lanczosmult(P,projected_Ritz_vectors);
      Integer k=P.coldim();
      Symmatrix tmpsym; tmpsym.newsize(k);
      chk_set_init(tmpsym,1);
      tmpsym.xetriu_yza(P,projected_Ritz_vectors);
      
      //---- determine maximum eigenvalue of tmpsym
      tmpsym.eig(projected_Ritz_vectors,projected_Ritz_values,false);
      if ((out)&&(print_level>0)){
	(*out).precision(8);
	(*out)<<"  eigmod="<<max(projected_Ritz_values);
      }
      return 0;
    }
    
    int compute_projection (const Matrix& P, Symmatrix& S)
    {
      assert(bigmat_init);
      //----  compute tmpsym=P'*bigm*P   for P=bundlevecs (same span as bundle)
      Matrix tmpmat;
      bigmat.lanczosmult(P,tmpmat);
      Integer k=P.coldim();
      S.newsize(k);
      chk_set_init(S,1);
      S.xetriu_yza(P,tmpmat);
      
      return 0;
    }
    
    int evaluate(const Matrix& bundlevecs, 
		 const double relprec, const double Ritz_bound,
		 Matrix& Ritz_vectors, Matrix&  Ritz_values)
    {
      assert(bigmat_init);
      /* for testing:
      (bigmat.get_symrep()).eig(Ritz_vectors,Ritz_values,false);
      std::cout<<"eval: lmax="<<Ritz_values(0)<<std::endl;
      */
      int status=0;
      Integer nreig=0;
      if (bigmat.lanczosdim()<dense_limit){
	(bigmat.get_symrep()).eig(Ritz_vectors,Ritz_values,false);
	nreig=dim;
	if ((out)&&(print_level>0)){
	  (*out)<<" dense: ";
	  (*out)<<" eigmax=";out->precision(8);(*out)<<Ritz_values(0)<<std::endl;
	}
      }
      //-------------------- sparse routine for lmax computation
      else {
	lanczos->set_relprec(relprec);
	Integer blocksz=1;
	if (exacteigs==0) {
          /*
	  if (relprec<1e-7){
	    lanczos->enable_stop_above(CB_plus_infinity);
	    nreig=5;
	    blocksz=5;
	  }
	  else{
	  */
	  lanczos->enable_stop_above(Ritz_bound);
	  nreig=1;
	  /* } */
	}
	else {
	  lanczos->enable_stop_above(CB_plus_infinity);
	  nreig=exacteigs;
	}
	
	bigmat.reset_nmult();
	
	
	Ritz_values.init(0,0,0.);
        if (bundlevecs.dim()!=0) {
	  Ritz_vectors=bundlevecs;
	}
	else {
	  Ritz_vectors.init(0,0,0.);
	}
	status=lanczos->compute(&bigmat,Ritz_values,Ritz_vectors,nreig,blocksz);
       

	/*
	Ritz_values.init(0,0,0.);
	if ((Ritz_vectors.rowdim()!=bigmat.lanczosdim())||(Ritz_vectors.coldim()==0))
	  Ritz_vectors.init(0,0,0.);
	Integer blocksz=1;
	status=lanczos->compute(&bigmat,Ritz_values,Ritz_vectors,nreig,blocksz);
	*/
	  

	nreig=Ritz_values.dim();
	if ((out)&&(print_level>0)){
	  (*out)<<"  Lanczos {"<<lanczos->get_iter()<<",";
	  (*out)<<bigmat.get_nmult()<<","<<nreig<<"} ";
	}
	
	lanczos->get_lanczosvecs(Ritz_values,Ritz_vectors);
	
	if ((out)&&(print_level>0)){
	  if (nreig>0) {(*out)<<" eigmax=";out->precision(8);(*out)<<Ritz_values(0);}
	  else {(*out)<<" Ritz_val=";out->precision(8);(*out)<<Ritz_values(0);}
	  (*out)<<" rp=";out->precision(2);(*out)<<relprec;
	  (*out)<<" vecs="<<Ritz_values.dim();
	  out->precision(8);(*out)<<" bd="<<Ritz_bound<<std::endl;
	}
      }
      return status;
    }
    
    
    void  set_out(std::ostream* o=0,int pril=1)
    {out=o;print_level=pril;if (lanczos) lanczos->set_out(out,pril-1);}
    
  };  //end of class MaxEigSolver
  
  
  
  
  
  
  




    
  MatrixSDPfunction::MatrixSDPfunction()
  {
    trace_value=-1;
    constant_trace=SDPtrace_fixed;
    generating_primal=0;
    maxvecs=5;
    
    out=0;
    print_level=0;
  }  
  
  MatrixSDPfunction::~MatrixSDPfunction()
  {
    for(std::map<Integer,SparseSDPCoeffVector*>::iterator mapj=opA_colrep.begin();
	mapj!=opA_colrep.end();++mapj){
      SparseSDPCoeffVector* col=(*mapj).second;
      
      for(std::map<Integer,Coeffmat*>::iterator mapi=col->begin();
	  mapi!=col->end();++mapi){
	delete (*mapi).second;
      }
      
      delete col;
    }
    opA_colrep.clear();
    for(std::map<Integer,SparseSDPCoeffVector*>::iterator mapi=opA_rowrep.begin();
	mapi!=opA_rowrep.end();++mapi){
      delete (*mapi).second;
    }
    opA_rowrep.clear();

    {for(SparseSDPCoeffVector::iterator mapj=C.begin();mapj!=C.end();++mapj){
      delete (*mapj).second;
    }}
    C.clear();
     
    for(unsigned int i=0;i<maxeigsolver.size();++i){
      delete maxeigsolver[i];
    }
    delete generating_primal;
  }
  
  void MatrixSDPfunction::get_trace_constraint(Real& tv,SDPtrace& ct)
  {
    tv=trace_value;
    ct=constant_trace;
  }
 
  int MatrixSDPfunction::gramip(const Matrix& P,Real& ip_C,Matrix& ip_opA)
  {

    if (Xdim.dim()==1){
      if (C.begin()!=C.end()) {
	ip_C=C.begin()->second->gramip(P);
      }
      else ip_C=0.;
      gramip_opA(P,ip_opA);
      //if (out) {
      //(*out)<<"P="; P.display(*out);
      //(*out)<<"ip_C="<<ip_C<<" ip_opA=";ip_opA.display(*out);
      //}
      return 0;
    }
    ip_C=0.;
    ip_opA.init(b.dim(),1,0.);
    SparseSDPCoeffVector::iterator cit=C.begin();
    SparseSDPCoeffMatrix::iterator coli=opA_colrep.begin();
    Integer start=0;
    for(Integer i=0;i<Xdim.dim();i++){
      if ((cit!=C.end())&&(cit->first==i)) {
	ip_C+=cit->second->gramip(P,start);
	++cit;
      }

      if ((coli!=opA_colrep.end())&&(coli->first==i)){
	SparseSDPCoeffVector& col=*coli->second;
	for(SparseSDPCoeffVector::iterator rowi=col.begin();rowi!=col.end(); ++rowi){
	  ip_opA(rowi->first)+=rowi->second->gramip(P,start);
	}
	++coli;
      }
      start+=Xdim(i);
    }
    return 0;
  }
  
  int MatrixSDPfunction::gramip_opA(const Matrix& P,Matrix& vec)
  {
    if (Xdim.dim()==1){
      vec.init(b.dim(),1,0.);
      if (opA_colrep.begin()==opA_colrep.end()) return 0;
      SparseSDPCoeffVector& col=*opA_colrep.begin()->second;
      for(SparseSDPCoeffVector::iterator rowi=col.begin();rowi!=col.end(); ++rowi){
	vec(rowi->first)+=rowi->second->gramip(P);
      }
      // if (out) {
      //(*out)<<"P="; P.display(*out);
      //(*out)<<" ip_opA=";vec.display(*out);
      //}
      return 0;
    }
    vec.init(b.dim(),1,0.);
    SparseSDPCoeffMatrix::iterator coli=opA_colrep.begin();
    Integer start=0;
    for(Integer i=0;i<Xdim.dim();i++){
      if ((coli!=opA_colrep.end())&&(coli->first==i)){
	SparseSDPCoeffVector& col=*coli->second;
	for(SparseSDPCoeffVector::iterator rowi=col.begin();rowi!=col.end(); ++rowi){
	  vec(rowi->first)+=rowi->second->gramip(P,start);
	}
	++coli;
      }
      start+=Xdim(i);
    }
    return 0;
  }
  
  
  const Matrix& MatrixSDPfunction::rhs() const
  {
    return b;
  }
  
  int MatrixSDPfunction::project_C(const Matrix& P,Symmatrix& S)
  {
    if (Xdim.dim()==1){
      if (C.begin()!=C.end()) {
	C.begin()->second->project(S,P);
	return 0;
      }
      S.init(P.coldim(),0.);
      return 0;
    }
    S.init(P.coldim(),0.);
    Integer start=0;
    Integer i=0;
    for (SparseSDPCoeffVector::iterator cit=C.begin();cit!=C.end();++cit){
      while(i<cit->first){
	start+=Xdim(i);
	++i;
      }
      cit->second->add_projection(S,P,start);
    }
    return 0;
  }

  int MatrixSDPfunction::project(const Integer i, const Matrix& P,Symmatrix& S)
  {

    SparseSDPCoeffMatrix::iterator rowi=opA_rowrep.find(i);
    if (rowi==opA_rowrep.end()){ //no nonzero matrix in row i
      S.init(P.coldim(),0.);
      return  0;
    }
    SparseSDPCoeffVector& row=*rowi->second;
    if (Xdim.dim()==1){
      row.begin()->second->project(S,P);
      return 0;
    }
    S.init(P.coldim(),0.);
    Integer start=0;
    Integer j=0;
    for (SparseSDPCoeffVector::iterator coli=row.begin();coli!=row.end();++coli){
      while(j<coli->first){
	start+=Xdim(j);
	++j;
      }
      coli->second->add_projection(S,P,start);
    }
    return 0;
  }

  SDPPrimal* MatrixSDPfunction::init_primal(const Matrix& P)
  {
    if (generating_primal){
      PrimalData* sp=generating_primal->clone_primal_data();
      SDPPrimal* sdpp=dynamic_cast<SDPPrimal*>(sp);
      if (sdpp==0) return 0; 
      sdpp->assign_Gram_matrix(P);
      return sdpp;
    }
    return 0;
  }
    
  int MatrixSDPfunction::primalip_opA(const SDPPrimal* primal,
				      const Indexmatrix& ind,
				      Matrix& ip_opA)
  {
    if (primal==0){ //is ok if opA has only zero rows corresponding to ind
      ip_opA.init(ind.dim(),1,0.);
      for(Integer i=0;i<ind.dim();i++){
	if (opA_rowrep.find(ind(i))!=opA_rowrep.end()) return 1;
      }
      return 0;
    }
    //--- check for dense symmetric primal
    const DenseSDPPrimal* dp=dynamic_cast<const DenseSDPPrimal*>(primal);
    if (dp){
      if (Xdim.dim()!=1) return 1;
      ip_opA.init(ind.dim(),1,0.);
      if (opA_colrep.begin()==opA_colrep.end()) return 0;
      SparseSDPCoeffVector& col=*(opA_colrep.begin()->second);
      for (Integer i=0;i<ind.dim();i++){
	SparseSDPCoeffVector::iterator rowi=col.find(ind(i));
        if (rowi==col.end()) continue;
	ip_opA(i)=rowi->second->ip(*dp);
      }
      return 0;
    } 
    //--- check for sparse symmetric primal
    const SparseSDPPrimal* sp=dynamic_cast<const SparseSDPPrimal*>(primal);
    if (sp){
      if (Xdim.dim()!=1) return 1;
      ip_opA.init(ind.dim(),1,0.);
      if (opA_colrep.begin()==opA_colrep.end()) return 0;
      SparseSDPCoeffVector& col=*(opA_colrep.begin()->second);
      for (Integer i=0;i<ind.dim();i++){
	SparseSDPCoeffVector::iterator rowi=col.find(ind(i));
	if (rowi==col.end()) continue;
	if (rowi->second->support_in(*sp)==0){
	  return 1;
	}
	ip_opA(i)=rowi->second->ip(*sp);
      }
      return 0;
    }
    //--- check for sparse symmetric primal
    const GramSparseSDPPrimal* gsp=dynamic_cast<const GramSparseSDPPrimal*>(primal);
    if (gsp){
      if (Xdim.dim()!=1) return 1;
      ip_opA.init(ind.dim(),1,0.);
      if (opA_colrep.begin()==opA_colrep.end()) return 0;
      SparseSDPCoeffVector& col=*(opA_colrep.begin()->second);
      for (Integer i=0;i<ind.dim();i++){
	SparseSDPCoeffVector::iterator rowi=col.find(ind(i));
	if (rowi==col.end()) continue;
	if (rowi->second->support_in(*gsp)==0){
	  return 1;
	}
	ip_opA(i)=rowi->second->ip(*gsp);
	if (gsp->get_grammatrix().dim()>0)
	  ip_opA(i)+=rowi->second->gramip(gsp->get_grammatrix());
      }
      return 0;
    }
    //--- check for block primal (must be of the same structure as opA)
    const BlockSDPPrimal* bp=dynamic_cast<const BlockSDPPrimal*>(primal);
    if (bp){
      if (bp->get_nblocks()!=Xdim.dim()) return 1;
      ip_opA.init(ind.dim(),1,0.);
      SparseSDPCoeffMatrix::iterator coli=opA_colrep.begin();
      for (Integer j=0;(coli!=opA_colrep.end())&&(j<Xdim.dim());j++){
	if (bp->block_dim(j)!=Xdim(j)) return 1;
	if (coli->first!=j) continue;
	SparseSDPCoeffVector& col=*(coli->second);
	coli++;
	const DenseSDPPrimal* dp=dynamic_cast<const DenseSDPPrimal*>(bp->block(j));
	if (dp){
	  for (Integer i=0;i<ind.dim();i++){
	    SparseSDPCoeffVector::iterator rowi=col.find(ind(i));
	    if (rowi==col.end()) continue;
	    ip_opA(i)+=rowi->second->ip(*dp);
	  }
	  continue;
	}
	const SparseSDPPrimal* sp=dynamic_cast<const SparseSDPPrimal*>(bp->block(j));
	if (sp){
	  for (Integer i=0;i<ind.dim();i++){
	    SparseSDPCoeffVector::iterator rowi=col.find(ind(i));
	    if (rowi==col.end()) continue;
	    if (rowi->second->support_in(*sp)==0){
	      return 1;
	    }
	    ip_opA(i)+=rowi->second->ip(*sp);
	  }
	  continue;
	}
	return 1;
      }
      return 0;
    }
    //--- type not recognized
    return 1;
  }
  
  int MatrixSDPfunction::primalip_C(const SDPPrimal* primal,Real& value)
  {
    if (primal==0){ //is ok if C is zero
      value=0.;
      if (C.begin()==C.end()) return 1;
      return 0;
    }
    //--- check for dense symmetric primal
    const DenseSDPPrimal* dp=dynamic_cast<const DenseSDPPrimal*>(primal);
    if (dp){
      if (Xdim.dim()!=1) return 1;
      value=0.;
      if (C.begin()==C.end()) return 0;
      value=C.begin()->second->ip(*dp);
      return 0;
    } 
    //--- check for sparse symmetric primal
    const SparseSDPPrimal* sp=dynamic_cast<const SparseSDPPrimal*>(primal);
    if (sp){
      if (Xdim.dim()!=1) return 1;
      value=0.;
      if (C.begin()==C.end()) return 0;
      if (C.begin()->second->support_in(*sp)==0) return 1;
      value=C.begin()->second->ip(*sp);
      return 0;
    }
    //--- check for block primal (must be of the same structure as opA)
    const BlockSDPPrimal* bp=dynamic_cast<const BlockSDPPrimal*>(primal);
    if (bp){
      if (bp->get_nblocks()!=Xdim.dim()) return 1;
      value=0.;
      SparseSDPCoeffVector::iterator coli=C.begin();
      for (Integer j=0;(coli!=C.end())&&(j<Xdim.dim());j++){
	if (bp->block_dim(j)!=Xdim(j)) return 1;
	if (coli->first!=j) continue;
	const DenseSDPPrimal* dp=dynamic_cast<const DenseSDPPrimal*>(bp->block(j));
	if (dp){
	  value+=coli->second->ip(*dp);
	  coli++;
	  continue;
	}
	const SparseSDPPrimal* sp=dynamic_cast<const SparseSDPPrimal*>(bp->block(j));
	if (sp){
	  if (coli->second->support_in(*sp)==0) return 1;
	  value+=coli->second->ip(*sp);
	  coli++;
	  continue;
	}
	return 1;
      }
      return 0;
    }
    //--- type not recognized
    return 1;
  }
  
  
  int MatrixSDPfunction::evaluate(const Matrix& current_point,
				  const Matrix& bundlevecs, 
				  const double relprec,
				  const double Ritz_bound,
				  Matrix& Ritz_vectors,
				  Matrix&  Ritz_values)
  {
    if (current_point.dim()!=b.dim()){
      if (out){
	(*out)<<"**** ERROR: MatrixSDPfunction::evaluate(...): argument dimension "<<current_point.dim()<<" does not equal number of constraints "<<b.dim()<<std::endl;
      }
      return 1;
    }
    
    //initialize the maxeigsolver if this is necessary  
    bool same_y=(last_bigmat_y.coldim()==1);
    if ((same_y)&&(current_point.dim()==last_bigmat_y.dim())){
      for(Integer i=0;i<current_point.dim();++i){      
	if (fabs(current_point(i)-last_bigmat_y(i))>eps_Real*max(fabs(last_bigmat_y(i)),1.)) {
	  same_y=false;
	  break;
	}
      }
    }
    else {
      same_y=false;
    }
    bool must_init=!same_y;
    if (!must_init){
      for (Integer i=0;i<Xdim.dim();++i){
	if (!maxeigsolver[i]->is_init()){
	  must_init=true;
	  break;
	}
      }
    }
    
    if (must_init){
      SparseSDPCoeffMatrix::iterator coli=opA_colrep.begin();
      SparseSDPCoeffVector::iterator cit=C.begin();
      for (Integer i=0;i<Xdim.dim();++i){
	if ((same_y)&&(maxeigsolver[i]->is_init())) continue;
	Coeffmat *cp;
	if ((cit!=C.end())&&(cit->first==i)){
	  cp=cit->second;
	  ++cit;
	}
	else {
	  cp=0;
	}
	SparseSDPCoeffVector* opAp;
	if ((coli!=opA_colrep.end())&&(coli->first==i)){
	  opAp=coli->second;
	  ++coli;
	}
	else {
	  opAp=0;
	}
	if (maxeigsolver[i]->init(current_point,Xdim(i),cp,opAp,(n_dense_coeffmat[i]!=0))){
	  if (out) {
	    (*out)<<"**** ERROR: MatrixSDPfunction::evaluate(...): initialization of Eigenvaluesolver failed for matrix "<<i<<std::endl;
	  }
	  return 1;
	}
      }
      last_bigmat_y=current_point;
    } 
    
    //evaluate
    int retval=0;
    if (Xdim.dim()==1){
      retval=maxeigsolver[0]->evaluate(bundlevecs,relprec,Ritz_bound,
				       Ritz_vectors,Ritz_values);
      if (retval){
	if (out) {
	  (*out)<<"**** ERROR: MatrixSDPfunction::evaluate(...): Eigenvaluesolver failed with code "<<retval<<std::endl;
	}
      }
    }
    else {
      typedef std::pair<Real,Integer> Sol;
      std::priority_queue<Sol> priq;
      Integer start_row=0;
      Integer dim=sum(Xdim);
      Matrix tmp_sol_vecs(dim,maxvecs); chk_set_init(tmp_sol_vecs,1);
      Matrix tmp_sol_vals(maxvecs,1); chk_set_init(tmp_sol_vals,1);
      for (Integer i=0;i<Xdim.dim();++i){
	Indexmatrix tmp_ind(Range(start_row,start_row+Xdim(i)-1));
	Matrix tmp_bundle;
        if (bundlevecs.rowdim()>=dim){
	  tmp_bundle=bundlevecs.rows(tmp_ind);
	}
	Matrix tmp_Ritz_vec;
        if (Ritz_vectors.rowdim()>=dim){
	  tmp_Ritz_vec=Ritz_vectors.rows(tmp_ind);
	}
	Matrix tmp_Ritz_val;
	int lretval=maxeigsolver[i]->evaluate(tmp_bundle,relprec,Ritz_bound,
					      tmp_Ritz_vec,tmp_Ritz_val);
	if (lretval){
	  retval|=lretval;
	  if (out) {
	    (*out)<<"**** ERROR: MatrixSDPfunction::evaluate(...): Eigenvaluesolver failed for matrix "<<i<<" with code"<<lretval<<std::endl;
	  }
	}
	else {
	  Indexmatrix ind;
	  sortindex(-tmp_Ritz_val,ind);
	  Integer j=0;
	  while(j<ind.dim()){
	    Integer jj=ind(j);
	    j++;
	    Real val=tmp_Ritz_val(jj);
	    Integer ii;
	    if (Integer(priq.size())<maxvecs){
	      ii=Integer(priq.size());
	    }
	    else if (-priq.top().first<val){
	      ii=priq.top().second;
	      priq.pop();
	    }
	    else {
	      break;
	    }
	    mat_xea(start_row,tmp_sol_vecs.get_store()+ii*dim,0.);
	    mat_xey(Xdim(i),tmp_sol_vecs.get_store()+ii*dim+start_row,tmp_Ritz_vec.get_store()+jj*Xdim(i));
	    mat_xea(dim-Xdim(i)-start_row,tmp_sol_vecs.get_store()+ii*dim+start_row+Xdim(i),0.);
	    tmp_sol_vals(ii)=val;
	    priq.push(Sol(-val,ii));
	  }  
	}
	start_row+=Xdim(i);
      } //end for
      tmp_sol_vecs.delete_cols(Range(Integer(priq.size()),maxvecs-1));
      tmp_sol_vals.delete_rows(Range(Integer(priq.size()),maxvecs-1));      
      Indexmatrix sind;
      sortindex(-tmp_sol_vals,sind);
      Ritz_vectors=tmp_sol_vecs.cols(sind);
      Ritz_values=tmp_sol_vals.rows(sind);  
      if ((out)&&(print_level>1)){
	(*out)<<" lmax="; out->precision(6);
	for(Integer i=0;i<min(Integer(5),Ritz_values.dim());i++){
	  (*out)<<" "<<Ritz_values(i);
	}
	(*out)<<std::endl;
      }
    }
    
    return retval;
  }

  
  int MatrixSDPfunction::evaluate_projection(const Matrix& current_point,
					     const Matrix& P, 
					     const double relprec,
					     const double Ritz_bound,
					     Matrix& projected_Ritz_vectors,
					     Matrix& projected_Ritz_values)
  {
    if (current_point.dim()!=b.dim()){
      if (out){
	(*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): argument dimension "<<current_point.dim()<<" does not equal number of constraints "<<b.dim()<<std::endl;
      }
      return 1;
    }
    

    //initialize the maxeigsolver if this is necessary  
    bool same_y=(last_bigmat_y.coldim()==1);
    if ((same_y)&&(current_point.dim()==last_bigmat_y.dim())){
      for(Integer i=0;i<current_point.dim();++i){
	if (fabs(current_point(i)-last_bigmat_y(i))>eps_Real*max(fabs(last_bigmat_y(i)),1.)) {
	  same_y=false;
	  break;
	}
      }
    }
    else {
      same_y=false;
    }
    bool must_init=!same_y;
    if (!must_init){
      for (Integer i=0;i<Xdim.dim();++i){
	if (!maxeigsolver[i]->is_init()){
	  must_init=true;
	  break;
	}
      }
    }
    
    if (must_init){
      SparseSDPCoeffMatrix::iterator coli=opA_colrep.begin();
      SparseSDPCoeffVector::iterator cit=C.begin();
      for (Integer i=0;i<Xdim.dim();++i){
	if ((same_y)&&(maxeigsolver[i]->is_init())) continue;
	Coeffmat *cp;
	if ((cit!=C.end())&&(cit->first==i)){
	  cp=cit->second;
	  ++cit;
	}
	else {
	  cp=0;
	}
	SparseSDPCoeffVector* opAp;
	if ((coli!=opA_colrep.end())&&(coli->first==i)){
	  opAp=coli->second;
	  ++coli;
	}
	else {
	  opAp=0;
	}
	if (maxeigsolver[i]->init(current_point,Xdim(i),cp,opAp,(n_dense_coeffmat[i]!=0))){
	  if (out) {
	    (*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): initialization of Eigenvaluesolver failed for matrix "<<i<<std::endl;
	  }
	  return 1;
	}
      }
      last_bigmat_y=current_point;
    } 
    
    
    
    //evaluate
    int retval=0;
    if (Xdim.dim()==1){
      retval=maxeigsolver[0]->evaluate_projection(P,relprec,Ritz_bound,
				    projected_Ritz_vectors,projected_Ritz_values);
      if (retval){
	if (out) {
	  (*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): Eigenvaluesolver failed with code "<<retval<<std::endl;
	}
      }
    }
    else {
      Integer start_row=0;
      Real maxval=CB_minus_infinity;
      Matrix tmp_sol_vecs;
      Matrix tmp_sol_vals;
      for (Integer i=0;i<Xdim.dim();++i){
	Indexmatrix tmp_ind(Range(start_row,start_row+Xdim(i)-1));
	Matrix tmp_P=P.rows(tmp_ind);
	Matrix tmp_Ritz_vec=projected_Ritz_vectors;
	Matrix tmp_Ritz_val=projected_Ritz_values;
	int lretval=maxeigsolver[i]->evaluate_projection(tmp_P,relprec,Ritz_bound,
					   tmp_Ritz_vec,tmp_Ritz_val);
	if (lretval){
	  retval|=lretval;
	  if (out) {
	    (*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): Eigenvaluesolver failed for matrix "<<i<<" with code"<<lretval<<std::endl;
	  }
	  continue;
	}
	if (maxval<max(tmp_Ritz_val)){
	  maxval=max(tmp_Ritz_val);
	  tmp_sol_vecs=tmp_Ritz_vec;
	  tmp_sol_vals=tmp_Ritz_val;
	}
      } //end for
      
      projected_Ritz_vectors=tmp_sol_vecs;
      projected_Ritz_values=tmp_sol_vals;      
    }
    
    return retval;
  }


  int MatrixSDPfunction::compute_projection(const Matrix& current_point,
					     const Matrix& P, 
					     Symmatrix& S)
  {
    if (current_point.dim()!=b.dim()){
      if (out){
	(*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): argument dimension "<<current_point.dim()<<" does not equal number of constraints "<<b.dim()<<std::endl;
      }
      return 1;
    }
    

    //initialize the maxeigsolver if this is necessary  
    bool same_y=(last_bigmat_y.coldim()==1);
    if ((same_y)&&(current_point.dim()==last_bigmat_y.dim())){
      for(Integer i=0;i<current_point.dim();++i){
	if (fabs(current_point(i)-last_bigmat_y(i))>eps_Real*max(fabs(last_bigmat_y(i)),1.)) {
	  same_y=false;
	  break;
	}
      }
    }
    else {
      same_y=false;
    }
    bool must_init=!same_y;
    if (!must_init){
      for (Integer i=0;i<Xdim.dim();++i){
	if (!maxeigsolver[i]->is_init()){
	  must_init=true;
	  break;
	}
      }
    }
    
    if (must_init){
      SparseSDPCoeffMatrix::iterator coli=opA_colrep.begin();
      SparseSDPCoeffVector::iterator cit=C.begin();
      for (Integer i=0;i<Xdim.dim();++i){
	if ((same_y)&&(maxeigsolver[i]->is_init())) continue;
	Coeffmat *cp;
	if ((cit!=C.end())&&(cit->first==i)){
	  cp=cit->second;
	  ++cit;
	}
	else {
	  cp=0;
	}
	SparseSDPCoeffVector* opAp;
	if ((coli!=opA_colrep.end())&&(coli->first==i)){
	  opAp=coli->second;
	  ++coli;
	}
	else {
	  opAp=0;
	}
	if (maxeigsolver[i]->init(current_point,Xdim(i),cp,opAp,(n_dense_coeffmat[i]!=0))){
	  if (out) {
	    (*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): initialization of Eigenvaluesolver failed for matrix "<<i<<std::endl;
	  }
	  return 1;
	}
      }
      last_bigmat_y=current_point;
    } 
    
    
    
    //evaluate
    int retval=0;
    if (Xdim.dim()==1){
      retval=maxeigsolver[0]->compute_projection(P,S);
      if (retval){
	if (out) {
	  (*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): Eigenvaluesolver failed with code "<<retval<<std::endl;
	}
      }
    }
    else {
      Integer start_row=0;
      S.init(P.coldim(),0.);
      for (Integer i=0;i<Xdim.dim();++i){
	Indexmatrix tmp_ind(Range(start_row,start_row+Xdim(i)-1));
	Matrix tmp_P=P.rows(tmp_ind);
	Symmatrix tmp_S;
	int lretval=maxeigsolver[i]->compute_projection(tmp_P,tmp_S);
	if (lretval){
	  retval|=lretval;
	  if (out) {
	    (*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): Eigenvaluesolver failed for matrix "<<i<<" with code"<<lretval<<std::endl;
	  }
	  continue;
	}
	S+=tmp_S;
      } //end for
    }
    
    return retval;
  }

  int MatrixSDPfunction::eig(const Matrix& current_point,
			     Matrix& P, 
			     Matrix& d)
  {
    if (current_point.dim()!=b.dim()){
      if (out){
	(*out)<<"**** ERROR: MatrixSDPfunction::evaluate_projection(...): argument dimension "<<current_point.dim()<<" does not equal number of constraints "<<b.dim()<<std::endl;
      }
      return 1;
    }
    
    Integer err=0;
    if (Xdim.dim()==1){
      const Bigmatrix& B=get_bigmat(current_point,0);
      const Symmatrix& symrep=B.get_symrep();
      if (symrep.eig(P,d,false))
	err=1;
    }
    else {
      Integer n=sum(Xdim);
      P.init(n,n,0.);
      d.init(n,1,0.);
      if (n==0) return 0;
      Integer r=0;
      for (Integer i=0;i<Xdim.dim();i++){
	const Bigmatrix& B=get_bigmat(current_point,i);
	const Symmatrix& symrep=B.get_symrep();
	Matrix tmpP,tmpd;
	if (symrep.eig(tmpP,tmpd,false))
	  err=i+1;
	Integer tmpn=tmpd.dim();
	for (Integer j=0;tmpn;j++){
	  d(j+r)=tmpd(j);
	  for(Integer k=0;tmpn;k++){
	    P(k+r,j+r)=tmpP(k,j);
	  }
	}
	r+=tmpn;   
      }
      Indexmatrix sind=sortindex(-d);
      d=d(sind);
      P=P.cols(sind);
    }

    return err;
  }


  //------------------  routines for full size computations ----------------
  
  /// computes the eigenvectors P and eigenvalues d of block ind of C-opAt(y) 
  int MatrixSDPfunction::eig(const Matrix& y,Integer ind,Matrix& P,Matrix& d) const
  {
    chk_range(ind,0,Xdim.dim(),1);
    SparseSDPCoeffMatrix::const_iterator coli=opA_colrep.find(ind);
    SparseSDPCoeffVector::const_iterator cit=C.find(ind);
    const Coeffmat *cp;
    if (cit!=C.end()){
      cp=cit->second;
    }
    else {
      cp=0;
    }
    const SparseSDPCoeffVector* opAp;
    if (coli!=opA_colrep.end()){
	opAp=coli->second;
    }
    else {
      opAp=0;
    }
    MaxEigSolver mes;
    if (mes.init(y,Xdim(ind),cp,opAp,(n_dense_coeffmat[ind]!=0))){
      if (out) {
	(*out)<<"**** ERROR: MatrixSDPfunction::eig(...): initialization of Eigenvaluesolver failed for matrix "<<ind<<std::endl;
      }
      return 1;
    }
    Matrix dummy_bundle;
    double relprec=1e-10;
    double Ritz_bound=1e20;
    mes.set_dense_limit(Xdim(ind)+1);
    if (mes.evaluate(dummy_bundle,relprec,Ritz_bound,P,d)){
      if (out) {
	(*out)<<"**** ERROR: MatrixSDPfunction::eig(...): evaluate of Eigenvaluesolver failed for matrix "<<ind<<std::endl;
      }
      return 1;
    }
      
    return 0;
  }
  
  /// returns the internal Bigmatrix representation of block ind of C-opAt(y) 
  const Bigmatrix& MatrixSDPfunction::get_bigmat(const CH_Matrix_Classes::Matrix& current_point,CH_Matrix_Classes::Integer ind)
  {
    chk_range(ind,0,Xdim.dim(),1);
    if (current_point.dim()!=b.dim()){
      if (out){
	(*out)<<"**** ERROR: MatrixSDPfunction::get_bigmat(...): argument dimension "<<current_point.dim()<<" does not equal number of constraints "<<b.dim()<<std::endl;
      }
      MEmessage(MatrixError(ME_unspec," exiting from MatrixSDPfunction::get_bigmat(...)",MTglobalfun));
      return maxeigsolver[ind]->get_bigmat();
    }
    
    //initialize the maxeigsolver if this is necessary  
    bool same_y=(last_bigmat_y.coldim()==1);
    if ((same_y)&&(current_point.dim()==last_bigmat_y.dim())){
      for(Integer i=0;i<current_point.dim();++i){      
	if (fabs(current_point(i)-last_bigmat_y(i))>eps_Real*max(fabs(last_bigmat_y(i)),1.)) {
	  same_y=false;
	  break;
	}
      }
    }
    else {
      same_y=false;
    }
    bool must_init=!same_y;
    if ((!must_init)&& (!maxeigsolver[ind]->is_init())){
      must_init=true;
    }
  
    if (must_init){
      SparseSDPCoeffMatrix::const_iterator coli=opA_colrep.find(ind);
      SparseSDPCoeffVector::const_iterator cit=C.find(ind);
      const Coeffmat *cp;
      if (cit!=C.end()){
	cp=cit->second;
      }
      else {
	cp=0;
      }
      const SparseSDPCoeffVector* opAp;
      if (coli!=opA_colrep.end()){
	opAp=coli->second;
      }
      else {
	opAp=0;
      }
      if (maxeigsolver[ind]->init(current_point,Xdim(ind),cp,opAp,(n_dense_coeffmat[ind]!=0))){
	if (out) {
	  (*out)<<"**** ERROR: MatrixSDPfunction::get_bigmat(...): initialization of Eigenvaluesolver failed for matrix "<<ind<<std::endl;
	}
	MEmessage(MatrixError(ME_unspec," exiting from MatrixSDPfunction::get_bigmat(...)",MTglobalfun));
	return maxeigsolver[ind]->get_bigmat();
      }
      last_bigmat_y=current_point;
    }
    return maxeigsolver[ind]->get_bigmat();
  }
  
 

  //------------------  routines for modifying the problem ----------------
  
  /// set the basic constraint <I,X> <=/= trace; #trace# must be nonnegative
  int MatrixSDPfunction::set_trace(const Real tr, SDPtrace constant_tr)
  {
    if (tr<0.) return 1;
    trace_value=tr;
    constant_trace=constant_tr;
    return 0;
  }
  
  
  /** defines in what form primal matrices should be aggregated.
      If the argument is NULL then no primal aggregation will take place.
      The control over the generating primal is
      passed over to this. This will delete an existing generating primal 
      whenever a new generating primal is set or upon destruction. */
  int MatrixSDPfunction::set_generating_primal(SDPPrimal* gen_primal)
  {
    delete generating_primal;
    generating_primal=gen_primal;
    return 0;
  }
  
  /** add new variables and their coefficients for existing constraints.
      It is not possible to introduce new constraints by this routine
      The control over the coefficient matrices pointed to by C and opA
      is passed over to this. This will delete C whenever a new
      cost matrix is set or upon destruction. */
  int MatrixSDPfunction::append_variables(const Indexmatrix append_Xdim,
		       const SparseSDPCoeffVector& append_C,
		       const std::map<Integer,SparseSDPCoeffVector>& opA_columns)
  {
    if (append_Xdim.dim()==0) return 0;
    assert(min(append_Xdim)>0);
    last_bigmat_y.init(0,0,0.);
    Integer add_dim=append_Xdim.dim();
    Integer old_dim=Xdim.dim();
    n_dense_coeffmat.concat_below(Indexmatrix(add_dim,1,Integer(0)));
    
    for(SparseSDPCoeffVector::const_iterator cit=append_C.begin();
	cit!=append_C.end();++cit){ 
      if ((cit->first<0)||(cit->first>=add_dim)){
	if (out){
	  (*out)<<"**** ERROR in MatrixSDPfunction::append_variables(...): index "<<cit->first<<" of a column of append_C exceeds number of added variables="<<add_dim<<std::endl;
	}
	return 1;
      }
      if (cit->second==0){
	if (out){
	  (*out)<<"**** ERROR in MatrixSDPfunction::append_variables(...): cost matrix of column "<<cit->first<<" holds a null pointer"<<std::endl;
	}
	return 1;
      }
      if (cit->second->dim()!=append_Xdim(cit->first)){
	if (out){
	  (*out)<<"**** ERROR in MatrixSDPfunction::append_variables(...): for matrix "<<cit->first<<" dimension of variable ("<<append_Xdim(cit->first)<<") and cost matrix ("<<cit->second->dim()<<") do not match"<<std::endl;
	}
	return 1;
      }	
      C[old_dim+cit->first]=cit->second;
      if (cit->second->dense()) n_dense_coeffmat(old_dim+cit->first)++;
    }

    for(std::map<Integer,SparseSDPCoeffVector>::const_iterator coli=opA_columns.begin();
	coli!=opA_columns.end();++coli){
      if ((coli->first<0)||(coli->first>=add_dim)){
	if (out){
	  (*out)<<"**** ERROR in MatrixSDPfunction::append_variables(...): index "<<coli->first<<" of a column of opA_columns exceeds number of added variables="<<add_dim<<std::endl;
	}
	return 1;
      }
      SparseSDPCoeffVector* newcol=new SparseSDPCoeffVector;
      assert(newcol);
      opA_colrep[coli->first+old_dim]=newcol;
      for(SparseSDPCoeffVector::const_iterator rowi=(coli->second).begin();
	  rowi!=(coli->second).end();++rowi){
	if ((rowi->first<0)||(rowi->first>=b.dim())){
	  if (out){
	    (*out)<<"**** ERROR in MatrixSDPfunction::append_variables(...): row index "<<rowi->first<<" of column "<<coli->first<<" of opA_columns exceeds constraints range="<<b.dim()<<std::endl;
	  }
	  return 1;
	}
	if (rowi->second==0){
	  if (out){
	    (*out)<<"**** ERROR in MatrixSDPfunction::append_variables(...): coefficient matrix in row "<<rowi->first<<", column "<<coli->first<<" is a null pointer"<<std::endl;
	  }
	  return 1;
	}
	if (rowi->second->dim()!=append_Xdim(coli->first)){
	  if (out){
	    (*out)<<"**** ERROR in MatrixSDPfunction::append_variables(...): dimension of coefficient matrix in row "<<rowi->first<<", column "<<coli->first<<" and its matrix variable ("<<append_Xdim(coli->first)<<") do not match"<<std::endl;
	  }
	  return 1;
	}
	(*newcol)[rowi->first]=rowi->second;
	if (rowi->second->dense()) n_dense_coeffmat(coli->first)++;
	SparseSDPCoeffMatrix::iterator rrow=opA_rowrep.find(rowi->first);
	if (rrow==opA_rowrep.end()){
	  SparseSDPCoeffVector* row=new SparseSDPCoeffVector;
	  assert(row);
	  opA_rowrep[rowi->first]=row;
	  (*row)[coli->first+old_dim]=rowi->second;
	}
	else {
	  (*rrow->second)[coli->first+old_dim]=rowi->second;
	}
      } //for all elements of the column
    } //for all new columns
    
    Xdim.concat_below(append_Xdim);
    maxeigsolver.resize(Xdim.dim(),0);
    for(unsigned int i=old_dim;i<maxeigsolver.size();i++){
      maxeigsolver[i]=new MaxEigSolver;
      assert(maxeigsolver[i]);
      maxeigsolver[i]->set_out(out,print_level-1);
    }
    
    return 0;
  }
    
  
  
  /** append constraints to the problem. The #SparseSDPCoeffVector#
      of #append_opA[i]# belongs to matrix variable i. 
      The control over all coefficient matrices pointed to in 
      #append_opA# is passed over to this. This will delete them
      when delete constraints is called or on destruction.
      The other values are copied. */
  int MatrixSDPfunction::append_constraints(const std::map<Integer,SparseSDPCoeffVector>& opA_rows,
			 const Matrix& append_b)
  {
    if (append_b.dim()==0) return 0;
    last_bigmat_y.init(0,0,0.);
    Integer add_dim=append_b.dim();
    Integer old_dim=b.dim();
    b.concat_below(append_b);
    

    for(std::map<Integer,SparseSDPCoeffVector>::const_iterator rowi=opA_rows.begin();
	rowi!=opA_rows.end();++rowi){
      if ((rowi->first<0)||(rowi->first>=add_dim)){
	if (out){
	  (*out)<<"**** ERROR in MatrixSDPfunction::append_constraints(...): index "<<rowi->first<<" of a row of opA_rows exceeds number of added constraints="<<add_dim<<std::endl;
	}
	return 1;
      }
      SparseSDPCoeffVector* newrow=new SparseSDPCoeffVector;
      assert(newrow);
      opA_rowrep[rowi->first+old_dim]=newrow;
      for(SparseSDPCoeffVector::const_iterator coli=(rowi->second).begin();
	  coli!=(rowi->second).end();++coli){
	if ((coli->first<0)||(coli->first>=Xdim.dim())){
	  if (out){
	    (*out)<<"**** ERROR in MatrixSDPfunction::append_constraints(...): column index "<<coli->first<<" of row "<<rowi->first<<" of opA_rows exceeds variables range="<<Xdim.dim()<<std::endl;
	  }
	  return 1;
	}
	if (coli->second==0){
	  if (out){
	    (*out)<<"**** ERROR in MatrixSDPfunction::append_constraints(...): coefficient matrix in row "<<rowi->first<<", column "<<coli->first<<" is a null pointer"<<std::endl;
	  }
	  return 1;
	}
	if (coli->second->dim()!=Xdim(coli->first)){
	  if (out){
	    (*out)<<"**** ERROR in MatrixSDPfunction::append_constraints(...): dimension of coefficient matrix in row "<<rowi->first<<", column "<<coli->first<<"  and its matrix variable ("<<Xdim(coli->first)<<") do not match"<<std::endl;
	  }
	  return 1;
	}
	(*newrow)[coli->first]=coli->second;
	if (coli->second->dense()) n_dense_coeffmat(coli->first)++;
	SparseSDPCoeffMatrix::iterator ccol=opA_colrep.find(coli->first);
	if (ccol==opA_colrep.end()){
	  SparseSDPCoeffVector* col=new SparseSDPCoeffVector;
	  assert(col);
	  opA_colrep[coli->first]=col;
	  (*col)[old_dim+rowi->first]=coli->second;
	}
	else {
	  (*ccol->second)[old_dim+rowi->first]=coli->second;
	}
      } //for all elements of the column
    } //for all new columns

    return 0;
}  


  int MatrixSDPfunction::reassign_variables(const Indexmatrix& map_to_old,
					    Indexmatrix* deleted_Xdim,
					    SparseSDPCoeffVector* deleted_C,
					    SparseSDPCoeffMatrix* deleted_cols)
  {
    Indexmatrix ncopies(Xdim.dim(),1,Integer(0));
    for(Integer i=0;i<map_to_old.dim();i++){
      ncopies(map_to_old(i))++;
    }
   
    //delete entries with ncopies(i)==0 in C and opA_colrep
    {
      if (deleted_Xdim) deleted_Xdim->init(0,0,Integer(0));
      if (deleted_C) {
	if (deleted_C->size()>0){
	  if (out) {
	    (*out)<<"*** WARNING: MatrixSDPfunction::reassign_variables(): *deleted_C is not empty, clearing now"<<std::endl;
	  }
	  deleted_C->clear();
	}
      }
      if (deleted_cols) {
	if (deleted_cols->size()>0){
	  if (out) {
	    (*out)<<"*** WARNING: MatrixSDPfunction::reassign_variables(): *deleted_cols is not empty, clearing now"<<std::endl;
	  }
	  deleted_cols->clear();
	}
      }

      int cnt=0;     
      for(Integer i=0;i<ncopies.dim();i++){
	if (ncopies(i)>0) continue;
	if (deleted_Xdim){
	  deleted_Xdim->concat_below(Xdim(i));
	}
	SparseSDPCoeffVector::iterator ci=C.find(i);
	if (ci!=C.end()){
	  if (deleted_C){
	    (*deleted_C)[cnt]=ci->second;
	  }
	  else {
	    delete ci->second;
	  }
	  C.erase(ci);
	}
	SparseSDPCoeffMatrix::iterator coli=opA_colrep.find(i);
	if (coli!=opA_colrep.end()){
	  if (deleted_cols){
	    (*deleted_cols)[cnt]=coli->second;
	  }
	  else {
	    SparseSDPCoeffVector* colp=coli->second;
	    for(SparseSDPCoeffVector::iterator ci=colp->begin();ci!=colp->end();ci++){
	      delete ci->second;
	    }
	    delete colp;
	  }
	  opA_colrep.erase(coli);
	}
	cnt++;
      }
    }

    //reassign matrices
    last_bigmat_y.init(0,0,0.);
    n_dense_coeffmat=n_dense_coeffmat(map_to_old);
    Xdim=Xdim(map_to_old);
    maxeigsolver.resize(map_to_old.dim());

    //construct new cost vector and  column representation
    SparseSDPCoeffMatrix newcols;
    SparseSDPCoeffVector newcost;
    //these are used to rebuild opA_rowrep; discard old opA_rowrep
    for(SparseSDPCoeffMatrix::iterator rowi=opA_rowrep.begin();rowi!=opA_rowrep.end();rowi++){
      delete rowi->second;
    }
    opA_rowrep.clear();

    //build up new columns columnwise
    {for (Integer i=0;i<map_to_old.dim();i++){

      Integer ind=map_to_old(i);
      ncopies(ind)--;
      bool clone_it=(ncopies(ind)!=0);

      //cost vector
      SparseSDPCoeffVector::iterator Cit=C.find(ind);
      if (Cit!=C.end()){
        if (clone_it) {
	  newcost[i]=Cit->second->clone();
	}
	else {
          newcost[i]=Cit->second;
	  C.erase(Cit);
	}
      }

      //column representation             
      SparseSDPCoeffMatrix::iterator coli=opA_colrep.find(ind);
      if ((coli==opA_colrep.end())||(coli->second->size()==0)) continue;
      if (clone_it) {
        SparseSDPCoeffVector* newcol=new SparseSDPCoeffVector(*(coli->second));
        newcols[i]=newcol;
	for(SparseSDPCoeffVector::iterator nit=newcol->begin();nit!=newcol->end();nit++){
	  nit->second=nit->second->clone();
	}
      }
      else {
	newcols[i]=coli->second;
	opA_colrep.erase(coli);
      }
 
      //row representation
      SparseSDPCoeffVector* colp=newcols[i];
      for(SparseSDPCoeffVector::iterator ci=colp->begin();ci!=colp->end();ci++){
	SparseSDPCoeffMatrix::iterator rowi=opA_rowrep.find(ci->first);
	SparseSDPCoeffVector* rowp;
	if (rowi==opA_rowrep.end()){
	  rowp=new SparseSDPCoeffVector;
	  opA_rowrep[ci->first]=rowp;
	}
	else {
	  rowp=rowi->second;
	}
	(*rowp)[i]=ci->second;
      }
    }}

	    
    C.clear();
    opA_colrep.clear();

    //swap new maps to C and opA_colrep
    C.swap(newcost);
    opA_colrep.swap(newcols);
    
    return 0;
  }
    
  
  /** delete variables. If #map_to_old# is not null then
      the correspondonce of the new indices to the old indices of
      the remaining variables is described in #map_to_old# so 
      that now position i corresponds to the variable previously 
      indexed by #map_to_old(i)#. */        
  int MatrixSDPfunction::delete_variables(const Indexmatrix& delete_indices,
					  Indexmatrix* map_to_old,
					  Indexmatrix* deleted_Xdim,
					  SparseSDPCoeffVector* deleted_C,
					  SparseSDPCoeffMatrix* deleted_cols)
  {
    Indexmatrix mapvec(Range(0,Xdim.dim()-1));
    mapvec.delete_rows(delete_indices);
    if (map_to_old) {
      *map_to_old=mapvec;
    }
    return reassign_variables(mapvec,deleted_Xdim,deleted_C,deleted_cols);
  }
  
  /** delete constraints. The new assignment of the remaining 
      constraints is described by #map_to_old# so that aftwards
      position i is assigned the constraint previously indexed
      by #map_to_old(i)#. */        
  int MatrixSDPfunction::reassign_constraints(const Indexmatrix& map_to_old,
					      SparseSDPCoeffMatrix* deleted_rows,
					      Matrix* deleted_rhs)
  {
    Indexmatrix ncopies(b.dim(),1,Integer(0));
    {for(Integer i=0;i<map_to_old.dim();i++){
      ncopies(map_to_old(i))++;
    }}

    //remove rows with ncopies==0 in opA_rowrep
    {
      if (deleted_rhs) deleted_rhs->init(0,0,0.);
      if (deleted_rows) {
	if (deleted_rows->size()>0){
	  if (out) {
	    (*out)<<"*** WARNING: MatrixSDPfunction::reassign_constraints(): *deleted_rows is not empty, clearing now"<<std::endl;
	  }
	  deleted_rows->clear();
	}
      }
      Integer cnt=0;
      for(Integer i=0;i<b.dim();i++){
	if (ncopies(i)>0) continue;
	if (deleted_rhs){
	  deleted_rhs->concat_below(b(i));
	}
	SparseSDPCoeffMatrix::iterator rowi=opA_rowrep.find(i);
	if (rowi!=opA_rowrep.end()){
	  SparseSDPCoeffVector* rowp=rowi->second;
	  if (deleted_rows){
	    (*deleted_rows)[cnt]=rowp;
	  }
	  else {
	    for(SparseSDPCoeffVector::iterator ri=rowp->begin();ri!=rowp->end();ri++){
	      delete ri->second;
	    }
	    delete rowp;
	  }
	  opA_rowrep.erase(rowi);
	}
	cnt++;
      }
    }

    //construct new row representation
    SparseSDPCoeffMatrix newrows;
    //this is then used to rebuild opA_colrep; discard old opA_colrep
    for(SparseSDPCoeffMatrix::iterator coli=opA_colrep.begin();coli!=opA_colrep.end();coli++){
      delete coli->second;
    }
    opA_colrep.clear();
    n_dense_coeffmat.init(Xdim.dim(),1,Integer(0));


    //reassign the constraint data and evaluation point
    last_bigmat_y.init(0,0,0.);
    b=b(map_to_old);

    //build up new rows rowwise
    for (Integer i=0;i<map_to_old.dim();i++){

      Integer ind=map_to_old(i);
      ncopies(ind)--;
      bool clone_it=(ncopies(ind)!=0);

      //row representation             
      SparseSDPCoeffMatrix::iterator rowi=opA_rowrep.find(ind);
      if ((rowi==opA_rowrep.end())||(rowi->second->size()==0)) continue;
      if (clone_it) {
        SparseSDPCoeffVector* newrow=new SparseSDPCoeffVector(*(rowi->second));
        newrows[i]=newrow;
	for(SparseSDPCoeffVector::iterator nit=newrow->begin();nit!=newrow->end();nit++){
	  nit->second=nit->second->clone();
	}
      }
      else {
	newrows[i]=rowi->second;
	opA_rowrep.erase(rowi);
      }
 
      //column representation
      SparseSDPCoeffVector* rowp=newrows[i];
      for(SparseSDPCoeffVector::iterator ri=rowp->begin();ri!=rowp->end();ri++){
	n_dense_coeffmat(ri->first)+=Integer(ri->second->dense());
	SparseSDPCoeffMatrix::iterator coli=opA_colrep.find(ri->first);
	SparseSDPCoeffVector* colp;
	if (coli==opA_colrep.end()){
	  colp=new SparseSDPCoeffVector;
	  opA_colrep[ri->first]=colp;
	}
	else {
	  colp=coli->second;
	}
	(*colp)[i]=ri->second;
      }
    }

    opA_rowrep.clear();


    //swap new map to opA_rowrep
    opA_rowrep.swap(newrows);
    
    return 0;
  }
  
  /** delete constraints. If #map_to_old# is not null then
      the correspondonce of the new indices to the old indices of
      the remaining constraints is described in #map_to_old# so 
      that now position i corresponds to the constraint previously 
      indexed by #map_to_old(i)#. */        
  int MatrixSDPfunction::delete_constraints(const Indexmatrix& delete_indices,
					    Indexmatrix* map_to_old,
					    SparseSDPCoeffMatrix* deleted_rows,
					    Matrix* deleted_rhs
					    )
  {
    Indexmatrix mapvec(Range(0,b.dim()-1));
    mapvec.delete_rows(delete_indices);
    if (map_to_old) {
      *map_to_old=mapvec;
    }
    return reassign_constraints(mapvec,deleted_rows,deleted_rhs);
  }
  
  const Coeffmat* MatrixSDPfunction::get_coeffmat(Integer constr_nr,
						  Integer block_nr) const
  {
    const Coeffmat* cmp=0;
    if (constr_nr==-1){
      SparseSDPCoeffVector::const_iterator itcol=C.find(block_nr);
      if (itcol!=C.end())
	cmp=itcol->second;
      return cmp;
    }
    SparseSDPCoeffMatrix::const_iterator itrow=opA_rowrep.find(constr_nr);
    if (itrow!=opA_rowrep.end()){
      SparseSDPCoeffVector::const_iterator itcol=itrow->second->find(block_nr);
      if (itcol!=itrow->second->end()){
	cmp=itcol->second;
      }
    }
    return cmp;
  } 

  void  MatrixSDPfunction::set_out(std::ostream* o,int pril)
  {
    out=o;
    print_level=pril;
    for(unsigned int i=0;i<maxeigsolver.size();i++){
      maxeigsolver[i]->set_out(o,pril-1);
    }
  }

  std::ostream& MatrixSDPfunction::print_problem_data(std::ostream& o)
  {
    o<<"\nBEGIN_SDP_PROBLEM\n";
    if (constant_trace==SDPtrace_fixed){
      o<<"\nTRACE = "<<trace_value;
    }
    else {
      o<<"\nTRACE <= "<<trace_value;
    }
    o<<"\nBLOCKS\n";
    for(Integer i=0;i<Xdim.dim();i++){
      o<<" "<<Xdim(i);
    }
    o<<"\n";
    o<<"\nCOST_MATRICES\n";
    for(SparseSDPCoeffVector::iterator coli=C.begin();coli!=C.end();coli++){
      o<<"\n"<<coli->first<<"\n";
      coli->second->out(o);
    }
    o<<"\nCONSTRAINT_MATRICES\n"; 
    for(SparseSDPCoeffMatrix::iterator rowi=opA_rowrep.begin();
	rowi!=opA_rowrep.end();rowi++){
      SparseSDPCoeffVector& row=*rowi->second;
      for(SparseSDPCoeffVector::iterator coli=row.begin();coli!=row.end();coli++){
	o<<"\n"<<rowi->first<<" "<<coli->first<<"\n";
	coli->second->out(o);
      }
    }
    o<<"\nRIGHT_HAND_SIDE\n";
    o<<b;
    o<<"\nEND_SDP_PROBLEM"<<std::endl;
    return o;
  }

  std::istream& MatrixSDPfunction::read_problem_data(std::istream& in)
  {
    Indexmatrix ind;
    reassign_constraints(ind);
    reassign_variables(ind);
    char next_word[80];
    char next_char;
    if (! in.good()){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<" instream is not good";
      }
      return in;
    }
    in>>next_word;
    if (strcmp(next_word,"BEGIN_SDP_PROBLEM")!=0){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected BEGIN_SDP_PROBLEM but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    in>>next_word;
    if (strcmp(next_word,"TRACE")!=0){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected TRACE but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    in>>std::ws;
    in.get(next_char);
    if ((next_char!='=')&&(next_char!='<')){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected = or <= but got "<<next_char<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    SDPtrace sdptrace;
    if (next_char=='='){
      sdptrace=SDPtrace_fixed;
    }
    else {
      in.get(next_char);
      sdptrace=SDPtrace_bounded;
    }
    Real traceval;
    in>>traceval;
    if (traceval<0.){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"trace_value is negative : "<<traceval<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    set_trace(traceval,sdptrace);
    in>>next_word;
    if (strcmp(next_word,"BLOCKS")!=0){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected BLOCKS but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    do {
      in>>std::ws;
      next_char=char(in.peek());
      if ((next_char<'0')||(next_char>'9'))
	break;
      Integer next_dim;
      in>>next_dim;
      ind.concat_below(next_dim);
    } while (in.good());
    SparseSDPCoeffVector append_C;
    in>>next_word;
    if (strcmp(next_word,"COST_MATRICES")!=0){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected COST_MATRICES but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    do {
      in>>std::ws;
      next_char=char(in.peek());
      if ((next_char<'0')||(next_char>'9'))
	break;
      Integer next_i;
      in>>next_i;
      Coeffmat* cmp=coeffmat_read(in);
      if (cmp==0){
	if (out){
	  (*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
	  (*out)<<"coeffmat_read failed in reading the cost coefficient matrix for "<<next_i<<std::endl;
	}
	in.clear(in.rdstate()|std::ios::failbit);
	return in;
      }
      SparseSDPCoeffVector::iterator coli=append_C.find(next_i);
      if (coli==append_C.end()){
	append_C[next_i]=cmp;
      }
      else {
	if (out){
	  (*out)<<"*** WARNING in MatrixSDPfunction::read_problem_data(): ";
	  (*out)<<" the cost coefficient matrix for "<<next_i<<" is given more than once, using the most recent one!"<<std::endl;
	}
	delete coli->second;
	coli->second=cmp;
      }
    } while (in.good());
    std::map<CH_Matrix_Classes::Integer,SparseSDPCoeffVector> opA_rows;
    if (append_variables(ind,append_C,opA_rows)){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"append_variables failed"<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }      
    in>>next_word;
    if (strcmp(next_word,"CONSTRAINT_MATRICES")!=0){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected CONSTRAINT_MATRICES but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    // for the constraint matrices we use the routine
    // append constraints to avoid to work of building both
    // representations
    do {
      in>>std::ws;
      next_char=char(in.peek());
      if ((next_char<'0')||(next_char>'9'))
	break;
      Integer next_i;
      in>>next_i;
      Integer next_j;
      in>>next_j;
      if (next_j<0) {
	if (out){
	  (*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
	  (*out)<<"coeffmat_read failed in reading the column index to row "<<next_i<<": "<<next_j<<std::endl;
	}
	in.clear(in.rdstate()|std::ios::failbit);
	return in;
      }
      Coeffmat* cmp=coeffmat_read(in);
      if (cmp==0){
	if (out){
	  (*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
	  (*out)<<"coeffmat_read failed in reading the constraint coefficient matrix for ("<<next_i<<","<<next_j<<")"<<std::endl;
	}
	in.clear(in.rdstate()|std::ios::failbit);
	return in;
      }
      SparseSDPCoeffVector* cvp=&(opA_rows[next_i]);
      SparseSDPCoeffVector::iterator colj=cvp->find(next_j);
      if (colj==cvp->end()){
	(*cvp)[next_j]=cmp;
      }
      else {
	if (out){
	  (*out)<<"*** WARNING in MatrixSDPfunction::read_problem_data(): ";
	  (*out)<<" the constraint coefficient matrix for ("<<next_i<<","<<next_j<<") is given more than once, using the most recent one!"<<std::endl;
	}
	delete colj->second;
	colj->second=cmp;
      }
    } while (in.good());
    in>>next_word;
    if (strcmp(next_word,"RIGHT_HAND_SIDE")!=0){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected RIGHT_HAND_SIDE but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    Matrix rhsb;
    in>>rhsb;
    if (append_constraints(opA_rows,rhsb)){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"append_constraints failed"<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    in>>next_word;
    if (strcmp(next_word,"END_SDP_PROBLEM")!=0){
      if (out){
	(*out)<<"*** ERROR in MatrixSDPfunction::read_problem_data(): ";
        (*out)<<"expected END_SDP_PROBLEM but got "<<next_word<<std::endl;
      }
      in.clear(in.rdstate()|std::ios::failbit);
      return in;
    }
    return in;
  }

  int convert_row(const Indexmatrix& sdpdim,const Sparsemat& row,SparseSDPCoeffVector& opA_row)
  {
    assert(min(sdpdim)>0);
    Integer offset=0;
    for(Integer j=0;j<sdpdim.dim();j++){
      Integer n=sdpdim(j);
      Sparsemat Aij=row.cols(Range(offset,offset+n*n-1));
      offset+=n*n;
      if (Aij.nonzeros()==0) continue;
      Indexmatrix indi;
      Indexmatrix indj;
      Matrix val;
      Aij.get_edge_rep(indi,indj,val);
      for (Integer k=0;k<indi.dim();k++){
	Integer ii=indj(k)/n;
	Integer jj=indj(k)%n;
	indi(k)=ii;
	indj(k)=jj;
	val(k)=val(k);
      }
      Aij.init(n,n,indi.dim(),indi,indj,val);
      opA_row[j]=new CMsymsparse(Aij);
    }
    
    return 0;
  }

  int convert_to_rows(const Indexmatrix& sdpdim,const Sparsemat& A,std::map<Integer,SparseSDPCoeffVector>& opA_rows)
  {
    opA_rows.clear();
    for(Integer i=0;i<A.rowdim();i++){
      if (A.row_nonzeros(i)==0) continue;
      convert_row(sdpdim,A.row(i),opA_rows[i]);
    }
    
    return 0;
  }
  



  int convert_to_cols(const Indexmatrix& sdpdim,const Sparsemat& A,std::map<Integer,SparseSDPCoeffVector>& opA_cols)
  {
    opA_cols.clear();
    Integer offset=0;
    for(Integer j=0;j<sdpdim.dim();j++){
      Integer n=sdpdim(j);
      Sparsemat col=A.cols(Range(offset,offset+n*n-1));
      offset+=n*n;
      if (A.nonzeros()==0) continue;
      for(Integer i=0;i<A.rowdim();i++){
	Sparsemat Aij=col.row(i);
	if (Aij.nonzeros()==0) continue;
	Indexmatrix indi;
	Indexmatrix indj;
	Matrix val;
	Aij.get_edge_rep(indi,indj,val);
	for (Integer k=0;k<indi.dim();k++){
	  Integer ii=indj(k)/n;
	  Integer jj=indj(k)%n;
	  indi(k)=ii;
	  indj(k)=jj;
	  val(k)=val(k);
	}
	Aij.init(n,n,indi.dim(),indi,indj,val);
	(opA_cols[j])[i]=new CMsymsparse(Aij);
      }
    }
	  	    
    return 0;
  }



}

