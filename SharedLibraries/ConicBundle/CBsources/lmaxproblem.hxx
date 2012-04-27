/* ****************************************************************************

    Copyright (C) 2004-2011  Christoph Helmberg

    ConicBundle, Version 0.3.10
    File:  CBsources/lmaxproblem.hxx

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



#ifndef CONICBUNDLE_LMAXPROBLEM_HXX
#define CONICBUNDLE_LMAXPROBLEM_HXX

#include "qp_solver.hxx"
#include "qp_sdpblock.hxx"
#include "problem.hxx"
#include "BaseSDPOracle.hxx"

namespace ConicBundle {

  /** @brief allows to read the internal bundle data 
  */

  class SDPBundleInformation: public SDPBundleParameters
  {
  public:
    /** #n_new_subgradients# gives an upper bound on the number of new
        Ritz-vectors added in each bundle update
  
        //int n_new_subgradients;  //inherited
  
        #n_bundle_size# gives an upper bound on the number of 
	old Ritz-vectors that may be included in the new bundle
     
        //int n_bundle_size;       //inherited
       */

    /** minimum number of old Ritz-vectors that should be kept in the bundle
	no matter whether they fall below the aggregation tolerance or not. 
        
        //int n_keep;             //inherited
	*/
    
    /** maximum number of aggregate matrices allowed
        //int n_aggregates;       //inherited
	*/

    /** bundle eigenvectors (with respect to the solution matrix of the
        quadratci semidefinite subproblem) whose eigenvalues are below 
        this aggregation tolerance times the maximum eigenvalue may be
        aggregated into an aggregate matrix 
       
        //double aggregation_tolerance; //inherited
	*/

    /** the current implementation supplies two update rules, 
	the default rule (0) and one without aggregation (1). 
        When switching to (1), all linear aggregates are lost (so the entire
        model should be reinitialized ,otherwise the bundle method may report
        augmented modle failures), the parameters n_bundle_size,n_aggregates 
        are not used, aggregation_tolerance should be set very small and specifies
        the relative size of QSDP-eigenvalues to be considered as zero, 
	n_keep indicates the number of vectors that may and should be kept on top 
	of the active and the new vectors.
           
        //int update_rule;   //inherited

    */

    /** bundle vectors for the semidefinite part of the model */ 
    const CH_Matrix_Classes::Matrix* bundlevecs;
    
    /** topvecs and skippedvecs hold the vectors with the highest Ritz values
        as computed in the routine lmaxproblem::update_model_top_spectrum();
        after null steps, currently, topvecs are identical to bundlevecs,
        after descent steps they are identical to primalvecs */ 
    const CH_Matrix_Classes::Matrix* topvecs;
    
    /**  topvecs and skippedvecs hold the vectors with the highest Ritz values
        as computed in the routine lmaxproblem::update_model_top_spectrum(),
	skippedvecs are typically not included in the bundle */ 
    const CH_Matrix_Classes::Matrix* skippedvecs;
    
    /** Ritz value of topvecs und skippedvecs consecutively */ 
    const CH_Matrix_Classes::Matrix* Ritz_values;
    
    /** function values of the linear aggregate subgradients at zero */
    const CH_Matrix_Classes::Matrix* subgCvalues;    

    /** the linear aggregate subgradients forming the linear part of the model */
    const CH_Matrix_Classes::Matrix* aggrsubgrads;   

    /** the eigenvalues of the primal QP-solution X */
    const CH_Matrix_Classes::Matrix* primaleigs;   

    /** the eigenvectors of the primal QP-solution X */
    const CH_Matrix_Classes::Matrix* primalvecs;   

    /** the Ritz values of the eigenvectors of X with respect to the dual Z */
    const CH_Matrix_Classes::Matrix* primalZval;   

    /** for some update rules this gives the number of active eigenvalues */
    int activedim;   



    SDPBundleInformation()
    { n_new_subgradients=1;n_bundle_size=20;n_keep=5;n_aggregates=1;
      aggregation_tolerance=0.01;update_rule=0;
      bundlevecs=0;topvecs=0;skippedvecs=0;Ritz_values=0;subgCvalues=0;
      aggrsubgrads=0;primaleigs=0;primalZval=0;
      activedim=0;
    }
  };




class LmaxProblem: public ConvexProblem
{
private:

  //--- problem description
  BaseSDPOracle* oracle;
  
  CH_Matrix_Classes::Real lmax_multiplier;
  SDPtrace trace_stat; 
  
  CH_Matrix_Classes::Indexmatrix bounds_index;  //indices i with y(i) bounded (sorted increasingly)
  CH_Matrix_Classes::Matrix lby;                //lower bounds on y, same dim as y
  CH_Matrix_Classes::Matrix uby;                //upper bounds on y, same dim as y

  int problem_modified;   //is set to true if the problem description has 
                          //been changed

  //--- data of center point
  int center_available;   //is set to true if the following four items 
                          //   are available
  CH_Matrix_Classes::Matrix y;               //current center of stability
  CH_Matrix_Classes::Matrix subg;     //a minorant in the current model: subg_val+ip(subg,.-y)
  CH_Matrix_Classes::Real subg_val;        
  CH_Matrix_Classes::Real ub_fun_val;        //upper bound on the function value in y

  CH_Matrix_Classes::Real ublmax;            //value used as upper bound on maximum eigenvalue
  CH_Matrix_Classes::Matrix eigval;          //approximate maximum eigenvalue
  CH_Matrix_Classes::Matrix eigvec;          //and a corresponding Ritz vector 

  int use_startheur;  //if 1 a heuristic is applied in starting_point to improve
                      //the starting point if the eigenvalues are well separated

  //--- data of last function evaluation 
  int cand_available;
  CH_Matrix_Classes::Matrix cand_y;
  CH_Matrix_Classes::Real cand_ublmax;
  CH_Matrix_Classes::Matrix cand_eigs;
  CH_Matrix_Classes::Matrix cand_vecs;
  CH_Matrix_Classes::Matrix cand_subg;     
  CH_Matrix_Classes::Real cand_subg_val;
  CH_Matrix_Classes::Real cand_ub_fun_val;
  int ret_code; 

  //--- data of last model evaluation
  int model_available;
  int model_subg_available;
  CH_Matrix_Classes::Real model_subg_val;
  CH_Matrix_Classes::Real model_ub_val;
  CH_Matrix_Classes::Matrix model_subg;  //this subgradient is usually not needed, lazy evaluation 

  int model_changed;  //if true and model_ind>=0, 
                      //then lazy subg evaluation is impossible
  CH_Matrix_Classes::Real    model_eig;   //max eigenvalue over all modelvalues
  CH_Matrix_Classes::Matrix  model_vec;   //for lazy evaluation if a bundle vector is maximizer
  CH_Matrix_Classes::Integer model_ind;   //if an aggregate subgradient is the maximizer, this
                       //is its index, otherwise its value is -1
  int model_ret_code; 

  //--- data describing the model
  int update_rule;                        //allows to choose the update rule  
  CH_Matrix_Classes::Integer maxkeepvecs; //max number of vectors kept in bundle from prev iteration
  CH_Matrix_Classes::Integer minkeepvecs; //min number of vectors kept in bundle from prev iteration
  CH_Matrix_Classes::Integer maxaddvecs;  //max number of new vectors added to old bundle >=1!
  CH_Matrix_Classes::Integer maxaggrcols; //max number of aggregate vectors >= 1!
  CH_Matrix_Classes::Real aggregtol;      //aggregate primalvector if primalvalue smaller than
                       //
  CH_Matrix_Classes::Integer n_new_vecs;
  CH_Matrix_Classes::Matrix bundlevecs;
  CH_Matrix_Classes::Matrix primalvecs;  //primalvecs and pimaleigs are eigenvalue decomposition of
  CH_Matrix_Classes::Matrix primaleigs;  //the qp-solution for P*V*P' with P=bundle
  CH_Matrix_Classes::Matrix primalZval;
  CH_Matrix_Classes::Matrix primal_tapia;
  CH_Matrix_Classes::Matrix dual_tapia;
  CH_Matrix_Classes::Real tapia_factor;

  CH_Matrix_Classes::Matrix subgCvalues;    //=ip(C,Wj) for  Wj\succeq 0, trace(Wj)=1; 
  CH_Matrix_Classes::Matrix aggrsubgrads;   //=A(Wj) for same Wj's
  CH_Matrix_Classes::Matrix aggrcoeff;      //the qp-solution for the aggregates
  
  CH_Matrix_Classes::Integer rankdefcnt;  //counts number of rank deficiencies  


  //optional: store/update primal aggregate matrices corresp. to aggrsubgrads
  std::vector<SDPPrimal*> primal;    

  //--- data of last augmented model evaluation
  int aug_available;                     //!=0 only if aug_* matches primal* and aggr* information
  int aug_subg_in_model;                 //!=0 only if aug_subg and aug_linconst are in the model
  CH_Matrix_Classes::Integer aug_xdim;
  CH_Matrix_Classes::Matrix aug_subg;    //this is the aggregate subgradient
  CH_Matrix_Classes::Real aug_linconst;  //aug_linconst+<aug_subg,.> is a minorant of the objective
   
  //--- augmented modle solver
  QP_SDPBlock block;

  //--- temporary matrices
  mutable CH_Matrix_Classes::Matrix tmpmat;
  mutable CH_Matrix_Classes::Matrix tmpvec;
  mutable CH_Matrix_Classes::Symmatrix tmpsym;
  mutable CH_Matrix_Classes::Indexmatrix tmpind;
  
  //--- model updating schemes
  int update_model_default(bool descent_step);
  int update_model_no_aggregate(bool descent_step);
  int update_model_top_spectrum(bool descent_step);
  int update_model_top_spectrum2(bool descent_step);
  int update_model_top_spectrum3(bool descent_step);
  int update_model_top_spectrum4(bool descent_step);
  int update_model_top_spectrum5(bool descent_step);
  int update_model_top_spectrum6(bool descent_step);
  int update_model_top_spectrum7(bool descent_step);
  int update_model_top_spectrum8(bool descent_step);
  int update_model_new_spectrum(bool descent_step);
  int update_model_tapia_scale(bool descent_step);
  int update_model_tapia_scale1(bool descent_step);
  int update_model_tapia_scale2(bool descent_step);
  int update_model_tapia_scale3(bool descent_step);
  int update_model_tapia_scale4(bool descent_step);
  int update_model_tapia_scale5(bool descent_step);
  int update_model_tapia_scale6(bool descent_step);
  int update_model_tapia_scale7(bool descent_step);
  int update_model_tapia_scale8(bool descent_step);
  int update_model_tapia_scale9(bool descent_step);
  int update_model_tapia_scale10(bool descent_step);
  int update_model_tapia_scale11(bool descent_step);
  int update_model_tapia_scale12(bool descent_step);

  CH_Matrix_Classes::Matrix topvecs;
  CH_Matrix_Classes::Matrix skippedvecs;
  CH_Matrix_Classes::Matrix Ritz_values;
  CH_Matrix_Classes::Integer activedim;
  CH_Matrix_Classes::Integer keepsize;
  CH_Matrix_Classes::Real curve_delta;

  std::ostream *out;
  int print_level;
 
  //--- basic operations

  void init_subgrad(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& d);

  void update_subgrad(const CH_Matrix_Classes::Matrix& P,const CH_Matrix_Classes::Matrix& d,CH_Matrix_Classes::Integer aggr_index);

  void delete_subgrads(const CH_Matrix_Classes::Indexmatrix& delind);

public:    
  LmaxProblem(BaseSDPOracle* oracle_function);
    
  ~LmaxProblem();

  //--- routines  of BundleProblem
  
  int init_center(int& init, CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix& yp, CH_Matrix_Classes::Matrix& subg,
                 CH_Matrix_Classes::Indexmatrix& bounds_index, CH_Matrix_Classes::Matrix& lby, CH_Matrix_Classes::Matrix& uby);

  int eval_function(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real nullstep_bound,CH_Matrix_Classes::Real relprec);

  int get_function_sol(CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp);

  int do_step(void);

  int eval_model(const CH_Matrix_Classes::Matrix& y,CH_Matrix_Classes::Real eval_model_bound,CH_Matrix_Classes::Real relprec);

  int get_model_sol(CH_Matrix_Classes::Real& lb, CH_Matrix_Classes::Real& ub, CH_Matrix_Classes::Matrix* subgp);

  int get_augmodel_sol(CH_Matrix_Classes::Real& linconst,CH_Matrix_Classes::Matrix& subg);

  int update_model(bool descent_step);

  //--- routines of ConvexProblem
  
  int intersect_box(CH_Matrix_Classes::Indexmatrix& bound_index,CH_Matrix_Classes::Matrix& lb,CH_Matrix_Classes::Matrix& ub);
  
  CH_Matrix_Classes::Real lb_function(const CH_Matrix_Classes::Matrix& y);
  
  CH_Matrix_Classes::Real lb_model(const CH_Matrix_Classes::Matrix& y);
  
  int start_augmodel(QP_Block*& blockp);
  
  int get_row(CH_Matrix_Classes::Integer index_y,CH_Matrix_Classes::Matrix& row,CH_Matrix_Classes::Real& b,CH_Matrix_Classes::Integer startindex) const;
  
  int make_aug_linmodel(CH_Matrix_Classes::Real* aug_linconst,CH_Matrix_Classes::Matrix* aug_subg,
			bool* conemult_increased,
				CH_Matrix_Classes::Real* function_value);  
  
   
  int adjust_multiplier();
  
  int change_variables(ChangeVarInfo* cvp);
  
  int recompute_center();

  //--- other

  int get_approximate_primal(PrimalData& primal) const;

  int get_center_primal(PrimalData& primal) const;

  int get_candidate_primal(PrimalData& primal) const;

  int set_bundle_parameters(const BundleParameters& bp)
  {
    maxkeepvecs=CH_Matrix_Classes::max(1,bp.n_bundle_size);
    maxaddvecs=CH_Matrix_Classes::max(1,bp.n_new_subgradients);
    const SDPBundleParameters* p=dynamic_cast<const SDPBundleParameters*>(&bp);
    if (p==0) return 0;
    minkeepvecs=CH_Matrix_Classes::max(0,p->n_keep);
    maxaggrcols=CH_Matrix_Classes::max(1,p->n_aggregates);
    aggregtol=CH_Matrix_Classes::max(0.,p->aggregation_tolerance);
    update_rule=p->update_rule;
    const SDPBundleInformation* pi=dynamic_cast<const SDPBundleInformation*>(&bp);
    if (pi==0) return 0;
    if ((pi->bundlevecs!=0)&&(pi->bundlevecs->dim()>0)){
      bundlevecs=*(pi->bundlevecs);
      aug_subg_in_model=0;
    }
    return 0;
  }

  int get_bundle_parameters(BundleParameters& bp) const
  {
    bp.n_bundle_size=maxkeepvecs;
    bp.n_new_subgradients=maxaddvecs;
    SDPBundleParameters* p=dynamic_cast<SDPBundleParameters*>(&bp);
    if (p==0) return 0;
    p->n_keep=minkeepvecs;
    p->n_aggregates=maxaggrcols;
    p->aggregation_tolerance=aggregtol;
    p->update_rule=update_rule;
    SDPBundleInformation* pi=dynamic_cast<SDPBundleInformation*>(&bp);
    if (pi==0) return 0;
    pi->bundlevecs=&bundlevecs;
    pi->topvecs=&topvecs;
    pi->skippedvecs=&skippedvecs;
    pi->Ritz_values=&Ritz_values;
    pi->subgCvalues=&subgCvalues;
    pi->aggrsubgrads=&aggrsubgrads;
    pi->primaleigs=&primaleigs;
    pi->primalvecs=&primalvecs;
    pi->primalZval=&primalZval;
    pi->activedim=activedim;
    return 0;
  }

  int get_bundle_values(BundleParameters& bp) const
  {
    bp.n_bundle_size=bundlevecs.coldim();
    bp.n_new_subgradients=n_new_vecs;
    SDPBundleParameters* p=dynamic_cast<SDPBundleParameters*>(&bp);
    if (p==0) return 0;
    p->n_keep=CH_Matrix_Classes::min(minkeepvecs,bundlevecs.coldim());
    p->n_aggregates=aggrcoeff.dim();
    p->aggregation_tolerance=aggregtol;
    p->update_rule=update_rule;
    SDPBundleInformation* pi=dynamic_cast<SDPBundleInformation*>(&bp);
    if (pi==0) return 0;
    pi->bundlevecs=&bundlevecs;
    pi->topvecs=&topvecs;
    pi->skippedvecs=&skippedvecs;
    pi->Ritz_values=&Ritz_values;
    pi->subgCvalues=&subgCvalues;
    pi->aggrsubgrads=&aggrsubgrads;
    pi->primaleigs=&primaleigs;
    pi->primalvecs=&primalvecs;
    pi->primalZval=&primalZval;
    pi->activedim=activedim;
    return 0;
  }

  void clear_model();  

  int clear_aggregates();

  int get_ret_code() const {return ret_code;}

  //--- other

  void get_y(CH_Matrix_Classes::Matrix& iny) const {iny=y;}
  const CH_Matrix_Classes::Matrix& get_y(void) const {return y;}

  void set_out(std::ostream* o=0,int pril=1)
  { 
    out=o;print_level=pril;     
    if (solverp) solverp->set_out(out,pril-1);
    block.set_out(out,pril-1);
  }

};

}

#endif

