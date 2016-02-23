/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef FILTERIPMALGORITHM_H
#define FILTERIPMALGORITHM_H

#include "Solver.h"
#include "OoqpVectorHandle.h"

class Data;
class Variables;
class ProblemFormulation;
class NlpInfo;

class RegularizationAlg;

class Filter2D;
class FilterIPMOption;
class PDRegularization;

/** 
 * Derived class of Solver implementing filter line search algorithm.
 * @ingroup NlpSolvers
 */
class FilterIPMSolver : public Solver
{
  protected:

	double *dualRegQuantities;
	double ALObj, ALPCon;
	
  	double xWx_0, xWx_done;
	double thd_0, thd_done;
	
  	int printlevel;
	
	double rho; // penalty parameter
	double mu, oldmu;
	double scal_commerr, scal_dualerr;	
	double alp_pri, alp_boundDual, alp_pri_max, alp_boundDual_max;
	double curr_ConNorm, curr_BarrObj, curr_BarrObjGradTimesD, curr_ConTimesD;
	double test_ConNorm, test_BarrObj;

	char* strPrt;
	char* strPrt_du;
	char strbuffer[8];
	char strbuffer2[8];
	char dashes[6];

	bool FirstTimeUpdateMu;
	
  	// storage for variables and residuals 
    Variables *mainstep, *SOCstep, *t_step, *n_step;	
    Variables *trial_iter, *KKT_sol;
	
    Residuals *trialIter_rhs_SOC;
	Residuals *KKT_Resid;	

    double tau_j,oldtau_j;
	double deltaRegular_pr,deltaRegular_du;	

	int NumBackCheck;

	bool DoFeasResto;
    bool SolfromSOC;
	bool IfTinyStep;
	bool DualStepIsTiny;
	bool LastIterRejectByFilter;

  	int LastIterAcceptByCondition;

	int cumulative_small_step;
	int cumulative_reject_filter;

  	int StepAcceptDueTo_SWC_AC;
  	int StepAcceptDueTo_SRC_obj;
  	int StepAcceptDueTo_SRC_con;
	int StepAcceptDueTo_SRC;

	char info_xWx_type;
	char info_PriAlpha_char;
	char info_PriAlpha_char_1st;
	
	Filter2D *Filter;
	FilterIPMOption *FilterIPMOpt;  
	PDRegularization *PD_Reg;
    ProblemFormulation * factory;  


	double AugLagVaryPart;
	
  public:

    FilterIPMSolver(ProblemFormulation * of, Data * prob );
    ~FilterIPMSolver();

	virtual int solve( Data *prob, Variables *iterate, Residuals * resid );

protected: 
  // modify var bounds
  virtual void ChangeBounds(Data *prob_in,Residuals *resid_in);

  // initialize var
  virtual void defaultStart( ProblemFormulation * formulation, Variables * iterate, Data * prob,
			 				Residuals * resid, Variables * rhs  );
  // initialize var  
  virtual void _initialVar(Variables * iterate_in, Data * prob_in, Residuals * resid_in);

  // get  optimal error scaling factor
  virtual void EvalErrScaling(Variables *iterate_in);

  // initialize filter and some parameters
  virtual void FilterInitializeAndPara(Data *prob_in, Variables *iterate_in, Residuals * resid_in);

  // update mu, reset filter if necessary
  virtual int updateBarrierParameter(Data *prob_in, Variables *iterate_in, Residuals * resid_in);

  // add damping term to where var/con only has single bound
  virtual void addDampingTermToKKT(Residuals *resid_in);

  // compute min step size
  virtual double getAlphaMin(const double priAlaph_in);

  // check if step is too tiny
  virtual bool isTinyStep(Variables *vars_in, Residuals * resid_in, Variables *step_in);

  // check if step is satisfied
  virtual bool testStep( const double trialAlpha, const double trial_ConsNorm_in, const double trial_Obj_in,
				const double curr_ConsNorm_in, const double curr_Obj_in, const double curr_BarrGradTimesD_in, const double curr_ConsTimesD_in, 
				FilterIPMOption *FilterIPMOpt);

  // test switching condition
  virtual bool testSwitchingCondition( const double trialAlpha,
				const double curr_filter_ConsNorm_in, const double augmented_BarrGradTimesD, FilterIPMOption *FilterIPMOpt);

  // test Armijo condition
  virtual bool testArmijo( const double trialAlpha, const double trial_Obj_in, const double curr_Obj_in,
				const double augmented_BarrGradTimesD, FilterIPMOption *FilterIPMOpt);

  // test sufficient reduction condition
  virtual bool testSufficientReductionCondition( const double trial_ConsNorm_in, const double trial_Obj_in,
				const double curr_ConsNorm_in, const double curr_Obj_in, FilterIPMOption *FilterIPMOpt);

  // accept trial step  
  virtual void accept_iterate(Data *data_in, Variables * vars_in, Residuals *resid_in, Variables * step_in);

  // compute xWx
  virtual double get_xWx( Data *prob_in, Residuals * resid_in, Variables *step_in);

  // fact/solve linear syste to get search direction
  virtual void compute_step_WithRegularization(Data * data_in, Variables * vars_in, 
			  Residuals * resids_in, Variables *step_in, bool inSOC = false);

  // evaluate barrobj
  virtual double evalBarrObj(Data *prob_in, Variables *iterate_in, Residuals * resid_in, const int isTrialIter=0);

  // evaluate constraint one norm
  virtual double evalConOneNorm(Data *prob_in, Variables *iterate_in, Residuals * resid_in, const int isTrialIter=0);
  
  // compute second order corrector
  virtual double computeSOCIter(Data *prob_in, Variables *iterate_in, Residuals * resid_in, 
				  const double AlphaStep);

  // check if trial step is accpted by filter and other conditions
  virtual bool ifStepAccepted(Data *prob_in, Variables *vars_in, Residuals * resid_in, 
				  Variables *step_in, Variables *trial_iter_in, const double trialAlpha, const bool inSOC = false);   

  // do filter line search to find step size
  virtual void line_search_Filter(Data * prob_in, Variables *vars_in, Residuals * resid_in, Variables *step_in);

  // do backtracking line search
  virtual int BackTrack( Data *prob_in, Variables *vars_in, Residuals * resid_in, Variables *step_in);

  // update filter if necessary
  virtual void UpdateFilter();

  // accept tiny step
  virtual void use_TinyStep(Variables *step_in);

  virtual void defaultMonitor( Data * data, Variables * vars, Residuals * resids,
							 double alpha, double sigma, int i, double mu, 
							 int status_code, int level ) ;	

  // print out solver stats
  virtual int defaultStatus(Data *  data_in, Variables * /* vars */, Residuals * resids_in,
			 int iterate, double mu, int /* level */);


  // compute quantities for Dual regularization
  virtual double* computeQuantitiesForDualReg(	Data *prob_in, Variables *iter_in, Residuals * resid_in, 
													Variables *step_in, PDRegularization* RegInfo);

};


#endif
