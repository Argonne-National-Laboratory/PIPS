/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include <cmath>
#include <cstdio>
#include <cassert>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>

#include "FilterIPMSolver.h"
#include "FilterIPMOption.h"
#include "Filter2D.h"
#include "PDRegularization.h"

#include "Variables.h"
#include "Residuals.h"
#include "LinearSystem.h"
#include "Status.h"
#include "Data.h"
#include "ProblemFormulation.h"

#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"
#include "NlpGenLinsys.h"

#include "../../global_var.h"
#ifdef TIMING
#include <mpi.h>
#include "../PIPS-NLP/Core/Utilities/PerfMetrics.h"
#endif

extern int gOoqpPrintLevel;
extern int gdWd_test;
extern int gUseFilter;
extern int gDoTinyStepTest;
extern int gFilterResetStep;
extern int gLineSearchMatStep;
extern int gdWd_test_soc;

extern double gkappa_tWt;
extern int gCheckSmallConstVio;
extern int gDoSOC;
extern int gkappaWithMu;

extern int gDoIR_Full;

extern int gUseDualRegAlg;


#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif


/** 
 * Derived class of Solver implementing filter line search algorithm.
 * @ingroup NlpSolvers
 */
FilterIPMSolver::FilterIPMSolver( ProblemFormulation * of, Data * prob )
	: 	KKT_Resid(NULL),
		KKT_sol(NULL),
		trialIter_rhs_SOC(NULL),
		SOCstep(NULL),
		dualRegQuantities(NULL)
	
{
  factory              = of;
  mainstep             = factory->makeVariables( prob );
  trial_iter           = factory->makeVariables( prob );

  if(gDoIR_Full==1){
    KKT_sol			 = factory->makeVariables( prob );
    KKT_Resid		 = factory->makeResiduals( prob );
  }
  
  FilterIPMOpt 	= new FilterIPMOption();
  Filter 		= new Filter2D();
  PD_Reg 		= new PDRegularization(FilterIPMOpt);

  maxit 	 = FilterIPMOpt->maxit;
  printlevel = FilterIPMOpt->printlevel;

  // allocate space to track mu and residual norms
  mu_history    = new double[maxit];
  rnorm_history = new double[maxit];

  // Use the defaultStatus method
  status   = 0; 

  FirstTimeUpdateMu 	 = true;
  DoFeasResto 			 = false;
  SolfromSOC 			 = false;
  DualStepIsTiny 		 = false;
  LastIterRejectByFilter = false;

  strcpy(dashes, "   - ");  
  info_PriAlpha_char 	 = 'n';
  info_PriAlpha_char_1st = '-';

  // for inertia-free test
  xWx_0  = 0.0; xWx_done=0.0;
  thd_0  = 0.0; thd_done=0.0;
  info_xWx_type = 'd';
  if(gdWd_test == 2 ){ 
	n_step = factory->makeVariables( prob );  	
    t_step = factory->makeVariables( prob );
    info_xWx_type = 't';	
  }

  NumBackCheck = 0;
  StepAcceptDueTo_SWC_AC 	= 0;
  StepAcceptDueTo_SRC_obj 	= 0;
  StepAcceptDueTo_SRC_con 	= 0;
  StepAcceptDueTo_SRC       = 0;
  
  LastIterAcceptByCondition = 0;

  cumulative_small_step 	= 0;
  cumulative_reject_filter	= 0;

  AugLagVaryPart = 0;
}

FilterIPMSolver::~FilterIPMSolver()
{
  if(trialIter_rhs_SOC) delete trialIter_rhs_SOC;
  if(SOCstep) delete SOCstep;

  if(gDoIR_Full==1){
    delete KKT_sol;
    delete KKT_Resid;
  }
  delete trial_iter;
  delete mainstep;

  delete PD_Reg;
  delete Filter;
  delete FilterIPMOpt; 

  delete [] mu_history;
  delete [] rnorm_history;
  
  if(sys) delete sys;

  if(dualRegQuantities) delete [] dualRegQuantities;
}


void
FilterIPMSolver::ChangeBounds(Data *prob_in, Residuals *resid_in)
{
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);  
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  
  double move_bound_scalar = FilterIPMOpt->move_bound_scalar;

  prob->moveBounds(resid->priWrk,resid->priWrk_S,move_bound_scalar);
}

void 
FilterIPMSolver::defaultStart( ProblemFormulation * formulation,
			   Variables * iterate_in, Data * prob_in,
			   Residuals * resid_in, Variables * rhs  )
{   
  _initialVar(iterate_in,prob_in,resid_in);
}


/*
	 Starting point is obtained by:

	only  lower bound:
	x0 <--	max{ x0, lb + k1*max{1,|lb|}

	only upper bound:
	x0 <--	min{ x0, ub - k1*max{1,|ub|}

	with both bounds:
	  pL(i) = min(k_1*max(1,abs(lb(i))),k_2*(ub(i)-lb(i)));
	  pU(i) = min(k_1*max(1,abs(ub(i))),k_2*(ub(i)-lb(i)));

	x is projected into the interval [lb+pL,ub-pU]

	The dual variables (y, z) are initialized to 0

	The bound multipliers (lambda pi gamma phi) are initialized to one
*/ 
void 
FilterIPMSolver::_initialVar(Variables * iterate_in, Data * prob_in, Residuals * resid_in)
{

  NlpGenVars *vars = (NlpGenVars *) iterate_in;
  NlpGenData *prob = dynamic_cast<NlpGenData*>(prob_in);
  NlpGenResiduals *resid = (NlpGenResiduals *) resid_in;  
  
  double k_1 = FilterIPMOpt->k_1;
  double k_2 = FilterIPMOpt->k_2;

  // initialize x and  push x from bounds 
  MESSAGE("_initalVar - before getting init x");
  prob->getInitX(vars->x);
  MESSAGE("_initalVar - after getting init x");
  vars->push_variables(vars->x, vars->v, vars->w, prob->blx, prob->bux,resid->priWrk,k_1,k_2,1);
  IF_VERBOSE_DO( vars->x->print(); );
  // initialize slack s, and t=s-lb, u=ub-s  note that Ci(x)-s =0, put primal var from bounds
  prob->evalConstraintBody(vars);
  prob->getInEqCons(*vars->s);
  vars->push_variables(vars->s, vars->t, vars->u, prob->bl, prob->bu,resid->priWrk_S,k_1,k_2,0);

  // initialize bound multipliers lambda pi gamma phi  
  vars->interiorBoundSlackDual(FilterIPMOpt->dual_var_init);

  // initialize cons dual var: y, z
  vars->interiorPointDualY(0.0);
  vars->interiorPointDualZ(0.0);
}


void
FilterIPMSolver::EvalErrScaling(Variables *iterate_in)
{
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (iterate_in);
  
  vars->getErrScaling(FilterIPMOpt->s_max,scal_commerr,scal_dualerr);
}



void
FilterIPMSolver::FilterInitializeAndPara(Data *prob_in, Variables *iterate_in, Residuals * resid_in)
{	
  double initConsOneNorm = evalConOneNorm(prob_in,iterate_in,resid_in);
  
  FilterIPMOpt->ConNorm_min = 1e-4*((1>initConsOneNorm)?1:initConsOneNorm);	
  FilterIPMOpt->ConNorm_max = 1e4 *((1>initConsOneNorm)?1:initConsOneNorm) ;	   
	
  if(gUseFilter==1)
  	Filter->Initialize(FilterIPMOpt);
}


int
FilterIPMSolver::updateBarrierParameter(Data *prob_in, Variables *iterate_in, Residuals * resid_in)
{
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (iterate_in);

  int barr_changed = 0; 
  double tempDouble =0.0;
  bool done = false, FastMuDecrease = true;

  double BarrProbRes = 0.0;

  double BarrKKTRes_Pri  = resid->priErr();
  double BarrKKTRes_Dual = resid->dualErr()/scal_dualerr;
  double BarrKKTRes_Comm = resid->comp_Err()/scal_commerr;
  
  BarrProbRes = MAX(BarrKKTRes_Pri,BarrKKTRes_Dual);
  BarrProbRes = MAX(BarrProbRes,BarrKKTRes_Comm);
  
  while(BarrProbRes <= FilterIPMOpt->k_ep*mu && !done ){
	oldmu = mu;
		
	mu = (FilterIPMOpt->k_mu)*oldmu;
	tempDouble = pow(oldmu,FilterIPMOpt->Q_mu);
	mu = (mu<tempDouble)?mu:tempDouble;
		
	mu = (mu>(FilterIPMOpt->opt_tol/11.))?mu:(FilterIPMOpt->opt_tol/11.);
	tau_j = ((FilterIPMOpt->tau_Min)>(1-mu))?(FilterIPMOpt->tau_Min):(1-mu);
	barr_changed++;

    BarrKKTRes_Comm = resid->getKKTError_Comp(prob, vars,mu,PIPS_NORM_INF) /scal_commerr;
	BarrProbRes = MAX(BarrKKTRes_Pri,BarrKKTRes_Dual);
	BarrProbRes = MAX(BarrProbRes,BarrKKTRes_Comm);

    if(FirstTimeUpdateMu || FastMuDecrease){
	  done = false;
    }else{
	  done = true;
	}
  }

  FirstTimeUpdateMu = false;

  prob->currMu = mu;

  // if mu has beend changed, update BarrObjValue with damping rate
  if(barr_changed>0)
	prob->BarrObj = prob->BarrObjValue(vars, prob->PriObj, FilterIPMOpt->k_d);

  // if mu has beend changed, reset filter
  if(barr_changed > 0 && gUseFilter==1){
  	Filter->Initialize(FilterIPMOpt);
	// FIXME_NY: Following two lines are new! 2015/09/23
    cumulative_reject_filter = 0;
	LastIterRejectByFilter=false;	
  }
  return barr_changed;
}



void
FilterIPMSolver::addDampingTermToKKT(Residuals *resid_in)
{
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);

  double dampingFact = FilterIPMOpt->k_d*mu;
  resid->addDampingTermToOneSidePart(dampingFact);
}


double
FilterIPMSolver::getAlphaMin(const double priAlaph_in)
{

  double alp_pri_min = FilterIPMOpt->r_Q;	
  double dtmp;

  double augmented_BarrGradTimesD;
//  augmented_BarrGradTimesD = curr_BarrObjGradTimesD + (1.0-priAlaph_in)*curr_ConTimesD;
  augmented_BarrGradTimesD = curr_BarrObjGradTimesD + priAlaph_in*AugLagVaryPart;
									
  if(augmented_BarrGradTimesD<0){
	dtmp =  -FilterIPMOpt->r_phi*curr_ConNorm/augmented_BarrGradTimesD;
	alp_pri_min = MIN( alp_pri_min, dtmp);
	  
	if(curr_ConNorm<=FilterIPMOpt->ConNorm_min){
	  dtmp = FilterIPMOpt->delta*pow(curr_ConNorm,FilterIPMOpt->powScalar_ConNorm)/(-augmented_BarrGradTimesD);
	  alp_pri_min = MIN( alp_pri_min, dtmp); 
	}
  }
  alp_pri_min = FilterIPMOpt->r_alpha*alp_pri_min;

  return alp_pri_min;
}


bool 
FilterIPMSolver::isTinyStep(Variables *vars_in, Residuals * resid_in, Variables *step_in)
{
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (vars_in);
  NlpGenVars *steps = dynamic_cast<NlpGenVars *> (step_in);   
  
  bool result = resid->findSmallStep(vars,steps,FilterIPMOpt->tol_mach);
  
  return result;
}


bool
FilterIPMSolver::testStep( const double trialAlpha, const double trial_ConsNorm_in, const double trial_Obj_in,
				const double curr_ConsNorm_in, const double curr_Obj_in, const double curr_BarrGradTimesD_in, const double curr_ConsTimesD_in,
				FilterIPMOption *FilterIPMOpt)
{
  bool acceptStep=false;

//  double augmented_BarrGradTimesD = curr_BarrGradTimesD_in - (1.0-trialAlpha)*curr_ConsTimesD_in;
  double augmented_BarrGradTimesD = curr_BarrGradTimesD_in + trialAlpha*AugLagVaryPart;

  // test switching condition
  if( testSwitchingCondition( trialAlpha, curr_ConsNorm_in, augmented_BarrGradTimesD, FilterIPMOpt) )
  {  
    // test Armijo condition
	acceptStep = testArmijo( trialAlpha, trial_Obj_in, curr_Obj_in, augmented_BarrGradTimesD, FilterIPMOpt);
  }
  else{
	// test sufficient reduction condition  	
	acceptStep = testSufficientReductionCondition(trial_ConsNorm_in, trial_Obj_in, 
									curr_ConsNorm_in, curr_Obj_in, FilterIPMOpt);
  }
  
  return acceptStep;
}


bool
FilterIPMSolver::testSwitchingCondition( const double trialAlpha,
				const double curr_ConsNorm_in, const double augmented_BarrGradTimesD, FilterIPMOption *FilterIPMOpt)
{
  bool SWChold = false;
  
  if(gCheckSmallConstVio == 1){
	if( curr_ConsNorm_in > FilterIPMOpt->ConNorm_min){
	  return SWChold;
    }
  }

  double valTemp1 = pow(-augmented_BarrGradTimesD,FilterIPMOpt->powScalar_Obj); 
  double valTemp2 = pow(curr_ConsNorm_in,FilterIPMOpt->powScalar_ConNorm);

  if( augmented_BarrGradTimesD < 0 && trialAlpha > FilterIPMOpt->delta*valTemp2/valTemp1)
  {
	SWChold = true;
  }
  return SWChold;
}


bool
FilterIPMSolver::testArmijo( const double trialAlpha,const double trial_Obj_in,
				const double curr_Obj_in, const double augmented_BarrGradTimesD, FilterIPMOption *FilterIPMOpt)
{
  bool AChold = false;

  double testLHS = trial_Obj_in - curr_Obj_in-10*(FilterIPMOpt->tol_mach)*fabs(curr_Obj_in);
  double testRHS = FilterIPMOpt->y_phi*trialAlpha*augmented_BarrGradTimesD;

  if( testLHS <= testRHS)
  {
    AChold = true;
	LastIterAcceptByCondition=1;
  }
  return AChold;
}


bool
FilterIPMSolver::testSufficientReductionCondition( const double trial_ConsNorm_in, const double trial_Obj_in,
				const double curr_ConsNorm_in, const double curr_Obj_in,FilterIPMOption *FilterIPMOpt)
{
  bool SRChold = false;
  bool SRC_obj = 	trial_Obj_in-10 * (FilterIPMOpt->tol_mach)*fabs(curr_Obj_in)
  				<= 	(curr_Obj_in-FilterIPMOpt->r_phi*curr_ConsNorm_in);
  bool SRC_con = 	trial_ConsNorm_in 
  				<=	(1-FilterIPMOpt->r_Q)*curr_ConsNorm_in;

  if( SRC_obj && SRC_con)
  {
   	SRChold = true;
	LastIterAcceptByCondition=4; 
  }else if( SRC_obj ){
	SRChold = true;
	LastIterAcceptByCondition=2;
  }else if( SRC_con ){
	SRChold = true;
	LastIterAcceptByCondition=3;
  }
  return SRChold;
}


void 
FilterIPMSolver::accept_iterate(Data *data_in, Variables * vars_in, Residuals *resid_in, Variables * step_in)
{
  NlpGenVars*  vars  = (NlpGenVars *) vars_in;
  NlpGenData * data  = (NlpGenData *) data_in;
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenVars*  steps  = (NlpGenVars *) step_in;

  // allow max step size for barrier multipliers
  alp_boundDual = alp_boundDual_max;
  
  // accept iterate
  vars->takeStep(step_in, alp_pri,alp_pri,alp_boundDual);
	  
  // update slack dual
  vars->updateSlackAndDual(resid->priWrk,resid->priWrk_S,FilterIPMOpt->k_sigma,mu); 

}


double  
FilterIPMSolver::get_xWx( Data *prob_in, Residuals * resid_in, Variables *step_in)
{
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenVars *steps = dynamic_cast<NlpGenVars *> (step_in);
  NlpGenLinsys *sysNLP = dynamic_cast<NlpGenLinsys *>(sys);
  
  double result;

  result = sysNLP->eval_xWx(prob,resid,steps); 

  return result;
}


void
FilterIPMSolver::compute_step_WithRegularization(Data *prob_in, Variables *iterate_in, 
					Residuals * resid_in, Variables *step_in, bool inSOC)
{
  NlpGenData *prob 		= (NlpGenData *) prob_in;
  NlpGenVars *vars 		= (NlpGenVars *) iterate_in;
  NlpGenVars *steps 	= (NlpGenVars *) step_in;
  NlpGenLinsys *sysNLP 	= dynamic_cast<NlpGenLinsys *>(sys);
	
#ifdef TIMING
	double stime0;
#endif
    
  if(inSOC==false){
#ifdef TIMING
	stime0=MPI_Wtime();
#endif
	// form rhs of linear system:
	resid_in->set_r3_xz_alpha( iterate_in, -mu );   

	// DoEvalReg =	1 -> when factorizing the matrix, add regularizations to correct inertia and singularity
	//			  	2 -> when factorizing the matrix, add regularizations to correct singularity only
	//			  	0 -> when factorizing the matrix, force to use primal regularizaion. called iff xWx tests fail
	if(gdWd_test!=0) 
	  PD_Reg->DoEvalReg=2;
	else 
	  PD_Reg->DoEvalReg=1;
#ifdef TIMING
  gprof.t_set_r3_xz_alpha+=MPI_Wtime()-stime0;
	stime0=MPI_Wtime();
#endif
	sysNLP->factor(prob, iterate_in,PD_Reg);

    double xWx, thd;
	double kappa_test = gkappa_tWt;
	if(gkappaWithMu) kappa_test *= mu;
	
#ifdef TIMING
  gprof.t_factor+=MPI_Wtime()-stime0;
	stime0=MPI_Wtime();
#endif
    if(gdWd_test<=1){
	  // solve d step, compute dWd and thd (this is the first fact/solve)
	  sysNLP->solve_IterRefine(prob, iterate_in, resid_in, step_in, KKT_Resid, KKT_sol);
      step_in->negate();
	  xWx_0 = get_xWx(prob_in, resid_in,step_in);
	  thd_0 = kappa_test*steps->computeXSDD(steps);  
    }
	else if(gdWd_test>=2){
	  // solve n&t steps, do tWt test
	  sysNLP->solve_NTsteps(prob, iterate_in, resid_in, n_step, t_step, step_in);
      step_in->negate();
	  xWx_0 = get_xWx(prob_in, resid_in,t_step);
	  thd_0 = kappa_test*steps->computeXSDD(t_step);	    
    }
#ifdef TIMING
				gprof.t_computeXSDD1+=MPI_Wtime()-stime0;
				stime0=MPI_Wtime();
#endif
	if(gUseDualRegAlg>0) computeQuantitiesForDualReg(prob_in,iterate_in,resid_in,step_in, PD_Reg);
#ifdef TIMING
	gprof.t_computeQuantitiesForDualReg+=MPI_Wtime()-stime0;
	stime0=MPI_Wtime();
#endif

	// do inertia-free test, if fails, add regularization
	xWx   = xWx_0;
	thd   = thd_0;
	if(gdWd_test!=0){
	  while( xWx<thd ){
	    PD_Reg->DoEvalReg=0;
		sysNLP->factor(prob, iterate_in,PD_Reg);

		if(gdWd_test==1){
		  // get d, do dWd test
		  sysNLP->solve_IterRefine(prob, iterate_in, resid_in, step_in, KKT_Resid, KKT_sol);
		  step_in->negate();
		  xWx = get_xWx(prob_in, resid_in,step_in);
		  thd = kappa_test*steps->computeXSDD(steps);
		}
		else if(gdWd_test==2){
		  // get t, do tWt test
		  sysNLP->solve_NTsteps(prob, iterate_in, resid_in, n_step, t_step, step_in);
		  step_in->negate();
		  xWx = get_xWx(prob_in, resid_in,t_step);
		  thd = kappa_test*steps->computeXSDD(t_step);
		}

		if( xWx >= thd ){
		  break;	
		}
	  }
	}
	xWx_done = xWx;
	thd_done = thd;
  }
  else{  
    // here we have already got correct reg values
    // rhs for second order correction is saved in resid_in, which is computed in computeSOCIter
	// no SOC for n_step and t_step, use full step
	
	sysNLP->solve_IterRefine(prob, iterate_in, resid_in, step_in, KKT_Resid, KKT_sol);	
	step_in->negate();

	if(gdWd_test && gdWd_test_soc){
	  double kappa_test = gkappa_tWt;
	  if(gkappaWithMu) kappa_test *= mu;
      xWx_done 	= get_xWx(prob_in, resid_in,step_in);  
	  thd_done 	= kappa_test*steps->computeDD();
	}
  }  
	#ifdef TIMING
	gprof.t_computeXSDD2+=MPI_Wtime()-stime0;
	#endif
}


double
FilterIPMSolver::evalBarrObj(Data *prob_in, Variables *iterate_in, Residuals * resid_in, const int isTrialIter)
{	
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenVars*vars = dynamic_cast<NlpGenVars *> (iterate_in);
  
  double barrobj;
  double trial_Obj, BarrObj_tmp;
  if(gUseDualRegAlg == 0){
    if(!isTrialIter) 
	  barrobj = prob->BarrObj;
	else{
      trial_Obj      = prob->objectiveValue(iterate_in);
	  barrobj = prob->BarrObjValue(vars, trial_Obj, FilterIPMOpt->k_d);	  
	}
  }else if(gUseDualRegAlg == 1){
	if(!isTrialIter) 
	  BarrObj_tmp 	 = prob->BarrObj;
	else{
	  trial_Obj   	 = prob->objectiveValue(iterate_in);
	  BarrObj_tmp 	 = prob->BarrObjValue(vars, trial_Obj, FilterIPMOpt->k_d);	 	  
	}
	barrobj 	 = prob->evalMeritFunc(BarrObj_tmp,iterate_in,resid_in);
  }else if(gUseDualRegAlg == 10){
    if(!isTrialIter) 
	  barrobj = prob->BarrObj;
	else{
      trial_Obj 	 = prob->objectiveValue(iterate_in);
      barrobj = prob->BarrObjValue(vars,trial_Obj,FilterIPMOpt->k_d);	
	}
  }
  return barrobj;
}


double
FilterIPMSolver::evalConOneNorm(Data *prob_in, Variables *iterate_in, Residuals * resid_in, const int isTrialIter)
{	
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenVars*vars = dynamic_cast<NlpGenVars *> (iterate_in);
  
  double cons_norm;
  if(gUseDualRegAlg == 0){
    cons_norm = resid->getKKTRhsNorm_Primal(prob,vars,PIPS_NORM_ONE,isTrialIter);
  }else if(gUseDualRegAlg >= 1){
	cons_norm = prob->evalScaledConstraintNorm(iterate_in,resid_in,isTrialIter);
  }
  return cons_norm;
}


double 
FilterIPMSolver::computeSOCIter(Data *prob_in, Variables *vars_in, Residuals * resid_in, const double alp_soc_init)
{
  if(!trialIter_rhs_SOC) trialIter_rhs_SOC = factory->makeResiduals(prob_in);
  if(!SOCstep) SOCstep = factory->makeVariables( prob_in );

  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenResiduals *rhs_SOC = dynamic_cast<NlpGenResiduals *>(trialIter_rhs_SOC);
  
  NlpGenVars *trialIter = dynamic_cast<NlpGenVars *> (trial_iter);
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (vars_in);

  int iter_SOC = 0;
  bool IfAccepted = false;
  
  double ConNorm_old_soc = test_ConNorm;
  double alp_pri_soc;

  // Compute second order corrector rhs
  rhs_SOC->copyFrom(resid); 
  rhs_SOC->updateSOCRhs(alp_soc_init,trialIter, prob);
		
  while (1)
  {
	iter_SOC++;
	
	compute_step_WithRegularization( prob, vars, rhs_SOC, SOCstep,true);
	  
	alp_pri_soc = vars->stepMax_Pri(SOCstep,tau_j);
	
	IfAccepted = ifStepAccepted(prob_in,vars_in,resid_in,SOCstep,trial_iter,alp_pri_soc,true);
	  
	if(IfAccepted==true){
	  if(gdWd_test_soc==1){  	
		if(xWx_done > thd_done){
		  SolfromSOC = true;		
		  break;
		}
	  }else if(gdWd_test_soc==0){
		SolfromSOC = true;
		break;
	  }
	}

	// note that trialIter and test_ConNorm have been updated in ifStepAccepted		
	if( test_ConNorm > FilterIPMOpt->k_soc*ConNorm_old_soc || iter_SOC >= FilterIPMOpt->socIter_max){
	  break;
	}else{ 
      rhs_SOC->updateSOCRhs(alp_pri_soc,trialIter, prob);	
	  ConNorm_old_soc = test_ConNorm;
	}
  }

  return alp_pri_soc;
}


bool 
FilterIPMSolver::ifStepAccepted(Data *prob_in, Variables *vars_in, Residuals * resid_in, Variables *step_in,
					Variables *trial_iter_in, const double trialAlpha, const bool inSOC)
{
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (vars_in);
  NlpGenVars *step = dynamic_cast<NlpGenVars *> (step_in);
	
  NlpGenVars *trialIter = dynamic_cast<NlpGenVars *> (trial_iter_in);

  int isTrialIter=1;
  bool IfInFilter=false, IfAccepted=false;
  double testAlpha = trialAlpha;

  // compute trial point
  trialIter->copy(vars_in);
  trialIter->takePrimalStep(step_in,trialAlpha,trialAlpha);

  // function evaluations at  trial point
  prob->evalConstraintBody(trialIter,isTrialIter);
  test_ConNorm = evalConOneNorm(prob_in,trialIter,resid_in,isTrialIter);
  test_BarrObj = evalBarrObj(prob_in,trialIter,resid_in,isTrialIter);

  // FIXME_NY: Following two lines are new! 2015/09/23.
  // if SOC step, use original step size to do test, which is alp_pri_max. 
  if(inSOC) testAlpha = alp_pri_max;
  
  // test conditions SWC, AC and SRC
  IfAccepted = testStep(testAlpha, test_ConNorm, test_BarrObj,
					  	curr_ConNorm, curr_BarrObj, curr_BarrObjGradTimesD, curr_ConTimesD,
					  	FilterIPMOpt );

  if(IfAccepted==false){
	LastIterRejectByFilter=false;	 
  }else{
    if(gUseFilter==1){  
  	  IfInFilter = Filter->WithinFilter(test_ConNorm,test_BarrObj);
	  // if in filter, means this step is rejected 
	  if(IfInFilter == true) {
	    info_PriAlpha_char_1st = 'r';
		IfAccepted=false;
	    LastIterRejectByFilter = true;
	  }
    }
  }

  if(IfAccepted){
	if(gFilterResetStep>0){
	  if (LastIterRejectByFilter){
		cumulative_reject_filter++;
		if (cumulative_reject_filter>=gFilterResetStep) {
		  if(FilterIPMOpt->ConNorm_max > 0.1*test_ConNorm)
			FilterIPMOpt->ConNorm_max *= 0.1;
		  Filter->Initialize(FilterIPMOpt);
    	  cumulative_reject_filter = 0;
        }
      }
      else {
        cumulative_reject_filter = 0;
      }
	}
	LastIterRejectByFilter=false;

    if(LastIterAcceptByCondition==1)
	  StepAcceptDueTo_SWC_AC++;
    else if(LastIterAcceptByCondition==2)
	  StepAcceptDueTo_SRC_obj++;
    else if(LastIterAcceptByCondition==3)
	  StepAcceptDueTo_SRC_con++;
    else if(LastIterAcceptByCondition==4)
	  StepAcceptDueTo_SRC++;

  }
  
  return IfAccepted;
}


void 
FilterIPMSolver::line_search_Filter( Data *prob_in, Variables *vars_in, Residuals * resid_in, Variables *step_in)
{
  SolfromSOC = false;
  info_PriAlpha_char_1st = '-';
  LastIterAcceptByCondition = -1;
			
  NumBackCheck = BackTrack(prob_in, vars_in,resid_in, step_in);
}


int
FilterIPMSolver::BackTrack( Data *prob_in, Variables *vars_in, Residuals * resid_in, Variables *step_in)
{
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenVars *vars = dynamic_cast<NlpGenVars *> (vars_in);
  NlpGenVars *step = dynamic_cast<NlpGenVars *> (step_in);
  NlpGenVars *trialIter = dynamic_cast<NlpGenVars *> (trial_iter);

  int BackSolveIter = 1, isTrialIter=1;
  double OldAlpha, alp_pri_test, alpha_k_soc;
  bool IfAccepted=false;

  // compute quantities from current point
  curr_BarrObj 			 = evalBarrObj(prob_in,vars_in,resid_in);
  curr_ConNorm 		 	 = evalConOneNorm(prob_in,vars_in,resid_in); 
  curr_BarrObjGradTimesD = prob->getBarrGradTimesD(vars_in,step_in,0,FilterIPMOpt->k_d);
  curr_ConTimesD 	 	 = prob->getConTimesD(vars_in,step_in,resid_in);

  //compute min step size
  double alp_pri_min = getAlphaMin(alp_pri_max);

  OldAlpha 		= alp_pri_max;
  alp_pri_test 	= alp_pri_max;

  while (BackSolveIter >= 0){	 

	// check if trial point is accepted
    IfAccepted = ifStepAccepted(prob_in,vars_in,resid_in,step_in,trial_iter,alp_pri_test);

    // line search is successful! 
    if(IfAccepted==true ) break;

	// check if we need to do second order correction
	if( BackSolveIter == 1 && gDoSOC ){
	  trialIter->copy(vars_in);
	  trialIter->takePrimalStep(step_in,alp_pri_max,alp_pri_max);
      prob->evalConstraintBody(trialIter,isTrialIter);
      test_ConNorm = evalConOneNorm(prob_in,trialIter,resid_in,isTrialIter);

	  if( test_ConNorm >= curr_ConNorm ){
		// Compute second order correction
		alpha_k_soc = computeSOCIter(prob_in, vars_in, resid_in, alp_pri_max);

		// seconde order correction is successful! 
		if (SolfromSOC == true) break;
	  }
	}

	// choose new trial step size
	OldAlpha	 = alp_pri_test;
	alp_pri_test = alp_pri_test/2;	
	
	if( alp_pri_test < alp_pri_min )
	{	
	  DoFeasResto = true;
	  break;
	}										
	if(BackSolveIter>gLineSearchMatStep) break; 
	
	BackSolveIter++;	
  }

  if(SolfromSOC == true)
 	alp_pri = alpha_k_soc;
  else
    alp_pri = alp_pri_test;

  return BackSolveIter;
}


void 
FilterIPMSolver::UpdateFilter()
{
  info_PriAlpha_char = 'f';
  
  if( LastIterAcceptByCondition>1 )			 
  {
    double addConNorm = (1-FilterIPMOpt->r_Q)*curr_ConNorm;
	double addBarrObj = curr_BarrObj-FilterIPMOpt->r_phi*curr_ConNorm;
	Filter->UpdateFilter(addConNorm,addBarrObj);
	if( LastIterAcceptByCondition==2) info_PriAlpha_char = 'o';
	else if( LastIterAcceptByCondition==3) info_PriAlpha_char = 'c';  
	else if( LastIterAcceptByCondition==4) info_PriAlpha_char = 'b';  		
  }
}


void 
FilterIPMSolver::use_TinyStep( Variables *step_in)
{
  NlpGenVars*  steps = (NlpGenVars *) step_in;

  LastIterAcceptByCondition = -1; 	
  SolfromSOC = false;
		
  if ( steps->dual_YZ_InfNorm() < FilterIPMOpt->tiny_step_dual ) {
    DualStepIsTiny = true;
	info_PriAlpha_char = 'T'; 
	cumulative_small_step++;
  }else{
	DualStepIsTiny = false;
	info_PriAlpha_char = ' '; 
	cumulative_small_step = 0;
  }
	
  if(cumulative_small_step == 5){
	DoFeasResto = true;
  }	
}


void FilterIPMSolver::defaultMonitor( Data * /* data */, Variables * /* vars */,
									Residuals * resids,
									double alpha, double sigma,
									int i, double mu, 
									int status_code,
									int level )
{
  assert("not done" && 0);
}


int 
FilterIPMSolver::defaultStatus(Data *  data_in, Variables * /* vars */,
				 Residuals * resids_in,
				 int iterate, double mu, 
				 int /* level */)
{
  NlpGenResiduals* resids = dynamic_cast<NlpGenResiduals*> (resids_in);
  NlpGenData* prob = (NlpGenData*) data_in;
  NlpGenVars*  steps = (NlpGenVars *) mainstep;

  if(DoFeasResto){
  	printf("OOPS! Feasiblity Restoration Phase. STOP! \n");
    return NEED_FEASIBILITY_RESTORATION;
  }

  int stop_code = NOT_FINISHED;
  int idx;

  double pObj  = prob->PriObj;

  double fullErr;
  double pErr = resids->priErr();
  double dErr = resids->dualErr();
  double err_Comple = resids->comp_Err();  
  double err_Comple_0 = resids->comp_Err_0();  

  int numberOfPrimalReg = PD_Reg->num_PrimReg;
  
  fullErr = MAX(pErr,dErr/scal_dualerr);
  fullErr = MAX(fullErr,err_Comple_0/scal_commerr);

  double pri_step_InfNorm = steps->primal_XS_InfNorm();

  deltaRegular_pr = PD_Reg->prim_reg_curr;
  if (deltaRegular_pr > 0){
	snprintf(strbuffer, sizeof(strbuffer), "  %1.1f", log10(deltaRegular_pr));
	strPrt = strbuffer;
  }else{
	strPrt = dashes;
  }
  
  deltaRegular_du = PD_Reg->dual_reg_curr;
  if (deltaRegular_du > 0){
	snprintf(strbuffer2, sizeof(strbuffer2), "	%1.1f", log10(deltaRegular_du));
	strPrt_du = strbuffer2;
  }else{
	strPrt_du = dashes;
  } 

  double LinSysErr 	= prob->linsysRes;
  double LinSysErr_KKT = prob->linsysRes_Full;

  idx = iterate-1;
  if(idx <  0     ) 
  	idx=0;
  if(idx >= maxit ) 
  	idx=maxit-1;

  // store the historical record
  mu_history[idx] = mu;
  rnorm_history[idx] = fullErr;

  if(printlevel>0 && (0==(iterate-1)%10) ){
    printf(" Iter\t  Objective     Inf_pr\t  Inf_du     Comp    lg(mu)   "
			"||d||   lg(PrRg)\t alpha_du   alpha_pr    LS   lg(DuRg)\tAugSysErr KKTSysErr  "
			"  dWd0      dWd       thd0       thd      #RegFact   #KryIt\n");
  }
  if(printlevel>0 && (1==iterate) ){
	printf(" %d\t% 1.7e  %1.2e  %1.2e  %1.2e  %1.1f\n",
		iterate-1,pObj,pErr,dErr,err_Comple,log10(mu)); 
  }else{
    if(printlevel>0){
	  printf(" %d\t% 1.7e  %1.2e  %1.2e  %1.2e  %1.1f"  
	      "  %1.2e   %s\t %1.2e  %1.2e%c%c   %d  %s\t%1.2e  %1.2e"
	      "  % 9.2e%c % 9.2e % 9.2e % 9.2e      %d        %d\n",
		  iterate-1,pObj,pErr,dErr,err_Comple,log10(mu),
		  pri_step_InfNorm,strPrt, alp_boundDual,alp_pri,info_PriAlpha_char_1st,info_PriAlpha_char,NumBackCheck, 
		  strPrt_du,  LinSysErr, LinSysErr_KKT, 
		  xWx_0,info_xWx_type,xWx_done,thd_0,thd_done, numberOfPrimalReg,prob->KryIter);  
	}
  }

  if ( fullErr <= FilterIPMOpt->opt_tol){ 
	stop_code = SUCCESSFUL_TERMINATION;
	if(printlevel>0){
	  printf("\n\n  Find Optimal solution! In Iter: %d",iterate-1);
	  printf("\n  Optimal solution is: %1.7e",pObj);
	  printf("\n  Addition Fact due to reg: %d \n",numberOfPrimalReg);
	}
  } else if (
    iterate-1 >= maxit ) {
    stop_code = MAX_ITS_EXCEEDED;
	if(printlevel>0){
	  printf("\n\n  EXIT: Max Iter = %d \n", FilterIPMOpt->maxit);
	  printf("\n  Last objective is: %1.7e",pObj);
	  printf("\n  Addition Fact due to reg: %d \n",numberOfPrimalReg);	
	}
  } else if ( pObj/FilterIPMOpt->df <= -1e+15) { 
	stop_code = UNKNOWN;
	if(printlevel>0){
	  printf("\n\n  EXIT: Iterations Diverging. Problem might be unbounded.\n");
	}	
  } else if ( pObj/FilterIPMOpt->df >= 1e+15   ) { 
    stop_code = UNKNOWN;
	if(printlevel>0){
	  printf("\n\n  EXIT: Iterations Diverging.\n");
	}			
  }

  if ( stop_code != NOT_FINISHED){ 
  	if(printlevel>0){
		printf("  Iter is accepted due to: \n");
		printf("							SWC and AC: %d \n", StepAcceptDueTo_SWC_AC);
		printf("							SRC __ Obj: %d \n", StepAcceptDueTo_SRC_obj);
		printf("							SRC __ Con: %d \n", StepAcceptDueTo_SRC_con);		
		printf("							SRC (both): %d \n", StepAcceptDueTo_SRC);	

	}
  }

  return stop_code;
  
}


double*  
FilterIPMSolver::computeQuantitiesForDualReg( Data *prob_in, Variables *iter_in, Residuals * resid_in, Variables *step_in,
												 PDRegularization* RegInfo)
{
  NlpGenResiduals *resid = dynamic_cast<NlpGenResiduals *>(resid_in);
  NlpGenData *prob = dynamic_cast<NlpGenData *> (prob_in);
  NlpGenVars *steps = dynamic_cast<NlpGenVars *> (step_in);
  NlpGenVars *iters = dynamic_cast<NlpGenVars *> (iter_in);
  NlpGenLinsys *sysNLP = dynamic_cast<NlpGenLinsys *>(sys);

  sysNLP->computeQuantitiesForDualReg(prob,iters,resid,steps,RegInfo->quantitiesForReg); 

  return dualRegQuantities;
}


int FilterIPMSolver::solve( Data *prob_in, Variables *iterate, Residuals * resid )
{
  NlpGenData * prob = dynamic_cast<NlpGenData *>(prob_in);
  NlpGenVars * vars = dynamic_cast<NlpGenVars *>(iterate);
  NlpGenVars*  steps = (NlpGenVars *) mainstep;

  int done, muChanged;
  int status_code;

  ChangeBounds(prob,resid);
  
  dnorm = prob->datanorm();
  
  // initialization of linear system
  sys = factory->makeLinsys( prob );

  // compute initial point 
  this->start( factory, iterate, prob, resid, steps );

  iter = 0;
  done = 0;
  mu =  FilterIPMOpt->mu_init;
  prob->currMu = mu;
  tau_j = (1-mu)>FilterIPMOpt->tau_Min ? 1-mu : FilterIPMOpt->tau_Min;

  /* ----------------------------- Start main loop ------------------------- */
  do
  {
	iter++;	

	EvalErrScaling(iterate);

    // evaluate functions (obj, con, jac, hes)
	prob->evalData(iterate);
	// update BarrObjValue with damping rate
	prob->BarrObj = prob->BarrObjValue(vars, prob->PriObj, FilterIPMOpt->k_d);

    // evaluate residuals
    resid->calcresids(prob, iterate);
	
	// initialize the filter  and some parameters
    if(iter==1 )
	  FilterInitializeAndPara(prob, iterate,resid);

    //  termination test:
    status_code = this->doStatus( prob, iterate, resid, iter, mu, 0 );
	if( status_code != NOT_FINISHED ) 
	  break;

	// update barrier parameter mu, reset filter if necessary
	muChanged=updateBarrierParameter(prob, iterate, resid);

    //add damping term to where var/con only has single bound, see filter line search IPM paper section 3.7 
    addDampingTermToKKT(resid);

    // compute search direction and update regularization if necessary
	compute_step_WithRegularization( prob, iterate, resid,steps);

	//apply frac-to-boundary to get max primal step length
	alp_pri_max = vars->stepMax_Pri(steps,tau_j);  

	// check if full step is too tiny
	IfTinyStep = isTinyStep(iterate,resid,steps);

	NumBackCheck = 0;
	if(!(IfTinyStep && gDoTinyStepTest==1) ){
	  // backtracking line search
	  line_search_Filter(prob,iterate,resid,steps);  
	}else{
	  // use tiny step
	  use_TinyStep(steps);
	}

	// terminate the process due to calling feasibility restoration phase
	if(DoFeasResto){
	  status_code = this->doStatus( prob, iterate, resid, iter, mu, 0 );
	  break;
	}
	
	if(SolfromSOC == true)
	{
	  steps->copy(SOCstep);
	}

	//apply frac-to-boundary to get max dual step length
	alp_boundDual_max = vars->stepMax_BoundDual(steps,tau_j);  

	// accept step
    accept_iterate(prob,iterate,resid,steps);
  
	// update filter if necessary 
	if(gUseFilter==1 && !(IfTinyStep && gDoTinyStepTest==1))
	  UpdateFilter();

	if(SolfromSOC == true)
	{
	  info_PriAlpha_char = (char)toupper(info_PriAlpha_char);
	}  
  }while(!done);

  return status_code;
}

