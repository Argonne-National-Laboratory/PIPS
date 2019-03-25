/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include "FilterIPMStochSolver.h"

#include "Status.h"
#include "Data.h"
#include "ProblemFormulation.h"

#include "OoqpVector.h"
#include "DoubleMatrix.h"

#include "StochTree.h"
//#include "NlpGenStoch.h"
#include "StochResourcesMonitor.h"

#include "sFactory.h"


#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;

#include <cstdio>
#include <cassert>
#include <cmath>

//#include "StochVector.h"
//#include "mpi.h"
#include "NlpGenVars.h"
#include "NlpGenData.h"
#include "NlpGenResiduals.h"


#include "FilterIPMOption.h"
#include "PDRegularization.h"

#include "global_var.h"

#ifdef TIMING
#include "../PIPS-NLP/Core/Utilities/PerfMetrics.h"
#endif

// gmu is used in EmtlStochSymIndefSolver to decide if the system should be perturbed
double gmu;
// double grnorm;
extern int gOoqpPrintLevel;
extern int gUseFilter;
extern int gDoTinyStepTest;


double g_iterNumber;


static int sleepFlag=0;




#ifndef MIN
#define MIN(a,b) ((a > b) ? b : a)
#endif

#ifndef MAX
#define MAX(a,b) ((a > b) ? a : b)
#endif



/** 
 * @ingroup NlpSolvers
 */
FilterIPMStochSolver::FilterIPMStochSolver( ProblemFormulation * of, Data * prob )
	: FilterIPMSolver(of, prob)

{}

FilterIPMStochSolver::~FilterIPMStochSolver()
{}




int FilterIPMStochSolver::solve( Data *prob_in, Variables *iterate, Residuals * resid )
{
  NlpGenData * prob = dynamic_cast<NlpGenData *>(prob_in);
  NlpGenVars * vars = dynamic_cast<NlpGenVars *>(iterate);
  NlpGenVars*  steps = (NlpGenVars *) mainstep;

  int done, muChanged;
  int status_code;

  sFactory* stochFactory = reinterpret_cast<sFactory*>(factory);

  gmu = 1000;
  g_iterNumber=0.0;
  giterNum = 0;

  ChangeBounds(prob,resid);
  
  dnorm = prob->datanorm();

  // initialization of linear system
  sys = factory->makeLinsys( prob );

  // compute initial point 
  stochFactory->iterateStarted();
  this->start( factory, iterate, prob, resid, steps );
  stochFactory->iterateEnded();

  iter = 0;
  done = 0;
  mu =  FilterIPMOpt->mu_init;
  prob->currMu = mu;
  tau_j = (1-mu)>FilterIPMOpt->tau_Min ? 1-mu : FilterIPMOpt->tau_Min;

  gmu = mu;

  /* ----------------------------- Start main loop ------------------------- */
  do
  {
	iter++; g_iterNumber=iter; giterNum = iter;
#ifdef TIMING
  double stime=MPI_Wtime();
	double stimet=stime;
#endif
    stochFactory->iterateStarted();

	EvalErrScaling(iterate);

    // evaluate functions (obj, con, jac, hes)
	prob->evalData(iterate);
#ifdef TIMING
	gprof.t_evalData+=MPI_Wtime()-stime;
#endif
	// update BarrObjValue with damping rate
#ifdef TIMING
	stime=MPI_Wtime();
#endif
	prob->BarrObj = prob->BarrObjValue(vars, prob->PriObj, FilterIPMOpt->k_d);
#ifdef TIMING
	gprof.t_BarrObj+=MPI_Wtime()-stime;
#endif


    // evaluate residuals
#ifdef TIMING
		stime=MPI_Wtime();
#endif
    resid->calcresids(prob, iterate);
#ifdef TIMING
		gprof.t_calcresids+=MPI_Wtime()-stime;
#endif

	// initialize the filter  and some parameters
#ifdef TIMING
	  stime=MPI_Wtime();
#endif
    if(iter==1 )
	  FilterInitializeAndPara(prob, iterate,resid);

    //  termination test:
    status_code = this->doStatus( prob, iterate, resid, iter, mu, 0 );
  //  if(iter==10) status_code = SUCCESSFUL_TERMINATION;
	if( status_code != NOT_FINISHED )
	  break;

	// update barrier parameter mu, reset filter if necessary
	muChanged=updateBarrierParameter(prob, iterate, resid);
#ifdef TIMING
	gprof.t_updateBarrierParameter+=MPI_Wtime()-stime;
#endif

    //add damping term to where var/con only has single bound, see filter line search IPM paper section 3.7
#ifdef TIMING
		stime=MPI_Wtime();
#endif
    addDampingTermToKKT(resid);
#ifdef TIMING
		gprof.t_addDampingTermToKKT+=MPI_Wtime()-stime;
#endif

    // compute search direction and update regularization if necessary
#ifdef TIMING
	stime=MPI_Wtime();
#endif
	compute_step_WithRegularization( prob, iterate, resid,steps);
#ifdef TIMING
	gprof.t_compute_step_WithRegularization+=MPI_Wtime()-stime;
  stime=MPI_Wtime();
#endif
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
#ifdef TIMING
 gprof.t_line_search+=MPI_Wtime()-stime;
 stime=MPI_Wtime();
#endif
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

    gmu = mu;
#ifdef TIMING
		gprof.t_rest+=MPI_Wtime()-stime;
		gprof.t_total+=MPI_Wtime()-stimet;
#endif
    stochFactory->iterateEnded();
  }while(!done);

  return status_code;
}


int 
FilterIPMStochSolver::defaultStatus(Data *  data_in, Variables * /* vars */,
				 Residuals * resids_in,
				 int iterate, double mu, 
				 int /* level */)
{
  NlpGenResiduals* resids = dynamic_cast<NlpGenResiduals*> (resids_in);
  NlpGenData* prob = (NlpGenData*) data_in;
  NlpGenVars*  steps = (NlpGenVars *) mainstep;

  if(DoFeasResto){
    //printf("OOPS! Feasiblity Restoration Phase. STOP! \n");
    int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    if(0==myRank)
      printf("Problem may be locally infeasible. Feasibility Restoration Phase is needed, but it is not yet implemented.\n");
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


  int myRank; MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  if(0==myRank){
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
		  pri_step_InfNorm, strPrt, alp_boundDual,alp_pri,info_PriAlpha_char_1st,info_PriAlpha_char,NumBackCheck, 
		  strPrt_du,  LinSysErr, LinSysErr_KKT, 
		  xWx_0,info_xWx_type,xWx_done,thd_0,thd_done, numberOfPrimalReg,prob->KryIter);  
	  }
    }
  }

  if ( fullErr <= FilterIPMOpt->opt_tol){ 
	stop_code = SUCCESSFUL_TERMINATION;
	if(0==myRank && printlevel>0){
	  printf("\n\n  Found optimal solution! In Iter: %d",iterate-1);
	  printf("\n  Optimal solution is: %1.7e",pObj);
	  printf("\n  Addition fact due to reg.: %d \n",numberOfPrimalReg);
	}
  } else if (
    iterate-1 >= maxit ) {
    stop_code = MAX_ITS_EXCEEDED;
	if(0==myRank && printlevel>0){
	  printf("\n\n  EXIT: Max Iter = %d \n", FilterIPMOpt->maxit);
	  printf("\n  Last objective is: %1.7e",pObj);
	  printf("\n  Addition Fact due to reg: %d \n",numberOfPrimalReg);	
	}
  } else if ( pObj/FilterIPMOpt->df <= -1e+15) { 
	stop_code = UNKNOWN;
	if(0==myRank && printlevel>0){
	  printf("\n\n  EXIT: Iterations Diverging. Problem might be unbounded.\n");
	}	
  } else if ( pObj/FilterIPMOpt->df >= 1e+15   ) { 
    stop_code = UNKNOWN;
	if(0==myRank && printlevel>0){
	  printf("\n\n  EXIT: Iterations Diverging.\n");
	}			
  }

  if ( stop_code != NOT_FINISHED){ 
  	if(0==myRank && printlevel>0){
		printf("  Iter is accepted due to: \n");
		printf("							SWC and AC: %d \n", StepAcceptDueTo_SWC_AC);
		printf("							SRC __ Obj: %d \n", StepAcceptDueTo_SRC_obj);
		printf("							SRC __ Con: %d \n", StepAcceptDueTo_SRC_con);		
		printf("							SRC (both): %d \n", StepAcceptDueTo_SRC);	

	}
  }

  
  return stop_code;
  
}

