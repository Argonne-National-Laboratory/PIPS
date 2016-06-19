/* PIPS-NLP                                             *
 * Authors: Feng Qiang and Cosmin G. Petra              *
 * (C) 2016 Argonne National Laboratory			*/

#include <cstdlib>

#include "parallelPipsNlp_C_Callback.h"

#include "../Core/NlpInfo/StructJuMPsInfo.h"
#include "../../Input/StructJuMPInput.h"
#include "../Core/NlpStoch/NlpPIPSIpmInterface.h"
#include "../Core/NlpStoch/sFactoryAug.h"
#include "../Core/NlpSolvers/FilterIPMStochSolver.h"


#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"
#include "NlpInfoCallBack.h"


#include "FilterIPMSolver.h"


#include "NlpGenSparseWithSolver.h"
#include "cNlpGenSparseNLP.h"

#include "Status.h"

#include "pipsOptions.h"

#include "../global_var.h"
#include "../PIPS-NLP/Core/Utilities/PerfMetrics.h"

extern "C"
PipsNlpProblemStructPtr 
CreatePipsNlpProblemStruct(
	MPI_Comm comm,
	int nscen,
	str_init_x0_cb init_x0,
	str_prob_info_cb prob_info,
	str_eval_f_cb eval_f,
	str_eval_g_cb eval_g,
	str_eval_grad_f_cb eval_grad_f,
	str_eval_jac_g_cb eval_jac_g,
	str_eval_h_cb eval_h,
	str_write_solution_cb write_solution,
	UserDataPtr userdata,
	str_get_link_matrix_cb get_link_matrix,
	str_link_info_cb link_info)
{
  MPI_Comm_rank(comm, &gmyid);
  MPI_Comm_size(comm, &gnprocs);
  //MESSAGE("on proc ["<<gmyid<<"] of ["<< gnprocs << "] MPI processes.");
  MESSAGE("CreatePipsNlpProblemStruct - C");

	PipsNlpProblemStructPtr retval = new PipsNlpProblemStruct;

	retval->comm = comm;
	retval->nscen = nscen;
	retval->init_x0 = init_x0;
	retval->prob_info = prob_info;
	retval->eval_f = eval_f;
	retval->eval_g = eval_g;
	retval->eval_grad_f = eval_grad_f;
	retval->eval_jac_g = eval_jac_g;
	retval->eval_h = eval_h;
	retval->write_solution = write_solution;
	retval->userdata = userdata;
	retval->get_link_matrix = get_link_matrix;
	retval->link_info = link_info;
	retval->objective = 0.0;
	retval->nvars = 0;
	retval->ncons = 0;
	return retval;
}


extern int gInnerSCsolve;
extern int gNP_Alg;
extern int gAddSlackParallelSetting;
extern int gSymLinearSolver;
extern int gUseReducedSpace;

extern "C"
int PipsNlpSolveStruct(PipsNlpProblemStruct* prob)
{
#ifdef NLPTIMING
  double stime = MPI_Wtime();
#endif
  MESSAGE("PipsNlpSolveStruct - for " << gnprocs << " processors");
  MPI_Comm comm = prob->comm;
  
  pipsOptions *pipsOpt = new pipsOptions();
  pipsOpt->readFile();
  pipsOpt->defGloOpt();
  gInnerSCsolve=0;
  
  StructJuMPInput *s = new StructJuMPInput(prob);
  MESSAGE("comm is "<<comm);
  assert(comm == MPI_COMM_WORLD);
  
  MESSAGE("before PIPSIpmInterface created .." );
  NlpPIPSIpmInterface<sFactoryAug, FilterIPMStochSolver, StructJuMPsInfo> pipsIpm(*s,comm);
  MESSAGE("PIPSIpmInterface created .." );
  
  if (gmyid == 0) {
    std::cout << "  \n  -----------------------------------------------\n"
	      << "  NLP Solver \n"
	      << "  Argonne National Laboratory, 2016\n"
	      << "  -----------------------------------------------\n" <<std::endl;
    
    if(gUseReducedSpace>0) {
      std::cout << "\n  Reduced Space Solver ------	 ";
      std::cout << "Reduced Space Solver with Umfpack and following linear solver.\n";
    }
    if(0==gSymLinearSolver)
      std::cout << "\n  Linear system solver ------	 Ma27.\n\n";
    else if(1==gSymLinearSolver)
      std::cout << "\n  Linear system solver ------	 Ma57.\n\n";
    else if(2==gSymLinearSolver)
      std::cout << "\n  Linear system solver ------	 Pardiso.\n\n";
    else if(3==gSymLinearSolver)
      std::cout << "\n  Linear system solver ------	 Umfpack.\n\n";
  }
  
  pipsIpm.computeProblemSize(prob->nvars,prob->ncons);

  int ret = pipsIpm.go();

  prob->objective = pipsIpm.getObjective();
  
  delete pipsOpt;
  delete s;

#ifdef NLPTIMING
  gprof.t_solver_lifetime = MPI_Wtime() - stime;
  gprof.report_timing();
#endif

  return ret;
}


extern "C"
void FreePipsNlpProblemStruct(PipsNlpProblemStruct* prob){
  //! shouldn't the signature be PipsNlpProblemStruct**
  delete prob;
}

extern "C"
int get_x(CallBackDataPtr data,double* x, double* lam_eq, double* lam_ieq)
{
  //to be implemented
  return 0;
}

extern "C"
double PipsNlpProblemStructGetObjective(PipsNlpProblemStruct* prob)
{
  if(prob)
    return prob->objective;
  else
    return 0.0;
}

extern "C"
int PipsNlpProblemStructGetTotalVars(PipsNlpProblemStruct* prob)
{
  if(prob)
      return prob->nvars;
    else
      return 0;
}

extern "C"
int PipsNlpProblemStructGetTotalCons(PipsNlpProblemStruct* prob)
{
  if(prob)
      return prob->ncons;
    else
      return 0;
}

