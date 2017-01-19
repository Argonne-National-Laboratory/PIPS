#ifdef NLPTIMING

#include "PerfMetrics.h"
#include "../../global_var.h"
#include <sys/stat.h>
#include <sstream>
#include <iostream>
#include <fstream>


PerfMetrics* PerfMetrics::_prof = NULL;
PerfMetrics& PerfMetrics::getPerfMetrics()
{
	if(PerfMetrics::_prof == NULL)
		PerfMetrics::_prof = new PerfMetrics();
	return *(PerfMetrics::_prof);
}

std::ostream& operator<<(std::ostream& os, const PerfMetrics& p)
{
	return os<<"t_struct_building "<<p.t_struct_building<< std::endl
			<<"t_model_evaluation "<<p.t_model_evaluation<< std::endl
			<<"t_solver_lifetime "<<p.t_solver_lifetime<< std::endl
			<<"t_solver_go "<<p.t_solver_go<< std::endl
			<<"------------------"<< std::endl
			<<"In solver"<< std::endl
			<<"------------------"<< std::endl
			<<"t_compute_step_WithRegularization "<<p.t_compute_step_WithRegularization << std::endl
			<<"t_evalData "<<p.t_evalData << std::endl
			<<"t_BarrObj "<<p.t_BarrObj << std::endl
			<<"t_calcresids "<<p.t_calcresids << std::endl
		  <<"t_updateBarrierParameter "<<p.t_updateBarrierParameter << std::endl
			<<"t_addDampingTermToKKT "<<p.t_addDampingTermToKKT << std::endl
			<<"t_line_search "<<p.t_line_search << std::endl
			<<"t_rest "<<p.t_rest << std::endl
			<<"t_total "<<p.t_total << std::endl
			<<"t_total_sum "<<p.t_compute_step_WithRegularization+p.t_evalData+
			p.t_BarrObj+p.t_calcresids+p.t_updateBarrierParameter
			+p.t_addDampingTermToKKT+p.t_line_search+p.t_rest << std::endl
			<<"------------------"<< std::endl
			<<"In compute_step_WithRegularization"<< std::endl
			<<"------------------"<< std::endl
			<<"t_set_r3_xz_alpha "<< p.t_set_r3_xz_alpha << std::endl
			<<"t_factor "<< p.t_factor << std::endl
			<<"t_computeXSDD1 "<< p.t_computeXSDD1 << std::endl
			<<"t_computeQuantitiesForDualRegp "<< p.t_computeQuantitiesForDualReg << std::endl
			<<"t_computeXSDD2 "<< p.t_computeXSDD2 << std::endl
			<<"------------------"<< std::endl
			<<"In factor()"<< std::endl
			<<"------------------"<< std::endl
			<<"t_factorNoMatChange "<< p.t_factorNoMatChange << std::endl
			<<"t_factor_rest "<< p.t_factor_rest << std::endl
			<<"t_computeRegularization "<< p.t_computeRegularization << std::endl
			<<"t_factorNoMatChange2 "<< p.t_factorNoMatChange2 << std::endl
			<<"t_factor2 "<< p.t_factor2 << std::endl
			<<"------------------"<< std::endl
			<<"In factor2()"<< std::endl
			<<"------------------"<< std::endl
			<<"t_initializeKKT "<< p.t_initializeKKT << std::endl
			<<"t_reduceKKT "<< p.t_reduceKKT << std::endl
			<<"t_finalizeKKT "<< p.t_finalizeKKT << std::endl
			<<"t_factorizeKKT "<< p.t_factorizeKKT << std::endl
			<<"t_factor2_total "<< p.t_factor2_total << std::endl
			<<"------------------"<< std::endl
			<<"In factorizeKKT()"<< std::endl
			<<"------------------"<< std::endl
			<<"t_matrixChanged "<< p.t_matrixChanged << std::endl
			<<"------------------"<< std::endl
			<<"n_prob_info "<<p.n_prob_info<< std::endl
			<<"n_init_x0 "<<p.n_init_x0<< std::endl
			<<"n_feval "<<p.n_feval<< std::endl
			<<"n_eval_g "<<p.n_eval_g<< std::endl
			<<"n_grad_f "<<p.n_grad_f<< std::endl
			<<"n_jac_g "<<p.n_jac_g<< std::endl
			<<"n_laghess "<<p.n_laghess<< std::endl
			<<"n_write_solution "<<p.n_write_solution<<" ";
}

void PerfMetrics::report_timing()
{
	std::ostringstream oss;
	oss<<"["<<gmyid<<"/"<<gnprocs<<"] " << " [ " <<(*this)<<" ]"<< std::endl;

	mkdir("./out", S_IRWXU);
	std::fstream fs;
	std::ostringstream fn;
	fn<<"./out/"<<gnprocs<<"."<<gmyid<<".c.txt";
	fs.open(fn.str().c_str(),std::fstream::out);
	fs<<oss.str();
	fs.close();

	if(gmyid == 0)
		std::cout<<oss.str();

}

#endif
