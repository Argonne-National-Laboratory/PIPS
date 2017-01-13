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
