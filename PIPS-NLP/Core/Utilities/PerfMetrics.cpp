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
	return os<<"t_struct_building "<<p.t_struct_building<<" "
			<<"t_model_evaluation "<<p.t_model_evaluation<<" "
			<<"t_solver_lifetime "<<p.t_solver_lifetime<<" "
			<<"n_prob_info "<<p.n_prob_info<<" "
			<<"n_init_x0 "<<p.n_init_x0<<" "
			<<"n_feval "<<p.n_feval<<" "
			<<"n_eval_g "<<p.n_eval_g<<" "
			<<"n_grad_f "<<p.n_grad_f<<" "
			<<"n_jac_g "<<p.n_jac_g<<" "
			<<"n_laghess "<<p.n_laghess<<" "
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
