#ifndef _PERFMETRICS_H_
#define _PERFMETRICS_H_

#ifdef NLPTIMING


#include <string>
#include <sstream>
#include <iostream>


class PerfMetrics
{

private:
	PerfMetrics():
		t_struct_building(0.0), t_model_evaluation(0.0), t_solver_lifetime(0.0),
		t_solver_go(0.0), t_compute_step_WithRegularization(0.0), t_evalData(0.0),
		t_BarrObj(0.0), t_calcresids(0.0), t_updateBarrierParameter(0.0),
		t_addDampingTermToKKT(0.0), t_line_search(0.0), t_rest(0.0), t_total(0.0),
		t_set_r3_xz_alpha(0.0), t_factor(0.0), t_computeXSDD1(0.0),
		t_computeQuantitiesForDualReg(0.0), t_computeXSDD2(0.0), 
		t_factorNoMatChange(0.0), t_computeRegularization(0.0),
		t_factorNoMatChange2(0.0), t_factor2(0.0), t_factor2_total(0.0),
		t_factor_rest(0.0), t_matrixChanged(0.0), t_ma57ad(0.0), t_ma57bd(0.0),
		t_initializeKKT(0.0), t_reduceKKT(0.0), t_factorizeKKT(0.0), t_finalizeKKT(0.0),
		n_prob_info(0), n_init_x0(0), n_feval(0), n_eval_g(0), n_grad_f(0),
		n_jac_g(0), n_laghess(0), n_write_solution(0)
	{
		PerfMetrics::_prof = NULL;
	};

public:
	static PerfMetrics& getPerfMetrics();
	~PerfMetrics() {};

	static PerfMetrics* _prof;

	double t_struct_building;
	double t_model_evaluation;
	double t_solver_lifetime;
	double t_solver_go;
	double t_compute_step_WithRegularization;
	double t_evalData;
	double t_BarrObj;
	double t_calcresids;
  double t_updateBarrierParameter;
	double t_addDampingTermToKKT;
	double t_line_search;
	double t_rest;
	double t_total;
	double t_set_r3_xz_alpha;
	double t_factor;
	double t_computeXSDD1;
	double t_computeQuantitiesForDualReg;
	double t_computeXSDD2;
	double t_factorNoMatChange;
	double t_computeRegularization;
	double t_factorNoMatChange2;
	double t_factor2;
	double t_initializeKKT;
	double t_reduceKKT;
	double t_factorizeKKT;
	double t_finalizeKKT;
	double t_factor_rest;
	double t_factor2_total;
	double t_matrixChanged;
	double t_ma57ad;
	double t_ma57bd;
	int n_prob_info;
	int n_init_x0;
	int n_feval;
	int n_eval_g;
	int n_grad_f;
	int n_jac_g;
	int n_laghess;
	int n_write_solution;

public:
	friend std::ostream& operator<<(std::ostream& os, const PerfMetrics& p);
	void report_timing();
};


#endif

#endif
