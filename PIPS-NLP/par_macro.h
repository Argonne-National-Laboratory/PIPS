#ifndef PAR_MACRO_H_
#define PAR_MACRO_H_

#include <string>
#include <sstream>
#include <iostream>

extern int gmyid;
extern int gnprocs;

#ifdef FQ_DEBUG
#define DEBUG_ENABLE 1
#else
#define DEBUG_ENABLE 0
#endif

#define DEBUG(X) do { \
	if (DEBUG_ENABLE) { X } \
} while (0)

#define PAR_DEBUG(x) do { \
  if (DEBUG_ENABLE) { std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< x << std::endl; } \
} while (0)


template <class T>
void print_array(const std::string& msg, T* data, size_t len)
{
	std::ostringstream oss;
	for(size_t i=0;i<len;i++){
		oss<<data[i]<<", ";
	}
	PAR_DEBUG(""<<msg<<" - "<<"Array [ "<<oss.str()<<"  ]");
}

#define PRINT_ARRAY(M, DATA, LEN) do { \
  if (0) { 	std::ostringstream oss; 	\
				for(size_t i=0;i<LEN;i++){ 	\
					oss<<DATA[i]<<", ";    	\
				} \
				std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] "<< M << "Array [ " <<oss.str() <<" ]"<< std::endl; \
			} \
} while (0)


extern void convert_to_csr(int m, int n, int* rowidx, int* colptr, double* elts, int nz, double* ret);

#ifdef PROF
struct Profile
{
	Profile():
		t_struct_building(0.0), t_model_evaluation(0.0), t_solver_lifetime(0.0),
		n_feval(0),n_eval_g(0),n_grad_f(0), n_jac_g(0), n_laghess(0)
	{ }

	double t_struct_building;
	double t_model_evaluation;
	double t_solver_lifetime;

	int n_prob_info;
	int n_init_x0;
	int n_feval;
	int n_eval_g;
	int n_grad_f;
	int n_jac_g;
	int n_laghess;
	int n_write_solution;
};

extern struct Profile gprof;
extern std::ostream& operator<<(std::ostream& os, const Profile& prof);
extern void report_timing(Profile& p);

#endif

#endif
