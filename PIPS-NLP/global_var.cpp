
#include "par_macro.h"
#include "CoinPackedMatrix.hpp"
#include <sstream>
#include <iostream>

int gmyid;
int gnprocs;
Profile gprof;

void convert_to_csr(int m, int n, int* rowidx, int* colptr, double* elts, int nz, double* ret)
{
	if(nz!=0){
		PRINT_ARRAY("rowidx",rowidx,nz);
		PRINT_ARRAY("colptr",colptr,n+1);
		PRINT_ARRAY("elts",elts,nz);
		CoinPackedMatrix mat;
		mat.copyOf(true,m,n,nz,&elts[0],&rowidx[0],&colptr[0],0);
		assert(nz == mat.getNumElements());
//		mat.dumpMatrix();
		mat.reverseOrdering();
//		mat.dumpMatrix();
		const double* csr_elts = mat.getElements();
		PRINT_ARRAY("csr_elts",csr_elts,nz);
		for(int i=0;i<nz;i++){
			ret[i] = csr_elts[i];
		}
	}
}

#ifdef PROF

std::ostream& operator<<(std::ostream& os, const Profile& prof)
{
	return os<<"t_struct_building "<<prof.t_struct_building<<" "
			<<"t_model_evaluation "<<prof.t_model_evaluation<<" "
			<<"t_solver_lifetime "<<prof.t_solver_lifetime<<" "
			<<"n_prob_info "<<prof.n_prob_info<<" "
			<<"n_init_x0 "<<prof.n_init_x0<<" "
			<<"n_feval "<<prof.n_feval<<" "
			<<"n_eval_g "<<prof.n_eval_g<<" "
			<<"n_grad_f "<<prof.n_grad_f<<" "
			<<"n_jac_g "<<prof.n_jac_g<<" "
			<<"n_laghess "<<prof.n_laghess<<" "
			<<"n_write_solution "<<prof.n_write_solution<<" ";
}

void report_timing(Profile& p)
{
	std::cout<<"["<<gmyid<<"/"<<gnprocs<<"] " << " [ " <<p<<" ]"<< std::endl;
}

#endif
