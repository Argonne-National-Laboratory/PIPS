#include "./Drivers/parallelPipsNlp_C_Callback.h"

#include "mpi.h"
#include "par_macro.h"
#include <iostream>
#include <cassert>
#include <math.h>

//# min (x1+x2)^2+(x1+x2)*x3 + (x1+x2)*x4
//# st.
//#     x1 * x2 = 10
//#     x2^2 + x3*x1 < 5
//#     x2^2 + x4*x1 < 6
//# x1, x2 , x3, x4 free variables


int str_init_x0(double* x0, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_init_x0 -- row " << row <<" col "<<col);
	assert(row == col);
	if(row == 0)
	{
		x0[0] = 1.0;
		x0[1] = 1.0;
	}
	else if(row == 1)
	{
		x0[0] = 1.0;
	}
	else if(row == 2)
	{
		x0[0] = 1.0;
	}
	else
		assert(false);

	return 1;
}

int str_prob_info(int* n, double* col_lb, double* col_ub, int* m,
		double* row_lb, double* row_ub, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_prob_info -- row " << row <<" col "<<col);
	assert(row == col);
	if(col_lb == NULL)
	{
		assert(row_lb == NULL);
		assert(row_ub == NULL);
		assert(col_ub == NULL);
		if(row==0)
		{
			*n = 2;
			*m = 1;
		}
		else if(row ==1 || row == 2)
		{
			*n = 1;
			*m = 1;
		}
		else
			assert(false);
	}
	else
	{
		if(row==0)
		{
			assert(*n==2 && *m == 1);
			col_lb[0] = -INFINITY;
			col_lb[1] = -INFINITY;
			col_ub[0] = INFINITY;
			col_ub[1] = INFINITY;
			row_lb[0] = 10;
			row_ub[0] = 10;
		}
		else if(row ==1)
		{
			assert(*n==1 && *m ==1 );
			col_lb[0] = -INFINITY;
			col_ub[0] = INFINITY;
			row_lb[0] = -INFINITY;
			row_ub[0] = 5;
		}
		else if(row ==2)
		{
			assert(*n==1 && *m ==1 );
			col_lb[0] = -INFINITY;
			col_ub[0] = INFINITY;
			row_lb[0] = -INFINITY;
			row_ub[0] = 6;
		}
		else
			assert(false);
	}

	return 1;
}

int str_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_prob_info -- row " << row <<" col "<<col );
	assert(row == col);
	if(row == 0 )
	{   // (x1+x2)^2
		*obj = (x0[0] + x0[1]) * (x0[0] + x0[1]) ;
	}
	else if(row == 1 || row == 2)
	{   //(x1+x2)*x3
		*obj = (x1[0]+x1[0])*x1[0];
	}
	else
		assert(false);
	return 1;
}

int str_eval_g(double* x0, double* x1, double* eq_g, double* inq_g,
		CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_eval_g  -- row " << row <<" col "<<col);
	assert(row == col);
	if(row == 0)
	{	//x1 + x2 = 100
		eq_g[0] = x0[0] * x0[1];
	}
	else if(row == 1 || row == 2)
	{   //x2^2 + x3|x4 * x1
		inq_g[0] = x0[1]*x0[1] + x1[0]*x0[0];
	}
	else
		assert(false);

	return 1;
}

int str_eval_grad_f(double* x0, double* x1, double* grad, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_eval_grad_f -- row " << row <<" col "<<col );

	if(row == 0 && col == 0)
	{
		grad[0] = 2.0 * (x0[0]+x0[1]);
		grad[1] = 2.0 * (x0[0]+x0[1]);
	}
	else if(row == 1 && col == 1)
	{
		grad[0] = (x0[0]+x0[1]);
	}
	else if(row == 2 && col == 2)
	{
		grad[0] = (x0[0]+x0[1]);
	}
	else if(row == 1 && col == 0)
	{
		grad[0] = x1[0];
		grad[1] = x1[0];
	}
	else if(row == 2 && col == 0)
	{
		grad[0] = x1[0];
		grad[1] = x1[0];
	}
	else
		assert(false);

	return 1;
}

int str_eval_jac_g(double* x0, double* x1, int* e_nz, double* e_elts,
		int* e_rowidx, int* e_colptr, int* i_nz, double* i_elts, int* i_rowidx,
		int* i_colptr, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_eval_jac_g  -- row " << row <<" col "<<col);
	if(e_colptr==NULL && i_colptr == NULL)
	{
		assert(e_elts == NULL && e_rowidx == NULL && e_colptr == NULL);
		assert(i_elts == NULL && i_rowidx == NULL && i_colptr == NULL);
		if (row == 0 && col == 0) {
			*e_nz = 2;
			*i_nz = 0;
		} else if (row == 1 && col == 1) {
			*e_nz = 0;
			*i_nz = 1;
		} else if (row == 2 && col == 2) {
			*e_nz = 0;
			*i_nz = 1;
		} else if (row == 1 && col == 0) {
			*e_nz = 0;
			*i_nz = 2;
		} else if (row == 2 && col == 0) {
			*e_nz = 0;
			*i_nz = 2;
		} else
			assert(false);
	}
	else
	{
		if (row == 0 && col == 0) {
			assert(*i_nz == 0 && *e_nz == 2);
			e_rowidx[0] = 0;
			e_rowidx[1] = 0;
			e_colptr[0] = 0;
			e_colptr[1] = 1;
			e_colptr[2] = 2;
			e_elts[0] = x0[1];
			e_elts[1] = x0[0];
		} else if (row == 1 && col == 1) {
			assert(*i_nz == 1 && *e_nz == 0);
			i_rowidx[0] = 0;
			i_colptr[0] = 0;
			i_colptr[1] = 1;
			i_elts[0] = x0[0];
		} else if (row == 2 && col == 2) {
			assert(*i_nz == 1 && *e_nz == 0);
			i_rowidx[0] = 0;
			i_colptr[0] = 0;
			i_colptr[1] = 1;
			i_elts[0] = x0[0];
		} else if (row == 1 && col == 0) {
			assert(*i_nz == 2 && *e_nz == 0);
			i_rowidx[0] = 0;
			i_rowidx[1] = 0;
			i_colptr[0] = 0;
			i_colptr[1] = 1;
			i_colptr[2] = 2;
			i_elts[0] = x1[0];
			i_elts[1] = 2.0*x0[1];
		} else if (row == 2 && col == 0) {
			assert(*i_nz == 2 && *e_nz == 0);
			i_rowidx[0] = 0;
			i_rowidx[1] = 0;
			i_colptr[0] = 0;
			i_colptr[1] = 1;
			i_colptr[2] = 2;
			i_elts[0] = x1[0];
			i_elts[1] = 2.0*x0[1];
		} else
			assert(false);
	}


	return 1;
}

int str_eval_h(double* x0, double* x1, double* lambda, int* nz, double* elts,
		int* rowidx, int* colptr, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	double obj_factor = 1.0;
	PAR_DEBUG("str_eval_h  -- row " << row <<" col "<<col);
	if(colptr==NULL)
	{
		assert(rowidx == NULL);
		assert(colptr == NULL);
		if (row == 0 && col == 0) {
			*nz = 3;
		} else if (row == 1 && col == 1) {
			*nz = 0;
		} else if (row == 2 && col == 2) {
			*nz = 0;
		} else if (row == 1 && col == 0) {
			*nz = 3;
		} else if (row == 2 && col == 0) {
			*nz = 3;
		} else if (row == 0 && col == 1) {
			*nz = 2;
		} else if (row == 0 && col == 2) {
			*nz = 2;
		} else
			assert(false);
	}
	else{
		if (row == 0 && col == 0) {
			rowidx[0] = 0;
			rowidx[1] = 1;
			rowidx[2] = 1;
			colptr[0] = 0;
			colptr[1] = 2;
			colptr[2] = 3;
			elts[0] = 2.0 * obj_factor;
			elts[1] = 2.0 * obj_factor + lambda[0]*1.0;
			elts[2] = 2.0 * obj_factor;
		} else if (row == 1 && col == 1) {
			assert(*nz == 0);
		} else if (row == 2 && col == 2) {
			assert(*nz == 0);
		} else if (row == 1 && col == 0) { //parent contribution
			rowidx[0] = 0;
			rowidx[1] = 1;
			rowidx[2] = 1;
			colptr[0] = 0;
			colptr[1] = 2;
			colptr[2] = 3;
			elts[0] = 0.0;
			elts[1] = 0.0;
			elts[2] = 2.0 * lambda[0];
		} else if (row == 2 && col == 0) { //parent contribution
			rowidx[0] = 0;
			rowidx[1] = 1;
			rowidx[2] = 1;
			colptr[0] = 0;
			colptr[1] = 2;
			colptr[2] = 3;
			elts[0] = 0.0;
			elts[1] = 0.0;
			elts[2] = 2.0 * lambda[0];
		} else if (row == 0 && col == 1) {
			rowidx[0] = 0;
			rowidx[1] = 0;
			colptr[0] = 0;
			colptr[1] = 1;
			colptr[2] = 2;
			elts[0] = obj_factor*1.0 + lambda[0]*1.0;
			elts[1] = obj_factor*1.0;
		} else if (row == 0 && col == 2) {
			rowidx[0] = 0;
			rowidx[1] = 0;
			colptr[0] = 0;
			colptr[1] = 1;
			colptr[2] = 2;
			elts[0] = obj_factor*1.0 + lambda[0]*1.0;
			elts[1] = obj_factor*1.0;
		} else
			assert(false);
	}

	return 1;
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	PAR_DEBUG("start");
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Comm_rank(comm, &gmyid);
	MPI_Comm_size(comm, &gnprocs);

	str_init_x0_cb init_x0 = &str_init_x0;
	str_prob_info_cb prob_info = &str_prob_info;
	str_eval_f_cb eval_f = &str_eval_f;
	str_eval_g_cb eval_g = &str_eval_g;
	str_eval_grad_f_cb eval_grad_f = &str_eval_grad_f;
	str_eval_jac_g_cb eval_jac_g = &str_eval_jac_g;
	str_eval_h_cb eval_h = &str_eval_h;

	PipsNlpProblemStructPtr prob = CreatePipsNlpProblemStruct(MPI_COMM_WORLD, 2,
			4, 3, init_x0, prob_info, eval_f, eval_g, eval_grad_f, eval_jac_g,
			eval_h, NULL);

	PAR_DEBUG("problem created");

	PipsNlpSolveStruct(prob);

	PAR_DEBUG("end solve ");
	MPI_Barrier(comm);
}
