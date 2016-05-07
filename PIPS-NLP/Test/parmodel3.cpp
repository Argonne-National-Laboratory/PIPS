#include "./Drivers/parallelPipsNlp_C_Callback.h"

#include "mpi.h"
#include "par_macro.h"
#include <iostream>
#include <cassert>
#include <math.h>


//# min 1.*x[1] + 3*x[3] + 4.*x[4] + 12*y1[2] + 13*y1[3];
//# st.
//s.t. con1: -100<= 0.1*x[1] + 	0.2*x[2] 																<= 100; 
//s.t. con2: -101<=  			1.1*x[2] + 	1.2*x[3] 													<= 101;
//s.t. con3: -102<=  						2.1*x[3] + 	2.2*x[4] 										<= 102;
//s.t. con4: -201<= 10.1*x[1] + 10.2*x[2] 						- 10.3*y1[1] 	+ 10.4*y1[2]				<= 201;
//s.t. con5: -202<= 			11.1*x[2] + 	11.2*x[3] 			- 11.3*y1[1] 	+ 11.4*y1[2] 	+ 11.5*y1[3]	<= 202;
//s.t. con6: -203<= 						12.1*x[3] + 	12.2*x[4]  		     	+ 12.4*y1[2] 	-  12.5*y1[3]	<= 203;
//s.t. con7: 								1*x[3] 	 + 	2*x[4] 										 == 107 ; 
//# x free variables

static double scalObj = 0.1;

int str_init_x0(double* x0, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_init_x0 -- row " << row <<" col "<<col);
	assert(row == col);
	if(row == 0)
	{
		x0[0] = 1.0;
		x0[1] = 1.0;
		x0[2] = 1.0;
		x0[3] = 1.0;		
	}
	else if(row == 1)
	{
		x0[0] = 1.0;
		x0[1] = 1.0;
		x0[2] = 1.0;
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
			*n = 4;
			*m = 4;
		}
		else if(row ==1)
		{
			*n = 3;
			*m = 3;
		}
		else
			assert(false);
	}
	else
	{
		if(row==0)
		{
			assert(*n==4 && *m == 4);
			col_lb[0] = -INFINITY;
			col_lb[1] = -INFINITY;
			col_lb[2] = -INFINITY;
			col_lb[3] = -INFINITY;			
			col_ub[0] = INFINITY;
			col_ub[1] = INFINITY;
			col_ub[2] = INFINITY;
			col_ub[3] = INFINITY;			
			row_lb[1] = -10.0;
			row_ub[1] = 10.0;
			row_lb[2] = -10.1;
			row_ub[2] = 10.1;
			row_lb[3] = -10.2;
			row_ub[3] = 10.2;
			row_lb[0] = 10.7;
			row_ub[0] = 10.7;			
		}
		else if(row ==1)
		{
			assert(*n==3 && *m ==3 );
			col_lb[0] = -INFINITY;
			col_lb[1] = -INFINITY;
			col_lb[2] = -INFINITY;
			col_ub[0] = INFINITY;
			col_ub[1] = INFINITY;
			col_ub[2] = INFINITY;
			row_lb[0] = -20.1;
			row_ub[0] = 20.1;
			row_lb[1] = -20.2;
			row_ub[1] = 20.2;
			row_lb[2] = -20.3;
			row_ub[2] = 20.3;
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
	{   
		*obj = scalObj*(1* x0[0] +  3 * x0[2] + 4*x0[3]);
	}
	else if(row == 1)
	{   
		*obj = scalObj*(12*x1[1] + 13*x1[2]);
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
		inq_g[0] = 0.1*x0[0] + 0.2*x0[1];
		inq_g[1] = 1.1*x0[1] + 1.2*x0[2];
		inq_g[2] = 2.1*x0[2] + 2.2*x0[3];
		eq_g[0] = 1.0*x0[2] + 2.0*x0[3];
	}
	else if(row == 1)
	{   //x2 + x3 + x4
		inq_g[0] = 10.1*x0[0] + 10.2*x0[1] - 10.3*x1[0] + 10.4*x1[1];
		inq_g[1] = 11.1*x0[1] + 11.2*x0[2] - 11.3*x1[0] + 11.4*x1[1] + 11.5*x1[2];
		inq_g[2] = 12.1*x0[2] + 12.2*x0[3] + 12.4*x1[1] - 12.5*x1[2];
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
		grad[0] = scalObj*1.0;
		grad[1] = scalObj*0.0;
		grad[2] = scalObj*3.0;
		grad[3] = scalObj*4.0;
	}
	else if(row == 1 && col == 1)
	{
		grad[0] = scalObj*0.0;
		grad[1] = scalObj*12.0;
		grad[2] = scalObj*13.0;
	}
	else if(row == 1 && col == 0)
	{
		grad[0] = scalObj*0.0;
		grad[1] = scalObj*0.0;
		grad[2] = scalObj*0.0;
		grad[3] = scalObj*0.0;		
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
			*i_nz = 6;
		} else if (row == 1 && col == 1) {
			*e_nz = 0;
			*i_nz = 7;
		} else if (row == 1 && col == 0) {
			*e_nz = 0;
			*i_nz = 6;
		} else
			assert(false);
	}
	else
	{
		if (row == 0 && col == 0) {
			assert(*i_nz == 6 && *e_nz == 2);
			e_rowidx[0] = 0;
			e_rowidx[1] = 0;
			e_colptr[0] = 0;
			e_colptr[1] = 0;
			e_colptr[2] = 0;
			e_colptr[3] = 1;
			e_colptr[4] = 2;			
			e_elts[0] = 1.0;
			e_elts[1] = 2.0;

		i_rowidx[0] = 0;
		i_rowidx[1] = 0;
		i_rowidx[2] = 1;
		i_rowidx[3] = 1;
		i_rowidx[4] = 2;
		i_rowidx[5] = 2;		
		i_colptr[0] = 0;
		i_colptr[1] = 1;
		i_colptr[2] = 3;
		i_colptr[3] = 5;
		i_colptr[4] = 6;
		i_elts[0] = 0.1;
		i_elts[1] = 0.2;
		i_elts[2] = 1.1;
		i_elts[3] = 1.2;
		i_elts[4] = 2.1;
		i_elts[5] = 2.2;

			
		} else if (row == 1 && col == 1) {
			assert(*i_nz == 7 && *e_nz == 0);
			
			i_rowidx[0] = 0;
			i_rowidx[1] = 1;
			i_rowidx[2] = 0;
			i_rowidx[3] = 1;
			i_rowidx[4] = 2;
			i_rowidx[5] = 1;		
			i_rowidx[6] = 2;		
			i_colptr[0] = 0;
			i_colptr[1] = 2;
			i_colptr[2] = 5;
			i_colptr[3] = 7;
			i_elts[0] = -10.3;
			i_elts[1] = -11.3;
			i_elts[2] = 10.4;
			i_elts[3] = 11.4;
			i_elts[4] = 12.4;
			i_elts[5] = 11.5;
			i_elts[6] = -12.5;
		} else if (row == 1 && col == 0) {
			assert(*i_nz == 6 && *e_nz == 0);
			i_rowidx[0] = 0;
			i_rowidx[1] = 0;
			i_rowidx[2] = 1;
			i_rowidx[3] = 1;
			i_rowidx[4] = 2;
			i_rowidx[5] = 2;		
			i_colptr[0] = 0;
			i_colptr[1] = 1;
			i_colptr[2] = 3;
			i_colptr[3] = 5;
			i_colptr[4] = 6;
			i_elts[0] = 10.1;
			i_elts[1] = 10.2;
			i_elts[2] = 11.1;
			i_elts[3] = 11.2;
			i_elts[4] = 12.1;
			i_elts[5] = 12.2;

		} else
			assert(false);
	}


	return 1;
}

int str_eval_h(double* x0, double* x1, double* lambda, int* nz, double* elts,
		int* rowidx, int* colptr, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	PAR_DEBUG("str_eval_h  -- row " << row <<" col "<<col);
	if(colptr==NULL)
	{
		assert(rowidx == NULL);
		assert(colptr == NULL);
		if (row == 0 && col == 0) {
			*nz = 0;
		} else if (row == 1 && col == 1) {
			*nz = 0;
		} else if (row == 1 && col == 0) {
			*nz = 0;
		} else if (row == 0 && col == 1) {
			*nz = 0;
		} else
			assert(false);
	}
	else{
		if (row == 0 && col == 0) {
			assert(*nz == 0);
		} else if (row == 1 && col == 1) {
			assert(*nz == 0);
		} else if (row == 1 && col == 0) { //parent contribution
			assert(*nz == 0);
		} else if (row == 0 && col == 1) {
			assert(*nz == 0);
//			assert(*nz == 2);
//			rowidx[0] = 0;
//			rowidx[1] = 1;
//			colptr[0] = 0;
//			colptr[1] = 1;
//			colptr[2] = 2;
//			elts[0] = 1.0;
//			elts[1] = 1.0;
		} else if (row == 0 && col == 2) {
			assert(*nz == 0);
//			assert(*nz == 2);
//			rowidx[0] = 0;
//			rowidx[1] = 1;
//			colptr[0] = 0;
//			colptr[1] = 1;
//			colptr[2] = 2;
//			elts[0] = 1.0;
//			elts[1] = 1.0;
		} else
			assert(false);
	}

	return 1;
}


int str_write_solution(double* x, double* lam_eq, double* lam_ieq, CallBackDataPtr cbd) {

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
	str_write_solution_cb write_solution = &str_write_solution;

	PipsNlpProblemStructPtr prob = CreatePipsNlpProblemStruct(MPI_COMM_WORLD, 1,
			init_x0, prob_info, eval_f, eval_g, eval_grad_f, eval_jac_g,
			eval_h, write_solution, NULL);


	PAR_DEBUG("problem created");

	PipsNlpSolveStruct(prob);

	PAR_DEBUG("end solve ");
	MPI_Barrier(comm);
}
