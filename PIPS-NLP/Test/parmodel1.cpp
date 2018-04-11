#include "./Drivers/parallelPipsNlp_C_Callback.h"

#include "mpi.h"
#include "global_var.h"
#include <iostream>
#include <cassert>
#include <math.h>

#ifdef TIMING
  double timeFromAMPL;
  double probGenTime;
  double PartSolver_GenTime;
  double PartSolver_SolTime;  
  double PartSolver_FactTime;
  int call_sol_Times;
  int call_fact_Times;

  int call_sol_Times_MA57;
  int call_fact_Times_MA57;  
  double genTime_localAmpl;
#endif

//# min (x1+x2)^2+(x1+x2)*x3 + (x1+x2)*x4
//# st.
//#     x1 * x2 = 10
//#     x2^2 + x3*x1 < 5
//#     x2^2 + x4*x1 < 6
//# x1, x2 , x3, x4 free variables

#define S  10    //S scenarios
#define N0 20    // number of first-stage vars
#define N1 2*N0  // number of second-stage vars for scen 1
#define NS 4*N0  // number of second-stage vars for scen > 1

/* min f0(x0) + 1/S sum fs(xs) 
 * s.t.  x0>=0.25, xs>=0.5, x0_1<=100, xs_1<=100
 *       x0_i + x_s_i =1, i=1,N0, s=1,S
 *       sum(x0)+sum(x1)<=N0+N1
 *       sum(xs) <= NS-2, s=1,S
 *
 * here f0(x) = 0.5*||x||^2, fs(x)=0.5*||x-1||^2
 * this problem has a different number of variables and constraints per scenario
 */


int str_init_x0(double* x0, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	MESSAGE("str_init_x0 -- row " << row <<" col "<<col);
	assert(row == col);
	if(row == 0)
	{
	  for(int i=0; i<N0; i++) x0[i]=0.75;
	}
	else if(row == 1)
	{
	  for(int i=0; i<N1; i++) x0[i]=0.75;
	}
	else //if(row >= 2)
	{
	  assert(row<=S);
	  for(int i=0; i<NS; i++) x0[i]=0.75;
	}
	return 1;
}

int str_prob_info(int* n, double* col_lb, double* col_ub, int* m,
		double* row_lb, double* row_ub, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	MESSAGE("str_prob_info -- row " << row <<" col "<<col);
	assert(row == col);
        int type = cbd->typeflag;
        if(type == 1){
	  if(row_lb == NULL){
	    assert(row_ub == NULL);
	    *m = 0;
	  }
	  return 1;
        }
	if(col_lb == NULL)
	{
	  assert(row_lb == NULL);
	  assert(row_ub == NULL);
	  assert(col_ub == NULL);
	  if(row==0) {
	    *n = N0;
	    *m = 0;
	  } else if(row ==1 ) {
	    *n = N1;
	    *m = N0+2;
	  } else {
	    assert(row<=S);
	    *n = NS;
	    *m = N0+1;
	  }
	} else {
	  if(row==0)
	    {
	      assert(*n==N0 && *m == 0);
	      for(int i=0; i<N0; i++) {
		col_lb[i] = 0.25;
		col_ub[i] = INFINITY;
	      }
	      col_ub[0] = 100.;
	    }
	  else if(row ==1) {
	    assert(*n==N1 && *m == N0+2 );
	    for(int i=0; i<N1; i++) {
		col_lb[i] = 0.5;
		col_ub[i] = INFINITY;
	    }
	    col_ub[0] = 100.;

	    for(int i=0; i<N0; i++) {row_lb[i] = 1.;  row_ub[i] = 1.;}
	    row_lb[N0]   = -INFINITY;  row_ub[N0]   = N0+N1;
	    row_lb[N0+1] = -INFINITY;  row_ub[N0+0] = N1-2;


	  }
	  else if(row ==2)  {
	    assert(row<=NS);
	    assert(*n==NS && *m ==NS+1 );
	    for(int i=0; i<NS; i++) { col_lb[i] = 0.5; col_ub[i]= INFINITY; }
	    col_ub[0] = 100.;

	    for(int i=0; i<N0; i++) { row_lb[i] = 1.; row_ub[i] = 1.;} //eq constraints
	    row_lb[N0] = -INFINITY; row_ub[N0] = NS-2;
	  }
	}

	return 1;
}

int str_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	MESSAGE("str_prob_info -- row " << row <<" col "<<col );
	assert(row == col);
	if(row == 0 ) { 
	  *obj = 0.; //1/2*||x||^2
	  for(int i=0; i<N0; i++) *obj += (x0[i]*x0[i]);
	  *obj *= 0.5;
	} else if(row == 1) {
	  double aux; *obj = 0.; //1/2*||x||^2
	   for(int i=0; i<N1; i++) {aux=x1[i]-1; *obj += (aux*aux); }
	   //*obj = (x0[0]+x0[1])*x1[0];
	   *obj *=0.5;
	} else 	{  
	  assert(row<=S);
	  double aux; *obj = 0.; //1/2*||x||^2
	   for(int i=0; i<NS; i++) {aux=x1[i]-1; *obj += (aux*aux); }
	   //*obj = (x0[0]+x0[1])*x1[0];
	   *obj *=0.5;
	}
	return 1;
}

int str_eval_g(double* x0, double* x1, double* eq_g, double* inq_g,
		CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	MESSAGE("str_eval_g  -- row " << row <<" col "<<col);
	assert(row == col);
	if(row == 0) {
	  ////x1 + x2 = 100
	  //	eq_g[0] = x0[0] * x0[1];
	  //no constraints
	} else if(row == 1) {
	  //x2^2 + x3|x4 * x1
	  //inq_g[0] = x0[1]*x0[1] + x1[0]*x0[0];
	  for(int i=0; i<N0; i++) eq_g[i] = x1[i]+x0[i];

	  inq_g[0]=0;
	  for(int i=0; i<N0; i++) inq_g[0] += x0[i];
	  for(int i=0; i<N1; i++) inq_g[0] += x1[i];

	  inq_g[1]=0;
	  for(int i=0; i<N1; i++) inq_g[1] += x1[i];
	} else {
	  assert(row <= S);
	  for(int i=0; i<N0; i++) eq_g[i] = x1[i]+x0[i];
	  inq_g[0]=0;
	  for(int i=0; i<NS; i++) inq_g[0] += x0[i];
	}
	return 1;
}

int str_eval_grad_f(double* x0, double* x1, double* grad, CallBackDataPtr cbd) {
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	MESSAGE("str_eval_grad_f -- row " << row <<" col "<<col );

	if(row == 0 && col == 0)
	{
	  //	grad[0] = 2.0 * (x0[0]+x0[1]);
	  //	grad[1] = 2.0 * (x0[0]+x0[1]);
	  for(int i=0; i<N0; i++) grad[i] = x0[i];
	}
	else if(row == 1 && col == 1)
	{
	  //	grad[0] = (x0[0]+x0[1]);
	  for(int i=0; i<N1; i++)  grad[i] = x1[i];

	} else if(row == col )
	{
	  assert(row<=NS);
	  //	grad[0] = (x0[0]+x0[1]);
	  for(int i=0; i<NS; i++) grad[i] = x1[i];
	}
	else if(row == 1 && col == 0)
	{
	  //row == 1 && col == 0
	  //	grad[0] = x1[0];
	  //	grad[1] = x1[0];
	  for(int i=0; i<N0; i++) grad[i]=0.;
	}
	else if(row >= 2 && col == 0)
	{
	  assert(row<=NS);
	  //row == 2 && col == 0
	  //	grad[0] = x1[0];
	  //	grad[1] = x1[0];
	  for(int i=0; i<N0; i++) grad[i]=0.;
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
	MESSAGE("str_eval_jac_g  -- row " << row <<" col "<<col);
	if(e_colptr==NULL && i_colptr == NULL)
	{
		assert(e_elts == NULL && e_rowidx == NULL && e_colptr == NULL);
		assert(i_elts == NULL && i_rowidx == NULL && i_colptr == NULL);
		if (row == 0 && col == 0) {
		  //	*e_nz = 2;
		  //	*i_nz = 0;
		  *e_nz=0; *i_nz=0;
		} else if (row == 1 && col == 1) {
		  //	*e_nz = 0;
		  //	*i_nz = 1;
		  *e_nz=N0; *i_nz=N1+N1;
		} else if (row == col) {
		  assert(row<=NS);
		  //	*e_nz = 0;
		  //	*i_nz = 1;
		  *e_nz=N0; *i_nz=NS;
		} else if (row == 1 && col == 0) {
		  //	*e_nz = 0;
		  //	*i_nz = 2;
		  *e_nz=N0; *i_nz=N0;
		} else if (row >= 2 && col == 0) {
		  //	*e_nz = 0;
		  //	*i_nz = 2;
		  assert(row<=S);
		  *e_nz=N0; *i_nz=0;
		} else
			assert(false);
	}
	else
	{
		if (row == 0 && col == 0) {
			assert(*i_nz == 0 && *e_nz == 0);

		} else if (row == 1 && col == 1) {
		  assert((*i_nz == N1+N1) && (*e_nz == N0));
		  //equalities
		  for(int i=0; i<N0; i++) {
		    e_rowidx[i] = i;
		    e_colptr[i] = i;
		    e_elts[i]=1.;
		  }
		  e_colptr[N1] = N0;

		  //two inequalities
		  i_colptr[0] = 0; i_colptr[1] = N1;
		  for(int i=0; i<N1; i++) {
		    i_rowidx[i] = 0;
		    i_elts[i] = 1.;
		  }
		  //second inequality 
		  for(int i=N1; i<N1+N1; i++) {
		    i_rowidx[i] = 0;
		    i_elts[i] = 1.;
		  }
		  i_colptr[2] = N1+N1;
		  
		  //	i_rowidx[0] = 0;
		  //	i_colptr[0] = 0;
		  //	i_colptr[1] = 1;
		  //	i_elts[0] = x0[0];
		  
		} else if (row == col ) {
		  assert(*i_nz == NS && *e_nz == N0);
		  //equalities
		  for(int i=0; i<N0; i++) {
		    e_rowidx[i] = i;
		    e_colptr[i] = i;
		    e_elts[i]=1.;
		  }
		  e_colptr[N1] = N0;

		  //one inequality
		  i_colptr[0]=0; i_colptr[1]=NS;
		  for(int i=0; i<NS; i++) {
		    i_rowidx[i]=0;
		    i_elts[i]=1.;
		  }

		  //	i_rowidx[0] = 0;
		  //	i_colptr[0] = 0;
		  //	i_colptr[1] = 1;
		  //	i_elts[0] = x0[0];

		} else if (row == 1 && col == 0) {
		  assert(*i_nz == N0 && *e_nz == N0);
		  //equalities
		  for(int i=0; i<N0; i++) {
		    e_rowidx[i] = i;
		    e_colptr[i] = i;
		    e_elts[i]=1.;
		  }
		  e_colptr[N0]=N0;
		  //the two inequalities
		  i_colptr[0]=0; i_colptr[1]=N0; i_colptr[2]=N0;
		  for(int i=0; i<N0; i++) {
		    i_rowidx[i]=i;
		    i_elts[i]=1.;
		  }
		  //	i_rowidx[0] = 0;
		  //	i_rowidx[1] = 0;
		  //	i_colptr[0] = 0;
		  //	i_colptr[1] = 1;
		  //	i_colptr[2] = 2;
		  //	i_elts[0] = x1[0];
		  //	i_elts[1] = 2.0*x0[1];
		} else if (row >= 2 && col == 0) {
		  assert(*i_nz == 0 && *e_nz == N0);
		  //the equalities
		  //equalities
		  for(int i=0; i<N0; i++) {
		    e_rowidx[i] = i;
		    e_colptr[i] = i;
		    e_elts[i]=1.;
		  }
		  e_colptr[N0]=N0;
		  //the inequality
		  //the matrix is empty
		  //i_colptr[0]=0; i_colptr[1]=0; 		  i_rowidx[0]=0;
		  
		  //	i_rowidx[0] = 0;
		  //	i_rowidx[1] = 0;
		  //	i_colptr[0] = 0;
		  //	i_colptr[1] = 1;
		  //	i_colptr[2] = 2;
		  //	i_elts[0] = x1[0];
		  //	i_elts[1] = 2.0*x0[1];
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
	MESSAGE("str_eval_h  -- row " << row <<" col "<<col);
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

int str_write_solution(double* x, double* lam_eq, double* lam_ieq,CallBackDataPtr cbd)
{
	int row = cbd->row_node_id;
	int col = cbd->col_node_id;
	assert(row == col);
	MESSAGE("write_solution  -- row " << row <<" col "<<col);
	if(row == 0)
	{
		PRINT_ARRAY("node = 0 - x ", x, 2);
		PRINT_ARRAY("node = 0 - eq ", lam_eq, 1);
		PRINT_ARRAY("node = 0 - ieq ", lam_ieq , 0);
	}
	else if(row == 1 || row == 2)
	{
		PRINT_ARRAY("node = "<<row<<" - x ", x, 1);
		PRINT_ARRAY("node = "<<row<<" - eq ", lam_eq, 0);
		PRINT_ARRAY("node = "<<row<<" - ieq ", lam_ieq , 1);
	}
	else
	{
		assert(false);
	}
	return 1;
}

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MESSAGE("start");
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

	PipsNlpProblemStructPtr prob = CreatePipsNlpProblemStruct(MPI_COMM_WORLD, 2,
			init_x0, prob_info, eval_f, eval_g, eval_grad_f, eval_jac_g,
			eval_h, write_solution, NULL);

	MESSAGE("problem created");

	PipsNlpSolveStruct(prob);

	MESSAGE("end solve ");
	MPI_Barrier(comm);
    MPI_Finalize();
}
