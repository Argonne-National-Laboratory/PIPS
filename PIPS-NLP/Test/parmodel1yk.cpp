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

// min x5
// st.  x3*x4/x2 - x5 =0
//      x1 -  x3 = 0
//     x1, first stage ,:=10
//     x2, [1, 100], :=50 
//     x3, [0, 100], :=10
//     x4, [1, 500], :=100
//     x5, [0, 3  ], :=1

int str_init_x0(double* x0, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  MESSAGE("str_init_x0 -- row " << row <<" col "<<col);
  assert(row == col);
  if (row == 0) {
    x0[0] = 10;
  }
  else if (row == 1) {
    x0[0] = 50;
    x0[1] = 10;
    x0[2] = 100;
    x0[3] = 1;
  }
  else
    assert(false);

  return 1;
}

int str_prob_info(int* n, double* col_lb, double* col_ub, int* m, double* row_lb, double* row_ub, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  MESSAGE("str_prob_info -- row " << row <<" col "<<col);
  assert(row == col);
  int type = cbd->typeflag;
  if (type == 1) {
    if (row_lb == NULL) {
      assert(row_ub == NULL);
      *m = 0;
    }
    return 1;
  }

  if (col_lb == NULL) {
    assert(row_lb == NULL);
    assert(row_ub == NULL);
    assert(col_ub == NULL);
    if (row == 0) {
      *n = 1;
      *m = 0;
    }
    else if (row == 1) {
      *n = 4;
      *m = 2;
    }
    else
      assert(false);
  }
  else {
    if (row == 0) {
      assert(*n == 1 && *m == 0);
      col_lb[0] = -INFINITY;
      col_ub[0] = INFINITY;
    }
    else if (row == 1) {
      assert(*n == 4 && *m == 2);
      col_lb[0] = 1;
      col_ub[0] = 100;
      col_lb[1] = 0;
      col_ub[1] = 100;
      col_lb[2] = 1;
      col_ub[2] = 500;
      col_lb[3] = 0;
      col_ub[3] = 3;

      row_lb[0] = 0;
      row_ub[0] = 0;
      row_lb[1] = 0;
      row_ub[1] = 0;
    }
    else
      assert(false);
  }

  return 1;
}

int str_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  MESSAGE("str_prob_info -- row " << row <<" col "<<col);
  assert(row == col);
  if (row == 0) {
    *obj = 0;
  }
  else if (row == 1) {
    *obj = x1[3];
  }
  else
    assert(false);
  MESSAGE(" f: "<<*obj);
  return 1;
}

int str_eval_g(double* x0, double* x1, double* eq_g, double* inq_g, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  MESSAGE("str_eval_g  -- row " << row <<" col "<<col);
  assert(row == col);
  if (row == 0) {
  }
  else if (row == 1) {
    eq_g[0] = -x1[2] * x1[1] / x1[0] + x1[3];
    eq_g[1] = -x0[0] + x1[1];
    PRINT_ARRAY("node = 1  - eq_g ", eq_g, 2);
  }
  else
    assert(false);

  return 1;
}

int str_eval_grad_f(double* x0, double* x1, double* grad, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  MESSAGE("str_eval_grad_f -- row " << row <<" col "<<col);
  if (row == 0 && col == 0) {
    grad[0] = 0;
    PRINT_ARRAY("grad =  ", grad, 1);
  }
  else if (row == 1 && col == 1) {
    grad[0] = 0;
    grad[1] = 0;
    grad[2] = 0;
    grad[3] = 1;
    PRINT_ARRAY("grad =  ", grad, 4);
  }
  else if (row == 1 && col == 0) {
    PRINT_ARRAY(" - x1 ", x1, 4);
    grad[0] = 0;
    PRINT_ARRAY("grad =  ", grad, 1);
  }
  else
    assert(false);

  return 1;
}

int str_eval_jac_g(double* x0, double* x1, int* e_nz, double* e_elts, int* e_rowidx, int* e_colptr, int* i_nz, double* i_elts, int* i_rowidx, int* i_colptr, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  MESSAGE("str_eval_jac_g  -- row " << row <<" col "<<col);
  if (e_colptr == NULL && i_colptr == NULL) {
    assert(e_elts == NULL && e_rowidx == NULL && e_colptr == NULL);
    assert(i_elts == NULL && i_rowidx == NULL && i_colptr == NULL);
    if (row == 0 && col == 0) {
      PRINT_ARRAY("x0 - ", x0, 1);
      *e_nz = 0;
      *i_nz = 0;
    }
    else if (row == 1 && col == 1) {
      PRINT_ARRAY("x0 - ", x0, 1);
      PRINT_ARRAY("x1 - ", x1, 4);
      *e_nz = 5;
      *i_nz = 0;
    }
    else if (row == 1 && col == 0) {
      PRINT_ARRAY("x0 - ", x0, 1);
      PRINT_ARRAY("x1 - ", x1, 4);
      *e_nz = 1;
      *i_nz = 0;
    }
    else
      assert(false);
  }
  else {
    if (row == 0 && col == 0) {
      PRINT_ARRAY("x0 - ", x0, 1);
    }
    else if (row == 1 && col == 1) {
      PRINT_ARRAY("x0 - ", x0, 1);
      PRINT_ARRAY("x1 - ", x1, 4);
      e_rowidx[0] = 0;
      e_rowidx[1] = 0;
      e_rowidx[2] = 1;
      e_rowidx[3] = 0;
      e_rowidx[4] = 0;

      e_colptr[0] = 0;
      e_colptr[1] = 1;
      e_colptr[2] = 3;
      e_colptr[3] = 4;
      e_colptr[4] = 5;

      e_elts[0] = x1[2] * x1[1] / x1[0] / x1[0];
      e_elts[1] = -x1[2] / x1[0];
      e_elts[2] = 1;
      e_elts[3] = -x1[1] / x1[0];
      e_elts[4] = 1;

      //PRINT_ARRAY("node = 1 1 - x ", x1, 4);
      //PRINT_ARRAY("node = 1 1 - rowidx ", e_rowidx , 5);
      //PRINT_ARRAY("node = 1 1 - colptr ", e_colptr , 5);
      //PRINT_ARRAY("node = 1 1 - elts ", e_elts , 5);

    }
    else if (row == 1 && col == 0) {
      PRINT_ARRAY("x0 - ", x0, 1);
      PRINT_ARRAY("x1 - ", x1, 4);
      e_rowidx[0] = 1;
      e_colptr[0] = 0;
      e_colptr[1] = 1;
      e_elts[0] = -1;

      //PRINT_ARRAY("node = 1 0 - rowidx ", e_rowidx , 1);
      //PRINT_ARRAY("node = 1 0 - colptr ", e_colptr , 2);
      //PRINT_ARRAY("node = 1 0 - elts ", e_elts , 1);
    }
    else
      assert(false);
  }

  return 1;
}

int str_eval_h(double* x0, double* x1, double* lambda, int* nz, double* elts, int* rowidx, int* colptr, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  double obj_factor = 1.0;
  MESSAGE("str_eval_h  -- row " << row <<" col "<<col);
  if (colptr == NULL) {
    assert(rowidx == NULL);
    assert(colptr == NULL);
    if (row == 0 && col == 0) {
      *nz = 0;
    }
    else if (row == 1 && col == 1) {
      *nz = 6;
    }
    else if (row == 1 && col == 0) {
      *nz = 0;
    }
    else if (row == 0 && col == 1) {
      *nz = 0;
    }
    else
      assert(false);
  }
  else {
    if (row == 0 && col == 0) {
    }
    else if (row == 1 && col == 1) {
      assert(*nz == 6);

      rowidx[0] = 0;
      rowidx[1] = 0;
      rowidx[2] = 1;
      rowidx[3] = 0;
      rowidx[4] = 1;
      rowidx[5] = 2;

      colptr[0] = 0;
      colptr[1] = 1;
      colptr[2] = 3;
      colptr[3] = 6;
      colptr[4] = 6;

      elts[0] = -2 * lambda[0] * x1[1] * x1[2] / x1[0] / x1[0] / x1[0];
      elts[1] = lambda[0] * x1[2] / x1[0] / x1[0];
      elts[3] = lambda[0] * x1[1] / x1[0] / x1[0];
      elts[2] = 0;
      elts[4] = -lambda[0] / x1[0];
      elts[5] = 0;
      rowidx[0] = 0;

      //PRINT_ARRAY("node = 1 - x ", x1, 4);
      //PRINT_ARRAY("node = 1 - lam ", lambda, 2);
      //PRINT_ARRAY("node = 1 - rowidx ", rowidx , 6);
      //PRINT_ARRAY("node = 1 - colptr ", colptr , 5);
      //PRINT_ARRAY("node = 1 - elts ", elts , 6);
    }
    else if (row == 1 && col == 0) { //parent contribution
    }
    else if (row == 0 && col == 1) {
    }
    else
      assert(false);
  }

  return 1;
}

int str_write_solution(double* x, double* lam_eq, double* lam_ieq, CallBackDataPtr cbd) {
  int row = cbd->row_node_id;
  int col = cbd->col_node_id;
  assert(row == col);
  MESSAGE("write_solution  -- row " << row <<" col "<<col);
  if (row == 0) {
    PRINT_ARRAY("node = 0 - x ", x, 1);
  }
  else if (row == 1 || row == 2) {
    PRINT_ARRAY("node = "<<row<<" - x ", x, 4);
  }
  else {
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

  PipsNlpProblemStructPtr prob = CreatePipsNlpProblemStruct(MPI_COMM_WORLD, 1, init_x0, prob_info, eval_f, eval_g, eval_grad_f, eval_jac_g, eval_h, write_solution, NULL);

  MESSAGE("problem created");

  PipsNlpSolveStruct(prob);

  MESSAGE("end solve ");
  MPI_Barrier(comm);
  MPI_Finalize();
}
