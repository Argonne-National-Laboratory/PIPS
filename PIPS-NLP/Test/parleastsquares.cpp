#include "./Drivers/parallelPipsNlp_C_Callback.h"

#include "mpi.h"
#include "global_var.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>

// min .5 * (x - x^*)^{T} (x - x^*)
// s.t. x \in [0, 1]^{nx}, given x^* \in [0,1]^{nx},
// where x, x^* are discrete representations of 1D functions on a grid of nx
// equispaced points, indexed by k, where 0 <= k < nx.
//
// Parallel decomposition of unknowns is as follows: x[k] (global) is
// stored in node/row (k % nblocks), nblocks = nscen + 1 x[k] is
// stored in index l such that k = l * nblocks + (k % nblocks) that
// is, l = k / nblocks. This parallel decomposition is a cyclic
// decomposition over the nodes/rows, ensuring automatic load
// balancing, because the difference between the maximum number of
// unknowns in a row and the minimum number of unknowns in a row is at
// most 1.

// NOTE: no error checking done for inline functions; asserts could be added
// Utility functions used for parallel data decomposition:
inline int minNumBlockVars(int nx, int nblocks) { return nx / nblocks; }
inline int lastRowWithExtraVar(int nx, int nblocks) { return nx % nblocks; }
inline double sqr(double x) { return x * x; }
inline double hat_function(double x, double xcen, double r)
{ return ((fabs(x - xcen) < r) ? 1 : 0); }
inline double get_pt(int idx, double lb, double ub, int n)
{ return ((idx * (ub - lb) / (n - 1)) + lb); }
inline int get_k(int l, int row, int nblocks) { return l * nblocks + row; }
inline double x_g(int l, int row, int nblocks, int nx, double xlb, double xub)
{
  const double xcen = 0.0; // center pt of hat function
  const double r = 0.6; // radius of circle
  const int k = get_k(l, row, nblocks);
  const double x = get_pt(k, xlb, xub, nx);
  const double f = hat_function(x, xcen, r);
  MESSAGE("Set x-coordinate for index " << k << " to " << x << " and f(x) = " << f);
  return f;
}
// End utility functions

struct ProbData // Problem data
{
  int nblocks; // # blocks = 1 (for 1st stage) + nscen (2nd stage blocks)
  int v; // minimum number of variables per row/scenario
  int lre; // last row with extra variable; if row < lre, # vars is v + 1
  double *x_given; // local scenario values to be fit
  int nx; // # mesh points in x direction = total # mesh points
  double xlb, xub; // x bounds
};

typedef struct ProbData ProbData;

void setMesh(ProbData* p, int nscen, int nx)
{ p->nblocks = nscen + 1; p->nx = nx; p->xlb = -1; p->xub = 1;}

void setProbData(ProbData* p, int row)
{
  assert(p->nblocks <= p->nx);
  p->lre = lastRowWithExtraVar(p->nx, p->nblocks);
  MESSAGE("Rows with index < " << p->lre << " have an extra variable");
  p->v = minNumBlockVars(p->nx, p->nblocks);
  MESSAGE("Each row must have at least " << p->v << " variables");
  const int v = (row < p->lre) ? (p->v + 1) : p->v;
  MESSAGE("Row " << row << " has " << p->v << " variables");
  p->x_given = new double[v];
  for (int i = 0; i < v; i++)
    { (p->x_given)[i] = x_g(i, row, p->nblocks, p->nx, p->xlb, p->xub); }
}

void delProbData(ProbData* const p) { delete[] p->x_given; }

int str_init_x0(double* x0, CallBackDataPtr cbd) {
  const int row = cbd->row_node_id;
  const int col = cbd->col_node_id;
  assert(row == col);
  MESSAGE("str_init_x0 -- row " << row <<" col "<<col);
  const ProbData* const p = static_cast<ProbData*>(cbd->prob);
  const int v = (row < p->lre) ? (p->v + 1) : p->v;

  // initialize initial guess to random value
  const double f = 1.0 / (static_cast<double>(RAND_MAX));
  for (int i=0; i < v; i++) { x0[i] = f * static_cast<double>(rand()); }
  return 1;
}

int str_prob_info(int* n, double* col_lb, double* col_ub, int* m,
		  double* row_lb, double* row_ub, CallBackDataPtr cbd) {
  const int row = cbd->row_node_id;
  const int col = cbd->col_node_id;
  MESSAGE("str_prob_info -- row " << row <<" col "<<col);
  assert(row == col);

  setProbData(static_cast<ProbData*>(cbd->prob), row);
  const ProbData* const p = static_cast<ProbData*>(cbd->prob);
  const int v = (row < p->lre) ? (p->v + 1) : p->v;
  MESSAGE("str_prob_info -- row " << row << " col " << col <<
	  " v = " << v);

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
    assert(row_lb == NULL); assert(row_ub == NULL); assert(col_ub == NULL);
    *n = v; *m = 0; // every scenario has only bound constraints
  }
  else
  {
    for (int i = 0; i < v; i++)
    { col_lb[i] = 0; col_ub[i] = 1; assert((0 == *m) && (v == *n)); }
  }

  return 1;
}

int str_eval_f(double* x0, double* x1, double* obj, CallBackDataPtr cbd) {
  const int row = cbd->row_node_id;
  const	int col = cbd->col_node_id;
  MESSAGE("str_prob_info -- row " << row <<" col "<<col );

  const ProbData* const p = static_cast<ProbData*>(cbd->prob);
  const int v = (row < p->lre) ? (p->v + 1) : p->v;
  const double * const x_given = static_cast<ProbData*>(cbd->prob)->x_given;

  assert(row == col);
  *obj = 0;

  if (row == 0) { for (int i=0;i<v;i++) { *obj += sqr(x0[i] - x_given[i]); } }
  else { for (int i = 0; i < v; i++) { *obj += sqr(x1[i] - x_given[i]); } }
  *obj *= 0.5; // Scale objective at end to minimize # of multiplications

  return 1;
}

int str_eval_g(double* x0, double* x1, double* eq_g, double* inq_g,
		CallBackDataPtr cbd) {
  const int row = cbd->row_node_id;
  const int col = cbd->col_node_id;
  MESSAGE("str_eval_g -- row " << row << " col " << col);
  assert(row == col);

  return 1;
}

int str_eval_grad_f(double* x0, double* x1, double* grad, CallBackDataPtr cbd) {
  const int row = cbd->row_node_id;
  const	int col = cbd->col_node_id;
  MESSAGE("str_eval_grad_f -- row " << row <<" col "<<col );

  const ProbData* const p = static_cast<ProbData*>(cbd->prob);
  const int v = (row < p->lre) ? (p->v + 1) : p->v;
  const double * const x_given = static_cast<ProbData*>(cbd->prob)->x_given;

  if(row==0) for(int i=0;i<v;i++) grad[i]=((row==col) ? (x0[i]-x_given[i]) : 0);
  else for(int i=0;i<v;i++) grad[i]=((row == col) ? (x1[i]-x_given[i]) : 0);
  
  return 1;
}

int str_eval_jac_g(double* x0, double* x1, int* e_nz, double* e_elts,
		int* e_rowidx, int* e_colptr, int* i_nz, double* i_elts, int* i_rowidx,
		int* i_colptr, CallBackDataPtr cbd) {
  // This problem has only bound constraints
  const int row = cbd->row_node_id;
  const int col = cbd->col_node_id;
  MESSAGE("str_eval_jac_g -- row " << row << " col " << col);
  if ((NULL == e_colptr) && (NULL == i_colptr))
  {
    assert((NULL == e_elts) && (NULL == e_rowidx) && (NULL == e_colptr));
    assert((NULL == i_elts) && (NULL == i_rowidx) && (NULL == i_colptr));
    *e_nz = 0; *i_nz = 0;
  }
  else
  {
    assert(*i_nz == 0 && *e_nz == 0);
  }

  return 1;
}

int str_eval_h(double* x0, double* x1, double* lambda, int* nz, double* elts,
		int* rowidx, int* colptr, CallBackDataPtr cbd) {
  const int row = cbd->row_node_id;
  const	int col = cbd->col_node_id;
  MESSAGE("str_eval_h  -- row " << row <<" col "<<col);

  const ProbData* const p = static_cast<ProbData*>(cbd->prob);
  const int v = (row < p->lre) ? (p->v + 1) : p->v;

  MESSAGE("str_eval_h in row " << row << " col " << col <<
	  " has " << ((row == col) ? v : 0) << " nonzero entries");

  if(colptr==NULL)
    {
      MESSAGE("assigning number of entries in (" << row << ", " << col <<
	      ") block of Hessian...");
      assert(rowidx == NULL);
      assert(colptr == NULL);
      *nz = ((row == col) || (0 == col)) ? v : 0;
      MESSAGE("number of entries in (" << row << ", " << col <<
	      ") block of Hessian assigned!");
    }
  else
 { // Diagonal block is identity; col == 0 blocks are zero element
   // "child blocks" of first stage variables and must be explicitly
   // allocated and assigned values (even if those values are zero!)
    if ((row == col) || (0 == col))
    {
      assert(*nz > 0);
      for (int i = 0; i < v; i++)
      {
	rowidx[i] = i; colptr[i] = i; elts[i] = (row == col) ? 1 : 0;
	MESSAGE("Assigned rowidx[" << i << "] = " << i << "; " <<
		"colptr[" << i << "] = " << i << "; " <<
		"elts[" << i << "] = " << i << " in " <<
		"(" << row << ", " << col << ") block of Hessian");
      }
      colptr[v] = v;
    }
    else if ((row != col) && (col > 0)) { assert(*nz == 0); }
    else assert(false);
  }

  return 1;
}

int str_write_solution(double* x, double* lam_eq, double* lam_ieq,CallBackDataPtr cbd)
{
  const int row = cbd->row_node_id;
  const int col = cbd->col_node_id;

  const ProbData* const p = static_cast<ProbData*>(cbd->prob);
  const int v = (row < p->lre) ? (p->v + 1) : p->v;
  const double * const x_given = static_cast<ProbData*>(cbd->prob)->x_given;
  const int nblocks = static_cast<ProbData*>(cbd->prob)->nblocks;

  assert(row == col);
  MESSAGE("write_solution  -- row " << row <<" col "<<col);

  for (int i = 0; i < v; i++)
  {
    assert(row < nblocks);
    const int k = get_k(i, row, nblocks);
    const double x_coord = get_pt(k, p->xlb, p->xub, p->nx);
    MESSAGE("node = " << row << ", x-coordinate = " << x_coord <<", x[" << k
	    << "] = " << x[i] << ", x_given[" << k << "] = " << x_given[i]
	    << ", lb = " << p->xlb << ", ub = " << p->xub << ", nx = " << p->nx);
  }
  return 1;
}

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  MESSAGE("start");
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &gmyid);
  MPI_Comm_size(comm, &gnprocs);

  MESSAGE("This process is " << gmyid << " of " << gnprocs);
  if (argc == 3)
  {
    const int nscen = atoi(argv[1]); const int nx = atoi(argv[2]);
    assert(nx > 1);
    ProbData* pd = new ProbData;
    setMesh(pd, nscen, nx);
    MESSAGE("# of scenarios is " << nscen << ", so # blocks is " << pd->nblocks);
    
    str_init_x0_cb init_x0 = &str_init_x0;
    str_prob_info_cb prob_info = &str_prob_info;
    str_eval_f_cb eval_f = &str_eval_f;
    str_eval_g_cb eval_g = &str_eval_g;
    str_eval_grad_f_cb eval_grad_f = &str_eval_grad_f;
    str_eval_jac_g_cb eval_jac_g = &str_eval_jac_g;
    str_eval_h_cb eval_h = &str_eval_h;
    str_write_solution_cb write_solution = &str_write_solution;
    
    PipsNlpProblemStructPtr prob = CreatePipsNlpProblemStruct(MPI_COMM_WORLD,
      nscen, init_x0, prob_info, eval_f, eval_g, eval_grad_f, eval_jac_g,
      eval_h, write_solution, static_cast<void*>(pd));

    MESSAGE("problem created");

    PipsNlpSolveStruct(prob);

    MESSAGE("end solve ");

    delProbData(pd);
    delete pd;
  }
  else
  {
    if(0 == gmyid)
    {
      std::cout << "parleastsquares takes 2 arguments" << std::endl;
      std::cout << "parleastsquares [# of scenarios] [# pts in x dir]" << std::endl;
      std::cout << "Example: parleastsquares 4 10" << std::endl;
    }
  }
  MPI_Barrier(comm);
  MPI_Finalize();
}
