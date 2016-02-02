/* PIPS-NLP                                                         	*
 * Authors: Feng Qiang                      		*
 * (C) 2016 Argonne National Laboratory			*/

#ifndef _PARALLELPIPSNLP_C_CALLBACK_H__
#define _PARALLELPIPSNLP_C_CALLBACK_H__


#include "mpi.h"

typedef void * UserDataPtr;

/*
 * prob points to the userdata field of the PipsNlpProblemStruct.
 * row_node_id and col_node_id are the node ids in the stochastic tree
 *
 * Note:
 * When the requested data are vectors, the row and col node should be equal and are used to identify the sub-vectors to be evaluated, where
 * 	eval_f_grad is an exception. See below for more detail.
 * When the requested data are sub-matrices, the row and col node ids are used to identify the block.
 */
typedef struct CallBackData
{
	UserDataPtr prob;
	int row_node_id;
	int col_node_id;
} CallBackData;

typedef CallBackData * CallBackDataPtr;

/*
 * This function is used to get the initial variable values to PIPS-NLP at the node Id specified in the CallBackData.
 * x0 should be already allocated for the size of declared variable in this node.
 * cbd is the call back data pointer.
 *
 * note: row and col node ids from the cbd should always be equal.
 */
extern "C" typedef int (*str_init_x0_cb)(double* x0, CallBackDataPtr cbd);

/*
 * This function is used to request the values and sizes of the col and row's lower and upper bound vectors.
 * The sizes of col and row are represented by n and m respectively. n and m should be written when setting the lower and upper
 * bounds pointers to NULL.
 * The vectors of the corresponding node id is given in the CallBackData.
 *
 *  note: when providing the row_lb and row_ub vectors, it assumes that the equality constraints should be at the top of these vectors.
 *  	Therefore, it is wise to order the constraint in a way that always keep the equality constraints in front of the
 *  	inequality constraints to keep everything else consistent.
 */
extern "C" typedef int (*str_prob_info_cb)(int* n, double* col_lb, double* col_up, int* m, double* row_lb, double* row_up, CallBackDataPtr cbd);

/*
 * This function evaluates the value of the objective function declared at the node whose Id is specified in the CallBackData.
 * x0 and x1 corresponds the values for first stage and second stage variables.
 * If the objective function on the root node (i.e. the only node in first stage) is requested, x0 and x1 are the same and both point to
 * the first stage variable values.
 */
extern "C" typedef int (*str_eval_f_cb)(double* x0, double* x1, double *obj, CallBackDataPtr cbd);

/*
 * This function computes the constraint value at current x0 and x1, where x0 and x1 are the first and second stage variable values respectively.
 * The results should be written in eq_g (for equality constraints) and inq_g (for inequality constraints).
 * If the constraints to be evaluated are in the root node, x0 and x1 are the same and both point to the first stage variable values.
 *
 */
extern "C" typedef int (*str_eval_g_cb)(double* x0, double* x1, double* eq_g,double* inq_g, CallBackDataPtr cbd);

/*
 * This function evaluates the objective gradient sub-vector. The objective function to be considered is declared in node with node id equals to row_node_id
 * from the CallBackData. The gradient sub-vector is the one with respect to the variables declared in node id equals to col_node_id from the CallBackData.
 *
 * Note:
 * 	evaluating 1st stage objective gradient, x0 and x1 are the same and both point to the first stage variable values.
 * 	evaluating 2st stage objective gradient, x0 and x1 represents the variables values of first and second stage variable values respectively.
 * 	the objective declared in the second stage node could contains a cross product term of first and second stage variables. This term can contribute
 * 		two gradient sub-vectors, where one is with respect to the first stage variables and the other is with respect to the second stage variables.
 * 		For example, let us look at the objective function declared in the node 1. By setting (1,0) in the CallBackData for row_node_id and col_node_id,
 * 		this function should compute the gradient sub-vector with respect to the root stage variables (where node id is 0).
 * 		The gradient sub-vector with respect to the second stage variables should be computed when the row and col node id pair (in the CallBackData)
 * 		equals to (1,1).
 *
 */
extern "C" typedef int (*str_eval_grad_f_cb)(double* x0, double* x1, double* vec_grad_f, CallBackDataPtr cbd);

/*
 * This function evaluates the structure or values of the sub-matrix of the Jacobian in two separate CCS constructs. One is for the
 * equality constraints and the other is for the inequality constraints.
 * The Jacobian sub-block is identified by the row_node_id and col_node_id given in the CallBackData.
 * When setting e_elts and i_elts to NULL, the function should write back the number of non-zeros to e_nz (for equality constraints)
 * and i_nz (for inequality constraints).
 * When e_elts or i_elts is not NULL, e_nz and i_nz should have the correct values set already by the solver. The function should
 * now computes the two Jacobian sub-blocks and write them to the values points triples whose memory is already allocated by the solver.
 *
 * (e_elts, e_rowidx, e_colptr) represents the Jacobian block of the equality constraints
 * (i_elts, i_rowidx, i_colptr) represents the Jacobian block of the inequality constraints
 */
extern "C" typedef int (*str_eval_jac_g_cb)(double* x0, double* x1,
		int* e_nz, double* e_elts, int* e_rowidx, int *e_colptr,
		int* i_nz, double* i_elts, int* i_rowidx, int *i_colptr,
		CallBackDataPtr cbd);

/*
 * This function evaluates the structure or values of the Hessian of the Lagrangian matrix in CCS format.
 * Since PIPS-NLP is not support for linking constraint, the root stage node only produces the diagonal sub-block for the Hessian
 * of the Lagrangian matrix.
 * For example, by setting (0,0) pair for row and col node id in the CallBackData, this function should compute above mentioned diagonal block.
 *
 * Each second stage constraints can produce three sub-blocks for the Hessian of the Lagrangian matrix, i.e. the diagonal block,
 * the border block and the block that contributes to the root diagonal block.
 * For example, let us consider the Hessian of the Lagrangian matrix in following structure (for a two stage n scenario stochastic programming problem),
 * where H_00 is called the root diagonal block and H_0i is called the border blocks for i \in {1...n}.
 *
 * H_00
 * H_01		H_11
 * H_02				H_22
 * H_03						H33
 * 	.							.
 * 	.								.
 * H_0n									H_nn
 *
 * 		The diagonal block: this can be requested by setting (1,1), (2,2) and etc. in the CallBackData for row and col node id.
 * 		The border block  : this can be requested by setting (0,1), (0,2) and etc. in the CallBackData for row and col node id.
 * 		The root diagonal contribution: this can be requested by setting (1,0), (2,0) and etc. in the CallBackData for row and col node id.
 * 		Then, PIPS-NLP will use MPI_Allreduce to compute the final H_00 blocks.
 *
 * Of course, the nz values of the above mentioned blocks should be computed by setting elts to NULL pointer.
 */
extern "C" typedef int (*str_eval_h_cb)(double* x0, double* x1, double* lamdba,
		int* nz, double* elts, int* rowidx, int *colptr,
		CallBackDataPtr cbd);

extern "C"
{

	struct PipsNlpProblemStruct
	{
	  MPI_Comm comm;
	  int nnodes;
	  int n;
	  int m;
	  str_init_x0_cb init_x0;
	  str_prob_info_cb prob_info;
	  str_eval_f_cb eval_f;
	  str_eval_g_cb eval_g;
	  str_eval_grad_f_cb eval_grad_f;
	  str_eval_jac_g_cb eval_jac_g;
	  str_eval_h_cb eval_h;
	  UserDataPtr userdata;
	};

	typedef struct PipsNlpProblemStruct* PipsNlpProblemStructPtr; 	/** Pointer to a pips_nlp Problem. **/

	PipsNlpProblemStructPtr CreatePipsNlpProblemStruct(
		  MPI_Comm comm
		, int nnodes /** number of nodes **/
		, int n /** Number of variables **/
		, int m /** Number of constraints. **/
		, str_init_x0_cb init_x0
		, str_prob_info_cb prob_info
		, str_eval_f_cb eval_f /** Callback function of objective function **/
		, str_eval_g_cb eval_g /** Callback function of constraint body **/
		, str_eval_grad_f_cb eval_grad_f /** Callback function of objective gradient **/
		, str_eval_jac_g_cb eval_jac_g /** Callback function of constraint Jacobian **/
		, str_eval_h_cb eval_h /** Callback function of Lagrangian Hessian **/
		, UserDataPtr userdata
	);

	void FreePipsNlpProblemStruct(PipsNlpProblemStruct* prob);

	int PipsNlpSolveStruct(PipsNlpProblemStruct* prob);
};
#endif

