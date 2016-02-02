/* PIPS-NLP                                                         	*
 * Authors: Feng Qiang                      		*
 * (C) 2016 Argonne National Laboratory			*/

#ifndef _PARALLELPIPSNLP_C_CALLBACK_H__
#define _PARALLELPIPSNLP_C_CALLBACK_H__


#include "mpi.h"

typedef void * UserDataPtr;

/*
 * prob points to the userdata field of the PipsNlpProblemStruct
 * row_node_id and col_node_id are of of the node id in the stochastic tree
 * When the requested data is vectors, the row and col node must be equal.
 * When the requested data is a sub-matrix, the row and col node ids are used to identify the block.
 */
typedef struct CallBackData
{
	UserDataPtr prob;
	int row_node_id;
	int col_node_id;
} CallBackData;

typedef CallBackData * CallBackDataPtr;

/*
 * getting initial variable values to PIPS-NLP at the node Id specified in the CallBackData.
 * x0 points the memory for variables. x0 should be allocated for the size of declared variable in this node.
 * cbd is the call back data pointer.
 * 	note: row and col node is in the cbd are always equals.
 */
extern "C" typedef int (*str_init_x0_cb)(double* x0, CallBackDataPtr cbd);

/*
 * requesting the col and row's lower and upper bound vectors.
 * The size of col and row represented by n and m are requested by setting these double pointers to NULL
 * The vectors of the corresponding node id is given in the CallBackData.
 *  note: when providing the row_lb and row_ub, it assumes that the equality constraints must be at the top of these vectors.
 *  	Therefore, it is wise to order the constraint in a way that always keep the equality constraints in front of the
 *  	inequality constraints.
 */
extern "C" typedef int (*str_prob_info_cb)(int* n, double* col_lb, double* col_up, int* m, double* row_lb, double* row_up, CallBackData* cbd);

/*
 * evaluate the objective function values for the objective function at the node whose Id is specified in the CallBackData.
 * x0 and x1 corresponds the values for first stage and second stage variables.
 * If the objective function on the root stage node (the only first stage node) is requested, x0 and x1 are the same and pointing to
 * the first stage variable values.
 */
extern "C" typedef int (*str_eval_f_cb)(double* x0, double* x1, double *obj, CallBackDataPtr user_data);

/*
 * computes the constraint value at current x0 and x1, where x0 is the first stage variable values and x1 is the second stage
 * variable values. The results should be written in eq_g (for equality constraints) and inq_g (for inequlity constraints).
 * If the constraints to be evaluated is in the root node, x0 and x1 are the same and pointing to the root node variable values.
 *
 */
extern "C" typedef int (*str_eval_g_cb)(double* x0, double* x1, double* eq_g,double* inq_g, CallBackDataPtr user_data);

/*
 * evaluate the objective gradient sub-vector that declared in the row_node_id in the CallBackData respect to the variables
 * declared in the col_node_id in the CallBackData.
 * For evaluating 1st stage objective gradient, x0 and x1 are the same and pointing to the root node variable values.
 * For evaluating 2st stage objective gradient, x0 and x1 represents the variables values of root stage and second stage variable
 * values. The second stage node is given in the col_node_id.
 *
 */
extern "C" typedef int (*str_eval_grad_f_cb)(double* x0, double* x1, double* vec_grad_f, CallBackDataPtr user_data);

/*
 * This function is to evaluate the structure or values of the sub-matrix of the Jacobian in two separate CCS constructs. One is for the
 * equality constraints and the other is for the inequality constraints.
 * The Jacobian sub-block is identified by the row_node_id and col_node_id given in the CallBackData.
 * When setting e_elts and i_elts to NULL, the function will write back the number of non-zeros to e_nz (for equality constraints)
 * and i_nz (for inequality constraints).
 * When e_elts and i_elts are not NULL, the function should computes the two Jacobian sub-blocks and write them to the values points,
 * (e_elts, e_rowidx, e_colptr) are for the Jacobian block of the equality constraints
 * (i_elts, i_rowidx, i_colptr) are for the Jacobian block of the inequality constraints
 */
extern "C" typedef int (*str_eval_jac_g_cb)(double* x0, double* x1,
		int* e_nz, double* e_elts, int* e_rowidx, int *e_colptr,
		int* i_nz, double* i_elts, int* i_rowidx, int *i_colptr,
		CallBackDataPtr user_data);

/*
 * This function is to evaluate the structure of values of the Hessian of the Lagrangian matrix in CCS format.
 * Since PIPS-NLP is not support for linking constraint, the root stage node only produces the diagonal sub-block for the Hessian
 * of the Lagrangian matrix.
 * 		By setting (0,0) pair for row and col node id in the CallBackData, it should enable this function to compute above mentioned diagonal block.
 *
 * Each second stage constraints can produce three sub-blocks for the Hessian of the Lagrangian matrix, i.e. the diagonal block,
 * the border block and the block that contributes the root diagonal block.
 * 		The diagonal block: this can be requested by setting (1,1), (2,2) and etc. in the CallBackData for row and col node id.
 * 		The border block: this can be requested by setting (0,1), (0,2) and etc. in the CallBackData for row and col node id.
 * 		The root diagonal contribution: this can be requested by setting (1,0), (2,0) and etc. in the CallBackData for row and col node id.
 *
 * Of course for above mentioned blocks, the nz value can be computed by also setting elts to NULL pointer.
 */
extern "C" typedef int (*str_eval_h_cb)(double* x0, double* x1, double* lamdba,
		int* nz, double* elts, int* rowidx, int *colptr,
		CallBackDataPtr user_data);

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

