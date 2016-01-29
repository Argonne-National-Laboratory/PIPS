/* PIPS-NLP                                                         	*
 * Authors: Feng Qiang                      		*
 * (C) 2016 Argonne National Laboratory			*/

#ifndef _PARALLELPIPSNLP_C_CALLBACK_H__
#define _PARALLELPIPSNLP_C_CALLBACK_H__


#include "mpi.h"

typedef void * UserDataPtr;
typedef struct CallBackData
{
	UserDataPtr prob;
	int row_node_id;
	int col_node_id;
} CallBackData;

typedef CallBackData * CallBackDataPtr;


extern "C" typedef int (*str_prob_info_cb)(int* n, double* col_lb, double* col_up, int* m, double* row_lb, double* row_up, CallBackData* cbd);
extern "C" typedef int (*str_eval_f_cb)(double* x0, double* x1, double *obj, CallBackDataPtr user_data);
extern "C" typedef int (*str_eval_g_cb)(double* x0, double* x1, double* eq_g,double* inq_g, CallBackDataPtr user_data);
extern "C" typedef int (*str_eval_grad_f_cb)(double* x0, double* x1, double* vec_grad_f, CallBackDataPtr user_data);
extern "C" typedef int (*str_eval_jac_g_cb)(double* x0, double* x1, double* vec_Jac, int* iRows, int *kCols, CallBackDataPtr user_data);
extern "C" typedef int (*str_eval_h_cb)(double* x0, double* x1, double* vec_lambda, double* vec_Hes, int* iRows, int *kCols, CallBackDataPtr user_data);


extern "C"
{

	struct PipsNlpProblemStruct
	{
	  MPI_Comm comm;
	  int nnodes;
	  int n;
	  int m;
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

