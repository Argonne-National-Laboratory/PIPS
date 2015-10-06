/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef _PIPSNLPCINTERFACE_H__
#define _PIPSNLPCINTERFACE_H__

#include "NlpInfoCallBack.h"

extern "C"
{

  struct PipsNlpProblemInfo;

  typedef struct PipsNlpProblemInfo* PipsNlpProblem; 	/** Pointer to a pips_nlp Problem. **/

  PipsNlpProblem CreatePipsNlpProblem(
      int n             	/** Number of variables **/
	, int m				  	/** Number of constraints. **/
    , double* x_L         	/** Variables lower bounds  **/
    , double* x_U         	/** Variables upper bounds **/
    , double* g_L          	/** Constraints lower bounds  **/
    , double* g_U         	/** Constraints upper bounds  **/
    , int nele_jac      	/** Number of Jacobian non-zeros **/
    , int nele_hess     	/** Number of Hessian non-zeros **/
    , eval_f_cb eval_f    	/** Callback function of objective function **/
    , eval_g_cb eval_g    	/** Callback function of constraint body **/
    , eval_grad_f_cb eval_grad_f 	/** Callback function of objective gradient **/
    , eval_jac_g_cb eval_jac_g		/** Callback function of constraint Jacobian **/
    , eval_h_cb eval_h    /** Callback function of Lagrangian Hessian **/
  );

  void FreePipsNlpProblem(PipsNlpProblem pipsnlp_problem);

  int PipsNlpSolve(
      PipsNlpProblem pipsnlp_problem
    , double* opt_obj   
	, double* x		   			/**  Input: initial value; Output: opt solution **/
    , UserDataPtr user_data		/**  Pointer to user data, for the latter use of the callback functions **/
  );

};
#endif

