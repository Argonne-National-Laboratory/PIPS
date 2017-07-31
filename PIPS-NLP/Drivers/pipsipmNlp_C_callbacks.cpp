/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#include <cstdlib>

#include "pipsipmNlp_C_callbacks.h"

#include "NlpGenData.h"
#include "NlpGenVars.h"
#include "NlpGenResiduals.h"
#include "NlpInfoCallBack.h"


#include "FilterIPMSolver.h"


#include "NlpGenSparseWithSolver.h"
#include "cNlpGenSparseNLP.h"

#include "Status.h"

#include "pipsOptions.h"

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

extern "C"
struct PipsNlpProblemInfo
{
	int n;
	double* x_L;
	double* x_U;
	int m;
	double* g_L;
	double* g_U;
	int nele_jac;
	int nele_hess;
	eval_f_cb eval_f;
	eval_g_cb eval_g;
	eval_grad_f_cb eval_grad_f;
	eval_jac_g_cb eval_jac_g;
	eval_h_cb eval_h;

	pipsOptions *pipsOpt;
	NlpInfoCallBack *pipsInfoCB;

	double	*Objgrad,  *xlow,  *xupp;
	int 	 nx,my,mz;
	int 	 nnzQ,nnzA,nnzC;
	int 	*irowQ,	*jcolQ;
	int 	*irowA,	*jcolA;
	int 	*irowC,	*jcolC;

	int nxL, nxU, nsL, nsU;
	int nnzCL, nnzCU;

	double	*dQ,	*dA,  *dC;
	double	*bA,  *clow,  *cupp;

	char	*ixlow,	*ixupp;
	char	*iclow,	*icupp;

};



extern "C"
PipsNlpProblem CreatePipsNlpProblem(
		int n,   		int m,
		double* x_L, 	double* x_U,
		double* g_L,  double* g_U,
		int nele_jac, int nele_hess,
		eval_f_cb eval_f,
		eval_g_cb eval_g,
		eval_grad_f_cb eval_grad_f,
		eval_jac_g_cb eval_jac_g,
		eval_h_cb eval_h)
{
	// make sure input is Ok
	if (n<1 || m<0 || !x_L || !x_U || (m>0 && (!g_L || !g_U)) ||
			(m==0 && nele_jac != 0) || (m>0 && nele_jac < 1) || nele_hess < 0 ||
			!eval_f || !eval_grad_f || (m>0 && (!eval_g || !eval_jac_g))) {
		return 0;
	}

	PipsNlpProblem retval = new PipsNlpProblemInfo;

	retval->n = n;
	retval->x_L = new double[n];
	for (int i=0; i<n; i++) {
		retval->x_L[i] = x_L[i];
	}
	retval->x_U = new double[n];
	for (int i=0; i<n; i++) {
		retval->x_U[i] = x_U[i];
	}

	retval->m = m;
	if (m>0) {
		retval->g_L = new double[m];
		for (int i=0; i<m; i++) {
			retval->g_L[i] = g_L[i];
		}
		retval->g_U = new double[m];
		for (int i=0; i<m; i++) {
			retval->g_U[i] = g_U[i];
		}
	}
	else {
		retval->g_L = NULL;
		retval->g_U = NULL;
	}

	retval->nele_jac = nele_jac;
	retval->nele_hess = nele_hess;
	retval->eval_f = eval_f;
	retval->eval_g = eval_g;
	retval->eval_grad_f = eval_grad_f;
	retval->eval_jac_g = eval_jac_g;
	retval->eval_h = eval_h;

	pipsOptions *pipsOpt = new pipsOptions();
	pipsOpt->readFile();
	pipsOpt->defGloOpt();
	if(pipsOpt->UseReducedSpace!=0)
		pipsOpt->UseReducedSpace=1;

	retval->pipsOpt 		= pipsOpt;

	return retval;
}



extern "C"
int PipsNlpSolve( PipsNlpProblem retval, double* obj_val, double* sol_x, UserDataPtr user_data)
{   
	pipsOptions *pipsOpt = retval->pipsOpt;

	double  *Objgrad     = NULL,  *xlow   = NULL,  *xupp = NULL;
	int      nx;
	int     *irowQ = NULL,  *jcolQ  = NULL;
	int      nnzQ;
	int     *irowA = NULL,  *jcolA  = NULL;
	int      nnzA, my;
	int     *irowC = NULL,  *jcolC  = NULL;
	int      nnzC, mz;
	int 	   result;

	int nxL=0, nxU=0, nsL=0, nsU=0;
	int nnzCL, nnzCU;

	double  *dQ    = NULL,  *dA     = NULL,  *dC   = NULL;
	double  *bA     = NULL,  *clow   = NULL,  *cupp = NULL;

	char    *ixlow = NULL,  *ixupp  = NULL;
	char    *iclow = NULL,  *icupp  = NULL;

	int nnzJ;
	int ierr;


	assert(1==pipsOpt->BuildSchurComp && 0==pipsOpt->UseReducedSpace);

	nnzJ = retval->nele_jac;
	nnzQ = retval->nele_hess;

	NlpInfoCallBack *updateNlpCB = new NlpInfoCallBack(retval->eval_f,retval->eval_g,retval->eval_grad_f,
			retval->eval_jac_g,retval->eval_h, user_data);



	// Count the problem sizes, and create a map between the rows of the
	// Jacobian and A and C (A is the Jacobian of the equality constraints,
	// C is the Jacobian of the inequality constraints.)
	// rowmap is moved into getAmplFunction as global variable



	updateNlpCB->_FindRowMap_AddSlack_NY( retval->n, retval->x_L, retval->x_U,
			retval->m, retval->g_L, retval->g_U,
			nnzJ, nx, nnzQ, my, nnzA, mz, nnzC,
			nnzCL, nnzCU, nxL, nxU, nsL, nsU);

	/*
	printf( "CP: my=%d\n",	my);
	printf( "CP: mz=%d\n",  mz);
	printf( "CP: retval->m=%d\n",  retval->m);
	printf( "CP: nnzA=%d\n",	nnzA);
	printf( "CP: nnzC=%d\n",  nnzC);
	printf( "CP: nnzJ=%d\n",  nnzJ);
	printf( "CP: nnzQ=%d\n",  nnzQ);
	*/

	// Note that n_con and n_var are MACROS the translate to
	// asl->n_con and asl->n_var. Arghh....
	assert( my + mz == retval->m );
	assert( nnzA + nnzC == nnzJ);

	updateNlpCB->setBaseInfo(nx,my,mz,nnzQ,nnzA,nnzC,nxL, nxU, nsL, nsU);

	// Allocate the space of the NLP
	newNlpGenSparse( &Objgrad, 	nx,
			&irowQ, nnzQ,  &jcolQ,	&dQ,
			&xlow,		   &ixlow,	&xupp, &ixupp,
			&irowA, nnzA,  &jcolA,	&dA,
			&bA, 	my,
			&irowC, nnzC,  &jcolC,	&dC,
			&clow,	mz,    &iclow,	&cupp, &icupp,
			&ierr );

	if( ierr != 0 ) {
		printf( "Couldn't allocate enough memory\n" );
		exit( 1 );
	}


	updateNlpCB->_get_bounds( retval->x_L, retval->x_U, retval->g_L, retval->g_U,
			xlow, nx, ixlow, xupp, ixupp,
			bA, my,
			clow, mz, iclow, cupp, icupp );

	updateNlpCB->_get_matrices_map(
			nx, nnzQ, my, nnzA, mz, nnzC,
			irowQ, jcolQ, dQ,
			irowA, jcolA, dA,
			irowC, jcolC, dC);


	retval->Objgrad 		= Objgrad;
	retval->nx 		= nx;
	retval->irowQ 		= irowQ;
	retval->nnzQ 		= nnzQ;
	retval->jcolQ 		= jcolQ;
	retval->dQ 		= dQ;
	retval->xlow 		= xlow;
	retval->ixlow 		= ixlow;
	retval->xupp 		= xupp;
	retval->ixupp 		= ixupp;
	retval->irowA 		= irowA;
	retval->nnzA 		= nnzA;
	retval->jcolA 		= jcolA;
	retval->dA 		= dA;
	retval->bA 		= bA;
	retval->my 		= my;
	retval->irowC 		= irowC;
	retval->nnzC 		= nnzC;
	retval->jcolC 		= jcolC;
	retval->dC 		= dC;
	retval->clow 		= clow;
	retval->mz 		= mz;
	retval->iclow 		= iclow;
	retval->cupp 		= cupp;
	retval->icupp 		= icupp;

	retval->nxL = nxL;
	retval->nxU = nxU;
	retval->nsL = nsL;
	retval->nsU = nsU;
	retval->nnzCL = nnzCL;
	retval->nnzCU = nnzCU;

	retval->pipsInfoCB 	= updateNlpCB;



	// add nickname
	PipsNlpProblem pipsnlp = retval;

	if(pipsOpt->prtLvl>0){
		printf("  \n  -----------------------------------------------\n");
		printf("  NLP Solver \n");
		printf("  Nai-Yuan Chiang & V.M. Zavala, Argonne National Laboratory, 2013\n");
		printf("  -----------------------------------------------\n");
		{
			printf("  Variables ================ %-5d \n",pipsnlp->nx);
			printf("  Equality Constraints ===== %-5d \n",pipsnlp->my);
			printf("  Inequality Constraints === %-5d \n",pipsnlp->mz);
			if(pipsnlp->xlow)
				printf("  Variable Lower Bounds ==== %-5d \n",pipsnlp->nxL);
			if(pipsnlp->xupp)
				printf("  Variable Upper Bounds ==== %-5d \n",pipsnlp->nxU);
			if(pipsnlp->clow)
				printf("  Inequality Lower Bounds == %-5d \n",pipsnlp->nsL);
			if(pipsnlp->cupp)
				printf("  Inequality Upper Bounds == %-5d \n\n",pipsnlp->nsU);
		}
		if(pipsnlp->pipsOpt->SymLinearSolver==0)
			std::cout << "\n  Linear system solver ------	 Ma27.\n\n";
		else if(pipsnlp->pipsOpt->SymLinearSolver==1)
			std::cout << "\n  Linear system solver ------	 Ma57.\n\n";
		else
			std::cout << "\n  Give me a linear system solver! Only support Ma27 or Ma57 \n\n";
	}


	//create the NLP formulation factory
	NlpGenSparseWithSolver* nlp  = NULL;
	NlpGenData	  * prob = NULL;
	NlpGenVars     	* vars;
	FilterIPMSolver  * s;
	NlpGenResiduals  	* resid;
	int *rowMap =NULL;




	nlp = new NlpGenSparseWithSolver( pipsnlp->nx, pipsnlp->my, pipsnlp->mz, pipsnlp->nnzQ, pipsnlp->nnzA, pipsnlp->nnzC);
	prob = (NlpGenData * ) nlp->copyDataFromSparseTriple(
			pipsnlp->Objgrad,
			pipsnlp->irowQ,  pipsnlp->nnzQ,   pipsnlp->jcolQ,  pipsnlp->dQ,
			pipsnlp->xlow,   pipsnlp->ixlow,  pipsnlp->xupp,   pipsnlp->ixupp,
			pipsnlp->irowA,  pipsnlp->nnzA,   pipsnlp->jcolA,  pipsnlp->dA,	  pipsnlp->bA,
			pipsnlp->irowC,  pipsnlp->nnzC,   pipsnlp->jcolC,  pipsnlp->dC,
			pipsnlp->clow,   pipsnlp->iclow,  pipsnlp->cupp,   pipsnlp->icupp , rowMap,
			pipsnlp->nxL, pipsnlp->nxU, pipsnlp->nsL, pipsnlp->nsU, pipsnlp->pipsInfoCB);


	vars	= (NlpGenVars * )  nlp->makeVariables( prob);
	resid = (NlpGenResiduals* ) nlp->makeResiduals( prob );
	s 	= new FilterIPMSolver( nlp, prob );

	prob->nxLOri= pipsnlp->nxL;
	prob->nxUOri = pipsnlp->nxU;
	prob->nsLOri = pipsnlp->nsL;
	prob->nsUOri = pipsnlp->nsU;

	s->monitorSelf();

	vars->x->copyFromArray(sol_x);
	//  vars->x->print();
	result = s->solve(prob,vars, resid);

	if(result !=0 ){
		std::cout << "Could not solve the problem.\n";
	} else {

		obj_val[0] = prob->objectiveValue(vars);

//		std::cout << "OK: " << prob->finalIter << " Iters. \n";
		vars->x->copyIntoArray(sol_x);
	}



	// Allocate the space of the NLP
	freeNlpGenSparse( &Objgrad,
			&irowQ,  &jcolQ,	&dQ,
			&xlow,		   &ixlow,	&xupp, &ixupp,
			&irowA,  &jcolA,	&dA,
			&bA,
			&irowC,  &jcolC,	&dC,
			&clow, &iclow,	&cupp, &icupp);


	delete s;
	delete vars;
	delete resid;
	delete prob;
	delete nlp;

	delete updateNlpCB;
	delete pipsOpt;
	return result;
}


extern "C"
void FreePipsNlpProblem(PipsNlpProblem pipsnlp_problem){
	if(pipsnlp_problem->x_L) delete [] pipsnlp_problem->x_L;
	if(pipsnlp_problem->x_U) delete [] pipsnlp_problem->x_U;
	if(pipsnlp_problem->g_L) delete [] pipsnlp_problem->g_L;
	if(pipsnlp_problem->g_U) delete [] pipsnlp_problem->g_U;

	delete pipsnlp_problem;
	pipsnlp_problem=NULL;
}

