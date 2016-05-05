/* PIPS-NLP                                                             
 * Author: Nai-Yuan Chiang
 * (C) 2015 Argonne National Laboratory
 */

#include "mpi.h"
#include "sNlpInfoFromNL.h"
#include "NlpGenVars.h"
#include "OoqpVector.h"
#include <cmath>

#include "SimpleVector.h"
#include "DoubleMatrixHandle.h"
#include "DoubleMatrix.h"

#include "sData.h"
#include "sVars.h"
#include "sTree.h"

#include <fstream>
#include <sstream>

#include "AmplData_NL.hpp"
#include "amplGenStochInput.hpp"
#include "amplGenStochInput_AddSlack.hpp"
#include "getAmplFunctionNew.h"
#include "../../par_macro.h"

using namespace std;

#include <map>
#include <vector>

#ifdef TIMING
#include "mpi.h"
extern double timeFromAMPL;
#endif

extern int gUseReducedSpace;
extern int gNP_Alg;
extern int gAddSlackParallelSetting;

void getMat_ValOnly(int nnz, double *MatElt, double *dataMatFull,
		std::map<int, int> *GoffMap) {
	if (GoffMap == NULL)
		return;

	map<int, int>::iterator it;

	for (it = GoffMap->begin(); it != GoffMap->end(); it++) {
		MatElt[it->second] = dataMatFull[it->first];
	}
}

sNlpInfoFromNL::sNlpInfoFromNL() {

}

sNlpInfoFromNL::~sNlpInfoFromNL() {

}

sNlpInfoFromNL::sNlpInfoFromNL(int nx_in, int my_in, int mz_in, int nzH_in,
		int nzA_in, int nzC_in) :
		sInfo(nx_in, my_in, mz_in, nzH_in, nzA_in, nzC_in) {
}

sNlpInfoFromNL::sNlpInfoFromNL(int nx_in, int my_in, int mz_in, int nzH_in,
		int nzA_in, int nzC_in, int nxL_in, int nxU_in, int nsL_in, int nsU_in) :
		sInfo(nx_in, my_in, mz_in, nzH_in, nzA_in, nzC_in, nxL_in, nxU_in,
				nsL_in, nsU_in) {
}

sNlpInfoFromNL::sNlpInfoFromNL(sData *data_in, stochasticInput& in) :
		sInfo(data_in) {
	amplGenStochInput & NL_Input = dynamic_cast<amplGenStochInput &>(in);
//  amplGenStochInputV2 & NL_Input = dynamic_cast<amplGenStochInputV2 &>(in);

	data_in->inputNlp = this;
	datarootname = data_in->datarootname;

	asl_local = NL_Input.firstStageData.locASL;
	ObjScal = NL_Input.ObjScale;

	amplRowMap = NL_Input.amplRowMap1st;

	nnzAeqLink = NL_Input.nnzALink1st;
	nnzBeqLoc = NL_Input.nnzALoc1st;
	nnzCineqLink = NL_Input.nnzCLink1st;
	nnzDineqLoc = NL_Input.nnzCLoc1st;

	LocAeqLinkJacGoffMap = NULL;
	LocBeqLocJacGoffMap = NULL;
	LocCineqLinkJacGoffMap = NULL;
	LocDineqLocJacGoffMap = NULL;

	LocQAmatHesGoffMap = &(NL_Input.LocQAmatHesGoffMap[0]);

	nnzQDiag = NL_Input.nnzQ1st;
	nnzQParent = 0; // root node
	nnzQCross = 0; // root node

	iAmDistrib = 0;
	if (MPI_COMM_NULL != mpiComm) {
		int size;
		MPI_Comm_size(mpiComm, &size);
		iAmDistrib = size == 1 ? 0 : 1;
	}

	createChildren(data_in, NL_Input);
}

sNlpInfoFromNL::sNlpInfoFromNL(sData *data_in, stochasticInput& in0,
		const int ChildIdx) :
		sInfo(data_in) {
	data_in->inputNlp = this;
	datarootname = data_in->datarootname;
	amplGenStochInput& in = dynamic_cast<amplGenStochInput&>(in0);

	asl_local = in.localData[ChildIdx].locASL;
	ObjScal = in.ObjScale;

	amplRowMap = in.amplRowMap2nd[ChildIdx];

	LocGloVarMap = &(in.LocGloVarMap[ChildIdx]);
	LocLocVarMap = &(in.LocLocVarMap[ChildIdx]);

	LocWmatJacGoffMap = &(in.LocWmatJacGoffMap[ChildIdx]);
	LocTmatJacGoffMap = &(in.LocTmatJacGoffMap[ChildIdx]);

	LocQAmatHesGoffMap = &(in.LocQAmatHesGoffMap[ChildIdx]);
	LocQWmatHesGoffMap = &(in.LocQWmatHesGoffMap[ChildIdx]);
	LocQTmatHesGoffMap = &(in.LocQTmatHesGoffMap[ChildIdx]);

	nnzQDiag = in.nnzQ2nd[ChildIdx];
	nnzQParent = in.nnzQ1st;
	nnzQCross = in.nnzQCross2nd[ChildIdx];

	nnzAeqLink = in.nnzALink2nd[ChildIdx];
	nnzBeqLoc = in.nnzALoc2nd[ChildIdx];
	nnzCineqLink = in.nnzCLink2nd[ChildIdx];
	nnzDineqLoc = in.nnzCLoc2nd[ChildIdx];

	LocAeqLinkJacGoffMap = &(in.JacALinkGoff2nd[ChildIdx]);
	LocBeqLocJacGoffMap = &(in.JacALocGoff2nd[ChildIdx]);
	LocCineqLinkJacGoffMap = &(in.JacCLinkGoff2nd[ChildIdx]);
	LocDineqLocJacGoffMap = &(in.JacCLocGoff2nd[ChildIdx]);

	if (2 == gUseReducedSpace) {
		data_in->schurVarConID = in.schurVarConIDinNL[ChildIdx];
		data_in->schurSize = in.decisionVarDim[ChildIdx];
	}

	if (gNP_Alg > 0) {
		data_in->var_Part_idx_in = in.NR_partIDX_var[ChildIdx];
		data_in->con_Part_idx_in = in.NR_partIDX_con[ChildIdx];
	}

	iAmDistrib = 0;
	if (MPI_COMM_NULL != mpiComm) {
		int size;
		MPI_Comm_size(mpiComm, &size);
		iAmDistrib = size == 1 ? 0 : 1;
	}

	createChildren(data_in, in);
}

void sNlpInfoFromNL::createChildren(sData *data_in, stochasticInput& in) {
	int mype_;
	MPI_Comm_rank(mpiComm/* MPI_COMM_WORLD*/, &mype_);

	for (size_t it = 0; it < data_in->children.size(); it++) {
		if (stochNode->children[it]->commWrkrs != MPI_COMM_NULL) {
			AddChild(new sNlpInfoFromNL(data_in->children[it], in, it));
		} else {
			AddChild(new sInfoDummy());
		}
		children[it]->parent = this;
//	this->iAmDistrib = 0;
	}
}

sNlpInfoFromNL::sNlpInfoFromNL(sData *data_in, amplGenStochInput_AddSlack& in,
		const int ChildIdx) :
		sInfo(data_in) {
	data_in->inputNlp = this;
	datarootname = data_in->datarootname;

	asl_local = in.localData[ChildIdx].locASL;
	ObjScal = in.ObjScale;

	amplRowMap = in.amplRowMap2nd[ChildIdx];

	LocGloVarMap = &(in.LocGloVarMap[ChildIdx]);
	LocLocVarMap = &(in.LocLocVarMap[ChildIdx]);

	LocWmatJacGoffMap = &(in.LocWmatJacGoffMap[ChildIdx]);
	LocTmatJacGoffMap = NULL;

	LocQAmatHesGoffMap = NULL;
	LocQWmatHesGoffMap = &(in.LocQWmatHesGoffMap[ChildIdx]);
	LocQTmatHesGoffMap = NULL;

	nnzQDiag = in.nnzQ2nd[ChildIdx];
	nnzQParent = in.nnzQ1st;
	nnzQCross = in.nnzQCross2nd[ChildIdx];

	nnzAeqLink = in.nnzALink2nd[ChildIdx];
	nnzBeqLoc = in.nnzALoc2nd[ChildIdx];
	nnzCineqLink = in.nnzCLink2nd[ChildIdx];
	nnzDineqLoc = in.nnzCLoc2nd[ChildIdx];

	LocAeqLinkJacGoffMap = &(in.JacALinkGoff2nd[ChildIdx]);
	LocBeqLocJacGoffMap = &(in.JacALocGoff2nd[ChildIdx]);
	LocCineqLinkJacGoffMap = &(in.JacCLinkGoff2nd[ChildIdx]);
	LocDineqLocJacGoffMap = &(in.JacCLocGoff2nd[ChildIdx]);

	if (2 == gUseReducedSpace) {
		data_in->schurVarConID = in.schurVarConIDinNL[ChildIdx];
		data_in->schurSize = in.decisionVarDim[ChildIdx];
	}

	iAmDistrib = 0;
	if (MPI_COMM_NULL != mpiComm) {
		int size;
		MPI_Comm_size(mpiComm, &size);
		iAmDistrib = size == 1 ? 0 : 1;
	}

	createChildren(data_in, in);
}

void sNlpInfoFromNL::createChildren(sData *data_in,
		amplGenStochInput_AddSlack& in) {
	int mype_;
	MPI_Comm_rank(mpiComm/* MPI_COMM_WORLD*/, &mype_);

	for (size_t it = 0; it < data_in->children.size(); it++) {
		if (stochNode->children[it]->commWrkrs != MPI_COMM_NULL) {
			AddChild(new sNlpInfoFromNL(data_in->children[it], in, it));
		} else {
			AddChild(new sInfoDummy());
		}
		children[it]->parent = this;
//	this->iAmDistrib = 0;
	}
}

double sNlpInfoFromNL::ObjValue(NlpGenVars * vars_) {
	double reVal = 0.0;
	if (gAddSlackParallelSetting == 0)
		reVal = ObjValue_General(vars_);
	else if (gAddSlackParallelSetting == 1)
		reVal = ObjValue_DummyCon(vars_);
	return reVal;
}

double sNlpInfoFromNL::ObjValue_General(NlpGenVars * vars_) {

#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*vars->x).vec);
	double locObj = 0;
	double objWrk;
	double result;

	for (size_t it = 0; it < children.size(); it++)
		locObj += children[it]->ObjValue(vars->children[it]);

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(asl_local);

		double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));

		// form x in AMPL order
		assert(parent->locNx + locNx == n_var);
		OoqpVector* parent_X = (vars_X.parent->vec);
		double *tempParX = (double*) malloc(parent->locNx * sizeof(double));
		double *tempLocX = (double*) malloc(locNx * sizeof(double));
		parent_X->copyIntoArray(tempParX);
		local_X.copyIntoArray(tempLocX);

		map<int, int>::iterator it;
		for (it = LocGloVarMap->begin(); it != LocGloVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempParX[it->second];
		}
		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempLocX[it->second];
		}

		locObj += Ampl_Eval_Obj(asl, tempX_Ampl) * ObjScal;

		free(tempLocX);
		free(tempParX);
		free(tempX_Ampl);
	}

	if (iAmDistrib) {
		MPI_Allreduce(&locObj, &objWrk, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
		locObj = objWrk;
	}

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

	PAR_DEBUG("ObjValue_General - "<<locObj);

	return locObj;

}

// not support seqSparseMatrix yet
double sNlpInfoFromNL::ObjValue_DummyCon(NlpGenVars * vars_) {

#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);

	double locObj = 0;
	double objWrk;

	for (size_t it = 0; it < children.size(); it++)
		locObj += children[it]->ObjValue(vars->children[it]);

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		SimpleVector& simplelocal_X = dynamic_cast<SimpleVector&>(*vars_X.vec);
		locObj += Ampl_Eval_Obj(asl, &simplelocal_X[0]) * ObjScal;
	}

	if (iAmDistrib) {
		MPI_Allreduce(&locObj, &objWrk, 1, MPI_DOUBLE, MPI_SUM, mpiComm);
		locObj = objWrk;
	}

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

	return locObj;

}

void sNlpInfoFromNL::ConstraintBody(NlpGenVars * vars_, OoqpVector *conEq,
		OoqpVector *conIneq) {

	if (gAddSlackParallelSetting == 0)
		ConstraintBody_General(vars_, conEq, conIneq);
	else if (gAddSlackParallelSetting == 1)
		ConstraintBody_DummyCon(vars_, conEq, conIneq);
}

void sNlpInfoFromNL::ConstraintBody_General(NlpGenVars * vars_,
		OoqpVector *conEq, OoqpVector *conIneq) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	int i;

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*vars->x).vec);
	StochVector* local_conEq = dynamic_cast<StochVector*>(conEq);
	StochVector* local_conIneq = dynamic_cast<StochVector*>(conIneq);

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(asl_local);

		double *tempConEq = (double*) malloc(locMy * sizeof(double));
		double *tempConInEq = (double*) malloc(locMz * sizeof(double));
		double *tempCon_Ampl = (double*) malloc(n_con * sizeof(double));
		double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));

		// form x in AMPL order
		assert(parent->locNx + locNx == n_var);
		OoqpVector* parent_X = (vars_X.parent->vec);
		double *tempParX = (double*) malloc(parent->locNx * sizeof(double));
		double *tempLocX = (double*) malloc(locNx * sizeof(double));
		parent_X->copyIntoArray(tempParX);
		local_X.copyIntoArray(tempLocX);

		map<int, int>::iterator it;
		for (it = LocGloVarMap->begin(); it != LocGloVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempParX[it->second];
		}
		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempLocX[it->second];
		}

		Ampl_Eval_Cons(asl, tempX_Ampl, tempCon_Ampl);

		int EqConIndex = 0, IneqIndex = 0, findEq = 0, findIneq = 0;
		for (int iampl = 0; iampl < n_con; iampl++) {
			if (amplRowMap[iampl] < 0) {
				EqConIndex = -(amplRowMap[iampl] + 1); // Recover i from the negative value: equality constraint

				tempConEq[EqConIndex] = tempCon_Ampl[iampl];
				findEq++;
			} else {
				IneqIndex = amplRowMap[iampl]; 			// Inequality constraint
				tempConInEq[IneqIndex] = tempCon_Ampl[iampl];
				findIneq++;
			}
		}
		assert(findIneq == locMz && findEq == locMy);

		local_conEq->vec->copyFromArray(tempConEq);
		local_conIneq->vec->copyFromArray(tempConInEq);

		PRINT_ARRAY("ConstraintBody_General - tempConEq ", tempConEq, locMy);
		PRINT_ARRAY("ConstraintBody_General - tempConInEq", tempConInEq, locMz);

		free(tempLocX);
		free(tempParX);
		free(tempX_Ampl);
		free(tempCon_Ampl);
		free(tempConInEq);
		free(tempConEq);

	}

	for (size_t it = 0; it < children.size(); it++)
		children[it]->ConstraintBody(vars->children[it],
				local_conEq->children[it], local_conIneq->children[it]);

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

}

void sNlpInfoFromNL::ConstraintBody_DummyCon(NlpGenVars * vars_,
		OoqpVector *conEq, OoqpVector *conIneq) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	int i;

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	SimpleVector& simplelocal_X = dynamic_cast<SimpleVector&>(*vars_X.vec);
	StochVector* conEq_st = dynamic_cast<StochVector*>(conEq);
	StochVector* conIneq_st = dynamic_cast<StochVector*>(conIneq);
	SimpleVector& local_conEq = dynamic_cast<SimpleVector&>(*conEq_st->vec);
	SimpleVector& local_conIneq = dynamic_cast<SimpleVector&>(*conIneq_st->vec);

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(asl_local);

//    double *tempConEq   	= (double*) malloc(locMy*sizeof(double));
//    double *tempConInEq 	= (double*) malloc(locMz*sizeof(double));
		double *tempCon_Ampl = (double*) malloc(n_con * sizeof(double));

		SimpleVector& parent_X =
				dynamic_cast<SimpleVector&>(*(vars_X.parent->vec));

		// form x in AMPL order
		assert(locNx == n_var);
		assert(locMy + locMz == n_con + parent->locNx);

		Ampl_Eval_Cons(asl, &simplelocal_X[0], tempCon_Ampl);

		int EqConIndex = 0, IneqIndex = 0, findEq = 0, findIneq = 0;
		for (int iampl = 0; iampl < n_con; iampl++) {
			if (amplRowMap[iampl] < 0) {
				EqConIndex = -(amplRowMap[iampl] + 1); // Recover i from the negative value: equality constraint

				local_conEq[EqConIndex] = tempCon_Ampl[iampl];
				findEq++;
			} else {
				IneqIndex = amplRowMap[iampl]; 			// Inequality constraint
				local_conIneq[IneqIndex] = tempCon_Ampl[iampl];
				findIneq++;
			}
		}
		assert(findIneq == locMz && findEq == locMy - parent->locNx);

		// for dummy constraint
		map<int, int>::iterator itVar;
		for (itVar = LocGloVarMap->begin(); itVar != LocGloVarMap->end();
				itVar++) {
			local_conEq[findEq] = simplelocal_X[itVar->first]
					- parent_X[itVar->second];
			findEq++;
		}
		assert(findEq == locMy);

//    local_conEq->vec->copyFromArray(tempConEq);
//    local_conIneq->vec->copyFromArray(tempConInEq);

//	free(tempLocX);
//    free(tempParX);
//    free(tempX_Ampl);
		free(tempCon_Ampl);
//    free(tempConInEq);	
//    free(tempConEq);  

	}

	for (size_t it = 0; it < children.size(); it++)
		children[it]->ConstraintBody(vars->children[it], conEq_st->children[it],
				conIneq_st->children[it]);

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

}

int sNlpInfoFromNL::ObjGrad(NlpGenVars * vars_, OoqpVector *grad_) {

	if (gAddSlackParallelSetting == 0)
		ObjGrad_General(vars_, grad_);
	else if (gAddSlackParallelSetting == 1)
		ObjGrad_DummyCon(vars_, grad_);
	return 1;
}

void sNlpInfoFromNL::ObjGrad_FromSon(NlpGenVars * vars_, OoqpVector *grad_,
		double* tempfromPar) {
	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*vars->x).vec);
	StochVector* sGrad = dynamic_cast<StochVector*>(grad_);

	assert(children.size() == vars->children.size());

	double *tempGrad = (double*) malloc(locNx * sizeof(double));
	for (int j = 0; j < locNx; j++) {
		tempGrad[j] = 0.;
	}

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);
	assert(mC + mA > 0);

	if (mC + mA > 0) {
		//not the root
		assert(asl_local);

		double *tempGrad_Ampl = (double*) malloc(n_var * sizeof(double));
		double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));

		assert(parent->locNx + locNx == n_var);
		OoqpVector* parent_X = (vars_X.parent->vec);
		double *tempParX = (double*) malloc(parent->locNx * sizeof(double));
		double *tempParGrad = (double*) malloc(parent->locNx * sizeof(double));
		double *tempLocX = (double*) malloc(locNx * sizeof(double));
		parent_X->copyIntoArray(tempParX);
		local_X.copyIntoArray(tempLocX);
		for (int j = 0; j < parent->locNx; j++) {
			tempParGrad[j] = 0.;
		}

		// form x in AMPL order
		map<int, int>::iterator it;
		for (it = LocGloVarMap->begin(); it != LocGloVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempParX[it->second];
		}
		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempLocX[it->second];
		}

		Ampl_Eval_ObjGrad(asl_local, tempX_Ampl, tempGrad_Ampl);

		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempGrad[it->second] = tempGrad_Ampl[it->first];
		}
		for (it = LocGloVarMap->begin(); it != LocGloVarMap->end(); it++) {
			tempParGrad[it->second] = tempGrad_Ampl[it->first];
		}

		for (int j = 0; j < parent->locNx; j++)
			tempfromPar[j] += tempParGrad[j];

		PRINT_ARRAY("ObjGrad_FromSon - tempGrad ", tempGrad, locNx);
		PRINT_ARRAY("ObjGrad_FromSon - tempParGrad ", tempParGrad,
				parent->locNx);

		free(tempParGrad);
		free(tempLocX);
		free(tempParX);
		free(tempX_Ampl);
		free(tempGrad_Ampl);
	}
#if 0
	else {
		// this is root node
		assert(asl_local);
		int mype_;
		MPI_Comm_rank(mpiComm/* MPI_COMM_WORLD*/, &mype_);
		if(mype_== 0) {
			double *tempGrad_Ampl = (double*) malloc(n_var*sizeof(double));
			double *tempX_Ampl = (double*) malloc(n_var*sizeof(double));

			sNlpInfoFromNL* childOne = dynamic_cast<sNlpInfoFromNL*>(children[0]);

			// form x in AMPL order
			assert(locNx + childOne->locNx == n_var);
			OoqpVector* son_X = (vars_X.children[0]->vec);
			double *tempSonX = (double*) malloc(childOne->locNx*sizeof(double));
			double *tempLocX = (double*) malloc(locNx*sizeof(double));
			son_X->copyIntoArray(tempSonX);
			local_X.copyIntoArray(tempLocX);

			map<int,int>::iterator it;
			for(it=childOne->LocGloVarMap->begin(); it!=childOne->LocGloVarMap->end(); it++) {
				tempX_Ampl[it->first] = tempLocX[it->second];
			}
			for(it=childOne->LocLocVarMap->begin(); it!=childOne->LocLocVarMap->end(); it++) {
				tempX_Ampl[it->first] = tempSonX[it->second];
			}

			Ampl_Eval_ObjGrad(asl_local, tempX_Ampl, tempGrad_Ampl);

			free(tempLocX);
			free(tempSonX);
			free(tempX_Ampl);
			free(tempGrad_Ampl);
		}

		MPI_Bcast(&tempGrad[0],locNx,MPI_DOUBLE,0,mpiComm);

	}
#endif

	sGrad->vec->copyFromArray(tempGrad);

	sGrad->vec->scale(ObjScal);

	free(tempGrad);

	assert(children.size() == 0);
	for (size_t it = 0; it < children.size(); it++) {
		PAR_DEBUG("it - "<<it);
		children[it]->ObjGrad(vars->children[it], sGrad->children[it]);
	}

}

void sNlpInfoFromNL::ObjGrad_General(NlpGenVars * vars_, OoqpVector *grad_) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*vars->x).vec);
	StochVector* sGrad = dynamic_cast<StochVector*>(grad_);

	assert(children.size() == vars->children.size());

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);
	assert(nA == 0 && nC == 0); // this is root!

	double *tempGrad = (double*) malloc(locNx * sizeof(double));
	double *tempFromSonGrad = (double*) malloc(locNx * sizeof(double));
	for (int j = 0; j < locNx; j++) {
		tempGrad[j] = 0.;
		tempFromSonGrad[j] = 0.;
	}

	for (size_t it = 0; it < children.size(); it++) {
		children[it]->ObjGrad_FromSon(vars->children[it], sGrad->children[it],
				tempFromSonGrad);
	}

	MPI_Allreduce(tempFromSonGrad, tempGrad, locNx, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

	PRINT_ARRAY("ObjGrad_General - tempGrad - ", tempGrad, locNx);
	PRINT_ARRAY("ObjGrad_General - tempFromSonGrad - ", tempFromSonGrad, locNx);

	sGrad->vec->copyFromArray(tempGrad);
	sGrad->vec->scale(ObjScal);

	free(tempFromSonGrad);
	free(tempGrad);

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

}

void sNlpInfoFromNL::ObjGrad_DummyCon(NlpGenVars * vars_, OoqpVector *grad_) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	SimpleVector& simplelocal_X = dynamic_cast<SimpleVector&>(*vars_X.vec);

	StochVector* sGrad = dynamic_cast<StochVector*>(grad_);
	SimpleVector& simplelocal_Grad = dynamic_cast<SimpleVector&>(*sGrad->vec);

	assert(children.size() == vars->children.size());

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(asl_local);
		Ampl_Eval_ObjGrad(asl_local, &simplelocal_X[0], &simplelocal_Grad[0]);
		simplelocal_Grad.scale(ObjScal);
	}

	for (size_t it = 0; it < children.size(); it++) {
		children[it]->ObjGrad(vars->children[it], sGrad->children[it]);
	}

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

}

void sNlpInfoFromNL::JacFull(NlpGenVars * vars_, GenMatrix* JacA,
		GenMatrix* JacC) {
	if (gAddSlackParallelSetting == 0)
		JacFull_General(vars_, JacA, JacC);
	else if (gAddSlackParallelSetting == 1)
		JacFull_DummyCon(vars_, JacA, JacC);
}

// Note that we assume 1st B&D are empty
void sNlpInfoFromNL::JacFull_General(NlpGenVars * vars_, GenMatrix* JacA,
		GenMatrix* JacC) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	OoqpVector* local_X = (vars_X.vec);
	assert(children.size() == vars->children.size());

	ASL_pfgh *asl = asl_local;

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		double *tempFullJac_Ampl = (double*) malloc(nzc * sizeof(double));
		double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));

		// form x in AMPL order
		assert(parent->locNx + locNx == n_var);
		OoqpVector* parent_X = (vars_X.parent->vec);
		double *tempParX = (double*) malloc(parent->locNx * sizeof(double));
		double *tempLocX = (double*) malloc(locNx * sizeof(double));
		parent_X->copyIntoArray(tempParX);
		local_X->copyIntoArray(tempLocX);

		map<int, int>::iterator it;
		for (it = LocGloVarMap->begin(); it != LocGloVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempParX[it->second];
		}
		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempLocX[it->second];
		}
		PRINT_ARRAY("JacFull_General - tempParX ", tempParX, parent->locNx);
		PRINT_ARRAY("JacFull_General - tempLocX ", tempLocX, locNx);

		Ampl_Eval_Jac(asl, tempX_Ampl, tempFullJac_Ampl);

		assert(Amat->numberOfNonZeros() == nnzAeqLink);
		assert(Bmat->numberOfNonZeros() == nnzBeqLoc);
		assert(Cmat->numberOfNonZeros() == nnzCineqLink);
		assert(Dmat->numberOfNonZeros() == nnzDineqLoc);

		assert(nnzAeqLink + nnzBeqLoc + nnzCineqLink + nnzDineqLoc == nzc);

		// A&C are linking Jac for Equallity constraint, B&D are loc Jac for InEquallity constraint
		double *tempAmat = (double*) malloc(
				Amat->numberOfNonZeros() * sizeof(double));
		double *tempBmat = (double*) malloc(
				Bmat->numberOfNonZeros() * sizeof(double));
		double *tempCmat = (double*) malloc(
				Cmat->numberOfNonZeros() * sizeof(double));
		double *tempDmat = (double*) malloc(
				Dmat->numberOfNonZeros() * sizeof(double));

		getMat_ValOnly(nnzAeqLink, tempAmat, tempFullJac_Ampl,
				LocAeqLinkJacGoffMap);
		getMat_ValOnly(nnzBeqLoc, tempBmat, tempFullJac_Ampl,
				LocBeqLocJacGoffMap);
		getMat_ValOnly(nnzCineqLink, tempCmat, tempFullJac_Ampl,
				LocCineqLinkJacGoffMap);
		getMat_ValOnly(nnzDineqLoc, tempDmat, tempFullJac_Ampl,
				LocDineqLocJacGoffMap);

		Amat->copyMtxFromDouble(Amat->numberOfNonZeros(), tempAmat);
		Bmat->copyMtxFromDouble(Bmat->numberOfNonZeros(), tempBmat);
		Cmat->copyMtxFromDouble(Cmat->numberOfNonZeros(), tempCmat);
		Dmat->copyMtxFromDouble(Dmat->numberOfNonZeros(), tempDmat);

		PRINT_ARRAY("JacFull_General - tempAmat ", tempAmat, nnzAeqLink);
		PRINT_ARRAY("JacFull_General - tempBmat ", tempBmat, nnzBeqLoc);
		PRINT_ARRAY("JacFull_General - tempCmat ", tempCmat, nnzCineqLink);
		PRINT_ARRAY("JacFull_General - tempDmat ", tempDmat, nnzDineqLoc);

		free(tempDmat);
		free(tempCmat);
		free(tempBmat);
		free(tempAmat);
		free(tempLocX);
		free(tempParX);
		free(tempX_Ampl);
		free(tempFullJac_Ampl);

	} else {

	}

	for (size_t it = 0; it < children.size(); it++) {
		PAR_DEBUG("it - "<<it);
		children[it]->JacFull(vars->children[it], NULL, NULL);
	}
#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

}

// Note that we assume 1st B&D are empty
void sNlpInfoFromNL::JacFull_DummyCon(NlpGenVars * vars_, GenMatrix* JacA,
		GenMatrix* JacC) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	SimpleVector& simplelocal_X = dynamic_cast<SimpleVector&>(*vars_X.vec);

	assert(children.size() == vars->children.size());

	ASL_pfgh *asl = asl_local;

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(locNx == n_var);
		assert(locMy + locMz == n_con + parent->locNx);

		double *tempFullJac_Ampl = (double*) malloc(nzc * sizeof(double));
		Ampl_Eval_Jac(asl, &simplelocal_X[0], tempFullJac_Ampl);

		assert(Amat->numberOfNonZeros() == nnzAeqLink);
		assert(Bmat->numberOfNonZeros() == nnzBeqLoc);
		assert(Cmat->numberOfNonZeros() == nnzCineqLink);
		assert(Dmat->numberOfNonZeros() == nnzDineqLoc);

		assert(nnzBeqLoc - parent->locNx + nnzDineqLoc == nzc);
		assert(nnzAeqLink == parent->locNx);
		assert(nnzCineqLink == 0);

		// A&C are linking Jac for Equallity constraint, B&D are loc Jac for InEquallity constraint
//	double *tempAmat = (double*) malloc(nnzAeqLink*sizeof(double));
//    double *tempBmat = (double*) malloc((nnzBeqLoc-parent->locNx)*sizeof(double));
//    double *tempCmat = (double*) malloc(Cmat->numberOfNonZeros()*sizeof(double));
//	double *tempDmat = (double*) malloc(0*sizeof(double));

		getMat_ValOnly(nnzBeqLoc - parent->locNx, Bmat->M(), tempFullJac_Ampl,
				LocBeqLocJacGoffMap);
		getMat_ValOnly(nnzDineqLoc, Dmat->M(), tempFullJac_Ampl,
				LocDineqLocJacGoffMap);

//	getMat_ValOnly(nnzAeqLink,tempAmat,tempFullJac_Ampl,LocAeqLinkJacGoffMap);
//	getMat_ValOnly(nnzBeqLoc,tempBmat,tempFullJac_Ampl,LocBeqLocJacGoffMap);
//	getMat_ValOnly(nnzCineqLink,tempCmat,tempFullJac_Ampl,LocCineqLinkJacGoffMap);	
//	getMat_ValOnly(nnzDineqLoc,tempDmat,tempFullJac_Ampl,LocDineqLocJacGoffMap);

//	Amat->copyMtxFromDouble(Amat->numberOfNonZeros(),tempAmat);
//	Bmat->copyMtxFromDouble(Bmat->numberOfNonZeros(),tempBmat);
//	Cmat->copyMtxFromDouble(Cmat->numberOfNonZeros(), tempCmat);
//	Dmat->copyMtxFromDouble(Dmat->numberOfNonZeros(), tempDmat);

//	free(tempDmat);	
//	free(tempCmat);
//	free(tempBmat);
//	free(tempAmat);
//	free(tempLocX);	
//	free(tempParX);	
//	free(tempX_Ampl);
		free(tempFullJac_Ampl);

	}

	for (size_t it = 0; it < children.size(); it++)
		children[it]->JacFull(vars->children[it], NULL, NULL);

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 

}

void sNlpInfoFromNL::Hessian(NlpGenVars * vars_, SymMatrix *Hess) {
	if (gAddSlackParallelSetting == 0)
		Hessian_General(vars_, Hess);
	else if (gAddSlackParallelSetting == 1)
		Hessian_DummyCon(vars_, Hess);
}

void sNlpInfoFromNL::Hessian_FromSon(NlpGenVars * vars_, double* tempfromPar) {

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	OoqpVector* local_X = (vars_X.vec);
	StochVector& vars_Y = dynamic_cast<StochVector&>(*vars->y);
	OoqpVector* local_Y = vars_Y.vec;
	StochVector& vars_Z = dynamic_cast<StochVector&>(*vars->z);
	OoqpVector* local_Z = vars_Z.vec;

	assert(children.size() == vars->children.size());

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mA + mC > 0) {
		//not the root
		double *tempY = (double*) malloc(locMy * sizeof(double));
		double *tempZ = (double*) malloc(locMz * sizeof(double));

		double *tempDual_Ampl = (double*) malloc(n_con * sizeof(double));
		double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));
		double *tempH_Ampl = (double*) malloc(
				(nnzQDiag + nnzQCross + nnzQParent) * sizeof(double));

		// form x in AMPL order
		assert(parent->locNx + locNx == n_var);
		OoqpVector* parent_X = (vars_X.parent->vec);
		double *tempParX = (double*) malloc(parent->locNx * sizeof(double));
		double *tempLocX = (double*) malloc(locNx * sizeof(double));
		parent_X->copyIntoArray(tempParX);
		local_X->copyIntoArray(tempLocX);

		map<int, int>::iterator it;
		for (it = LocGloVarMap->begin(); it != LocGloVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempParX[it->second];
		}
		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempLocX[it->second];
		}

		local_Y->negate();
		local_Y->copyIntoArray(tempY);
		local_Y->negate();
		local_Z->negate();
		local_Z->copyIntoArray(tempZ);
		local_Z->negate();

		int EqConIndex = 0, IneqIndex = 0;
		// form dualVar in AMPL order
		for (int iampl = 0; iampl < n_con; iampl++) {
			if (amplRowMap[iampl] < 0) {
				EqConIndex = -(amplRowMap[iampl] + 1); // Recover i from the negative value: equality constraint
				tempDual_Ampl[iampl] = tempY[EqConIndex];
			} else {
				IneqIndex = amplRowMap[iampl];		  	// Inequality constraint
				tempDual_Ampl[iampl] = tempZ[IneqIndex];
			}
		}

		Ampl_Eval_Hessian_Tri(asl, tempDual_Ampl, tempH_Ampl, ObjScal);

		double *tempQDiag = (double*) malloc(nnzQDiag * sizeof(double));
		double *tempQCross = (double*) malloc(nnzQCross * sizeof(double));
		double *tempQPar = (double*) malloc(nnzQParent * sizeof(double));

		getMat_ValOnly(nnzQDiag, tempQDiag, tempH_Ampl, LocQWmatHesGoffMap);
		getMat_ValOnly(nnzQCross, tempQCross, tempH_Ampl, LocQTmatHesGoffMap);
		getMat_ValOnly(nnzQParent, tempQPar, tempH_Ampl, LocQAmatHesGoffMap);

		assert(
				Qdiag->numberOfNonZeros() == nnzQDiag
						&& Qborder->numberOfNonZeros() == nnzQCross);
		Qdiag->copyMtxFromDouble(Qdiag->numberOfNonZeros(), tempQDiag);
		Qborder->copyMtxFromDouble(Qborder->numberOfNonZeros(), tempQCross);
		for (int j = 0; j < nnzQParent; j++)
			tempfromPar[j] += tempQPar[j];

		PRINT_ARRAY("Hessian_FromSon - tempQDiag - ", tempQDiag, nnzQDiag);
		PRINT_ARRAY("Hessian_FromSon - tempQCross - ", tempQCross, nnzQCross);
		PRINT_ARRAY("Hessian_FromSon - tempQPar - ", tempQPar, nnzQParent);

		free(tempQPar);
		free(tempQCross);
		free(tempQDiag);
		free(tempLocX);
		free(tempParX);
		free(tempX_Ampl);
		free(tempDual_Ampl);
		free(tempH_Ampl);
		free(tempZ);
		free(tempY);

	}

	// assume only 1 level
	assert(children.size() == 0);
//  for(size_t it=0; it<children.size(); it++){
//    children[it]->Hessian_FromSon(vars->children[it],NULL);
//  }

}

void sNlpInfoFromNL::Hessian_General(NlpGenVars * vars_, SymMatrix *Hess) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	StochVector& vars_Y = dynamic_cast<StochVector&>(*vars->y);
	StochVector& vars_Z = dynamic_cast<StochVector&>(*vars->z);
	OoqpVector* local_X = vars_X.vec;
	OoqpVector* local_Y = vars_Y.vec;
	OoqpVector* local_Z = vars_Z.vec;

	assert(children.size() == vars->children.size());

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);
	assert(nA == 0 && nC == 0); // this is root!

	assert(Qdiag->numberOfNonZeros() == nnzQDiag);
	double *tempQDiag = (double*) malloc(nnzQDiag * sizeof(double));
	double *tempFromSonH = (double*) malloc(nnzQDiag * sizeof(double));

	for (int j = 0; j < nnzQDiag; j++) {
		tempQDiag[j] = 0;
		tempFromSonH[j] = 0;
	}

	for (size_t it = 0; it < children.size(); it++) {
		PAR_DEBUG("it - "<<it);
		children[it]->Hessian_FromSon(vars->children[it], tempFromSonH);
	}

	MPI_Allreduce(tempFromSonH, tempQDiag, nnzQDiag, MPI_DOUBLE, MPI_SUM,
			MPI_COMM_WORLD);

#if 0
	int mype_;
	MPI_Comm_rank(mpiComm/* MPI_COMM_WORLD*/, &mype_);

	if(mype_== 0) {
		assert(locMy==0 && locMz==0); // this is root!
		//use 1st scenario to get info
		sNlpInfoFromNL* childOne = dynamic_cast<sNlpInfoFromNL*>(children[0]);

		assert(locNx + childOne->locNx == n_var);
		assert(childOne->nnzQParent = nnzQDiag);
		double *tempDual_Ampl = (double*) malloc(n_con*sizeof(double));
		double *tempH_Ampl = (double*) malloc((childOne->nnzQDiag+childOne->nnzQCross+childOne->nnzQParent)*sizeof(double));

		for(int j=0; j<n_con;j++) {
			tempDual_Ampl[j]=0.;
		}

		/*
		 OoqpVector* son_Y = (vars_Y.children[0]->vec);
		 OoqpVector* son_Z = (vars_Z.children[0]->vec);
		 double *tempSonY = (double*) malloc(childOne->locMy*sizeof(double));
		 double *tempSonZ = (double*) malloc(childOne->locMz*sizeof(double));

		 son_Y->negate();    son_Y->copyIntoArray(tempSonY);    son_Y->negate();
		 son_Z->negate();    son_Z->copyIntoArray(tempSonZ);    son_Z->negate();

		 int EqConIndex=0, IneqIndex=0;
		 // form dualVar in AMPL order
		 for( int iampl = 0; iampl < n_con; iampl++ ) {
		 if( childOne->amplRowMap[iampl]<0 ) {
		 EqConIndex = - ( childOne->amplRowMap[iampl] + 1 );   // Recover i from the negative value: equality constraint
		 tempDual_Ampl[iampl] = tempSonY[EqConIndex];
		 } else {
		 IneqIndex = amplRowMap[iampl];		  		// Inequality constraint
		 tempDual_Ampl[iampl] = tempSonZ[IneqIndex];
		 }
		 }
		 */
		Ampl_Eval_Hessian_Tri(asl, tempDual_Ampl, tempH_Ampl, 1);

		double *tempQPar = (double*) malloc(childOne->nnzQParent * sizeof(double));
		for(int j=0;j<nnzQDiag;j++)
		tempQPar[j] = 0.0;

		getMat_ValOnly(nnzQDiag, tempQPar, tempH_Ampl,childOne->LocQAmatHesGoffMap);

		for(int j=0;j<nnzQDiag;j++)
		tempQDiag[j] += tempQPar[j];

		free(tempQPar);
		free(tempH_Ampl);
		free(tempDual_Ampl);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&tempQDiag[0],nnzQDiag,MPI_DOUBLE,0,mpiComm);
	MPI_Barrier(MPI_COMM_WORLD);

#endif

	assert(Qdiag->numberOfNonZeros() == nnzQDiag);
	Qdiag->copyMtxFromDouble(Qdiag->numberOfNonZeros(), tempQDiag);

	PRINT_ARRAY("Hessian_General - tempFromSonH - ", tempFromSonH, nnzQDiag);
	PRINT_ARRAY("Hessian_General - tempQDiag - ", tempQDiag, nnzQDiag);

	free(tempFromSonH);
	free(tempQDiag);

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 
}

void sNlpInfoFromNL::Hessian_DummyCon(NlpGenVars * vars_, SymMatrix *Hess) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	StochVector& vars_Y = dynamic_cast<StochVector&>(*vars->y);
	StochVector& vars_Z = dynamic_cast<StochVector&>(*vars->z);

	SimpleVector& simplelocal_X = dynamic_cast<SimpleVector&>(*vars_X.vec);
	SimpleVector& simplelocal_Y = dynamic_cast<SimpleVector&>(*vars_Y.vec);
	SimpleVector& simplelocal_Z = dynamic_cast<SimpleVector&>(*vars_Z.vec);

	assert(children.size() == vars->children.size());

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mA + mC > 0) {
		//not the root
		assert(locNx == n_var);
		assert(locMy + locMz == n_con + parent->locNx);
		assert(nnzQCross + nnzQParent == 0);

		double *tempDual_Ampl = (double*) malloc((n_con) * sizeof(double));
		double *tempH_Ampl = (double*) malloc((nnzQDiag) * sizeof(double));

		// form x in AMPL order

		int EqConIndex = 0, IneqIndex = 0;
		// form dualVar in AMPL order
		for (int iampl = 0; iampl < n_con; iampl++) {
			if (amplRowMap[iampl] < 0) {
				EqConIndex = -(amplRowMap[iampl] + 1); // Recover i from the negative value: equality constraint
				tempDual_Ampl[iampl] = -simplelocal_Y[EqConIndex];
			} else {
				IneqIndex = amplRowMap[iampl];		  	// Inequality constraint
				tempDual_Ampl[iampl] = -simplelocal_Z[IneqIndex];
			}
		}

		Ampl_Eval_Hessian_Tri(asl, tempDual_Ampl, Qdiag->M(), ObjScal);

//	double *tempQDiag  = (double*) malloc(nnzQDiag   * sizeof(double));	
//    double *tempQCross = (double*) malloc(nnzQCross  * sizeof(double));
//    double *tempQPar   = (double*) malloc(nnzQParent * sizeof(double));
//	getMat_ValOnly(nnzQDiag,Qdiag->M(),tempH_Ampl,LocQWmatHesGoffMap); 

//	Qdiag->copyMtxFromDouble(Qdiag->numberOfNonZeros(),tempH);

		assert(
				Qdiag->numberOfNonZeros() == nnzQDiag
						&& Qborder->numberOfNonZeros() == 0);

//    free(tempQPar);
//	free(tempQCross);
//    free(tempQDiag);
//    free(tempLocX);
//	free(tempParX);
		free(tempDual_Ampl);
		free(tempH_Ampl);
//    free(tempZ);
//    free(tempY);

	}

//  assert(Qdiag->numberOfNonZeros()==nnzQDiag);
//  Qdiag->copyMtxFromDouble(Qdiag->numberOfNonZeros(),tempQDiag);

//  free(tempFromSonH);
//  free(tempQDiag);

	for (size_t it = 0; it < children.size(); it++) {
		children[it]->Hessian(vars->children[it], Hess);
	}

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 
}

void sNlpInfoFromNL::get_InitX0(OoqpVector* vX) {

	if (gAddSlackParallelSetting == 0)
		get_InitX0_General(vX);
	else if (gAddSlackParallelSetting == 1)
		get_InitX0_DummyCon(vX);
}

void sNlpInfoFromNL::get_InitX0_General(OoqpVector* vX) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	StochVector* vars_X = dynamic_cast<StochVector*>(vX);
	OoqpVector* local_X = vars_X->vec;
	double *tempX = (double*) malloc(locNx * sizeof(double));

	assert(children.size() == vars_X->children.size());

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(asl_local);

		double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));
		Ampl_Eval_InitX0(asl_local, tempX_Ampl);

		map<int, int>::iterator it;
		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempX[it->second] = tempX_Ampl[it->first];
		}
		free(tempX_Ampl);
	} else {
		// this is root node
		assert(asl_local);
		int mype_;
		MPI_Comm_rank(mpiComm/* MPI_COMM_WORLD*/, &mype_);
		if (mype_ == 0) {
			double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));
			sNlpInfoFromNL* childOne =
					dynamic_cast<sNlpInfoFromNL*>(children[0]);

			Ampl_Eval_InitX0(asl_local, tempX_Ampl);

			map<int, int>::iterator it;
			for (it = childOne->LocGloVarMap->begin();
					it != childOne->LocGloVarMap->end(); it++) {
				tempX[it->second] = tempX_Ampl[it->first];
			}
			free(tempX_Ampl);
		}

		MPI_Bcast(&tempX[0], locNx, MPI_DOUBLE, 0, mpiComm);

	}

	local_X->copyFromArray(tempX);
	free(tempX);

	for (size_t it = 0; it < children.size(); it++)
		children[it]->get_InitX0(vars_X->children[it]);

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 	

}

void sNlpInfoFromNL::get_InitX0_DummyCon(OoqpVector* vX) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif

	ASL_pfgh *asl = asl_local;

	StochVector* vars_X = dynamic_cast<StochVector*>(vX);
	OoqpVector* local_X = vars_X->vec;
	double *tempX = (double*) malloc(locNx * sizeof(double));
	;

	assert(children.size() == vars_X->children.size());

	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(locNx == n_var);
		Ampl_Eval_InitX0(asl_local, tempX);
		local_X->copyFromArray(tempX);
	} else {
		// this is root node, use the 1st scenatio info to get initial point
		assert(asl != NULL);
//	asl = ((sNlpInfoFromNL*)children[0])->asl_local;

		int mype_;
		MPI_Comm_rank(mpiComm/* MPI_COMM_WORLD*/, &mype_);
		if (mype_ == 0) {
			double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));
			sNlpInfoFromNL* childOne =
					dynamic_cast<sNlpInfoFromNL*>(children[0]);

			Ampl_Eval_InitX0(asl, tempX_Ampl);

			map<int, int>::iterator it;
			for (it = childOne->LocGloVarMap->begin();
					it != childOne->LocGloVarMap->end(); it++) {
				tempX[it->second] = tempX_Ampl[it->first];
			}
			free(tempX_Ampl);
		}

		MPI_Bcast(&tempX[0], locNx, MPI_DOUBLE, 0, mpiComm);

	}

	local_X->copyFromArray(tempX);
	free(tempX);

	for (size_t it = 0; it < children.size(); it++)
		children[it]->get_InitX0(vars_X->children[it]);

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 	

}

void sNlpInfoFromNL::writeSolution(NlpGenVars * vars_) {
#ifdef TIMING
	double tTot=MPI_Wtime();
#endif
	ASL_pfgh *asl = asl_local;
	sVars * vars = dynamic_cast<sVars*>(vars_);
	StochVector& vars_X = dynamic_cast<StochVector&>(*vars->x);
	OoqpVector& local_X = *(dynamic_cast<StochVector&>(*vars->x).vec);
	OoqpVector& local_Y = *(dynamic_cast<StochVector&>(*vars->y).vec);
	OoqpVector& local_Z = *(dynamic_cast<StochVector&>(*vars->z).vec);

	for (size_t it = 0; it < children.size(); it++) {
		children[it]->writeSolution(vars->children[it]);
	}
	long long mA, nA, mC, nC;
	Amat->getSize(mA, nA);
	Cmat->getSize(mC, nC);

	if (mC + mA > 0) {
		//not the root
		assert(asl_local);

		double *tempX_Ampl = (double*) malloc(n_var * sizeof(double));

		// form x in AMPL order
		assert(parent->locNx + locNx == n_var);
		OoqpVector* parent_X = (vars_X.parent->vec);
		double *tempParX = (double*) malloc(parent->locNx * sizeof(double));
		double *tempLocX = (double*) malloc(locNx * sizeof(double));

		parent_X->copyIntoArray(tempParX);
		local_X.copyIntoArray(tempLocX);

		map<int, int>::iterator it;
		for (it = LocGloVarMap->begin(); it != LocGloVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempParX[it->second];
		}
		for (it = LocLocVarMap->begin(); it != LocLocVarMap->end(); it++) {
			tempX_Ampl[it->first] = tempLocX[it->second];
		}

		double *tempY = (double*) malloc(locMy * sizeof(double));
		double *tempZ = (double*) malloc(locMz * sizeof(double));
		double *dualWrk = (double*) malloc(n_con * sizeof(double));
		local_Y.copyIntoArray(tempY);
		local_Z.copyIntoArray(tempZ);
		int EqConIndex = 0, IneqIndex = 0, findEq = 0, findIneq = 0;
		for (int iampl = 0; iampl < n_con; iampl++) {
			if (amplRowMap[iampl] < 0) {
				EqConIndex = -(amplRowMap[iampl] + 1); // Recover i from the negative value: equality constraint
				if (tempY)
					dualWrk[iampl] = tempY[EqConIndex];
				else
					dualWrk[iampl] = 0;
				findEq++;
			} else {
				IneqIndex = amplRowMap[iampl]; // Inequality constraint
				if (tempZ)
					dualWrk[iampl] = tempZ[IneqIndex];
				else
					dualWrk[iampl] = 0;
				findIneq++;
			}
		}
		assert(findIneq == locMz && findEq == locMy);
		ampl_write_solution(asl, tempX_Ampl, dualWrk);

		free(dualWrk);
		free(tempZ);
		free(tempY);
		free(tempLocX);
		free(tempParX);
		free(tempX_Ampl);
	}

#ifdef TIMING
	timeFromAMPL += MPI_Wtime()-tTot;
#endif 
}

