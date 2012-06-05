#include "ScipLagrangeSolver.hpp"
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <../examples/Queens/src/scip_exception.hpp> // why isn't this in the base?

using namespace std;


ScipLagrangeSolver::ScipLagrangeSolver(stochasticInput &input, int scenarioNumber, const vector<double> &lagrangeMults) {

	absgap = ratio = 0.;
	nvar1 = input.nFirstStageVars();
	int nvar2 = input.nSecondStageVars(scenarioNumber);
	int ncons1 = input.nFirstStageCons();
	int ncons2 = input.nSecondStageCons(scenarioNumber);

	

	SCIP_CALL_EXC( SCIPcreate(&scip) );
	SCIP_CALL_EXC( SCIPincludeDefaultPlugins(scip) );

	// comment to enable output
	SCIP_CALL_EXC( SCIPsetMessagehdlr(NULL) );
	
	SCIP_CALL_EXC( SCIPcreateProb(scip, "lagrange", NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

	vector<double> const obj1 = input.getFirstStageObj(),
		l1 = input.getFirstStageColLB(),
		u1 = input.getFirstStageColUB();
	vector<string> const names1 = input.getFirstStageColNames();
	
	// rescale first-stage objective
	double scale = input.scenarioProbability(scenarioNumber);
	assert(lagrangeMults.size() == static_cast<unsigned>(nvar1));

	for (int i = 0; i < nvar1; i++) {
		SCIP_Var *var;
		SCIP_VARTYPE t;
		if (input.isFirstStageColInteger(i)) {
			if (l1[i] == 0. && u1[i] == 1.) {
				t = SCIP_VARTYPE_BINARY;
			} else {
				t = SCIP_VARTYPE_INTEGER;
			}
		} else {
			t = SCIP_VARTYPE_CONTINUOUS;
		}
		SCIP_CALL_EXC( SCIPcreateVar(scip, &var, names1[i].c_str(), l1[i], u1[i], scale*obj1[i]+lagrangeMults[i], t,
			TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
		SCIP_CALL_EXC( SCIPaddVar(scip, var) );
		vars1.push_back(var);
	}

	vector<double> const obj2 = input.getSecondStageObj(scenarioNumber),
		l2 = input.getSecondStageColLB(scenarioNumber),
		u2 = input.getSecondStageColUB(scenarioNumber);
	vector<string> const names2 = input.getSecondStageColNames(scenarioNumber);
	
	for (int i = 0; i < nvar2; i++) {
		SCIP_Var *var;
		SCIP_VARTYPE t;
		if (input.isSecondStageColInteger(scenarioNumber,i)) {
			if (l2[i] == 0. && u2[i] == 1.) {
				t = SCIP_VARTYPE_BINARY;
			} else {
				t = SCIP_VARTYPE_INTEGER;
			}
		} else {
			t = SCIP_VARTYPE_CONTINUOUS;
		}
		SCIP_CALL_EXC( SCIPcreateVar(scip, &var, names2[i].c_str(), l2[i], u2[i], obj2[i], t,
			TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
		SCIP_CALL_EXC( SCIPaddVar(scip, var) );
		vars2.push_back(var);
	}




	CoinPackedMatrix Arow, Trow, Wrow;
	Arow.reverseOrderedCopyOf(input.getFirstStageConstraints());
	Trow.reverseOrderedCopyOf(input.getLinkingConstraints(scenarioNumber));
	Wrow.reverseOrderedCopyOf(input.getSecondStageConstraints(scenarioNumber));

	vector<double> const lb1 = input.getFirstStageRowLB(),
		ub1 = input.getFirstStageRowUB();
	vector<string> const rownames1 = input.getFirstStageRowNames();

	for (int i = 0; i < ncons1; i++) {
		SCIP_CONS *cons;
		SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &cons, rownames1[i].c_str(), 0, 
			NULL, NULL, lb1[i], ub1[i], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
		
		const double *elts = Arow.getElements();
		const int *idx = Arow.getIndices();
		for (CoinBigIndex k = Arow.getVectorFirst(i); k < Arow.getVectorLast(i); k++) {
			SCIP_CALL_EXC( SCIPaddCoefLinear(scip, cons, vars1[idx[k]], elts[k]) );
		}
		SCIP_CALL_EXC( SCIPaddCons(scip,cons) );
		this->cons.push_back(cons);
	}

	vector<double> const lb2 = input.getSecondStageRowLB(scenarioNumber),
		ub2 = input.getSecondStageRowUB(scenarioNumber);
	vector<string> const rownames2 = input.getSecondStageRowNames(scenarioNumber);

	for (int i = 0; i < ncons2; i++) {
		SCIP_CONS *cons;
		SCIP_CALL_EXC( SCIPcreateConsLinear(scip, &cons, rownames2[i].c_str(), 0, 
			NULL, NULL, lb2[i], ub2[i], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
		
		const double *eltsT = Trow.getElements();
		const int *idxT = Trow.getIndices();
		const double *eltsW = Wrow.getElements();
		const int *idxW = Wrow.getIndices();
		for (CoinBigIndex k = Trow.getVectorFirst(i); k < Trow.getVectorLast(i); k++) {
			SCIP_CALL_EXC( SCIPaddCoefLinear(scip, cons, vars1[idxT[k]], eltsT[k]) );
		}
		for (CoinBigIndex k = Wrow.getVectorFirst(i); k < Wrow.getVectorLast(i); k++) {
			SCIP_CALL_EXC( SCIPaddCoefLinear(scip, cons, vars2[idxW[k]], eltsW[k]) );
		}
		SCIP_CALL_EXC( SCIPaddCons(scip,cons) );
		this->cons.push_back(cons);
	}

}

ScipLagrangeSolver::~ScipLagrangeSolver() {
	
	for (unsigned i = 0; i < vars1.size(); i++) {
		SCIP_CALL_EXC( SCIPreleaseVar(scip, &vars1[i]) );
	}

	for (unsigned i = 0; i < vars2.size(); i++) {
		SCIP_CALL_EXC( SCIPreleaseVar(scip, &vars2[i]) );
	}

	for (unsigned i = 0; i < cons.size(); i++) {
		SCIP_CALL_EXC( SCIPreleaseCons(scip, &cons[i]) );
	}

	SCIP_CALL_EXC( SCIPfree(&scip) );

}

void ScipLagrangeSolver::go() {

	SCIP_CALL_EXC( SCIPsetRealParam(scip, "limits/time", 600) ); // 10 minute time limit
	//SCIP_CALL_EXC( SCIPsetEmphasis(scip,SCIP_PARAMEMPHASIS_HARDLP,true) );
	SCIP_CALL_EXC( SCIPsolve(scip) );
   
}

double ScipLagrangeSolver::getBestPossibleObjective() const {
	return SCIPgetDualbound(scip);
}

double ScipLagrangeSolver::getBestFeasibleObjective() const {
	assert(SCIPisPrimalboundSol(scip));
	return SCIPgetPrimalbound(scip);
}

solverState ScipLagrangeSolver::getStatus() const {
	
	SCIP_STATUS s = SCIPgetStatus(scip);

	switch (s) {
		case SCIP_STATUS_OPTIMAL:
		case SCIP_STATUS_NODELIMIT:
		case SCIP_STATUS_GAPLIMIT:
		case SCIP_STATUS_TIMELIMIT:
		   return Optimal;
		case SCIP_STATUS_INFEASIBLE:
		   return ProvenInfeasible;
		default:
		   assert(0 && "Status unknown");
		   return Initialized;
	}
	
}


vector<double> ScipLagrangeSolver::getBestFirstStageSolution() const {
	SCIP_SOL *sol = SCIPgetBestSol(scip);
	assert(sol);

	vector<double> s(nvar1);
	for (int i = 0; i < nvar1; i++) {
		s[i] = SCIPgetSolVal(scip, sol, vars1[i]);
	}
	return s;
}


