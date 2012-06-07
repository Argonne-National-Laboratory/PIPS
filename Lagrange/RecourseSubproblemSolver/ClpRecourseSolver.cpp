#include "ClpRecourseSolver.hpp"

using namespace std;

ClpRecourseSolver::ClpRecourseSolver(stochasticInput &input, int scen, const vector<double>& firstStageSolution) {

	int nvar1 = input.nFirstStageVars();
	//int nvar2 = input.nSecondStageVars(scen);
	//int ncons1 = input.nFirstStageCons();
	int ncons2 = input.nSecondStageCons(scen);
	assert(firstStageSolution.size() == static_cast<unsigned>(nvar1));

	const CoinPackedMatrix &Tmat = input.getLinkingConstraints(scen),
		&Wmat = input.getSecondStageConstraints(scen);

	vector<double> Tx(ncons2);
	Tmat.times(&firstStageSolution[0],&Tx[0]);
	const vector<double> &collb = input.getSecondStageColLB(scen),
		&colub = input.getSecondStageColUB(scen),
		&obj = input.getSecondStageObj(scen);
	vector<double> rowlb = input.getSecondStageRowLB(scen),
		rowub = input.getSecondStageRowUB(scen);

	for (int k = 0; k < ncons2; k++) {
		if (rowub[k] < 1e20) {
			rowub[k] -= Tx[k];
		}
		if (rowlb[k] >-1e20) {
			rowlb[k] -= Tx[k];
		}
	}

	solver.loadProblem(Wmat,&collb[0],&colub[0],&obj[0],&rowlb[0],&rowub[0]);
	solver.copyNames(input.getSecondStageRowNames(scen),input.getSecondStageColNames(scen));
	solver.createStatus();

}

void ClpRecourseSolver::go() {
	ClpSolve solvectl;
	// disable presolve so we can re-use optimal bases
	solvectl.setPresolveType(ClpSolve::presolveOff);
	solvectl.setSolveType(ClpSolve::useDual);
	solver.setMaximumSeconds(300);
	solver.initialSolve(solvectl);

}	

double ClpRecourseSolver::getObjective() const {
	return solver.getObjValue();
}

solverState ClpRecourseSolver::getStatus() const {
	switch (solver.status()) {
		case -1:
			return Initialized; // unknown
		case 0:
			return Optimal;
		case 1:
			return ProvenInfeasible;
		case 2:
			return ProvenUnbounded;
		case 3:
			return ProvenInfeasible;
		case 4:
		case 5:
			return Stopped;
		default:
			assert(0 && "Unknown clp status\n");
			return Uninitialized;
	}
}

static variableState clpStatusToState(ClpSimplex::Status s) {
	switch (s) {
		case ClpSimplex::basic:
			return Basic;
		case ClpSimplex::atLowerBound:
			return AtLower;
		case ClpSimplex::atUpperBound:
			return AtUpper;
		case ClpSimplex::isFixed:
			return AtLower;
		case ClpSimplex::isFree:
			return AtLower;
		case ClpSimplex::superBasic:
			printf("Super basic?\n");
			return AtLower;
	}
	assert(0);
	return AtLower;
	
}

static ClpSimplex::Status stateToClpStatus(variableState s) {
	switch (s) {
		case Basic:
			return ClpSimplex::basic;
		case AtLower:
			return ClpSimplex::atLowerBound;
		case AtUpper:
			return ClpSimplex::atUpperBound;
	}
	assert(0);
	return ClpSimplex::atLowerBound;
	
}

void ClpRecourseSolver::setSecondStageColState(int idx,variableState s) {
	solver.setColumnStatus(idx,stateToClpStatus(s));
}

void ClpRecourseSolver::setSecondStageRowState(int idx,variableState s) {
	solver.setRowStatus(idx,stateToClpStatus(s));
}

variableState ClpRecourseSolver::getSecondStageColState(int idx) const {
	return clpStatusToState(solver.getColumnStatus(idx));
}

variableState ClpRecourseSolver::getSecondStageRowState(int idx) const {
	return clpStatusToState(solver.getRowStatus(idx));
}


void ClpRecourseSolver::setDualObjectiveLimit(double d) { 
	solver.setDualObjectiveLimit(d); 
}
