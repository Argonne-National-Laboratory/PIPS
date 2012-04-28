#include "CbcRecourseSolver.hpp"

using namespace std;

CbcRecourseSolver::CbcRecourseSolver(stochasticInput &input, int scen, const vector<double>& firstStageSolution) {

	int nvar1 = input.nFirstStageVars();
	int nvar2 = input.nSecondStageVars(scen);
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

	model.loadProblem(Wmat,&collb[0],&colub[0],&obj[0],&rowlb[0],&rowub[0]);

	for (int i = 0; i < nvar2; i++) {
		if (input.isSecondStageColInteger(scen,i)) {
			model.setInteger(i);
		}
	}
	
	model.messageHandler()->setLogLevel(0);
	cbcm.reset(new CbcModel(model));
	CbcMain0(*cbcm);
	cbcm->messageHandler()->setLogLevel(0);

}

void CbcRecourseSolver::go() {
	const char * argv2[]={"","-solve","-quit"};
	//cbcm->setMaximumNodes(1);
	//cbcm->setNumberThreads(2);
	//CbcMain1(3,argv2,*cbcm);
	cbcm->branchAndBound();

}	

double CbcRecourseSolver::getObjective() const {
	return cbcm->getObjValue();
}


solverState CbcRecourseSolver::getStatus() const {

	if (cbcm->isProvenOptimal()) {
		return Optimal;
	} 
	if (cbcm->isProvenInfeasible()) {
		return ProvenInfeasible;
	}
	assert(!cbcm->isAbandoned());

	return Initialized;
}

void CbcRecourseSolver::setSecondStageColState(int idx,variableState s) {
	assert(0);
}

void CbcRecourseSolver::setSecondStageRowState(int idx,variableState s) {
	assert(0);
}

variableState CbcRecourseSolver::getSecondStageColState(int idx) const {
	assert(0); return AtLower;
}

variableState CbcRecourseSolver::getSecondStageRowState(int idx) const {
	assert(0); return AtLower;
}


void CbcRecourseSolver::setDualObjectiveLimit(double d) { 
	model.setDblParam(OsiDualObjectiveLimit,d);
}
