#include "CbcLagrangeSolver.hpp"

using namespace std;

namespace{
template<typename T1, typename T2> void concatenateOne(stochasticInput &data, 
			T1 &out, 
			T2 (stochasticInput::*first)(),
			T2 (stochasticInput::*second)(int), int scen) {

	const T2& arr1 = (data.*first)();
	for (unsigned k = 0; k < arr1.size(); k++) {
		out[k] = arr1[k];
	}
	unsigned r = arr1.size();
	const T2& arr2 = (data.*second)(scen);
	for (unsigned k = 0; k < arr2.size(); k++) {
		out[r++] = arr2[k];
	}

}
}

CbcLagrangeSolver::CbcLagrangeSolver(stochasticInput &input, int scenarioNumber, const vector<double> &lagrangeDiff) {

	ratio = 0.;
	nvar1 = input.nFirstStageVars();
	int nvar2 = input.nSecondStageVars(scenarioNumber);
	int ncons1 = input.nFirstStageCons();
	int ncons2 = input.nSecondStageCons(scenarioNumber);

	const CoinPackedMatrix &Amat = input.getFirstStageConstraints(),
		&Tmat = input.getLinkingConstraints(scenarioNumber),
		&Wmat = input.getSecondStageConstraints(scenarioNumber);
	
	
	int totalVar = nvar1 + nvar2;
	int totalCons = ncons1 + ncons2;

	CoinBigIndex totalNnz = Amat.getNumElements() + Tmat.getNumElements() + Wmat.getNumElements();

	// CoinPackedMatrix takes ownership of these, so we don't free them
	CoinBigIndex *starts = new CoinBigIndex[totalVar+1];
	double *elts = new double[totalNnz];
	int *rowIdx = new int[totalNnz];
	
	CoinBigIndex nnz = 0, start, end;
	// put first-stage variables first, as is customary
	int const *Aidx = Amat.getIndices();
	double const *Aelts = Amat.getElements();
	int const *Tidx = Tmat.getIndices();
	double const *Telts = Tmat.getElements();

	for (int c = 0; c < nvar1; c++) {
		starts[c] = nnz;
		start = Amat.getVectorFirst(c);
		end = Amat.getVectorLast(c);
		for (CoinBigIndex j = start; j < end; j++) {
			elts[nnz] = Aelts[j];
			rowIdx[nnz++] = Aidx[j];
		}
		start = Tmat.getVectorFirst(c);
		end = Tmat.getVectorLast(c);
		for (CoinBigIndex j = start; j < end; j++) {
			elts[nnz] = Telts[j];
			rowIdx[nnz++] = Tidx[j]+ncons1;
		}
	
	}

	// now W blocks
	int rowOffset = ncons1;
	int colOffset = nvar1;
	int const *Widx = Wmat.getIndices();
	double const *Welts = Wmat.getElements();

	for (int c = 0; c < nvar2; c++) {
		starts[colOffset++] = nnz;
		start = Wmat.getVectorFirst(c);
		end = Wmat.getVectorLast(c);
		for (CoinBigIndex j = start; j < end; j++) {
			elts[nnz] = Welts[j];
			rowIdx[nnz++] = Widx[j]+rowOffset;
		}
	}
	starts[totalVar] = nnz;
	assert(nnz == totalNnz);

	CoinPackedMatrix *constr = new CoinPackedMatrix();
	int *lens = 0;
	constr->assignMatrix(true,totalCons,totalVar,totalNnz,elts,rowIdx,starts,lens);
	constr->verifyMtx(); // debugging

	// OsiClpSolverInterface takes ownership of these
	double *collb = new double[totalVar];
	double *colub = new double[totalVar];
	double *obj = new double[totalVar];
	double *rowlb = new double[totalCons];
	double *rowub = new double[totalCons];

	concatenateOne(input, collb, &stochasticInput::getFirstStageColLB,
			&stochasticInput::getSecondStageColLB, scenarioNumber);
	
	concatenateOne(input, colub, &stochasticInput::getFirstStageColUB,
			&stochasticInput::getSecondStageColUB, scenarioNumber);

	concatenateOne(input, obj, &stochasticInput::getFirstStageObj,
			&stochasticInput::getSecondStageObj, scenarioNumber);
	
	concatenateOne(input, rowlb, &stochasticInput::getFirstStageRowLB,
			&stochasticInput::getSecondStageRowLB, scenarioNumber);

	concatenateOne(input, rowub, &stochasticInput::getFirstStageRowUB,
			&stochasticInput::getSecondStageRowUB, scenarioNumber);

	// rescale first-stage objective
	double scale = input.scenarioProbability(scenarioNumber);
	assert(lagrangeDiff.size() == static_cast<unsigned>(nvar1));
	for (int i = 0; i < nvar1; i++) {
		obj[i] = scale*obj[i] + lagrangeDiff[i];
	}

	m.assignProblem(constr, collb, colub, obj, rowlb, rowub);

	for (int i = 0; i < nvar1; i++) {
		if (input.isFirstStageColInteger(i)) m.setInteger(i);
	}
	for (int i = 0; i < nvar2; i++) {
		if (input.isSecondStageColInteger(scenarioNumber,i)) m.setInteger(i+nvar1);
	}

	cbcm.reset(new CbcModel(m));
	CbcMain0(*cbcm);
	//m.messageHandler()->setLogLevel(0);
	cbcm->messageHandler()->setLogLevel(0);

	//assert(CbcModel::haveMultiThreadSupport());
	//cbcm->setNumberThreads(4); // this should be adjusted to the system	


}


void CbcLagrangeSolver::go() {
	// TODO: are all the settings correct?
	// preprocessing, cuts, heuristics all enabled?

	const char * argv2[]={"","-solve","-log","0","-quit"};
	if (ratio != 0.) cbcm->setDblParam(CbcModel::CbcAllowableFractionGap,ratio);
	//cbcm->setMaximumNodes(10);
	CbcMain1(5,argv2,*cbcm);
    	//cbcm->branchAndBound();
   
}

double CbcLagrangeSolver::getBestPossibleObjective() const {
	return cbcm->getBestPossibleObjValue();
}

double CbcLagrangeSolver::getBestFeasibleObjective() const {
	return cbcm->getObjValue();
}

solverState CbcLagrangeSolver::getStatus() const {
	
	if (cbcm->isProvenOptimal()) {
		return Optimal;
	} 
	if (cbcm->isProvenInfeasible()) {
		return ProvenInfeasible;
	}
	assert(!cbcm->isAbandoned());

	return Initialized;

}


vector<double> CbcLagrangeSolver::getBestFirstStageSolution() const {
	//const double *s = cbcm->solver()->getColSolution();
	const double *s = cbcm->bestSolution();
	return vector<double>(s,s+nvar1);
}
