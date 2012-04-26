#include "combinedInput.hpp"

using namespace std;

combinedInput::combinedInput(stochasticInput &inner, std::vector<std::vector<int> > const& scenarioMap) :
	scenarioMap(scenarioMap), inner(inner) {

	assert(scenarioMap.size() > 0);
	equalScenarios = true;
	unsigned siz = scenarioMap[0].size();
	for (unsigned i = 1; i < scenarioMap.size(); i++) {
		if (scenarioMap[i].size() != siz) {
			equalScenarios = false;
			break;
		}
	}

}


int combinedInput::nSecondStageVars(int scen) {
	int nvar = 0;
	for (unsigned i = 0; i < scenarioMap[scen].size(); i++) {
		nvar += inner.nSecondStageVars(scenarioMap[scen][i]);
	}
	return nvar;
}

int combinedInput::nSecondStageCons(int scen) {
	int ncons = 0;
	for (unsigned i = 0; i < scenarioMap[scen].size(); i++) {
		ncons += inner.nSecondStageCons(scenarioMap[scen][i]);
	}
	return ncons;
}

namespace{
template <typename T> vector<T> concatenateSubset(stochasticInput &data,  
			vector<T>(stochasticInput::*second)(int),
			int (stochasticInput::*secondDims)(int),
			vector<int> const &scens) {
	
	vector<T> out;

	for (unsigned i = 0; i < scens.size(); i++) {
		int scen = scens[i];
		int len2 = (data.*secondDims)(scen);
		vector<T> const &secondArray = (data.*second)(scen);
		assert(len2 == (int)secondArray.size());
		for (int k = 0; k < len2; k++) {
			out.push_back(secondArray[k]);
		}
	}
	return out;

}
}



vector<double> combinedInput::getSecondStageColLB(int scen) {
	return concatenateSubset(inner,&stochasticInput::getSecondStageColLB,
		&stochasticInput::nSecondStageVars,scenarioMap[scen]);
}

vector<double> combinedInput::getSecondStageColUB(int scen) {
	return concatenateSubset(inner,&stochasticInput::getSecondStageColUB,
		&stochasticInput::nSecondStageVars,scenarioMap[scen]);
}

vector<double> combinedInput::getSecondStageObj(int scen) {
	return concatenateSubset(inner,&stochasticInput::getSecondStageObj,
		&stochasticInput::nSecondStageVars,scenarioMap[scen]);
}

vector<string> combinedInput::getSecondStageColNames(int scen) {
	return concatenateSubset(inner,&stochasticInput::getSecondStageColNames,
		&stochasticInput::nSecondStageVars,scenarioMap[scen]);
}

vector<double> combinedInput::getSecondStageRowLB(int scen) {
	return concatenateSubset(inner,&stochasticInput::getSecondStageRowLB,
		&stochasticInput::nSecondStageCons,scenarioMap[scen]);
}

vector<double> combinedInput::getSecondStageRowUB(int scen) {
	return concatenateSubset(inner,&stochasticInput::getSecondStageRowUB,
		&stochasticInput::nSecondStageCons,scenarioMap[scen]);
}

vector<string> combinedInput::getSecondStageRowNames(int scen) {
	return concatenateSubset(inner,&stochasticInput::getSecondStageRowNames,
		&stochasticInput::nSecondStageCons,scenarioMap[scen]);
}


double combinedInput::scenarioProbability(int scen) {
	double p = 0.;
	for (unsigned i = 0; i < scenarioMap[scen].size(); i++) 
		p += inner.scenarioProbability(scenarioMap[scen][i]);
	return p;
}

bool combinedInput::isSecondStageColInteger(int scen, int col) {
	int i = 0, colMax = 0, nvarThis;
	while (col >= colMax + (nvarThis=inner.nSecondStageVars(scenarioMap[scen][i])) ) {
		colMax += nvarThis; i++;
	}
	return inner.isSecondStageColInteger(i,col-colMax);
}


CoinPackedMatrix combinedInput::getSecondStageConstraints(int scen) {
	int const totalVar = nSecondStageVars(scen);
	int const totalCons = nSecondStageCons(scen);
	int const nScenarios = scenarioMap[scen].size();
	
	CoinBigIndex totalNnz = 0;
	for (int k = 0; k < nScenarios; k++) {
		int s = scenarioMap[scen][k];
		totalNnz += inner.getSecondStageConstraints(s).getNumElements();
	}

	// CoinPackedMatrix takes ownership of these, so we don't free them
	CoinBigIndex *starts = new CoinBigIndex[totalVar+1];
	double *elts = new double[totalNnz];
	int *rowIdx = new int[totalNnz];
	
	CoinBigIndex nnz = 0, start, end;
	
	int rowOffset = 0;
	int colOffset = 0;
	for (int k = 0; k < nScenarios; k++) {
		int s = scenarioMap[scen][k];
		int nSecondStageVars = inner.nSecondStageVars(s);
		CoinPackedMatrix const &Wmat = inner.getSecondStageConstraints(s);
		int const *Widx = Wmat.getIndices();
		double const *Welts = Wmat.getElements();

		for (int c = 0; c < nSecondStageVars; c++) {
			starts[colOffset++] = nnz;
			start = Wmat.getVectorFirst(c);
			end = Wmat.getVectorLast(c);
			for (CoinBigIndex j = start; j < end; j++) {
				elts[nnz] = Welts[j];
				rowIdx[nnz++] = Widx[j]+rowOffset;
			}
		}
		rowOffset += inner.nSecondStageCons(s);
	}
	starts[totalVar] = nnz;
	assert(nnz == totalNnz);

	CoinPackedMatrix constr;
	int *lens = 0;
	constr.assignMatrix(true,totalCons,totalVar,totalNnz,
		elts, rowIdx, starts, lens);

	return constr;
}

CoinPackedMatrix combinedInput::getLinkingConstraints(int scen) {
	int const nScenarios = scenarioMap[scen].size();
	
	
	CoinPackedMatrix mat = inner.getLinkingConstraints(scenarioMap[scen][0]);

	for (int k = 1; k < nScenarios; k++) {
		mat.bottomAppendPackedMatrix(inner.getLinkingConstraints(scenarioMap[scen][k]));
	}

	return mat;
}
