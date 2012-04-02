#include "OsiSubproblemWrapper.hpp"
#include <cmath>

using namespace std;

namespace {

vector<double> concat(const vector<double>&a, const vector<double>&b) {
	vector<double> out(a);
	out.insert(out.end(),b.begin(),b.end());
	return out;
}

}

OsiSubproblemWrapper::OsiSubproblemWrapper(stochasticInput& in, int scen) {

	nvar1 = in.nFirstStageVars();
	nvar2 = in.nSecondStageVars(scen);
	ncons2 = in.nSecondStageCons(scen);

	// first-stage variables first
	cMat = in.getLinkingConstraints(scen);
	cMat.rightAppendPackedMatrix(in.getSecondStageConstraints(scen));
	rMat.reverseOrderedCopyOf(cMat);

	collb = concat(in.getFirstStageColLB(),in.getSecondStageColLB(scen));
	colub = concat(in.getFirstStageColUB(),in.getSecondStageColUB(scen));
	rowlb = concat(in.getFirstStageRowLB(),in.getSecondStageRowLB(scen));
	rowub = concat(in.getFirstStageRowUB(),in.getSecondStageRowUB(scen));

	for (int i = 0; i < nvar1; i++) {
		isInteger.push_back(in.isFirstStageColInteger(i));
	}
	for (int i = 0; i < nvar2; i++) {
		isInteger.push_back(in.isSecondStageColInteger(scen,i));
	}

	rhs.resize(nvar1+nvar2);

	for (int i = 0; i < nvar1 + nvar2; i++) {
		if (rowlb[i] < -1e-20) {
			rowlb[i] = -COIN_DBL_MAX;
			if (rowub[i] > 1e20) {
				rowub[i] = COIN_DBL_MAX;
				rowSense.push_back('N');
				rhs[i] = 0.0;
			} else {
				rowSense.push_back('L');
				rhs[i] = rowub[i];
			}
		} else {
			if (rowub[i] > 1e20){
				rowub[i] = COIN_DBL_MAX;
				rowSense.push_back('G');
				rhs[i] = rowlb[i];
			} else if (fabs(rowlb[i] - rowub[i]) < 1e-7) {
				rowSense.push_back('E');
				rhs[i] = rowlb[i];
			} else {
				rowSense.push_back('R');
				rhs[i] = rowub[i];
			}
		}
	}
				

}



void OsiSubproblemWrapper::setCurrentSolution(const std::vector<double> stg1, const std::vector<double> stg2) {
	colsol = concat(stg1, stg2);
	rowActivity.resize(nvar1+nvar2);
	cMat.times(&colsol[0], &rowActivity[0]);
}

void OsiSubproblemWrapper::setCurrentReducedCosts(const std::vector<double> stg1, const std::vector<double> stg2) {
	coldualsol = concat(stg1, stg2);
}
