#include "stochasticInput.hpp"



CoinPackedMatrix stochasticInput::getFirstStageHessian() {
	int nvar1 = nFirstStageVars();
	std::vector<CoinBigIndex> starts(nvar1+1,0.);
	return CoinPackedMatrix(true,nvar1,nvar1,0,0,0,&starts[0],0);
}

CoinPackedMatrix stochasticInput::getSecondStageHessian(int scen) {
	int nvar2 = nSecondStageVars(scen);
	std::vector<CoinBigIndex> starts(nvar2+1,0.);
	return CoinPackedMatrix(true,nvar2,nvar2,0,0,0,&starts[0],0);
}

CoinPackedMatrix stochasticInput::getSecondStageCrossHessian(int scen) {
	int nvar1 = nFirstStageVars();
	int nvar2 = nSecondStageVars(scen);
	std::vector<CoinBigIndex> starts(nvar1+1,0.);
	return CoinPackedMatrix(true,nvar2,nvar1,0,0,0,&starts[0],0);
}
