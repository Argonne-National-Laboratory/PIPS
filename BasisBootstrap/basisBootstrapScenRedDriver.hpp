#ifndef BASISBOOTSTRAPSCENRED_HPP
#define BASISBOOTSTRAPSCENRED_HPP

#include "basisBootstrapDriver.hpp"
#include "ScenarioReduction/scenarioReduction.hpp"

// select the first group of scenarios by scenario reduction
template<typename BALPSolver, typename RecourseSolver>
class basisBootstrapScenRedDriver : public basisBootstrapDriver<BALPSolver,RecourseSolver> {
public:
	basisBootstrapScenRedDriver(stochasticInput &input, BAContext& ctx) : basisBootstrapDriver<BALPSolver,RecourseSolver>(input,ctx) {}

	void goFromScratch(int startingScenarios, int addedPer) {
		// get starting scenarios from scenario reduction
		std::vector<int> s = fastForwardSelection(this->input,startingScenarios);

		this->initializeMaster(s);
		this->solveMaster();
		this->countSquare();
		
		while(this->scenariosNotInMaster.size()) {
			
			for (int i = 0; i < addedPer && this->scenariosNotInMaster.size(); i++) {
				this->addScenarioToMaster(this->scenariosNotInMaster[0]);
			}
			this->solveMaster();
			this->countSquare();

		}

	}


};

#endif
