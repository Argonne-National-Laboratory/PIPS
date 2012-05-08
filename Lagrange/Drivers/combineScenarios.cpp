#include "scenarioReduction.hpp"
#include "combineScenarios.hpp"
#include <algorithm>

// dummy wrapper so we can send successively smaller problems to scenario reduction

class stochasticInputSubsetWrapper : public stochasticInput {
public:
	stochasticInputSubsetWrapper(stochasticInput &inner) :
		inner(inner) {
		int nscen = inner.nScenarios();
		for (int i = 0; i < nscen; i++) {
			realScenarios.push_back(i);
		}
		rescale = 1.0;
	}
	virtual int nScenarios() { return realScenarios.size(); }
	virtual int nFirstStageVars() { return inner.nFirstStageVars(); }
	virtual int nFirstStageCons() { return inner.nFirstStageCons(); }
	virtual int nSecondStageVars(int scen) { return inner.nSecondStageVars(realScenarios[scen]); }
	virtual int nSecondStageCons(int scen) { return inner.nSecondStageCons(realScenarios[scen]); }

	virtual std::vector<double> getFirstStageColLB() { return inner.getFirstStageColLB(); }
	virtual std::vector<double> getFirstStageColUB() { return inner.getFirstStageColUB(); }
	virtual std::vector<double> getFirstStageObj() { return inner.getFirstStageObj(); }
	virtual std::vector<std::string> getFirstStageColNames() { return inner.getFirstStageColNames(); }
	virtual std::vector<double> getFirstStageRowLB() { return inner.getFirstStageRowLB(); }
	virtual std::vector<double> getFirstStageRowUB() { return inner.getFirstStageRowUB(); }
	virtual std::vector<std::string> getFirstStageRowNames() { return inner.getFirstStageRowNames(); }
	virtual bool isFirstStageColInteger(int col) { return inner.isFirstStageColInteger(col); }

	virtual std::vector<double> getSecondStageColLB(int scen) { return inner.getSecondStageColLB(realScenarios[scen]); }
	virtual std::vector<double> getSecondStageColUB(int scen) { return inner.getSecondStageColUB(realScenarios[scen]); }
	// should really rescale objective, but this isn't needed for scenario reduction
	virtual std::vector<double> getSecondStageObj(int scen) { return inner.getSecondStageObj(realScenarios[scen]); }
	virtual std::vector<std::string> getSecondStageColNames(int scen) { return inner.getSecondStageColNames(realScenarios[scen]); }
	virtual std::vector<double> getSecondStageRowUB(int scen) { return inner.getSecondStageRowUB(realScenarios[scen]); }
	virtual std::vector<double> getSecondStageRowLB(int scen) { return inner.getSecondStageRowLB(realScenarios[scen]); }
	virtual std::vector<std::string> getSecondStageRowNames(int scen) { return inner.getSecondStageRowNames(realScenarios[scen]); }
	virtual double scenarioProbability(int scen) { return inner.scenarioProbability(realScenarios[scen])/rescale; }
	virtual bool isSecondStageColInteger(int scen, int col) { return inner.isSecondStageColInteger(realScenarios[scen],col); }

	virtual CoinPackedMatrix getFirstStageConstraints() { return inner.getFirstStageConstraints(); }
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) { return inner.getSecondStageConstraints(realScenarios[scen]); }
	virtual CoinPackedMatrix getLinkingConstraints(int scen) { return inner.getLinkingConstraints(realScenarios[scen]); }

	void removeScenario(int idx) { rescale -= inner.scenarioProbability(idx);
		realScenarios.erase(std::find(realScenarios.begin(),realScenarios.end(),idx)); }
	int realScenarioIdx(int idx) { return realScenarios[idx]; }

	virtual bool scenarioDimensionsEqual() { return inner.scenarioDimensionsEqual(); }
	virtual bool onlyBoundsVary() { return inner.onlyBoundsVary(); }
	virtual bool allProbabilitiesEqual() { return inner.allProbabilitiesEqual(); }
	virtual bool continuousRecourse() { return inner.continuousRecourse(); }
private:
	stochasticInput &inner;
	std::vector<int> realScenarios;
	double rescale;
};

combinedInput combineScenarios(stochasticInput &input, 
	int nper, bool scenred) {

	int nscen = input.nScenarios(); 
	stochasticInputSubsetWrapper wrapper(input);

	std::vector<std::vector<int> > scenarioMap;
	
	int nsubproblems = 0;
	for (int i = 0; i < nscen; i += nper) {
		scenarioMap.resize(nsubproblems+1);
		int scenthis = std::min(nper,nscen-i);
		if (scenred && scenthis == nper) {
			std::vector<int> subproblemScen = fastForwardSelection(wrapper,scenthis);
			for (int k = 0; k < nper; k++) {
				scenarioMap[nsubproblems].push_back(wrapper.realScenarioIdx(subproblemScen[k]));
			}
		} else {
			for (int k = 0; k < scenthis; k++) {
				scenarioMap[nsubproblems].push_back(wrapper.realScenarioIdx(k));
			}
		}
		for (int k = 0; k < scenthis; k++) {
			wrapper.removeScenario(scenarioMap[nsubproblems][k]);
			std::cout << scenarioMap[nsubproblems][k] << " ";
		}
		std::cout << std::endl;
		nsubproblems++;
	}

	return combinedInput(input,scenarioMap);

}

