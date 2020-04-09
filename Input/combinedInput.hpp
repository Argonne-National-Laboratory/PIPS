#ifndef COMBINEDINPUTHPP
#define COMBINEDINPUTHPP

#include "stochasticInput.hpp"
#include <cmath>

// wrapper for combining scenarios for lagrangian subproblems

class combinedInput : public stochasticInput {
public:
	combinedInput(stochasticInput &inner, std::vector<std::vector<int> > const& scenarioMap);
	virtual int nScenarios() { return scenarioMap.size(); }
	virtual int nFirstStageVars() { return inner.nFirstStageVars(); }
	virtual int nFirstStageCons() { return inner.nFirstStageCons(); }
	virtual int nSecondStageVars(int scen);
	virtual int nSecondStageCons(int scen);

	virtual std::vector<double> getFirstStageColLB() { return inner.getFirstStageColLB(); }
	virtual std::vector<double> getFirstStageColUB() { return inner.getFirstStageColUB(); }
	virtual std::vector<double> getFirstStageObj() { return inner.getFirstStageObj(); }
	virtual std::vector<std::string> getFirstStageColNames() { return inner.getFirstStageColNames(); }
	virtual std::vector<double> getFirstStageRowLB() { return inner.getFirstStageRowLB(); }
	virtual std::vector<double> getFirstStageRowUB() { return inner.getFirstStageRowUB(); }
	virtual std::vector<std::string> getFirstStageRowNames() { return inner.getFirstStageRowNames(); }
	virtual bool isFirstStageColInteger(int col) { return inner.isFirstStageColInteger(col); }
	virtual bool isFirstStageColBinary(int col) {
		bool isInteger = this->isFirstStageColInteger(col);
		// CoinMpsIO has no isBinary member function, but some preprocessing features require
		// knowledge of binary variables, so kludge in an "isBinary" member function by
		// relying on CoinMpsIO setting lower and upper bounds to zero and one, respectively.
		// Also note: CoinMpsIO uses a default tolerance of 1.0e-8 on integrality comparisons.
		const double intTol = 1.0e-8;
		bool isLBzero = (fabs(this->getFirstStageColLB().at(col)) < intTol);
		bool isUBone = (fabs(this->getFirstStageColUB().at(col) - 1.0) < intTol);
		return (isInteger && isLBzero && isUBone);
	}

	virtual std::vector<double> getSecondStageColLB(int scen);
	virtual std::vector<double> getSecondStageColUB(int scen);
	virtual std::vector<double> getSecondStageObj(int scen);
	virtual std::vector<std::string> getSecondStageColNames(int scen);
	virtual std::vector<double> getSecondStageRowUB(int scen);
	virtual std::vector<double> getSecondStageRowLB(int scen);
	virtual std::vector<std::string> getSecondStageRowNames(int scen);
	virtual double scenarioProbability(int scen);
	virtual bool isSecondStageColInteger(int scen, int col);
	virtual bool isSecondStageColBinary(int scen, int col) {
		bool isInteger = this->isSecondStageColInteger(scen, col);
		// CoinMpsIO has no isBinary member function, but some preprocessing features require
		// knowledge of binary variables, so kludge in an "isBinary" member function by
		// relying on CoinMpsIO setting lower and upper bounds to zero and one, respectively.
		// Also note: CoinMpsIO uses a default tolerance of 1.0e-8 on integrality comparisons.
		const double intTol = 1.0e-8;
		bool isLBzero = (fabs(this->getSecondStageColLB(scen).at(col)) < intTol);
		bool isUBone = (fabs(this->getSecondStageColUB(scen).at(col) - 1.0) < intTol);
		return (isInteger && isLBzero && isUBone);
	}

	virtual CoinPackedMatrix getFirstStageConstraints() { return inner.getFirstStageConstraints(); }
	virtual CoinPackedMatrix getSecondStageConstraints(int scen);
	virtual CoinPackedMatrix getLinkingConstraints(int scen);

	

	
	virtual bool scenarioDimensionsEqual() { return inner.scenarioDimensionsEqual() && equalScenarios; }
	virtual bool onlyBoundsVary() { return inner.onlyBoundsVary(); }
	virtual bool allProbabilitiesEqual() { return equalScenarios && inner.allProbabilitiesEqual(); }
	virtual bool continuousRecourse() { return inner.continuousRecourse(); }

private:
	// map from "fake" scenario index to group of scenarios it represents
	std::vector<std::vector<int> > scenarioMap;
	bool equalScenarios; // equal number of scenarios in each combined scenario
	stochasticInput &inner;

};



#endif
