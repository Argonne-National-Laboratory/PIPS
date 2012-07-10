#ifndef CBCLAGRANGE_HPP
#define CBCLAGRANGE_HPP

#include "LagrangeSubproblemInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

class CbcLagrangeSolver : public LagrangeSubproblemInterface {
public:
	// scenario number is scenario which will be included in the subproblem.
	// if want to have multiple scenarios, combine them at the stochasticInput level
	// lagrangeMults is \lambda^i - \lambda^{i-1}. will be added to first-stage objective vector
	CbcLagrangeSolver(stochasticInput &input, int scenarioNumber, const std::vector<double>& lagrangeMults);

	void go();
	double getBestPossibleObjective() const;
	double getBestFeasibleObjective() const;
	solverState getStatus() const;

	// for setting states to warm-start the root node 
	void setFirstStageColRootState(int idx,variableState);
	void setFirstStageRowRootState(int idx,variableState);
	void setSecondStageColRootState(int idx,variableState);
	void setSecondStageRowRootState(int idx,variableState);
	void commitStates(); 

	variableState getFirstStageColRootState(int idx) const;
	variableState getFirstStageRowRootState(int idx) const;
	variableState getSecondStageColRootState(int idx) const;
	variableState getSecondStageRowRootState(int idx) const;

	// optimality gap ratio to terminate branch and bound, zero to disable
	void setRatio(double rel) { ratio = rel; }
	// absolute optimality gap for termination, zero to disable
	void setAbsoluteGap(double g) { absgap = g; }


	std::vector<double> getBestFirstStageSolution() const;
	std::vector<PrimalSolution> getBestFirstStageSolutions(double relcutoff) const { assert(0 && "not implemented"); return std::vector<PrimalSolution>(); } 

	struct WarmStart{
		WarmStart(std::vector<double> const& bestSol, CoinWarmStart* basis) :
			bestSol(bestSol), basis(basis) {}
		std::vector<double> bestSol; boost::shared_ptr<CoinWarmStart> basis;
	};

	// user must free!
	WarmStart* getWarmStart() const;
	void setWarmStart(const WarmStart&);

protected:
	OsiClpSolverInterface m;
	boost::scoped_ptr<CbcModel> cbcm;
	int nvar1;
	double ratio;
	double absgap;

};


#endif
