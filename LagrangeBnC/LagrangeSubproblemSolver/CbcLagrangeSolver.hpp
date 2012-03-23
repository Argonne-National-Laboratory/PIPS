#ifndef CBCLAGRANGE_HPP
#define CBCLAGRANGE_HPP

#include "LagrangeSubproblemInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include <boost/scoped_ptr.hpp>

class CbcLagrangeSolver : public LagrangeSubproblemInterface {
public:
	// scenario number is scenario which will be included in the subproblem.
	// if want to have multiple scenarios, combine them at the stochasticInput level
	// lagrangeMults is \lambda^i - \lambda^{i-1}. will be added to first-stage objective vector
	CbcLagrangeSolver(stochasticInput &input, int scenarioNumber, const std::vector<double>& lagrangeMults);

	void go();
	double getBestPossibleObjective() const;
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


	std::vector<double> getBestFirstStageSolution() const;


protected:
	OsiClpSolverInterface m;
	boost::scoped_ptr<CbcModel> cbcm;
	int nvar1;

};


#endif
