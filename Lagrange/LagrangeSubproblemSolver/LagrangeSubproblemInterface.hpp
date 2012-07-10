#ifndef LAGRANGESUBINTERFACE_HPP
#define LAGRANGESUBINTERFACE_HPP

#include "BALPSolverInterface.hpp"


// this is basically a "concept" interface.
// implementing classes will be used as template parameters, so functions not virtual.
// no real need to actually derive from this class
class LagrangeSubproblemInterface {

public:
	
	void go();
	double getObjective() const;
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

	
	std::vector<double> getFirstStagePrimalColSolution() const;

};


struct PrimalSolution {
	std::vector<double> sol;
	double objval;
};

#endif
