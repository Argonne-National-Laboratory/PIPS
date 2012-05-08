#ifndef RECOURSESUBINTERFACE_HPP
#define RECOURSESUBINTERFACE_HPP

#include "BALPSolverInterface.hpp"

// this is basically a "concept" interface.
// implementing classes will be used as template parameters, so functions not virtual.
// no real need to actually derive from this class
class RecourseSubproblemInterface {

public:
	
	void go();
	double getObjective() const;
	solverState getStatus() const;

	// for setting states to warm-start the solve 
	void setSecondStageColState(int idx,variableState);
	void setSecondStageRowState(int idx,variableState);
	void commitStates(); 

	variableState getSecondStageColState(int idx) const;
	variableState getSecondStageRowState(int idx) const;

	// objective limit above which problem is considered infeasible when
	// using dual simplex
	void setDualObjectiveLimit(double);

	
};


#endif
