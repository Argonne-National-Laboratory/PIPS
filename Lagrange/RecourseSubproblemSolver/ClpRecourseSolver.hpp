#ifndef CLPRECOURSE_HPP
#define CLPRECOURSE_HPP

#include "RecourseSubproblemInterface.hpp"
#include "ClpSimplex.hpp"

class ClpRecourseSolver : public RecourseSubproblemInterface {
public:
	ClpRecourseSolver(stochasticInput &input, int scenarioNumber, const std::vector<double> &firstStageSolution);

	void go();
	double getObjective() const;
	solverState getStatus() const;

	// for setting states to warm-start the solve 
	void setSecondStageColState(int idx,variableState);
	void setSecondStageRowState(int idx,variableState);
	void commitStates() {} 

	variableState getSecondStageColState(int idx) const;
	variableState getSecondStageRowState(int idx) const;

	void setDualObjectiveLimit(double d);

protected:
	ClpSimplex solver;
};



#endif
