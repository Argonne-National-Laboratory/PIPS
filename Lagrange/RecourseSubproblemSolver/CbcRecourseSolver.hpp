#ifndef CBCRECOURSE_HPP
#define CBCRECOURSE_HPP

#include "RecourseSubproblemInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include <boost/scoped_ptr.hpp>

class CbcRecourseSolver : public RecourseSubproblemInterface {
public:
	CbcRecourseSolver(stochasticInput &input, int scenarioNumber, const std::vector<double> &firstStageSolution);

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
	OsiClpSolverInterface model;
	boost::scoped_ptr<CbcModel> cbcm;
};



#endif
