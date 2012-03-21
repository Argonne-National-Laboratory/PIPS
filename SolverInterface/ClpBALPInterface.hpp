#ifndef CLPBALPINTERFACE_HPP
#define CLPBALPINTERFACE_HPP


#include "stochasticInput.hpp"
#include "BALPSolverInterface.hpp"
#include "ClpSimplex.hpp"
#include "BA.hpp"

class ClpBALPInterface : public BALPSolverInterface<ClpBALPInterface> {
public:

	ClpBALPInterface(stochasticInput &in, BAContext &,solveType t = useDual);
	//ClpBALPInterface(const BAData &d, solveType t);
	~ClpBALPInterface() {}

	void writeStatus(const std::string &filebase);
	void loadStatus(const std::string &filebase);

	void go();

	void setPrimalTolerance(double val) { model.setCurrentPrimalTolerance(val); }
	void setDualTolerance(double val) { model.setCurrentDualTolerance(val); }
	double getObjective() const { return model.objectiveValue(); }
	solverState getStatus() const;

	std::vector<double> getFirstStagePrimalColSolution() const;

	void setFirstStageColState(int idx,variableState); 
	void setFirstStageRowState(int idx,variableState);
	void setSecondStageColState(int scen, int idx,variableState);
	void setSecondStageRowState(int scen, int idx,variableState);
	void commitStates() {} 

	variableState getFirstStageColState(int idx) const;
	variableState getFirstStageRowState(int idx) const;
	variableState getSecondStageColState(int scen, int idx) const;
	variableState getSecondStageRowState(int scen, int idx) const;

	// add row that preserves block-angular structure
	// e.g. lb <= elts1*x + elts2*y_i <= ub
	void addRow(const std::vector<double>& elts1, const std::vector<double> &elts2, int scen, double lb = -COIN_DBL_MAX, double ub = COIN_DBL_MAX);

protected:
	ClpSimplex model;
	BADimensions dims;
	solveType t;

	const BADimensions& getDims() const { return dims; }

};


#endif
