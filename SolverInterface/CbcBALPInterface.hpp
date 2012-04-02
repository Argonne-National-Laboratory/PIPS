#ifndef CBCBALPINTERFACE_HPP
#define CBCBALPINTERFACE_HPP


#include "stochasticInput.hpp"
#include "BALPSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "BA.hpp"

#include <boost/scoped_ptr.hpp>

class CbcBALPInterface : public BALPSolverInterface<CbcBALPInterface> {
public:

	CbcBALPInterface(stochasticInput &in, BAContext &);
	//ClpBALPInterface(const BAData &d, solveType t);
	~CbcBALPInterface() {}

	void writeStatus(const std::string &filebase);
	void loadStatus(const std::string &filebase);

	void go();

	void setPrimalTolerance(double val) { assert(0); }
	void setDualTolerance(double val) { assert(0); }
	void setDualObjectiveLimit(double val) { assert(0); }
	double getObjective() const { return model.getObjValue(); }
	solverState getStatus() const;

	std::vector<double> getFirstStagePrimalColSolution() const;
	std::vector<double> getSecondStagePrimalColSolution(int scen) const;
	std::vector<double> getFirstStageDualColSolution() const;
	std::vector<double> getSecondStageDualColSolution(int scen) const;

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

	void setFirstStageColLB(int idx, double newLb);
	void setFirstStageColUB(int idx, double newUb);

protected:
	OsiClpSolverInterface model;
	
	BADimensions dims;
	boost::scoped_ptr<CbcModel> cbcm;

	const BADimensions& getDims() const { return dims; }

};


#endif
