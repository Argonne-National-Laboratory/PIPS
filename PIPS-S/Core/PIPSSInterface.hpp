#ifndef PIPSSINTERFACE_HPP
#define PIPSSINTERFACE_HPP

#include "BALPSolverBase.hpp"
#include "stochasticInput.hpp"
#include "BALPSolverInterface.hpp"

class BALPSolverDual;

// public interface to the solver
// handles switching between primal/dual solvers to clean up infeasibilities
class PIPSSInterface : public BALPSolverInterface<PIPSSInterface> {
public:
	PIPSSInterface(stochasticInput &in, BAContext &ctx, solveType t);
	PIPSSInterface(const BAData &d, solveType t);
	~PIPSSInterface();


	// if writeBasisOnly is false, we also write edge weights and perturbed objective/lb/ub's.
	// this is for resuming in the middle of a solve
	void writeStatus(const std::string &filebase, bool appendIterateNumber = false,
			 bool writeBasisOnly = true) { solver->writeStatus(filebase,appendIterateNumber,writeBasisOnly); }
	void loadStatus(const std::string &filebase) { solver->loadStatus(filebase); }

	void go();

	void setPrimalTolerance(double val) { solver->setPrimalTolerance(val); }
	void setDualTolerance(double val) { solver->setDualTolerance(val); }
	double getObjective() const { return solver->objval; }
	solverState getStatus() const { return solver->status; }

	std::vector<double> getFirstStagePrimalColSolution() const;
	std::vector<double> getSecondStagePrimalColSolution(int scen) const;
	std::vector<double> getFirstStageDualColSolution() const;
	std::vector<double> getSecondStageDualColSolution(int scen) const;

	void setFirstStageColState(int idx,variableState s); 
	void setFirstStageRowState(int idx,variableState);
	void setSecondStageColState(int scen, int idx,variableState);
	void setSecondStageRowState(int scen, int idx,variableState);
	void commitStates(); 

	variableState getFirstStageColState(int idx) const;
	variableState getFirstStageRowState(int idx) const;
	variableState getSecondStageColState(int scen, int idx) const;
	variableState getSecondStageRowState(int scen, int idx) const;
	
	// override default
	void setStates(const BAFlagVector<variableState> &s) { solver->setStates(s); }

	// will dump current status every d iterations. zero to disable.
	void setDumpFrequency(int d, const std::string &outputname) { solver->setDumpFrequency(d,outputname); }

	void setFirstStageColLB(int idx, double newLb);
	void setFirstStageColUB(int idx, double newUb);

protected:
	const BAFlagVector<variableState>& getStatesRef() const { return solver->getStates(); }
	const BADimensions& getDims() const { return d.dims.inner; }
	
	BALPSolverBase *solver;
	BAData d;
	bool boundsChanged;

friend class BALPSolverDual;

};


#endif
