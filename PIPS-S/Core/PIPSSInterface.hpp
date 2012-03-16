#ifndef PIPSSINTERFACE_HPP
#define PIPSSINTERFACE_HPP

#include "BALPSolverBase.hpp"
#include "stochasticInput.hpp"
#include "BALPSolverInterface.hpp"

class BALPSolverDual;

// public interface to the solver
// handles switching between primal/dual solvers to clean up infeasibilities
class PIPSSInterface : public BALPSolverInterface {
public:
	enum solveType { usePrimal, useDual };
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

	//const denseBAVector& getPrimalSolution() const { return solver->getPrimalSolution(); }
	std::vector<double> getFirstStagePrimalColSolution() const;

	virtual void setFirstStageColState(int idx,variableState s); 
	virtual void setFirstStageRowState(int idx,variableState);
	virtual void setSecondStageColState(int scen, int idx,variableState);
	virtual void setSecondStageRowState(int scen, int idx,variableState);
	virtual void commitStates(); 

	virtual variableState getFirstStageColState(int idx) const;
	virtual variableState getFirstStageRowState(int idx) const;
	virtual variableState getSecondStageColState(int scen, int idx) const;
	virtual variableState getSecondStageRowState(int scen, int idx) const;

	// will dump current status every d iterations. zero to disable.
	void setDumpFrequency(int d, const std::string &outputname) { solver->setDumpFrequency(d,outputname); }

protected:
	void setStates(const BAFlagVector<variableState> &s) { solver->setStates(s); }
	const BAFlagVector<variableState>& getStates() const { return solver->getStates(); }
	
	BALPSolverBase *solver;
	BAData d;

friend class BALPSolverDual;

};


#endif
