#ifndef BALPSOLVERINTERFACE_HPP
#define BALPSOLVERINTERFACE_HPP

enum constraintType { Free, LB, UB, Range, Fixed };

enum variableState { Basic, AtLower, AtUpper };

// unbounded/infeasible are in the sense of the primal problem 
enum solverState { Uninitialized, LoadedFromFile, Initialized, PrimalFeasible, DualFeasible, Optimal, ProvenUnbounded, ProvenInfeasible };

// abstract interface to a BALP (block-angular linear program) solver
class BALPSolverInterface {
public:
	virtual ~BALPSolverInterface() {}
	
	virtual void go() = 0;

	virtual void setPrimalTolerance(double val) = 0;
	virtual void setDualTolerance(double val) = 0;
	virtual double getObjective() const = 0;
	virtual solverState getStatus() const = 0;

	virtual std::vector<double> getFirstStagePrimalColSolution() const = 0;

	virtual void setFirstStageColState(int idx,variableState) = 0;
	virtual void setFirstStageRowState(int idx,variableState) = 0;
	virtual void setSecondStageColState(int scen, int idx,variableState) = 0;
	virtual void setSecondStageRowState(int scen, int idx,variableState) = 0;
	virtual void commitStates() = 0; // must call this after setting states

	virtual variableState getFirstStageColState(int idx) const = 0;
	virtual variableState getFirstStageRowState(int idx) const = 0;
	virtual variableState getSecondStageColState(int scen, int idx) const = 0;
	virtual variableState getSecondStageRowState(int scen, int idx) const = 0;




};


#endif

