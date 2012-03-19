#ifndef BALPSOLVERINTERFACE_HPP
#define BALPSOLVERINTERFACE_HPP
#include "BAVector.hpp"

enum constraintType { Free, LB, UB, Range, Fixed };

enum variableState { Basic, AtLower, AtUpper };

// unbounded/infeasible are in the sense of the primal problem 
enum solverState { Uninitialized, LoadedFromFile, Initialized, PrimalFeasible, DualFeasible, Optimal, ProvenUnbounded, ProvenInfeasible, Stopped };

// abstract interface to a BALP (block-angular linear program) solver
// We use CRTP (curiously recurring template pattern)
template <typename Derived> class BALPSolverInterface {
public:
	enum solveType { usePrimal, useDual };
	~BALPSolverInterface() {}
	
	// derived classes should implement the following

	void go();

	void setPrimalTolerance(double val);
	void setDualTolerance(double val);
	double getObjective() const;
	solverState getStatus() const;

	std::vector<double> getFirstStagePrimalColSolution() const;

	void setFirstStageColState(int idx,variableState);
	void setFirstStageRowState(int idx,variableState);
	void setSecondStageColState(int scen, int idx,variableState);
	void setSecondStageRowState(int scen, int idx,variableState);
	void commitStates(); // must call this after setting states

	variableState getFirstStageColState(int idx) const;
	variableState getFirstStageRowState(int idx) const;
	variableState getSecondStageColState(int scen, int idx) const;
	variableState getSecondStageRowState(int scen, int idx) const;

	// shared functions
	// derived classes need not implement the following
	 	
	void setStates(const BAFlagVector<variableState> &s);
	void getStates(BAFlagVector<variableState> &s);

protected:
	// derived classes should implement the following if the shared functions will be used
	
	// don't we like C++... need to use virtual even inside CRTP
	virtual const BADimensions& getDims() const = 0;
	

};

template <typename Derived> void BALPSolverInterface<Derived>::setStates(const BAFlagVector<variableState> &s) {
	const BADimensions &dims = getDims();

	int nscen = dims.numScenarios();
	int nvar1 = dims.numFirstStageVars();
	int ncons1 = dims.numFirstStageCons();
	const denseFlagVector<variableState> &s1 = s.getFirstStageVec();
	for (int i = 0; i < nvar1; i++) {
		static_cast<Derived*>(this)->setFirstStageColState(i,s1[i]);
	}
	for (int i = 0; i < ncons1; i++) {
		static_cast<Derived*>(this)->setFirstStageRowState(i,s1[i+nvar1]);
	}
	for (int scen = 0; scen < nscen; scen++) {
		if (!s.hasScenario(scen)) continue;
		int nvar2 = dims.numSecondStageVars(scen);
		int ncons2 = dims.numSecondStageCons(scen);
		const denseFlagVector<variableState> &s2 = s.getSecondStageVec(scen);
		for (int i = 0; i < nvar2; i++) {
			static_cast<Derived*>(this)->setSecondStageColState(scen,i,s2[i]);
		}
		for (int i = 0; i < ncons2; i++) {
			static_cast<Derived*>(this)->setSecondStageRowState(scen,i,s2[i+nvar2]);
		}
	}

	static_cast<Derived*>(this)->commitStates();
}

template <typename Derived> void BALPSolverInterface<Derived>::getStates(BAFlagVector<variableState> &s) {
	const BADimensions &dims = getDims();

	int nscen = dims.numScenarios();
	int nvar1 = dims.numFirstStageVars();
	int ncons1 = dims.numFirstStageCons();
	denseFlagVector<variableState> &s1 = s.getFirstStageVec();
	for (int i = 0; i < nvar1; i++) {
		s1[i] = static_cast<Derived*>(this)->getFirstStageColState(i);
	}
	for (int i = 0; i < ncons1; i++) {
		printf("var: %d cons: %d acc: %d\n",nvar1, ncons1, i+nvar1);
		s1[i+nvar1] = static_cast<Derived*>(this)->getFirstStageRowState(i);
	}
	for (int scen = 0; scen < nscen; scen++) {
		if (!s.hasScenario(scen)) continue;
		int nvar2 = dims.numSecondStageVars(scen);
		int ncons2 = dims.numSecondStageCons(scen);
		denseFlagVector<variableState> &s2 = s.getSecondStageVec(scen);
		for (int i = 0; i < nvar2; i++) {
			s2[i] = static_cast<Derived*>(this)->getSecondStageColState(scen,i);
		}
		for (int i = 0; i < ncons2; i++) {
			s2[i+nvar2] = static_cast<Derived*>(this)->getSecondStageRowState(scen,i);
		}
	}

}

#endif

