#ifndef SCIPLAGRANGE_HPP
#define SCIPLAGRANGE_HPP

#include "LagrangeSubproblemInterface.hpp"

struct Scip;
struct SCIP_Var;
struct SCIP_Cons;

class ScipLagrangeSolver : public LagrangeSubproblemInterface {
public:
	// scenario number is scenario which will be included in the subproblem.
	// if want to have multiple scenarios, combine them at the stochasticInput level
	// lagrangeMults will be added to first-stage objective vector
	ScipLagrangeSolver(stochasticInput &input, int scenarioNumber, const std::vector<double>& lagrangeMults);
	~ScipLagrangeSolver();

	void go();
	double getBestPossibleObjective() const;
	double getBestFeasibleObjective() const;
	solverState getStatus() const;


	// optimality gap ratio to terminate branch and bound, zero to disable
	void setRatio(double rel) { ratio = rel; }
	// absolute optimality gap for termination, zero to disable
	void setAbsoluteGap(double g) { absgap = g; }

	void setFirstStageColLB(int idx, double newLb);


	std::vector<double> getBestFirstStageSolution() const;
	/*
	struct WarmStart{
		WarmStart(std::vector<double> const& bestSol, CoinWarmStart* basis) :
			bestSol(bestSol), basis(basis) {}
		std::vector<double> bestSol; boost::shared_ptr<CoinWarmStart> basis;
	};*/

	// user must free!
	typedef int WarmStart; // dummy
	WarmStart* getWarmStart() const { return new int; }
	void setWarmStart(const WarmStart&) {}
	
protected:
	Scip *scip;
	std::vector<SCIP_Var*> vars1, vars2; // first and second-stage vars
	std::vector<SCIP_Cons*> cons;
	int nvar1;
	double ratio;
	double absgap;

private:
	ScipLagrangeSolver(const ScipLagrangeSolver&);

};


#endif
