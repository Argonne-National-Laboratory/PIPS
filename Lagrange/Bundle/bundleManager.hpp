#ifndef BUNDLEMANAGER_HPP
#define BUNDLEMANAGER_HPP

#include "BA.hpp"
#include "BALPSolverInterface.hpp"
#include "CoinFinite.hpp"


struct cutInfo {
	std::vector<double> evaluatedAt, subgradient, primalSol;
	double objval;
	double computeC() const; // compute component of c as defined above
	// could include error estimates here later
};


// outer dimension is scenario, inner is list of cuts.
// this allows bundle to be distributed per scenario, i.e,
// if a scenario s isn't assigned, bundle[s] will be empty. 
typedef std::vector<std::vector<cutInfo> > bundle_t; 


// this is specialized particularly for lagrangian relaxation of nonanticipativity constraints, not general bundle solver
template<typename BALPSolver, typename LagrangeSolver, typename RecourseSolver> class bundleManager {
public:
	bundleManager(stochasticInput &input, BAContext & ctx) : ctx(ctx), input(input) {
		int nscen = input.nScenarios();
		// use zero as initial iterate
		currentSolution.resize(nscen,std::vector<double>(input.nFirstStageVars(),0.));
		nIter = -1;
		bundle.resize(nscen);
		bestPrimalObj = COIN_DBL_MAX;
		relativeConvergenceTol = 1e-6;
		terminated_ = false;
	}


	void setRelConvergenceTol(double t) { relativeConvergenceTol = t; }
	bool terminated() { return terminated_; }
	void iterate();

protected:

	virtual void doStep() = 0;
	cutInfo solveSubproblem(std::vector<double> const& at, int scen); 
	double testPrimal(std::vector<double> const& primal); // test a primal solution and return objective (COIN_DBL_MAX) if infeasible
	// evaluates trial solution and updates the bundle
	double evaluateSolution(std::vector<std::vector<double> > const& sol);

	void evaluateAndUpdate();
	void checkLastPrimals();
	
	bundle_t bundle;
	std::vector<std::vector<double> > currentSolution; // dual solution
	double currentObj;
	std::vector<double> bestPrimal;
	double bestPrimalObj;
	BAContext &ctx;
	stochasticInput &input;
	int nIter;
	double relativeConvergenceTol;
	bool terminated_;
	


};


template<typename B,typename L, typename R> cutInfo bundleManager<B,L,R>::solveSubproblem(std::vector<double> const& at, int scen) {
	using namespace std;	
	int nvar1 = input.nFirstStageVars();
	//double t = MPI_Wtime();
	

	L lsol(input,scen,at);
	//lsol.setRatio(100*fabs(relprec));
	lsol.go();
	//cout << "SCEN " << scen << " DONE, " << MPI_Wtime() - t << " SEC" << endl;

	//assert(lsol.getStatus() == Optimal);
	assert(lsol.getStatus() != ProvenInfeasible);
	
	cutInfo cut;
	cut.objval = -lsol.getBestPossibleObjective();
	cut.primalSol = lsol.getBestFirstStageSolution(); 
	cut.evaluatedAt = at;
	cut.subgradient.resize(nvar1);
	for (int i = 0; i < nvar1; i++) {
		cut.subgradient[i] = -cut.primalSol[i];
	}

	return cut;

}


// this won't work with distributed scenarios
template<typename B, typename L, typename R> double bundleManager<B,L,R>::testPrimal(std::vector<double> const& primal) {

	const std::vector<double> &obj1 = input.getFirstStageObj();
	int nscen = input.nScenarios();
	int nvar1 = input.nFirstStageVars();
	double obj = 0.;
	bool infeas = false;
	for (int s = 0; s < nscen; s++) {
		R rsol(input,s,primal);
		rsol.go();
		obj += rsol.getObjective();
		
		if (rsol.getStatus() == ProvenInfeasible) {
			printf("got infeasible 1st stage\n");
			infeas = true; break;
		}
	}
	if (infeas) return COIN_DBL_MAX;
	for (int k = 0; k < nvar1; k++) {
		obj += primal[k]*obj1[k];
	}
	
	return obj;

}

template<typename B, typename L, typename R> double bundleManager<B,L,R>::evaluateSolution(std::vector<std::vector<double> > const& sol) {

	int nscen = input.nScenarios();
	double obj = 0.;
	for (int i = 0; i < nscen; i++) {
		cutInfo cut = solveSubproblem(sol[i],i);
		bundle[i].push_back(cut);
		obj -= cut.objval; // note cut has obj flipped
	}
	return obj;
}

template<typename B, typename L, typename R> void bundleManager<B,L,R>::evaluateAndUpdate() {

	currentObj = evaluateSolution(currentSolution);
}

template<typename B, typename L, typename R> void bundleManager<B,L,R>::checkLastPrimals() {

	int nscen = input.nScenarios();
	for (int i = 0; i < 5; i++) {
	//for (int i = 0; i < nscen; i++) {
		assert(bundle[i].size());
		std::vector<double> const& p = bundle[i][bundle[i].size()-1].primalSol;
		double o = testPrimal(p);
		if (o < bestPrimalObj) {
			bestPrimalObj = o;
			bestPrimal = p;
		}
	}
}


template<typename B, typename L, typename R> void bundleManager<B,L,R>::iterate() {

	assert(!terminated_);
	

	if (nIter++ == -1) { // just evaluate and generate subgradients
		evaluateAndUpdate();
		checkLastPrimals();
		return;
	}

	doStep();
	
}

#endif
