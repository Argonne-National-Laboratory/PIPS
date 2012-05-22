#ifndef BUNDLEMANAGER_HPP
#define BUNDLEMANAGER_HPP

#include "BA.hpp"
#include "BALPSolverInterface.hpp"
#include "CoinFinite.hpp"
#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

struct cutInfo {
	std::vector<double> evaluatedAt, subgradient, primalSol;
	double objval; // underestimate of convex objective
	double objmax; // overestimate of convex objective
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
		hotstarts.resize(nscen);
		recourseRowStates.resize(nscen);
		recourseColStates.resize(nscen);
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
	std::vector<shared_ptr<typename LagrangeSolver::WarmStart> > hotstarts;
	std::vector<std::vector<variableState> > recourseRowStates, recourseColStates;
	


};


template<typename B,typename L, typename R> cutInfo bundleManager<B,L,R>::solveSubproblem(std::vector<double> const& at, int scen) {
	using namespace std;	
	int nvar1 = input.nFirstStageVars();
	//double t = MPI_Wtime();
	

	L lsol(input,scen,at);
	// hot start for root node LP
	if (hotstarts[scen]) {
		lsol.setWarmStart(hotstarts[scen].get());
	} else if (hotstarts[ctx.localScenarios().at(1)]) {
		lsol.setWarmStart(hotstarts[ctx.localScenarios()[1]].get());
	}
	//lsol.setRatio(100*fabs(relprec));
	lsol.go();
	hotstarts[scen].reset(lsol.getWarmStart());
	//cout << "SCEN " << scen << " DONE, " << MPI_Wtime() - t << " SEC" << endl;

	//assert(lsol.getStatus() == Optimal);
	assert(lsol.getStatus() != ProvenInfeasible);
	
	cutInfo cut;
	cut.objval = -lsol.getBestFeasibleObjective();
	cut.objmax = -lsol.getBestPossibleObjective();
	cut.primalSol = lsol.getBestFirstStageSolution(); 
	cut.evaluatedAt = at;
	cut.subgradient.resize(nvar1);
	for (int i = 0; i < nvar1; i++) {
		cut.subgradient[i] = -cut.primalSol[i];
	}

	return cut;

}


template<typename B, typename L, typename R> double bundleManager<B,L,R>::testPrimal(std::vector<double> const& primal) {

	const std::vector<double> &obj1 = input.getFirstStageObj();
	const std::vector<int> &localScen = ctx.localScenarios();
	int nvar1 = input.nFirstStageVars();
	double obj = 0.;
	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		int nvar2 = input.nSecondStageVars(scen);
		int ncons2 = input.nSecondStageCons(scen);
		R rsol(input,scen,primal);
		rsol.setDualObjectiveLimit(1e10);
		if (input.continuousRecourse() && recourseRowStates[scen].size()) {
			for (int k = 0; k < nvar2; k++) {
				rsol.setSecondStageColState(k,recourseColStates[scen][k]);
			}
			for (int k = 0; k < ncons2; k++) {
				rsol.setSecondStageRowState(k,recourseRowStates[scen][k]);
			}
		}
		
		rsol.go();
		if (input.continuousRecourse()) {
			recourseColStates[scen].resize(nvar2);
			for (int k = 0; k < nvar2; k++) {
				recourseColStates[scen][k] = rsol.getSecondStageColState(k);
			}
			recourseRowStates[scen].resize(ncons2);
			for (int k = 0; k < ncons2; k++) {
				recourseRowStates[scen][k] = rsol.getSecondStageRowState(k);
			}
		}

		obj += rsol.getObjective();
		
		if (rsol.getStatus() == ProvenInfeasible) {
			printf("got infeasible 1st stage\n");
			obj = COIN_DBL_MAX;
			break;
		}
	}
	double allobj;
	MPI_Allreduce(&obj,&allobj,1,MPI_DOUBLE,MPI_SUM,ctx.comm());
	for (int k = 0; k < nvar1; k++) {
		allobj += primal[k]*obj1[k];
	}
	
	return allobj;

}

template<typename B, typename L, typename R> double bundleManager<B,L,R>::evaluateSolution(std::vector<std::vector<double> > const& sol) {

	std::vector<int> const& localScen = ctx.localScenarios();
	int nscen = input.nScenarios();
	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		cutInfo cut = solveSubproblem(sol[scen],scen);
		bundle[scen].push_back(cut);
	}


	int nvar1 = input.nFirstStageVars();
	double obj = 0.;
	for (int scen = 0; scen < nscen; scen++) {
		int proc = ctx.owner(scen);
		int i;
		if (ctx.mype() != proc) {
			i = bundle[scen].size();
			bundle[scen].push_back(cutInfo());
			bundle[scen][i].primalSol.resize(nvar1);
			bundle[scen][i].evaluatedAt.resize(nvar1);
			bundle[scen][i].subgradient.resize(nvar1);
		} else {
			i = bundle[scen].size()-1;
		}
		MPI_Bcast(&bundle[scen][i].primalSol[0],nvar1,MPI_DOUBLE,proc,ctx.comm());
		MPI_Bcast(&bundle[scen][i].evaluatedAt[0],nvar1,MPI_DOUBLE,proc,ctx.comm());
		MPI_Bcast(&bundle[scen][i].subgradient[0],nvar1,MPI_DOUBLE,proc,ctx.comm());
		MPI_Bcast(&bundle[scen][i].objval,1,MPI_DOUBLE,proc,ctx.comm());
		MPI_Bcast(&bundle[scen][i].objmax,1,MPI_DOUBLE,proc,ctx.comm());
		obj -= bundle[scen][i].objmax; // note cut has obj flipped
	}
	return obj;
}

template<typename B, typename L, typename R> void bundleManager<B,L,R>::evaluateAndUpdate() {

	currentObj = evaluateSolution(currentSolution);
}

template<typename B, typename L, typename R> void bundleManager<B,L,R>::checkLastPrimals() {

	int nscen = input.nScenarios();
	//for (int i = 0; i < 5; i++) {
	for (int i = 0; i < nscen; i++) {
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
		return;
	}

	checkLastPrimals();
	doStep();
	
}

#endif
