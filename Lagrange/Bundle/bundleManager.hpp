#ifndef BUNDLEMANAGER_HPP
#define BUNDLEMANAGER_HPP

#include "BA.hpp"
#include "BALPSolverInterface.hpp"

// TODO: remove?
#include "cuttingPlaneBALP.hpp"



// this is specialized particularly for lagrangian relaxation of nonanticipativity constraints, not general bundle solver
template<typename BALPSolver, typename LagrangeSolver, typename RecourseSolver> class bundleManager {
public:
	bundleManager(stochasticInput &input, BAContext & ctx) : ctx(ctx), input(input) {
		// use zero as initial iterate
		currentSolution.resize(input.nScenarios(),std::vector<double>(input.nFirstStageVars(),0.));
		nIter = -1;
		bundle.resize(input.nScenarios());
		bestPrimalObj = COIN_DBL_MAX;
	}

	cutInfo solveSubproblem(std::vector<double> const& at, int scen); 
	
	double testPrimal(std::vector<double> const& primal); // test a primal solution and return objective (COIN_DBL_MAX) if infeasible

	virtual void iterate();

protected:

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


};


template<typename B,typename L, typename R> cutInfo bundleManager<B,L,R>::solveSubproblem(std::vector<double> const& at, int scen) {
	using namespace std;	
	int nvar1 = input.nFirstStageVars();
	double t = MPI_Wtime();
	

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


template<typename B, typename L, typename R> void bundleManager<B,L,R>::evaluateAndUpdate() {

	int nscen = input.nScenarios();
	currentObj = 0.;
	for (int i = 0; i < nscen; i++) {
		cutInfo cut = solveSubproblem(currentSolution[i],i);
		bundle[i].push_back(cut);
		currentObj -= cut.objval; // note cut has obj flipped
	}
}

template<typename B, typename L, typename R> void bundleManager<B,L,R>::checkLastPrimals() {

	int nscen = input.nScenarios();
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
	using namespace std;

	int nscen = input.nScenarios();
	int nvar1 = input.nFirstStageVars();

	if (nIter++ == -1) { // just evaluate and generate subgradients
		evaluateAndUpdate();
		checkLastPrimals();
	}

	printf("Iter %d Current Objective: %f Best Primal: %f\n",nIter-1,currentObj,bestPrimalObj);
	
	// use cutting plane model for now
	cuttingPlaneModel cpm(nvar1,bundle,-bestPrimalObj);
	
	B solver(cpm,ctx);
	solver.go();
	cout << "Model objective: " << solver.getObjective() << endl;

	for (int i = 0; i < nscen; i++) {
		std::vector<double> const& iterate = solver.getSecondStageDualRowSolution(i);
		for (int k = 0; k < nvar1; k++) {
			currentSolution[i][k] = -iterate[k+1];
		}
	}
	evaluateAndUpdate();
}

#endif
