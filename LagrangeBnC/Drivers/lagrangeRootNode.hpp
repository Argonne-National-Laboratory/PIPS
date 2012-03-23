#ifndef LAGRANGEROOTNODE_HPP
#define LAGRANGEROOTNODE_HPP

#include "stochasticInput.hpp"
#include <boost/scoped_ptr.hpp>

template <typename LagrangeSolver, typename RecourseSolver> void lagrangeRootNode(stochasticInput &input, 
	MPI_Comm comm = MPI_COMM_WORLD) {
	using namespace std;
	using boost::scoped_ptr; 

	BAContext ctx(comm);
	ctx.initializeAssignment(input.nScenarios());

	
	int nvar1 = input.nFirstStageVars();
	int nscen = input.nScenarios();

	const vector<int> &localScen = ctx.localScenarios();
	vector<double> lagrangeObjs;
	vector<vector<double> > lagrangeSolutions;
	double objsum = 0.0;
	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		LagrangeSolver lsol(input, scen, vector<double>(nvar1,0.0));
		lsol.go();
		printf("Objective from scen %d: %f\n",scen,lsol.getBestPossibleObjective());
		objsum += lsol.getBestPossibleObjective();
		lagrangeObjs.push_back(lsol.getBestPossibleObjective());
		lagrangeSolutions.push_back(lsol.getBestFirstStageSolution());
	}
	const vector<double> &obj1 = input.getFirstStageObj();


	
	// evaluate solutions
	vector<double> bestSolution;
	double bestObj = COIN_DBL_MAX;
	double sumObj = 0.0;

	assert(input.scenarioDimensionsEqual());
	int nvar2 = input.nSecondStageVars(0);
	int ncons2 = input.nSecondStageCons(0);

	vector<variableState> rowSave(ncons2), colSave(nvar2);
	bool havesave = false;
	for (unsigned i = 0; i < lagrangeSolutions.size(); i++) {
		double sum = 0.0;
		bool infeas = false;
		const vector<double> &curSolution = lagrangeSolutions[i]; 
		for (int k = 0; k < nvar1; k++) sum += curSolution[k]*obj1[k];
		for (int scen = 0; scen < nscen; scen++) {
			RecourseSolver rsol(input, scen, curSolution);
			rsol.setDualObjectiveLimit(1e7);
			
			if (havesave) {
				for (int r = 0; r < nvar2; r++) {
					rsol.setSecondStageColState(r,colSave[r]);
				}
				for (int r = 0; r < ncons2; r++) {
					rsol.setSecondStageRowState(r,rowSave[r]);
				}
			}
			rsol.go();
			sum += rsol.getObjective();
			
			if (rsol.getStatus() == ProvenInfeasible) {
				printf("got infeasible 1st stage\n");
				infeas = true; break;
			}
			assert(rsol.getStatus() == Optimal);
			if (!havesave) {
				for (int r = 0; r < nvar2; r++) {
					colSave[r] = rsol.getSecondStageColState(r);
				}
				for (int r = 0; r < ncons2; r++) {
					rowSave[r] = rsol.getSecondStageRowState(r);
				}
				havesave = true;
			}

		}
		if (infeas) { sumObj += COIN_DBL_MAX; continue; }
		sumObj += sum;

		if (sum < bestObj) {
			bestObj = sum;
			bestSolution = curSolution;
		}

	}

	printf("Lagrange LB: %f\n",objsum);
	printf("Best UB: %f\n",bestObj);
	printf("Average UB: %f\n",sumObj/lagrangeSolutions.size());
	// formula for positive values
	printf("Gap: %f%%\n",100*(bestObj-objsum)/objsum);


}


#endif
