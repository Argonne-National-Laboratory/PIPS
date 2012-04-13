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
	int mype = ctx.mype();

	const vector<int> &localScen = ctx.localScenarios();
	vector<double> lagrangeObjs;
	vector<vector<double> > lagrangeSolutions;
	double objsum_local = 0.0;
	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		LagrangeSolver lsol(input, scen, vector<double>(nvar1,0.0));
		lsol.go();
		printf("Objective from scen %d: %f\n",scen,lsol.getBestPossibleObjective());
		objsum_local += lsol.getBestPossibleObjective();
		lagrangeObjs.push_back(lsol.getBestPossibleObjective());
		lagrangeSolutions.push_back(lsol.getBestFirstStageSolution());
	}
	double lagrangelb;
	MPI_Allreduce(&objsum_local,&lagrangelb,1,MPI_DOUBLE,MPI_SUM,comm);

	const vector<double> &obj1 = input.getFirstStageObj();


	
	// evaluate solutions
	vector<double> bestSolution;
	double bestObj = COIN_DBL_MAX;

	assert(input.scenarioDimensionsEqual());
	int nvar2 = input.nSecondStageVars(0);
	int ncons2 = input.nSecondStageCons(0);

	vector<variableState> rowSave(ncons2), colSave(nvar2);
	bool havesave = false;

	vector<vector<double> > uniqueSolutions;
	vector<double> objs;
	vector<int> freqCount;
	for (int scen_ = 0; scen_ < nscen; scen_++) {
		vector<double> curSolution(nvar1);
		
		if (mype == ctx.owner(scen_)) {
			int idx = find(localScen.begin(),localScen.end(),scen_)-localScen.begin()-1;
			curSolution = lagrangeSolutions[idx];
		}
		MPI_Bcast(&curSolution[0],nvar1,MPI_DOUBLE,ctx.owner(scen_),comm);

		bool diff = true;
		for (unsigned r = 0; r < uniqueSolutions.size(); r++) {
			const vector<double>& v = uniqueSolutions[r];
			diff = false;
			for (int i = 0; i < nvar1; i++) {
				if (fabs(v[i]-curSolution[i]) > 1e-5) {
					diff = true;
					break;
				}
			}
			if (!diff) {
				freqCount[r]++;
				break;
			} else {
				diff = true;
			}
		}
		if (!diff) continue;

		uniqueSolutions.push_back(curSolution);
		freqCount.push_back(1);


		double sum = 0.0;
		bool infeas = false;
		for (unsigned q = 1; q < localScen.size(); q++) {
			int scen = localScen[q];
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
		if (infeas) { sum += COIN_DBL_MAX; }
		double sum_all;
		MPI_Allreduce(&sum,&sum_all,1,MPI_DOUBLE,MPI_SUM,comm);
		for (int k = 0; k < nvar1; k++) sum_all += curSolution[k]*obj1[k];
		
		objs.push_back(sum_all);


		if (sum_all < bestObj) {
			bestObj = sum_all;
			bestSolution = curSolution;
		}

	}
	if (mype == 0) {
		printf("%lu unique solutions\nCount\tObj\n", uniqueSolutions.size());
		for (unsigned r = 0; r < uniqueSolutions.size(); r++) {
			printf("%d\t", freqCount[r]);
			if (objs[r] > 1e20) {
				printf("inf\n");
			} else {
				printf("%f\n",objs[r]);
			}
		}
		printf("Lagrange LB: %f\n",lagrangelb);
		printf("Best UB: %f\n",bestObj);
		// formula for positive values
		printf("Gap: %f%%\n",100*(bestObj-lagrangelb)/lagrangelb);
	}


}


#endif
