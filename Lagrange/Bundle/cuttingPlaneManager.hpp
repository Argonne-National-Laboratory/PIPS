#ifndef CUTTINGPLANEMANGAGER_HPP
#define CUTTINGPLANEMANGAGER_HPP

#include "bundleManager.hpp"

#include "cuttingPlaneBALP.hpp"


// simple cutting plane approach
// next iterate is the solution of the cutting plane model, no regularization
template<typename BALPSolver, typename LagrangeSolver, typename RecourseSolver>
class cuttingPlaneManager : public bundleManager<BALPSolver,LagrangeSolver,RecourseSolver> {
public:
	cuttingPlaneManager(stochasticInput &input, BAContext & ctx) : bundleManager<BALPSolver,LagrangeSolver,RecourseSolver>(input,ctx) {
		int nscen = input.nScenarios();
		rows2.resize(nscen);
		cols2.resize(nscen);
		t = MPI_Wtime();
		t2 = 0;
	}


protected: 
	virtual void doStep() {
		using namespace std;
		
		vector<int> const &localScen = this->ctx.localScenarios();
		int nvar1 = this->input.nFirstStageVars();
		if ((lastModelObj-this->currentObj)/(1.+fabs(this->currentObj)) < this->relativeConvergenceTol) {
			this->terminated_ = true;
		}
		if (this->ctx.mype() == 0) printf("Iter %d Current Objective: %f Model Objective: %f Elapsed: %f (%f in LP solve)\n",this->nIter-1,this->currentObj,lastModelObj,MPI_Wtime()-t,t2);
		if (this->terminated_) return;	
		
		cuttingPlaneModel cpm(nvar1,this->bundle,-this->bestPrimalObj);

		

		BALPSolver solver(cpm,this->ctx, (this->nIter == 1) ? BALPSolver::useDual : BALPSolver::usePrimal);

		if (this->nIter > 1) {
			for (int k = 0; k < cpm.nFirstStageVars(); k++) {
				solver.setFirstStageColState(k,cols1[k]);
			}
			for (unsigned r = 1; r < localScen.size(); r++) {
				int scen = localScen[r];
				for (int k = 0; k < cpm.nSecondStageCons(scen); k++) {
					solver.setSecondStageRowState(scen,k,rows2[scen][k]);
				}
				cols2[scen].push_back(AtLower);
				for (int k = 0; k < cpm.nSecondStageVars(scen); k++) {
					solver.setSecondStageColState(scen,k,cols2[scen][k]);
				}
			}
			solver.commitStates();
		}
		double tstart = MPI_Wtime();
		solver.go();
		t2 += MPI_Wtime() - tstart;
		lastModelObj = solver.getObjective();

		cols1.resize(nvar1+1);
		for (int k = 0; k < nvar1+1; k++) {
			cols1[k] = solver.getFirstStageColState(k);
		}
		for (unsigned r = 1; r < localScen.size(); r++) {
			int scen = localScen[r];
			rows2[scen].resize(cpm.nSecondStageCons(scen));
			for (int k = 0; k < cpm.nSecondStageCons(scen); k++) {
				rows2[scen][k] = solver.getSecondStageRowState(scen,k);
			}
			cols2[scen].resize(cpm.nSecondStageVars(scen));
			for (int k = 0; k < cpm.nSecondStageVars(scen); k++) {
				cols2[scen][k] = solver.getSecondStageColState(scen,k);
			}
		}

		for (unsigned r = 1; r < localScen.size(); r++) {
			int scen = localScen[r];
			std::vector<double> const& iterate = solver.getSecondStageDualRowSolution(scen);
			for (int k = 0; k < nvar1; k++) {
				this->currentSolution[scen][k] = -iterate[k+1];
			}
		}
		
		this->evaluateAndUpdate();


	}

private:
	double lastModelObj;
	double t, t2;
	// saved states for cutting plane lp
	std::vector<std::vector<variableState> > rows2, cols2;
	std::vector<variableState> cols1;

};


#endif
