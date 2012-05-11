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
	}


protected: 
	virtual void doStep() {
		using namespace std;
		
		int nscen = this->input.nScenarios();
		int nvar1 = this->input.nFirstStageVars();
		if (fabs(lastModelObj-this->currentObj)/fabs(this->currentObj) < this->relativeConvergenceTol) {
			this->terminated_ = true;
			//checkLastPrimals();
		}
		printf("Iter %d Current Objective: %f Model Objective: %f Best Primal: %f\n",this->nIter-1,this->currentObj,lastModelObj,this->bestPrimalObj);
		if (this->terminated_) return;	
		
		cuttingPlaneModel cpm(nvar1,this->bundle,-this->bestPrimalObj);

		

		BALPSolver solver(cpm,this->ctx, (this->nIter == 1) ? BALPSolver::useDual : BALPSolver::usePrimal);

		if (this->nIter > 1) {
			for (int k = 0; k < nvar1+1; k++) {
				solver.setFirstStageColState(k,cols1[k]);
			}
			for (int i = 0; i < nscen; i++) {
				for (int k = 0; k < nvar1+1; k++) {
					solver.setSecondStageRowState(i,k,rows2[i][k]);
				}
				cols2[i].push_back(AtLower);
				for (unsigned k = 0; k < this->bundle[i].size(); k++) {
					solver.setSecondStageColState(i,k,cols2[i][k]);
				}
			}
			solver.commitStates();
		}

		solver.go();
		lastModelObj = solver.getObjective();

		cols1.resize(nvar1+1);
		for (int k = 0; k < nvar1+1; k++) {
			cols1[k] = solver.getFirstStageColState(k);
		}
		for (int i = 0; i < nscen; i++) {
			rows2[i].resize(nvar1+1);
			for (int k = 0; k < nvar1+1; k++) {
				rows2[i][k] = solver.getSecondStageRowState(i,k);
			}
			cols2[i].resize(this->bundle[i].size());
			for (unsigned k = 0; k < this->bundle[i].size(); k++) {
				cols2[i][k] = solver.getSecondStageColState(i,k);
			}
		}

		for (int i = 0; i < nscen; i++) {
			std::vector<double> const& iterate = solver.getSecondStageDualRowSolution(i);
			for (int k = 0; k < nvar1; k++) {
				this->currentSolution[i][k] = -iterate[k+1];
			}
		}
		this->evaluateAndUpdate();


	}

private:
	double lastModelObj;
	// saved states for cutting plane lp
	std::vector<std::vector<variableState> > rows2, cols2;
	std::vector<variableState> cols1;

};


#endif
