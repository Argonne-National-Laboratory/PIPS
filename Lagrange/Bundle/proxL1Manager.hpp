#ifndef PROXLEVELMANGAGER_HPP
#define PROXLEVELMANGAGER_HPP

#include "bundleManager.hpp"

#include "cuttingPlaneBALP.hpp"
#include "l1bundleBALP.hpp"

// like levelManager, but only accept a step if it's in the right direction
template<typename BALPSolver, typename LagrangeSolver, typename RecourseSolver>
class proxL1Manager : public bundleManager<BALPSolver,LagrangeSolver,RecourseSolver> {
public:
	proxL1Manager(stochasticInput &input, BAContext & ctx) : 
		bundleManager<BALPSolver,LagrangeSolver,RecourseSolver>(input,ctx) {
		int nscen = input.nScenarios();
		trialSolution.resize(nscen,std::vector<double>(input.nFirstStageVars(),0.));
		rows2l.resize(nscen);
		cols2l.resize(nscen);
	}


protected: 
	virtual void doStep() {
		using namespace std;
		
		int nscen = this->input.nScenarios();
		int nvar1 = this->input.nFirstStageVars();
		/*if (fabs(lastModelObj-this->currentObj)/fabs(this->currentObj) < this->relativeConvergenceTol) {
			this->terminated_ = true;
			//checkLastPrimals();
		}*/
		printf("Iter %d Current Objective: %f\n",this->nIter-1,this->currentObj);
		if (this->terminated_) return;	
		

		{
			l1bundleModel lm(nvar1,this->bundle,100,this->currentSolution);
			BALPSolver solver(lm,this->ctx, BALPSolver::useDual);
			
			if (this->nIter > 1) {
				for (int k = 0; k < nvar1+1; k++) {
					solver.setFirstStageColState(k,cols1l[k]);
				}
				for (int i = 0; i < nscen; i++) {
					for (int k = 0; k < 2*nvar1+1; k++) {
						solver.setSecondStageRowState(i,k,rows2l[i][k]);
					}
					cols2l[i].push_back(AtLower);
					for (unsigned k = 0; k < 2*nvar1+this->bundle[i].size(); k++) {
						solver.setSecondStageColState(i,k,cols2l[i][k]);
					}
				}
				solver.commitStates();
			}



			solver.go();

			cols1l.resize(nvar1+1);
			for (int k = 0; k < nvar1+1; k++) {
				cols1l[k] = solver.getFirstStageColState(k);
			}
			for (int i = 0; i < nscen; i++) {
				rows2l[i].resize(2*nvar1+1);
				for (int k = 0; k < 2*nvar1+1; k++) {
					rows2l[i][k] = solver.getSecondStageRowState(i,k);
				}
				cols2l[i].resize(2*nvar1+this->bundle[i].size());
				for (unsigned k = 0; k < 2*nvar1+this->bundle[i].size(); k++) {
					cols2l[i][k] = solver.getSecondStageColState(i,k);
				}
			}


			for (int i = 0; i < nscen; i++) {
				std::vector<double> const& iterate = solver.getSecondStageDualRowSolution(i);
				for (int k = 0; k < nvar1; k++) {
					trialSolution[i][k] = -iterate[k+1];
				}
			}

			double newObj = this->evaluateSolution(trialSolution);
			cout << "Trial solution has obj = " << newObj;
			if (newObj > this->currentObj) {
				cout << ", accepting\n";
				swap(this->currentSolution,trialSolution);
				this->currentObj = newObj;
			} else {
				cout << ", null step\n";
			}
		}


	}

private:
	double lastModelObj;
	// saved states for regularized cutting plane lp
	std::vector<std::vector<variableState> > rows2l, cols2l;
	std::vector<variableState> cols1l;

	std::vector<std::vector<double> > trialSolution; 


};


#endif
