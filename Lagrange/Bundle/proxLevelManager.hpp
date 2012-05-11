#ifndef PROXLEVELMANGAGER_HPP
#define PROXLEVELMANGAGER_HPP

#include "bundleManager.hpp"

#include "cuttingPlaneBALP.hpp"
#include "levelBALP.hpp"

// like levelManager, but only accept a step if it's in the right direction
template<typename BALPSolver, typename LagrangeSolver, typename RecourseSolver>
class proxLevelManager : public bundleManager<BALPSolver,LagrangeSolver,RecourseSolver> {
public:
	// levelParameter is \lambda as defined on p.112 "New variants of bundle methods", lemarechal et al.
	proxLevelManager(stochasticInput &input, BAContext & ctx, double levelParam) : 
		bundleManager<BALPSolver,LagrangeSolver,RecourseSolver>(input,ctx), levelParam(levelParam) {
		int nscen = input.nScenarios();
		trialSolution.resize(nscen,std::vector<double>(input.nFirstStageVars(),0.));
		rows2.resize(nscen);
		rows2l.resize(nscen);
		cols2.resize(nscen);
		cols2l.resize(nscen);
		assert(levelParam >= 0. && levelParam <= 1.);
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
		

		
		{
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
		}

		if (lastModelObj < this->bestPrimalObj - 1e-5) levelParam = 0.5;

		double newLevel = levelParam*(-this->currentObj)+(1.-levelParam)*(-lastModelObj);

		cout << "Target Level = " << -newLevel << endl;
		{
			levelModel lm(nvar1,this->bundle,newLevel,this->currentSolution);
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
	double levelParam;
	// saved states for cutting plane lp
	std::vector<std::vector<variableState> > rows2, cols2;
	std::vector<variableState> cols1;
	// saved states for level set lp
	std::vector<std::vector<variableState> > rows2l, cols2l;
	std::vector<variableState> cols1l;

	std::vector<std::vector<double> > trialSolution; 


};


#endif
