#ifndef PROXLINFTRUSTMANGAGER_HPP
#define PROXLINFTRUSTMANGAGER_HPP

#include "bundleManager.hpp"

#include "proximalBAQP.hpp"

// like levelManager, but only accept a step if it's in the right direction
template<typename BAQPSolver, typename LagrangeSolver, typename RecourseSolver>
class proxQPManager : public bundleManager<BAQPSolver,LagrangeSolver,RecourseSolver> {
public:
	proxQPManager(stochasticInput &input, BAContext & ctx) : 
		bundleManager<BAQPSolver,LagrangeSolver,RecourseSolver>(input,ctx) {
		int nscen = input.nScenarios();
		trialSolution.resize(nscen,std::vector<double>(input.nFirstStageVars(),0.));
		tau = 1e-2;
		t = MPI_Wtime();
	}


protected: 
	virtual void doStep() {
		using namespace std;
		
		int nvar1 = this->input.nFirstStageVars();
		int nscen = this->input.nScenarios();
		vector<int> const &localScen = this->ctx.localScenarios();
		
		if (this->ctx.mype() == 0) printf("Iter %d Current Objective: %f Best Primal: %f, Relerr: %g Elapsed: %f\n",this->nIter-1,this->currentObj,this->bestPrimalObj,fabs(lastModelObj-this->currentObj)/fabs(this->currentObj),MPI_Wtime()-t);
		if (this->terminated_) return;	
		

		proximalQPModel lm(nvar1,this->bundle,this->currentSolution,tau);
		BAQPSolver solver(lm);
		

		solver.go();
		lastModelObj = solver.getObjective();
	
		std::vector<double> const& y = solver.getFirstStagePrimalColSolution();
		std::vector<double> solSum(nvar1,0.);
		for (unsigned r = 1; r < localScen.size(); r++) {
			int scen = localScen[r];
			std::vector<double> const& z = solver.getSecondStagePrimalColSolution(scen);
			for (int k = 0; k < nvar1; k++) {
				double sum = y[k]/sqrt((double)nscen);
				double zsum = 0.;
				for (unsigned i = 0; i < this->bundle[scen].size(); i++) {
					sum -= this->bundle[scen][i].subgradient[k]*z[i];
					zsum += z[i];
					assert(z[i] > -1e-7);
				}
				assert(fabs(zsum-1.) < 1e-7);
				trialSolution[scen][k] = this->currentSolution[scen][k]+sum/tau;
				solSum[k] += trialSolution[scen][k];
			}
		}
		for (int k = 0; k < nvar1; k++) {
			cout << trialSolution[0][k] << " ";
			assert(fabs(solSum[k]) < 1e-7);
		}
		cout << endl;

		double newObj = this->evaluateSolution(trialSolution);
		if (this->ctx.mype() == 0) cout << "Trial solution has obj = " << newObj;
		if (newObj > this->currentObj) {
			if (this->ctx.mype() == 0) cout << ", accepting\n";
			swap(this->currentSolution,trialSolution);
			this->currentObj = newObj;
			if (fabs(lastModelObj-this->currentObj)/fabs(this->currentObj) < this->relativeConvergenceTol) {
				this->terminated_ = true;
				//checkLastPrimals();
			}

		} else {
			if (this->ctx.mype() == 0) cout << ", null step\n";
		}


	}

private:
	double lastModelObj;
	double tau;
	double t;

	std::vector<std::vector<double> > trialSolution; 


};


#endif
