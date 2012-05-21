#ifndef PROXLINFTRUSTMANGAGER_HPP
#define PROXLINFTRUSTMANGAGER_HPP

#include "bundleManager.hpp"

#include "cuttingPlaneBALP.hpp"
#include "lInfTrustBALP.hpp"
#include "lInfTrustBALP2.hpp"

// like levelManager, but only accept a step if it's in the right direction
template<typename BALPSolver, typename LagrangeSolver, typename RecourseSolver>
class proxLinfTrustManager : public bundleManager<BALPSolver,LagrangeSolver,RecourseSolver> {
public:
	proxLinfTrustManager(stochasticInput &input, BAContext & ctx) : 
		bundleManager<BALPSolver,LagrangeSolver,RecourseSolver>(input,ctx) {
		int nscen = input.nScenarios();
		trialSolution.resize(nscen,std::vector<double>(input.nFirstStageVars(),0.));
		rows2l.resize(nscen);
		cols2l.resize(nscen);
		maxRadius = 10.;
		curRadius = maxRadius/10.;
		counter = 0;
		t = MPI_Wtime();
	}


protected: 
	virtual void doStep() {
		using namespace std;
		
		int nscen = this->input.nScenarios();
		int nvar1 = this->input.nFirstStageVars();
		vector<int> const &localScen = this->ctx.localScenarios();
		
		if (this->ctx.mype() == 0) printf("Iter %d Current Objective: %f Best Primal: %f, Relerr: %g Elapsed: %f\n",this->nIter-1,this->currentObj,this->bestPrimalObj,fabs(lastModelObj-this->currentObj)/fabs(this->currentObj),MPI_Wtime()-t);
		if (this->terminated_) return;	
		

		{
			lInfTrustModel lm(nvar1,this->bundle,curRadius,this->currentSolution);
			BALPSolver solver(lm,this->ctx, (this->nIter > 1) ? BALPSolver::usePrimal : BALPSolver::useDual);
			
			if (this->nIter > 1) {
				for (int k = 0; k < lm.nFirstStageVars(); k++) {
					solver.setFirstStageColState(k,cols1l[k]);
				}
				for (unsigned r = 1; r < localScen.size(); r++) {
					int scen = localScen[r];
					for (int k = 0; k < lm.nSecondStageCons(scen); k++) {
						solver.setSecondStageRowState(scen,k,rows2l[scen][k]);
					}
					cols2l[scen].push_back(AtLower);
					for (int k = 0; k < lm.nSecondStageVars(scen); k++) {
						solver.setSecondStageColState(scen,k,cols2l[scen][k]);
					}
				}
				solver.commitStates();
			}



			solver.go();
			lastModelObj = solver.getObjective();

			cols1l.resize(lm.nFirstStageVars());
			for (int k = 0; k < lm.nFirstStageVars(); k++) {
				cols1l[k] = solver.getFirstStageColState(k);
			}
			for (unsigned r = 1; r < localScen.size(); r++) {
				int scen = localScen[r];
				rows2l[scen].resize(lm.nSecondStageCons(scen));
				for (int k = 0; k < lm.nSecondStageCons(scen); k++) {
					rows2l[scen][k] = solver.getSecondStageRowState(scen,k);
				}
				cols2l[scen].resize(lm.nSecondStageVars(scen));
				for (int k = 0; k < lm.nSecondStageVars(scen); k++) {
					cols2l[scen][k] = solver.getSecondStageColState(scen,k);
				}
			}

			double maxdifftemp = 0.;
			for (unsigned r = 1; r < localScen.size(); r++) {
				int scen = localScen[r];
				std::vector<double> const& iterate = solver.getSecondStageDualRowSolution(scen);
				for (int k = 0; k < nvar1; k++) {
					trialSolution[scen][k] = -iterate[k+1];
					maxdifftemp = max(fabs(trialSolution[scen][k]-this->currentSolution[scen][k]),maxdifftemp);
				}
			}
			double maxdiff;
			MPI_Allreduce(&maxdifftemp,&maxdiff,1,MPI_DOUBLE,MPI_MAX,this->ctx.comm());

			double newObj = this->evaluateSolution(trialSolution);
			if (this->ctx.mype() == 0) cout << "Trial solution has obj = " << newObj;
			if (newObj > this->currentObj) {
				if (this->ctx.mype() == 0) cout << ", accepting\n";
				if (maxdiff > curRadius - 1e-5 && newObj > this->currentObj + 0.5*(lastModelObj-this->currentObj)) {
					curRadius = min(2.*curRadius,maxRadius);
					if (this->ctx.mype() == 0) cout << "Increased trust region radius to " << curRadius << endl;
				}
				swap(this->currentSolution,trialSolution);
				this->currentObj = newObj;
				counter = 0;
				if (fabs(lastModelObj-this->currentObj)/fabs(this->currentObj) < this->relativeConvergenceTol) {
					this->terminated_ = true;
					//checkLastPrimals();
				}

			} else {
				if (this->ctx.mype() == 0) cout << ", null step\n";
				double rho = min(1.,curRadius)*(newObj-this->currentObj)/(this->currentObj-lastModelObj);
				if (this->ctx.mype() == 0) cout << "Rho: " << rho << " counter: " << counter << endl;
				if (rho > 0) counter++;
				if (rho > 3 || (counter >= 3 && 1. < rho && rho <= 3.)) {
					curRadius *= 1./min(rho,4.);
					if (this->ctx.mype() == 0) cout << "Decreased trust region radius to " << curRadius << endl;
					counter = 0;
				}
			}
		}


	}

private:
	double lastModelObj;
	double maxRadius, curRadius;
	double t;
	int counter;
	// saved states for regularized cutting plane lp
	std::vector<std::vector<variableState> > rows2l, cols2l;
	std::vector<variableState> cols1l;

	std::vector<std::vector<double> > trialSolution; 


};


#endif
