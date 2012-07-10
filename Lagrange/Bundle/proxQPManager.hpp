#ifndef PROXLINFTRUSTMANGAGER_HPP
#define PROXLINFTRUSTMANGAGER_HPP

#include "bundleManager.hpp"

#include "proximalBAQP.hpp"

class filter {
public:
	bool operator()(cutInfo const&c) { return c.z < 1e-7; }
};

// like levelManager, but only accept a step if it's in the right direction
template<typename BAQPSolver, typename LagrangeSolver, typename RecourseSolver>
class proxQPManager : public bundleManager<BAQPSolver,LagrangeSolver,RecourseSolver> {
public:
	proxQPManager(stochasticInput &input, BAContext & ctx) : 
		bundleManager<BAQPSolver,LagrangeSolver,RecourseSolver>(input,ctx) {
		int nscen = input.nScenarios();
		trialSolution.resize(nscen);//,std::vector<double>(input.nFirstStageVars(),0.));
		lastModelObj.resize(nscen);
		//proxCenterModelObj.resize(nscen);
		if (BAQPSolver::isDistributed()) {
			localScen = ctx.localScenarios();
		} else {
			localScen.resize(nscen+1);
			for (int i = -1; i < nscen; i++) localScen[i+1] = i;
		}
		for (unsigned i = 1; i < localScen.size(); i++) {
			int scen = localScen[i];
			trialSolution[scen].resize(input.nFirstStageVars());
		}
		u = 1.;
		t = MPI_Wtime();
		t2 = 0.;
		mL = 0.1;
		mR = 0.1;
		eps_sol = 0.; // absolute scale, need to adjust
	}


protected: 
	virtual void doStep() {
		using namespace std;
		
		int nvar1 = this->input.nFirstStageVars();
		int nscen = this->input.nScenarios();

		/*if (this->nIter == 1) {
			for (unsigned r = 1; r < localScen.size(); r++) {
				int scen = localScen[r];
				proxCenterModelObj[scen] = this->bundle[scen].at(0).objmax;
			}
		}*/
		
		if (this->ctx.mype() == 0) printf("Iter %d Current Objective: %f Relerr: %g Elapsed: %f (%f in QP solve)\n",this->nIter-1,this->currentObj,fabs(lastModelObjSum-this->currentObj)/(1.+fabs(this->currentObj)),MPI_Wtime()-t,t2);
		if (this->terminated_) return;	
		

		proximalQPModel lm(nvar1,this->bundle,this->currentSolution,u);
		BAQPSolver solver(lm);
		
		double tstart = MPI_Wtime();
		solver.go();
		t2 += MPI_Wtime() - tstart;
	
		std::vector<double> const& y = solver.getFirstStagePrimalColSolution();
		double lastModelObjSum_this = 0.;
		//double v_this = 0.;
		for (unsigned r = 1; r < localScen.size(); r++) { // for distributed solver
			int scen = localScen[r];
			std::vector<double> const& z = solver.getSecondStagePrimalColSolution(scen);
			assert(trialSolution[scen].size());
			for (int k = 0; k < nvar1; k++) {
				double sum = y[k]/sqrt((double)nscen);
				for (unsigned i = 0; i < this->bundle[scen].size(); i++) {
					sum -= this->bundle[scen][i].subgradient[k]*z[i];
				}
				trialSolution[scen][k] = this->currentSolution[scen][k]+sum/u;
			}
			lastModelObj[scen] = -solver.getSecondStageDualRowSolution(scen)[0];
			lastModelObjSum_this += -lastModelObj[scen];
			//v_this += lastModelObj[scen] - proxCenterModelObj[scen];// - eps_sol/nscen;
			for (unsigned i = 0; i < this->bundle[scen].size();i++) {
				this->bundle[scen][i].z = z[i];
			}
		
			//this->bundle[scen].erase(std::remove_if(this->bundle[scen].begin(),this->bundle[scen].end(),filter()),this->bundle[scen].end());
				
		}
		//double v;
		if (localScen.size() == (size_t)nscen+1) {
			//v = v_this;
			lastModelObjSum = lastModelObjSum_this;
		} else {
			MPI_Allreduce(&lastModelObjSum_this,&lastModelObjSum,1,MPI_DOUBLE,MPI_SUM,this->ctx.comm());
			//MPI_Allreduce(&v_this,&v,1,MPI_DOUBLE,MPI_SUM,this->ctx.comm());
		}
		double v = this->currentObj-lastModelObjSum; 
		eps_sol = -(mR-mL)*v; 
		if (this->ctx.mype() == 0) {
			//printf("v^k = %g, eps_sol = %g\n",v, eps_sol);
			printf("v^k = %g\n",v);
			stringstream ss;
			ss << "primalconv" << this->nIter;
			ofstream f(ss.str().c_str());
			std::vector<double> y = solver.getFirstStagePrimalColSolution();
			for (unsigned k = 0; k < y.size(); k++) {
				y[k] /= -sqrt((double)nscen);
				f << y[k] << " ";
			}
			f << endl;
		}
		if (-v/(1.+fabs(this->currentObj)) < this->relativeConvergenceTol) {
			this->terminated_ = true; this->nIter++; doStep(); // for printout
			/*
						double val = this->testPrimal(y);
			printf("Primal Obj: %.10g\n",val);*/
			//checkLastPrimals();
			return;
		}
		

		double newObj = this->evaluateSolution(trialSolution, eps_sol/nscen);
		
		// update tau according to p.536 of kiwiel's 1995 paper
		double u_int = 2.*u*(1-(this->currentObj-newObj)/v);
		u = min(max(max(u_int,u/10.),1e-4),10.*u);
		if (this->ctx.mype() == 0) printf("updated u to %g\n",u);

		if (this->ctx.mype() == 0) cout << "Trial solution has obj = " << newObj;
		if (newObj - this->currentObj > -mL*v) {
			if (this->ctx.mype() == 0) cout << ", accepting\n";
			swap(this->currentSolution,trialSolution);
			/*for (unsigned r = 1; r < localScen.size(); r++) {
				int scen = localScen[r];
				proxCenterModelObj[scen] = this->bundle[scen][this->bundle[scen].size()-1].objmax;
			}*/
			this->currentObj = newObj;
			

		} else {
			if (this->ctx.mype() == 0) cout << ", null step\n";
		}
		


	}

private:
	std::vector<double> lastModelObj;
	//std::vector<double> proxCenterModelObj;
	std::vector<int> localScen;
	double lastModelObjSum;
	double u;
	double t, t2;
	double mL, mR;
	double eps_sol; // allowable imprecision in the solution 

	std::vector<std::vector<double> > trialSolution; 


};


#endif
