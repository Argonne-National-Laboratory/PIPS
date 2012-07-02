#ifndef CONICBUNDLEPARDRIVER1HPP
#define CONICBUNDLEPARDRIVER1HPP

#include "conicBundleDriver.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>

// use static assignment of scenarios

template<typename LagrangeSolver> class lagrangeSubproblemMaster : public FunctionOracle {
public:
	lagrangeSubproblemMaster(stochasticInput &input, BAContext const &ctx) :
		input(input),ctx(ctx) {
		int nscen = input.nScenarios();
		int nvar1 = input.nFirstStageVars();
		ndual = (nscen-1)*nvar1;
		dummies.resize(nscen);
		
		const std::vector<int> &localScen = ctx.localScenarios();
		for (unsigned i = 1; i < localScen.size(); i++) {
			int scen = localScen[i];
			inner.push_back(lagrangeSubproblem<LagrangeSolver>(&input,scen));
		}
		for (int i = 0; i < nscen; i++) {
			dummies[i].subgrad.resize(ndual);
		}
		nIter = 0;
		bestobj = -1e20;
		startTime = MPI_Wtime();
		subproblemTime = 0.;
			
	}

	virtual int evaluate(  const DVector& dual, /* argument/Lagrange multipliers */
                         double relprec,
		         double&   objective_value,
			 DVector&  cut_vals,
			 std::vector<DVector>&    subgradients,
			 std::vector<PrimalData*>&     primal_solutions,
			 PrimalExtender*& e
		      ) {
		nIter++;	
		double t = MPI_Wtime();
		const std::vector<int> &localScen = ctx.localScenarios();
		for (unsigned i = 1; i < localScen.size(); i++) {
			int scen = localScen[i];
			DVector cut_vals;
			std::vector<DVector> subgrads;
			std::vector<PrimalData*> primalsol;

			inner[i-1].evaluate(dual,relprec,dummies[scen].objval,cut_vals,
				subgrads,primalsol,e);
			std::swap(dummies[scen].subgrad,subgrads[0]);
			dummies[scen].cutval = cut_vals[0];
			dummies[scen].primalsol = inner[i-1].getSavedSolutions().at(nIter-1);
		}
		int nscen = input.nScenarios();
		int nvar1 = input.nFirstStageVars();
		/*std::ofstream f;
		std::stringstream ss;
		ss << "cbsols" << nIter;
		if (ctx.mype() == 0) f.open(ss.str().c_str());*/
		double allobj = 0.;
		for (int scen = 0; scen < nscen; scen++) {
			MPI_Bcast(&dummies[scen].subgrad[0],ndual,MPI_DOUBLE,
				ctx.owner(scen),ctx.comm());
			MPI_Bcast(&dummies[scen].cutval,1,MPI_DOUBLE,
				ctx.owner(scen),ctx.comm());
			MPI_Bcast(&dummies[scen].objval,1,MPI_DOUBLE,
				ctx.owner(scen),ctx.comm());
			allobj -= dummies[scen].objval;
			std::vector<double> sol(nvar1);
			if (ctx.mype() == ctx.owner(scen)) {
				sol = dummies[scen].primalsol;
			}
			MPI_Bcast(&sol[0],nvar1,MPI_DOUBLE,ctx.owner(scen),ctx.comm());
			/*if (ctx.mype() == 0) {
				for (int k = 0; k < nvar1; k++) {
					f << sol[k] << " ";
				}
				f << std::endl;
			}*/
		}
		subproblemTime += MPI_Wtime()-t;
		bestobj = std::max(allobj,bestobj);
		if (ctx.mype() == 0) {
			double elap = MPI_Wtime()-startTime;
			printf("Iter %d Current Objective: %g Elapsed: %f (%f in QP Solve)\n",
				nIter-1,bestobj,elap,elap-subproblemTime);
		}
		return dummies[0].evaluate(dual,relprec,objective_value,cut_vals,
			subgradients,primal_solutions,e);

	}

	class lagrangeSubproblemDummy : public FunctionOracle {
	public:
		virtual int evaluate(  const DVector& dual, /* argument/Lagrange multipliers */
                         double relprec,
		         double&   objective_value,
			 DVector&  cut_vals,
			 std::vector<DVector>&    subgradients,
			 std::vector<PrimalData*>&     primal_solutions,
			 PrimalExtender*&
		      ) {
		      	objective_value = objval;
			subgradients.push_back(subgrad);
			cut_vals.push_back(cutval);

			return 0;

		}


		DVector subgrad;
		double cutval;
		double objval;
		std::vector<double> primalsol;
	
	
	};

	std::vector<lagrangeSubproblemDummy> dummies;
	std::vector<lagrangeSubproblem<LagrangeSolver> > inner;
private:

	stochasticInput &input;
	BAContext const &ctx;
	int ndual;
	int nIter;
	double bestobj;
	double startTime;
	double subproblemTime;

};

template <typename LagrangeSolver, typename RecourseSolver> void conicBundleParDriver1(stochasticInput &input, 
	MPI_Comm comm = MPI_COMM_WORLD) {
	using namespace std;
	using boost::scoped_ptr; 

	BAContext ctx(comm);
	ctx.initializeAssignment(input.nScenarios());

	
	int nvar1 = input.nFirstStageVars();
	int nscen = input.nScenarios();
	int mype = ctx.mype();
	
	CBSolver solver;

	solver.init_problem(nvar1*(nscen-1));
	
	lagrangeSubproblemMaster<LagrangeSolver> master(input,ctx);
	solver.add_function(master);
	//solver.set_max_bundlesize(master,5);
	for (int i = 1; i < nscen; i++) {
		solver.add_function(master.dummies[i]);
		//solver.set_max_bundlesize(master.dummies[i],5);
	}
	if (mype == 0) solver.set_out(&cout,1);
	solver.set_term_relprec(1e-7);
	do {
		solver.do_descent_step();
	} while (!solver.termination_code());

	if (mype == 0) solver.print_termination_code(cout);
	
	/*	
	
	const vector<double> &obj1 = input.getFirstStageObj();
	// test last solutions (only)
	double best = COIN_DBL_MAX;

	const vector<int> &localScen = ctx.localScenarios();
	for (int i = 0; i < nscen; i++) {
			
		vector<double> sol(nvar1);
		if (mype == ctx.owner(i)) {
			int idx = *std::find(localScen.begin(),localScen.end(),i);
			sol = *(master.inner[idx].getSavedSolutions().end()-1);
		}
		MPI_Bcast(&sol[0],nvar1,MPI_DOUBLE,ctx.owner(i),ctx.comm());
		
		double obj = 0.;
		bool infeas = false;
		for (int s = 0; s < nscen; s++) {
			RecourseSolver rsol(input,s,sol);
			rsol.go();
			obj += rsol.getObjective();
			
			if (rsol.getStatus() == ProvenInfeasible) {
				printf("got infeasible 1st stage\n");
				infeas = true; break;
			}
		}
		if (infeas) obj += COIN_DBL_MAX;
		for (int k = 0; k < nvar1; k++) {
			obj += sol[k]*obj1[k];
		}
		best = min(obj,best);
	}

	cout << "Best LB: " << -solver.get_objval() << " Best UB: " << best << endl;
	*/

}


#endif
