#ifndef CONICBUNDLEDRIVER_HPP
#define CONICBUNDLEDRIVER_HPP

#include "stochasticInput.hpp"
#include <boost/scoped_ptr.hpp>
#include "CBSolver.hxx"

using namespace ConicBundle;

template<typename LagrangeSolver> class lagrangeSubproblem : public FunctionOracle {
public:
	lagrangeSubproblem(stochasticInput *input, int scen) : input(input), scen(scen) {}

	virtual int evaluate(  const DVector& dual, /* argument/Lagrange multipliers */
                         double relprec,
		         double&   objective_value,
			 DVector&  cut_vals,
			 std::vector<DVector>&    subgradients,
			 std::vector<PrimalData*>&     primal_solutions,
			 PrimalExtender*&
		      ) {
		using namespace std;
		int nvar1 = input->nFirstStageVars();
		int nscen = input->nScenarios();
		std::vector<double> lagrangeDiff(nvar1, 0.0);
		if (scen > 0) {
			// add -\lambda^{scen-1}
			int offset = nvar1*(scen-1);
			for (int i = 0; i < nvar1; i++) {
				lagrangeDiff[i] -= dual[i+offset];
			}
		}
		if (scen < nscen-1) {
			// add \lambda^{scen}
			int offset = nvar1*scen;
			for (int i = 0; i < nvar1; i++) {
				lagrangeDiff[i] += dual[i+offset];
			}
		}
		double t = MPI_Wtime();
		/*for (unsigned i = 0; i < dual.size(); i++) {
			cout << dual[i] << " ";
		}
		cout << endl;*/

		LagrangeSolver lsol(*input,scen,lagrangeDiff);
		//lsol.setRatio(100*fabs(relprec));
		lsol.go();
		//cout << "SCEN " << scen << " DONE, " << MPI_Wtime() - t << " SEC" << endl;

		//assert(lsol.getStatus() == Optimal);
		assert(lsol.getStatus() != ProvenInfeasible);
		objective_value = -lsol.getBestPossibleObjective();

		std::vector<double> subgrad(nvar1*(nscen-1),0.0);
		std::vector<double> sol = lsol.getBestFirstStageSolution();
		sols.push_back(sol);
		if (scen > 0) {
			int offset = nvar1*(scen-1);
			for (int i = 0; i < nvar1; i++) {
				subgrad[i+offset] = sol[i];
			}
		}
		if (scen < nscen-1) {
			int offset = nvar1*scen;
			for (int i = 0; i < nvar1; i++) {
				subgrad[i+offset] = -sol[i];
			}
		}
		subgradients.push_back(subgrad);
		cut_vals.push_back(objective_value);

		return 0;
	}
	std::vector<std::vector<double> > const& getSavedSolutions() const { return sols; }
			
private:
	stochasticInput *input;
	int scen;
	std::vector<std::vector<double> > sols;
};

template <typename LagrangeSolver, typename RecourseSolver> void conicBundleDriver(stochasticInput &input, 
	MPI_Comm comm = MPI_COMM_WORLD) {
	using namespace std;
	using boost::scoped_ptr; 

	BAContext ctx(comm);
	ctx.initializeAssignment(input.nScenarios());

	
	int nvar1 = input.nFirstStageVars();
	int nscen = input.nScenarios();
	int mype = ctx.mype();
	
	assert(mype == 0); // serial only for now

	CBSolver solver;

	solver.init_problem(nvar1*(nscen-1));

	std::vector<lagrangeSubproblem<LagrangeSolver> > funcs;
	funcs.reserve(nscen);

	for (int i = 0; i < nscen; i++) {
		funcs.push_back(lagrangeSubproblem<LagrangeSolver>(&input,i));
		solver.add_function(funcs[i]);
	}
	if (mype == 0) solver.set_out(&cout,1);
	solver.set_term_relprec(1e-6);
	double t = MPI_Wtime();
	do {
		solver.do_descent_step();
		if (mype == 0) cout << "Lagrange Objective: " << -solver.get_objval() << " Elapsed: " << MPI_Wtime()-t << endl;
	} while (!solver.termination_code());

	solver.print_termination_code(cout);

	
	const vector<double> &obj1 = input.getFirstStageObj();
	// test last solutions (only)
	double best = COIN_DBL_MAX;
	for (int i = 0; i < nscen; i++) {
		
		const vector<double> &sol = *(funcs[i].getSavedSolutions().end()-1);
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

	if (mype == 0) cout << "Best LB: " << -solver.get_objval() << " Best UB: " << best << endl;

}


#endif
