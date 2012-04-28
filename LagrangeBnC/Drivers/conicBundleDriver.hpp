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

		LagrangeSolver lsol(*input,scen,lagrangeDiff);
		lsol.go();
		assert(lsol.getStatus() == Optimal);
		objective_value = -lsol.getBestPossibleObjective();

		std::vector<double> subgrad(nvar1*(nscen-1),0.0);
		std::vector<double> sol = lsol.getBestFirstStageSolution();
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
			
private:
	stochasticInput *input;
	int scen;
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
	solver.set_out(&cout,1);
	solver.set_term_relprec(1e-4);
	double t = MPI_Wtime();
	do {
		solver.do_descent_step();
		cout << "Lagrange Objective: " << -solver.get_objval() << " Elapsed: " << MPI_Wtime()-t << endl;
	} while (!solver.termination_code());

	solver.print_termination_code(cout);



}


#endif
