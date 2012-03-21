#ifndef LAGRANGEROOTNODE_HPP
#define LAGRANGEROOTNODE_HPP

#include "stochasticInput.hpp"
#include <boost/scoped_ptr.hpp>

template <typename BALPSolver, typename LagrangeSolver, typename RecourseSolver> void lagrangeRootNode(stochasticInput &input, 
	const std::string &LPBasis, MPI_Comm comm = MPI_COMM_WORLD) {
	using namespace std;
	using boost::scoped_ptr; 

	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(input.nScenarios());

	/*
	BALPSolver solver(input, ctx, BALPSolver::useDual);
	solver.loadStatus(LPBasis);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();
*/
	int nvar1 = input.nFirstStageVars();

	const vector<int> &localScen = ctx.localScenarios();
	double objsum = 0.0;
	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		LagrangeSolver lsol(input, scen, vector<double>(nvar1,0.0));
		lsol.go();
		printf("Objective from scen %d: %f\n",scen,lsol.getObjective());
		objsum += lsol.getObjective();
	}

	printf("Lagrange LB: %f\n",objsum);




}


#endif
