#ifndef BOOTSTRAPDRIVER_HPP
#define BOOTSTRAPDRIVER_HPP

#include "rawInput.hpp"
#include "ClpSimplex.hpp"
#include <boost/scoped_ptr.hpp>
#include "BALPSolverInterface.hpp"

#include <cstdlib>

static variableState convertClpStatus(ClpSimplex::Status s) {
	if (s == ClpSimplex::basic) return Basic;
	else if (s == ClpSimplex::atLowerBound || s == ClpSimplex::isFixed) return AtLower;
	else if (s == ClpSimplex::atUpperBound) return AtUpper;
	else { assert(0); return AtLower; }
}


template <typename SolverInterface> void bootstrapDriver(const std::string &datarootname, const std::string &basein, 
	const std::string &baseout, int nscenIn, int nscen, MPI_Comm comm = MPI_COMM_WORLD) {
	using namespace std;
	using boost::scoped_ptr; // replace with unique_ptr for C++11

	int mype;
	MPI_Comm_rank(comm,&mype);

	if (mype == 0) printf("Initializing data interface\n");
	scoped_ptr<rawInput> s(new rawInput(datarootname,nscen,comm));
	assert(s->onlyBoundsVary());
	BAContext ctx(comm);
	ctx.initializeAssignment(nscen);
	/*
	if (mype == 0) printf("Loading data\n");
	BAData d(*s,ctx);
	s.reset(0);
	
	if (mype == 0) printf("Data loaded\n");
	*/
	BADimensions dims(*s,ctx);
	BADimensionsSlacks dimsSlacks(dims);

	int nvar1real = dims.numFirstStageVars();
	//int ncons1 = dims.numFirstStageCons();
	int nvar2real = dims.numSecondStageVars(0);
	int ncons2 = dims.numSecondStageCons(0);

	BAFlagVector<variableState> states;
	
	states.allocate(dimsSlacks,ctx, PrimalVector);

	assert(nscenIn <= nscen);
	vector<double> firstStageSol;
	{
		scoped_ptr<rawInput> s2(new rawInput(datarootname, nscenIn));

		SolverInterface sub(*s2,ctx,SolverInterface::useDual);
		s2.reset(0);

		sub.loadStatus(basein);
		sub.setPrimalTolerance(1e-6);
		sub.setDualTolerance(1e-6);
		sub.go();

		firstStageSol = sub.getFirstStagePrimalColSolution();

		sub.getStates(states);
	}
	/*for (int i = 0; i < nvar1real; i++) {
		cout << firstStageSol[i] << endl;
	}
	return 0;*/

	// solve second-stage subproblems
	// store basis from first feasible 2nd stage subproblem
	vector<ClpSimplex::Status> saveState(nvar2real + ncons2);
	bool gotBasis = false;
	vector<int> wasInfeasible;
	const CoinPackedMatrix &Wmat = s->getSecondStageConstraints(0), &Tmat = s->getLinkingConstraints(0);
	assert(Wmat.getNumCols() == nvar2real);
	assert(Tmat.getNumCols() == nvar1real);
	for (int scen = nscenIn; scen < nscen; scen++) {
		if (!ctx.assignedScenario(scen)) continue;
		vector<double> Tx(ncons2);
		Tmat.times(&firstStageSol[0],&Tx[0]);
		const vector<double> &collb = s->getSecondStageColLB(scen),
			&colub = s->getSecondStageColUB(scen),
			&obj = s->getSecondStageObj(scen);
		vector<double> rowlb = s->getSecondStageRowLB(scen),
			rowub = s->getSecondStageRowUB(scen);

		for (int k = 0; k < ncons2; k++) {
			if (rowub[k] < 1e20) {
				rowub[k] -= Tx[k];
			}
			if (rowlb[k] >-1e20) {
				rowlb[k] -= Tx[k];
			}
		}
		
		ClpSimplex sub;
		sub.loadProblem(Wmat,&collb[0],&colub[0],&obj[0],&rowlb[0],&rowub[0]);
		sub.createStatus();
		if (gotBasis) {
			for (int i = 0; i < nvar2real; i++) 
				sub.setColumnStatus(i,saveState[i]);
			for (int i = 0; i < ncons2; i++)
				sub.setRowStatus(i,saveState[i+nvar2real]);
		}
		sub.setDualObjectiveLimit(1e6); // depends on problem!
		ClpSolve solvectl;
		solvectl.setPresolveType(ClpSolve::presolveOff);
		solvectl.setSolveType(ClpSolve::useDual);
		sub.initialSolve(solvectl);
		
		if (!sub.primalFeasible()) { wasInfeasible.push_back(scen); continue; }
		bool saveBasis = !gotBasis;
		for (int i = 0; i < nvar2real; i++) {
			states.getSecondStageVec(scen)[i] = convertClpStatus(sub.getColumnStatus(i));
			if (saveBasis) saveState[i] = sub.getColumnStatus(i);
		}
		for (int i = 0; i < ncons2; i++) {
			states.getSecondStageVec(scen)[i+nvar2real] = convertClpStatus(sub.getRowStatus(i));
			if (saveBasis) saveState[i+nvar2real] = sub.getRowStatus(i);
		}
		if (saveBasis) gotBasis = true;


	}

	struct intint { int gotBasis; int proc; } ii, ii2;
	ii.gotBasis = gotBasis;
	ii.proc = mype;
	MPI_Allreduce(&ii,&ii2,1,MPI_2INT,MPI_MAXLOC,MPI_COMM_WORLD);
	if (!ii2.gotBasis) {
		if (mype == 0) printf("Nobody got a feasible subproblem.\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	MPI_Bcast(&saveState[0],(nvar2real+ncons2)*sizeof(ClpSimplex::Status),MPI_BYTE,ii2.proc,MPI_COMM_WORLD);
	int nInfeasible = wasInfeasible.size();
	// infeasible subproblems can give very poorly conditioned starting bases,
	// so copy from one of the feasible subproblems
	for (int j = 0; j < nInfeasible; j++) {
		int scen = wasInfeasible[j];
		denseFlagVector<variableState> &s = states.getSecondStageVec(scen);
		for (int i = 0; i < nvar2real+ncons2; i++) {
			s[i] = convertClpStatus(saveState[i]);
		}
	}
	int nInfeasibleAll;
	MPI_Allreduce(&nInfeasible,&nInfeasibleAll,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	



	ctx.initializeAssignment(nscen);
	double t = MPI_Wtime();
	// no infeasible subproblems -> we have primal feasible basis
	// otherwise, primal and dual infeasible
	SolverInterface extensive(*s,ctx, (nInfeasibleAll == 0) ? SolverInterface::usePrimal : SolverInterface::useDual);

	s.reset(0);

	extensive.setStates(states);
	extensive.setPrimalTolerance(1e-6);
	extensive.setDualTolerance(1e-6);
	extensive.go();

	t = MPI_Wtime() - t;

	if (mype == 0) printf("Solve took %f seconds\n",t);
	// - indicates no output
	if (baseout != "-") {
		if (mype == 0) printf("Writing solution\n");
		extensive.writeStatus(baseout);
		if (mype == 0) printf("Finished writing solution\n");
	}


}

#endif
