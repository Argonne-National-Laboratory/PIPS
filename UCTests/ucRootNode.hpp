#ifndef UCROOTNODE_HPP
#define UCROOTNODE_HPP


#include "OsiSubproblemWrapper.hpp"
#include "stochasticInput.hpp"
#include <cmath>
#include "CglMixedIntegerRounding2.hpp"
#include "CglProbing.hpp"
#include "ClpRecourseSolver.hpp"


template<typename BALPSolver> void ucRootNode(stochasticInput &input, 
	std::string const& LPBasis, MPI_Comm comm = MPI_COMM_WORLD) {

	using namespace std;

	BAContext ctx(comm);
	ctx.initializeAssignment(input.nScenarios());
	
	int mype;
	MPI_Comm_rank(comm,&mype);

	vector<int> const &localScen = ctx.localScenarios();

	int nvar1 = input.nFirstStageVars();
	
	BALPSolver solver(input, ctx, BALPSolver::useDual);
	solver.loadStatus(LPBasis);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();

	double lprelax = solver.getObjective();

	vector<int> fixstg1(nvar1,0), fixstg1_all(nvar1);

	for(unsigned i = 1; i < localScen.size(); i++) {
		//CglMixedIntegerRounding2 round;
		int scen = localScen[i];
		CglProbing probing;
		probing.setUsingObjective(-1);
		OsiCuts cuts;
		OsiSubproblemWrapper wrap(input,scen);
		wrap.setCurrentSolution(solver.getFirstStagePrimalColSolution(),solver.getSecondStagePrimalColSolution(scen));
		wrap.setCurrentReducedCosts(solver.getFirstStageDualColSolution(),solver.getSecondStageDualColSolution(scen));
		//round.generateCuts(wrap,cuts);
		probing.generateCuts(wrap,cuts);
		
		//printf("%d column cuts, %d row cuts\n",cuts.sizeColCuts(),cuts.sizeRowCuts());

		for (int i = 0; i < cuts.sizeColCuts(); i++) {
			const OsiColCut& cut = cuts.colCut(i);
			
			const CoinPackedVector &lbs = cut.lbs();
			{
				int nnz = lbs.getNumElements();
				const int* idx = lbs.getIndices();
				const double* elts = lbs.getElements();
				for (int r = 0; r < nnz; r++) {
					if (idx[r] >= nvar1) {
						printf("got second-stage column cut\n");
					} else {
						assert(fabs(elts[r]-1.0) < 1e-7);
						fixstg1[idx[r]] = 1;
					}
				}
			}
		}
		/*
		int nvar2 = input.nSecondStageVars(scen);
		for (int i = 0; i < cuts.sizeRowCuts(); i++) {
			vector<double> elts1(nvar1,0.0), elts2(nvar2,0.0);
			const OsiRowCut& cut = cuts.rowCut(i);
			const CoinPackedVector &v = cut.row();
			
			int nnz = v.getNumElements();
			const int* idx = v.getIndices();
			const double* elts = v.getElements();
			for (int r = 0; r < nnz; r++) {
				if (idx[r] >= nvar1) {
					elts2[idx[r]-nvar1] = elts[r];
				} else {
					elts1[idx[r]] = elts[r];
				}
			}
			solver.addRow(elts1,elts2,scen,cut.lb(),cut.ub());
		}*/
	}
	//solver.commitNewRows();

	MPI_Allreduce(&fixstg1[0],&fixstg1_all[0],nvar1,MPI_INT,MPI_SUM,comm);
	for (int i = 0; i < nvar1; i++) {
		if (fixstg1_all[i]) solver.setFirstStageColLB(i,1.0);
	}
	solver.commitStates();
	solver.go();

	double lpcuts = solver.getObjective();

	
	vector<double> firstSol = solver.getFirstStagePrimalColSolution();
	for (int i = 0; i < nvar1; i++) {
		if (firstSol[i] < 0.01) {
			firstSol[i] = 0.0;
		} else {
			firstSol[i] = 1.0;
		}
	}

	int nvar2 = input.nSecondStageVars(0);
	int ncons2 = input.nSecondStageCons(0);
	
	const vector<double> &obj1 = input.getFirstStageObj();

	vector<variableState> rowSave(ncons2), colSave(nvar2);
	bool havesave = false;
	double sum = 0.0;
	bool infeas = false;
	typedef ClpRecourseSolver RecourseSolver;
	for (unsigned i = 1; i < localScen.size(); i++) {
		int scen = localScen[i];
		RecourseSolver rsol(input, scen, firstSol);
		rsol.setDualObjectiveLimit(1e7);
		
		if (havesave) {
			for (int r = 0; r < nvar2; r++) {
				rsol.setSecondStageColState(r,colSave[r]);
			}
			for (int r = 0; r < ncons2; r++) {
				rsol.setSecondStageRowState(r,rowSave[r]);
			}
		}
		rsol.go();
		sum += rsol.getObjective();
		
		if (rsol.getStatus() == ProvenInfeasible) {
			printf("got infeasible 1st stage\n");
			sum = COIN_DBL_MAX;
			infeas = true; break;
		}
		assert(rsol.getStatus() == Optimal);
		if (!havesave) {
			for (int r = 0; r < nvar2; r++) {
				colSave[r] = rsol.getSecondStageColState(r);
			}
			for (int r = 0; r < ncons2; r++) {
				rowSave[r] = rsol.getSecondStageRowState(r);
			}
			havesave = true;
		}

	}
	double allsum;
	MPI_Allreduce(&sum,&allsum,1,MPI_DOUBLE,MPI_SUM,comm);
	for (int k = 0; k < nvar1; k++) allsum += firstSol[k]*obj1[k];


	if (mype == 0) printf("LP Relaxation: %f\n",lprelax);
	if (mype == 0) printf("LP Relaxation + Probing Cuts: %f\n",lpcuts);
	if (mype == 0) printf("Rounding solution: %f\n", allsum);
	if (mype == 0) printf("Gap: %f%%\n",100*(allsum-lpcuts)/lpcuts);



}


#endif
