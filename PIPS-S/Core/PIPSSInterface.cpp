#include "PIPSSInterface.hpp"
#include "BALPSolverDual.hpp"
#include "BALPSolverPrimal.hpp"

PIPSSInterface::PIPSSInterface(stochasticInput &in, BAContext &ctx, solveType t) : d(in,ctx) {

	if (t == usePrimal) {
		solver = new BALPSolverPrimal(d);
	} else {
		solver = new BALPSolverDual(d);
	}
	if (d.ctx.mype() == 0) {
		printf("First stage: %d cons %d vars\n", d.dims.numFirstStageCons(), d.dims.inner.numFirstStageVars());
		printf("Second stage: %d cons %d vars %d scenarios\n", d.dims.numSecondStageCons(0), d.dims.inner.numSecondStageVars(0),
			d.dims.numScenarios());
		printf("reinvert every %d\n",solver->reinvertFrequency);
	}

}

PIPSSInterface::PIPSSInterface(const BAData& d, solveType t) : d(d) {

	if (t == usePrimal) {
		solver = new BALPSolverPrimal(d);
	} else {
		solver = new BALPSolverDual(d);
	}
}


PIPSSInterface::~PIPSSInterface() {
	delete solver;
}


void PIPSSInterface::go() {
	double t = MPI_Wtime();
	int mype = d.ctx.mype();

	solver->go();

	// clean up infeasibilities
	int count = 0;
	while (solver->getStatus() != Optimal && count < 10) {
		if (solver->getStatus() == ProvenInfeasible ||
		    solver->getStatus() == ProvenUnbounded) break;
		BALPSolverBase *solver2;
		if (solver->getStatus() == DualFeasible) {
			solver2 = new BALPSolverDual(d);
			if (mype == 0) printf("Switching to dual\n");
		} else {
			assert(solver->getStatus() == PrimalFeasible);
			solver2 = new BALPSolverPrimal(d);
			if (mype == 0) printf("Switching to primal\n");
		}
		solver2->setStates(solver->getStates());
		solver2->setPrimalTolerance(solver->getPrimalTolerance());
		solver2->setDualTolerance(solver->getDualTolerance());
		solver2->setReinversionFrequency(5); // be very careful

		solver2->ftranTime = solver->ftranTime;
		solver2->btranTime = solver->btranTime;
		solver2->updateIteratesTime = solver->updateIteratesTime;
		solver2->ftranDSETime = solver->ftranDSETime;
		solver2->invertTime = solver->invertTime;
		solver2->selectEnteringTime = solver->selectEnteringTime;
		solver2->selectLeavingTime = solver->selectLeavingTime;
		solver2->priceTime = solver->priceTime;
		solver2->updateColumnTime = solver->updateColumnTime;
		solver2->startTime = solver->startTime;
		
		solver2->replaceFirst = solver->replaceFirst;
		solver2->firstReplaceSecond = solver->firstReplaceSecond;
		solver2->secondReplaceFirst = solver->secondReplaceFirst;
		solver2->replaceSecondSelf = solver->replaceSecondSelf;
		solver2->replaceSecondOther = solver->replaceSecondOther;

		solver2->nIter = solver->nIter;

		delete solver; // free memory
		solver2->go();
		solver = solver2;
	}

 	if (solver->getStatus() != Optimal) {
		if (mype == 0) printf("Switched between primal and dual %d times, and still not optimal!\n", count);
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	t = MPI_Wtime() - t;

	if (mype == 0) printf("Solve took %f seconds\n",t);


}

std::vector<double> PIPSSInterface::getFirstStagePrimalColSolution() const {
	const denseVector &x = solver->getPrimalSolution().getFirstStageVec();
	int nvar1real = d.dims.inner.numFirstStageVars();
	return std::vector<double>(&x[0],&x[nvar1real]);
}

void PIPSSInterface::setFirstStageColState(int idx,variableState s) {
	solver->states.getFirstStageVec()[idx] = s;
}

void PIPSSInterface::setFirstStageRowState(int idx,variableState s) {
	int nvarreal = d.dims.inner.numFirstStageVars();
	solver->states.getFirstStageVec()[idx+nvarreal] = s;
}

void PIPSSInterface::setSecondStageColState(int scen, int idx,variableState s) {
	solver->states.getSecondStageVec(scen)[idx] = s;
}

void PIPSSInterface::setSecondStageRowState(int scen, int idx,variableState s) {
	int nvarreal = d.dims.inner.numSecondStageVars(scen);
	solver->states.getSecondStageVec(scen)[idx+nvarreal] = s;
}
void PIPSSInterface::commitStates() {
	solver->status = LoadedFromFile;
	solver->setupIndices();
}

variableState PIPSSInterface::getFirstStageColState(int idx) const {
	return solver->states.getFirstStageVec()[idx];
}

variableState PIPSSInterface::getFirstStageRowState(int idx) const {
	int nvarreal = d.dims.inner.numFirstStageVars();
	return solver->states.getFirstStageVec()[nvarreal+idx];
}

variableState PIPSSInterface::getSecondStageColState(int scen, int idx) const {
	return solver->states.getSecondStageVec(scen)[idx];
}

variableState PIPSSInterface::getSecondStageRowState(int scen, int idx) const {
	int nvarreal = d.dims.inner.numSecondStageVars(scen);
	return solver->states.getSecondStageVec(scen)[nvarreal+idx];
}	
