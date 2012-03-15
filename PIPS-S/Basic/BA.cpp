#include "BA.hpp"

#include <cassert>

#define MIN(a,b) ( (a>b) ? b : a )

BAContext::BAContext(MPI_Comm comm) : mpicomm(comm) {
	int flag;
	MPI_Initialized(&flag);
	assert(flag);

	MPI_Comm_size(comm, &_nprocs);
	MPI_Comm_rank(comm, &_mype);

}

// Assign scenarios monotonically to processes, 
// putting odd scenarios at the end.
// Note that other code assumes this assignment,
// Arbitrary assignments WILL NOT work.
// Must relabel scenarios to achieve this.
void BAContext::initializeAssignment(int nscen) {
	if (!assignedScen.size()) { // this check is useful in balpAdvancedAndDumpBasis
		assignedScen.resize(nscen);
		nper = nscen / _nprocs;
		for (int i = 0; i < nscen; i++) {
			assignedScen[i] = MIN(i/nper, _nprocs-1);
		}
	}
	localScen.clear();
	localScen.push_back(-1);
	for (int i = 0; i < nscen; i++) {
		if (assignedScen[i] == _mype) localScen.push_back(i);
	}	

}

double BAContext::reduce(double d) const {
	double out;
	MPI_Allreduce(&d,&out,1,MPI_DOUBLE,MPI_SUM,mpicomm);
	return out;
}

int BAContext::reduce(int d) const {
	int out;
	MPI_Allreduce(&d,&out,1,MPI_INT,MPI_SUM,mpicomm);
	return out;
}


int BAContext::owner(int scen) const {
	if (scen == -1) return 0;
	else return (assignedScen.at(scen));
}


BADimensions::BADimensions(stochasticInput &in, const BAContext &ctx) {
	
	same = in.scenarioDimensionsEqual();
	nScen = in.nScenarios();
	nFirstStageVars_ = in.nFirstStageVars();
	nFirstStageCons_ = in.nFirstStageCons();

	const std::vector<int> &localScen = ctx.localScenarios();

	if (same) {
		nSecondStageVars_ = in.nSecondStageVars(localScen.at(1));
		nSecondStageCons_ = in.nSecondStageCons(localScen.at(1));
	} else {
		nSecondStageVars_v.resize(nScen);
		nSecondStageCons_v.resize(nScen);
		// we don't know dimensions of scenarios not assigned to us
		for (unsigned k = 1; k < localScen.size(); k++) {
			int i = localScen[k];
			nSecondStageVars_v[i] = in.nSecondStageVars(i);
			nSecondStageCons_v[i] = in.nSecondStageCons(i);
		}
	}

}

