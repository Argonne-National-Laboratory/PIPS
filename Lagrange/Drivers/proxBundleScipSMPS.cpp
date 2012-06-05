#include "SMPSInput.hpp"
#include "ScipLagrangeSolver.hpp"
#include "CbcRecourseSolver.hpp"
#include "ClpRecourseSolver.hpp"
#include "OOQPInterface.hpp"
#include "proxQPManager.hpp"
#include "QpGenSparseMa57.h"
#include "MehrotraSolver.h"

using namespace std;

// example:
// ./proxBundleSMPS /homes/mlubin/sslp/sslp_5_25_50
// compare with:
// ./conicBundleSMPS /homes/mlubin/sslp/sslp_5_25_50
// optimal value is -121.6

int main(int argc, char **argv) {

	
	MPI_Init(&argc, &argv);

	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc != 2) {
		if (mype == 0) printf("Usage: %s [SMPS root name]\n",argv[0]);
		return 1;
	}

	string smpsrootname(argv[1]);

	SMPSInput input(smpsrootname+".cor",smpsrootname+".tim",smpsrootname+".sto");
	
	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(input.nScenarios());

	proxQPManager<OOQPInterface<MehrotraSolver,QpGenSparseMa57>,ScipLagrangeSolver,ClpRecourseSolver> manager(input,ctx);

	while(!manager.terminated()) {
		manager.iterate();
	}

	MPI_Finalize();

	return 0;

}
