#include "ucRollingModel.hpp"
#include "ScipLagrangeSolver.hpp"
#include "CbcRecourseSolver.hpp"
#include "ClpRecourseSolver.hpp"
#include "OOQPInterface.hpp"
#include "proxQPManager.hpp"
#include "QpGenSparseMa57.h"
#include "MehrotraSolver.h"

using namespace std;

int main(int argc, char **argv) {
	MPI_Init(&argc,&argv);

	if (argc != 3) {
		printf("Usage: %s [num scenarios] [time horizon]\n",argv[0]);
		return 1;
	}

	int nscen = atoi(argv[1]);
	int T = atoi(argv[2]);
	ucRollingModel input("../../../apps/unitcommitment_rolling/Illinois",nscen,0,T);

	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(input.nScenarios());

	proxQPManager<OOQPInterface<MehrotraSolver,QpGenSparseMa57>,ScipLagrangeSolver,ClpRecourseSolver> manager(input,ctx);

	while(!manager.terminated()) {
		manager.iterate();
	}

	MPI_Finalize();

}
