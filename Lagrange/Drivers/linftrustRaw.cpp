#include "rawInput.hpp"
#include "CbcLagrangeSolver.hpp"
#include "CbcRecourseSolver.hpp"
#include "ClpRecourseSolver.hpp"
#include "ClpBALPInterface.hpp"
#include "PIPSSInterface.hpp"
#include "proxLinfTrustManager.hpp"

using namespace std;

int main(int argc, char **argv) {

	
	MPI_Init(&argc, &argv);

	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc != 3) {
		if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios]\n",argv[0]);
		return 1;
	}
	string datarootname(argv[1]);
	int nscen = atoi(argv[2]);

	rawInput input(datarootname,nscen);
	
	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(input.nScenarios());

	proxLinfTrustManager<ClpBALPInterface,CbcLagrangeSolver,ClpRecourseSolver> manager(input,ctx);

	while(!manager.terminated()) {
		manager.iterate();
	}

	MPI_Finalize();

	return 0;

}
