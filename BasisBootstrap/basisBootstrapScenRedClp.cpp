#include "basisBootstrapScenRedDriver.hpp"
#include "ClpRecourseSolver.hpp"
#include "ClpBALPInterface.hpp"
#include "PIPSSInterface.hpp"
#include "rawInput.hpp"

int main(int argc, char **argv) {
	MPI_Init(&argc,&argv);

	if (argc != 5) {
		printf("Usage: %s [raw problem data] [total scenarios] [initial scenario] [added per iteration]\n",argv[0]);
		return 1;
	}

	std::string basename(argv[1]);
	int nscen = atoi(argv[2]);
	int initial = atoi(argv[3]);
	int nper = atoi(argv[4]);

	rawInput input(basename,nscen);
	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(nscen);
	basisBootstrapScenRedDriver<ClpBALPInterface,ClpRecourseSolver> s(input,ctx);
	//simpleBasisBootstrapDriver<PIPSSInterface,ClpRecourseSolver> s(input,ctx);

	s.goFromScratch(initial,nper);



	MPI_Finalize();
}
