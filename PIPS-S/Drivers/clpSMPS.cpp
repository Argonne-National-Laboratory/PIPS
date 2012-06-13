#include "SMPSInput.hpp"
#include "ClpBALPInterface.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>

using boost::scoped_ptr; // replace with unique_ptr for C++11
using namespace std;

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

	ClpBALPInterface solver(input, ctx, ClpBALPInterface::useDual);
	if (argc == 5) {
		solver.loadStatus(argv[4]);
	}
	//solver.setDumpFrequency(5000,argv[3]);	
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();


	if (argc >= 4 && argv[3][0] != '-') {
		if (mype == 0) printf("Writing solution\n");
		solver.writeStatus(argv[3]);
		if (mype == 0) printf("Finished writing solution\n");
	}

	MPI_Finalize();

	return 0;
}

