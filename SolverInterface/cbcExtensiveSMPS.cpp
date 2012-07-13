#include <cmath>
#include "CbcBALPInterface.hpp"
#include "SMPSInput.hpp"

using namespace std;


// solve extensive form with Cbc
// assumes UC model

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

	CbcBALPInterface solver(input, ctx);
	solver.go();

	MPI_Finalize();

}

