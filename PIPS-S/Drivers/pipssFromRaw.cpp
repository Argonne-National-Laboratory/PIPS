#include "BAData.hpp"
#include "rawInput.hpp"
#include "PIPSSInterface.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>

using boost::scoped_ptr; // replace with unique_ptr
using namespace std;

int main(int argc, char **argv) {
	
	MPI_Init(&argc, &argv);

	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc < 3) {
		if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [solution output root name] [starting basis root name]\n",argv[0]);
		return 1;
	}

	string datarootname(argv[1]);
	int nscen = atoi(argv[2]);

	scoped_ptr<rawInput> s(new rawInput(datarootname,nscen));
	
	BAContext ctx(MPI_COMM_WORLD);

	PIPSSInterface solver(*s, ctx, PIPSSInterface::useDual);
	s.reset(0); // free memory
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

