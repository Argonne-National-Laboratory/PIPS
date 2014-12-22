#include "SMPSInput.hpp"
#include "BAData.hpp"
#include "PIPSSInterface.hpp"
#include <boost/scoped_ptr.hpp>
#include <cstdlib>

using boost::scoped_ptr; // replace with unique_ptr for C++11
using namespace std;

int main(int argc, char **argv) {

	MPI_Init(&argc, &argv);

	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc < 2) {
	  if (mype == 0) printf("Usage: %s [SMPS root name] [solution output root name] [starting basis root name]\n",argv[0]);
	  return 1;
	}

	PIPSLogging::init_logging(1);

	string smpsrootname(argv[1]);

	SMPSInput input(smpsrootname+".cor",smpsrootname+".tim",smpsrootname+".sto");

	//scoped_ptr<SMPSInput> s(new SMPSInput(datarootname,nscen));

	BAContext ctx(MPI_COMM_WORLD);

	PIPSSInterface solver(input, ctx, PIPSSInterface::useDual);
	//s.reset(0); // free memory
	if (argc == 4) {
		solver.loadStatus(argv[4]);
	}
	//solver.setDumpFrequency(5000,argv[3]);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();


	if (argc >= 3 && argv[2][0] != '-') {
	  if (mype == 0) PIPS_APP_LOG_SEV(info) << "Writing solution";
	  solver.writeStatus(argv[2]);
	  if (mype == 0) PIPS_APP_LOG_SEV(info) << "Finished writing solution";
	}

	MPI_Finalize();

	return 0;
}

