#include "bootstrapDriver.hpp"
#include "ClpBALPInterface.hpp"

using namespace std;

int main(int argc, char **argv) {
	
	MPI_Init(&argc,&argv);
	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc != 6) {
		if (mype == 0) printf("Usage: %s [raw data input basename] [basis input basename] [basis output basename] [num scenarios in] [num scenarios out]\n", argv[0]);
		return 1;
	}

	string datarootname(argv[1]);
	string basein(argv[2]);
	string baseout(argv[3]);
	int nscenIn = atoi(argv[4]);
	int nscen = atoi(argv[5]);

	bootstrapDriver<ClpBALPInterface>(datarootname, basein, baseout, nscenIn, nscen);

	MPI_Finalize();

	return 0;
}	
