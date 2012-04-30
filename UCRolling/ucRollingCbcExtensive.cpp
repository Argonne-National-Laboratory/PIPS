#include "ucRollingModel.hpp"
#include "CbcBALPInterface.hpp"

using namespace std;

int main(int argc, char **argv) {
	MPI_Init(&argc,&argv);

	if (argc != 3) {
		printf("Usage: %s [num scenarios] [time horizon]\n",argv[0]);
		return 1;
	}

	int nscen = atoi(argv[1]);
	int T = atoi(argv[2]);
	ucRollingModel model("../../../apps/unitcommitment_rolling/Illinois",nscen,0,T);

	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(nscen);

	CbcBALPInterface cbc(model, ctx);

	
	cbc.go();

	MPI_Finalize();

}
