#include "ucRollingModel.hpp"
#include "CbcBALPInterface.hpp"
#include "ClpBALPInterface.hpp"

int main(int argc, char **argv) {
	MPI_Init(&argc,&argv);

	int nscen = 3;
	ucRollingModel model("/sandbox/vzavala/Illinois",nscen,0,3);

	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(nscen);

	ClpBALPInterface cbc(model, ctx);

	cbc.go();

	MPI_Finalize();

}
