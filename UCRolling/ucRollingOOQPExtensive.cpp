#include "ucRollingModel.hpp"
#include "OOQPInterface.hpp"
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
	ucRollingModel model("../../../apps/unitcommitment_rolling/Illinois",nscen,0,T);

	OOQPInterface<MehrotraSolver,QpGenSparseMa57> ooqp(model);
	
	ooqp.go();

	MPI_Finalize();

}
