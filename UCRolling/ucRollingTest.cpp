#include "ucRollingModel.hpp"


int main(int argc, char **argv) {
	MPI_Init(&argc,&argv);

	ucRollingModel model("/sandbox/vzavala/Illinois",10,0,4);



	MPI_Finalize();

}
