#include "ucRootNode.hpp"
#include "PIPSSInterface.hpp"
#include "rawInput.hpp"

using namespace std;

class fakeIntegerWrapper : public rawInput {
public:
	fakeIntegerWrapper(const std::string &datarootname, int overrideScenarioNumber = 0, MPI_Comm comm = MPI_COMM_WORLD) :
		rawInput(datarootname, overrideScenarioNumber, comm) {}
	
	virtual bool isFirstStageColInteger(int col) { return true; }

};


int main(int argc, char **argv) {

	
	MPI_Init(&argc, &argv);

	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc != 4) {
		if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [LP relaxation basis]\n",argv[0]);
		return 1;
	}

	string datarootname(argv[1]);
	int nscen = atoi(argv[2]);
	string LPBasis(argv[3]);

	if (mype == 0) printf("Initializing data interface\n");
	fakeIntegerWrapper s(datarootname,nscen);

	ucRootNode<PIPSSInterface>(s,LPBasis);


	MPI_Finalize();
	return 0;
}
