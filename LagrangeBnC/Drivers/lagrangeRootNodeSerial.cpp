
#include "rawInput.hpp"
#include "ClpBALPInterface.hpp"
#include "CbcLagrangeSolver.hpp"
#include "ClpRecourseSolver.hpp"
#include "lagrangeRootNode.hpp"
#include <boost/scoped_ptr.hpp>

using namespace std;
using boost::scoped_ptr;


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
		if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [starting basis root name]\n",argv[0]);
		return 1;
	}

	string datarootname(argv[1]);
	int nscen = atoi(argv[2]);
	string lpbasisname(argv[3]);

	if (mype == 0) printf("Initializing data interface\n");
	scoped_ptr<fakeIntegerWrapper> s(new fakeIntegerWrapper(datarootname,nscen));

	lagrangeRootNode<ClpBALPInterface,CbcLagrangeSolver,ClpRecourseSolver>(*s,argv[1]);

	MPI_Finalize();

	return 0;

}
