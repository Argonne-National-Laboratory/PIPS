#include "ClpRecourseSolver.hpp"
#include "rawInput.hpp"
#include "ucRootNode.hpp"

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
		if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [fractional solution file]\n",argv[0]);
		return 1;
	}

	string datarootname(argv[1]);
	int nscen = atoi(argv[2]);

	fakeIntegerWrapper s(datarootname,nscen);
	int nvar1 = s.nFirstStageVars();
	vector<double> sol(nvar1);

	{
		ifstream f(argv[3]);
		for (int i = 0; i < nvar1; i++) {
			f >> sol[i];
		}
	}
	
	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(nscen);
	double obj = roundSolution<ClpRecourseSolver>(s,ctx,sol,0.01);

	printf("Primal objective: %g\n",obj);

	MPI_Finalize();
	return 0;
}
