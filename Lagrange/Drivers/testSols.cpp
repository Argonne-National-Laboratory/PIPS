#include "rawInput.hpp"
#include "CbcLagrangeSolver.hpp"
#include "CbcRecourseSolver.hpp"
#include "ClpRecourseSolver.hpp"
#include "ClpBALPInterface.hpp"
#include "PIPSSInterface.hpp"
#include "bundleManager.hpp"

using namespace std;

class solutionTester : public bundleManager<PIPSSInterface,CbcLagrangeSolver,ClpRecourseSolver> 
{
public:
	solutionTester(stochasticInput &input, BAContext & ctx, string const &solbase, int niter) :
		bundleManager<PIPSSInterface,CbcLagrangeSolver,ClpRecourseSolver>(input,ctx),
		solbase(solbase), niter(niter) {}

	void go() {
		int nscen = input.nScenarios();
		int nvar1 = input.nFirstStageVars();
		double obj = COIN_DBL_MAX;
		for (int i = 1; i <= niter; i++) {
			stringstream ss;
			ss << solbase << i;
			ifstream f(ss.str().c_str());
			for (int s = 0; s < nscen; s++) {
				vector<double> sol(nvar1);
				for (int k = 0; k < nvar1; k++) {
					f >> sol[k];
				}
				obj = min(testPrimal(sol),obj);
			}
			if (ctx.mype() == 0) printf("Iter %d Best Solution %f\n",i,obj);
		}
	}

	void doStep() {}
	


private:
	string solbase;
	int niter;

};

int main(int argc, char **argv) {

	
	MPI_Init(&argc, &argv);

	int mype;
	MPI_Comm_rank(MPI_COMM_WORLD,&mype);

	if (argc != 5) {
		if (mype == 0) printf("Usage: %s [rawdump root name] [num scenarios] [solution filename base] [number of iterations]\n",argv[0]);
		return 1;
	}
	string datarootname(argv[1]);
	int nscen = atoi(argv[2]);

	rawInput input(datarootname,nscen);
	
	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(input.nScenarios());

	solutionTester s(input, ctx, argv[3], atoi(argv[4]));

	s.go();

	MPI_Finalize();

}
