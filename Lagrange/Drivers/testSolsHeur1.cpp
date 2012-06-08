#include "rawInput.hpp"
#include "ScipLagrangeSolver.hpp"
#include "CbcRecourseSolver.hpp"
#include "ClpRecourseSolver.hpp"
#include "ClpBALPInterface.hpp"
#include "bundleManager.hpp"

using namespace std;

class fakeIntegerWrapper : public rawInput {
public:
	fakeIntegerWrapper(const std::string &datarootname, int overrideScenarioNumber = 0, MPI_Comm comm = MPI_COMM_WORLD) :
		rawInput(datarootname, overrideScenarioNumber, comm) {}
	
	virtual bool isFirstStageColInteger(int col) { return true; }

};

// heuristic for uc model:
// given a candidate solution, solve IP subproblem for each scenario where
// units are fixed to ON if they are ON in the candidate solution.
// this gives us a feasible solution for each scenario. Take the logical union/bitwise AND of the
// solutions to get a feasible solution for the extensive form IP
class solutionTester : public bundleManager<ClpBALPInterface,ScipLagrangeSolver,ClpRecourseSolver> 
{
public:
	solutionTester(stochasticInput &input, BAContext & ctx, string const &solbase, int niter) :
		bundleManager<ClpBALPInterface,ScipLagrangeSolver,ClpRecourseSolver>(input,ctx),
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
				vector<double> sol(nvar1), at(nvar1);
				for (int k = 0; k < nvar1; k++) {
					f >> sol[k];
					f >> at[k];
					if (fabs(sol[k]) < 1e-5) sol[k] = 0.;
					else {
						assert(fabs(sol[k]-1.) < 1e-5);
						sol[k] = 1.;
					}
				}
				obj = min(primalHeur(sol,at),obj);
			}
			if (ctx.mype() == 0) printf("Iter %d Best Solution %f\n",i,obj);
		}
	}

	void doStep() {}

	double primalHeur(vector<double> &sol, vector<double> const &at) {
		int nscen = input.nScenarios();
		int nvar1 = input.nFirstStageVars();
		for (int s = 0; s < nscen; s++) {
			ScipLagrangeSolver solver(input,s,at);
			for (int k = 0; k < nvar1; k++) {
				if (sol[k] == 1.) {
					solver.setFirstStageColLB(k,1.);
				}
			}
			solver.go();
			vector<double> solThis = solver.getBestFirstStageSolution();
			int cnt = 0;
			for (int k = 0; k < nvar1; k++) {
				if (fabs(solThis[k]-1.) < 1e-7) {
					if (sol[k] != 1.) cnt++;
					sol[k] = 1.;
				}
			}
			printf("set %d new to 1\n",cnt);
		}
		double newobj = testPrimal(sol);
		printf("Got feasible objective val %g\n",newobj);
		return newobj;

	}


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

	fakeIntegerWrapper input(datarootname,nscen);
	
	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(input.nScenarios());

	solutionTester s(input, ctx, argv[3], atoi(argv[4]));

	s.go();

	MPI_Finalize();

}
