#include "stochasticInput.hpp"
#include <boost/scoped_ptr.hpp>
#include <cmath>
#include "ClpBALPInterface.hpp"
#include "ClpRecourseSolver.hpp"
#include "rawInput.hpp"

using namespace std;
using boost::scoped_ptr;


template <typename BALPSolver> void ucCoefficientDiving(stochasticInput &input, const std::string &LPBasis,
	const vector<int> &fixed,
	MPI_Comm comm = MPI_COMM_WORLD) {

	using namespace std;

	BAContext ctx(comm);
	ctx.initializeAssignment(input.nScenarios());

	int nvar1 = input.nFirstStageVars();
	
	BALPSolver solver(input, ctx, BALPSolver::useDual);
	solver.loadStatus(LPBasis);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	for (unsigned i = 0; i < fixed.size(); i++) {
		solver.setFirstStageColLB(fixed[i],1.0);
	}
	solver.go();


	int nIter = 0;
	while (true) {
		const vector<double>& firstSol = solver.getFirstStagePrimalColSolution();
		int atZero = 0, atOne = 0, fractional = 0;
		double bestDist = 1.;
		int bestIdx = -1;
		for (int i = 0; i < nvar1; i++) {
			if (fabs(firstSol[i]-0.0) < 1e-5) {
				atZero++;
			} else if (fabs(firstSol[i]-1.0) < 1e-5) {
				atOne++;
			} else {
				assert(solver.getFirstStageColState(i) == Basic);
				fractional++;
				if (fabs(firstSol[i]-1.0) < bestDist) {
					bestDist = fabs(firstSol[i]-1.0);
					bestIdx = i;
				}
			}
		}
		printf("Iteration %d\n",nIter++);
		printf("%d at zero, %d at one, %d fractional\n",atZero,atOne,fractional);
		if (fractional == 0) break;
		printf("Candidate for diving: %d %f\n",bestIdx,firstSol.at(bestIdx));
		
		solver.setFirstStageColLB(bestIdx,1.0);
		solver.go();
		assert(solver.getStatus() == Optimal);

	}
	printf("Objective from coefficient diving: %f\n", solver.getObjective());



}

template <typename BALPSolver, typename RecourseSolver> void ucRounding(stochasticInput &input, const std::string &LPBasis,
	const vector<int> &fixed,
	MPI_Comm comm = MPI_COMM_WORLD) {

	using namespace std;

	BAContext ctx(comm);
	ctx.initializeAssignment(input.nScenarios());

	int nvar1 = input.nFirstStageVars();
	
	BALPSolver solver(input, ctx, BALPSolver::useDual);
	solver.loadStatus(LPBasis);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);

	for (unsigned i = 0; i < fixed.size(); i++) {
		solver.setFirstStageColLB(fixed[i],1.0);
	}
	solver.go();

	vector<double> firstSol = solver.getFirstStagePrimalColSolution();
	for (int i = 0; i < nvar1; i++) {
		if (fabs(firstSol[i]-0.0) < 0.15) {
			firstSol[i] = 0.0;
		} else {
			firstSol[i] = 1.0;
		}
		/*
		if (fabs(firstSol[i]-0.0) < fabs(firstSol[i]-1.0)) {
			firstSol[i] = 0.0;
		} else {
			firstSol[i] = 1.0;
		}*/	
	}
	
	int nvar2 = input.nSecondStageVars(0);
	int nscen = input.nScenarios();
	int ncons2 = input.nSecondStageCons(0);
	
	const vector<double> &obj1 = input.getFirstStageObj();

	vector<variableState> rowSave(ncons2), colSave(nvar2);
	bool havesave = false;
	double sum = 0.0;
	bool infeas = false;
	for (int k = 0; k < nvar1; k++) sum += firstSol[k]*obj1[k];
	for (int scen = 0; scen < nscen; scen++) {
		RecourseSolver rsol(input, scen, firstSol);
		rsol.setDualObjectiveLimit(1e7);
		
		if (havesave) {
			for (int r = 0; r < nvar2; r++) {
				rsol.setSecondStageColState(r,colSave[r]);
			}
			for (int r = 0; r < ncons2; r++) {
				rsol.setSecondStageRowState(r,rowSave[r]);
			}
		}
		rsol.go();
		sum += rsol.getObjective();
		
		if (rsol.getStatus() == ProvenInfeasible) {
			printf("got infeasible 1st stage\n");
			sum = COIN_DBL_MAX;
			infeas = true; break;
		}
		assert(rsol.getStatus() == Optimal);
		if (!havesave) {
			for (int r = 0; r < nvar2; r++) {
				colSave[r] = rsol.getSecondStageColState(r);
			}
			for (int r = 0; r < ncons2; r++) {
				rowSave[r] = rsol.getSecondStageRowState(r);
			}
			havesave = true;
		}

	}


	printf("Rounding solution: %f\n", sum);

}


// improve the lower bound by fixing generators to on if they are infeasible off
template <typename BALPSolver> void ucFakeStrongBranching(stochasticInput &input, const std::string &LPBasis,
	MPI_Comm comm = MPI_COMM_WORLD) {

	using namespace std;

	BAContext ctx(comm);
	ctx.initializeAssignment(input.nScenarios());

	int nvar1 = input.nFirstStageVars();
	

	BALPSolver solver(input, ctx, BALPSolver::useDual);
	solver.loadStatus(LPBasis);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();

	solver.setDualObjectiveLimit(10.*solver.getObjective());

	BADimensions dims(input,ctx);
	BADimensionsSlacks dimsSlacks(dims);
	BAFlagVector<variableState> states;
	states.allocate(dimsSlacks,ctx, PrimalVector);
	solver.getStates(states);

	const vector<double>& firstSol = solver.getFirstStagePrimalColSolution();
	vector<int> fixed;
	for (int i = 0; i < nvar1; i++) {
		if (fabs(firstSol[i]-0.0) < 1e-5) {
		} else if (fabs(firstSol[i]-1.0) < 1e-5) {
		} else {
			solver.setFirstStageColUB(i,0.0);
			solver.go();
			
			solver.setFirstStageColUB(i,1.0);
			solver.setStates(states);
			if (solver.getStatus() == ProvenInfeasible) {
				solver.setFirstStageColLB(i,1.0);
				fixed.push_back(i);
				solver.go();
				solver.getStates(states);
				printf("New Lower Bound: %f\n", solver.getObjective());
			} else {
			}
		}
	}
	cout << fixed.size() << " fixed at UB\n";
	for (unsigned i = 0; i < fixed.size(); i++) {
		printf("%d,",fixed[i]);
	}
	printf("\n");

}

// 194 fixed
int fixed4h10s[] = {36,63,64,93,96,116,117,120,121,122,123,143,144,145,147,151,157,158,159,164,167,168,169,170,171,172,191,192,193,194,195,196,197,198,204,205,206,207,212,213,214,215,218,219,220,221,238,241,257,296,297,298,299,325,354,355,357,377,378,381,382,383,384,408,412,418,419,420,424,425,428,429,430,431,432,433,452,453,454,455,456,457,458,459,465,466,467,468,473,474,475,476,479,480,481,482,499,518,557,558,559,560,586,615,616,618,638,639,642,643,644,645,669,673,679,680,681,685,686,689,690,691,692,693,694,713,714,715,716,717,718,719,720,726,727,728,729,734,735,736,737,740,741,742,743,760,779,818,819,820,821,847,876,879,899,900,903,904,905,906,930,934,940,941,942,946,947,950,951,952,953,954,955,974,975,976,977,978,979,980,981,987,988,989,990,995,996,997,998,1001,1002,1003,1004,1040};

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
	scoped_ptr<fakeIntegerWrapper> s(new fakeIntegerWrapper(datarootname,nscen));

	//ucFakeStrongBranching<ClpBALPInterface>(*s,LPBasis);
	
	//ucCoefficientDiving<ClpBALPInterface>(*s,LPBasis, vector<int>(fixed4h10s,fixed4h10s+194));
	ucRounding<ClpBALPInterface,ClpRecourseSolver>(*s,LPBasis, vector<int>(fixed4h10s,fixed4h10s+194));

	MPI_Finalize();

	return 0;

}
