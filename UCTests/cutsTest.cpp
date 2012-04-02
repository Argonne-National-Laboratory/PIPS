#include "OsiSubproblemWrapper.hpp"
#include "stochasticInput.hpp"
#include <boost/scoped_ptr.hpp>
#include <cmath>
#include "ClpBALPInterface.hpp"
#include "ClpRecourseSolver.hpp"
#include "rawInput.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglProbing.hpp"
#include "CbcLagrangeSolver.hpp"

using namespace std;
using boost::scoped_ptr;

class fakeIntegerWrapper : public rawInput {
public:
	fakeIntegerWrapper(const std::string &datarootname, int overrideScenarioNumber = 0, MPI_Comm comm = MPI_COMM_WORLD) :
		rawInput(datarootname, overrideScenarioNumber, comm) {}
	
	virtual bool isFirstStageColInteger(int col) { return true; }

};

int fixed4h10s[] = {36,63,64,93,96,116,117,120,121,122,123,143,144,145,147,151,157,158,159,164,167,168,169,170,171,172,191,192,193,194,195,196,197,198,204,205,206,207,212,213,214,215,218,219,220,221,238,241,257,296,297,298,299,325,354,355,357,377,378,381,382,383,384,408,412,418,419,420,424,425,428,429,430,431,432,433,452,453,454,455,456,457,458,459,465,466,467,468,473,474,475,476,479,480,481,482,499,518,557,558,559,560,586,615,616,618,638,639,642,643,644,645,669,673,679,680,681,685,686,689,690,691,692,693,694,713,714,715,716,717,718,719,720,726,727,728,729,734,735,736,737,740,741,742,743,760,779,818,819,820,821,847,876,879,899,900,903,904,905,906,930,934,940,941,942,946,947,950,951,952,953,954,955,974,975,976,977,978,979,980,981,987,988,989,990,995,996,997,998,1001,1002,1003,1004,1040};

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

	BAContext ctx(MPI_COMM_WORLD);
	ctx.initializeAssignment(s->nScenarios());

	int nvar1 = s->nFirstStageVars();
	
	ClpBALPInterface solver(*s, ctx, ClpBALPInterface::useDual);
	solver.loadStatus(LPBasis);
	solver.setPrimalTolerance(1e-6);
	solver.setDualTolerance(1e-6);
	solver.go();


	printf("LP Relaxation LB: %f\n",solver.getObjective());
	
	vector<double> fixed(fixed4h10s,fixed4h10s+194);
	for (unsigned i = 0; i < fixed.size(); i++) {
		solver.setFirstStageColLB(fixed[i],1.0);
	}
	solver.go();
	

	for(int scen = 0; scen < s->nScenarios(); scen++) {
		CglMixedIntegerRounding2 round;
		CglProbing probing;
		probing.setUsingObjective(-1);
		OsiCuts cuts;
		OsiSubproblemWrapper wrap(*s,scen);
		wrap.setCurrentSolution(solver.getFirstStagePrimalColSolution(),solver.getSecondStagePrimalColSolution(scen));
		wrap.setCurrentReducedCosts(solver.getFirstStageDualColSolution(),solver.getSecondStageDualColSolution(scen));
		//round.generateCuts(wrap,cuts);
		probing.generateCuts(wrap,cuts);
		
		printf("%d column cuts, %d row cuts\n",cuts.sizeColCuts(),cuts.sizeRowCuts());

		int nvar2 = s->nSecondStageVars(scen);

		for (int i = 0; i < cuts.sizeRowCuts(); i++) {
			vector<double> elts1(nvar1,0.0), elts2(nvar2,0.0);
			const OsiRowCut& cut = cuts.rowCut(i);
			const CoinPackedVector &v = cut.row();
			
			int nnz = v.getNumElements();
			const int* idx = v.getIndices();
			const double* elts = v.getElements();
			for (int r = 0; r < nnz; r++) {
				if (idx[r] >= nvar1) {
					elts2[idx[r]-nvar1] = elts[r];
				} else {
					elts1[idx[r]] = elts[r];
				}
			}
			solver.addRow(elts1,elts2,scen,cut.lb(),cut.ub());
		}

		for (int i = 0; i < cuts.sizeColCuts(); i++) {
			const OsiColCut& cut = cuts.colCut(i);
			
			const CoinPackedVector &lbs = cut.lbs();
			{
				int nnz = lbs.getNumElements();
				const int* idx = lbs.getIndices();
				const double* elts = lbs.getElements();
				for (int r = 0; r < nnz; r++) {
					if (idx[r] >= nvar1) {
						printf("got second-stage column cut\n");
					} else {
						solver.setFirstStageColLB(idx[r],elts[r]);
					}
				}
			}
			
		}
	}
	solver.go();

	for(int scen = 0; scen < s->nScenarios(); scen++) {
		CglMixedIntegerRounding2 round;
		CglProbing probing;
		probing.setUsingObjective(-1);
		OsiCuts cuts;
		OsiSubproblemWrapper wrap(*s,scen);
		wrap.setCurrentSolution(solver.getFirstStagePrimalColSolution(),solver.getSecondStagePrimalColSolution(scen));
		wrap.setCurrentReducedCosts(solver.getFirstStageDualColSolution(),solver.getSecondStageDualColSolution(scen));
		round.generateCuts(wrap,cuts);
		//probing.generateCuts(wrap,cuts);
		
		printf("%d column cuts, %d row cuts\n",cuts.sizeColCuts(),cuts.sizeRowCuts());

		int nvar2 = s->nSecondStageVars(scen);

		for (int i = 0; i < cuts.sizeRowCuts(); i++) {
			vector<double> elts1(nvar1,0.0), elts2(nvar2,0.0);
			const OsiRowCut& cut = cuts.rowCut(i);
			const CoinPackedVector &v = cut.row();
			
			int nnz = v.getNumElements();
			const int* idx = v.getIndices();
			const double* elts = v.getElements();
			for (int r = 0; r < nnz; r++) {
				if (idx[r] >= nvar1) {
					elts2[idx[r]-nvar1] = elts[r];
				} else {
					elts1[idx[r]] = elts[r];
				}
			}
			solver.addRow(elts1,elts2,scen,cut.lb(),cut.ub());
		}

		for (int i = 0; i < cuts.sizeColCuts(); i++) {
			const OsiColCut& cut = cuts.colCut(i);
			
			const CoinPackedVector &lbs = cut.lbs();
			{
				int nnz = lbs.getNumElements();
				const int* idx = lbs.getIndices();
				const double* elts = lbs.getElements();
				for (int r = 0; r < nnz; r++) {
					if (idx[r] >= nvar1) {
						printf("got second-stage column cut\n");
					} else {
						solver.setFirstStageColLB(idx[r],elts[r]);
					}
				}
			}
			
		}
	}

	solver.go();
	

	/*
	vector<double> firstSol = solver.getFirstStagePrimalColSolution();
	for (int i = 0; i < nvar1; i++) {
		if (firstSol[i] < 0.01) {
			firstSol[i] = 0.0;
		} else {
			firstSol[i] = 1.0;
		}
	}

	stochasticInput &input = *s;
	int nvar2 = input.nSecondStageVars(0);
	int ncons2 = input.nSecondStageCons(0);
	
	const vector<double> &obj1 = input.getFirstStageObj();

	vector<variableState> rowSave(ncons2), colSave(nvar2);
	bool havesave = false;
	double sum = 0.0;
	bool infeas = false;
	typedef ClpRecourseSolver RecourseSolver;
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
	*/

	MPI_Finalize();

	return 0;

}
