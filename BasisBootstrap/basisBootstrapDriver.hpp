#ifndef BASISBOOTSTRAPDRIVER_HPP
#define BASISBOOTSTRAPDRIVER_HPP

#include "BALPSolverInterface.hpp"
#include "BAVector.hpp"
#include <boost/scoped_ptr.hpp>

class subsetWrapper;

// serial-only code for now
template<typename BALPSolver, typename RecourseSolver>
class basisBootstrapDriver {
public:
	
	basisBootstrapDriver(stochasticInput &input, BAContext& ctx) : input(input), ctx(ctx) {
		masterIsValid = false;
		feasible = false;
		int nscen = input.nScenarios();
		rowStates2.resize(nscen);
		colStates2.resize(nscen);
	}

	// start from scratch
	void initializeMaster(std::vector<int> const& scenariosInMaster) {
		int nscen = input.nScenarios();
		this->scenariosInMaster = scenariosInMaster;
		scenariosNotInMaster.clear();
		std::vector<bool> in(nscen,false);
		for (unsigned i = 0; i < scenariosInMaster.size(); i++) {
			in[scenariosInMaster[i]] = true;
		}
		for (int i = 0; i < nscen; i++) {
			if (in[i] == false) scenariosNotInMaster.push_back(i);
		}
		firstSolve = true;
	}
	void solveMaster() {
		subsetWrapper wrapper(input, scenariosInMaster);
		ctx.initializeAssignment(scenariosInMaster.size());
		solver.reset(new BALPSolver(wrapper, ctx, feasible ? BALPSolver::usePrimal : BALPSolver::useDual));

		int nvar1 = input.nFirstStageVars();
		int ncons1 = input.nFirstStageCons();
		int nscen = input.nScenarios();
		
		if (!firstSolve) { // copy states
			for (int i = 0; i < nvar1; i++) {
				solver->setFirstStageColState(i,colStates1[i]);
			}
			for (int i = 0; i < ncons1; i++) {
				solver->setFirstStageRowState(i,rowStates1[i]);
			}
			for (unsigned k = 0; k < scenariosInMaster.size(); k++) {
				int scen = scenariosInMaster[k];
				int nvar2 = input.nSecondStageVars(scen);
				int ncons2 = input.nSecondStageCons(scen);
				for (int i = 0; i < nvar2; i++) {
					solver->setSecondStageColState(k,i,colStates2[scen][i]);
				}
				for (int i = 0; i < ncons2; i++) {
					solver->setSecondStageRowState(k,i,rowStates2[scen][i]);
				}
			}
			solver->commitStates();
		}
		solver->go();

		colStates1.resize(nvar1); // won't cause reallocation if already done
		rowStates1.resize(ncons1);
		for (int i = 0; i < nvar1; i++) {
			colStates1[i] = solver->getFirstStageColState(i);
		}
		for (int i = 0; i < ncons1; i++) {
			rowStates1[i] = solver->getFirstStageRowState(i);
		}
		for (int k = 0; k < scenariosInMaster.size(); k++) {
			int scen = scenariosInMaster[k];
			int nvar2 = input.nSecondStageVars(scen);
			int ncons2 = input.nSecondStageCons(scen);
			colStates2[scen].resize(nvar2);
			rowStates2[scen].resize(ncons2);
			for (int i = 0; i < nvar2; i++) {
				colStates2[scen][i] = solver->getSecondStageColState(k,i);
			}
			for (int i = 0; i < ncons2; i++) {
				rowStates2[scen][i] = solver->getSecondStageRowState(k,i);
			}
		}


		firstSolve = false;
		masterIsValid = true;
		feasible = true;
	}
	std::vector<double> getFirstStageSolution() const { assert(solver.get()); return solver->getFirstStagePrimalColSolution(); }

	void addScenarioToMaster(int scen) {
		RecourseSolver rsolver(input, scen, getFirstStageSolution());

		// may need to change this
		rsolver.setDualObjectiveLimit(1e10);

		int nvar2 = input.nSecondStageVars(scen);
		int ncons2 = input.nSecondStageCons(scen);
		if (goodRowStates.size() && input.onlyBoundsVary()) {
			for (int i = 0; i < nvar2; i++) {
				rsolver.setSecondStageColState(i,goodColStates[i]);
			}
			for (int i = 0; i < ncons2; i++) {
				rsolver.setSecondStageRowState(i,goodRowStates[i]);
			}
			rsolver.commitStates();
		}
		rsolver.go();

		std::vector<variableState> cols(nvar2), rows(ncons2);
		for (int i = 0; i < nvar2; i++) {
			cols[i] = rsolver.getSecondStageColState(i);
		}
		for (int i = 0; i < ncons2; i++) {
			rows[i] = rsolver.getSecondStageRowState(i);
		}

		addScenarioToMaster(scen, cols, rows, rsolver.getStatus() == Optimal);

	}

	void addScenarioToMaster(int scen, std::vector<variableState> const& colstates, 
				 std::vector<variableState> const& rowstates, bool isFeasible) {
		feasible = feasible && isFeasible;
		assert(feasible);
		if (feasible) {
			if (input.onlyBoundsVary()) {
				goodRowStates = rowstates;
				goodColStates = colstates;
			}
			rowStates2[scen] = rowstates;
			colStates2[scen] = colstates;
		} else {
			if (input.onlyBoundsVary()) {
				// reuse a good basis
				rowStates2[scen] = goodRowStates;
				colStates2[scen] = goodColStates;
			} else {
				// slack basis
				colStates2[scen].clear();
				colStates2[scen].resize(input.nSecondStageVars(scen),AtLower);
				rowStates2[scen].clear();
				rowStates2[scen].resize(input.nSecondStageCons(scen),Basic);
			}
		}
		scenariosInMaster.push_back(scen);
		scenariosNotInMaster.erase(std::find(scenariosNotInMaster.begin(),scenariosNotInMaster.end(),scen));
		assert(scenariosInMaster.size() + scenariosNotInMaster.size() == input.nScenarios());
		masterIsValid = false;
		
	}

	void writeStatus(std::string const &filebase);
	void loadStatus(std::string const &filebase);

	double getObjective() {
		assert(masterIsValid && scenariosInMaster.size() == input.nScenarios());
		return solver->getObjective();
	}

protected:
	std::vector<variableState> rowStates1, colStates1;
	std::vector<variableState> goodRowStates, goodColStates; // basis from a feasible recourse problem, for hot starting
	std::vector<std::vector<variableState> > rowStates2, colStates2;
	std::vector<int> scenariosInMaster, scenariosNotInMaster;
	BAContext &ctx;
	stochasticInput &input;
	bool masterIsValid; // false if we've added scenarios and not reformulated/resolved master
	bool feasible; // if all scenarios are feasible (so we can use primal simplex)
	bool firstSolve; // we haven't solved a master problem before
	

	boost::scoped_ptr<BALPSolver> solver;


};

// wrap stochasticInput so that only a subset of scenarios are presented
class subsetWrapper : public stochasticInput {
public:
        subsetWrapper(stochasticInput &inner, std::vector<int> const& includedScenarios) :
		realScenarios(includedScenarios), probabilitySum(0.), inner(inner) {
			for (unsigned i = 0; i < realScenarios.size(); i++) {
				probabilitySum += inner.scenarioProbability(realScenarios[i]);
			}
	}
	virtual int nScenarios() { return realScenarios.size(); }
        virtual int nFirstStageVars() { return inner.nFirstStageVars(); }
        virtual int nFirstStageCons() { return inner.nFirstStageCons(); }
        virtual int nSecondStageVars(int scen) { return inner.nSecondStageVars(realScenarios[scen]); }
        virtual int nSecondStageCons(int scen) { return inner.nSecondStageCons(realScenarios[scen]); }

        virtual std::vector<double> getFirstStageColLB() { return inner.getFirstStageColLB(); }
        virtual std::vector<double> getFirstStageColUB() { return inner.getFirstStageColUB(); }
        virtual std::vector<double> getFirstStageObj() { return inner.getFirstStageObj(); }
        virtual std::vector<std::string> getFirstStageColNames() { return inner.getFirstStageColNames(); }
        virtual std::vector<double> getFirstStageRowLB() { return inner.getFirstStageRowLB(); }
        virtual std::vector<double> getFirstStageRowUB() { return inner.getFirstStageRowUB(); }
        virtual std::vector<std::string> getFirstStageRowNames() { return inner.getFirstStageRowNames(); }
        virtual bool isFirstStageColInteger(int col) { return inner.isFirstStageColInteger(col); }

        virtual std::vector<double> getSecondStageColLB(int scen) { return inner.getSecondStageColLB(realScenarios[scen]); }
        virtual std::vector<double> getSecondStageColUB(int scen) { return inner.getSecondStageColUB(realScenarios[scen]); }
        virtual std::vector<double> getSecondStageObj(int scen) {
		std::vector<double> v = inner.getSecondStageObj(realScenarios[scen]);
		for (unsigned i = 0; i < v.size(); i++) v[i] /= probabilitySum;
		return v;
	}
        virtual std::vector<std::string> getSecondStageColNames(int scen) { return inner.getSecondStageColNames(realScenarios[scen]); }
        virtual std::vector<double> getSecondStageRowUB(int scen) { return inner.getSecondStageRowUB(realScenarios[scen]); }
        virtual std::vector<double> getSecondStageRowLB(int scen) { return inner.getSecondStageRowLB(realScenarios[scen]); }
        virtual std::vector<std::string> getSecondStageRowNames(int scen) { return inner.getSecondStageRowNames(realScenarios[scen]); }
        virtual double scenarioProbability(int scen) { return inner.scenarioProbability(realScenarios[scen])/probabilitySum; }
        virtual bool isSecondStageColInteger(int scen, int col) { return inner.isSecondStageColInteger(realScenarios[scen],col); }

        virtual CoinPackedMatrix getFirstStageConstraints() { return inner.getFirstStageConstraints(); }
        virtual CoinPackedMatrix getSecondStageConstraints(int scen) { return inner.getSecondStageConstraints(realScenarios[scen]); }
        virtual CoinPackedMatrix getLinkingConstraints(int scen) { return inner.getLinkingConstraints(realScenarios[scen]); }

        

        virtual bool scenarioDimensionsEqual() { return inner.scenarioDimensionsEqual(); }
        virtual bool onlyBoundsVary() { return inner.onlyBoundsVary(); }
        virtual bool allProbabilitiesEqual() { return inner.allProbabilitiesEqual(); }
        virtual bool continuousRecourse() { return inner.continuousRecourse(); }

private:
	std::vector<int> const realScenarios;
	double probabilitySum;
	stochasticInput &inner;	
};

// select scenarios by number only
template<typename BALPSolver, typename RecourseSolver>
class simpleBasisBootstrapDriver : public basisBootstrapDriver<BALPSolver,RecourseSolver> {
public:
	simpleBasisBootstrapDriver(stochasticInput &input, BAContext& ctx) : basisBootstrapDriver<BALPSolver,RecourseSolver>(input,ctx) {}

	void goFromScratch(int startingScenarios, int addedPer) {
		std::vector<int> s;
		for (int i = 0; i < startingScenarios; i++) s.push_back(i);
		this->initializeMaster(s);
		this->solveMaster();
		
		while(this->scenariosNotInMaster.size()) {
			for (int i = 0; i < addedPer && this->scenariosNotInMaster.size(); i++) {
				this->addScenarioToMaster(this->scenariosNotInMaster[0]);
			}
			this->solveMaster();
		}

	}


};
#endif
