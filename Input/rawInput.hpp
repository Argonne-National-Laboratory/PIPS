#ifndef RAWINPUT_HPP
#define RAWINPUT_HPP

#include "stochasticInput.hpp"
#include "mpi.h"

// reads format written by dumpSmlModel.cpp
class rawInput : public stochasticInput {
public:

	rawInput(const std::string &datarootname, int overrideScenarioNumber = 0, MPI_Comm comm = MPI_COMM_WORLD);

	virtual int nScenarios() { return nScenarios_; }
	virtual int nFirstStageVars() { return nFirstStageVars_; }
	virtual int nFirstStageCons() { return nFirstStageCons_; }
	virtual int nSecondStageVars(int scen) { return nSecondStageVars_; }
	virtual int nSecondStageCons(int scen) { return nSecondStageCons_; }

	virtual std::vector<double> getFirstStageColLB() { return firstStageData.collb; }
	virtual std::vector<double> getFirstStageColUB() { return firstStageData.colub; }
	virtual std::vector<double> getFirstStageObj() { return firstStageData.obj; }
	virtual std::vector<std::string> getFirstStageColNames() { return firstStageData.colnames; }
	virtual std::vector<double> getFirstStageRowLB() { return firstStageData.rowlb; }
	virtual std::vector<double> getFirstStageRowUB() { return firstStageData.rowub; }
	virtual std::vector<std::string> getFirstStageRowNames() { return firstStageData.rownames; }
	virtual bool isFirstStageColInteger(int col) { return false; }

	virtual std::vector<double> getSecondStageColLB(int scen);
	virtual std::vector<double> getSecondStageColUB(int scen);
	virtual std::vector<double> getSecondStageObj(int scen);
	virtual std::vector<std::string> getSecondStageColNames(int scen);
	virtual std::vector<double> getSecondStageRowUB(int scen);
	virtual std::vector<double> getSecondStageRowLB(int scen);
	virtual std::vector<std::string> getSecondStageRowNames(int scen);
	virtual double scenarioProbability(int scen) { return 1.0/nScenarios_; }
	virtual bool isSecondStageColInteger(int scen, int col) { return false; }

	virtual CoinPackedMatrix getFirstStageConstraints() { return Amat; }
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) { return Wmat; }
	virtual CoinPackedMatrix getLinkingConstraints(int scen) { return Tmat; }

	

	virtual bool scenarioDimensionsEqual() { return true; }
	virtual bool onlyBoundsVary() { return true; }
	virtual bool allProbabilitiesEqual() { return true; }
	virtual bool continuousRecourse() { return true; }
	

protected:
	struct scenData {
		std::vector<double> collb, colub, rowlb, rowub, obj;
		std::vector<std::string> rownames, colnames;
		bool didLoad;
		scenData() : didLoad(false) {}
		void initialize(int nvar, int ncons);
	};
	void loadLocalScenData(int scen);
	CoinPackedMatrix Amat, Tmat, Wmat;
	std::vector<scenData> localData;
	const std::string datarootname;
	scenData firstStageData;
	int nScenariosTrue;
	int mype_;

	int nScenarios_, nFirstStageVars_, nFirstStageCons_, nSecondStageVars_, nSecondStageCons_;

};


#endif
