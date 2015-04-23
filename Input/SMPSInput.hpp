#ifndef SMPSINPUT_HPP
#define SMPSINPUT_HPP

#include "stochasticInput.hpp"
#include "CoinMpsIO.hpp"
#include <cmath>

// Only works with SMPS files with SCENARIO format!
class SMPSInput : public stochasticInput {
public:
	virtual ~SMPSInput() {}
	SMPSInput(std::string const& cor, std::string const& tim, std::string const& sto);
	virtual int nScenarios() { return nscen; }
	virtual int nFirstStageVars() { return nvar1; }
	virtual int nFirstStageCons() { return ncons1; }
	virtual int nSecondStageVars(int scen) { return nvar2; }
	virtual int nSecondStageCons(int scen) { return ncons2; }

	virtual std::vector<double> getFirstStageColLB() { return firstStageData.collb; }
	virtual std::vector<double> getFirstStageColUB() { return firstStageData.colub; }
	virtual std::vector<double> getFirstStageObj() { return firstStageData.obj; }
	virtual std::vector<std::string> getFirstStageColNames() { return firstStageData.colname; }
	virtual std::vector<double> getFirstStageRowLB() { return firstStageData.rowlb; }
	virtual std::vector<double> getFirstStageRowUB() { return firstStageData.rowub; }
	virtual std::vector<std::string> getFirstStageRowNames() { return firstStageData.rowname; }
	virtual bool isFirstStageColInteger(int col) { return firstStageData.isColInteger.at(col); }
        virtual bool isFirstStageColBinary(int col) {
	  bool isInteger = this->isFirstStageColInteger(col);
	  // CoinMpsIO has no isBinary member function, but some preprocessing features require
	  // knowledge of binary variables, so kludge in an "isBinary" member function by
	  // relying on CoinMpsIO setting lower and upper bounds to zero and one, respectively.
	  // Also note: CoinMpsIO uses a default tolerance of 1.0e-8 on integrality comparisons.
	  const double intTol = 1.0e-8;
	  bool isLBzero = (fabs(this->getFirstStageColLB().at(col)) < intTol);
	  bool isUBone = (fabs(this->getFirstStageColUB().at(col) - 1.0) < intTol);
	  return (isInteger && isLBzero && isUBone);
	}

	virtual std::vector<double> getSecondStageColLB(int scen);
	virtual std::vector<double> getSecondStageColUB(int scen);
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen);
	virtual std::vector<std::string> getSecondStageColNames(int scen);
	virtual std::vector<double> getSecondStageRowUB(int scen);
	virtual std::vector<double> getSecondStageRowLB(int scen);
	virtual std::vector<std::string> getSecondStageRowNames(int scen);
	virtual double scenarioProbability(int scen) { return (probabilitiesequal) ? 1.0/nscen : probabilities.at(scen); }
	virtual bool isSecondStageColInteger(int scen, int col) { return secondStageTemplate.isColInteger.at(col); }
        virtual bool isSecondStageColBinary(int scen, int col) {
	  bool isInteger = this->isSecondStageColInteger(scen, col);
	  // CoinMpsIO has no isBinary member function, but some preprocessing features require
	  // knowledge of binary variables, so kludge in an "isBinary" member function by
	  // relying on CoinMpsIO setting lower and upper bounds to zero and one, respectively.
	  // Also note: CoinMpsIO uses a default tolerance of 1.0e-8 on integrality comparisons.
	  const double intTol = 1.0e-8;
	  bool isLBzero = (fabs(this->getSecondStageColLB(scen).at(col)) < intTol);
	  bool isUBone = (fabs(this->getSecondStageColUB(scen).at(col) - 1.0) < intTol);
	  return (isInteger && isLBzero && isUBone);
	}

	// returns the column-oriented first-stage constraint matrix (A matrix)
	virtual CoinPackedMatrix getFirstStageConstraints() { return firstStageData.mat; }
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen);
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen);



	virtual bool scenarioDimensionsEqual() { return true; }
	virtual bool onlyBoundsVary() { return onlyboundsvary; }
	virtual bool allProbabilitiesEqual() { return probabilitiesequal; }
	virtual bool continuousRecourse() { return continuousrecourse; }


private:
	struct problemData {

		problemData() : ncol(-1), nrow(-1) {}
		bool isInitialized() { return (ncol >= 0 && nrow >= 0); }

		int ncol, nrow;
		std::vector<double> collb, colub, rowlb, rowub, obj;
		std::vector<bool> isColInteger;
		CoinPackedMatrix mat;
		std::vector<std::string> colname, rowname;

		void initialize(int ncol, int nrow) {
			this->ncol = ncol;
			this->nrow = nrow;
			colub.resize(ncol);
			collb.resize(ncol);
			colname.resize(ncol);
			obj.resize(ncol);
			rowub.resize(nrow);
			rowlb.resize(nrow);
			rowname.resize(nrow);
			isColInteger.resize(ncol);
		}

	};

	void cacheScenario(int scen);

	int nscen, nvar1, ncons1, nvar2, ncons2;
	int nvar, ncons; // total variables
	std::vector<problemData> scenarioData;
	problemData firstStageData, secondStageTemplate;
	CoinPackedMatrix TmatTemplate;
	std::vector<CoinPackedMatrix> Tmats;
	std::vector<double> probabilities;
	std::string const corfile, timfile, stofile;
	CoinMpsIO reader;
	bool onlyboundsvary;
	bool probabilitiesequal;
	bool continuousrecourse;

	// save locations in sto file for later reading
	std::vector<std::streampos> scenarioStarts;
	std::vector<int> scenarioLens; // number of lines in sto file for each scenario, not including "SC" line



};

#endif
