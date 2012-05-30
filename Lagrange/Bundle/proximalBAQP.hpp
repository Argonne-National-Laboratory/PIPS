#ifndef PROXIMALBAQP_HPP
#define PROXIMALBAQP_HPP

#include "stochasticInput.hpp"
#include "bundleManager.hpp"

// See doc/PROJECTS/LAGRANGE/NOTES/proximalqp.tex for the formulation
// This is the version with y rescaled.


class proximalQPModel : public stochasticInput {
public:
	proximalQPModel(int nvar1, bundle_t const &cuts, std::vector<std::vector<double> > const& proxCenter, double tau);
	virtual int nScenarios() { return cuts.size(); }
	virtual int nFirstStageVars() { return nvar1; }
	virtual int nFirstStageCons() { return 0; }
	virtual int nSecondStageVars(int scen) { return cuts[scen].size(); }
	virtual int nSecondStageCons(int scen) { return 1; }

	virtual std::vector<double> getFirstStageColLB();
	virtual std::vector<double> getFirstStageColUB();
	virtual std::vector<double> getFirstStageObj();
	virtual std::vector<std::string> getFirstStageColNames();
	virtual std::vector<double> getFirstStageRowLB();
	virtual std::vector<double> getFirstStageRowUB();
	virtual std::vector<std::string> getFirstStageRowNames();
	virtual bool isFirstStageColInteger(int col) { return false; }

	virtual std::vector<double> getSecondStageColLB(int scen);
	virtual std::vector<double> getSecondStageColUB(int scen);
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen);
	virtual std::vector<std::string> getSecondStageColNames(int scen);
	virtual std::vector<double> getSecondStageRowUB(int scen);
	virtual std::vector<double> getSecondStageRowLB(int scen);
	virtual std::vector<std::string> getSecondStageRowNames(int scen);
	virtual double scenarioProbability(int scen) { return 1.0/cuts.size(); }
	virtual bool isSecondStageColInteger(int scen, int col) { return false; }

	virtual CoinPackedMatrix getFirstStageConstraints();
	virtual CoinPackedMatrix getSecondStageConstraints(int scen);
	virtual CoinPackedMatrix getLinkingConstraints(int scen);

	

	virtual bool scenarioDimensionsEqual() { return false; } // not worth checking
	virtual bool onlyBoundsVary() { return false; }
	virtual bool allProbabilitiesEqual() { return true; }
	virtual bool continuousRecourse() { return true; }

	virtual CoinPackedMatrix getFirstStageHessian();
	virtual CoinPackedMatrix getSecondStageHessian(int scen);
	virtual CoinPackedMatrix getSecondStageCrossHessian(int scen);

private:
	int nvar1; // dimension of each \gamma
	bundle_t const & cuts;
	std::vector<std::vector<double> > const& proxCenter;
	double tau;


};




#endif
