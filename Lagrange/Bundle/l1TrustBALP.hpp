#ifndef L1BUNDLEBALP_HPP
#define L1BUNDLEBALP_HPP

#include "stochasticInput.hpp"
#include "bundleManager.hpp"

/* This is the l_1 trust region problem:

min \sum_i \theta_i 
s.t.\sum_i \gamma_i = 0
    \theta_i e_K - G^i\gamma_i >= c_i, i = 1, ..., N
    z_i >= \gamma_i - \gamma_i^+
    z_i >= \gamma_i^+ - \gamma_i
    \sum e^Tz_i <= \Delta

\gamma_i^+ and \Delta are given.

This LP has a primal block angular structure:


      (t)(g)(z)
min [e          ][e          ]   ... [e           ]
s.t.[     I      [     I         ... [      I        = 0
            -e^T]        -e^T]                -e^T] >= -\Delta
    [     I   I                                     >= \gamma_1^+
         -I   I                                     >= -\gamma_1^+
      e -G_1    ]                                   >= c_1
                 [     I   I                        >= \gamma_2^+
		      -I   I                        >= -\gamma_2^+
                   e -G_2    ]                      >= c_2
		                 ...
			            [      I   I    >= \gamma_N^+
				          -I   I    >= -\gamma_N^+
				      e  -G_N    ]  >= c_N

We formulate and solve the *dual* of the above as a dual block-angular LP.

*/

// note that the dual is a maximization problem, but input format assumes minimization.
// so, sign of objective is flipped!!
class l1TrustModel : public stochasticInput {
public:
	l1TrustModel(int nvar1, bundle_t const &cuts, double trustRadius, std::vector<std::vector<double> > const& center);
	virtual int nScenarios() { return cuts.size(); }
	virtual int nFirstStageVars() { return nvar1+1; }
	virtual int nFirstStageCons() { return 0; }
	virtual int nSecondStageVars(int scen) { return cuts[scen].size() + 2*nvar1; }
	virtual int nSecondStageCons(int scen) { return 2*nvar1 + 1; }

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

private:
	int nvar1; // dimension of each \gamma
	bundle_t const & cuts;
	std::vector<std::vector<double> > const& center;
	double trustRadius;


};


#endif
