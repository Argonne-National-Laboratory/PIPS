#ifndef LINFTRUSTBALP_HPP
#define LINFTRUSTBALP_HPP

#include "stochasticInput.hpp"
#include "bundleManager.hpp"

/*

This is the cutting plane model with an L_infinity trust region of radius \Delta around
the \gamma's. That is:

min (wrt. \theta_i,\gamma_i) \sum_i \theta_i
s.t. 
\sum_i \gamma_i = 0,
\theta_i e_K - G^i\gamma_i >= c_i, i = 1, ..., N
\gamma_i^+-\Delta e <= \gamma_i <= \gamma_i^+ + \Delta e

\theta_i's are free.
\gamma_i^+ and \Delta are given.


where G^i is a matrix with subgradients of f_i on the rows.
This comes from the subgradient inequality:
f_i(\gamma_i) >= f_i(\gamma_i^r) + g^T(\gamma_i - \gamma_i^r),
where g is a subgradient of f_i at the point \gamma_i^r.
From this we derive that the rth component of c_i is
f_i(\gamma_i^r) - g^T\gamma_i^r.

This LP has a primal block angular structure:


     (t) (g)  (t) (g)         (t) (g)
min [ 1     ][ 1     ]       [ 1     ]
s.t.[     I ][     I ]       [     I ]   = 0    (l)
    [ e -G_1]                           >= c_1  (u_1)
             [ e -G_2]                  >= c_2  (u_2)
                        ...
                             [ e -G_N]  >= c_N  (u_N)

We formulate and solve the *dual* of the LP as a dual block-angular problem, 
because we have a solver for those!

*/



// note that the dual is a maximization problem, but input format assumes minimization.
// so, sign of objective is flipped!!
class lInfTrustModel : public stochasticInput {
public:
	lInfTrustModel(int nvar1, bundle_t const &cuts, double trustRadius, std::vector<std::vector<double> > const& center);
	virtual int nScenarios() { return cuts.size(); }
	virtual int nFirstStageVars() { return nvar1 ; }
	virtual int nFirstStageCons() { return 0; }
	virtual int nSecondStageVars(int scen) { return cuts[scen].size()+2*nvar1; }
	virtual int nSecondStageCons(int scen) { return nvar1 + 1; }

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
	double trustRadius;
	std::vector<std::vector<double> > const& center;


};



#endif
