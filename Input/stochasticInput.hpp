#ifndef STOCHASTICINPUTHPP
#define STOCHASTICINPUTHPP

#include <vector>
#include <string>
#include "CoinPackedMatrix.hpp"

// Stochastic (MI)LP readers will implement this virtual class
// Assume 2-stage stochastic (MI)LP with recourse, with a discrete distribution

/* 
constraints are of form:
A  x
T_1x + W_1y1
T_2x          + W_2y_2
...
T_Nx                            + W_Ny_N

where x are first-stage variables, and y second-stage variables
each row (constraint) and column (variable) can have an upper and lower bound
if these are equal for a row, then you have an equality constraint
if these are equal for a column, then you have a fixed variable
*/


// Meant to be a "nice" interface that can be used from other languages,
// e.g. with SWIG, so no pointers involved. 
// If worried about copying objects on return, use move semantics in 
// the implementation.
// Accessors not marked const because implementing class may need to
// change internal data structures when called, e.g. for caching the answer.
class stochasticInput {
public:
	virtual ~stochasticInput() {}
	virtual int nScenarios() = 0;
	virtual int nFirstStageVars() = 0;
	virtual int nFirstStageCons() = 0;
	virtual int nSecondStageVars(int scen) = 0;
	virtual int nSecondStageCons(int scen) = 0;

	virtual std::vector<double> getFirstStageColLB() = 0;
	virtual std::vector<double> getFirstStageColUB() = 0;
	virtual std::vector<double> getFirstStageObj() = 0;
	virtual std::vector<std::string> getFirstStageColNames() = 0;
	virtual std::vector<double> getFirstStageRowLB() = 0;
	virtual std::vector<double> getFirstStageRowUB() = 0;
	virtual std::vector<std::string> getFirstStageRowNames() = 0;
	virtual bool isFirstStageColInteger(int col) = 0;

	virtual std::vector<double> getSecondStageColLB(int scen) = 0;
	virtual std::vector<double> getSecondStageColUB(int scen) = 0;
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen) = 0;
	virtual std::vector<std::string> getSecondStageColNames(int scen) = 0;
	virtual std::vector<double> getSecondStageRowUB(int scen) = 0;
	virtual std::vector<double> getSecondStageRowLB(int scen) = 0;
	virtual std::vector<std::string> getSecondStageRowNames(int scen) = 0;
	virtual double scenarioProbability(int scen) = 0;
	virtual bool isSecondStageColInteger(int scen, int col) = 0;

	// returns the column-oriented first-stage constraint matrix (A matrix) 
	virtual CoinPackedMatrix getFirstStageConstraints() = 0;
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) = 0;
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen) = 0;

	

	// some problem characteristics that could be helpful to know
	
	// all scenarios have the same number of variables and constraints
	virtual bool scenarioDimensionsEqual() = 0;
	// constraint matrices are identical for each scenario,
	// column and row bounds and objective are allowed to vary
	virtual bool onlyBoundsVary() = 0;
	// all scenarios equally likely
	virtual bool allProbabilitiesEqual() = 0;
	// all second-stage variables continuous
	virtual bool continuousRecourse() = 0;

};

#endif

