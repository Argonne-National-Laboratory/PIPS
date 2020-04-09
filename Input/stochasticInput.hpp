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
        virtual std::vector<double> getLinkRowLB(){ return std::vector<double>(); }
        virtual std::vector<double> getLinkRowUB(){ return std::vector<double>(); }
	virtual std::vector<std::string> getFirstStageRowNames() = 0;
	virtual bool isFirstStageColInteger(int col) = 0;
	virtual bool isFirstStageColBinary(int col) = 0;

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
	virtual bool isSecondStageColBinary(int scen, int col) = 0;

	// returns the column-oriented first-stage constraint matrix (A matrix) 
	virtual CoinPackedMatrix getFirstStageConstraints() = 0;
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) = 0;
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen) = 0;

	

	// some problem characteristics that could be helpful to know
	
	// all scenarios have the same number of variables and constraints
	virtual bool scenarioDimensionsEqual() = 0;
	// constraint (and hessian) matrices are identical for each scenario,
	// column and row bounds and objective are allowed to vary
	virtual bool onlyBoundsVary() = 0;
	// all scenarios equally likely
	virtual bool allProbabilitiesEqual() = 0;
	// all second-stage variables continuous
	virtual bool continuousRecourse() = 0;

	/* Quadratic terms:
	We allow the input to specify quadratic objective in the form of:
	(1/2)x^TQx + c^T + \sum_{i=1}^N p_i ((1/2)y_i^TQ_iy_i + x_i^T{\hat Q_i^T}y_i + c_i^Ty)

	Q is the first-stage hessian
	Q_i is the second-stage hessian
	\hat Q_i is the second-stage cross hessian

	Default implementations are provided so that these do not need to be implemented if not used
	*/

	// column-oriented *lower triangle only*
	// Q
	virtual CoinPackedMatrix getFirstStageHessian();
	// Q_i
	virtual CoinPackedMatrix getSecondStageHessian(int scen);
	// column-oriented, \hat Q_i
	// Note: this has the second-stage variables on the rows and first-stage on the columns
	virtual CoinPackedMatrix getSecondStageCrossHessian(int scen);


        virtual int nLinkCons(){ return 0; }
        virtual int nLinkECons(){ return 0; }
        virtual int nLinkICons(){ return 0; }
        virtual CoinPackedMatrix getLinkMatrix(int nodeid){ return CoinPackedMatrix(); }
	std::string datarootname;
    int useInputDate;

};



#endif

