#ifndef BALPSOLVERDUAL_HPP
#define BALPSOLVERDUAL_HPP

#include "BALPSolverBase.hpp"

class BALPSolverPrimal;

class BALPSolverDual : public BALPSolverBase {

public:
	BALPSolverDual(const BAData &);
	//~BALPSolverDual();
	
	void iterate();

	// use DSE prices if priceDSE, returns (basic) index of leaving variable
	BAIndex price() const;



	// main driver
	void go();

	void writeStatus(const std::string &filebase, bool appendIterateNumber = false,
			 bool writeBasisOnly = true);

protected:
	// use index of infeasibilities
	bool quickCheckPrimalFeasible() const;

	void makeFeasible();

	// rebuild basic indices before reinverting,
	// copying over DSE weights if appropriate
	void rebuildIndices();
	// calculate DSE weights from scratch
	void initializeDSE(bool slackbasis = false);
	// flip bounds on boxed dual infeasible variables to make them feasible
	void flipBounds();

	
	virtual void initializeIndexedInfeasibilities();

	void removePerturbation();
	void updateDSE(sparseBAVector &rho, const sparseBAVector &aq, BAIndex leave, double pivot);
	void updateDuals(const sparseBAVector &alpha, BAIndex leaveIdx, BAIndex enterIdx, const double thetad);
	void updatePrimals(const sparseBAVector &aq, BAIndex enterIdx, BAIndex enter, BAIndex leave, double thetap);

	// bound-flipping ratio test.
	BAIndex ratioBFRT(const sparseBAVector &alpha2, double delta0);
	BAIndex ratioHarris(const sparseBAVector &alpha2, double delta0);

	void forceDualFeasible(); // add perturbations to objective to force dual feasibility

	bool DSEPricing;
	
	// following 8.2.2.2, maintain an indexed vector of (squared) primal infeasibilities
	BAContainer<CoinIndexedVector> primalInfeas;

	denseBAVector dse; // dse weights

	// objective with our perturbations 
	denseBAVector cPerturb;
	bool didperturb;

	
	virtual const denseBAVector & myObjective() const { return cPerturb; }

// indicates variable is actually feasible
#define FEASIBLE_TINY 1e-100

	// given index, check if primal infeasible
	inline infeasType checkInfeas(const BAIndex idx) const {
			if (x[idx] < data.l[idx] - primalTol) {
				return Below;
			}
			if (x[idx] > data.u[idx] + primalTol) {
				return Above;
			}
		return NotInfeasible;
	}

	// given index in basis and global index, update primalInfeas
	inline void updateInfeas(const int scen, const int basicIdx, const int globalIdx) {
		BAIndex global = {scen, globalIdx};
		infeasType t = checkInfeas(global);
		double &infeasEntry = primalInfeas.getVec(scen).denseVector()[basicIdx];
		if (t != NotInfeasible) {
			double infeas;
			if (t == Below) {
				infeas = data.l[global]-x[global];
			} else {
				infeas = data.u[global]-x[global];
			}
			infeas *= infeas;
			if (infeasEntry) {
				infeasEntry = infeas;
			} else {
				primalInfeas.getVec(scen).insert(basicIdx,infeas);
			}
		} else {
			if (infeasEntry) {
				infeasEntry = FEASIBLE_TINY;
			}
		}
	}
	void loadStatus2(int scen, std::istream&);
	void setSlackBasis();

};


#endif

