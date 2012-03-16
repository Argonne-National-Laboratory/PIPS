#ifndef BALPSOLVERPRIMAL_HPP
#define BALPSOLVERPRIMAL_HPP

#include "BALPSolverBase.hpp"

class BALPSolverPrimal : public BALPSolverBase {

public:
	BALPSolverPrimal(const BAData &);
	//~BALPSolverPrimal();
	
	void iterate();

	// use DSE prices if priceDSE, returns (basic) index of leaving variable
	BAIndex price() const;



	// main driver
	void go();

	void writeStatus(const std::string &filebase, bool appendIterateNumber = false,
			 bool writeBasisOnly = true);


protected:
	bool quickCheckDualFeasible() const;
	void makeFeasible(bool &basischange);

	void rebuildIndices() { setupIndices(); }

	
	virtual void initializeIndexedInfeasibilities();

	void removePerturbation();
	void updateDuals(const sparseBAVector &alpha, BAIndex leaveIdx, BAIndex enterIdx, const double thetad);
	template<bool doDevex> void updatePrimals(const sparseBAVector &aq, BAIndex enterIdx, double thetap);

	template<infeasType> BAIndex ratioHarris(const sparseBAVector &pivotCol);

	void forcePrimalFeasible(); // add perturbations to LB/UBs to force primal feasibility

	// maintain an indexed vector of (squared) dual infeasibilities.
	// uses global indices (we don't have non-basic indices).
	BAContainer<CoinIndexedVector> dualInfeas;

	// temporary vector used during ratio test
	packedBAVector packed; 

	// primal upper/lower bounds, with local perturbations
	denseBAVector lPerturb, uPerturb;
	bool didperturb;

	// edge weights for devex
	denseBAVector edgeWeights; 
	bool devexPricing;
	int devexIter; // iteration number of last reset
	double pivotDevex; // local part of pivot column devex weight
	
	BAFlagVector<bool> reference;

	void resetDevex();
	
	void updateDevexWeights(const sparseBAVector &pivotRow, double pivotValue,
		BAIndex enterIdx, BAIndex leaveIdx);

	virtual const denseBAVector & myLB() const { return lPerturb; }
	virtual const denseBAVector & myUB() const { return uPerturb; }

// indicates variable is actually feasible
#define FEASIBLE_TINY 1e-100

	// given global index, check if dual infeasible
	inline infeasType checkInfeas(const BAIndex idx) const {
		if (data.vartype[idx] == Fixed) return NotInfeasible;
		if (data.vartype[idx] != Free) {
			if (states[idx] == AtLower && d[idx] < -dualTol) {
				return Below;
			}
			if (states[idx] == AtUpper && d[idx] > dualTol) {
				return Above;
			}
		} else {
			if (d[idx] < -dualTol) return Below;
			if (d[idx] > dualTol) return Above;
		}
		return NotInfeasible;
	}

	// given global index, update dualInfeas
	inline void updateInfeas(const int scen, const int globalIdx) {
		BAIndex global = {scen, globalIdx};
		infeasType t = checkInfeas(global);
		double &infeasEntry = dualInfeas.getVec(scen).denseVector()[globalIdx];
		if (t != NotInfeasible) {
			double infeas = d[global];
			infeas *= infeas;
			if (infeasEntry) {
				infeasEntry = infeas;
			} else {
				dualInfeas.getVec(scen).insert(globalIdx,infeas);
			}
		} else {
			if (infeasEntry) {
				infeasEntry = FEASIBLE_TINY;
			}
		}
	}

	void loadStatus2(int scen, std::istream&);

};


#endif

