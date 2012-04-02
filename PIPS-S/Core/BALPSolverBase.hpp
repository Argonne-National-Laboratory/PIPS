#ifndef BALPSOLVERBASE_HPP
#define BALPSOLVERBASE_HPP

#include "BAData.hpp"
#include "BALinearAlgebra.hpp"
#include "BAPFIPar.hpp"
#include "BAPFIParWrapper.hpp"

#include <fstream>


// kind of infeasibility
enum infeasType { Below, Above, NotInfeasible };

class PIPSSInterface;

class BALPSolverBase {

public:
	BALPSolverBase(const BAData &);
	virtual ~BALPSolverBase();
	
	// set variable states 
	void setStates(const BAFlagVector<variableState> &);
	const BAFlagVector<variableState>& getStates() const { return states; }

	//void iterate();
	bool checkDualFeasible(BAContainer<std::vector<int> > &infeas);
	bool checkPrimalFeasible() const;

	solverState getStatus() const { return status; }
	
	void setLogIterFile(const char *f) { if (logIters) delete logIters; logIters = new std::ofstream(f); }
	void writeBasis(const std::string &filebase);

	const denseBAVector& getPrimalSolution() const { return x; }
	const denseBAVector& getDualColSolution() const { return d; }
	void setPrimalTolerance(double val) { primalTol = val; }
	void setDualTolerance(double val) { dualTol = val; }
	double getPrimalTolerance() const { return primalTol; }
	double getDualTolerance() const { return dualTol; }

	// write status (basis and edge weights)
	// format is specific to BALP structure, not interchangeable
	virtual void writeStatus(const std::string &filebase, bool appendIterateNumber = false,
			 bool writeBasisOnly = true) = 0;
	
	virtual void loadStatus(const std::string &filebase);



	virtual void go() = 0;

	void setReinversionFrequency(int r) { reinvertFrequency = reinvertFrequency_backup = r; }

	// will dump current status every d iterations. zero to disable.
	void setDumpFrequency(int d, const std::string &outputname) { dumpEvery = d; outputName = outputname; }

	friend class PIPSSInterface;

protected:

	int nIter;
	double startTime;

	// initial set-up for basic indices
	void setupIndices();
	// rebuild basic indices before reinverting,
	// copying over edge weights if appropriate
	virtual void rebuildIndices() = 0;
	// reinvert basis, initialize vectors, etc
	void initialize(bool reinvert = true);

	// dual has perturbed objective and primal has perturbed
	// LB and UB, so let these be overridden.
	virtual const denseBAVector & myObjective() const { return data.c; }
	virtual const denseBAVector & myLB() const { return data.l; }
	virtual const denseBAVector & myUB() const { return data.u; }

	virtual void initializeIndexedInfeasibilities() = 0;

	BAPFIParWrapper<BALinearAlgebra,BAPFIPar> *la;
	//BAPFI<BALinearAlgebra2> *la;
	
	const BAData &data;

	solverState status;

	// from data
	//void initializeBounds();

	// iterates
	denseBAVector x, d;

	// working vectors
	sparseBAVector rowVec; // for output of PRICE
	sparseBAVector ftranVec;
	sparseBAVector ftranVec2;
	sparseBAVector btranVec;


	// states completely describes current solution (except for edge weights)
	BAFlagVector<variableState> states;
	BAContainer<std::vector<int> > basicIdx;
	double objval;

	// tolerences. see section 6.2.1 in koberstein thesis 
	double dualTol, primalTol, zeroTol, primalRelTol, pivotTol;
	
	
	bool phase1;
	int doreport;
	int reportFrequency, reinvertFrequency, reinvertFrequency_backup;
	int lastGoodReinvert, lastBadReinvert, lastReinvert;
	int dumpEvery; std::string outputName;

	// record of densities
	double DSEDensity,ftranDensity,btranDensity;
	// record of pivots
	int replaceFirst, firstReplaceSecond, secondReplaceFirst, replaceSecondSelf, replaceSecondOther;

	void recordPivot(int enterScen, int leaveScen);

	// can be used to store a list of infeasible indices.
	// this isn't the index of infeasibilities for pricing.
	BAContainer<std::vector<int> > infeasList;


	double largestPrimalError;
	double calculateLargestPrimalError() const;
	int wasOptimal; // iteration number when appeared optimal

	std::ofstream *logIters;


	double ftranTime, btranTime, updateIteratesTime, ftranDSETime, invertTime,
		selectEnteringTime, selectLeavingTime, priceTime, updateColumnTime;

	virtual void loadStatus2(int scen, std::istream&) = 0;	

};


#endif

