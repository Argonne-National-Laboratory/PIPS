#ifndef BA_HPP
#define BA_HPP

#include <vector>
#include <cstdio>
#include "stochasticInput.hpp"
#include "mpi.h"

// BA == block angular
// the language of stochastic programming is used, e.g. stage, scenario,
// however, the code works with any dual block-angular problem
// first stage refers to linking variables
// second stage/scenario refers to blocks

// define basic structures

// struct for indexing vectors
// scen == -1 means first stage
// this should not be used for looping through vectors,
// but for specifying columns/rows of the constraint matrix
struct BAIndex {
	int scen;
	int idx;
	bool operator==(const BAIndex&r) const { return (idx == r.idx && scen == r.scen); }
};

// communication context class
// contains MPI communicator 
// handles logic for assigning scenarios
class BAContext {
public:
	BAContext(MPI_Comm comm);
	
	// returns true if current MPI process is assigned scenario idx
	bool assignedScenario(int idx) const { 
		if (idx == -1) return true; 
		else return (assignedScen.at(idx) == _mype);
	}
	// returns MPI process that "owns" scenario
	int owner(int scen) const;
	// reduce (sum) a single value (convenience functions)
	double reduce(double in) const;
	int reduce(int in) const;

	int mype() const { return _mype; }
	int nprocs() const { return _nprocs; }
	MPI_Comm comm() const { return mpicomm; }
	void initializeAssignment(int nscen);
	// localScenarios is list of scenarios for which a vector has been allocated (except for basic vectors)
	// includes -1 as first element, indicating first stage
	// IF YOU WANT TO ONLY ITERATE OVER SECOND STAGE, START AT localScen[1]! 
	const std::vector<int>& localScenarios() const { return localScen; } 

protected:
	int _mype;
	int _nprocs;
	int nper;
	MPI_Comm mpicomm;
	std::vector<int> assignedScen; // map from scenario to assigned proc
	std::vector<int> localScen;


};

// for accessing dimensions of problem
class BADimensions {

public:
	BADimensions() {}
	BADimensions(stochasticInput &i, const BAContext&);
	BADimensions(const BADimensions& d) { *this = d; } 
	virtual ~BADimensions() {}
	virtual int numFirstStageVars() const { return nFirstStageVars_; }
	virtual int numFirstStageCons() const { return nFirstStageCons_; }
	virtual int numSecondStageVars(int idx) const { if (same) return nSecondStageVars_; else return nSecondStageVars_v[idx]; }
	virtual int numSecondStageCons(int idx) const { if (same) return nSecondStageCons_; else return nSecondStageCons_v[idx]; }
	virtual int numScenarios() const { return nScen; }
	virtual int totalVars() const {
		int v = this->numFirstStageVars();
		for (int i = 0; i < numScenarios(); i++) {
			v += this->numSecondStageVars(i);
		}
		return v;
	}
	virtual int totalCons() const {
		int v = this->numFirstStageCons();
		
		for (int i = 0; i < numScenarios(); i++) {
			v += this->numSecondStageCons(i);
		}
		return v;
	}
	virtual int numVars(int idx) const {
		if (idx == -1) return this->numFirstStageVars();
		else return this->numSecondStageVars(idx);
	}
	virtual int numCons(int idx) const {
		if (idx == -1) return this->numFirstStageCons();
		else return this->numSecondStageCons(idx);
	}

protected:
	int nFirstStageVars_, nFirstStageCons_, nScen;
	int nSecondStageVars_, nSecondStageCons_;
	std::vector<int> nSecondStageVars_v, nSecondStageCons_v;
	bool same;

};

// add slacks
class BADimensionsSlacks : public BADimensions {
public:
	BADimensionsSlacks(const BADimensions &inner) : inner(inner) { }
	BADimensionsSlacks() {}

	virtual int numScenarios() const { return inner.numScenarios(); }
	virtual int numFirstStageCons() const { return inner.numFirstStageCons(); }
	virtual int numFirstStageVars() const { return inner.numFirstStageVars() + inner.numFirstStageCons(); }
	virtual int numSecondStageVars(int idx) const { return inner.numSecondStageVars(idx) + inner.numSecondStageCons(idx); }
	virtual int numSecondStageCons(int idx) const { return inner.numSecondStageCons(idx); }

	BADimensions inner;
};



#endif
