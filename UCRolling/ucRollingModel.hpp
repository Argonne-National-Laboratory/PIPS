#ifndef UCROLLINGMODEL_HPP
#define UCROLLINGMODEL_HPP

#include "stochasticInput.hpp"
#include "mpi.h"

#include <map>

// simple object that represents a set of variables
// all it does it keep a map of the indices
// assumes variables are indexed by a zero-based index and by a template parameter T,
// which is a key to a given map
template <typename T, typename Tmap> class variableSet {
public:
	void initialize(int nIndex1, int offset, Tmap const& map) {
		idxmap = &map; this->offset = offset; this->nIndex1 = nIndex1;
	}

	// given zero-based index
	int operator()(int i1, int i2) {
		return i1+ nIndex1*i2 + offset;
	}

	// given key
	int lookup(int i1, T i2) {
		return i1+ nIndex1*idxmap->find(i2)->second + offset;
	}

	int totalVars() {
		return nIndex1*idxmap->size();
	}

private:
	int nIndex1; // first index is 0,1,..,nIndex1-1 (doesn't include nIndex1)
	int offset;
	Tmap const *idxmap;
};

class ucRollingModel : public stochasticInput {
public:
	// when initial conditions are unknown
	ucRollingModel(std::string const& dataRoot, int nscen, int tOffset, int tHorizon,
		MPI_Comm = MPI_COMM_WORLD);
	// given initial conditions
	ucRollingModel(std::string const& dataRoot, int nscen, int tOffset, int tHorizon,
		std::vector<double> const& Pgen_init, std::vector<double> const& SOC_init,
		MPI_Comm = MPI_COMM_WORLD);
	
	virtual int nScenarios() { return nscen; }
	virtual int nFirstStageVars() { return nvar1; }
	virtual int nFirstStageCons() { return 0; }
	virtual int nSecondStageVars(int scen) { return nvar2; }
	virtual int nSecondStageCons(int scen) { return ncons2; }

	virtual std::vector<double> getFirstStageColLB();
	virtual std::vector<double> getFirstStageColUB();
	virtual std::vector<double> getFirstStageObj();
	virtual std::vector<std::string> getFirstStageColNames();
	virtual std::vector<double> getFirstStageRowLB();
	virtual std::vector<double> getFirstStageRowUB();
	virtual std::vector<std::string> getFirstStageRowNames();
	virtual bool isFirstStageColInteger(int col) { return true; }

	virtual std::vector<double> getSecondStageColLB(int scen);
	virtual std::vector<double> getSecondStageColUB(int scen);
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen);
	virtual std::vector<std::string> getSecondStageColNames(int scen);
	virtual std::vector<double> getSecondStageRowUB(int scen);
	virtual std::vector<double> getSecondStageRowLB(int scen);
	virtual std::vector<std::string> getSecondStageRowNames(int scen);
	virtual double scenarioProbability(int scen) { return 1.0/nscen; }
	virtual bool isSecondStageColInteger(int scen, int col) { return false; }

	// returns the column-oriented first-stage constraint matrix (A matrix) 
	virtual CoinPackedMatrix getFirstStageConstraints();
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen);
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen);

	

	virtual bool scenarioDimensionsEqual() { return true; }
	virtual bool onlyBoundsVary() { return true; }
	virtual bool allProbabilitiesEqual() { return true; }
	virtual bool continuousRecourse() { return true; }

private:

	void readData(std::string const& dataRoot, MPI_Comm);
	void generateWind(int scen, double sigma);
	void initializeVariables();

	/* PROBLEM DATA */
	int nscen;
	int timeOffset; // offset from time zero
	int horizon; // horizon length (number of time periods in the future, not including now)
	double sigma;

	struct genStruct {
		int bus_gen;
		double np_cap;
		double sum_cap;
		double win_cap;
		double min_cap;
		std::string fuel;
		double min_hrate;
		double min_power;
		double max_ur;
		double max_dr;
		double gen_cost;
		double Pgen_init;
	};

	std::vector<genStruct> genData;
	std::map<int,int> genMap;

	struct linStruct {
		int idx;
		int snd_bus;
		int rec_bus;
		double X;
		double V;
		double Pmax;
	};

	std::vector<linStruct> linData;
	std::map<std::string,int> linMap;

	struct fuelStruct {
		double HV;
		double Unitprice;
	};
	
	std::vector<fuelStruct> fuelData;
	std::map<std::string,int> fuelMap;
	
	/*struct busStruct {
		int idx;
	};*/
	std::map<int,int> busMap;
	//std::map<int,busStruct> busData;
	int ref_bus;

	struct loadStruct {
		int bus_load;
		double SOCcap;
		double SOC_init;
	};

	std::map<int,int> loadMap;
	std::vector<loadStruct> loadData;

	// outer index is time, inner is LOAD index
	std::vector<std::vector<double> > loads;

	struct windStruct {
		int bus_wind;
		double wind_share;
	};
	
	std::map<int,int> windMap;
	std::vector<windStruct> windData;
	
	// indexed by time
	std::vector<double> wind_total_determ;
	std::vector<double> SOCprof_determ;

	/* VARIABLES */
	variableSet<int,std::map<int,int> > yuc; // these are the only first stage variables
	variableSet<int,std::map<int,int> > Pgen, dPgen;
	variableSet<int,std::map<int,int> > Pwind;
	variableSet<int,std::map<int,int> > theta, slackp, slackm;
	variableSet<std::string,std::map<std::string,int> > P;
	variableSet<int,std::map<int,int> > SOC, Psto;
	

	/* Instance data */
	int nvar1, nvar2, ncons2;
	bool givenInitial;
	// randomly generated winds for each scenario
	// outer index is scenario, inner index is time (global index)
	std::vector<std::vector<double> > wind_total;

};



#endif
