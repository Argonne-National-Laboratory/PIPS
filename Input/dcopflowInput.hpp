/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/
 

#ifndef DCOPFINPUT_HPP
#define DCOPFINPUT_HPP

#include "stochasticInput.hpp"
#include "mpi.h"

#include "ps.h"
#include "dcopflow.hpp"


//  solve problem as general stochastic programming problem
//  min f(x_0) + \sum f(x_i)
//  st    g(x_i)  = 0
//

// Jacobian looks like 
//  	W_1           		T_1	
//		W_2           	T_1	
//			...
//				W_n T_n
//					  A

// Hessian looks like 
//  	QW_1           			QT_1	
//			QW_2           	QT_1	
//				      ...
//				     QW_n	QT_n
//					  	  QA



// reads format written by dumpSmlModel.cpp
class dcopflowInput : public stochasticInput {
public:
    dcopflowInput(){};
	~dcopflowInput();
    dcopflowInput(const DCPS* powersys, MPI_Comm comm = MPI_COMM_WORLD);
  
	virtual int nScenarios() { return nScenarios_; }
	virtual int nFirstStageVars() { return nFirstStageVars_; }
	virtual int nFirstStageCons() { return nFirstStageCons_; }
	virtual int nSecondStageVars(int scen) { return nSecondStageVars_[scen]; }
	virtual int nSecondStageCons(int scen) { return nSecondStageCons_[scen]; }

	virtual bool isFirstStageColInteger(int col) { return false; }	
	virtual bool isSecondStageColInteger(int scen, int col) { return false; }
	
	virtual bool scenarioDimensionsEqual() { return false; }
	virtual bool onlyBoundsVary() { return false; }
	virtual bool allProbabilitiesEqual() { return false; }
	virtual bool continuousRecourse() { return false; }

    virtual double scenarioProbability(int scen) { return 1.0; }

  	// First stage obj and bounds   
  	virtual std::vector<double> getFirstStageObj(){ return firstStageData.objLin; }
  	virtual std::vector<double> getFirstStageColLB(){ return firstStageData.collb; }
  	virtual std::vector<double> getFirstStageColUB(){ return firstStageData.colub; } 
  	virtual std::vector<std::string> getFirstStageColNames() { return firstStageData.colnames; }
  	virtual std::vector<double> getFirstStageRowLB(){ return firstStageData.rowlb; }
  	virtual std::vector<double> getFirstStageRowUB(){ return firstStageData.rowub; }  
  	virtual std::vector<std::string> getFirstStageRowNames() { return firstStageData.rownames; }

  	// Second stage obj and bounds  
  	virtual std::vector<double> getSecondStageObj(int scen){loadLocalScenData(scen);return localData[scen].objLin;}  
  	virtual std::vector<double> getSecondStageColLB(int scen){loadLocalScenData(scen);return localData[scen].collb;}
  	virtual std::vector<double> getSecondStageColUB(int scen){loadLocalScenData(scen);return localData[scen].colub;}
  	virtual std::vector<std::string> getSecondStageColNames(int scen){loadLocalScenData(scen);return localData[scen].colnames;}
  	virtual std::vector<double> getSecondStageRowLB(int scen){loadLocalScenData(scen);return localData[scen].rowlb;}
  	virtual std::vector<double> getSecondStageRowUB(int scen){loadLocalScenData(scen);return localData[scen].rowub;}  
  	virtual std::vector<std::string> getSecondStageRowNames(int scen){loadLocalScenData(scen);return localData[scen].rownames;}

  	// A
  	virtual CoinPackedMatrix getFirstStageConstraints() { return Amat; }
  	// W
  	virtual CoinPackedMatrix getSecondStageConstraints(int scen) {loadLocalScenData(scen); return Wmat[scen]; }
  	// T
  	virtual CoinPackedMatrix getLinkingConstraints(int scen) {loadLocalScenData(scen); return Tmat[scen]; }
  
  	// QA
//  	virtual CoinPackedMatrix getFirstStageHessian(){return QAmat;}
  	// QW
  	virtual CoinPackedMatrix getSecondStageHessian(int scen){loadLocalScenData(scen);return QWmat[scen];}
  	// QT
//  	virtual CoinPackedMatrix getSecondStageCrossHessian(int scen){loadLocalScenData(scen);return QTmat[scen];}

protected:
	struct scenData {
		std::vector<double> collb, colub, rowlb, rowub, objLin;		
		std::vector<std::string> rownames, colnames;
		bool didLoad;
		scenData() : didLoad(false) {}
		void initialize(int nvar, int ncons);
		void copyFrom(scenData s2);
	};
	
	virtual void loadLocalScenData(int scen);
	
	int mype_;

	double ObjScale;
	int nScenarios_, nFirstStageVars_, nFirstStageCons_;
	std::vector<int> nSecondStageVars_, nSecondStageCons_;

	scenData firstStageData;
    std::vector<scenData> localData;
	
	CoinPackedMatrix Amat, QAmat; 
	std::vector<CoinPackedMatrix> Tmat, QTmat, Wmat, QWmat;

	
private:
  void parseZeroData();

  DCOPFLOW *dcopf;
    
};


#endif

