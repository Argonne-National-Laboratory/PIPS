/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLCROSSHESSIANINPUT__HPP
#define NLCROSSHESSIANINPUT__HPP

#include <map>
#include <vector>
#include <string>

#include "stochasticInput.hpp"
#include "AmplData_NL.hpp"
#include "global_var.h"

struct ASL_pfgh;
struct SufDecl;
struct SufDesc;

#include "mpi.h"
#include "../PIPS-NLP/global_var.h"
//  solve problem as general stochastic programming problem
//  min f(x_0) + \sum f(x_i)
//  st    g(x_i)  = 0
//
//  we do split/reorder the Hessian and Jacobian in this case.

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

class amplGenStochInput: public stochasticInput {
public:
	amplGenStochInput() {
	}
	amplGenStochInput(const std::string &datarootname,
			int overrideScenarioNumber = 0, MPI_Comm comm = MPI_COMM_WORLD);
	~amplGenStochInput();

	virtual int nScenarios() {
		return nScenarios_;
	}
	virtual int nFirstStageVars() {
		return nFirstStageVars_;
	}
	virtual int nFirstStageCons() {
		return nFirstStageCons_;
	}
	virtual int nSecondStageVars(int scen) {
		return nSecondStageVars_[scen];
	}
	virtual int nSecondStageCons(int scen) {
		return nSecondStageCons_[scen];
	}

	virtual std::vector<double> getFirstStageColLB() {
		return firstStageData.collb;
	}
	virtual std::vector<double> getFirstStageColUB() {
		return firstStageData.colub;
	}
	virtual std::vector<double> getFirstStageObj() {
		return firstStageData.objGrad;
	}
	virtual std::vector<std::string> getFirstStageColNames() {
		return firstStageData.colnames;
	}
	virtual std::vector<double> getFirstStageRowLB() {
		return firstStageData.rowlb;
	}
	virtual std::vector<double> getFirstStageRowUB() {
		return firstStageData.rowub;
	}
	virtual std::vector<std::string> getFirstStageRowNames() {
		return firstStageData.rownames;
	}

	virtual bool isFirstStageColInteger(int col) {
		return false;
	}

	virtual std::vector<double> getSecondStageColLB(int scen) {
		loadLocalNLdata(scen);
		return localData[scen].collb;
	}
	virtual std::vector<double> getSecondStageColUB(int scen) {
		loadLocalNLdata(scen);
		return localData[scen].colub;
	}
	virtual std::vector<double> getSecondStageObj(int scen) {
		loadLocalNLdata(scen);
		return localData[scen].objGrad;
	}
	virtual std::vector<std::string> getSecondStageColNames(int scen) {
		loadLocalNLdata(scen);
		return localData[scen].colnames;
	}
	virtual std::vector<double> getSecondStageRowLB(int scen) {
		loadLocalNLdata(scen);
		return localData[scen].rowlb;
	}
	virtual std::vector<double> getSecondStageRowUB(int scen) {
		loadLocalNLdata(scen);
		return localData[scen].rowub;
	}
	virtual std::vector<std::string> getSecondStageRowNames(int scen) {
		loadLocalNLdata(scen);
		return localData[scen].rownames;
	}
	virtual double scenarioProbability(int scen) {
		return 1.0 / nScenarios_;
	}
	virtual bool isSecondStageColInteger(int scen, int col) {
		return false;
	}

	virtual CoinPackedMatrix getFirstStageConstraints() {
		MESSAGE("getFirstStageConstraints");
		IF_VERBOSE_DO( Amat.dumpMatrix(); );
		return Amat;
	}
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) {
		loadLocalNLdata(scen);
		MESSAGE("getSecondStageConstraints scen" << scen);
		IF_VERBOSE_DO( Wmat[scen].dumpMatrix(); );
		return Wmat[scen];
	}
	virtual CoinPackedMatrix getLinkingConstraints(int scen) {
		loadLocalNLdata(scen);
		MESSAGE("getLinkingConstraints scen" << scen);
		IF_VERBOSE_DO( Tmat[scen].dumpMatrix(); );
		return Tmat[scen];
	}

	virtual CoinPackedMatrix getFirstStageHessian() {
		MESSAGE("getFirstStageHessian");
		IF_VERBOSE_DO( QAmat.dumpMatrix(); );
		return QAmat;
	}
	// Q_i
	virtual CoinPackedMatrix getSecondStageHessian(int scen) {
		loadLocalNLdata(scen);
		MESSAGE("getSecondStageHessian scen" << scen);
		IF_VERBOSE_DO( QWmat[scen].dumpMatrix(); );
		return QWmat[scen];
	}
	// column-oriented, \hat Q_i
	// Note: this has the second-stage variables on the rows and first-stage on the columns
	virtual CoinPackedMatrix getSecondStageCrossHessian(int scen) {
		loadLocalNLdata(scen);
		MESSAGE("getSecondStageCrossHessian scen" << scen);
		IF_VERBOSE_DO( QTmat[scen].dumpMatrix(); );
		return QTmat[scen];
	}

	virtual bool scenarioDimensionsEqual() {
		return dimsEqual;
	}
	virtual bool onlyBoundsVary() {
		return false;
	} // no easy way to check
	virtual bool allProbabilitiesEqual() {
		return true;
	}
	virtual bool continuousRecourse() {
		return true;
	}

	virtual void doNetworkPart(const int scen, AmplSuffix* amplSuffix);

protected:
	bool dimsEqual;
	bool haveRootNL;

	double ObjScale;

	int nScenarios_, nFirstStageVars_, nFirstStageCons_;
	std::vector<int> nSecondStageVars_, nSecondStageCons_;
	std::vector<int*> LocGloVarIdx;

	std::vector<std::map<int, int> > LocLocVarMap;
	std::vector<std::map<int, int> > LocGloVarMap;

	std::vector<std::map<int, int> > LocWmatJacGoffMap;
	std::vector<std::map<int, int> > LocTmatJacGoffMap;

	std::vector<std::map<int, int> > LocQAmatHesGoffMap;
	std::vector<std::map<int, int> > LocQWmatHesGoffMap;
	std::vector<std::map<int, int> > LocQTmatHesGoffMap;

	std::vector<int> decisionVarDim;
	std::vector<int*> schurVarConIDinNL;

	std::vector<int*> NR_partIDX_var;
	std::vector<int*> NR_partIDX_con;

	int mype_;

	std::vector<AmplData_NL> localData;
	AmplData_NL firstStageData;

	std::vector<ASL_pfgh*> asl_i;
	ASL_pfgh* asl_0;

	CoinPackedMatrix Amat, QAmat;
	std::vector<CoinPackedMatrix> Tmat, QTmat, Wmat, QWmat;

	virtual void loadLocalNLdata(int scen);
	virtual void splitMatrices(const int scen);

	// no copying
	amplGenStochInput(const amplGenStochInput&);
	amplGenStochInput& operator=(const amplGenStochInput&);

////////////////////////////////////////////////////////////////////////////////////

	virtual void getRowMap(const int scen);

	int nnzEqJac1st, nnzIneqJac1st;
	int* amplRowMap1st;

	std::vector<int> nnzEqJac2nd, nnzIneqJac2nd;
	std::vector<int*> amplRowMap2nd;

////////////////////////////////////////////////////////////////////////////////////

	virtual void getJacGoffMap(const int scen);

	int nnzA1st, nnzC1st;
	std::vector<int> nnzA2nd, nnzC2nd;


virtual int nLinkCons(){return 0;};
virtual int nLinkECons(){return 0;};
virtual int nLinkICons(){return 0;};


////////////////////////////////////////////////////////////////////////////////////

	int nnzALink1st, nnzCLink1st, nnzALoc1st, nnzCLoc1st;
	std::vector<int> nnzALink2nd, nnzCLink2nd, nnzALoc2nd, nnzCLoc2nd;

	std::map<int, int> JacALinkGoff1st, JacCLinkGoff1st;
	std::map<int, int> JacALocGoff1st, JacCLocGoff1st;

	std::vector<std::map<int, int> > JacALinkGoff2nd, JacCLinkGoff2nd;
	std::vector<std::map<int, int> > JacALocGoff2nd, JacCLocGoff2nd;

	int nnzQ1st;
	std::vector<int> nnzQ2nd, nnzQCross2nd;

////////////////////////////////////////////////////////////////////////////////////	

	std::vector<int> starts1stSt;
	friend class sNlpInfoFromNL;
};

#endif

