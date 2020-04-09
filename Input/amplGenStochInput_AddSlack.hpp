/* PIPS-NLP                                                         	*
 * Authors: Nai-Yuan Chiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef NLCROSSHESSIANINPUT_C__HPP
#define NLCROSSHESSIANINPUT_C__HPP

#include <map>
#include <vector>
#include <string>

#include "stochasticInput.hpp"
#include "amplGenStochInput.hpp"
#include "AmplData_NL.hpp"

#include "mpi.h"


struct ASL_pfgh;
struct SufDecl;
struct SufDesc;


//  solve problem as 
//  min 1/s f(x_s)
//  st    g(x_s)  = 0
//         M_sx_s = x_0
//
//  this input is differenent from the amplGenStochInput one. 
//  In amplGenStochInput, we do split/reorder the Hessian and Jacobian in this case.
//  Here we need to add dummy variable/constraint
//
//
// Jacobian looks like 
//  	W_1           		T_1	
//		W_2           	T_1	
//			...
//				W_n T_n
//					  A
//
//
// Hessian looks like 
//  	QW_1           			QT_1	
//			QW_2           	QT_1	
//				      ...
//				     QW_n	QT_n
//					  	  QA


class amplGenStochInput_AddSlack : public amplGenStochInput {
public:
	amplGenStochInput_AddSlack(){}
	amplGenStochInput_AddSlack(const std::string &datarootname, int overrideScenarioNumber = 0, 
						MPI_Comm comm = MPI_COMM_WORLD);
	~amplGenStochInput_AddSlack() { }

	virtual int nScenarios() { return nScenarios_; }
	virtual int nFirstStageVars() { return nFirstStageVars_; }
	virtual int nFirstStageCons() { return nFirstStageCons_; }
	virtual int nSecondStageVars(int scen) { return nSecondStageVars_[scen]; }
	virtual int nSecondStageCons(int scen) { return nSecondStageCons_[scen]; }

	virtual std::vector<double> getFirstStageColLB(){ return firstStageData.collb; }
	virtual std::vector<double> getFirstStageColUB(){ return firstStageData.colub; } 
	virtual std::vector<double> getFirstStageObj(){ return firstStageData.objGrad; }
	virtual std::vector<std::string> getFirstStageColNames(){return firstStageData.colnames;}
	virtual std::vector<double> getFirstStageRowLB(){ return firstStageData.rowlb; }
	virtual std::vector<double> getFirstStageRowUB(){ return firstStageData.rowub; }
	virtual std::vector<std::string> getFirstStageRowNames(){ return firstStageData.rownames; }
	
	virtual bool isFirstStageColInteger(int col) { return false; }
	virtual bool isFirstStageColBinary(int col) { return false; }

	virtual std::vector<double> getSecondStageColLB(int scen){loadLocalNLdata(scen);return localData[scen].collb;}
	virtual std::vector<double> getSecondStageColUB(int scen){loadLocalNLdata(scen);return localData[scen].colub;}
	virtual std::vector<double> getSecondStageObj(int scen){loadLocalNLdata(scen);return localData[scen].objGrad;}
	virtual std::vector<std::string> getSecondStageColNames(int scen){loadLocalNLdata(scen);return localData[scen].colnames;}
	virtual std::vector<double> getSecondStageRowLB(int scen){loadLocalNLdata(scen);return localData[scen].rowlb;}
	virtual std::vector<double> getSecondStageRowUB(int scen){loadLocalNLdata(scen);return localData[scen].rowub;}
	virtual std::vector<std::string> getSecondStageRowNames(int scen){loadLocalNLdata(scen);return localData[scen].rownames;}
	virtual double scenarioProbability(int scen) { return 1.0/nScenarios_; }
	virtual bool isSecondStageColInteger(int scen, int col) { return false; }
	virtual bool isSecondStageColBinary(int scen, int col) { return false; }


	virtual CoinPackedMatrix getFirstStageConstraints() { return Amat; }
	virtual CoinPackedMatrix getSecondStageConstraints(int scen) {loadLocalNLdata(scen); return Wmat[scen]; }
	virtual CoinPackedMatrix getLinkingConstraints(int scen) {loadLocalNLdata(scen); return Tmat[scen]; }



	virtual CoinPackedMatrix getFirstStageHessian(){return QAmat;}
	// Q_i
	virtual CoinPackedMatrix getSecondStageHessian(int scen){loadLocalNLdata(scen);return QWmat[scen];}
	// column-oriented, \hat Q_i
	// Note: this has the second-stage variables on the rows and first-stage on the columns
	virtual CoinPackedMatrix getSecondStageCrossHessian(int scen){loadLocalNLdata(scen);return QTmat[scen];}


	

	virtual bool scenarioDimensionsEqual() { return dimsEqual; }
	virtual bool onlyBoundsVary() { return false; } // no easy way to check
	virtual bool allProbabilitiesEqual() { return true; }
	virtual bool continuousRecourse() { return true; }

	virtual void doNetworkPart(const int scen, AmplSuffix* amplSuffix);

	virtual void loadLocalNLdata(int scen);
	virtual void splitMatrices(const int scen);

	virtual void getRowMap(const int scen);

	virtual void getJacGoffMap(const int scen);


	friend class sNlpInfoFromNL;


	
};


#endif


