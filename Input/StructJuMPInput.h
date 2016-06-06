/* PIPS-NLP                                                         	*
 * Authors: Feng Qiang                      		*
 * (C) 2015 Argonne National Laboratory			*/

#ifndef _STRUCTJUMP_INPUT_H_
#define _STRUCTJUMP_INPUT_H_

#include <map>
#include <vector>
#include <string>

#include "stochasticInput.hpp"
#include "../PIPS-NLP/Drivers/parallelPipsNlp_C_Callback.h"
#include "mpi.h"

class StructJuMPInput : public stochasticInput {

protected:
	std::map<int, int> nvar_map;
	std::map<int, int> ncon_map;
	std::map<int, std::vector<double> > collb_map;
	std::map<int, std::vector<double> > colub_map;
	std::map<int, std::vector<double> > rowlb_map;
	std::map<int, std::vector<double> > rowub_map;
	std::map<int, int> e_ncon_map;
	std::map<int, int> i_ncon_map;

	CoinPackedMatrix amat;
	std::map<int, CoinPackedMatrix> wmat_map;
	std::map<int, CoinPackedMatrix> tmat_map;

	CoinPackedMatrix qamat;
	std::map<int, CoinPackedMatrix> qwmat_map;
	std::map<int, CoinPackedMatrix> qtmat_map;

	std::vector<double> linklb;
	std::vector<double> linkub;
        int mlink;
        int e_ml;
        int i_ml;
	std::map<int,CoinPackedMatrix> Emat_map; // 0 is corresponding to the first stage, 1 to scenario 1. linking jacobian 
public:
	PipsNlpProblemStruct* prob;

public:
	StructJuMPInput(PipsNlpProblemStruct* prob);
	~StructJuMPInput();

	void get_prob_info(int nodeid);
	virtual int nScenarios();
	virtual int nFirstStageVars();
	virtual int nFirstStageCons();
	virtual int nSecondStageVars(int scen);
	virtual int nSecondStageCons(int scen);

	virtual std::vector<double> getFirstStageColLB();
	virtual std::vector<double> getFirstStageColUB();
	virtual std::vector<double> getFirstStageObj();
	virtual std::vector<std::string> getFirstStageColNames();
	virtual std::vector<double> getFirstStageRowLB();
	virtual std::vector<double> getFirstStageRowUB();
	virtual std::vector<std::string> getFirstStageRowNames();
	virtual bool isFirstStageColInteger(int col);

	virtual std::vector<double> getSecondStageColLB(int scen);
	virtual std::vector<double> getSecondStageColUB(int scen);
	// objective vector, already multiplied by probability
	virtual std::vector<double> getSecondStageObj(int scen);
	virtual std::vector<std::string> getSecondStageColNames(int scen);
	virtual std::vector<double> getSecondStageRowUB(int scen);
	virtual std::vector<double> getSecondStageRowLB(int scen);
	virtual std::vector<std::string> getSecondStageRowNames(int scen);
	virtual double scenarioProbability(int scen);
	virtual bool isSecondStageColInteger(int scen, int col);


	//to deal with constraints links first stage variables and second stage from scenario 1, 2 ...s 
	virtual int nLinkCons(){ return mlink;};
        virtual int nLinkECons(){ return e_ml;};
        virtual int nLinkICons(){ return i_ml;};
        virtual std::vector<double> getLinkRowLB();
        virtual std::vector<double> getLinkRowUB();
	virtual CoinPackedMatrix getLinkMatrix(int nodeid);

	// returns the column-oriented first-stage constraint matrix (A matrix)
	virtual CoinPackedMatrix getFirstStageConstraints();
	// returns the column-oriented second-stage constraint matrix (W matrix)
	virtual CoinPackedMatrix getSecondStageConstraints(int scen);
	// returns the column-oriented matrix linking the first-stage to the second (T matrix)
	virtual CoinPackedMatrix getLinkingConstraints(int scen);

	virtual CoinPackedMatrix getFirstStageHessian();
	// Q_i
	virtual CoinPackedMatrix getSecondStageHessian(int scen);
	// column-oriented, \hat Q_i
	// Note: this has the second-stage variables on the rows and first-stage on the columns
	virtual CoinPackedMatrix getSecondStageCrossHessian(int scen);

	// some problem characteristics that could be helpful to know

	// all scenarios have the same number of variables and constraints
	virtual bool scenarioDimensionsEqual();
	// constraint (and hessian) matrices are identical for each scenario,
	// column and row bounds and objective are allowed to vary
	virtual bool onlyBoundsVary();
	// all scenarios equally likely
	virtual bool allProbabilitiesEqual();
	// all second-stage variables continuous
	virtual bool continuousRecourse();
};

#endif
