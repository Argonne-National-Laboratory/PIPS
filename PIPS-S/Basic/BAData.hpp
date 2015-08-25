#ifndef BADATA_HPP
#define BADATA_HPP

#include <boost/shared_ptr.hpp>

#include "BA.hpp"
#include "BAVector.hpp"
#include "stochasticInput.hpp"
#include "BALPSolverInterface.hpp"
#include "CoinPackedVector.hpp"

// stores actual problem data
// collect data into standard format for solver (introducing slacks),
// and perform linear algebra (PRICE) here

class BAData {
public:
	BAData(stochasticInput &input, BAContext &ctx);
	BAData(const BAData&);
	~BAData();

	void getCol(sparseBAVector &v, BAIndex i) const;
	void addColToVec(sparseBAVector &v, BAIndex i, double mult) const;
	
	// pick out nonbasic elements from "in" and multiply them with nonbasic columns
	// effectively same as multiplying with whole matrix assuming basic "in" are zero
	// a bit strange for dense input and sparse output, but this is the most convenient presently
	void multiply(const denseBAVector &in, sparseBAVector &out, const BAFlagVector<variableState>&) const; 

	// multiply entire constraint matrix transpose by "in"
	// taking linear combinations of the rows
	void multiplyT(const sparseBAVector &in, sparseBAVector &out) const;


	void addRow(const CoinPackedVectorBase& elts1, const CoinPackedVectorBase &elts2, int scen, double lb = -COIN_DBL_MAX, double ub = COIN_DBL_MAX);

	void addFirstStageRows(const std::vector<CoinPackedVector*> &v1, std::vector<double> &lb, std::vector<double> &ub, int nRows);
	
	void addSecondStageConsecutiveRows(const std::vector<CoinPackedVector*> &elts1, const std::vector<CoinPackedVector*> &elts2, int scenario, std::vector<double> &lb, std::vector<double> &ub, int nRows);

	int addFirstStageColumn(double lb, double ub, double c);

	void addFirstStageRow(const CoinPackedVectorBase& elts1, double lb, double ub);

	int addSecondStageColumn(int scen,double lb, double ub, double cobj);

	void deleteLastFirstStageRows(int nRows);

	void deleteLastSecondStageConsecutiveRows(int scenario, int nRows);

	void deleteLastFirstStageColumns(int nCols);

	const CoinShallowPackedVector retrieveARow(int index) const;

	const CoinShallowPackedVector retrieveWRow(int index,int scen) const;

	const CoinShallowPackedVector retrieveTRow(int index,int scen) const;

	const CoinShallowPackedVector retrieveACol(int index) const;

	const CoinShallowPackedVector retrieveWCol(int index,int scen) const;

	const CoinShallowPackedVector retrieveTCol (int index,int scen) const;


	// make members easily accessable for use in solver
	BADimensionsSlacks dims;
	denseBAVector l,u,c;
	BAFlagVector<constraintType> vartype;
	BAFlagVector<std::string> names;
	boost::shared_ptr<CoinPackedMatrix> Acol; 
	boost::shared_ptr<CoinPackedMatrix> Arow;
	std::vector<boost::shared_ptr<CoinPackedMatrix> > Tcol;
	std::vector<boost::shared_ptr<CoinPackedMatrix> > Trow;
	std::vector<boost::shared_ptr<CoinPackedMatrix> > Wcol; 
	std::vector<boost::shared_ptr<CoinPackedMatrix> > Wrow; 
	BAContext &ctx;

protected:
	bool onlyBoundsVary;
	mutable CoinIndexedVector out1Send; // buffer for multiplyT
	//BAFlagVector<variableState> stateCol, stateRow;

};

#endif
