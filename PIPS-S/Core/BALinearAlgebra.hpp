#ifndef BALPLINALG_HPP
#define BALPLINALG_HPP

#include "BAData.hpp"
#include "CoinBALPFactorization.hpp"

// treat first stage as sparse and solve with CoinUtils

class BALinearAlgebra  {
public:
	BALinearAlgebra(const BAData& d);
	~BALinearAlgebra();

	void reinvert(const BAFlagVector<variableState> &);
	void ftran(sparseBAVector &);
	void btran(sparseBAVector &);

protected:
	virtual void reinvertFirstStage(const std::vector<int> &basicCols1);
	int nbasic1;

	std::vector<CoinBALPFactorization*> f;
	const BAData &data;

	// first-stage matrix
	CoinPackedMatrix m;
	CoinIndexedVector region1;
	CoinBALPFactorization f0;
	// triplets for first stage
	std::vector<double> myElements, allElements;
	std::vector<int> myIndicesRow, allIndicesRow, myIndicesColumn, allIndicesColumn;
	std::vector<int> rowsPerProc; // number of rows that come into the first stage from each process
	int rowsIn; // total number of rows that come into the first stage
	int maxRowsIn; // max number of rows from any process (used for packing)
	int myOffset; // row index where this process's scenarios start

	std::vector<double> ftranSend, ftranRecv; // send recv buffers for Allgather during FTRAN
	std::vector<double> btranSend;

	std::vector<CoinIndexedVector> regions;

	double tmp_t, tmp2_t;

};

#endif
