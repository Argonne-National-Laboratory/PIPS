#include "CbcBALPInterface.hpp"
#include <sstream>
#include <fstream>

using namespace std;


namespace{
template<typename T1, typename T2> void concatenateAll(stochasticInput &data, 
			T1 &out, 
			T2 (stochasticInput::*first)(),
			T2 (stochasticInput::*second)(int)) {

	const T2& arr1 = (data.*first)();
	for (unsigned k = 0; k < arr1.size(); k++) {
		out[k] = arr1[k];
	}
	unsigned r = arr1.size();
	int nscen = data.nScenarios();
	for (int i = 0; i < nscen; i++) {
		const T2& arr2 = (data.*second)(i);
		for (unsigned k = 0; k < arr2.size(); k++) {
			out[r++] = arr2[k];
		}
	}

}
}

CbcBALPInterface::CbcBALPInterface(stochasticInput &input, BAContext &ctx) : dims(input,ctx) {
	assert(ctx.nprocs() == 1);

	int const totalVar = dims.totalVars();
	int const totalCons = dims.totalCons();
	int const nScenarios = dims.numScenarios();
	
	CoinBigIndex totalNnz = input.getFirstStageConstraints().getNumElements();
	for (int scen = 0; scen < nScenarios; scen++) {
		totalNnz += input.getSecondStageConstraints(scen).getNumElements();
		totalNnz += input.getLinkingConstraints(scen).getNumElements();
	}

	// CoinPackedMatrix takes ownership of these, so we don't free them
	CoinBigIndex *starts = new CoinBigIndex[totalVar+1];
	double *elts = new double[totalNnz];
	int *rowIdx = new int[totalNnz];
	
	CoinBigIndex nnz = 0, start, end;
	// put first-stage variables first, as is customary
	CoinPackedMatrix const &Amat = input.getFirstStageConstraints();
	int const *Aidx = Amat.getIndices();
	double const *Aelts = Amat.getElements();
	int const nFirstStageVars = input.nFirstStageVars();
	int const nFirstStageCons = input.nFirstStageCons();
	CoinPackedMatrix const &Tmat0 = input.getLinkingConstraints(ctx.localScenarios().at(1));

	for (int c = 0; c < nFirstStageVars; c++) {
		starts[c] = nnz;
		start = Amat.getVectorFirst(c);
		end = Amat.getVectorLast(c);
		for (CoinBigIndex j = start; j < end; j++) {
			elts[nnz] = Aelts[j];
			rowIdx[nnz++] = Aidx[j];
		}
		int rowOffset = nFirstStageCons;
		for (int scen = 0; scen < nScenarios; scen++) {
			// this could be very inefficient if we make another copy of the T matrix for each column
			CoinPackedMatrix const &Tmat = (input.onlyBoundsVary()) ? Tmat0 : input.getLinkingConstraints(scen);
			int const *Tidx = Tmat.getIndices();
			double const *Telts = Tmat.getElements();
			start = Tmat.getVectorFirst(c);
			end = Tmat.getVectorLast(c);
			for (CoinBigIndex j = start; j < end; j++) {
				elts[nnz] = Telts[j];
				rowIdx[nnz++] = Tidx[j]+rowOffset;
			}
			rowOffset += input.nSecondStageCons(scen);
		}
	
	}

	// now W blocks
	int rowOffset = nFirstStageCons;
	int colOffset = nFirstStageVars;
	for (int scen = 0; scen < nScenarios; scen++) {
		int nSecondStageVars = input.nSecondStageVars(scen);
		CoinPackedMatrix const &Wmat = input.getSecondStageConstraints(scen);
		int const *Widx = Wmat.getIndices();
		double const *Welts = Wmat.getElements();

		for (int c = 0; c < nSecondStageVars; c++) {
			starts[colOffset++] = nnz;
			start = Wmat.getVectorFirst(c);
			end = Wmat.getVectorLast(c);
			for (CoinBigIndex j = start; j < end; j++) {
				elts[nnz] = Welts[j];
				rowIdx[nnz++] = Widx[j]+rowOffset;
			}
		}
		rowOffset += input.nSecondStageCons(scen);
	}
	starts[totalVar] = nnz;
	assert(nnz == totalNnz);

	CoinPackedMatrix constr;
	int *lens = 0;
	constr.assignMatrix(true,totalCons,totalVar,totalNnz,
		elts, rowIdx, starts, lens);

	double *collb = new double[totalVar], *colub = new double[totalVar];
	double *obj = new double[totalVar];
	double *rowlb = new double[totalCons], *rowub = new double[totalCons];

	concatenateAll(input, collb, &stochasticInput::getFirstStageColLB,
			&stochasticInput::getSecondStageColLB);
	
	concatenateAll(input, colub, &stochasticInput::getFirstStageColUB,
			&stochasticInput::getSecondStageColUB);

	concatenateAll(input, obj, &stochasticInput::getFirstStageObj,
			&stochasticInput::getSecondStageObj);
	
	concatenateAll(input, rowlb, &stochasticInput::getFirstStageRowLB,
			&stochasticInput::getSecondStageRowLB);

	concatenateAll(input, rowub, &stochasticInput::getFirstStageRowUB,
			&stochasticInput::getSecondStageRowUB);
	
	// don't copy names for now

	model.loadProblem(constr,collb,colub,obj,rowlb,rowub);

	for (int i = 0; i < nFirstStageVars; i++) {
		if (input.isFirstStageColInteger(i)) {
			model.setInteger(i);
		}
	}
	int offset = nFirstStageVars;
	for (int scen = 0; scen < nScenarios; scen++) {
		int nSecondStageVars = input.nSecondStageVars(scen);
		for (int i = 0; i < nSecondStageVars; i++) {
			if (input.isSecondStageColInteger(scen,i)) {
				model.setInteger(i+offset);
			}
		}
		offset += nSecondStageVars;
	}

	delete [] collb;
	delete [] colub;
	delete [] obj;
	delete [] rowlb;
	delete [] rowub;

	model.messageHandler()->setLogLevel(0);
	cbcm.reset(new CbcModel(model));
	CbcMain0(*cbcm);
	cbcm->messageHandler()->setLogLevel(0);

}

void CbcBALPInterface::go() {
	
	const char * argv2[]={"","-solve","-quit"};
	//cbcm->setMaximumNodes(1);
	cbcm->setNumberThreads(2);
	CbcMain1(3,argv2,*cbcm);
	//cbcm->branchAndBound();


}

solverState CbcBALPInterface::getStatus() const {

	if (cbcm->isProvenOptimal()) {
		return Optimal;
	} 
	if (cbcm->isProvenInfeasible()) {
		return ProvenInfeasible;
	}
	assert(!cbcm->isAbandoned());

	return Initialized;
}

/*
static variableState clpStatusToState(ClpSimplex::Status s) {
	switch (s) {
		case ClpSimplex::basic:
			return Basic;
		case ClpSimplex::atLowerBound:
			return AtLower;
		case ClpSimplex::atUpperBound:
			return AtUpper;
		case ClpSimplex::isFixed:
			return AtLower;
		case ClpSimplex::isFree:
			return AtLower;
		case ClpSimplex::superBasic:
			printf("Super basic?\n");
			return AtLower;
	}
	assert(0);
	return AtLower;
	
}

static ClpSimplex::Status stateToClpStatus(variableState s) {
	switch (s) {
		case Basic:
			return ClpSimplex::basic;
		case AtLower:
			return ClpSimplex::atLowerBound;
		case AtUpper:
			return ClpSimplex::atUpperBound;
	}
	assert(0);
	return ClpSimplex::atLowerBound;
	
}

variableState CbcBALPInterface::getFirstStageColState(int idx) const {
	return clpStatusToState(model.getColumnStatus(idx));
}

variableState CbcBALPInterface::getFirstStageRowState(int idx) const {
	return clpStatusToState(model.getRowStatus(idx));
}

variableState CbcBALPInterface::getSecondStageColState(int scen, int idx) const {
	int offset = dims.numFirstStageVars();
	for (int i = 0; i < scen; i++) offset += dims.numSecondStageVars(i);
	return clpStatusToState(model.getColumnStatus(offset+idx));
}

variableState CbcBALPInterface::getSecondStageRowState(int scen, int idx) const {
	int offset = dims.numFirstStageCons();
	for (int i = 0; i < scen; i++) offset += dims.numSecondStageCons(i);
	return clpStatusToState(model.getRowStatus(offset+idx));
}

void CbcBALPInterface::setFirstStageColState(int idx,variableState s) {
	model.setColumnStatus(idx, stateToClpStatus(s));
}

void CbcBALPInterface::setFirstStageRowState(int idx,variableState s) {
	model.setRowStatus(idx,stateToClpStatus(s));
}

void CbcBALPInterface::setSecondStageColState(int scen, int idx,variableState s) {
	int offset = dims.numFirstStageVars();
	for (int i = 0; i < scen; i++) offset += dims.numSecondStageVars(i);
	model.setColumnStatus(offset+idx,stateToClpStatus(s));
}

void CbcBALPInterface::setSecondStageRowState(int scen, int idx,variableState s) {
	int offset = dims.numFirstStageCons();
	for (int i = 0; i < scen; i++) offset += dims.numSecondStageCons(i);
	model.setRowStatus(offset+idx,stateToClpStatus(s));
}


vector<double> CbcBALPInterface::getFirstStagePrimalColSolution() const {
	const double *sol = model.primalColumnSolution();
	int nvar1 = dims.numFirstStageVars();
	return vector<double>(sol,sol+nvar1);
}

vector<double> CbcBALPInterface::getSecondStagePrimalColSolution(int scen) const {
	const double *sol = model.primalColumnSolution();
	int offset = dims.numFirstStageVars();
	for (int i = 0; i < scen; i++) {
		offset += dims.numSecondStageVars(i);
	}
	return vector<double>(sol+offset,sol+offset+dims.numSecondStageVars(scen));
}

vector<double> CbcBALPInterface::getFirstStageDualColSolution() const {
	const double *sol = model.dualColumnSolution();
	int nvar1 = dims.numFirstStageVars();
	return vector<double>(sol,sol+nvar1);
}

vector<double> CbcBALPInterface::getSecondStageDualColSolution(int scen) const {
	const double *sol = model.dualColumnSolution();
	int offset = dims.numFirstStageVars();
	for (int i = 0; i < scen; i++) {
		offset += dims.numSecondStageVars(i);
	}
	return vector<double>(sol+offset,sol+offset+dims.numSecondStageVars(scen));
}

static const char* statusString(ClpSimplex::Status s) {
	if (s == ClpSimplex::basic) {
		return "Basic";
	} else if (s == ClpSimplex::atLowerBound || s == ClpSimplex::isFixed) {
		return "AtLower";
	} else if (s == ClpSimplex::atUpperBound) {
		return "AtUpper";
	} else if (s == ClpSimplex::isFree) {
		printf("Got a free variable, why?\n");
		return "AtLower";
	} else if (s == ClpSimplex::superBasic) {
		printf("Super basic??\n");
		return "AtLower";
	} else {
		printf("bad value: %d\n", (int)s);
		assert(0);
		return "error";
	}
}

void CbcBALPInterface::writeStatus(const std::string &filebase) {
	
	int nscen = dims.numScenarios();
	int nvar1real = dims.numFirstStageVars();
	int ncons1 = dims.numFirstStageCons();
	//int nvar2real = dims.numSecondStageVars(0);
	//int ncons2 = dims.numSecondStageCons(0);
	int totalRows = dims.totalCons();

	const double *sol = model.primalColumnSolution();
	const double *rowSol = model.primalRowSolution();

	// TODO: think about extracting DSE weights

	int nbasic = 0;
	int roffset = ncons1;
	int coffset = nvar1real;
	for (int k = 0; k < nscen; k++) {
		stringstream fname;
		fname << filebase << k+1;
		ofstream f(fname.str().c_str());
		int nbasicThis = 0;
		int nvar2real = dims.numSecondStageVars(k);
		int ncons2 = dims.numSecondStageCons(k);
		for (int i = 0; i < nvar2real; i++) {
			int idx = i + coffset;
			if (model.getColumnStatus(idx) == ClpSimplex::basic) {
				nbasicThis++; nbasic++;
			} else if (model.getColumnStatus(idx) == ClpSimplex::superBasic) {
				double s = sol[idx], l = model.getColLower()[idx], u = model.getColUpper()[idx];
				printf("Got superbasic: at: %e lb: %e ub: %e\n",s,l,u);
				if (abs(s-l) < abs(s-u)) {
					printf("set to lower\n");
					model.setColumnStatus(idx,ClpSimplex::atLowerBound);
				} else {
					printf("set to upper\n");
					model.setColumnStatus(idx,ClpSimplex::atUpperBound);
				}
			}	
		}
		for (int i = 0; i < ncons2; i++) {
			int idx = i+roffset;
			if (model.getRowStatus(idx) == ClpSimplex::basic) {
				nbasicThis++; nbasic++;
			} else if (model.getRowStatus(idx) == ClpSimplex::superBasic) {
				double s = rowSol[idx], l = model.getRowLower()[idx], u = model.getRowLower()[idx];
				printf("Got superbasic: at: %e lb: %e ub: %e\n",s,l,u);
				if (abs(s-l) < abs(s-u)) {
					printf("set to lower\n");
					model.setRowStatus(idx,ClpSimplex::atLowerBound);
				} else {
					printf("set to upper\n");
					model.setRowStatus(idx,ClpSimplex::atUpperBound);
				}
			}
		}
		f << nbasicThis << " BasisOnly\n";
		for (int i = 0; i < nvar2real; i++) {
			f << i << " " << statusString(model.getColumnStatus(i+coffset)) << "\n";	
		}
		for (int i = 0; i < ncons2; i++) {
			f << i + nvar2real << " " << statusString(model.getRowStatus(i+roffset)) << "\n";
		}

		f.close();
		roffset += ncons2;
		coffset += nvar2real;
	}
	stringstream fname;
	fname << filebase << 0;
	ofstream f(fname.str().c_str());

	int nbasicThis = 0;
	for (int i = 0; i < nvar1real; i++) {
		int idx = i;
		if (model.getColumnStatus(idx) == ClpSimplex::basic) {
			nbasicThis++; nbasic++;
		} else if (model.getColumnStatus(idx) == ClpSimplex::superBasic) {
			double s = sol[idx], l = model.getColLower()[idx], u = model.getColUpper()[idx];
			printf("Got superbasic: at: %e lb: %e ub: %e\n",s,l,u);
			if (abs(s-l) < abs(s-u)) {
				printf("set to lower\n");
				model.setColumnStatus(idx,ClpSimplex::atLowerBound);
			} else {
				printf("set to upper\n");
				model.setColumnStatus(idx,ClpSimplex::atUpperBound);
			}
		}	
	}
	for (int i = 0; i < ncons1; i++) {
		int idx = i;
		if (model.getRowStatus(idx) == ClpSimplex::basic) {
			nbasicThis++; nbasic++;
		} else if (model.getRowStatus(idx) == ClpSimplex::superBasic) {
			double s = rowSol[idx], l = model.getColLower()[idx], u = model.getColUpper()[idx];
			printf("Got superbasic: at: %e lb: %e ub: %e\n",s,l,u);
			if (abs(s-l) < abs(s-u)) {
				printf("set to lower\n");
				model.setRowStatus(idx,ClpSimplex::atLowerBound);
			} else {
				printf("set to upper\n");
				model.setRowStatus(idx,ClpSimplex::atUpperBound);
			}
		}
	}
	f << nbasicThis << " BasisOnly\n";
	for (int i = 0; i < nvar1real; i++) {
		f << i << " " << statusString(model.getColumnStatus(i)) << "\n";
	}
	for (int i = 0; i < ncons1; i++) {
		f << i+nvar1real <<  " " << statusString(model.getRowStatus(i)) << "\n";
	}
	f.close();
	assert(nbasic == totalRows);

}

static ClpSimplex::Status statusFromString(const string& s) {
	if (s == "Basic") {
		return ClpSimplex::basic;
	} else if (s == "AtLower") {
		return ClpSimplex::atLowerBound;
	} else if (s == "AtUpper") {
		return ClpSimplex::atUpperBound;
	} else {
		assert(0); return ClpSimplex::atLowerBound;
	}
}

void CbcBALPInterface::loadStatus(const std::string &filebase) {

	int nscen = dims.numScenarios();
	int nvar1real = dims.numFirstStageVars();
	int ncons1 = dims.numFirstStageCons();
	
	int roffset = ncons1;
	int coffset = nvar1real;
	string line;

	for (int k = 0; k < nscen; k++) {
		stringstream fname;
		fname << filebase << k+1;
		ifstream f(fname.str().c_str());
		f.exceptions(ifstream::failbit | ifstream::badbit);
		int nbasicThis, r;
		getline(f,line);
		assert(line.find("BasisOnly") != string::npos);
		istringstream iss(line);
		iss >> nbasicThis;
		string status;
		int nvar2real = dims.numSecondStageVars(k);
		int ncons2 = dims.numSecondStageCons(k);
		for (int i = 0; i < nvar2real; i++) {
			f >> r;
			f >> status;
			assert(r == i);
			model.setColumnStatus(r+coffset, statusFromString(status));

		}
		for (int i = 0; i < ncons2; i++) {
			f >> r;
			f >> status;
			assert(r == i + nvar2real);
			model.setRowStatus(i+roffset, statusFromString(status));
		}

		f.close();
		roffset += ncons2;
		coffset += nvar2real;
	}
	stringstream fname;
	fname << filebase << 0;
	ifstream f(fname.str().c_str());
	f.exceptions(ifstream::failbit | ifstream::badbit);

	int nbasicThis, r;
	getline(f,line);
	assert(line.find("BasisOnly") != string::npos);
	istringstream iss(line);
	iss >> nbasicThis;

	string status;
	for (int i = 0; i < nvar1real; i++) {
		f >> r;
		f >> status;
		assert(r == i);
		model.setColumnStatus(i, statusFromString(status));

	}
	for (int i = 0; i < ncons1; i++) {
		f >> r;
		f >> status;
		assert(r == i + nvar1real);
		model.setRowStatus(i, statusFromString(status));
	}

}


void CbcBALPInterface::addRow(const std::vector<double>& elts1, const std::vector<double> &elts2, int scen, double lb, double ub) {
	
	vector<double> elts;
	vector<int> idx;
	
	int nvar1 = dims.numFirstStageVars();
	assert(elts1.size() == static_cast<unsigned>(nvar1));

	for (int i = 0; i < nvar1; i++) {
		if (elts1[i]) {
			elts.push_back(elts1[i]);
			idx.push_back(i);
		}
	}
	int nvar2 = dims.numSecondStageVars(scen);
	assert(elts2.size() == static_cast<unsigned>(nvar2));

	int offset = nvar1;
	for (int i = 0; i < scen; i++) offset += dims.numSecondStageVars(i);

	for (int i = 0; i < nvar2; i++) {
		if (elts2[i]) {
			elts.push_back(elts2[i]);
			idx.push_back(offset+i);
		}
	}

	model.addRow(elts.size(),&idx[0],&elts[0],lb,ub);

}



void CbcBALPInterface::setFirstStageColLB(int idx, double newLb) {
	model.setColLower(idx, newLb);
}

void CbcBALPInterface::setFirstStageColUB(int idx, double newUb) {
	model.setColUpper(idx, newUb);
}
*/
