#include "cuttingPlaneBALP.hpp"
#include <cmath>

using namespace std;


cuttingPlaneModel::cuttingPlaneModel(int nvar1, std::vector<std::vector<cutInfo> > const &cuts, double LB) :
	nvar1(nvar1), cuts(cuts), LB(LB) {}

vector<double> cuttingPlaneModel::getFirstStageColLB() {
	// v >= 0, l free
	vector<double> l(nvar1 + 1);
	l[0] = 0.;
	std::fill(l.begin()+1,l.end(),-COIN_DBL_MAX);
	//std::fill(l.begin()+1,l.end(),0.);

	return l;
}

vector<double> cuttingPlaneModel::getFirstStageColUB() {

	return vector<double>(nvar1+1,COIN_DBL_MAX);
}

vector<double> cuttingPlaneModel::getFirstStageObj() {

	vector<double> obj(nvar1+1,0.);
	obj[0] = -LB;
	return obj;
}
vector<string> cuttingPlaneModel::getFirstStageColNames() {

	vector<string> names(nvar1+1);
	names[0] = "v";
	std::fill(names.begin()+1,names.end(),"l"); // don't worry about numbers for now
	return names;
}

vector<double> cuttingPlaneModel::getFirstStageRowLB() {
	return vector<double>(0);
}
vector<double> cuttingPlaneModel::getFirstStageRowUB() {
	return vector<double>(0);
}

vector<string> cuttingPlaneModel::getFirstStageRowNames() {
	return vector<string>(0);
}


vector<double> cuttingPlaneModel::getSecondStageColLB(int scen) {
	return vector<double>(cuts[scen].size(),0.);
}

vector<double> cuttingPlaneModel::getSecondStageColUB(int scen) {
	return vector<double>(cuts[scen].size(),COIN_DBL_MAX);
}

// f_i(\gamma_i^r) - g^T\gamma_i^r
double cutInfo::computeC() const {
	assert(evaluatedAt.size() == subgradient.size());
	double out = objval;
	for (unsigned i = 0; i < subgradient.size(); i++) {
		out -= subgradient[i]*evaluatedAt[i];
	}
	return out;
}

vector<double> cuttingPlaneModel::getSecondStageObj(int scen) {
	vector<double> obj;
	obj.reserve(cuts[scen].size());
	for (unsigned i = 0; i < cuts[scen].size(); i++) {
		obj.push_back(-cuts[scen][i].computeC());
	}
	return obj;
}

vector<string> cuttingPlaneModel::getSecondStageColNames(int scen) {
	vector<string> names(cuts[scen].size(),"u");
	return names;
}

vector<double> cuttingPlaneModel::getSecondStageRowUB(int scen) {
	vector<double> ub(nvar1+1,0.);
	ub[0] = 1.;
	return ub;
}

vector<double> cuttingPlaneModel::getSecondStageRowLB(int scen) {
	return getSecondStageRowUB(scen);
}

vector<string> cuttingPlaneModel::getSecondStageRowNames(int scen) {
	vector<string> names(nvar1+1);
	names[0] = "t";
	std::fill(names.begin()+1,names.end(),"g");
	return names;

}

CoinPackedMatrix cuttingPlaneModel::getFirstStageConstraints() {
	vector<CoinBigIndex> starts(nvar1+2,0);
	return CoinPackedMatrix(true, 0, nvar1+1,0,0,0,&starts[0],0);
}

CoinPackedMatrix cuttingPlaneModel::getSecondStageConstraints(int scen) {
	/* diagonal blocks are of the form:
	   [  e^T ]
	   [-G_i^T].
	   So each column starts with 1 and then has a (negative) subgradient. */
	
	CoinBigIndex nnz = 0;
	vector<double> elts;
	vector<int> idx;
	unsigned ncol = cuts[scen].size();
	vector<CoinBigIndex> starts(ncol+1);
	for (unsigned col = 0; col < ncol; col++) {
		elts.push_back(1.);
		idx.push_back(0);
		starts[col] = nnz++;
		vector<double> const& subgrad = cuts[scen][col].subgradient;
		for (unsigned k = 0; k < subgrad.size(); k++) {
			if (fabs(subgrad[k]) > 1e-7) {
				elts.push_back(-subgrad[k]);
				idx.push_back(k+1);
				nnz++;
			}
		}
	}
	starts[ncol] = nnz;

	return CoinPackedMatrix(true,nvar1+1,ncol,nnz,&elts[0],&idx[0],&starts[0],0);

}

CoinPackedMatrix cuttingPlaneModel::getLinkingConstraints(int scen) {
	// this is easy, (nvar+1)x(nvar+1) identity matrix!
	
	vector<double> elts(nvar1+1,1.);
	vector<int> idx;
	vector<CoinBigIndex> starts;
	for (int i = 0; i < nvar1+1; i++) {
		starts.push_back(i);
		idx.push_back(i);
	}
	starts.push_back(nvar1+1);

	return CoinPackedMatrix(true,nvar1+1,nvar1+1,nvar1+1,&elts[0],&idx[0],&starts[0],0);

}


