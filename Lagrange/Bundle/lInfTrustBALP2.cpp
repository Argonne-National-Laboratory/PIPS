#include "lInfTrustBALP2.hpp"
#include <cmath>

using namespace std;


lInfTrustModel2::lInfTrustModel2(int nvar1, bundle_t const &cuts, double trustRadius, std::vector<std::vector<double> > const& center) :
	nvar1(nvar1), cuts(cuts), trustRadius(trustRadius), center(center) {}

vector<double> lInfTrustModel2::getFirstStageColLB() {
	vector<double> l(nvar1,-COIN_DBL_MAX);

	return l;
}

vector<double> lInfTrustModel2::getFirstStageColUB() {

	return vector<double>(nvar1,COIN_DBL_MAX);
}

vector<double> lInfTrustModel2::getFirstStageObj() {

	vector<double> obj(nvar1,0.);
	return obj;
}
vector<string> lInfTrustModel2::getFirstStageColNames() {

	vector<string> names(nvar1,"l");
	return names;
}

vector<double> lInfTrustModel2::getFirstStageRowLB() {
	return vector<double>(0);
}
vector<double> lInfTrustModel2::getFirstStageRowUB() {
	return vector<double>(0);
}

vector<string> lInfTrustModel2::getFirstStageRowNames() {
	return vector<string>(0);
}


vector<double> lInfTrustModel2::getSecondStageColLB(int scen) {
	vector<double> lb(cuts[scen].size()+3*nvar1,0.);
	return lb;
}

vector<double> lInfTrustModel2::getSecondStageColUB(int scen) {
	vector<double> ub(cuts[scen].size()+3*nvar1,COIN_DBL_MAX);
	return ub;
}


vector<double> lInfTrustModel2::getSecondStageObj(int scen) {
	vector<double> obj;
	obj.reserve(cuts[scen].size()+3*nvar1);
	for (int i = 0; i < nvar1; i++) {
		obj.push_back(-center[scen][i]);
	}
	for (int i = 0; i < nvar1; i++) {
		obj.push_back(center[scen][i]);
	}
	for (int i = 0; i < nvar1; i++) {
		obj.push_back(trustRadius);
	}
	for (unsigned i = 0; i < cuts[scen].size(); i++) {
		obj.push_back(-cuts[scen][i].computeC());
	}
	
	return obj;
}

vector<string> lInfTrustModel2::getSecondStageColNames(int scen) {
	vector<string> names(cuts[scen].size()+3*nvar1,"");
	return names;
}

vector<double> lInfTrustModel2::getSecondStageRowUB(int scen) {
	vector<double> ub(2*nvar1+1,0.);
	ub[0] = 1.;
	return ub;
}

vector<double> lInfTrustModel2::getSecondStageRowLB(int scen) {
	return getSecondStageRowUB(scen);
}

vector<string> lInfTrustModel2::getSecondStageRowNames(int scen) {
	vector<string> names(2*nvar1+1);
	return names;

}

CoinPackedMatrix lInfTrustModel2::getFirstStageConstraints() {
	vector<CoinBigIndex> starts(nvar1+1,0);
	return CoinPackedMatrix(true, 0, nvar1,0,0,0,&starts[0],0);
}

CoinPackedMatrix lInfTrustModel2::getSecondStageConstraints(int scen) {
	/* diagonal blocks are of the form:
	   [            e^T   ] 
	   [ I  -I     -G_i^T ]
	   [ I   I  -I        ]
	   This order is convenient for updating the basis when subgradients are added
	 */
	
	CoinBigIndex nnz = 0;
	unsigned ncol = cuts[scen].size() + 3*nvar1;
	vector<double> elts; elts.reserve(cuts[scen].size()*(nvar1+1)+3*nvar1);
	vector<int> idx; idx.reserve(cuts[scen].size()*(nvar1+1)+3*nvar1);
	vector<CoinBigIndex> starts(ncol+1);
	for (int i = 0; i < nvar1; i++) {
		elts.push_back(1.);
		elts.push_back(1.);
		idx.push_back(i+1);
		idx.push_back(i+1+nvar1);
		starts[i] = nnz;
		nnz += 2;
	}
	for (int i = 0; i < nvar1; i++) {
		elts.push_back(-1.);
		elts.push_back(1.);
		idx.push_back(i+1);
		idx.push_back(i+1+nvar1);
		starts[i+nvar1] = nnz;
		nnz += 2;
	}
	for (int i = 0; i < nvar1; i++) {
		elts.push_back(-1.);
		idx.push_back(i+1+nvar1);
		starts[i+2*nvar1] = nnz++;
	}
	for (unsigned col = 0; col < cuts[scen].size(); col++) {
		elts.push_back(1.);
		idx.push_back(0);
		starts[col+3*nvar1] = nnz++;
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

	return CoinPackedMatrix(true,nSecondStageCons(scen),ncol,nnz,&elts[0],&idx[0],&starts[0],0);

}

CoinPackedMatrix lInfTrustModel2::getLinkingConstraints(int scen) {
	// [   ]
	// [ I ]
	// [   ]

	vector<double> elts(nvar1,1.);
	vector<int> idx; idx.reserve(nvar1);
	vector<CoinBigIndex> starts; starts.reserve(nvar1+1);
	for (int i = 0; i < nvar1; i++) {
		starts.push_back(i);
		idx.push_back(i+1);
	}
	starts.push_back(nvar1);

	return CoinPackedMatrix(true,nSecondStageCons(scen),nvar1,nvar1,&elts[0],&idx[0],&starts[0],0);

}


